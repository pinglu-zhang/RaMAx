#include "external_tool.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <limits.h>
#include <signal.h>
#include <spawn.h>
#include <sstream>
#include <stdexcept>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>

extern char** environ;

#ifndef RAMAX_TOOL_BIN_CONFIGURED_PATH
#define RAMAX_TOOL_BIN_CONFIGURED_PATH ""
#endif

namespace RaMAxExternalTool {

namespace {

using Clock = std::chrono::steady_clock;

CommandResult decodeStatus(
    int status, std::chrono::milliseconds elapsed, bool timed_out) {
    CommandResult result;
    result.timed_out = timed_out;
    result.elapsed = elapsed;
    if (WIFSIGNALED(status)) {
        result.termination_signal = WTERMSIG(status);
    }
    if (timed_out) {
        // Keep timeout distinct from a command that naturally exits 124.
        result.exit_code = 124;
    } else if (WIFEXITED(status)) {
        result.exit_code = WEXITSTATUS(status);
    } else if (WIFSIGNALED(status)) {
        result.exit_code = 128 + WTERMSIG(status);
    } else {
        result.exit_code = 255;
    }
    return result;
}

pid_t waitWithoutBlocking(pid_t pid, int& status) {
    while (true) {
        const pid_t waited = ::waitpid(pid, &status, WNOHANG);
        if (waited >= 0) return waited;
        if (errno == EINTR) continue;
        throw std::runtime_error(
            "Cannot wait for child process: " + std::string(std::strerror(errno)));
    }
}

void waitBlocking(pid_t pid, int& status, const std::filesystem::path& executable) {
    while (::waitpid(pid, &status, 0) < 0) {
        if (errno == EINTR) continue;
        throw std::runtime_error("Cannot wait for " + executable.string() +
                                 ": " + std::strerror(errno));
    }
}

void sendSignal(pid_t pid, int signal_number, bool process_group) {
    if (process_group) {
        if (::kill(-pid, signal_number) != 0 && errno != ESRCH) {
            throw std::runtime_error(
                "Cannot signal child process group: " +
                std::string(std::strerror(errno)));
        }
        return;
    }
    if (::kill(pid, signal_number) != 0 && errno != ESRCH) {
        throw std::runtime_error(
            "Cannot signal child process: " + std::string(std::strerror(errno)));
    }
}

bool processGroupExists(pid_t pid) {
    if (::kill(-pid, 0) == 0) return true;
    return errno == EPERM;
}

void sleepUntilNextPoll(
    Clock::time_point deadline, std::chrono::milliseconds poll_interval) {
    const auto now = Clock::now();
    if (now >= deadline) return;
    const auto remaining = std::chrono::duration_cast<std::chrono::milliseconds>(
        deadline - now);
    const auto normalized_poll = std::max(
        std::chrono::milliseconds(1), poll_interval);
    std::this_thread::sleep_for(std::min(normalized_poll, remaining));
}

}  // namespace

bool isExecutable(const std::filesystem::path& candidate) {
    std::error_code error;
    return !candidate.empty() &&
           std::filesystem::is_regular_file(candidate, error) && !error &&
           ::access(candidate.c_str(), X_OK) == 0;
}

std::filesystem::path executableDirectory() {
    std::array<char, PATH_MAX + 1> buffer{};
    const ssize_t length = ::readlink("/proc/self/exe", buffer.data(), PATH_MAX);
    if (length <= 0) return {};
    buffer[static_cast<size_t>(length)] = '\0';
    return std::filesystem::path(buffer.data()).parent_path();
}

std::filesystem::path searchPath(std::string_view name) {
    const char* raw_path = std::getenv("PATH");
    if (!raw_path) return {};
    std::string_view path(raw_path);
    size_t offset = 0;
    while (offset <= path.size()) {
        const size_t separator = path.find(':', offset);
        const size_t end = separator == std::string_view::npos
            ? path.size() : separator;
        if (end > offset) {
            const auto candidate =
                std::filesystem::path(path.substr(offset, end - offset)) / name;
            if (isExecutable(candidate)) {
                return std::filesystem::absolute(candidate);
            }
        }
        if (separator == std::string_view::npos) break;
        offset = separator + 1;
    }
    return {};
}

std::filesystem::path locateExecutable(
    std::string_view name,
    const std::filesystem::path& configured_path) {
    if (isExecutable(configured_path)) {
        return std::filesystem::absolute(configured_path);
    }
    const std::filesystem::path configured_directory(
        RAMAX_TOOL_BIN_CONFIGURED_PATH);
    if (!configured_directory.empty()) {
        const auto candidate = configured_directory / name;
        if (isExecutable(candidate)) {
            return std::filesystem::absolute(candidate);
        }
    }
    const auto sibling = executableDirectory() / name;
    if (isExecutable(sibling)) return std::filesystem::absolute(sibling);
    return searchPath(name);
}

CommandResult run(
    const std::filesystem::path& executable,
    const std::vector<std::string>& arguments,
    const std::filesystem::path& stdout_path,
    const std::filesystem::path& stderr_path) {
    return run(executable, arguments, stdout_path, stderr_path, RunOptions{});
}

CommandResult run(
    const std::filesystem::path& executable,
    const std::vector<std::string>& arguments,
    const std::filesystem::path& stdout_path,
    const std::filesystem::path& stderr_path,
    const RunOptions& options) {
    const auto start = Clock::now();
    const int stdout_fd = ::open(stdout_path.c_str(),
        O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (stdout_fd < 0) {
        throw std::runtime_error("Cannot open command stdout file " +
                                 stdout_path.string() + ": " +
                                 std::strerror(errno));
    }
    const int stderr_fd = ::open(stderr_path.c_str(),
        O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (stderr_fd < 0) {
        ::close(stdout_fd);
        throw std::runtime_error("Cannot open command stderr file " +
                                 stderr_path.string() + ": " +
                                 std::strerror(errno));
    }

    posix_spawn_file_actions_t actions;
    posix_spawn_file_actions_init(&actions);
    posix_spawn_file_actions_adddup2(&actions, stdout_fd, STDOUT_FILENO);
    posix_spawn_file_actions_adddup2(&actions, stderr_fd, STDERR_FILENO);
    posix_spawn_file_actions_addclose(&actions, stdout_fd);
    posix_spawn_file_actions_addclose(&actions, stderr_fd);

    std::vector<std::string> storage;
    storage.reserve(arguments.size() + 1);
    storage.push_back(executable.string());
    storage.insert(storage.end(), arguments.begin(), arguments.end());
    std::vector<char*> argv;
    argv.reserve(storage.size() + 1);
    for (auto& value : storage) argv.push_back(value.data());
    argv.push_back(nullptr);

    posix_spawnattr_t attributes;
    bool attributes_initialized = false;
    posix_spawnattr_t* attributes_pointer = nullptr;
    if (options.create_process_group) {
        int attribute_error = ::posix_spawnattr_init(&attributes);
        attributes_initialized = attribute_error == 0;
        if (attribute_error == 0) {
            attribute_error = ::posix_spawnattr_setpgroup(&attributes, 0);
        }
        if (attribute_error == 0) {
            attribute_error = ::posix_spawnattr_setflags(
                &attributes, POSIX_SPAWN_SETPGROUP);
        }
        if (attribute_error != 0) {
            if (attributes_initialized) {
                ::posix_spawnattr_destroy(&attributes);
            }
            ::posix_spawn_file_actions_destroy(&actions);
            ::close(stdout_fd);
            ::close(stderr_fd);
            throw std::runtime_error(
                "Cannot configure child process group for " +
                executable.string() + ": " + std::strerror(attribute_error));
        }
        attributes_pointer = &attributes;
    }

    pid_t pid = -1;
    const int spawn_error = ::posix_spawn(
        &pid, executable.c_str(), &actions, attributes_pointer,
        argv.data(), environ);
    if (attributes_pointer != nullptr) {
        ::posix_spawnattr_destroy(&attributes);
    }
    posix_spawn_file_actions_destroy(&actions);
    ::close(stdout_fd);
    ::close(stderr_fd);
    if (spawn_error != 0) {
        throw std::runtime_error("Cannot start " + executable.string() +
                                 ": " + std::strerror(spawn_error));
    }

    int status = 0;
    if (options.timeout <= std::chrono::milliseconds::zero()) {
        waitBlocking(pid, status, executable);
        return decodeStatus(
            status,
            std::chrono::duration_cast<std::chrono::milliseconds>(
                Clock::now() - start),
            false);
    }

    const auto deadline = start + options.timeout;
    while (Clock::now() < deadline) {
        const pid_t waited = waitWithoutBlocking(pid, status);
        if (waited == pid) {
            return decodeStatus(
                status,
                std::chrono::duration_cast<std::chrono::milliseconds>(
                    Clock::now() - start),
                false);
        }
        sleepUntilNextPoll(deadline, options.poll_interval);
    }

    // Close the race where the child exits between the final poll and timeout.
    if (waitWithoutBlocking(pid, status) == pid) {
        return decodeStatus(
            status,
            std::chrono::duration_cast<std::chrono::milliseconds>(
                Clock::now() - start),
            false);
    }

    sendSignal(pid, SIGTERM, options.create_process_group);
    bool child_reaped = false;
    const auto grace_deadline = Clock::now() +
        std::max(std::chrono::milliseconds::zero(), options.termination_grace);
    while (Clock::now() < grace_deadline) {
        if (!child_reaped && waitWithoutBlocking(pid, status) == pid) {
            child_reaped = true;
        }
        const bool process_tree_alive = options.create_process_group
            ? processGroupExists(pid) : !child_reaped;
        if (child_reaped && !process_tree_alive) {
            return decodeStatus(
                status,
                std::chrono::duration_cast<std::chrono::milliseconds>(
                    Clock::now() - start),
                true);
        }
        sleepUntilNextPoll(grace_deadline, options.poll_interval);
    }

    if (!child_reaped && waitWithoutBlocking(pid, status) == pid) {
        child_reaped = true;
    }
    const bool process_tree_alive = options.create_process_group
        ? processGroupExists(pid) : !child_reaped;
    if (process_tree_alive) {
        sendSignal(pid, SIGKILL, options.create_process_group);
    }
    if (!child_reaped) {
        waitBlocking(pid, status, executable);
        child_reaped = true;
    }

    if (options.create_process_group && process_tree_alive) {
        const auto group_cleanup_deadline = Clock::now() + std::chrono::seconds(2);
        while (processGroupExists(pid) && Clock::now() < group_cleanup_deadline) {
            sleepUntilNextPoll(group_cleanup_deadline, options.poll_interval);
        }
    }
    return decodeStatus(
        status,
        std::chrono::duration_cast<std::chrono::milliseconds>(
            Clock::now() - start),
        true);
}

std::string readText(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("Cannot read file: " + path.string());
    }
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

}  // namespace RaMAxExternalTool
