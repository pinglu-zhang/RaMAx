#include "external_tool.h"

#include <array>
#include <chrono>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <limits.h>
#include <spawn.h>
#include <sstream>
#include <stdexcept>
#include <sys/wait.h>
#include <sys/resource.h>
#include <unistd.h>

extern char** environ;

namespace RaMAxExternalTool {

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
    const auto sibling = executableDirectory() / name;
    if (isExecutable(sibling)) return std::filesystem::absolute(sibling);
    return searchPath(name);
}

CommandResult run(
    const std::filesystem::path& executable,
    const std::vector<std::string>& arguments,
    const std::filesystem::path& stdout_path,
    const std::filesystem::path& stderr_path) {
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

    const auto started = std::chrono::steady_clock::now();
    pid_t pid = -1;
    const int spawn_error = ::posix_spawn(
        &pid, executable.c_str(), &actions, nullptr, argv.data(), environ);
    posix_spawn_file_actions_destroy(&actions);
    ::close(stdout_fd);
    ::close(stderr_fd);
    if (spawn_error != 0) {
        throw std::runtime_error("Cannot start " + executable.string() +
                                 ": " + std::strerror(spawn_error));
    }

    int status = 0;
    struct rusage usage {};
    while (::wait4(pid, &status, 0, &usage) < 0) {
        if (errno == EINTR) continue;
        throw std::runtime_error("Cannot wait for " + executable.string() +
                                 ": " + std::strerror(errno));
    }
    const double wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    const double user_seconds = static_cast<double>(usage.ru_utime.tv_sec) +
        static_cast<double>(usage.ru_utime.tv_usec) / 1'000'000.0;
    const double system_seconds = static_cast<double>(usage.ru_stime.tv_sec) +
        static_cast<double>(usage.ru_stime.tv_usec) / 1'000'000.0;
    const int exit_code = WIFEXITED(status)
        ? WEXITSTATUS(status)
        : (WIFSIGNALED(status) ? 128 + WTERMSIG(status) : 255);
    return {exit_code, wall_seconds, user_seconds, system_seconds,
            usage.ru_maxrss};
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
