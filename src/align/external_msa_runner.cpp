#include "external_msa_runner.h"

#include <spdlog/spdlog.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cctype>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <spawn.h>
#include <sstream>
#include <system_error>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <utility>
#include <vector>

#if defined(__linux__)
#include <linux/memfd.h>
#include <sys/syscall.h>
#endif

extern char** environ;

namespace RaMesh::Alignment {
namespace {

std::mutex configuration_mutex;
std::string default_executable;
std::atomic<uint64_t> file_counter{0};
std::atomic<uint64_t> completed{0};
std::atomic<uint64_t> input_nanoseconds{0};
std::atomic<uint64_t> process_nanoseconds{0};
std::atomic<uint64_t> parse_nanoseconds{0};
std::atomic<uint64_t> memfd_inputs{0};
std::atomic<uint64_t> file_inputs{0};
std::atomic<uint64_t> failures{0};
std::atomic<uint32_t> progress_index{0};
std::atomic_flag spawn_warning = ATOMIC_FLAG_INIT;
std::atomic_flag command_warning = ATOMIC_FLAG_INIT;
std::atomic_flag parse_warning = ATOMIC_FLAG_INIT;
constexpr std::array<uint64_t, 10> progress_milestones{
    1000, 2000, 4000, 8000, 12000, 16000, 20000, 24000, 28000, 32000};

uint64_t elapsedNanoseconds(
    std::chrono::steady_clock::time_point start,
    std::chrono::steady_clock::time_point finish) {
    return static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            finish - start).count());
}

bool sameUngappedSequence(const std::string& aligned,
                          const std::string& original) {
    size_t left = 0;
    size_t right = 0;
    while (true) {
        while (left < aligned.size() &&
               (aligned[left] == '-' ||
                std::isspace(static_cast<unsigned char>(aligned[left])))) {
            ++left;
        }
        while (right < original.size() &&
               (original[right] == '-' ||
                std::isspace(static_cast<unsigned char>(original[right])))) {
            ++right;
        }
        if (left == aligned.size() || right == original.size()) {
            return left == aligned.size() && right == original.size();
        }
        if (std::toupper(static_cast<unsigned char>(aligned[left])) !=
            std::toupper(static_cast<unsigned char>(original[right]))) {
            return false;
        }
        ++left;
        ++right;
    }
}

void warnParseOnce(const char* message) {
    failures.fetch_add(1, std::memory_order_relaxed);
    if (!parse_warning.test_and_set(std::memory_order_relaxed)) {
        spdlog::warn("[external-msa] {}; further parse failures are aggregated",
                     message);
    }
}

class ScopedFd {
public:
    ScopedFd() = default;
    explicit ScopedFd(int fd) : fd_(fd) {}
    ScopedFd(const ScopedFd&) = delete;
    ScopedFd& operator=(const ScopedFd&) = delete;
    ScopedFd(ScopedFd&& other) noexcept : fd_(other.release()) {}
    ScopedFd& operator=(ScopedFd&& other) noexcept {
        if (this != &other) reset(other.release());
        return *this;
    }
    ~ScopedFd() { reset(); }

    int get() const { return fd_; }
    explicit operator bool() const { return fd_ >= 0; }
    int release() { return std::exchange(fd_, -1); }
    void reset(int fd = -1) {
        if (fd_ >= 0) ::close(fd_);
        fd_ = fd;
    }

private:
    int fd_ = -1;
};

class TemporaryInput {
public:
    std::filesystem::path path;
    ~TemporaryInput() {
        std::error_code error;
        if (!path.empty()) std::filesystem::remove(path, error);
    }
};

bool writeAll(int fd, const std::string& contents) {
    size_t offset = 0;
    while (offset < contents.size()) {
        const ssize_t written = ::write(
            fd, contents.data() + offset, contents.size() - offset);
        if (written < 0) {
            if (errno == EINTR) continue;
            return false;
        }
        offset += static_cast<size_t>(written);
    }
    return true;
}

ScopedFd createMemfd(const std::string& contents) {
#if defined(__linux__) && defined(SYS_memfd_create)
    const int fd = static_cast<int>(
        ::syscall(SYS_memfd_create, "ramax-minipoa-input", MFD_CLOEXEC));
    if (fd >= 0) {
        ScopedFd result(fd);
        if (writeAll(fd, contents) && ::lseek(fd, 0, SEEK_SET) == 0) {
            return result;
        }
    }
#else
    (void)contents;
#endif
    return {};
}

bool readAll(int fd, std::string& contents) {
    std::array<char, 8192> buffer{};
    while (true) {
        const ssize_t count = ::read(fd, buffer.data(), buffer.size());
        if (count == 0) return true;
        if (count < 0) {
            if (errno == EINTR) continue;
            return false;
        }
        contents.append(buffer.data(), static_cast<size_t>(count));
    }
}

bool waitForChild(pid_t pid, int& status) {
    while (::waitpid(pid, &status, 0) < 0) {
        if (errno != EINTR) return false;
    }
    return true;
}

bool parseOutput(
    const std::string& text,
    const std::vector<ChrName>& keys,
    std::unordered_map<ChrName, std::string>& sequences) {
    std::istringstream output(text);
    std::unordered_map<std::string, std::string> aligned_by_id;
    std::string current_id;
    std::string line;
    while (std::getline(output, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty() && line.front() == '>') {
            std::istringstream header(line.substr(1));
            header >> current_id;
            if (current_id.empty() ||
                !aligned_by_id.emplace(current_id, std::string{}).second) {
                warnParseOnce("invalid FASTA header");
                return false;
            }
            continue;
        }
        if (current_id.empty()) {
            if (line.empty()) continue;
            warnParseOnce("sequence before FASTA header");
            return false;
        }
        for (const unsigned char c : line) {
            if (!std::isspace(c)) {
                aligned_by_id[current_id].push_back(static_cast<char>(c));
            }
        }
    }
    if (aligned_by_id.size() != keys.size()) {
        warnParseOnce("output row count mismatch");
        return false;
    }

    size_t aligned_length = 0;
    std::vector<std::string> aligned(keys.size());
    for (size_t index = 0; index < keys.size(); ++index) {
        const auto found = aligned_by_id.find("s" + std::to_string(index));
        if (found == aligned_by_id.end() ||
            !sameUngappedSequence(found->second,
                                  sequences.at(keys[index]))) {
            warnParseOnce("output sequence validation failed");
            return false;
        }
        if (index == 0) aligned_length = found->second.size();
        else if (found->second.size() != aligned_length) {
            warnParseOnce("output rows have unequal lengths");
            return false;
        }
        aligned[index] = found->second;
    }
    for (size_t index = 0; index < keys.size(); ++index) {
        sequences[keys[index]] = std::move(aligned[index]);
    }
    return true;
}

}  // namespace

ExternalMsaRunner& ExternalMsaRunner::instance() {
    static ExternalMsaRunner runner;
    return runner;
}

void ExternalMsaRunner::configureDefaultExecutable(std::string executable) {
    std::lock_guard lock(configuration_mutex);
    default_executable = std::move(executable);
}

bool ExternalMsaRunner::alignWithDefault(
    std::unordered_map<ChrName, std::string>& sequences) {
    std::string executable;
    {
        std::lock_guard lock(configuration_mutex);
        executable = default_executable;
    }
    return align(executable, sequences);
}

bool ExternalMsaRunner::align(
    const std::string& executable,
    std::unordered_map<ChrName, std::string>& sequences) {
    if (executable.empty() || sequences.size() < 2) return false;
    const auto input_start = std::chrono::steady_clock::now();

    std::vector<ChrName> keys;
    keys.reserve(sequences.size());
    for (const auto& [key, unused] : sequences) {
        (void)unused;
        keys.push_back(key);
    }
    std::sort(keys.begin(), keys.end());
    size_t input_capacity = keys.size() * 8;
    for (const auto& key : keys) input_capacity += sequences.at(key).size();
    std::string input_text;
    input_text.reserve(input_capacity);
    for (size_t index = 0; index < keys.size(); ++index) {
        input_text += ">s";
        input_text += std::to_string(index);
        input_text.push_back('\n');
        input_text += sequences.at(keys[index]);
        input_text.push_back('\n');
    }

    ScopedFd input_fd = createMemfd(input_text);
    TemporaryInput temporary;
    std::string input_path;
    if (!input_fd) {
        file_inputs.fetch_add(1, std::memory_order_relaxed);
        const uint64_t serial = file_counter.fetch_add(1);
        const uint64_t stamp = static_cast<uint64_t>(
            std::chrono::steady_clock::now().time_since_epoch().count());
        temporary.path = std::filesystem::temp_directory_path() /
            ("ramax-insertion-msa-" + std::to_string(stamp) + "-" +
             std::to_string(serial) + ".input.fa");
        std::ofstream input(temporary.path, std::ios::binary);
        if (!input || !(input << input_text)) {
            throw std::runtime_error(
                "Cannot create external MSA input: " +
                temporary.path.string());
        }
        input_path = temporary.path.string();
    } else {
        memfd_inputs.fetch_add(1, std::memory_order_relaxed);
    }

    int pipe_fds[2] = {-1, -1};
#if defined(__linux__)
    const int pipe_result = ::pipe2(pipe_fds, O_CLOEXEC);
#else
    const int pipe_result = ::pipe(pipe_fds);
#endif
    if (pipe_result != 0) {
        throw std::system_error(
            errno, std::generic_category(), "Cannot create MSA output pipe");
    }
    ScopedFd output_read(pipe_fds[0]);
    ScopedFd output_write(pipe_fds[1]);
#if !defined(__linux__)
    for (const int fd : pipe_fds) {
        const int flags = ::fcntl(fd, F_GETFD);
        if (flags < 0 || ::fcntl(fd, F_SETFD, flags | FD_CLOEXEC) < 0) {
            throw std::system_error(
                errno, std::generic_category(),
                "Cannot mark MSA output pipe close-on-exec");
        }
    }
#endif

    int child_input_fd = 100;
    while (child_input_fd == input_fd.get() ||
           child_input_fd == output_read.get() ||
           child_input_fd == output_write.get()) ++child_input_fd;
    if (input_fd) {
        input_path = "/proc/self/fd/" + std::to_string(child_input_fd);
    }

    posix_spawn_file_actions_t actions;
    const int init_error = posix_spawn_file_actions_init(&actions);
    if (init_error != 0) {
        throw std::system_error(init_error, std::generic_category(),
                                "Cannot initialize MSA spawn actions");
    }
    int action_error = 0;
    if (input_fd) action_error = posix_spawn_file_actions_adddup2(
        &actions, input_fd.get(), child_input_fd);
    if (action_error == 0) action_error = posix_spawn_file_actions_adddup2(
        &actions, output_write.get(), STDOUT_FILENO);
    if (action_error == 0) action_error = posix_spawn_file_actions_addclose(
        &actions, output_read.get());
    if (action_error == 0 && output_write.get() != STDOUT_FILENO) {
        action_error = posix_spawn_file_actions_addclose(
            &actions, output_write.get());
    }
    if (action_error != 0) {
        posix_spawn_file_actions_destroy(&actions);
        throw std::system_error(action_error, std::generic_category(),
                                "Cannot configure MSA spawn actions");
    }

    std::vector<char*> arguments{
        const_cast<char*>(executable.c_str()), const_cast<char*>("-r"),
        const_cast<char*>("1"), const_cast<char*>("-t"),
        const_cast<char*>("1"), input_path.data(), nullptr};
    pid_t pid = -1;
    const auto process_start = std::chrono::steady_clock::now();
    input_nanoseconds.fetch_add(
        elapsedNanoseconds(input_start, process_start));
    const int spawn_error = posix_spawnp(
        &pid, executable.c_str(), &actions, nullptr,
        arguments.data(), environ);
    posix_spawn_file_actions_destroy(&actions);
    if (spawn_error != 0) {
        failures.fetch_add(1, std::memory_order_relaxed);
        if (!spawn_warning.test_and_set(std::memory_order_relaxed)) {
            spdlog::warn(
                "[external-msa] cannot spawn {}: {}; further spawn failures "
                "are aggregated", executable,
                std::generic_category().message(spawn_error));
        }
        return false;
    }
    output_write.reset();

    std::string output_text;
    const bool read_ok = readAll(output_read.get(), output_text);
    output_read.reset();
    int child_status = 0;
    const bool waited = waitForChild(pid, child_status);
    if (!read_ok || !waited || !WIFEXITED(child_status) ||
        WEXITSTATUS(child_status) != 0) {
        failures.fetch_add(1, std::memory_order_relaxed);
        if (!command_warning.test_and_set(std::memory_order_relaxed)) {
            spdlog::warn(
                "[external-msa] command failed: {}; further command failures "
                "are aggregated", executable);
        }
        return false;
    }

    const auto parse_start = std::chrono::steady_clock::now();
    process_nanoseconds.fetch_add(
        elapsedNanoseconds(process_start, parse_start));
    if (!parseOutput(output_text, keys, sequences)) return false;
    const auto parse_finish = std::chrono::steady_clock::now();
    parse_nanoseconds.fetch_add(elapsedNanoseconds(parse_start, parse_finish));

    const uint64_t current = completed.fetch_add(1) + 1;
    uint32_t milestone = progress_index.load(std::memory_order_relaxed);
    while (milestone < progress_milestones.size() &&
           current >= progress_milestones[milestone]) {
        if (progress_index.compare_exchange_weak(
                milestone, milestone + 1, std::memory_order_relaxed)) {
            spdlog::info("[external-msa] progress completed={}", current);
            break;
        }
    }
    return true;
}

void ExternalMsaRunner::logSummary() const {
    constexpr double kNsPerSecond = 1.0e9;
    spdlog::info(
        "[external-msa] completed={} failures={} input/process/parse="
        "{:.3f}/{:.3f}/{:.3f}s inputs(memfd/file)={}/{}",
        completed.load(), failures.load(),
        input_nanoseconds.load() / kNsPerSecond,
        process_nanoseconds.load() / kNsPerSecond,
        parse_nanoseconds.load() / kNsPerSecond, memfd_inputs.load(),
        file_inputs.load());
}

}  // namespace RaMesh::Alignment
