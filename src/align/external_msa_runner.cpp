#include "external_msa_runner.h"

#include <spdlog/spdlog.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cctype>
#include <condition_variable>
#include <cstdlib>
#include <deque>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <memory>
#include <spawn.h>
#include <string_view>
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
std::atomic<uint64_t> cache_hits{0};
std::atomic<uint64_t> cache_misses{0};
std::atomic<uint64_t> cache_single_flight_hits{0};
std::atomic<uint32_t> progress_index{0};
std::atomic_flag spawn_warning = ATOMIC_FLAG_INIT;
std::atomic_flag command_warning = ATOMIC_FLAG_INIT;
std::atomic_flag parse_warning = ATOMIC_FLAG_INIT;
constexpr std::array<uint64_t, 10> progress_milestones{
    1000, 2000, 4000, 8000, 12000, 16000, 20000, 24000, 28000, 32000};
constexpr size_t kCacheMaximumBytes = 64ULL * 1024ULL * 1024ULL;
constexpr uint64_t kCacheDisableProbeRequests = 4096;

struct CacheEntry {
    size_t hash = 0;
    std::string executable;
    std::string input;
    std::string output;
    size_t bytes = 0;
};

struct InFlight {
    size_t hash = 0;
    std::string executable;
    std::string input;
    bool finished = false;
    bool success = false;
    std::string output;
    std::condition_variable ready;
};

std::mutex cache_mutex;
std::deque<std::shared_ptr<CacheEntry>> msa_cache;
std::unordered_multimap<size_t, std::shared_ptr<CacheEntry>> cache_index;
std::vector<std::shared_ptr<InFlight>> in_flight;
size_t cache_bytes = 0;
uint64_t cache_requests = 0;
uint64_t cache_successful_hits = 0;
std::atomic<bool> cache_enabled{true};

struct ThreadBuffers {
    std::vector<ChrName> keys;
    std::string input;
    std::string output;
};

thread_local ThreadBuffers thread_buffers;

bool cacheEligible(const std::string& executable) {
    const char* disabled = std::getenv("RAMAX_DISABLE_EXTERNAL_MSA_CACHE");
    return executable == "/usr/local/bin/minipoa" &&
           !(disabled && *disabled == '1');
}

size_t cacheHash(const std::string& executable, const std::string& input) {
    size_t hash = std::hash<std::string>{}(executable);
    const size_t input_hash = std::hash<std::string>{}(input);
    hash ^= input_hash + 0x9e3779b97f4a7c15ULL + (hash << 6U) +
            (hash >> 2U);
    // The fixed invocation contract is part of the key even though it is not
    // currently configurable.
    hash ^= std::hash<std::string_view>{}("-r 1 -t 1") +
            (hash << 6U) + (hash >> 2U);
    return hash;
}

bool exactKey(const CacheEntry& entry, size_t hash,
              const std::string& executable, const std::string& input) {
    return entry.hash == hash && entry.executable == executable &&
           entry.input == input;
}

bool exactKey(const InFlight& entry, size_t hash,
              const std::string& executable, const std::string& input) {
    return entry.hash == hash && entry.executable == executable &&
           entry.input == input;
}

struct CacheLookup {
    bool hit = false;
    std::string output;
    std::shared_ptr<InFlight> owner;
};

CacheLookup lookupCache(const std::string& executable,
                        const std::string& input) {
    CacheLookup result;
    if (!cacheEligible(executable)) return result;
    const size_t hash = cacheHash(executable, input);
    std::unique_lock lock(cache_mutex);
    if (!cache_enabled) return result;
    ++cache_requests;
    const auto [begin, end] = cache_index.equal_range(hash);
    for (auto it = begin; it != end; ++it) {
        if (exactKey(*it->second, hash, executable, input)) {
            ++cache_successful_hits;
            cache_hits.fetch_add(1, std::memory_order_relaxed);
            result.hit = true;
            result.output = it->second->output;
            return result;
        }
    }
    for (const auto& flight : in_flight) {
        if (!exactKey(*flight, hash, executable, input)) continue;
        flight->ready.wait(lock, [&] { return flight->finished; });
        if (flight->success) {
            ++cache_successful_hits;
            cache_hits.fetch_add(1, std::memory_order_relaxed);
            cache_single_flight_hits.fetch_add(1, std::memory_order_relaxed);
            result.hit = true;
            result.output = flight->output;
        }
        return result;
    }
    cache_misses.fetch_add(1, std::memory_order_relaxed);
    if (cache_requests >= kCacheDisableProbeRequests &&
        cache_successful_hits == 0) {
        cache_enabled = false;
        msa_cache.clear();
        cache_index.clear();
        cache_bytes = 0;
        return result;
    }
    result.owner = std::make_shared<InFlight>();
    result.owner->hash = hash;
    result.owner->executable = executable;
    result.owner->input = input;
    in_flight.push_back(result.owner);
    return result;
}

void finishCache(const std::shared_ptr<InFlight>& flight, bool success,
                 const std::string& output) {
    if (!flight) return;
    std::lock_guard lock(cache_mutex);
    flight->success = success;
    if (success) flight->output = output;
    flight->finished = true;
    if (success && cache_enabled) {
        auto entry = std::make_shared<CacheEntry>();
        entry->hash = flight->hash;
        entry->executable = flight->executable;
        entry->input = flight->input;
        entry->output = output;
        entry->bytes = entry->executable.size() + entry->input.size() +
                       entry->output.size();
        while (!msa_cache.empty() &&
               cache_bytes + entry->bytes > kCacheMaximumBytes) {
            const auto evicted = msa_cache.front();
            cache_bytes -= evicted->bytes;
            const auto [first, last] = cache_index.equal_range(evicted->hash);
            for (auto it = first; it != last; ++it) {
                if (it->second == evicted) {
                    cache_index.erase(it);
                    break;
                }
            }
            msa_cache.pop_front();
        }
        if (entry->bytes <= kCacheMaximumBytes) {
            cache_bytes += entry->bytes;
            cache_index.emplace(entry->hash, entry);
            msa_cache.push_back(std::move(entry));
        }
    }
    flight->ready.notify_all();
    in_flight.erase(
        std::remove(in_flight.begin(), in_flight.end(), flight),
        in_flight.end());
}

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

void recordCompleted() {
    const uint64_t current = completed.fetch_add(1) + 1;
    uint32_t milestone = progress_index.load(std::memory_order_relaxed);
    while (milestone < progress_milestones.size() &&
           current >= progress_milestones[milestone]) {
        if (progress_index.compare_exchange_weak(
                milestone, milestone + 1, std::memory_order_relaxed)) {
            spdlog::debug("[external-msa] progress completed={}", current);
            break;
        }
    }
}

class CacheFlightGuard {
public:
    explicit CacheFlightGuard(std::shared_ptr<InFlight> flight)
        : flight_(std::move(flight)) {}
    ~CacheFlightGuard() {
        if (flight_) finishCache(flight_, false, {});
    }
    void succeed(const std::string& output) {
        finishCache(flight_, true, output);
        flight_.reset();
    }

private:
    std::shared_ptr<InFlight> flight_;
};

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
    std::vector<std::string> aligned(keys.size());
    std::vector<bool> seen(keys.size(), false);
    size_t current = keys.size();
    size_t offset = 0;
    while (offset < text.size()) {
        size_t end = text.find('\n', offset);
        if (end == std::string::npos) end = text.size();
        std::string_view line(text.data() + offset, end - offset);
        if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
        offset = end == text.size() ? end : end + 1;
        if (!line.empty() && line.front() == '>') {
            line.remove_prefix(1);
            const size_t token_end = line.find_first_of(" \t");
            const std::string_view id = line.substr(0, token_end);
            if (id.size() < 2 || id.front() != 's' ||
                (id.size() > 2 && id[1] == '0')) {
                warnParseOnce("invalid FASTA header");
                return false;
            }
            size_t index = 0;
            for (size_t i = 1; i < id.size(); ++i) {
                const unsigned char c = static_cast<unsigned char>(id[i]);
                if (!std::isdigit(c) ||
                    index > (keys.size() - 1) / 10) {
                    warnParseOnce("invalid FASTA header");
                    return false;
                }
                index = index * 10 + static_cast<size_t>(c - '0');
            }
            if (index >= keys.size() || seen[index]) {
                warnParseOnce("invalid FASTA header");
                return false;
            }
            current = index;
            seen[index] = true;
            continue;
        }
        if (current == keys.size()) {
            if (line.empty()) continue;
            warnParseOnce("sequence before FASTA header");
            return false;
        }
        for (const unsigned char c : line) {
            if (!std::isspace(c)) {
                aligned[current].push_back(static_cast<char>(c));
            }
        }
    }
    if (std::find(seen.begin(), seen.end(), false) != seen.end()) {
        warnParseOnce("output row count mismatch");
        return false;
    }

    size_t aligned_length = 0;
    for (size_t index = 0; index < keys.size(); ++index) {
        if (!sameUngappedSequence(aligned[index], sequences.at(keys[index]))) {
            warnParseOnce("output sequence validation failed");
            return false;
        }
        if (index == 0) aligned_length = aligned[index].size();
        else if (aligned[index].size() != aligned_length) {
            warnParseOnce("output rows have unequal lengths");
            return false;
        }
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

    auto& keys = thread_buffers.keys;
    keys.clear();
    keys.reserve(sequences.size());
    for (const auto& [key, unused] : sequences) {
        (void)unused;
        keys.push_back(key);
    }
    std::sort(keys.begin(), keys.end());
    size_t input_capacity = keys.size() * 8;
    for (const auto& key : keys) input_capacity += sequences.at(key).size();
    auto& input_text = thread_buffers.input;
    input_text.clear();
    input_text.reserve(input_capacity);
    for (size_t index = 0; index < keys.size(); ++index) {
        input_text += ">s";
        input_text += std::to_string(index);
        input_text.push_back('\n');
        input_text += sequences.at(keys[index]);
        input_text.push_back('\n');
    }

    CacheLookup cache_lookup = lookupCache(executable, input_text);
    if (cache_lookup.hit) {
        const auto parse_start = std::chrono::steady_clock::now();
        input_nanoseconds.fetch_add(
            elapsedNanoseconds(input_start, parse_start));
        if (!parseOutput(cache_lookup.output, keys, sequences)) return false;
        const auto parse_finish = std::chrono::steady_clock::now();
        parse_nanoseconds.fetch_add(
            elapsedNanoseconds(parse_start, parse_finish));
        recordCompleted();
        return true;
    }
    CacheFlightGuard cache_guard(std::move(cache_lookup.owner));

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
    const int spawn_error = executable.find('/') != std::string::npos
        ? posix_spawn(&pid, executable.c_str(), &actions, nullptr,
                      arguments.data(), environ)
        : posix_spawnp(&pid, executable.c_str(), &actions, nullptr,
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

    auto& output_text = thread_buffers.output;
    output_text.clear();
    output_text.reserve(std::max<size_t>(input_text.size(), 4096));
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

    cache_guard.succeed(output_text);
    recordCompleted();
    return true;
}

void ExternalMsaRunner::logSummary() const {
    constexpr double kNsPerSecond = 1.0e9;
    spdlog::info(
        "[external-msa] completed={} failures={} input/process/parse="
        "{:.3f}/{:.3f}/{:.3f}s inputs(memfd/file)={}/{} "
        "cache(hit/miss/single_flight)={}/{}/{} enabled={}",
        completed.load(), failures.load(),
        input_nanoseconds.load() / kNsPerSecond,
        process_nanoseconds.load() / kNsPerSecond,
        parse_nanoseconds.load() / kNsPerSecond, memfd_inputs.load(),
        file_inputs.load(), cache_hits.load(), cache_misses.load(),
        cache_single_flight_hits.load(), cache_enabled.load());
}

}  // namespace RaMesh::Alignment
