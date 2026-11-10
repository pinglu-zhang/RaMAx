// Small launcher regression. The child echoes FASTA; this is not an MSA test.
#include "external_msa_runner.h"
#include "spawn_compat.h"
#include <atomic>
#include <cassert>
#include <csignal>
#include <fstream>
#include <spawn.h>
#include <thread>
#include <unistd.h>
#include <sys/syscall.h>
#include <sys/vfs.h>
#include <sys/wait.h>
#include <cstdarg>
#include <fcntl.h>

namespace {
std::atomic<int> fail_init{0}, fail_flags{0}, fail_action{0}, live_attributes{0};
std::atomic<int> launches{0}, path_launches{0};
bool expect_legacy = false;
std::string backend;
void checkAttributes(const posix_spawnattr_t* attributes) {
    if (!expect_legacy) { assert(attributes == nullptr); return; }
    assert(attributes != nullptr);
    short flags = 0;
    assert(posix_spawnattr_getflags(attributes, &flags) == 0);
    assert(flags == POSIX_SPAWN_USEVFORK);
}
size_t fdCount() {
    return std::distance(std::filesystem::directory_iterator("/proc/self/fd"),
                         std::filesystem::directory_iterator{});
}
}

extern "C" {
long __real_syscall(long, ...);
ssize_t __real_write(int, const void*, size_t);
int __real_mkostemp(char*, int);
int __real_posix_spawn_file_actions_adddup2(posix_spawn_file_actions_t*, int, int);
ssize_t __wrap_write(int fd, const void* data, size_t size) {
    struct statfs fs{};
    if (backend == "fullshm" && fstatfs(fd, &fs) == 0 && fs.f_type == 0x01021994) {
        errno = ENOSPC; return -1;
    }
    return __real_write(fd, data, size);
}
long __wrap_syscall(long number, ...) {
    va_list args; va_start(args, number);
    const char* name = va_arg(args, const char*);
    unsigned flags = va_arg(args, unsigned);
    va_end(args);
    if (backend != "memfd") { errno = ENOSYS; return -1; }
    return __real_syscall(number, name, flags);
}
int __wrap_mkostemp(char* name, int flags) {
    if (backend == "file" && std::string(name).starts_with("/dev/shm/")) {
        errno = EACCES; return -1;
    }
    return __real_mkostemp(name, flags);
}
int __wrap_posix_spawn_file_actions_adddup2(posix_spawn_file_actions_t* a, int f, int t) {
    if (fail_action.exchange(0)) return EBADF;
    return __real_posix_spawn_file_actions_adddup2(a, f, t);
}
int __real_posix_spawnattr_init(posix_spawnattr_t*);
int __real_posix_spawnattr_setflags(posix_spawnattr_t*, short);
int __real_posix_spawnattr_destroy(posix_spawnattr_t*);
int __real_posix_spawn(pid_t*, const char*, const posix_spawn_file_actions_t*,
                       const posix_spawnattr_t*, char* const[], char* const[]);
int __real_posix_spawnp(pid_t*, const char*, const posix_spawn_file_actions_t*,
                        const posix_spawnattr_t*, char* const[], char* const[]);
int __wrap_posix_spawnattr_init(posix_spawnattr_t* a) {
    if (fail_init.exchange(0)) return ENOMEM;
    int result = __real_posix_spawnattr_init(a);
    if (!result) ++live_attributes;
    return result;
}
int __wrap_posix_spawnattr_setflags(posix_spawnattr_t* a, short flags) {
    if (fail_flags.exchange(0)) return EINVAL;
    return __real_posix_spawnattr_setflags(a, flags);
}
int __wrap_posix_spawnattr_destroy(posix_spawnattr_t* a) {
    --live_attributes;
    return __real_posix_spawnattr_destroy(a);
}
int __wrap_posix_spawn(pid_t* p, const char* e,
                      const posix_spawn_file_actions_t* a,
                      const posix_spawnattr_t* s, char* const v[], char* const n[]) {
    checkAttributes(s); ++launches;
    return __real_posix_spawn(p, e, a, s, v, n);
}
int __wrap_posix_spawnp(pid_t* p, const char* e,
                       const posix_spawn_file_actions_t* a,
                       const posix_spawnattr_t* s, char* const v[], char* const n[]) {
    checkAttributes(s); ++path_launches;
    return __real_posix_spawnp(p, e, a, s, v, n);
}
}

int main(int argc, char** argv) {
    if (argc == 6 && std::string(argv[2]) == "-S") {
        assert(std::filesystem::current_path() == std::filesystem::path(getenv("TEST_EXPECT_CWD")));
        std::ofstream("minipoa_paths_" + std::to_string(getpid()) + ".tmp") << "fixture";
        assert(std::string(argv[3]) == "-f" && std::string(argv[4]) == "0");
        assert(std::string(argv[5]) == "-r1");
        const auto name = std::filesystem::path(argv[0]).filename().string();
        if (name == "exit7") return 7;
        if (name == "signal") { raise(SIGTERM); return 9; }
        if (name == "invalid") { std::cout << ">invalid\nAC\n"; return 0; }
        std::ifstream input(argv[1], std::ios::binary);
        if (!input) return 8;
        std::cout << input.rdbuf();
        return 0;
    }
    assert(argc == 4);
#ifdef RAMAX_LEGACY_SPAWN
    expect_legacy = true;
#endif
    backend = argv[3];
    const auto directory = std::filesystem::absolute(argv[1]);
    std::filesystem::create_directories(directory);
    const auto self = std::filesystem::canonical("/proc/self/exe");
    for (const auto* name : {"echo_msa", "exit7", "signal", "invalid"})
        std::filesystem::create_symlink(self, directory / name);
    const auto path = directory.string() + ":" + (getenv("PATH") ? getenv("PATH") : "");
    assert(setenv("PATH", path.c_str(), 1) == 0);
    auto& runner = RaMesh::Alignment::ExternalMsaRunner::instance();
    const auto scratch = directory / "local";
    assert(setenv("RAMAX_MINIPOA_TMPDIR", scratch.c_str(), 1) == 0);
    assert(setenv("TEST_EXPECT_CWD", scratch.c_str(), 1) == 0);
    runner.configureScratchDirectory(directory / "unused");
    assert(!std::filesystem::exists(directory / "unused"));
    const size_t original_fds = fdCount();
    using Rows = std::unordered_map<ChrName, std::string>;
    const Rows rows{{"z", "ACGTNacgtn"}, {"a", "ACGTTacgtn"}};
    auto aligned = rows;
    assert(runner.align((directory / "echo_msa").string(), aligned));
    assert(aligned == rows);
    aligned = rows;
    assert(runner.align("echo_msa", aligned) && aligned == rows);
#ifdef RAMAX_SPAWN_HAS_CHDIR
    assert(path_launches == 1);
#else
    assert(path_launches == 0);
#endif
    const int before_cache = launches.load();
    assert(runner.align((directory / "echo_msa").string(), aligned));
    assert(launches == before_cache);
    for (const auto* name : {"exit7", "signal", "invalid", "missing"}) {
        aligned = rows;
        assert(!runner.align((directory / name).string(), aligned));
        assert(aligned == rows);
    }
    if (expect_legacy) {
        for (auto* failure : {&fail_init, &fail_flags}) {
            *failure = 1;
            aligned = {{"a", "AT"}, {"z", "AC"}};
            bool threw = false;
            try { runner.align((directory / "echo_msa").string(), aligned); }
            catch (const std::system_error&) { threw = true; }
            assert(threw && live_attributes == 0 && fdCount() == original_fds);
        }
    }
    fail_action = 1;
    bool action_threw = false;
    aligned = {{"a", "AG"}, {"z", "AC"}};
    try { runner.align((directory / "echo_msa").string(), aligned); }
    catch (const std::system_error&) { action_threw = true; }
    assert(action_threw && fdCount() == original_fds);
#ifndef RAMAX_SPAWN_HAS_CHDIR
    assert(setenv("RAMAX_SPAWN_HELPER", "/missing/ramax-spawn-helper", 1) == 0);
    bool helper_threw = false;
    try { runner.align((directory / "echo_msa").string(), aligned); }
    catch (const std::runtime_error&) { helper_threw = true; }
    assert(helper_threw && fdCount() == original_fds);
    unsetenv("RAMAX_SPAWN_HELPER");
#endif
    std::vector<std::thread> threads;
    for (int i = 0; i < 8; ++i) threads.emplace_back([&, i] {
        for (int j = 0; j < 4; ++j) {
            // Exceeds a pipe buffer, with unique input to exercise real launches.
            std::string sequence(70000 + i * 4 + j, 'A');
            Rows input{{"b", sequence}, {"a", sequence}};
            const Rows expected = input;
            assert(runner.align((directory / "echo_msa").string(), input));
            assert(input == expected);
        }
    });
    for (auto& t : threads) t.join();
    pid_t child = fork();
    assert(child >= 0);
    if (!child) {
        close(STDIN_FILENO); close(STDOUT_FILENO); close(STDERR_FILENO);
        Rows closed{{"a", "AACT"}, {"b", "AATT"}};
        const auto expected = closed;
        assert(runner.align((directory / "echo_msa").string(), closed) && closed == expected);
        _exit(0);
    }
    int status;
    assert(waitpid(child, &status, 0) == child && WIFEXITED(status) && WEXITSTATUS(status) == 0);
    if (const char* real = getenv("RAMAX_TEST_MINIPOA")) {
        std::string reference;
        for (int i = 0; i < 16; ++i) reference += "ACGTTGCAAGTCCGAT";
        Rows real_rows{{"a", reference}, {"b", reference.substr(0, 91) + reference.substr(94)},
                       {"c", reference.substr(0, 77) + "GGA" + reference.substr(77)}};
        const auto expected = real_rows;
        assert(runner.align(real, real_rows));
        std::vector<std::string> ordered;
        for (const auto& key : {"a", "b", "c"}) {
            auto ungapped = real_rows.at(key);
            std::erase(ungapped, '-');
            assert(ungapped == expected.at(key));
            ordered.push_back(real_rows.at(key));
        }
        assert(ordered[0].size() == ordered[1].size() && ordered[0].size() == ordered[2].size());
        std::cout << "REAL_OUTPUT " << ordered[0] << '|' << ordered[1] << '|' << ordered[2] << '\n';
    }
    assert(live_attributes == 0 && fdCount() == original_fds);
    assert(std::filesystem::is_empty(scratch));
    for (const auto& entry : std::filesystem::directory_iterator(std::filesystem::current_path()))
        assert(!entry.path().filename().string().starts_with("minipoa_paths_"));
    runner.logSummary();
    std::cout << "PASS: path/PATH, cache, exit 7, signal, invalid output, "
                 "missing executable, 32 concurrent launches, descriptors, cwd, temp cleanup, input backend=" << backend;
    if (expect_legacy) std::cout << ", attribute failure cleanup (legacy branch)";
    std::cout << '\n';
}
