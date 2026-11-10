// Small launcher regression. The child echoes FASTA; this is not an MSA test.
#include "external_msa_runner.h"
#include <atomic>
#include <cassert>
#include <csignal>
#include <fstream>
#include <spawn.h>
#include <thread>
#include <unistd.h>

namespace {
std::atomic<int> fail_init{0}, fail_flags{0}, live_attributes{0};
std::atomic<int> launches{0}, path_launches{0};
bool expect_legacy = false;
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
    assert(argc == 3);
    expect_legacy = std::string(argv[2]) == "legacy";
    const auto directory = std::filesystem::absolute(argv[1]);
    std::filesystem::create_directories(directory);
    const auto self = std::filesystem::canonical("/proc/self/exe");
    for (const auto* name : {"echo_msa", "exit7", "signal", "invalid"})
        std::filesystem::create_symlink(self, directory / name);
    const auto path = directory.string() + ":" + (getenv("PATH") ? getenv("PATH") : "");
    assert(setenv("PATH", path.c_str(), 1) == 0);
    auto& runner = RaMesh::Alignment::ExternalMsaRunner::instance();
    runner.configureScratchDirectory(directory);
    const size_t original_fds = fdCount();
    using Rows = std::unordered_map<ChrName, std::string>;
    const Rows rows{{"z", "ACGTNacgtn"}, {"a", "ACGTTacgtn"}};
    auto aligned = rows;
    assert(runner.align((directory / "echo_msa").string(), aligned));
    assert(aligned == rows);
    aligned = rows;
    assert(runner.align("echo_msa", aligned) && aligned == rows);
    assert(path_launches == 1);
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
    assert(live_attributes == 0 && fdCount() == original_fds);
    std::cout << "PASS: path/PATH, cache, exit 7, signal, invalid output, "
                 "missing executable, 32 concurrent launches, descriptors";
    if (expect_legacy) std::cout << ", attribute failure cleanup (simulated headers)";
    std::cout << '\n';
}
