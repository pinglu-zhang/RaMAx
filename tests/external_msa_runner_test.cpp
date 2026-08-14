#include "external_msa_runner.h"

#include <atomic>
#include <barrier>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <sys/stat.h>
#include <unistd.h>

namespace {

using RaMesh::Alignment::ExternalMsaRunner;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void writeExecutable(const std::filesystem::path& path,
                     const std::string& body) {
    std::ofstream output(path);
    require(static_cast<bool>(output), "cannot create fake MSA executable");
    output << "#!/bin/sh\n" << body;
    output.close();
    require(::chmod(path.c_str(), 0755) == 0,
            "cannot mark fake MSA executable as executable");
}

size_t countMatching(const std::filesystem::path& directory,
                     const std::string& prefix) {
    size_t count = 0;
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (entry.path().filename().string().starts_with(prefix)) ++count;
    }
    return count;
}

size_t countLines(const std::filesystem::path& path) {
    std::ifstream input(path);
    size_t count = 0;
    std::string line;
    while (std::getline(input, line)) ++count;
    return count;
}

std::unordered_map<ChrName, std::string> sequences(size_t suffix = 0) {
    std::string query = "ACGTACGTACGT";
    if (suffix % 2 != 0) query.back() = 'A';
    return {{"ref", "ACGTACGTACGT"}, {"query", std::move(query)}};
}

}  // namespace

int main() {
    const auto root = std::filesystem::temp_directory_path() /
        ("ramax-external-msa-test-" + std::to_string(::getpid()));
    const auto scratch = root / "work" / "minipoa_tmp";
    const auto launcher = root / "launcher";
    std::filesystem::remove_all(root);
    std::filesystem::create_directories(scratch);
    std::filesystem::create_directories(launcher);

    const auto stale = scratch / "minipoa_paths_123.tmp";
    std::ofstream(stale) << "stale";
    const auto unrelated = scratch / "keep.txt";
    std::ofstream(unrelated) << "keep";

    auto& runner = ExternalMsaRunner::instance();
    runner.configureScratchDirectory(scratch);
    require(!std::filesystem::exists(stale),
            "stale internal minipoa file was not removed");
    require(std::filesystem::exists(unrelated),
            "scratch cleanup removed an unrelated file");

    const auto success = root / "success-minipoa";
    writeExecutable(success,
        "touch \"minipoa_paths_$$.tmp\"\ncat \"$5\"\n");
    const auto failure = root / "failure-minipoa";
    writeExecutable(failure,
        "touch \"minipoa_paths_$$.tmp\"\nexit 7\n");
    const auto malformed = root / "malformed-minipoa";
    writeExecutable(malformed,
        "touch \"minipoa_paths_$$.tmp\"\nprintf '>s0\\nACGT\\n'\n");
    const auto cleanup_failure = root / "cleanup-failure-minipoa";
    writeExecutable(cleanup_failure,
        "mkdir \"minipoa_paths_$$.tmp\"\ncat \"$5\"\n");

    const auto original_cwd = std::filesystem::current_path();
    std::filesystem::current_path(launcher);
    auto aligned = sequences();
    require(runner.align(success.string(), aligned),
            "successful fake MSA failed");
    require(countMatching(scratch, "minipoa_paths_") == 0,
            "successful invocation left an internal temp file");
    require(countMatching(launcher, "minipoa_paths_") == 0,
            "child polluted the launcher directory");

    auto failed = sequences();
    require(!runner.align(failure.string(), failed),
            "non-zero fake MSA unexpectedly succeeded");
    require(countMatching(scratch, "minipoa_paths_") == 0,
            "failed invocation left an internal temp file");

    auto malformed_rows = sequences();
    require(!runner.align(malformed.string(), malformed_rows),
            "malformed fake MSA unexpectedly succeeded");
    require(countMatching(scratch, "minipoa_paths_") == 0,
            "malformed invocation left an internal temp file");

    require(::setenv("RAMAX_FORCE_EXTERNAL_MSA_FILE_INPUT", "1", 1) == 0,
            "cannot force external MSA input-file fallback");
    auto fallback = sequences();
    require(runner.align(success.string(), fallback),
            "input-file fallback failed");
    ::unsetenv("RAMAX_FORCE_EXTERNAL_MSA_FILE_INPUT");
    require(countMatching(scratch, "ramax-minipoa-input-") == 0,
            "input-file fallback was not cleaned");
    require(countMatching(launcher, "ramax-minipoa-input-") == 0,
            "input-file fallback polluted the launcher directory");

    auto cleanup_failed = sequences();
    require(runner.align(cleanup_failure.string(), cleanup_failed),
            "cleanup failure changed successful MSA result");
    require(countMatching(scratch, "minipoa_paths_") == 1,
            "cleanup failure was not preserved for inspection");
    for (const auto& entry : std::filesystem::directory_iterator(scratch)) {
        if (entry.path().filename().string().starts_with("minipoa_paths_")) {
            std::filesystem::remove_all(entry.path());
        }
    }

    std::atomic<size_t> parallel_failures{0};
    std::vector<std::thread> workers;
    for (size_t index = 0; index < 16; ++index) {
        workers.emplace_back([&, index] {
            auto rows = sequences(index);
            if (!runner.align(success.string(), rows)) ++parallel_failures;
        });
    }
    for (auto& worker : workers) worker.join();
    require(parallel_failures.load() == 0,
            "parallel fake MSA invocation failed");
    require(countMatching(scratch, "minipoa_paths_") == 0,
            "parallel invocations left internal temp files");
    require(countMatching(launcher, "minipoa_paths_") == 0,
            "parallel invocations polluted the launcher directory");

    const auto cached_success = root / "cached-success-minipoa";
    const auto invocation_log = root / "cached-invocations.txt";
    writeExecutable(cached_success,
        "printf 'invoked\\n' >> \"$RAMAX_TEST_EXTERNAL_MSA_INVOCATIONS\"\n"
        "sleep 0.2\n"
        "touch \"minipoa_paths_$$.tmp\"\n"
        "cat \"$5\"\n");
    require(::setenv("RAMAX_TEST_EXTERNAL_MSA_CACHE_EXECUTABLE",
                     cached_success.c_str(), 1) == 0,
            "cannot enable cache regression test executable");
    require(::setenv("RAMAX_TEST_EXTERNAL_MSA_INVOCATIONS",
                     invocation_log.c_str(), 1) == 0,
            "cannot configure cache regression invocation log");

    std::atomic<size_t> cached_failures{0};
    std::barrier cache_start(16);
    workers.clear();
    for (size_t index = 0; index < 16; ++index) {
        workers.emplace_back([&] {
            auto rows = sequences();
            cache_start.arrive_and_wait();
            if (!runner.align(cached_success.string(), rows)) {
                ++cached_failures;
            }
        });
    }
    for (auto& worker : workers) worker.join();
    require(cached_failures.load() == 0,
            "single-flight cache request failed or stalled");
    require(countLines(invocation_log) == 1,
            "single-flight cache launched duplicate MSA processes");
    auto cached_again = sequences();
    require(runner.align(cached_success.string(), cached_again),
            "completed cache entry could not be reused");
    require(countLines(invocation_log) == 1,
            "cache hit unexpectedly launched another MSA process");
    require(countMatching(scratch, "minipoa_paths_") == 0,
            "cached invocation left an internal temp file");
    ::unsetenv("RAMAX_TEST_EXTERNAL_MSA_INVOCATIONS");
    ::unsetenv("RAMAX_TEST_EXTERNAL_MSA_CACHE_EXECUTABLE");

    std::filesystem::current_path(original_cwd);
    std::filesystem::remove_all(root);
    return 0;
}
