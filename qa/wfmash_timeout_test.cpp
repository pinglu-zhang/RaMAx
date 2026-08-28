#include "external_tool.h"
#include "wfmash_router.h"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <csignal>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <sys/types.h>
#include <unistd.h>

namespace {

using namespace std::chrono_literals;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

std::filesystem::path testRoot() {
    return std::filesystem::temp_directory_path() /
        ("ramax-wfmash-timeout-" + std::to_string(::getpid()));
}

void writeText(const std::filesystem::path& path, const std::string& text) {
    std::filesystem::create_directories(path.parent_path());
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    require(static_cast<bool>(output), "cannot create " + path.string());
    output << text;
    output.flush();
    require(static_cast<bool>(output), "cannot finish " + path.string());
}

void makeExecutable(const std::filesystem::path& path) {
    std::filesystem::permissions(
        path,
        std::filesystem::perms::owner_exec |
            std::filesystem::perms::group_exec |
            std::filesystem::perms::others_exec,
        std::filesystem::perm_options::add);
}

RaMAxExternalTool::CommandResult runShell(
    const std::filesystem::path& root,
    const std::string& command,
    const RaMAxExternalTool::RunOptions& options,
    std::string_view name) {
    return RaMAxExternalTool::run(
        "/bin/sh", {"-c", command},
        root / (std::string(name) + ".stdout"),
        root / (std::string(name) + ".stderr"), options);
}

void testPairThreadSchedule() {
    struct Case {
        uint_t threads;
        size_t tasks;
        std::vector<uint_t> expected;
    };
    const std::vector<Case> cases{
        {1, 10, {1}},
        {2, 10, {2}},
        {3, 10, {3}},
        {4, 10, {4}},
        {8, 10, {4, 4}},
        {18, 10, {5, 5, 4, 4}},
        {32, 3, {11, 11, 10}},
        {32, 20, {4, 4, 4, 4, 4, 4, 4, 4}},
        {64, 20, {4, 4, 4, 4, 4, 4, 4, 4,
                  4, 4, 4, 4, 4, 4, 4, 4}},
    };
    for (const auto& test : cases) {
        const auto schedule = WfmashRouterDetail::pairThreadSchedule(
            test.tasks, test.threads, 4);
        require(schedule.threads_per_worker == test.expected,
                "unexpected pair thread schedule");
        const uint_t allocated = std::accumulate(
            schedule.threads_per_worker.begin(),
            schedule.threads_per_worker.end(), 0U);
        require(allocated <= test.threads, "thread schedule oversubscribes");
        const auto [minimum, maximum] = std::minmax_element(
            schedule.threads_per_worker.begin(),
            schedule.threads_per_worker.end());
        require(*maximum - *minimum <= 1, "worker budgets are unbalanced");
        if (test.threads >= 4) {
            require(*minimum >= 4, "worker received fewer than four threads");
        }
    }
    require(WfmashRouterDetail::pairThreadSchedule(0, 32, 4).workers() == 0,
            "empty task set must not create workers");
}

void testExternalToolTimeouts(const std::filesystem::path& root) {
    RaMAxExternalTool::RunOptions blocking;
    const auto normal = runShell(root, "exit 0", blocking, "normal");
    require(normal.exit_code == 0 && !normal.timed_out,
            "normal command result");

    const auto natural_124 = runShell(root, "exit 124", blocking, "exit124");
    require(natural_124.exit_code == 124 && !natural_124.timed_out,
            "natural exit 124 must not be a timeout");

    RaMAxExternalTool::RunOptions timeout;
    timeout.timeout = 100ms;
    timeout.termination_grace = 50ms;
    timeout.poll_interval = 10ms;
    timeout.create_process_group = true;
    const auto terminated = runShell(root, "sleep 5", timeout, "term");
    require(terminated.exit_code == 124 && terminated.timed_out,
            "sleep must time out");
    require(terminated.termination_signal == SIGTERM,
            "cooperative process must stop on SIGTERM");

    const auto pid_path = root / "ignored-term.pid";
    const std::string ignored_command =
        "echo $$ > " + pid_path.string() +
        "; trap '' TERM; sleep 5 & wait";
    const auto killed = runShell(root, ignored_command, timeout, "kill");
    require(killed.exit_code == 124 && killed.timed_out,
            "TERM-ignoring process must time out");
    require(killed.termination_signal == SIGKILL,
            "TERM-ignoring process must be killed with SIGKILL");

    std::ifstream pid_input(pid_path);
    pid_t process_group = -1;
    pid_input >> process_group;
    require(process_group > 0, "timeout test did not record process group");
    errno = 0;
    require(::kill(-process_group, 0) != 0 && errno == ESRCH,
            "timed-out process group still exists");
}

void writeFakeSamtools(const std::filesystem::path& path) {
    writeText(path,
        "#!/bin/sh\n"
        "if [ \"$1\" = \"--version\" ]; then\n"
        "  printf 'samtools 1.23.1\\nUsing htslib 1.23.1\\n'\n"
        "  exit 0\n"
        "fi\n"
        "if [ \"$1\" = \"faidx\" ] && [ \"$2\" = \"--fai-idx\" ]; then\n"
        "  out=$3\n"
        "  fasta=$4\n"
        "  name=$(sed -n '1s/^>//p' \"$fasta\")\n"
        "  length=$(tail -n +2 \"$fasta\" | tr -d '\\r\\n' | wc -c)\n"
        "  printf '%s\\t%s\\t0\\t%s\\t%s\\n' \"$name\" \"$length\" "
        "\"$length\" \"$((length + 1))\" > \"$out\"\n"
        "  exit 0\n"
        "fi\n"
        "exit 2\n");
    makeExecutable(path);
}

void writeAlignmentTimeoutWfmash(const std::filesystem::path& path) {
    writeText(path,
        "#!/bin/sh\n"
        "if [ \"$1\" = \"--version\" ]; then\n"
        "  printf 'v0.14.0-0-g517e1bc\\n'\n"
        "  exit 0\n"
        "fi\n"
        "for arg in \"$@\"; do\n"
        "  if [ \"$arg\" = \"--approx-map\" ]; then\n"
        "    printf 'qrychr\\t100\\t0\\t100\\t+\\trefchr\\t100\\t0\\t100\\t100\\t100\\t60\\n'\n"
        "    exit 0\n"
        "  fi\n"
        "done\n"
        "trap '' TERM\n"
        "sleep 5 & wait\n");
    makeExecutable(path);
}

void writeRetrySuccessWfmash(const std::filesystem::path& path) {
    writeText(path,
        "#!/bin/sh\n"
        "if [ \"$1\" = \"--version\" ]; then\n"
        "  printf 'v0.14.0-0-g517e1bc\\n'\n"
        "  exit 0\n"
        "fi\n"
        "for arg in \"$@\"; do\n"
        "  if [ \"$arg\" = \"--approx-map\" ]; then\n"
        "    printf 'qrychr\\t100\\t0\\t100\\t+\\trefchr\\t100\\t0\\t100\\t100\\t100\\t60\\n'\n"
        "    exit 0\n"
        "  fi\n"
        "done\n"
        "for arg in \"$@\"; do\n"
        "  case \"$arg\" in\n"
        "    *alignment.attempt1*) exit 7 ;;\n"
        "  esac\n"
        "done\n"
        "printf 'qrychr\\t100\\t0\\t100\\t+\\trefchr\\t100\\t0\\t100\\t100\\t100\\t60\\tcg:Z:100=\\n'\n");
    makeExecutable(path);
}

SeqPro::SharedManagerVariant managerFor(const std::filesystem::path& path) {
    SeqPro::ManagerVariant manager{
        std::make_unique<SeqPro::SequenceManager>(path)};
    return std::make_shared<SeqPro::ManagerVariant>(std::move(manager));
}

std::filesystem::path singleSpeciesDirectory(
    const std::filesystem::path& work_directory) {
    std::vector<std::filesystem::path> directories;
    for (const auto& entry :
         std::filesystem::directory_iterator(work_directory)) {
        if (entry.is_directory() && entry.path().filename() != "views") {
            directories.push_back(entry.path());
        }
    }
    require(directories.size() == 1,
            "expected exactly one species work directory");
    return directories.front();
}

void testRouterAlignmentTimeoutFallback(const std::filesystem::path& root) {
    const auto integration_root = root / "router";
    const auto samtools = integration_root / "samtools";
    const auto wfmash = integration_root / "wfmash";
    writeFakeSamtools(samtools);
    writeAlignmentTimeoutWfmash(wfmash);

    const auto reference = integration_root / "reference.fa";
    const auto query = integration_root / "query.fa";
    writeText(reference, ">refchr\n" + std::string(100, 'A') + "\n");
    writeText(query, ">qrychr\n" + std::string(100, 'A') + "\n");

    std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
    managers.emplace("ref", managerFor(reference));
    managers.emplace("query", managerFor(query));
    std::vector<MashDistanceRecord> distances{
        {"ref", "query", 0.001, 0.0, 19900, 20000,
         reference, query}
    };

    WfmashRouterDetail::ExecutionPolicy policy;
    policy.pair_timeout = 250ms;
    policy.termination_grace = 50ms;
    policy.poll_interval = 10ms;
    FirstRoundWfmashRouter router(
        samtools, wfmash, integration_root / "work", 8, policy);
    const auto result = router.run("ref", distances, 0.01, 0.02, managers);
    require(result.successful_species.empty(),
            "timed-out query must not be marked successful");
    require(result.anchors_by_species.empty(),
            "timed-out query must not publish anchors");

    const auto query_work = singleSpeciesDirectory(integration_root / "work");
    require(!std::filesystem::exists(query_work / "alignment.paf"),
            "timed-out alignment must not be published");
    require(!std::filesystem::exists(query_work / "alignment.stderr.log"),
            "timed-out first alignment must not run the retry");
    const std::string routing = RaMAxExternalTool::readText(
        integration_root / "work" / "routing.tsv");
    require(routing.find("timeout_fallback") != std::string::npos &&
            routing.find("alignment-attempt1 timeout") != std::string::npos,
            "routing table must identify alignment timeout fallback");
    const std::string tools = RaMAxExternalTool::readText(
        integration_root / "work" / "tools.tsv");
    require(tools.find("minimum_threads_per_process=4;pair_timeout_ms=250") !=
                std::string::npos,
            "tools table must record the effective execution policy");
}

void testRouterOrdinaryFailureStillRetries(const std::filesystem::path& root) {
    const auto integration_root = root / "router-retry";
    const auto samtools = integration_root / "samtools";
    const auto wfmash = integration_root / "wfmash";
    writeFakeSamtools(samtools);
    writeRetrySuccessWfmash(wfmash);

    const auto reference = integration_root / "reference.fa";
    const auto query = integration_root / "query.fa";
    writeText(reference, ">refchr\n" + std::string(100, 'A') + "\n");
    writeText(query, ">qrychr\n" + std::string(100, 'A') + "\n");

    std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
    managers.emplace("ref", managerFor(reference));
    managers.emplace("query", managerFor(query));
    std::vector<MashDistanceRecord> distances{
        {"ref", "query", 0.001, 0.0, 19900, 20000,
         reference, query}
    };

    WfmashRouterDetail::ExecutionPolicy policy;
    policy.pair_timeout = 2s;
    policy.termination_grace = 50ms;
    policy.poll_interval = 10ms;
    FirstRoundWfmashRouter router(
        samtools, wfmash, integration_root / "work", 8, policy);
    const auto result = router.run("ref", distances, 0.01, 0.02, managers);
    require(result.successful_species.contains("query"),
            "ordinary alignment failure must retain one retry");
    require(result.anchors_by_species.at("query").size() == 1,
            "successful retry must publish its anchor");

    const auto query_work = singleSpeciesDirectory(integration_root / "work");
    require(std::filesystem::is_regular_file(query_work / "mappings.paf") &&
            std::filesystem::is_regular_file(query_work / "alignment.paf"),
            "successful retry must publish final PAF files");
    require(std::filesystem::is_regular_file(
                query_work / "alignment.attempt1.stderr.log") &&
            std::filesystem::is_regular_file(
                query_work / "alignment.stderr.log"),
            "both alignment attempt logs must be retained");
    const std::string routing = RaMAxExternalTool::readText(
        integration_root / "work" / "routing.tsv");
    require(routing.find("wfmash\tsuccess\t1 anchors") != std::string::npos,
            "successful retry must remain a wfmash route");
}

}  // namespace

int main() {
    const auto root = testRoot();
    std::error_code ignored;
    std::filesystem::remove_all(root, ignored);
    std::filesystem::create_directories(root);
    try {
        testPairThreadSchedule();
        testExternalToolTimeouts(root);
        testRouterAlignmentTimeoutFallback(root);
        testRouterOrdinaryFailureStillRetries(root);
    } catch (...) {
        throw;
    }
    std::filesystem::remove_all(root, ignored);
    return 0;
}
