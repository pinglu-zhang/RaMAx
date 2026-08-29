#include "external_tool.h"
#include "wfmash_router.h"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <csignal>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/types.h>
#include <thread>
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

void testMappingChunkPlanning() {
    using WfmashRouterDetail::ParsedPafRecord;

    ParsedPafRecord ordinary;
    ordinary.query_start = 0;
    ordinary.query_end = 1000;
    ordinary.target_start = 0;
    ordinary.target_end = 800;
    ordinary.matches = 90;
    ordinary.block_length = 100;
    require(WfmashRouterDetail::mappingAlignmentCost(ordinary) == 100,
            "mapping alignment cost formula changed");

    ParsedPafRecord saturated;
    saturated.query_end = std::numeric_limits<uint64_t>::max();
    saturated.target_end = std::numeric_limits<uint64_t>::max();
    saturated.matches = 0;
    saturated.block_length = std::numeric_limits<uint64_t>::max();
    require(WfmashRouterDetail::mappingAlignmentCost(saturated) ==
                std::numeric_limits<uint64_t>::max(),
            "mapping alignment cost must saturate");

    require(WfmashRouterDetail::makeMappingChunkPlan({}, 4)
                .chunk_by_record.empty(),
            "empty mapping input must not create chunks");
    const std::vector<uint64_t> costs{10, 10, 10, 10, 9, 1};
    const auto plan = WfmashRouterDetail::makeMappingChunkPlan(costs, 2);
    const auto repeated = WfmashRouterDetail::makeMappingChunkPlan(costs, 2);
    require(plan.chunk_by_record == repeated.chunk_by_record &&
                plan.estimated_cost == repeated.estimated_cost &&
                plan.record_count == repeated.record_count,
            "mapping chunk plan must be deterministic");
    require(plan.chunk_by_record.size() == costs.size() &&
                plan.estimated_cost.size() == 2 &&
                plan.record_count.size() == 2,
            "mapping chunk plan has unexpected dimensions");
    require(std::accumulate(plan.record_count.begin(),
                            plan.record_count.end(), size_t{0}) == costs.size(),
            "mapping chunk plan lost or duplicated records");
    require(std::accumulate(plan.estimated_cost.begin(),
                            plan.estimated_cost.end(), uint64_t{0}) ==
                std::accumulate(costs.begin(), costs.end(), uint64_t{0}),
            "mapping chunk cost accounting changed");
    require(plan.chunk_by_record[0] == 0 && plan.chunk_by_record[1] == 1,
            "LPT tie breaking must prefer the lowest chunk index");

    const auto seqnum_schedule = WfmashRouterDetail::allocateMappingChunks(
        {100, 100, 100}, {10, 10, 10}, 12);
    require(seqnum_schedule == std::vector<size_t>({4, 4, 4}),
            "equal-cost seqnum-4 allocation must be balanced");
    const auto bounded = WfmashRouterDetail::allocateMappingChunks(
        {1000, 1, 1}, {2, 20, 20}, 1000);
    require(std::accumulate(bounded.begin(), bounded.end(), size_t{0}) == 42 &&
                bounded[0] <= 2 && bounded[1] <= 20 && bounded[2] <= 20,
            "chunk allocation must be bounded by mapping record counts");
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

    const auto cancelled_pid_path = root / "cancelled.pid";
    std::atomic<bool> cancellation_requested{false};
    RaMAxExternalTool::RunOptions cancellation = timeout;
    cancellation.timeout = 5s;
    cancellation.cancellation_requested = &cancellation_requested;
    std::thread request_cancellation([&]() {
        std::this_thread::sleep_for(100ms);
        cancellation_requested.store(true, std::memory_order_relaxed);
    });
    const auto cancelled = runShell(
        root,
        "echo $$ > " + cancelled_pid_path.string() +
            "; trap '' TERM; sleep 5 & wait",
        cancellation, "cancelled");
    request_cancellation.join();
    require(cancelled.exit_code == 125 && cancelled.cancelled &&
                !cancelled.timed_out,
            "explicit cancellation must remain distinct from timeout");
    std::ifstream cancelled_pid_input(cancelled_pid_path);
    pid_t cancelled_group = -1;
    cancelled_pid_input >> cancelled_group;
    require(cancelled_group > 0,
            "cancellation test did not record process group");
    errno = 0;
    require(::kill(-cancelled_group, 0) != 0 && errno == ESRCH,
            "cancelled process group still exists");
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

void writeConcurrentChunkWfmash(
    const std::filesystem::path& path,
    const std::filesystem::path& state_directory) {
    const auto lock = state_directory / "lock";
    const auto processes = state_directory / "processes.current";
    const auto maximum_processes = state_directory / "processes.maximum";
    const auto threads = state_directory / "threads.current";
    const auto maximum_threads = state_directory / "threads.maximum";
    std::ostringstream script;
    script
        << "#!/bin/sh\n"
        << "if [ \"$1\" = \"--version\" ]; then\n"
        << "  printf 'v0.14.0-0-g517e1bc\\n'\n"
        << "  exit 0\n"
        << "fi\n"
        << "approx=false\n"
        << "mapping=\n"
        << "declared_threads=0\n"
        << "previous=\n"
        << "for arg in \"$@\"; do\n"
        << "  if [ \"$previous\" = \"-i\" ]; then mapping=$arg; fi\n"
        << "  if [ \"$previous\" = \"-t\" ]; then declared_threads=$arg; fi\n"
        << "  if [ \"$arg\" = \"--approx-map\" ]; then approx=true; fi\n"
        << "  previous=$arg\n"
        << "done\n"
        << "if $approx; then\n"
        << "  index=0\n"
        << "  while [ \"$index\" -lt 4 ]; do\n"
        << "    start=$((index * 100))\n"
        << "    end=$((start + 100))\n"
        << "    printf 'qrychr\\t400\\t%s\\t%s\\t+\\trefchr\\t400\\t%s\\t%s\\t99\\t100\\t60\\tzz:Z:raw-%s\\n' \"$start\" \"$end\" \"$start\" \"$end\" \"$index\"\n"
        << "    index=$((index + 1))\n"
        << "  done\n"
        << "  exit 0\n"
        << "fi\n"
        << "[ -n \"$mapping\" ] || exit 3\n"
        << "exec 9>\"" << lock.string() << "\"\n"
        << "/usr/bin/flock 9\n"
        << "current=$(cat \"" << processes.string() << "\")\n"
        << "current=$((current + 1))\n"
        << "printf '%s\\n' \"$current\" > \"" << processes.string() << "\"\n"
        << "maximum=$(cat \"" << maximum_processes.string() << "\")\n"
        << "if [ \"$current\" -gt \"$maximum\" ]; then printf '%s\\n' \"$current\" > \"" << maximum_processes.string() << "\"; fi\n"
        << "thread_count=$(cat \"" << threads.string() << "\")\n"
        << "thread_count=$((thread_count + declared_threads))\n"
        << "printf '%s\\n' \"$thread_count\" > \"" << threads.string() << "\"\n"
        << "thread_maximum=$(cat \"" << maximum_threads.string() << "\")\n"
        << "if [ \"$thread_count\" -gt \"$thread_maximum\" ]; then printf '%s\\n' \"$thread_count\" > \"" << maximum_threads.string() << "\"; fi\n"
        << "/usr/bin/flock -u 9\n"
        << "sleep 0.30\n"
        << "while IFS= read -r line; do\n"
        << "  [ -n \"$line\" ] && printf '%s\\tcg:Z:100=\\n' \"$line\"\n"
        << "done < \"$mapping\"\n"
        << "/usr/bin/flock 9\n"
        << "current=$(cat \"" << processes.string() << "\")\n"
        << "printf '%s\\n' \"$((current - 1))\" > \"" << processes.string() << "\"\n"
        << "thread_count=$(cat \"" << threads.string() << "\")\n"
        << "printf '%s\\n' \"$((thread_count - declared_threads))\" > \"" << threads.string() << "\"\n"
        << "/usr/bin/flock -u 9\n";
    writeText(path, script.str());
    makeExecutable(path);
}

void writeSiblingFailureWfmash(
    const std::filesystem::path& path,
    const std::filesystem::path& state_directory) {
    const auto slow_pid = state_directory / "query1-slow.pid";
    std::ostringstream script;
    script
        << "#!/bin/sh\n"
        << "if [ \"$1\" = \"--version\" ]; then\n"
        << "  printf 'v0.14.0-0-g517e1bc\\n'\n"
        << "  exit 0\n"
        << "fi\n"
        << "approx=false\n"
        << "mapping=\n"
        << "previous=\n"
        << "for arg in \"$@\"; do\n"
        << "  if [ \"$previous\" = \"-i\" ]; then mapping=$arg; fi\n"
        << "  if [ \"$arg\" = \"--approx-map\" ]; then approx=true; fi\n"
        << "  previous=$arg\n"
        << "done\n"
        << "if $approx; then\n"
        << "  index=0\n"
        << "  while [ \"$index\" -lt 4 ]; do\n"
        << "    start=$((index * 100))\n"
        << "    end=$((start + 100))\n"
        << "    printf 'qrychr\\t400\\t%s\\t%s\\t+\\trefchr\\t400\\t%s\\t%s\\t99\\t100\\t60\\tzz:Z:raw-%s\\n' \"$start\" \"$end\" \"$start\" \"$end\" \"$index\"\n"
        << "    index=$((index + 1))\n"
        << "  done\n"
        << "  exit 0\n"
        << "fi\n"
        << "[ -n \"$mapping\" ] || exit 3\n"
        << "case \"$mapping\" in\n"
        << "  */query1_*/*chunk-00000.mapping.paf)\n"
        << "    printf '%s\\n' \"$$\" > \"" << slow_pid.string() << "\"\n"
        << "    trap '' TERM\n"
        << "    sleep 5 & wait\n"
        << "    ;;\n"
        << "  */query1_*/*chunk-00001.mapping.paf)\n"
        << "    sleep 0.05\n"
        << "    exit 7\n"
        << "    ;;\n"
        << "esac\n"
        << "sleep 0.10\n"
        << "while IFS= read -r line; do\n"
        << "  [ -n \"$line\" ] && printf '%s\\tcg:Z:100=\\n' \"$line\"\n"
        << "done < \"$mapping\"\n";
    writeText(path, script.str());
    makeExecutable(path);
}

std::multiset<std::string> nonemptyLines(const std::string& text) {
    std::multiset<std::string> lines;
    std::istringstream input(text);
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty()) lines.emplace(std::move(line));
    }
    return lines;
}

int readInteger(const std::filesystem::path& path) {
    std::ifstream input(path);
    int value = -1;
    input >> value;
    require(static_cast<bool>(input), "cannot read integer from " + path.string());
    return value;
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

std::filesystem::path speciesDirectory(
    const std::filesystem::path& work_directory,
    std::string_view species) {
    const std::string prefix = std::string(species) + '_';
    std::filesystem::path selected;
    for (const auto& entry :
         std::filesystem::directory_iterator(work_directory)) {
        if (!entry.is_directory() ||
            !entry.path().filename().string().starts_with(prefix)) {
            continue;
        }
        require(selected.empty(),
                "multiple work directories found for " + std::string(species));
        selected = entry.path();
    }
    require(!selected.empty(),
            "work directory is absent for " + std::string(species));
    return selected;
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
            "timed-out alignment must not publish a consolidated success log");
    require(std::filesystem::is_regular_file(
                query_work / "alignment_chunks.tsv"),
            "timed-out alignment must retain its chunk manifest");
    const std::string manifest = RaMAxExternalTool::readText(
        query_work / "alignment_chunks.tsv");
    require(manifest.find("\ttrue\tfalse\t0\t0\tfailed") !=
                std::string::npos,
            "timed-out chunk must not run the ordinary-error retry");
    const std::string routing = RaMAxExternalTool::readText(
        integration_root / "work" / "routing.tsv");
    require(routing.find("timeout_fallback") != std::string::npos &&
            routing.find("alignment-chunk timeout") != std::string::npos &&
            routing.find("failed_chunk=0") != std::string::npos,
            "routing table must identify alignment timeout fallback");
    const std::string tools = RaMAxExternalTool::readText(
        integration_root / "work" / "tools.tsv");
    require(tools.find("minimum_threads_per_process=4") != std::string::npos &&
                tools.find("pair_timeout_ms=250") != std::string::npos &&
                tools.find("maximum_alignment_processes=4") != std::string::npos &&
                tools.find("alignment_chunks_per_worker=4") != std::string::npos,
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
                query_work / "alignment.stderr.log") &&
                std::filesystem::is_regular_file(
                    query_work / "alignment_chunks.tsv"),
            "successful retry must retain consolidated logs and manifest");
    const std::string alignment_log = RaMAxExternalTool::readText(
        query_work / "alignment.stderr.log");
    require(alignment_log.find("attempt1.stderr.log") != std::string::npos &&
                alignment_log.find("attempt2.stderr.log") != std::string::npos,
            "consolidated log must identify both alignment attempts");
    const std::string manifest = RaMAxExternalTool::readText(
        query_work / "alignment_chunks.tsv");
    require(manifest.find("\tfalse\tfalse\t1\t1\tsuccess") !=
                std::string::npos,
            "successful retry must be recorded exactly once in the manifest");
    const std::string routing = RaMAxExternalTool::readText(
        integration_root / "work" / "routing.tsv");
    require(routing.find("wfmash\tsuccess\t1 anchors;mapping_records=1;chunks=1") !=
                std::string::npos,
            "successful retry must remain a wfmash route");
}

void testRouterChunkConcurrencyAndThreadBudget(
    const std::filesystem::path& root) {
    const auto integration_root = root / "router-chunks";
    const auto state = integration_root / "state";
    std::filesystem::create_directories(state);
    for (const auto& name : {"processes.current", "processes.maximum",
                             "threads.current", "threads.maximum"}) {
        writeText(state / name, "0\n");
    }
    const auto samtools = integration_root / "samtools";
    const auto wfmash = integration_root / "wfmash";
    writeFakeSamtools(samtools);
    writeConcurrentChunkWfmash(wfmash, state);

    const auto reference = integration_root / "reference.fa";
    writeText(reference, ">refchr\n" + std::string(400, 'A') + "\n");
    std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
    managers.emplace("ref", managerFor(reference));
    std::vector<MashDistanceRecord> distances;
    for (size_t index = 0; index < 3; ++index) {
        const SpeciesName species = "query" + std::to_string(index + 1);
        const auto query = integration_root / (species + ".fa");
        writeText(query, ">qrychr\n" + std::string(400, 'A') + "\n");
        managers.emplace(species, managerFor(query));
        distances.push_back(
            {"ref", species, 0.001, 0.0, 19900, 20000,
             reference, query});
    }

    WfmashRouterDetail::ExecutionPolicy policy;
    policy.pair_timeout = 5s;
    policy.termination_grace = 50ms;
    policy.poll_interval = 10ms;
    policy.alignment_chunks_per_worker = 2;
    FirstRoundWfmashRouter router(
        samtools, wfmash, integration_root / "work", 32, policy);
    const auto result = router.run(
        "ref", distances, 0.01, 0.02, managers);
    require(result.successful_species.size() == 3 &&
                result.anchors_by_species.size() == 3,
            "all chunked species must complete successfully");
    for (size_t index = 0; index < 3; ++index) {
        const SpeciesName species = "query" + std::to_string(index + 1);
        require(result.successful_species.contains(species) &&
                    result.anchors_by_species.at(species).size() == 4,
                "chunked species lost canonical anchors");
    }
    require(readInteger(state / "processes.maximum") == 3,
            "seqnum-4 policy must expose exactly three alignment slots");
    require(readInteger(state / "threads.maximum") <= 32,
            "concurrent wfmash children exceeded the user thread budget");
    require(readInteger(state / "processes.current") == 0 &&
                readInteger(state / "threads.current") == 0,
            "fake wfmash accounting did not return to zero");

    size_t species_directories = 0;
    for (const auto& entry : std::filesystem::directory_iterator(
             integration_root / "work")) {
        if (!entry.is_directory() || entry.path().filename() == "views") {
            continue;
        }
        ++species_directories;
        const auto manifest = RaMAxExternalTool::readText(
            entry.path() / "alignment_chunks.tsv");
        require(std::count(manifest.begin(), manifest.end(), '\n') == 3,
                "each species must have exactly two chunk manifest rows");
        const auto mappings = RaMAxExternalTool::readText(
            entry.path() / "mappings.paf");
        require(mappings.find("zz:Z:raw-0") != std::string::npos &&
                    mappings.find("zz:Z:raw-3") != std::string::npos,
                "mapping PAF optional tags were not preserved");
        for (const auto& child :
             std::filesystem::directory_iterator(entry.path())) {
            require(child.path().filename().string().find(
                        "alignment-chunks.") != 0,
                    "successful chunk scratch must be removed");
        }
    }
    require(species_directories == 3,
            "unexpected number of chunked species directories");
}

void testChunkFailureCancelsOnlySiblingChunks(
    const std::filesystem::path& root) {
    const auto integration_root = root / "router-sibling-failure";
    const auto state = integration_root / "state";
    std::filesystem::create_directories(state);
    const auto samtools = integration_root / "samtools";
    const auto wfmash = integration_root / "wfmash";
    writeFakeSamtools(samtools);
    writeSiblingFailureWfmash(wfmash, state);

    const auto reference = integration_root / "reference.fa";
    writeText(reference, ">refchr\n" + std::string(400, 'A') + "\n");
    std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
    managers.emplace("ref", managerFor(reference));
    std::vector<MashDistanceRecord> distances;
    for (size_t index = 0; index < 3; ++index) {
        const SpeciesName species = "query" + std::to_string(index + 1);
        const auto query = integration_root / (species + ".fa");
        writeText(query, ">qrychr\n" + std::string(400, 'A') + "\n");
        managers.emplace(species, managerFor(query));
        distances.push_back(
            {"ref", species, 0.001, 0.0, 19900, 20000,
             reference, query});
    }

    WfmashRouterDetail::ExecutionPolicy policy;
    policy.pair_timeout = 5s;
    policy.termination_grace = 50ms;
    policy.poll_interval = 10ms;
    policy.alignment_chunks_per_worker = 2;
    FirstRoundWfmashRouter router(
        samtools, wfmash, integration_root / "work", 32, policy);
    const auto result = router.run(
        "ref", distances, 0.01, 0.02, managers);
    require(!result.successful_species.contains("query1") &&
                !result.anchors_by_species.contains("query1"),
            "failed species must not publish anchors");
    for (const SpeciesName species : {"query2", "query3"}) {
        require(result.successful_species.contains(species) &&
                    result.anchors_by_species.at(species).size() == 4,
                "another species was incorrectly cancelled");
    }

    const auto query_work = speciesDirectory(
        integration_root / "work", "query1");
    require(std::filesystem::is_regular_file(query_work / "mappings.paf"),
            "failed species must retain its mapping PAF");
    require(!std::filesystem::exists(query_work / "alignment.paf"),
            "failed species must not publish a final alignment PAF");
    const std::string manifest = RaMAxExternalTool::readText(
        query_work / "alignment_chunks.tsv");
    require(manifest.find("\t1\t0\tfailed") != std::string::npos &&
                manifest.find("cancelled") != std::string::npos,
            "manifest must record one retrying failure and sibling cancellation");

    std::filesystem::path chunk_root;
    for (const auto& entry : std::filesystem::directory_iterator(query_work)) {
        if (entry.is_directory() && entry.path().filename().string().starts_with(
                "alignment-chunks.")) {
            require(chunk_root.empty(),
                    "failed species retained multiple chunk roots");
            chunk_root = entry.path();
        }
    }
    require(!chunk_root.empty(),
            "failed species must retain chunk inputs and partial outputs");
    std::multiset<std::string> chunk_lines;
    size_t chunk_files = 0;
    for (const auto& entry : std::filesystem::directory_iterator(chunk_root)) {
        const std::string name = entry.path().filename().string();
        if (!entry.is_regular_file() ||
            !name.ends_with(".mapping.paf")) {
            continue;
        }
        ++chunk_files;
        const auto current = nonemptyLines(
            RaMAxExternalTool::readText(entry.path()));
        chunk_lines.insert(current.begin(), current.end());
    }
    require(chunk_files == 2,
            "failed species must retain both deterministic chunk inputs");
    require(chunk_lines == nonemptyLines(
                RaMAxExternalTool::readText(query_work / "mappings.paf")),
            "chunk inputs changed, duplicated, or dropped raw mapping lines");

    std::ifstream slow_pid_input(state / "query1-slow.pid");
    pid_t slow_group = -1;
    slow_pid_input >> slow_group;
    require(slow_group > 0,
            "slow sibling chunk did not start before failure");
    errno = 0;
    require(::kill(-slow_group, 0) != 0 && errno == ESRCH,
            "failed species left a sibling wfmash process group alive");

    const std::string routing = RaMAxExternalTool::readText(
        integration_root / "work" / "routing.tsv");
    require(routing.find("query1\t0.001\tlegacy\tfallback") !=
                std::string::npos &&
                routing.find("query2\t0.001\twfmash\tsuccess") !=
                    std::string::npos &&
                routing.find("query3\t0.001\twfmash\tsuccess") !=
                    std::string::npos,
            "routing did not isolate the failed species fallback");
}

}  // namespace

int main() {
    const auto root = testRoot();
    std::error_code ignored;
    std::filesystem::remove_all(root, ignored);
    std::filesystem::create_directories(root);
    try {
        testPairThreadSchedule();
        testMappingChunkPlanning();
        testExternalToolTimeouts(root);
        testRouterAlignmentTimeoutFallback(root);
        testRouterOrdinaryFailureStillRetries(root);
        testRouterChunkConcurrencyAndThreadBudget(root);
        testChunkFailureCancelsOnlySiblingChunks(root);
    } catch (...) {
        throw;
    }
    std::filesystem::remove_all(root, ignored);
    return 0;
}
