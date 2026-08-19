#include "mm2plus_router.h"
#include "external_tool.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template<class Function>
void expectFailure(Function&& function) {
    bool failed = false;
    try {
        function();
    } catch (const std::exception&) {
        failed = true;
    }
    require(failed, "expected operation to fail");
}

void testVersion() {
    require(Mm2plusRouterDetail::validateMm2plusVersion(
        "1.3\n", "Usage: mm2plus [options] <target.fa> [query.fa]\n") ==
        "mm2plus 1.3; upstream Minimap2 2.31-r1302", "mm2plus version");
    expectFailure([] {
        Mm2plusRouterDetail::validateMm2plusVersion(
            "1.2\n", "Usage: mm2plus [options]\n");
    });
    expectFailure([] {
        Mm2plusRouterDetail::validateMm2plusVersion(
            "1.3\n", "Usage: minimap2 [options]\n");
    });
}

void testArguments() {
    const auto index = Mm2plusRouterDetail::indexArguments(
        3, 1000, "/tmp/reference.asm20.mmi.part", "/tmp/reference.fa");
    const std::vector<std::string> expected_index{
        "-x", "asm20", "-I", "1001", "-t", "3", "-d",
        "/tmp/reference.asm20.mmi.part", "/tmp/reference.fa"};
    require(index == expected_index, "single-part index arguments");

    const auto alignment = Mm2plusRouterDetail::alignmentArguments(
        16, "/tmp/reference.asm20.mmi", "/tmp/query with spaces.fa");
    const std::vector<std::string> expected_alignment{
        "-x", "asm20", "-c", "--eqx", "--secondary=no", "-t", "16",
        "/tmp/reference.asm20.mmi", "/tmp/query with spaces.fa"};
    require(alignment == expected_alignment, "alignment arguments");
    require(std::find(alignment.begin(), alignment.end(), "-N") == alignment.end(),
            "must not use -N 0");
    require(std::find(alignment.begin(), alignment.end(), "--cs") == alignment.end(),
            "must not emit cs");
}

void testScheduling() {
    require(Mm2plusRouterDetail::threadsPerPair(1) == 1, "one thread");
    require(Mm2plusRouterDetail::threadsPerPair(12) == 12, "twelve threads");
    require(Mm2plusRouterDetail::threadsPerPair(32) == 16, "thread cap");
    require(Mm2plusRouterDetail::workerCount(4, 12) == 1, "12 thread wave");
    require(Mm2plusRouterDetail::workerCount(4, 32) == 2, "32/4 waves");
    require(Mm2plusRouterDetail::workerCount(2, 64) == 2, "64/2 cap");
}

void testStrictPafAndNormalization() {
    auto make = [](std::string_view line) {
        return WfmashRouterDetail::parsePafLine(line, true);
    };
    std::vector<WfmashRouterDetail::ParsedPafRecord> records{
        make("q\t1000\t0\t500\t+\tt\t1000\t0\t500\t490\t500\t60"
             "\ttp:A:P\tAS:i:900\tdv:f:0.02\tcg:Z:490=10X"),
        make("q\t1000\t0\t500\t+\tt\t1000\t0\t500\t490\t500\t60"
             "\ttp:A:P\tAS:i:900\tdv:f:0.02\tcg:Z:490=10X"),
        make("q\t1000\t490\t990\t+\tt\t1000\t490\t990\t500\t500\t40"
             "\ttp:A:P\tAS:i:800\tcg:Z:500=")
    };
    const auto stats = Mm2plusRouterDetail::normalizeForGraph(records, 65);
    require(stats.raw_records == 3 && stats.exact_duplicates == 1,
            "exact duplicate accounting");
    require(stats.trimmed_records == 1 && stats.overlap_discarded == 0,
            "boundary overlap trimming");
    require(records.size() == 2, "canonical PAF record count");
    const auto trimmed = std::find_if(records.begin(), records.end(),
        [](const auto& record) { return record.query_end == 990; });
    require(trimmed != records.end() && trimmed->query_start == 500 &&
            trimmed->target_start == 500 && trimmed->cigar_text == "490=",
            "score-aware lower record trim");

    const auto parsed = make(
        "q\t1000\t0\t100\t-\tt\t1000\t10\t110\t90\t100\t30"
        "\ttp:A:I\tAS:i:100\tcg:Z:90=10X");
    require(parsed.alignment_type == 'I' && parsed.alignment_score == 100 &&
            parsed.equal_bases == 90 && parsed.mismatch_bases == 10,
            "mm2plus optional tags and eqx accounting");
    require(Mm2plusRouterDetail::passesGraphQualityFilters(parsed),
            "positive-score primary passes quality filter");
    const auto negative_score = make(
        "q\t1000\t0\t100\t+\tt\t1000\t10\t110\t90\t100\t60"
        "\ttp:A:P\tAS:i:-10\tcg:Z:90=10X");
    require(!Mm2plusRouterDetail::passesGraphQualityFilters(negative_score),
            "negative alignment score is filtered per record");
    const auto low_identity = make(
        "q\t1000\t0\t100\t+\tt\t1000\t10\t110\t49\t100\t60"
        "\ttp:A:P\tAS:i:10\tcg:Z:49=51X");
    require(!Mm2plusRouterDetail::passesGraphQualityFilters(low_identity),
            "sub-50-percent identity is filtered per record");
}

std::string makeSequence(size_t length, uint64_t seed) {
    static constexpr char bases[] = "ACGT";
    std::string sequence;
    sequence.reserve(length);
    uint64_t state = seed;
    for (size_t i = 0; i < length; ++i) {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        sequence.push_back(bases[state & 3U]);
    }
    return sequence;
}

void mutate(std::string& sequence, size_t stride, size_t phase) {
    for (size_t i = phase; i < sequence.size(); i += stride) {
        sequence[i] = sequence[i] == 'A' ? 'C' : 'A';
    }
}

void writeFasta(const std::filesystem::path& path,
                const std::vector<std::pair<std::string, std::string>>& records) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    require(static_cast<bool>(output), "write FASTA");
    for (const auto& [name, sequence] : records) {
        output << '>' << name << '\n';
        for (size_t offset = 0; offset < sequence.size(); offset += 80) {
            output.write(sequence.data() + offset,
                static_cast<std::streamsize>(
                    std::min<size_t>(80, sequence.size() - offset)));
            output.put('\n');
        }
    }
}

SeqPro::SharedManagerVariant managerFor(const std::filesystem::path& path) {
    SeqPro::ManagerVariant variant{
        std::make_unique<SeqPro::SequenceManager>(path)};
    return std::make_shared<SeqPro::ManagerVariant>(std::move(variant));
}

void testRealRoutingIfAvailable() {
    const auto executable = RaMAxExternalTool::searchPath("mm2plus");
    if (executable.empty()) return;

    const auto root = std::filesystem::current_path() /
                      "mm2plus_router_integration_artifacts";
    std::error_code ignored;
    std::filesystem::remove_all(root, ignored);
    std::filesystem::create_directories(root);

    const std::string ref1 = makeSequence(50000, 0x123456789ULL);
    const std::string ref2 = makeSequence(52000, 0x987654321ULL);
    std::string query1 = ref1;
    std::string query2 = ref2;
    mutate(query1, 29, 7);
    mutate(query2, 31, 11);
    const auto ref_path = root / "reference with spaces.fa";
    const auto query_path = root / "query with spaces.fa";
    writeFasta(ref_path, {{"chr1", ref1}, {"chr2", ref2}});
    writeFasta(query_path, {{"chr1", query1}, {"private", query2}});

    std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
    managers.emplace("ref", managerFor(ref_path));
    managers.emplace("query", managerFor(query_path));
    const std::vector<MashDistanceRecord> distances{
        {"ref", "query", 0.03, 0.0, 8000, 20000,
         ref_path, query_path}
    };
    FirstRoundMm2plusRouter router(
        executable, root / "work" / "mm2plus" / "round_0", 4);
    const auto result = router.run("ref", distances, 0.02, 65, managers);
    require(result.successful_species.contains("query"),
            "real mm2plus route");
    require(!result.anchors_by_species.at("query").empty(),
            "real mm2plus anchors");
    std::filesystem::path canonical_path;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(
             root / "work" / "mm2plus" / "round_0")) {
        if (entry.path().filename() == "alignment.paf") {
            canonical_path = entry.path();
            break;
        }
    }
    require(!canonical_path.empty(), "canonical alignment PAF");
    const auto canonical = RaMAxExternalTool::readText(canonical_path);
    require(canonical.find("query.chr1") != std::string::npos &&
            canonical.find("ref.chr1") != std::string::npos,
            "canonical species.contig names");
}

}  // namespace

int main() {
    testVersion();
    testArguments();
    testScheduling();
    testStrictPafAndNormalization();
    testRealRoutingIfAvailable();
    return 0;
}
