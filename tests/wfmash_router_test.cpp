#include "wfmash_router.h"
#include "external_tool.h"

#include <sstream>
#include <stdexcept>
#include <string>
#include <filesystem>
#include <fstream>
#include <map>

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

void testVersions() {
    require(WfmashRouterDetail::validateSamtoolsVersion(
        "samtools 1.24\nUsing htslib 1.24\n") ==
        "samtools 1.24; HTSlib 1.24", "samtools version");
    expectFailure([] {
        WfmashRouterDetail::validateSamtoolsVersion(
            "samtools 1.24\nUsing htslib 1.23\n");
    });
    require(WfmashRouterDetail::validateWfmashVersion(
        "v0.24.2-0-g774c01ff\n") == "v0.24.2-0-g774c01ff",
        "wfmash version");
    expectFailure([] {
        WfmashRouterDetail::validateWfmashVersion("v0.24.1\n");
    });
    require(RaMAxExternalTool::locateExecutable(
                "ramax-tool-that-does-not-exist", "/missing/ramax/tool").empty(),
            "missing external tool must remain missing");
}

void testFai() {
    std::istringstream valid(
        "chr1\t12\t6\t12\t13\nchr2\t4\t25\t4\t5\r\n");
    const auto records = WfmashRouterDetail::parseFai(valid);
    require(records.size() == 2, "FAI record count");
    require(records[0].name == "chr1" && records[0].length == 12,
            "FAI first record");
    require(records[1].line_width == 5, "FAI CRLF record");

    expectFailure([] {
        std::istringstream malformed("chr1\t12\t6\t12\n");
        WfmashRouterDetail::parseFai(malformed);
    });
    expectFailure([] {
        std::istringstream malformed("chr1\t12\t6\t13\t12\n");
        WfmashRouterDetail::parseFai(malformed);
    });
}

void testPaf() {
    const std::string forward =
        "q1\t100\t10\t42\t+\tt1\t200\t20\t50\t28\t32\t60\tcg:Z:10M2I20M";
    const auto parsed = WfmashRouterDetail::parsePafLine(forward, true);
    require(parsed.strand == Strand::FORWARD, "forward PAF strand");
    require(parsed.cigar.size() == 3, "forward PAF CIGAR");

    const std::string reverse =
        "q2\t80\t5\t37\t-\tt2\t90\t11\t43\t25\t34\t20\tcg:Z:12M2D18M2I";
    const auto reverse_parsed =
        WfmashRouterDetail::parsePafLine(reverse, true);
    require(reverse_parsed.strand == Strand::REVERSE, "reverse PAF strand");

    // wfmash 0.24.2 reports field 11 as max(query span, target span), not
    // necessarily the number of CIGAR columns when both insertions and
    // deletions occur. Coordinate consumption remains the strict invariant.
    const std::string wfmash_indel =
        "q3\t58411566\t270000\t315691\t+\tChr09\t59416394\t271283\t317466"
        "\t45491\t46183\t17\tcg:Z:23011=692D22412=200I68=";
    const auto indel_parsed =
        WfmashRouterDetail::parsePafLine(wfmash_indel, true);
    require(indel_parsed.block_length == 46183,
            "preserve wfmash PAF block length");
    require(indel_parsed.cigar_columns == 46383,
            "derive alignment columns from CIGAR");

    std::vector<WfmashRouterDetail::ParsedPafRecord> overlapping{
        WfmashRouterDetail::parsePafLine(
            "q4\t1000\t0\t500\t+\tt4\t1000\t0\t500\t500\t500\t60\tcg:Z:500=", true),
        WfmashRouterDetail::parsePafLine(
            "q4\t1000\t490\t990\t+\tt4\t1000\t490\t990\t500\t500\t60\tcg:Z:500=", true)
    };
    const auto normalization =
        WfmashRouterDetail::normalizePafForGraph(overlapping);
    require(normalization.trimmed_records == 1 &&
            normalization.skipped_records == 0 && overlapping.size() == 2,
            "overlapping PAF records must be trimmed, not discarded");
    require(overlapping[1].query_start == 500 &&
            overlapping[1].target_start == 500 &&
            overlapping[1].cigar_text == "490=",
            "PAF prefix trim must update coordinates and CIGAR");

    std::vector<WfmashRouterDetail::ParsedPafRecord> reverse_overlapping{
        WfmashRouterDetail::parsePafLine(
            "q5\t1000\t500\t1000\t-\tt5\t1000\t0\t500\t500\t500\t60\tcg:Z:500=", true),
        WfmashRouterDetail::parsePafLine(
            "q5\t1000\t10\t510\t-\tt5\t1000\t490\t990\t500\t500\t60\tcg:Z:500=", true)
    };
    const auto reverse_normalization =
        WfmashRouterDetail::normalizePafForGraph(reverse_overlapping);
    require(reverse_normalization.trimmed_records == 1 &&
            reverse_normalization.skipped_records == 0,
            "reverse PAF overlap must be normalized");
    require(reverse_overlapping[1].query_end == 500 &&
            reverse_overlapping[1].target_start == 500 &&
            reverse_overlapping[1].cigar_text == "490=",
            "reverse PAF prefix trim must update query end");

    expectFailure([&] {
        WfmashRouterDetail::parsePafLine(
            "q1\t100\t10\t40\t+\tt1\t200\t20\t50\t28\t30\t60", true);
    });
    expectFailure([&] {
        WfmashRouterDetail::parsePafLine(
            "q1\t100\t10\t41\t+\tt1\t200\t20\t50\t28\t31\t60\tcg:Z:30M", true);
    });
}

void testScheduling() {
    require(WfmashRouterDetail::threadsPerTask(4, 16) == 4, "16/4 task threads");
    require(WfmashRouterDetail::workerCount(4, 16) == 4, "16/4 workers");
    require(WfmashRouterDetail::threadsPerTask(16, 16) == 1, "16/16 task threads");
    require(WfmashRouterDetail::workerCount(16, 16) == 16, "16/16 workers");
    require(WfmashRouterDetail::threadsPerTask(32, 16) == 1, "16/32 task threads");
    require(WfmashRouterDetail::workerCount(32, 16) == 16, "16/32 workers");
}

std::string makeSequence(size_t length, uint64_t seed) {
    std::string sequence;
    sequence.reserve(length);
    static constexpr char bases[] = "ACGT";
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
        const char current = sequence[i];
        sequence[i] = current == 'A' ? 'C' : 'A';
    }
}

void writeFasta(const std::filesystem::path& path,
                const std::vector<std::pair<std::string, std::string>>& records) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    require(static_cast<bool>(output), "create integration FASTA");
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

void testRealMultiFastaRouting() {
    const auto root = std::filesystem::current_path() /
                      "wfmash_router_integration_artifacts";
    std::error_code ignored;
    std::filesystem::remove_all(root, ignored);
    std::filesystem::create_directories(root);

    const std::string ref1 = makeSequence(40000, 0x123456789ULL);
    const std::string ref2 = makeSequence(42000, 0x987654321ULL);
    std::string q11 = ref1;
    std::string q12 = ref2;
    std::string q21 = ref1;
    std::string q22 = ref2;
    const std::string unhit = makeSequence(30000, 0xabcdef123ULL);
    const std::string fallback1 = makeSequence(41000, 0x13579bdfULL);
    const std::string fallback2 = makeSequence(43000, 0x2468ace0ULL);
    mutate(q11, 251, 17);
    mutate(q12, 263, 29);
    mutate(q21, 277, 31);
    mutate(q22, 283, 43);

    const auto ref_path = root / "reference with spaces.fa";
    const auto q1_path = root / "query one.fa";
    const auto q2_path = root / "query two.fa";
    const auto fallback_path = root / "query fallback.fa";
    writeFasta(ref_path, {{"chr1", ref1}, {"chr2", ref2}});
    writeFasta(q1_path, {{"chr1", q11}, {"private1", q12},
                         {"unhit", unhit}});
    writeFasta(q2_path, {{"chr1", q21}, {"private2", q22}});
    writeFasta(fallback_path, {{"fallback1", fallback1},
                               {"fallback2", fallback2}});

    std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
    managers.emplace("ref", managerFor(ref_path));
    managers.emplace("query one", managerFor(q1_path));
    managers.emplace("query two", managerFor(q2_path));
    managers.emplace("query fallback", managerFor(fallback_path));
    std::vector<MashDistanceRecord> distances{
        {"ref", "query one", 0.001, 0.0, 19900, 20000,
         ref_path, q1_path},
        {"ref", "query two", 0.002, 0.0, 19800, 20000,
         ref_path, q2_path},
        {"ref", "query fallback", 0.003, 0.0, 19700, 20000,
         ref_path, fallback_path}
    };

    FirstRoundWfmashRouter router(
        locateSamtoolsExecutable(), locateWfmashExecutable(),
        root / "work" / "wfmash" / "round_0", 4);
    const auto routed = router.run("ref", distances, 0.01, 0.02, managers);
    require(routed.successful_species.size() == 2,
            "two real wfmash routes must succeed");
    require(!routed.anchors_by_species.at("query one").empty(),
            "query one anchors");
    require(!routed.anchors_by_species.at("query two").empty(),
            "query two anchors");
    require(!routed.successful_species.contains("query fallback"),
            "no-hit query must fall back to legacy");
    require(routed.anchors_by_species.at("query one").size() == 2,
            "unhit contig must not create a false anchor");
    require(std::filesystem::is_regular_file(ref_path.string() + ".fai"),
            "reference samtools FAI");
    require(std::filesystem::is_regular_file(
                root / "work" / "wfmash" / "round_0" /
                "reference.wfmash.index"),
            "shared wfmash index");
    const auto index_path = root / "work" / "wfmash" / "round_0" /
                            "reference.wfmash.index";
    const auto index_mtime = std::filesystem::last_write_time(index_path);
    const auto rerouted = router.run("ref", distances, 0.01, 0.02, managers);
    require(rerouted.successful_species.size() == 2,
            "cached-index reroute must succeed");
    require(std::filesystem::last_write_time(index_path) == index_mtime,
            "shared wfmash index must be reused");

    std::vector<MashDistanceRecord> legacy_distances{
        {"ref", "query one", 0.01, 0.0, 19000, 20000,
         ref_path, q1_path},
        {"ref", "query two", 0.02, 0.0, 18000, 20000,
         ref_path, q2_path}
    };
    const auto legacy = router.run(
        "ref", legacy_distances, 0.01, 0.02, managers);
    require(legacy.successful_species.empty(),
            "equal/greater near threshold must stay legacy");
}

}  // namespace

int main() {
    testVersions();
    testFai();
    testPaf();
    testScheduling();
    testRealMultiFastaRouting();
    return 0;
}
