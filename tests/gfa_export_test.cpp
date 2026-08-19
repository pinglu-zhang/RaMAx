#include "align.h"
#include "ramesh.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

namespace fs = std::filesystem;

void require(bool condition, const std::string& message) {
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

std::vector<std::string> split(const std::string& value, char delimiter) {
    std::vector<std::string> fields;
    std::stringstream stream(value);
    std::string field;
    while (std::getline(stream, field, delimiter)) fields.push_back(field);
    return fields;
}

void writeFasta(
    const fs::path& path,
    const std::vector<std::pair<std::string, std::string>>& records) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    require(static_cast<bool>(output), "create test FASTA");
    for (const auto& [name, sequence] : records) {
        output << '>' << name << '\n' << sequence << '\n';
    }
}

SeqPro::SharedManagerVariant managerFor(const fs::path& path) {
    SeqPro::ManagerVariant variant{
        std::make_unique<SeqPro::SequenceManager>(path)};
    return std::make_shared<SeqPro::ManagerVariant>(std::move(variant));
}

std::string readFile(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    require(static_cast<bool>(input), "open test output");
    return std::string(std::istreambuf_iterator<char>(input),
                       std::istreambuf_iterator<char>());
}

char complement(char base) {
    switch (base) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    case 'N': return 'N';
    default: throw std::runtime_error("invalid test base");
    }
}

struct ParsedGfa {
    std::string header;
    std::vector<std::string> segments_and_links;
    std::map<std::string, std::string> nodes;
    std::map<std::string, std::vector<std::pair<std::string, bool>>> walks;
    std::size_t paths{0};
    std::size_t wlines{0};
};

ParsedGfa parseGfa(const fs::path& path) {
    ParsedGfa parsed;
    std::ifstream input(path);
    require(static_cast<bool>(input), "open test GFA");
    std::string line;
    while (std::getline(input, line)) {
        const auto fields = split(line, '\t');
        require(!fields.empty(), "non-empty GFA line");
        if (fields[0] == "H") {
            parsed.header = line;
        } else if (fields[0] == "S") {
            require(fields.size() >= 3, "valid S-line");
            parsed.segments_and_links.push_back(line);
            require(parsed.nodes.emplace(fields[1], fields[2]).second,
                    "unique node id");
        } else if (fields[0] == "L") {
            require(fields.size() >= 6, "valid L-line");
            parsed.segments_and_links.push_back(line);
            require(fields[5] == "0M", "zero-overlap link");
        } else if (fields[0] == "P") {
            require(fields.size() >= 4, "valid P-line");
            ++parsed.paths;
            auto& walk = parsed.walks[fields[1]];
            for (const auto& token : split(fields[2], ',')) {
                require(token.size() >= 2, "valid P step");
                const char orientation = token.back();
                require(orientation == '+' || orientation == '-',
                        "valid P orientation");
                walk.emplace_back(token.substr(0, token.size() - 1),
                                  orientation == '-');
            }
        } else if (fields[0] == "W") {
            require(fields.size() >= 7, "valid W-line");
            ++parsed.wlines;
            const std::string key = fields[1] + "." + fields[3];
            auto& walk = parsed.walks[key];
            const std::string& encoded = fields[6];
            std::size_t cursor = 0;
            while (cursor < encoded.size()) {
                const bool reverse = encoded[cursor] == '<';
                require(reverse || encoded[cursor] == '>',
                        "valid W orientation");
                const std::size_t start = ++cursor;
                while (cursor < encoded.size() &&
                       encoded[cursor] != '<' && encoded[cursor] != '>') {
                    ++cursor;
                }
                require(cursor > start, "valid W node id");
                walk.emplace_back(encoded.substr(start, cursor - start),
                                  reverse);
            }
        }
    }
    return parsed;
}

std::string spell(
    const ParsedGfa& parsed,
    const std::vector<std::pair<std::string, bool>>& walk) {
    std::string sequence;
    for (const auto& [node, reverse] : walk) {
        const auto found = parsed.nodes.find(node);
        require(found != parsed.nodes.end(), "walk node exists");
        if (!reverse) {
            sequence += found->second;
        } else {
            for (auto it = found->second.rbegin();
                 it != found->second.rend(); ++it) {
                sequence.push_back(complement(*it));
            }
        }
    }
    return sequence;
}

void testVersionsAndLosslessWalks(const fs::path& root) {
    const std::string ref_chr1 = "AACCGGTTNNACGT";
    const std::string qry_chr1 = "AACCTGTTNNACGT";
    const fs::path reference = root / "reference multi.fa";
    const fs::path query = root / "query multi.fa";
    writeFasta(reference, {{"chr1", ref_chr1}, {"private", "TTTT"}});
    writeFasta(query, {{"chr1", qry_chr1}, {"private", "GGGG"}});

    std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
    managers.emplace("ref", managerFor(reference));
    managers.emplace("query", managerFor(query));
    RaMesh::RaMeshMultiGenomeGraph graph(managers);
    graph.reference_order = {"ref", "query"};

    Cigar_t cigar;
    parseCigarString(std::to_string(ref_chr1.size()) + "M", cigar);
    Anchor anchor(
        0, 0, static_cast<Length_t>(ref_chr1.size()),
        0, 0, static_cast<Length_t>(qry_chr1.size()),
        Strand::FORWARD, static_cast<uint_t>(ref_chr1.size()),
        static_cast<uint_t>(ref_chr1.size()), cigar);
    graph.insertAnchorIntoGraph(
        *managers.at("ref"), *managers.at("query"),
        "ref", "query", anchor, true);

    const fs::path maf_before = root / "before.maf";
    const fs::path paf_before = root / "before.paf";
    graph.exportToMaf(maf_before, managers, true, false);
    graph.exportToPaf(paf_before, managers);

    RaMesh::Gfa::GfaExportOptions options;
    options.threads = 2;
    options.version = RaMesh::Gfa::Version::V1_0;
    const fs::path gfa10 = root / "graph-1.0.gfa";
    const auto stats10 = graph.exportToGfa(gfa10, managers, options);
    options.version = RaMesh::Gfa::Version::V1_1;
    const fs::path gfa11 = root / "graph-1.1.gfa";
    const auto stats11 = graph.exportToGfa(gfa11, managers, options);
    const fs::path maf_after = root / "after.maf";
    const fs::path paf_after = root / "after.paf";
    graph.exportToMaf(maf_after, managers, true, false);
    graph.exportToPaf(paf_after, managers);

    const ParsedGfa parsed10 = parseGfa(gfa10);
    const ParsedGfa parsed11 = parseGfa(gfa11);
    require(parsed10.header == "H\tVN:Z:1.0", "GFA 1.0 header");
    require(parsed11.header == "H\tVN:Z:1.1\tRS:Z:ref",
            "GFA 1.1 header and reference sample");
    require(parsed10.paths == 4 && parsed10.wlines == 0,
            "GFA 1.0 uses only P-lines");
    require(parsed11.paths == 0 && parsed11.wlines == 4,
            "GFA 1.1 uses only W-lines");
    require(parsed10.segments_and_links == parsed11.segments_and_links,
            "GFA versions share byte-identical S/L records");
    require(parsed10.walks == parsed11.walks,
            "P-lines and W-lines represent identical logical walks");
    require(spell(parsed10, parsed10.walks.at("ref.chr1")) == ref_chr1,
            "reference chr1 lossless");
    require(spell(parsed10, parsed10.walks.at("query.chr1")) == qry_chr1,
            "query chr1 lossless");
    require(spell(parsed11, parsed11.walks.at("ref.private")) == "TTTT",
            "reference private contig lossless");
    require(spell(parsed11, parsed11.walks.at("query.private")) == "GGGG",
            "query private contig lossless");
    require(stats10.walks == 4 && stats11.walks == 4,
            "one graph walk per input contig");
    require(stats10.nodes == stats11.nodes &&
            stats10.edges == stats11.edges &&
            stats10.exact_runs == stats11.exact_runs,
            "version does not alter graph IR");
    require(stats10.shared_sequence_bp > 0,
            "aligned exact ACGT runs are shared");
    require(readFile(maf_before) == readFile(maf_after),
            "GFA export does not mutate MAF output");
    require(readFile(paf_before) == readFile(paf_after),
            "GFA export does not mutate PAF output");
}

void testPathNameCollision(const fs::path& root) {
    const fs::path first = root / "collision-first.fa";
    const fs::path second = root / "collision-second.fa";
    writeFasta(first, {{"c", "ACGT"}});
    writeFasta(second, {{"b.c", "ACGT"}});
    std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
    managers.emplace("a.b", managerFor(first));
    managers.emplace("a", managerFor(second));
    RaMesh::RaMeshMultiGenomeGraph graph(managers);
    graph.reference_order = {"a.b", "a"};

    RaMesh::Gfa::GfaExportOptions options;
    options.version = RaMesh::Gfa::Version::V1_0;
    expectFailure([&] {
        graph.exportToGfa(root / "collision-1.0.gfa", managers, options);
    });
    require(!fs::exists(root / "collision-1.0.gfa"),
            "failed GFA 1.0 is not published");

    options.version = RaMesh::Gfa::Version::V1_1;
    graph.exportToGfa(root / "collision-1.1.gfa", managers, options);
    const ParsedGfa parsed = parseGfa(root / "collision-1.1.gfa");
    require(parsed.wlines == 2, "W metadata disambiguates P-name collision");
}

struct CompactFixtureResult {
    RaMesh::Gfa::GfaExportStats exact_stats;
    RaMesh::Gfa::GfaExportStats compact_stats;
    fs::path exact;
    fs::path compact;
    fs::path work;
    std::string reference;
    std::string query;
};

CompactFixtureResult runCompactFixture(const fs::path& root,
                                       std::size_t short_run,
                                       int threads,
                                       const std::string& label,
                                       std::size_t minimum_exact_run = 25) {
    const std::string shared_prefix(10000, 'A');
    CompactFixtureResult result;
    result.reference = shared_prefix + "C" + std::string(short_run, 'G');
    result.query = shared_prefix + "T" + std::string(short_run, 'G');
    const fs::path reference = root / (label + "-reference.fa");
    const fs::path query = root / (label + "-query.fa");
    writeFasta(reference, {{"chr1", result.reference}});
    writeFasta(query, {{"chr1", result.query}});

    std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
    managers.emplace("ref", managerFor(reference));
    managers.emplace("query", managerFor(query));
    RaMesh::RaMeshMultiGenomeGraph graph(managers);
    graph.reference_order = {"ref", "query"};
    Cigar_t cigar;
    parseCigarString(std::to_string(result.reference.size()) + "M", cigar);
    Anchor anchor(
        0, 0, static_cast<Length_t>(result.reference.size()),
        0, 0, static_cast<Length_t>(result.query.size()),
        Strand::FORWARD, static_cast<uint_t>(result.reference.size()),
        static_cast<uint_t>(result.query.size()), cigar);
    graph.insertAnchorIntoGraph(
        *managers.at("ref"), *managers.at("query"),
        "ref", "query", anchor, true);

    RaMesh::Gfa::GfaExportOptions options;
    options.version = RaMesh::Gfa::Version::V1_1;
    options.threads = threads;
    result.exact = root / (label + "-exact.gfa");
    result.exact_stats = graph.exportToGfa(result.exact, managers, options);
    options.profile = RaMesh::Gfa::Profile::COMPACT;
    options.compact.minimum_exact_run_bp = minimum_exact_run;
    result.work = root / (label + "-work");
    options.work_dir = result.work;
    result.compact = root / (label + "-compact.gfa");
    result.compact_stats = graph.exportToGfa(
        result.compact, managers, options);
    return result;
}

void testCompactThresholdAndDeterminism(const fs::path& root) {
    const auto below = runCompactFixture(root, 24, 1, "threshold-24-t1");
    require(below.compact_stats.suppressed_exact_runs == 1,
            "24 bp exact relation is suppressed");
    require(below.compact_stats.nodes < below.exact_stats.nodes,
            "short relation suppression reduces nodes");
    require(readFile(below.exact) ==
                readFile(below.work / "gfa" / "exact.gfa"),
            "compact exact shadow is byte-identical to exact output");
    require(fs::exists(below.work / "gfa" / "compact_transform.tsv") &&
            fs::exists(below.work / "gfa" / "compact_stats.tsv") &&
            fs::exists(below.work / "gfa" / "compact_rejections.tsv") &&
            fs::exists(below.work / "gfa" / "compact_parameters.tsv"),
            "compact audit reports are published");
    const std::string stage_stats = readFile(
        below.work / "gfa" / "compact_stats.tsv");
    require(stage_stats.find(
                "metric\texact\tafter_filter\tafter_unitig\t"
                "after_allele_islands\tfinal") != std::string::npos,
            "compact report contains all transform stages");
    const std::string parameters = readFile(
        below.work / "gfa" / "compact_parameters.tsv");
    require(parameters.find("profile\tcompact-v2-balanced") !=
                std::string::npos &&
            parameters.find("minimum_exact_run_bp\t25") !=
                std::string::npos,
            "compact report records the profile and effective configuration");
    const auto parsed = parseGfa(below.compact);
    require(spell(parsed, parsed.walks.at("ref.chr1")) == below.reference,
            "compact reference walk is lossless");
    require(spell(parsed, parsed.walks.at("query.chr1")) == below.query,
            "compact query walk is lossless");

    const auto boundary = runCompactFixture(root, 25, 1, "threshold-25");
    require(boundary.compact_stats.suppressed_exact_runs == 0,
            "25 bp exact relation is retained");

    const auto parallel12 = runCompactFixture(
        root, 24, 12, "threshold-24-t12");
    const auto parallel32 = runCompactFixture(
        root, 24, 32, "threshold-24-t32");
    require(readFile(below.compact) == readFile(parallel12.compact) &&
            readFile(below.compact) == readFile(parallel32.compact),
            "compact output is byte-identical at 1, 12, and 32 threads");
}

void testUnitigAndCompoundAllele(const fs::path& root) {
    {
        const std::string sequence(200, 'A');
        const fs::path reference = root / "unitig-reference.fa";
        const fs::path query = root / "unitig-query.fa";
        writeFasta(reference, {{"chr1", sequence}});
        writeFasta(query, {{"chr1", sequence}});
        std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
        managers.emplace("ref", managerFor(reference));
        managers.emplace("query", managerFor(query));
        RaMesh::RaMeshMultiGenomeGraph graph(managers);
        graph.reference_order = {"ref", "query"};
        for (std::size_t start : {std::size_t{0}, std::size_t{100}}) {
            Cigar_t cigar;
            parseCigarString("100M", cigar);
            Anchor anchor(
                0, static_cast<Coord_t>(start), 100,
                0, static_cast<Coord_t>(start), 100,
                Strand::FORWARD, 100, 100, cigar);
            graph.insertAnchorIntoGraph(
                *managers.at("ref"), *managers.at("query"),
                "ref", "query", anchor, true);
        }
        RaMesh::Gfa::GfaExportOptions options;
        options.version = RaMesh::Gfa::Version::V1_1;
        options.profile = RaMesh::Gfa::Profile::COMPACT;
        options.work_dir = root / "unitig-work";
        options.threads = 2;
        const auto stats = graph.exportToGfa(
            root / "unitig.gfa", managers, options);
        require(stats.unitig_merges >= 1,
                "occurrence-aware unitig compaction merges a forced chain");
        require(stats.nodes == 1, "forced shared chain becomes one node");
    }

    {
        const std::string prefix(10000, 'A');
        const std::string middle(30, 'C');
        const std::string suffix(10000, 'G');
        const std::string reference_sequence =
            prefix + "A" + middle + "T" + suffix;
        const std::string query_sequence =
            prefix + "T" + middle + "A" + suffix;
        const fs::path reference = root / "allele-reference.fa";
        const fs::path query = root / "allele-query.fa";
        writeFasta(reference, {{"chr1", reference_sequence}});
        writeFasta(query, {{"chr1", query_sequence}});
        std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
        managers.emplace("ref", managerFor(reference));
        managers.emplace("query", managerFor(query));
        RaMesh::RaMeshMultiGenomeGraph graph(managers);
        graph.reference_order = {"ref", "query"};
        Cigar_t cigar;
        parseCigarString(std::to_string(reference_sequence.size()) + "M",
                         cigar);
        Anchor anchor(
            0, 0, static_cast<Length_t>(reference_sequence.size()),
            0, 0, static_cast<Length_t>(query_sequence.size()),
            Strand::FORWARD,
            static_cast<uint_t>(reference_sequence.size()),
            static_cast<uint_t>(query_sequence.size()), cigar);
        graph.insertAnchorIntoGraph(
            *managers.at("ref"), *managers.at("query"),
            "ref", "query", anchor, true);
        RaMesh::Gfa::GfaExportOptions options;
        options.version = RaMesh::Gfa::Version::V1_1;
        options.profile = RaMesh::Gfa::Profile::COMPACT;
        options.work_dir = root / "allele-work";
        options.threads = 4;
        const auto stats = graph.exportToGfa(
            root / "allele.gfa", managers, options);
        require(stats.allele_islands == 1,
                "dense small variants become an observed compound allele");
        const auto parsed = parseGfa(root / "allele.gfa");
        require(spell(parsed, parsed.walks.at("ref.chr1")) ==
                    reference_sequence,
                "compound reference allele is lossless");
        require(spell(parsed, parsed.walks.at("query.chr1")) ==
                    query_sequence,
                "compound query allele is lossless");

        options.compact.enable_compound_alleles = false;
        options.work_dir = root / "allele-disabled-work";
        const auto disabled_stats = graph.exportToGfa(
            root / "allele-disabled.gfa", managers, options);
        require(disabled_stats.allele_islands == 0,
                "compound allele rewriting can be disabled internally");
        const auto disabled = parseGfa(root / "allele-disabled.gfa");
        require(spell(disabled, disabled.walks.at("ref.chr1")) ==
                    reference_sequence &&
                spell(disabled, disabled.walks.at("query.chr1")) ==
                    query_sequence,
                "disabling compound alleles remains lossless");
    }
}

}  // namespace

int main(int argc, char** argv) {
    const bool preserve_outputs = argc == 2;
    const fs::path root = preserve_outputs
        ? fs::path(argv[1])
        : fs::current_path() / "gfa_export_test_tmp";
    fs::remove_all(root);
    fs::create_directories(root);
    try {
        require(RaMesh::Gfa::parseVersion("1.0") ==
                    RaMesh::Gfa::Version::V1_0,
                "parse GFA 1.0");
        require(RaMesh::Gfa::parseVersion("1.1") ==
                    RaMesh::Gfa::Version::V1_1,
                "parse GFA 1.1");
        expectFailure([] { RaMesh::Gfa::parseVersion("1"); });
        expectFailure([] { RaMesh::Gfa::parseVersion("1.2"); });
        require(RaMesh::Gfa::parseProfile("exact") ==
                    RaMesh::Gfa::Profile::EXACT,
                "parse exact GFA profile");
        require(RaMesh::Gfa::parseProfile("compact") ==
                    RaMesh::Gfa::Profile::COMPACT,
                "parse compact GFA profile");
        expectFailure([] { RaMesh::Gfa::parseProfile("compact-v1"); });
        testVersionsAndLosslessWalks(root);
        testPathNameCollision(root);
        testCompactThresholdAndDeterminism(root);
        testUnitigAndCompoundAllele(root);
        if (!preserve_outputs) fs::remove_all(root);
        return 0;
    } catch (...) {
        if (!preserve_outputs) fs::remove_all(root);
        throw;
    }
}
