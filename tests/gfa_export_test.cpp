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
        testVersionsAndLosslessWalks(root);
        testPathNameCollision(root);
        if (!preserve_outputs) fs::remove_all(root);
        return 0;
    } catch (...) {
        if (!preserve_outputs) fs::remove_all(root);
        throw;
    }
}
