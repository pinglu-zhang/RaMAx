#include "ramesh.h"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

namespace fs = std::filesystem;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void writeFasta(const fs::path& path, const std::string& name,
                const std::string& sequence) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("cannot create test FASTA");
    output << '>' << name << '\n' << sequence << '\n';
    if (!output) throw std::runtime_error("cannot write test FASTA");
}

SeqPro::SharedManagerVariant makeManager(const fs::path& fasta) {
    SeqPro::ManagerVariant manager{
        std::make_unique<SeqPro::SequenceManager>(fasta)};
    return std::make_shared<SeqPro::ManagerVariant>(std::move(manager));
}

RaMesh::BlockPtr makeBlock(
    const SpeciesName& block_reference_species,
    const SpeciesName& first_species,
    const SpeciesName& second_species,
    uint_t first_start,
    uint_t first_length,
    Strand first_strand,
    Cigar_t first_cigar,
    uint_t second_start,
    uint_t second_length,
    Strand second_strand,
    Cigar_t second_cigar) {
    auto block = RaMesh::Block::createEmpty(
        block_reference_species, "chr", 2);
    auto first = RaMesh::Segment::create(
        first_start, first_length, first_strand, std::move(first_cigar),
        RaMesh::AlignRole::PRIMARY, RaMesh::SegmentRole::SEGMENT, block);
    auto second = RaMesh::Segment::create(
        second_start, second_length, second_strand, std::move(second_cigar),
        RaMesh::AlignRole::PRIMARY, RaMesh::SegmentRole::SEGMENT, block);
    block->anchors[{first_species, "chr"}] = first;
    block->anchors[{second_species, "chr"}] = second;
    return block;
}

std::size_t lineCount(const fs::path& path) {
    std::ifstream input(path);
    std::size_t count = 0;
    std::string line;
    while (std::getline(input, line)) ++count;
    return count;
}

std::string readFile(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(input),
                       std::istreambuf_iterator<char>());
}

void testMalformedBlocksAndReverseRecord(const fs::path& directory) {
    writeFasta(directory / "ref.fa", "chr", "AGTCAAAA");
    writeFasta(directory / "rev.fa", "chr", "GACTAAAA");
    writeFasta(directory / "bad.fa", "chr", "AGTCAAAA");

    std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
    managers.emplace("ref", makeManager(directory / "ref.fa"));
    managers.emplace("rev", makeManager(directory / "rev.fa"));
    managers.emplace("bad", makeManager(directory / "bad.fa"));
    RaMesh::RaMeshMultiGenomeGraph graph(managers);

    std::vector<RaMesh::BlockPtr> storage;
    storage.push_back(makeBlock(
        "ref", "ref", "rev", 0, 4, Strand::FORWARD, {},
        0, 4, Strand::REVERSE, Cigar_t{cigarToInt('M', 4)}));
    storage.push_back(makeBlock(
        "ref", "ref", "bad", 0, 4, Strand::FORWARD, {},
        0, 4, Strand::FORWARD, Cigar_t{cigarToInt('M', 3)}));
    storage.push_back(makeBlock(
        "ghost", "ref", "bad", 0, 4, Strand::FORWARD, {},
        0, 4, Strand::FORWARD, Cigar_t{cigarToInt('M', 4)}));
    storage.push_back(makeBlock(
        "ref", "ref", "bad", 0, 4, Strand::FORWARD, {},
        7, 4, Strand::FORWARD, Cigar_t{cigarToInt('M', 4)}));
    for (const auto& block : storage) graph.blocks.push_back(block);

    const fs::path output = directory / "malformed.paf";
    const auto stats = graph.exportToPaf(output, managers);
    require(stats.blocks_seen == 4, "exporter Block count");
    require(stats.blocks_exported == 1, "only valid Block should export");
    require(stats.invalid_blocks_skipped == 3,
            "three malformed Blocks should be skipped");
    require(stats.records_written == 1, "one valid PAF record expected");
    require(stats.reverse_records == 1, "reverse record count");
    require(lineCount(output) == 1, "malformed Blocks must not be partial");

    const std::string record = readFile(output);
    require(record.find("rev.chr\t8\t0\t4\t-\tref.chr\t8\t0\t4\t4\t4\t255\t") == 0,
            "reverse PAF coordinates or strand are incorrect");
    require(record.find("cg:Z:4=\ttp:A:P\tNM:i:0\n") !=
                std::string::npos,
            "reverse PAF tags are incorrect");
}

void testQualifiedNameCollisionIsAtomic(const fs::path& directory) {
    writeFasta(directory / "collision-a.fa", "b.c", "AGTCAAAA");
    writeFasta(directory / "collision-ab.fa", "c", "AGTCAAAA");

    std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
    managers.emplace("a", makeManager(directory / "collision-a.fa"));
    managers.emplace("a.b", makeManager(directory / "collision-ab.fa"));
    RaMesh::RaMeshMultiGenomeGraph graph(managers);

    const fs::path output = directory / "atomic.paf";
    {
        std::ofstream existing(output, std::ios::binary | std::ios::trunc);
        existing << "sentinel-existing-output\n";
    }

    bool threw = false;
    try {
        (void)graph.exportToPaf(output, managers);
    } catch (const std::runtime_error& error) {
        threw = std::string(error.what()).find(
                    "Duplicate qualified PAF name") != std::string::npos;
    }
    require(threw, "qualified-name collision must fail globally");
    require(readFile(output) == "sentinel-existing-output\n",
            "global error must preserve existing output");
    for (const auto& entry : fs::directory_iterator(directory)) {
        require(entry.path().filename().string().find("atomic.paf.tmp.") != 0,
                "global error left a temporary PAF");
    }
}

void testMafAtomicPublish(const fs::path& directory) {
    writeFasta(directory / "maf-ref.fa", "chr", "AGTCAAAA");
    writeFasta(directory / "maf-query.fa", "chr", "AGTCAAAA");

    std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
    managers.emplace("ref", makeManager(directory / "maf-ref.fa"));
    managers.emplace("query", makeManager(directory / "maf-query.fa"));
    RaMesh::RaMeshMultiGenomeGraph graph(managers);

    auto block = makeBlock(
        "ref", "ref", "query", 0, 4, Strand::FORWARD, {},
        0, 4, Strand::FORWARD, Cigar_t{cigarToInt('M', 4)});
    graph.blocks.push_back(block);

    const fs::path output = directory / "atomic.maf";
    {
        std::ofstream existing(output, std::ios::binary | std::ios::trunc);
        existing << "sentinel-existing-output\n";
    }

    graph.exportToMaf(output, managers, false);
    const std::string maf = readFile(output);
    require(maf.starts_with("##maf version=1 scoring=none\n"),
            "MAF atomic publish did not replace existing output");
    require(maf.find("s ref.chr") != std::string::npos,
            "MAF output is missing reference row");
    require(maf.find("s query.chr") != std::string::npos,
            "MAF output is missing query row");
    for (const auto& entry : fs::directory_iterator(directory)) {
        require(entry.path().filename().string().find("atomic.maf.tmp.") != 0,
                "successful MAF export left a temporary file");
    }
}

}  // namespace

int main() {
    const auto stamp = std::chrono::steady_clock::now()
                           .time_since_epoch().count();
    const fs::path directory = fs::temp_directory_path() /
        ("ramax-paf-export-test-" + std::to_string(stamp));
    struct Cleanup {
        fs::path path;
        ~Cleanup() {
            std::error_code error;
            fs::remove_all(path, error);
        }
    } cleanup{directory};

    try {
        fs::create_directories(directory);
        testMalformedBlocksAndReverseRecord(directory);
        testQualifiedNameCollisionIsAtomic(directory);
        testMafAtomicPublish(directory);
    } catch (const std::exception& error) {
        std::cerr << "ramax_paf_export_tests: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    std::cout << "ramax_paf_export_tests: all checks passed\n";
    return EXIT_SUCCESS;
}
