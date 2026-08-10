#include "ramesh.h"

#include "SeqPro.h"
#include "align.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

using RaMesh::AlignRole;
using RaMesh::Block;
using RaMesh::BlockPtr;
using RaMesh::GenomeEnd;
using RaMesh::Segment;
using RaMesh::SegmentRole;
using RaMesh::SegPtr;
using RaMesh::SpeciesChrPair;

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void writeInput(const std::filesystem::path& path) {
    std::ofstream output(path);
    output << ">chr1\n" << std::string(10000, 'A') << '\n';
    if (!output) {
        throw std::runtime_error("cannot write zero-gap test FASTA");
    }
}

SegPtr addSegment(const BlockPtr& block,
                  const std::string& species,
                  uint_t start,
                  uint_t length,
                  Strand strand) {
    auto segment = Segment::create(
        start, length, strand,
        Cigar_t{cigarToInt('M', length)},
        AlignRole::PRIMARY, SegmentRole::SEGMENT, block);
    require(block->anchors.emplace(
                SpeciesChrPair{species, "chr1"}, segment).second,
            "duplicate test Block anchor");
    return segment;
}

void linkPath(GenomeEnd& genome_end,
              const std::vector<SegPtr>& segments) {
    auto previous = genome_end.head;
    for (const auto& segment : segments) {
        previous->primary_path.next.store(
            segment, std::memory_order_release);
        segment->primary_path.prev.store(
            previous, std::memory_order_release);
        previous = segment;
    }
    previous->primary_path.next.store(
        genome_end.tail, std::memory_order_release);
    genome_end.tail->primary_path.prev.store(
        previous, std::memory_order_release);
    genome_end.sample_vec.assign(1, genome_end.head);
    for (const auto& segment : segments) {
        genome_end.setToSampling(segment);
    }
}

void populateZeroGapWindow(RaMesh::RaMeshMultiGenomeGraph& graph) {
    auto left = Block::createEmpty("chr1", 4);
    auto middle = Block::createEmpty("chr1", 3);
    auto right = Block::createEmpty("chr1", 4);

    const auto chimp_left = addSegment(
        left, "simChimp", 0, 10, Strand::FORWARD);
    const auto chimp_middle = addSegment(
        middle, "simChimp", 10, 23, Strand::FORWARD);
    const auto chimp_right = addSegment(
        right, "simChimp", 33, 10, Strand::FORWARD);

    const auto gorilla_left = addSegment(
        left, "simGorilla", 133, 10, Strand::REVERSE);
    const auto gorilla_middle = addSegment(
        middle, "simGorilla", 110, 23, Strand::REVERSE);
    const auto gorilla_right = addSegment(
        right, "simGorilla", 100, 10, Strand::REVERSE);
    gorilla_left->cigar = Cigar_t{
        cigarToInt('M', 5), cigarToInt('I', 1),
        cigarToInt('D', 1), cigarToInt('M', 4)};
    gorilla_middle->cigar = Cigar_t{
        cigarToInt('M', 10), cigarToInt('D', 1),
        cigarToInt('I', 1), cigarToInt('M', 12)};
    gorilla_right->cigar = Cigar_t{
        cigarToInt('M', 4), cigarToInt('I', 1),
        cigarToInt('D', 1), cigarToInt('M', 5)};

    const auto human_left = addSegment(
        left, "simHuman", 200, 10, Strand::FORWARD);
    const auto human_middle = addSegment(
        middle, "simHuman", 210, 23, Strand::FORWARD);
    const auto human_right = addSegment(
        right, "simHuman", 233, 10, Strand::FORWARD);

    const auto orang_left = addSegment(
        left, "simOrang", 310, 10, Strand::REVERSE);
    const auto orang_right = addSegment(
        right, "simOrang", 300, 10, Strand::REVERSE);

    linkPath(graph.species_graphs.at("simChimp").chr2end.at("chr1"),
             {chimp_left, chimp_middle, chimp_right});
    linkPath(graph.species_graphs.at("simGorilla").chr2end.at("chr1"),
             {gorilla_right, gorilla_middle, gorilla_left});
    linkPath(graph.species_graphs.at("simHuman").chr2end.at("chr1"),
             {human_left, human_middle, human_right});
    linkPath(graph.species_graphs.at("simOrang").chr2end.at("chr1"),
             {orang_right, orang_left});

    graph.blocks = {left, middle, right};
}

void populateSubsetZeroGapWindow(
    RaMesh::RaMeshMultiGenomeGraph& graph,
    const std::vector<std::string>& middle_species,
    const std::string& nonzero_missing_species = "") {
    auto left = Block::createEmpty("chr1", 4);
    auto middle = Block::createEmpty("chr1", middle_species.size());
    auto right = Block::createEmpty("chr1", 4);
    const std::vector<std::string> species_names = {
        "simChimp", "simGorilla", "simHuman", "simOrang"};
    for (size_t species_index = 0;
         species_index < species_names.size(); ++species_index) {
        const auto& species = species_names[species_index];
        const uint_t base = static_cast<uint_t>(species_index * 100);
        const bool present =
            std::find(middle_species.begin(), middle_species.end(), species) !=
            middle_species.end();
        const auto left_segment = addSegment(
            left, species, base, 10, Strand::FORWARD);
        std::vector<SegPtr> path{left_segment};
        if (present) {
            path.push_back(addSegment(
                middle, species, base + 10, 20, Strand::FORWARD));
        }
        path.push_back(addSegment(
            right, species,
            present ? base + 30
                    : base + 10 +
                          (species == nonzero_missing_species ? 5 : 0),
            10, Strand::FORWARD));
        linkPath(graph.species_graphs.at(species).chr2end.at("chr1"), path);
    }
    graph.blocks = {left, middle, right};
}

void populateInteriorChain(
    RaMesh::RaMeshMultiGenomeGraph& graph,
    const std::vector<std::vector<std::string>>& participants,
    const std::string& reverse_species = "",
    uint_t interval_length = 30) {
    auto left = Block::createEmpty("chr1", 4);
    auto right = Block::createEmpty("chr1", 4);
    std::vector<BlockPtr> interiors;
    for (const auto& participant_set : participants) {
        interiors.push_back(
            Block::createEmpty("chr1", participant_set.size()));
    }
    const std::vector<std::string> species_names = {
        "simChimp", "simGorilla", "simHuman", "simOrang"};
    for (size_t species_index = 0;
         species_index < species_names.size(); ++species_index) {
        const auto& species = species_names[species_index];
        const uint_t base = static_cast<uint_t>(species_index * 100);
        const bool reverse = species == reverse_species;
        std::vector<SegPtr> reference_order;
        reference_order.push_back(addSegment(
            left, species,
            reverse ? base + 10 + interval_length : base,
            10, reverse ? Strand::REVERSE : Strand::FORWARD));
        for (size_t index = 0; index < participants.size(); ++index) {
            if (std::find(
                    participants[index].begin(), participants[index].end(),
                    species) == participants[index].end()) {
                continue;
            }
            const uint_t start =
                reverse ? base + interval_length -
                              static_cast<uint_t>(index * 10)
                        : base + 10 + static_cast<uint_t>(index * 10);
            reference_order.push_back(addSegment(
                interiors[index], species, start, 10,
                reverse ? Strand::REVERSE : Strand::FORWARD));
        }
        reference_order.push_back(addSegment(
            right, species,
            reverse ? base : base + 10 + interval_length,
            10, reverse ? Strand::REVERSE : Strand::FORWARD));
        if (reverse) {
            std::reverse(reference_order.begin(), reference_order.end());
        }
        linkPath(
            graph.species_graphs.at(species).chr2end.at("chr1"),
            reference_order);
    }
    graph.blocks = {left};
    graph.blocks.insert(
        graph.blocks.end(), interiors.begin(), interiors.end());
    graph.blocks.push_back(right);
}

void populateConflictingOrdinaryHybridWindows(
    RaMesh::RaMeshMultiGenomeGraph& graph) {
    auto a = Block::createEmpty("chr1", 4);
    auto b1 = Block::createEmpty("chr1", 3);
    auto b2 = Block::createEmpty("chr1", 3);
    auto c = Block::createEmpty("chr1", 4);
    auto d1 = Block::createEmpty("chr1", 3);
    auto d2 = Block::createEmpty("chr1", 3);
    auto e = Block::createEmpty("chr1", 4);
    const std::vector<std::string> species_names = {
        "simChimp", "simGorilla", "simHuman", "simOrang"};
    for (size_t species_index = 0;
         species_index < species_names.size(); ++species_index) {
        const auto& species = species_names[species_index];
        const uint_t base = static_cast<uint_t>(species_index * 200);
        std::vector<SegPtr> path;
        path.push_back(addSegment(
            a, species, base, 10, Strand::FORWARD));
        if (species != "simOrang") {
            path.push_back(addSegment(
                b1, species, base + 10, 10, Strand::FORWARD));
            path.push_back(addSegment(
                b2, species, base + 20, 10, Strand::FORWARD));
        }
        path.push_back(addSegment(
            c, species, base + 40, 10, Strand::FORWARD));
        if (species != "simHuman") {
            path.push_back(addSegment(
                d1, species, base + 50, 10, Strand::FORWARD));
            path.push_back(addSegment(
                d2, species, base + 60, 10, Strand::FORWARD));
        }
        path.push_back(addSegment(
            e, species,
            species == "simHuman" ? base + 50 : base + 80,
            10, Strand::FORWARD));
        linkPath(
            graph.species_graphs.at(species).chr2end.at("chr1"), path);
    }
    graph.blocks = {a, b1, b2, c, d1, d2, e};
}

void populateFiveBlockWindows(
    RaMesh::RaMeshMultiGenomeGraph& graph,
    bool invalidate_first_window,
    uint_t orang_second_gap) {
    auto a = Block::createEmpty("chr1", 4);
    auto b = Block::createEmpty("chr1", 3);
    auto c = Block::createEmpty("chr1", 4);
    auto d = Block::createEmpty("chr1", 3);
    auto e = Block::createEmpty("chr1", 4);

    std::map<std::string, std::vector<SegPtr>> paths;
    const std::map<std::string, uint_t> starts = {
        {"simChimp", 0},
        {"simGorilla", 100},
        {"simHuman", 200},
        {"simOrang", 300}};
    for (const auto& [species, start] : starts) {
        auto a_segment = addSegment(
            a, species, start, 10, Strand::FORWARD);
        paths[species].push_back(a_segment);
        if (species != "simOrang") {
            paths[species].push_back(addSegment(
                b, species, start + 10, 23, Strand::FORWARD));
        }
        paths[species].push_back(addSegment(
            c, species,
            species == "simOrang" ? start + 10 : start + 33,
            10, Strand::FORWARD));
        if (species != "simOrang") {
            paths[species].push_back(addSegment(
                d, species, start + 43, 17, Strand::FORWARD));
        }
        paths[species].push_back(addSegment(
            e, species,
            species == "simOrang"
                ? start + 20 + orang_second_gap
                : start + 60,
            10, Strand::FORWARD));

        if (species == "simGorilla" && invalidate_first_window) {
            a_segment->cigar = Cigar_t{cigarToInt('M', 9)};
        }
    }

    for (auto& [species, segments] : paths) {
        linkPath(
            graph.species_graphs.at(species).chr2end.at("chr1"),
            segments);
    }
    graph.blocks = {a, b, c, d, e};
}

void writePassthroughMsa(const std::filesystem::path& executable,
                         const std::filesystem::path& counter) {
    std::ofstream script(executable);
    script << "#!/bin/sh\n"
           << "printf x >> '" << counter.string() << "'\n"
           << "exec /bin/cat \"$5\"\n";
    if (!script) {
        throw std::runtime_error("cannot write passthrough MSA script");
    }
    script.close();
    std::filesystem::permissions(
        executable,
        std::filesystem::perms::owner_read |
            std::filesystem::perms::owner_write |
            std::filesystem::perms::owner_exec,
        std::filesystem::perm_options::replace);
}

void writeCapturingMsa(const std::filesystem::path& executable,
                       const std::filesystem::path& counter,
                       const std::filesystem::path& captured_input) {
    std::ofstream script(executable);
    script << "#!/bin/sh\n"
           << "printf x >> '" << counter.string() << "'\n"
           << "/bin/cp \"$5\" '" << captured_input.string() << "'\n"
           << "exec /bin/cat \"$5\"\n";
    if (!script) {
        throw std::runtime_error("cannot write capturing MSA script");
    }
    script.close();
    std::filesystem::permissions(
        executable,
        std::filesystem::perms::owner_read |
            std::filesystem::perms::owner_write |
            std::filesystem::perms::owner_exec,
        std::filesystem::perm_options::replace);
}

void writeFailOnceMsa(const std::filesystem::path& executable,
                      const std::filesystem::path& counter,
                      const std::filesystem::path& marker) {
    std::ofstream script(executable);
    script << "#!/bin/sh\n"
           << "printf x >> '" << counter.string() << "'\n"
           << "if [ ! -e '" << marker.string() << "' ]; then\n"
           << "  : > '" << marker.string() << "'\n"
           << "  exit 1\n"
           << "fi\n"
           << "exec /bin/cat \"$5\"\n";
    if (!script) {
        throw std::runtime_error("cannot write fail-once MSA script");
    }
    script.close();
    std::filesystem::permissions(
        executable,
        std::filesystem::perms::owner_read |
            std::filesystem::perms::owner_write |
            std::filesystem::perms::owner_exec,
        std::filesystem::perm_options::replace);
}

void writeReferenceInsertionMsa(const std::filesystem::path& executable) {
    std::ofstream script(executable);
    script << "#!/bin/sh\n"
           << "cat <<'EOF'\n"
           << ">s0\n"
           << "AAAAAAAAAAAAAAA-AAAAAAAAAAAAAAA\n"
           << ">s1\n"
           << "AAAAAAAAAAAAAAAA-AAAAAAAAAAAAAA\n"
           << ">s2\n"
           << "AAAAAAAAAAAAAAAA-AAAAAAAAAAAAAA\n"
           << "EOF\n";
    if (!script) {
        throw std::runtime_error("cannot write insertion MSA script");
    }
    script.close();
    std::filesystem::permissions(
        executable,
        std::filesystem::perms::owner_read |
            std::filesystem::perms::owner_write |
            std::filesystem::perms::owner_exec,
        std::filesystem::perm_options::replace);
}

BlockPtr onlyActiveBlock(RaMesh::RaMeshMultiGenomeGraph& graph) {
    BlockPtr active;
    for (const auto& weak_block : graph.blocks) {
        if (const auto block = weak_block.lock()) {
            require(!active, "more than one active Block remains");
            active = block;
        }
    }
    require(active != nullptr, "no active Block remains");
    return active;
}

std::string graphSignature(RaMesh::RaMeshMultiGenomeGraph& graph) {
    std::string signature;
    auto segment = graph.species_graphs.at("simChimp")
                       .chr2end.at("chr1")
                       .head->primary_path.next.load(
                           std::memory_order_acquire);
    while (segment && !segment->isTail()) {
        const auto block = segment->parent_block;
        require(block != nullptr, "signature path contains detached Block");
        signature += "B" + std::to_string(block->anchors.size()) + ":";
        for (const auto& [key, anchor] : block->anchors) {
            signature += key.first + "." + key.second + "@" +
                         std::to_string(anchor->start) + "+" +
                         std::to_string(anchor->length) + "/" +
                         std::to_string(static_cast<int>(anchor->strand)) +
                         "/" + cigarToString(anchor->cigar) + ";";
        }
        signature += "|";
        segment = segment->primary_path.next.load(
            std::memory_order_acquire);
    }
    return signature;
}

}  // namespace

int main() {
    const auto temp = std::filesystem::path("/tmp") /
        ("ramax-zero-gap-merge-" + std::to_string(getpid()));
    std::filesystem::remove_all(temp);
    std::filesystem::create_directories(temp);

    try {
        std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
        for (const std::string species :
             {"simChimp", "simGorilla", "simHuman", "simOrang"}) {
            const auto fasta = temp / (species + ".fa");
            writeInput(fasta);
            SeqPro::ManagerVariant manager =
                std::make_unique<SeqPro::SequenceManager>(fasta);
            managers[species] =
                std::make_shared<SeqPro::ManagerVariant>(
                    std::move(manager));
        }

        RaMesh::RaMeshMultiGenomeGraph graph(managers);
        populateZeroGapWindow(graph);
        const size_t replaced =
            graph.realignSingleMissingSpeciesWindows(
                "simChimp", managers, "/not-used-for-zero-gap",
                3000, 1, 200);
        require(replaced == 1, "zero-gap 4-3-4 window was not merged");
        require(graph.blocks.size() == 1,
                "three old Blocks were not replaced by one Block");

        const auto merged = onlyActiveBlock(graph);
        require(merged->anchors.size() == 4,
                "merged Block does not contain all four species");
        const auto orang = merged->anchors.at({"simOrang", "chr1"});
        require(orang->start == 300 && orang->length == 20 &&
                    orang->strand == Strand::REVERSE,
                "reverse-strand missing-species coordinates changed");
        require(cigarToString(orang->cigar) == "10M23D10M",
                "missing species did not receive the expected deletion CIGAR");
        require(countRefLength(orang->cigar) == 43 &&
                    countQryLength(orang->cigar) == 20,
                "missing-species CIGAR consumption is inconsistent");
        const auto gorilla = merged->anchors.at({"simGorilla", "chr1"});
        const auto gorilla_cigar = cigarToString(gorilla->cigar);
        require(gorilla->strand == Strand::REVERSE &&
                    gorilla_cigar.find('I') != std::string::npos &&
                    gorilla_cigar.find('D') != std::string::npos &&
                    countRefLength(gorilla->cigar) == 43 &&
                    countQryLength(gorilla->cigar) == 43,
                "complex reverse-strand CIGAR was not preserved");
        require(graph.verifyGraphCorrectness(false),
                "graph verification failed after zero-gap merge");

        RaMesh::RaMeshMultiGenomeGraph two_missing_graph(managers);
        populateSubsetZeroGapWindow(
            two_missing_graph, {"simChimp", "simGorilla"});
        require(two_missing_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers, "/not-used-for-zero-gap",
                    3000, 2, 200) == 1,
                "zero-gap 4-2-4 window was not merged");
        const auto two_missing_block = onlyActiveBlock(two_missing_graph);
        require(cigarToString(two_missing_block->anchors.at(
                    {"simHuman", "chr1"})->cigar) == "10M20D10M" &&
                    cigarToString(two_missing_block->anchors.at(
                    {"simOrang", "chr1"})->cigar) == "10M20D10M",
                "4-2-4 missing species did not receive deletion CIGARs");
        require(two_missing_graph.verifyGraphCorrectness(false),
                "graph verification failed after 4-2-4 merge");

        RaMesh::RaMeshMultiGenomeGraph three_missing_graph(managers);
        populateSubsetZeroGapWindow(three_missing_graph, {"simChimp"});
        require(three_missing_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers, "/not-used-for-zero-gap",
                    3000, 3, 200) == 1,
                "zero-gap 4-1-4 window was not merged");
        const auto three_missing_block = onlyActiveBlock(three_missing_graph);
        for (const std::string species :
             {"simGorilla", "simHuman", "simOrang"}) {
            require(cigarToString(three_missing_block->anchors.at(
                        {species, "chr1"})->cigar) == "10M20D10M",
                    "4-1-4 missing species deletion CIGAR is incorrect");
        }
        require(three_missing_graph.verifyGraphCorrectness(false),
                "graph verification failed after 4-1-4 merge");

        RaMesh::RaMeshMultiGenomeGraph mixed_gap_graph(managers);
        populateSubsetZeroGapWindow(
            mixed_gap_graph, {"simChimp", "simGorilla"}, "simHuman");
        require(mixed_gap_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers, "/bin/false",
                    3000, 2, 200) == 0,
                "mixed zero/nonzero missing intervals were accepted");
        require(mixed_gap_graph.blocks.size() == 3,
                "mixed-gap rejection modified the Block pool");

        RaMesh::RaMeshMultiGenomeGraph rejected_graph(managers);
        populateZeroGapWindow(rejected_graph);
        require(rejected_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers, "/not-used-for-zero-gap",
                    3000, 1, 22) == 0,
                "zero-gap window exceeded the configured threshold");
        require(rejected_graph.blocks.size() == 3,
                "threshold rejection modified the original Block pool");

        RaMesh::RaMeshMultiGenomeGraph fallback_graph(managers);
        populateFiveBlockWindows(fallback_graph, true, 0);
        require(fallback_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers, "/not-used-for-zero-gap",
                    3000, 1, 200) == 1,
                "a failed zero-gap preparation hid the following window");
        require(fallback_graph.blocks.size() == 3,
                "fallback selection did not preserve A/B and merge C/D/E");
        const auto fallback_orang =
            fallback_graph.species_graphs.at("simOrang")
                .chr2end.at("chr1").head->primary_path.next.load(
                    std::memory_order_acquire)
                ->primary_path.next.load(std::memory_order_acquire);
        require(fallback_orang && !fallback_orang->isTail() &&
                    cigarToString(fallback_orang->cigar) == "10M17D10M",
                "the valid following zero-gap window was not committed");
        require(fallback_graph.verifyGraphCorrectness(false),
                "fallback graph verification failed after later-window merge");

        RaMesh::RaMeshMultiGenomeGraph cascade_graph(managers);
        populateFiveBlockWindows(cascade_graph, false, 0);
        require(cascade_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers, "/not-used-for-zero-gap",
                    3000, 4, 200) == 2,
                "cascading zero-gap windows did not reach a fixed point");
        require(cascade_graph.blocks.size() == 1,
                "fixed-point zero-gap merge did not reduce five Blocks to one");
        require(cascade_graph.verifyGraphCorrectness(false),
                "graph verification failed after cascading zero-gap merge");

        const auto passthrough_msa = temp / "passthrough-minipoa.sh";
        const auto msa_counter = temp / "passthrough-minipoa.calls";
        writePassthroughMsa(passthrough_msa, msa_counter);

        const std::vector<std::vector<std::string>> pattern_4324 = {
            {"simChimp", "simGorilla", "simHuman"},
            {"simChimp", "simGorilla"}};
        const auto multiblock_msa = temp / "multiblock-minipoa.sh";
        const auto multiblock_counter = temp / "multiblock-minipoa.calls";
        writePassthroughMsa(multiblock_msa, multiblock_counter);
        RaMesh::RaMeshMultiGenomeGraph multiblock_graph(managers);
        populateInteriorChain(multiblock_graph, pattern_4324);
        require(multiblock_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers, multiblock_msa.string(),
                    3000, 4, 200) == 1,
                "4-3-2-4 window was not realigned");
        require(multiblock_graph.blocks.size() == 3,
                "4-3-2-4 interior Blocks were not replaced together");
        require(multiblock_graph.verifyGraphCorrectness(false),
                "graph verification failed after 4-3-2-4 realignment");

        RaMesh::RaMeshMultiGenomeGraph span_rejected_graph(managers);
        populateInteriorChain(
            span_rejected_graph, pattern_4324, "", 3001);
        require(span_rejected_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers, "/bin/false",
                    3000, 4, 200) == 0,
                "multi-Block window exceeding 3 kb was realigned");

        const auto capturing_msa = temp / "capturing-minipoa.sh";
        const auto capturing_counter = temp / "capturing-minipoa.calls";
        const auto captured_input = temp / "captured-minipoa.input.fa";
        writeCapturingMsa(
            capturing_msa, capturing_counter, captured_input);
        const std::vector<std::vector<std::string>> pattern_42324 = {
            {"simChimp", "simGorilla"},
            {"simChimp", "simHuman", "simOrang"},
            {"simChimp", "simGorilla"}};

        const auto hybrid_capturing_msa =
            temp / "hybrid-capturing-minipoa.sh";
        const auto hybrid_capturing_counter =
            temp / "hybrid-capturing-minipoa.calls";
        const auto hybrid_captured_input =
            temp / "hybrid-captured-minipoa.input.fa";
        writeCapturingMsa(
            hybrid_capturing_msa, hybrid_capturing_counter,
            hybrid_captured_input);
        const std::vector<std::vector<std::string>> pattern_4334 = {
            {"simChimp", "simGorilla", "simHuman"},
            {"simChimp", "simGorilla", "simHuman"}};
        RaMesh::RaMeshMultiGenomeGraph hybrid_graph(managers);
        populateInteriorChain(hybrid_graph, pattern_4334);
        const auto hybrid_right = hybrid_graph.blocks.back().lock();
        require(hybrid_right != nullptr,
                "hybrid right boundary is unavailable");
        hybrid_right->anchors.at({"simOrang", "chr1"})->start = 310;
        require(hybrid_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers,
                    hybrid_capturing_msa.string(),
                    3000, 4, 200) == 1,
                "hybrid empty-sequence window was not realigned");
        std::ifstream hybrid_captured(hybrid_captured_input);
        const std::string hybrid_captured_text(
            (std::istreambuf_iterator<char>(hybrid_captured)),
            std::istreambuf_iterator<char>());
        require(hybrid_captured_text.find(">simOrang") ==
                    std::string::npos,
                "empty species was written to the minipoa FASTA");
        const auto hybrid_block = onlyActiveBlock(hybrid_graph);
        const auto hybrid_orang = hybrid_block->anchors.at(
            {"simOrang", "chr1"});
        require(hybrid_orang->start == 300 &&
                    hybrid_orang->length == 20 &&
                    cigarToString(hybrid_orang->cigar) == "10M30D10M",
                "hybrid empty species did not receive the expected gaps");
        for (const auto& [key, segment] : hybrid_block->anchors) {
            (void)key;
            require(segment && segment->length != 0,
                    "hybrid preparation created a zero-length Segment");
        }
        require(hybrid_graph.verifyGraphCorrectness(false),
                "hybrid empty-sequence graph verification failed");

        RaMesh::RaMeshMultiGenomeGraph serial_hybrid_graph(managers);
        populateInteriorChain(serial_hybrid_graph, pattern_4334);
        const auto serial_hybrid_right =
            serial_hybrid_graph.blocks.back().lock();
        require(serial_hybrid_right != nullptr,
                "serial hybrid right boundary is unavailable");
        serial_hybrid_right->anchors.at({"simOrang", "chr1"})->start = 310;
        require(serial_hybrid_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers,
                    hybrid_capturing_msa.string(),
                    3000, 1, 200) == 1,
                "serial hybrid window was not realigned");
        require(graphSignature(serial_hybrid_graph) ==
                    graphSignature(hybrid_graph),
                "serial and parallel hybrid graph signatures diverged");

        const auto insertion_msa = temp / "reference-insertion-minipoa.sh";
        writeReferenceInsertionMsa(insertion_msa);
        RaMesh::RaMeshMultiGenomeGraph insertion_graph(managers);
        populateInteriorChain(insertion_graph, pattern_4334);
        const auto insertion_right = insertion_graph.blocks.back().lock();
        require(insertion_right != nullptr,
                "reference-insertion right boundary is unavailable");
        insertion_right->anchors.at({"simOrang", "chr1"})->start = 310;
        require(insertion_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers, insertion_msa.string(),
                    3000, 4, 200) == 1,
                "hybrid reference-insertion window was not realigned");
        const auto insertion_maf = temp / "hybrid-reference-insertion.maf";
        insertion_graph.exportToMaf(
            insertion_maf, managers, true, false);
        std::ifstream insertion_maf_input(insertion_maf);
        const std::string insertion_maf_text(
            (std::istreambuf_iterator<char>(insertion_maf_input)),
            std::istreambuf_iterator<char>());
        require(insertion_maf_text.find(
                    "AAAAAAAAAA-------------------------------AAAAAAAAAA") !=
                    std::string::npos,
                "empty species did not remain gap-only across an MSA "
                "reference-insertion column");

        const std::vector<std::vector<std::string>>
            reverse_empty_pattern = {
                {"simChimp", "simHuman", "simOrang"},
                {"simChimp", "simHuman", "simOrang"}};
        RaMesh::RaMeshMultiGenomeGraph reverse_empty_graph(managers);
        populateInteriorChain(
            reverse_empty_graph, reverse_empty_pattern,
            "simGorilla");
        const auto reverse_empty_left =
            reverse_empty_graph.blocks.front().lock();
        require(reverse_empty_left != nullptr,
                "reverse hybrid left boundary is unavailable");
        reverse_empty_left->anchors.at(
            {"simGorilla", "chr1"})->start = 110;
        require(reverse_empty_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers,
                    hybrid_capturing_msa.string(),
                    3000, 4, 200) == 1,
                "reverse empty-sequence hybrid window was not realigned");
        const auto reverse_empty_block =
            onlyActiveBlock(reverse_empty_graph);
        const auto reverse_empty_gorilla =
            reverse_empty_block->anchors.at(
                {"simGorilla", "chr1"});
        require(reverse_empty_gorilla->start == 100 &&
                    reverse_empty_gorilla->length == 20 &&
                    reverse_empty_gorilla->strand == Strand::REVERSE &&
                    cigarToString(reverse_empty_gorilla->cigar) ==
                        "10M30D10M",
                "reverse empty species coordinates or CIGAR changed");
        require(reverse_empty_graph.verifyGraphCorrectness(false),
                "reverse hybrid graph verification failed");

        const std::vector<std::vector<std::string>> pattern_4114 = {
            {"simChimp"}, {"simChimp"}};
        RaMesh::RaMeshMultiGenomeGraph reference_only_graph(managers);
        populateInteriorChain(reference_only_graph, pattern_4114);
        const auto reference_only_right =
            reference_only_graph.blocks.back().lock();
        require(reference_only_right != nullptr,
                "reference-only right boundary is unavailable");
        reference_only_right->anchors.at({"simGorilla", "chr1"})->start =
            110;
        reference_only_right->anchors.at({"simHuman", "chr1"})->start =
            210;
        reference_only_right->anchors.at({"simOrang", "chr1"})->start =
            310;
        require(reference_only_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers, "/bin/false",
                    3000, 4, 200) == 1,
                "reference-only hybrid window did not bypass minipoa");
        const auto reference_only_block =
            onlyActiveBlock(reference_only_graph);
        for (const std::string species :
             {"simGorilla", "simHuman", "simOrang"}) {
            require(cigarToString(reference_only_block->anchors.at(
                        {species, "chr1"})->cigar) == "10M30D10M",
                    "reference-only empty species CIGAR is incorrect");
        }
        require(reference_only_graph.verifyGraphCorrectness(false),
                "reference-only hybrid graph verification failed");

        const auto fail_once_msa = temp / "fail-once-minipoa.sh";
        const auto fail_once_counter = temp / "fail-once-minipoa.calls";
        const auto fail_once_marker = temp / "fail-once-minipoa.marker";
        writeFailOnceMsa(
            fail_once_msa, fail_once_counter, fail_once_marker);
        RaMesh::RaMeshMultiGenomeGraph conflict_graph(managers);
        populateConflictingOrdinaryHybridWindows(conflict_graph);
        require(conflict_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers, fail_once_msa.string(),
                    3000, 4, 200) == 2,
                "unified fallback did not recover the conflicting windows");
        std::ifstream fail_once_counter_input(fail_once_counter);
        std::string fail_once_calls;
        fail_once_counter_input >> fail_once_calls;
        require(fail_once_calls == "xxx",
                "ordinary failure did not trigger same-snapshot hybrid "
                "fallback followed by fixed-point ordinary retry");
        require(conflict_graph.blocks.size() == 3 &&
                    conflict_graph.verifyGraphCorrectness(false),
                "unified ordinary/hybrid fixed-point graph is invalid");

        RaMesh::RaMeshMultiGenomeGraph reverse_multiblock_graph(managers);
        populateInteriorChain(
            reverse_multiblock_graph, pattern_42324, "simGorilla");
        require(reverse_multiblock_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers, capturing_msa.string(),
                    3000, 4, 200) == 1,
                "negative-strand 4-2-3-2-4 window was not realigned");
        std::ifstream captured(captured_input);
        const std::string captured_text(
            (std::istreambuf_iterator<char>(captured)),
            std::istreambuf_iterator<char>());
        require(captured_text.find(std::string(30, 'T')) !=
                    std::string::npos,
                "negative-strand sequence was not reverse-complemented");
        auto replacement = reverse_multiblock_graph.species_graphs
                               .at("simChimp").chr2end.at("chr1")
                               .head->primary_path.next.load(
                                   std::memory_order_acquire)
                               ->primary_path.next.load(
                                   std::memory_order_acquire)
                               ->parent_block;
        const auto reverse_anchor = replacement->anchors.at(
            {"simGorilla", "chr1"});
        require(reverse_anchor->strand == Strand::REVERSE &&
                    reverse_anchor->start == 110 &&
                    reverse_anchor->length == 30,
                "negative-strand replacement coordinates changed");
        require(reverse_multiblock_graph.blocks.size() == 3 &&
                    reverse_multiblock_graph.verifyGraphCorrectness(false),
                "negative-strand multiblock graph is invalid");

        RaMesh::RaMeshMultiGenomeGraph strand_mismatch_graph(managers);
        populateInteriorChain(
            strand_mismatch_graph, pattern_42324, "simGorilla");
        const auto mismatch_right = strand_mismatch_graph.blocks.back().lock();
        require(mismatch_right != nullptr,
                "strand mismatch right boundary is unavailable");
        mismatch_right->anchors.at({"simGorilla", "chr1"})->strand =
            Strand::FORWARD;
        require(strand_mismatch_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers, "/bin/false",
                    3000, 4, 200) == 0,
                "window crossing a strand transition was accepted");

        RaMesh::RaMeshMultiGenomeGraph order_mismatch_graph(managers);
        populateInteriorChain(
            order_mismatch_graph, pattern_42324);
        const auto order_mismatch_right =
            order_mismatch_graph.blocks.back().lock();
        require(order_mismatch_right != nullptr,
                "anchor-order mismatch right boundary is unavailable");
        order_mismatch_right->anchors.at({"simOrang", "chr1"})->start =
            100;
        require(order_mismatch_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers, "/bin/false",
                    3000, 4, 200) == 0,
                "backward anchor-order window was accepted");
        require(order_mismatch_graph.blocks.size() == 5,
                "backward anchor-order rejection modified the Block pool");

        RaMesh::RaMeshMultiGenomeGraph staged_graph(managers);
        populateFiveBlockWindows(staged_graph, false, 17);
        require(staged_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers, passthrough_msa.string(),
                    3000, 4, 200) == 2,
                "minipoa rescan did not run after zero-gap commit");
        std::ifstream counter_input(msa_counter);
        std::string counter_contents;
        counter_input >> counter_contents;
        require(counter_contents == "x",
                "staged workflow did not invoke minipoa exactly once");
        require(staged_graph.blocks.size() == 3,
                "ordinary minipoa replacement unexpectedly changed Block count");
        require(staged_graph.verifyGraphCorrectness(false),
                "graph verification failed after staged minipoa commit");

        RaMesh::RaMeshMultiGenomeGraph serial_staged_graph(managers);
        populateFiveBlockWindows(serial_staged_graph, false, 17);
        require(serial_staged_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers, passthrough_msa.string(),
                    3000, 1, 200) == 2,
                "serial staged workflow changed the replacement count");
        require(graphSignature(serial_staged_graph) ==
                    graphSignature(staged_graph),
                "serial and parallel staged workflows diverged");
        std::ifstream second_counter_input(msa_counter);
        counter_contents.clear();
        second_counter_input >> counter_contents;
        require(counter_contents == "xx",
                "serial and parallel workflows did not each invoke one MSA");

        RaMesh::RaMeshMultiGenomeGraph retained_zero_graph(managers);
        populateFiveBlockWindows(retained_zero_graph, false, 17);
        require(retained_zero_graph.realignSingleMissingSpeciesWindows(
                    "simChimp", managers, "/bin/false",
                    3000, 4, 200) == 1,
                "minipoa failure did not preserve committed zero-gap work");
        require(retained_zero_graph.blocks.size() == 3,
                "failed minipoa phase changed the committed zero-gap graph");
        require(retained_zero_graph.verifyGraphCorrectness(false),
                "graph verification failed after retained zero-gap commit");

        std::filesystem::remove_all(temp);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "zero_gap_missing_species_merge_test: "
                  << error.what() << '\n';
        std::filesystem::remove_all(temp);
        return 1;
    }
}
