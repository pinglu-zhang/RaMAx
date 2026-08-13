#include "short_block_repair.h"

#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

#include "align.h"
#include "ramesh.h"

namespace {

using namespace RaMesh;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void writeFasta(const std::filesystem::path& path) {
    std::ofstream output(path);
    output << ">chr1\n" << std::string(1000, 'A') << '\n';
    require(static_cast<bool>(output), "cannot write short-Block FASTA");
}

SegPtr add(const BlockPtr& block,
           const SpeciesChrPair& key,
           uint_t start,
           uint_t length) {
    auto segment = Segment::create(
        start, length, Strand::FORWARD,
        Cigar_t{cigarToInt('M', length)}, AlignRole::PRIMARY,
        SegmentRole::SEGMENT, block);
    require(block->anchors.emplace(key, segment).second,
            "duplicate short-Block test anchor");
    return segment;
}

void link(GenomeEnd& end, const std::vector<SegPtr>& segments) {
    auto previous = end.head;
    for (const auto& segment : segments) {
        previous->primary_path.next.store(segment);
        segment->primary_path.prev.store(previous);
        previous = segment;
    }
    previous->primary_path.next.store(end.tail);
    end.tail->primary_path.prev.store(previous);
    end.sample_vec.assign(1, end.head);
    for (const auto& segment : segments) end.setToSampling(segment);
}

size_t activeBlocks(const RaMeshMultiGenomeGraph& graph) {
    size_t result = 0;
    for (const auto& weak : graph.blocks) result += !weak.expired();
    return result;
}

struct Fixture {
    BlockPtr left;
    BlockPtr small;
    BlockPtr right;
    SegPtr ref_left;
    SegPtr ref_small;
    SegPtr ref_right;
    SegPtr q1_left;
    SegPtr q1_small;
    SegPtr q1_right;
};

Fixture baseFixture(RaMeshMultiGenomeGraph& graph) {
    Fixture fixture;
    fixture.left = Block::createEmpty("chr1", 2);
    fixture.small = Block::createEmpty("chr1", 2);
    fixture.right = Block::createEmpty("chr1", 2);
    fixture.ref_left = add(fixture.left, {"ref", "chr1"}, 0, 100);
    fixture.ref_small = add(fixture.small, {"ref", "chr1"}, 100, 20);
    fixture.ref_right = add(fixture.right, {"ref", "chr1"}, 120, 80);
    fixture.q1_left = add(fixture.left, {"q1", "chr1"}, 0, 100);
    fixture.q1_small = add(fixture.small, {"q1", "chr1"}, 100, 20);
    fixture.q1_right = add(fixture.right, {"q1", "chr1"}, 120, 80);
    link(graph.species_graphs.at("ref").chr2end.at("chr1"),
         {fixture.ref_left, fixture.ref_small, fixture.ref_right});
    link(graph.species_graphs.at("q1").chr2end.at("chr1"),
         {fixture.q1_left, fixture.q1_small, fixture.q1_right});
    graph.blocks = {fixture.left, fixture.small, fixture.right};
    return fixture;
}

}  // namespace

int main() {
    const auto temporary = std::filesystem::path("/tmp") /
        ("ramax-short-block-" + std::to_string(getpid()));
    std::filesystem::remove_all(temporary);
    std::filesystem::create_directories(temporary);
    try {
        std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
        for (const auto& species : {std::string("ref"), std::string("q1"),
                                    std::string("q2")}) {
            const auto fasta = temporary / (species + ".fa");
            writeFasta(fasta);
            auto manager = std::make_unique<SeqPro::SequenceManager>(fasta);
            managers.emplace(
                species,
                std::make_shared<SeqPro::ManagerVariant>(std::move(manager)));
        }

        {
            RaMeshMultiGenomeGraph graph(managers);
            const auto fixture = baseFixture(graph);
            ShortBlockRepair::Options options;
            const auto result = ShortBlockRepair::repairFinalShortBlocks(
                graph, {"ref"}, managers, options);
            require(result.left_merged == 0 && result.deleted_blocks == 0 &&
                        fixture.ref_small->parent_block == fixture.small,
                    "disabled short-Block repair changed graph");
        }

        {
            RaMeshMultiGenomeGraph graph(managers);
            const auto fixture = baseFixture(graph);
            ShortBlockRepair::Options options;
            options.enabled = true;
            const auto result = ShortBlockRepair::repairFinalShortBlocks(
                graph, {"ref"}, managers, options);
            require(result.left_merged == 1 && result.ksw2_calls == 0 &&
                        result.deleted_blocks == 0 && activeBlocks(graph) == 2,
                    "zero-KSW left merge failed");
            const auto merged = graph.species_graphs.at("ref")
                                    .chr2end.at("chr1")
                                    .head->primary_path.next.load();
            require(merged && merged->length == 120 &&
                        merged->parent_block->anchors.size() == 2,
                    "zero-KSW merged Block is invalid");
        }

        {
            RaMeshMultiGenomeGraph graph(managers);
            auto fixture = baseFixture(graph);
            const auto q2_left = add(fixture.left, {"q2", "chr1"}, 0, 100);
            const auto q2_right = add(fixture.right, {"q2", "chr1"}, 120, 80);
            link(graph.species_graphs.at("q2").chr2end.at("chr1"),
                 {q2_left, q2_right});
            ShortBlockRepair::Options options;
            options.enabled = true;
            const auto result = ShortBlockRepair::repairFinalShortBlocks(
                graph, {"ref"}, managers, options);
            require(result.left_merged == 1 && result.ksw2_calls >= 1 &&
                        result.ksw2_passed >= 1 &&
                        result.deleted_blocks == 0,
                    "missing species was not filled by KSW2");
            const auto q2_merged = graph.species_graphs.at("q2")
                                       .chr2end.at("chr1")
                                       .head->primary_path.next.load();
            require(q2_merged && q2_merged->length == 120 &&
                        q2_merged->primary_path.next.load()->length == 80,
                    "KSW2-filled query interval was not inserted");
        }

        {
            RaMeshMultiGenomeGraph graph(managers);
            BlockPtr left = Block::createEmpty("chr1", 1);
            BlockPtr first = Block::createEmpty("chr1", 2);
            BlockPtr second = Block::createEmpty("chr1", 2);
            BlockPtr right = Block::createEmpty("chr1", 1);
            const auto ref_left = add(left, {"ref", "chr1"}, 0, 100);
            const auto ref_first = add(first, {"ref", "chr1"}, 100, 20);
            const auto ref_second = add(second, {"ref", "chr1"}, 120, 20);
            const auto ref_right = add(right, {"ref", "chr1"}, 140, 100);
            const auto q_first = add(first, {"q1", "chr1"}, 100, 20);
            const auto q_second = add(second, {"q1", "chr1"}, 120, 20);
            link(graph.species_graphs.at("ref").chr2end.at("chr1"),
                 {ref_left, ref_first, ref_second, ref_right});
            link(graph.species_graphs.at("q1").chr2end.at("chr1"),
                 {q_first, q_second});
            graph.blocks = {left, first, second, right};
            ShortBlockRepair::Options options;
            options.enabled = true;
            const auto result = ShortBlockRepair::repairFinalShortBlocks(
                graph, {"ref"}, managers, options);
            require(result.right_merged == 1 &&
                        result.deleted_blocks == 1 &&
                        result.deleted_segments == 2 && activeBlocks(graph) == 2,
                    "consecutive failed short Blocks were not deleted");
            require(ref_left->primary_path.next.load() == ref_right &&
                        ref_right->primary_path.prev.load() == ref_left &&
                        graph.species_graphs.at("q1").chr2end.at("chr1")
                                .head->primary_path.next.load()
                            ->isTail(),
                    "consecutive deletion corrupted paths");
            require(ref_first->parent_block == nullptr &&
                        ref_second->parent_block == nullptr &&
                        q_first->parent_block == nullptr &&
                        q_second->parent_block == nullptr,
                    "original consecutive short Segments survived deletion");
        }
    } catch (...) {
        std::filesystem::remove_all(temporary);
        throw;
    }
    std::filesystem::remove_all(temporary);
    return 0;
}
