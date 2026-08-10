#include "window_detector.h"

#include <algorithm>
#include <cassert>
#include <iostream>

namespace WD = RaMesh::WindowDetection;

namespace {

WD::SegmentSnapshot makeSegment(const std::string& block_id,
                                const std::string& species,
                                uint64_t start,
                                uint64_t length) {
    WD::SegmentSnapshot segment;
    segment.segment_id = block_id + "_" + species;
    segment.block_id = block_id;
    segment.species = species;
    segment.chromosome = "chr1";
    segment.graph_start = start;
    segment.graph_end = start + length;
    segment.original_start = start;
    segment.original_end = start + length;
    segment.original_coordinates_resolved = true;
    return segment;
}

WD::BlockSnapshot makeBlock(const std::string& block_id,
                            uint64_t start,
                            uint64_t length,
                            const std::vector<std::string>& species) {
    WD::BlockSnapshot block;
    block.block_id = block_id;
    block.canonical_signature = block_id;
    block.min_segment_length = length;
    block.max_segment_length = length;
    block.participating_genomes = species;
    for (const auto& name : species) {
        block.segments.push_back(makeSegment(block_id, name, start, length));
    }
    return block;
}

WD::GraphSnapshot makeMicroSnapshot(uint64_t right_start = 505) {
    WD::GraphSnapshot snapshot;
    snapshot.round_id = 0;
    snapshot.current_reference = "g1";
    snapshot.input_genomes = {"g1", "g2", "g3"};
    snapshot.blocks = {
        makeBlock("A", 0, 500, snapshot.input_genomes),
        makeBlock("M", 500, 5, snapshot.input_genomes),
        makeBlock("C", right_start, 500, snapshot.input_genomes),
    };
    for (const auto& species : snapshot.input_genomes) {
        WD::PathSnapshot path;
        path.species = species;
        path.chromosome = "chr1";
        path.segments = {
            makeSegment("A", species, 0, 500),
            makeSegment("M", species, 500, 5),
            makeSegment("C", species, right_start, 500),
        };
        snapshot.paths.push_back(std::move(path));
    }
    snapshot.active_segment_count = 9;
    return snapshot;
}

void testMicroWindow() {
    WD::Options options;
    auto result = WD::detectProblemWindows(makeMicroSnapshot(), options);
    assert(result.windows.size() == 1);
    const auto& window = result.windows.front();
    assert(window.priority_tier == "A");
    assert(window.has_two_sided_anchors);
    assert(window.max_possible_k == 3);
    assert(std::find(window.signals.begin(), window.signals.end(),
                     "MICRO_BLOCK") != window.signals.end());
    assert(result.boundaries.size() == 2);
    for (const auto& boundary : result.boundaries) {
        assert(boundary.detector_recommends_bridge);
    }
}

void testHardBoundaryProtection() {
    WD::Options options;
    auto snapshot = makeMicroSnapshot(2006);
    auto result = WD::detectProblemWindows(std::move(snapshot), options);
    assert(!result.windows.empty());
    bool saw_protected = false;
    bool saw_hard_boundary = false;
    for (const auto& window : result.windows) {
        if (window.priority_tier == "PROTECTED") {
            saw_protected = true;
        }
    }
    for (const auto& boundary : result.boundaries) {
        if (boundary.hard_boundary) {
            saw_hard_boundary = true;
            assert(!boundary.detector_recommends_bridge);
        }
    }
    assert(saw_protected);
    assert(saw_hard_boundary);
}

void testNToKCompatibilitySearch() {
    WD::Options options;
    auto snapshot = makeMicroSnapshot();
    for (auto& block : snapshot.blocks) {
        if (block.block_id != "C") {
            continue;
        }
        for (auto& segment : block.segments) {
            if (segment.species == "g3") {
                segment.reverse = true;
            }
        }
    }
    for (auto& path : snapshot.paths) {
        if (path.species == "g3") {
            path.segments.back().reverse = true;
        }
    }
    auto result = WD::detectProblemWindows(std::move(snapshot), options);
    assert(result.windows.size() == 1);
    const auto& window = result.windows.front();
    assert(window.max_possible_k == 2);
    assert(window.max_k_genomes == std::vector<std::string>({"g1", "g2"}));
    assert(window.excluded_genomes == std::vector<std::string>({"g3"}));
}

}  // namespace

int main() {
    testMicroWindow();
    testHardBoundaryProtection();
    testNToKCompatibilitySearch();
    std::cout << "window detector tests passed\n";
    return 0;
}
