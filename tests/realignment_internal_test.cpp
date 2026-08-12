#include "realignment_internal.h"

#include <stdexcept>
#include <string>

namespace {

using namespace RaMesh;
using namespace RaMesh::Realignment;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

SegPtr addSegment(const BlockPtr& block,
                  const SpeciesName& species,
                  uint_t start) {
    auto segment = Segment::create(
        start, 10, Strand::FORWARD,
        Cigar_t{cigarToInt('M', 10)}, AlignRole::PRIMARY,
        SegmentRole::SEGMENT, block);
    block->anchors.emplace(SpeciesChrPair{species, "chr1"}, segment);
    return segment;
}

void testBlockViewProfilesAndCache() {
    auto block = Block::createEmpty("chr1", 2);
    addSegment(block, "ref", 0);
    addSegment(block, "query", 10);

    BlockViewBuilder builder("ref", 1);
    BlockView view;
    require(builder.build(block, BlockViewProfile::ExactMerge, view),
            "exact profile rejected a valid block");
    require(view.species_count == 2 && view.reference_segment,
            "valid Block view lost participants or reference");
    require(builder.build(block, BlockViewProfile::MissingWindow, view),
            "missing-window profile rejected a valid block");
    require(builder.build(block, BlockViewProfile::Diagnostics, view, "chr1"),
            "diagnostic profile rejected a valid block");

    builder.clear(2);
    require(builder.build(block, BlockViewProfile::ExactMerge, view),
            "cache generation reset changed a valid Block view");

    auto singleton = Block::createEmpty("chr1", 1);
    addSegment(singleton, "ref", 0);
    require(!builder.build(singleton, BlockViewProfile::ExactMerge, view),
            "exact profile accepted a singleton block");
    require(builder.build(singleton, BlockViewProfile::MissingWindow, view),
            "missing-window profile rejected a singleton block");
}

void testPlannerConflictSelection() {
    auto a = Block::createEmpty("chr1", 0);
    auto b = Block::createEmpty("chr1", 0);
    auto c = Block::createEmpty("chr1", 0);
    const std::vector<PlannerConflictFootprint> footprints{
        {{a.get(), b.get()}, {b.get()}},
        {{b.get(), c.get()}, {c.get()}},
        {{c.get()}, {a.get()}}};
    const auto selected = MissingWindowPlanner::selectConflictFreeBatch(
        {0, 1, 2}, footprints);
    require(selected.size() == 1 && selected.front() == 0,
            "planner changed stable conflict precedence");

    LocalRepairCandidate candidate;
    candidate.reference_chromosome = "chr1";
    candidate.read_blocks = {a, b};
    candidate.replaced_blocks = {b};
    candidate.affected_paths = {{"ref", "chr1"}};
    PreparedGraphReplacement prepared;
    prepared.provisional_block = c;
    prepared.replaced_blocks = candidate.replaced_blocks;
    prepared.progress_before = 2;
    prepared.progress_after = 1;
    require(prepared.progress_after < prepared.progress_before,
            "prepared replacement does not encode strict progress");
}

}  // namespace

int main() {
    testBlockViewProfilesAndCache();
    testPlannerConflictSelection();
    return 0;
}
