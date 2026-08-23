#include "hal/export.h"

#include <cstdint>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using RaMesh::hal_export::BinaryInferenceResult;
using RaMesh::hal_export::BucketedDonor;
using RaMesh::hal_export::EdgeSupport;
using RaMesh::hal_export::ExportBlockOrderKey;
using RaMesh::hal_export::ExportStats;
using RaMesh::hal_export::GenomeSequenceName;
using RaMesh::hal_export::LeafInterval;
using RaMesh::hal_export::RunOrderKey;
using RaMesh::hal_export::TreeMeta;

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void expect(bool condition, const std::string& message) {
    if (!condition) {
        fail(message);
    }
}

uint8_t stateByName(const TreeMeta& tree,
                    const BinaryInferenceResult& result,
                    const std::string& node_name) {
    auto node_it = tree.name_to_id.find(node_name);
    if (node_it == tree.name_to_id.end()) {
        fail("missing node in TreeMeta: " + node_name);
    }
    return result.present_by_node.at(static_cast<size_t>(node_it->second));
}

TreeMeta parseTree(const std::string& newick) {
    NewickParser parser(newick);
    return RaMesh::hal_export::buildTreeMeta(parser);
}

void testTreeMetaBuild() {
    TreeMeta tree = parseTree("((A:1,B:1)AB:1,C:1)root;");
    expect(tree.root_id != -1, "root_id should be valid");
    expect(tree.nodes.at(static_cast<size_t>(tree.root_id)).name == "root", "root name mismatch");
    expect(tree.leaf_ids.size() == 3, "expected exactly three leaves");
    expect(tree.internal_postorder.size() == 2, "expected exactly two internal nodes");
    expect(tree.name_to_id.contains("AB"), "missing AB internal node");
}

void testSubtreeSharedInference() {
    TreeMeta tree = parseTree("((A:1,B:1)AB:1,C:0.05)root;");
    std::unordered_map<std::string, bool> leaf_presence{
        {"A", true},
        {"B", true},
        {"C", false},
    };

    BinaryInferenceResult result = RaMesh::hal_export::inferDescendantUnion(tree, leaf_presence);
    expect(stateByName(tree, result, "AB") == 1, "AB should be present for subtree-shared run");
    expect(stateByName(tree, result, "root") == 1, "root must contain the descendant union");
}

void testSingleLeafInsertionInference() {
    TreeMeta tree = parseTree("((A:1,B:0.05)AB:0.05,C:0.05)root;");
    std::unordered_map<std::string, bool> leaf_presence{
        {"A", true},
        {"B", false},
        {"C", false},
    };

    BinaryInferenceResult result = RaMesh::hal_export::inferDescendantUnion(tree, leaf_presence);
    expect(stateByName(tree, result, "AB") == 1, "AB must retain a run present in any descendant");
    expect(stateByName(tree, result, "root") == 1, "root must retain a run present in any descendant");
}
void testDescendantUnionIgnoresBranchLength() {
    std::unordered_map<std::string, bool> leaf_presence{
        {"A", true},
        {"B", false},
        {"C", false},
    };

    TreeMeta long_b_tree = parseTree("((A:1,B:10)AB:10,C:0.05)root;");
    TreeMeta short_b_tree = parseTree("((A:1,B:0.05)AB:10,C:0.05)root;");
    BinaryInferenceResult long_b =
        RaMesh::hal_export::inferDescendantUnion(long_b_tree, leaf_presence);
    BinaryInferenceResult short_b =
        RaMesh::hal_export::inferDescendantUnion(short_b_tree, leaf_presence);

    expect(long_b.present_by_node == short_b.present_by_node,
           "descendant-union segment presence must not depend on branch-length tuning");
}

void testLeafIntervalProjection() {
    LeafInterval forward = RaMesh::hal_export::projectLeafInterval(100, 20, false, 3, 5);
    expect(forward.start == 103, "forward projection start mismatch");
    expect(forward.length == 5, "forward projection length mismatch");
    expect(forward.forward_to_parent, "forward projection orientation mismatch");

    LeafInterval reverse = RaMesh::hal_export::projectLeafInterval(100, 20, true, 3, 5);
    expect(reverse.start == 112, "reverse projection start mismatch");
    expect(reverse.length == 5, "reverse projection length mismatch");
    expect(!reverse.forward_to_parent, "reverse projection orientation mismatch");
}

void testElementaryRunsPreserveAllOccurrences() {
    using RaMesh::hal_export::AlignedOccurrence;
    std::vector<AlignedOccurrence> rows{
        {"A.chr1", "A", "chr1", 100, 4, true, "AC-GT"},
        {"A.chr2", "A", "chr2", 200, 5, false, "ACCGT"},
        {"B.chr1", "B", "chr1", 300, 4, false, "A-CGT"},
    };

    auto runs = RaMesh::hal_export::projectElementaryRuns(rows);
    expect(runs.size() == 4, "participant-mask changes should create four elementary runs");
    expect(runs[0].col_beg == 0 && runs[0].col_end == 1,
           "first elementary run column range mismatch");
    expect(runs[0].occurrences.size() == 3,
           "same-genome occurrences on distinct sequences must both survive projection");

    size_t genome_a_count = 0;
    for (const auto& occurrence : runs[0].occurrences) {
        if (occurrence.genome_name == "A") {
            ++genome_a_count;
        }
        if (occurrence.row_id == "A.chr1") {
            expect(occurrence.start == 103,
                   "reverse occurrence must project from the segment's right edge");
            expect(occurrence.reversed, "reverse occurrence orientation was lost");
        }
    }
    expect(genome_a_count == 2, "occurrence projection collapsed rows by genome name");
}

void testElementaryRunsRejectDuplicateOccurrenceIds() {
    using RaMesh::hal_export::AlignedOccurrence;
    std::vector<AlignedOccurrence> rows{
        {"A.chr1", "A", "chr1", 0, 2, false, "AC"},
        {"A.chr1", "A", "chr2", 0, 2, false, "AC"},
    };
    bool threw = false;
    try {
        (void)RaMesh::hal_export::projectElementaryRuns(rows);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    expect(threw, "duplicate occurrence row ids must fail instead of overwriting a row");
}

void testC2HBottomFlagEncoding() {
    expect(RaMesh::hal_export::c2hHasBottomFlag(0) == 0, "top-only sequence should encode as 0");
    expect(RaMesh::hal_export::c2hHasBottomFlag(1) == 1, "bottom-bearing sequence should encode as 1");
    expect(RaMesh::hal_export::c2hHasBottomFlag(34) == 1,
           "multiple bottom segments must still encode as boolean 1");
}

void testCactusZScoreWeightsLengthAndGap() {
    expect(
        std::abs(
            RaMesh::hal_export::calculateCactusZScore(
                20, 30, 7, 0.0L) -
            600.0L) <
            1e-12L,
        "theta zero must reduce the Cactus Z-score to the segment length product");
    const long double adjacent =
        RaMesh::hal_export::calculateCactusZScore(
            100, 100, 0, 0.01L);
    const long double separated =
        RaMesh::hal_export::calculateCactusZScore(
            100, 100, 100, 0.01L);
    expect(
        adjacent > separated && separated > 0.0L,
        "Cactus Z-score must decay with observed gap length");
}

void testOutputSequenceOrderMatchesWriters() {
    std::vector<std::string> genome_order{"root", "simOrang", "simHuman"};
    std::vector<GenomeSequenceName> unordered{
        {"simOrang", "simOrang.simOrang.chrE"},
        {"root", "root.scf000002"},
        {"simHuman", "simHuman.simHuman.chrA"},
        {"simOrang", "simOrang.simOrang.chrB"},
        {"root", "root.scf000001"},
        {"simOrang", "simOrang.simOrang.chrD"},
    };

    auto ordered = RaMesh::hal_export::buildOutputSequenceOrder(genome_order, unordered);
    std::vector<GenomeSequenceName> expected{
        {"root", "root.scf000001"},
        {"root", "root.scf000002"},
        {"simOrang", "simOrang.simOrang.chrB"},
        {"simOrang", "simOrang.simOrang.chrD"},
        {"simOrang", "simOrang.simOrang.chrE"},
        {"simHuman", "simHuman.simHuman.chrA"},
    };
    expect(ordered == expected,
           "output sequence ordering should follow genome_order first, then seq_name lexicographically");
}

void testInternalChildOrientationEmission() {
    expect(!RaMesh::hal_export::computeForwardToParent(false, true),
           "reverse child vs forward parent should emit reverse top orientation");
}

void testLeafOrientationRespectsParentPlacement() {
    expect(RaMesh::hal_export::computeForwardToParent(false, false),
           "reverse leaf vs reverse parent should still be forward_to_parent");
    expect(!RaMesh::hal_export::computeForwardToParent(true, false),
           "forward leaf vs reverse parent should emit reverse top orientation");
}

void testReversePlacementChangesAncestorDNA() {
    expect(RaMesh::hal_export::orientRunDNAForPlacement("AACG", true) == "AACG",
           "forward placement should keep canonical DNA");
    expect(RaMesh::hal_export::orientRunDNAForPlacement("AACG", false) == "CGTT",
           "reverse placement should reverse-complement canonical DNA");
}

void testMixedOrientationPreservesAdjacency() {
    auto mixed = RaMesh::hal_export::orientAdjacencyVote(10, true, 20, false);
    expect(mixed.has_value(), "mixed orientation adjacency must be retained");
    expect(mixed->left_run_id == 10 &&
               mixed->left_forward_to_canonical &&
               mixed->right_run_id == 20 &&
               !mixed->right_forward_to_canonical,
           "mixed orientation vote must preserve both segment orientations");

    std::vector<uint64_t> run_ids{10, 20};
    std::vector<EdgeSupport> edges{
        EdgeSupport{10, 20, 1, {1}, true, false},
    };
    std::unordered_map<uint64_t, RunOrderKey> run_keys{
        {10, RunOrderKey{1, 10}},
        {20, RunOrderKey{1, 20}},
    };
    auto paths = RaMesh::hal_export::buildMaximumCardinalityWeightPathCover(
        run_ids, edges, run_keys, nullptr);
    expect(paths == std::vector<std::vector<uint64_t>>({{10, 20}}),
           "mixed orientation extremities must remain in one ancestral path");

    auto self = RaMesh::hal_export::orientAdjacencyVote(10, true, 10, false);
    expect(!self.has_value(), "collapsed paralogous self-adjacency must remain explicit");
}

void testMaximumWeightChoosesSupportedJoin() {
    std::vector<uint64_t> run_ids{1, 2, 3};
    std::vector<EdgeSupport> edges{
        EdgeSupport{1, 2, 3, {10, 11}},
        EdgeSupport{1, 3, 2, {10, 11}},
    };
    std::unordered_map<uint64_t, RunOrderKey> run_keys{
        {1, RunOrderKey{1, 10}},
        {2, RunOrderKey{1, 20}},
        {3, RunOrderKey{1, 30}},
    };
    auto paths = RaMesh::hal_export::buildMaximumCardinalityWeightPathCover(
        run_ids, edges, run_keys, nullptr);
    expect(paths.size() == 2, "one outgoing endpoint may select exactly one supported join");
    expect(paths[0] == std::vector<uint64_t>({1, 2}),
           "maximum-weight matching should select the stronger observed adjacency");
}

void testMaximumWeightMatchingIsGlobal() {
    std::vector<uint64_t> run_ids{1, 2, 3, 4};
    std::vector<EdgeSupport> edges{
        EdgeSupport{1, 3, 100, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}},
        EdgeSupport{1, 4, 0, {1, 2, 3, 4, 5, 6, 7, 8, 9}},
        EdgeSupport{2, 3, 0, {1, 2, 3, 4, 5, 6, 7, 8, 9}},
    };
    std::unordered_map<uint64_t, RunOrderKey> run_keys{
        {1, RunOrderKey{1, 10}},
        {2, RunOrderKey{1, 20}},
        {3, RunOrderKey{1, 30}},
        {4, RunOrderKey{1, 40}},
    };
    auto paths = RaMesh::hal_export::buildMaximumCardinalityWeightPathCover(
        run_ids, edges, run_keys, nullptr);
    expect(paths == std::vector<std::vector<uint64_t>>({{1, 4}, {2, 3}}),
           "matching must maximize total child support instead of selecting the strongest edge greedily");
}

void testEqualCardinalityAndWeightTieBreakIsDeterministic() {
    std::vector<uint64_t> run_ids{1, 2, 3, 4};
    std::vector<EdgeSupport> edges{
        EdgeSupport{1, 3, 1, {1}},
        EdgeSupport{1, 4, 1, {1}},
        EdgeSupport{2, 3, 1, {2}},
        EdgeSupport{2, 4, 1, {2}},
    };
    std::unordered_map<uint64_t, RunOrderKey> run_keys{
        {1, RunOrderKey{1, 10}},
        {2, RunOrderKey{1, 20}},
        {3, RunOrderKey{1, 30}},
        {4, RunOrderKey{1, 40}},
    };
    auto paths = RaMesh::hal_export::buildMaximumCardinalityWeightPathCover(
        run_ids, edges, run_keys, nullptr);
    expect(
        paths ==
            std::vector<std::vector<uint64_t>>(
                {{1, 3}, {2, 4}}),
        "equal-cardinality and equal-weight matchings must use a stable endpoint-order tie break");
}

void testSupportedJoinCardinalityPrecedesWeight() {
    std::vector<uint64_t> run_ids{1, 2, 3, 4};
    std::vector<EdgeSupport> edges{
        EdgeSupport{1, 3, 100, {1, 2}},
        EdgeSupport{1, 4, 1, {1}},
        EdgeSupport{2, 3, 1, {2}},
    };
    std::unordered_map<uint64_t, RunOrderKey> run_keys{
        {1, RunOrderKey{1, 10}},
        {2, RunOrderKey{1, 20}},
        {3, RunOrderKey{1, 30}},
        {4, RunOrderKey{1, 40}},
    };
    auto paths = RaMesh::hal_export::buildMaximumCardinalityWeightPathCover(
        run_ids, edges, run_keys, nullptr);
    expect(
        paths ==
            std::vector<std::vector<uint64_t>>(
                {{1, 4}, {2, 3}}),
        "two endpoint-compatible observed joins must precede one stronger join");
}

void testObservedEdgeIsRetained() {
    std::vector<uint64_t> run_ids{1, 2};
    std::vector<EdgeSupport> edges{
        EdgeSupport{1, 2, 2, {10, 11}},
    };
    std::unordered_map<uint64_t, RunOrderKey> run_keys{
        {1, RunOrderKey{1, 10}},
        {2, RunOrderKey{1, 20}},
    };
    auto paths = RaMesh::hal_export::buildMaximumCardinalityWeightPathCover(
        run_ids, edges, run_keys, nullptr);
    expect(paths == std::vector<std::vector<uint64_t>>({{1, 2}}),
           "an observed adjacency must be retained when its endpoints have no competing evidence");
}

void testDisconnectedComponentsStaySeparate() {
    std::vector<uint64_t> run_ids{1, 2, 3, 4};
    std::vector<EdgeSupport> edges{
        EdgeSupport{1, 2, 2, {10, 11}},
        EdgeSupport{3, 4, 2, {12, 13}},
    };
    std::unordered_map<uint64_t, RunOrderKey> run_keys{
        {1, RunOrderKey{1, 10}},
        {2, RunOrderKey{1, 20}},
        {3, RunOrderKey{2, 10}},
        {4, RunOrderKey{2, 20}},
    };
    auto paths = RaMesh::hal_export::buildMaximumCardinalityWeightPathCover(run_ids, edges, run_keys, nullptr);
    expect(paths.size() == 2, "disconnected components should remain as two paths");
    expect(paths[0] == std::vector<uint64_t>({1, 2}), "first component path mismatch");
    expect(paths[1] == std::vector<uint64_t>({3, 4}), "second component path mismatch");
}

void testReferencePackingScaffoldsWithinEvidenceComponent() {
    std::vector<uint64_t> occurrence_ids{1, 2, 3};
    std::vector<EdgeSupport> edges{
        EdgeSupport{1, 2, 2, {10, 11}, true, true, 10.0L, 0},
        EdgeSupport{1, 3, 1, {10}, true, true, 1.0L, 0},
    };
    std::unordered_map<uint64_t, RunOrderKey> occurrence_keys{
        {1, RunOrderKey{1, 10}},
        {2, RunOrderKey{1, 20}},
        {3, RunOrderKey{1, 30}},
    };
    RaMesh::hal_export::ExportStats stats;
    auto assembly = RaMesh::hal_export::buildAncestralSequenceAssembly(
        occurrence_ids,
        edges,
        occurrence_keys,
        {},
        10,
        &stats);
    expect(
        assembly.sequences.size() == 1,
        "unmatched paths in one evidence component must be scaffolded together");
    std::vector<uint64_t> observed;
    size_t direct_joins = 0;
    size_t scaffold_joins = 0;
    for (const auto& sequence : assembly.sequences) {
        for (const auto& fragment :
             sequence.supported_fragments) {
            observed.insert(
                observed.end(),
                fragment.begin(),
                fragment.end());
        }
        for (const auto& join : sequence.joins) {
            if (join.kind ==
                RaMesh::hal_export::ReferenceJoinKind::DIRECT) {
                expect(
                    join.gap_length == 0,
                    "direct evidence must remain zero gap");
                ++direct_joins;
            } else {
                expect(
                    join.kind ==
                            RaMesh::hal_export::ReferenceJoinKind::SCAFFOLD &&
                        join.gap_length == 10,
                    "component scaffold must use one explicit 10N join");
                ++scaffold_joins;
            }
        }
    }
    std::sort(observed.begin(), observed.end());
    expect(
        observed == occurrence_ids &&
            direct_joins == 1 &&
            scaffold_joins == 1 &&
            stats.path_vertex_count == 3 &&
            stats.supported_join_count == 1 &&
            stats.direct_join_count == 1 &&
            stats.indirect_join_count == 0 &&
            stats.scaffold_join_count == 1 &&
            stats.scaffold_gap_bases == 10,
        "component-scaffold reference packing statistics mismatch");
}
void testDisconnectedEvidenceComponentsShareReferenceInterval() {
    const std::vector<uint64_t> occurrence_ids{1, 2, 3, 4, 5, 6};
    const std::vector<EdgeSupport> edges{
        EdgeSupport{1, 2, 10, {10}, true, true, 10.0L, 0},
        EdgeSupport{1, 3, 9, {10}, true, true, 9.0L, 0},
        EdgeSupport{4, 5, 10, {11}, true, true, 10.0L, 0},
        EdgeSupport{4, 6, 9, {11}, true, true, 9.0L, 0},
    };
    const std::unordered_map<uint64_t, RunOrderKey>
        occurrence_keys{
            {1, RunOrderKey{1, 10}},
            {2, RunOrderKey{1, 20}},
            {3, RunOrderKey{1, 30}},
            {4, RunOrderKey{2, 10}},
            {5, RunOrderKey{2, 20}},
            {6, RunOrderKey{2, 30}},
        };
    RaMesh::hal_export::ExportStats stats;
    const auto assembly =
        RaMesh::hal_export::buildAncestralSequenceAssembly(
            occurrence_ids,
            edges,
            occurrence_keys,
            {},
            10,
            &stats);
    expect(
        assembly.sequences.size() == 1 &&
            stats.supported_join_count == 2 &&
            stats.direct_join_count == 2 &&
            stats.scaffold_join_count == 3 &&
            stats.scaffold_gap_bases == 30,
        "top-level stub matching must pack disconnected evidence components into one reference interval");
    std::set<uint64_t> occurrences;
    for (const auto& fragment :
         assembly.sequences.front().supported_fragments) {
        occurrences.insert(
            fragment.begin(),
            fragment.end());
    }
    expect(
        occurrences ==
            std::set<uint64_t>({1, 2, 3, 4, 5, 6}),
        "reference scaffolding must preserve every occurrence exactly once");
}
void testDuplicateHomologyOccurrencesRemainDistinct() {
    std::vector<uint64_t> occurrence_ids{11, 12, 21, 22};
    std::vector<EdgeSupport> edges{
        EdgeSupport{
            11, 21, 2, {10, 11}, true, true, 10.0L, 0},
        EdgeSupport{
            12, 22, 2, {10, 11}, true, true, 10.0L, 0},
    };
    std::unordered_map<uint64_t, RunOrderKey> occurrence_keys{
        {11, RunOrderKey{1, 10}},
        {12, RunOrderKey{1, 10}},
        {21, RunOrderKey{2, 20}},
        {22, RunOrderKey{2, 20}},
    };
    auto assembly =
        RaMesh::hal_export::buildAncestralSequenceAssembly(
            occurrence_ids,
            edges,
            occurrence_keys,
            {},
            10,
            nullptr);
    std::vector<uint64_t> observed;
    for (const auto& sequence : assembly.sequences) {
        for (const auto& fragment :
             sequence.supported_fragments) {
            observed.insert(
                observed.end(),
                fragment.begin(),
                fragment.end());
        }
    }
    std::sort(observed.begin(), observed.end());
    expect(
        observed == occurrence_ids &&
            assembly.forward_by_occurrence.size() ==
                occurrence_ids.size(),
        "homologous copies sharing one block order key must retain distinct occurrence identities");
}


void testConfirmedTerminalEndsPreserveReferenceIntervals() {
    std::vector<uint64_t> occurrence_ids{1, 2, 3, 4};
    std::vector<EdgeSupport> edges{
        EdgeSupport{
            1, 2, 2, {10, 11}, true, true, 10.0L, 0},
        EdgeSupport{
            3, 4, 2, {12, 13}, true, true, 10.0L, 0},
    };
    std::unordered_map<uint64_t, RunOrderKey> occurrence_keys{
        {1, RunOrderKey{1, 10}},
        {2, RunOrderKey{1, 20}},
        {3, RunOrderKey{2, 10}},
        {4, RunOrderKey{2, 20}},
    };
    std::vector<RaMesh::hal_export::TerminalEndSupport> terminal_ends{
        {{1, RaMesh::hal_export::OccurrenceEndSide::LEFT},
         2,
         {10, 11},
         10.0L},
        {{2, RaMesh::hal_export::OccurrenceEndSide::RIGHT},
         2,
         {10, 11},
         10.0L},
        {{3, RaMesh::hal_export::OccurrenceEndSide::LEFT},
         2,
         {12, 13},
         10.0L},
        {{4, RaMesh::hal_export::OccurrenceEndSide::RIGHT},
         2,
         {12, 13},
         10.0L},
    };
    RaMesh::hal_export::ExportStats stats;
    auto assembly =
        RaMesh::hal_export::buildAncestralSequenceAssembly(
            occurrence_ids,
            edges,
            occurrence_keys,
            terminal_ends,
            10,
            &stats);
    expect(
        assembly.sequences.size() == 2,
        "confirmed terminal ends must preserve two supported reference paths");
    expect(
        assembly.sequences[0].supported_fragments ==
                std::vector<std::vector<uint64_t>>(
                    {{1}, {2}}) &&
            assembly.sequences[1].supported_fragments ==
                std::vector<std::vector<uint64_t>>(
                    {{3}, {4}}),
        "terminal-constrained reference path ordering mismatch");
    expect(
        stats.terminal_end_candidate_count == 4 &&
            stats.confirmed_terminal_end_count == 4 &&
            stats.reference_interval_count == 2 &&
            stats.supported_join_count == 2,
        "terminal-end and reference interval statistics mismatch");
}

void testSingleLineageTerminalEvidenceDoesNotCreateBoundaries() {
    const std::vector<uint64_t> occurrence_ids{1, 2, 3, 4};
    const std::unordered_map<uint64_t, RunOrderKey>
        occurrence_keys{
            {1, RunOrderKey{1, 10}},
            {2, RunOrderKey{1, 20}},
            {3, RunOrderKey{2, 10}},
            {4, RunOrderKey{2, 20}},
        };
    const std::vector<RaMesh::hal_export::TerminalEndSupport>
        terminal_ends{
            {{1, RaMesh::hal_export::OccurrenceEndSide::LEFT},
             1,
             {10},
             1.0L},
            {{2, RaMesh::hal_export::OccurrenceEndSide::RIGHT},
             1,
             {10},
             1.0L},
            {{3, RaMesh::hal_export::OccurrenceEndSide::LEFT},
             1,
             {11},
             1.0L},
            {{4, RaMesh::hal_export::OccurrenceEndSide::RIGHT},
             1,
             {11},
             1.0L},
        };
    RaMesh::hal_export::ExportStats stats;
    const auto assembly =
        RaMesh::hal_export::buildAncestralSequenceAssembly(
            occurrence_ids,
            {},
            occurrence_keys,
            terminal_ends,
            10,
            &stats);
    std::vector<uint64_t> observed;
    size_t scaffold_joins = 0;
    for (const auto& sequence : assembly.sequences) {
        for (const auto& join : sequence.joins) {
            expect(
                join.kind ==
                        RaMesh::hal_export::ReferenceJoinKind::SCAFFOLD &&
                    join.gap_length == 10,
                "zero-support reference packing must use explicit 10N scaffold joins");
            ++scaffold_joins;
        }
        for (const auto& fragment :
             sequence.supported_fragments) {
            observed.insert(
                observed.end(),
                fragment.begin(),
                fragment.end());
        }
    }
    std::sort(observed.begin(), observed.end());
    expect(
        assembly.sequences.size() == 1 &&
            observed == occurrence_ids &&
            scaffold_joins == 3 &&
            stats.terminal_end_candidate_count == 4 &&
            stats.confirmed_terminal_end_count == 0 &&
            stats.supported_join_count == 0 &&
            stats.scaffold_join_count == 3 &&
            stats.scaffold_gap_bases == 30 &&
            stats.reference_interval_count == 1,
        "single-lineage terminal evidence must not create ancestral reference boundaries");
}

void testConfirmedTerminalEndRemainsExternalDuringScaffolding() {
    const std::vector<uint64_t> occurrence_ids{1, 2, 3};
    const std::vector<EdgeSupport> edges{
        EdgeSupport{
            1, 2, 10, {10, 11}, true, true, 10.0L, 0},
        EdgeSupport{
            1, 3, 9, {10}, true, true, 9.0L, 0},
    };
    const std::unordered_map<uint64_t, RunOrderKey>
        occurrence_keys{
            {1, RunOrderKey{1, 10}},
            {2, RunOrderKey{1, 20}},
            {3, RunOrderKey{1, 30}},
        };
    const std::vector<RaMesh::hal_export::TerminalEndSupport>
        terminal_ends{
            {{3, RaMesh::hal_export::OccurrenceEndSide::RIGHT},
             2,
             {10, 11},
             2.0L},
        };
    RaMesh::hal_export::ExportStats stats;
    const auto assembly =
        RaMesh::hal_export::buildAncestralSequenceAssembly(
            occurrence_ids,
            edges,
            occurrence_keys,
            terminal_ends,
            10,
            &stats);
    expect(
        assembly.sequences.size() == 1 &&
            stats.confirmed_terminal_end_count == 1 &&
            stats.supported_join_count == 1 &&
            stats.scaffold_join_count == 1,
        "one confirmed terminal must allow component scaffolding through the opposite end");
    const auto& fragments =
        assembly.sequences.front().supported_fragments;
    const bool terminal_is_external =
        (fragments.front().front() == 3 &&
         !assembly.forward_by_occurrence.at(3)) ||
        (fragments.back().back() == 3 &&
         assembly.forward_by_occurrence.at(3));
    expect(
        terminal_is_external,
        "confirmed physical right end must remain an external sequence end");
}
void testMinimumScaffoldCoverPairsOppositeTerminals() {
    const std::vector<uint64_t> occurrence_ids{1, 2, 3, 4};
    const std::vector<EdgeSupport> edges{
        EdgeSupport{
            1, 2, 10, {10, 11}, true, true, 10.0L, 0},
        EdgeSupport{
            1, 3, 9, {10}, true, true, 9.0L, 0},
        EdgeSupport{
            1, 4, 8, {11}, true, false, 8.0L, 0},
    };
    const std::unordered_map<uint64_t, RunOrderKey>
        occurrence_keys{
            {1, RunOrderKey{1, 10}},
            {2, RunOrderKey{1, 20}},
            {3, RunOrderKey{1, 30}},
            {4, RunOrderKey{1, 40}},
        };
    const std::vector<RaMesh::hal_export::TerminalEndSupport>
        terminal_ends{
            {{3, RaMesh::hal_export::OccurrenceEndSide::RIGHT},
             2,
             {10, 11},
             2.0L},
            {{4, RaMesh::hal_export::OccurrenceEndSide::LEFT},
             2,
             {10, 11},
             2.0L},
        };
    RaMesh::hal_export::ExportStats stats;
    const auto assembly =
        RaMesh::hal_export::buildAncestralSequenceAssembly(
            occurrence_ids,
            edges,
            occurrence_keys,
            terminal_ends,
            10,
            &stats);
    expect(
        assembly.sequences.size() == 1 &&
            stats.confirmed_terminal_end_count == 2 &&
            stats.supported_join_count == 1 &&
            stats.scaffold_join_count == 2 &&
            stats.reference_interval_count == 1,
        "opposite confirmed terminals must bound one minimum scaffold path");
    const auto& fragments =
        assembly.sequences.front().supported_fragments;
    expect(
        fragments.front().front() == 3 &&
            fragments.back().back() == 4 &&
            !assembly.forward_by_occurrence.at(3) &&
            !assembly.forward_by_occurrence.at(4),
        "minimum scaffold path must preserve both physical terminal ends externally");
}

void testSingleLineageTerminalDoesNotOverrideSupportedJoin() {
    std::vector<uint64_t> occurrence_ids{1, 2};
    std::vector<EdgeSupport> edges{
        EdgeSupport{
            1, 2, 2, {10, 11}, true, true, 10.0L, 0},
    };
    std::unordered_map<uint64_t, RunOrderKey> occurrence_keys{
        {1, RunOrderKey{1, 10}},
        {2, RunOrderKey{1, 20}},
    };
    std::vector<RaMesh::hal_export::TerminalEndSupport>
        terminal_ends{
            {{1, RaMesh::hal_export::OccurrenceEndSide::LEFT},
             1,
             {10},
             1.0L},
            {{1, RaMesh::hal_export::OccurrenceEndSide::RIGHT},
             1,
             {10},
             1.0L},
        };
    RaMesh::hal_export::ExportStats stats;
    auto assembly =
        RaMesh::hal_export::buildAncestralSequenceAssembly(
            occurrence_ids,
            edges,
            occurrence_keys,
            terminal_ends,
            10,
            &stats);
    expect(
        assembly.sequences.size() == 1 &&
            assembly.sequences[0].supported_fragments ==
                std::vector<std::vector<uint64_t>>(
                    {{1}, {2}}),
        "one-child terminal evidence must not split an observed direct chain");
    expect(
        assembly.sequences[0].joins.size() == 1 &&
            stats.path_vertex_count == 2 &&
            stats.terminal_end_candidate_count == 2 &&
            stats.confirmed_terminal_end_count == 0 &&
            stats.supported_join_count == 1 &&
            stats.direct_join_count == 1 &&
            stats.indirect_join_count == 0,
        "one-child terminal evidence must not change supported reference threading");
}

void testSingleChildOvervoteCannotOverrideMultiChildSupport() {
    std::vector<uint64_t> run_ids{1, 2, 3};
    std::vector<EdgeSupport> edges{
        EdgeSupport{1, 2, 2, {10, 11}},
        EdgeSupport{1, 3, 5, {10}},
    };
    std::unordered_map<uint64_t, RunOrderKey> run_keys{
        {1, RunOrderKey{1, 10}},
        {2, RunOrderKey{1, 20}},
        {3, RunOrderKey{1, 30}},
    };
    auto paths = RaMesh::hal_export::buildMaximumCardinalityWeightPathCover(run_ids, edges, run_keys, nullptr);
    expect(paths.size() == 2, "multi-child support should win over single-child overvote");
    expect(paths[0] == std::vector<uint64_t>({1, 2}), "expected stable multi-child path");
}

void testAcyclicMatchingRejectsCycleClosingEdge() {
    std::vector<uint64_t> run_ids{1, 2, 3};
    std::vector<EdgeSupport> edges{
        EdgeSupport{1, 2, 4, {10, 11, 12, 13}},
        EdgeSupport{2, 3, 4, {10, 11, 12, 13}},
        EdgeSupport{3, 1, 1, {10}},
    };
    std::unordered_map<uint64_t, RunOrderKey> run_keys{
        {1, RunOrderKey{1, 10}},
        {2, RunOrderKey{1, 20}},
        {3, RunOrderKey{1, 30}},
    };
    auto paths = RaMesh::hal_export::buildMaximumCardinalityWeightPathCover(
        run_ids, edges, run_keys, nullptr);
    expect(
        paths ==
            std::vector<std::vector<uint64_t>>(
                {{1, 2, 3}}),
        "acyclic matching must reject the weakest cycle-closing edge during selection");
}

void testAcyclicMatchingPreservesStrongerEndpointEdges() {
    std::vector<uint64_t> run_ids{1, 2, 3, 4};
    std::vector<EdgeSupport> edges{
        EdgeSupport{1, 2, 1, {10, 11, 12}},
        EdgeSupport{2, 3, 1, {10, 11, 12}},
        EdgeSupport{3, 1, 1, {10, 11, 12}},
        EdgeSupport{4, 2, 1, {10, 11}},
    };
    std::unordered_map<uint64_t, RunOrderKey> run_keys{
        {1, RunOrderKey{1, 10}},
        {2, RunOrderKey{1, 20}},
        {3, RunOrderKey{1, 30}},
        {4, RunOrderKey{1, 40}},
    };
    auto paths =
        RaMesh::hal_export::buildMaximumCardinalityWeightPathCover(
            run_ids,
            edges,
            run_keys,
            nullptr);
    std::vector<uint64_t> observed;
    for (const auto& path : paths) {
        observed.insert(
            observed.end(),
            path.begin(),
            path.end());
    }
    std::sort(observed.begin(), observed.end());
    expect(
        paths.size() == 2 &&
            paths[0].size() + paths[1].size() == 4 &&
            (paths[0].size() == 3 ||
             paths[1].size() == 3) &&
            observed == run_ids,
        "acyclic matching must retain three compatible endpoint edges without losing a run");
}

void testPresenceMarginExposed() {
    TreeMeta tree = parseTree("((A:1,B:1)AB:1,C:1)root;");
    std::unordered_map<std::string, bool> leaf_presence{
        {"A", true},
        {"B", true},
        {"C", false},
    };
    BinaryInferenceResult result = RaMesh::hal_export::inferDescendantUnion(tree, leaf_presence);
    expect(result.score0.size() == tree.nodes.size(), "score0 size mismatch");
    expect(result.score1.size() == tree.nodes.size(), "score1 size mismatch");
    expect(result.margin.size() == tree.nodes.size(), "margin size mismatch");
    expect(result.margin.at(static_cast<size_t>(tree.name_to_id.at("AB"))) > 0.0, "AB margin should be positive");
}

void testDonorBucketDeduplication() {
    std::vector<BucketedDonor> donors{
        BucketedDonor{1, false, "AAAA", 10.0},
        BucketedDonor{1, true,  "CCCC", 1.0},
        BucketedDonor{2, false, "GGGG", 2.0},
    };
    auto selected = RaMesh::hal_export::selectBestDonorsByBucket(donors);
    expect(selected.size() == 2, "bucket dedup should keep one donor per bucket");
    expect(selected[0].first == "CCCC", "internal donor should win within same bucket");
    expect(selected[1].first == "GGGG", "second bucket donor mismatch");
}

void testConsensusThresholdControlsNRate() {
    std::vector<std::pair<std::string, double>> donors{
        {"A", 0.55},
        {"C", 0.45},
    };
    expect(RaMesh::hal_export::buildConsensusDNA(donors, 1, 0.60) == "N",
           "higher threshold should mask low-confidence column");
    expect(RaMesh::hal_export::buildConsensusDNA(donors, 1, 0.55) == "A",
           "lower threshold should keep best-supported base");
    std::vector<std::pair<std::string, double>> softmasked_donors{
        {"a", 0.7},
        {"A", 0.3},
    };
    expect(RaMesh::hal_export::buildConsensusDNA(softmasked_donors, 1, 0.60) == "a",
           "ancestor consensus must preserve weighted-majority softmask state");
}

void testBlockAnchorsPreserveDuplicateOccurrences() {
    auto block = RaMesh::Block::create(2);
    auto first = RaMesh::Segment::create(
        10,
        5,
        Strand::FORWARD,
        Cigar_t{cigarToInt('M', 5)},
        RaMesh::AlignRole::PRIMARY,
        RaMesh::SegmentRole::SEGMENT,
        block);
    auto second = RaMesh::Segment::create(
        30,
        5,
        Strand::FORWARD,
        Cigar_t{cigarToInt('M', 5)},
        RaMesh::AlignRole::PRIMARY,
        RaMesh::SegmentRole::SEGMENT,
        block);
    const RaMesh::SpeciesChrPair occurrence_key{"leafA", "chr1"};
    block->anchors.emplace(occurrence_key, first);
    block->anchors.emplace(occurrence_key, second);
    expect(block->anchors.count(occurrence_key) == 2,
           "block anchor storage must preserve repeated occurrences on one leaf sequence");
}
void testSparseDegreeTwoMatchingScalesLinearly() {
    constexpr uint64_t occurrence_count = 50000;
    std::vector<uint64_t> occurrence_ids;
    std::vector<EdgeSupport> edges;
    std::unordered_map<uint64_t, RunOrderKey>
        occurrence_keys;
    occurrence_ids.reserve(occurrence_count);
    edges.reserve(occurrence_count - 1);
    occurrence_keys.reserve(occurrence_count);
    for (uint64_t occurrence_id = 1;
         occurrence_id <= occurrence_count;
         ++occurrence_id) {
        occurrence_ids.push_back(occurrence_id);
        occurrence_keys.emplace(
            occurrence_id,
            RunOrderKey{occurrence_id, 0});
        if (occurrence_id != occurrence_count) {
            edges.push_back(
                EdgeSupport{
                    occurrence_id,
                    occurrence_id + 1,
                    2,
                    {10, 11},
                    false,
                    true,
                    1.0L,
                    0});
        }
    }
    auto paths =
        RaMesh::hal_export::buildMaximumCardinalityWeightPathCover(
            occurrence_ids,
            edges,
            occurrence_keys,
            nullptr);
    std::vector<uint64_t> observed;
    observed.reserve(occurrence_count);
    for (const auto& path : paths) {
        observed.insert(
            observed.end(),
            path.begin(),
            path.end());
    }
    std::sort(observed.begin(), observed.end());
    expect(
        paths.size() == occurrence_count / 2 &&
            observed == occurrence_ids,
        "large degree-two matching must remain exact without quadratic blossom storage");
}

void testReferenceSequenceCountDoesNotScaleWithEvidenceComponents() {
    constexpr uint64_t occurrence_count = 256;
    std::vector<uint64_t> occurrence_ids;
    std::unordered_map<uint64_t, RunOrderKey>
        occurrence_keys;
    occurrence_ids.reserve(occurrence_count);
    occurrence_keys.reserve(occurrence_count);
    for (uint64_t occurrence_id = 1;
         occurrence_id <= occurrence_count;
         ++occurrence_id) {
        occurrence_ids.push_back(occurrence_id);
        occurrence_keys.emplace(
            occurrence_id,
            RunOrderKey{occurrence_id, 0});
    }
    RaMesh::hal_export::ExportStats stats;
    const auto assembly =
        RaMesh::hal_export::buildAncestralSequenceAssembly(
            occurrence_ids,
            {},
            occurrence_keys,
            {},
            10,
            &stats);
    expect(
        assembly.sequences.size() == 1 &&
            stats.reference_interval_count == 1 &&
            stats.supported_join_count == 0 &&
            stats.scaffold_join_count ==
                occurrence_count - 1 &&
            stats.scaffold_gap_bases ==
                (occurrence_count - 1) * 10,
        "zero-weight top-level stub matching must prevent sequence count from scaling with disconnected evidence components");
}


void testReferenceThreadingTiePerturbationIsDeterministic() {
    const std::vector<uint64_t> occurrence_ids{
        1, 2, 3, 4, 5, 6};
    const std::vector<EdgeSupport> edges{
        EdgeSupport{
            1, 2, 2, {10, 11}, true, true, 8.0L, 0},
        EdgeSupport{
            3, 4, 2, {10, 11}, true, true, 8.0L, 0},
        EdgeSupport{
            5, 6, 2, {10, 11}, true, true, 8.0L, 0},
        EdgeSupport{
            2, 3, 2, {10, 11}, true, true, 2.0L, 5},
        EdgeSupport{
            2, 5, 2, {10, 11}, true, true, 2.0L, 5},
        EdgeSupport{
            4, 5, 2, {10, 11}, true, true, 2.0L, 5},
    };
    const std::unordered_map<uint64_t, RunOrderKey>
        occurrence_keys{
            {1, RunOrderKey{1, 10}},
            {2, RunOrderKey{1, 20}},
            {3, RunOrderKey{2, 10}},
            {4, RunOrderKey{2, 20}},
            {5, RunOrderKey{3, 10}},
            {6, RunOrderKey{3, 20}},
        };
    const std::vector<RaMesh::hal_export::TerminalEndSupport>
        terminal_ends{
            {{1, RaMesh::hal_export::OccurrenceEndSide::LEFT},
             2,
             {10, 11},
             8.0L},
            {{2, RaMesh::hal_export::OccurrenceEndSide::RIGHT},
             2,
             {10, 11},
             8.0L},
        };
    auto fingerprint = [&](const auto& assembly) {
        std::string value;
        for (const auto& sequence : assembly.sequences) {
            value += '[';
            for (const auto& fragment :
                 sequence.supported_fragments) {
                value += '{';
                for (uint64_t occurrence_id : fragment) {
                    value += std::to_string(occurrence_id);
                    value += ',';
                }
                value += '}';
            }
            for (const auto& join : sequence.joins) {
                value += ':';
                value += std::to_string(
                    static_cast<unsigned>(join.kind));
                value += ':';
                value += std::to_string(join.gap_length);
            }
            value += ']';
        }
        for (uint64_t occurrence_id : occurrence_ids) {
            value += assembly.forward_by_occurrence.at(
                         occurrence_id)
                         ? '+'
                         : '-';
        }
        return value;
    };

    std::string expected;
    for (size_t iteration = 0; iteration < 8; ++iteration) {
        auto assembly =
            RaMesh::hal_export::buildAncestralSequenceAssembly(
                occurrence_ids,
                edges,
                occurrence_keys,
                terminal_ends,
                10,
                nullptr);
        const std::string observed =
            fingerprint(assembly);
        if (iteration == 0) {
            expected = observed;
        } else {
            expect(
                observed == expected,
                "fixed-seed equal-score threading perturbations must be reproducible");
        }
    }
}


void testExportOrderingDeterministic() {
    std::vector<ExportBlockOrderKey> keys{
        {"simHuman", "chr2", 30, 9},
        {"simChimp", "chr1", 10, 7},
        {"simChimp", "chr1", 5, 8},
    };
    std::sort(keys.begin(), keys.end(), RaMesh::hal_export::exportBlockOrderLess);
    expect(keys[0].ref_start == 5 && keys[1].ref_start == 10 && keys[2].ref_species == "simHuman",
           "export block ordering should be lexicographic and stable");
}



} // namespace

int main() {
    testTreeMetaBuild();
    testSubtreeSharedInference();
    testSingleLeafInsertionInference();
    testDescendantUnionIgnoresBranchLength();
    testLeafIntervalProjection();
    testElementaryRunsPreserveAllOccurrences();
    testElementaryRunsRejectDuplicateOccurrenceIds();
    testC2HBottomFlagEncoding();
    testCactusZScoreWeightsLengthAndGap();
    testOutputSequenceOrderMatchesWriters();
    testInternalChildOrientationEmission();
    testLeafOrientationRespectsParentPlacement();
    testReversePlacementChangesAncestorDNA();
    testMixedOrientationPreservesAdjacency();
    testMaximumWeightChoosesSupportedJoin();
    testMaximumWeightMatchingIsGlobal();
    testEqualCardinalityAndWeightTieBreakIsDeterministic();
    testSupportedJoinCardinalityPrecedesWeight();
    testObservedEdgeIsRetained();
    testDisconnectedComponentsStaySeparate();
    testReferencePackingScaffoldsWithinEvidenceComponent();
    testDisconnectedEvidenceComponentsShareReferenceInterval();
    testDuplicateHomologyOccurrencesRemainDistinct();
    testConfirmedTerminalEndsPreserveReferenceIntervals();
    testSingleLineageTerminalEvidenceDoesNotCreateBoundaries();
    testConfirmedTerminalEndRemainsExternalDuringScaffolding();
    testMinimumScaffoldCoverPairsOppositeTerminals();
    testSingleLineageTerminalDoesNotOverrideSupportedJoin();
    testSingleChildOvervoteCannotOverrideMultiChildSupport();
    testAcyclicMatchingRejectsCycleClosingEdge();
    testAcyclicMatchingPreservesStrongerEndpointEdges();
    testPresenceMarginExposed();
    testDonorBucketDeduplication();
    testConsensusThresholdControlsNRate();
    testBlockAnchorsPreserveDuplicateOccurrences();
    testReferenceThreadingTiePerturbationIsDeterministic();
    testSparseDegreeTwoMatchingScalesLinearly();
    testReferenceSequenceCountDoesNotScaleWithEvidenceComponents();
    testExportOrderingDeterministic();
    return 0;
}
