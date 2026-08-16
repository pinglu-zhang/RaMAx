#include "paf_export.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using RaMesh::Paf::Mode;
using RaMesh::Paf::SequencePair;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

bool hasPair(const std::vector<SequencePair>& pairs,
             std::size_t left, std::size_t right) {
    if (left > right) std::swap(left, right);
    return std::find(pairs.begin(), pairs.end(),
                     SequencePair{left, right}) != pairs.end();
}

void testReferenceStar() {
    const std::vector<std::string> rows{"ACGT", "ACGT", "ACGT"};
    const std::vector<std::string> names{"ref.c", "a.c", "b.c"};
    const auto selected = RaMesh::Paf::selectPairs(
        rows, names, 0, Mode::CONNECTED);
    require(selected.verified, "reference star must verify");
    require(selected.base_pairs == 2, "reference star base count");
    require(selected.supplemental_pairs == 0,
            "reference star must not add edges");
    require(selected.pairs.size() == 2, "reference star pair count");
}

void testReferenceGapInsertion() {
    const std::vector<std::string> rows{
        "AA---TT", "AACGGTT", "AACGGTT", "AACGGTT"};
    const std::vector<std::string> names{
        "ref.c", "a.c", "b.c", "c.c"};
    const auto selected = RaMesh::Paf::selectPairs(
        rows, names, 0, Mode::CONNECTED);
    require(selected.verified, "reference-gap insertion must verify");
    require(selected.base_pairs == 3, "reference-gap star count");
    require(selected.supplemental_pairs == 2,
            "three insertion rows require two supplemental pairs");
    require(selected.pairs.size() == 5, "reference-gap total pair count");
    const auto all = RaMesh::Paf::selectPairs(rows, names, 0, Mode::ALL);
    require(all.verified, "reference-gap all mode must verify");
    require(all.pairs.size() == 6, "reference-gap all pair count");
    require(selected.pairs.size() < all.pairs.size(),
            "connected must be strictly sparser in this fixture");
}

void testSharedNonReferenceAllele() {
    const std::vector<std::string> rows{"A", "G", "G"};
    const std::vector<std::string> names{"ref.c", "a.c", "b.c"};
    const auto selected = RaMesh::Paf::selectPairs(
        rows, names, 0, Mode::CONNECTED);
    require(selected.verified, "shared non-reference allele must verify");
    require(selected.base_pairs == 2, "allele star count");
    require(selected.supplemental_pairs == 1,
            "G/G requires a supplemental edge");
    require(hasPair(selected.pairs, 1, 2), "missing G/G edge");
}

void testPairReuseAndChangingParticipants() {
    const std::vector<std::string> rows{
        "--AC-GT", "TTAC-GT", "TTACCGT", "--ACCGT"};
    const std::vector<std::string> names{
        "ref.c", "a.c", "b.c", "c.c"};
    const auto selected = RaMesh::Paf::selectPairs(
        rows, names, 0, Mode::CONNECTED);
    require(selected.verified, "nested gaps must verify");
    require(std::adjacent_find(selected.pairs.begin(), selected.pairs.end()) ==
                selected.pairs.end(),
            "a pair must be emitted once per Block");
    const auto repeated = RaMesh::Paf::selectPairs(
        rows, names, 0, Mode::CONNECTED);
    require(selected.pairs == repeated.pairs,
            "connected selection must be deterministic");
}

void testHubSelection() {
    const std::vector<std::string> rows{
        "--AA", "TTAA", "T-AA", "TTAA"};
    const std::vector<std::string> names{
        "ref.c", "z.c", "m.c", "a.c"};
    const auto selected = RaMesh::Paf::selectPairs(
        rows, names, 0, Mode::CONNECTED);
    require(selected.verified, "hub-selection fixture must verify");
    require(selected.supplemental_pairs == 2,
            "hub selection supplemental count");
    require(hasPair(selected.pairs, 1, 3),
            "longest/name-minimum hub must connect row 1");
    require(hasPair(selected.pairs, 2, 3),
            "longest/name-minimum hub must connect row 2");
}

void testAllModeEligibility() {
    const std::vector<std::string> rows{"AA--", "--TT", "A-T-"};
    const std::vector<std::string> names{"ref.c", "a.c", "b.c"};
    const auto selected = RaMesh::Paf::selectPairs(rows, names, 0, Mode::ALL);
    require(selected.theoretical_pairs == 3, "all theoretical pair count");
    require(selected.eligible_pairs == 2, "all eligible pair count");
    require(selected.pairs.size() == 2, "all output pair count");
    require(!hasPair(selected.pairs, 0, 1),
            "disjoint rows must not produce PAF");
    require(selected.verified, "eligible all mode must verify");
}

void testCaseAndXConnectivity() {
    const std::vector<std::string> rows{"aX", "AX", "gX"};
    const std::vector<SequencePair> only_a{{0, 1}};
    require(RaMesh::Paf::verifyColumnConnectivity(rows, only_a),
            "case-normalized A rows should connect and X should be ignored");
    require(!RaMesh::Paf::verifyColumnConnectivity(
                std::vector<std::string>{"a", "A"}, {}),
            "unconnected case-normalized bases must fail");
}

void testProjection() {
    const auto forward = RaMesh::Paf::projectPair(
        "AC-GTX", "AcT-TX", false);
    require(forward.valid, "forward projection must be valid");
    require(forward.cigar == "2=1I1D1=1X", "forward extended CIGAR");
    require(forward.matching_bases == 3, "forward match count");
    require(forward.block_length == 6, "forward block length");
    require(forward.target_consumed == 5, "forward target consumption");
    require(forward.query_consumed == 5, "forward query consumption");
    require(forward.edit_distance == 3, "forward edit distance");

    const auto reverse = RaMesh::Paf::projectPair(
        "AC-GTX", "AcT-TX", true);
    require(reverse.valid, "reverse projection must be valid");
    require(reverse.cigar == "1X1=1D1I2=", "reverse extended CIGAR");
    require(reverse.matching_bases == forward.matching_bases,
            "reverse match count");
    require(reverse.edit_distance == forward.edit_distance,
            "reverse edit distance");
}

void testInvalidInputs() {
    bool threw = false;
    try {
        (void)RaMesh::Paf::selectPairs(
            std::vector<std::string>{"AA", "A"},
            std::vector<std::string>{"a", "b"}, 0, Mode::CONNECTED);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    require(threw, "unequal MSA rows must be rejected");

    threw = false;
    try {
        (void)RaMesh::Paf::projectPair("AA", "A", false);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    require(threw, "unequal projection rows must be rejected");
}

}  // namespace

int main() {
    try {
        testReferenceStar();
        testReferenceGapInsertion();
        testSharedNonReferenceAllele();
        testPairReuseAndChangingParticipants();
        testHubSelection();
        testAllModeEligibility();
        testCaseAndXConnectivity();
        testProjection();
        testInvalidInputs();
    } catch (const std::exception& error) {
        std::cerr << "ramax_paf_tests: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    std::cout << "ramax_paf_tests: all checks passed\n";
    return EXIT_SUCCESS;
}
