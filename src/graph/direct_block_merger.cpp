#include "ramesh.h"

#include <spdlog/spdlog.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <numeric>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace RaMesh {
namespace {

using OrderedAnchors = std::map<SpeciesChrPair, SegPtr>;

struct BlockView {
    BlockPtr block;
    OrderedAnchors anchors;
    SegPtr reference_segment;
    ChrName declared_reference_chromosome;
    size_t species_count = 0;
};

struct CandidateChain {
    std::vector<BlockView*> views;
    std::vector<BlockPtr> blocks;
    ChrName reference_chromosome;
    uint_t reference_start = 0;
    uint_t reference_length = 0;
    size_t species_count = 0;
};

struct PathReplacement {
    SpeciesChrPair key;
    std::vector<SegPtr> old_segments;
    SegPtr merged_segment;
};

struct PreparedChain {
    CandidateChain candidate;
    BlockPtr merged_block;
    std::vector<PathReplacement> paths;
};

struct AppliedSplice {
    SegPtr old_first;
    SegPtr old_last;
    SegPtr previous;
    SegPtr next;
    SegPtr replacement;
};

struct SamplingSnapshot {
    GenomeEnd* genome_end = nullptr;
    std::vector<SegPtr> sample_vec;
};

constexpr uint64_t kMaximumCigarOperationLength = (1ULL << 28U) - 1U;

enum class BoundaryReason : size_t {
    MERGEABLE = 0,
    LEFT_BLOCK_INVALID,
    RIGHT_BLOCK_INVALID,
    BLOCK_REFERENCE_MISMATCH,
    PARTICIPANT_SPECIES_MISMATCH,
    PARTICIPANT_CHROMOSOME_MISMATCH,
    REFERENCE_GAP,
    REFERENCE_OVERLAP,
    REFERENCE_PATH_DISCONTINUITY,
    QUERY_STRAND_MISMATCH,
    QUERY_GAP,
    QUERY_OVERLAP,
    QUERY_PATH_DISCONTINUITY,
    CIGAR_INVALID,
    MAX_REFERENCE_SPAN,
    COUNT
};

constexpr size_t kBoundaryReasonCount =
    static_cast<size_t>(BoundaryReason::COUNT);
constexpr size_t kGapBinCount = 5;

struct BoundaryDiagnostics {
    size_t total_boundaries = 0;
    std::array<size_t, kBoundaryReasonCount> reasons{};
    std::array<size_t, kGapBinCount> reference_gap_bins{};
    std::array<size_t, kGapBinCount> query_gap_bins{};
    std::map<std::pair<size_t, size_t>, size_t> participant_transitions;
};

const char* boundaryReasonName(BoundaryReason reason) {
    switch (reason) {
        case BoundaryReason::MERGEABLE:
            return "mergeable";
        case BoundaryReason::LEFT_BLOCK_INVALID:
            return "left_block_invalid";
        case BoundaryReason::RIGHT_BLOCK_INVALID:
            return "right_block_invalid";
        case BoundaryReason::BLOCK_REFERENCE_MISMATCH:
            return "block_reference_mismatch";
        case BoundaryReason::PARTICIPANT_SPECIES_MISMATCH:
            return "participant_species_mismatch";
        case BoundaryReason::PARTICIPANT_CHROMOSOME_MISMATCH:
            return "participant_chromosome_mismatch";
        case BoundaryReason::REFERENCE_GAP:
            return "reference_gap";
        case BoundaryReason::REFERENCE_OVERLAP:
            return "reference_overlap";
        case BoundaryReason::REFERENCE_PATH_DISCONTINUITY:
            return "reference_path_discontinuity";
        case BoundaryReason::QUERY_STRAND_MISMATCH:
            return "query_strand_mismatch";
        case BoundaryReason::QUERY_GAP:
            return "query_gap";
        case BoundaryReason::QUERY_OVERLAP:
            return "query_overlap";
        case BoundaryReason::QUERY_PATH_DISCONTINUITY:
            return "query_path_discontinuity";
        case BoundaryReason::CIGAR_INVALID:
            return "cigar_invalid";
        case BoundaryReason::MAX_REFERENCE_SPAN:
            return "max_reference_span";
        case BoundaryReason::COUNT:
            break;
    }
    return "unknown";
}

size_t gapBin(uint64_t gap) {
    if (gap <= 10) {
        return 0;
    }
    if (gap <= 50) {
        return 1;
    }
    if (gap <= 100) {
        return 2;
    }
    if (gap <= 500) {
        return 3;
    }
    return 4;
}

const char* gapBinName(size_t bin) {
    static constexpr const char* names[kGapBinCount] = {
        "1_10", "11_50", "51_100", "101_500", "gt500"};
    return names[bin];
}

uint64_t segmentEnd(const SegPtr& segment) {
    return static_cast<uint64_t>(segment->start) +
           static_cast<uint64_t>(segment->length);
}

bool normalizeAndValidateCigar(const SegPtr& segment,
                               uint_t reference_length,
                               Cigar_t& normalized) {
    normalized.clear();
    if (!segment || !segment->isSegment() || !segment->isPrimary() ||
        segment->length == 0 || reference_length == 0) {
        return false;
    }

    if (segment->cigar.empty()) {
        if (segment->length != reference_length) {
            return false;
        }
        normalized.push_back(cigarToInt('M', reference_length));
        return true;
    }

    uint64_t reference_consumed = 0;
    uint64_t query_consumed = 0;
    for (const auto unit : segment->cigar) {
        char operation = '?';
        uint32_t length = 0;
        intToCigar(unit, operation, length);
        if (length == 0 ||
            (operation != 'M' && operation != 'I' && operation != 'D' &&
             operation != '=' && operation != 'X')) {
            return false;
        }

        if (operation == 'M' || operation == '=' || operation == 'X') {
            reference_consumed += length;
            query_consumed += length;
        } else if (operation == 'I') {
            query_consumed += length;
        } else {
            reference_consumed += length;
        }

        appendCigarOp(normalized, operation, length);
    }

    return reference_consumed == reference_length &&
           query_consumed == segment->length;
}

bool appendCigarChecked(Cigar_t& destination, const Cigar_t& source) {
    for (const auto unit : source) {
        char operation = '?';
        uint32_t length = 0;
        intToCigar(unit, operation, length);
        if (length == 0 ||
            (operation != 'M' && operation != 'I' && operation != 'D' &&
             operation != '=' && operation != 'X')) {
            return false;
        }

        if (!destination.empty()) {
            char previous_operation = '?';
            uint32_t previous_length = 0;
            intToCigar(destination.back(), previous_operation, previous_length);
            if (previous_operation == operation) {
                const uint64_t combined =
                    static_cast<uint64_t>(previous_length) + length;
                if (combined > kMaximumCigarOperationLength) {
                    return false;
                }
            }
        }
        appendCigarOp(destination, operation, length);
    }
    return true;
}

bool buildBlockView(const BlockPtr& block,
                    const SpeciesName& reference_species,
                    BlockView& view) {
    if (!block) {
        return false;
    }

    view = BlockView{};
    view.block = block;

    std::set<SpeciesName> species;
    {
        std::shared_lock block_lock(block->rw);
        if (block->ref_chr.empty() || block->anchors.size() < 2) {
            return false;
        }

        view.declared_reference_chromosome = block->ref_chr;
        for (const auto& [key, segment] : block->anchors) {
            if (!segment || !segment->isSegment() || !segment->isPrimary() ||
                segment->length == 0 || segment->parent_block.get() != block.get()) {
                return false;
            }
            if (!species.insert(key.first).second) {
                return false;
            }
            view.anchors.emplace(key, segment);
        }

        const SpeciesChrPair reference_key{reference_species, block->ref_chr};
        const auto reference_it = view.anchors.find(reference_key);
        if (reference_it == view.anchors.end() ||
            reference_it->second->strand != Strand::FORWARD) {
            return false;
        }
        view.reference_segment = reference_it->second;
    }

    view.species_count = species.size();
    return view.species_count >= 2;
}

bool sameAnchorKeys(const BlockView& left, const BlockView& right) {
    if (left.anchors.size() != right.anchors.size() ||
        left.species_count != right.species_count) {
        return false;
    }

    auto left_it = left.anchors.begin();
    auto right_it = right.anchors.begin();
    for (; left_it != left.anchors.end(); ++left_it, ++right_it) {
        if (left_it->first != right_it->first) {
            return false;
        }
    }
    return true;
}

int64_t segmentGapInReferenceOrder(const SegPtr& left,
                                   const SegPtr& right) {
    if (!left || !right || left->strand != right->strand) {
        return std::numeric_limits<int64_t>::min();
    }
    if (left->strand == Strand::FORWARD) {
        return static_cast<int64_t>(right->start) -
               static_cast<int64_t>(segmentEnd(left));
    }
    return static_cast<int64_t>(left->start) -
           static_cast<int64_t>(segmentEnd(right));
}

bool segmentPathsAreDirectlyAdjacent(const SegPtr& left,
                                     const SegPtr& right) {
    if (!left || !right || left->strand != right->strand) {
        return false;
    }
    if (left->strand == Strand::FORWARD) {
        return left->primary_path.next.load(std::memory_order_acquire) ==
                   right &&
               right->primary_path.prev.load(std::memory_order_acquire) ==
                   left;
    }
    return right->primary_path.next.load(std::memory_order_acquire) ==
               left &&
           left->primary_path.prev.load(std::memory_order_acquire) ==
               right;
}

bool segmentsAreMergeable(const SegPtr& left,
                          const SegPtr& right,
                          uint_t maximum_gap) {
    const int64_t gap = segmentGapInReferenceOrder(left, right);
    return gap >= 0 &&
           static_cast<uint64_t>(gap) <= maximum_gap &&
           segmentPathsAreDirectlyAdjacent(left, right);
}

bool segmentsAreDirectlyAdjacent(const SegPtr& left, const SegPtr& right) {
    return segmentsAreMergeable(left, right, 0);
}

bool canMergePair(const BlockView& left,
                  const BlockView& right,
                  const SpeciesName& reference_species,
                  uint_t maximum_query_gap) {
    if (!left.block || !right.block || left.block == right.block ||
        left.block->ref_chr != right.block->ref_chr ||
        !sameAnchorKeys(left, right)) {
        return false;
    }

    const SpeciesChrPair reference_key{
        reference_species, left.block->ref_chr};
    const auto left_reference_it = left.anchors.find(reference_key);
    const auto right_reference_it = right.anchors.find(reference_key);
    if (left_reference_it == left.anchors.end() ||
        right_reference_it == right.anchors.end() ||
        left.reference_segment != left_reference_it->second ||
        right.reference_segment != right_reference_it->second ||
        !segmentsAreDirectlyAdjacent(left.reference_segment,
                                     right.reference_segment)) {
        return false;
    }

    for (const auto& [key, left_segment] : left.anchors) {
        const auto right_it = right.anchors.find(key);
        const uint_t allowed_gap =
            key == reference_key ? 0 : maximum_query_gap;
        if (right_it == right.anchors.end() ||
            !segmentsAreMergeable(
                left_segment, right_it->second, allowed_gap)) {
            return false;
        }

        Cigar_t normalized;
        if (!normalizeAndValidateCigar(
                left_segment, left.reference_segment->length, normalized) ||
            !normalizeAndValidateCigar(
                right_it->second, right.reference_segment->length, normalized)) {
            return false;
        }
    }

    return true;
}

GenomeEnd& genomeEndFor(RaMeshMultiGenomeGraph& graph,
                        const SpeciesChrPair& key) {
    auto species_it = graph.species_graphs.find(key.first);
    if (species_it == graph.species_graphs.end()) {
        throw std::runtime_error(
            "Missing species graph for " + key.first);
    }
    auto chromosome_it = species_it->second.chr2end.find(key.second);
    if (chromosome_it == species_it->second.chr2end.end()) {
        throw std::runtime_error(
            "Missing chromosome graph for " + key.first + "." + key.second);
    }
    return chromosome_it->second;
}

PreparedChain prepareChain(const CandidateChain& candidate,
                           const SpeciesName& reference_species) {
    if (candidate.views.size() < 2 ||
        candidate.views.size() != candidate.blocks.size()) {
        throw std::runtime_error("Invalid exact-contiguous candidate chain");
    }

    PreparedChain prepared;
    prepared.candidate = candidate;
    prepared.merged_block =
        Block::createEmpty(candidate.reference_chromosome,
                           candidate.species_count);

    const auto& keys = candidate.views.front()->anchors;
    for (const auto& [key, unused] : keys) {
        (void)unused;
        PathReplacement path;
        path.key = key;
        path.old_segments.reserve(candidate.views.size());

        uint64_t merged_start = std::numeric_limits<uint64_t>::max();
        uint64_t merged_end = 0;
        uint64_t merged_length = 0;
        Cigar_t merged_cigar;
        Strand strand = Strand::FORWARD;

        for (size_t index = 0; index < candidate.views.size(); ++index) {
            const auto* view = candidate.views[index];
            const auto segment_it = view->anchors.find(key);
            if (segment_it == view->anchors.end()) {
                throw std::runtime_error(
                    "Candidate chain changed its participant set");
            }
            const auto& segment = segment_it->second;
            if (index == 0) {
                strand = segment->strand;
            } else if (segment->strand != strand) {
                throw std::runtime_error(
                    "Candidate chain changed strand while being prepared");
            }

            if (index > 0) {
                const auto previous_it =
                    candidate.views[index - 1]->anchors.find(key);
                if (previous_it ==
                    candidate.views[index - 1]->anchors.end()) {
                    throw std::runtime_error(
                        "Candidate chain lost its previous Segment");
                }
                const int64_t gap =
                    segmentGapInReferenceOrder(previous_it->second, segment);
                if (gap < 0 ||
                    static_cast<uint64_t>(gap) >
                        kMaximumCigarOperationLength) {
                    throw std::runtime_error(
                        "Candidate chain contains an invalid query gap");
                }
                appendCigarOp(
                    merged_cigar, 'I', static_cast<uint32_t>(gap));
                merged_length += static_cast<uint64_t>(gap);
            }
            Cigar_t normalized;
            if (!normalizeAndValidateCigar(
                    segment, view->reference_segment->length, normalized) ||
                !appendCigarChecked(merged_cigar, normalized)) {
                throw std::runtime_error(
                    "Candidate chain contains an invalid CIGAR");
            }

            path.old_segments.push_back(segment);
            merged_start = std::min<uint64_t>(merged_start, segment->start);
            merged_end = std::max<uint64_t>(
                merged_end, segmentEnd(segment));
            merged_length += segment->length;
        }

        if (merged_start > std::numeric_limits<uint_t>::max() ||
            merged_length == 0 ||
            merged_length > std::numeric_limits<uint_t>::max() ||
            merged_end - merged_start != merged_length) {
            throw std::runtime_error(
                "Merged Segment coordinates exceed the graph coordinate type");
        }

        path.merged_segment = Segment::create(
            static_cast<uint_t>(merged_start),
            static_cast<uint_t>(merged_length),
            strand,
            std::move(merged_cigar),
            AlignRole::PRIMARY,
            SegmentRole::SEGMENT,
            prepared.merged_block);
        path.merged_segment->left_extend =
            path.old_segments.front()->left_extend;
        path.merged_segment->right_extend =
            path.old_segments.back()->right_extend;

        prepared.merged_block->anchors.emplace(key, path.merged_segment);
        prepared.paths.push_back(std::move(path));
    }

    const SpeciesChrPair reference_key{
        reference_species, candidate.reference_chromosome};
    const auto merged_reference =
        prepared.merged_block->anchors.find(reference_key);
    if (merged_reference == prepared.merged_block->anchors.end() ||
        merged_reference->second->strand != Strand::FORWARD ||
        merged_reference->second->start != candidate.reference_start ||
        merged_reference->second->length != candidate.reference_length) {
        throw std::runtime_error(
            "Prepared Block has an inconsistent reference Segment");
    }

    return prepared;
}

void rebuildSampling(GenomeEnd& genome_end) {
    genome_end.sample_vec.clear();
    genome_end.sample_vec.resize(1, genome_end.head);
    auto current =
        genome_end.head->primary_path.next.load(std::memory_order_acquire);
    while (current && !current->isTail()) {
        genome_end.setToSampling(current);
        current =
            current->primary_path.next.load(std::memory_order_acquire);
    }
}

bool auditAffectedGraph(
    RaMeshMultiGenomeGraph& graph,
    const std::map<SpeciesChrPair, SamplingSnapshot>& affected,
    const std::vector<PreparedChain>& prepared,
    const std::unordered_set<const Block*>& old_blocks) {
    std::unordered_set<const Block*> pool_blocks;
    for (const auto& weak_block : graph.blocks) {
        const auto block = weak_block.lock();
        if (!block) {
            continue;
        }
        if (!pool_blocks.insert(block.get()).second ||
            old_blocks.count(block.get()) != 0) {
            return false;
        }
    }

    std::unordered_set<const Segment*> expected_new_segments;
    for (const auto& chain : prepared) {
        if (!chain.merged_block ||
            pool_blocks.count(chain.merged_block.get()) != 1 ||
            chain.merged_block->anchors.size() != chain.paths.size()) {
            return false;
        }
        for (const auto& path : chain.paths) {
            if (!path.merged_segment ||
                path.merged_segment->parent_block != chain.merged_block ||
                chain.merged_block->anchors.find(path.key) ==
                    chain.merged_block->anchors.end() ||
                chain.merged_block->anchors.at(path.key) !=
                    path.merged_segment ||
                !expected_new_segments.insert(
                    path.merged_segment.get()).second) {
                return false;
            }
        }
    }

    std::unordered_set<const Segment*> seen_new_segments;
    for (const auto& [key, snapshot] : affected) {
        const auto* genome_end = snapshot.genome_end;
        if (!genome_end || !genome_end->head || !genome_end->tail) {
            return false;
        }

        std::unordered_set<const Segment*> path_segments;
        auto previous = genome_end->head;
        auto current =
            genome_end->head->primary_path.next.load(std::memory_order_acquire);
        uint64_t previous_start = 0;
        bool have_previous_segment = false;

        while (current && !current->isTail()) {
            if (!current->isSegment() ||
                !path_segments.insert(current.get()).second ||
                current->primary_path.prev.load(std::memory_order_acquire) !=
                    previous ||
                previous->primary_path.next.load(std::memory_order_acquire) !=
                    current ||
                !current->parent_block) {
                return false;
            }
            if (have_previous_segment && current->start < previous_start) {
                return false;
            }
            previous_start = current->start;
            have_previous_segment = true;

            const auto parent = current->parent_block;
            const auto anchor_it = parent->anchors.find(key);
            if (anchor_it == parent->anchors.end() ||
                anchor_it->second != current ||
                pool_blocks.count(parent.get()) == 0) {
                return false;
            }

            if (expected_new_segments.count(current.get()) != 0) {
                if (!seen_new_segments.insert(current.get()).second) {
                    return false;
                }
            }

            previous = current;
            current =
                current->primary_path.next.load(std::memory_order_acquire);
        }

        if (!current || current != genome_end->tail ||
            genome_end->tail->primary_path.prev.load(
                std::memory_order_acquire) != previous ||
            previous->primary_path.next.load(std::memory_order_acquire) !=
                genome_end->tail) {
            return false;
        }
    }

    return seen_new_segments == expected_new_segments;
}

void detachPreparedBlocks(std::vector<PreparedChain>& prepared) {
    for (auto& chain : prepared) {
        for (auto& path : chain.paths) {
            if (path.merged_segment) {
                path.merged_segment->primary_path.prev.store(
                    nullptr, std::memory_order_release);
                path.merged_segment->primary_path.next.store(
                    nullptr, std::memory_order_release);
                path.merged_segment->parent_block.reset();
            }
        }
        if (chain.merged_block) {
            chain.merged_block->anchors.clear();
        }
    }
}

bool buildDiagnosticBlockView(const BlockPtr& block,
                              const SpeciesName& reference_species,
                              const ChrName& reference_chromosome,
                              BlockView& view) {
    if (!block) {
        return false;
    }

    view = BlockView{};
    view.block = block;

    std::set<SpeciesName> species;
    std::shared_lock block_lock(block->rw);
    if (block->anchors.size() < 2) {
        return false;
    }

    view.declared_reference_chromosome = block->ref_chr;
    for (const auto& [key, segment] : block->anchors) {
        if (!segment || !segment->isSegment() || !segment->isPrimary() ||
            segment->length == 0 ||
            segment->parent_block.get() != block.get() ||
            !species.insert(key.first).second) {
            return false;
        }
        view.anchors.emplace(key, segment);
    }

    const SpeciesChrPair reference_key{
        reference_species, reference_chromosome};
    const auto reference_it = view.anchors.find(reference_key);
    if (reference_it == view.anchors.end() ||
        reference_it->second->strand != Strand::FORWARD) {
        return false;
    }

    view.reference_segment = reference_it->second;
    view.species_count = species.size();
    return view.species_count >= 2;
}

bool sameParticipantSpecies(const BlockView& left,
                            const BlockView& right) {
    std::set<SpeciesName> left_species;
    std::set<SpeciesName> right_species;
    for (const auto& [key, unused] : left.anchors) {
        (void)unused;
        left_species.insert(key.first);
    }
    for (const auto& [key, unused] : right.anchors) {
        (void)unused;
        right_species.insert(key.first);
    }
    return left_species == right_species;
}

bool pathLinksAreDirectlyAdjacent(const SegPtr& left,
                                  const SegPtr& right) {
    if (!left || !right || left->strand != right->strand) {
        return false;
    }
    if (left->strand == Strand::FORWARD) {
        return left->primary_path.next.load(std::memory_order_acquire) ==
                   right &&
               right->primary_path.prev.load(std::memory_order_acquire) ==
                   left;
    }
    return right->primary_path.next.load(std::memory_order_acquire) ==
               left &&
           left->primary_path.prev.load(std::memory_order_acquire) ==
               right;
}

int64_t signedGapInReferenceOrder(const SegPtr& left,
                                  const SegPtr& right) {
    if (left->strand == Strand::FORWARD) {
        return static_cast<int64_t>(right->start) -
               static_cast<int64_t>(segmentEnd(left));
    }
    return static_cast<int64_t>(left->start) -
           static_cast<int64_t>(segmentEnd(right));
}

BoundaryReason classifyBoundary(
    const BlockView& left,
    const BlockView& right,
    const SpeciesName& reference_species,
    const ChrName& reference_chromosome,
    uint_t maximum_reference_span,
    BoundaryDiagnostics& diagnostics) {
    diagnostics.participant_transitions[
        {left.species_count, right.species_count}]++;

    if (left.block == right.block ||
        left.declared_reference_chromosome != reference_chromosome ||
        right.declared_reference_chromosome != reference_chromosome) {
        return BoundaryReason::BLOCK_REFERENCE_MISMATCH;
    }

    if (!sameParticipantSpecies(left, right)) {
        return BoundaryReason::PARTICIPANT_SPECIES_MISMATCH;
    }
    if (!sameAnchorKeys(left, right)) {
        return BoundaryReason::PARTICIPANT_CHROMOSOME_MISMATCH;
    }

    const int64_t reference_gap =
        signedGapInReferenceOrder(
            left.reference_segment, right.reference_segment);
    if (reference_gap > 0) {
        diagnostics.reference_gap_bins[
            gapBin(static_cast<uint64_t>(reference_gap))]++;
        return BoundaryReason::REFERENCE_GAP;
    }
    if (reference_gap < 0) {
        return BoundaryReason::REFERENCE_OVERLAP;
    }
    if (!pathLinksAreDirectlyAdjacent(
            left.reference_segment, right.reference_segment)) {
        return BoundaryReason::REFERENCE_PATH_DISCONTINUITY;
    }

    const SpeciesChrPair reference_key{
        reference_species, reference_chromosome};
    bool query_strand_mismatch = false;
    bool query_path_discontinuity = false;
    uint64_t maximum_query_gap = 0;
    uint64_t maximum_query_overlap = 0;

    for (const auto& [key, left_segment] : left.anchors) {
        if (key == reference_key) {
            continue;
        }
        const auto right_it = right.anchors.find(key);
        if (right_it == right.anchors.end()) {
            return BoundaryReason::PARTICIPANT_CHROMOSOME_MISMATCH;
        }
        const auto& right_segment = right_it->second;
        if (left_segment->strand != right_segment->strand) {
            query_strand_mismatch = true;
            continue;
        }

        const int64_t gap =
            signedGapInReferenceOrder(left_segment, right_segment);
        if (gap > 0) {
            maximum_query_gap = std::max(
                maximum_query_gap, static_cast<uint64_t>(gap));
        } else if (gap < 0) {
            maximum_query_overlap = std::max(
                maximum_query_overlap, static_cast<uint64_t>(-gap));
        }
        if (!pathLinksAreDirectlyAdjacent(left_segment, right_segment)) {
            query_path_discontinuity = true;
        }
    }

    if (query_strand_mismatch) {
        return BoundaryReason::QUERY_STRAND_MISMATCH;
    }
    if (maximum_query_gap > 0) {
        diagnostics.query_gap_bins[gapBin(maximum_query_gap)]++;
        return BoundaryReason::QUERY_GAP;
    }
    if (maximum_query_overlap > 0) {
        return BoundaryReason::QUERY_OVERLAP;
    }
    if (query_path_discontinuity) {
        return BoundaryReason::QUERY_PATH_DISCONTINUITY;
    }

    for (const auto& [key, left_segment] : left.anchors) {
        const auto right_it = right.anchors.find(key);
        Cigar_t normalized;
        if (right_it == right.anchors.end() ||
            !normalizeAndValidateCigar(
                left_segment, left.reference_segment->length, normalized) ||
            !normalizeAndValidateCigar(
                right_it->second, right.reference_segment->length,
                normalized)) {
            return BoundaryReason::CIGAR_INVALID;
        }
    }

    const uint64_t proposed_span =
        segmentEnd(right.reference_segment) -
        static_cast<uint64_t>(left.reference_segment->start);
    if (proposed_span > maximum_reference_span ||
        proposed_span > std::numeric_limits<uint_t>::max()) {
        return BoundaryReason::MAX_REFERENCE_SPAN;
    }

    return BoundaryReason::MERGEABLE;
}

BoundaryDiagnostics collectBoundaryDiagnostics(
    const RaMeshMultiGenomeGraph& graph,
    const SpeciesName& reference_species,
    uint_t maximum_reference_span) {
    BoundaryDiagnostics diagnostics;
    const auto reference_graph_it =
        graph.species_graphs.find(reference_species);
    if (reference_graph_it == graph.species_graphs.end()) {
        throw std::runtime_error(
            "Exact Block scan reference is absent from the graph: " +
            reference_species);
    }

    for (const auto& [chromosome, genome_end] :
         reference_graph_it->second.chr2end) {
        auto left_segment =
            genome_end.head->primary_path.next.load(
                std::memory_order_acquire);
        while (left_segment && !left_segment->isTail()) {
            const auto right_segment =
                left_segment->primary_path.next.load(
                    std::memory_order_acquire);
            if (!right_segment || right_segment->isTail()) {
                break;
            }

            ++diagnostics.total_boundaries;
            BlockView left;
            BlockView right;
            BoundaryReason reason = BoundaryReason::LEFT_BLOCK_INVALID;
            if (!buildDiagnosticBlockView(
                    left_segment->parent_block, reference_species,
                    chromosome, left)) {
                reason = BoundaryReason::LEFT_BLOCK_INVALID;
            } else if (!buildDiagnosticBlockView(
                           right_segment->parent_block, reference_species,
                           chromosome, right)) {
                reason = BoundaryReason::RIGHT_BLOCK_INVALID;
            } else {
                reason = classifyBoundary(
                    left, right, reference_species, chromosome,
                    maximum_reference_span, diagnostics);
            }
            diagnostics.reasons[static_cast<size_t>(reason)]++;
            left_segment = right_segment;
        }
    }

    return diagnostics;
}

void logBoundaryDiagnostics(const BoundaryDiagnostics& diagnostics,
                            const SpeciesName& reference_species,
                            const std::string& stage,
                            uint_t maximum_reference_span) {
    const size_t mergeable =
        diagnostics.reasons[
            static_cast<size_t>(BoundaryReason::MERGEABLE)];
    const size_t rejected =
        diagnostics.total_boundaries >= mergeable
            ? diagnostics.total_boundaries - mergeable
            : 0;
    spdlog::info(
        "[exact-block-scan] stage={} reference={} total_boundaries={} "
        "mergeable={} rejected={} max_reference_span={}",
        stage, reference_species, diagnostics.total_boundaries,
        mergeable, rejected, maximum_reference_span);

    size_t classified = 0;
    for (size_t index = 0; index < kBoundaryReasonCount; ++index) {
        const size_t count = diagnostics.reasons[index];
        classified += count;
        if (count == 0) {
            continue;
        }
        const double percent =
            diagnostics.total_boundaries == 0
                ? 0.0
                : 100.0 * static_cast<double>(count) /
                      static_cast<double>(diagnostics.total_boundaries);
        spdlog::info(
            "[exact-block-scan] stage={} reason={} count={} percent={:.6f}",
            stage,
            boundaryReasonName(static_cast<BoundaryReason>(index)),
            count, percent);
    }

    for (size_t bin = 0; bin < kGapBinCount; ++bin) {
        if (diagnostics.reference_gap_bins[bin] != 0) {
            spdlog::info(
                "[exact-block-scan] stage={} reference_gap_bin={} count={}",
                stage, gapBinName(bin),
                diagnostics.reference_gap_bins[bin]);
        }
        if (diagnostics.query_gap_bins[bin] != 0) {
            spdlog::info(
                "[exact-block-scan] stage={} query_max_gap_bin={} count={}",
                stage, gapBinName(bin),
                diagnostics.query_gap_bins[bin]);
        }
    }

    for (const auto& [transition, count] :
         diagnostics.participant_transitions) {
        spdlog::info(
            "[exact-block-scan] stage={} participant_transition={}_to_{} "
            "count={}",
            stage, transition.first, transition.second, count);
    }

    if (classified != diagnostics.total_boundaries) {
        throw std::runtime_error(
            "Exact Block scan classification count mismatch");
    }
}

}  // namespace

void RaMeshMultiGenomeGraph::inspectExactContiguousBlockBoundaries(
    const SpeciesName& reference_species,
    const std::string& stage,
    uint_t maximum_reference_span) const {
    if (maximum_reference_span == 0) {
        throw std::invalid_argument(
            "maximum_reference_span must be greater than zero");
    }

    std::shared_lock graph_lock(rw);
    const auto diagnostics = collectBoundaryDiagnostics(
        *this, reference_species, maximum_reference_span);
    logBoundaryDiagnostics(
        diagnostics, reference_species, stage, maximum_reference_span);
}

size_t RaMeshMultiGenomeGraph::mergeExactContiguousBlocks(
    const SpeciesName& reference_species,
    uint_t maximum_reference_span,
    uint_t maximum_query_gap) {
    if (maximum_reference_span == 0) {
        throw std::invalid_argument(
            "maximum_reference_span must be greater than zero");
    }

    std::unique_lock graph_lock(rw);

    const auto reference_graph_it = species_graphs.find(reference_species);
    if (reference_graph_it == species_graphs.end()) {
        throw std::runtime_error(
            "Exact Block merge reference is absent from the graph: " +
            reference_species);
    }

    const auto diagnostics = collectBoundaryDiagnostics(
        *this, reference_species, maximum_reference_span);
    logBoundaryDiagnostics(
        diagnostics, reference_species, "pre-merge-current-reference",
        maximum_reference_span);

    std::unordered_map<const Block*, BlockView> view_cache;
    std::unordered_set<const Block*> invalid_views;
    auto get_view = [&](const BlockPtr& block) -> BlockView* {
        if (!block) {
            return nullptr;
        }
        const auto cached = view_cache.find(block.get());
        if (cached != view_cache.end()) {
            return &cached->second;
        }
        if (invalid_views.count(block.get()) != 0) {
            return nullptr;
        }

        BlockView view;
        if (!buildBlockView(block, reference_species, view)) {
            invalid_views.insert(block.get());
            return nullptr;
        }
        return &view_cache.emplace(block.get(), std::move(view)).first->second;
    };

    std::vector<CandidateChain> candidates;
    for (auto& [chromosome, genome_end] :
         reference_graph_it->second.chr2end) {
        auto current =
            genome_end.head->primary_path.next.load(std::memory_order_acquire);
        while (current && !current->isTail()) {
            const auto current_block = current->parent_block;
            auto* first_view = get_view(current_block);
            if (!first_view || first_view->reference_segment != current ||
                first_view->block->ref_chr != chromosome) {
                current =
                    current->primary_path.next.load(std::memory_order_acquire);
                continue;
            }

            CandidateChain candidate;
            candidate.views.push_back(first_view);
            candidate.blocks.push_back(first_view->block);
            candidate.reference_chromosome = chromosome;
            candidate.reference_start = current->start;
            candidate.reference_length = current->length;
            candidate.species_count = first_view->species_count;

            auto last_reference = current;
            auto next_reference =
                current->primary_path.next.load(std::memory_order_acquire);
            while (next_reference && !next_reference->isTail()) {
                auto* next_view = get_view(next_reference->parent_block);
                if (!next_view ||
                    next_view->reference_segment != next_reference ||
                    next_view->species_count != candidate.species_count ||
                    !canMergePair(
                        *candidate.views.back(), *next_view,
                        reference_species, maximum_query_gap)) {
                    break;
                }

                const uint64_t proposed_end = segmentEnd(next_reference);
                const uint64_t proposed_span =
                    proposed_end - candidate.reference_start;
                if (proposed_span > maximum_reference_span ||
                    proposed_span > std::numeric_limits<uint_t>::max()) {
                    break;
                }

                candidate.views.push_back(next_view);
                candidate.blocks.push_back(next_view->block);
                candidate.reference_length =
                    static_cast<uint_t>(proposed_span);
                last_reference = next_reference;
                next_reference =
                    next_reference->primary_path.next.load(
                        std::memory_order_acquire);
            }

            if (candidate.blocks.size() > 1) {
                candidates.push_back(std::move(candidate));
                current =
                    last_reference->primary_path.next.load(
                        std::memory_order_acquire);
            } else {
                current = next_reference;
            }
        }
    }

    if (candidates.empty()) {
        spdlog::info(
            "[exact-block-merge] reference={} candidates=0 eliminated_boundaries=0",
            reference_species);
        return 0;
    }

    std::sort(
        candidates.begin(), candidates.end(),
        [](const CandidateChain& left, const CandidateChain& right) {
            if (left.species_count != right.species_count) {
                return left.species_count > right.species_count;
            }
            if (left.reference_chromosome != right.reference_chromosome) {
                return left.reference_chromosome <
                       right.reference_chromosome;
            }
            return left.reference_start < right.reference_start;
        });

    std::unordered_map<const Block*, size_t> old_block_to_candidate;
    std::unordered_set<const Block*> old_block_set;
    size_t eliminated_boundaries = 0;
    for (size_t candidate_index = 0;
         candidate_index < candidates.size(); ++candidate_index) {
        const auto& candidate = candidates[candidate_index];
        eliminated_boundaries += candidate.blocks.size() - 1;
        for (const auto& block : candidate.blocks) {
            if (!old_block_set.insert(block.get()).second) {
                throw std::runtime_error(
                    "Exact Block merge produced overlapping candidate chains");
            }
            old_block_to_candidate.emplace(block.get(), candidate_index);
        }
    }

    std::vector<PreparedChain> prepared;
    prepared.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        prepared.push_back(
            prepareChain(candidate, reference_species));
    }

    std::vector<WeakBlock> replacement_pool;
    replacement_pool.reserve(
        blocks.size() >= eliminated_boundaries
            ? blocks.size() - eliminated_boundaries
            : blocks.size());
    std::vector<bool> replacement_emitted(candidates.size(), false);
    std::vector<size_t> old_block_occurrences(candidates.size(), 0);
    std::unordered_set<const Block*> active_pool_blocks;
    size_t active_blocks_before = 0;

    for (const auto& weak_block : blocks) {
        const auto block = weak_block.lock();
        if (!block) {
            replacement_pool.push_back(weak_block);
            continue;
        }
        if (!active_pool_blocks.insert(block.get()).second) {
            throw std::runtime_error(
                "Exact Block merge found a duplicate Block pool entry");
        }
        ++active_blocks_before;

        const auto candidate_it =
            old_block_to_candidate.find(block.get());
        if (candidate_it == old_block_to_candidate.end()) {
            replacement_pool.push_back(weak_block);
            continue;
        }

        const size_t candidate_index = candidate_it->second;
        ++old_block_occurrences[candidate_index];
        if (!replacement_emitted[candidate_index]) {
            replacement_pool.emplace_back(
                prepared[candidate_index].merged_block);
            replacement_emitted[candidate_index] = true;
        }
    }

    for (size_t index = 0; index < candidates.size(); ++index) {
        if (!replacement_emitted[index] ||
            old_block_occurrences[index] != candidates[index].blocks.size()) {
            throw std::runtime_error(
                "Exact Block merge candidate is inconsistent with the Block pool");
        }
    }

    std::map<SpeciesChrPair, SamplingSnapshot> sampling_snapshots;
    for (const auto& chain : prepared) {
        for (const auto& path : chain.paths) {
            if (sampling_snapshots.count(path.key) != 0) {
                continue;
            }
            auto& genome_end = genomeEndFor(*this, path.key);
            sampling_snapshots.emplace(
                path.key,
                SamplingSnapshot{&genome_end, genome_end.sample_vec});
        }
    }

    std::vector<AppliedSplice> applied_splices;
    applied_splices.reserve(
        std::accumulate(
            prepared.begin(), prepared.end(), size_t{0},
            [](size_t total, const PreparedChain& chain) {
                return total + chain.paths.size();
            }));

    bool pool_replaced = false;
    try {
        for (auto& chain : prepared) {
            for (auto& path : chain.paths) {
                const bool reverse =
                    path.merged_segment->strand == Strand::REVERSE;
                const auto old_first =
                    reverse ? path.old_segments.back()
                            : path.old_segments.front();
                const auto old_last =
                    reverse ? path.old_segments.front()
                            : path.old_segments.back();

                const auto previous =
                    old_first->primary_path.prev.load(
                        std::memory_order_acquire);
                const auto next =
                    old_last->primary_path.next.load(
                        std::memory_order_acquire);
                if (!previous || !next ||
                    previous->primary_path.next.load(
                        std::memory_order_acquire) != old_first ||
                    next->primary_path.prev.load(
                        std::memory_order_acquire) != old_last) {
                    throw std::runtime_error(
                        "Exact Block merge path changed before commit");
                }

                path.merged_segment->primary_path.prev.store(
                    previous, std::memory_order_release);
                path.merged_segment->primary_path.next.store(
                    next, std::memory_order_release);
                previous->primary_path.next.store(
                    path.merged_segment, std::memory_order_release);
                next->primary_path.prev.store(
                    path.merged_segment, std::memory_order_release);

                applied_splices.push_back(
                    {old_first, old_last, previous, next,
                     path.merged_segment});
            }
        }

        blocks.swap(replacement_pool);
        pool_replaced = true;

        for (auto& [key, snapshot] : sampling_snapshots) {
            (void)key;
            rebuildSampling(*snapshot.genome_end);
        }

        std::atomic_thread_fence(std::memory_order_seq_cst);
        if (!auditAffectedGraph(
                *this, sampling_snapshots, prepared, old_block_set)) {
            throw std::runtime_error(
                "Exact Block merge post-commit graph audit failed");
        }
    } catch (...) {
        if (pool_replaced) {
            blocks.swap(replacement_pool);
        }

        for (auto splice = applied_splices.rbegin();
             splice != applied_splices.rend(); ++splice) {
            splice->previous->primary_path.next.store(
                splice->old_first, std::memory_order_release);
            splice->old_first->primary_path.prev.store(
                splice->previous, std::memory_order_release);
            splice->old_last->primary_path.next.store(
                splice->next, std::memory_order_release);
            splice->next->primary_path.prev.store(
                splice->old_last, std::memory_order_release);
            splice->replacement->primary_path.prev.store(
                nullptr, std::memory_order_release);
            splice->replacement->primary_path.next.store(
                nullptr, std::memory_order_release);
        }

        for (auto& [key, snapshot] : sampling_snapshots) {
            (void)key;
            snapshot.genome_end->sample_vec.swap(snapshot.sample_vec);
        }

        detachPreparedBlocks(prepared);
        std::atomic_thread_fence(std::memory_order_seq_cst);
        spdlog::error(
            "[exact-block-merge] commit failed; the original graph was restored");
        throw;
    }

    for (auto& chain : prepared) {
        for (auto& path : chain.paths) {
            for (auto& old_segment : path.old_segments) {
                old_segment->primary_path.prev.store(
                    nullptr, std::memory_order_release);
                old_segment->primary_path.next.store(
                    nullptr, std::memory_order_release);
                old_segment->parent_block.reset();
            }
        }
    }

    std::map<size_t, std::pair<size_t, size_t>> by_species_count;
    size_t merged_old_blocks = 0;
    size_t longest_chain = 0;
    for (const auto& candidate : candidates) {
        auto& stats = by_species_count[candidate.species_count];
        ++stats.first;
        stats.second += candidate.blocks.size() - 1;
        merged_old_blocks += candidate.blocks.size();
        longest_chain =
            std::max(longest_chain, candidate.blocks.size());
    }

    const size_t active_blocks_after =
        active_blocks_before - eliminated_boundaries;
    spdlog::info(
        "[exact-block-merge] reference={} candidates={} old_blocks={} "
        "eliminated_boundaries={} blocks_before={} blocks_after={} "
        "longest_chain={} max_reference_span={} max_query_gap={}",
        reference_species, candidates.size(), merged_old_blocks,
        eliminated_boundaries, active_blocks_before, active_blocks_after,
        longest_chain, maximum_reference_span, maximum_query_gap);
    for (const auto& [species_count, stats] : by_species_count) {
        spdlog::info(
            "[exact-block-merge] participants={} chains={} eliminated_boundaries={}",
            species_count, stats.first, stats.second);
    }

    return eliminated_boundaries;
}

}  // namespace RaMesh
