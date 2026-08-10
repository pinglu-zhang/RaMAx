#include "ramesh.h"
#include "align.h"

#include <spdlog/spdlog.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <exception>
#include <numeric>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
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

void detachPreparedChain(PreparedChain& chain) {
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

void detachPreparedBlocks(std::vector<PreparedChain>& prepared) {
    for (auto& chain : prepared) {
        detachPreparedChain(chain);
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

enum class MissingWindowReject : size_t {
    INVALID_BLOCK = 0,
    PARTICIPANT_PATTERN,
    ANCHOR_MISMATCH,
    PATH_MISMATCH,
    INTERVAL_INVALID,
    ANCHOR_ORDER_MISMATCH,
    MIXED_ZERO_NONZERO_INTERVAL,
    ZERO_GAP_SPAN_EXCEEDED,
    ZERO_GAP_NOT_EXACT_CONTIGUOUS,
    ZERO_GAP_PREPARATION_INVALID,
    OVERLAPPING_PREPARED_WINDOW,
    SEQUENCE_UNAVAILABLE,
    MINIPOA_FAILED,
    MSA_INVALID,
    COUNT
};

const char* missingWindowRejectName(MissingWindowReject reason) {
    switch (reason) {
        case MissingWindowReject::INVALID_BLOCK:
            return "invalid_block";
        case MissingWindowReject::PARTICIPANT_PATTERN:
            return "participant_pattern";
        case MissingWindowReject::ANCHOR_MISMATCH:
            return "anchor_mismatch";
        case MissingWindowReject::PATH_MISMATCH:
            return "path_mismatch";
        case MissingWindowReject::INTERVAL_INVALID:
            return "interval_invalid";
        case MissingWindowReject::ANCHOR_ORDER_MISMATCH:
            return "anchor_order_mismatch";
        case MissingWindowReject::MIXED_ZERO_NONZERO_INTERVAL:
            return "mixed_zero_nonzero_interval";
        case MissingWindowReject::ZERO_GAP_SPAN_EXCEEDED:
            return "zero_gap_span_exceeded";
        case MissingWindowReject::ZERO_GAP_NOT_EXACT_CONTIGUOUS:
            return "zero_gap_not_exact_contiguous";
        case MissingWindowReject::ZERO_GAP_PREPARATION_INVALID:
            return "zero_gap_preparation_invalid";
        case MissingWindowReject::OVERLAPPING_PREPARED_WINDOW:
            return "overlapping_prepared_window";
        case MissingWindowReject::SEQUENCE_UNAVAILABLE:
            return "sequence_unavailable";
        case MissingWindowReject::MINIPOA_FAILED:
            return "minipoa_failed";
        case MissingWindowReject::MSA_INVALID:
            return "msa_invalid";
        case MissingWindowReject::COUNT:
            break;
    }
    return "unknown";
}

struct MissingWindowPathContext {
    SpeciesChrPair key;
    SegPtr left_anchor;
    std::vector<SegPtr> interior_segments;
    SegPtr right_anchor;
    uint_t interval_start = 0;
    uint_t interval_length = 0;
    Strand strand = Strand::FORWARD;
};

struct MissingWindowCandidate {
    BlockView left;
    std::vector<BlockView> interiors;
    BlockView right;
    std::set<SpeciesName> missing_species;
    ChrName reference_chromosome;
    std::vector<MissingWindowPathContext> paths;
    std::string participant_pattern;
    bool zero_gap_deletion = false;
};

bool buildRealignBlockView(const BlockPtr& block,
                           const SpeciesName& reference_species,
                           BlockView& view) {
    if (!block) {
        return false;
    }
    view = BlockView{};
    view.block = block;
    std::set<SpeciesName> species;
    std::shared_lock block_lock(block->rw);
    if (block->ref_chr.empty() || block->anchors.empty()) {
        return false;
    }
    view.declared_reference_chromosome = block->ref_chr;
    for (const auto& [key, segment] : block->anchors) {
        if (!segment || !segment->isSegment() || !segment->isPrimary() ||
            segment->length == 0 || segment->parent_block.get() != block.get() ||
            !species.insert(key.first).second) {
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
    view.species_count = species.size();
    return true;
}

bool participantSpeciesAreSubset(const BlockView& subset,
                                 const BlockView& superset) {
    for (const auto& [key, unused] : subset.anchors) {
        (void)unused;
        if (superset.anchors.find(key) == superset.anchors.end()) {
            return false;
        }
    }
    return true;
}

MissingWindowReject buildMissingWindowCandidate(
    const BlockPtr& left_block,
    const std::vector<BlockPtr>& interior_blocks,
    const BlockPtr& right_block,
    const SpeciesName& reference_species,
    const ChrName& reference_chromosome,
    size_t expected_species_count,
    uint_t maximum_span,
    uint_t zero_gap_maximum_span,
    bool zero_gap_phase,
    MissingWindowCandidate& candidate) {
    candidate = MissingWindowCandidate{};
    if (!buildBlockView(left_block, reference_species, candidate.left) ||
        !buildBlockView(right_block, reference_species, candidate.right) ||
        interior_blocks.empty()) {
        return MissingWindowReject::INVALID_BLOCK;
    }
    candidate.interiors.reserve(interior_blocks.size());
    for (const auto& block : interior_blocks) {
        BlockView view;
        if (!buildRealignBlockView(block, reference_species, view)) {
            return MissingWindowReject::INVALID_BLOCK;
        }
        candidate.interiors.push_back(std::move(view));
    }
    if (candidate.left.declared_reference_chromosome !=
            reference_chromosome ||
        candidate.right.declared_reference_chromosome !=
            reference_chromosome ||
        candidate.left.species_count != expected_species_count ||
        candidate.right.species_count != expected_species_count ||
        !sameAnchorKeys(candidate.left, candidate.right)) {
        return MissingWindowReject::PARTICIPANT_PATTERN;
    }
    if (zero_gap_phase && candidate.interiors.size() != 1) {
        return MissingWindowReject::PARTICIPANT_PATTERN;
    }
    candidate.participant_pattern =
        std::to_string(expected_species_count);
    for (const auto& interior : candidate.interiors) {
        if (interior.declared_reference_chromosome != reference_chromosome ||
            interior.species_count >= expected_species_count ||
            !participantSpeciesAreSubset(interior, candidate.left)) {
            return MissingWindowReject::PARTICIPANT_PATTERN;
        }
        candidate.participant_pattern +=
            "-" + std::to_string(interior.species_count);
    }
    candidate.participant_pattern +=
        "-" + std::to_string(expected_species_count);
    candidate.reference_chromosome = reference_chromosome;
    candidate.paths.reserve(expected_species_count);

    for (const auto& [key, left_segment] : candidate.left.anchors) {
        const auto right_it = candidate.right.anchors.find(key);
        if (right_it == candidate.right.anchors.end() ||
            !left_segment || !right_it->second ||
            left_segment->strand != right_it->second->strand) {
            return MissingWindowReject::ANCHOR_MISMATCH;
        }
        const auto& right_segment = right_it->second;
        const int64_t interval_length =
            segmentGapInReferenceOrder(left_segment, right_segment);
        if (interval_length < 0) {
            return MissingWindowReject::ANCHOR_ORDER_MISMATCH;
        }
        if (static_cast<uint64_t>(interval_length) >
            std::numeric_limits<uint_t>::max()) {
            return MissingWindowReject::INTERVAL_INVALID;
        }
        MissingWindowPathContext path;
        path.key = key;
        path.left_anchor = left_segment;
        path.right_anchor = right_segment;
        path.interval_length =
            static_cast<uint_t>(interval_length);
        path.strand = left_segment->strand;
        path.interval_start =
            path.strand == Strand::FORWARD
                ? static_cast<uint_t>(segmentEnd(left_segment))
                : static_cast<uint_t>(segmentEnd(right_segment));

        SegPtr previous = left_segment;
        for (const auto& interior : candidate.interiors) {
            const auto interior_it = interior.anchors.find(key);
            if (interior_it == interior.anchors.end()) {
                candidate.missing_species.insert(key.first);
                continue;
            }
            const auto& segment = interior_it->second;
            if (!segment || segment->strand != path.strand ||
                segmentGapInReferenceOrder(previous, segment) < 0 ||
                !segmentPathsAreDirectlyAdjacent(previous, segment)) {
                return MissingWindowReject::PATH_MISMATCH;
            }
            path.interior_segments.push_back(segment);
            previous = segment;
        }
        if (segmentGapInReferenceOrder(previous, right_segment) < 0 ||
            !segmentPathsAreDirectlyAdjacent(previous, right_segment)) {
            return MissingWindowReject::PATH_MISMATCH;
        }
        candidate.paths.push_back(std::move(path));
    }
    if (candidate.missing_species.empty()) {
        return MissingWindowReject::PARTICIPANT_PATTERN;
    }

    if (!zero_gap_phase) {
        bool any_zero = false;
        bool any_nonzero = false;
        for (const auto& path : candidate.paths) {
            any_zero = any_zero || path.interval_length == 0;
            any_nonzero = any_nonzero || path.interval_length != 0;
            if (path.interval_length > maximum_span) {
                return MissingWindowReject::INTERVAL_INVALID;
            }
        }
        if (any_zero) {
            return any_nonzero
                       ? MissingWindowReject::MIXED_ZERO_NONZERO_INTERVAL
                       : MissingWindowReject::INTERVAL_INVALID;
        }
        return MissingWindowReject::COUNT;
    }

    candidate.zero_gap_deletion = true;
    bool missing_zero = false;
    bool missing_nonzero = false;
    for (const auto& path : candidate.paths) {
        if (path.interior_segments.empty()) {
            missing_zero = missing_zero || path.interval_length == 0;
            missing_nonzero = missing_nonzero || path.interval_length != 0;
            if (path.interval_length == 0 &&
                !segmentsAreDirectlyAdjacent(
                    path.left_anchor, path.right_anchor)) {
                return MissingWindowReject::ZERO_GAP_NOT_EXACT_CONTIGUOUS;
            }
            continue;
        }

        if (path.interval_length == 0 ||
            path.interval_length > zero_gap_maximum_span) {
            return MissingWindowReject::ZERO_GAP_SPAN_EXCEEDED;
        }
        if (path.interior_segments.size() != 1 ||
            !segmentsAreDirectlyAdjacent(
                path.left_anchor, path.interior_segments.front()) ||
            !segmentsAreDirectlyAdjacent(
                path.interior_segments.front(), path.right_anchor)) {
            return MissingWindowReject::ZERO_GAP_NOT_EXACT_CONTIGUOUS;
        }
    }

    if (missing_nonzero) {
        return missing_zero
                   ? MissingWindowReject::MIXED_ZERO_NONZERO_INTERVAL
                   : MissingWindowReject::INTERVAL_INVALID;
    }

    return MissingWindowReject::COUNT;
}

bool fetchWindowSequence(
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    const MissingWindowPathContext& path,
    std::string& sequence) {
    const auto manager_it = managers.find(path.key.first);
    if (manager_it == managers.end() || !manager_it->second) {
        return false;
    }
    try {
        sequence = std::visit(
            [&](auto& manager) -> std::string {
                using T = std::decay_t<decltype(manager)>;
                if constexpr (std::is_same_v<
                                  T,
                                  std::unique_ptr<
                                      SeqPro::SequenceManager>>) {
                    return manager->getSubSequence(
                        path.key.second, path.interval_start,
                        path.interval_length);
                } else {
                    return manager->getOriginalManager().getSubSequence(
                        path.key.second, path.interval_start,
                        path.interval_length);
                }
            },
            *manager_it->second);
    } catch (const std::exception&) {
        return false;
    }
    if (sequence.size() != path.interval_length) {
        return false;
    }
    if (path.strand == Strand::REVERSE) {
        reverseComplement(sequence);
    }
    return true;
}

bool cigarFromMsaRows(const std::string& reference,
                      const std::string& query,
                      Cigar_t& cigar) {
    cigar.clear();
    if (reference.empty() || reference.size() != query.size()) {
        return false;
    }
    for (size_t index = 0; index < reference.size(); ++index) {
        const bool reference_gap = reference[index] == '-';
        const bool query_gap = query[index] == '-';
        if (reference_gap && query_gap) {
            continue;
        }
        const char operation =
            reference_gap ? 'I' : (query_gap ? 'D' : 'M');
        appendCigarOp(cigar, operation, 1);
    }
    return !cigar.empty();
}

std::optional<PreparedChain> prepareZeroGapDeletionWindow(
    const MissingWindowCandidate& candidate,
    const SpeciesName& reference_species,
    MissingWindowReject& rejection) {
    rejection = MissingWindowReject::ZERO_GAP_PREPARATION_INVALID;
    if (!candidate.zero_gap_deletion || candidate.interiors.size() != 1 ||
        candidate.paths.empty() ||
        !candidate.left.reference_segment ||
        !candidate.interiors.front().reference_segment ||
        !candidate.right.reference_segment) {
        return std::nullopt;
    }

    const auto& middle = candidate.interiors.front();

    const uint64_t reference_length =
        static_cast<uint64_t>(candidate.left.reference_segment->length) +
        middle.reference_segment->length +
        candidate.right.reference_segment->length;
    if (reference_length == 0 ||
        reference_length > std::numeric_limits<uint_t>::max() ||
        candidate.left.reference_segment->strand != Strand::FORWARD ||
        middle.reference_segment->strand != Strand::FORWARD ||
        candidate.right.reference_segment->strand != Strand::FORWARD ||
        segmentEnd(candidate.right.reference_segment) -
                candidate.left.reference_segment->start !=
            reference_length) {
        return std::nullopt;
    }

    PreparedChain prepared;
    prepared.candidate.blocks = {
        candidate.left.block,
        middle.block,
        candidate.right.block};
    prepared.candidate.reference_chromosome =
        candidate.reference_chromosome;
    prepared.candidate.reference_start =
        candidate.left.reference_segment->start;
    prepared.candidate.reference_length =
        static_cast<uint_t>(reference_length);
    prepared.candidate.species_count = candidate.paths.size();
    prepared.merged_block = Block::createEmpty(
        candidate.reference_chromosome, candidate.paths.size());
    prepared.paths.reserve(candidate.paths.size());
    const auto reject_prepared = [&]() -> std::optional<PreparedChain> {
        detachPreparedChain(prepared);
        return std::nullopt;
    };

    for (const auto& context : candidate.paths) {
        PathReplacement path;
        path.key = context.key;
        path.old_segments.push_back(context.left_anchor);
        if (!context.interior_segments.empty()) {
            path.old_segments.push_back(context.interior_segments.front());
        }
        path.old_segments.push_back(context.right_anchor);

        Cigar_t merged_cigar;
        Cigar_t normalized;
        if (!normalizeAndValidateCigar(
                context.left_anchor,
                candidate.left.reference_segment->length,
                normalized) ||
            !appendCigarChecked(merged_cigar, normalized)) {
            return reject_prepared();
        }

        if (!context.interior_segments.empty()) {
            if (!normalizeAndValidateCigar(
                    context.interior_segments.front(),
                    middle.reference_segment->length,
                    normalized) ||
                !appendCigarChecked(merged_cigar, normalized)) {
                return reject_prepared();
            }
        } else {
            appendCigarOp(
                merged_cigar, 'D',
                middle.reference_segment->length);
        }

        if (!normalizeAndValidateCigar(
                context.right_anchor,
                candidate.right.reference_segment->length,
                normalized) ||
            !appendCigarChecked(merged_cigar, normalized)) {
            return reject_prepared();
        }

        uint64_t merged_start =
            std::numeric_limits<uint64_t>::max();
        uint64_t merged_end = 0;
        uint64_t merged_length = 0;
        for (const auto& segment : path.old_segments) {
            if (!segment || segment->strand != context.strand) {
                return reject_prepared();
            }
            merged_start = std::min<uint64_t>(
                merged_start, segment->start);
            merged_end = std::max<uint64_t>(
                merged_end, segmentEnd(segment));
            merged_length += segment->length;
        }
        if (merged_start > std::numeric_limits<uint_t>::max() ||
            merged_length == 0 ||
            merged_length > std::numeric_limits<uint_t>::max() ||
            merged_end - merged_start != merged_length ||
            countRefLength(merged_cigar) != reference_length ||
            countQryLength(merged_cigar) != merged_length) {
            return reject_prepared();
        }

        path.merged_segment = Segment::create(
            static_cast<uint_t>(merged_start),
            static_cast<uint_t>(merged_length), context.strand,
            std::move(merged_cigar), AlignRole::PRIMARY,
            SegmentRole::SEGMENT, prepared.merged_block);
        path.merged_segment->left_extend =
            path.old_segments.front()->left_extend;
        path.merged_segment->right_extend =
            path.old_segments.back()->right_extend;
        if (!prepared.merged_block->anchors.emplace(
                context.key, path.merged_segment).second) {
            return reject_prepared();
        }
        prepared.paths.push_back(std::move(path));
    }

    const SpeciesChrPair reference_key{
        reference_species, candidate.reference_chromosome};
    const auto reference_it =
        prepared.merged_block->anchors.find(reference_key);
    if (reference_it == prepared.merged_block->anchors.end() ||
        reference_it->second->strand != Strand::FORWARD ||
        reference_it->second->start != prepared.candidate.reference_start ||
        reference_it->second->length !=
            prepared.candidate.reference_length) {
        return reject_prepared();
    }

    rejection = MissingWindowReject::COUNT;
    return prepared;
}

std::optional<PreparedChain> prepareMissingWindow(
    const MissingWindowCandidate& candidate,
    const SpeciesName& reference_species,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    const std::string& msa_executable,
    MissingWindowReject& rejection) {
    std::unordered_map<ChrName, std::string> raw_sequences;
    raw_sequences.reserve(candidate.paths.size());
    for (const auto& path : candidate.paths) {
        std::string sequence;
        if (!fetchWindowSequence(managers, path, sequence) ||
            !raw_sequences.emplace(
                path.key.first, std::move(sequence)).second) {
            rejection = MissingWindowReject::SEQUENCE_UNAVAILABLE;
            return std::nullopt;
        }
    }

    auto aligned_sequences = raw_sequences;
    bool msa_succeeded = false;
    try {
        msa_succeeded = alignSequencesWithExternalMsa(
            msa_executable, aligned_sequences);
    } catch (const std::exception& error) {
        spdlog::warn(
            "[species-mismatch-realign] minipoa invocation failed: {}",
            error.what());
    }
    if (!msa_succeeded) {
        rejection = MissingWindowReject::MINIPOA_FAILED;
        return std::nullopt;
    }
    const auto reference_it =
        aligned_sequences.find(reference_species);
    if (reference_it == aligned_sequences.end()) {
        rejection = MissingWindowReject::MSA_INVALID;
        return std::nullopt;
    }

    PreparedChain prepared;
    for (const auto& interior : candidate.interiors) {
        prepared.candidate.blocks.push_back(interior.block);
    }
    prepared.candidate.reference_chromosome =
        candidate.reference_chromosome;
    prepared.candidate.species_count = candidate.paths.size();
    prepared.merged_block = Block::createEmpty(
        candidate.reference_chromosome, candidate.paths.size());
    prepared.paths.reserve(candidate.paths.size());

    for (const auto& context : candidate.paths) {
        const auto aligned_it =
            aligned_sequences.find(context.key.first);
        if (aligned_it == aligned_sequences.end()) {
            rejection = MissingWindowReject::MSA_INVALID;
            return std::nullopt;
        }

        Cigar_t cigar;
        if (context.key.first == reference_species) {
            cigar.push_back(
                cigarToInt('M', context.interval_length));
        } else if (!cigarFromMsaRows(
                       reference_it->second, aligned_it->second,
                       cigar)) {
            rejection = MissingWindowReject::MSA_INVALID;
            return std::nullopt;
        }
        if (countRefLength(cigar) !=
                raw_sequences.at(reference_species).size() ||
            countQryLength(cigar) != context.interval_length) {
            rejection = MissingWindowReject::MSA_INVALID;
            return std::nullopt;
        }

        PathReplacement path;
        path.key = context.key;
        path.old_segments = context.interior_segments;
        path.merged_segment = Segment::create(
            context.interval_start, context.interval_length,
            context.strand, std::move(cigar),
            AlignRole::PRIMARY, SegmentRole::SEGMENT,
            prepared.merged_block);
        if (!context.interior_segments.empty()) {
            path.merged_segment->left_extend =
                context.interior_segments.front()->left_extend;
            path.merged_segment->right_extend =
                context.interior_segments.back()->right_extend;
        }
        prepared.merged_block->anchors.emplace(
            context.key, path.merged_segment);
        prepared.paths.push_back(std::move(path));
    }

    rejection = MissingWindowReject::COUNT;
    return prepared;
}

struct MissingWindowCommitResult {
    size_t replaced_windows = 0;
    size_t replaced_old_blocks = 0;
    size_t blocks_before = 0;
    size_t blocks_after = 0;
};

MissingWindowCommitResult commitPreparedMissingWindows(
    RaMeshMultiGenomeGraph& graph,
    std::vector<MissingWindowCandidate>& accepted_candidates,
    std::vector<PreparedChain>& prepared,
    const std::unordered_set<const Block*>& accepted_old_blocks,
    size_t species_count,
    const char* phase,
    bool require_block_reduction) {
    if (prepared.empty()) {
        return {};
    }
    if (prepared.size() != accepted_candidates.size()) {
        throw std::runtime_error(
            "Species-mismatch prepared candidate vectors changed size");
    }

    size_t replaced_old_blocks = 0;
    std::unordered_map<const Block*, size_t> old_to_prepared;
    for (size_t index = 0; index < prepared.size(); ++index) {
        const auto& blocks_in_candidate =
            prepared[index].candidate.blocks;
        if (blocks_in_candidate.empty()) {
            throw std::runtime_error(
                "Species-mismatch prepared candidate has no old Blocks");
        }
        replaced_old_blocks += blocks_in_candidate.size();
        for (const auto& old_block : blocks_in_candidate) {
            if (!old_block ||
                !old_to_prepared.emplace(
                    old_block.get(), index).second) {
                throw std::runtime_error(
                    "Species-mismatch prepared candidate has invalid or "
                    "overlapping old Blocks");
            }
        }
    }
    if (replaced_old_blocks < prepared.size()) {
        throw std::runtime_error(
            "Species-mismatch replacement Block accounting underflow");
    }

    const size_t blocks_before = graph.blocks.size();
    const size_t expected_reduction =
        replaced_old_blocks - prepared.size();
    if (require_block_reduction && expected_reduction == 0) {
        throw std::runtime_error(
            "Zero-gap replacement did not reduce the Block pool");
    }

    std::vector<WeakBlock> replacement_pool;
    replacement_pool.reserve(graph.blocks.size());
    std::vector<size_t> old_occurrences(prepared.size(), 0);
    std::vector<bool> replacement_emitted(prepared.size(), false);
    std::unordered_set<const Block*> pool_blocks;
    for (const auto& weak_block : graph.blocks) {
        const auto block = weak_block.lock();
        if (!block) {
            replacement_pool.push_back(weak_block);
            continue;
        }
        if (!pool_blocks.insert(block.get()).second) {
            throw std::runtime_error(
                "Species-mismatch realignment found duplicate Block "
                "pool entries");
        }
        const auto old_it = old_to_prepared.find(block.get());
        if (old_it == old_to_prepared.end()) {
            replacement_pool.push_back(weak_block);
            continue;
        }
        ++old_occurrences[old_it->second];
        if (!replacement_emitted[old_it->second]) {
            replacement_pool.emplace_back(
                prepared[old_it->second].merged_block);
            replacement_emitted[old_it->second] = true;
        }
    }
    for (size_t index = 0; index < old_occurrences.size(); ++index) {
        if (!replacement_emitted[index] ||
            old_occurrences[index] !=
                prepared[index].candidate.blocks.size()) {
            throw std::runtime_error(
                "Species-mismatch old Blocks are inconsistent with "
                "the Block pool");
        }
    }

    std::map<SpeciesChrPair, SamplingSnapshot> sampling_snapshots;
    for (const auto& chain : prepared) {
        for (const auto& path : chain.paths) {
            if (sampling_snapshots.count(path.key) != 0) {
                continue;
            }
            auto& genome_end = genomeEndFor(graph, path.key);
            sampling_snapshots.emplace(
                path.key,
                SamplingSnapshot{&genome_end, genome_end.sample_vec});
        }
    }

    std::vector<AppliedSplice> applied_splices;
    applied_splices.reserve(prepared.size() * species_count);
    bool pool_replaced = false;
    try {
        for (size_t candidate_index = 0;
             candidate_index < prepared.size();
             ++candidate_index) {
            auto& chain = prepared[candidate_index];
            const auto& candidate =
                accepted_candidates[candidate_index];
            if (chain.paths.size() != candidate.paths.size()) {
                throw std::runtime_error(
                    "Species-mismatch path preparation changed size");
            }
            for (size_t path_index = 0;
                 path_index < chain.paths.size();
                 ++path_index) {
                auto& path = chain.paths[path_index];
                const auto& context = candidate.paths[path_index];
                const bool reverse =
                    path.merged_segment->strand == Strand::REVERSE;
                SegPtr old_first;
                SegPtr old_last;
                SegPtr previous;
                SegPtr next;
                if (path.old_segments.empty()) {
                    previous = reverse ? context.right_anchor
                                       : context.left_anchor;
                    next = reverse ? context.left_anchor
                                   : context.right_anchor;
                } else {
                    old_first = reverse
                                    ? path.old_segments.back()
                                    : path.old_segments.front();
                    old_last = reverse
                                   ? path.old_segments.front()
                                   : path.old_segments.back();
                    previous = old_first->primary_path.prev.load(
                        std::memory_order_acquire);
                    next = old_last->primary_path.next.load(
                        std::memory_order_acquire);
                    for (size_t old_index = 1;
                         old_index < path.old_segments.size();
                         ++old_index) {
                        if (!segmentPathsAreDirectlyAdjacent(
                                path.old_segments[old_index - 1],
                                path.old_segments[old_index])) {
                            throw std::runtime_error(
                                "Species-mismatch old path is no longer "
                                "contiguous");
                        }
                    }
                }

                if (!previous || !next ||
                    previous->primary_path.next.load(
                        std::memory_order_acquire) !=
                        (old_first ? old_first : next) ||
                    next->primary_path.prev.load(
                        std::memory_order_acquire) !=
                        (old_last ? old_last : previous)) {
                    throw std::runtime_error(
                        "Species-mismatch path changed before commit");
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

        graph.blocks.swap(replacement_pool);
        pool_replaced = true;
        if (graph.blocks.size() + expected_reduction != blocks_before) {
            throw std::runtime_error(
                "Species-mismatch replacement changed the Block pool by "
                "an unexpected amount");
        }
        for (auto& [key, snapshot] : sampling_snapshots) {
            (void)key;
            rebuildSampling(*snapshot.genome_end);
        }

        std::atomic_thread_fence(std::memory_order_seq_cst);
        if (!auditAffectedGraph(
                graph, sampling_snapshots, prepared,
                accepted_old_blocks)) {
            throw std::runtime_error(
                "Species-mismatch post-commit graph audit failed");
        }
    } catch (...) {
        if (pool_replaced) {
            graph.blocks.swap(replacement_pool);
        }
        for (auto splice = applied_splices.rbegin();
             splice != applied_splices.rend();
             ++splice) {
            if (splice->old_first) {
                splice->previous->primary_path.next.store(
                    splice->old_first, std::memory_order_release);
                splice->old_first->primary_path.prev.store(
                    splice->previous, std::memory_order_release);
                splice->old_last->primary_path.next.store(
                    splice->next, std::memory_order_release);
                splice->next->primary_path.prev.store(
                    splice->old_last, std::memory_order_release);
            } else {
                splice->previous->primary_path.next.store(
                    splice->next, std::memory_order_release);
                splice->next->primary_path.prev.store(
                    splice->previous, std::memory_order_release);
            }
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
            "[species-mismatch-realign][{}] commit failed; phase "
            "changes rolled back",
            phase);
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

    return MissingWindowCommitResult{
        prepared.size(), replaced_old_blocks,
        blocks_before, graph.blocks.size()};
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

size_t RaMeshMultiGenomeGraph::realignSingleMissingSpeciesWindows(
    const SpeciesName& reference_species,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>&
        seqpro_managers,
    const std::string& msa_executable,
    uint_t maximum_span,
    uint_t parallel_threads,
    uint_t zero_gap_maximum_span) {
    if (maximum_span == 0) {
        throw std::invalid_argument(
            "species-mismatch realignment maximum_span must be positive");
    }
    if (zero_gap_maximum_span == 0) {
        throw std::invalid_argument(
            "species-mismatch zero-gap maximum span must be positive");
    }
    if (seqpro_managers.size() < 3) {
        spdlog::info(
            "[species-mismatch-realign] reference={} skipped: "
            "at least three species are required",
            reference_species);
        return 0;
    }

    std::unique_lock graph_lock(rw);
    const auto reference_graph_it =
        species_graphs.find(reference_species);
    if (reference_graph_it == species_graphs.end()) {
        throw std::runtime_error(
            "Species-mismatch realignment reference is absent: " +
            reference_species);
    }

    constexpr size_t rejection_count =
        static_cast<size_t>(MissingWindowReject::COUNT);
    using RejectionCounts = std::array<size_t, rejection_count>;

    RejectionCounts zero_gap_rejections{};
    std::map<SpeciesName, size_t> zero_gap_accepted_by_missing_species;
    std::map<std::string, size_t> zero_gap_pattern_counts;
    size_t zero_gap_scan_count = 0;
    size_t zero_gap_candidate_events = 0;
    size_t zero_gap_prepared_events = 0;
    size_t zero_gap_overlap_events = 0;
    size_t zero_gap_replaced_windows = 0;
    size_t zero_gap_replaced_old_blocks = 0;
    double zero_gap_prepare_seconds = 0.0;
    const size_t zero_gap_iteration_limit = blocks.size() + 1;

    while (true) {
        ++zero_gap_scan_count;
        if (zero_gap_scan_count > zero_gap_iteration_limit) {
            throw std::runtime_error(
                "Zero-gap fixed-point iteration exceeded the initial "
                "Block-pool bound");
        }

        size_t scanned_this_iteration = 0;
        size_t structural_this_iteration = 0;
        std::vector<MissingWindowCandidate> zero_gap_candidates;
        for (auto& [chromosome, genome_end] :
             reference_graph_it->second.chr2end) {
            auto left_reference =
                genome_end.head->primary_path.next.load(
                    std::memory_order_acquire);
            while (left_reference && !left_reference->isTail()) {
                const auto middle_reference =
                    left_reference->primary_path.next.load(
                        std::memory_order_acquire);
                if (!middle_reference || middle_reference->isTail()) {
                    break;
                }
                const auto right_reference =
                    middle_reference->primary_path.next.load(
                        std::memory_order_acquire);
                if (!right_reference || right_reference->isTail()) {
                    break;
                }

                ++scanned_this_iteration;
                MissingWindowCandidate candidate;
                const auto candidate_result =
                    buildMissingWindowCandidate(
                        left_reference->parent_block,
                        std::vector<BlockPtr>{
                            middle_reference->parent_block},
                        right_reference->parent_block,
                        reference_species, chromosome,
                        seqpro_managers.size(), maximum_span,
                        zero_gap_maximum_span, true, candidate);
                if (candidate_result != MissingWindowReject::COUNT) {
                    ++zero_gap_rejections[
                        static_cast<size_t>(candidate_result)];
                    left_reference = middle_reference;
                    continue;
                }

                ++structural_this_iteration;
                if (candidate.zero_gap_deletion) {
                    zero_gap_candidates.push_back(std::move(candidate));
                }
                left_reference = right_reference;
            }
        }
        zero_gap_candidate_events += zero_gap_candidates.size();

        if (zero_gap_candidates.empty()) {
            spdlog::info(
                "[species-mismatch-realign][zero-gap] iteration={} "
                "scanned_triples={} structural_candidates={} candidates=0 "
                "status=stable",
                zero_gap_scan_count, scanned_this_iteration,
                structural_this_iteration);
            break;
        }

        const uint_t zero_gap_threads = static_cast<uint_t>(
            std::max<size_t>(
                1,
                std::min<size_t>(
                    std::max<uint_t>(1, parallel_threads),
                    zero_gap_candidates.size())));
        std::vector<std::optional<PreparedChain>> prepared_slots(
            zero_gap_candidates.size());
        std::vector<MissingWindowReject> prepare_results(
            zero_gap_candidates.size(), MissingWindowReject::COUNT);
        std::vector<std::exception_ptr> prepare_exceptions(
            zero_gap_candidates.size());

        const auto prepare_start = std::chrono::steady_clock::now();
#pragma omp parallel for schedule(dynamic, 1) num_threads(zero_gap_threads)
        for (std::int64_t candidate_index = 0;
             candidate_index < static_cast<std::int64_t>(
                                   zero_gap_candidates.size());
             ++candidate_index) {
            const size_t index = static_cast<size_t>(candidate_index);
            try {
                MissingWindowReject prepare_result =
                    MissingWindowReject::COUNT;
                auto prepared_candidate =
                    prepareZeroGapDeletionWindow(
                        zero_gap_candidates[index], reference_species,
                        prepare_result);
                if (!prepared_candidate.has_value() &&
                    prepare_result == MissingWindowReject::COUNT) {
                    prepare_result =
                        MissingWindowReject::ZERO_GAP_PREPARATION_INVALID;
                }
                prepare_results[index] = prepare_result;
                prepared_slots[index] = std::move(prepared_candidate);
            } catch (...) {
                prepare_exceptions[index] = std::current_exception();
            }
        }
        const double prepare_seconds =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - prepare_start)
                .count();
        zero_gap_prepare_seconds += prepare_seconds;

        for (size_t index = 0;
             index < prepare_exceptions.size();
             ++index) {
            if (!prepare_exceptions[index]) {
                continue;
            }
            for (auto& slot : prepared_slots) {
                if (slot.has_value()) {
                    detachPreparedChain(*slot);
                }
            }
            spdlog::error(
                "[species-mismatch-realign][zero-gap] parallel "
                "preparation failed at candidate_index={}",
                index);
            std::rethrow_exception(prepare_exceptions[index]);
        }

        size_t prepared_this_iteration = 0;
        size_t overlap_this_iteration = 0;
        std::vector<MissingWindowCandidate> accepted_candidates;
        std::vector<PreparedChain> prepared;
        std::unordered_set<const Block*> accepted_old_blocks;
        std::unordered_set<const Block*> reserved_read_blocks;
        accepted_candidates.reserve(zero_gap_candidates.size());
        prepared.reserve(zero_gap_candidates.size());
        accepted_old_blocks.reserve(zero_gap_candidates.size() * 3);
        reserved_read_blocks.reserve(zero_gap_candidates.size() * 3);

        for (size_t index = 0;
             index < zero_gap_candidates.size();
             ++index) {
            if (!prepared_slots[index].has_value()) {
                auto prepare_result = prepare_results[index];
                if (prepare_result == MissingWindowReject::COUNT) {
                    prepare_result =
                        MissingWindowReject::ZERO_GAP_PREPARATION_INVALID;
                }
                ++zero_gap_rejections[
                    static_cast<size_t>(prepare_result)];
                continue;
            }
            ++prepared_this_iteration;

            auto& candidate = zero_gap_candidates[index];
            const std::array<const Block*, 3> read_blocks = {
                candidate.left.block.get(),
                candidate.interiors.front().block.get(),
                candidate.right.block.get()};
            const bool overlaps = std::any_of(
                read_blocks.begin(), read_blocks.end(),
                [&](const Block* block) {
                    return block && reserved_read_blocks.count(block) != 0;
                });
            if (overlaps) {
                ++overlap_this_iteration;
                ++zero_gap_rejections[static_cast<size_t>(
                    MissingWindowReject::OVERLAPPING_PREPARED_WINDOW)];
                detachPreparedChain(*prepared_slots[index]);
                prepared_slots[index].reset();
                continue;
            }

            const auto& replaced_blocks =
                prepared_slots[index]->candidate.blocks;
            for (const auto& replaced_block : replaced_blocks) {
                if (!replaced_block ||
                    !accepted_old_blocks.insert(
                        replaced_block.get()).second) {
                    for (auto& slot : prepared_slots) {
                        if (slot.has_value()) {
                            detachPreparedChain(*slot);
                        }
                    }
                    detachPreparedBlocks(prepared);
                    throw std::runtime_error(
                        "Zero-gap selection produced overlapping "
                        "replacement Blocks");
                }
            }
            for (const auto* block : read_blocks) {
                if (block) {
                    reserved_read_blocks.insert(block);
                }
            }
            for (const auto& species : candidate.missing_species) {
                ++zero_gap_accepted_by_missing_species[species];
            }
            ++zero_gap_pattern_counts[candidate.participant_pattern];
            accepted_candidates.push_back(std::move(candidate));
            prepared.push_back(std::move(*prepared_slots[index]));
        }

        zero_gap_prepared_events += prepared_this_iteration;
        zero_gap_overlap_events += overlap_this_iteration;
        if (prepared.empty()) {
            spdlog::warn(
                "[species-mismatch-realign][zero-gap] iteration={} "
                "scanned_triples={} structural_candidates={} candidates={} "
                "prepared={} selected=0 failed={} "
                "overlap={} parallel_threads={} wall_seconds={:.3f} "
                "status=stalled",
                zero_gap_scan_count, scanned_this_iteration,
                structural_this_iteration, zero_gap_candidates.size(),
                prepared_this_iteration,
                zero_gap_candidates.size() - prepared_this_iteration,
                overlap_this_iteration, zero_gap_threads,
                prepare_seconds);
            break;
        }

        const auto commit = commitPreparedMissingWindows(
            *this, accepted_candidates, prepared, accepted_old_blocks,
            seqpro_managers.size(), "zero-gap", true);
        zero_gap_replaced_windows += commit.replaced_windows;
        zero_gap_replaced_old_blocks += commit.replaced_old_blocks;
        spdlog::info(
            "[species-mismatch-realign][zero-gap] iteration={} "
            "scanned_triples={} structural_candidates={} candidates={} "
            "prepared={} selected={} failed={} "
            "overlap={} parallel_threads={} wall_seconds={:.3f} "
            "replaced_old_blocks={} blocks_before={} blocks_after={} "
            "status=committed",
            zero_gap_scan_count, scanned_this_iteration,
            structural_this_iteration, zero_gap_candidates.size(),
            prepared_this_iteration, commit.replaced_windows,
            zero_gap_candidates.size() - prepared_this_iteration,
            overlap_this_iteration, zero_gap_threads, prepare_seconds,
            commit.replaced_old_blocks, commit.blocks_before,
            commit.blocks_after);
    }

    std::array<size_t, rejection_count> rejections{};
    size_t scanned_triples = 0;
    size_t structural_candidates = 0;
    std::map<SpeciesName, size_t> accepted_by_missing_species =
        zero_gap_accepted_by_missing_species;
    std::map<std::string, size_t> minipoa_pattern_counts;
    std::vector<MissingWindowCandidate> candidate_slots;
    std::vector<MissingWindowCandidate> accepted_candidates;
    std::vector<PreparedChain> prepared;
    std::unordered_set<const Block*> accepted_old_blocks;

    for (auto& [chromosome, genome_end] :
         reference_graph_it->second.chr2end) {
        auto left_reference =
            genome_end.head->primary_path.next.load(
                std::memory_order_acquire);
        while (left_reference && !left_reference->isTail()) {
            BlockView left_view;
            if (!buildBlockView(
                    left_reference->parent_block, reference_species,
                    left_view) ||
                left_view.species_count != seqpro_managers.size()) {
                left_reference = left_reference->primary_path.next.load(
                    std::memory_order_acquire);
                continue;
            }

            std::vector<BlockPtr> interior_blocks;
            auto cursor = left_reference->primary_path.next.load(
                std::memory_order_acquire);
            SegPtr right_reference;
            while (cursor && !cursor->isTail()) {
                BlockView cursor_view;
                if (!buildRealignBlockView(
                        cursor->parent_block, reference_species,
                        cursor_view)) {
                    interior_blocks.push_back(cursor->parent_block);
                    cursor = cursor->primary_path.next.load(
                        std::memory_order_acquire);
                    continue;
                }
                if (cursor_view.species_count ==
                    seqpro_managers.size()) {
                    right_reference = cursor;
                    break;
                }
                interior_blocks.push_back(cursor->parent_block);
                cursor = cursor->primary_path.next.load(
                    std::memory_order_acquire);
            }
            if (!right_reference) {
                break;
            }
            if (interior_blocks.empty()) {
                left_reference = right_reference;
                continue;
            }

            ++scanned_triples;
            MissingWindowCandidate candidate;
            const auto candidate_result =
                buildMissingWindowCandidate(
                    left_reference->parent_block, interior_blocks,
                    right_reference->parent_block,
                    reference_species, chromosome,
                    seqpro_managers.size(), maximum_span,
                    zero_gap_maximum_span, false, candidate);
            if (candidate_result != MissingWindowReject::COUNT) {
                ++rejections[static_cast<size_t>(candidate_result)];
                left_reference = right_reference;
                continue;
            }

            ++structural_candidates;
            candidate_slots.push_back(std::move(candidate));
            left_reference = right_reference;
        }
    }

    const uint_t effective_threads = static_cast<uint_t>(
        std::max<size_t>(
            1,
            std::min<size_t>(
                std::max<uint_t>(1, parallel_threads),
                candidate_slots.empty() ? 1 : candidate_slots.size())));
    std::vector<std::optional<PreparedChain>> prepared_slots(
        candidate_slots.size());
    std::vector<MissingWindowReject> prepare_results(
        candidate_slots.size(), MissingWindowReject::COUNT);
    std::vector<std::exception_ptr> prepare_exceptions(
        candidate_slots.size());

    const auto prepare_start = std::chrono::steady_clock::now();
#pragma omp parallel for schedule(dynamic, 1) num_threads(effective_threads)
    for (std::int64_t candidate_index = 0;
         candidate_index <
             static_cast<std::int64_t>(candidate_slots.size());
         ++candidate_index) {
        const size_t index = static_cast<size_t>(candidate_index);
        try {
            MissingWindowReject prepare_result =
                MissingWindowReject::COUNT;
            auto prepared_candidate = prepareMissingWindow(
                candidate_slots[index], reference_species,
                seqpro_managers, msa_executable, prepare_result);
            if (!prepared_candidate.has_value() &&
                prepare_result == MissingWindowReject::COUNT) {
                prepare_result = MissingWindowReject::MSA_INVALID;
            }
            prepare_results[index] = prepare_result;
            prepared_slots[index] = std::move(prepared_candidate);
        } catch (...) {
            prepare_exceptions[index] = std::current_exception();
        }
    }
    const double prepare_seconds =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - prepare_start)
            .count();

    for (size_t index = 0; index < prepare_exceptions.size(); ++index) {
        if (!prepare_exceptions[index]) {
            continue;
        }
        for (auto& slot : prepared_slots) {
            if (slot.has_value()) {
                detachPreparedChain(*slot);
            }
        }
        spdlog::error(
            "[species-mismatch-realign][minipoa] parallel preparation "
            "failed at candidate_index={}",
            index);
        std::rethrow_exception(prepare_exceptions[index]);
    }

    const size_t attempted_msa = candidate_slots.size();
    size_t successful_msa = 0;
    for (size_t index = 0; index < candidate_slots.size(); ++index) {
        if (prepared_slots[index].has_value()) {
            ++successful_msa;
        }
    }

    accepted_candidates.reserve(candidate_slots.size());
    prepared.reserve(candidate_slots.size());
    for (size_t index = 0; index < candidate_slots.size(); ++index) {
        if (!prepared_slots[index].has_value()) {
            auto prepare_result = prepare_results[index];
            if (prepare_result == MissingWindowReject::COUNT) {
                prepare_result = MissingWindowReject::MSA_INVALID;
            }
            ++rejections[static_cast<size_t>(prepare_result)];
            continue;
        }

        auto& candidate = candidate_slots[index];
        const auto& replaced_blocks =
            prepared_slots[index]->candidate.blocks;
        for (const auto& replaced_block : replaced_blocks) {
            if (!replaced_block ||
                !accepted_old_blocks.insert(
                    replaced_block.get()).second) {
                for (auto& slot : prepared_slots) {
                    if (slot.has_value()) {
                        detachPreparedChain(*slot);
                    }
                }
                detachPreparedBlocks(prepared);
                throw std::runtime_error(
                    "Species-mismatch realignment produced overlapping "
                    "replacement Blocks");
            }
        }
        for (const auto& species : candidate.missing_species) {
            ++accepted_by_missing_species[species];
        }
        ++minipoa_pattern_counts[candidate.participant_pattern];
        accepted_candidates.push_back(std::move(candidate));
        prepared.push_back(std::move(*prepared_slots[index]));
    }

    const size_t failed_msa = attempted_msa - successful_msa;
    const double windows_per_second =
        prepare_seconds > 0.0
            ? static_cast<double>(attempted_msa) / prepare_seconds
            : 0.0;
    spdlog::info(
        "[species-mismatch-realign][minipoa] reference={} "
        "rescanned_triples={} structural_candidates={} "
        "attempted={} succeeded={} failed={} "
        "selected={} parallel_threads={} wall_seconds={:.3f} "
        "windows_per_second={:.3f} maximum_span={}",
        reference_species, scanned_triples, structural_candidates,
        attempted_msa, successful_msa, failed_msa, prepared.size(),
        effective_threads, prepare_seconds, windows_per_second,
        maximum_span);
    spdlog::info(
        "[species-mismatch-realign][zero-gap] reference={} scans={} "
        "candidate_events={} prepared_events={} selected={} "
        "overlap_events={} replaced_old_blocks={} "
        "prepare_wall_seconds={:.3f} zero_gap_maximum_span={}",
        reference_species, zero_gap_scan_count,
        zero_gap_candidate_events, zero_gap_prepared_events,
        zero_gap_replaced_windows, zero_gap_overlap_events,
        zero_gap_replaced_old_blocks, zero_gap_prepare_seconds,
        zero_gap_maximum_span);

    for (size_t index = 0; index < rejection_count; ++index) {
        if (zero_gap_rejections[index] != 0) {
            spdlog::info(
                "[species-mismatch-realign][zero-gap] "
                "reject_reason={} count={}",
                missingWindowRejectName(
                    static_cast<MissingWindowReject>(index)),
                zero_gap_rejections[index]);
        }
        if (rejections[index] != 0) {
            spdlog::info(
                "[species-mismatch-realign][minipoa] "
                "reject_reason={} count={}",
                missingWindowRejectName(
                    static_cast<MissingWindowReject>(index)),
                rejections[index]);
        }
    }
    for (const auto& [species, count] :
         accepted_by_missing_species) {
        spdlog::info(
            "[species-mismatch-realign] missing_species={} committed={}",
            species, count);
    }
    for (const auto& [pattern, count] : zero_gap_pattern_counts) {
        spdlog::info(
            "[species-mismatch-realign][zero-gap] pattern={} committed={}",
            pattern, count);
    }
    for (const auto& [pattern, count] : minipoa_pattern_counts) {
        spdlog::info(
            "[species-mismatch-realign][minipoa] pattern={} committed={}",
            pattern, count);
    }

    if (prepared.empty()) {
        spdlog::info(
            "[species-mismatch-realign] reference={} replaced_windows={} "
            "zero_gap_deletion_windows={} minipoa_windows=0 "
            "replaced_old_blocks={} new_full_species_blocks={}",
            reference_species, zero_gap_replaced_windows,
            zero_gap_replaced_windows, zero_gap_replaced_old_blocks,
            zero_gap_replaced_windows);
        return zero_gap_replaced_windows;
    }

    const auto minipoa_commit = commitPreparedMissingWindows(
        *this, accepted_candidates, prepared, accepted_old_blocks,
        seqpro_managers.size(), "minipoa", false);
    const size_t total_replaced_windows =
        zero_gap_replaced_windows + minipoa_commit.replaced_windows;
    const size_t total_replaced_old_blocks =
        zero_gap_replaced_old_blocks +
        minipoa_commit.replaced_old_blocks;
    spdlog::info(
        "[species-mismatch-realign] reference={} replaced_windows={} "
        "zero_gap_deletion_windows={} minipoa_windows={} "
        "replaced_old_blocks={} new_full_species_blocks={}",
        reference_species, total_replaced_windows,
        zero_gap_replaced_windows, minipoa_commit.replaced_windows,
        total_replaced_old_blocks, total_replaced_windows);
    return total_replaced_windows;
}

}  // namespace RaMesh
