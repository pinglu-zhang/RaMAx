#include "ramesh.h"
#include "align.h"
#include "realignment_internal.h"

#include <spdlog/spdlog.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <exception>
#include <fstream>
#include <numeric>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
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

using BlockView = Realignment::BlockView;

using BlockViewCache = Realignment::BlockViewBuilder;

struct CandidateChain {
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

constexpr size_t kExactBatchMaximumChains = 4096;
constexpr size_t kExactBatchMaximumAnchors = 262144;
constexpr size_t kZeroGapBatchMaximumWindows = 4096;
constexpr size_t kZeroGapBatchMaximumAnchors = 262144;
constexpr size_t kMinipoaBatchMaximumWindows = 1024;
constexpr size_t kMinipoaBatchMaximumAnchors = 131072;
constexpr size_t kMaximumPreparationThreads = 32;

thread_local const char* g_allocation_phase = "idle";
thread_local const char* g_last_allocation_phase = "idle";
thread_local size_t g_allocation_candidates = 0;
thread_local size_t g_allocation_prepared = 0;
thread_local size_t g_allocation_batch_size = 0;

class AllocationPhaseScope {
public:
    explicit AllocationPhaseScope(const char* phase) noexcept
        : previous_(g_allocation_phase) {
        g_allocation_phase = phase;
        g_last_allocation_phase = phase;
    }
    AllocationPhaseScope(const AllocationPhaseScope&) = delete;
    AllocationPhaseScope& operator=(const AllocationPhaseScope&) = delete;
    ~AllocationPhaseScope() { g_allocation_phase = previous_; }

private:
    const char* previous_;
};

struct ProcessMemorySnapshot {
    uint64_t rss_kib = 0;
    uint64_t peak_kib = 0;
    uint64_t virtual_kib = 0;
    uint64_t cgroup_limit_bytes = 0;
};

ProcessMemorySnapshot readProcessMemorySnapshot() noexcept {
    ProcessMemorySnapshot result;
    try {
        std::ifstream status("/proc/self/status");
        std::string line;
        while (std::getline(status, line)) {
            std::istringstream fields(line);
            std::string key;
            uint64_t value = 0;
            std::string unit;
            if (!(fields >> key >> value >> unit)) continue;
            if (key == "VmRSS:") result.rss_kib = value;
            else if (key == "VmPeak:") result.peak_kib = value;
            else if (key == "VmSize:") result.virtual_kib = value;
        }
        for (const char* path : {
                 "/sys/fs/cgroup/memory.max",
                 "/sys/fs/cgroup/memory/memory.limit_in_bytes"}) {
            std::ifstream limit_file(path);
            std::string limit;
            if (limit_file >> limit && limit != "max") {
                result.cgroup_limit_bytes = std::stoull(limit);
                if (result.cgroup_limit_bytes >= (1ULL << 60U)) {
                    result.cgroup_limit_bytes = 0;
                }
                break;
            }
        }
    } catch (...) {
    }
    return result;
}

void logAllocationFailure(const char* module,
                          const SpeciesName& reference_species,
                          size_t blocks,
                          size_t candidates,
                          size_t prepared,
                          size_t batch_size) noexcept {
    const auto memory = readProcessMemorySnapshot();
    try {
        spdlog::error(
            "[{}][bad-alloc] reference={} phase={} blocks={} "
            "candidates={} prepared={} batch_size={} vm_rss_kib={} "
            "vm_peak_kib={} vm_size_kib={} cgroup_limit_bytes={}",
            module, reference_species, g_last_allocation_phase, blocks,
            candidates, prepared, batch_size, memory.rss_kib,
            memory.peak_kib, memory.virtual_kib,
            memory.cgroup_limit_bytes);
    } catch (...) {
    }
}

class BlockPoolEditor {
public:
    using AddedSlotMap = std::unordered_map<const Block*, size_t>;

    struct SlotChange {
        size_t slot = 0;
        WeakBlock block;
        bool removed = false;
    };

    struct Change {
        std::vector<SlotChange> slots;
        std::vector<AddedSlotMap::node_type> restored_added;
        std::vector<const Block*> added;
        size_t active_before = 0;
        size_t removed_before = 0;
    };

    explicit BlockPoolEditor(std::vector<WeakBlock>& pool)
        : pool_(pool), removed_slots_(pool.size(), false) {
        initial_slots_.reserve(pool_.size());
        for (size_t slot = 0; slot < pool_.size(); ++slot) {
            if (const auto block = pool_[slot].lock()) {
                initial_slots_.emplace_back(block.get(), slot);
            }
        }
        std::sort(
            initial_slots_.begin(), initial_slots_.end(),
            [](const InitialSlot& left, const InitialSlot& right) {
                return std::less<const Block*>{}(
                    left.first, right.first);
            });
        for (size_t index = 1; index < initial_slots_.size(); ++index) {
            if (initial_slots_[index - 1].first ==
                initial_slots_[index].first) {
                throw std::runtime_error(
                    "Block pool contains a duplicate active Block");
            }
        }
        active_count_ = initial_slots_.size();
    }

    [[nodiscard]] size_t activeCount() const noexcept {
        return active_count_;
    }

    [[nodiscard]] bool contains(const Block* block) const {
        if (!block) return false;
        if (added_slots_.count(block) != 0) return true;
        const auto initial = findInitial(block);
        if (initial == initial_slots_.end() || initial->first != block) {
            return false;
        }
        const auto active = pool_[initial->second].lock();
        return active && active.get() == block;
    }

    Change replace(const std::vector<PreparedChain>& prepared) {
        Change change;
        change.active_before = active_count_;
        change.removed_before = removed_count_;
        size_t changed_slots = 0;
        for (const auto& chain : prepared) {
            changed_slots += chain.candidate.blocks.size();
        }
        change.slots.reserve(changed_slots);
        change.restored_added.reserve(changed_slots);
        change.added.reserve(prepared.size());
        added_slots_.reserve(added_slots_.size() + prepared.size());

        try {
            for (const auto& chain : prepared) {
                if (!chain.merged_block || chain.candidate.blocks.empty()) {
                    throw std::runtime_error(
                        "Prepared replacement has no Block payload");
                }
                std::vector<size_t> slots;
                slots.reserve(chain.candidate.blocks.size());
                for (const auto& old_block : chain.candidate.blocks) {
                    const auto slot = activeSlot(old_block.get());
                    if (!slot.has_value()) {
                        throw std::runtime_error(
                            "Prepared replacement Block is absent from pool");
                    }
                    slots.push_back(*slot);
                }
                std::sort(slots.begin(), slots.end());
                if (std::adjacent_find(slots.begin(), slots.end()) !=
                    slots.end()) {
                    throw std::runtime_error(
                        "Prepared replacement repeats a pool slot");
                }
                const size_t replacement_slot = slots.front();
                for (const size_t slot : slots) {
                    change.slots.push_back(
                        SlotChange{slot, pool_[slot], removed_slots_[slot]});
                    pool_[slot].reset();
                    if (!removed_slots_[slot]) {
                        removed_slots_[slot] = true;
                        ++removed_count_;
                    }
                }
                pool_[replacement_slot] = chain.merged_block;
                if (removed_slots_[replacement_slot]) {
                    removed_slots_[replacement_slot] = false;
                    --removed_count_;
                }

                for (size_t index = 0;
                     index < chain.candidate.blocks.size(); ++index) {
                    const Block* old =
                        chain.candidate.blocks[index].get();
                    const auto added = added_slots_.find(old);
                    if (added != added_slots_.end()) {
                        change.restored_added.push_back(
                            added_slots_.extract(added));
                    }
                }
                change.added.push_back(chain.merged_block.get());
                added_slots_.emplace(
                    chain.merged_block.get(), replacement_slot);
                active_count_ -= chain.candidate.blocks.size() - 1;
            }
        } catch (...) {
            rollback(change);
            throw;
        }
        return change;
    }

    void rollback(Change& change) noexcept {
        for (const Block* block : change.added) {
            added_slots_.erase(block);
        }
        for (auto& entry : change.restored_added) {
            added_slots_.insert(std::move(entry));
        }
        for (auto iterator = change.slots.rbegin();
             iterator != change.slots.rend(); ++iterator) {
            pool_[iterator->slot] = iterator->block;
            removed_slots_[iterator->slot] = iterator->removed;
        }
        active_count_ = change.active_before;
        removed_count_ = change.removed_before;
        change = {};
    }

    void compact() {
        const size_t original_size = pool_.size();
        size_t write = 0;
        for (size_t read = 0; read < original_size; ++read) {
            if (removed_slots_[read]) continue;
            if (write != read) pool_[write] = std::move(pool_[read]);
            ++write;
        }
        pool_.resize(write);
        if (removed_count_ > original_size ||
            write + removed_count_ != original_size) {
            throw std::runtime_error(
                "Block pool compaction changed the slot count");
        }
        removed_slots_.clear();
        removed_count_ = 0;
    }

    [[nodiscard]] bool auditPool() const {
        std::vector<const Block*> active;
        active.reserve(active_count_);
        for (const auto& weak_block : pool_) {
            if (const auto block = weak_block.lock()) {
                active.push_back(block.get());
            }
        }
        std::sort(
            active.begin(), active.end(), std::less<const Block*>{});
        return active.size() == active_count_ &&
            std::adjacent_find(active.begin(), active.end()) == active.end();
    }

private:
    using InitialSlot = std::pair<const Block*, size_t>;

    [[nodiscard]] std::vector<InitialSlot>::const_iterator findInitial(
        const Block* block) const {
        return std::lower_bound(
            initial_slots_.begin(), initial_slots_.end(),
            InitialSlot{block, 0},
            [](const InitialSlot& left, const InitialSlot& right) {
                return std::less<const Block*>{}(
                    left.first, right.first);
            });
    }

    [[nodiscard]] std::optional<size_t> activeSlot(
        const Block* block) const {
        if (!block) return std::nullopt;
        const auto added = added_slots_.find(block);
        if (added != added_slots_.end()) return added->second;
        const auto initial = findInitial(block);
        if (initial == initial_slots_.end() || initial->first != block) {
            return std::nullopt;
        }
        const auto active = pool_[initial->second].lock();
        if (!active || active.get() != block) return std::nullopt;
        return initial->second;
    }

    std::vector<WeakBlock>& pool_;
    std::vector<InitialSlot> initial_slots_;
    AddedSlotMap added_slots_;
    std::vector<bool> removed_slots_;
    size_t active_count_ = 0;
    size_t removed_count_ = 0;
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
    Realignment::BlockViewBuilder builder(reference_species);
    return builder.build(block, Realignment::BlockViewProfile::ExactMerge,
                         view);
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

    size_t positive_query_gap_count = 0;
    for (const auto& [key, left_segment] : left.anchors) {
        const auto right_it = right.anchors.find(key);
        const uint_t allowed_gap =
            key == reference_key ? 0 : maximum_query_gap;
        if (right_it != right.anchors.end() && key != reference_key) {
            const int64_t gap = segmentGapInReferenceOrder(
                left_segment, right_it->second);
            if (gap > 0) {
                ++positive_query_gap_count;
            }
        }
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

    return positive_query_gap_count <= 1;
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
    if (candidate.blocks.size() < 2) {
        throw std::runtime_error("Invalid exact-contiguous candidate chain");
    }
    std::vector<BlockView> views;
    views.reserve(candidate.blocks.size());
    for (const auto& block : candidate.blocks) {
        BlockView view;
        if (!buildBlockView(block, reference_species, view)) {
            throw std::runtime_error(
                "Exact-contiguous candidate changed before preparation");
        }
        views.push_back(std::move(view));
    }

    PreparedChain prepared;
    prepared.candidate.blocks = candidate.blocks;
    prepared.candidate.reference_chromosome =
        candidate.reference_chromosome;
    prepared.candidate.reference_start = candidate.reference_start;
    prepared.candidate.reference_length = candidate.reference_length;
    prepared.candidate.species_count = candidate.species_count;
    prepared.merged_block =
        Block::createEmpty(candidate.reference_chromosome,
                           candidate.species_count);

    const auto& keys = views.front().anchors;
    for (const auto& [key, unused] : keys) {
        (void)unused;
        PathReplacement path;
        path.key = key;
        path.old_segments.reserve(views.size());

        uint64_t merged_start = std::numeric_limits<uint64_t>::max();
        uint64_t merged_end = 0;
        uint64_t merged_length = 0;
        Cigar_t merged_cigar;
        Strand strand = Strand::FORWARD;

        for (size_t index = 0; index < views.size(); ++index) {
            const auto& view = views[index];
            const auto segment_it = view.anchors.find(key);
            if (segment_it == view.anchors.end()) {
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
                    views[index - 1].anchors.find(key);
                if (previous_it ==
                    views[index - 1].anchors.end()) {
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
                    segment, view.reference_segment->length, normalized) ||
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

bool rebuildSamplingAndAuditAffectedGraph(
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
        auto* genome_end = snapshot.genome_end;
        if (!genome_end || !genome_end->head || !genome_end->tail) {
            return false;
        }

        genome_end->sample_vec.clear();
        genome_end->sample_vec.resize(1, genome_end->head);

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

            genome_end->setToSampling(current);

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

bool rebuildSamplingAndAuditAffectedGraphBounded(
    RaMeshMultiGenomeGraph& graph,
    const std::map<SpeciesChrPair, SamplingSnapshot>& affected,
    const std::vector<PreparedChain>& prepared,
    const BlockPoolEditor& pool_editor) {
    std::unordered_set<const Segment*> expected_new_segments;
    expected_new_segments.reserve(
        std::accumulate(
            prepared.begin(), prepared.end(), size_t{0},
            [](size_t total, const PreparedChain& chain) {
                return total + chain.paths.size();
            }));
    for (const auto& chain : prepared) {
        if (!chain.merged_block ||
            !pool_editor.contains(chain.merged_block.get()) ||
            chain.merged_block->anchors.size() != chain.paths.size()) {
            return false;
        }
        for (const auto& path : chain.paths) {
            const auto anchor = chain.merged_block->anchors.find(path.key);
            if (!path.merged_segment ||
                path.merged_segment->parent_block != chain.merged_block ||
                anchor == chain.merged_block->anchors.end() ||
                anchor->second != path.merged_segment ||
                !expected_new_segments.insert(
                    path.merged_segment.get()).second) {
                return false;
            }
        }
    }

    std::unordered_set<const Segment*> seen_new_segments;
    seen_new_segments.reserve(expected_new_segments.size());
    for (const auto& [key, snapshot] : affected) {
        auto* genome_end = snapshot.genome_end;
        if (!genome_end || !genome_end->head || !genome_end->tail) {
            return false;
        }

        genome_end->sample_vec.clear();
        genome_end->sample_vec.resize(1, genome_end->head);
        std::unordered_set<const Segment*> path_segments;
        auto previous = genome_end->head;
        auto current = genome_end->head->primary_path.next.load(
            std::memory_order_acquire);
        uint64_t previous_start = 0;
        bool have_previous_segment = false;

        while (current && !current->isTail()) {
            if (!current->isSegment() ||
                !path_segments.insert(current.get()).second ||
                current->primary_path.prev.load(
                    std::memory_order_acquire) != previous ||
                previous->primary_path.next.load(
                    std::memory_order_acquire) != current ||
                !current->parent_block) {
                return false;
            }
            if (have_previous_segment && current->start < previous_start) {
                return false;
            }
            previous_start = current->start;
            have_previous_segment = true;

            const auto parent = current->parent_block;
            const auto anchor = parent->anchors.find(key);
            if (anchor == parent->anchors.end() ||
                anchor->second != current ||
                !pool_editor.contains(parent.get())) {
                return false;
            }
            if (expected_new_segments.count(current.get()) != 0 &&
                !seen_new_segments.insert(current.get()).second) {
                return false;
            }
            genome_end->setToSampling(current);
            previous = current;
            current = current->primary_path.next.load(
                std::memory_order_acquire);
        }

        if (!current || current != genome_end->tail ||
            genome_end->tail->primary_path.prev.load(
                std::memory_order_acquire) != previous ||
            previous->primary_path.next.load(
                std::memory_order_acquire) != genome_end->tail) {
            return false;
        }
    }
    return seen_new_segments == expected_new_segments;
}

bool auditPreparedReplacementsLocally(
    const std::vector<PreparedChain>& prepared,
    const BlockPoolEditor& pool_editor) {
    for (const auto& chain : prepared) {
        if (!chain.merged_block ||
            !pool_editor.contains(chain.merged_block.get()) ||
            chain.merged_block->anchors.size() != chain.paths.size()) {
            return false;
        }
        for (const auto& path : chain.paths) {
            const auto anchor = chain.merged_block->anchors.find(path.key);
            if (!path.merged_segment ||
                path.merged_segment->parent_block != chain.merged_block ||
                anchor == chain.merged_block->anchors.end() ||
                anchor->second != path.merged_segment) {
                return false;
            }
            const auto previous = path.merged_segment->primary_path.prev.load(
                std::memory_order_acquire);
            const auto next = path.merged_segment->primary_path.next.load(
                std::memory_order_acquire);
            if (!previous || !next ||
                previous->primary_path.next.load(
                    std::memory_order_acquire) != path.merged_segment ||
                next->primary_path.prev.load(
                    std::memory_order_acquire) != path.merged_segment) {
                return false;
            }
        }
    }
    return true;
}

bool rebuildSamplingAndAuditPaths(
    RaMeshMultiGenomeGraph& graph,
    const std::set<SpeciesChrPair>& affected_paths,
    const BlockPoolEditor& pool_editor) {
    std::unordered_set<const Segment*> globally_seen_segments;
    for (const auto& key : affected_paths) {
        auto& genome_end = genomeEndFor(graph, key);
        if (!genome_end.head || !genome_end.tail) return false;
        genome_end.sample_vec.clear();
        genome_end.sample_vec.resize(1, genome_end.head);

        auto previous = genome_end.head;
        auto current = genome_end.head->primary_path.next.load(
            std::memory_order_acquire);
        uint64_t previous_start = 0;
        bool have_previous = false;
        while (current && !current->isTail()) {
            if (!current->isSegment() ||
                !globally_seen_segments.insert(current.get()).second ||
                current->primary_path.prev.load(
                    std::memory_order_acquire) != previous ||
                previous->primary_path.next.load(
                    std::memory_order_acquire) != current ||
                !current->parent_block ||
                !pool_editor.contains(current->parent_block.get())) {
                return false;
            }
            if (have_previous && current->start < previous_start) {
                return false;
            }
            previous_start = current->start;
            have_previous = true;
            const auto anchor = current->parent_block->anchors.find(key);
            if (anchor == current->parent_block->anchors.end() ||
                anchor->second != current) {
                return false;
            }
            genome_end.setToSampling(current);
            previous = current;
            current = current->primary_path.next.load(
                std::memory_order_acquire);
        }
        if (!current || current != genome_end.tail ||
            genome_end.tail->primary_path.prev.load(
                std::memory_order_acquire) != previous ||
            previous->primary_path.next.load(
                std::memory_order_acquire) != genome_end.tail) {
            return false;
        }
    }
    return true;
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

class PreparedReplacementOwner {
public:
    explicit PreparedReplacementOwner(std::vector<PreparedChain>& prepared)
        : prepared_(&prepared) {}
    PreparedReplacementOwner(const PreparedReplacementOwner&) = delete;
    PreparedReplacementOwner& operator=(const PreparedReplacementOwner&) =
        delete;
    ~PreparedReplacementOwner() {
        if (prepared_) {
            detachPreparedBlocks(*prepared_);
        }
    }

    void release() noexcept { prepared_ = nullptr; }

private:
    std::vector<PreparedChain>* prepared_;
};

void commitPreparedExactBatch(
    RaMeshMultiGenomeGraph& graph,
    std::vector<PreparedChain>& prepared,
    BlockPoolEditor& pool_editor,
    std::set<SpeciesChrPair>& affected_paths) {
    if (prepared.empty()) return;
    PreparedReplacementOwner provisional_owner(prepared);

    for (const auto& chain : prepared) {
        for (const auto& path : chain.paths) {
            affected_paths.insert(path.key);
        }
    }

    std::vector<AppliedSplice> applied_splices;
    applied_splices.reserve(
        std::accumulate(
            prepared.begin(), prepared.end(), size_t{0},
            [](size_t total, const PreparedChain& chain) {
                return total + chain.paths.size();
            }));
    std::optional<BlockPoolEditor::Change> pool_change;
    try {
        for (auto& chain : prepared) {
            for (auto& path : chain.paths) {
                const bool reverse =
                    path.merged_segment->strand == Strand::REVERSE;
                const auto old_first = reverse
                    ? path.old_segments.back()
                    : path.old_segments.front();
                const auto old_last = reverse
                    ? path.old_segments.front()
                    : path.old_segments.back();
                const auto previous = old_first->primary_path.prev.load(
                    std::memory_order_acquire);
                const auto next = old_last->primary_path.next.load(
                    std::memory_order_acquire);
                if (!previous || !next ||
                    previous->primary_path.next.load(
                        std::memory_order_acquire) != old_first ||
                    next->primary_path.prev.load(
                        std::memory_order_acquire) != old_last) {
                    throw std::runtime_error(
                        "Exact Block merge path changed before batch commit");
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

        pool_change.emplace(pool_editor.replace(prepared));
        std::atomic_thread_fence(std::memory_order_seq_cst);
        if (!auditPreparedReplacementsLocally(prepared, pool_editor)) {
            throw std::runtime_error(
                "Exact Block merge batch local audit failed");
        }
    } catch (...) {
        if (pool_change.has_value()) {
            pool_editor.rollback(*pool_change);
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
        std::atomic_thread_fence(std::memory_order_seq_cst);
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
    provisional_owner.release();
}

bool buildDiagnosticBlockView(const BlockPtr& block,
                              const SpeciesName& reference_species,
                              const ChrName& reference_chromosome,
                              BlockView& view) {
    Realignment::BlockViewBuilder builder(reference_species);
    return builder.build(
        block, Realignment::BlockViewProfile::Diagnostics, view,
        reference_chromosome);
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
    spdlog::debug(
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
        spdlog::debug(
            "[exact-block-scan] stage={} reason={} count={} percent={:.6f}",
            stage,
            boundaryReasonName(static_cast<BoundaryReason>(index)),
            count, percent);
    }

    for (size_t bin = 0; bin < kGapBinCount; ++bin) {
        if (diagnostics.reference_gap_bins[bin] != 0) {
            spdlog::debug(
                "[exact-block-scan] stage={} reference_gap_bin={} count={}",
                stage, gapBinName(bin),
                diagnostics.reference_gap_bins[bin]);
        }
        if (diagnostics.query_gap_bins[bin] != 0) {
            spdlog::debug(
                "[exact-block-scan] stage={} query_max_gap_bin={} count={}",
                stage, gapBinName(bin),
                diagnostics.query_gap_bins[bin]);
        }
    }

    for (const auto& [transition, count] :
         diagnostics.participant_transitions) {
        spdlog::debug(
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
    EMPTY_INTERVAL_NOT_ADJACENT,
    EMPTY_REFERENCE_INTERVAL,
    ADJACENT_PAIR_ALL_ZERO,
    ADJACENT_PAIR_REQUIRES_MSA,
    ADJACENT_PAIR_SINGLE_GAP_EXCEEDED,
    INCOMPATIBLE_BOUNDARY_SIGNATURE,
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
        case MissingWindowReject::EMPTY_INTERVAL_NOT_ADJACENT:
            return "empty_interval_not_adjacent";
        case MissingWindowReject::EMPTY_REFERENCE_INTERVAL:
            return "empty_reference_interval";
        case MissingWindowReject::ADJACENT_PAIR_ALL_ZERO:
            return "adjacent_pair_all_zero";
        case MissingWindowReject::ADJACENT_PAIR_REQUIRES_MSA:
            return "adjacent_pair_requires_msa";
        case MissingWindowReject::ADJACENT_PAIR_SINGLE_GAP_EXCEEDED:
            return "adjacent_pair_single_gap_exceeded";
        case MissingWindowReject::INCOMPATIBLE_BOUNDARY_SIGNATURE:
            return "incompatible_boundary_signature";
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

enum class MissingWindowKind : uint8_t {
    SUBSET_CHAIN = 0,
    ADJACENT_PAIR,
    FULL_K_TRIPLE
};

struct MissingWindowCandidate {
    BlockView left;
    std::vector<BlockView> interiors;
    BlockView right;
    std::set<SpeciesName> missing_species;
    ChrName reference_chromosome;
    std::vector<MissingWindowPathContext> paths;
    std::string participant_pattern;
    MissingWindowKind kind = MissingWindowKind::SUBSET_CHAIN;
    size_t boundary_species_count = 0;
    bool zero_gap_deletion = false;
    bool hybrid_empty = false;
    bool adjacent_pair = false;
    bool direct_adjacent_pair = false;
    size_t empty_species_count = 0;
    size_t nonempty_species_count = 0;
    uint64_t gap_burden = 0;
};

struct MissingWindowPlan {
    BlockPtr left;
    std::vector<BlockPtr> interiors;
    BlockPtr right;
    ChrName reference_chromosome;
    MissingWindowKind kind = MissingWindowKind::SUBSET_CHAIN;
    size_t boundary_species_count = 0;
    bool hybrid_empty = false;
    bool adjacent_pair = false;
    bool reference_empty = false;
};

MissingWindowPlan makeMissingWindowPlan(
    const MissingWindowCandidate& candidate,
    const SpeciesName& reference_species) {
    MissingWindowPlan plan;
    plan.left = candidate.left.block;
    plan.interiors.reserve(candidate.interiors.size());
    for (const auto& interior : candidate.interiors) {
        plan.interiors.push_back(interior.block);
    }
    plan.right = candidate.right.block;
    plan.reference_chromosome = candidate.reference_chromosome;
    plan.kind = candidate.kind;
    plan.boundary_species_count = candidate.boundary_species_count;
    plan.hybrid_empty = candidate.hybrid_empty;
    plan.adjacent_pair = candidate.adjacent_pair;
    const auto reference_path = std::find_if(
        candidate.paths.begin(), candidate.paths.end(),
        [&](const MissingWindowPathContext& path) {
            return path.key.first == reference_species;
        });
    plan.reference_empty =
        reference_path != candidate.paths.end() &&
        reference_path->interval_length == 0;
    return plan;
}

struct MissingWindowPreparationTiming {
    double total_seconds = 0.0;
    double sequence_fetch_seconds = 0.0;
    double msa_seconds = 0.0;
};

class ScopedSeconds {
public:
    explicit ScopedSeconds(double* destination)
        : destination_(destination),
          start_(std::chrono::steady_clock::now()) {}
    ScopedSeconds(const ScopedSeconds&) = delete;
    ScopedSeconds& operator=(const ScopedSeconds&) = delete;
    ~ScopedSeconds() {
        if (destination_) {
            *destination_ += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start_).count();
        }
    }

private:
    double* destination_;
    std::chrono::steady_clock::time_point start_;
};

bool buildRealignBlockView(const BlockPtr& block,
                           const SpeciesName& reference_species,
                           BlockView& view) {
    Realignment::BlockViewBuilder builder(reference_species);
    return builder.build(
        block, Realignment::BlockViewProfile::MissingWindow, view);
}

bool getCachedRealignBlockView(
    const BlockPtr& block,
    const SpeciesName& reference_species,
    BlockViewCache* cache,
    BlockView& view) {
    if (!cache) {
        return buildRealignBlockView(block, reference_species, view);
    }
    return cache->build(
        block, Realignment::BlockViewProfile::MissingWindow, view);
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
    uint_t adjacent_pair_gap_max,
    bool zero_gap_phase,
    MissingWindowCandidate& candidate,
    bool full_k_triple = false,
    BlockViewCache* view_cache = nullptr) {
    candidate = MissingWindowCandidate{};
    if (!getCachedRealignBlockView(
            left_block, reference_species, view_cache, candidate.left) ||
        !getCachedRealignBlockView(
            right_block, reference_species, view_cache, candidate.right) ||
        candidate.left.species_count < 2 ||
        candidate.right.species_count < 2) {
        return MissingWindowReject::INVALID_BLOCK;
    }
    candidate.adjacent_pair = interior_blocks.empty();
    candidate.kind = full_k_triple
                         ? MissingWindowKind::FULL_K_TRIPLE
                         : (candidate.adjacent_pair
                                ? MissingWindowKind::ADJACENT_PAIR
                                : MissingWindowKind::SUBSET_CHAIN);
    if (full_k_triple && interior_blocks.size() != 1) {
        return MissingWindowReject::PARTICIPANT_PATTERN;
    }
    candidate.interiors.reserve(interior_blocks.size());
    for (const auto& block : interior_blocks) {
        BlockView view;
        if (!getCachedRealignBlockView(
                block, reference_species, view_cache, view)) {
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
    if (zero_gap_phase && !candidate.adjacent_pair &&
        candidate.interiors.size() != 1) {
        return MissingWindowReject::PARTICIPANT_PATTERN;
    }
    candidate.participant_pattern =
        std::to_string(expected_species_count);
    candidate.boundary_species_count = expected_species_count;
    if (candidate.adjacent_pair) {
        candidate.participant_pattern += "-0";
    }
    for (const auto& interior : candidate.interiors) {
        const bool participant_pattern_valid =
            full_k_triple
                ? interior.species_count == expected_species_count &&
                      sameAnchorKeys(interior, candidate.left)
                : interior.species_count < expected_species_count &&
                      participantSpeciesAreSubset(
                          interior, candidate.left);
        if (interior.declared_reference_chromosome != reference_chromosome ||
            !participant_pattern_valid) {
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
    if (candidate.adjacent_pair) {
        for (const auto& path : candidate.paths) {
            if (path.interval_length == 0) {
                candidate.missing_species.insert(path.key.first);
            }
        }
    }
    if (!candidate.adjacent_pair && candidate.missing_species.empty() &&
        !full_k_triple) {
        return MissingWindowReject::PARTICIPANT_PATTERN;
    }

    if (full_k_triple) {
        for (const auto& path : candidate.paths) {
            if (path.interior_segments.size() != 1) {
                return MissingWindowReject::PARTICIPANT_PATTERN;
            }
            const int64_t left_gap = segmentGapInReferenceOrder(
                path.left_anchor, path.interior_segments.front());
            const int64_t right_gap = segmentGapInReferenceOrder(
                path.interior_segments.front(), path.right_anchor);
            if (left_gap < 0 || right_gap < 0) {
                return MissingWindowReject::ANCHOR_ORDER_MISMATCH;
            }
            candidate.gap_burden +=
                static_cast<uint64_t>(left_gap) +
                static_cast<uint64_t>(right_gap);
        }
        if (candidate.gap_burden == 0) {
            return MissingWindowReject::PARTICIPANT_PATTERN;
        }
    }

    if (zero_gap_phase && !candidate.adjacent_pair &&
        (candidate.interiors.front().species_count == 0 ||
         candidate.interiors.front().species_count >=
             expected_species_count ||
         candidate.missing_species.size() !=
             expected_species_count -
                 candidate.interiors.front().species_count)) {
        return MissingWindowReject::PARTICIPANT_PATTERN;
    }

    if (!zero_gap_phase) {
        bool any_zero = false;
        bool any_nonzero = false;
        bool reference_nonzero = false;
        for (const auto& path : candidate.paths) {
            any_zero = any_zero || path.interval_length == 0;
            any_nonzero = any_nonzero || path.interval_length != 0;
            candidate.nonempty_species_count +=
                path.interval_length != 0 ? 1 : 0;
            if (path.key.first == reference_species) {
                reference_nonzero = path.interval_length != 0;
            }
            if (path.interval_length > maximum_span) {
                return MissingWindowReject::INTERVAL_INVALID;
            }
        }
        if (!any_nonzero) {
            return candidate.adjacent_pair
                       ? MissingWindowReject::ADJACENT_PAIR_ALL_ZERO
                       : MissingWindowReject::INTERVAL_INVALID;
        }
        if (candidate.adjacent_pair &&
            candidate.nonempty_species_count == 1) {
            const auto nonempty_it = std::find_if(
                candidate.paths.begin(), candidate.paths.end(),
                [](const MissingWindowPathContext& path) {
                    return path.interval_length != 0;
                });
            if (nonempty_it == candidate.paths.end() ||
                nonempty_it->interval_length > adjacent_pair_gap_max) {
                return MissingWindowReject::
                    ADJACENT_PAIR_SINGLE_GAP_EXCEEDED;
            }
            return MissingWindowReject::ADJACENT_PAIR_REQUIRES_MSA;
        }
        if (any_zero) {
            if (!reference_nonzero && !candidate.adjacent_pair) {
                return MissingWindowReject::EMPTY_REFERENCE_INTERVAL;
            }
            candidate.hybrid_empty = true;
            for (const auto& path : candidate.paths) {
                if (path.interval_length != 0) {
                    continue;
                }
                ++candidate.empty_species_count;
                if (!path.interior_segments.empty() ||
                    !segmentsAreDirectlyAdjacent(
                        path.left_anchor, path.right_anchor)) {
                    return MissingWindowReject::EMPTY_INTERVAL_NOT_ADJACENT;
                }
            }
        }
        return MissingWindowReject::COUNT;
    }

    if (candidate.adjacent_pair) {
        for (const auto& path : candidate.paths) {
            candidate.nonempty_species_count +=
                path.interval_length != 0 ? 1 : 0;
        }
        if (candidate.nonempty_species_count == 0) {
            return MissingWindowReject::ADJACENT_PAIR_ALL_ZERO;
        }
        if (candidate.nonempty_species_count != 1) {
            return MissingWindowReject::ADJACENT_PAIR_REQUIRES_MSA;
        }
        const auto nonempty_it = std::find_if(
            candidate.paths.begin(), candidate.paths.end(),
            [](const MissingWindowPathContext& path) {
                return path.interval_length != 0;
            });
        if (nonempty_it == candidate.paths.end() ||
            nonempty_it->interval_length > adjacent_pair_gap_max) {
            return MissingWindowReject::ADJACENT_PAIR_SINGLE_GAP_EXCEEDED;
        }
        candidate.zero_gap_deletion = true;
        candidate.direct_adjacent_pair = true;
        candidate.empty_species_count =
            candidate.paths.size() - 1;
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

MissingWindowReject materializeMissingWindowPlan(
    const MissingWindowPlan& plan,
    const SpeciesName& reference_species,
    uint_t maximum_span,
    uint_t zero_gap_maximum_span,
    uint_t adjacent_pair_gap_max,
    bool zero_gap_phase,
    MissingWindowCandidate& candidate) {
    return buildMissingWindowCandidate(
        plan.left, plan.interiors, plan.right, reference_species,
        plan.reference_chromosome, plan.boundary_species_count,
        maximum_span, zero_gap_maximum_span, adjacent_pair_gap_max,
        zero_gap_phase, candidate,
        plan.kind == MissingWindowKind::FULL_K_TRIPLE, nullptr);
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

std::optional<PreparedChain> prepareAdjacentPairDirect(
    const MissingWindowCandidate& candidate,
    const SpeciesName& reference_species,
    MissingWindowReject& rejection) {
    rejection = MissingWindowReject::ZERO_GAP_PREPARATION_INVALID;
    if (!candidate.adjacent_pair || !candidate.interiors.empty() ||
        candidate.paths.empty() || candidate.nonempty_species_count == 0 ||
        !candidate.left.reference_segment ||
        !candidate.right.reference_segment ||
        candidate.left.reference_segment->strand != Strand::FORWARD ||
        candidate.right.reference_segment->strand != Strand::FORWARD) {
        return std::nullopt;
    }

    const auto reference_path_it = std::find_if(
        candidate.paths.begin(), candidate.paths.end(),
        [&](const MissingWindowPathContext& path) {
            return path.key.first == reference_species;
        });
    if (reference_path_it == candidate.paths.end()) {
        return std::nullopt;
    }
    const uint64_t reference_interval_length =
        reference_path_it->interval_length;
    const bool direct_single_gap =
        candidate.direct_adjacent_pair &&
        candidate.nonempty_species_count == 1;
    const bool reference_empty = reference_interval_length == 0;
    if (!direct_single_gap && !reference_empty) {
        return std::nullopt;
    }
    const uint64_t full_reference_length =
        static_cast<uint64_t>(candidate.left.reference_segment->length) +
        reference_interval_length +
        candidate.right.reference_segment->length;
    if (full_reference_length == 0 ||
        full_reference_length > std::numeric_limits<uint_t>::max() ||
        segmentEnd(candidate.right.reference_segment) -
                candidate.left.reference_segment->start !=
            full_reference_length) {
        return std::nullopt;
    }

    PreparedChain prepared;
    prepared.candidate.blocks = {
        candidate.left.block, candidate.right.block};
    prepared.candidate.reference_chromosome =
        candidate.reference_chromosome;
    prepared.candidate.reference_start =
        candidate.left.reference_segment->start;
    prepared.candidate.reference_length =
        static_cast<uint_t>(full_reference_length);
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
        path.old_segments = {
            context.left_anchor, context.right_anchor};

        Cigar_t merged_cigar;
        Cigar_t normalized;
        if (!normalizeAndValidateCigar(
                context.left_anchor,
                candidate.left.reference_segment->length,
                normalized) ||
            !appendCigarChecked(merged_cigar, normalized)) {
            return reject_prepared();
        }

        if (reference_interval_length != 0) {
            appendCigarOp(
                merged_cigar,
                context.key.first == reference_species ? 'M' : 'D',
                static_cast<uint32_t>(reference_interval_length));
        } else if (context.interval_length != 0) {
            appendCigarOp(
                merged_cigar, 'I', context.interval_length);
        }

        if (!normalizeAndValidateCigar(
                context.right_anchor,
                candidate.right.reference_segment->length,
                normalized) ||
            !appendCigarChecked(merged_cigar, normalized)) {
            return reject_prepared();
        }

        const uint64_t merged_start = std::min<uint64_t>(
            context.left_anchor->start,
            context.right_anchor->start);
        const uint64_t merged_end = std::max<uint64_t>(
            segmentEnd(context.left_anchor),
            segmentEnd(context.right_anchor));
        const uint64_t merged_length =
            static_cast<uint64_t>(context.left_anchor->length) +
            context.interval_length + context.right_anchor->length;
        if (merged_length == 0 ||
            merged_length > std::numeric_limits<uint_t>::max() ||
            merged_start > std::numeric_limits<uint_t>::max() ||
            merged_end - merged_start != merged_length ||
            countRefLength(merged_cigar) != full_reference_length ||
            countQryLength(merged_cigar) != merged_length) {
            return reject_prepared();
        }

        path.merged_segment = Segment::create(
            static_cast<uint_t>(merged_start),
            static_cast<uint_t>(merged_length), context.strand,
            std::move(merged_cigar), AlignRole::PRIMARY,
            SegmentRole::SEGMENT, prepared.merged_block);
        path.merged_segment->left_extend =
            context.left_anchor->left_extend;
        path.merged_segment->right_extend =
            context.right_anchor->right_extend;
        if (!prepared.merged_block->anchors.emplace(
                context.key, path.merged_segment).second) {
            return reject_prepared();
        }
        prepared.paths.push_back(std::move(path));
    }

    const SpeciesChrPair reference_key{
        reference_species, candidate.reference_chromosome};
    const auto prepared_reference =
        prepared.merged_block->anchors.find(reference_key);
    if (prepared_reference == prepared.merged_block->anchors.end() ||
        prepared_reference->second->strand != Strand::FORWARD ||
        prepared_reference->second->start !=
            prepared.candidate.reference_start ||
        prepared_reference->second->length !=
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
    MissingWindowReject& rejection,
    bool& msa_invoked,
    bool& reference_only_bypass,
    bool& reference_empty_bypass,
    MissingWindowPreparationTiming* timing = nullptr) {
    MissingWindowPreparationTiming ignored_timing;
    auto& measured = timing ? *timing : ignored_timing;
    ScopedSeconds total_timer(&measured.total_seconds);
    msa_invoked = false;
    reference_only_bypass = false;
    reference_empty_bypass = false;
    if (candidate.adjacent_pair) {
        const auto reference_path_it = std::find_if(
            candidate.paths.begin(), candidate.paths.end(),
            [&](const MissingWindowPathContext& path) {
                return path.key.first == reference_species;
            });
        if (reference_path_it != candidate.paths.end() &&
            reference_path_it->interval_length == 0) {
            auto prepared = prepareAdjacentPairDirect(
                candidate, reference_species, rejection);
            reference_empty_bypass = prepared.has_value();
            return prepared;
        }
    }
    std::unordered_map<ChrName, std::string> raw_sequences;
    raw_sequences.reserve(candidate.paths.size());
    {
        ScopedSeconds fetch_timer(&measured.sequence_fetch_seconds);
        for (const auto& path : candidate.paths) {
            if (candidate.hybrid_empty && path.interval_length == 0) {
                continue;
            }
            std::string sequence;
            if (!fetchWindowSequence(managers, path, sequence) ||
                !raw_sequences.emplace(
                    path.key.first, std::move(sequence)).second) {
                rejection = MissingWindowReject::SEQUENCE_UNAVAILABLE;
                return std::nullopt;
            }
        }
    }

    const auto raw_reference_it = raw_sequences.find(reference_species);
    const size_t raw_reference_length =
        raw_reference_it == raw_sequences.end()
            ? 0
            : raw_reference_it->second.size();
    auto aligned_sequences = std::move(raw_sequences);
    bool msa_succeeded = false;
    if (aligned_sequences.size() == 1 &&
        aligned_sequences.count(reference_species) == 1) {
        reference_only_bypass = true;
        msa_succeeded = true;
    } else {
        msa_invoked = true;
        try {
            ScopedSeconds msa_timer(&measured.msa_seconds);
            msa_succeeded = alignSequencesWithExternalMsa(
                msa_executable, aligned_sequences);
        } catch (const std::exception& error) {
            spdlog::warn(
                "[species-mismatch-realign] minipoa invocation failed: {}",
                error.what());
        }
    }
    if (!msa_succeeded) {
        rejection = MissingWindowReject::MINIPOA_FAILED;
        return std::nullopt;
    }
    const auto reference_it =
        aligned_sequences.find(reference_species);
    std::string synthetic_reference_alignment;
    const std::string* aligned_reference = nullptr;
    if (reference_it != aligned_sequences.end()) {
        aligned_reference = &reference_it->second;
    } else if (candidate.adjacent_pair) {
        if (aligned_sequences.empty()) {
            rejection = MissingWindowReject::MSA_INVALID;
            return std::nullopt;
        }
        const size_t aligned_length =
            aligned_sequences.begin()->second.size();
        if (aligned_length == 0 ||
            std::any_of(
                aligned_sequences.begin(), aligned_sequences.end(),
                [&](const auto& item) {
                    return item.second.size() != aligned_length;
                })) {
            rejection = MissingWindowReject::MSA_INVALID;
            return std::nullopt;
        }
        synthetic_reference_alignment.assign(
            aligned_length, '-');
        aligned_reference = &synthetic_reference_alignment;
    } else {
        rejection = MissingWindowReject::MSA_INVALID;
        return std::nullopt;
    }

    if (candidate.hybrid_empty || candidate.adjacent_pair) {
        if (!candidate.left.reference_segment ||
            !candidate.right.reference_segment ||
            candidate.left.reference_segment->strand != Strand::FORWARD ||
            candidate.right.reference_segment->strand != Strand::FORWARD) {
            rejection = MissingWindowReject::MSA_INVALID;
            return std::nullopt;
        }
        const auto reference_path_it = std::find_if(
            candidate.paths.begin(), candidate.paths.end(),
            [&](const MissingWindowPathContext& path) {
                return path.key.first == reference_species;
            });
        if (reference_path_it == candidate.paths.end()) {
            rejection = MissingWindowReject::MSA_INVALID;
            return std::nullopt;
        }
        const uint64_t reference_interval_length =
            reference_path_it->interval_length;
        const uint64_t full_reference_length =
            static_cast<uint64_t>(
                candidate.left.reference_segment->length) +
            reference_interval_length +
            candidate.right.reference_segment->length;
        if ((!candidate.adjacent_pair &&
             reference_interval_length == 0) ||
            full_reference_length >
                std::numeric_limits<uint_t>::max() ||
            segmentEnd(candidate.right.reference_segment) -
                    candidate.left.reference_segment->start !=
                full_reference_length) {
            rejection = MissingWindowReject::MSA_INVALID;
            return std::nullopt;
        }

        PreparedChain prepared;
        prepared.candidate.blocks.push_back(candidate.left.block);
        for (const auto& interior : candidate.interiors) {
            prepared.candidate.blocks.push_back(interior.block);
        }
        prepared.candidate.blocks.push_back(candidate.right.block);
        prepared.candidate.reference_chromosome =
            candidate.reference_chromosome;
        prepared.candidate.reference_start =
            candidate.left.reference_segment->start;
        prepared.candidate.reference_length =
            static_cast<uint_t>(full_reference_length);
        prepared.candidate.species_count = candidate.paths.size();
        prepared.merged_block = Block::createEmpty(
            candidate.reference_chromosome, candidate.paths.size());
        prepared.paths.reserve(candidate.paths.size());
        const auto reject_prepared =
            [&]() -> std::optional<PreparedChain> {
            detachPreparedChain(prepared);
            return std::nullopt;
        };

        for (const auto& context : candidate.paths) {
            Cigar_t middle_cigar;
            if (context.interval_length == 0) {
                if (reference_interval_length == 0) {
                    middle_cigar.clear();
                } else {
                appendCigarOp(
                    middle_cigar, 'D',
                    static_cast<uint32_t>(reference_interval_length));
                }
            } else if (context.key.first == reference_species) {
                appendCigarOp(
                    middle_cigar, 'M', context.interval_length);
            } else {
                const auto aligned_it =
                    aligned_sequences.find(context.key.first);
                if (aligned_it == aligned_sequences.end() ||
                    !cigarFromMsaRows(
                        *aligned_reference, aligned_it->second,
                        middle_cigar)) {
                    rejection = MissingWindowReject::MSA_INVALID;
                    return reject_prepared();
                }
            }
            if (countRefLength(middle_cigar) !=
                    reference_interval_length ||
                countQryLength(middle_cigar) !=
                    context.interval_length) {
                rejection = MissingWindowReject::MSA_INVALID;
                return reject_prepared();
            }

            PathReplacement path;
            path.key = context.key;
            path.old_segments.push_back(context.left_anchor);
            path.old_segments.insert(
                path.old_segments.end(),
                context.interior_segments.begin(),
                context.interior_segments.end());
            path.old_segments.push_back(context.right_anchor);

            Cigar_t merged_cigar;
            Cigar_t normalized;
            if (!normalizeAndValidateCigar(
                    context.left_anchor,
                    candidate.left.reference_segment->length,
                    normalized) ||
                !appendCigarChecked(merged_cigar, normalized) ||
                !appendCigarChecked(merged_cigar, middle_cigar) ||
                !normalizeAndValidateCigar(
                    context.right_anchor,
                    candidate.right.reference_segment->length,
                    normalized) ||
                !appendCigarChecked(merged_cigar, normalized)) {
                rejection = MissingWindowReject::MSA_INVALID;
                return reject_prepared();
            }

            const uint64_t merged_length =
                static_cast<uint64_t>(context.left_anchor->length) +
                context.interval_length +
                context.right_anchor->length;
            const uint64_t merged_start = std::min<uint64_t>(
                context.left_anchor->start,
                context.right_anchor->start);
            const uint64_t merged_end = std::max<uint64_t>(
                segmentEnd(context.left_anchor),
                segmentEnd(context.right_anchor));
            if (merged_length == 0 ||
                merged_length > std::numeric_limits<uint_t>::max() ||
                merged_start > std::numeric_limits<uint_t>::max() ||
                merged_end - merged_start != merged_length ||
                countRefLength(merged_cigar) !=
                    full_reference_length ||
                countQryLength(merged_cigar) != merged_length) {
                rejection = MissingWindowReject::MSA_INVALID;
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
                rejection = MissingWindowReject::MSA_INVALID;
                return reject_prepared();
            }
            prepared.paths.push_back(std::move(path));
        }

        const SpeciesChrPair reference_key{
            reference_species, candidate.reference_chromosome};
        const auto prepared_reference =
            prepared.merged_block->anchors.find(reference_key);
        if (prepared_reference == prepared.merged_block->anchors.end() ||
            prepared_reference->second->strand != Strand::FORWARD ||
            prepared_reference->second->start !=
                prepared.candidate.reference_start ||
            prepared_reference->second->length !=
                prepared.candidate.reference_length) {
            rejection = MissingWindowReject::MSA_INVALID;
            return reject_prepared();
        }

        rejection = MissingWindowReject::COUNT;
        return prepared;
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
    const auto reject_prepared = [&]() -> std::optional<PreparedChain> {
        detachPreparedChain(prepared);
        return std::nullopt;
    };

    for (const auto& context : candidate.paths) {
        const auto aligned_it =
            aligned_sequences.find(context.key.first);
        if (aligned_it == aligned_sequences.end()) {
            rejection = MissingWindowReject::MSA_INVALID;
            return reject_prepared();
        }

        Cigar_t cigar;
        if (context.key.first == reference_species) {
            cigar.push_back(
                cigarToInt('M', context.interval_length));
        } else if (!cigarFromMsaRows(
                       reference_it->second, aligned_it->second,
                       cigar)) {
            rejection = MissingWindowReject::MSA_INVALID;
            return reject_prepared();
        }
        if (countRefLength(cigar) !=
                raw_reference_length ||
            countQryLength(cigar) != context.interval_length) {
            rejection = MissingWindowReject::MSA_INVALID;
            return reject_prepared();
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
        if (!prepared.merged_block->anchors.emplace(
                context.key, path.merged_segment).second) {
            rejection = MissingWindowReject::MSA_INVALID;
            return reject_prepared();
        }
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
    BlockPoolEditor& pool_editor,
    std::set<SpeciesChrPair>& affected_paths,
    size_t species_count,
    const char* phase,
    bool require_block_reduction,
    bool require_participant_deficit_reduction = false) {
    if (prepared.empty()) {
        return {};
    }
    if (prepared.size() != accepted_candidates.size()) {
        throw std::runtime_error(
            "Species-mismatch prepared candidate vectors changed size");
    }
    PreparedReplacementOwner provisional_owner(prepared);

    size_t replaced_old_blocks = 0;
    std::unordered_set<const Block*> batch_old_blocks;
    batch_old_blocks.reserve(prepared.size() * 3);
    bool requires_participant_deficit_progress = false;
    uint64_t old_participant_deficit = 0;
    uint64_t new_participant_deficit = 0;
    for (size_t index = 0; index < prepared.size(); ++index) {
        const auto& boundary = accepted_candidates[index].left;
        const auto& merged_block = prepared[index].merged_block;
        if (!merged_block) {
            throw std::runtime_error(
                "Species-mismatch prepared candidate has no merged Block");
        }
        {
            std::shared_lock merged_lock(merged_block->rw);
            if (merged_block->anchors.size() !=
                boundary.anchors.size()) {
                throw std::runtime_error(
                    "Species-mismatch merged Block changed the boundary "
                    "participant count");
            }
            for (const auto& [key, unused] : boundary.anchors) {
                (void)unused;
                if (merged_block->anchors.count(key) == 0) {
                    throw std::runtime_error(
                        "Species-mismatch merged Block changed the "
                        "boundary participant keys");
                }
            }
        }
        const auto& blocks_in_candidate =
            prepared[index].candidate.blocks;
        if (blocks_in_candidate.empty()) {
            throw std::runtime_error(
                "Species-mismatch prepared candidate has no old Blocks");
        }
        replaced_old_blocks += blocks_in_candidate.size();
        if (merged_block->anchors.size() > species_count) {
            throw std::runtime_error(
                "Species-mismatch merged Block exceeds species count");
        }
        new_participant_deficit +=
            species_count - merged_block->anchors.size();
        if (require_participant_deficit_reduction &&
            blocks_in_candidate.size() == 1) {
            if (accepted_candidates[index].kind ==
                MissingWindowKind::FULL_K_TRIPLE) {
                if (accepted_candidates[index].gap_burden == 0 ||
                    blocks_in_candidate.front()->anchors.size() !=
                        merged_block->anchors.size()) {
                    throw std::runtime_error(
                        "Full-K triple does not reduce its boundary gap "
                        "burden");
                }
            } else {
                requires_participant_deficit_progress = true;
                if (blocks_in_candidate.front()->anchors.size() >=
                    merged_block->anchors.size()) {
                    throw std::runtime_error(
                        "Species-mismatch candidate does not reduce either "
                        "the Block count or its participant deficit");
                }
            }
        }
        for (const auto& old_block : blocks_in_candidate) {
            if (!old_block ||
                !batch_old_blocks.insert(old_block.get()).second ||
                old_block->anchors.size() > species_count) {
                throw std::runtime_error(
                    "Species-mismatch prepared candidate has invalid or "
                    "overlapping old Blocks");
            }
            old_participant_deficit +=
                species_count - old_block->anchors.size();
        }
    }
    if (replaced_old_blocks < prepared.size()) {
        throw std::runtime_error(
            "Species-mismatch replacement Block accounting underflow");
    }

    const size_t blocks_before = pool_editor.activeCount();
    const size_t expected_reduction =
        replaced_old_blocks - prepared.size();
    if (require_block_reduction && expected_reduction == 0) {
        throw std::runtime_error(
            "Zero-gap replacement did not reduce the Block pool");
    }

    for (const auto& chain : prepared) {
        for (const auto& path : chain.paths) {
            affected_paths.insert(path.key);
        }
    }

    std::vector<AppliedSplice> applied_splices;
    applied_splices.reserve(prepared.size() * species_count);
    std::optional<BlockPoolEditor::Change> pool_change;
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

        pool_change.emplace(pool_editor.replace(prepared));
        if (pool_editor.activeCount() + expected_reduction != blocks_before) {
            throw std::runtime_error(
                "Species-mismatch replacement changed the Block pool by "
                "an unexpected amount");
        }
        if (require_participant_deficit_reduction &&
            requires_participant_deficit_progress &&
            expected_reduction == 0 &&
            new_participant_deficit >= old_participant_deficit) {
            throw std::runtime_error(
                "Species-mismatch fixed-point replacement did not reduce "
                "the global participant deficit");
        }
        std::atomic_thread_fence(std::memory_order_seq_cst);
        if (!auditPreparedReplacementsLocally(prepared, pool_editor)) {
            throw std::runtime_error(
                "Species-mismatch post-commit local audit failed");
        }
    } catch (...) {
        if (pool_change.has_value()) {
            pool_editor.rollback(*pool_change);
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

    provisional_owner.release();

    return MissingWindowCommitResult{
        prepared.size(), replaced_old_blocks,
        blocks_before, pool_editor.activeCount()};
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
    uint_t maximum_query_gap) try {
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
    if (spdlog::should_log(spdlog::level::debug)) {
        const auto diagnostics = collectBoundaryDiagnostics(
            *this, reference_species, maximum_reference_span);
        logBoundaryDiagnostics(
            diagnostics, reference_species, "pre-merge-current-reference",
            maximum_reference_span);
    }
    AllocationPhaseScope scan_phase("scan");

    std::vector<CandidateChain> candidates;
    candidates.reserve(std::min<size_t>(blocks.size() / 8 + 1, 65536));
    for (auto& [chromosome, genome_end] :
         reference_graph_it->second.chr2end) {
        auto current =
            genome_end.head->primary_path.next.load(std::memory_order_acquire);
        while (current && !current->isTail()) {
            const auto current_block = current->parent_block;
            BlockView first_view;
            if (!buildBlockView(
                    current_block, reference_species, first_view) ||
                first_view.reference_segment != current ||
                first_view.block->ref_chr != chromosome) {
                current =
                    current->primary_path.next.load(std::memory_order_acquire);
                continue;
            }

            CandidateChain candidate;
            candidate.blocks.push_back(first_view.block);
            candidate.reference_chromosome = chromosome;
            candidate.reference_start = current->start;
            candidate.reference_length = current->length;
            candidate.species_count = first_view.species_count;
            BlockView last_view = std::move(first_view);

            auto last_reference = current;
            auto next_reference =
                current->primary_path.next.load(std::memory_order_acquire);
            while (next_reference && !next_reference->isTail()) {
                BlockView next_view;
                if (!buildBlockView(
                        next_reference->parent_block,
                        reference_species, next_view) ||
                    next_view.reference_segment != next_reference ||
                    next_view.species_count != candidate.species_count ||
                    !canMergePair(
                        last_view, next_view,
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

                candidate.blocks.push_back(next_view.block);
                last_view = std::move(next_view);
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
        spdlog::debug(
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

    size_t eliminated_boundaries = 0;
    size_t merged_old_blocks = 0;
    size_t longest_chain = 0;
    std::map<size_t, std::pair<size_t, size_t>> by_species_count;
    const bool detailed_stats = spdlog::should_log(spdlog::level::debug);
    std::vector<const Block*> candidate_blocks;
    candidate_blocks.reserve(
        std::min<size_t>(blocks.size(), candidates.size() * 2));
    for (const auto& candidate : candidates) {
        eliminated_boundaries += candidate.blocks.size() - 1;
        merged_old_blocks += candidate.blocks.size();
        longest_chain = std::max(longest_chain, candidate.blocks.size());
        if (detailed_stats) {
            auto& stats = by_species_count[candidate.species_count];
            ++stats.first;
            stats.second += candidate.blocks.size() - 1;
        }
        for (const auto& block : candidate.blocks) {
            candidate_blocks.push_back(block.get());
        }
    }
    std::sort(
        candidate_blocks.begin(), candidate_blocks.end(),
        std::less<const Block*>{});
    if (std::adjacent_find(
            candidate_blocks.begin(), candidate_blocks.end()) !=
        candidate_blocks.end()) {
        throw std::runtime_error(
            "Exact Block merge produced overlapping candidate chains");
    }
    candidate_blocks.clear();
    candidate_blocks.shrink_to_fit();

    AllocationPhaseScope pool_index_phase("pool-index");
    BlockPoolEditor pool_editor(blocks);
    std::set<SpeciesChrPair> affected_paths;
    const size_t active_blocks_before = pool_editor.activeCount();
    size_t batch_begin = 0;
    size_t committed_chains = 0;
    bool warned_batch_reduction = false;
    g_allocation_candidates = candidates.size();
    try {
        while (batch_begin < candidates.size()) {
            size_t target_end = batch_begin;
            size_t anchor_records = 0;
            while (target_end < candidates.size() &&
                   target_end - batch_begin < kExactBatchMaximumChains) {
                const size_t next_anchors =
                    candidates[target_end].species_count;
                if (target_end != batch_begin &&
                    anchor_records + next_anchors >
                        kExactBatchMaximumAnchors) {
                    break;
                }
                anchor_records += next_anchors;
                ++target_end;
            }
            size_t attempt_end = target_end;
            while (true) {
                std::vector<PreparedChain> prepared;
                g_allocation_batch_size = attempt_end - batch_begin;
                g_allocation_prepared = 0;
                try {
                    AllocationPhaseScope prepare_phase("preparation");
                    prepared.reserve(attempt_end - batch_begin);
                    for (size_t index = batch_begin;
                         index < attempt_end; ++index) {
                        prepared.push_back(prepareChain(
                            candidates[index], reference_species));
                        g_allocation_prepared = prepared.size();
                    }
                    AllocationPhaseScope commit_phase("commit-audit");
                    commitPreparedExactBatch(
                        *this, prepared, pool_editor, affected_paths);
                } catch (const std::bad_alloc&) {
                    detachPreparedBlocks(prepared);
                    const size_t attempt_size =
                        attempt_end - batch_begin;
                    if (attempt_size <= 1) throw;
                    attempt_end = batch_begin + attempt_size / 2;
                    if (!warned_batch_reduction) {
                        spdlog::warn(
                            "[exact-block-merge] reference={} allocation "
                            "pressure detected; reducing batch size "
                            "(subsequent reductions are debug-only)",
                            reference_species);
                        warned_batch_reduction = true;
                    }
                    spdlog::debug(
                        "[exact-block-merge] reference={} reducing "
                        "allocation batch from {} to {} chains",
                        reference_species, attempt_size,
                        attempt_end - batch_begin);
                    continue;
                }
                for (size_t index = batch_begin;
                     index < attempt_end; ++index) {
                    candidates[index].blocks.clear();
                }
                committed_chains += attempt_end - batch_begin;
                batch_begin = attempt_end;
                break;
            }
        }
        AllocationPhaseScope audit_phase("final-pool-audit");
        if (!rebuildSamplingAndAuditPaths(
                *this, affected_paths, pool_editor)) {
            throw std::runtime_error(
                "Exact Block merge final path audit failed");
        }
        pool_editor.compact();
        if (!pool_editor.auditPool()) {
            throw std::runtime_error(
                "Exact Block merge final Block pool audit failed");
        }
    } catch (...) {
        spdlog::error(
            "[exact-block-merge] reference={} failed after {} committed "
            "chains; the current batch was restored",
            reference_species, committed_chains);
        throw;
    }
    const size_t active_blocks_after = pool_editor.activeCount();
    spdlog::info(
        "[exact-block-merge] reference={} candidates={} old_blocks={} "
        "eliminated_boundaries={} blocks_before={} blocks_after={} "
        "longest_chain={} max_reference_span={} max_query_gap={}",
        reference_species, candidates.size(), merged_old_blocks,
        eliminated_boundaries, active_blocks_before, active_blocks_after,
        longest_chain, maximum_reference_span, maximum_query_gap);
    for (const auto& [species_count, stats] : by_species_count) {
        spdlog::debug(
            "[exact-block-merge] participants={} chains={} eliminated_boundaries={}",
            species_count, stats.first, stats.second);
    }

    return eliminated_boundaries;
} catch (const std::bad_alloc&) {
    logAllocationFailure(
        "exact-block-merge", reference_species, blocks.size(),
        g_allocation_candidates, g_allocation_prepared,
        g_allocation_batch_size);
    throw;
}

size_t RaMeshMultiGenomeGraph::realignSingleMissingSpeciesWindows(
    const SpeciesName& reference_species,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>&
        seqpro_managers,
    const std::string& msa_executable,
    uint_t maximum_span,
    uint_t parallel_threads,
    uint_t zero_gap_maximum_span,
    uint_t adjacent_pair_gap_max) try {
    if (maximum_span == 0) {
        throw std::invalid_argument(
            "species-mismatch realignment maximum_span must be positive");
    }
    if (zero_gap_maximum_span == 0) {
        throw std::invalid_argument(
            "species-mismatch zero-gap maximum span must be positive");
    }
    if (seqpro_managers.size() < 2) {
        spdlog::debug(
            "[species-mismatch-realign] reference={} skipped: "
            "at least two species are required",
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
    AllocationPhaseScope pool_index_phase("pool-index");
    BlockPoolEditor pool_editor(blocks);
    std::set<SpeciesChrPair> affected_paths;

    constexpr size_t rejection_count =
        static_cast<size_t>(MissingWindowReject::COUNT);
    const bool detailed_stats =
        spdlog::should_log(spdlog::level::debug);
    using RejectionCounts = std::array<size_t, rejection_count>;

    RejectionCounts zero_gap_rejections{};
    std::map<SpeciesName, size_t> zero_gap_accepted_by_missing_species;
    std::map<std::string, size_t> zero_gap_pattern_counts;
    std::map<size_t, size_t> zero_gap_boundary_k_counts;
    size_t zero_gap_scan_count = 0;
    size_t zero_gap_candidate_events = 0;
    size_t zero_gap_prepared_events = 0;
    size_t zero_gap_overlap_events = 0;
    size_t zero_gap_replaced_windows = 0;
    size_t zero_gap_replaced_old_blocks = 0;
    size_t zero_gap_adjacent_pair_windows = 0;
    double zero_gap_scan_seconds = 0.0;
    double zero_gap_prepare_seconds = 0.0;
    double zero_gap_commit_seconds = 0.0;
    const size_t zero_gap_iteration_limit = blocks.size() + 1;
    bool zero_gap_full_scan = true;
    std::vector<std::pair<ChrName, SegPtr>> zero_gap_dirty_starts;

    while (true) {
        AllocationPhaseScope iteration_phase("zero-gap-scan");
        ++zero_gap_scan_count;
        if (zero_gap_scan_count > zero_gap_iteration_limit) {
            throw std::runtime_error(
                "Zero-gap fixed-point iteration exceeded the initial "
                "Block-pool bound");
        }

        size_t scanned_this_iteration = 0;
        size_t structural_this_iteration = 0;
        std::vector<MissingWindowPlan> zero_gap_candidates;
        zero_gap_candidates.reserve(
            zero_gap_full_scan ? blocks.size() / 8 + 1
                               : zero_gap_dirty_starts.size() * 4 + 1);
        BlockViewCache zero_gap_view_cache(
            reference_species, zero_gap_scan_count);
        const auto scan_zero_gap_start =
            [&](const ChrName& chromosome,
                const SegPtr& left_reference) {
                if (!left_reference || left_reference->isTail()) {
                    return;
                }
                const auto next_reference =
                    left_reference->primary_path.next.load(
                        std::memory_order_acquire);
                if (!next_reference || next_reference->isTail()) {
                    return;
                }
                BlockView left_view;
                if (!getCachedRealignBlockView(
                        left_reference->parent_block,
                        reference_species, &zero_gap_view_cache,
                        left_view) ||
                    left_view.species_count < 2 ||
                    left_view.species_count > seqpro_managers.size()) {
                    ++zero_gap_rejections[static_cast<size_t>(
                        MissingWindowReject::INVALID_BLOCK)];
                    return;
                }

                const auto right_reference =
                    next_reference->primary_path.next.load(
                        std::memory_order_acquire);
                if (right_reference && !right_reference->isTail()) {
                    ++scanned_this_iteration;
                    MissingWindowCandidate candidate;
                    const auto candidate_result =
                        buildMissingWindowCandidate(
                            left_reference->parent_block,
                            std::vector<BlockPtr>{
                                next_reference->parent_block},
                            right_reference->parent_block,
                            reference_species, chromosome,
                            left_view.species_count, maximum_span,
                            zero_gap_maximum_span,
                            adjacent_pair_gap_max, true, candidate,
                            false, &zero_gap_view_cache);
                    if (candidate_result == MissingWindowReject::COUNT) {
                        ++structural_this_iteration;
                        zero_gap_candidates.push_back(
                            makeMissingWindowPlan(
                                candidate, reference_species));
                    } else {
                        ++zero_gap_rejections[
                            static_cast<size_t>(candidate_result)];
                    }
                }

                MissingWindowCandidate pair_candidate;
                const auto pair_result = buildMissingWindowCandidate(
                    left_reference->parent_block, {},
                    next_reference->parent_block,
                    reference_species, chromosome,
                    left_view.species_count, maximum_span,
                    zero_gap_maximum_span,
                    adjacent_pair_gap_max, true, pair_candidate,
                    false, &zero_gap_view_cache);
                if (pair_result == MissingWindowReject::COUNT) {
                    ++structural_this_iteration;
                    zero_gap_candidates.push_back(
                        makeMissingWindowPlan(
                            pair_candidate, reference_species));
                } else if (
                    pair_result !=
                        MissingWindowReject::ADJACENT_PAIR_REQUIRES_MSA &&
                    pair_result !=
                        MissingWindowReject::ADJACENT_PAIR_ALL_ZERO &&
                    pair_result !=
                        MissingWindowReject::PARTICIPANT_PATTERN) {
                    ++zero_gap_rejections[
                        static_cast<size_t>(pair_result)];
                }
            };

        const auto zero_gap_scan_start =
            std::chrono::steady_clock::now();
        if (zero_gap_full_scan) {
            for (auto& [chromosome, genome_end] :
                 reference_graph_it->second.chr2end) {
                auto left_reference =
                    genome_end.head->primary_path.next.load(
                        std::memory_order_acquire);
                while (left_reference && !left_reference->isTail()) {
                    scan_zero_gap_start(chromosome, left_reference);
                    left_reference =
                        left_reference->primary_path.next.load(
                            std::memory_order_acquire);
                }
            }
        } else {
            std::sort(
                zero_gap_dirty_starts.begin(),
                zero_gap_dirty_starts.end(),
                [](const auto& lhs, const auto& rhs) {
                    if (lhs.first != rhs.first) {
                        return lhs.first < rhs.first;
                    }
                    if (lhs.second->start != rhs.second->start) {
                        return lhs.second->start < rhs.second->start;
                    }
                    return lhs.second.get() < rhs.second.get();
                });
            std::unordered_set<const Segment*> seen_dirty_starts;
            for (const auto& [chromosome, left_reference] :
                 zero_gap_dirty_starts) {
                if (left_reference &&
                    seen_dirty_starts.insert(
                        left_reference.get()).second) {
                    scan_zero_gap_start(chromosome, left_reference);
                }
            }
        }
        zero_gap_scan_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() -
            zero_gap_scan_start).count();
        zero_gap_candidate_events += zero_gap_candidates.size();
        g_allocation_candidates = zero_gap_candidates.size();

        std::stable_sort(
            zero_gap_candidates.begin(), zero_gap_candidates.end(),
            [](const MissingWindowPlan& lhs,
               const MissingWindowPlan& rhs) {
                return lhs.boundary_species_count >
                       rhs.boundary_species_count;
            });

        if (zero_gap_candidates.empty()) {
            spdlog::debug(
                "[species-mismatch-realign][zero-gap] iteration={} "
                "scanned_triples={} structural_candidates={} candidates=0 "
                "status=stable",
                zero_gap_scan_count, scanned_this_iteration,
                structural_this_iteration);
            break;
        }

        size_t prepared_this_iteration = 0;
        size_t overlap_this_iteration = 0;
        std::unordered_set<const Block*> reserved_read_blocks;
        reserved_read_blocks.reserve(
            std::min<size_t>(zero_gap_candidates.size() * 3, 262144));
        std::vector<size_t> remaining(zero_gap_candidates.size());
        std::iota(remaining.begin(), remaining.end(), 0);
        size_t committed_this_iteration = 0;
        size_t replaced_old_this_iteration = 0;
        size_t blocks_before_iteration = pool_editor.activeCount();
        size_t maximum_threads_used = 1;
        double prepare_seconds = 0.0;
        zero_gap_dirty_starts.clear();

        const auto read_blocks = [](const MissingWindowPlan& candidate) {
            std::vector<const Block*> result;
            result.reserve(candidate.interiors.size() + 2);
            result.push_back(candidate.left.get());
            for (const auto& interior : candidate.interiors) {
                result.push_back(interior.get());
            }
            result.push_back(candidate.right.get());
            return result;
        };

        while (!remaining.empty()) {
            std::vector<size_t> batch;
            std::vector<size_t> deferred;
            std::unordered_set<const Block*> batch_reads;
            batch.reserve(std::min<size_t>(
                kZeroGapBatchMaximumWindows, remaining.size()));
            deferred.reserve(remaining.size());
            size_t anchor_records = 0;
            for (const size_t index : remaining) {
                const auto reads = read_blocks(zero_gap_candidates[index]);
                const bool conflicts_committed = std::any_of(
                    reads.begin(), reads.end(), [&](const Block* block) {
                        return block &&
                            reserved_read_blocks.count(block) != 0;
                    });
                if (conflicts_committed) {
                    ++overlap_this_iteration;
                    ++zero_gap_rejections[static_cast<size_t>(
                        MissingWindowReject::OVERLAPPING_PREPARED_WINDOW)];
                    continue;
                }
                const bool conflicts_batch = std::any_of(
                    reads.begin(), reads.end(), [&](const Block* block) {
                        return block && batch_reads.count(block) != 0;
                    });
                const size_t next_anchors =
                    zero_gap_candidates[index].boundary_species_count;
                const bool batch_full =
                    batch.size() >= kZeroGapBatchMaximumWindows ||
                    (!batch.empty() &&
                     anchor_records + next_anchors >
                         kZeroGapBatchMaximumAnchors);
                if (conflicts_batch || batch_full) {
                    deferred.push_back(index);
                    continue;
                }
                batch.push_back(index);
                anchor_records += next_anchors;
                batch_reads.insert(reads.begin(), reads.end());
            }
            remaining = std::move(deferred);
            if (batch.empty()) break;

            const uint_t effective_threads = static_cast<uint_t>(
                std::max<size_t>(
                    1, std::min<size_t>(
                        {std::max<uint_t>(1, parallel_threads),
                         kMaximumPreparationThreads, batch.size()})));
            g_allocation_batch_size = batch.size();
            g_allocation_prepared = 0;
            maximum_threads_used = std::max<size_t>(
                maximum_threads_used, effective_threads);
            std::vector<std::optional<PreparedChain>> prepared_slots(
                batch.size());
            std::vector<std::optional<MissingWindowCandidate>>
                materialized_candidates(batch.size());
            std::vector<MissingWindowReject> prepare_results(
                batch.size(), MissingWindowReject::COUNT);
            std::vector<std::exception_ptr> prepare_exceptions(batch.size());
            const auto prepare_start = std::chrono::steady_clock::now();
            {
                AllocationPhaseScope prepare_phase("zero-gap-preparation");
#pragma omp parallel for schedule(dynamic, 1) num_threads(effective_threads)
                for (std::int64_t batch_index = 0;
                     batch_index < static_cast<std::int64_t>(batch.size());
                     ++batch_index) {
                    const size_t slot = static_cast<size_t>(batch_index);
                    const size_t candidate_index = batch[slot];
                    try {
                        MissingWindowReject result =
                            MissingWindowReject::COUNT;
                        MissingWindowCandidate candidate;
                        result = materializeMissingWindowPlan(
                            zero_gap_candidates[candidate_index],
                            reference_species, maximum_span,
                            zero_gap_maximum_span, adjacent_pair_gap_max,
                            true, candidate);
                        std::optional<PreparedChain> prepared_candidate;
                        if (result == MissingWindowReject::COUNT) {
                            prepared_candidate = candidate.adjacent_pair
                                ? prepareAdjacentPairDirect(
                                      candidate,
                                      reference_species, result)
                                : prepareZeroGapDeletionWindow(
                                      candidate,
                                      reference_species, result);
                        }
                        if (!prepared_candidate.has_value() &&
                            result == MissingWindowReject::COUNT) {
                            result = MissingWindowReject::
                                ZERO_GAP_PREPARATION_INVALID;
                        }
                        prepare_results[slot] = result;
                        prepared_slots[slot] =
                            std::move(prepared_candidate);
                        materialized_candidates[slot] =
                            std::move(candidate);
                    } catch (...) {
                        prepare_exceptions[slot] =
                            std::current_exception();
                    }
                }
            }
            prepare_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - prepare_start).count();
            for (size_t slot = 0; slot < batch.size(); ++slot) {
                if (!prepare_exceptions[slot]) continue;
                for (auto& prepared_slot : prepared_slots) {
                    if (prepared_slot) detachPreparedChain(*prepared_slot);
                }
                std::rethrow_exception(prepare_exceptions[slot]);
            }

            std::vector<MissingWindowCandidate> accepted_candidates;
            std::vector<PreparedChain> prepared;
            accepted_candidates.reserve(batch.size());
            prepared.reserve(batch.size());
            for (size_t slot = 0; slot < batch.size(); ++slot) {
                const size_t candidate_index = batch[slot];
                if (!prepared_slots[slot]) {
                    auto result = prepare_results[slot];
                    if (result == MissingWindowReject::COUNT) {
                        result = MissingWindowReject::
                            ZERO_GAP_PREPARATION_INVALID;
                    }
                    ++zero_gap_rejections[static_cast<size_t>(result)];
                    continue;
                }
                ++prepared_this_iteration;
                ++g_allocation_prepared;
                auto& candidate = *materialized_candidates[slot];
                const auto reads =
                    read_blocks(zero_gap_candidates[candidate_index]);
                reserved_read_blocks.insert(reads.begin(), reads.end());
                if (detailed_stats) {
                    for (const auto& species : candidate.missing_species) {
                        ++zero_gap_accepted_by_missing_species[species];
                    }
                    ++zero_gap_pattern_counts[candidate.participant_pattern];
                    ++zero_gap_boundary_k_counts[
                        candidate.boundary_species_count];
                }
                if (candidate.adjacent_pair) {
                    ++zero_gap_adjacent_pair_windows;
                }
                accepted_candidates.push_back(std::move(candidate));
                prepared.push_back(std::move(*prepared_slots[slot]));
            }
            if (prepared.empty()) continue;

            const auto commit_start = std::chrono::steady_clock::now();
            MissingWindowCommitResult commit;
            {
                AllocationPhaseScope commit_phase("zero-gap-commit-audit");
                commit = commitPreparedMissingWindows(
                    *this, accepted_candidates, prepared, pool_editor,
                    affected_paths,
                    seqpro_managers.size(), "zero-gap", true);
            }
            zero_gap_commit_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - commit_start).count();
            committed_this_iteration += commit.replaced_windows;
            replaced_old_this_iteration += commit.replaced_old_blocks;
            for (const auto& chain : prepared) {
                const auto reference_path = std::find_if(
                    chain.paths.begin(), chain.paths.end(),
                    [&](const PathReplacement& path) {
                        return path.key.first == reference_species;
                    });
                if (reference_path == chain.paths.end() ||
                    !reference_path->merged_segment) {
                    continue;
                }
                auto start = reference_path->merged_segment;
                for (size_t distance = 0; distance < 3 && start;
                     ++distance) {
                    if (!start->isHead() && !start->isTail()) {
                        zero_gap_dirty_starts.emplace_back(
                            chain.candidate.reference_chromosome, start);
                    }
                    start = start->primary_path.prev.load(
                        std::memory_order_acquire);
                }
            }
        }

        zero_gap_prepare_seconds += prepare_seconds;
        zero_gap_prepared_events += prepared_this_iteration;
        zero_gap_overlap_events += overlap_this_iteration;
        if (committed_this_iteration == 0) {
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
                overlap_this_iteration, maximum_threads_used,
                prepare_seconds);
            break;
        }
        zero_gap_full_scan = false;
        zero_gap_replaced_windows += committed_this_iteration;
        zero_gap_replaced_old_blocks += replaced_old_this_iteration;
        spdlog::debug(
            "[species-mismatch-realign][zero-gap] iteration={} "
            "scanned_triples={} structural_candidates={} candidates={} "
            "prepared={} selected={} failed={} "
            "overlap={} parallel_threads={} wall_seconds={:.3f} "
            "replaced_old_blocks={} blocks_before={} blocks_after={} "
            "status=committed",
            zero_gap_scan_count, scanned_this_iteration,
            structural_this_iteration, zero_gap_candidates.size(),
            prepared_this_iteration, committed_this_iteration,
            zero_gap_candidates.size() - prepared_this_iteration,
            overlap_this_iteration, maximum_threads_used, prepare_seconds,
            replaced_old_this_iteration, blocks_before_iteration,
            pool_editor.activeCount());
    }

    RejectionCounts rejections{};
    std::map<SpeciesName, size_t> accepted_by_missing_species;
    if (detailed_stats) {
        accepted_by_missing_species = zero_gap_accepted_by_missing_species;
    }
    std::map<std::string, size_t> ordinary_pattern_counts;
    std::map<std::string, size_t> hybrid_pattern_counts;
    std::map<size_t, size_t> ordinary_candidates_by_k;
    std::map<size_t, size_t> hybrid_candidates_by_k;
    std::map<size_t, size_t> ordinary_committed_by_k;
    std::map<size_t, size_t> hybrid_committed_by_k;
    size_t minipoa_rounds = 0;
    size_t scanned_windows_total = 0;
    size_t ordinary_candidates_total = 0;
    size_t hybrid_candidates_total = 0;
    size_t ordinary_committed_total = 0;
    size_t hybrid_committed_total = 0;
    size_t adjacent_pair_candidates_total = 0;
    size_t adjacent_pair_committed_total = 0;
    size_t adjacent_pair_reference_empty_total = 0;
    size_t full_k_triple_candidates_total = 0;
    size_t full_k_triple_committed_total = 0;
    size_t full_k_triple_fallback_total = 0;
    size_t minipoa_invocations_total = 0;
    size_t reference_only_bypass_total = 0;
    size_t reference_empty_bypass_total = 0;
    size_t preparation_failures_total = 0;
    size_t fallback_batches_total = 0;
    size_t conflict_deferred_total = 0;
    size_t minipoa_replaced_old_blocks = 0;
    double minipoa_scan_seconds = 0.0;
    double minipoa_prepare_seconds = 0.0;
    double minipoa_commit_seconds = 0.0;
    double minipoa_task_seconds = 0.0;
    double minipoa_sequence_fetch_seconds = 0.0;
    double minipoa_msa_task_seconds = 0.0;
    double legacy_minipoa_task_seconds = 0.0;
    double full_k_triple_task_seconds = 0.0;
    uint64_t initial_participant_deficit = 0;
    for (const auto& weak_block : blocks) {
        const auto block = weak_block.lock();
        if (!block || block->anchors.size() > seqpro_managers.size()) {
            continue;
        }
        initial_participant_deficit +=
            seqpro_managers.size() - block->anchors.size();
    }
    uint64_t minipoa_iteration_bound =
        static_cast<uint64_t>(blocks.size());
    if (std::numeric_limits<uint64_t>::max() -
            minipoa_iteration_bound <=
        initial_participant_deficit) {
        minipoa_iteration_bound =
            std::numeric_limits<uint64_t>::max();
    } else {
        minipoa_iteration_bound +=
            initial_participant_deficit + 1;
    }
    const size_t minipoa_iteration_limit =
        minipoa_iteration_bound >
                std::numeric_limits<size_t>::max()
            ? std::numeric_limits<size_t>::max()
            : static_cast<size_t>(minipoa_iteration_bound);
    bool minipoa_full_scan = true;
    std::vector<std::pair<ChrName, SegPtr>> minipoa_dirty_starts;

    while (true) {
        AllocationPhaseScope round_phase("minipoa-scan");
        ++minipoa_rounds;
        if (minipoa_rounds > minipoa_iteration_limit) {
            throw std::runtime_error(
                "Unified minipoa fixed-point iteration exceeded the "
                "initial participant-deficit bound");
        }

        size_t scanned_this_round = 0;
        size_t ordinary_this_round = 0;
        size_t hybrid_this_round = 0;
        std::vector<MissingWindowPlan> candidate_slots;
        candidate_slots.reserve(
            std::min<size_t>(
                minipoa_full_scan ? blocks.size() / 8 + 1
                                  : minipoa_dirty_starts.size() * 4 + 1,
                65536));
        BlockViewCache minipoa_view_cache(
            reference_species, minipoa_rounds);
        const auto scan_minipoa_start =
            [&](const ChrName& chromosome,
                const SegPtr& left_reference) {
                if (!left_reference || left_reference->isTail()) {
                    return;
                }
                const auto next_reference =
                    left_reference->primary_path.next.load(
                        std::memory_order_acquire);
                BlockView left_view;
                if (!getCachedRealignBlockView(
                        left_reference->parent_block,
                        reference_species, &minipoa_view_cache,
                        left_view) ||
                    left_view.species_count < 2 ||
                    left_view.species_count > seqpro_managers.size()) {
                    return;
                }

                if (next_reference && !next_reference->isTail()) {
                    const auto triple_right =
                        next_reference->primary_path.next.load(
                            std::memory_order_acquire);
                    if (triple_right && !triple_right->isTail()) {
                        BlockView middle_view;
                        BlockView right_view;
                        if (getCachedRealignBlockView(
                                next_reference->parent_block,
                                reference_species, &minipoa_view_cache,
                                middle_view) &&
                            getCachedRealignBlockView(
                                triple_right->parent_block,
                                reference_species, &minipoa_view_cache,
                                right_view) &&
                            sameAnchorKeys(left_view, middle_view) &&
                            sameAnchorKeys(left_view, right_view)) {
                            MissingWindowCandidate triple_candidate;
                            const auto triple_result =
                                buildMissingWindowCandidate(
                                    left_reference->parent_block,
                                    std::vector<BlockPtr>{
                                        next_reference->parent_block},
                                    triple_right->parent_block,
                                    reference_species, chromosome,
                                    left_view.species_count, maximum_span,
                                    zero_gap_maximum_span,
                                    adjacent_pair_gap_max, false,
                                    triple_candidate, true,
                                    &minipoa_view_cache);
                            if (triple_result ==
                                MissingWindowReject::COUNT) {
                                ++full_k_triple_candidates_total;
                                ++ordinary_this_round;
                                if (detailed_stats) {
                                    ++ordinary_candidates_by_k[
                                        triple_candidate
                                            .boundary_species_count];
                                }
                                candidate_slots.push_back(
                                    makeMissingWindowPlan(
                                        triple_candidate,
                                        reference_species));
                            } else {
                                ++rejections[static_cast<size_t>(
                                    triple_result)];
                            }
                        }
                    }
                }

                std::vector<BlockPtr> interior_blocks;
                auto cursor = next_reference;
                SegPtr right_reference;
                bool incompatible_boundary = false;
                while (cursor && !cursor->isTail()) {
                    BlockView cursor_view;
                    if (!getCachedRealignBlockView(
                            cursor->parent_block, reference_species,
                            &minipoa_view_cache, cursor_view)) {
                        incompatible_boundary = true;
                        break;
                    }
                    if (sameAnchorKeys(left_view, cursor_view)) {
                        right_reference = cursor;
                        break;
                    }
                    if (cursor_view.species_count >=
                            left_view.species_count ||
                        !participantSpeciesAreSubset(
                            cursor_view, left_view)) {
                        incompatible_boundary = true;
                        break;
                    }
                    interior_blocks.push_back(cursor->parent_block);
                    cursor = cursor->primary_path.next.load(
                        std::memory_order_acquire);
                }
                if (!right_reference) {
                    if (incompatible_boundary) {
                        ++rejections[static_cast<size_t>(
                            MissingWindowReject::
                                INCOMPATIBLE_BOUNDARY_SIGNATURE)];
                    }
                    return;
                }
                ++scanned_this_round;
                MissingWindowCandidate candidate;
                const auto candidate_result =
                    buildMissingWindowCandidate(
                        left_reference->parent_block, interior_blocks,
                        right_reference->parent_block,
                        reference_species, chromosome,
                        left_view.species_count, maximum_span,
                        zero_gap_maximum_span,
                        adjacent_pair_gap_max, false, candidate,
                        false, &minipoa_view_cache);
                if (candidate_result != MissingWindowReject::COUNT) {
                    ++rejections[static_cast<size_t>(candidate_result)];
                    return;
                }
                if (candidate.hybrid_empty) {
                    ++hybrid_this_round;
                    if (detailed_stats) {
                        ++hybrid_candidates_by_k[
                            candidate.boundary_species_count];
                    }
                } else {
                    ++ordinary_this_round;
                    if (detailed_stats) {
                        ++ordinary_candidates_by_k[
                            candidate.boundary_species_count];
                    }
                }
                if (candidate.adjacent_pair) {
                    ++adjacent_pair_candidates_total;
                    const auto reference_path_it = std::find_if(
                        candidate.paths.begin(), candidate.paths.end(),
                        [&](const MissingWindowPathContext& path) {
                            return path.key.first == reference_species;
                        });
                    if (reference_path_it != candidate.paths.end() &&
                        reference_path_it->interval_length == 0) {
                        ++adjacent_pair_reference_empty_total;
                    }
                }
                candidate_slots.push_back(
                    makeMissingWindowPlan(candidate, reference_species));
            };

        const auto minipoa_scan_start =
            std::chrono::steady_clock::now();
        if (minipoa_full_scan) {
            for (auto& [chromosome, genome_end] :
                 reference_graph_it->second.chr2end) {
                auto left_reference =
                    genome_end.head->primary_path.next.load(
                        std::memory_order_acquire);
                while (left_reference && !left_reference->isTail()) {
                    scan_minipoa_start(chromosome, left_reference);
                    left_reference =
                        left_reference->primary_path.next.load(
                            std::memory_order_acquire);
                }
            }
        } else {
            std::sort(
                minipoa_dirty_starts.begin(),
                minipoa_dirty_starts.end(),
                [](const auto& lhs, const auto& rhs) {
                    if (lhs.first != rhs.first) {
                        return lhs.first < rhs.first;
                    }
                    if (lhs.second->start != rhs.second->start) {
                        return lhs.second->start < rhs.second->start;
                    }
                    return lhs.second.get() < rhs.second.get();
                });
            std::unordered_set<const Segment*> seen_dirty_starts;
            for (const auto& [chromosome, left_reference] :
                 minipoa_dirty_starts) {
                if (left_reference &&
                    seen_dirty_starts.insert(
                        left_reference.get()).second) {
                    scan_minipoa_start(chromosome, left_reference);
                }
            }
        }
        minipoa_scan_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() -
            minipoa_scan_start).count();
        scanned_windows_total += scanned_this_round;
        ordinary_candidates_total += ordinary_this_round;
        hybrid_candidates_total += hybrid_this_round;

        if (candidate_slots.empty()) {
            spdlog::debug(
                "[species-mismatch-realign][minipoa-unified] round={} "
                "scanned_windows={} ordinary_candidates=0 "
                "hybrid_candidates=0 status=stable",
                minipoa_rounds, scanned_this_round);
            break;
        }
        g_allocation_candidates = candidate_slots.size();

        std::vector<size_t> remaining(candidate_slots.size());
        std::iota(remaining.begin(), remaining.end(), 0);
        std::stable_sort(
            remaining.begin(), remaining.end(),
            [&](size_t lhs, size_t rhs) {
                const auto priority = [](const MissingWindowPlan& item) {
                    if (item.kind == MissingWindowKind::FULL_K_TRIPLE) {
                        return 0;
                    }
                    if (!item.hybrid_empty && !item.adjacent_pair) {
                        return 1;
                    }
                    if (!item.hybrid_empty && item.adjacent_pair) {
                        return 2;
                    }
                    return 3;
                };
                const int lhs_priority = priority(candidate_slots[lhs]);
                const int rhs_priority = priority(candidate_slots[rhs]);
                if (lhs_priority != rhs_priority) {
                    return lhs_priority < rhs_priority;
                }
                return candidate_slots[lhs].boundary_species_count >
                       candidate_slots[rhs].boundary_species_count;
            });

        std::unordered_set<const Block*> reserved_read_blocks;
        std::unordered_set<const Block*> reserved_write_blocks;
        size_t round_batches = 0;
        size_t round_invocations = 0;
        size_t round_bypasses = 0;
        size_t round_failures = 0;
        size_t round_ordinary_committed = 0;
        size_t round_hybrid_committed = 0;
        size_t round_replaced_old_blocks = 0;
        size_t round_blocks_before = pool_editor.activeCount();
        size_t round_blocks_after = round_blocks_before;
        double round_prepare_seconds = 0.0;
        std::unordered_set<const Block*>
            failed_full_k_triple_blocks;
        minipoa_dirty_starts.clear();

        const auto read_blocks = [](const MissingWindowPlan& candidate) {
            std::vector<const Block*> result;
            result.reserve(candidate.interiors.size() + 2);
            result.push_back(candidate.left.get());
            for (const auto& interior : candidate.interiors) {
                result.push_back(interior.get());
            }
            result.push_back(candidate.right.get());
            return result;
        };
        const auto write_blocks = [](const MissingWindowPlan& candidate) {
            std::vector<const Block*> result;
            const bool replaces_boundaries =
                candidate.hybrid_empty || candidate.adjacent_pair;
            result.reserve(
                candidate.interiors.size() +
                (replaces_boundaries ? 2 : 0));
            if (replaces_boundaries) {
                result.push_back(candidate.left.get());
            }
            for (const auto& interior : candidate.interiors) {
                result.push_back(interior.get());
            }
            if (replaces_boundaries) {
                result.push_back(candidate.right.get());
            }
            return result;
        };

        std::vector<Realignment::PlannerConflictFootprint> footprints;
        footprints.reserve(candidate_slots.size());
        for (const auto& candidate : candidate_slots) {
            footprints.push_back(
                {read_blocks(candidate), write_blocks(candidate)});
        }

        while (!remaining.empty()) {
            const std::vector<const Block*> reserved_reads(
                reserved_read_blocks.begin(), reserved_read_blocks.end());
            const std::vector<const Block*> reserved_writes(
                reserved_write_blocks.begin(), reserved_write_blocks.end());
            const std::vector<size_t> selected =
                Realignment::MissingWindowPlanner::selectConflictFreeBatch(
                    remaining, footprints, reserved_reads, reserved_writes);
            std::vector<size_t> batch;
            batch.reserve(std::min<size_t>(
                kMinipoaBatchMaximumWindows, selected.size()));
            size_t anchor_records = 0;
            for (const size_t index : selected) {
                const size_t next_anchors =
                    candidate_slots[index].boundary_species_count;
                if (batch.size() >= kMinipoaBatchMaximumWindows ||
                    (!batch.empty() &&
                     anchor_records + next_anchors >
                         kMinipoaBatchMaximumAnchors)) {
                    break;
                }
                batch.push_back(index);
                anchor_records += next_anchors;
            }
            const std::unordered_set<size_t> selected_indices(
                batch.begin(), batch.end());
            std::vector<size_t> deferred;
            for (const size_t index : remaining) {
                if (selected_indices.count(index) == 0) {
                    deferred.push_back(index);
                }
            }
            if (batch.empty()) {
                conflict_deferred_total += deferred.size();
                break;
            }
            remaining = std::move(deferred);
            ++round_batches;

            const uint_t effective_threads = static_cast<uint_t>(
                std::max<size_t>(
                    1, std::min<size_t>(
                        {std::max<uint_t>(1, parallel_threads),
                         kMaximumPreparationThreads, batch.size()})));
            AllocationPhaseScope prepare_phase("minipoa-preparation");
            g_allocation_batch_size = batch.size();
            g_allocation_prepared = 0;
            std::vector<std::optional<PreparedChain>> prepared_slots(
                batch.size());
            std::vector<std::optional<MissingWindowCandidate>>
                materialized_candidates(batch.size());
            std::vector<MissingWindowReject> prepare_results(
                batch.size(), MissingWindowReject::COUNT);
            std::vector<std::exception_ptr> prepare_exceptions(
                batch.size());
            std::vector<bool> msa_invoked(batch.size(), false);
            std::vector<bool> reference_only_bypass(
                batch.size(), false);
            std::vector<bool> reference_empty_bypass(
                batch.size(), false);
            std::vector<MissingWindowPreparationTiming>
                preparation_timings(batch.size());
            const auto prepare_start =
                std::chrono::steady_clock::now();
#pragma omp parallel for schedule(dynamic, 1) num_threads(effective_threads)
            for (std::int64_t batch_index = 0;
                 batch_index <
                     static_cast<std::int64_t>(batch.size());
                 ++batch_index) {
                const size_t slot = static_cast<size_t>(batch_index);
                try {
                    MissingWindowReject prepare_result =
                        MissingWindowReject::COUNT;
                    bool invoked = false;
                    bool bypassed = false;
                    bool empty_reference_bypassed = false;
                    MissingWindowCandidate candidate;
                    prepare_result = materializeMissingWindowPlan(
                        candidate_slots[batch[slot]], reference_species,
                        maximum_span, zero_gap_maximum_span,
                        adjacent_pair_gap_max, false, candidate);
                    std::optional<PreparedChain> prepared_candidate;
                    if (prepare_result == MissingWindowReject::COUNT) {
                        prepared_candidate = prepareMissingWindow(
                            candidate, reference_species,
                            seqpro_managers, msa_executable,
                            prepare_result, invoked, bypassed,
                            empty_reference_bypassed,
                            &preparation_timings[slot]);
                    }
                    if (!prepared_candidate.has_value() &&
                        prepare_result == MissingWindowReject::COUNT) {
                        prepare_result = MissingWindowReject::MSA_INVALID;
                    }
                    prepare_results[slot] = prepare_result;
                    prepared_slots[slot] =
                        std::move(prepared_candidate);
                    materialized_candidates[slot] =
                        std::move(candidate);
                    msa_invoked[slot] = invoked;
                    reference_only_bypass[slot] = bypassed;
                    reference_empty_bypass[slot] =
                        empty_reference_bypassed;
                } catch (...) {
                    prepare_exceptions[slot] =
                        std::current_exception();
                }
            }
            const double batch_seconds =
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - prepare_start)
                    .count();
            round_prepare_seconds += batch_seconds;

            for (size_t slot = 0; slot < batch.size(); ++slot) {
                if (!prepare_exceptions[slot]) {
                    continue;
                }
                for (auto& prepared_slot : prepared_slots) {
                    if (prepared_slot.has_value()) {
                        detachPreparedChain(*prepared_slot);
                    }
                }
                spdlog::error(
                    "[species-mismatch-realign][minipoa-unified] "
                    "parallel preparation failed at candidate_index={}",
                    batch[slot]);
                std::rethrow_exception(prepare_exceptions[slot]);
            }

            std::vector<MissingWindowCandidate> accepted_candidates;
            std::vector<PreparedChain> prepared;
            accepted_candidates.reserve(batch.size());
            prepared.reserve(batch.size());
            for (size_t slot = 0; slot < batch.size(); ++slot) {
                minipoa_task_seconds +=
                    preparation_timings[slot].total_seconds;
                minipoa_sequence_fetch_seconds +=
                    preparation_timings[slot].sequence_fetch_seconds;
                minipoa_msa_task_seconds +=
                    preparation_timings[slot].msa_seconds;
                round_invocations += msa_invoked[slot] ? 1 : 0;
                round_bypasses += reference_only_bypass[slot] ? 1 : 0;
                reference_empty_bypass_total +=
                    reference_empty_bypass[slot] ? 1 : 0;
                const size_t candidate_index = batch[slot];
                if (candidate_slots[candidate_index].kind ==
                    MissingWindowKind::FULL_K_TRIPLE) {
                    full_k_triple_task_seconds +=
                        preparation_timings[slot].total_seconds;
                } else {
                    legacy_minipoa_task_seconds +=
                        preparation_timings[slot].total_seconds;
                }
                if (!prepared_slots[slot].has_value()) {
                    auto result = prepare_results[slot];
                    if (result == MissingWindowReject::COUNT) {
                        result = MissingWindowReject::MSA_INVALID;
                    }
                    ++rejections[static_cast<size_t>(result)];
                    ++round_failures;
                    if (candidate_slots[candidate_index].kind ==
                        MissingWindowKind::FULL_K_TRIPLE) {
                        const auto failed_reads = read_blocks(
                            candidate_slots[candidate_index]);
                        failed_full_k_triple_blocks.insert(
                            failed_reads.begin(), failed_reads.end());
                    }
                    continue;
                }

                auto& candidate = *materialized_candidates[slot];
                const auto reads =
                    read_blocks(candidate_slots[candidate_index]);
                const auto writes =
                    write_blocks(candidate_slots[candidate_index]);
                std::unordered_set<const Block*> expected_writes(
                    writes.begin(), writes.end());
                std::unordered_set<const Block*> prepared_writes;
                for (const auto& block :
                     prepared_slots[slot]->candidate.blocks) {
                    if (block) {
                        prepared_writes.insert(block.get());
                    }
                }
                if (expected_writes != prepared_writes) {
                    for (auto& prepared_slot : prepared_slots) {
                        if (prepared_slot.has_value()) {
                            detachPreparedChain(*prepared_slot);
                        }
                    }
                    throw std::runtime_error(
                        "Unified minipoa candidate write set changed "
                        "during preparation");
                }
                reserved_read_blocks.insert(reads.begin(), reads.end());
                reserved_write_blocks.insert(writes.begin(), writes.end());
                if (detailed_stats) {
                    for (const auto& species : candidate.missing_species) {
                        ++accepted_by_missing_species[species];
                    }
                }
                if (candidate.hybrid_empty) {
                    ++round_hybrid_committed;
                    if (detailed_stats) {
                        ++hybrid_committed_by_k[
                            candidate.boundary_species_count];
                        ++hybrid_pattern_counts[
                            candidate.participant_pattern];
                    }
                } else {
                    ++round_ordinary_committed;
                    if (detailed_stats) {
                        ++ordinary_committed_by_k[
                            candidate.boundary_species_count];
                        ++ordinary_pattern_counts[
                            candidate.participant_pattern];
                    }
                }
                if (candidate.kind == MissingWindowKind::FULL_K_TRIPLE) {
                    ++full_k_triple_committed_total;
                }
                if (candidate.adjacent_pair) {
                    if (std::any_of(
                            reads.begin(), reads.end(),
                            [&](const Block* block) {
                                return block &&
                                       failed_full_k_triple_blocks.count(
                                           block) != 0;
                            })) {
                        ++full_k_triple_fallback_total;
                    }
                    ++adjacent_pair_committed_total;
                }
                accepted_candidates.push_back(std::move(candidate));
                prepared.push_back(
                    std::move(*prepared_slots[slot]));
                prepared_slots[slot].reset();
                ++g_allocation_prepared;
            }

            if (!prepared.empty()) {
                const auto commit_start =
                    std::chrono::steady_clock::now();
                MissingWindowCommitResult commit;
                {
                    AllocationPhaseScope commit_phase(
                        "minipoa-commit-audit");
                    commit = commitPreparedMissingWindows(
                        *this, accepted_candidates, prepared, pool_editor,
                        affected_paths,
                        seqpro_managers.size(), "minipoa-unified",
                        false, true);
                }
                minipoa_commit_seconds +=
                    std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - commit_start)
                        .count();
                round_replaced_old_blocks +=
                    commit.replaced_old_blocks;
                round_blocks_after = commit.blocks_after;
                for (const auto& chain : prepared) {
                    const auto reference_path = std::find_if(
                        chain.paths.begin(), chain.paths.end(),
                        [&](const PathReplacement& path) {
                            return path.key.first == reference_species;
                        });
                    if (reference_path == chain.paths.end() ||
                        !reference_path->merged_segment) {
                        continue;
                    }
                    auto start = reference_path->merged_segment;
                    bool is_replacement = true;
                    while (start && !start->isHead()) {
                        if (!start->isTail()) {
                            minipoa_dirty_starts.emplace_back(
                                chain.candidate.reference_chromosome,
                                start);
                        }
                        if (!is_replacement) {
                            BlockView barrier_view;
                            if (buildBlockView(
                                    start->parent_block,
                                    reference_species, barrier_view) &&
                                barrier_view.species_count ==
                                    seqpro_managers.size()) {
                                break;
                            }
                        }
                        start = start->primary_path.prev.load(
                            std::memory_order_acquire);
                        is_replacement = false;
                    }
                }
            }
        }

        minipoa_invocations_total += round_invocations;
        reference_only_bypass_total += round_bypasses;
        preparation_failures_total += round_failures;
        fallback_batches_total += round_batches > 0
                                      ? round_batches - 1
                                      : 0;
        minipoa_prepare_seconds += round_prepare_seconds;

        if (round_ordinary_committed + round_hybrid_committed == 0) {
            spdlog::warn(
                "[species-mismatch-realign][minipoa-unified] round={} "
                "scanned_windows={} ordinary_candidates={} "
                "hybrid_candidates={} batches={} minipoa_calls={} "
                "reference_only_bypass={} failed={} selected=0 "
                "status=stalled",
                minipoa_rounds, scanned_this_round,
                ordinary_this_round, hybrid_this_round,
                round_batches, round_invocations, round_bypasses,
                round_failures);
            break;
        }

        minipoa_full_scan = false;
        ordinary_committed_total += round_ordinary_committed;
        hybrid_committed_total += round_hybrid_committed;
        minipoa_replaced_old_blocks += round_replaced_old_blocks;
        spdlog::debug(
            "[species-mismatch-realign][minipoa-unified] round={} "
            "scanned_windows={} ordinary_candidates={} "
            "hybrid_candidates={} batches={} minipoa_calls={} "
            "reference_only_bypass={} failed={} ordinary_committed={} "
            "hybrid_committed={} replaced_old_blocks={} "
            "blocks_before={} blocks_after={} dirty_starts={} "
            "wall_seconds={:.3f} "
            "status=committed",
            minipoa_rounds, scanned_this_round,
            ordinary_this_round, hybrid_this_round, round_batches,
            round_invocations, round_bypasses, round_failures,
            round_ordinary_committed, round_hybrid_committed,
            round_replaced_old_blocks, round_blocks_before,
            round_blocks_after, minipoa_dirty_starts.size(),
            round_prepare_seconds);
    }

    {
        AllocationPhaseScope audit_phase("final-pool-audit");
        if (!rebuildSamplingAndAuditPaths(
                *this, affected_paths, pool_editor)) {
            throw std::runtime_error(
                "Species-mismatch final path audit failed");
        }
        pool_editor.compact();
        if (!pool_editor.auditPool()) {
            throw std::runtime_error(
                "Species-mismatch final Block pool audit failed");
        }
    }

    const size_t minipoa_committed_total =
        ordinary_committed_total + hybrid_committed_total;
    const double windows_per_second =
        minipoa_prepare_seconds > 0.0
            ? static_cast<double>(minipoa_invocations_total) /
                  minipoa_prepare_seconds
            : 0.0;
    const double minipoa_cigar_task_seconds = std::max(
        0.0, minipoa_task_seconds -
                 minipoa_sequence_fetch_seconds -
                 minipoa_msa_task_seconds);
    spdlog::debug(
        "[species-mismatch-realign][minipoa-unified] reference={} "
        "rounds={} scanned_windows={} ordinary_candidates={} "
        "hybrid_candidates={} minipoa_calls={} "
        "reference_only_bypass={} reference_empty_bypass={} "
        "preparation_failures={} "
        "fallback_batches={} conflict_deferred={} "
        "ordinary_committed={} hybrid_committed={} "
        "scan_seconds={:.3f} prepare_wall_seconds={:.3f} "
        "commit_seconds={:.3f} task_seconds={:.3f} "
        "sequence_fetch_task_seconds={:.3f} msa_task_seconds={:.3f} "
        "cigar_task_seconds={:.3f} legacy_task_seconds={:.3f} "
        "full_k_triple_task_seconds={:.3f} "
        "windows_per_second={:.3f} "
        "maximum_span={}",
        reference_species, minipoa_rounds, scanned_windows_total,
        ordinary_candidates_total, hybrid_candidates_total,
        minipoa_invocations_total, reference_only_bypass_total,
        reference_empty_bypass_total, preparation_failures_total,
        fallback_batches_total,
        conflict_deferred_total, ordinary_committed_total,
        hybrid_committed_total, minipoa_scan_seconds,
        minipoa_prepare_seconds, minipoa_commit_seconds,
        minipoa_task_seconds, minipoa_sequence_fetch_seconds,
        minipoa_msa_task_seconds, minipoa_cigar_task_seconds,
        legacy_minipoa_task_seconds, full_k_triple_task_seconds,
        windows_per_second, maximum_span);
    spdlog::debug(
        "[species-mismatch-realign][zero-gap] reference={} scans={} "
        "candidate_events={} prepared_events={} selected={} "
        "adjacent_pair_direct={} overlap_events={} replaced_old_blocks={} "
        "scan_seconds={:.3f} prepare_wall_seconds={:.3f} "
        "commit_seconds={:.3f} zero_gap_maximum_span={} "
        "adjacent_pair_gap_max={}",
        reference_species, zero_gap_scan_count,
        zero_gap_candidate_events, zero_gap_prepared_events,
        zero_gap_replaced_windows, zero_gap_adjacent_pair_windows,
        zero_gap_overlap_events,
        zero_gap_replaced_old_blocks, zero_gap_scan_seconds,
        zero_gap_prepare_seconds, zero_gap_commit_seconds,
        zero_gap_maximum_span, adjacent_pair_gap_max);
    spdlog::info(
        "[species-mismatch-realign] reference={} replaced={} "
        "zero_gap={} minipoa={} calls={} bypass={} "
        "time(scan/prepare/commit)={:.3f}/{:.3f}/{:.3f}s",
        reference_species,
        zero_gap_replaced_windows + minipoa_committed_total,
        zero_gap_replaced_windows, minipoa_committed_total,
        minipoa_invocations_total,
        reference_only_bypass_total + reference_empty_bypass_total,
        zero_gap_scan_seconds + minipoa_scan_seconds,
        zero_gap_prepare_seconds + minipoa_prepare_seconds,
        zero_gap_commit_seconds + minipoa_commit_seconds);
    spdlog::debug(
        "[species-mismatch-realign][adjacent-pair] reference={} "
        "minipoa_candidates={} minipoa_committed={} "
        "reference_empty_candidates={} maximum_span={}",
        reference_species, adjacent_pair_candidates_total,
        adjacent_pair_committed_total,
        adjacent_pair_reference_empty_total, maximum_span);
    spdlog::debug(
        "[species-mismatch-realign][full-k-triple] reference={} "
        "candidates={} committed={} fallback_to_pair={} maximum_span={}",
        reference_species, full_k_triple_candidates_total,
        full_k_triple_committed_total, full_k_triple_fallback_total,
        maximum_span);

    for (size_t index = 0; index < rejection_count; ++index) {
        if (zero_gap_rejections[index] != 0) {
            spdlog::debug(
                "[species-mismatch-realign][zero-gap] "
                "reject_reason={} count={}",
                missingWindowRejectName(
                    static_cast<MissingWindowReject>(index)),
                zero_gap_rejections[index]);
        }
        if (rejections[index] != 0) {
            spdlog::debug(
                "[species-mismatch-realign][minipoa] "
                "reject_reason={} count={}",
                missingWindowRejectName(
                    static_cast<MissingWindowReject>(index)),
                rejections[index]);
        }
    }
    for (const auto& [species, count] :
         accepted_by_missing_species) {
        spdlog::debug(
            "[species-mismatch-realign] missing_species={} committed={}",
            species, count);
    }
    for (const auto& [pattern, count] : zero_gap_pattern_counts) {
        spdlog::debug(
            "[species-mismatch-realign][zero-gap] pattern={} committed={}",
            pattern, count);
    }
    for (const auto& [boundary_k, count] :
         zero_gap_boundary_k_counts) {
        spdlog::debug(
            "[species-mismatch-realign][zero-gap] boundary_k={} "
            "committed={}",
            boundary_k, count);
    }
    for (const auto& [pattern, count] : ordinary_pattern_counts) {
        spdlog::debug(
            "[species-mismatch-realign][minipoa-ordinary] "
            "pattern={} committed={}",
            pattern, count);
    }
    for (const auto& [pattern, count] : hybrid_pattern_counts) {
        spdlog::debug(
            "[species-mismatch-realign][minipoa-hybrid] "
            "pattern={} committed={}",
            pattern, count);
    }
    for (const auto& [boundary_k, candidates] :
         ordinary_candidates_by_k) {
        spdlog::debug(
            "[species-mismatch-realign][minipoa-ordinary] "
            "boundary_k={} candidates={} committed={}",
            boundary_k, candidates,
            ordinary_committed_by_k[boundary_k]);
    }
    for (const auto& [boundary_k, candidates] :
         hybrid_candidates_by_k) {
        spdlog::debug(
            "[species-mismatch-realign][minipoa-hybrid] "
            "boundary_k={} candidates={} committed={}",
            boundary_k, candidates,
            hybrid_committed_by_k[boundary_k]);
    }

    const size_t total_replaced_windows =
        zero_gap_replaced_windows + minipoa_committed_total;
    const size_t total_replaced_old_blocks =
        zero_gap_replaced_old_blocks +
        minipoa_replaced_old_blocks;
    spdlog::debug(
        "[species-mismatch-realign] reference={} replaced_windows={} "
        "zero_gap_deletion_windows={} ordinary_minipoa_windows={} "
        "hybrid_minipoa_windows={} minipoa_windows={} "
        "replaced_old_blocks={} new_boundary_complete_blocks={}",
        reference_species, total_replaced_windows,
        zero_gap_replaced_windows, ordinary_committed_total,
        hybrid_committed_total, minipoa_committed_total,
        total_replaced_old_blocks, total_replaced_windows);
    return total_replaced_windows;
} catch (const std::bad_alloc&) {
    logAllocationFailure(
        "species-mismatch-realign", reference_species,
        blocks.size(), g_allocation_candidates,
        g_allocation_prepared, g_allocation_batch_size);
    throw;
}

}  // namespace RaMesh
