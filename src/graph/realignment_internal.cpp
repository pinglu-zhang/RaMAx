#include "realignment_internal.h"

#include <set>
#include <shared_mutex>
#include <stdexcept>
#include <unordered_set>

namespace RaMesh::Realignment {
namespace {

constexpr size_t kMaximumBlockViewCacheEntries = 65536;

}  // namespace

BlockView::OrderedAnchorRefs::OrderedAnchorRefs()
    : entries_(std::make_shared<std::vector<const AnchorEntry*>>()) {}

void BlockView::OrderedAnchorRefs::assign(const ChrHeadMap& anchors) {
    auto entries = std::make_shared<std::vector<const AnchorEntry*>>();
    entries->reserve(anchors.size());
    for (const auto& entry : anchors) {
        entries->push_back(&entry);
    }
    std::sort(
        entries->begin(), entries->end(),
        [](const AnchorEntry* left, const AnchorEntry* right) {
            return left->first < right->first;
        });
    entries_ = std::move(entries);
}

void BlockView::OrderedAnchorRefs::clear() {
    entries_ = std::make_shared<std::vector<const AnchorEntry*>>();
}

BlockView::OrderedAnchorRefs::const_iterator
BlockView::OrderedAnchorRefs::find(const SpeciesChrPair& key) const {
    const auto found = std::lower_bound(
        entries_->begin(), entries_->end(), key,
        [](const AnchorEntry* entry, const SpeciesChrPair& value) {
            return entry->first < value;
        });
    if (found == entries_->end() || (*found)->first != key) {
        return end();
    }
    return const_iterator(found);
}

const SegPtr& BlockView::OrderedAnchorRefs::at(
    const SpeciesChrPair& key) const {
    const auto found = find(key);
    if (found == end()) {
        throw std::out_of_range("BlockView anchor key is absent");
    }
    return found->second;
}

BlockViewBuilder::BlockViewBuilder(SpeciesName reference_species,
                                   uint64_t graph_version)
    : reference_species_(std::move(reference_species)),
      graph_version_(graph_version) {}

size_t BlockViewBuilder::CacheKeyHash::operator()(
    const CacheKey& key) const noexcept {
    size_t result = std::hash<const Block*>{}(key.block);
    result ^= static_cast<size_t>(key.profile) + 0x9e3779b9U +
        (result << 6U) + (result >> 2U);
    const size_t chromosome =
        std::hash<ChrName>{}(key.diagnostic_reference_chromosome);
    return result ^ (chromosome + 0x9e3779b9U +
                     (result << 6U) + (result >> 2U));
}

void BlockViewBuilder::clear(uint64_t graph_version) {
    graph_version_ = graph_version;
    cache_.clear();
}

bool BlockViewBuilder::build(
    const BlockPtr& block,
    BlockViewProfile profile,
    BlockView& view,
    const ChrName& diagnostic_reference_chromosome) {
    (void)graph_version_;
    if (!block) return false;
    CacheKey key{block.get(), profile, diagnostic_reference_chromosome};
    const auto cached = cache_.find(key);
    if (cached != cache_.end()) {
        if (!cached->second) return false;
        view = *cached->second;
        return true;
    }
    if (cache_.size() >= kMaximumBlockViewCacheEntries) {
        cache_.clear();
    }
    BlockView built;
    if (!buildUncached(block, profile, built,
                       diagnostic_reference_chromosome)) {
        cache_.emplace(std::move(key), std::nullopt);
        return false;
    }
    view = built;
    cache_.emplace(std::move(key), std::move(built));
    return true;
}

bool BlockViewBuilder::buildUncached(
    const BlockPtr& block,
    BlockViewProfile profile,
    BlockView& view,
    const ChrName& diagnostic_reference_chromosome) const {
    view = {};
    view.block = block;
    std::set<SpeciesName> species;
    std::shared_lock lock(block->rw);

    const size_t minimum_anchors =
        profile == BlockViewProfile::MissingWindow ? 1 : 2;
    if (block->anchors.size() < minimum_anchors) return false;
    if (profile != BlockViewProfile::Diagnostics && block->ref_chr.empty()) {
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
    }
    view.anchors.assign(block->anchors);

    const ChrName& reference_chromosome =
        profile == BlockViewProfile::Diagnostics
            ? diagnostic_reference_chromosome
            : block->ref_chr;
    const auto reference = view.anchors.find(
        SpeciesChrPair{reference_species_, reference_chromosome});
    if (reference == view.anchors.end() ||
        reference->second->strand != Strand::FORWARD) {
        return false;
    }
    view.reference_segment = reference->second;
    view.species_count = species.size();
    return profile == BlockViewProfile::MissingWindow ||
        view.species_count >= 2;
}

std::vector<size_t> MissingWindowPlanner::selectConflictFreeBatch(
    const std::vector<size_t>& ordered_candidates,
    const std::vector<PlannerConflictFootprint>& footprints,
    const std::vector<const Block*>& reserved_reads,
    const std::vector<const Block*>& reserved_writes) {
    std::unordered_set<const Block*> reads(
        reserved_reads.begin(), reserved_reads.end());
    std::unordered_set<const Block*> writes(
        reserved_writes.begin(), reserved_writes.end());
    std::vector<size_t> selected;
    selected.reserve(ordered_candidates.size());
    for (const size_t index : ordered_candidates) {
        if (index >= footprints.size()) continue;
        const auto& footprint = footprints[index];
        bool conflict = false;
        for (const Block* block : footprint.writes) {
            if (!block || reads.count(block) != 0) {
                conflict = true;
                break;
            }
        }
        if (!conflict) {
            for (const Block* block : footprint.reads) {
                if (!block || writes.count(block) != 0) {
                    conflict = true;
                    break;
                }
            }
        }
        if (conflict) continue;
        selected.push_back(index);
        reads.insert(footprint.reads.begin(), footprint.reads.end());
        writes.insert(footprint.writes.begin(), footprint.writes.end());
    }
    return selected;
}

MissingWindowBatchScheduler::MissingWindowBatchScheduler(
    std::vector<size_t> ordered_candidates,
    const std::vector<PlannerConflictFootprint>& footprints,
    const std::vector<size_t>& anchor_counts)
    : footprints_(&footprints), anchor_counts_(&anchor_counts),
      pending_(ordered_candidates.begin(), ordered_candidates.end()) {
    if (anchor_counts.size() != footprints.size()) {
        throw std::invalid_argument(
            "Missing-window scheduler footprint/count size mismatch");
    }
}

bool MissingWindowBatchScheduler::conflicts(
    const PlannerConflictFootprint& footprint,
    const std::unordered_set<const Block*>& reads,
    const std::unordered_set<const Block*>& writes) {
    for (const Block* block : footprint.writes) {
        if (!block || reads.count(block) != 0) return true;
    }
    for (const Block* block : footprint.reads) {
        if (!block || writes.count(block) != 0) return true;
    }
    return false;
}

std::vector<size_t> MissingWindowBatchScheduler::nextBatch(
    size_t maximum_windows,
    size_t maximum_anchors,
    bool defer_anchor_overflow) {
    if (batch_active_) {
        throw std::logic_error(
            "Missing-window scheduler batch was not completed");
    }
    if (maximum_windows == 0 || maximum_anchors == 0) {
        throw std::invalid_argument(
            "Missing-window scheduler limits must be positive");
    }

    batch_deferred_.clear();
    std::unordered_set<const Block*> batch_reads;
    std::unordered_set<const Block*> batch_writes;
    std::vector<size_t> batch;
    batch.reserve(std::min(maximum_windows, pending_.size()));
    size_t anchors = 0;

    while (!pending_.empty()) {
        const size_t index = pending_.front();
        pending_.pop_front();
        if (index >= footprints_->size()) {
            ++dropped_conflicts_;
            continue;
        }
        const auto& footprint = (*footprints_)[index];
        if (conflicts(footprint, committed_reads_, committed_writes_)) {
            ++dropped_conflicts_;
            continue;
        }
        if (conflicts(footprint, batch_reads, batch_writes)) {
            batch_deferred_.push_back(index);
            continue;
        }

        const size_t next_anchors = (*anchor_counts_)[index];
        if (!batch.empty() && batch.size() >= maximum_windows) {
            pending_.push_front(index);
            break;
        }
        if (!batch.empty() &&
            anchors + next_anchors > maximum_anchors) {
            if (defer_anchor_overflow) {
                batch_deferred_.push_back(index);
                continue;
            }
            pending_.push_front(index);
            break;
        }

        batch.push_back(index);
        anchors += next_anchors;
        batch_reads.insert(footprint.reads.begin(), footprint.reads.end());
        batch_writes.insert(footprint.writes.begin(), footprint.writes.end());
    }

    batch_active_ = !batch.empty();
    if (!batch_active_ && !batch_deferred_.empty()) {
        throw std::logic_error(
            "Missing-window scheduler deferred candidates without a blocker");
    }
    return batch;
}

void MissingWindowBatchScheduler::completeBatch(
    const std::vector<size_t>& batch,
    const std::vector<bool>& succeeded) {
    if (!batch_active_ || batch.size() != succeeded.size()) {
        throw std::logic_error(
            "Missing-window scheduler completed an invalid batch");
    }
    for (size_t slot = 0; slot < batch.size(); ++slot) {
        if (!succeeded[slot]) continue;
        const size_t index = batch[slot];
        if (index >= footprints_->size()) {
            throw std::logic_error(
                "Missing-window scheduler success index is invalid");
        }
        const auto& footprint = (*footprints_)[index];
        committed_reads_.insert(
            footprint.reads.begin(), footprint.reads.end());
        committed_writes_.insert(
            footprint.writes.begin(), footprint.writes.end());
    }

    retried_candidates_ += batch_deferred_.size();
    for (auto iterator = batch_deferred_.rbegin();
         iterator != batch_deferred_.rend(); ++iterator) {
        pending_.push_front(*iterator);
    }
    batch_deferred_.clear();
    batch_active_ = false;
}

bool MissingWindowBatchScheduler::empty() const noexcept {
    return pending_.empty() && batch_deferred_.empty() && !batch_active_;
}

size_t MissingWindowBatchScheduler::droppedConflicts() const noexcept {
    return dropped_conflicts_;
}

size_t MissingWindowBatchScheduler::retriedCandidates() const noexcept {
    return retried_candidates_;
}

size_t MissingWindowBatchScheduler::pendingCandidates() const noexcept {
    return pending_.size() + batch_deferred_.size();
}

}  // namespace RaMesh::Realignment
