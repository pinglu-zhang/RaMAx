#include "realignment_internal.h"

#include <set>
#include <shared_mutex>
#include <stdexcept>
#include <unordered_set>

namespace RaMesh::Realignment {

void BlockView::OrderedAnchorRefs::assign(const ChrHeadMap& anchors) {
    entries_.clear();
    entries_.reserve(anchors.size());
    for (const auto& entry : anchors) {
        entries_.push_back(&entry);
    }
    std::sort(
        entries_.begin(), entries_.end(),
        [](const AnchorEntry* left, const AnchorEntry* right) {
            return left->first < right->first;
        });
}

BlockView::OrderedAnchorRefs::const_iterator
BlockView::OrderedAnchorRefs::find(const SpeciesChrPair& key) const {
    const auto found = std::lower_bound(
        entries_.begin(), entries_.end(), key,
        [](const AnchorEntry* entry, const SpeciesChrPair& value) {
            return entry->first < value;
        });
    if (found == entries_.end() || (*found)->first != key) {
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

}  // namespace RaMesh::Realignment
