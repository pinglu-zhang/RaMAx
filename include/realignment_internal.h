#ifndef RAMAX_REALIGNMENT_INTERNAL_H
#define RAMAX_REALIGNMENT_INTERNAL_H

#include "ramesh.h"

#include <compare>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <iterator>
#include <memory>
#include <optional>
#include <unordered_set>
#include <unordered_map>
#include <vector>

namespace RaMesh::Realignment {

enum class BlockViewProfile : uint8_t {
    ExactMerge,
    Diagnostics,
    MissingWindow
};

struct BlockView {
    using AnchorEntry = ChrHeadMap::value_type;

    class OrderedAnchorRefs {
    public:
        class const_iterator {
        public:
            using iterator_category = std::random_access_iterator_tag;
            using value_type = AnchorEntry;
            using difference_type = std::ptrdiff_t;
            using pointer = const AnchorEntry*;
            using reference = const AnchorEntry&;

            const_iterator() = default;
            explicit const_iterator(
                std::vector<const AnchorEntry*>::const_iterator iterator)
                : iterator_(iterator) {}

            reference operator*() const { return **iterator_; }
            pointer operator->() const { return *iterator_; }
            const_iterator& operator++() { ++iterator_; return *this; }
            const_iterator operator++(int) {
                auto copy = *this;
                ++*this;
                return copy;
            }
            const_iterator& operator--() { --iterator_; return *this; }
            const_iterator operator--(int) {
                auto copy = *this;
                --*this;
                return copy;
            }
            const_iterator& operator+=(difference_type offset) {
                iterator_ += offset;
                return *this;
            }
            const_iterator& operator-=(difference_type offset) {
                iterator_ -= offset;
                return *this;
            }
            friend const_iterator operator+(
                const_iterator iterator, difference_type offset) {
                iterator += offset;
                return iterator;
            }
            friend const_iterator operator-(
                const_iterator iterator, difference_type offset) {
                iterator -= offset;
                return iterator;
            }
            friend difference_type operator-(
                const const_iterator& left, const const_iterator& right) {
                return left.iterator_ - right.iterator_;
            }
            friend bool operator==(
                const const_iterator&, const const_iterator&) = default;
            friend auto operator<=> (
                const const_iterator&, const const_iterator&) = default;

        private:
            std::vector<const AnchorEntry*>::const_iterator iterator_{};
        };

        OrderedAnchorRefs();
        void assign(const ChrHeadMap& anchors);
        void clear();
        [[nodiscard]] size_t size() const noexcept {
            return entries_->size();
        }
        [[nodiscard]] bool empty() const noexcept {
            return entries_->empty();
        }
        [[nodiscard]] const_iterator begin() const {
            return const_iterator(entries_->begin());
        }
        [[nodiscard]] const_iterator end() const {
            return const_iterator(entries_->end());
        }
        [[nodiscard]] const_iterator find(const SpeciesChrPair& key) const;
        [[nodiscard]] size_t count(const SpeciesChrPair& key) const {
            return find(key) == end() ? 0U : 1U;
        }
        [[nodiscard]] const SegPtr& at(const SpeciesChrPair& key) const;

    private:
        std::shared_ptr<std::vector<const AnchorEntry*>> entries_;
    };

    BlockPtr block;
    OrderedAnchorRefs anchors;
    SegPtr reference_segment;
    ChrName declared_reference_chromosome;
    size_t species_count = 0;
};

class BlockViewBuilder {
public:
    explicit BlockViewBuilder(SpeciesName reference_species,
                              uint64_t graph_version = 0);

    bool build(const BlockPtr& block,
               BlockViewProfile profile,
               BlockView& view,
               const ChrName& diagnostic_reference_chromosome = {});
    void clear(uint64_t graph_version);

private:
    struct CacheKey {
        const Block* block = nullptr;
        BlockViewProfile profile = BlockViewProfile::MissingWindow;
        ChrName diagnostic_reference_chromosome;

        bool operator==(const CacheKey&) const = default;
    };

    struct CacheKeyHash {
        size_t operator()(const CacheKey& key) const noexcept;
    };

    bool buildUncached(const BlockPtr& block,
                       BlockViewProfile profile,
                       BlockView& view,
                       const ChrName& diagnostic_reference_chromosome) const;

    SpeciesName reference_species_;
    uint64_t graph_version_ = 0;
    std::unordered_map<CacheKey, std::optional<BlockView>, CacheKeyHash>
        cache_;
};

struct PlannerConflictFootprint {
    std::vector<const Block*> reads;
    std::vector<const Block*> writes;
};

class MissingWindowPlanner {
public:
    static std::vector<size_t> selectConflictFreeBatch(
        const std::vector<size_t>& ordered_candidates,
        const std::vector<PlannerConflictFootprint>& footprints,
        const std::vector<const Block*>& reserved_reads = {},
        const std::vector<const Block*>& reserved_writes = {});
};

// Stateful equivalent of repeatedly calling selectConflictFreeBatch().  The
// old fixed-size batching code rebuilt and rescanned every remaining candidate
// after each batch, which becomes quadratic on multi-million-Block graphs.
// This scheduler retains the same stable greedy order while visiting a
// candidate again only when an earlier conflicting candidate failed.
class MissingWindowBatchScheduler {
public:
    MissingWindowBatchScheduler(
        std::vector<size_t> ordered_candidates,
        const std::vector<PlannerConflictFootprint>& footprints,
        const std::vector<size_t>& anchor_counts);

    [[nodiscard]] std::vector<size_t> nextBatch(
        size_t maximum_windows,
        size_t maximum_anchors,
        bool defer_anchor_overflow = false);

    void completeBatch(
        const std::vector<size_t>& batch,
        const std::vector<bool>& succeeded);

    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] size_t droppedConflicts() const noexcept;
    [[nodiscard]] size_t retriedCandidates() const noexcept;
    [[nodiscard]] size_t pendingCandidates() const noexcept;

private:
    [[nodiscard]] static bool conflicts(
        const PlannerConflictFootprint& footprint,
        const std::unordered_set<const Block*>& reads,
        const std::unordered_set<const Block*>& writes);

    const std::vector<PlannerConflictFootprint>* footprints_ = nullptr;
    const std::vector<size_t>* anchor_counts_ = nullptr;
    std::deque<size_t> pending_;
    std::vector<size_t> batch_deferred_;
    std::unordered_set<const Block*> committed_reads_;
    std::unordered_set<const Block*> committed_writes_;
    bool batch_active_ = false;
    size_t dropped_conflicts_ = 0;
    size_t retried_candidates_ = 0;
};

}  // namespace RaMesh::Realignment

#endif
