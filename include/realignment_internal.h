#ifndef RAMAX_REALIGNMENT_INTERNAL_H
#define RAMAX_REALIGNMENT_INTERNAL_H

#include "ramesh.h"

#include <cstdint>
#include <optional>
#include <unordered_map>
#include <vector>

namespace RaMesh::Realignment {

enum class BlockViewProfile : uint8_t {
    ExactMerge,
    Diagnostics,
    MissingWindow
};

struct BlockView {
    BlockPtr block;
    std::map<SpeciesChrPair, SegPtr> anchors;
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

}  // namespace RaMesh::Realignment

#endif
