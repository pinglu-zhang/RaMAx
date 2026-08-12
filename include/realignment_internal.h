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

enum class LocalRepairKind : uint8_t {
    ExistingMissingWindow,
    StructuralBreak
};

struct StablePriorityKey {
    size_t participant_count = 0;
    uint_t reference_start = 0;
    uint64_t serial = 0;
};

struct LocalRepairCandidate {
    LocalRepairKind kind = LocalRepairKind::ExistingMissingWindow;
    ChrName reference_chromosome;
    uint_t reference_start = 0;
    uint_t reference_length = 0;
    StablePriorityKey priority;
    std::vector<BlockPtr> read_blocks;
    std::vector<BlockPtr> replaced_blocks;
    std::vector<SpeciesChrPair> affected_paths;
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

struct PreparedPathReplacement {
    SpeciesChrPair key;
    std::vector<SegPtr> old_segments;
    SegPtr replacement;
};

struct PreparedGraphReplacement {
    BlockPtr provisional_block;
    // Structural-break repair can detach an invalid homology edge while
    // retaining the original genomic Segment in a residual singleton Block.
    // Missing-window transactions leave this vector empty.
    std::vector<BlockPtr> residual_blocks;
    std::vector<PreparedPathReplacement> paths;
    std::vector<BlockPtr> replaced_blocks;
    std::vector<SpeciesChrPair> audit_scope;
    uint64_t progress_before = 0;
    uint64_t progress_after = 0;
};

}  // namespace RaMesh::Realignment

#endif
