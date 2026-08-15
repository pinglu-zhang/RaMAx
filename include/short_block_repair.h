#ifndef RAMAX_SHORT_BLOCK_REPAIR_H
#define RAMAX_SHORT_BLOCK_REPAIR_H

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "SeqPro.h"
#include "config.hpp"

namespace RaMesh {
class RaMeshMultiGenomeGraph;
}

namespace RaMesh::ShortBlockRepair {

struct Options {
    bool enabled = false;
    uint_t maximum_short_length = 500;
    // Zero disables deletion of Blocks that cannot be merged.
    uint_t maximum_delete_length = 0;
    uint_t maximum_missing_span = 200;
    double minimum_coverage = 0.70;
    double minimum_identity = 0.60;
    uint_t maximum_query_gap = 100;
    uint_t parallel_threads = 1;
};

struct Result {
    uint64_t reference_paths = 0;
    uint64_t scanned_blocks = 0;
    uint64_t unique_short_blocks = 0;
    uint64_t left_candidates = 0;
    uint64_t right_candidates = 0;
    uint64_t participant_rejected = 0;
    uint64_t path_rejected = 0;
    uint64_t strand_rejected = 0;
    uint64_t span_rejected = 0;
    uint64_t cigar_rejected = 0;
    uint64_t ksw2_calls = 0;
    uint64_t ksw2_passed = 0;
    uint64_t ksw2_failed = 0;
    uint64_t left_merged = 0;
    uint64_t right_merged = 0;
    uint64_t fixed_point_generations = 0;
    uint64_t deleted_blocks = 0;
    uint64_t deleted_segments = 0;
    uint64_t transaction_rollbacks = 0;
    std::map<SpeciesName, uint64_t> deleted_bases_by_species;
    uint64_t blocks_before = 0;
    uint64_t blocks_after = 0;
    uint64_t blocks_le_10_before = 0;
    uint64_t blocks_le_50_before = 0;
    uint64_t blocks_le_100_before = 0;
    uint64_t blocks_le_500_before = 0;
    uint64_t blocks_le_10_after = 0;
    uint64_t blocks_le_50_after = 0;
    uint64_t blocks_le_100_after = 0;
    uint64_t blocks_le_500_after = 0;
    double scan_seconds = 0.0;
    double ksw2_seconds = 0.0;
    double transaction_seconds = 0.0;
    double total_seconds = 0.0;
};

Result repairFinalShortBlocks(
    RaMeshMultiGenomeGraph& graph,
    const std::vector<SpeciesName>& reference_order,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    const Options& options);

}  // namespace RaMesh::ShortBlockRepair

#endif
