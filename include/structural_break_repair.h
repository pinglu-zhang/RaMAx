#ifndef RAMAX_STRUCTURAL_BREAK_REPAIR_H
#define RAMAX_STRUCTURAL_BREAK_REPAIR_H

#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>

#include "SeqPro.h"
#include "config.hpp"

namespace RaMesh {
class RaMeshMultiGenomeGraph;
}

namespace RaMesh::StructuralBreakRepair {

class FailureCache {
public:
    bool contains(const std::string& signature) const;
    void remember(std::string signature, std::string reason);
    void clear();

private:
    mutable std::mutex mutex_;
    std::unordered_map<std::string, std::string> failures_;
};

// Precision-first structural repair settings are independent of missing-window
// settings: a structural discontinuity is only a recall signal and never
// sufficient by itself to edit the graph.
struct Options {
    bool enabled = false;
    uint_t maximum_span = 1000;
    uint_t minimum_outer_anchor = 100;
    uint_t strong_outer_anchor = 500;
    uint_t maximum_interior_blocks = 5;
    double minimum_coverage = 0.70;
    double minimum_identity = 0.60;
    double maximum_anchor_identity_drop = 0.15;
    uint_t parallel_threads = 1;
    std::string msa_executable;
    std::shared_ptr<FailureCache> failure_cache =
        std::make_shared<FailureCache>();
};

struct Result {
    uint64_t scanned_windows = 0;
    uint64_t structural_candidates = 0;
    uint64_t conflict_deferred = 0;
    uint64_t sequence_prepared = 0;
    uint64_t msa_calls = 0;
    uint64_t msa_failed = 0;
    uint64_t quality_rejected = 0;
    uint64_t transaction_rejected = 0;
    uint64_t committed = 0;
    uint64_t failure_cache_hits = 0;

    uint64_t outer_anchor_invalid = 0;
    uint64_t span_exceeded = 0;
    uint64_t reference_empty = 0;
    uint64_t too_many_interior_blocks = 0;
    uint64_t large_gap_only = 0;
    uint64_t unsupported_participant_count = 0;

    uint64_t target_switch = 0;
    uint64_t strand_switch = 0;
    uint64_t order_break = 0;
    std::map<size_t, uint64_t> candidates_by_k;
    std::map<size_t, uint64_t> prepared_by_k;
    std::map<size_t, uint64_t> committed_by_k;

    double scan_seconds = 0.0;
    double sequence_seconds = 0.0;
    double msa_seconds = 0.0;
    double quality_seconds = 0.0;
    double transaction_seconds = 0.0;
    std::map<std::string, uint64_t> rejection_reasons;
};

// Runs after exact/missing-window cleanup and before masking. The graph is
// unchanged for every rejected or failed candidate.
Result repairAnchorBoundedStructuralBreaks(
    RaMeshMultiGenomeGraph& graph,
    const SpeciesName& reference_species,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    const Options& options);

}  // namespace RaMesh::StructuralBreakRepair

#endif  // RAMAX_STRUCTURAL_BREAK_REPAIR_H
