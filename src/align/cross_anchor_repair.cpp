#include "cross_anchor_repair.h"

#include <spdlog/spdlog.h>

#include <algorithm>

namespace RaMesh::Alignment {

void CrossAnchorRepairCounters::reset() {
    blocks_scanned.store(0);
    cigar_candidates.store(0);
    short_insertions_skipped.store(0);
    same_anchor_skipped.store(0);
    distance_skipped.store(0);
    similarity_pairs.store(0);
    ksw2_repairs.store(0);
    minipoa_calls.store(0);
    cache_hits.store(0);
    accepted.store(0);
    fallback.store(0);
    nanoseconds.store(0);
}

CrossAnchorInsertionRepairSession&
CrossAnchorInsertionRepairSession::instance() {
    static CrossAnchorInsertionRepairSession session;
    return session;
}

void CrossAnchorInsertionRepairSession::configure(
    std::string executable, uint_t maximum_window_span) {
    {
        std::lock_guard lock(configuration_mutex_);
        configuration_.executable = std::move(executable);
        configuration_.maximum_window_span =
            std::max<uint_t>(1, maximum_window_span);
    }
    {
        std::lock_guard lock(cache_mutex);
        msa_cache.clear();
    }
    counters.reset();
}

CrossAnchorRepairConfiguration
CrossAnchorInsertionRepairSession::configurationSnapshot() {
    std::lock_guard lock(configuration_mutex_);
    return configuration_;
}

void CrossAnchorInsertionRepairSession::logStats() const {
    constexpr double kNanosecondsPerSecond = 1.0e9;
    spdlog::info(
        "[cross-anchor-insertion] blocks_scanned={} cigar_candidates={} "
        "short_insertions_skipped={} same_anchor_skipped={} "
        "distance_skipped={} similar_pairs={} ksw2_repaired={} "
        "minipoa_calls={} cache_hits={} accepted={} fallback={} "
        "wall_seconds={:.3f}",
        counters.blocks_scanned.load(), counters.cigar_candidates.load(),
        counters.short_insertions_skipped.load(),
        counters.same_anchor_skipped.load(), counters.distance_skipped.load(),
        counters.similarity_pairs.load(), counters.ksw2_repairs.load(),
        counters.minipoa_calls.load(), counters.cache_hits.load(),
        counters.accepted.load(), counters.fallback.load(),
        counters.nanoseconds.load() / kNanosecondsPerSecond);
}

}  // namespace RaMesh::Alignment
