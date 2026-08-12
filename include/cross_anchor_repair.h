#ifndef RAMAX_CROSS_ANCHOR_REPAIR_H
#define RAMAX_CROSS_ANCHOR_REPAIR_H

#include "align.h"

#include <atomic>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace RaMesh::Alignment {

struct CrossAnchorRepairConfiguration {
    std::string executable;
    uint_t maximum_window_span = 3000;
};

struct CrossAnchorRepairCounters {
    std::atomic<uint64_t> blocks_scanned{0};
    std::atomic<uint64_t> cigar_candidates{0};
    std::atomic<uint64_t> short_insertions_skipped{0};
    std::atomic<uint64_t> same_anchor_skipped{0};
    std::atomic<uint64_t> distance_skipped{0};
    std::atomic<uint64_t> similarity_pairs{0};
    std::atomic<uint64_t> ksw2_repairs{0};
    std::atomic<uint64_t> minipoa_calls{0};
    std::atomic<uint64_t> cache_hits{0};
    std::atomic<uint64_t> accepted{0};
    std::atomic<uint64_t> fallback{0};
    std::atomic<uint64_t> nanoseconds{0};

    void reset();
};

class CrossAnchorInsertionRepairSession {
public:
    static CrossAnchorInsertionRepairSession& instance();

    void configure(std::string executable, uint_t maximum_window_span);
    CrossAnchorRepairConfiguration configurationSnapshot();
    void logStats() const;

    std::mutex cache_mutex;
    std::unordered_map<
        std::string,
        std::vector<std::pair<ChrName, std::string>>>
        msa_cache;
    CrossAnchorRepairCounters counters;

private:
    std::mutex configuration_mutex_;
    CrossAnchorRepairConfiguration configuration_;
};

}  // namespace RaMesh::Alignment

#endif
