#ifndef RAMAX_PAF_EXPORT_H
#define RAMAX_PAF_EXPORT_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace RaMesh::Paf {

enum class Mode {
    CONNECTED,
    ALL
};

struct PafExportOptions {
    Mode mode = Mode::CONNECTED;
    bool only_primary = true;
};

struct PafExportStats {
    std::size_t blocks_seen = 0;
    std::size_t blocks_exported = 0;
    std::size_t blocks_expired = 0;
    std::size_t blocks_too_small = 0;
    std::size_t invalid_blocks_skipped = 0;
    std::size_t fallback_blocks = 0;
    std::size_t theoretical_all_pairs = 0;
    std::size_t eligible_all_pairs = 0;
    std::size_t base_pairs = 0;
    std::size_t supplemental_pairs = 0;
    std::size_t records_written = 0;
    std::size_t reverse_records = 0;
    std::uint64_t matching_bases = 0;
    std::uint64_t alignment_columns = 0;
    double elapsed_seconds = 0.0;
    std::string first_invalid_reference;
    std::string first_invalid_reason;
};

using SequencePair = std::pair<std::size_t, std::size_t>;

struct PairSelectionResult {
    std::vector<SequencePair> pairs;
    std::size_t base_pairs = 0;
    std::size_t supplemental_pairs = 0;
    std::size_t theoretical_pairs = 0;
    std::size_t eligible_pairs = 0;
    bool verified = false;
};

struct PairProjection {
    std::string cigar;
    std::uint64_t matching_bases = 0;
    std::uint64_t block_length = 0;
    std::uint64_t query_consumed = 0;
    std::uint64_t target_consumed = 0;
    std::uint64_t edit_distance = 0;
    std::uint64_t overlap_columns = 0;
    bool valid = false;
};

PairSelectionResult selectPairs(
    const std::vector<std::string>& aligned_rows,
    const std::vector<std::string>& names,
    std::size_t reference_index,
    Mode mode);

bool verifyColumnConnectivity(
    const std::vector<std::string>& aligned_rows,
    const std::vector<SequencePair>& pairs);

PairProjection projectPair(
    const std::string& target_aligned,
    const std::string& query_aligned,
    bool reverse_columns);

}  // namespace RaMesh::Paf

#endif
