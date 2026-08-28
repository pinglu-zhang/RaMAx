#ifndef RAMAX_GRAPH_MERGE_INTERNAL_H
#define RAMAX_GRAPH_MERGE_INTERNAL_H

#include "align.h"

#include <cstdint>

namespace RaMesh::detail {

struct MergeCigarPiece {
  Cigar_t cigar;
  uint32_t reference_length = 0;
  uint32_t query_length = 0;
};

struct MergeCigarSplit {
  bool has_prefix = false;
  bool has_suffix = false;
  MergeCigarPiece prefix;
  MergeCigarPiece overlap;
  MergeCigarPiece suffix;
  uint64_t source_units_scanned = 0;
};

// Split a query CIGAR according to reference-consuming lengths without a
// text round trip. This intentionally preserves the legacy merge rules:
// I does not consume reference, every non-D operation consumes query, and
// a terminal query-only remainder is appended to the final piece with the
// same boundary coalescing used by appendCigar().
MergeCigarSplit splitCigarForOverlapMerge(
    const Cigar_t& source,
    bool has_prefix,
    uint32_t prefix_reference_length,
    uint32_t overlap_reference_length,
    bool has_suffix,
    uint32_t suffix_reference_length);

}  // namespace RaMesh::detail

#endif  // RAMAX_GRAPH_MERGE_INTERNAL_H
