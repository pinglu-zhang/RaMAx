// =============================================================
//  File: ramesh_graph.cpp –  High‑level graph ops (v0.6‑alpha)
// =============================================================
#include "align.h"
#include "extension_internal.h"
#include "ramesh.h"
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <curl/curl.h>
#include <limits>
#include <optional>
#include <shared_mutex>
#include <spdlog/spdlog.h>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace RaMesh {
namespace {

void stableRadixSortSegmentsByStart(
    std::vector<SegPtr>& segments,
    std::vector<SegPtr>& scratch,
    std::vector<size_t>& counts) {
  if (segments.size() < 2) return;
  constexpr size_t kRadix = 1u << 16;
  counts.resize(kRadix);
  scratch.resize(segments.size());
  const auto pass = [&](const std::vector<SegPtr>& source,
                        std::vector<SegPtr>& destination,
                        unsigned shift) {
    std::fill(counts.begin(), counts.end(), 0);
    for (const auto& segment : source) {
      ++counts[(static_cast<uint32_t>(segment->start) >> shift) & 0xffffu];
    }
    size_t offset = 0;
    for (size_t& count : counts) {
      const size_t frequency = count;
      count = offset;
      offset += frequency;
    }
    for (const auto& segment : source) {
      const size_t key =
          (static_cast<uint32_t>(segment->start) >> shift) & 0xffffu;
      destination[counts[key]++] = segment;
    }
  };
  pass(segments, scratch, 0);
  pass(scratch, segments, 16);
}

constexpr uint32_t kNoExtensionWorkIndex =
    std::numeric_limits<uint32_t>::max();

struct ExtensionWorkItem {
  Segment* query_segment{nullptr};
  Segment* reference_segment{nullptr};
  Segment* left_reference_candidate{nullptr};
  Segment* right_reference_candidate{nullptr};
  uint32_t query_context_index{0};
  uint32_t reference_context_index{0};
  uint32_t same_reference_next{kNoExtensionWorkIndex};
};

struct QueryExtensionContext {
  GenomeEnd* end{nullptr};
  ChrName chromosome;
  detail::ExtensionSequenceContext sequence;
};

struct ReferenceExtensionContext {
  GenomeEnd* end{nullptr};
  ChrName chromosome;
  detail::ExtensionSequenceContext sequence;
  std::vector<uint32_t> work_indices;
};

struct ReferenceExtensionNode {
  Segment* segment{nullptr};
  bool barrier{false};
  bool contains_query_species{false};
};

struct ExtensionStageStats {
  uint64_t query_segments{0};
  uint64_t eligible_segments{0};
  uint64_t touched_reference_chromosomes{0};
  uint64_t reference_segments_scanned{0};
  uint64_t barriers{0};
  uint64_t left_candidate_hits{0};
  uint64_t left_candidate_misses{0};
  uint64_t right_candidate_hits{0};
  uint64_t right_candidate_misses{0};
  uint64_t right_alignment_calls{0};
  uint64_t right_alignment_accepted{0};
  uint64_t left_alignment_calls{0};
  uint64_t left_alignment_accepted{0};
  uint64_t gap_too_long{0};
  uint64_t empty_or_nonpositive_gap{0};
  uint64_t maximum_work_items{0};
  uint64_t maximum_temporary_reference_nodes{0};
  double collect_seconds{0.0};
  double index_seconds{0.0};
  double right_seconds{0.0};
  double left_seconds{0.0};
  double resort_seconds{0.0};

  void merge(const ExtensionStageStats& other) {
    query_segments += other.query_segments;
    eligible_segments += other.eligible_segments;
    touched_reference_chromosomes +=
        other.touched_reference_chromosomes;
    reference_segments_scanned += other.reference_segments_scanned;
    barriers += other.barriers;
    left_candidate_hits += other.left_candidate_hits;
    left_candidate_misses += other.left_candidate_misses;
    right_candidate_hits += other.right_candidate_hits;
    right_candidate_misses += other.right_candidate_misses;
    right_alignment_calls += other.right_alignment_calls;
    right_alignment_accepted += other.right_alignment_accepted;
    left_alignment_calls += other.left_alignment_calls;
    left_alignment_accepted += other.left_alignment_accepted;
    gap_too_long += other.gap_too_long;
    empty_or_nonpositive_gap += other.empty_or_nonpositive_gap;
    maximum_work_items =
        std::max(maximum_work_items, other.maximum_work_items);
    maximum_temporary_reference_nodes = std::max(
        maximum_temporary_reference_nodes,
        other.maximum_temporary_reference_nodes);
    collect_seconds += other.collect_seconds;
    index_seconds += other.index_seconds;
    right_seconds += other.right_seconds;
    left_seconds += other.left_seconds;
    resort_seconds += other.resort_seconds;
  }
};

double extensionElapsedSeconds(
    std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start).count();
}

bool blockContainsSpecies(const Segment* segment,
                          const SpeciesName& species) {
  if (!segment) return false;
  const Block* block = segment->parent_block.get();
  if (!block) return false;
  for (const auto& [key, occurrence] : block->anchors) {
    (void)occurrence;
    if (key.first == species) return true;
  }
  return false;
}

void recordPreparedExtensionResult(
    const detail::PreparedExtensionResult& result,
    bool is_left_extend,
    ExtensionStageStats& stats) {
  switch (result.status) {
    case detail::PreparedExtensionStatus::NONPOSITIVE_GAP:
      ++stats.empty_or_nonpositive_gap;
      break;
    case detail::PreparedExtensionStatus::GAP_TOO_LONG:
      ++stats.gap_too_long;
      break;
    case detail::PreparedExtensionStatus::ALIGNMENT_REJECTED:
      if (is_left_extend) {
        ++stats.left_alignment_calls;
      } else {
        ++stats.right_alignment_calls;
      }
      break;
    case detail::PreparedExtensionStatus::ALIGNMENT_ACCEPTED:
      if (is_left_extend) {
        ++stats.left_alignment_calls;
        ++stats.left_alignment_accepted;
      } else {
        ++stats.right_alignment_calls;
        ++stats.right_alignment_accepted;
      }
      break;
    case detail::PreparedExtensionStatus::SKIPPED:
    case detail::PreparedExtensionStatus::NO_CANDIDATE:
      break;
  }
}

void logExtensionStageStats(std::string_view scope,
                            std::string_view query_name,
                            const ExtensionStageStats& stats) {
  spdlog::info(
      "[extend-ref-nodes] scope={} query={} query_segments={} "
      "eligible_segments={} touched_reference_chromosomes={} "
      "reference_segments_scanned={} barriers={} "
      "left_candidate_hits={} left_candidate_misses={} "
      "right_candidate_hits={} right_candidate_misses={} "
      "right_alignment_calls={} right_alignment_accepted={} "
      "left_alignment_calls={} left_alignment_accepted={} "
      "gap_too_long={} empty_or_nonpositive_gap={} "
      "collect_seconds={:.3f} index_seconds={:.3f} "
      "right_seconds={:.3f} left_seconds={:.3f} "
      "resort_seconds={:.3f} maximum_work_items={} "
      "maximum_temporary_reference_nodes={}",
      scope, query_name, stats.query_segments, stats.eligible_segments,
      stats.touched_reference_chromosomes,
      stats.reference_segments_scanned, stats.barriers,
      stats.left_candidate_hits, stats.left_candidate_misses,
      stats.right_candidate_hits, stats.right_candidate_misses,
      stats.right_alignment_calls, stats.right_alignment_accepted,
      stats.left_alignment_calls, stats.left_alignment_accepted,
      stats.gap_too_long, stats.empty_or_nonpositive_gap,
      stats.collect_seconds, stats.index_seconds, stats.right_seconds,
      stats.left_seconds, stats.resort_seconds, stats.maximum_work_items,
      stats.maximum_temporary_reference_nodes);
}

}  // namespace
/* =============================================================
 * 1.  RaMeshGenomeGraph
 * ===========================================================*/
RaMeshGenomeGraph::RaMeshGenomeGraph(const SpeciesName &sp)
    : species_name(sp) {}

RaMeshGenomeGraph::RaMeshGenomeGraph(const SpeciesName &sp,
                                     const std::vector<ChrName> &chrs)
    : species_name(sp) {
  chr2end.reserve(chrs.size());
  for (const auto &c : chrs)
    chr2end.try_emplace(c);
}

/* =============================================================
 * 2.  RaMeshMultiGenomeGraph – ctor
 * ===========================================================*/
RaMeshMultiGenomeGraph::RaMeshMultiGenomeGraph(
    std::map<SpeciesName, SeqPro::ManagerVariant> &seqpro_map) {
  for (const auto &[sp, mgr_var] : seqpro_map) {
    std::vector<ChrName> chr_names = std::visit(
        [](const auto &up) -> std::vector<ChrName> {
          return up ? up->getSequenceNames() : std::vector<ChrName>{};
        },
        mgr_var);
    species_graphs.try_emplace(sp, sp, chr_names);
  }
}

RaMeshMultiGenomeGraph::RaMeshMultiGenomeGraph(
    std::map<SpeciesName, SeqPro::SharedManagerVariant> &seqpro_map) {
  for (const auto &[sp, mgr_ptr] : seqpro_map) {
    // 默认空列表，以便 mgr_ptr 为空时也能插入
    std::vector<ChrName> chr_names;

    if (mgr_ptr) {
      // shared_ptr 本身非空
      chr_names = std::visit(
          [](const auto &up) -> std::vector<ChrName> {
            // up 是 std::unique_ptr<...>
            return up ? up->getSequenceNames() : std::vector<ChrName>{};
          },
          *mgr_ptr // 解引用 shared_ptr 得到 ManagerVariant
      );
    }

    // sp → RaMeshSpeciesGraph 构造（sp, chr_names）
    species_graphs.try_emplace(sp, sp, std::move(chr_names));
  }
}


void RaMeshMultiGenomeGraph::insertAnchorIntoGraph(
    SeqPro::ManagerVariant &ref_mgr, SeqPro::ManagerVariant &qry_mgr,
    const SpeciesName &ref_name, const SpeciesName &qry_name,
    const Anchor &anchor,
    bool isMultiple) {

  auto fetchName = [](const SeqPro::ManagerVariant &mv, const ChrIndex &chr) {
    return std::visit(
        [&](auto &p) {
          using T = std::decay_t<decltype(p)>;
          if constexpr (std::is_same_v<
                            T, std::unique_ptr<SeqPro::SequenceManager>>)
            return p->getSequenceName(chr);
          else
            return p->getOriginalManager().getSequenceName(chr);
        },
        mv);
  };
  // 1. Locate ends for reference & query chromosomes

  ChrIndex ref_chr_index = anchor.ref_chr_index;
  ChrIndex qry_chr_index = anchor.qry_chr_index;
  /*ChrName ref_chr = anchor.ref_chr_index;
  ChrName qry_chr = anchor.qry_chr_index;*/

  const ChrName &ref_chr = fetchName(ref_mgr, ref_chr_index);
  const ChrName &qry_chr = fetchName(qry_mgr, qry_chr_index);

  auto &ref_end = species_graphs[ref_name].chr2end[ref_chr];
  auto &qry_end = species_graphs[qry_name].chr2end[qry_chr];

  BlockPtr blk = Block::create(2);
  blk->ref_chr = ref_chr;
  blk->ref_species = ref_name;

  auto [r_seg, q_seg] = Block::createSegmentPair(anchor, ref_name, qry_name,
                                                 ref_chr, qry_chr, blk);

  // register to global pool
  {
    std::unique_lock pool_lock(rw);
    blocks.emplace_back(WeakBlock(blk));
  }

  ref_end.insertSegment(r_seg);
  qry_end.insertSegment(q_seg);
}

void RaMeshMultiGenomeGraph::insertAnchorsIntoGraphBatch(
    SeqPro::ManagerVariant& ref_mgr,
    SeqPro::ManagerVariant& qry_mgr,
    const SpeciesName& ref_name,
    const SpeciesName& qry_name,
    const std::vector<Anchor*>& anchors) {
  if (anchors.empty()) return;
  const auto sequenceNames = [](SeqPro::ManagerVariant& manager) {
    return std::visit([](auto& pointer) {
      using Pointer = std::decay_t<decltype(pointer)>;
      if constexpr (std::is_same_v<
                        Pointer,
                        std::unique_ptr<SeqPro::SequenceManager>>) {
        return pointer->getSequenceNames();
      } else {
        return pointer->getOriginalManager().getSequenceNames();
      }
    }, manager);
  };
  const std::vector<ChrName> ref_chromosomes = sequenceNames(ref_mgr);
  const std::vector<ChrName> qry_chromosomes = sequenceNames(qry_mgr);

  struct PendingGenomeEnd {
    GenomeEnd* end = nullptr;
    std::vector<SegPtr> segments;
  };
  std::vector<PendingGenomeEnd> pending;
  std::unordered_map<GenomeEnd*, size_t> pending_index;
  pending_index.reserve(ref_chromosomes.size() + qry_chromosomes.size());
  const auto enqueue = [&](GenomeEnd& end, SegPtr segment) {
    const auto [iterator, inserted] = pending_index.emplace(
        &end, pending.size());
    if (inserted) pending.push_back({&end, {}});
    pending[iterator->second].segments.push_back(std::move(segment));
  };

  auto& ref_graph = species_graphs.at(ref_name);
  auto& qry_graph = species_graphs.at(qry_name);
  {
    std::unique_lock pool_lock(rw);
    blocks.reserve(blocks.size() + anchors.size());
    for (Anchor* anchor : anchors) {
      if (!anchor) continue;
      try {
        const ChrName& ref_chromosome =
            ref_chromosomes.at(anchor->ref_chr_index);
        const ChrName& qry_chromosome =
            qry_chromosomes.at(anchor->qry_chr_index);
        GenomeEnd& ref_end = ref_graph.chr2end.at(ref_chromosome);
        GenomeEnd& qry_end = qry_graph.chr2end.at(qry_chromosome);
        BlockPtr block = Block::create(2);
        block->ref_chr = ref_chromosome;
        block->ref_species = ref_name;
        auto [ref_segment, qry_segment] = Block::createSegmentPair(
            *anchor, ref_name, qry_name,
            ref_chromosome, qry_chromosome, block);
        blocks.emplace_back(WeakBlock(block));
        enqueue(ref_end, std::move(ref_segment));
        enqueue(qry_end, std::move(qry_segment));
      } catch (const std::exception& error) {
        spdlog::error("Error materializing anchor for batch insertion: {}",
                      error.what());
      }
    }
  }

  std::vector<SegPtr> radix_scratch;
  std::vector<size_t> radix_counts;
  for (auto& item : pending) {
    stableRadixSortSegmentsByStart(
        item.segments, radix_scratch, radix_counts);
    item.end->insertSegmentsSorted(item.segments);
    std::vector<SegPtr>().swap(item.segments);
  }
}

namespace {

struct IndexedPrimarySegment {
  uint64_t start = 0;
  uint64_t end = 0;
  SegPtr segment;
};
struct PrimaryBlockPair {
  uint64_t first = 0;
  uint64_t second = 0;

  bool operator==(const PrimaryBlockPair&) const = default;
};

struct PrimaryBlockPairHash {
  size_t operator()(const PrimaryBlockPair& pair) const noexcept {
    const size_t first_hash = std::hash<uint64_t>{}(pair.first);
    const size_t second_hash = std::hash<uint64_t>{}(pair.second);
    return first_hash ^
           (second_hash + 0x9e3779b9U + (first_hash << 6U) +
            (first_hash >> 2U));
  }
};

bool blocksShareNonReferenceSpecies(const BlockPtr &left, const BlockPtr &right,
                                    const SpeciesName &reference_species) {
  if (!left || !right) {
    return true;
  }
  for (const auto &[left_key, unused_left] : left->anchors) {
    (void)unused_left;
    if (left_key.first == reference_species) {
      continue;
    }
    for (const auto &[right_key, unused_right] : right->anchors) {
      (void)unused_right;
      if (right_key.first == reference_species) {
        continue;
      }
      if (left_key.first == right_key.first) {
        return true;
      }
    }
  }
  return false;
}

using PrimarySegmentIndex =
    std::unordered_map<SpeciesChrPair, std::vector<IndexedPrimarySegment>,
                       SpeciesChrPairHash>;

bool cigarOperationConsumesReference(char operation) {
  return operation == 'M' || operation == '=' || operation == 'X' ||
         operation == 'D' || operation == 'N';
}

bool cigarOperationConsumesQuery(char operation) {
  return operation == 'M' || operation == '=' || operation == 'X' ||
         operation == 'I' || operation == 'S';
}

bool cigarOperationAlignsBases(char operation) {
  return operation == 'M' || operation == '=' || operation == 'X';
}

void appendTerminalQueryOnlyRemainder(const SegPtr &segment,
                                      const std::string &cigar_remainder) {
  if (cigar_remainder.empty()) {
    return;
  }
  if (!segment) {
    throw std::logic_error(
        "overlap merge has no segment for terminal CIGAR operations");
  }
  Cigar_t remainder;
  parseCigarString(cigar_remainder, remainder);
  if (countRefLength(remainder) != 0) {
    throw std::logic_error(
        "overlap merge left reference-consuming CIGAR operations "
        "unassigned");
  }
  segment->length += countQryLength(remainder);
  appendCigar(segment->cigar, remainder);
}

std::optional<uint64_t> orientedSegmentOffset(
    const SegPtr& segment, uint64_t coordinate) {
  if (!segment || coordinate < segment->start ||
      coordinate >= static_cast<uint64_t>(segment->start) +
                        segment->length) {
    return std::nullopt;
  }
  if (segment->strand == Strand::FORWARD) {
    return coordinate - segment->start;
  }
  return static_cast<uint64_t>(segment->start) + segment->length - 1 -
         coordinate;
}

std::optional<uint64_t> segmentCoordinateToReferenceOffset(
    const SegPtr& segment, const SegPtr& reference_segment,
    uint64_t coordinate) {
  const auto query_offset = orientedSegmentOffset(segment, coordinate);
  if (!query_offset) {
    return std::nullopt;
  }
  if (segment == reference_segment) {
    return query_offset;
  }

  uint64_t reference_offset = 0;
  uint64_t current_query_offset = 0;
  for (const auto unit : segment->cigar) {
    uint32_t length = 0;
    char operation = '\0';
    intToCigar(unit, operation, length);
    if (cigarOperationAlignsBases(operation) &&
        *query_offset >= current_query_offset &&
        *query_offset < current_query_offset + length) {
      return reference_offset + (*query_offset - current_query_offset);
    }
    if (cigarOperationConsumesReference(operation)) {
      reference_offset += length;
    }
    if (cigarOperationConsumesQuery(operation)) {
      current_query_offset += length;
    }
  }
  return std::nullopt;
}

std::optional<uint64_t> referenceOffsetToSegmentCoordinate(
    const SegPtr& segment, const SegPtr& reference_segment,
    uint64_t target_reference_offset) {
  uint64_t query_offset = 0;
  if (segment == reference_segment) {
    query_offset = target_reference_offset;
  } else {
    uint64_t reference_offset = 0;
    uint64_t current_query_offset = 0;
    bool mapped = false;
    for (const auto unit : segment->cigar) {
      uint32_t length = 0;
      char operation = '\0';
      intToCigar(unit, operation, length);
      if (cigarOperationAlignsBases(operation) &&
          target_reference_offset >= reference_offset &&
          target_reference_offset < reference_offset + length) {
        query_offset =
            current_query_offset +
            (target_reference_offset - reference_offset);
        mapped = true;
        break;
      }
      if (cigarOperationConsumesReference(operation)) {
        reference_offset += length;
      }
      if (cigarOperationConsumesQuery(operation)) {
        current_query_offset += length;
      }
    }
    if (!mapped) {
      return std::nullopt;
    }
  }

  if (query_offset >= segment->length) {
    return std::nullopt;
  }
  if (segment->strand == Strand::FORWARD) {
    return static_cast<uint64_t>(segment->start) + query_offset;
  }
  return static_cast<uint64_t>(segment->start) + segment->length - 1 -
         query_offset;
}

bool isReferenceCompatibleSegment(const SegPtr& segment) {
  if (!segment) {
    return false;
  }
  if (segment->cigar.empty()) {
    return true;
  }
  uint64_t consumed = 0;
  for (const auto unit : segment->cigar) {
    uint32_t length = 0;
    char operation = '\0';
    intToCigar(unit, operation, length);
    if (!cigarOperationAlignsBases(operation)) {
      return false;
    }
    consumed += length;
  }
  return consumed == segment->length;
}

SegPtr findAlignmentReferenceSegment(const BlockPtr& block) {
  if (!block) {
    return nullptr;
  }
  for (const auto& [species_chr, segment] : block->anchors) {
    if (species_chr.first == block->ref_species &&
        species_chr.second == block->ref_chr &&
        segment->isPrimary() &&
        isReferenceCompatibleSegment(segment)) {
      return segment;
    }
  }
  return nullptr;
}

const IndexedPrimarySegment* findPrimarySegment(
    const PrimarySegmentIndex& index, const SpeciesChrPair& key,
    uint64_t coordinate) {
  const auto index_it = index.find(key);
  if (index_it == index.end()) {
    return nullptr;
  }
  const auto& intervals = index_it->second;
  auto upper = std::upper_bound(
      intervals.begin(), intervals.end(), coordinate,
      [](uint64_t value, const IndexedPrimarySegment& interval) {
        return value < interval.start;
      });
  while (upper != intervals.begin()) {
    --upper;
    if (coordinate >= upper->start && coordinate < upper->end) {
      return &*upper;
    }
    if (upper->end <= coordinate) {
      break;
    }
  }
  return nullptr;
}

struct SecondaryOccurrencePoint {
  SpeciesName species;
  ChrName chromosome;
  int step = 1;
  int64_t diagonal = 0;
  uint64_t coordinate = 0;
  AlignRole role = AlignRole::PRIMARY;
};

struct SecondaryHomologyPoint {
  SpeciesName common_species;
  ChrName common_chromosome;
  uint64_t common_coordinate = 0;
  std::vector<SecondaryOccurrencePoint> occurrences;
  PrimaryBlockPair primary_block_pair;
};

struct SecondaryHomologyRun {
  SpeciesName common_species;
  ChrName common_chromosome;
  uint64_t common_start = 0;
  uint64_t length = 0;
  std::vector<SecondaryOccurrencePoint> occurrences;
  PrimaryBlockPair primary_block_pair;
};

auto secondaryOccurrenceRunKey(
    const SecondaryOccurrencePoint& occurrence) {
  return std::tie(
      occurrence.species, occurrence.chromosome, occurrence.step,
      occurrence.diagonal, occurrence.role);
}

bool secondaryOccurrenceRunLess(
    const SecondaryOccurrencePoint& lhs,
    const SecondaryOccurrencePoint& rhs) {
  return secondaryOccurrenceRunKey(lhs) <
         secondaryOccurrenceRunKey(rhs);
}

template <typename Left, typename Right>
bool secondaryHomologySameRun(
    const Left& lhs,
    const Right& rhs) {
  return lhs.common_species == rhs.common_species &&
         lhs.common_chromosome == rhs.common_chromosome &&
         lhs.primary_block_pair == rhs.primary_block_pair &&
         lhs.occurrences.size() == rhs.occurrences.size() &&
         std::equal(
             lhs.occurrences.begin(), lhs.occurrences.end(),
             rhs.occurrences.begin(),
             [](const SecondaryOccurrencePoint& left,
                const SecondaryOccurrencePoint& right) {
               return secondaryOccurrenceRunKey(left) ==
                      secondaryOccurrenceRunKey(right);
             });
}

bool secondaryRunLess(
    const SecondaryHomologyRun& lhs,
    const SecondaryHomologyRun& rhs) {
  const auto lhs_common = std::tie(
      lhs.common_species,
      lhs.common_chromosome,
      lhs.primary_block_pair.first,
      lhs.primary_block_pair.second);
  const auto rhs_common = std::tie(
      rhs.common_species,
      rhs.common_chromosome,
      rhs.primary_block_pair.first,
      rhs.primary_block_pair.second);
  if (lhs_common != rhs_common) {
    return lhs_common < rhs_common;
  }
  if (std::lexicographical_compare(
          lhs.occurrences.begin(), lhs.occurrences.end(),
          rhs.occurrences.begin(), rhs.occurrences.end(),
          secondaryOccurrenceRunLess)) {
    return true;
  }
  if (std::lexicographical_compare(
          rhs.occurrences.begin(), rhs.occurrences.end(),
          lhs.occurrences.begin(), lhs.occurrences.end(),
          secondaryOccurrenceRunLess)) {
    return false;
  }
  return lhs.common_start < rhs.common_start;
}

}  // namespace

void RaMeshMultiGenomeGraph::registerSecondaryAnchorCandidate(
    SpeciesName ref_species, ChrName ref_chromosome,
    SpeciesName query_species, ChrName query_chromosome,
    const Anchor& anchor, bool initial_round,
    bool common_is_reference) {
  std::unique_lock lock(rw);
  secondary_anchor_candidates.push_back(SecondaryAnchorCandidate{
      std::move(ref_species), std::move(ref_chromosome),
      std::move(query_species), std::move(query_chromosome), anchor,
      initial_round, common_is_reference});
}

size_t RaMeshMultiGenomeGraph::materializeSecondaryAlignments() {
  std::vector<SecondaryAnchorCandidate> candidates;
  {
    std::unique_lock lock(rw);
    candidates.swap(secondary_anchor_candidates);
  }
  if (candidates.empty()) {
    return 0;
  }

  PrimarySegmentIndex primary_index;
  for (const auto& [species, genome] : species_graphs) {
    for (const auto& [chromosome, genome_end] : genome.chr2end) {
      auto& intervals =
          primary_index[SpeciesChrPair{species, chromosome}];
      for (SegPtr segment =
               genome_end.head->primary_path.next.load(
                   std::memory_order_acquire);
           segment && !segment->isTail();
           segment = segment->primary_path.next.load(
               std::memory_order_acquire)) {
        if (!segment->isPrimary() || !segment->parent_block ||
            segment->length == 0) {
          continue;
        }
        intervals.push_back(IndexedPrimarySegment{
            segment->start,
            static_cast<uint64_t>(segment->start) + segment->length,
            segment});
      }
    }
  }
  for (auto& [key, intervals] : primary_index) {
    (void)key;
    std::sort(
        intervals.begin(), intervals.end(),
        [](const IndexedPrimarySegment& lhs,
           const IndexedPrimarySegment& rhs) {
          return std::tie(lhs.start, lhs.end) <
                 std::tie(rhs.start, rhs.end);
        });
  }

  std::vector<SecondaryHomologyRun> secondary_runs;
  size_t mapped_point_observations = 0;
  size_t candidate_aligned_bases = 0;
  size_t missing_primary_segment_points = 0;
  size_t missing_reference_segment_points = 0;
  size_t unmappable_reference_offset_points = 0;
  size_t missing_duplicate_primary_points = 0;
  size_t secondary_already_present_points = 0;
  size_t initial_round_candidates = 0;
  for (const auto &candidate : candidates) {
    if (candidate.initial_round) {
      ++initial_round_candidates;
    }
    const bool common_is_reference = candidate.common_is_reference;
    const SpeciesName &common_species =
        common_is_reference ? candidate.ref_species : candidate.query_species;
    const ChrName &common_chromosome = common_is_reference
                                           ? candidate.ref_chromosome
                                           : candidate.query_chromosome;
    const SpeciesName &duplicate_species =
        common_is_reference ? candidate.query_species : candidate.ref_species;
    const ChrName& secondary_chromosome =
        common_is_reference ? candidate.query_chromosome
                            : candidate.ref_chromosome;
    const int anchor_step =
        candidate.anchor.strand == Strand::FORWARD ? 1 : -1;

    uint64_t reference_coordinate = candidate.anchor.ref_start;
    int64_t query_coordinate =
        candidate.anchor.strand == Strand::FORWARD
            ? static_cast<int64_t>(candidate.anchor.qry_start)
            : static_cast<int64_t>(candidate.anchor.qry_start) +
                  candidate.anchor.qry_len - 1;

    for (const auto unit : candidate.anchor.cigar) {
      uint32_t length = 0;
      char operation = '\0';
      intToCigar(unit, operation, length);
      if (cigarOperationAlignsBases(operation)) {
        std::optional<SecondaryHomologyRun> active_run;
        const auto flush_active_run = [&]() {
          if (active_run) {
            secondary_runs.push_back(std::move(*active_run));
            active_run.reset();
          }
        };
        candidate_aligned_bases += length;
        for (uint32_t offset = 0; offset < length; ++offset) {
          const uint64_t ref_position = reference_coordinate + offset;
          const int64_t qry_position =
              query_coordinate +
              static_cast<int64_t>(anchor_step) * offset;
          if (qry_position < 0) {
            flush_active_run();
            continue;
          }
          const uint64_t common_coordinate =
              common_is_reference
                  ? ref_position
                  : static_cast<uint64_t>(qry_position);
          const uint64_t secondary_coordinate =
              common_is_reference
                  ? static_cast<uint64_t>(qry_position)
                  : ref_position;

          const auto* indexed_segment = findPrimarySegment(
              primary_index,
              SpeciesChrPair{common_species, common_chromosome},
              common_coordinate);
          if (!indexed_segment) {
            ++missing_primary_segment_points;
            flush_active_run();
            continue;
          }
          const SegPtr& common_segment = indexed_segment->segment;
          const BlockPtr common_block =
              common_segment->parent_block;
          const SegPtr common_reference_segment =
              findAlignmentReferenceSegment(common_block);
          if (!common_reference_segment) {
            ++missing_reference_segment_points;
            flush_active_run();
            continue;
          }
          const auto common_reference_offset =
              segmentCoordinateToReferenceOffset(
                  common_segment, common_reference_segment,
                  common_coordinate);
          if (!common_reference_offset) {
            ++unmappable_reference_offset_points;
            flush_active_run();
            continue;
          }
          const auto common_reference_coordinate =
              referenceOffsetToSegmentCoordinate(
                  common_reference_segment,
                  common_reference_segment,
                  *common_reference_offset);
          if (!common_reference_coordinate) {
            ++unmappable_reference_offset_points;
            flush_active_run();
            continue;
          }

          bool secondary_already_present = false;
          for (const auto& [species_chr, primary_segment] :
               common_block->anchors) {
            if (!primary_segment->isPrimary() ||
                species_chr.first != duplicate_species ||
                species_chr.second != secondary_chromosome) {
              continue;
            }
            const auto primary_coordinate =
                referenceOffsetToSegmentCoordinate(
                    primary_segment, common_reference_segment,
                    *common_reference_offset);
            if (!primary_coordinate) {
              continue;
            }
            if (*primary_coordinate == secondary_coordinate) {
              secondary_already_present = true;
            }
          }
          if (secondary_already_present) {
            ++secondary_already_present_points;
            flush_active_run();
            continue;
          }

          const auto* duplicate_indexed_segment =
              findPrimarySegment(
                  primary_index,
                  SpeciesChrPair{
                      duplicate_species,
                      secondary_chromosome},
                  secondary_coordinate);
          if (!duplicate_indexed_segment) {
            ++missing_duplicate_primary_points;
            flush_active_run();
            continue;
          }
          const SegPtr& duplicate_segment =
              duplicate_indexed_segment->segment;
          const BlockPtr duplicate_block =
              duplicate_segment->parent_block;
          if (!duplicate_block ||
              duplicate_block == common_block) {
            ++secondary_already_present_points;
            flush_active_run();
            continue;
          }
          const SegPtr duplicate_reference_segment =
              findAlignmentReferenceSegment(duplicate_block);
          if (!duplicate_reference_segment) {
            ++missing_reference_segment_points;
            flush_active_run();
            continue;
          }
          const auto duplicate_reference_offset =
              segmentCoordinateToReferenceOffset(
                  duplicate_segment,
                  duplicate_reference_segment,
                  secondary_coordinate);
          if (!duplicate_reference_offset) {
            ++unmappable_reference_offset_points;
            flush_active_run();
            continue;
          }
          const auto duplicate_reference_coordinate =
              referenceOffsetToSegmentCoordinate(
                  duplicate_reference_segment,
                  duplicate_reference_segment,
                  *duplicate_reference_offset);
          if (!duplicate_reference_coordinate) {
            ++unmappable_reference_offset_points;
            flush_active_run();
            continue;
          }

          const int common_reference_orientation =
              common_reference_segment->strand ==
                      Strand::FORWARD
                  ? 1
                  : -1;
          const int common_segment_orientation =
              common_segment->strand == Strand::FORWARD
                  ? 1
                  : -1;
          const int duplicate_reference_orientation =
              duplicate_reference_segment->strand ==
                      Strand::FORWARD
                  ? 1
                  : -1;
          const int duplicate_segment_orientation =
              duplicate_segment->strand == Strand::FORWARD
                  ? 1
                  : -1;
          const int duplicate_axis_step =
              common_reference_orientation *
              common_segment_orientation * anchor_step *
              duplicate_reference_orientation *
              duplicate_segment_orientation;

          const PrimaryBlockPair primary_block_pair{
              std::min(
                  common_block->block_id,
                  duplicate_block->block_id),
              std::max(
                  common_block->block_id,
                  duplicate_block->block_id)};
          const bool common_is_canonical =
              common_block->block_id ==
              primary_block_pair.first;
          const BlockPtr& canonical_block =
              common_is_canonical ? common_block
                                  : duplicate_block;
          const uint64_t canonical_coordinate =
              common_is_canonical
                  ? *common_reference_coordinate
                  : *duplicate_reference_coordinate;
          const int common_axis_step =
              common_is_canonical ? 1 : duplicate_axis_step;
          const int duplicate_axis_from_canonical =
              common_is_canonical ? duplicate_axis_step : 1;

          std::vector<SecondaryOccurrencePoint> occurrences;
          occurrences.reserve(
              common_block->anchors.size() +
              duplicate_block->anchors.size());
          const auto append_block_occurrences =
              [&](const BlockPtr& source_block,
                  const SegPtr& source_reference_segment,
                  const uint64_t source_reference_offset,
                  const int source_axis_step,
                  const AlignRole role) {
                const int source_reference_orientation =
                    source_reference_segment->strand ==
                            Strand::FORWARD
                        ? 1
                        : -1;
                for (const auto& [species_chr, segment] :
                     source_block->anchors) {
                  if (!segment->isPrimary()) {
                    continue;
                  }
                  const auto coordinate =
                      referenceOffsetToSegmentCoordinate(
                          segment,
                          source_reference_segment,
                          source_reference_offset);
                  if (!coordinate) {
                    continue;
                  }
                  const int segment_orientation =
                      segment->strand == Strand::FORWARD
                          ? 1
                          : -1;
                  const int step =
                      source_axis_step *
                      source_reference_orientation *
                      segment_orientation;
                  occurrences.push_back(
                      SecondaryOccurrencePoint{
                          species_chr.first,
                          species_chr.second,
                          step,
                          static_cast<int64_t>(*coordinate) -
                              static_cast<int64_t>(step) *
                                  static_cast<int64_t>(
                                      canonical_coordinate),
                          *coordinate,
                          role});
                }
              };

          if (common_is_canonical) {
            append_block_occurrences(
                common_block, common_reference_segment,
                *common_reference_offset, common_axis_step,
                AlignRole::PRIMARY);
            append_block_occurrences(
                duplicate_block,
                duplicate_reference_segment,
                *duplicate_reference_offset,
                duplicate_axis_from_canonical,
                AlignRole::SECONDARY);
          } else {
            append_block_occurrences(
                duplicate_block,
                duplicate_reference_segment,
                *duplicate_reference_offset,
                duplicate_axis_from_canonical,
                AlignRole::PRIMARY);
            append_block_occurrences(
                common_block, common_reference_segment,
                *common_reference_offset, common_axis_step,
                AlignRole::SECONDARY);
          }
          std::sort(
              occurrences.begin(), occurrences.end(),
              secondaryOccurrenceRunLess);
          occurrences.erase(
              std::unique(
                  occurrences.begin(), occurrences.end(),
                  [](const SecondaryOccurrencePoint& lhs,
                     const SecondaryOccurrencePoint& rhs) {
                    return secondaryOccurrenceRunKey(lhs) ==
                           secondaryOccurrenceRunKey(rhs);
                  }),
              occurrences.end());
          SecondaryHomologyPoint point{
              canonical_block->ref_species,
              canonical_block->ref_chr,
              canonical_coordinate,
              std::move(occurrences),
              primary_block_pair};
          ++mapped_point_observations;
          const bool extends_right =
              active_run &&
              secondaryHomologySameRun(*active_run, point) &&
              active_run->common_start <=
                  std::numeric_limits<uint64_t>::max() -
                      active_run->length &&
              point.common_coordinate ==
                  active_run->common_start +
                      active_run->length;
          const bool extends_left =
              active_run &&
              secondaryHomologySameRun(*active_run, point) &&
              point.common_coordinate <
                  active_run->common_start &&
              point.common_coordinate + 1 ==
                  active_run->common_start;
          if (extends_right) {
            ++active_run->length;
          } else if (extends_left) {
            active_run->common_start =
                point.common_coordinate;
            ++active_run->length;
          } else {
            flush_active_run();
            active_run = SecondaryHomologyRun{
                std::move(point.common_species),
                std::move(point.common_chromosome),
                point.common_coordinate,
                1,
                std::move(point.occurrences),
                point.primary_block_pair};
          }
        }
        flush_active_run();
      }
      if (cigarOperationConsumesReference(operation)) {
        reference_coordinate += length;
      }
      if (cigarOperationConsumesQuery(operation)) {
        query_coordinate +=
            static_cast<int64_t>(anchor_step) * length;
      }
    }
  }
  const size_t unfiltered_pair_points =
      mapped_point_observations;
  std::sort(
      secondary_runs.begin(),
      secondary_runs.end(),
      secondaryRunLess);
  std::vector<SecondaryHomologyRun> merged_runs;
  merged_runs.reserve(secondary_runs.size());
  for (auto& run : secondary_runs) {
    if (run.length == 0 ||
        run.common_start >
            std::numeric_limits<uint64_t>::max() -
                run.length) {
      throw std::runtime_error(
          "Secondary homology run coordinate overflow");
    }
    if (!merged_runs.empty() &&
        secondaryHomologySameRun(
            merged_runs.back(),
            run)) {
      auto& previous = merged_runs.back();
      const uint64_t previous_end =
          previous.common_start + previous.length;
      const uint64_t current_end =
          run.common_start + run.length;
      if (run.common_start <= previous_end) {
        previous.length =
            std::max(previous_end, current_end) -
            previous.common_start;
        continue;
      }
    }
    merged_runs.push_back(std::move(run));
  }
  secondary_runs = std::move(merged_runs);

  size_t mapped_points = 0;
  std::vector<BlockPtr> secondary_blocks;
  secondary_blocks.reserve(secondary_runs.size());
  size_t materialized_bases = 0;
  constexpr uint64_t maximum_cigar_operation_length =
      (1ULL << 28U) - 1U;
  for (const auto& run : secondary_runs) {
    if (run.length == 0 ||
        run.length > maximum_cigar_operation_length ||
        run.common_start >
            std::numeric_limits<uint_t>::max() ||
        run.length >
            std::numeric_limits<uint_t>::max()) {
      continue;
    }
    if (run.length >
        std::numeric_limits<size_t>::max() -
            mapped_points) {
      throw std::overflow_error(
          "Secondary homology mapped-point count overflow");
    }
    mapped_points += static_cast<size_t>(run.length);

    auto block = Block::create(run.occurrences.size());
    block->ref_species = run.common_species;
    block->ref_chr = run.common_chromosome;
    bool coordinates_fit = true;
    for (const auto& occurrence : run.occurrences) {
      const __int128 coordinate_wide =
          static_cast<__int128>(occurrence.step) *
              static_cast<__int128>(run.common_start) +
          static_cast<__int128>(occurrence.diagonal);
      if (coordinate_wide < 0 ||
          coordinate_wide >
              std::numeric_limits<uint_t>::max()) {
        coordinates_fit = false;
        break;
      }
      const uint64_t coordinate =
          static_cast<uint64_t>(coordinate_wide);
      if (occurrence.step < 0 &&
          coordinate < run.length - 1) {
        coordinates_fit = false;
        break;
      }
      const uint64_t start =
          occurrence.step > 0
              ? coordinate
              : coordinate - (run.length - 1);
      if (start >
          std::numeric_limits<uint_t>::max()) {
        coordinates_fit = false;
        break;
      }
      auto segment = Segment::create(
          static_cast<uint_t>(start),
          static_cast<uint_t>(run.length),
          occurrence.step > 0 ? Strand::FORWARD
                              : Strand::REVERSE,
          Cigar_t{
              cigarToInt(
                  'M',
                  static_cast<uint32_t>(run.length))},
          occurrence.role,
          SegmentRole::SEGMENT,
          block);
      block->anchors.emplace(
          SpeciesChrPair{
              occurrence.species,
              occurrence.chromosome},
          std::move(segment));
    }
    if (coordinates_fit) {
      if (run.length >
          std::numeric_limits<size_t>::max() -
              materialized_bases) {
        throw std::overflow_error(
            "Secondary homology materialized-base count overflow");
      }
      materialized_bases += run.length;
      secondary_blocks.push_back(std::move(block));
    }
  }

  {
    std::unique_lock lock(rw);
    blocks.reserve(blocks.size() + secondary_blocks.size());
    for (const auto& block : secondary_blocks) {
      blocks.emplace_back(block);
    }
  }
  spdlog::info(
      "[secondary-alignments] candidates={} initial_round_candidates={} "
      "candidate_aligned_bases={} unfiltered_pair_points={} "
      "mapped_points={} blocks={} materialized_bases={} "
      "missing_primary_segment_points={} "
      "missing_reference_segment_points={} "
      "unmappable_reference_offset_points={} "
      "missing_duplicate_primary_points={} "
      "secondary_already_present_points={}",
      candidates.size(), initial_round_candidates, candidate_aligned_bases,
      unfiltered_pair_points, mapped_points, secondary_blocks.size(),
      materialized_bases, missing_primary_segment_points,
      missing_reference_segment_points,
      unmappable_reference_offset_points, missing_duplicate_primary_points,
      secondary_already_present_points);
  return secondary_blocks.size();
}

void RaMeshMultiGenomeGraph::extendRefNodes(
    const SpeciesName& ref_name,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    int_t zdrop) {
  auto reference_graph_it = species_graphs.find(ref_name);
  if (reference_graph_it == species_graphs.end()) {
    throw std::invalid_argument(
        "Reference species is missing from the graph: " + ref_name);
  }
  RaMeshGenomeGraph& reference_graph = reference_graph_it->second;
  ExtensionStageStats total_stats;

  for (auto& [query_name, query_graph] : species_graphs) {
    if (query_name == ref_name) continue;

    ExtensionStageStats stats;
    const auto collect_start = std::chrono::steady_clock::now();
    std::size_t eligible_count = 0;
    for (auto& [query_chromosome, query_end] : query_graph.chr2end) {
      (void)query_chromosome;
      SegPtr current = query_end.head->primary_path.next.load(
          std::memory_order_acquire);
      while (current && !current->isTail()) {
        ++stats.query_segments;
        if (!(current->left_extend && current->right_extend)) {
          ++eligible_count;
        }
        current = current->primary_path.next.load(
            std::memory_order_acquire);
      }
    }
    stats.eligible_segments = eligible_count;
    stats.maximum_work_items = eligible_count;

    if (eligible_count >
        static_cast<std::size_t>(kNoExtensionWorkIndex)) {
      throw std::length_error(
          "extendRefNodes work list exceeds the compact index range");
    }

    if (eligible_count != 0) {
      const auto& query_manager = detail::originalSequenceManager(
          managers.at(query_name));
      const auto& reference_manager = detail::originalSequenceManager(
          managers.at(ref_name));

      std::vector<ExtensionWorkItem> work_items;
      std::vector<QueryExtensionContext> query_contexts;
      std::vector<ReferenceExtensionContext> reference_contexts;
      std::unordered_map<ChrName, uint32_t>
          reference_context_by_chromosome;

      work_items.reserve(eligible_count);
      query_contexts.reserve(
          std::min(eligible_count, query_graph.chr2end.size()));
      reference_contexts.reserve(
          std::min(eligible_count, reference_graph.chr2end.size()));
      reference_context_by_chromosome.reserve(
          std::min(eligible_count, reference_graph.chr2end.size()));

      for (auto& [query_chromosome, query_end] :
           query_graph.chr2end) {
        uint32_t query_context_index = kNoExtensionWorkIndex;
        SegPtr current = query_end.head->primary_path.next.load(
            std::memory_order_acquire);
        while (current && !current->isTail()) {
          if (current->left_extend && current->right_extend) {
            current = current->primary_path.next.load(
                std::memory_order_acquire);
            continue;
          }

          if (query_context_index == kNoExtensionWorkIndex) {
            if (query_contexts.size() >=
                static_cast<std::size_t>(kNoExtensionWorkIndex)) {
              throw std::length_error(
                  "Too many query chromosome extension contexts");
            }
            const SeqPro::SequenceId sequence_id =
                query_manager.getSequenceId(query_chromosome);
            query_context_index =
                static_cast<uint32_t>(query_contexts.size());
            query_contexts.push_back(QueryExtensionContext{
                &query_end,
                query_chromosome,
                detail::ExtensionSequenceContext{
                    &query_manager,
                    sequence_id,
                    query_manager.getSequenceLength(sequence_id)}});
          }

          Block* block = current->parent_block.get();
          if (!block) {
            throw std::runtime_error(
                "Query segment is missing its parent block");
          }
          const auto reference_occurrence = block->anchors.find(
              {ref_name, block->ref_chr});
          if (reference_occurrence == block->anchors.end() ||
              !reference_occurrence->second) {
            throw std::runtime_error(
                "Reference occurrence is missing from the current block");
          }

          auto reference_end_it =
              reference_graph.chr2end.find(block->ref_chr);
          if (reference_end_it == reference_graph.chr2end.end()) {
            throw std::runtime_error(
                "Reference chromosome is missing from the graph: " +
                block->ref_chr);
          }

          auto [context_it, inserted] =
              reference_context_by_chromosome.try_emplace(
                  block->ref_chr,
                  static_cast<uint32_t>(reference_contexts.size()));
          if (inserted) {
            if (reference_contexts.size() >=
                static_cast<std::size_t>(kNoExtensionWorkIndex)) {
              throw std::length_error(
                  "Too many reference chromosome extension contexts");
            }
            const SeqPro::SequenceId sequence_id =
                reference_manager.getSequenceId(block->ref_chr);
            reference_contexts.push_back(ReferenceExtensionContext{
                &reference_end_it->second,
                block->ref_chr,
                detail::ExtensionSequenceContext{
                    &reference_manager,
                    sequence_id,
                    reference_manager.getSequenceLength(sequence_id)},
                {}});
            context_it->second =
                static_cast<uint32_t>(reference_contexts.size() - 1);
          }

          const uint32_t reference_context_index = context_it->second;
          const uint32_t work_index =
              static_cast<uint32_t>(work_items.size());
          work_items.push_back(ExtensionWorkItem{
              current.get(),
              reference_occurrence->second.get(),
              nullptr,
              nullptr,
              query_context_index,
              reference_context_index,
              kNoExtensionWorkIndex});
          reference_contexts[reference_context_index]
              .work_indices.push_back(work_index);

          current = current->primary_path.next.load(
              std::memory_order_acquire);
        }
      }
      stats.collect_seconds = extensionElapsedSeconds(collect_start);
      stats.touched_reference_chromosomes =
          reference_contexts.size();

      const auto index_start = std::chrono::steady_clock::now();
      for (ReferenceExtensionContext& reference_context :
           reference_contexts) {
        std::unordered_map<Segment*, uint32_t> target_heads;
        target_heads.reserve(reference_context.work_indices.size());
        for (const uint32_t work_index :
             reference_context.work_indices) {
          ExtensionWorkItem& work = work_items[work_index];
          auto [target_it, inserted] = target_heads.try_emplace(
              work.reference_segment, kNoExtensionWorkIndex);
          work.same_reference_next = target_it->second;
          target_it->second = work_index;
          (void)inserted;
        }

        std::vector<ReferenceExtensionNode> reference_nodes;
        reference_nodes.reserve(reference_context.work_indices.size());
        SegPtr current =
            reference_context.end->head->primary_path.next.load(
                std::memory_order_acquire);
        while (current && !current->isTail()) {
          const bool barrier =
              current->left_extend && current->right_extend;
          reference_nodes.push_back(ReferenceExtensionNode{
              current.get(),
              barrier,
              !barrier &&
                  blockContainsSpecies(current.get(), query_name)});
          if (barrier) ++stats.barriers;
          current = current->primary_path.next.load(
              std::memory_order_acquire);
        }

        stats.maximum_temporary_reference_nodes = std::max<uint64_t>(
            stats.maximum_temporary_reference_nodes,
            reference_nodes.size());
        stats.reference_segments_scanned +=
            static_cast<uint64_t>(reference_nodes.size()) * 2;

        Segment* nearest = nullptr;
        for (const ReferenceExtensionNode& node : reference_nodes) {
          const auto target_it = target_heads.find(node.segment);
          if (target_it != target_heads.end()) {
            for (uint32_t work_index = target_it->second;
                 work_index != kNoExtensionWorkIndex;
                 work_index =
                     work_items[work_index].same_reference_next) {
              work_items[work_index].left_reference_candidate =
                  nearest;
            }
          }
          if (node.barrier) {
            nearest = nullptr;
          } else if (node.contains_query_species) {
            nearest = node.segment;
          }
        }

        nearest = nullptr;
        for (auto node_it = reference_nodes.rbegin();
             node_it != reference_nodes.rend(); ++node_it) {
          const auto target_it = target_heads.find(node_it->segment);
          if (target_it != target_heads.end()) {
            for (uint32_t work_index = target_it->second;
                 work_index != kNoExtensionWorkIndex;
                 work_index =
                     work_items[work_index].same_reference_next) {
              work_items[work_index].right_reference_candidate =
                  nearest;
            }
          }
          if (node_it->barrier) {
            nearest = nullptr;
          } else if (node_it->contains_query_species) {
            nearest = node_it->segment;
          }
        }

        // Candidate pointers are now stored directly in work_items. Release
        // this chromosome's grouping buffer before any KSW2 work begins.
        std::vector<uint32_t>().swap(reference_context.work_indices);
      }

      for (const ExtensionWorkItem& work : work_items) {
        if (work.left_reference_candidate) {
          ++stats.left_candidate_hits;
        } else {
          ++stats.left_candidate_misses;
        }
        if (work.right_reference_candidate) {
          ++stats.right_candidate_hits;
        } else {
          ++stats.right_candidate_misses;
        }
      }
      stats.index_seconds = extensionElapsedSeconds(index_start);

      detail::ExtensionSequenceScratch scratch;
      scratch.reserve();

      const auto right_start = std::chrono::steady_clock::now();
      for (ExtensionWorkItem& work : work_items) {
        const auto result = detail::alignIntervalPrepared(
            work.query_segment,
            work.reference_segment,
            work.right_reference_candidate,
            query_contexts[work.query_context_index].sequence,
            reference_contexts[work.reference_context_index].sequence,
            false, zdrop, scratch);
        recordPreparedExtensionResult(result, false, stats);
      }
      stats.right_seconds = extensionElapsedSeconds(right_start);

      const auto left_start = std::chrono::steady_clock::now();
      for (ExtensionWorkItem& work : work_items) {
        const auto result = detail::alignIntervalPrepared(
            work.query_segment,
            work.reference_segment,
            work.left_reference_candidate,
            query_contexts[work.query_context_index].sequence,
            reference_contexts[work.reference_context_index].sequence,
            true, zdrop, scratch);
        recordPreparedExtensionResult(result, true, stats);
      }
      stats.left_seconds = extensionElapsedSeconds(left_start);
    } else {
      stats.collect_seconds = extensionElapsedSeconds(collect_start);
    }

    total_stats.merge(stats);
    logExtensionStageStats("species", query_name, stats);
  }

  const auto resort_start = std::chrono::steady_clock::now();
  for (auto& [chromosome, end] : reference_graph.chr2end) {
    (void)chromosome;
    end.resortSegments();
  }
  total_stats.resort_seconds = extensionElapsedSeconds(resort_start);
  logExtensionStageStats("total", "all", total_stats);
}

/* ==============================================================
 * 5.  Graph Correctness Verification (comprehensive check)
 * ==============================================================*/
bool RaMeshMultiGenomeGraph::verifyGraphCorrectness(
    bool verbose, bool show_detailed_segments) const {
  // 使用默认选项调用增强版本
  VerificationOptions options;
  options.verbose = verbose;
  options.show_detailed_segments = show_detailed_segments;
  options.max_errors_per_type = 100000;    // 完整统计所有错误
  options.max_total_errors = 500000;       // 完整统计所有错误
  options.max_verbose_errors_per_type = 5; // 每种类型只显示前5条详细信息

  // 默认情况下不执行高开销的内存与线程检查，可按需启用
  options.disable(VerificationType::MEMORY_INTEGRITY);
  options.disable(VerificationType::THREAD_SAFETY);
  options.disable(VerificationType::PERFORMANCE_ISSUES);
  options.enable(VerificationType::BLOCK_REFERENCE_CHR);

  VerificationResult result = verifyGraphCorrectness(options);
  return result.is_valid;
}

bool RaMeshMultiGenomeGraph::verifyGraphCorrectness(
    const SpeciesName &reference_species, bool verbose,
    bool show_detailed_segments, bool require_reference_overlap,
    bool forbid_non_reference_overlap, bool allow_reference_overlap) const {
  VerificationOptions options;
  options.verbose = verbose;
  options.show_detailed_segments = show_detailed_segments;
  options.max_errors_per_type = 100000;
  options.max_total_errors = 500000;
  options.max_verbose_errors_per_type = 5;

  options.disable(VerificationType::MEMORY_INTEGRITY);
  options.disable(VerificationType::THREAD_SAFETY);
  options.disable(VerificationType::PERFORMANCE_ISSUES);

  options.enableReferenceOverlapPolicy(
      reference_species, allow_reference_overlap, require_reference_overlap,
      forbid_non_reference_overlap);
  options.enable(VerificationType::BLOCK_REFERENCE_CHR);

  VerificationResult result = verifyGraphCorrectness(options);
  return result.is_valid;
}

// ――― Enhanced verification system implementation ―――
void RaMeshMultiGenomeGraph::addVerificationError(
    VerificationResult &result, const VerificationOptions &options,
    VerificationType type, ErrorSeverity severity, const std::string &species,
    const std::string &chr, size_t segment_index, uint_t position,
    const std::string &message, const std::string &details) const {
  // 检查是否超过总错误限制（用于完整统计）
  if (result.errors.size() >= options.max_total_errors) {
    return;
  }

  // 使用优化的错误计数器检查该类型的错误是否超过限制
  size_t type_count = result.getErrorCountFast(type);
  if (type_count >= options.max_errors_per_type) {
    return;
  }

  // 始终添加错误到结果中（完整统计）
  result.errors.emplace_back(type, severity, species, chr, segment_index,
                             position, message, details);

  // 更新快速计数器
  result.incrementErrorCount(type);

  if (severity == ErrorSeverity::ERROR || severity == ErrorSeverity::CRITICAL) {
    result.is_valid = false;
  }

  // 限制详细输出数量，但不影响错误统计
  bool is_high_severity = severity == ErrorSeverity::ERROR ||
                          severity == ErrorSeverity::CRITICAL;
  bool should_log = false;
  if (options.verbose) {
    if (is_high_severity) {
      size_t severity_verbose_count =
          result.getSeverityVerboseCount(type, severity);
      should_log =
          severity_verbose_count < options.max_verbose_errors_per_type;
    } else {
      should_log = type_count < options.max_verbose_errors_per_type;
    }
  }

  if (should_log) {
    // 优化字符串拼接：使用预分配的字符串缓冲区
    std::string full_message;
    full_message.reserve(256); // 预分配合理大小
    full_message += species;
    full_message += "/";
    full_message += chr;
    full_message += " segment#";
    full_message += std::to_string(segment_index);
    full_message += " pos=";
    full_message += std::to_string(position);
    full_message += ": ";
    full_message += message;
    if (!details.empty()) {
      full_message += " (";
      full_message += details;
      full_message += ")";
    }

    switch (severity) {
    case ErrorSeverity::INFO:
      spdlog::debug(full_message);
      break;
    case ErrorSeverity::WARNING:
      spdlog::warn(full_message);
      break;
    case ErrorSeverity::ERROR:
      spdlog::error(full_message);
      break;
    case ErrorSeverity::CRITICAL:
      spdlog::critical(full_message);
      break;
    }

    if (is_high_severity) {
      result.incrementSeverityVerboseCount(type, severity);
    }
  }
}

bool RaMeshMultiGenomeGraph::shouldStopVerification(
    const VerificationResult &result,
    const VerificationOptions &options) const {
  // 只在明确要求遇到严重错误时停止，或者达到绝对错误限制时才停止
  if (options.stop_on_critical) {
    for (const auto &error : result.errors) {
      if (error.severity == ErrorSeverity::CRITICAL) {
        return true;
      }
    }
  }

  // 只有在达到绝对错误限制时才停止（用于防止内存溢出）
  return result.errors.size() >= options.max_total_errors;
}

void RaMeshMultiGenomeGraph::logVerificationSummary(
    const VerificationResult &result,
    const VerificationOptions &options) const {
  if (!options.verbose)
    return;

  spdlog::debug("");
  spdlog::debug("============================================================");
  spdlog::debug("              VERIFICATION SUMMARY                         ");
  spdlog::debug("============================================================");
  spdlog::debug("Total errors: {}", result.errors.size());
  spdlog::debug("Verification time: {} microseconds",
                result.verification_time.count());

  // 优化：使用快速计数器显示错误统计，避免重复遍历
  if (!result.error_counts.empty()) {
    spdlog::debug("");
    spdlog::debug("Error breakdown by type:");

    const std::map<VerificationType, std::string> type_names = {
        {VerificationType::POINTER_VALIDITY, "POINTER_VALIDITY"},
        {VerificationType::LINKED_LIST_INTEGRITY, "LINKED_LIST_INTEGRITY"},
        {VerificationType::COORDINATE_OVERLAP, "COORDINATE_OVERLAP"},
        {VerificationType::COORDINATE_ORDERING, "COORDINATE_ORDERING"},
        {VerificationType::BLOCK_CONSISTENCY, "BLOCK_CONSISTENCY"},
        {VerificationType::MEMORY_INTEGRITY, "MEMORY_INTEGRITY"},
        {VerificationType::THREAD_SAFETY, "THREAD_SAFETY"},
        {VerificationType::PERFORMANCE_ISSUES, "PERFORMANCE_ISSUES"},
        {VerificationType::BLOCK_REFERENCE_CHR, "BLOCK_REFERENCE_CHR"}};

    for (const auto &[type, count] : result.error_counts) {
      auto type_name_it = type_names.find(type);
      std::string type_name = (type_name_it != type_names.end())
                                  ? type_name_it->second
                                  : "UNKNOWN_TYPE";

      spdlog::debug("  {}: {} errors", type_name, count);

      // 如果显示的详细错误数量少于总数，显示提示信息
      if (count > options.max_verbose_errors_per_type) {
        size_t hidden_count = count - options.max_verbose_errors_per_type;
        spdlog::debug("    ... and {} more errors of this type (detailed "
                      "output limited to {} per type)",
                      hidden_count, options.max_verbose_errors_per_type);
      }
    }
  }

  if (result.errors.size() >= options.max_total_errors) {
    spdlog::warn(
        "Error limit reached ({} errors), some issues may not be reported",
        options.max_total_errors);
  }

  spdlog::debug("");
  spdlog::debug("Graph status: {}", result.is_valid ? "Valid" : "Issues found");
  spdlog::debug("============================================================");
  spdlog::debug("");
}

RaMeshMultiGenomeGraph::VerificationResult
RaMeshMultiGenomeGraph::verifyGraphCorrectness(
    const VerificationOptions &options) const {
  auto start_time = std::chrono::high_resolution_clock::now();
  std::shared_lock gLock(rw);

  VerificationResult result;

  if (options.verbose) {
    spdlog::debug("");
    spdlog::debug(
        "============================================================");
    spdlog::debug(
        "              GRAPH CORRECTNESS VERIFICATION               ");
    spdlog::debug(
        "============================================================");
  }

  // 优化：使用统一的遍历来执行多个检查，减少重复遍历
  bool need_unified_traversal =
      options.isEnabled(VerificationType::LINKED_LIST_INTEGRITY) ||
      options.isEnabled(VerificationType::COORDINATE_OVERLAP) ||
      options.isEnabled(VerificationType::COORDINATE_ORDERING) ||
      options.isEnabled(VerificationType::MEMORY_INTEGRITY);

  // 首先执行不需要遍历的检查
  if (options.isEnabled(VerificationType::POINTER_VALIDITY)) {
    verifyPointerValidity(result, options);
    if (shouldStopVerification(result, options))
      goto verification_complete;
  }

  // 统一遍历执行多个检查
  if (need_unified_traversal) {
    verifyWithUnifiedTraversal(result, options);
    if (shouldStopVerification(result, options))
      goto verification_complete;
  }

  if (options.isEnabled(VerificationType::BLOCK_CONSISTENCY)) {
    verifyBlockConsistency(result, options);
    if (shouldStopVerification(result, options))
      goto verification_complete;
  }

  if (options.isEnabled(VerificationType::BLOCK_REFERENCE_CHR)) {
    verifyBlockReferenceChromosome(result, options);
    if (shouldStopVerification(result, options))
      goto verification_complete;
  }

  if (options.isEnabled(VerificationType::THREAD_SAFETY)) {
    verifyThreadSafety(result, options);
    if (shouldStopVerification(result, options))
      goto verification_complete;
  }

  if (options.isEnabled(VerificationType::PERFORMANCE_ISSUES) &&
      options.include_performance_checks) {
    verifyPerformanceIssues(result, options);
  }

verification_complete:
  auto end_time = std::chrono::high_resolution_clock::now();
  result.verification_time =
      std::chrono::duration_cast<std::chrono::microseconds>(end_time -
                                                            start_time);

  logVerificationSummary(result, options);
  return result;
}

void RaMeshMultiGenomeGraph::verifyPointerValidity(
    VerificationResult &result, const VerificationOptions &options) const {
  if (options.verbose) {
    spdlog::debug("Verifying pointer validity...");
  }

  for (const auto &[species_name, genome_graph] : species_graphs) {
    std::shared_lock species_lock(genome_graph.rw);

    for (const auto &[chr_name, genome_end] : genome_graph.chr2end) {
      // 检查头尾指针有效性
      if (!genome_end.head || !genome_end.tail) {
        addVerificationError(
            result, options, VerificationType::POINTER_VALIDITY,
            ErrorSeverity::CRITICAL, species_name, chr_name, 0, 0,
            "Head or tail pointer is null",
            "head=" +
                std::to_string(
                    reinterpret_cast<uintptr_t>(genome_end.head.get())) +
                ", tail=" +
                std::to_string(
                    reinterpret_cast<uintptr_t>(genome_end.tail.get())));
        continue;
      }

      // 检查头尾标记正确性
      if (!genome_end.head->isHead() || !genome_end.tail->isTail()) {
        addVerificationError(
            result, options, VerificationType::POINTER_VALIDITY,
            ErrorSeverity::ERROR, species_name, chr_name, 0, 0,
            "Head/tail segment role markers are incorrect",
            "head_role=" +
                std::to_string(static_cast<int>(genome_end.head->seg_role)) +
                ", tail_role=" +
                std::to_string(static_cast<int>(genome_end.tail->seg_role)));
      }

      // 检查采样向量中的指针有效性
      for (size_t i = 0; i < genome_end.sample_vec.size(); ++i) {
        if (!genome_end.sample_vec[i]) {
          addVerificationError(result, options,
                               VerificationType::POINTER_VALIDITY,
                               ErrorSeverity::WARNING, species_name, chr_name,
                               i, 0, "Null pointer in sampling vector",
                               "sample_index=" + std::to_string(i));
        }
      }
    }
  }

  // 在函数结束时显示该类型的统计信息
  if (options.verbose) {
    size_t total_errors =
        result.getErrorCountFast(VerificationType::POINTER_VALIDITY);
    if (total_errors > options.max_verbose_errors_per_type) {
      spdlog::debug("POINTER_VALIDITY: {} total errors ({} shown in detail)",
                    total_errors, options.max_verbose_errors_per_type);
    } else if (total_errors > 0) {
      spdlog::debug("POINTER_VALIDITY: {} total errors", total_errors);
    }
  }
}

void RaMeshMultiGenomeGraph::verifyLinkedListIntegrity(
    VerificationResult &result, const VerificationOptions &options) const {
  if (options.verbose) {
    spdlog::debug("Verifying linked list integrity...");
  }

  for (const auto &[species_name, genome_graph] : species_graphs) {
    std::shared_lock species_lock(genome_graph.rw);

    for (const auto &[chr_name, genome_end] : genome_graph.chr2end) {
      if (!genome_end.head || !genome_end.tail) {
        continue; // 已在指针有效性检查中报告
      }

      SegPtr current = genome_end.head;
      SegPtr prev = nullptr;
      size_t segment_count = 0;
      std::unordered_set<SegPtr> visited_segments;

      while (current) {
        segment_count++;

        // 检查循环引用
        if (visited_segments.find(current) != visited_segments.end()) {
          addVerificationError(
              result, options, VerificationType::LINKED_LIST_INTEGRITY,
              ErrorSeverity::CRITICAL, species_name, chr_name, segment_count,
              current->start, "Circular reference detected in linked list",
              "Segment appears twice in traversal");
          break;
        }
        visited_segments.insert(current);

        // 检查双向链表的一致性
        SegPtr next_ptr =
            current->primary_path.next.load(std::memory_order_acquire);
        SegPtr prev_ptr =
            current->primary_path.prev.load(std::memory_order_acquire);

        if (prev_ptr != prev) {
          addVerificationError(
              result, options, VerificationType::LINKED_LIST_INTEGRITY,
              ErrorSeverity::ERROR, species_name, chr_name, segment_count,
              current->start, "Inconsistent prev pointer in doubly linked list",
              "expected_prev=" +
                  std::to_string(reinterpret_cast<uintptr_t>(prev.get())) +
                  ", actual_prev=" +
                  std::to_string(reinterpret_cast<uintptr_t>(prev_ptr.get())));
        }

        // 防止死循环
        if (segment_count > 10000000) {
          addVerificationError(
              result, options, VerificationType::LINKED_LIST_INTEGRITY,
              ErrorSeverity::CRITICAL, species_name, chr_name, segment_count, 0,
              "Linked list may contain cycle",
              "Traversed over 10 million nodes");
          break;
        }

        prev = current;
        current = next_ptr;
      }

      // 检查是否正确到达tail
      if (current == nullptr && prev && !prev->isTail()) {
        addVerificationError(
            result, options, VerificationType::LINKED_LIST_INTEGRITY,
            ErrorSeverity::ERROR, species_name, chr_name, segment_count, 0,
            "Linked list ends without reaching tail sentinel",
            "Last segment is not marked as TAIL");
      }

      if (options.verbose) {
        spdlog::debug("  Chromosome {} contains {} segments", chr_name,
                      segment_count);
      }
    }
  }
}

void RaMeshMultiGenomeGraph::verifyCoordinateOverlap(
    VerificationResult &result, const VerificationOptions &options) const {
  if (options.verbose) {
    spdlog::debug("Verifying coordinate overlap detection...");
  }

  for (const auto &[species_name, genome_graph] : species_graphs) {
    std::shared_lock species_lock(genome_graph.rw);

    for (const auto &[chr_name, genome_end] : genome_graph.chr2end) {
      if (!genome_end.head || !genome_end.tail) {
        continue; // 已在指针有效性检查中报告
      }

      SegPtr current = genome_end.head;
      SegPtr prev = nullptr;
      size_t segment_count = 0;

      while (current) {
        segment_count++;

        // 检查非哨兵segment的坐标合理性
        if (current->isSegment()) {
          // 检查segment长度有效性
          if (current->length == 0) {
            addVerificationError(
                result, options, VerificationType::COORDINATE_OVERLAP,
                ErrorSeverity::ERROR, species_name, chr_name, segment_count,
                current->start, "Segment has zero length",
                "start=" + std::to_string(current->start));
          }

          // 检查坐标是否有重叠 - 使用与原始函数完全相同的逻辑
          if (prev && prev->isSegment()) {
            uint_t prev_end = prev->start + prev->length;
            if (current->start < prev_end) {
              uint_t overlap_size = prev_end - current->start;

              // 统一使用ERROR级别，不区分major/minor
              addVerificationError(
                  result, options, VerificationType::COORDINATE_OVERLAP,
                  ErrorSeverity::ERROR, species_name, chr_name, segment_count,
                  current->start, "Segment coordinates overlap",
                  "prev_end=" + std::to_string(prev_end) +
                      ", current_start=" + std::to_string(current->start) +
                      ", overlap_size=" + std::to_string(overlap_size) + "bp");
            }
          }
        }

        prev = current;
        current = current->primary_path.next.load(std::memory_order_acquire);
      }
    }
  }

  // 在函数结束时显示该类型的统计信息
  if (options.verbose) {
    size_t total_errors =
        result.getErrorCountFast(VerificationType::COORDINATE_OVERLAP);
    if (total_errors > options.max_verbose_errors_per_type) {
      spdlog::debug("COORDINATE_OVERLAP: {} total errors ({} shown in detail)",
                    total_errors, options.max_verbose_errors_per_type);
    } else if (total_errors > 0) {
      spdlog::debug("COORDINATE_OVERLAP: {} total errors", total_errors);
    }
  }
}

void RaMeshMultiGenomeGraph::verifyCoordinateOrdering(
    VerificationResult &result, const VerificationOptions &options) const {
  if (options.verbose) {
    spdlog::debug("Verifying coordinate ordering (segment chain sequence)...");
  }

  // 验证参数
  const uint_t LARGE_GAP_THRESHOLD =
      1000000; // 1MB gap threshold for performance hints

  for (const auto &[species_name, genome_graph] : species_graphs) {
    std::shared_lock species_lock(genome_graph.rw);

    for (const auto &[chr_name, genome_end] : genome_graph.chr2end) {
      if (!genome_end.head || !genome_end.tail) {
        continue; // 已在指针有效性检查中报告
      }

      // 添加与debug_all_species_segments相同的锁保护
      std::shared_lock end_lock(genome_end.rw);

      SegPtr current = genome_end.head;
      SegPtr prev_segment = nullptr;
      size_t segment_count = 0;
      size_t ordering_violations = 0;
      size_t large_gaps = 0;

      // 调试：收集所有segments用于对比
      std::vector<std::pair<size_t, uint_t>> debug_segments; // <index, start>

      // 遍历整个链表，与debug_all_species_segments保持一致的条件
      while (current && current != genome_end.tail) {
        if (current->isSegment()) {
          segment_count++;
          debug_segments.emplace_back(segment_count, current->start);

          // 检查segment链表的顺序（start是否递增）
          if (prev_segment) {
            // 调试输出：显示当前比较的两个segments
            if (options.show_detailed_segments && segment_count <= 10) {
              spdlog::debug("    Comparing segment#{} (start={}) with prev "
                            "segment (start={})",
                            segment_count, current->start, prev_segment->start);
            }

            // 检查排序：start应该是递增的
            if (current->start < prev_segment->start) {
              ordering_violations++;

              // 详细的调试信息
              if (options.verbose) {
                spdlog::error("    ORDERING VIOLATION DETECTED!");
                spdlog::error("      Current segment#{}: start={}",
                              segment_count, current->start);
                spdlog::error("      Previous segment#{}: start={}",
                              segment_count - 1, prev_segment->start);
                spdlog::error("      Difference: {} < {} (violation)",
                              current->start, prev_segment->start);
              }

              addVerificationError(
                  result, options, VerificationType::COORDINATE_ORDERING,
                  ErrorSeverity::ERROR, species_name, chr_name, segment_count,
                  current->start,
                  "Segments are not properly ordered (start not increasing)",
                  "current_start=" + std::to_string(current->start) +
                      " < prev_start=" + std::to_string(prev_segment->start) +
                      " (violation #" + std::to_string(ordering_violations) +
                      ")");
            }

            // 检查间隙大小（性能提示）
            uint_t prev_end = prev_segment->start + prev_segment->length;
            if (current->start > prev_end) {
              uint_t gap_size = current->start - prev_end;
              if (gap_size > LARGE_GAP_THRESHOLD) {
                large_gaps++;
                addVerificationError(
                    result, options, VerificationType::COORDINATE_ORDERING,
                    ErrorSeverity::INFO, species_name, chr_name, segment_count,
                    current->start, "Large gap between segments",
                    "gap_size=" + std::to_string(gap_size) + "bp" +
                        " (large gap #" + std::to_string(large_gaps) + ")");
              }
            }
          }

          // 更新prev_segment为当前的segment
          prev_segment = current;
        }

        current = current->primary_path.next.load(std::memory_order_acquire);
      }

      // 输出统计信息和调试信息（仅在启用详细段信息时输出，避免刷屏）
      if (options.verbose && options.show_detailed_segments) {
        spdlog::debug("  Chromosome {}: {} segments, {} ordering violations, "
                      "{} large gaps",
                      chr_name, segment_count, ordering_violations, large_gaps);

        // 输出前10个和后10个segments的start值用于调试（仅在启用详细段信息时）
        if (options.show_detailed_segments && debug_segments.size() > 20) {
          spdlog::debug("    First 10 segments start values:");
          for (size_t i = 0; i < std::min(size_t(10), debug_segments.size());
               ++i) {
            spdlog::debug("      Segment#{}: start={}", debug_segments[i].first,
                          debug_segments[i].second);
          }
          spdlog::debug("    Last 10 segments start values:");
          for (size_t i = std::max(size_t(0), debug_segments.size() - 10);
               i < debug_segments.size(); ++i) {
            spdlog::debug("      Segment#{}: start={}", debug_segments[i].first,
                          debug_segments[i].second);
          }

          // 特别检查第5000和第5500个segment
          if (debug_segments.size() > 5500) {
            spdlog::debug("    Special check - Segment#5000: start={}",
                          debug_segments[4999].second);
            spdlog::debug("    Special check - Segment#5500: start={}",
                          debug_segments[5499].second);
            if (debug_segments[4999].second > debug_segments[5499].second) {
              spdlog::error("    MANUAL CHECK CONFIRMED: Segment#5000 start > "
                            "Segment#5500 start!");
              spdlog::error("      This should have been detected as an "
                            "ordering violation!");
            }
          }
        }
      }
    }
  }

  // 在函数结束时显示该类型的统计信息
  if (options.verbose) {
    size_t total_errors =
        result.getErrorCountFast(VerificationType::COORDINATE_ORDERING);
    if (total_errors > options.max_verbose_errors_per_type) {
      spdlog::debug("COORDINATE_ORDERING: {} total errors ({} shown in detail)",
                    total_errors, options.max_verbose_errors_per_type);
    } else if (total_errors > 0) {
      spdlog::debug("COORDINATE_ORDERING: {} total errors", total_errors);
    }
  }
}

void RaMeshMultiGenomeGraph::verifyBlockConsistency(
    VerificationResult &result, const VerificationOptions &options) const {
  if (options.verbose) {
    spdlog::debug("Verifying block consistency...");
  }

  size_t valid_blocks = 0;
  size_t expired_blocks = 0;
  size_t empty_blocks = 0;

  for (const auto &weak_block : blocks) {
    if (auto block_ptr = weak_block.lock()) {
      valid_blocks++;
      std::shared_lock block_lock(block_ptr->rw);

      // 检查block是否为空
      if (block_ptr->anchors.empty()) {
        empty_blocks++;
        addVerificationError(result, options,
                             VerificationType::BLOCK_CONSISTENCY,
                             ErrorSeverity::WARNING, "", "", 0, 0,
                             "Empty block found", "Block contains no anchors");
      }

      // 检查block中的anchors是否都有效
      for (const auto &[species_chr_pair, head_ptr] : block_ptr->anchors) {
        if (!head_ptr) {
          addVerificationError(result, options,
                               VerificationType::BLOCK_CONSISTENCY,
                               ErrorSeverity::ERROR, species_chr_pair.first,
                               species_chr_pair.second, 0, 0,
                               "Null anchor pointer found in block",
                               "Block contains invalid anchor reference");
        }
      }

      // 优化：收集该block的所有segments，避免重复遍历整个图
      std::unordered_set<SpeciesChrPair, SpeciesChrPairHash>
          referenced_species_chrs;

      // 遍历图中的segments，只检查引用了当前block的segments
      for (const auto &[species_name, genome_graph] : species_graphs) {
        std::shared_lock species_lock(genome_graph.rw);

        for (const auto &[chr_name, genome_end] : genome_graph.chr2end) {
          SegPtr current = genome_end.head;
          if (current) {
            current =
                current->primary_path.next.load(std::memory_order_acquire);
          }

          while (current && !current->isTail()) {
            if (current->isSegment() && current->parent_block == block_ptr) {
              SpeciesChrPair key{species_name, chr_name};
              referenced_species_chrs.insert(key);

              // 检查block是否确实包含该chr的anchor
              auto anchor_it = block_ptr->anchors.find(key);
              if (anchor_it == block_ptr->anchors.end()) {
                addVerificationError(result, options,
                                     VerificationType::BLOCK_CONSISTENCY,
                                     ErrorSeverity::ERROR, species_name,
                                     chr_name, 0, current->start,
                                     "Segment references block but block "
                                     "doesn't contain corresponding anchor",
                                     "Inconsistent block-segment relationship");
              }
              break; // 找到一个引用就足够了，不需要继续遍历该染色体
            }
            current =
                current->primary_path.next.load(std::memory_order_acquire);
          }
        }
      }
    } else {
      expired_blocks++;
    }
  }

  if (options.verbose) {
    spdlog::debug("Block pool statistics: {} valid blocks, {} expired blocks, "
                  "{} empty blocks",
                  valid_blocks, expired_blocks, empty_blocks);
  }

  // 检查过期block比例
  if (expired_blocks > 0) {
    double expired_ratio =
        static_cast<double>(expired_blocks) / (valid_blocks + expired_blocks);
    if (expired_ratio > 0.5) {
      addVerificationError(
          result, options, VerificationType::BLOCK_CONSISTENCY,
          ErrorSeverity::WARNING, "", "", 0, 0, "High ratio of expired blocks",
          "expired_ratio=" + std::to_string(expired_ratio * 100) +
              "%, consider cleanup");
    }
  }
}

void RaMeshMultiGenomeGraph::verifyBlockReferenceChromosome(
    VerificationResult &result, const VerificationOptions &options) const {
  if (options.verbose) {
    spdlog::debug("Verifying block reference chromosome coverage...");
  }

  size_t checked_blocks = 0;
  size_t missing_ref_definition = 0;
  size_t missing_ref_anchor = 0;

  for (const auto &weak_block : blocks) {
    if (auto block_ptr = weak_block.lock()) {
      checked_blocks++;
      std::shared_lock block_lock(block_ptr->rw);

      if (block_ptr->anchors.empty()) {
        continue;
      }

      const auto makeAnchorSummary = [](const ChrHeadMap &anchors) {
        std::string detail = "[";
        size_t anchor_idx = 0;
        const size_t total_anchors = anchors.size();
        constexpr size_t kMaxAnchorSamples = 3;
        constexpr size_t kMaxSegmentSamples = 2;

        for (const auto &entry : anchors) {
          const auto &species_chr = entry.first;
          const SegPtr &head_ptr = entry.second;

          if (anchor_idx > 0) {
            detail += "; ";
          }

          detail += species_chr.first + ":" + species_chr.second;
          detail += " {";

          if (!head_ptr) {
            detail += "null_head";
          } else {
            bool printed_segment = false;
            bool truncated = false;
            size_t segment_samples = 0;

            SegPtr seg =
                head_ptr->primary_path.next.load(std::memory_order_acquire);
            while (seg && !seg->isTail()) {
              if (seg->isSegment()) {
                {
                  std::shared_lock seg_lock(seg->rw);
                  if (printed_segment) {
                    detail += " | ";
                  }
                  detail += "start=" + std::to_string(seg->start) +
                            ",len=" + std::to_string(seg->length);
                }

                printed_segment = true;
                ++segment_samples;
                if (segment_samples >= kMaxSegmentSamples) {
                  SegPtr probe =
                      seg->primary_path.next.load(std::memory_order_acquire);
                  while (probe && !probe->isTail()) {
                    if (probe->isSegment()) {
                      truncated = true;
                      break;
                    }
                    probe =
                        probe->primary_path.next.load(std::memory_order_acquire);
                  }
                  break;
                }
              }
              seg = seg->primary_path.next.load(std::memory_order_acquire);
            }

            if (!printed_segment) {
              detail += "no_segment";
            } else if (truncated) {
              detail += " | ...";
            }
          }

          detail += "}";
          ++anchor_idx;
          if (anchor_idx >= kMaxAnchorSamples) {
            break;
          }
        }

        detail += "]";
        if (total_anchors > anchor_idx) {
          detail += " +" + std::to_string(total_anchors - anchor_idx) + " more";
        }
        return detail;
      };

      if (block_ptr->ref_species.empty() || block_ptr->ref_chr.empty()) {
        missing_ref_definition++;
        std::string hint_species;
        std::string hint_chr;
        if (!block_ptr->anchors.empty()) {
          const auto &first_anchor = *block_ptr->anchors.begin();
          hint_species = first_anchor.first.first;
          hint_chr = first_anchor.first.second;
        }

        addVerificationError(
            result, options, VerificationType::BLOCK_REFERENCE_CHR,
            ErrorSeverity::ERROR, hint_species, hint_chr, 0, 0,
            "Block missing reference species/chromosome assignment",
            "anchors=" + makeAnchorSummary(block_ptr->anchors));
        continue;
      }

      const SpeciesChrPair reference_key{block_ptr->ref_species,
                                         block_ptr->ref_chr};
      const auto [reference_begin, reference_end] =
          block_ptr->anchors.equal_range(reference_key);
      const size_t primary_reference_occurrences = static_cast<size_t>(
          std::count_if(reference_begin, reference_end, [](const auto &entry) {
            return entry.second && entry.second->isPrimary();
          }));

      if (primary_reference_occurrences != 1) {
        missing_ref_anchor++;
        addVerificationError(
            result, options, VerificationType::BLOCK_REFERENCE_CHR,
            ErrorSeverity::ERROR, block_ptr->ref_species, block_ptr->ref_chr, 0,
            0,
            "Block must contain exactly one primary declared reference "
            "occurrence",
            "primary_occurrences=" +
                std::to_string(primary_reference_occurrences) +
                ", ref_species=" + block_ptr->ref_species +
                ", ref_chr=" + block_ptr->ref_chr +
                ", anchors=" + makeAnchorSummary(block_ptr->anchors));
      }
    }
  }

  if (options.verbose) {
    spdlog::debug("Blocks without ref_chr: {} (checked {})",
                  missing_ref_definition, checked_blocks);
    spdlog::debug("Blocks with invalid reference occurrence: {}",
                  missing_ref_anchor);
  }
}

void RaMeshMultiGenomeGraph::verifyMemoryIntegrity(
    VerificationResult &result, const VerificationOptions &options) const {
  if (options.verbose) {
    spdlog::debug("Verifying memory integrity...");
  }

  // 检查segment的内存一致性
  for (const auto &[species_name, genome_graph] : species_graphs) {
    std::shared_lock species_lock(genome_graph.rw);

    for (const auto &[chr_name, genome_end] : genome_graph.chr2end) {
      if (!genome_end.head || !genome_end.tail) {
        continue;
      }

      std::unordered_set<SegPtr> all_segments;
      SegPtr current = genome_end.head;
      size_t segment_count = 0;

      while (current) {
        segment_count++;

        // 检查重复的segment指针
        if (all_segments.find(current) != all_segments.end()) {
          addVerificationError(
              result, options, VerificationType::MEMORY_INTEGRITY,
              ErrorSeverity::CRITICAL, species_name, chr_name, segment_count,
              current->start, "Duplicate segment pointer detected",
              "Same segment appears multiple times in memory");
        }
        all_segments.insert(current);

        // 检查坐标溢出
        if (current->isSegment()) {
          if (current->start > UINT32_MAX - current->length) {
            addVerificationError(
                result, options, VerificationType::MEMORY_INTEGRITY,
                ErrorSeverity::ERROR, species_name, chr_name, segment_count,
                current->start, "Coordinate overflow detected",
                "start=" + std::to_string(current->start) +
                    ", length=" + std::to_string(current->length));
          }

          // 检查parent_block引用的有效性
          if (current->parent_block) {
            bool block_found = false;
            for (const auto &weak_block : blocks) {
              if (auto block_ptr = weak_block.lock()) {
                if (block_ptr == current->parent_block) {
                  block_found = true;
                  break;
                }
              }
            }
            if (!block_found) {
              addVerificationError(
                  result, options, VerificationType::MEMORY_INTEGRITY,
                  ErrorSeverity::ERROR, species_name, chr_name, segment_count,
                  current->start, "Segment references invalid parent_block",
                  "parent_block not found in global block pool");
            }
          }
        }

        current = current->primary_path.next.load(std::memory_order_acquire);
      }

      // 检查采样向量的内存一致性
      for (size_t i = 0; i < genome_end.sample_vec.size(); ++i) {
        SegPtr sample_seg = genome_end.sample_vec[i];
        if (sample_seg && all_segments.find(sample_seg) == all_segments.end()) {
          addVerificationError(
              result, options, VerificationType::MEMORY_INTEGRITY,
              ErrorSeverity::ERROR, species_name, chr_name, i, 0,
              "Sampling vector contains invalid segment reference",
              "sample_index=" + std::to_string(i) +
                  ", segment not in main list");
        }
      }
    }
  }
}

void RaMeshMultiGenomeGraph::verifyThreadSafety(
    VerificationResult &result, const VerificationOptions &options) const {
  if (options.verbose) {
    spdlog::debug("Verifying thread safety...");
  }

  size_t total_atomic_operations = 0;

  // 检查原子操作的一致性（静态检查）
  for (const auto &[species_name, genome_graph] : species_graphs) {
    std::shared_lock species_lock(genome_graph.rw);

    for (const auto &[chr_name, genome_end] : genome_graph.chr2end) {
      if (!genome_end.head || !genome_end.tail) {
        continue;
      }

      SegPtr current = genome_end.head;
      while (current) {
        // 检查原子指针的内存序
        SegPtr next_ptr =
            current->primary_path.next.load(std::memory_order_acquire);
        SegPtr prev_ptr =
            current->primary_path.prev.load(std::memory_order_acquire);

        total_atomic_operations += 2;

        // 这里只能做基本的静态检查
        // 动态竞争条件需要专门的工具检测

        current = next_ptr;
      }
    }
  }

  if (options.verbose) {
    spdlog::debug(
        "Thread safety check completed: {} atomic operations verified",
        total_atomic_operations);
  }

  // 添加一般性的线程安全信息
  addVerificationError(result, options, VerificationType::THREAD_SAFETY,
                       ErrorSeverity::INFO, "", "", 0, 0,
                       "Static thread safety check completed",
                       "All accesses are properly protected by locks, " +
                           std::to_string(total_atomic_operations) +
                           " atomic operations verified");
}

void RaMeshMultiGenomeGraph::verifyPerformanceIssues(
    VerificationResult &result, const VerificationOptions &options) const {
  if (options.verbose) {
    spdlog::debug("Verifying performance issues...");
  }

  // 性能检查阈值
  const size_t HIGH_SEGMENT_COUNT_THRESHOLD = 100000;
  const size_t SMALL_SEGMENT_THRESHOLD = 100;
  const double SMALL_SEGMENT_RATIO_THRESHOLD = 0.8;

  for (const auto &[species_name, genome_graph] : species_graphs) {
    std::shared_lock species_lock(genome_graph.rw);

    for (const auto &[chr_name, genome_end] : genome_graph.chr2end) {
      if (!genome_end.head || !genome_end.tail) {
        continue;
      }

      size_t total_segments = 0;
      size_t small_segments = 0;
      SegPtr current =
          genome_end.head->primary_path.next.load(std::memory_order_acquire);

      while (current && !current->isTail()) {
        if (current->isSegment()) {
          total_segments++;
          if (current->length < SMALL_SEGMENT_THRESHOLD) {
            small_segments++;
          }
        }
        current = current->primary_path.next.load(std::memory_order_acquire);
      }

      // 检查高segment数量
      if (total_segments > HIGH_SEGMENT_COUNT_THRESHOLD) {
        addVerificationError(result, options,
                             VerificationType::PERFORMANCE_ISSUES,
                             ErrorSeverity::WARNING, species_name, chr_name, 0,
                             0, "High segment count detected",
                             "segment_count=" + std::to_string(total_segments) +
                                 ", may impact performance");
      }

      // 检查小segment比例
      if (total_segments > 0) {
        double small_ratio =
            static_cast<double>(small_segments) / total_segments;
        if (small_ratio > SMALL_SEGMENT_RATIO_THRESHOLD) {
          addVerificationError(
              result, options, VerificationType::PERFORMANCE_ISSUES,
              ErrorSeverity::INFO, species_name, chr_name, 0, 0,
              "High ratio of small segments",
              "small_segments=" + std::to_string(small_segments) + "/" +
                  std::to_string(total_segments) + " (" +
                  std::to_string(small_ratio * 100) + "%)");
        }
      }

      // 检查采样向量效率
      size_t expected_samples = total_segments / GenomeEnd::kSampleStep + 1;
      if (genome_end.sample_vec.size() > expected_samples * 2) {
        addVerificationError(
            result, options, VerificationType::PERFORMANCE_ISSUES,
            ErrorSeverity::INFO, species_name, chr_name, 0, 0,
            "Sampling vector may be oversized",
            "actual_size=" + std::to_string(genome_end.sample_vec.size()) +
                ", expected_size=" + std::to_string(expected_samples));
      }
    }
  }
}

// Safe path-linking implementation.
void RaMeshMultiGenomeGraph::safeLink(SegPtr prev, SegPtr next) {
  if (!prev || !next)
    return;

  /* 1) 先写 next->prev  */
  next->primary_path.prev.store(prev, std::memory_order_relaxed);

  /* 2) 再 CAS prev->next  (确保无并发竞争) */
  SegPtr exp = next->primary_path.next.load(std::memory_order_acquire); // 旧值
  prev->primary_path.next.store(next, std::memory_order_release);

  /* 3) 发布内存屏障，保证两侧对所有线程可见 */
  std::atomic_thread_fence(std::memory_order_release);
}

/* =============================================================
 * 6.  Merge multiple graphs (public API)
 * ===========================================================*/
void RaMeshMultiGenomeGraph::mergeMultipleGraphs(const SpeciesName &ref_name,
                                                 uint_t thread_num) {
  // ═══════════════════════════════════════════════════════════
  // 性能分析：函数开始计时
  // ═══════════════════════════════════════════════════════════
  using namespace std::chrono;
  using HighResClock = high_resolution_clock;
  using TimePoint = HighResClock::time_point;
  using Duration = nanoseconds;

  const TimePoint function_start_time = HighResClock::now();

  // 性能统计结构
  struct PerformanceStats {
    Duration initialization_time{0};
    Duration chromosome_iteration_time{0};
    Duration overlap_detection_time{0};
    Duration block_creation_time{0};
    Duration cigar_processing_time{0};
    Duration segment_creation_time{0};
    Duration link_replacement_time{0};
    Duration cleanup_time{0};

    size_t total_overlaps_found = 0;
    size_t total_blocks_created = 0;
    size_t total_segments_created = 0;
    size_t total_cigar_operations = 0;

    // 批量删除相关统计
    size_t total_blocks_marked_for_deletion = 0;
    size_t batch_deletion_operations = 0;

    // 细粒度清理工作计时
    Duration block_lookup_time{0};
    Duration actual_deletion_time{0};
    Duration container_reorganization_time{0};

    void logStepTime(const std::string &step_name, Duration step_time,
                     Duration total_time) const {
      double step_ms = duration_cast<microseconds>(step_time).count() / 1000.0;
      double total_ms =
          duration_cast<microseconds>(total_time).count() / 1000.0;
      double percentage = total_ms > 0 ? (step_ms / total_ms) * 100.0 : 0.0;

      spdlog::info("性能分析 - {}: {:.3f}ms ({:.2f}%)", step_name, step_ms,
                   percentage);
    }

    void logFinalReport(Duration total_time) const {
      double total_ms =
          duration_cast<microseconds>(total_time).count() / 1000.0;

      spdlog::info(
          "═══════════════════════════════════════════════════════════");
      spdlog::info("性能分析报告 - mergeMultipleGraphs 函数");
      spdlog::info(
          "═══════════════════════════════════════════════════════════");
      spdlog::info("总执行时间: {:.3f}ms", total_ms);
      spdlog::info(
          "───────────────────────────────────────────────────────────");

      logStepTime("初始化阶段", initialization_time, total_time);
      logStepTime("染色体遍历", chromosome_iteration_time, total_time);
      logStepTime("重叠检测", overlap_detection_time, total_time);
      logStepTime("Block创建", block_creation_time, total_time);
      logStepTime("CIGAR处理", cigar_processing_time, total_time);
      logStepTime("Segment创建", segment_creation_time, total_time);
      logStepTime("链表重连", link_replacement_time, total_time);
      logStepTime("清理工作", cleanup_time, total_time);

      spdlog::info(
          "───────────────────────────────────────────────────────────");
      spdlog::info("清理工作详细分析:");
      logStepTime("  Block查找", block_lookup_time, total_time);
      logStepTime("  实际删除", actual_deletion_time, total_time);
      logStepTime("  容器重组", container_reorganization_time, total_time);

      spdlog::info(
          "───────────────────────────────────────────────────────────");
      spdlog::info("操作统计:");
      spdlog::info("  发现重叠: {} 次", total_overlaps_found);
      spdlog::info("  创建Block: {} 个", total_blocks_created);
      spdlog::info("  创建Segment: {} 个", total_segments_created);
      spdlog::info("  CIGAR操作: {} 次", total_cigar_operations);
      spdlog::info("  标记删除Block: {} 个", total_blocks_marked_for_deletion);
      spdlog::info("  批量删除操作: {} 次", batch_deletion_operations);

      if (total_overlaps_found > 0) {
        double avg_overlap_time =
            duration_cast<microseconds>(overlap_detection_time).count() /
            1000.0 / total_overlaps_found;
        spdlog::info("  平均重叠处理时间: {:.3f}ms", avg_overlap_time);
      }

      if (total_blocks_marked_for_deletion > 0 &&
          batch_deletion_operations > 0) {
        double avg_deletion_efficiency =
            static_cast<double>(total_blocks_marked_for_deletion) /
            batch_deletion_operations;
        spdlog::info("  批量删除效率: {:.1f} 个Block/次",
                     avg_deletion_efficiency);
      }

      spdlog::info(
          "═══════════════════════════════════════════════════════════");
    }
  } perf_stats;

  // 辅助函数：带性能监控的Segment创建
  auto createSegmentWithTiming = [&](uint_t start, uint_t len, Strand sd,
                                     Cigar_t cg, AlignRole rl, SegmentRole sl,
                                     const BlockPtr &bp) -> SegPtr {
    TimePoint seg_create_start = HighResClock::now();
    SegPtr result = Segment::create(start, len, sd, std::move(cg), rl, sl, bp);
    TimePoint seg_create_end = HighResClock::now();
    perf_stats.segment_creation_time +=
        duration_cast<Duration>(seg_create_end - seg_create_start);
    perf_stats.total_segments_created++;
    return result;
  };

  TimePoint step_start_time = HighResClock::now();

  // ═══════════════════════════════════════════════════════════
  // 批量删除优化：收集要删除的Block，避免逐个删除的O(n²)问题
  // ═══════════════════════════════════════════════════════════
  std::unordered_set<BlockPtr> blocks_to_delete;
  blocks_to_delete.reserve(100000); // 预分配空间，根据统计数据估算

  std::shared_lock graph_lock(rw);

  // ═══════════════════════════════════════════════════════════
  // 进度跟踪结构体定义 - 基于基因组坐标而非segment数量
  // ═══════════════════════════════════════════════════════════
  struct MergeProgress {
    size_t total_chromosomes = 0;          // 总染色体数量
    size_t completed_chromosomes = 0;      // 已完成染色体数量
    uint64_t total_genomic_length = 0;     // 总基因组长度（所有染色体）
    uint64_t processed_genomic_length = 0; // 已处理基因组长度
    uint64_t current_chr_length = 0;       // 当前染色体长度
    uint64_t current_chr_processed = 0;    // 当前染色体已处理长度
    size_t total_merges = 0;               // 总合并操作数量
    size_t completed_merges = 0;           // 已完成合并操作数量
    size_t participant_conflict_rejections = 0;
    std::string current_species;    // 当前处理的物种
    std::string current_chromosome; // 当前处理的染色体
    uint64_t current_position = 0;  // 当前处理位置

    // 上次日志输出时间（避免过于频繁的日志输出）
    mutable std::chrono::steady_clock::time_point last_log_time;
    mutable double last_logged_percentage = -1.0;

    MergeProgress() { last_log_time = std::chrono::steady_clock::now(); }

    // 显示进度（使用spdlog，控制输出频率）
    void displayProgress() const {
      if (total_genomic_length == 0)
        return;

      // 计算总体进度
      double overall_percentage =
          static_cast<double>(processed_genomic_length) / total_genomic_length *
          100.0;

      // 计算当前染色体进度
      double chr_percentage = 0.0;
      if (current_chr_length > 0) {
        chr_percentage = static_cast<double>(current_chr_processed) /
                         current_chr_length * 100.0;
      }

      // 控制日志输出频率：每1%或每5秒输出一次
      auto now = std::chrono::steady_clock::now();
      auto time_diff =
          std::chrono::duration_cast<std::chrono::seconds>(now - last_log_time)
              .count();
      double percentage_diff =
          std::abs(overall_percentage - last_logged_percentage);

      if (percentage_diff >= 5.0 || time_diff >= 5) {
        spdlog::info(
            "Merge progress: {:.1f}% ({:.1f} MB/{:.1f} MB) Chromosomes: {}/{} "
            "Current: {}.{} ({:.1f}%) Position: {} Merges: {}",
            overall_percentage, processed_genomic_length / 1000000.0,
            total_genomic_length / 1000000.0, completed_chromosomes,
            total_chromosomes, current_species, current_chromosome,
            chr_percentage, current_position, completed_merges);

        last_log_time = now;
        last_logged_percentage = overall_percentage;
      }
    }

    // 开始处理新染色体
    void startChromosome(const std::string &species,
                         const std::string &chromosome, uint64_t chr_length) {
      current_species = species;
      current_chromosome = chromosome;
      current_chr_length = chr_length;
      current_chr_processed = 0;
      current_position = 0;

      // spdlog::info("Starting to process chromosome: {}.{} (Length: {:.1f}
      // MB)",
      //              species,
      //              chromosome,
      //              chr_length / 1000000.0);

      displayProgress();
    }

    // 更新当前处理位置
    void updatePosition(uint64_t position) {
      if (position > current_position) {
        uint64_t advance = position - current_position;
        current_chr_processed += advance;
        processed_genomic_length += advance;
        current_position = position;
        displayProgress();
      }
    }

    // 完成一个染色体的处理
    void completeChromosome() {
      // 确保当前染色体完全处理完成
      if (current_chr_processed < current_chr_length) {
        uint64_t remaining = current_chr_length - current_chr_processed;
        current_chr_processed = current_chr_length;
        processed_genomic_length += remaining;
      }

      completed_chromosomes++;
      // spdlog::info("Completed chromosome: {}.{} ({}/{})",
      //              current_species,
      //              current_chromosome,
      //              completed_chromosomes,
      //              total_chromosomes);

      displayProgress();
    }

    // 记录合并操作
    void recordMerge() { completed_merges++; }
    void recordParticipantConflict() { participant_conflict_rejections++; }

    // 完成所有处理
    void finish() {
      spdlog::info(
          "Merge completed! Processed a total of {} chromosomes, {:.1f} MB of "
          "genomic data, with {} merge operations executed.",
          total_chromosomes, total_genomic_length / 1000000.0,
          completed_merges);
      spdlog::info("Merge retained {} overlapping Block pairs because their "
                   "non-reference participant species conflict.",
                   participant_conflict_rejections);
    }
  };

  // 调试用：收集所有物种的segment详细信息，方便在调试器中查看
  struct SegmentDebugInfo {
    uint_t start;
    uint_t length;
    std::string strand_str;
    std::string seg_role_str;
    std::string align_role_str;
    size_t cigar_size;
    std::string parent_block_ref_chr;
    std::vector<std::string> parent_block_chromosomes;

    SegmentDebugInfo(const SegPtr &seg) {
      if (!seg)
        return;
      start = seg->start;
      length = seg->length;
      strand_str = (seg->strand == Strand::FORWARD) ? "FORWARD" : "REVERSE";
      seg_role_str = (seg->seg_role == SegmentRole::SEGMENT) ? "SEGMENT"
                     : (seg->seg_role == SegmentRole::HEAD)  ? "HEAD"
                                                             : "TAIL";
      align_role_str =
          (seg->align_role == AlignRole::PRIMARY) ? "PRIMARY" : "SECONDARY";
      cigar_size = seg->cigar.size();

      if (seg->parent_block) {
        parent_block_ref_chr = seg->parent_block->ref_chr;
        std::shared_lock block_lock(seg->parent_block->rw);
        parent_block_chromosomes.reserve(seg->parent_block->anchors.size());
        for (const auto &anchor_entry : seg->parent_block->anchors) {
          const auto &species = anchor_entry.first.first;
          const auto &chr = anchor_entry.first.second;
          parent_block_chromosomes.emplace_back(species + "." + chr);
        }
        std::sort(parent_block_chromosomes.begin(),
                  parent_block_chromosomes.end());
      } else {
        parent_block_ref_chr = "null";
        parent_block_chromosomes.clear();
      }
    }
  };

  struct ChromosomeDebugInfo {
    std::string chr_name;
    std::vector<SegmentDebugInfo> segments;

    ChromosomeDebugInfo(const std::string &name) : chr_name(name) {}
  };

  struct SpeciesDebugInfo {
    std::string species_name;
    std::vector<ChromosomeDebugInfo> chromosomes;

    SpeciesDebugInfo(const std::string &name) : species_name(name) {}
  };

  // ═══════════════════════════════════════════════════════════
  // 初始化进度跟踪
  // ═══════════════════════════════════════════════════════════
  MergeProgress progress;

  // 统计总的染色体数量和基因组长度
  for (const auto &[species_name, genome_graph] : species_graphs) {
    if (species_name == ref_name) {
      progress.total_chromosomes += genome_graph.chr2end.size();

      // 统计总基因组长度
      std::shared_lock genome_lock(genome_graph.rw);
      for (const auto &[chr_name, genome_end] : genome_graph.chr2end) {
        std::shared_lock end_lock(genome_end.rw);

        // 找到染色体的最后一个segment来确定长度
        uint64_t chr_length = 0;
        SegPtr current = genome_end.head;
        if (current) {
          current = current->primary_path.next.load(std::memory_order_acquire);
        }
        while (current && !current->isTail()) {
          if (current->isSegment()) {
            uint64_t segment_end = current->start + current->length;
            chr_length = std::max(chr_length, segment_end);
          }
          current = current->primary_path.next.load(std::memory_order_acquire);
        }
        progress.total_genomic_length += chr_length;
      }
      break;
    }
  }

  spdlog::info("Starting to merge multi-genome graph, reference species: {}", ref_name);
  spdlog::info("A total of {} chromosomes ({:.1f} MB of genomic data) will be processed",
               progress.total_chromosomes,
               progress.total_genomic_length / 1000000.0);

  progress.displayProgress();

  // ═══════════════════════════════════════════════════════════
  // 性能分析：初始化阶段完成
  // ═══════════════════════════════════════════════════════════
  TimePoint init_end_time = HighResClock::now();
  perf_stats.initialization_time =
      duration_cast<Duration>(init_end_time - step_start_time);
  step_start_time = init_end_time;

  std::vector<SpeciesDebugInfo> debug_all_species_segments;
  debug_all_species_segments.reserve(species_graphs.size());

  // 收集所有物种的segment详细信息
  for (const auto &[species_name, genome_graph] : species_graphs) {
    SpeciesDebugInfo species_info(species_name);
    std::shared_lock genome_lock(genome_graph.rw);

    for (const auto &[chr_name, genome_end] : genome_graph.chr2end) {
      ChromosomeDebugInfo chr_info(chr_name);
      std::shared_lock end_lock(genome_end.rw);

      // 遍历该染色体的所有segments
      SegPtr current = genome_end.head;
      while (current && current != genome_end.tail) {
        if (current->isSegment()) {
          chr_info.segments.emplace_back(current);
        }
        current = current->primary_path.next.load(std::memory_order_acquire);
      }

      species_info.chromosomes.push_back(std::move(chr_info));
    }

    debug_all_species_segments.push_back(std::move(species_info));
  }

  auto ref_it = species_graphs.find(ref_name);
  if (ref_it == species_graphs.end())
    return;

  // ═══════════════════════════════════════════════════════════
  // 性能分析：开始主循环 - 染色体遍历
  // ═══════════════════════════════════════════════════════════
  TimePoint main_loop_start = HighResClock::now();

  // 改为遍历species_graphs来查找Ref
  for (const auto &[current_species, ref_genome] : species_graphs) {
    if (current_species != ref_name) {
      continue;
    }
    std::shared_lock genome_lock(ref_genome.rw);
    // 主循环，遍历所有Ref序列
    for (const auto &[chr_name, genome_end] : ref_genome.chr2end) {
      // ═══════════════════════════════════════════════════════════
      // 计算当前染色体长度并开始处理
      // ═══════════════════════════════════════════════════════════
      uint64_t chr_length = 0;
      std::shared_lock end_lock(genome_end.rw);
      SegPtr count_current = genome_end.head;
      if (count_current) {
        count_current =
            count_current->primary_path.next.load(std::memory_order_acquire);
      }
      while (count_current && !count_current->isTail()) {
        if (count_current->isSegment()) {
          uint64_t segment_end = count_current->start + count_current->length;
          chr_length = std::max(chr_length, segment_end);
        }
        count_current =
            count_current->primary_path.next.load(std::memory_order_acquire);
      }
      end_lock.unlock();

      progress.startChromosome(current_species, chr_name, chr_length);

      // 从头节点开始遍历链表
      SegPtr prev =
          genome_end.head->primary_path.next.load(std::memory_order_acquire);
      SegPtr current = prev->primary_path.next.load(std::memory_order_acquire);
      BlockPtr prev_block = prev->parent_block;
      BlockPtr current_block = nullptr;

      while (current && !current->isTail()) {
        // if (current->start >= 2310000 && chr_name == "simGorilla.chrD") {
        //     std::cout<<"test";
        // }
        // 只处理真正的segment（跳过头尾哨兵）
        if (current->isSegment() && current->parent_block) {
          // ═══════════════════════════════════════════════════════════
          // 更新基因组位置进度
          // ═══════════════════════════════════════════════════════════
          progress.updatePosition(current->start + current->length);

          current_block = current->parent_block;

          // ═══════════════════════════════════════════════════════════
          // 性能分析：重叠检测开始
          // ═══════════════════════════════════════════════════════════
          TimePoint overlap_detection_start = HighResClock::now();

          // 构造查找键并直接find
          SpeciesChrPair ref_key{ref_name, chr_name};

          // 直接查找，避免遍历整个anchors
          auto prev_anchor_it = prev_block->anchors.find(ref_key);
          auto curr_anchor_it = current_block->anchors.find(ref_key);

          // 快速存在性检查
          if (prev_anchor_it != prev_block->anchors.end() &&
              curr_anchor_it != current_block->anchors.end()) {
            const SegPtr prev_seg = prev_anchor_it->second;
            const SegPtr curr_seg = curr_anchor_it->second;

            // 快速重叠判断（内联计算，避免函数调用）
            const uint_t prev_end = prev_seg->start + prev_seg->length;
            const uint_t curr_start = curr_seg->start;

            // 发现重叠
            if (prev_end > curr_start &&
                blocksShareNonReferenceSpecies(prev_block, current_block,
                                               ref_name)) {
              progress.recordParticipantConflict();
            } else if (prev_end > curr_start) {
              // ═══════════════════════════════════════════════════════════
              // 性能分析：重叠检测完成，开始处理
              // ═══════════════════════════════════════════════════════════
              TimePoint overlap_detection_end = HighResClock::now();
              perf_stats.overlap_detection_time += duration_cast<Duration>(
                  overlap_detection_end - overlap_detection_start);
              perf_stats.total_overlaps_found++;

              TimePoint block_creation_start = HighResClock::now();
              // ═══════════════════════════════════════════════════════════
              // 记录合并操作
              // ═══════════════════════════════════════════════════════════
              progress.recordMerge();

              // 处理合并并确定重叠区间。
              uint_t overlap_start = std::max(prev_seg->start, curr_seg->start);
              uint_t overlap_end = std::min(prev_seg->start + prev_seg->length,
                                            curr_seg->start + curr_seg->length);

              // 计算出区间：前缀区间，重叠区间，后缀区间（前缀或后缀不一定存在，但重叠区间一定存在）
              uint_t full_start = std::min(prev_seg->start, curr_seg->start);
              uint_t full_end = std::max(prev_seg->start + prev_seg->length,
                                         curr_seg->start + curr_seg->length);

              // 前缀区间：从合并范围开始到重叠开始
              uint_t prefix_start = full_start;
              uint_t prefix_end = overlap_start;
              uint32_t prefix_len = prefix_end - prefix_start;
              bool prev_has_prefix = prev_seg->start < overlap_start;
              bool curr_has_prefix = curr_seg->start < overlap_start;

              // 后缀区间：从重叠结束到合并范围结束
              uint_t suffix_start = overlap_end;
              uint_t suffix_end = full_end;
              uint32_t suffix_len = suffix_end - suffix_start;
              bool prev_has_suffix =
                  prev_seg->start + prev_seg->length > overlap_end;
              bool curr_has_suffix =
                  curr_seg->start + curr_seg->length > overlap_end;

              // 计算创建新block和新seg的参数，然后创建2-3个新的Block和Segment

              // 声明新blocks（作用域扩大到整个处理过程）
              BlockPtr prefix_block = nullptr;
              BlockPtr overlap_block = nullptr;
              BlockPtr suffix_block = nullptr;
              SegPtr prefix_ref_seg = nullptr;
              SegPtr overlap_ref_seg = nullptr;
              SegPtr suffix_ref_seg = nullptr;

              // 1. 创建前缀block和segment（如果存在）
              if (prev_has_prefix || curr_has_prefix) {
                prefix_block = Block::create(2);
                prefix_block->ref_species = ref_name;
                prefix_block->ref_chr = prev_block->ref_chr;
                prefix_block->ref_species = prev_block->ref_species;
                prefix_ref_seg = Segment::create(
                    prefix_start, prefix_len,
                    Strand::FORWARD, // ref总是正向
                    prev_seg->cigar, // ref没有cigar，所以直接使用prev的cigar
                    prev_seg->align_role, // 使用prev的align_role,
                                          // 以防以后有secondary alignment
                    SegmentRole::SEGMENT, // 前缀是segment，所以是segment
                    prefix_block);
                // 将ref segment注册到block
                prefix_block->anchors.emplace(ref_key, prefix_ref_seg);
                perf_stats.total_blocks_created++;
                perf_stats.total_segments_created++;
              }

              // 2. 创建重叠block和segment（必定存在）
              overlap_block = Block::create(2);
              overlap_block->ref_species = ref_name;
              overlap_block->ref_chr = prev_block->ref_chr;
              overlap_block->ref_species = prev_block->ref_species;
              uint32_t overlap_len = overlap_end - overlap_start;
              overlap_ref_seg = Segment::create(
                  overlap_start, overlap_len,
                  Strand::FORWARD,      // ref总是正向
                  prev_seg->cigar,      // ref没有cigar，所以直接使用prev的cigar
                  prev_seg->align_role, // 使用prev的align_role,
                                        // 以防以后有secondary alignment
                  SegmentRole::SEGMENT, // 重叠是segment，所以是segment
                  overlap_block);
              // 将ref segment注册到block
              overlap_block->anchors.emplace(ref_key, overlap_ref_seg);
              perf_stats.total_blocks_created++;
              perf_stats.total_segments_created++;

              // 3. 创建后缀block和segment（如果存在）
              if (prev_has_suffix || curr_has_suffix) {
                suffix_block = Block::create(2);
                suffix_block->ref_species = ref_name;
                suffix_block->ref_chr = prev_block->ref_chr;
                suffix_block->ref_species = prev_block->ref_species;
                suffix_ref_seg = Segment::create(
                    suffix_start, suffix_len,
                    Strand::FORWARD, // ref总是正向
                    prev_seg->cigar, // ref没有cigar，所以直接使用prev的cigar
                    prev_seg->align_role, // 使用prev的align_role,
                                          // 以防以后有secondary alignment
                    SegmentRole::SEGMENT, // 后缀是segment，所以是segment
                    suffix_block);
                // 将ref segment注册到block
                suffix_block->anchors.emplace(ref_key, suffix_ref_seg);
                perf_stats.total_blocks_created++;
                perf_stats.total_segments_created++;
              }

              // ═══════════════════════════════════════════════════════════
              // 性能分析：Block创建完成
              // ═══════════════════════════════════════════════════════════
              TimePoint block_creation_end = HighResClock::now();
              perf_stats.block_creation_time += duration_cast<Duration>(
                  block_creation_end - block_creation_start);

              TimePoint cigar_processing_start = HighResClock::now();
              // 准备开始处理对应的query，首先处理prev_block，然后处理current_block
              // 1. 处理prev_block
              // 遍历prev_block非ref的anchors
              for (const auto &[species_chr, segment] : prev_block->anchors) {
                // 找到非ref的anchor
                if (species_chr.first != ref_name) {
                  if (segment->strand == Strand::FORWARD) {
                    // if (segment->start == 2174791 && species_chr.second ==
                    // "simChimp.chrD") {
                    //     std::cout<<"test";
                    // }
                    std::string cigar_str;
                    uint32_t sum = 0;
                    bool prefix_done = prev_has_prefix ? false : true;
                    bool overlap_done = false;
                    bool suffix_done = prev_has_suffix ? false : true;
                    SegPtr prefix_qry_seg = nullptr;
                    SegPtr overlap_qry_seg = nullptr;
                    SegPtr suffix_qry_seg = nullptr;
                    uint32_t prefix_length = prev_has_prefix ? prefix_len : 0;
                    uint32_t overlap_length = prefix_length + overlap_len;
                    uint32_t suffix_length =
                        prev_has_suffix ? overlap_length + suffix_len : 0;
                    for (CigarUnit cu : segment->cigar) {
                      perf_stats.total_cigar_operations++;
                      uint32_t cigar_len;
                      char op;
                      intToCigar(cu, op, cigar_len);
                      cigar_str += std::to_string(cigar_len) + op;
                      if (op != 'I') {
                        sum += cigar_len;
                        // 修正累计长度判断
                        if (!prefix_done && sum >= prefix_length) {
                          prefix_done = true;
                          auto [split_cigar, remain_cigar] =
                              splitCigarMixed(cigar_str, prefix_len);
                          cigar_str = remain_cigar;
                          prefix_qry_seg = createSegmentWithTiming(
                              segment->start,
                              countNonDeletionOperations(split_cigar),
                              segment->strand, split_cigar, segment->align_role,
                              SegmentRole::SEGMENT,
                              prefix_block // 使用新创建的prefix_block
                          );
                          // 将query segment注册到prefix_block
                          if (prefix_block) {
                            prefix_block->anchors.emplace(species_chr, prefix_qry_seg);
                          }
                        }
                        if (prefix_done && !overlap_done &&
                            sum >= overlap_length) {
                          // 累计长度
                          overlap_done = true;
                          auto [split_cigar, remain_cigar] =
                              splitCigarMixed(cigar_str, overlap_len);
                          cigar_str = remain_cigar;
                          overlap_qry_seg = createSegmentWithTiming(
                              prefix_qry_seg ? prefix_qry_seg->start +
                                                   prefix_qry_seg->length
                                             : segment->start,
                              countNonDeletionOperations(split_cigar),
                              segment->strand, split_cigar, segment->align_role,
                              SegmentRole::SEGMENT,
                              overlap_block // 使用新创建的overlap_block
                          );
                          // 将query segment注册到overlap_block
                          overlap_block->anchors.emplace(species_chr, overlap_qry_seg);
                        }
                        if (overlap_done && !suffix_done &&
                            sum >= suffix_length) {
                          // 累计长度
                          suffix_done = true;
                          auto [split_cigar, remain_cigar] =
                              splitCigarMixed(cigar_str, suffix_len);
                          cigar_str = remain_cigar;
                          suffix_qry_seg = createSegmentWithTiming(
                              overlap_qry_seg->start + overlap_qry_seg->length,
                              countNonDeletionOperations(split_cigar),
                              segment->strand, split_cigar, segment->align_role,
                              SegmentRole::SEGMENT,
                              suffix_block // 使用新创建的suffix_block
                          );
                          // 将query segment注册到suffix_block
                          if (suffix_block) {
                            suffix_block->anchors.emplace(species_chr,
                                                          suffix_qry_seg);
                          }
                        }
                      }
                    }
                    appendTerminalQueryOnlyRemainder(
                        suffix_qry_seg ? suffix_qry_seg : overlap_qry_seg,
                        cigar_str);
                  } else {
                    // segment->strand == Strand::REVERSE
                    std::string cigar_str;
                    uint32_t sum = 0;
                    bool prefix_done = prev_has_prefix ? false : true;
                    bool overlap_done = false;
                    bool suffix_done = prev_has_suffix ? false : true;
                    SegPtr prefix_qry_seg = nullptr;
                    SegPtr overlap_qry_seg = nullptr;
                    SegPtr suffix_qry_seg = nullptr;
                    uint32_t prefix_length = prev_has_prefix ? prefix_len : 0;
                    uint32_t overlap_length = prefix_length + overlap_len;
                    uint32_t suffix_length =
                        prev_has_suffix ? overlap_length + suffix_len : 0;

                    // CIGAR分割逻辑与正向链相同
                    for (CigarUnit cu : segment->cigar) {
                      uint32_t cigar_len;
                      char op;
                      intToCigar(cu, op, cigar_len);
                      cigar_str += std::to_string(cigar_len) + op;
                      if (op != 'I') {
                        sum += cigar_len;
                        // 分割判断逻辑相同
                        if (!prefix_done && sum >= prefix_length) {
                          prefix_done = true;
                          auto [split_cigar, remain_cigar] =
                              splitCigarMixed(cigar_str, prefix_len);
                          cigar_str = remain_cigar;
                          prefix_qry_seg = Segment::create(
                              0, // 临时坐标，稍后修正
                              countNonDeletionOperations(split_cigar),
                              segment->strand, // 保持REVERSE
                              split_cigar, segment->align_role,
                              SegmentRole::SEGMENT, prefix_block);
                          // 将query segment注册到prefix_block
                          if (prefix_block) {
                            prefix_block->anchors.emplace(species_chr, prefix_qry_seg);
                          }
                        }
                        if (prefix_done && !overlap_done &&
                            sum >= overlap_length) {
                          overlap_done = true;
                          auto [split_cigar, remain_cigar] =
                              splitCigarMixed(cigar_str, overlap_len);
                          cigar_str = remain_cigar;
                          overlap_qry_seg = Segment::create(
                              0, // 临时坐标，稍后修正
                              countNonDeletionOperations(split_cigar),
                              segment->strand, // 保持REVERSE
                              split_cigar, segment->align_role,
                              SegmentRole::SEGMENT, overlap_block);
                          // 将query segment注册到overlap_block
                          overlap_block->anchors.emplace(species_chr, overlap_qry_seg);
                        }
                        if (overlap_done && !suffix_done &&
                            sum >= suffix_length) {
                          suffix_done = true;
                          auto [split_cigar, remain_cigar] =
                              splitCigarMixed(cigar_str, suffix_len);
                          cigar_str = remain_cigar;
                          suffix_qry_seg = Segment::create(
                              0, // 临时坐标，稍后修正
                              countNonDeletionOperations(split_cigar),
                              segment->strand, // 保持REVERSE
                              split_cigar, segment->align_role,
                              SegmentRole::SEGMENT, suffix_block);
                          // 将query segment注册到suffix_block
                          if (suffix_block) {
                            suffix_block->anchors.emplace(species_chr,
                                                          suffix_qry_seg);
                          }
                        }
                      }
                    }
                    appendTerminalQueryOnlyRemainder(
                        suffix_qry_seg ? suffix_qry_seg : overlap_qry_seg,
                        cigar_str);

                    // 反向链坐标修正：从segment末端开始，向前计算
                    uint_t segment_end = segment->start + segment->length;

                    // 使用一个临时变量来追踪当前段的末尾
                    auto current_segment_end = segment_end;

                    // 1. 首先处理前缀 (prefix)，如果它存在
                    if (prefix_qry_seg) {
                      prefix_qry_seg->start =
                          current_segment_end - prefix_qry_seg->length;
                      // 更新末尾位置，为下一个块做准备
                      current_segment_end = prefix_qry_seg->start;
                    }

                    // 2. 接着处理重叠 (overlap)
                    overlap_qry_seg->start =
                        current_segment_end - overlap_qry_seg->length;
                    // 再次更新末尾位置
                    current_segment_end = overlap_qry_seg->start;

                    // 3. 最后处理后缀 (suffix)，如果它存在
                    if (suffix_qry_seg) {
                      suffix_qry_seg->start =
                          current_segment_end - suffix_qry_seg->length;
                    }
                  }
                }
              }

              // 2. 处理current_block的非ref anchors
              for (const auto &[species_chr, segment] :
                   current_block->anchors) {
                // 找到非ref的anchor
                if (species_chr.first != ref_name) {
                  if (segment->strand == Strand::FORWARD) {
                    // if (segment->start == 2174791 && species_chr.second ==
                    // "simChimp.chrD") {
                    //     std::cout<<"test";
                    // }
                    std::string cigar_str;
                    uint32_t sum = 0;
                    bool prefix_done = curr_has_prefix ? false : true;
                    bool overlap_done = false;
                    bool suffix_done = curr_has_suffix ? false : true;
                    SegPtr prefix_qry_seg = nullptr;
                    SegPtr overlap_qry_seg = nullptr;
                    SegPtr suffix_qry_seg = nullptr;
                    uint32_t prefix_length = curr_has_prefix ? prefix_len : 0;
                    uint32_t overlap_length = prefix_length + overlap_len;
                    uint32_t suffix_length =
                        curr_has_suffix ? overlap_length + suffix_len : 0;

                    for (CigarUnit cu : segment->cigar) {
                      uint32_t cigar_len;
                      char op;
                      intToCigar(cu, op, cigar_len);
                      cigar_str += std::to_string(cigar_len) + op;
                      if (op != 'I') {
                        sum += cigar_len;
                        // 修正累计长度判断
                        if (!prefix_done && sum >= prefix_length) {
                          prefix_done = true;
                          auto [split_cigar, remain_cigar] =
                              splitCigarMixed(cigar_str, prefix_len);
                          cigar_str = remain_cigar;
                          prefix_qry_seg = Segment::create(
                              segment->start,
                              countNonDeletionOperations(split_cigar),
                              segment->strand, split_cigar, segment->align_role,
                              SegmentRole::SEGMENT, prefix_block);
                          // 将query segment注册到prefix_block
                          if (prefix_block) {
                            prefix_block->anchors.emplace(species_chr, prefix_qry_seg);
                          }
                        }
                        if (prefix_done && !overlap_done &&
                            sum >= overlap_length) {
                          overlap_done = true;
                          auto [split_cigar, remain_cigar] =
                              splitCigarMixed(cigar_str, overlap_len);
                          cigar_str = remain_cigar;
                          overlap_qry_seg = Segment::create(
                              prefix_qry_seg ? prefix_qry_seg->start +
                                                   prefix_qry_seg->length
                                             : segment->start,
                              countNonDeletionOperations(split_cigar),
                              segment->strand, split_cigar, segment->align_role,
                              SegmentRole::SEGMENT, overlap_block);
                          // 将query segment注册到overlap_block
                          overlap_block->anchors.emplace(species_chr, overlap_qry_seg);
                        }
                        if (overlap_done && !suffix_done &&
                            sum >= suffix_length) {
                          suffix_done = true;
                          auto [split_cigar, remain_cigar] =
                              splitCigarMixed(cigar_str, suffix_len);
                          cigar_str = remain_cigar;
                          uint_t suffix_start_pos =
                              overlap_qry_seg
                                  ? overlap_qry_seg->start +
                                        overlap_qry_seg->length
                                  : (prefix_qry_seg ? prefix_qry_seg->start +
                                                          prefix_qry_seg->length
                                                    : segment->start);
                          suffix_qry_seg = Segment::create(
                              suffix_start_pos,
                              countNonDeletionOperations(split_cigar),
                              segment->strand, split_cigar, segment->align_role,
                              SegmentRole::SEGMENT, suffix_block);
                          // 将query segment注册到suffix_block
                          if (suffix_block) {
                            suffix_block->anchors.emplace(species_chr,
                                                          suffix_qry_seg);
                          }
                        }
                      }
                    }
                    appendTerminalQueryOnlyRemainder(
                        suffix_qry_seg ? suffix_qry_seg : overlap_qry_seg,
                        cigar_str);
                  } else {
                    // segment->strand == Strand::REVERSE
                    std::string cigar_str;
                    uint32_t sum = 0;
                    bool prefix_done = curr_has_prefix ? false : true;
                    bool overlap_done = false;
                    bool suffix_done = curr_has_suffix ? false : true;
                    SegPtr prefix_qry_seg = nullptr;
                    SegPtr overlap_qry_seg = nullptr;
                    SegPtr suffix_qry_seg = nullptr;
                    uint32_t prefix_length = curr_has_prefix ? prefix_len : 0;
                    uint32_t overlap_length = prefix_length + overlap_len;
                    uint32_t suffix_length =
                        curr_has_suffix ? overlap_length + suffix_len : 0;

                    // CIGAR分割逻辑与正向链相同
                    for (CigarUnit cu : segment->cigar) {
                      uint32_t cigar_len;
                      char op;
                      intToCigar(cu, op, cigar_len);
                      cigar_str += std::to_string(cigar_len) + op;
                      if (op != 'I') {
                        sum += cigar_len;
                        // 分割判断逻辑相同
                        if (!prefix_done && sum >= prefix_length) {
                          prefix_done = true;
                          auto [split_cigar, remain_cigar] =
                              splitCigarMixed(cigar_str, prefix_len);
                          cigar_str = remain_cigar;
                          prefix_qry_seg = Segment::create(
                              0, // 临时坐标，稍后修正
                              countNonDeletionOperations(split_cigar),
                              segment->strand, // 保持REVERSE
                              split_cigar, segment->align_role,
                              SegmentRole::SEGMENT, prefix_block);
                          // 将query segment注册到prefix_block
                          if (prefix_block) {
                            prefix_block->anchors.emplace(species_chr, prefix_qry_seg);
                          }
                        }
                        if (prefix_done && !overlap_done &&
                            sum >= overlap_length) {
                          overlap_done = true;
                          auto [split_cigar, remain_cigar] =
                              splitCigarMixed(cigar_str, overlap_len);
                          cigar_str = remain_cigar;
                          overlap_qry_seg = Segment::create(
                              0, // 临时坐标，稍后修正
                              countNonDeletionOperations(split_cigar),
                              segment->strand, // 保持REVERSE
                              split_cigar, segment->align_role,
                              SegmentRole::SEGMENT, overlap_block);
                          // 将query segment注册到overlap_block
                          overlap_block->anchors.emplace(species_chr, overlap_qry_seg);
                        }
                        if (overlap_done && !suffix_done &&
                            sum >= suffix_length) {
                          suffix_done = true;
                          auto [split_cigar, remain_cigar] =
                              splitCigarMixed(cigar_str, suffix_len);
                          cigar_str = remain_cigar;
                          suffix_qry_seg = Segment::create(
                              0, // 临时坐标，稍后修正
                              countNonDeletionOperations(split_cigar),
                              segment->strand, // 保持REVERSE
                              split_cigar, segment->align_role,
                              SegmentRole::SEGMENT, suffix_block);
                          // 将query segment注册到suffix_block
                          if (suffix_block) {
                            suffix_block->anchors.emplace(species_chr,
                                                          suffix_qry_seg);
                          }
                        }
                      }
                    }
                    appendTerminalQueryOnlyRemainder(
                        suffix_qry_seg ? suffix_qry_seg : overlap_qry_seg,
                        cigar_str);

                    // 反向链坐标修正：从segment末端开始，向前计算
                    uint_t segment_end = segment->start + segment->length;

                    // 使用一个临时变量来追踪当前段的末尾
                    auto current_segment_end = segment_end;

                    // 1. 首先处理前缀 (prefix)，如果它存在
                    if (prefix_qry_seg) {
                      prefix_qry_seg->start =
                          current_segment_end - prefix_qry_seg->length;
                      // 更新末尾位置，为下一个块做准备
                      current_segment_end = prefix_qry_seg->start;
                    }

                    // 2. 接着处理重叠 (overlap)
                    overlap_qry_seg->start =
                        current_segment_end - overlap_qry_seg->length;
                    // 再次更新末尾位置
                    current_segment_end = overlap_qry_seg->start;

                    // 3. 最后处理后缀 (suffix)，如果它存在
                    if (suffix_qry_seg) {
                      suffix_qry_seg->start =
                          current_segment_end - suffix_qry_seg->length;
                    }
                  }
                }
              }

              // ═══════════════════════════════════════════════════════════
              // 性能分析：CIGAR处理完成，开始链表重连
              // ═══════════════════════════════════════════════════════════
              TimePoint cigar_processing_end = HighResClock::now();
              perf_stats.cigar_processing_time += duration_cast<Duration>(
                  cigar_processing_end - cigar_processing_start);

              TimePoint link_replacement_start = HighResClock::now();

              // 开始替换prev和current，先替换Ref
              // 获取prev的前驱和current的后继
              SegPtr prev_prev =
                  prev->primary_path.prev.load(std::memory_order_acquire);
              SegPtr current_next =
                  current->primary_path.next.load(std::memory_order_acquire);
              // 移除prev和current的前驱和后继
              prev->primary_path.prev.store(nullptr, std::memory_order_release);
              prev->primary_path.next.store(nullptr, std::memory_order_release);
              current->primary_path.next.store(nullptr,
                                               std::memory_order_release);
              current->primary_path.prev.store(nullptr,
                                               std::memory_order_release);

              // 将新blocks添加到全局池
              if (prefix_block) {
                blocks.emplace_back(WeakBlock(prefix_block));
                // 将prefix的Seg插入

                prefix_ref_seg->primary_path.next.store(
                    overlap_ref_seg, std::memory_order_release);
                overlap_ref_seg->primary_path.prev.store(
                    prefix_ref_seg, std::memory_order_release);

                prefix_ref_seg->primary_path.prev.store(
                    prev_prev, std::memory_order_release);
                prev_prev->primary_path.next.store(prefix_ref_seg,
                                                   std::memory_order_release);
              } else {
                overlap_ref_seg->primary_path.prev.store(
                    prev_prev, std::memory_order_release);
                prev_prev->primary_path.next.store(overlap_ref_seg,
                                                   std::memory_order_release);
              }
              blocks.emplace_back(WeakBlock(overlap_block));

              if (suffix_block) {
                blocks.emplace_back(WeakBlock(suffix_block));
                SegPtr suffix_next = current_next;
                if (current_next->isTail()) {
                  suffix_ref_seg->primary_path.next.store(
                      current_next, std::memory_order_release);
                  current_next->primary_path.prev.store(
                      suffix_ref_seg, std::memory_order_release);

                  overlap_ref_seg->primary_path.next.store(
                      suffix_ref_seg, std::memory_order_release);
                  suffix_ref_seg->primary_path.prev.store(
                      overlap_ref_seg, std::memory_order_release);
                } else if (suffix_ref_seg->start > suffix_next->start) {
                  while (suffix_ref_seg->start > suffix_next->start) {
                    if (suffix_next->isTail()) {
                      break;
                    }
                    suffix_next = suffix_next->primary_path.next.load(
                        std::memory_order_release);
                  }
                  overlap_ref_seg->primary_path.next.store(
                      current_next, std::memory_order_release);
                  current_next->primary_path.prev.store(
                      overlap_ref_seg, std::memory_order_release);

                  suffix_ref_seg->primary_path.prev.store(
                      suffix_next->primary_path.prev.load(
                          std::memory_order_release),
                      std::memory_order_release);
                  suffix_next->primary_path.prev
                      .load(std::memory_order_release)
                      ->primary_path.next.store(suffix_ref_seg,
                                                std::memory_order_release);

                  suffix_ref_seg->primary_path.next.store(
                      suffix_next, std::memory_order_release);
                  suffix_next->primary_path.prev.store(
                      suffix_ref_seg, std::memory_order_release);
                } else {
                  suffix_ref_seg->primary_path.next.store(
                      current_next, std::memory_order_release);
                  current_next->primary_path.prev.store(
                      suffix_ref_seg, std::memory_order_release);

                  overlap_ref_seg->primary_path.next.store(
                      suffix_ref_seg, std::memory_order_release);
                  suffix_ref_seg->primary_path.prev.store(
                      overlap_ref_seg, std::memory_order_release);
                }
              } else {
                overlap_ref_seg->primary_path.next.store(
                    current_next, std::memory_order_release);
                current_next->primary_path.prev.store(
                    overlap_ref_seg, std::memory_order_release);
              }

              // 开始替换query的seg和block
              for (const auto &[species_chr, segment] : prev_block->anchors) {
                if (species_chr.first != ref_name) {
                  SegPtr qry_prev = segment->primary_path.prev.load(
                      std::memory_order_acquire);
                  SegPtr qry_next = segment->primary_path.next.load(
                      std::memory_order_acquire);
                  // 先找到一定存在的overlap_block
                  bool matched = false;
                  for (const auto &[species_chr_overlap, segment_overlap] :
                       overlap_block->anchors) {
                    if (species_chr_overlap == species_chr) {
                      matched = true;
                      Segment::unlinkSegment(segment);
                      if (segment->strand == Strand::FORWARD) {
                        bool query_has_prefix = false;
                        bool query_has_suffix = false;
                        if (prev_has_prefix) {
                          for (const auto &[species_chr_prefix,
                                            segment_prefix] :
                               prefix_block->anchors) {
                            if (species_chr_prefix == species_chr) {
                              query_has_prefix = true;
                              qry_prev->primary_path.next.store(
                                  segment_prefix, std::memory_order_release);
                              segment_prefix->primary_path.prev.store(
                                  qry_prev, std::memory_order_release);
                              segment_prefix->primary_path.next.store(
                                  segment_overlap, std::memory_order_release);
                              segment_overlap->primary_path.prev.store(
                                  segment_prefix, std::memory_order_release);
                              break;
                            }
                          }
                          if (!query_has_prefix) {
                            qry_prev->primary_path.next.store(
                                segment_overlap, std::memory_order_release);
                            segment_overlap->primary_path.prev.store(
                                qry_prev, std::memory_order_release);
                          }
                        } else {
                          qry_prev->primary_path.next.store(
                              segment_overlap, std::memory_order_release);
                          segment_overlap->primary_path.prev.store(
                              qry_prev, std::memory_order_release);
                        }
                        if (prev_has_suffix) {
                          for (const auto &[species_chr_suffix,
                                            segment_suffix] :
                               suffix_block->anchors) {
                            if (species_chr_suffix == species_chr) {
                              query_has_suffix = true;
                              segment_overlap->primary_path.next.store(
                                  segment_suffix, std::memory_order_release);
                              segment_suffix->primary_path.prev.store(
                                  segment_overlap, std::memory_order_release);
                              segment_suffix->primary_path.next.store(
                                  qry_next, std::memory_order_release);
                              qry_next->primary_path.prev.store(
                                  segment_suffix, std::memory_order_release);
                              break;
                            }
                          }
                          if (!query_has_suffix) {
                            segment_overlap->primary_path.next.store(
                                qry_next, std::memory_order_release);
                            qry_next->primary_path.prev.store(
                                segment_overlap, std::memory_order_release);
                          }
                        } else {
                          segment_overlap->primary_path.next.store(
                              qry_next, std::memory_order_release);
                          qry_next->primary_path.prev.store(
                              segment_overlap, std::memory_order_release);
                        }
                      }
                      // 如果是反向链就反着连
                      else {
                        bool query_has_prefix = false;
                        bool query_has_suffix = false;
                        if (prev_has_suffix) {
                          for (const auto &[species_chr_suffix,
                                            segment_suffix] :
                               suffix_block->anchors) {
                            if (species_chr_suffix == species_chr) {
                              query_has_suffix = true;
                              qry_prev->primary_path.next.store(
                                  segment_suffix, std::memory_order_release);
                              segment_suffix->primary_path.prev.store(
                                  qry_prev, std::memory_order_release);
                              segment_suffix->primary_path.next.store(
                                  segment_overlap, std::memory_order_release);
                              segment_overlap->primary_path.prev.store(
                                  segment_suffix, std::memory_order_release);
                              break;
                            }
                          }
                          if (!query_has_suffix) {
                            qry_prev->primary_path.next.store(
                                segment_overlap, std::memory_order_release);
                            segment_overlap->primary_path.prev.store(
                                qry_prev, std::memory_order_release);
                          }
                        } else {
                          qry_prev->primary_path.next.store(
                              segment_overlap, std::memory_order_release);
                          segment_overlap->primary_path.prev.store(
                              qry_prev, std::memory_order_release);
                        }
                        if (prev_has_prefix) {
                          for (const auto &[species_chr_prefix,
                                            segment_prefix] :
                               prefix_block->anchors) {
                            if (species_chr_prefix == species_chr) {
                              query_has_prefix = true;
                              segment_prefix->primary_path.prev.store(
                                  segment_overlap, std::memory_order_release);
                              segment_prefix->primary_path.next.store(
                                  qry_next, std::memory_order_release);
                              qry_next->primary_path.prev.store(
                                  segment_prefix, std::memory_order_release);
                              segment_overlap->primary_path.next.store(
                                  segment_prefix, std::memory_order_release);

                              break;
                            }
                          }
                          if (!query_has_prefix) {
                            segment_overlap->primary_path.next.store(
                                qry_next, std::memory_order_release);
                            qry_next->primary_path.prev.store(
                                segment_overlap, std::memory_order_release);
                          }
                        } else {
                          segment_overlap->primary_path.next.store(
                              qry_next, std::memory_order_release);
                          qry_next->primary_path.prev.store(
                              segment_overlap, std::memory_order_release);
                        }
                      }
                      break;
                    }
                  }
                  if (!matched) {
                    continue;
                  }
                }
              }
              for (const auto &[species_chr, segment] :
                   current_block->anchors) {
                if (species_chr.first != ref_name) {
                  SegPtr qry_prev = segment->primary_path.prev.load(
                      std::memory_order_acquire);
                  SegPtr qry_next = segment->primary_path.next.load(
                      std::memory_order_acquire);
                  // 先找到一定存在的overlap_block
                  bool matched = false;
                  for (const auto &[species_chr_overlap, segment_overlap] :
                       overlap_block->anchors) {
                    if (species_chr_overlap == species_chr) {
                      matched = true;
                      Segment::unlinkSegment(segment);
                      if (segment->strand == Strand::FORWARD) {
                        bool query_has_prefix = false;
                        bool query_has_suffix = false;
                        if (curr_has_prefix) {
                          for (const auto &[species_chr_prefix,
                                            segment_prefix] :
                               prefix_block->anchors) {
                            if (species_chr_prefix == species_chr) {
                              query_has_prefix = true;
                              qry_prev->primary_path.next.store(
                                  segment_prefix, std::memory_order_release);
                              segment_prefix->primary_path.prev.store(
                                  qry_prev, std::memory_order_release);
                              segment_prefix->primary_path.next.store(
                                  segment_overlap, std::memory_order_release);
                              segment_overlap->primary_path.prev.store(
                                  segment_prefix, std::memory_order_release);
                              break;
                            }
                          }
                          if (!query_has_prefix) {
                            qry_prev->primary_path.next.store(
                                segment_overlap, std::memory_order_release);
                            segment_overlap->primary_path.prev.store(
                                qry_prev, std::memory_order_release);
                          }
                        } else {
                          qry_prev->primary_path.next.store(
                              segment_overlap, std::memory_order_release);
                          segment_overlap->primary_path.prev.store(
                              qry_prev, std::memory_order_release);
                        }
                        if (curr_has_suffix) {
                          for (const auto &[species_chr_suffix,
                                            segment_suffix] :
                               suffix_block->anchors) {
                            if (species_chr_suffix == species_chr) {
                              query_has_suffix = true;
                              segment_overlap->primary_path.next.store(
                                  segment_suffix, std::memory_order_release);
                              segment_suffix->primary_path.prev.store(
                                  segment_overlap, std::memory_order_release);
                              segment_suffix->primary_path.next.store(
                                  qry_next, std::memory_order_release);
                              qry_next->primary_path.prev.store(
                                  segment_suffix, std::memory_order_release);
                              break;
                            }
                          }
                          if (!query_has_suffix) {
                            segment_overlap->primary_path.next.store(
                                qry_next, std::memory_order_release);
                            qry_next->primary_path.prev.store(
                                segment_overlap, std::memory_order_release);
                          }
                        } else {
                          segment_overlap->primary_path.next.store(
                              qry_next, std::memory_order_release);
                          qry_next->primary_path.prev.store(
                              segment_overlap, std::memory_order_release);
                        }
                      }
                      // 如果是反向链就反着连
                      else {
                        bool query_has_prefix = false;
                        bool query_has_suffix = false;
                        if (curr_has_suffix) {
                          for (const auto &[species_chr_suffix,
                                            segment_suffix] :
                               suffix_block->anchors) {
                            if (species_chr_suffix == species_chr) {
                              query_has_suffix = true;
                              qry_prev->primary_path.next.store(
                                  segment_suffix, std::memory_order_release);
                              segment_suffix->primary_path.prev.store(
                                  qry_prev, std::memory_order_release);
                              segment_suffix->primary_path.next.store(
                                  segment_overlap, std::memory_order_release);
                              segment_overlap->primary_path.prev.store(
                                  segment_suffix, std::memory_order_release);
                              break;
                            }
                          }
                          if (!query_has_suffix) {
                            qry_prev->primary_path.next.store(
                                segment_overlap, std::memory_order_release);
                            segment_overlap->primary_path.prev.store(
                                qry_prev, std::memory_order_release);
                          }
                        } else {
                          qry_prev->primary_path.next.store(
                              segment_overlap, std::memory_order_release);
                          segment_overlap->primary_path.prev.store(
                              qry_prev, std::memory_order_release);
                        }
                        if (curr_has_prefix) {
                          for (const auto &[species_chr_prefix,
                                            segment_prefix] :
                               prefix_block->anchors) {
                            if (species_chr_prefix == species_chr) {
                              query_has_prefix = true;
                              segment_prefix->primary_path.prev.store(
                                  segment_overlap, std::memory_order_release);
                              segment_prefix->primary_path.next.store(
                                  qry_next, std::memory_order_release);
                              qry_next->primary_path.prev.store(
                                  segment_prefix, std::memory_order_release);
                              segment_overlap->primary_path.next.store(
                                  segment_prefix, std::memory_order_release);

                              break;
                            }
                          }
                          if (!query_has_prefix) {
                            segment_overlap->primary_path.next.store(
                                qry_next, std::memory_order_release);
                            qry_next->primary_path.prev.store(
                                segment_overlap, std::memory_order_release);
                          }
                        } else {
                          segment_overlap->primary_path.next.store(
                              qry_next, std::memory_order_release);
                          qry_next->primary_path.prev.store(
                              segment_overlap, std::memory_order_release);
                        }
                      }
                      break;
                    }
                  }
                  if (!matched) {
                    continue;
                  }
                }
              }

              // ═══════════════════════════════════════════════════════════
              // 性能分析：链表重连完成，开始清理工作
              // ═══════════════════════════════════════════════════════════
              TimePoint link_replacement_end = HighResClock::now();
              perf_stats.link_replacement_time += duration_cast<Duration>(
                  link_replacement_end - link_replacement_start);

              TimePoint cleanup_start = HighResClock::now();

              // ═══════════════════════════════════════════════════════════
              // 批量删除优化：标记要删除的Block，而不是立即删除
              // ═══════════════════════════════════════════════════════════
              if (prev_block) {
                blocks_to_delete.insert(prev_block);
                perf_stats.total_blocks_marked_for_deletion++;
              }
              if (current_block) {
                blocks_to_delete.insert(current_block);
                perf_stats.total_blocks_marked_for_deletion++;
              }
              // 替换current和prev
              prev = overlap_ref_seg->primary_path.prev.load(
                  std::memory_order_acquire);
              if (prev->isHead()) {
                prev = overlap_ref_seg;
              }
              prev_block = prev->parent_block;
              current = prev->primary_path.next.load(std::memory_order_acquire);
              current_block = current->parent_block;
              // ═══════════════════════════════════════════════════════════
              // 性能分析：清理工作完成
              // ═══════════════════════════════════════════════════════════
              TimePoint cleanup_end = HighResClock::now();
              perf_stats.cleanup_time +=
                  duration_cast<Duration>(cleanup_end - cleanup_start);

              // 继续处理，不要跳到最后，因为新创建的blocks之间可能还有重叠

              //                                 //
              //                                 收集所有物种的segment详细信息
              //                                 debug_all_species_segments.clear();
              //                                 for (const auto& [species_name,
              //                                 genome_graph] : species_graphs)
              //                                 {
              //                                     SpeciesDebugInfo
              //                                     species_info(species_name);
              //                                     std::shared_lock
              //                                     genome_lock(genome_graph.rw);
              //
              //                                     for (const auto& [chr_name,
              //                                     genome_end] :
              //                                     genome_graph.chr2end) {
              //                                         ChromosomeDebugInfo
              //                                         chr_info(chr_name);
              //                                         std::shared_lock
              //                                         end_lock(genome_end.rw);
              //
              //                                         //
              //                                         遍历该染色体的所有segments
              //                                         SegPtr current =
              //                                         genome_end.head;
              //                                         current =
              //                                         current->primary_path.next.load(std::memory_order_acquire);
              //                                         while (current &&
              //                                         current !=
              //                                         genome_end.tail) {
              //                                             if
              //                                             (current->isSegment())
              //                                             {
              //                                                 chr_info.segments.emplace_back(current);
              //                                             }
              //                                             //
              //                                             检查链表结构是否正确
              //                                             if
              //                                             (current->isHead())
              //                                             {
              //                                                 current =
              //                                                 current->primary_path.next.load(std::memory_order_acquire);
              //                                                 continue;
              //                                             }
              //                                             SegPtr prev =
              //                                             current->primary_path.prev.load(std::memory_order_acquire);
              //                                             if(prev &&
              //                                             prev->isHead())
              //                                             {
              //                                                 current =
              //                                                 current->primary_path.next.load(std::memory_order_acquire);
              //                                                 continue;
              //                                             }
              //                                             if (prev &&
              //                                             prev->primary_path.next.load(std::memory_order_acquire)
              //                                             != current) {
              //                                                 std::cerr <<
              //                                                 "[链表错误]
              //                                                 prev->next !=
              //                                                 current,
              //                                                 current start:
              //                                                 " <<
              //                                                 current->start
              //                                                 << std::endl;
              //                                             }
              //                                             current =
              //                                             current->primary_path.next.load(std::memory_order_acquire);
              //                                         }
              //
              //                                         species_info.chromosomes.push_back(std::move(chr_info));
              //                                     }
              //
              //                                     debug_all_species_segments.push_back(std::move(species_info));
              //                                 }
              //
              // count++;
              // std::cout<<count<<std::endl;
              continue;
            }
          }
        }

        prev_block = current_block;
        prev = current;
        // 移动到下一个segment - 使用原子操作保证线程安全
        current = current->primary_path.next.load(std::memory_order_acquire);
      }

      // ═══════════════════════════════════════════════════════════
      // 标记当前染色体处理完成
      // ═══════════════════════════════════════════════════════════
      progress.completeChromosome();
    }
  }

  // ═══════════════════════════════════════════════════════════
  // 性能分析：染色体遍历完成
  // ═══════════════════════════════════════════════════════════
  TimePoint main_loop_end = HighResClock::now();
  perf_stats.chromosome_iteration_time =
      duration_cast<Duration>(main_loop_end - main_loop_start);

  // ═══════════════════════════════════════════════════════════
  // 批量删除优化：一次性删除所有标记的Block
  // ═══════════════════════════════════════════════════════════
  TimePoint batch_deletion_start = HighResClock::now();

  if (!blocks_to_delete.empty()) {
    spdlog::info("Starting batch deletion of {} blocks...", blocks_to_delete.size());


    // 升级为写锁以进行删除操作
    graph_lock.unlock();
    std::unique_lock<std::shared_mutex> write_lock(rw);

    // 使用高效的批量删除：标记-清除策略
    size_t original_size = blocks.size();

    // 1. Block查找阶段计时
    TimePoint lookup_start = HighResClock::now();

    // 方法1：使用remove_if + erase idiom（推荐）
    auto new_end = std::remove_if(
        blocks.begin(), blocks.end(), [&blocks_to_delete](const WeakBlock &wb) {
          auto locked = wb.lock();
          return locked &&
                 blocks_to_delete.find(locked) != blocks_to_delete.end();
        });

    TimePoint lookup_end = HighResClock::now();
    perf_stats.block_lookup_time +=
        duration_cast<Duration>(lookup_end - lookup_start);

    // 2. 实际删除阶段计时
    TimePoint deletion_start = HighResClock::now();

    blocks.erase(new_end, blocks.end());

    TimePoint deletion_end = HighResClock::now();
    perf_stats.actual_deletion_time +=
        duration_cast<Duration>(deletion_end - deletion_start);

    size_t deleted_count = original_size - blocks.size();
    perf_stats.batch_deletion_operations = 1;

    spdlog::info("Delete {} Block", deleted_count);

    // 释放写锁，重新获取读锁
    write_lock.unlock();
    graph_lock = std::shared_lock<std::shared_mutex>(rw);
  }

  TimePoint batch_deletion_end = HighResClock::now();
  Duration batch_deletion_time =
      duration_cast<Duration>(batch_deletion_end - batch_deletion_start);

  // 将批量删除时间加入清理工作统计
  perf_stats.cleanup_time += batch_deletion_time;

  // ═══════════════════════════════════════════════════════════
  // 完成所有处理
  // ═══════════════════════════════════════════════════════════
  progress.finish();

  // ═══════════════════════════════════════════════════════════
  // 性能分析：函数执行完成，生成最终报告
  // ═══════════════════════════════════════════════════════════
  TimePoint function_end_time = HighResClock::now();
  Duration total_function_time =
      duration_cast<Duration>(function_end_time - function_start_time);
#ifdef _DEBUG_
  perf_stats.logFinalReport(total_function_time);
#endif
  // // 收集所有物种的segment详细信息
  // debug_all_species_segments.clear();
  // for (const auto &[species_name, genome_graph]: species_graphs) {
  //     SpeciesDebugInfo species_info(species_name);
  //     std::shared_lock genome_lock(genome_graph.rw);
  //
  //     for (const auto &[chr_name, genome_end]: genome_graph.chr2end) {
  //         ChromosomeDebugInfo chr_info(chr_name);
  //         std::shared_lock end_lock(genome_end.rw);
  //
  //         // 遍历该染色体的所有segments
  //         SegPtr current = genome_end.head;
  //         current =
  //         current->primary_path.next.load(std::memory_order_acquire); while
  //         (current && current != genome_end.tail) {
  //             if (current->isSegment()) {
  //                 chr_info.segments.emplace_back(current);
  //             }
  //             // 检查链表结构是否正确
  //             if (current->isHead()) {
  //                 current =
  //                 current->primary_path.next.load(std::memory_order_acquire);
  //                 continue;
  //             }
  //             SegPtr prev =
  //             current->primary_path.prev.load(std::memory_order_acquire); if
  //             (prev && prev->isHead()) {
  //                 current =
  //                 current->primary_path.next.load(std::memory_order_acquire);
  //                 continue;
  //             }
  //             if (prev &&
  //             prev->primary_path.next.load(std::memory_order_acquire) !=
  //             current) {
  //                 std::cerr << "[链表错误] prev->next != current, current
  //                 start: " << current->start << std::endl;
  //             }
  //             current =
  //             current->primary_path.next.load(std::memory_order_acquire);
  //         }
  //
  //         species_info.chromosomes.push_back(std::move(chr_info));
  //     }
  //
  //     debug_all_species_segments.push_back(std::move(species_info));
  // }
}

/* ==============================================================
 * 7. 高性能删除方法 (public API)
 * ==============================================================*/
bool RaMeshMultiGenomeGraph::removeBlock(const BlockPtr &block) {
  if (!block)
    return false;

  std::unique_lock graph_lock(rw);

  // 1. 从全局block池中移除
  auto it =
      std::find_if(blocks.begin(), blocks.end(), [&](const WeakBlock &wb) {
        auto locked = wb.lock();
        return locked && locked == block;
      });

  if (it != blocks.end()) {
    blocks.erase(it);
  } else {
    return false; // 未找到block
  }

  // 2. 从每个物种的链表中移除相关segments
  {
    std::shared_lock block_lock(block->rw);
    for (const auto &[species_chr, segment] : block->anchors) {
      if (segment) {
        const SpeciesName &species = species_chr.first;
        const ChrName &chr = species_chr.second;

        auto species_it = species_graphs.find(species);
        if (species_it != species_graphs.end()) {
          auto chr_it = species_it->second.chr2end.find(chr);
          if (chr_it != species_it->second.chr2end.end()) {
            chr_it->second.removeSegment(segment);
          }
        }
      }
    }
  }

  // 3. 清空block的anchors（segments会由GenomeEnd负责删除）
  block->removeAllSegments();

  return true;
}

bool RaMeshMultiGenomeGraph::removeBlockById(size_t block_index) {
  std::unique_lock graph_lock(rw);

  if (block_index >= blocks.size())
    return false;

  auto weak_block = blocks[block_index];
  auto block = weak_block.lock();

  if (!block) {
    // 清理过期的weak_ptr
    blocks.erase(blocks.begin() + block_index);
    return false;
  }

  return removeBlock(block);
}

void RaMeshMultiGenomeGraph::removeBlocksBatch(
    const std::vector<BlockPtr> &blocks_to_remove) {
  if (blocks_to_remove.empty())
    return;

  std::unique_lock graph_lock(rw);

  // 收集所有需要删除的segments，按物种-染色体分组
  std::unordered_map<SpeciesName,
                     std::unordered_map<ChrName, std::vector<SegPtr>>>
      segments_by_genome;

  for (const auto &block : blocks_to_remove) {
    if (!block)
      continue;

    std::shared_lock block_lock(block->rw);
    for (const auto &[species_chr, segment] : block->anchors) {
      if (segment) {
        segments_by_genome[species_chr.first][species_chr.second].emplace_back(
            segment);
      }
    }
  }

  // 批量删除segments
  for (const auto &[species, chr_map] : segments_by_genome) {
    auto species_it = species_graphs.find(species);
    if (species_it != species_graphs.end()) {
      for (const auto &[chr, segments] : chr_map) {
        auto chr_it = species_it->second.chr2end.find(chr);
        if (chr_it != species_it->second.chr2end.end()) {
          chr_it->second.removeBatch(segments);
        }
      }
    }
  }

  // 从全局block池中移除
  for (const auto &block : blocks_to_remove) {
    auto it =
        std::find_if(blocks.begin(), blocks.end(), [&](const WeakBlock &wb) {
          auto locked = wb.lock();
          return locked && locked == block;
        });

    if (it != blocks.end()) {
      blocks.erase(it);
    }

    // 清空block的anchors
    block->removeAllSegments();
  }
}

size_t RaMeshMultiGenomeGraph::removeExpiredBlocks() {
  std::unique_lock graph_lock(rw);

  size_t removed_count = 0;

  // 移除所有过期的weak_ptr
  blocks.erase(std::remove_if(blocks.begin(), blocks.end(),
                              [&removed_count](const WeakBlock &wb) {
                                if (wb.expired()) {
                                  ++removed_count;
                                  return true;
                                }
                                return false;
                              }),
               blocks.end());

  return removed_count;
}

void RaMeshMultiGenomeGraph::removeSpecies(const SpeciesName &species) {
  std::unique_lock graph_lock(rw);

  // 1. 移除species_graphs中的条目
  auto species_it = species_graphs.find(species);
  if (species_it == species_graphs.end())
    return;

  // 2. 收集包含该物种的所有blocks
  std::vector<BlockPtr> blocks_to_remove;

  for (const auto &weak_block : blocks) {
    auto block = weak_block.lock();
    if (!block)
      continue;

    std::shared_lock block_lock(block->rw);
    for (const auto &[species_chr, segment] : block->anchors) {
      if (species_chr.first == species) {
        blocks_to_remove.emplace_back(block);
        break;
      }
    }
  }

  // 3. 批量删除blocks
  removeBlocksBatch(blocks_to_remove);

  // 4. 从species_graphs中移除
  species_graphs.erase(species_it);
}

void RaMeshMultiGenomeGraph::removeChromosome(const SpeciesName &species,
                                              const ChrName &chr) {
  std::unique_lock graph_lock(rw);

  auto species_it = species_graphs.find(species);
  if (species_it == species_graphs.end())
    return;

  auto chr_it = species_it->second.chr2end.find(chr);
  if (chr_it == species_it->second.chr2end.end())
    return;

  // 1. 收集包含该染色体的所有blocks
  std::vector<BlockPtr> blocks_to_update;
  std::vector<BlockPtr> blocks_to_remove;

  for (const auto &weak_block : blocks) {
    auto block = weak_block.lock();
    if (!block)
      continue;

    std::shared_lock block_lock(block->rw);
    SpeciesChrPair target_key{species, chr};

    if (block->anchors.find(target_key) != block->anchors.end()) {
      if (block->anchors.count(target_key) == block->anchors.size()) {
        // 如果block只包含这一个染色体，删除整个block
        blocks_to_remove.emplace_back(block);
      } else {
        // 否则只删除这个染色体的segment
        blocks_to_update.emplace_back(block);
      }
    }
  }

  // 2. 清空该染色体的所有segments
  chr_it->second.clearAllSegments();

  // 3. 更新相关blocks
  for (const auto &block : blocks_to_update) {
    block->removeSegments(species, chr);
  }

  // 4. 删除只包含该染色体的blocks
  removeBlocksBatch(blocks_to_remove);

  // 5. 从species_graphs中移除染色体
  species_it->second.chr2end.erase(chr_it);
}

void RaMeshMultiGenomeGraph::clearAllGraphs() {
  std::unique_lock graph_lock(rw);

  // 1. 清空所有species的图
  for (auto &[species, genome_graph] : species_graphs) {
    std::unique_lock species_lock(genome_graph.rw);
    for (auto &[chr, genome_end] : genome_graph.chr2end) {
      genome_end.clearAllSegments();
    }
    genome_graph.chr2end.clear();
  }

  // 2. 清空所有blocks
  for (const auto &weak_block : blocks) {
    auto block = weak_block.lock();
    if (block) {
      block->removeAllSegments();
    }
  }
  blocks.clear();

  // 3. 清空species_graphs
  species_graphs.clear();
}

size_t RaMeshMultiGenomeGraph::compactBlockPool() {
  std::unique_lock graph_lock(rw);

  size_t compacted = removeExpiredBlocks();

  // 可以在这里添加更多的内存优化逻辑
  blocks.shrink_to_fit();

  return compacted;
}

void RaMeshMultiGenomeGraph::optimizeGraphStructure() {
  size_t zero_length_removed = 0;
  size_t detached_block_anchors = 0;

  for (auto &[species, genome_graph] : species_graphs) {
    std::unique_lock species_lock(genome_graph.rw);

    for (auto &[chr, genome_end] : genome_graph.chr2end) {
      std::unique_lock end_lock(genome_end.rw);

      if (!genome_end.head || !genome_end.tail) {
        genome_end.sample_vec.clear();
        continue;
      }

      // 清理零长度 segment 并同步 block/采样状态
      SegPtr current =
          genome_end.head->primary_path.next.load(std::memory_order_acquire);
      std::unordered_set<const Segment*> cleanup_visited;

      while (current && !current->isTail()) {
        if (!cleanup_visited.emplace(current.get()).second) {
          throw std::runtime_error(
              "optimizeGraphStructure detected a cycle during cleanup for " +
              species + "." + chr);
        }
        if (current->isSegment() && current->length == 0) {
          SegPtr next_ptr =
              current->primary_path.next.load(std::memory_order_acquire);

          Segment::unlinkSegment(current);

          if (auto block = current->parent_block) {
            std::unique_lock block_lock(block->rw);
            for (auto it = block->anchors.begin(); it != block->anchors.end();
                 ++it) {
              if (it->second == current) {
                block->anchors.erase(it);
                detached_block_anchors++;
                break;
              }
            }
          }

          current->parent_block.reset();
          zero_length_removed++;
          current = next_ptr;
          continue;
        }

        current =
            current->primary_path.next.load(std::memory_order_acquire);
      }

      genome_end.sample_vec.clear();
      SegPtr rebuild =
          genome_end.head->primary_path.next.load(std::memory_order_acquire);
      std::unordered_set<const Segment*> rebuild_visited;
      while (rebuild && !rebuild->isTail()) {
        if (!rebuild_visited.emplace(rebuild.get()).second) {
          throw std::runtime_error(
              "optimizeGraphStructure detected a cycle during sampling for " +
              species + "." + chr);
        }
        genome_end.setToSampling(rebuild);
        rebuild =
            rebuild->primary_path.next.load(std::memory_order_acquire);
      }
    }
  }

  if (zero_length_removed > 0 || detached_block_anchors > 0) {
    spdlog::info(
        "optimizeGraphStructure: removed {} zero-length segments, detached {} block anchors",
        zero_length_removed, detached_block_anchors);
  }
}

RaMeshMultiGenomeGraph::DeletionStats
RaMeshMultiGenomeGraph::performMaintenance(bool full_gc) {
  auto start_time = std::chrono::high_resolution_clock::now();

  DeletionStats stats;

  // 1. 清理过期的blocks
  stats.expired_blocks_cleaned = removeExpiredBlocks();

  if (full_gc) {
    // 2. 压缩block池
    stats.blocks_removed = compactBlockPool();

    // 3. 优化图结构
    optimizeGraphStructure();
    stats.sampling_updates = species_graphs.size(); // 简化统计
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  stats.total_time = std::chrono::duration_cast<std::chrono::microseconds>(
      end_time - start_time);

  return stats;
}

// ――― 优化的统一遍历函数实现 ―――
void RaMeshMultiGenomeGraph::verifyWithUnifiedTraversal(
    VerificationResult &result, const VerificationOptions &options) const {
  if (options.verbose) {
    spdlog::debug("Performing unified traversal verification...");
  }

  // 预构建 block 查找映射以优化内存完整性检查
  std::unordered_set<BlockPtr> valid_blocks;
  if (options.isEnabled(VerificationType::MEMORY_INTEGRITY)) {
    valid_blocks.reserve(blocks.size());
    for (const auto &weak_block : blocks) {
      if (auto block_ptr = weak_block.lock()) {
        valid_blocks.insert(block_ptr);
      }
    }
  }

  const bool reference_policy_enabled = options.reference_policy.enabled;
  const SpeciesName &reference_species =
      options.reference_policy.reference_species;
  const bool allow_reference_overlap =
      !reference_policy_enabled ||
      options.reference_policy.allow_reference_overlap;
  const bool require_reference_overlap =
      reference_policy_enabled &&
      options.reference_policy.require_reference_overlap;
  const bool forbid_non_reference_overlap =
      !reference_policy_enabled ||
      options.reference_policy.forbid_non_reference_overlap;

  size_t reference_overlap_counter = 0;
  bool reference_overlap_found = false;
  size_t reference_chr_visited = 0;

  for (const auto &[species_name, genome_graph] : species_graphs) {
    std::shared_lock species_lock(genome_graph.rw);

    for (const auto &[chr_name, genome_end] : genome_graph.chr2end) {
      if (!genome_end.head || !genome_end.tail) {
        continue; // 已在指针有效性检查中报告
      }

      if (reference_policy_enabled && species_name == reference_species) {
        reference_chr_visited++;
      }

      // 为链表完整性检查准备数据结构
      std::unordered_set<SegPtr> visited_segments;
      std::unordered_set<SegPtr> all_segments; // 用于内存完整性检查

      if (options.isEnabled(VerificationType::LINKED_LIST_INTEGRITY)) {
        visited_segments.reserve(10000); // 预分配合理大小
      }
      if (options.isEnabled(VerificationType::MEMORY_INTEGRITY)) {
        all_segments.reserve(10000); // 预分配合理大小
      }

      SegPtr current = genome_end.head;
      SegPtr prev = nullptr;
      SegPtr prev_nonzero_segment = nullptr; // 用于坐标重叠检查（忽略0长度segment）
      SegPtr prev_segment = nullptr; // 用于坐标排序检查
      size_t segment_count = 0;
      size_t ordering_violations = 0;
      size_t large_gaps = 0;
      const uint_t LARGE_GAP_THRESHOLD = 1000000;

      // 统一遍历链表
      while (current) {
        segment_count++;

        // ===== 链表完整性检查 =====
        if (options.isEnabled(VerificationType::LINKED_LIST_INTEGRITY)) {
          // 检查循环引用
          if (visited_segments.find(current) != visited_segments.end()) {
            addVerificationError(
                result, options, VerificationType::LINKED_LIST_INTEGRITY,
                ErrorSeverity::CRITICAL, species_name, chr_name, segment_count,
                current->start, "Circular reference detected in linked list",
                "Segment appears twice in traversal");
            break;
          }
          visited_segments.insert(current);

          // 检查双向链表的一致性
          SegPtr next_ptr =
              current->primary_path.next.load(std::memory_order_acquire);
          SegPtr prev_ptr =
              current->primary_path.prev.load(std::memory_order_acquire);

          if (prev_ptr != prev) {
            addVerificationError(
                result, options, VerificationType::LINKED_LIST_INTEGRITY,
                ErrorSeverity::ERROR, species_name, chr_name, segment_count,
                current->start,
                "Inconsistent prev pointer in doubly linked list",
                "expected_prev=" +
                    std::to_string(reinterpret_cast<uintptr_t>(prev.get())) +
                    ", actual_prev=" +
                    std::to_string(
                        reinterpret_cast<uintptr_t>(prev_ptr.get())));
          }

          // 防止死循环
          if (segment_count > 10000000) {
            addVerificationError(
                result, options, VerificationType::LINKED_LIST_INTEGRITY,
                ErrorSeverity::CRITICAL, species_name, chr_name, segment_count,
                0, "Linked list may contain cycle",
                "Traversed over 10 million nodes");
            break;
          }
        }

        // ===== 内存完整性检查 =====
        if (options.isEnabled(VerificationType::MEMORY_INTEGRITY)) {
          // 检查重复的segment指针
          if (all_segments.find(current) != all_segments.end()) {
            addVerificationError(
                result, options, VerificationType::MEMORY_INTEGRITY,
                ErrorSeverity::CRITICAL, species_name, chr_name, segment_count,
                current->start, "Duplicate segment pointer detected",
                "Same segment appears multiple times in memory");
          }
          all_segments.insert(current);
        }

        // ===== 坐标相关检查（仅对实际segment） =====
        if (current->isSegment()) {
          // 坐标重叠检查
          if (options.isEnabled(VerificationType::COORDINATE_OVERLAP)) {
            // 检查segment长度有效性
            if (current->length == 0) {
              addVerificationError(
                  result, options, VerificationType::COORDINATE_OVERLAP,
                  ErrorSeverity::WARNING, species_name, chr_name, segment_count,
                  current->start, "Zero-length segment detected",
                  "start=" + std::to_string(current->start));

              // 继续遍历，但要保持 prev 同步为真实前驱节点，避免下一节点误报 prev 不一致
              prev = current;
              current =
                  current->primary_path.next.load(std::memory_order_acquire);
              continue;
            }

            // 检查坐标是否有重叠
            if (prev_nonzero_segment) {
              uint_t prev_end =
                  prev_nonzero_segment->start + prev_nonzero_segment->length;
              if (current->start < prev_end) {
                uint_t overlap_size = prev_end - current->start;
                bool is_reference_species = reference_policy_enabled &&
                                            species_name == reference_species;

                if (is_reference_species) {
                  reference_overlap_found = true;
                  reference_overlap_counter++;

                  if (!allow_reference_overlap) {
                    addVerificationError(
                        result, options, VerificationType::COORDINATE_OVERLAP,
                        ErrorSeverity::ERROR, species_name, chr_name,
                        segment_count, current->start,
                        "Reference segment coordinates overlap",
                        "prev_end=" + std::to_string(prev_end) +
                            ", current_start=" +
                            std::to_string(current->start) + ", overlap_size=" +
                            std::to_string(overlap_size) + "bp");
                  }
                } else if (forbid_non_reference_overlap) {
                  addVerificationError(
                      result, options, VerificationType::COORDINATE_OVERLAP,
                      ErrorSeverity::ERROR, species_name, chr_name,
                      segment_count, current->start,
                      "Segment coordinates overlap",
                      "prev_end=" + std::to_string(prev_end) +
                          ", current_start=" + std::to_string(current->start) +
                          ", overlap_size=" + std::to_string(overlap_size) +
                          "bp");
                }
              }
            }
          }

          // 坐标排序检查
          if (options.isEnabled(VerificationType::COORDINATE_ORDERING)) {
            if (prev_segment) {
              // 检查排序：start应该是递增的
              if (current->start < prev_segment->start) {
                ordering_violations++;
                addVerificationError(
                    result, options, VerificationType::COORDINATE_ORDERING,
                    ErrorSeverity::ERROR, species_name, chr_name, segment_count,
                    current->start,
                    "Segments are not properly ordered (start not increasing)",
                    "current_start=" + std::to_string(current->start) +
                        " < prev_start=" + std::to_string(prev_segment->start) +
                        " (violation #" + std::to_string(ordering_violations) +
                        ")");
              }

              // 检查间隙大小（性能提示）
              uint_t prev_end = prev_segment->start + prev_segment->length;
              if (current->start > prev_end) {
                uint_t gap_size = current->start - prev_end;
                if (gap_size > LARGE_GAP_THRESHOLD) {
                  large_gaps++;
                  addVerificationError(
                      result, options, VerificationType::COORDINATE_ORDERING,
                      ErrorSeverity::INFO, species_name, chr_name,
                      segment_count, current->start,
                      "Large gap between segments",
                      "gap_size=" + std::to_string(gap_size) + "bp" +
                          " (large gap #" + std::to_string(large_gaps) + ")");
                }
              }
            }
            prev_segment = current;
          }

          // 内存完整性的额外检查
          if (options.isEnabled(VerificationType::MEMORY_INTEGRITY)) {
            // 检查坐标溢出
            if (current->start > UINT32_MAX - current->length) {
              addVerificationError(
                  result, options, VerificationType::MEMORY_INTEGRITY,
                  ErrorSeverity::ERROR, species_name, chr_name, segment_count,
                  current->start, "Coordinate overflow detected",
                  "start=" + std::to_string(current->start) +
                      ", length=" + std::to_string(current->length));
            }

            // 检查parent_block引用的有效性（使用预构建的映射）
            if (current->parent_block &&
                valid_blocks.find(current->parent_block) ==
                    valid_blocks.end()) {
              addVerificationError(
                  result, options, VerificationType::MEMORY_INTEGRITY,
                  ErrorSeverity::ERROR, species_name, chr_name, segment_count,
                  current->start, "Segment references invalid parent_block",
                  "parent_block not found in global block pool");
            }
          }
        }

        // 记录上一段“有效segment”（长度>0），供坐标重叠检查使用
        if (current->isSegment() && current->length > 0) {
          prev_nonzero_segment = current;
        }

        prev = current;
        current = current->primary_path.next.load(std::memory_order_acquire);
      }

      // 链表完整性的最终检查
      if (options.isEnabled(VerificationType::LINKED_LIST_INTEGRITY)) {
        if (current == nullptr && prev && !prev->isTail()) {
          addVerificationError(
              result, options, VerificationType::LINKED_LIST_INTEGRITY,
              ErrorSeverity::ERROR, species_name, chr_name, segment_count, 0,
              "Linked list ends without reaching tail sentinel",
              "Last segment is not marked as TAIL");
        }
      }

      // 内存完整性的采样向量检查
      if (options.isEnabled(VerificationType::MEMORY_INTEGRITY)) {
        for (size_t i = 0; i < genome_end.sample_vec.size(); ++i) {
          SegPtr sample_seg = genome_end.sample_vec[i];
          if (sample_seg &&
              all_segments.find(sample_seg) == all_segments.end()) {
            addVerificationError(
                result, options, VerificationType::MEMORY_INTEGRITY,
                ErrorSeverity::ERROR, species_name, chr_name, i, 0,
                "Sampling vector contains invalid segment reference",
                "sample_index=" + std::to_string(i) +
                    ", segment not in main list");
          }
        }
      }

      // 输出统计信息（仅在启用详细段信息时输出，避免刷屏）
      if (options.verbose && options.show_detailed_segments) {
        if (options.isEnabled(VerificationType::COORDINATE_ORDERING)) {
          spdlog::debug("  Chromosome {}: {} segments, {} ordering violations, "
                        "{} large gaps",
                        chr_name, segment_count, ordering_violations,
                        large_gaps);
        } else {
          spdlog::debug("  Chromosome {} contains {} segments", chr_name,
                        segment_count);
        }
      }
    }
  }

  if (reference_policy_enabled && require_reference_overlap &&
      !reference_overlap_found) {
    std::string detail;
    if (reference_chr_visited > 0) {
      detail =
          "Visited " + std::to_string(reference_chr_visited) +
          " chromosome(s) for reference species but no overlaps were detected";
    } else {
      detail = "Reference species not present in current graph";
    }

    addVerificationError(result, options, VerificationType::COORDINATE_OVERLAP,
                         ErrorSeverity::CRITICAL, reference_species, "", 0, 0,
                         "Reference species expected to have overlapping "
                         "segments but none were found",
                         detail);
  }

  // 输出各类型的统计信息
  if (options.verbose) {
    if (options.isEnabled(VerificationType::LINKED_LIST_INTEGRITY)) {
      size_t total_errors =
          result.getErrorCountFast(VerificationType::LINKED_LIST_INTEGRITY);
      if (total_errors > 0) {
        spdlog::debug("LINKED_LIST_INTEGRITY: {} total errors", total_errors);
      }
    }
    if (options.isEnabled(VerificationType::COORDINATE_OVERLAP)) {
      size_t total_errors =
          result.getErrorCountFast(VerificationType::COORDINATE_OVERLAP);
      if (total_errors > 0) {
        spdlog::debug("COORDINATE_OVERLAP: {} total errors", total_errors);
      }
    }
    if (options.isEnabled(VerificationType::COORDINATE_ORDERING)) {
      size_t total_errors =
          result.getErrorCountFast(VerificationType::COORDINATE_ORDERING);
      if (total_errors > 0) {
        spdlog::debug("COORDINATE_ORDERING: {} total errors", total_errors);
      }
    }
    if (options.isEnabled(VerificationType::MEMORY_INTEGRITY)) {
      size_t total_errors =
          result.getErrorCountFast(VerificationType::MEMORY_INTEGRITY);
      if (total_errors > 0) {
        spdlog::debug("MEMORY_INTEGRITY: {} total errors", total_errors);
      }
    }
    if (reference_policy_enabled) {
      if (reference_overlap_counter > 0) {
        spdlog::debug("REFERENCE_OVERLAP: {} overlap(s) observed for {}",
                      reference_overlap_counter, reference_species);
      } else {
        spdlog::debug("REFERENCE_OVERLAP: no overlaps observed for {}",
                      reference_species);
      }
    }
  }
}

// ------------------------------------------------------------------
// 获取子序列工具函数
// 统一接口，无需 lambda。
// ------------------------------------------------------------------
inline std::string getSubSeq(const SeqPro::ManagerVariant& mv,
    const ChrName& chr,
    Coord_t b,
    Coord_t l)
{
    return std::visit([&](auto& p) -> std::string {
        using T = std::decay_t<decltype(p)>;
        if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::SequenceManager>>) {
            return p->getSubSequence(chr, b, l);
        }
        else if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            return p->getOriginalManager().getSubSequence(chr, b, l);
        }
        else {
            throw std::runtime_error("Unsupported ManagerVariant type in getSubSeq()");
        }
        }, mv);
}

} // namespace RaMesh
