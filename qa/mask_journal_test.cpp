#include "rare_aligner.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using Intervals = std::vector<SeqPro::MaskInterval>;
using ManagerMap =
    std::map<SpeciesName, SeqPro::SharedManagerVariant>;

void require(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

Intervals legacyNormalize(Intervals intervals) {
  std::sort(intervals.begin(), intervals.end(),
            [](const SeqPro::MaskInterval& left,
               const SeqPro::MaskInterval& right) {
              return left.start < right.start;
            });
  if (intervals.empty()) {
    return intervals;
  }

  Intervals merged;
  merged.reserve(intervals.size());
  merged.push_back(intervals.front());
  for (size_t index = 1; index < intervals.size(); ++index) {
    if (merged.back().end >= intervals[index].start) {
      merged.back().end =
          std::max(merged.back().end, intervals[index].end);
    } else {
      merged.push_back(intervals[index]);
    }
  }
  return merged;
}

void requireSameIntervals(const Intervals& actual,
                          const Intervals& expected,
                          const std::string& context) {
  require(actual.size() == expected.size(),
          context + ": interval count differs");
  for (size_t index = 0; index < actual.size(); ++index) {
    require(actual[index].start == expected[index].start &&
                actual[index].end == expected[index].end,
            context + ": interval differs at index " +
                std::to_string(index));
  }
}

void verifyIncrementalMerge(const Intervals& history,
                            const Intervals& delta,
                            const std::string& context) {
  constexpr SeqPro::SequenceId sequence_id = 7;
  SeqPro::MaskManager manager;
  manager.addMaskIntervals(sequence_id, history);
  manager.finalizeMaskIntervals(sequence_id);

  const Intervals normalized_history = legacyNormalize(history);
  const Intervals normalized_delta = legacyNormalize(delta);
  Intervals expected = history;
  expected.insert(expected.end(), delta.begin(), delta.end());
  expected = legacyNormalize(std::move(expected));

  const SeqPro::MaskBatchMergeStats stats =
      manager.mergeFinalizedMaskIntervals(sequence_id, delta);
  requireSameIntervals(manager.getMaskIntervals(sequence_id), expected,
                       context);
  require(stats.incoming_intervals == delta.size(),
          context + ": incoming interval count differs");
  require(stats.normalized_delta_intervals == normalized_delta.size(),
          context + ": normalized delta count differs");
  require(stats.previous_intervals == normalized_history.size(),
          context + ": previous interval count differs");
  require(stats.final_intervals == expected.size(),
          context + ": final interval count differs");
  require(stats.touched_sequences == (delta.empty() ? 0U : 1U),
          context + ": touched sequence count differs");

  constexpr SeqPro::Length original_length = 5000;
  const SeqPro::FinalizedMaskSummary summary =
      manager.summarizeFinalizedIntervals(sequence_id, original_length);
  require(summary.masked_sequence_length ==
              manager.getMaskedSequenceLength(sequence_id, original_length),
          context + ": one-pass masked length differs");
  require(summary.separator_count ==
              manager.getSeparatorCount(sequence_id, original_length),
          context + ": one-pass separator count differs");
  require(summary.masked_bases ==
              original_length - summary.masked_sequence_length,
          context + ": one-pass masked base count differs");
}

void testIncrementalMergeBoundaries() {
  struct Case {
    Intervals history;
    Intervals delta;
  };
  const std::vector<Case> cases{
      {{}, {}},
      {{}, {{10, 20}}},
      {{{10, 20}}, {}},
      {{{10, 20}}, {{20, 30}}},
      {{{10, 20}}, {{19, 30}}},
      {{{10, 30}}, {{15, 20}}},
      {{{20, 30}}, {{10, 20}}},
      {{{10, 20}, {30, 40}}, {{18, 32}}},
      {{{10, 20}, {40, 50}}, {{25, 30}}},
      {{{100, 120}}, {{1, 2}, {2, 3}, {120, 121}}},
      {{{0, 1}, {std::numeric_limits<uint64_t>::max() - 20,
                  std::numeric_limits<uint64_t>::max() - 10}},
       {{1, 2}, {std::numeric_limits<uint64_t>::max() - 12,
                  std::numeric_limits<uint64_t>::max()}}},
  };

  for (size_t index = 0; index < cases.size(); ++index) {
    verifyIncrementalMerge(cases[index].history, cases[index].delta,
                           "boundary case " + std::to_string(index));
  }
}

Intervals randomIntervals(std::mt19937_64& random, size_t count) {
  Intervals intervals;
  intervals.reserve(count);
  for (size_t index = 0; index < count; ++index) {
    SeqPro::Position start = random() % 4096;
    if ((random() & 31U) == 0U) {
      start += uint64_t{1} << 32;
    }
    const SeqPro::Length length = 1 + random() % 256;
    intervals.emplace_back(start, start + length);
  }
  return intervals;
}

void testRandomizedIncrementalMerge() {
  std::mt19937_64 random(0x4d41534b4a4f5552ULL);
  constexpr size_t trials = 100000;
  for (size_t trial = 0; trial < trials; ++trial) {
    const Intervals history =
        randomIntervals(random, static_cast<size_t>(random() % 16));
    const Intervals delta =
        randomIntervals(random, static_cast<size_t>(random() % 16));
    verifyIncrementalMerge(history, delta,
                           "random trial " + std::to_string(trial));
  }
}

struct TemporaryDirectory {
  std::filesystem::path path;

  TemporaryDirectory() {
    const auto stamp = std::chrono::steady_clock::now()
                           .time_since_epoch()
                           .count();
    path = std::filesystem::temp_directory_path() /
           ("ramax-mask-journal-test-" + std::to_string(stamp));
    std::filesystem::create_directories(path);
  }

  ~TemporaryDirectory() {
    std::error_code error;
    std::filesystem::remove_all(path, error);
  }
};

void writeFasta(const std::filesystem::path& path) {
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  require(static_cast<bool>(output),
          "cannot create test FASTA: " + path.string());
  output << ">chr1\n" << std::string(1024, 'A') << '\n';
  output << ">chr2\n" << std::string(257, 'C') << '\n';
  output.flush();
  require(static_cast<bool>(output),
          "cannot finish test FASTA: " + path.string());
}

std::unique_ptr<SeqPro::MaskedSequenceManager> makeMaskedManager(
    const std::filesystem::path& fasta) {
  auto original = std::make_unique<SeqPro::SequenceManager>(fasta);
  return std::make_unique<SeqPro::MaskedSequenceManager>(
      std::move(original));
}

void requireSameManagerState(
    SeqPro::MaskedSequenceManager& actual,
    SeqPro::MaskedSequenceManager& expected,
    const std::string& context) {
  require(actual.getSequenceNames() == expected.getSequenceNames(),
          context + ": sequence names differ");
  for (const auto& sequence : actual.getSequenceNames()) {
    const SeqPro::SequenceId actual_id = actual.getSequenceId(sequence);
    const SeqPro::SequenceId expected_id = expected.getSequenceId(sequence);
    requireSameIntervals(actual.getMaskIntervals(actual_id),
                         expected.getMaskIntervals(expected_id),
                         context + ":" + sequence);
    require(actual.getSequenceLength(actual_id) ==
                expected.getSequenceLength(expected_id),
            context + ": masked sequence length differs");
    require(actual.getSeparatorCount(actual_id) ==
                expected.getSeparatorCount(expected_id),
            context + ": separator count differs");
    require(actual.getSequenceLengthWithSeparators(actual_id) ==
                expected.getSequenceLengthWithSeparators(expected_id),
            context + ": length with separators differs");

    const SeqPro::SequenceInfo* actual_info =
        actual.getOriginalManager().getIndex().getSequenceInfo(actual_id);
    const SeqPro::SequenceInfo* expected_info =
        expected.getOriginalManager().getIndex().getSequenceInfo(expected_id);
    require(actual_info && expected_info,
            context + ": missing sequence metadata");
    require(actual_info->masked_global_start_pos ==
                expected_info->masked_global_start_pos &&
                actual_info->masked_length == expected_info->masked_length,
            context + ": finalized sequence metadata differs");

    if (actual.getSequenceLength(actual_id) > 0) {
      require(actual.localToGlobal(actual_id, 0) ==
                  expected.localToGlobal(expected_id, 0),
              context + ": global offset cache differs");
      require(actual.localToGlobalSeparated(actual_id, 0) ==
                  expected.localToGlobalSeparated(expected_id, 0),
              context + ": separated global offset differs");
    }
  }
  require(actual.getTotalLength() == expected.getTotalLength(),
          context + ": total masked length differs");
  require(actual.getTotalLengthWithSeparators() ==
              expected.getTotalLengthWithSeparators(),
          context + ": total separated length differs");
  require(actual.getTotalSeparatorCount() ==
              expected.getTotalSeparatorCount(),
          context + ": total separator count differs");
}

void testMetadataRefresh(const std::filesystem::path& fasta) {
  auto legacy = makeMaskedManager(fasta);
  auto incremental = makeMaskedManager(fasta);

  const Intervals chr1_history{{10, 20}, {30, 40}};
  const Intervals chr2_history{{0, 3}, {100, 110}};
  const Intervals chr1_delta{{18, 35}, {900, 1024}};
  const Intervals chr2_delta{{3, 9}, {200, 300}};

  Intervals chr1_all = chr1_history;
  chr1_all.insert(chr1_all.end(), chr1_delta.begin(), chr1_delta.end());
  Intervals chr2_all = chr2_history;
  chr2_all.insert(chr2_all.end(), chr2_delta.begin(), chr2_delta.end());
  legacy->addMaskIntervals("chr1", chr1_all);
  legacy->addMaskIntervals("chr2", chr2_all);
  legacy->finalizeMaskIntervals();

  incremental->addMaskIntervals("chr1", chr1_history);
  incremental->addMaskIntervals("chr2", chr2_history);
  incremental->finalizeMaskIntervals();
  std::vector<SeqPro::MaskIntervalDelta> deltas;
  deltas.push_back(
      {incremental->getSequenceId("chr1"), chr1_delta});
  deltas.push_back(
      {incremental->getSequenceId("chr2"), chr2_delta});
  const SeqPro::MaskBatchMergeStats stats =
      incremental->applyFinalizedMaskDeltas(std::move(deltas));
  require(stats.touched_sequences == 2,
          "metadata refresh touched sequence count differs");
  requireSameManagerState(*incremental, *legacy,
                          "incremental metadata refresh");

  const auto before_chr1 = incremental->getMaskIntervals("chr1");
  const auto before_chr2 = incremental->getMaskIntervals("chr2");
  const SeqPro::MaskBatchMergeStats empty =
      incremental->applyFinalizedMaskDeltas({});
  require(empty.touched_sequences == 0 &&
              empty.incoming_intervals == 0 &&
              empty.sort_seconds == 0.0 &&
              empty.merge_seconds == 0.0 &&
              empty.metadata_seconds == 0.0,
          "empty delta must not sort, merge, or refresh metadata");
  requireSameIntervals(incremental->getMaskIntervals("chr1"), before_chr1,
                       "empty delta chr1");
  requireSameIntervals(incremental->getMaskIntervals("chr2"), before_chr2,
                       "empty delta chr2");
}

ManagerMap makeManagers(const std::filesystem::path& reference,
                        const std::filesystem::path& query) {
  ManagerMap managers;
  for (const auto& [species, path] :
       std::vector<std::pair<std::string, std::filesystem::path>>{
           {"reference", reference}, {"query", query}}) {
    std::unique_ptr<SeqPro::MaskedSequenceManager> masked =
        makeMaskedManager(path);
    managers.emplace(
        species,
        std::make_shared<SeqPro::ManagerVariant>(std::move(masked)));
  }
  return managers;
}

SeqPro::MaskedSequenceManager& maskedManager(
    ManagerMap& managers, const SpeciesName& species) {
  auto& manager =
      std::get<std::unique_ptr<SeqPro::MaskedSequenceManager>>(
          *managers.at(species));
  require(static_cast<bool>(manager),
          "missing masked manager for " + species);
  return *manager;
}

void legacyFullMaskScan(const RaMesh::RaMeshMultiGenomeGraph& graph,
                        ManagerMap& managers) {
  using Collected = std::unordered_map<
      SpeciesName,
      std::unordered_map<ChrName, Intervals>>;
  Collected collected;
  for (const auto& weak_block : graph.blocks) {
    const auto block = weak_block.lock();
    if (!block) {
      continue;
    }
    for (const auto& [key, segment] : block->anchors) {
      if (!segment || !segment->isSegment() || segment->length == 0 ||
          !managers.contains(key.first)) {
        continue;
      }
      collected[key.first][key.second].emplace_back(
          segment->start, segment->start + segment->length);
    }
  }

  for (auto& [species, chromosomes] : collected) {
    auto& manager = maskedManager(managers, species);
    bool changed = false;
    for (auto& [chromosome, intervals] : chromosomes) {
      if (manager.getSequenceId(chromosome) ==
          SeqPro::SequenceIndex::INVALID_ID) {
        continue;
      }
      manager.addMaskIntervals(chromosome, intervals);
      changed = true;
    }
    if (changed) {
      manager.finalizeMaskIntervals();
    }
  }
}

struct AddedBlock {
  RaMesh::BlockPtr block;
  RaMesh::SegPtr reference_segment;
  RaMesh::SegPtr query_segment;
};

AddedBlock appendPairBlock(RaMesh::RaMeshMultiGenomeGraph& graph,
                           uint_t reference_start,
                           uint_t query_start,
                           uint_t length) {
  AddedBlock added;
  added.block =
      RaMesh::Block::createEmpty("reference", "chr1", 2);
  added.reference_segment = RaMesh::Segment::create(
      reference_start, length, FORWARD, {},
      RaMesh::AlignRole::PRIMARY, RaMesh::SegmentRole::SEGMENT,
      added.block);
  added.query_segment = RaMesh::Segment::create(
      query_start, length, FORWARD, {},
      RaMesh::AlignRole::PRIMARY, RaMesh::SegmentRole::SEGMENT,
      added.block);
  added.block->anchors.emplace(
      RaMesh::SpeciesChrPair{"reference", "chr1"},
      added.reference_segment);
  added.block->anchors.emplace(
      RaMesh::SpeciesChrPair{"query", "chr1"},
      added.query_segment);
  graph.blocks.emplace_back(added.block);
  return added;
}

RaMesh::SegPtr appendSingleSegmentBlock(
    RaMesh::RaMeshMultiGenomeGraph& graph,
    std::vector<RaMesh::BlockPtr>& owners,
    const SpeciesName& species,
    const ChrName& chromosome,
    uint_t start,
    uint_t length) {
  auto block =
      RaMesh::Block::createEmpty("reference", "chr1", 1);
  auto segment = RaMesh::Segment::create(
      start, length, FORWARD, {}, RaMesh::AlignRole::PRIMARY,
      RaMesh::SegmentRole::SEGMENT, block);
  block->anchors.emplace(
      RaMesh::SpeciesChrPair{species, chromosome}, segment);
  graph.blocks.emplace_back(block);
  owners.push_back(std::move(block));
  return segment;
}

void requireSameGraphMasks(ManagerMap& actual,
                           ManagerMap& expected,
                           const std::string& context) {
  for (const SpeciesName species : {"reference", "query"}) {
    requireSameManagerState(maskedManager(actual, species),
                            maskedManager(expected, species),
                            context + ":" + species);
  }
}

void testGraphMaskJournal(const std::filesystem::path& reference,
                          const std::filesystem::path& query) {
  auto graph_managers = makeManagers(reference, query);
  auto legacy_managers = makeManagers(reference, query);
  auto journal_managers = makeManagers(reference, query);
  RaMesh::RaMeshMultiGenomeGraph graph(graph_managers);

  AddedBlock first = appendPairBlock(graph, 10, 20, 30);
  AddedBlock second = appendPairBlock(graph, 100, 120, 25);

  legacyFullMaskScan(graph, legacy_managers);
  addAlignedRegionsAsMask(graph, journal_managers, "reference");
  requireSameGraphMasks(journal_managers, legacy_managers,
                        "first full mask round");
  require(first.reference_segment->mask_journal_start == 10 &&
              first.reference_segment->mask_journal_length == 30 &&
              first.query_segment->mask_journal_start == 20 &&
              first.query_segment->mask_journal_length == 30,
          "first round did not publish journal coordinates");

  const Intervals unchanged_reference =
      maskedManager(journal_managers, "reference")
          .getMaskIntervals("chr1");
  const Intervals unchanged_query =
      maskedManager(journal_managers, "query")
          .getMaskIntervals("chr1");
  legacyFullMaskScan(graph, legacy_managers);
  addAlignedRegionsAsMask(graph, journal_managers, "reference");
  requireSameGraphMasks(journal_managers, legacy_managers,
                        "no-change mask round");
  requireSameIntervals(
      maskedManager(journal_managers, "reference")
          .getMaskIntervals("chr1"),
      unchanged_reference, "no-change reference intervals");
  requireSameIntervals(
      maskedManager(journal_managers, "query")
          .getMaskIntervals("chr1"),
      unchanged_query, "no-change query intervals");

  first.query_segment->start = 60;
  first.query_segment->length = 35;
  legacyFullMaskScan(graph, legacy_managers);
  addAlignedRegionsAsMask(graph, journal_managers, "reference");
  requireSameGraphMasks(journal_managers, legacy_managers,
                        "moved segment round");
  require(first.query_segment->mask_journal_start == 60 &&
              first.query_segment->mask_journal_length == 35,
          "moved segment journal was not refreshed");

  AddedBlock replacement = appendPairBlock(graph, 300, 340, 40);
  legacyFullMaskScan(graph, legacy_managers);
  addAlignedRegionsAsMask(graph, journal_managers, "reference");
  requireSameGraphMasks(journal_managers, legacy_managers,
                        "new block round");

  replacement.block->anchors.clear();
  replacement.reference_segment->parent_block.reset();
  replacement.query_segment->parent_block.reset();
  replacement.block.reset();
  replacement.reference_segment.reset();
  replacement.query_segment.reset();
  legacyFullMaskScan(graph, legacy_managers);
  addAlignedRegionsAsMask(graph, journal_managers, "reference");
  requireSameGraphMasks(journal_managers, legacy_managers,
                        "removed block cumulative mask round");

  std::vector<RaMesh::BlockPtr> special_owners;
  const auto missing = appendSingleSegmentBlock(
      graph, special_owners, "query", "missing_chr", 11, 12);
  const auto zero = appendSingleSegmentBlock(
      graph, special_owners, "query", "chr1", 77, 0);
  const auto unknown = appendSingleSegmentBlock(
      graph, special_owners, "unknown", "chr1", 5, 6);
  (void)zero;
  (void)unknown;

  auto sentinel_block =
      RaMesh::Block::createEmpty("reference", "chr1", 3);
  sentinel_block->anchors.emplace(
      RaMesh::SpeciesChrPair{"query", "chr1"},
      RaMesh::Segment::createHead());
  sentinel_block->anchors.emplace(
      RaMesh::SpeciesChrPair{"reference", "chr1"},
      RaMesh::Segment::createTail());
  sentinel_block->anchors.emplace(
      RaMesh::SpeciesChrPair{"query", "chr2"}, nullptr);
  graph.blocks.emplace_back(sentinel_block);
  special_owners.push_back(std::move(sentinel_block));
  {
    auto expired =
        RaMesh::Block::createEmpty("reference", "chr1", 0);
    graph.blocks.emplace_back(expired);
  }

  legacyFullMaskScan(graph, legacy_managers);
  addAlignedRegionsAsMask(graph, journal_managers, "reference");
  requireSameGraphMasks(journal_managers, legacy_managers,
                        "ignored graph entries round");
  require(missing->mask_journal_length == 0,
          "missing chromosome must not publish a journal snapshot");

  std::unique_ptr<SeqPro::SequenceManager> null_original;
  journal_managers.emplace(
      "broken",
      std::make_shared<SeqPro::ManagerVariant>(
          std::move(null_original)));
  const auto broken = appendSingleSegmentBlock(
      graph, special_owners, "broken", "chr1", 50, 10);
  addAlignedRegionsAsMask(graph, journal_managers, "reference");
  require(broken->mask_journal_length == 0,
          "failed species must not publish a journal snapshot");

  if constexpr (sizeof(uint_t) >= sizeof(uint64_t)) {
    const uint_t large_start =
        static_cast<uint_t>((uint64_t{1} << 32) + 123);
    const uint_t large_length =
        static_cast<uint_t>((uint64_t{1} << 32) + 321);
    const auto large = appendSingleSegmentBlock(
        graph, special_owners, "query", "chr2",
        large_start, large_length);
    legacyFullMaskScan(graph, legacy_managers);
    addAlignedRegionsAsMask(graph, journal_managers, "reference");
    requireSameGraphMasks(journal_managers, legacy_managers,
                          "64-bit journal round");
    require(large->mask_journal_start == large_start &&
                large->mask_journal_length == large_length,
            "64-bit journal coordinates were truncated");
  }

  (void)second;
}

}  // namespace

int main() {
  try {
    testIncrementalMergeBoundaries();
    testRandomizedIncrementalMerge();

    TemporaryDirectory root;
    const auto reference = root.path / "reference.fa";
    const auto query = root.path / "query.fa";
    writeFasta(reference);
    writeFasta(query);
    testMetadataRefresh(reference);
    testGraphMaskJournal(reference, query);
    return 0;
  } catch (const std::exception& error) {
    std::fprintf(stderr, "mask journal QA failed: %s\n", error.what());
    return 1;
  }
}
