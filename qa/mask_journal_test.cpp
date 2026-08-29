#ifdef NDEBUG
#undef NDEBUG
#endif

#include "rare_aligner.h"

#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using ManagerMap =
    std::map<SpeciesName, SeqPro::SharedManagerVariant>;

ManagerMap makeManagers(const std::filesystem::path& reference,
                        const std::filesystem::path& query) {
  ManagerMap managers;
  for (const auto& [species, path] :
       std::vector<std::pair<std::string, std::filesystem::path>>{
           {"reference", reference}, {"query", query}}) {
    auto original = std::make_unique<SeqPro::SequenceManager>(path);
    auto masked =
        std::make_unique<SeqPro::MaskedSequenceManager>(std::move(original));
    managers.emplace(
        species,
        std::make_shared<SeqPro::ManagerVariant>(std::move(masked)));
  }
  return managers;
}

SeqPro::MaskedSequenceManager& maskedManager(
    ManagerMap& managers, const SpeciesName& species) {
  return *std::get<std::unique_ptr<SeqPro::MaskedSequenceManager>>(
      *managers.at(species));
}

void legacyFullMaskScan(const RaMesh::RaMeshMultiGenomeGraph& graph,
                        ManagerMap& managers) {
  using Intervals = std::unordered_map<
      SpeciesName,
      std::unordered_map<ChrName, std::vector<SeqPro::MaskInterval>>>;
  Intervals collected;
  for (const auto& weak_block : graph.blocks) {
    const auto block = weak_block.lock();
    if (!block) continue;
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
    for (auto& [chromosome, intervals] : chromosomes) {
      manager.addMaskIntervals(chromosome, intervals);
    }
    manager.finalizeMaskIntervals();
  }
}

void assertSameMasks(ManagerMap& left, ManagerMap& right) {
  for (const std::string species : {"reference", "query"}) {
    const auto& left_intervals =
        maskedManager(left, species).getMaskIntervals("chr1");
    const auto& right_intervals =
        maskedManager(right, species).getMaskIntervals("chr1");
    assert(left_intervals.size() == right_intervals.size());
    for (size_t index = 0; index < left_intervals.size(); ++index) {
      assert(left_intervals[index].start == right_intervals[index].start);
      assert(left_intervals[index].end == right_intervals[index].end);
    }
  }
}

std::pair<RaMesh::BlockPtr, RaMesh::SegPtr> appendBlock(
    RaMesh::RaMeshMultiGenomeGraph& graph,
    uint_t reference_start,
    uint_t query_start,
    uint_t length) {
  auto block = RaMesh::Block::createEmpty("reference", "chr1", 2);
  auto reference_segment = RaMesh::Segment::create(
      reference_start, length, FORWARD, {}, RaMesh::AlignRole::PRIMARY,
      RaMesh::SegmentRole::SEGMENT, block);
  auto query_segment = RaMesh::Segment::create(
      query_start, length, FORWARD, {}, RaMesh::AlignRole::PRIMARY,
      RaMesh::SegmentRole::SEGMENT, block);
  block->anchors.emplace(
      RaMesh::SpeciesChrPair{"reference", "chr1"}, reference_segment);
  block->anchors.emplace(
      RaMesh::SpeciesChrPair{"query", "chr1"}, query_segment);
  graph.blocks.emplace_back(block);
  return {std::move(block), std::move(query_segment)};
}

}  // namespace

int main() {
  const auto stamp = std::chrono::steady_clock::now()
                         .time_since_epoch().count();
  const auto root = std::filesystem::current_path() /
      ("ramax-mask-journal-test-" + std::to_string(stamp));
  std::filesystem::create_directories(root);
  const auto reference_path = root / "reference.fa";
  const auto query_path = root / "query.fa";
  for (const auto& path : {reference_path, query_path}) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    output << ">chr1\n" << std::string(1000, 'A') << '\n';
  }

  auto graph_managers = makeManagers(reference_path, query_path);
  auto legacy_managers = makeManagers(reference_path, query_path);
  auto journal_managers = makeManagers(reference_path, query_path);
  {
    RaMesh::RaMeshMultiGenomeGraph graph(graph_managers);
    std::vector<RaMesh::BlockPtr> owners;
    auto [first_block, moving_query] = appendBlock(graph, 10, 20, 30);
    owners.push_back(first_block);
    auto [second_block, unused_query] = appendBlock(graph, 100, 120, 25);
    owners.push_back(second_block);

    legacyFullMaskScan(graph, legacy_managers);
    addAlignedRegionsAsMask(graph, journal_managers, "reference");
    assertSameMasks(legacy_managers, journal_managers);

    // A complete no-change round must append nothing while producing the same
    // finalized interval union as the legacy full-graph scan.
    legacyFullMaskScan(graph, legacy_managers);
    addAlignedRegionsAsMask(graph, journal_managers, "reference");
    assertSameMasks(legacy_managers, journal_managers);

    // Moving an existing Segment preserves the legacy cumulative-union
    // behavior: the old interval remains and the new interval is appended.
    moving_query->start = 60;
    moving_query->length = 35;
    legacyFullMaskScan(graph, legacy_managers);
    addAlignedRegionsAsMask(graph, journal_managers, "reference");
    assertSameMasks(legacy_managers, journal_managers);

    auto [third_block, third_query] = appendBlock(graph, 300, 340, 40);
    (void)third_query;
    owners.push_back(third_block);
    legacyFullMaskScan(graph, legacy_managers);
    addAlignedRegionsAsMask(graph, journal_managers, "reference");
    assertSameMasks(legacy_managers, journal_managers);
  }
  std::filesystem::remove_all(root);
  return 0;
}
