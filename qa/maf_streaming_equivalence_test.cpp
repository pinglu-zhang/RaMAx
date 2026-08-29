#ifdef NDEBUG
#undef NDEBUG
#endif

#include "ramesh.h"

#include "SeqPro.h"
#include "hal/export.h"
#include "softmask_index.h"

#include "hal.h"

#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace {

using RaMesh::AlignRole;
using RaMesh::Block;
using RaMesh::BlockPtr;
using RaMesh::GenomeEnd;
using RaMesh::Segment;
using RaMesh::SegmentRole;
using RaMesh::SegPtr;
using RaMesh::SpeciesChrPair;

void writeFasta(const std::filesystem::path& path, char first, char second) {
  std::string sequence(256, first);
  for (size_t index = 0; index < sequence.size(); index += 11) {
    sequence[index] = second;
  }
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  output << ">chr1\n" << sequence << '\n';
  assert(output.good());
}

SegPtr addSegment(const BlockPtr& block, const std::string& species,
                  uint_t start, uint_t length, Cigar_t cigar) {
  auto segment = Segment::create(
      start, length, Strand::FORWARD, std::move(cigar),
      AlignRole::PRIMARY, SegmentRole::SEGMENT, block);
  block->anchors.emplace(SpeciesChrPair{species, "chr1"}, segment);
  return segment;
}

void linkPath(GenomeEnd& end, const std::vector<SegPtr>& segments) {
  SegPtr previous = end.head;
  for (const auto& segment : segments) {
    previous->primary_path.next.store(segment, std::memory_order_release);
    segment->primary_path.prev.store(previous, std::memory_order_release);
    previous = segment;
  }
  previous->primary_path.next.store(end.tail, std::memory_order_release);
  end.tail->primary_path.prev.store(previous, std::memory_order_release);
}

void verifyDensePresenceApi() {
  using RaMesh::hal_export::TreeMeta;
  using RaMesh::hal_export::TreeNodeMeta;

  TreeMeta tree;
  tree.nodes = {
      TreeNodeMeta{0, "root", -1, {1, 4}, 0.0, false, -1},
      TreeNodeMeta{1, "ancestor", 0, {2, 3}, 1.0, false, -1},
      TreeNodeMeta{2, "human", 1, {}, 1.0, true, 0},
      TreeNodeMeta{3, "chimp", 1, {}, 1.0, true, 1},
      TreeNodeMeta{4, "gorilla", 0, {}, 1.0, true, 2},
  };
  tree.name_to_id = {{"root", 0}, {"ancestor", 1}, {"human", 2},
                     {"chimp", 3}, {"gorilla", 4}};
  tree.leaf_ids = {2, 3, 4};
  tree.internal_postorder = {1, 0};
  tree.root_id = 0;

  const auto inference = RaMesh::hal_export::inferDescendantUnion(
      tree, {{"human", true}, {"chimp", false}, {"gorilla", true}});
  assert((inference.present_by_node ==
          std::vector<uint8_t>{1, 1, 1, 0, 1}));
  assert((inference.score0 ==
          std::vector<double>{1.0, 1.0, 1.0, 0.0, 1.0}));
  assert((inference.score1 ==
          std::vector<double>{0.0, 0.0, 0.0, 1.0, 0.0}));
  assert(inference.margin ==
         std::vector<double>(tree.nodes.size(), 1.0));
}

void writeHalSnapshot(const std::filesystem::path& hal_path,
                      const std::filesystem::path& snapshot_path) {
  const auto alignment = hal::openHalAlignment(hal_path.string());
  assert(alignment);
  hal::validateAlignment(alignment.get());

  std::set<std::string> pending{alignment->getRootName()};
  std::set<std::string> visited;
  std::ofstream output(snapshot_path, std::ios::binary | std::ios::trunc);
  output << "root\t" << alignment->getRootName() << '\n';
  while (!pending.empty()) {
    const std::string genome_name = *pending.begin();
    pending.erase(pending.begin());
    if (!visited.insert(genome_name).second) {
      continue;
    }
    const auto* genome = alignment->openGenomeCheck(genome_name);
    output << "genome\t" << genome_name << '\t'
           << alignment->getParentName(genome_name) << '\t'
           << genome->getNumSequences() << '\t'
           << genome->getSequenceLength() << '\t'
           << genome->getNumTopSegments() << '\t'
           << genome->getNumBottomSegments() << '\n';
    auto sequence = genome->getSequenceIterator();
    while (!sequence->atEnd()) {
      const auto* current = sequence->getSequence();
      output << "sequence\t" << genome_name << '\t'
             << current->getName() << '\t'
             << current->getSequenceLength() << '\t'
             << current->getNumTopSegments() << '\t'
             << current->getNumBottomSegments() << '\n';
      sequence->toNext();
    }
    for (const auto& child : alignment->getChildNames(genome_name)) {
      pending.insert(child);
    }
    sequence.reset();
    alignment->closeGenome(genome);
  }
  assert(output.good());
}

}  // namespace

int main(int argc, char** argv) {
  verifyDensePresenceApi();
  const auto stamp = std::chrono::steady_clock::now()
      .time_since_epoch().count();
  const bool keep_output = argc >= 2;
  const bool export_hal = argc >= 3;
  const std::filesystem::path output_path = keep_output
      ? std::filesystem::path(argv[1])
      : std::filesystem::current_path() /
            ("ramax-maf-streaming-" + std::to_string(stamp) + ".maf");
  const std::filesystem::path root = output_path.parent_path() /
      (output_path.stem().string() + "-inputs-" + std::to_string(stamp));
  std::filesystem::create_directories(root);

  const auto human_fasta = root / "human.fa";
  const auto chimp_fasta = root / "chimp.fa";
  writeFasta(human_fasta, 'A', 'C');
  writeFasta(chimp_fasta, 'A', 'G');

  const auto human_upper = root / "human.upper.fa";
  const auto chimp_upper = root / "chimp.upper.fa";
  const auto human_index = root / "human.softmask.idx";
  const auto chimp_index = root / "chimp.softmask.idx";
  assert(!SoftMask::ensureUppercaseFastaAndIndex(
      human_fasta, human_upper, human_index,
      root / "human.softmask.complete"));
  assert(!SoftMask::ensureUppercaseFastaAndIndex(
      chimp_fasta, chimp_upper, chimp_index,
      root / "chimp.softmask.complete"));

  std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
  for (const auto& [species, fasta] :
       std::vector<std::pair<std::string, std::filesystem::path>>{
           {"human", human_upper}, {"chimp", chimp_upper}}) {
    SeqPro::ManagerVariant manager =
        std::make_unique<SeqPro::SequenceManager>(fasta);
    managers.emplace(
        species,
        std::make_shared<SeqPro::ManagerVariant>(std::move(manager)));
  }

  RaMesh::RaMeshMultiGenomeGraph graph(managers);
  std::vector<BlockPtr> blocks;
  std::vector<SegPtr> human_segments;
  std::vector<SegPtr> chimp_segments;

  auto first = Block::createEmpty("human", "chr1", 2);
  human_segments.push_back(addSegment(
      first, "human", 0, 8, Cigar_t{}));
  chimp_segments.push_back(addSegment(
      first, "chimp", 100, 8, Cigar_t{cigarToInt('M', 8)}));
  blocks.push_back(first);

  auto second = Block::createEmpty("human", "chr1", 2);
  human_segments.push_back(addSegment(
      second, "human", 8, 7, Cigar_t{}));
  chimp_segments.push_back(addSegment(
      second, "chimp", 108, 8,
      Cigar_t{cigarToInt('M', 4), cigarToInt('I', 1),
              cigarToInt('M', 3)}));
  blocks.push_back(second);

  auto third = Block::createEmpty("human", "chr1", 2);
  human_segments.push_back(addSegment(
      third, "human", 15, 6, Cigar_t{}));
  chimp_segments.push_back(addSegment(
      third, "chimp", 116, 6, Cigar_t{cigarToInt('M', 6)}));
  blocks.push_back(third);

  linkPath(graph.species_graphs.at("human").chr2end.at("chr1"),
           human_segments);
  linkPath(graph.species_graphs.at("chimp").chr2end.at("chr1"),
           chimp_segments);
  // Deliberately scramble the weak Block pool.  Export order must continue
  // to be determined by the legacy reference path/start/block-id key.
  for (size_t index : std::vector<size_t>{2, 0, 1}) {
    graph.blocks.emplace_back(blocks[index]);
  }

  assert(graph.verifyGraphCorrectness(false));
  graph.exportToMaf(output_path, managers, false);
  std::ifstream input(output_path, std::ios::binary);
  const std::string maf((std::istreambuf_iterator<char>(input)),
                        std::istreambuf_iterator<char>());
  const std::string expected =
      "##maf version=1 scoring=none\n"
      "a score=0\n"
      "s chimp.chr1                   100           8 +         256 AAAAAAAA\n"
      "s human.chr1                     0           8 +         256 CAAAAAAA\n"
      "\n"
      "a score=0\n"
      "s chimp.chr1                   108           8 +         256 AAGAAAAA\n"
      "s human.chr1                     8           7 +         256 AAAC-AAA\n"
      "\n"
      "a score=0\n"
      "s chimp.chr1                   116           6 +         256 AAAAAG\n"
      "s human.chr1                    15           6 +         256 AAAAAA\n"
      "\n";
  assert(maf == expected);

  if (export_hal) {
    const std::filesystem::path hal_path(argv[2]);
    const SoftMask::IndexMap softmask_indexes = SoftMask::loadIndexes(
        {{"human", human_index}, {"chimp", chimp_index}});
    RaMesh::hal_export::ExportStats stats;
    RaMesh::hal_export::exportToHal(
        graph.blocks, hal_path, managers,
        NewickParser("(human:1,chimp:1)ancestor;"), "ancestor",
        softmask_indexes,
        RaMesh::hal_export::ExportConfig{.parallel_threads = 4},
        &stats);
    assert(stats.block_count == blocks.size());
    writeHalSnapshot(hal_path, hal_path.string() + ".snapshot");
  }

  std::filesystem::remove_all(root);
  if (!keep_output) {
    std::filesystem::remove(output_path);
  }
  return 0;
}
