#include "ramesh.h"
#include "data_process.h"
#include "hal/export.h"

#include "SeqPro.h"
#include "align.h"
#include "anchor.h"
#include "halAlignmentInstance.h"
#include "halColumnIterator.h"
#include "halGenome.h"
#include <algorithm>

#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <unistd.h>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

bool hasUpperAndLower(const std::string& dna) {
    bool uppercase = false;
    bool lowercase = false;
    for (char base : dna) {
        const auto byte = static_cast<unsigned char>(base);
        uppercase = uppercase || std::isupper(byte) != 0;
        lowercase = lowercase || std::islower(byte) != 0;
    }
    return uppercase && lowercase;
}

void writeInput(const std::filesystem::path& path, const std::string& dna) {
    std::ofstream output(path);
    output << ">chr1\n" << dna << '\n';
    if (!output) throw std::runtime_error("cannot write test FASTA");
}

void writeMultiInput(
    const std::filesystem::path& path,
    const std::vector<std::pair<std::string, std::string>>& sequences) {
    std::ofstream output(path);
    for (const auto& [name, dna] : sequences) {
        output << '>' << name << '\n' << dna << '\n';
    }
    if (!output) {
        throw std::runtime_error(
            "cannot write multi-sequence test FASTA");
    }
}

std::string readGenome(hal::AlignmentPtr alignment, const std::string& name) {
    hal::Genome* genome = alignment->openGenome(name);
    if (!genome) throw std::runtime_error("missing genome in HAL: " + name);
    std::string dna;
    genome->getString(dna);
    alignment->closeGenome(genome);
    return dna;
}

uint64_t halPairCoverage(
    const std::filesystem::path& hal_path,
    const std::string& reference_name,
    const std::string& target_name) {
    hal::AlignmentPtr alignment = hal::openHalAlignment(
        hal_path.string(), nullptr, hal::READ_ACCESS);
    const hal::Genome* reference =
        alignment->openGenome(reference_name);
    const hal::Genome* target =
        alignment->openGenome(target_name);
    require(
        reference != nullptr && target != nullptr,
        "HAL coverage contract is missing a leaf genome");
    std::set<const hal::Genome*> targets{target};
    auto columns = reference->getColumnIterator(
        &targets,
        0,
        0,
        hal::NULL_INDEX,
        false,
        true,
        false,
        true,
        false);
    uint64_t covered = 0;
    while (true) {
        const auto* column = columns->getColumnMap();
        const bool target_present = std::any_of(
            column->begin(),
            column->end(),
            [&](const auto& entry) {
                return entry.first->getGenome() == target &&
                       entry.second != nullptr &&
                       !entry.second->empty();
            });
        covered += static_cast<uint64_t>(target_present);
        if (columns->lastColumn()) {
            break;
        }
        columns->toRight();
    }
    alignment->closeGenome(target);
    alignment->closeGenome(reference);
    return covered;
}

uint64_t mafPairCoverage(
    const std::filesystem::path& maf_path,
    const std::string& reference_name,
    const std::string& target_name) {
    struct MafRow {
        std::string source;
        uint64_t start = 0;
        char strand = '+';
        uint64_t source_size = 0;
        std::string dna;
    };
    std::ifstream input(maf_path);
    if (!input) {
        throw std::runtime_error(
            "cannot open coverage-contract MAF");
    }
    std::set<
        std::pair<std::string, uint64_t>>
        covered_positions;
    std::map<std::string, std::vector<MafRow>>
        rows;
    auto flush = [&] {
        const auto reference =
            rows.find(reference_name);
        const auto target =
            rows.find(target_name);
        if (reference == rows.end() ||
            target == rows.end()) {
            rows.clear();
            return;
        }
        for (const auto& reference_row :
             reference->second) {
            uint64_t non_gap_before = 0;
            for (size_t column = 0;
                 column <
                     reference_row.dna.size();
                 ++column) {
                if (reference_row.dna[column] ==
                    '-') {
                    continue;
                }
                const bool target_present =
                    std::any_of(
                        target->second.begin(),
                        target->second.end(),
                        [&](const MafRow& row) {
                            require(
                                row.dna.size() ==
                                    reference_row
                                        .dna.size(),
                                "coverage-contract MAF "
                                "row lengths differ");
                            return row.dna[column] !=
                                   '-';
                        });
                if (target_present) {
                    const uint64_t position =
                        reference_row.strand == '+'
                        ? reference_row.start +
                              non_gap_before
                        : reference_row.source_size -
                              reference_row.start -
                              non_gap_before -
                              1;
                    covered_positions.emplace(
                        reference_row.source,
                        position);
                }
                ++non_gap_before;
            }
        }
        rows.clear();
    };
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) {
            flush();
            continue;
        }
        if (!line.starts_with("s ")) {
            continue;
        }
        std::istringstream fields(line);
        char tag = '\0';
        MafRow row;
        uint64_t size = 0;
        fields >> tag >> row.source >>
            row.start >> size >>
            row.strand >>
            row.source_size >>
            row.dna;
        require(
            fields && tag == 's',
            "cannot parse coverage-contract MAF row");
        const size_t delimiter =
            row.source.find('.');
        const std::string species =
            row.source.substr(
                0,
                delimiter);
        rows[species].push_back(
            std::move(row));
    }
    flush();
    return covered_positions.size();
}
void testNamespacedLeafSequenceIsNotDuplicated(
    const std::filesystem::path& temp) {
    std::map<
        SpeciesName,
        SeqPro::SharedManagerVariant>
        managers;
    for (const auto& species :
         std::array<std::string, 2>{
             "leafA", "leafB"}) {
        const auto fasta =
            temp /
            (species + ".namespaced.fa");
        std::ofstream output(fasta);
        output << '>' << species
               << ".chr1\nACGT\n";
        output.close();
        SeqPro::ManagerVariant manager =
            std::make_unique<
                SeqPro::SequenceManager>(
                fasta);
        managers[species] =
            std::make_shared<
                SeqPro::ManagerVariant>(
                std::move(manager));
    }

    auto block = RaMesh::Block::create(2);
    block->ref_species = "leafA";
    block->ref_chr = "leafA.chr1";
    const Cigar_t cigar{
        cigarToInt('M', 4)};
    for (const auto& species :
         std::array<std::string, 2>{
             "leafA", "leafB"}) {
        auto segment =
            RaMesh::Segment::create(
                0,
                4,
                Strand::FORWARD,
                cigar,
                RaMesh::AlignRole::PRIMARY,
                RaMesh::SegmentRole::SEGMENT,
                block);
        block->anchors.emplace(
            RaMesh::SpeciesChrPair{
                species,
                species + ".chr1"},
            std::move(segment));
    }

    auto graph_managers = managers;
    RaMesh::RaMeshMultiGenomeGraph graph(
        graph_managers);
    graph.blocks.emplace_back(block);
    const auto maf_path =
        temp / "namespaced.maf";
    graph.exportToMaf(maf_path, managers, false);
    std::ifstream maf_input(maf_path);
    const std::string maf{
        std::istreambuf_iterator<char>(
            maf_input),
        std::istreambuf_iterator<char>()};
    require(
        maf.find("leafA.chr1") !=
                std::string::npos &&
            maf.find(
                "leafA.leafA.chr1") ==
                std::string::npos &&
            maf.find(
                "leafB.leafB.chr1") ==
                std::string::npos,
        "MAF duplicated an existing species "
        "namespace in a leaf sequence name");
}

void testMafReassemblesGappedBlock(
    const std::filesystem::path& temp) {
    std::map<
        SpeciesName,
        SeqPro::SharedManagerVariant>
        managers;
    for (const auto& [species, dna] :
         std::map<std::string, std::string>{
             {"leafA", "AACCGGTT"},
             {"leafB", "AATTGGTT"}}) {
        const auto fasta =
            temp / (species + ".gapped.fa");
        writeInput(fasta, dna);
        SeqPro::ManagerVariant manager =
            std::make_unique<
                SeqPro::SequenceManager>(
                fasta);
        managers[species] =
            std::make_shared<
                SeqPro::ManagerVariant>(
                std::move(manager));
    }

    auto block = RaMesh::Block::create(2);
    block->ref_species = "leafA";
    block->ref_chr = "chr1";
    block->anchors.emplace(
        RaMesh::SpeciesChrPair{
            "leafA",
            "chr1"},
        RaMesh::Segment::create(
            0,
            8,
            Strand::FORWARD,
            Cigar_t{cigarToInt('M', 8)},
            RaMesh::AlignRole::PRIMARY,
            RaMesh::SegmentRole::SEGMENT,
            block));
    block->anchors.emplace(
        RaMesh::SpeciesChrPair{
            "leafB",
            "chr1"},
        RaMesh::Segment::create(
            0,
            8,
            Strand::FORWARD,
            Cigar_t{
                cigarToInt('M', 2),
                cigarToInt('D', 2),
                cigarToInt('I', 2),
                cigarToInt('M', 4)},
            RaMesh::AlignRole::PRIMARY,
            RaMesh::SegmentRole::SEGMENT,
            block));

    auto graph_managers = managers;
    RaMesh::RaMeshMultiGenomeGraph graph(
        graph_managers);
    graph.blocks.emplace_back(block);
    const auto maf_path =
        temp / "gapped-reassembled.maf";
    graph.exportToMaf(
        maf_path,
        managers,
        false);

    std::ifstream input(maf_path);
    size_t block_count = 0;
    std::vector<std::string> aligned_rows;
    std::string line;
    while (std::getline(input, line)) {
        if (line.starts_with("a ")) {
            ++block_count;
            continue;
        }
        if (!line.starts_with("s ")) {
            continue;
        }
        std::istringstream fields(line);
        std::string tag;
        std::string source;
        uint64_t start = 0;
        uint64_t size = 0;
        char strand = '+';
        uint64_t source_size = 0;
        std::string dna;
        fields >> tag >> source >> start >>
            size >> strand >> source_size >>
            dna;
        require(
            fields &&
                size ==
                    static_cast<uint64_t>(
                        std::count_if(
                            dna.begin(),
                            dna.end(),
                            [](char base) {
                                return base != '-';
                            })),
            "reassembled MAF row has an invalid size");
        aligned_rows.push_back(
            std::move(dna));
    }
    require(
        block_count == 1,
        "one gapped graph Block must remain one MAF block");
    require(
        aligned_rows.size() == 2 &&
            aligned_rows[0].size() == 10 &&
            aligned_rows[1].size() == 10 &&
            aligned_rows[0].find('-') !=
                std::string::npos &&
            aligned_rows[1].find('-') !=
                std::string::npos,
        "reassembled MAF did not preserve the gapped alignment");
    require(
        mafPairCoverage(
            maf_path,
            "leafA",
            "leafB") == 6 &&
            mafPairCoverage(
                maf_path,
                "leafB",
                "leafA") == 6,
        "reassembled MAF changed pairwise homology");
}

void testParalogousOccurrencesRemainThreaded(
    const std::filesystem::path &temp,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant> &managers,
    const SoftMask::PathMap &softmask_paths) {
  auto block = RaMesh::Block::create(4);
  block->ref_species = "leafA";
  block->ref_chr = "chr1";
  const Cigar_t cigar{cigarToInt('M', 4)};
  auto reference = RaMesh::Segment::create(0, 4, Strand::FORWARD, cigar,
                                           RaMesh::AlignRole::PRIMARY,
                                           RaMesh::SegmentRole::SEGMENT, block);
  auto secondary_reference = RaMesh::Segment::create(
      4, 4, Strand::FORWARD, cigar, RaMesh::AlignRole::SECONDARY,
      RaMesh::SegmentRole::SEGMENT, block);
  auto first_copy = RaMesh::Segment::create(
      0, 4, Strand::FORWARD, cigar, RaMesh::AlignRole::PRIMARY,
      RaMesh::SegmentRole::SEGMENT, block);
  auto second_copy = RaMesh::Segment::create(
      4, 4, Strand::FORWARD, cigar, RaMesh::AlignRole::PRIMARY,
      RaMesh::SegmentRole::SEGMENT, block);
  block->anchors.emplace(RaMesh::SpeciesChrPair{"leafA", "chr1"}, reference);
  block->anchors.emplace(RaMesh::SpeciesChrPair{"leafA", "chr1"},
                         secondary_reference);
  block->anchors.emplace(RaMesh::SpeciesChrPair{"leafB", "chr1"}, first_copy);
  block->anchors.emplace(RaMesh::SpeciesChrPair{"leafB", "chr1"}, second_copy);
  auto maf_managers = managers;
  RaMesh::RaMeshMultiGenomeGraph maf_graph(maf_managers);
  maf_graph.blocks.emplace_back(block);
  const auto maf_path = temp / "paralogy.maf";
  maf_graph.exportToMaf(maf_path, managers, false);
  std::ifstream maf_input(maf_path);
  size_t maf_rows = 0;
  std::vector<std::string> maf_sources;
  std::string maf_line;
  while (std::getline(maf_input, maf_line)) {
    if (!maf_line.starts_with("s ")) {
      continue;
    }
    ++maf_rows;
    std::istringstream fields(maf_line);
    std::string tag;
    std::string source;
    fields >> tag >> source;
    maf_sources.push_back(std::move(source));
  }
  require(maf_rows == 4, "MAF export dropped a primary or secondary occurrence "
                         "sharing the declared reference key");
  require(
      maf_sources.size() == 4 &&
          maf_sources[0].starts_with("leafA.") &&
          maf_sources[1].starts_with("leafB.") &&
          maf_sources[2].starts_with("leafA.") &&
          maf_sources[3].starts_with("leafB."),
      "MAF paralog rows are not interleaved by copy rank");

  const auto hal_path = temp / "paralogy.hal";
  std::vector<std::weak_ptr<RaMesh::Block>> blocks{block};
  RaMesh::hal_export::ExportStats stats;
  RaMesh::hal_export::exportToHal(
      blocks, hal_path, managers, NewickParser("(leafA:0.1,leafB:0.1)anc0;"),
      "anc0", SoftMask::loadIndexes(softmask_paths), {}, &stats);
  require(stats.observed_occurrence_count == 4 &&
              stats.ancestor_occurrence_count == 2,
          "copy-number reconciliation must keep primary and secondary "
          "occurrences distinct while inferring both ancestral copies");
  require(stats.paralogous_top_count == 0 &&
              stats.paralogy_self_adjacency_count == 0,
          "balanced copies in both children must not be mislabeled as child-only paralogy");

  hal::AlignmentPtr alignment =
      hal::openHalAlignment(hal_path.string(), nullptr, hal::READ_ACCESS);
  require(static_cast<bool>(alignment), "cannot reopen paralogy HAL");
  hal::Genome *ancestor = alignment->openGenome("anc0");
  hal::Genome *leaf_b = alignment->openGenome("leafB");
  require(ancestor != nullptr && leaf_b != nullptr,
          "paralogy HAL is missing expected genomes");
  require(ancestor->getNumSequences() == 1 &&
              ancestor->getNumBottomSegments() == 2 &&
              leaf_b->getNumTopSegments() == 3,
          "paralogy HAL segment cardinalities are inconsistent");
  alignment->closeGenome(leaf_b);
  alignment->closeGenome(ancestor);
  alignment.reset();
}
void testSecondaryLinkRunsUnifyPrimaryOccurrences(
    const std::filesystem::path& temp,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    const SoftMask::PathMap& softmask_paths) {
    const Cigar_t base_cigar{cigarToInt('M', 8)};
    const Cigar_t link_cigar{cigarToInt('M', 4)};
    auto make_leaf_block = [&](const SpeciesName& species) {
        auto block = RaMesh::Block::create(1);
        block->ref_species = species;
        block->ref_chr = "chr1";
        auto segment = RaMesh::Segment::create(
            0,
            8,
            Strand::FORWARD,
            base_cigar,
            RaMesh::AlignRole::PRIMARY,
            RaMesh::SegmentRole::SEGMENT,
            block);
        block->anchors.emplace(
            RaMesh::SpeciesChrPair{species, "chr1"},
            segment);
        return block;
    };

    auto leaf_a_block = make_leaf_block("leafA");
    auto leaf_b_block = make_leaf_block("leafB");
    auto make_link_block = [&](Strand leaf_b_strand) {
        auto block = RaMesh::Block::create(2);
        block->ref_species = "leafA";
        block->ref_chr = "chr1";
        auto leaf_a_link = RaMesh::Segment::create(
            2,
            4,
            Strand::FORWARD,
            link_cigar,
            RaMesh::AlignRole::PRIMARY,
            RaMesh::SegmentRole::SEGMENT,
            block);
        auto leaf_b_link = RaMesh::Segment::create(
            2,
            4,
            leaf_b_strand,
            link_cigar,
            RaMesh::AlignRole::SECONDARY,
            RaMesh::SegmentRole::SEGMENT,
            block);
        block->anchors.emplace(
            RaMesh::SpeciesChrPair{"leafA", "chr1"},
            leaf_a_link);
        block->anchors.emplace(
            RaMesh::SpeciesChrPair{"leafB", "chr1"},
            leaf_b_link);
        return block;
    };
    auto link_block =
        make_link_block(Strand::FORWARD);
    auto conflicting_link_block =
        make_link_block(Strand::REVERSE);

    const auto hal_path = temp / "secondary-link.hal";
    std::vector<std::weak_ptr<RaMesh::Block>> blocks{
        leaf_a_block,
        leaf_b_block,
        link_block,
        conflicting_link_block};
    const auto rejected_blocks =
        RaMesh::hal_export::
            findRejectedSecondaryHomologyBlocks(
                blocks,
                managers);
    require(
        rejected_blocks.size() == 1 &&
            rejected_blocks.contains(
                conflicting_link_block->block_id) &&
            !rejected_blocks.contains(
                link_block->block_id),
        "secondary homology selection did not reject the "
        "conflicting block atomically");

    auto graph_managers = managers;
    RaMesh::RaMeshMultiGenomeGraph graph(
        graph_managers);
    graph.blocks = blocks;
    const auto maf_path =
        temp / "secondary-link.maf";
    graph.exportToMaf(maf_path, managers, false);
    std::ifstream maf_stream(maf_path);
    require(
        static_cast<bool>(maf_stream),
        "failed to open filtered secondary-link MAF");
    size_t maf_block_count = 0;
    std::string maf_line;
    while (std::getline(
        maf_stream,
        maf_line)) {
        if (maf_line.starts_with("a ")) {
            ++maf_block_count;
        }
    }
    require(
        maf_block_count == 1,
        "direct MAF export did not apply the shared "
        "secondary-homology block selection");
    RaMesh::hal_export::ExportStats stats;
    RaMesh::hal_export::exportToHal(
        blocks,
        hal_path,
        managers,
        NewickParser("(leafA:0.1,leafB:0.1)anc0;"),
        "anc0",
        SoftMask::loadIndexes(softmask_paths),
        {},
        &stats);

    require(
        halPairCoverage(hal_path, "leafA", "leafB") == 4 &&
            halPairCoverage(hal_path, "leafB", "leafA") == 4,
        "partial secondary homology link was lost or "
        "incorrectly extrapolated into unsupported flanks");
    require(
        mafPairCoverage(
            maf_path,
            "leafA",
            "leafB") ==
                halPairCoverage(
                    hal_path,
                    "leafA",
                    "leafB") &&
            mafPairCoverage(
                maf_path,
                "leafB",
                "leafA") ==
                halPairCoverage(
                    hal_path,
                    "leafB",
                    "leafA"),
        "secondary homology coverage differs between "
        "direct MAF and HAL");
    require(
        stats.observed_occurrence_count == 6,
        "secondary homology normalization did not split and deduplicate "
        "the two primary leaf intervals at the link boundaries");
}
void testSecondaryMaterializationBridgesReferenceGap(
    const std::filesystem::path& temp,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers) {
    auto graph_managers = managers;
    RaMesh::RaMeshMultiGenomeGraph graph(graph_managers);
    const Cigar_t gapped_cigar{
        cigarToInt('M', 2),
        cigarToInt('D', 4),
        cigarToInt('M', 2)};
    Anchor gapped_primary(
        0,
        0,
        8,
        0,
        0,
        4,
        Strand::FORWARD,
        8,
        4,
        gapped_cigar);
    Anchor second_primary(
        0,
        8,
        4,
        0,
        8,
        4,
        Strand::FORWARD,
        4,
        4,
        Cigar_t{cigarToInt('M', 4)});
    graph.insertAnchorIntoGraph(
        *managers.at("leafA"),
        *managers.at("leafB"),
        "leafA",
        "leafB",
        gapped_primary);
    graph.insertAnchorIntoGraph(
        *managers.at("leafA"),
        *managers.at("leafB"),
        "leafA",
        "leafB",
        second_primary);

    Anchor secondary(
        0,
        2,
        4,
        0,
        8,
        4,
        Strand::FORWARD,
        4,
        4,
        Cigar_t{cigarToInt('M', 4)});
    graph.registerSecondaryAnchorCandidate(
        "leafA",
        "chr1",
        "leafB",
        "chr1",
        secondary,
        true,
        true);
    graph.registerSecondaryAnchorCandidate(
        "leafA",
        "chr1",
        "leafB",
        "chr1",
        secondary,
        true,
        false);
    require(
        graph.materializeSecondaryAlignments() == 1,
        "reciprocal secondary anchors did not bridge a primary reference gap");

    const auto maf_path =
        temp / "secondary-reference-gap.maf";
    graph.exportToMaf(maf_path, managers, false);
    require(
        mafPairCoverage(
            maf_path,
            "leafA",
            "leafB") == 12,
        "secondary reference-gap bridge did not restore all supported "
        "pairwise reference bases");
}
void testScaffoldGapSegmentsTileSequences(
    const std::filesystem::path& temp,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>&
        managers,
    const SoftMask::PathMap& softmask_paths) {
    const Cigar_t cigar{cigarToInt('M', 4)};
    auto make_block = [&](uint64_t start) {
        auto block = RaMesh::Block::create(2);
        block->ref_species = "leafA";
        block->ref_chr = "chr1";
        auto reference = RaMesh::Segment::create(
            start,
            4,
            Strand::FORWARD,
            cigar,
            RaMesh::AlignRole::PRIMARY,
            RaMesh::SegmentRole::SEGMENT,
            block);
        auto query = RaMesh::Segment::create(
            start,
            4,
            Strand::FORWARD,
            cigar,
            RaMesh::AlignRole::PRIMARY,
            RaMesh::SegmentRole::SEGMENT,
            block);
        block->anchors.emplace(
            RaMesh::SpeciesChrPair{"leafA", "chr1"},
            reference);
        block->anchors.emplace(
            RaMesh::SpeciesChrPair{"leafB", "chr1"},
            query);
        return block;
    };
    auto first = make_block(0);
    auto second = make_block(8);
    std::vector<std::weak_ptr<RaMesh::Block>> blocks{
        first,
        second};
    const auto hal_path = temp / "scaffold-gap.hal";
    RaMesh::hal_export::ExportStats stats;
    RaMesh::hal_export::exportToHal(
        blocks,
        hal_path,
        managers,
        NewickParser("(leafA:0.1,leafB:0.1)anc0;"),
        "anc0",
        SoftMask::loadIndexes(softmask_paths),
        {},
        &stats);
    require(
        stats.path_vertex_count == 2 &&
            stats.supported_join_count == 1 &&
            stats.indirect_join_count == 1 &&
            stats.scaffold_gap_bases == 10,
        "positive-gap descendant adjacency must produce one explicit 10N ancestor join");

    hal::AlignmentPtr alignment = hal::openHalAlignment(
        hal_path.string(), nullptr, hal::READ_ACCESS);
    require(
        static_cast<bool>(alignment),
        "cannot reopen scaffold-gap HAL");
    hal::Genome* ancestor = alignment->openGenome("anc0");
    hal::Genome* leaf_a = alignment->openGenome("leafA");
    hal::Genome* leaf_b = alignment->openGenome("leafB");
    require(
        ancestor != nullptr &&
            leaf_a != nullptr &&
            leaf_b != nullptr,
        "scaffold-gap HAL is missing expected genomes");
    std::string ancestor_dna;
    ancestor->getString(ancestor_dna);
    require(
        ancestor->getNumSequences() == 1 &&
            ancestor->getNumBottomSegments() == 3 &&
            ancestor_dna.size() == 18 &&
            ancestor_dna.substr(4, 10) ==
                std::string(10, 'N'),
        "ancestor aligned occurrences and scaffold gap must tile one 18-base sequence");
    require(
        leaf_a->getNumTopSegments() == 3 &&
            leaf_b->getNumTopSegments() == 3,
        "leaf top segments must continuously cover aligned and unaligned intervals");
    alignment->closeGenome(leaf_b);
    alignment->closeGenome(leaf_a);
    alignment->closeGenome(ancestor);
    alignment.reset();
    auto mutable_managers = managers;
    RaMesh::RaMeshMultiGenomeGraph graph(
        mutable_managers);
    graph.blocks = blocks;
    const auto maf_path =
        temp / "scaffold-gap.maf";
    graph.exportToMaf(maf_path, managers, false);
    const uint64_t hal_leaf_a_coverage =
        halPairCoverage(
            hal_path,
            "leafA",
            "leafB");
    const uint64_t hal_leaf_b_coverage =
        halPairCoverage(
            hal_path,
            "leafB",
            "leafA");
    const uint64_t maf_leaf_a_coverage =
        mafPairCoverage(
            maf_path,
            "leafA",
            "leafB");
    const uint64_t maf_leaf_b_coverage =
        mafPairCoverage(
            maf_path,
            "leafB",
            "leafA");
    require(
        hal_leaf_a_coverage == 8 &&
            hal_leaf_b_coverage == 8 &&
            hal_leaf_a_coverage ==
                maf_leaf_a_coverage &&
            hal_leaf_b_coverage ==
                maf_leaf_b_coverage,
        "HAL leaf-pair coverage must exactly equal direct MAF coverage in both reference directions");
}
void testBilateralContextPreservesAncestralCopies(
    const std::filesystem::path& temp,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>&
        base_managers,
    const SoftMask::PathMap& base_softmask_paths) {
    auto managers = base_managers;
    auto softmask_paths = base_softmask_paths;
    const std::string reference_species = "leafR";
    const auto input =
        temp / (reference_species + ".input.fa");
    const auto uppercase =
        temp / (reference_species + ".align-v2.fasta");
    const auto index =
        temp / (reference_species + ".softmask-v1.bin");
    const auto marker =
        temp / (reference_species + ".complete.json");
    writeInput(input, "AAAcccGGttTT");
    SoftMask::ensureUppercaseFastaAndIndex(
        input,
        uppercase,
        index,
        marker);
    softmask_paths[reference_species] = index;
    SeqPro::ManagerVariant manager =
        std::make_unique<SeqPro::SequenceManager>(
            uppercase);
    managers[reference_species] =
        std::make_shared<SeqPro::ManagerVariant>(
            std::move(manager));

    auto block = RaMesh::Block::create(5);
    block->ref_species = reference_species;
    block->ref_chr = "chr1";
    const Cigar_t cigar{cigarToInt('M', 4)};
    auto reference = RaMesh::Segment::create(
        0,
        4,
        Strand::FORWARD,
        cigar,
        RaMesh::AlignRole::PRIMARY,
        RaMesh::SegmentRole::SEGMENT,
        block);
    block->anchors.emplace(
        RaMesh::SpeciesChrPair{
            reference_species,
            "chr1"},
        reference);
    for (const std::string species :
         {"leafA", "leafB"}) {
        for (uint64_t start :
             {uint64_t{0}, uint64_t{4}}) {
            auto segment = RaMesh::Segment::create(
                start,
                4,
                Strand::FORWARD,
                cigar,
                RaMesh::AlignRole::PRIMARY,
                RaMesh::SegmentRole::SEGMENT,
                block);
            block->anchors.emplace(
                RaMesh::SpeciesChrPair{species, "chr1"},
                segment);
        }
    }
    const auto hal_path =
        temp / "bilateral-copies.hal";
    std::vector<std::weak_ptr<RaMesh::Block>> blocks{block};
    RaMesh::hal_export::ExportStats stats;
    RaMesh::hal_export::exportToHal(
        blocks,
        hal_path,
        managers,
        NewickParser(
            "((leafA:0.1,leafB:0.1)anc1:0.1,"
            "leafR:0.1)anc0;"),
        "anc0",
        SoftMask::loadIndexes(softmask_paths),
        {},
        &stats);
    require(
        stats.observed_occurrence_count == 5 &&
            stats.ancestor_occurrence_count == 3,
        "two child lineages with matching bilateral context must preserve two copies in their ancestor");

    hal::AlignmentPtr alignment = hal::openHalAlignment(
        hal_path.string(), nullptr, hal::READ_ACCESS);
    require(
        static_cast<bool>(alignment),
        "cannot reopen bilateral-copy HAL");
    hal::Genome* anc1 = alignment->openGenome("anc1");
    require(
        anc1 != nullptr &&
            anc1->getNumSequences() == 1 &&
            anc1->getNumBottomSegments() == 2,
        "bilateral child context must materialize two distinct aligned ancestor occurrences");
    alignment->closeGenome(anc1);
    alignment.reset();
}

void testRecursiveOccurrenceThreading(
    const std::filesystem::path& temp,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>&
        base_managers,
    const SoftMask::PathMap& base_softmask_paths) {
    auto managers = base_managers;
    auto softmask_paths = base_softmask_paths;
    for (const std::string species :
         {"leafC", "leafD"}) {
        const auto input = temp / (species + ".input.fa");
        const auto uppercase =
            temp / (species + ".align-v2.fasta");
        const auto index =
            temp / (species + ".softmask-v1.bin");
        const auto marker =
            temp / (species + ".complete.json");
        writeInput(input, "AAAcccGGttTT");
        SoftMask::ensureUppercaseFastaAndIndex(
            input,
            uppercase,
            index,
            marker);
        softmask_paths[species] = index;
        SeqPro::ManagerVariant manager =
            std::make_unique<SeqPro::SequenceManager>(
                uppercase);
        managers[species] =
            std::make_shared<SeqPro::ManagerVariant>(
                std::move(manager));
    }

    const Cigar_t cigar{cigarToInt('M', 4)};
    auto make_block = [&](uint64_t start) {
        auto block = RaMesh::Block::create(4);
        block->ref_species = "leafA";
        block->ref_chr = "chr1";
        for (const std::string species :
             {"leafA", "leafB", "leafC", "leafD"}) {
            auto segment = RaMesh::Segment::create(
                start,
                4,
                Strand::FORWARD,
                cigar,
                RaMesh::AlignRole::PRIMARY,
                RaMesh::SegmentRole::SEGMENT,
                block);
            block->anchors.emplace(
                RaMesh::SpeciesChrPair{species, "chr1"},
                segment);
        }
        return block;
    };
    auto first = make_block(0);
    auto second = make_block(8);
    std::vector<std::weak_ptr<RaMesh::Block>> blocks{
        first,
        second};
    const auto hal_path = temp / "recursive.hal";
    RaMesh::hal_export::ExportStats stats;
    RaMesh::hal_export::exportToHal(
        blocks,
        hal_path,
        managers,
        NewickParser(
            "((leafA:0.1,leafB:0.1)anc1:0.1,"
            "(leafC:0.1,leafD:0.1)anc2:0.1)anc0;"),
        "anc0",
        SoftMask::loadIndexes(softmask_paths),
        {},
        &stats);
    require(
        stats.root_sequence_count == 1,
        "recursive occurrence threading must pack the root into one sequence");

    hal::AlignmentPtr alignment = hal::openHalAlignment(
        hal_path.string(), nullptr, hal::READ_ACCESS);
    require(
        static_cast<bool>(alignment),
        "cannot reopen recursive HAL");
    hal::Genome* root = alignment->openGenome("anc0");
    hal::Genome* anc1 = alignment->openGenome("anc1");
    hal::Genome* anc2 = alignment->openGenome("anc2");
    hal::Genome* leaf_a = alignment->openGenome("leafA");
    hal::Genome* leaf_c = alignment->openGenome("leafC");
    require(
        root != nullptr &&
            anc1 != nullptr &&
            anc2 != nullptr &&
            leaf_a != nullptr &&
            leaf_c != nullptr,
        "recursive HAL is missing expected genomes");
    std::string root_dna;
    root->getString(root_dna);
    require(
        root->getNumSequences() == 1 &&
            root->getNumBottomSegments() == 3 &&
            root_dna.size() == 18 &&
            root_dna.substr(4, 10) ==
                std::string(10, 'N'),
        "recursive root occurrences and gap must tile one sequence");
    require(
        anc1->getNumTopSegments() == 3 &&
            anc1->getNumBottomSegments() == 3 &&
            anc2->getNumTopSegments() == 3 &&
            anc2->getNumBottomSegments() == 3 &&
            leaf_a->getNumTopSegments() == 3 &&
            leaf_c->getNumTopSegments() == 3,
        "recursive parent bottoms and child tops must continuously tile every sequence");
    alignment->closeGenome(leaf_c);
    alignment->closeGenome(leaf_a);
    alignment->closeGenome(anc2);
    alignment->closeGenome(anc1);
    alignment->closeGenome(root);
    alignment.reset();
}

void testParentReferenceProjectsChildContainers(
    const std::filesystem::path& temp,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>&
        base_managers,
    const SoftMask::PathMap& base_softmask_paths) {
    SoftMask::PathMap softmask_paths;
    std::map<SpeciesName, SeqPro::SharedManagerVariant>
        managers;
    for (const std::string species :
         {"leafA", "leafB", "leafC", "leafD"}) {
        const auto input =
            temp / (species + ".projection.input.fa");
        const auto uppercase =
            temp /
            (species + ".projection.align-v2.fasta");
        const auto index =
            temp /
            (species + ".projection.softmask-v1.bin");
        const auto marker =
            temp /
            (species + ".projection.complete.json");
        writeMultiInput(
            input,
            {{"chr1", "AAAACCCC"},
             {"chr2", "AAAACCCC"}});
        SoftMask::ensureUppercaseFastaAndIndex(
            input,
            uppercase,
            index,
            marker);
        softmask_paths[species] = index;
        SeqPro::ManagerVariant manager =
            std::make_unique<SeqPro::SequenceManager>(
                uppercase);
        managers[species] =
            std::make_shared<SeqPro::ManagerVariant>(
                std::move(manager));
    }

    const Cigar_t cigar{cigarToInt('M', 4)};
    auto make_block =
        [&](const std::string& child_one_chr,
            uint64_t child_one_start,
            const std::string& child_two_chr,
            uint64_t child_two_start) {
            auto block = RaMesh::Block::create(4);
            block->ref_species = "leafA";
            block->ref_chr = child_one_chr;
            for (const std::string species :
                 {"leafA", "leafB"}) {
                auto segment =
                    RaMesh::Segment::create(
                        child_one_start,
                        4,
                        Strand::FORWARD,
                        cigar,
                        RaMesh::AlignRole::PRIMARY,
                        RaMesh::SegmentRole::SEGMENT,
                        block);
                block->anchors.emplace(
                    RaMesh::SpeciesChrPair{
                        species,
                        child_one_chr},
                    segment);
            }
            for (const std::string species :
                 {"leafC", "leafD"}) {
                auto segment =
                    RaMesh::Segment::create(
                        child_two_start,
                        4,
                        Strand::FORWARD,
                        cigar,
                        RaMesh::AlignRole::PRIMARY,
                        RaMesh::SegmentRole::SEGMENT,
                        block);
                block->anchors.emplace(
                    RaMesh::SpeciesChrPair{
                        species,
                        child_two_chr},
                    segment);
            }
            return block;
        };
    auto first = make_block(
        "chr1", 0, "chr1", 4);
    auto second = make_block(
        "chr2", 0, "chr1", 0);
    std::vector<std::weak_ptr<RaMesh::Block>>
        blocks{first, second};
    const auto hal_path =
        temp / "parent-projection.hal";
    RaMesh::hal_export::ExportStats stats;
    RaMesh::hal_export::exportToHal(
        blocks,
        hal_path,
        managers,
        NewickParser(
            "((leafA:0.1,leafB:0.1)anc1:0.1,"
            "(leafC:0.1,leafD:0.1)anc2:0.1)anc0;"),
        "anc0",
        SoftMask::loadIndexes(softmask_paths),
        {},
        &stats);

    hal::AlignmentPtr alignment =
        hal::openHalAlignment(
            hal_path.string(),
            nullptr,
            hal::READ_ACCESS);
    require(
        static_cast<bool>(alignment),
        "cannot reopen parent-projection HAL");
    hal::Genome* root =
        alignment->openGenome("anc0");
    hal::Genome* anc1 =
        alignment->openGenome("anc1");
    hal::Genome* anc2 =
        alignment->openGenome("anc2");
    require(
        root != nullptr &&
            anc1 != nullptr &&
            anc2 != nullptr,
        "parent-projection HAL is missing ancestor genomes");
    require(
        root->getNumSequences() == 1 &&
            anc1->getNumSequences() == 1 &&
            anc2->getNumSequences() == 1,
        "parent reference must replace each internal child's independent temporary container boundaries");
    require(
        root->getNumBottomSegments() == 2 &&
            anc1->getNumTopSegments() == 3 &&
            anc1->getNumBottomSegments() == 3 &&
            anc2->getNumTopSegments() == 2 &&
            anc2->getNumBottomSegments() == 2,
        "parent-projected ancestor top and bottom segments must continuously tile");
    require(
        stats.root_sequence_count == 1 &&
            stats.root_singleton_sequence_count == 0,
        "root reference must remain evidence-driven while child containers follow parent projection");
    alignment->closeGenome(anc2);
    alignment->closeGenome(anc1);
    alignment->closeGenome(root);
    alignment.reset();
}

void testNonMonotonicChildReferenceFragmentIsPreserved(
    const std::filesystem::path& temp) {
    SoftMask::PathMap softmask_paths;
    std::map<SpeciesName, SeqPro::SharedManagerVariant>
        managers;
    for (const std::string species :
         {"leafA", "leafB", "leafC", "leafD"}) {
        const auto input =
            temp /
            (species + ".rearrangement.input.fa");
        const auto uppercase =
            temp /
            (species + ".rearrangement.align-v2.fasta");
        const auto index =
            temp /
            (species + ".rearrangement.softmask-v1.bin");
        const auto marker =
            temp /
            (species + ".rearrangement.complete.json");
        writeInput(input, "AAAACCCCGGGG");
        SoftMask::ensureUppercaseFastaAndIndex(
            input,
            uppercase,
            index,
            marker);
        softmask_paths[species] = index;
        SeqPro::ManagerVariant manager =
            std::make_unique<SeqPro::SequenceManager>(
                uppercase);
        managers[species] =
            std::make_shared<SeqPro::ManagerVariant>(
                std::move(manager));
    }

    const Cigar_t cigar{cigarToInt('M', 4)};
    auto make_block =
        [&](uint64_t child_one_start,
            uint64_t child_two_start) {
            auto block = RaMesh::Block::create(4);
            block->ref_species = "leafA";
            block->ref_chr = "chr1";
            for (const std::string species :
                 {"leafA", "leafB"}) {
                auto segment =
                    RaMesh::Segment::create(
                        child_one_start,
                        4,
                        Strand::FORWARD,
                        cigar,
                        RaMesh::AlignRole::PRIMARY,
                        RaMesh::SegmentRole::SEGMENT,
                        block);
                block->anchors.emplace(
                    RaMesh::SpeciesChrPair{
                        species,
                        "chr1"},
                    segment);
            }
            for (const std::string species :
                 {"leafC", "leafD"}) {
                auto segment =
                    RaMesh::Segment::create(
                        child_two_start,
                        4,
                        Strand::FORWARD,
                        cigar,
                        RaMesh::AlignRole::PRIMARY,
                        RaMesh::SegmentRole::SEGMENT,
                        block);
                block->anchors.emplace(
                    RaMesh::SpeciesChrPair{
                        species,
                        "chr1"},
                    segment);
            }
            return block;
        };
    auto first = make_block(0, 0);
    auto second = make_block(4, 8);
    auto third = make_block(8, 4);
    std::vector<std::weak_ptr<RaMesh::Block>>
        blocks{first, second, third};
    const auto hal_path =
        temp / "parent-rearrangement.hal";
    RaMesh::hal_export::exportToHal(
        blocks,
        hal_path,
        managers,
        NewickParser(
            "((leafA:0.1,leafB:0.1)anc1:0.1,"
            "(leafC:0.1,leafD:0.1)anc2:0.1)anc0;"),
        "anc0",
        SoftMask::loadIndexes(softmask_paths));

    hal::AlignmentPtr alignment =
        hal::openHalAlignment(
            hal_path.string(),
            nullptr,
            hal::READ_ACCESS);
    require(
        static_cast<bool>(alignment),
        "cannot reopen parent-rearrangement HAL");
    hal::Genome* root =
        alignment->openGenome("anc0");
    hal::Genome* anc1 =
        alignment->openGenome("anc1");
    hal::Genome* anc2 =
        alignment->openGenome("anc2");
    require(
        root != nullptr &&
            anc1 != nullptr &&
            anc2 != nullptr,
        "parent-rearrangement HAL is missing ancestor genomes");
    std::string root_dna;
    std::string anc1_dna;
    std::string anc2_dna;
    root->getString(root_dna);
    anc1->getString(anc1_dna);
    anc2->getString(anc2_dna);
    require(
        root->getNumSequences() == 1 &&
            anc1->getNumSequences() == 1 &&
            anc2->getNumSequences() == 1 &&
            root_dna.size() == 12 &&
            anc1_dna.size() == 12 &&
            anc2_dna.size() == 12,
        "parent projection must preserve every supported rearranged child adjacency without scaffold gaps");
    require(
        anc2_dna == "AAAACCCCGGGG",
        "a child rearrangement must retain its own supported reference order");
    alignment->closeGenome(anc2);
    alignment->closeGenome(anc1);
    alignment->closeGenome(root);
    alignment.reset();
}

void testParentProjectionPreservesScaffoldedParentContainer(
    const std::filesystem::path& temp) {
    SoftMask::PathMap softmask_paths;
    std::map<SpeciesName, SeqPro::SharedManagerVariant>
        managers;
    for (const std::string species :
         {"leafA", "leafB", "leafC", "leafD"}) {
        const auto input =
            temp /
            (species + ".parent-boundary.input.fa");
        const auto uppercase =
            temp /
            (species + ".parent-boundary.align-v2.fasta");
        const auto index =
            temp /
            (species + ".parent-boundary.softmask-v1.bin");
        const auto marker =
            temp /
            (species + ".parent-boundary.complete.json");
        const std::vector<std::pair<std::string, std::string>>
            sequences =
                species == "leafC" ||
                        species == "leafD"
                    ? std::vector<std::pair<std::string, std::string>>{
                          {"chr1", "AAAACCCC"}}
                    : std::vector<std::pair<std::string, std::string>>{
                          {"chr1", "AAAA"},
                          {"chr2", "CCCC"}};
        writeMultiInput(
            input,
            sequences);
        SoftMask::ensureUppercaseFastaAndIndex(
            input,
            uppercase,
            index,
            marker);
        softmask_paths[species] = index;
        SeqPro::ManagerVariant manager =
            std::make_unique<SeqPro::SequenceManager>(
                uppercase);
        managers[species] =
            std::make_shared<SeqPro::ManagerVariant>(
                std::move(manager));
    }

    const Cigar_t cigar{cigarToInt('M', 4)};
    auto make_block =
        [&](const std::string& parent_chr,
            uint64_t child_start) {
            auto block = RaMesh::Block::create(4);
            block->ref_species = "leafA";
            block->ref_chr = parent_chr;
            for (const std::string species :
                 {"leafA", "leafB"}) {
                auto segment =
                    RaMesh::Segment::create(
                        0,
                        4,
                        Strand::FORWARD,
                        cigar,
                        RaMesh::AlignRole::PRIMARY,
                        RaMesh::SegmentRole::SEGMENT,
                        block);
                block->anchors.emplace(
                    RaMesh::SpeciesChrPair{
                        species,
                        parent_chr},
                    segment);
            }
            for (const std::string species :
                 {"leafC", "leafD"}) {
                auto segment =
                    RaMesh::Segment::create(
                        child_start,
                        4,
                        Strand::FORWARD,
                        cigar,
                        RaMesh::AlignRole::PRIMARY,
                        RaMesh::SegmentRole::SEGMENT,
                        block);
                block->anchors.emplace(
                    RaMesh::SpeciesChrPair{
                        species,
                        "chr1"},
                    segment);
            }
            return block;
        };
    auto first = make_block("chr1", 0);
    auto second = make_block("chr2", 4);
    std::vector<std::weak_ptr<RaMesh::Block>>
        blocks{first, second};
    const auto hal_path =
        temp / "parent-boundary.hal";
    RaMesh::hal_export::exportToHal(
        blocks,
        hal_path,
        managers,
        NewickParser(
            "((leafA:0.1,leafB:0.1)anc1:0.1,"
            "(leafC:0.1,leafD:0.1)anc2:0.1)anc0;"),
        "anc0",
        SoftMask::loadIndexes(softmask_paths));

    hal::AlignmentPtr alignment =
        hal::openHalAlignment(
            hal_path.string(),
            nullptr,
            hal::READ_ACCESS);
    require(
        static_cast<bool>(alignment),
        "cannot reopen parent-scaffold HAL");
    hal::Genome* root =
        alignment->openGenome("anc0");
    hal::Genome* anc2 =
        alignment->openGenome("anc2");
    require(
        root != nullptr &&
            anc2 != nullptr,
        "parent-scaffold HAL is missing ancestor genomes");
    require(
        root->getNumSequences() == 1 &&
            anc2->getNumSequences() == 1,
        "parent top-level scaffolding must project one inherited container into its child");
    alignment->closeGenome(anc2);
    alignment->closeGenome(root);
    alignment.reset();
}




}  // namespace

int main() {
    const auto temp = std::filesystem::path("/tmp") /
        ("ramax-hal-softmask-export-" + std::to_string(getpid()));
    std::filesystem::remove_all(temp);
    std::filesystem::create_directories(temp);

    try {
        const std::string leaf_a_expected = "AAAcccGGttTT";
        const std::string leaf_b_expected = "AAAcccGGTTtt";
        const std::string ancestor_expected = "AAAcccGGtttt";

        SoftMask::PathMap softmask_paths;
        std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
        for (const auto& [species, dna] :
             std::map<std::string, std::string>{{"leafA", leaf_a_expected},
                                                 {"leafB", leaf_b_expected}}) {
            const auto input = temp / (species + ".input.fa");
            const auto uppercase = temp / (species + ".align-v2.fasta");
            const auto index = temp / (species + ".softmask-v1.bin");
            const auto marker = temp / (species + ".complete.json");
            writeInput(input, dna);
            SoftMask::ensureUppercaseFastaAndIndex(input, uppercase, index, marker);
            softmask_paths[species] = index;

            SeqPro::ManagerVariant manager =
                std::make_unique<SeqPro::SequenceManager>(uppercase);
            managers[species] =
                std::make_shared<SeqPro::ManagerVariant>(std::move(manager));
        }

        RaMesh::RaMeshMultiGenomeGraph graph(managers);
        constexpr uint32_t length = 12;
        Anchor anchor(0, 0, length, 0, 0, length, Strand::FORWARD,
                      length, length, Cigar_t{cigarToInt('M', length)});
        graph.insertAnchorIntoGraph(*managers.at("leafA"), *managers.at("leafB"),
                                    "leafA", "leafB", anchor);

        const auto hal_path = temp / "exported.hal";
        NewickParser parser("(leafA:0.1,leafB:0.1)anc0;");
        graph.exportToHal(
            hal_path, managers, parser, "anc0", 1, softmask_paths);

        hal::AlignmentPtr alignment = hal::openHalAlignment(
            hal_path.string(), nullptr, hal::READ_ACCESS);
        require(static_cast<bool>(alignment), "cannot reopen exported HAL");
        const std::string leaf_a = readGenome(alignment, "leafA");
        const std::string leaf_b = readGenome(alignment, "leafB");
        const std::string ancestor = readGenome(alignment, "anc0");

        require(leaf_a == leaf_a_expected, "exported leafA mask differs from input");
        require(leaf_b == leaf_b_expected, "exported leafB mask differs from input");
        require(ancestor == ancestor_expected,
                "exported ancestor mask differs from lowercase vote");
        require(hasUpperAndLower(leaf_a), "exported leafA lacks mixed case");
        require(hasUpperAndLower(leaf_b), "exported leafB lacks mixed case");
        require(hasUpperAndLower(ancestor), "exported ancestor lacks mixed case");
        alignment.reset();
        testNamespacedLeafSequenceIsNotDuplicated(
            temp);
        testMafReassemblesGappedBlock(
            temp);
        testParalogousOccurrencesRemainThreaded(
            temp,
            managers,
            softmask_paths);
        testSecondaryLinkRunsUnifyPrimaryOccurrences(
            temp,
            managers,
            softmask_paths);
        testSecondaryMaterializationBridgesReferenceGap(
            temp,
            managers);
        testScaffoldGapSegmentsTileSequences(
            temp,
            managers,
            softmask_paths);
        testBilateralContextPreservesAncestralCopies(
            temp,
            managers,
            softmask_paths);
        testRecursiveOccurrenceThreading(
            temp,
            managers,
            softmask_paths);
        testParentReferenceProjectsChildContainers(
            temp,
            managers,
            softmask_paths);
        testNonMonotonicChildReferenceFragmentIsPreserved(
            temp);
        testParentProjectionPreservesScaffoldedParentContainer(
            temp);

        std::cout << "genome\tdna\n"
                  << "leafA\t" << leaf_a << '\n'
                  << "leafB\t" << leaf_b << '\n'
                  << "anc0\t" << ancestor << '\n';
        std::filesystem::remove_all(temp);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "hal_softmask_export_test: " << error.what() << '\n';
        std::filesystem::remove_all(temp);
        return 1;
    }
}
