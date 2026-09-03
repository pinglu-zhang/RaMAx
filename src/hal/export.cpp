#include "hal/export.h"

#include "align.h"
#include "dependency_preflight.h"
#include "external_tool.h"

#include "halAlignmentInstance.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/maximum_weighted_matching.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iterator>
#include <list>
#include <limits>
#include <map>
#include <optional>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <set>
#include <unordered_set>
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace RaMesh::hal_export {


struct PathCoverResult {
    std::vector<std::vector<OccurrenceId>> paths;
    std::unordered_map<OccurrenceId, bool> forward_by_occurrence;
};


PathCoverResult buildMaximumCardinalityWeightPathCoverDetailedImpl(
    const std::vector<uint64_t>& run_ids,
    const std::vector<EdgeSupport>& edges,
    const std::unordered_map<uint64_t, RunOrderKey>& run_order_keys,
    ExportStats* stats);



namespace {


struct LeafRow {
    std::string row_id;
    std::string leaf_name;
    std::string chr_name;
    // Internal species-qualified key.  The public HAL leaf sequence name is
    // the original chr_name and is assigned only when emissions are written.
    std::string hal_sequence_name;
    uint64_t segment_start = 0;
    uint32_t segment_length = 0;
    bool reversed = false;
    std::string aligned;
};

struct BlockMSA {
    BlockPtr block;
    uint64_t block_id = 0;
    std::string ref_row_id;
    size_t alignment_length = 0;
    bool secondary_homology = false;
    ExportBlockOrderKey order_key;
    std::vector<LeafRow> leaf_rows;
};

struct LeafRunSpan {
    std::string row_id;
    std::string leaf_name;
    std::string chr_name;
    // Internal species-qualified key; never emit this as the HAL sequence
    // name because HAL sequence names are scoped by genome.
    std::string hal_sequence_name;
    uint64_t start = 0;
    uint32_t length = 0;
    bool reversed = false;
    std::string dna;
};

struct LeafOccurrence {
    uint64_t run_id = 0;
    LeafRunSpan span;
    bool forward_to_canonical = true;
    uint32_t copy_index = 0;
};

struct ColumnRun {
    uint64_t run_id = 0;
    uint64_t block_id = 0;
    uint32_t col_beg = 0;
    uint32_t col_end = 0;
    bool secondary_homology = false;
    std::vector<uint8_t> leaf_present;
    std::vector<LeafRunSpan> leaf_spans;
    std::vector<uint64_t> source_block_ids;
    std::vector<uint8_t> present_by_node;
    std::vector<double> presence_margin;
};

struct OrientedOccurrence {
    OccurrenceId occurrence_id = 0;
    uint64_t run_id = 0;
    uint32_t copy_index = 0;
    bool forward = true;
};
struct SequenceGap {
    uint64_t start = 0;
    uint32_t length = 0;
};


struct SequenceModel {
    std::string seq_name;
    std::string dna;
    std::vector<OrientedOccurrence> path;
    std::vector<ReferenceJoin> joins;
    std::vector<SequenceGap> gaps;
};

struct NodePlacement {
    std::string seq_name;
    uint64_t start = 0;
    uint32_t length = 0;
    bool forward = true;
};

using ChildRunCopyMap = std::unordered_map<
    int,
    std::unordered_map<
        uint64_t,
        std::vector<uint32_t>>>;


struct NodeModel {
    std::string genome_name;
    std::vector<SequenceModel> sequences;
    std::unordered_map<uint64_t, std::string> run_dna;
    std::unordered_map<OccurrenceId, uint64_t> run_by_occurrence;
    std::unordered_map<uint64_t, std::vector<OccurrenceId>> occurrences_by_run;
    std::unordered_map<OccurrenceId, uint32_t> copy_index_by_occurrence;
    ChildRunCopyMap parent_copy_by_child_run;
    std::unordered_map<OccurrenceId, NodePlacement> placements;
    std::vector<TerminalEndSupport> terminal_ends;
};

struct TopSegmentLine {
    OccurrenceId occurrence_id = 0;
    uint64_t start = 0;
    uint32_t length = 0;
    std::optional<uint64_t> parent_bottom_name;
    bool forward_to_parent = true;
};

struct BottomSegmentLine {
    OccurrenceId occurrence_id = 0;
    uint64_t name = 0;
    uint64_t start = 0;
    uint32_t length = 0;
};

struct SequenceEmission {
    std::string genome_name;
    std::string seq_name;
    size_t bottom_count = 0;
    std::vector<BottomSegmentLine> bottoms;
    std::vector<TopSegmentLine> tops;
    std::optional<std::string> dna;
    std::optional<std::pair<std::string, std::string>> leaf_source;
};


struct ChildEdgeContribution {
    uint64_t from = 0;
    uint64_t to = 0;
    bool from_forward_to_canonical = true;
    bool to_forward_to_canonical = true;
    uint32_t occurrence_support = 0;
    uint32_t left_length = 0;
    uint32_t right_length = 0;
    uint64_t minimum_gap = std::numeric_limits<uint64_t>::max();
    long double weighted_support = 0.0L;
    int child_id = -1;
};

struct FastInferenceResult {
    std::vector<uint8_t> present_by_node;
    std::vector<double> margin;
};

struct NodeModelBuildTimings {
    size_t candidate_occurrence_count = 0;
    size_t edge_count = 0;
    size_t path_count = 0;
    uint64_t candidate_filter_ms = 0;
    uint64_t run_dna_ms = 0;
    uint64_t edge_collect_ms = 0;
    uint64_t path_decompose_ms = 0;
    uint64_t sequence_materialize_ms = 0;
    uint64_t total_ms = 0;
};

struct RunPairKey {
    uint64_t from = 0;
    uint64_t to = 0;
    bool from_forward_to_canonical = true;
    bool to_forward_to_canonical = true;

    bool operator==(const RunPairKey& other) const {
        return from == other.from &&
               to == other.to &&
               from_forward_to_canonical == other.from_forward_to_canonical &&
               to_forward_to_canonical == other.to_forward_to_canonical;
    }
};

struct RunPairKeyHash {
    size_t operator()(const RunPairKey& key) const noexcept {
        size_t seed = std::hash<uint64_t>{}(key.from);
        seed ^= std::hash<uint64_t>{}(key.to) + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
        seed ^= std::hash<bool>{}(key.from_forward_to_canonical) +
                0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
        seed ^= std::hash<bool>{}(key.to_forward_to_canonical) +
                0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
        return seed;
    }
};


template <typename VariantLike>
std::string fetchSubSequence(const VariantLike& shared_mv,
                             const std::string& chr,
                             uint64_t start,
                             uint64_t length) {
    return std::visit([&](const auto& mgr) -> std::string {
        using PtrType = std::decay_t<decltype(mgr)>;
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
            auto chr_id = mgr->getSequenceId(chr);
            return mgr->getSubSequence(chr_id, start, length);
        } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            auto chr_id = mgr->getSequenceId(chr);
            return mgr->getOriginalManager().getSubSequence(chr_id, start, length);
        } else {
            throw std::runtime_error("Unsupported SeqPro manager variant");
        }
    }, *shared_mv);
}

template <typename VariantLike>
uint64_t fetchSequenceLength(const VariantLike& shared_mv, const std::string& chr) {
    return std::visit([&](const auto& mgr) -> uint64_t {
        using PtrType = std::decay_t<decltype(mgr)>;
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
            return mgr->getSequenceLength(chr);
        } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            return mgr->getOriginalManager().getSequenceLength(chr);
        } else {
            throw std::runtime_error("Unsupported SeqPro manager variant");
        }
    }, *shared_mv);
}

template <typename VariantLike>
std::vector<std::string> fetchSequenceNames(const VariantLike& shared_mv) {
    return std::visit([&](const auto& mgr) -> std::vector<std::string> {
        return mgr->getSequenceNames();
    }, *shared_mv);
}

std::string makeQualifiedLeafSequenceName(
    const std::string& species,
    const std::string& chr) {
    const std::string prefix =
        species + ".";
    if (chr.starts_with(prefix)) {
        return chr;
    }
    return prefix + chr;
}

std::string makeLeafRowId(const std::string& species, const std::string& chr) {
    return species + '\x1f' + chr;
}

bool normalizeRootNode(NewickParser& parser, const std::string& requested_root_name) {
    const auto& nodes = parser.getNodes();
    int root_id = -1;
    for (const auto& node : nodes) {
        if (node.father != -1) {
            continue;
        }
        if (root_id != -1) {
            throw std::runtime_error("HAL export phylogeny contains multiple roots");
        }
        root_id = node.id;
    }
    if (root_id == -1) {
        throw std::runtime_error("HAL export phylogeny has no root");
    }

    if (nodes[static_cast<size_t>(root_id)].isLeaf) {
        NewickTreeNode artificial_root;
        artificial_root.id = parser.currentIndex_++;
        artificial_root.name =
            requested_root_name.empty() ? std::string("ancestor") : requested_root_name;
        artificial_root.father = -1;
        artificial_root.isLeaf = false;
        artificial_root.branchLength = 0.0;
        artificial_root.leftChild = root_id;
        artificial_root.rightChild = -1;
        parser.nodes_[static_cast<size_t>(root_id)].father = artificial_root.id;
        parser.nodes_[static_cast<size_t>(root_id)].branchLength = 1.0;
        parser.nodes_.push_back(std::move(artificial_root));
        return true;
    }

    if (nodes[static_cast<size_t>(root_id)].name.empty()) {
        parser.nodes_[static_cast<size_t>(root_id)].name =
            requested_root_name.empty() ? std::string("ancestor") : requested_root_name;
        return true;
    }
    return false;
}

void validateLeafNamesExact(
    const NewickParser& parser,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {
    std::vector<std::string> tree_leaves = parser.getLeafNames();
    std::vector<std::string> manager_leaves;
    manager_leaves.reserve(seqpro_managers.size());
    for (const auto& [species_name, manager] : seqpro_managers) {
        (void)manager;
        manager_leaves.push_back(species_name);
    }
    std::sort(tree_leaves.begin(), tree_leaves.end());
    std::sort(manager_leaves.begin(), manager_leaves.end());
    if (tree_leaves == manager_leaves) {
        return;
    }

    std::ostringstream message;
    message << "HAL export leaf names differ between the phylogeny and sequence managers";
    std::vector<std::string> missing_in_tree;
    std::vector<std::string> missing_in_managers;
    std::set_difference(
        manager_leaves.begin(),
        manager_leaves.end(),
        tree_leaves.begin(),
        tree_leaves.end(),
        std::back_inserter(missing_in_tree));
    std::set_difference(
        tree_leaves.begin(),
        tree_leaves.end(),
        manager_leaves.begin(),
        manager_leaves.end(),
        std::back_inserter(missing_in_managers));
    for (const auto& name : missing_in_tree) {
        message << "; missing in tree: " << name;
    }
    for (const auto& name : missing_in_managers) {
        message << "; missing sequence manager: " << name;
    }
    throw std::runtime_error(message.str());
}

std::string stripGaps(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (c != '-') {
            out.push_back(c);
        }
    }
    return out;
}

std::string formatScaffoldName(const std::string& genome_name, size_t index) {
    return genome_name + "refChr" + std::to_string(index);
}

std::vector<int> computePostorderInternal(const TreeMeta& tree, int node_id) {
    std::vector<int> out;
    std::function<void(int)> dfs = [&](int cur) {
        for (int child : tree.nodes[cur].children) {
            dfs(child);
        }
        if (!tree.nodes[cur].is_leaf) {
            out.push_back(cur);
        }
    };
    dfs(node_id);
    return out;
}

FastInferenceResult inferDescendantUnionFast(
    const TreeMeta& tree,
    const std::vector<uint8_t>& leaf_present) {
    if (leaf_present.size() != tree.leaf_ids.size()) {
        throw std::runtime_error("Descendant-union inference received a malformed leaf state vector");
    }

    FastInferenceResult result;
    result.present_by_node.assign(tree.nodes.size(), 0);
    result.margin.assign(tree.nodes.size(), 1.0);
    for (int leaf_id : tree.leaf_ids) {
        const int leaf_index = tree.nodes[leaf_id].leaf_index;
        result.present_by_node[leaf_id] = leaf_present[static_cast<size_t>(leaf_index)];
    }
    for (int node_id : tree.internal_postorder) {
        const auto& node = tree.nodes[node_id];
        result.present_by_node[node_id] = std::any_of(
            node.children.begin(), node.children.end(), [&](int child_id) {
                return result.present_by_node[child_id] != 0;
            });
    }
    return result;
}

BlockMSA buildBlockMSA(
    const BlockPtr& block,
    const TreeMeta& tree,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
    const SoftMask::IndexMap* softmask_indexes) {
    if (!block) {
        throw std::runtime_error("Null block passed to buildBlockMSA");
    }

    std::unordered_map<std::string, std::string> sequences;
    // Alignment decisions must not depend on output-only soft-mask case.
    // Keep an oriented case-preserving copy and overlay it only after the
    // uppercase rows have been aligned.
    std::unordered_map<std::string, std::string> output_sequences;
    std::unordered_map<ChrName, Cigar_t> cigars;
    std::unordered_map<std::string, LeafRow> rows;
    std::vector<std::pair<SpeciesChrPair, SegPtr>> occurrences;
    bool secondary_homology = false;

    {
        std::shared_lock block_lock(block->rw);
        occurrences.reserve(block->anchors.size());
        for (const auto& [species_chr, segment] : block->anchors) {
            const auto node_it = tree.name_to_id.find(species_chr.first);
            if (node_it == tree.name_to_id.end() || !tree.nodes[node_it->second].is_leaf) {
                continue;
            }
            secondary_homology =
                secondary_homology ||
                !segment->isPrimary();
            occurrences.emplace_back(species_chr, segment);
        }
    }
    std::sort(occurrences.begin(), occurrences.end(),
              [](const auto &lhs, const auto &rhs) {
                const auto &lhs_key = lhs.first;
                const auto &rhs_key = rhs.first;
                const auto &lhs_segment = lhs.second;
                const auto &rhs_segment = rhs.second;
                return std::tie(lhs_key.first, lhs_key.second,
                                lhs_segment->start, lhs_segment->length,
                                lhs_segment->strand, lhs_segment->align_role) <
                       std::tie(rhs_key.first, rhs_key.second,
                                rhs_segment->start, rhs_segment->length,
                                rhs_segment->strand, rhs_segment->align_role);
              });

    std::unordered_map<std::string, uint32_t> next_occurrence;
    std::vector<std::string> reference_rows;
    for (const auto &[species_chr, segment] : occurrences) {
      auto mgr_it = seqpro_managers.find(species_chr.first);
      if (mgr_it == seqpro_managers.end()) {
        continue;
      }

      const std::string sequence_key =
          makeLeafRowId(species_chr.first, species_chr.second);
      const uint32_t occurrence = next_occurrence[sequence_key]++;
      const std::string row_id =
          sequence_key + '\1' + std::to_string(occurrence);
      std::string dna = fetchSubSequence(mgr_it->second, species_chr.second,
                                         segment->start, segment->length);
      std::string output_dna;
      if (softmask_indexes != nullptr) {
        const auto softmask_it =
            softmask_indexes->find(species_chr.first);
        if (softmask_it == softmask_indexes->end() ||
            !softmask_it->second) {
          throw std::runtime_error("Missing soft-mask index for leaf genome: " +
                                   species_chr.first);
        }
        output_dna = dna;
        softmask_it->second->restore(
            species_chr.second, segment->start, output_dna);
      }
      if (segment->strand == Strand::REVERSE) {
        reverseComplement(dna);
        if (!output_dna.empty()) {
          reverseComplement(output_dna);
        }
      }
      sequences.emplace(row_id, std::move(dna));
      if (!output_dna.empty()) {
        output_sequences.emplace(
            row_id, std::move(output_dna));
      }
      cigars.emplace(row_id, segment->cigar);

      LeafRow row;
      row.row_id = row_id;
      row.leaf_name = species_chr.first;
      row.chr_name = species_chr.second;
      row.hal_sequence_name =
          makeQualifiedLeafSequenceName(species_chr.first, species_chr.second);
      row.segment_start = segment->start;
      row.segment_length = segment->length;
      row.reversed = (segment->strand == Strand::REVERSE);
      rows.emplace(row_id, std::move(row));
      if (species_chr.first == block->ref_species &&
          species_chr.second == block->ref_chr && segment->isPrimary()) {
        reference_rows.push_back(row_id);
      }
    }

    if (sequences.empty()) {
      return {};
    }
    if (reference_rows.size() != 1) {
      throw std::runtime_error(
          "HAL export failed: block " + std::to_string(block->block_id) +
          " must contain exactly one primary declared reference occurrence '" +
          block->ref_species + "." + block->ref_chr + "', found " +
          std::to_string(reference_rows.size()));
    }
    const std::string &ref_row_id = reference_rows.front();
    const auto ref_row_it = rows.find(ref_row_id);

    try {
      mergeAlignmentByRef(ref_row_id, sequences, cigars);
    } catch (const std::exception &e) {
      std::ostringstream oss;
      oss << "HAL export failed: mergeAlignmentByRef on block "
          << block->block_id << " (" << block->ref_species << ","
          << block->ref_chr << ") threw: " << e.what();
      throw std::runtime_error(oss.str());
    }

    BlockMSA msa;
    msa.block = block;
    msa.block_id = block->block_id;
    msa.ref_row_id = ref_row_id;
    msa.alignment_length = sequences.at(ref_row_id).size();
    msa.secondary_homology = secondary_homology;
    msa.order_key.ref_species = block->ref_species;
    msa.order_key.ref_chr = block->ref_chr;
    msa.order_key.ref_start = ref_row_it->second.segment_start;
    msa.order_key.block_id = block->block_id;
    msa.leaf_rows.reserve(rows.size());

    for (auto& [row_id, aligned] : sequences) {
        if (aligned.size() != msa.alignment_length) {
            throw std::runtime_error(
                "HAL export failed: inconsistent aligned row length in block " +
                std::to_string(block->block_id));
        }
        const auto output_it =
            output_sequences.find(row_id);
        if (output_it != output_sequences.end()) {
            size_t source_index = 0;
            for (char& base : aligned) {
                if (base == '-') {
                    continue;
                }
                if (source_index >=
                    output_it->second.size()) {
                    throw std::runtime_error(
                        "HAL export failed: aligned occurrence exceeds "
                        "its soft-mask source");
                }
                base = output_it->second[
                    source_index++];
            }
            if (source_index !=
                output_it->second.size()) {
                throw std::runtime_error(
                    "HAL export failed: aligned occurrence did not consume "
                    "its soft-mask source");
            }
        }
        auto row_it = rows.find(row_id);
        if (row_it == rows.end()) {
            throw std::runtime_error("HAL export failed: aligned occurrence lost its row metadata");
        }
        row_it->second.aligned = std::move(aligned);
        msa.leaf_rows.push_back(std::move(row_it->second));
    }
    std::sort(msa.leaf_rows.begin(), msa.leaf_rows.end(), [](const LeafRow& lhs, const LeafRow& rhs) {
        return lhs.row_id < rhs.row_id;
    });
    return msa;
}

std::vector<ColumnRun> buildColumnRuns(const std::vector<BlockMSA>& block_msas, const TreeMeta& tree) {
    std::vector<ColumnRun> runs;
    uint64_t next_run_id = 1;
    for (const auto& msa : block_msas) {
        if (!msa.block || msa.alignment_length == 0) {
            continue;
        }

        std::vector<AlignedOccurrence> aligned_rows;
        aligned_rows.reserve(msa.leaf_rows.size());
        for (const auto& row : msa.leaf_rows) {
            aligned_rows.push_back(AlignedOccurrence{
                row.row_id,
                row.leaf_name,
                row.chr_name,
                row.segment_start,
                row.segment_length,
                row.reversed,
                row.aligned});
        }

        for (auto& projection : projectElementaryRuns(aligned_rows)) {
            ColumnRun run;
            run.run_id = next_run_id++;
            run.block_id = msa.block_id;
            run.source_block_ids.push_back(
                msa.block_id);
            run.col_beg = projection.col_beg;
            run.col_end = projection.col_end;
            run.secondary_homology =
                msa.secondary_homology;
            run.leaf_present.assign(tree.leaf_ids.size(), 0);
            run.leaf_spans.reserve(projection.occurrences.size());

            for (auto& occurrence : projection.occurrences) {
                const auto node_it = tree.name_to_id.find(occurrence.genome_name);
                if (node_it == tree.name_to_id.end()) {
                    throw std::runtime_error("HAL export failed: projected leaf is absent from the phylogeny");
                }
                const auto& leaf_node = tree.nodes[node_it->second];
                LeafRunSpan span;
                span.row_id = std::move(occurrence.row_id);
                span.leaf_name = std::move(occurrence.genome_name);
                span.chr_name = std::move(occurrence.sequence_name);
                span.hal_sequence_name = makeQualifiedLeafSequenceName(span.leaf_name, span.chr_name);
                span.start = occurrence.start;
                span.length = occurrence.length;
                span.reversed = occurrence.reversed;
                span.dna = std::move(occurrence.dna);
                run.leaf_spans.push_back(std::move(span));
                run.leaf_present[static_cast<size_t>(leaf_node.leaf_index)] = 1;
            }
            runs.push_back(std::move(run));
        }
    }
    return runs;
}

struct RunCoordinateTransform {
    uint32_t root = 0;
    int8_t sign = 1;
    int64_t shift = 0;
};

class RunCoordinateUnion {
public:
    explicit RunCoordinateUnion(size_t count)
        : parent_(count),
          rank_(count, 0),
          component_size_(count, 1),
          sign_to_parent_(count, 1),
          shift_to_parent_(count, 0) {
        if (count > static_cast<size_t>(
                        std::numeric_limits<uint32_t>::max())) {
            throw std::overflow_error(
                "Too many HAL column runs for coordinate unification");
        }
        std::iota(parent_.begin(), parent_.end(), 0U);
    }

    RunCoordinateTransform find(uint32_t node) {
        const uint32_t parent = parent_.at(node);
        if (parent == node) {
            return {node, 1, 0};
        }
        const auto parent_transform = find(parent);
        const int8_t old_sign = sign_to_parent_[node];
        const int64_t old_shift = shift_to_parent_[node];
        parent_[node] = parent_transform.root;
        sign_to_parent_[node] = static_cast<int8_t>(
            parent_transform.sign * old_sign);
        shift_to_parent_[node] =
            static_cast<int64_t>(parent_transform.sign) *
                old_shift +
            parent_transform.shift;
        return {
            parent_[node],
            sign_to_parent_[node],
            shift_to_parent_[node]};
    }

    bool unite(
        uint32_t first,
        uint32_t second,
        int8_t second_sign_from_first,
        int64_t second_shift_from_first) {
        const auto first_transform = find(first);
        const auto second_transform = find(second);
        if (first_transform.root == second_transform.root) {
            const int8_t expected_sign = static_cast<int8_t>(
                second_transform.sign * second_sign_from_first);
            const int64_t expected_shift =
                static_cast<int64_t>(second_transform.sign) *
                    second_shift_from_first +
                second_transform.shift;
            if (first_transform.sign != expected_sign ||
                first_transform.shift != expected_shift) {
                throw std::runtime_error(
                    "HAL run coordinate unification found an "
                    "inconsistent homology cycle");
            }
            return false;
        }

        const uint32_t first_root = first_transform.root;
        const uint32_t second_root = second_transform.root;
        if (rank_[first_root] >= rank_[second_root]) {
            const int8_t root_sign = static_cast<int8_t>(
                first_transform.sign *
                second_sign_from_first *
                second_transform.sign);
            const int64_t root_shift =
                first_transform.shift -
                static_cast<int64_t>(
                    first_transform.sign *
                    second_sign_from_first) *
                    second_shift_from_first -
                static_cast<int64_t>(root_sign) *
                    second_transform.shift;
            parent_[second_root] = first_root;
            sign_to_parent_[second_root] = root_sign;
            shift_to_parent_[second_root] = root_shift;
            component_size_[first_root] +=
                component_size_[second_root];
            if (rank_[first_root] == rank_[second_root]) {
                ++rank_[first_root];
            }
        } else {
            const int8_t root_sign = static_cast<int8_t>(
                second_transform.sign *
                second_sign_from_first *
                first_transform.sign);
            const int64_t root_shift =
                static_cast<int64_t>(second_transform.sign) *
                    second_shift_from_first +
                second_transform.shift -
                static_cast<int64_t>(root_sign) *
                    first_transform.shift;
            parent_[first_root] = second_root;
            sign_to_parent_[first_root] = root_sign;
            shift_to_parent_[first_root] = root_shift;
            component_size_[second_root] +=
                component_size_[first_root];
        }
        return true;
    }

    uint32_t componentSize(uint32_t root) const {
        return component_size_.at(root);
    }

private:
    std::vector<uint32_t> parent_;
    std::vector<uint8_t> rank_;
    std::vector<uint32_t> component_size_;
    std::vector<int8_t> sign_to_parent_;
    std::vector<int64_t> shift_to_parent_;
};


class RollbackRunConnectivity {
public:
    explicit RollbackRunConnectivity(size_t count)
        : parent_(count),
          component_size_(count, 1) {
        std::iota(parent_.begin(), parent_.end(), 0U);
    }

    size_t snapshot() const noexcept {
        return history_.size();
    }

    uint32_t find(uint32_t node) const {
        while (parent_.at(node) != node) {
            node = parent_[node];
        }
        return node;
    }

    void unite(uint32_t first, uint32_t second) {
        first = find(first);
        second = find(second);
        if (first == second) {
            return;
        }
        if (component_size_[first] <
            component_size_[second]) {
            std::swap(first, second);
        }
        history_.push_back(
            Change{
                second,
                first,
                component_size_[first]});
        parent_[second] = first;
        component_size_[first] +=
            component_size_[second];
    }

    void rollback(size_t snapshot) {
        while (history_.size() > snapshot) {
            const Change change = history_.back();
            history_.pop_back();
            parent_[change.child] = change.child;
            component_size_[change.parent] =
                change.parent_size;
        }
    }

private:
    struct Change {
        uint32_t child = 0;
        uint32_t parent = 0;
        uint32_t parent_size = 0;
    };

    std::vector<uint32_t> parent_;
    std::vector<uint32_t> component_size_;
    std::vector<Change> history_;
};

class RollbackRunCoordinateUnion {
public:
    explicit RollbackRunCoordinateUnion(size_t count)
        : parent_(count),
          component_size_(count, 1),
          sign_to_parent_(count, 1),
          shift_to_parent_(count, 0) {
        if (count > static_cast<size_t>(
                        std::numeric_limits<uint32_t>::max())) {
            throw std::overflow_error(
                "Too many HAL runs for coordinate selection");
        }
        std::iota(parent_.begin(), parent_.end(), 0U);
    }

    size_t snapshot() const noexcept {
        return history_.size();
    }

    RunCoordinateTransform find(uint32_t node) const {
        int8_t sign = 1;
        int64_t shift = 0;
        while (parent_.at(node) != node) {
            const int8_t edge_sign =
                sign_to_parent_[node];
            shift =
                static_cast<int64_t>(edge_sign) *
                    shift +
                shift_to_parent_[node];
            sign = static_cast<int8_t>(
                edge_sign * sign);
            node = parent_[node];
        }
        return {node, sign, shift};
    }

    bool unite(
        uint32_t first,
        uint32_t second,
        int8_t second_sign_from_first,
        int64_t second_shift_from_first) {
        const auto first_transform = find(first);
        const auto second_transform = find(second);
        if (first_transform.root ==
            second_transform.root) {
            return first_transform.sign ==
                       static_cast<int8_t>(
                           second_transform.sign *
                           second_sign_from_first) &&
                   first_transform.shift ==
                       static_cast<int64_t>(
                           second_transform.sign) *
                           second_shift_from_first +
                           second_transform.shift;
        }

        uint32_t first_root =
            first_transform.root;
        uint32_t second_root =
            second_transform.root;
        if (component_size_[first_root] >=
            component_size_[second_root]) {
            const int8_t root_sign =
                static_cast<int8_t>(
                    first_transform.sign *
                    second_sign_from_first *
                    second_transform.sign);
            const int64_t root_shift =
                first_transform.shift -
                static_cast<int64_t>(
                    first_transform.sign *
                    second_sign_from_first) *
                    second_shift_from_first -
                static_cast<int64_t>(root_sign) *
                    second_transform.shift;
            history_.push_back(
                Change{
                    second_root,
                    first_root,
                    component_size_[first_root],
                    sign_to_parent_[second_root],
                    shift_to_parent_[second_root]});
            parent_[second_root] = first_root;
            sign_to_parent_[second_root] =
                root_sign;
            shift_to_parent_[second_root] =
                root_shift;
            component_size_[first_root] +=
                component_size_[second_root];
        } else {
            const int8_t root_sign =
                static_cast<int8_t>(
                    second_transform.sign *
                    second_sign_from_first *
                    first_transform.sign);
            const int64_t root_shift =
                static_cast<int64_t>(
                    second_transform.sign) *
                    second_shift_from_first +
                second_transform.shift -
                static_cast<int64_t>(root_sign) *
                    first_transform.shift;
            history_.push_back(
                Change{
                    first_root,
                    second_root,
                    component_size_[second_root],
                    sign_to_parent_[first_root],
                    shift_to_parent_[first_root]});
            parent_[first_root] = second_root;
            sign_to_parent_[first_root] =
                root_sign;
            shift_to_parent_[first_root] =
                root_shift;
            component_size_[second_root] +=
                component_size_[first_root];
        }
        return true;
    }

    void rollback(size_t snapshot) {
        while (history_.size() > snapshot) {
            const Change change = history_.back();
            history_.pop_back();
            parent_[change.child] =
                change.child;
            sign_to_parent_[change.child] =
                change.child_sign;
            shift_to_parent_[change.child] =
                change.child_shift;
            component_size_[change.parent] =
                change.parent_size;
        }
    }

private:
    struct Change {
        uint32_t child = 0;
        uint32_t parent = 0;
        uint32_t parent_size = 0;
        int8_t child_sign = 1;
        int64_t child_shift = 0;
    };

    std::vector<uint32_t> parent_;
    std::vector<uint32_t> component_size_;
    std::vector<int8_t> sign_to_parent_;
    std::vector<int64_t> shift_to_parent_;
    std::vector<Change> history_;
};
struct RunNormalizationStats {
    size_t input_runs = 0;
    size_t output_runs = 0;
    size_t overlap_constraints = 0;
    size_t connected_components = 0;
    size_t duplicate_leaf_occurrences = 0;
};

int64_t checkedSignedCoordinate(uint64_t coordinate) {
    if (coordinate >
        static_cast<uint64_t>(
            std::numeric_limits<int64_t>::max())) {
        throw std::overflow_error(
            "HAL leaf coordinate exceeds signed normalization range");
    }
    return static_cast<int64_t>(coordinate);
}

int64_t runColumnLength(const ColumnRun& run) {
    if (run.col_end <= run.col_beg) {
        throw std::runtime_error(
            "HAL column run has an empty column interval");
    }
    return static_cast<int64_t>(run.col_end - run.col_beg);
}

int64_t leafBoundaryOffset(
    const LeafRunSpan& span,
    uint64_t coordinate) {
    const uint64_t end = span.start + span.length;
    if (end < span.start ||
        coordinate < span.start ||
        coordinate > end) {
        throw std::runtime_error(
            "HAL leaf overlap boundary lies outside its run span");
    }
    const uint64_t offset = span.reversed
        ? end - coordinate
        : coordinate - span.start;
    return checkedSignedCoordinate(offset);
}

struct SecondaryRunSelectionStats {
    size_t candidate_runs = 0;
    size_t accepted_runs = 0;
    size_t rejected_runs = 0;
    size_t conflict_rejected_runs = 0;
    std::unordered_set<uint64_t>
        rejected_block_ids;
    size_t redundant_runs = 0;
    uint64_t rejected_bases = 0;
    size_t primary_overlap_constraints = 0;
};

SecondaryRunSelectionStats
selectCoordinateConsistentSecondaryRuns(
    std::vector<ColumnRun>& runs) {
    SecondaryRunSelectionStats stats;
    if (runs.empty()) {
        return stats;
    }

    struct PrimaryWindow {
        uint64_t start = 0;
        uint64_t end = 0;
        uint32_t run_index = 0;
        uint32_t span_index = 0;
    };
    struct PrimaryConstraint {
        uint32_t primary_run = 0;
        uint32_t secondary_run = 0;
        int8_t secondary_sign_from_primary = 1;
        int64_t secondary_shift_from_primary = 0;
    };

    std::unordered_map<
        std::string,
        std::vector<PrimaryWindow>>
        primary_windows_by_sequence;
    std::vector<uint32_t> secondary_runs;
    for (uint32_t run_index = 0;
         run_index < runs.size();
         ++run_index) {
        const auto& run = runs[run_index];
        const int64_t length = runColumnLength(run);
        if (run.secondary_homology) {
            secondary_runs.push_back(run_index);
            continue;
        }
        if (run.leaf_spans.size() >
            static_cast<size_t>(
                std::numeric_limits<uint32_t>::max())) {
            throw std::overflow_error(
                "Too many leaf spans in a primary HAL run");
        }
        for (uint32_t span_index = 0;
             span_index < run.leaf_spans.size();
             ++span_index) {
            const auto& span =
                run.leaf_spans[span_index];
            if (span.length != length) {
                throw std::runtime_error(
                    "Primary HAL run span length is inconsistent");
            }
            const uint64_t end =
                span.start + span.length;
            if (end < span.start) {
                throw std::overflow_error(
                    "Primary HAL leaf span coordinate overflow");
            }
            primary_windows_by_sequence[
                span.hal_sequence_name]
                .push_back(
                    PrimaryWindow{
                        span.start,
                        end,
                        run_index,
                        span_index});
        }
    }
    stats.candidate_runs = secondary_runs.size();
    if (secondary_runs.empty()) {
        return stats;
    }

    for (auto& [sequence_name, windows] :
         primary_windows_by_sequence) {
        std::sort(
            windows.begin(),
            windows.end(),
            [](const PrimaryWindow& lhs,
               const PrimaryWindow& rhs) {
                return std::tie(
                           lhs.start,
                           lhs.end,
                           lhs.run_index,
                           lhs.span_index) <
                       std::tie(
                           rhs.start,
                           rhs.end,
                           rhs.run_index,
                           rhs.span_index);
            });
        uint64_t covered_until = 0;
        bool have_window = false;
        for (const auto& window : windows) {
            if (have_window &&
                window.start < covered_until) {
                throw std::runtime_error(
                    "HAL primary runs overlap on leaf sequence " +
                    sequence_name);
            }
            covered_until = window.end;
            have_window = true;
        }
    }

    std::unordered_map<
        uint32_t,
        std::vector<PrimaryConstraint>>
        constraints_by_secondary;
    constraints_by_secondary.reserve(
        secondary_runs.size());
    std::vector<uint8_t> secondary_has_novel_span(
        runs.size(),
        0);
    for (uint32_t secondary_run :
         secondary_runs) {
        const auto& run = runs[secondary_run];
        const int64_t length = runColumnLength(run);
        for (const auto& secondary_span :
             run.leaf_spans) {
            bool overlaps_primary = false;
            if (secondary_span.length != length) {
                throw std::runtime_error(
                    "Secondary HAL run span length is inconsistent");
            }
            const auto windows_it =
                primary_windows_by_sequence.find(
                    secondary_span
                        .hal_sequence_name);
            if (windows_it ==
                primary_windows_by_sequence.end()) {
                secondary_has_novel_span[
                    secondary_run] = 1;
                continue;
            }
            const uint64_t secondary_end =
                secondary_span.start +
                secondary_span.length;
            const auto& windows =
                windows_it->second;
            auto window_it = std::partition_point(
                windows.begin(),
                windows.end(),
                [&](const PrimaryWindow& window) {
                    return window.end <=
                           secondary_span.start;
                });
            while (window_it != windows.end() &&
                   window_it->start <
                       secondary_end) {
                overlaps_primary = true;
                const auto& primary_span =
                    runs[window_it->run_index]
                        .leaf_spans[
                            window_it->span_index];
                const uint64_t shared_coordinate =
                    std::max(
                        window_it->start,
                        secondary_span.start);
                const int64_t primary_offset =
                    leafBoundaryOffset(
                        primary_span,
                        shared_coordinate);
                const int64_t secondary_offset =
                    leafBoundaryOffset(
                        secondary_span,
                        shared_coordinate);
                const int8_t sign =
                    primary_span.reversed ==
                            secondary_span.reversed
                        ? 1
                        : -1;
                constraints_by_secondary[
                    secondary_run]
                    .push_back(
                        PrimaryConstraint{
                            window_it->run_index,
                            secondary_run,
                            sign,
                            secondary_offset -
                                static_cast<int64_t>(
                                    sign) *
                                    primary_offset});
                ++stats.primary_overlap_constraints;
                ++window_it;
            }
            if (!overlaps_primary) {
                secondary_has_novel_span[
                    secondary_run] = 1;
            }
        }
    }

    for (auto& [secondary_run, constraints] :
         constraints_by_secondary) {
        (void)secondary_run;
        std::sort(
            constraints.begin(),
            constraints.end(),
            [](const PrimaryConstraint& lhs,
               const PrimaryConstraint& rhs) {
                return std::tie(
                           lhs.primary_run,
                           lhs.secondary_sign_from_primary,
                           lhs.secondary_shift_from_primary) <
                       std::tie(
                           rhs.primary_run,
                           rhs.secondary_sign_from_primary,
                           rhs.secondary_shift_from_primary);
            });
        constraints.erase(
            std::unique(
                constraints.begin(),
                constraints.end(),
                [](const PrimaryConstraint& lhs,
                   const PrimaryConstraint& rhs) {
                    return lhs.primary_run ==
                               rhs.primary_run &&
                           lhs.secondary_sign_from_primary ==
                               rhs.secondary_sign_from_primary &&
                           lhs.secondary_shift_from_primary ==
                               rhs.secondary_shift_from_primary;
                }),
            constraints.end());
    }

    const auto support_weight =
        [&](uint32_t run_index) {
            const auto& run = runs[run_index];
            const unsigned __int128 occurrence_count =
                run.leaf_spans.size();
            const unsigned __int128 pair_count =
                occurrence_count > 1
                ? occurrence_count *
                      (occurrence_count - 1) /
                      2
                : 1;
            return pair_count *
                   static_cast<unsigned __int128>(
                       run.col_end - run.col_beg);
        };
    const auto run_content_less =
        [&](uint32_t lhs_index,
            uint32_t rhs_index) {
            const auto& lhs = runs[lhs_index];
            const auto& rhs = runs[rhs_index];
            const auto span_less =
                [](const LeafRunSpan& left,
                   const LeafRunSpan& right) {
                    return std::tie(
                               left.leaf_name,
                               left.chr_name,
                               left.start,
                               left.length,
                               left.reversed,
                               left.row_id) <
                           std::tie(
                               right.leaf_name,
                               right.chr_name,
                               right.start,
                               right.length,
                               right.reversed,
                               right.row_id);
                };
            if (std::lexicographical_compare(
                    lhs.leaf_spans.begin(),
                    lhs.leaf_spans.end(),
                    rhs.leaf_spans.begin(),
                    rhs.leaf_spans.end(),
                    span_less)) {
                return true;
            }
            if (std::lexicographical_compare(
                    rhs.leaf_spans.begin(),
                    rhs.leaf_spans.end(),
                    lhs.leaf_spans.begin(),
                    lhs.leaf_spans.end(),
                    span_less)) {
                return false;
            }
            return std::tie(
                       lhs.col_beg,
                       lhs.col_end) <
                   std::tie(
                       rhs.col_beg,
                       rhs.col_end);
        };
    struct SecondaryRunGroup {
        uint64_t block_id = 0;
        unsigned __int128 support = 0;
        std::vector<uint32_t> runs;
    };
    std::map<uint64_t, SecondaryRunGroup>
        group_by_block;
    for (uint32_t secondary_run :
         secondary_runs) {
        const uint64_t block_id =
            runs[secondary_run].block_id;
        auto& group =
            group_by_block[block_id];
        group.block_id = block_id;
        const auto weight =
            support_weight(secondary_run);
        if (std::numeric_limits<
                unsigned __int128>::max() -
                group.support <
            weight) {
            throw std::overflow_error(
                "HAL secondary block support weight overflow");
        }
        group.support += weight;
        group.runs.push_back(secondary_run);
    }
    std::vector<SecondaryRunGroup> groups;
    groups.reserve(group_by_block.size());
    for (auto& [block_id, group] :
         group_by_block) {
        (void)block_id;
        std::sort(
            group.runs.begin(),
            group.runs.end(),
            run_content_less);
        groups.push_back(std::move(group));
    }
    std::sort(
        groups.begin(),
        groups.end(),
        [&](const SecondaryRunGroup& lhs,
            const SecondaryRunGroup& rhs) {
            if (lhs.support != rhs.support) {
                return lhs.support > rhs.support;
            }
            if (std::lexicographical_compare(
                    lhs.runs.begin(),
                    lhs.runs.end(),
                    rhs.runs.begin(),
                    rhs.runs.end(),
                    run_content_less)) {
                return true;
            }
            if (std::lexicographical_compare(
                    rhs.runs.begin(),
                    rhs.runs.end(),
                    lhs.runs.begin(),
                    lhs.runs.end(),
                    run_content_less)) {
                return false;
            }
            return lhs.block_id <
                   rhs.block_id;
        });

    RollbackRunCoordinateUnion coordinates(
        runs.size());
    std::vector<uint8_t> rejected(
        runs.size(),
        0);
    for (const auto& group : groups) {
        const size_t snapshot =
            coordinates.snapshot();
        bool consistent = true;
        bool contributes_homology = false;
        for (uint32_t secondary_run :
             group.runs) {
            contributes_homology =
                contributes_homology ||
                secondary_has_novel_span[
                    secondary_run] != 0;
            const auto constraints_it =
                constraints_by_secondary.find(
                    secondary_run);
            if (constraints_it ==
                constraints_by_secondary.end()) {
                continue;
            }
            const auto& constraints =
                constraints_it->second;
            const auto& anchor =
                constraints.front();
            for (size_t constraint_index = 1;
                 constraint_index <
                     constraints.size();
                 ++constraint_index) {
                const auto& constraint =
                    constraints[constraint_index];
                const int8_t sign =
                    static_cast<int8_t>(
                        constraint
                            .secondary_sign_from_primary *
                        anchor
                            .secondary_sign_from_primary);
                const __int128 shift_wide =
                    static_cast<__int128>(
                        constraint
                            .secondary_sign_from_primary) *
                    (static_cast<__int128>(
                         anchor
                             .secondary_shift_from_primary) -
                     static_cast<__int128>(
                         constraint
                             .secondary_shift_from_primary));
                if (shift_wide <
                        std::numeric_limits<int64_t>::min() ||
                    shift_wide >
                        std::numeric_limits<int64_t>::max()) {
                    throw std::overflow_error(
                        "HAL secondary homology relation "
                        "exceeds signed coordinate range");
                }
                const bool separate_components =
                    coordinates.find(
                        anchor.primary_run).root !=
                    coordinates.find(
                        constraint.primary_run).root;
                if (!coordinates.unite(
                        anchor.primary_run,
                        constraint.primary_run,
                        sign,
                        static_cast<int64_t>(
                            shift_wide))) {
                    consistent = false;
                    break;
                }
                contributes_homology =
                    contributes_homology ||
                    separate_components;
            }
            if (!consistent) {
                break;
            }
        }
        if (consistent &&
            contributes_homology) {
            stats.accepted_runs +=
                group.runs.size();
            continue;
        }

        coordinates.rollback(snapshot);
        stats.rejected_block_ids.insert(
            group.block_id);
        stats.rejected_runs +=
            group.runs.size();
        if (consistent) {
            stats.redundant_runs +=
                group.runs.size();
        } else {
            stats.conflict_rejected_runs +=
                group.runs.size();
        }
        for (uint32_t secondary_run :
             group.runs) {
            rejected[secondary_run] = 1;
            const uint64_t run_length =
                runs[secondary_run].col_end -
                runs[secondary_run].col_beg;
            if (std::numeric_limits<uint64_t>::max() -
                    stats.rejected_bases <
                run_length) {
                throw std::overflow_error(
                    "HAL rejected secondary base count overflow");
            }
            stats.rejected_bases +=
                run_length;
        }
    }

    size_t output_index = 0;
    for (size_t input_index = 0;
         input_index < runs.size();
         ++input_index) {
        if (rejected[input_index] != 0) {
            continue;
        }
        if (output_index != input_index) {
            runs[output_index] =
                std::move(runs[input_index]);
        }
        ++output_index;
    }
    runs.resize(output_index);
    return stats;
}

TreeMeta buildLeafOnlyTree(
    const std::map<
        SpeciesName,
        SeqPro::SharedManagerVariant>&
        seqpro_managers) {
    TreeMeta leaf_tree;
    leaf_tree.nodes.reserve(
        seqpro_managers.size());
    leaf_tree.leaf_ids.reserve(
        seqpro_managers.size());
    int leaf_index = 0;
    for (const auto& [species, unused_manager] :
         seqpro_managers) {
        (void)unused_manager;
        const int node_id =
            static_cast<int>(
                leaf_tree.nodes.size());
        leaf_tree.nodes.push_back(
            TreeNodeMeta{
                .id = node_id,
                .name = species,
                .parent = -1,
                .children = {},
                .branch_length_to_parent = 0.0,
                .is_leaf = true,
                .leaf_index = leaf_index++});
        leaf_tree.name_to_id.emplace(
            species,
            node_id);
        leaf_tree.leaf_ids.push_back(
            node_id);
    }
    return leaf_tree;
}

std::vector<BlockMSA> buildLeafBlockMSAs(
    const std::vector<std::weak_ptr<Block>>& blocks,
    const TreeMeta& leaf_tree,
    const std::map<
        SpeciesName,
        SeqPro::SharedManagerVariant>&
        seqpro_managers) {
    std::vector<BlockMSA> block_msas;
    block_msas.reserve(blocks.size());
    for (const auto& weak_block : blocks) {
        const auto block =
            weak_block.lock();
        if (!block) {
            continue;
        }
        auto msa = buildBlockMSA(
            block,
            leaf_tree,
            seqpro_managers,
            nullptr);
        if (msa.block &&
            !msa.leaf_rows.empty()) {
            block_msas.push_back(
                std::move(msa));
        }
    }
    std::sort(
        block_msas.begin(),
        block_msas.end(),
        [](const BlockMSA& lhs,
           const BlockMSA& rhs) {
            return exportBlockOrderLess(
                lhs.order_key,
                rhs.order_key);
        });
    return block_msas;
}

std::unordered_set<uint64_t>
findRejectedSecondaryHomologyBlocksImpl(
    const std::vector<std::weak_ptr<Block>>& blocks,
    const std::map<
        SpeciesName,
        SeqPro::SharedManagerVariant>&
        seqpro_managers) {
    const TreeMeta leaf_tree =
        buildLeafOnlyTree(
            seqpro_managers);
    const auto block_msas =
        buildLeafBlockMSAs(
            blocks,
            leaf_tree,
            seqpro_managers);
    auto runs =
        buildColumnRuns(
            block_msas,
            leaf_tree);
    auto selection =
        selectCoordinateConsistentSecondaryRuns(
            runs);
    return std::move(
        selection.rejected_block_ids);
}

ColumnRun sliceColumnRun(
    const ColumnRun& run,
    uint32_t local_begin,
    uint32_t local_end) {
    const uint32_t run_length = run.col_end - run.col_beg;
    if (local_begin >= local_end ||
        local_end > run_length) {
        throw std::runtime_error(
            "HAL column run slice lies outside its source run");
    }

    ColumnRun result;
    result.run_id = run.run_id;
    result.block_id = run.block_id;
    result.source_block_ids =
        run.source_block_ids;
    result.col_beg = run.col_beg + local_begin;
    result.col_end = run.col_beg + local_end;
    result.secondary_homology =
        run.secondary_homology;
    result.leaf_spans.reserve(run.leaf_spans.size());
    const uint32_t slice_length = local_end - local_begin;
    for (const auto& span : run.leaf_spans) {
        if (span.length != run_length ||
            span.dna.size() != run_length) {
            throw std::runtime_error(
                "HAL elementary run span length is inconsistent");
        }
        LeafRunSpan sliced = span;
        sliced.start += span.reversed
            ? span.length - local_end
            : local_begin;
        sliced.length = slice_length;
        sliced.dna =
            span.dna.substr(local_begin, slice_length);
        result.leaf_spans.push_back(std::move(sliced));
    }
    return result;
}

void reverseColumnRun(ColumnRun& run) {
    for (auto& span : run.leaf_spans) {
        reverseComplement(span.dna);
        span.reversed = !span.reversed;
    }
}

void finalizeNormalizedRun(
    ColumnRun& run,
    const TreeMeta& tree,
    RunNormalizationStats& stats) {
    std::sort(
        run.source_block_ids.begin(),
        run.source_block_ids.end());
    run.source_block_ids.erase(
        std::unique(
            run.source_block_ids.begin(),
            run.source_block_ids.end()),
        run.source_block_ids.end());
    if (run.source_block_ids.empty()) {
        throw std::runtime_error(
            "HAL normalized run has no source Block provenance");
    }
    std::sort(
        run.leaf_spans.begin(),
        run.leaf_spans.end(),
        [](const LeafRunSpan& lhs, const LeafRunSpan& rhs) {
            return std::tie(
                       lhs.hal_sequence_name,
                       lhs.start,
                       lhs.length,
                       lhs.row_id) <
                   std::tie(
                       rhs.hal_sequence_name,
                       rhs.start,
                       rhs.length,
                       rhs.row_id);
        });

    std::vector<LeafRunSpan> unique_spans;
    unique_spans.reserve(run.leaf_spans.size());
    for (auto& span : run.leaf_spans) {
        if (span.length != run.col_end - run.col_beg ||
            span.dna.size() != span.length) {
            throw std::runtime_error(
                "HAL normalized leaf span length is inconsistent");
        }
        if (!unique_spans.empty()) {
            const auto& previous = unique_spans.back();
            const bool same_occurrence =
                previous.hal_sequence_name ==
                    span.hal_sequence_name &&
                previous.start == span.start &&
                previous.length == span.length;
            if (same_occurrence) {
                if (previous.reversed != span.reversed ||
                    previous.dna != span.dna) {
                    throw std::runtime_error(
                        "HAL homology links assign conflicting "
                        "orientations or DNA to one leaf occurrence");
                }
                ++stats.duplicate_leaf_occurrences;
                continue;
            }
        }
        unique_spans.push_back(std::move(span));
    }
    run.leaf_spans = std::move(unique_spans);
    run.leaf_present.assign(tree.leaf_ids.size(), 0);
    for (const auto& span : run.leaf_spans) {
        const auto node_it =
            tree.name_to_id.find(span.leaf_name);
        if (node_it == tree.name_to_id.end()) {
            throw std::runtime_error(
                "HAL normalized leaf is absent from the phylogeny");
        }
        const auto& node = tree.nodes[node_it->second];
        if (!node.is_leaf || node.leaf_index < 0) {
            throw std::runtime_error(
                "HAL normalized occurrence belongs to a non-leaf node");
        }
        run.leaf_present[
            static_cast<size_t>(node.leaf_index)] = 1;
    }
}

RunNormalizationStats normalizeOverlappingColumnRuns(
    std::vector<ColumnRun>& runs,
    const TreeMeta& tree) {
    RunNormalizationStats stats;
    stats.input_runs = runs.size();
    if (runs.empty()) {
        return stats;
    }

    struct LeafWindow {
        uint64_t start = 0;
        uint64_t end = 0;
        uint32_t run_index = 0;
        uint32_t span_index = 0;
    };
    struct OverlapConstraint {
        uint32_t first_run = 0;
        uint32_t second_run = 0;
        int64_t first_local_begin = 0;
        int64_t first_local_end = 0;
    };

    std::vector<OverlapConstraint> overlap_constraints;

    RunCoordinateUnion coordinates(runs.size());
    std::unordered_map<
        std::string,
        std::vector<LeafWindow>>
        windows_by_sequence;
    for (uint32_t run_index = 0;
         run_index < runs.size();
         ++run_index) {
        const auto& run = runs[run_index];
        const int64_t run_length = runColumnLength(run);
        if (run.leaf_spans.size() >
            static_cast<size_t>(
                std::numeric_limits<uint32_t>::max())) {
            throw std::overflow_error(
                "Too many leaf spans in one HAL column run");
        }
        for (uint32_t span_index = 0;
             span_index < run.leaf_spans.size();
             ++span_index) {
            const auto& span = run.leaf_spans[span_index];
            if (span.length != run_length ||
                span.length == 0) {
                throw std::runtime_error(
                    "HAL elementary run span length is inconsistent");
            }
            const uint64_t end = span.start + span.length;
            if (end < span.start) {
                throw std::overflow_error(
                    "HAL leaf span coordinate overflow");
            }
            windows_by_sequence[
                span.hal_sequence_name]
                .push_back(
                    LeafWindow{
                        span.start,
                        end,
                        run_index,
                        span_index});
        }
    }

    for (auto& [sequence_name, windows] :
         windows_by_sequence) {
        std::sort(
            windows.begin(),
            windows.end(),
            [](const LeafWindow& lhs,
               const LeafWindow& rhs) {
                return std::tie(
                           lhs.start,
                           lhs.end,
                           lhs.run_index,
                           lhs.span_index) <
                       std::tie(
                           rhs.start,
                           rhs.end,
                           rhs.run_index,
                           rhs.span_index);
            });
        std::vector<size_t> active;
        for (size_t window_index = 0;
             window_index < windows.size();
             ++window_index) {
            const auto& current = windows[window_index];
            std::erase_if(
                active,
                [&](size_t active_index) {
                    return windows[active_index].end <=
                           current.start;
                });
            for (size_t active_index : active) {
                const auto& previous =
                    windows[active_index];
                if (previous.end <= current.start) {
                    continue;
                }
                if (previous.run_index ==
                    current.run_index) {
                    throw std::runtime_error(
                        "HAL column run contains overlapping "
                        "copies of leaf sequence " +
                        sequence_name);
                }

                const uint64_t shared_coordinate =
                    std::max(
                        previous.start,
                        current.start);
                const uint64_t shared_end =
                    std::min(
                        previous.end,
                        current.end);
                const auto& previous_span =
                    runs[previous.run_index]
                        .leaf_spans[
                            previous.span_index];
                const auto& current_span =
                    runs[current.run_index]
                        .leaf_spans[
                            current.span_index];
                const int64_t previous_offset =
                    leafBoundaryOffset(
                        previous_span,
                        shared_coordinate);
                const int64_t current_offset =
                    leafBoundaryOffset(
                        current_span,
                        shared_coordinate);
                const int8_t sign =
                    previous_span.reversed ==
                            current_span.reversed
                        ? 1
                        : -1;
                const int64_t shift =
                    current_offset -
                    static_cast<int64_t>(sign) *
                        previous_offset;
                coordinates.unite(
                    previous.run_index,
                    current.run_index,
                    sign,
                    shift);
                const int64_t previous_end_offset =
                    leafBoundaryOffset(
                        previous_span,
                        shared_end);
                overlap_constraints.push_back(
                    OverlapConstraint{
                        previous.run_index,
                        current.run_index,
                        std::min(
                            previous_offset,
                            previous_end_offset),
                        std::max(
                            previous_offset,
                            previous_end_offset)});
                ++stats.overlap_constraints;
            }
            active.push_back(window_index);
        }
    }
    windows_by_sequence.clear();

    std::map<uint32_t, std::vector<uint32_t>>
        members_by_root;
    std::vector<RunCoordinateTransform> transforms(
        runs.size());
    for (uint32_t run_index = 0;
         run_index < runs.size();
         ++run_index) {
        transforms[run_index] =
            coordinates.find(run_index);
        if (coordinates.componentSize(
                transforms[run_index].root) > 1) {
            members_by_root[
                transforms[run_index].root]
                .push_back(run_index);
        }
    }
    stats.connected_components =
        members_by_root.size();

    std::map<
        uint32_t,
        std::vector<const OverlapConstraint*>>
        constraints_by_root;
    for (const auto& constraint :
         overlap_constraints) {
        const uint32_t first_root =
            transforms[constraint.first_run].root;
        const uint32_t second_root =
            transforms[constraint.second_run].root;
        if (first_root != second_root) {
            throw std::runtime_error(
                "HAL overlap constraint endpoints were not unified");
        }
        constraints_by_root[first_root].push_back(
            &constraint);
    }

    std::vector<ColumnRun> normalized;
    normalized.reserve(runs.size());
    for (uint32_t run_index = 0;
         run_index < runs.size();
         ++run_index) {
        if (coordinates.componentSize(
                transforms[run_index].root) == 1) {
            normalized.push_back(
                std::move(runs[run_index]));
        }
    }

    for (auto& [root, members] : members_by_root) {
        std::sort(members.begin(), members.end());
        std::unordered_map<uint32_t, uint32_t>
            local_index_by_run;
        local_index_by_run.reserve(members.size());
        for (uint32_t local_index = 0;
             local_index < members.size();
             ++local_index) {
            local_index_by_run.emplace(
                members[local_index],
                local_index);
        }

        const auto constraint_it =
            constraints_by_root.find(root);
        if (constraint_it ==
            constraints_by_root.end()) {
            throw std::runtime_error(
                "HAL connected run component has no overlap constraints");
        }
        const auto& component_constraints =
            constraint_it->second;
        std::vector<int64_t> boundaries;
        boundaries.reserve(
            2 * (members.size() +
                 component_constraints.size()));
        for (uint32_t run_index : members) {
            const auto& transform =
                transforms[run_index];
            const int64_t length =
                runColumnLength(runs[run_index]);
            const int64_t first =
                transform.shift;
            const int64_t second =
                static_cast<int64_t>(
                    transform.sign) *
                    length +
                transform.shift;
            boundaries.push_back(
                std::min(first, second));
            boundaries.push_back(
                std::max(first, second));
        }
        for (const auto* constraint :
             component_constraints) {
            const auto& transform =
                transforms[
                    constraint->first_run];
            const int64_t first =
                static_cast<int64_t>(
                    transform.sign) *
                    constraint->first_local_begin +
                transform.shift;
            const int64_t second =
                static_cast<int64_t>(
                    transform.sign) *
                    constraint->first_local_end +
                transform.shift;
            boundaries.push_back(
                std::min(first, second));
            boundaries.push_back(
                std::max(first, second));
        }
        std::sort(boundaries.begin(), boundaries.end());
        boundaries.erase(
            std::unique(
                boundaries.begin(),
                boundaries.end()),
            boundaries.end());
        if (boundaries.size() < 2) {
            throw std::runtime_error(
                "HAL connected run component has no positive span");
        }

        const size_t cell_count =
            boundaries.size() - 1;
        size_t tree_size = 1;
        while (tree_size < cell_count) {
            if (tree_size >
                std::numeric_limits<size_t>::max() /
                    2) {
                throw std::overflow_error(
                    "HAL overlap event tree size overflow");
            }
            tree_size *= 2;
        }
        struct LocalEdge {
            uint32_t first = 0;
            uint32_t second = 0;
        };
        struct IntervalEvents {
            std::vector<uint32_t> active_runs;
            std::vector<LocalEdge> edges;
        };
        std::vector<IntervalEvents> events(
            tree_size * 2);

        auto boundaryRange = [&](
            int64_t low,
            int64_t high) {
            const size_t begin =
                static_cast<size_t>(
                    std::lower_bound(
                        boundaries.begin(),
                        boundaries.end(),
                        low) -
                    boundaries.begin());
            const size_t end =
                static_cast<size_t>(
                    std::lower_bound(
                        boundaries.begin(),
                        boundaries.end(),
                        high) -
                    boundaries.begin());
            if (begin >= end ||
                begin >= cell_count ||
                end > cell_count ||
                boundaries[begin] != low ||
                boundaries[end] != high) {
                throw std::runtime_error(
                    "HAL overlap event boundary is absent "
                    "from its coordinate component");
            }
            return std::pair{begin, end};
        };
        auto addRunInterval = [&](
            size_t begin,
            size_t end,
            uint32_t local_run) {
            begin += tree_size;
            end += tree_size;
            while (begin < end) {
                if ((begin & 1U) != 0) {
                    events[begin++].active_runs
                        .push_back(local_run);
                }
                if ((end & 1U) != 0) {
                    events[--end].active_runs
                        .push_back(local_run);
                }
                begin /= 2;
                end /= 2;
            }
        };
        auto addEdgeInterval = [&](
            size_t begin,
            size_t end,
            LocalEdge edge) {
            begin += tree_size;
            end += tree_size;
            while (begin < end) {
                if ((begin & 1U) != 0) {
                    events[begin++].edges
                        .push_back(edge);
                }
                if ((end & 1U) != 0) {
                    events[--end].edges
                        .push_back(edge);
                }
                begin /= 2;
                end /= 2;
            }
        };

        for (uint32_t local_index = 0;
             local_index < members.size();
             ++local_index) {
            const uint32_t run_index =
                members[local_index];
            const auto& transform =
                transforms[run_index];
            const int64_t length =
                runColumnLength(runs[run_index]);
            const int64_t first =
                transform.shift;
            const int64_t second =
                static_cast<int64_t>(
                    transform.sign) *
                    length +
                transform.shift;
            const auto [begin, end] =
                boundaryRange(
                    std::min(first, second),
                    std::max(first, second));
            addRunInterval(
                begin,
                end,
                local_index);
        }
        for (const auto* constraint :
             component_constraints) {
            const auto& transform =
                transforms[
                    constraint->first_run];
            const int64_t first =
                static_cast<int64_t>(
                    transform.sign) *
                    constraint->first_local_begin +
                transform.shift;
            const int64_t second =
                static_cast<int64_t>(
                    transform.sign) *
                    constraint->first_local_end +
                transform.shift;
            const auto [begin, end] =
                boundaryRange(
                    std::min(first, second),
                    std::max(first, second));
            addEdgeInterval(
                begin,
                end,
                LocalEdge{
                    local_index_by_run.at(
                        constraint->first_run),
                    local_index_by_run.at(
                        constraint->second_run)});
        }

        auto emitGroup = [&](
            const std::vector<uint32_t>&
                local_members,
            int64_t global_begin,
            int64_t global_end) {
            if (local_members.empty() ||
                global_begin >= global_end) {
                throw std::runtime_error(
                    "HAL normalized homology group is empty");
            }
            ColumnRun combined;
            bool initialized = false;
            for (uint32_t local_index :
                 local_members) {
                const uint32_t run_index =
                    members.at(local_index);
                const auto& source =
                    runs[run_index];
                const auto& transform =
                    transforms[run_index];
                const int64_t length =
                    runColumnLength(source);
                const int64_t local_begin =
                    transform.sign > 0
                    ? global_begin -
                          transform.shift
                    : transform.shift -
                          global_end;
                const int64_t local_end =
                    transform.sign > 0
                    ? global_end -
                          transform.shift
                    : transform.shift -
                          global_begin;
                if (local_begin < 0 ||
                    local_end > length ||
                    local_begin >= local_end ||
                    local_end >
                        std::numeric_limits<
                            uint32_t>::max()) {
                    throw std::runtime_error(
                        "HAL normalized run boundary is invalid");
                }
                ColumnRun fragment =
                    sliceColumnRun(
                        source,
                        static_cast<uint32_t>(
                            local_begin),
                        static_cast<uint32_t>(
                            local_end));
                if (transform.sign < 0) {
                    reverseColumnRun(fragment);
                }
                combined.secondary_homology =
                    combined.secondary_homology ||
                    fragment.secondary_homology;
                combined.source_block_ids.insert(
                    combined.source_block_ids.end(),
                    fragment.source_block_ids.begin(),
                    fragment.source_block_ids.end());
                if (!initialized ||
                    std::tie(
                        fragment.block_id,
                        fragment.col_beg,
                        fragment.run_id) <
                    std::tie(
                        combined.block_id,
                        combined.col_beg,
                        combined.run_id)) {
                    combined.run_id =
                        fragment.run_id;
                    combined.block_id =
                        fragment.block_id;
                    combined.col_beg =
                        fragment.col_beg;
                    combined.col_end =
                        fragment.col_end;
                    initialized = true;
                }
                combined.leaf_spans.insert(
                    combined.leaf_spans.end(),
                    std::make_move_iterator(
                        fragment.leaf_spans.begin()),
                    std::make_move_iterator(
                        fragment.leaf_spans.end()));
            }
            finalizeNormalizedRun(
                combined,
                tree,
                stats);
            normalized.push_back(
                std::move(combined));
        };

        RollbackRunConnectivity connectivity(
            members.size());
        std::vector<uint32_t> active_runs;
        std::map<std::vector<uint32_t>, int64_t>
            open_groups;
        auto consumeCell = [&](size_t cell_index) {
            std::map<
                uint32_t,
                std::vector<uint32_t>>
                members_by_local_root;
            for (uint32_t local_run : active_runs) {
                members_by_local_root[
                    connectivity.find(local_run)]
                    .push_back(local_run);
            }
            std::vector<std::vector<uint32_t>>
                current_groups;
            current_groups.reserve(
                members_by_local_root.size());
            for (auto& [local_root, group] :
                 members_by_local_root) {
                (void)local_root;
                std::sort(group.begin(), group.end());
                current_groups.push_back(
                    std::move(group));
            }
            std::sort(
                current_groups.begin(),
                current_groups.end());

            const int64_t global_begin =
                boundaries[cell_index];
            for (auto open_it = open_groups.begin();
                 open_it != open_groups.end();) {
                if (!std::binary_search(
                        current_groups.begin(),
                        current_groups.end(),
                        open_it->first)) {
                    emitGroup(
                        open_it->first,
                        open_it->second,
                        global_begin);
                    open_it =
                        open_groups.erase(open_it);
                } else {
                    ++open_it;
                }
            }
            for (const auto& group :
                 current_groups) {
                open_groups.try_emplace(
                    group,
                    global_begin);
            }
        };

        std::function<void(size_t, size_t, size_t)>
            visitEventTree;
        visitEventTree = [&](
            size_t node,
            size_t cell_begin,
            size_t cell_end) {
            if (cell_begin >= cell_count) {
                return;
            }
            const size_t connectivity_snapshot =
                connectivity.snapshot();
            const size_t active_snapshot =
                active_runs.size();
            active_runs.insert(
                active_runs.end(),
                events[node].active_runs.begin(),
                events[node].active_runs.end());
            for (const auto& edge :
                 events[node].edges) {
                connectivity.unite(
                    edge.first,
                    edge.second);
            }
            if (cell_end - cell_begin == 1) {
                consumeCell(cell_begin);
            } else {
                const size_t middle =
                    cell_begin +
                    (cell_end - cell_begin) / 2;
                visitEventTree(
                    node * 2,
                    cell_begin,
                    middle);
                visitEventTree(
                    node * 2 + 1,
                    middle,
                    cell_end);
            }
            active_runs.resize(active_snapshot);
            connectivity.rollback(
                connectivity_snapshot);
        };
        visitEventTree(1, 0, tree_size);
        for (const auto& [group, global_begin] :
             open_groups) {
            emitGroup(
                group,
                global_begin,
                boundaries.back());
        }
    }

    std::sort(
        normalized.begin(),
        normalized.end(),
        [](const ColumnRun& lhs,
           const ColumnRun& rhs) {
            return std::tie(
                       lhs.block_id,
                       lhs.col_beg,
                       lhs.col_end,
                       lhs.run_id) <
                   std::tie(
                       rhs.block_id,
                       rhs.col_beg,
                       rhs.col_end,
                       rhs.run_id);
        });
    for (size_t run_index = 0;
         run_index < normalized.size();
         ++run_index) {
        normalized[run_index].run_id =
            run_index + 1;
    }
    runs = std::move(normalized);
    stats.output_runs = runs.size();
    return stats;
}

void sanitizeLeafCoverage(
    const std::vector<ColumnRun>& runs,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {

    struct LeafWindow {
        uint64_t start = 0;
        uint64_t end = 0;
        uint64_t run_id = 0;
    };

    std::unordered_map<std::string, std::vector<LeafWindow>> by_sequence;
    for (const auto& run : runs) {
        for (const auto& span : run.leaf_spans) {
            by_sequence[span.hal_sequence_name].push_back(
                LeafWindow{span.start, span.start + span.length, run.run_id});
        }
    }

    for (auto& [seq_name, windows] : by_sequence) {
        std::sort(windows.begin(), windows.end(), [](const LeafWindow& a, const LeafWindow& b) {
            if (a.start != b.start) {
                return a.start < b.start;
            }
            return a.end < b.end;
        });

        for (size_t i = 1; i < windows.size(); ++i) {
            if (windows[i].start < windows[i - 1].end) {
                std::ostringstream oss;
                oss << "HAL export sanitizer failed: overlapping leaf intervals on " << seq_name
                    << " between run " << windows[i - 1].run_id
                    << " and run " << windows[i].run_id;
                throw std::runtime_error(oss.str());
            }
        }
    }

    for (const auto& [species_name, shared_mgr] : seqpro_managers) {
        for (const auto& chr_name : fetchSequenceNames(shared_mgr)) {
            (void)fetchSequenceLength(shared_mgr, chr_name);
        }
    }
}

struct MafReblockedRow {
    std::string row_id;
    std::string leaf_name;
    std::string chr_name;
    std::string hal_sequence_name;
    uint64_t interval_start = 0;
    uint64_t interval_end = 0;
    uint64_t ungapped_length = 0;
    bool reversed = false;
    std::string aligned_dna;
};

struct MafReblockedBlock {
    bool initialized = false;
    bool has_multi_occurrence_column = false;
    uint32_t column_end = 0;
    size_t alignment_width = 0;
    std::vector<MafReblockedRow> rows;
    std::unordered_map<
        std::string,
        std::vector<size_t>>
        row_indices;
};

bool canAppendMafRun(
    const MafReblockedBlock& block,
    const ColumnRun& run) {
    return !block.initialized ||
        run.col_beg ==
            block.column_end;
}

bool mafRowContinues(
    const MafReblockedRow& row,
    const LeafRunSpan& span,
    uint64_t span_end) {
    if (row.leaf_name !=
            span.leaf_name ||
        row.chr_name !=
            span.chr_name ||
        row.hal_sequence_name !=
            span.hal_sequence_name ||
        row.reversed !=
            span.reversed) {
        return false;
    }
    return span.reversed
        ? span_end ==
              row.interval_start
        : span.start ==
              row.interval_end;
}

void appendMafRun(
    MafReblockedBlock& block,
    const ColumnRun& run,
    size_t reserve_width) {
    const uint32_t run_length =
        run.col_end - run.col_beg;
    if (run_length == 0) {
        throw std::runtime_error(
            "Cannot reblock an empty MAF run");
    }
    if (block.initialized &&
        !canAppendMafRun(block, run)) {
        throw std::logic_error(
            "Cannot append a discontinuous MAF run");
    }

    const size_t old_width =
        block.alignment_width;
    for (auto& row : block.rows) {
        row.aligned_dna.append(
            run_length,
            '-');
    }
    std::unordered_set<size_t>
        used_rows;
    used_rows.reserve(
        run.leaf_spans.size());
    for (const auto& span :
         run.leaf_spans) {
        if (span.length !=
                run_length ||
            span.dna.size() !=
                run_length) {
            throw std::runtime_error(
                "MAF run span length is inconsistent");
        }
        const uint64_t span_end =
            span.start + span.length;
        if (span_end < span.start) {
            throw std::overflow_error(
                "MAF run coordinate overflow");
        }

        size_t row_index =
            std::numeric_limits<
                size_t>::max();
        const auto row_it =
            block.row_indices.find(
                span.row_id);
        if (row_it !=
            block.row_indices.end()) {
            for (const size_t candidate :
                 row_it->second) {
                if (used_rows.count(
                        candidate) == 0 &&
                    mafRowContinues(
                        block.rows[candidate],
                        span,
                        span_end)) {
                    row_index = candidate;
                    break;
                }
            }
        }
        if (row_index ==
            std::numeric_limits<
                size_t>::max()) {
            MafReblockedRow row;
            row.row_id = span.row_id;
            row.leaf_name =
                span.leaf_name;
            row.chr_name =
                span.chr_name;
            row.hal_sequence_name =
                span.hal_sequence_name;
            row.interval_start =
                span.start;
            row.interval_end =
                span_end;
            row.ungapped_length =
                span.length;
            row.reversed =
                span.reversed;
            row.aligned_dna.reserve(
                reserve_width);
            row.aligned_dna.assign(
                old_width,
                '-');
            row.aligned_dna.append(
                span.dna);
            row_index =
                block.rows.size();
            block.rows.push_back(
                std::move(row));
            block.row_indices[
                span.row_id].push_back(
                    row_index);
            used_rows.insert(
                row_index);
            continue;
        }
        used_rows.insert(
            row_index);

        auto& row =
            block.rows[row_index];
        std::copy(
            span.dna.begin(),
            span.dna.end(),
            row.aligned_dna.end() -
                run_length);
        if (span.reversed) {
            row.interval_start =
                span.start;
        } else {
            row.interval_end =
                span_end;
        }
        row.ungapped_length +=
            span.length;
    }

    block.initialized = true;
    block.has_multi_occurrence_column =
        block.has_multi_occurrence_column ||
        run.leaf_spans.size() >= 2;
    block.column_end = run.col_end;
    block.alignment_width +=
        run_length;
}

bool emitReblockedMafBlock(
    std::ostream& output,
    const MafReblockedBlock& block,
    const std::map<
        SpeciesName,
        SeqPro::SharedManagerVariant>&
        seqpro_managers,
    bool pairwise_mode) {
    if (!block.initialized ||
        !block.has_multi_occurrence_column ||
        block.rows.size() < 2) {
        return false;
    }

    std::map<
        SpeciesName,
        std::vector<
            const MafReblockedRow*>>
        rows_by_leaf;
    size_t maximum_copy_count = 0;
    for (const auto& row :
         block.rows) {
        auto& copies =
            rows_by_leaf[row.leaf_name];
        copies.push_back(&row);
        maximum_copy_count =
            std::max(
                maximum_copy_count,
                copies.size());
    }
    for (auto& [leaf_name, copies] :
         rows_by_leaf) {
        (void)leaf_name;
        std::sort(
            copies.begin(),
            copies.end(),
            [](const MafReblockedRow* lhs,
               const MafReblockedRow* rhs) {
                return std::tie(
                           lhs->hal_sequence_name,
                           lhs->interval_start,
                           lhs->interval_end,
                           lhs->row_id) <
                       std::tie(
                           rhs->hal_sequence_name,
                           rhs->interval_start,
                           rhs->interval_end,
                           rhs->row_id);
            });
    }

    output << "a score=0\n";
    for (size_t copy_index = 0;
         copy_index <
             maximum_copy_count;
         ++copy_index) {
        for (const auto& [
                 leaf_name,
                 copies] :
             rows_by_leaf) {
            (void)leaf_name;
            if (copy_index >=
                copies.size()) {
                continue;
            }
            const auto& row =
                *copies[copy_index];
            const auto manager_it =
                seqpro_managers.find(
                    row.leaf_name);
            if (manager_it ==
                seqpro_managers.end()) {
                throw std::runtime_error(
                    "MAF reblocked row references "
                    "an unknown leaf genome: " +
                    row.leaf_name);
            }
            const uint64_t sequence_length =
                fetchSequenceLength(
                    manager_it->second,
                    row.chr_name);
            if (row.interval_end <
                    row.interval_start ||
                row.interval_end >
                    sequence_length ||
                row.ungapped_length !=
                    row.interval_end -
                        row.interval_start ||
                row.aligned_dna.size() !=
                    block.alignment_width) {
                throw std::runtime_error(
                    "MAF reblocked row has invalid "
                    "coordinates or alignment width");
            }
            const uint64_t maf_start =
                row.reversed
                ? sequence_length -
                      row.interval_end
                : row.interval_start;
            const std::string& source =
                pairwise_mode
                ? row.chr_name
                : row.hal_sequence_name;
            output << "s " <<
                std::left <<
                std::setw(20) <<
                source <<
                std::right <<
                std::setw(12) <<
                maf_start <<
                std::setw(12) <<
                row.ungapped_length <<
                ' ' <<
                (row.reversed
                     ? '-'
                     : '+') <<
                std::setw(12) <<
                sequence_length <<
                ' ' <<
                row.aligned_dna <<
                '\n';
        }
    }
    output << '\n';
    return true;
}

size_t emitReblockedMafRuns(
    std::ostream& output,
    std::vector<ColumnRun>::const_iterator begin,
    std::vector<ColumnRun>::const_iterator end,
    const std::map<
        SpeciesName,
        SeqPro::SharedManagerVariant>&
        seqpro_managers,
    bool pairwise_mode) {
    size_t reserve_width = 0;
    for (auto run = begin;
         run != end;
         ++run) {
        reserve_width +=
            run->col_end -
            run->col_beg;
    }

    size_t emitted_blocks = 0;
    MafReblockedBlock block;
    for (auto run = begin;
         run != end;
         ++run) {
        if (block.initialized &&
            !canAppendMafRun(
                block,
                *run)) {
            emitted_blocks +=
                emitReblockedMafBlock(
                    output,
                    block,
                    seqpro_managers,
                    pairwise_mode)
                ? 1
                : 0;
            block = {};
        }
        appendMafRun(
            block,
            *run,
            reserve_width);
    }
    emitted_blocks +=
        emitReblockedMafBlock(
            output,
            block,
            seqpro_managers,
            pairwise_mode)
        ? 1
        : 0;
    return emitted_blocks;
}

void exportCanonicalMafImpl(
    const std::vector<std::weak_ptr<Block>>& blocks,
    const std::filesystem::path& maf_path,
    const std::map<
        SpeciesName,
        SeqPro::SharedManagerVariant>&
        seqpro_managers,
    bool pairwise_mode) {
    const TreeMeta leaf_tree =
        buildLeafOnlyTree(
            seqpro_managers);
    const auto block_msas =
        buildLeafBlockMSAs(
            blocks,
            leaf_tree,
            seqpro_managers);
    auto runs =
        buildColumnRuns(
            block_msas,
            leaf_tree);
    const size_t projected_run_count =
        runs.size();
    const auto selection =
        selectCoordinateConsistentSecondaryRuns(
            runs);
    const auto normalization =
        normalizeOverlappingColumnRuns(
            runs,
            leaf_tree);
    sanitizeLeafCoverage(
        runs,
        seqpro_managers);

    if (!maf_path.parent_path().empty()) {
        std::filesystem::create_directories(
            maf_path.parent_path());
    }
    std::ofstream output(
        maf_path,
        std::ios::binary |
            std::ios::trunc);
    if (!output) {
        throw std::runtime_error(
            "Cannot open: " +
            maf_path.string());
    }
    output << "##maf version=1 scoring=none\n";
    size_t emitted_blocks = 0;
    size_t reblocked_groups = 0;
    size_t isolated_overlap_runs = 0;
    auto run = runs.cbegin();
    while (run != runs.cend()) {
        auto group_end = run;
        ++group_end;
        if (run->source_block_ids.size() ==
            1) {
            while (
                group_end != runs.cend() &&
                group_end->source_block_ids ==
                    run->source_block_ids) {
                ++group_end;
            }
            ++reblocked_groups;
        } else {
            ++isolated_overlap_runs;
        }
        emitted_blocks +=
            emitReblockedMafRuns(
                output,
                run,
                group_end,
                seqpro_managers,
                pairwise_mode);
        run = group_end;
    }
    spdlog::info(
        "MAF export split {} elementary column runs; "
        "selected {}/{} secondary runs "
        "({} conflicting and {} redundant runs rejected, "
        "{} bp), normalized to {} runs across {} "
        "homology components, reblocked {} single-source "
        "groups, kept {} overlap-normalized runs isolated, "
        "and emitted {} MAF blocks",
        projected_run_count,
        selection.accepted_runs,
        selection.candidate_runs,
        selection.conflict_rejected_runs,
        selection.redundant_runs,
        selection.rejected_bases,
        normalization.output_runs,
        normalization.connected_components,
        reblocked_groups,
        isolated_overlap_runs,
        emitted_blocks);
}

std::string buildConsensusDNAImpl(
    const std::vector<std::pair<std::string, double>>& donors,
    size_t expected_length,
    double consensus_threshold) {

    if (donors.empty()) {
        return {};
    }
    std::string consensus(expected_length, 'N');
    for (size_t col = 0; col < expected_length; ++col) {
        std::array<double, 5> weights{0.0, 0.0, 0.0, 0.0, 0.0};
        std::array<double, 5> masked_weights{0.0, 0.0, 0.0, 0.0, 0.0};
        double total = 0.0;
        for (const auto& [dna, weight] : donors) {
            if (col >= dna.size()) {
                continue;
            }
            const unsigned char raw_base = static_cast<unsigned char>(dna[col]);
            char base = static_cast<char>(std::toupper(raw_base));
            size_t idx = 4;
            switch (base) {
            case 'A': idx = 0; break;
            case 'C': idx = 1; break;
            case 'G': idx = 2; break;
            case 'T': idx = 3; break;
            default: idx = 4; break;
            }
            weights[idx] += weight;
            if (std::islower(raw_base)) {
                masked_weights[idx] += weight;
            }
            total += weight;
        }

        size_t best_idx = 4;
        double best_weight = -1.0;
        for (size_t i = 0; i < weights.size(); ++i) {
            if (weights[i] > best_weight) {
                best_weight = weights[i];
                best_idx = i;
            }
        }

        if (total <= 0.0 || best_weight / total < consensus_threshold) {
            consensus[col] = 'N';
        } else {
            static const char kBaseMap[5] = {'A', 'C', 'G', 'T', 'N'};
            const char base = kBaseMap[best_idx];
            consensus[col] =
                masked_weights[best_idx] * 2.0 >= best_weight
                    ? static_cast<char>(std::tolower(
                          static_cast<unsigned char>(base)))
                    : base;
        }
    }
    return consensus;
}


std::vector<EdgeSupport> collectAdjacencySupport(
    int node_id,
    const TreeMeta& tree,
    const std::unordered_map<uint64_t, std::vector<OccurrenceId>>&
        occurrences_by_run,
    const ChildRunCopyMap& parent_copy_by_child_run,
    const std::unordered_map<OccurrenceId, uint32_t>& candidate_index,
    const std::unordered_map<int, NodeModel>& prior_models,
    const std::unordered_map<std::string, std::vector<LeafOccurrence>>&
        leaf_paths,
    const ExportConfig& config,
    ExportStats* stats) {

    auto canonical_key = [](const AdjacencyVote& vote) {
        std::pair<uint64_t, bool> first_end{
            vote.left_run_id,
            vote.left_forward_to_canonical};
        std::pair<uint64_t, bool> second_end{
            vote.right_run_id,
            !vote.right_forward_to_canonical};
        if (second_end < first_end) {
            std::swap(first_end, second_end);
        }
        return RunPairKey{
            first_end.first,
            second_end.first,
            first_end.second,
            !second_end.second};
    };

    struct ProjectedOccurrence {
        OccurrenceId occurrence_id = 0;
        bool forward = true;
        uint64_t start = 0;
        uint32_t length = 0;
    };

    std::vector<ChildEdgeContribution> contributions;
    for (int child_id :
         tree.nodes[static_cast<size_t>(node_id)].children) {
        const auto child_mapping_it =
            parent_copy_by_child_run.find(child_id);
        if (child_mapping_it ==
            parent_copy_by_child_run.end()) {
            continue;
        }
        const long double branch_length =
            tree.nodes[static_cast<size_t>(child_id)]
                .branch_length_to_parent;
        const long double child_weight =
            std::exp(
                -static_cast<long double>(
                    config.phylogenetic_phi) *
                branch_length);

        std::unordered_map<
            RunPairKey,
            ChildEdgeContribution,
            RunPairKeyHash>
            child_edges;
        auto parent_occurrence =
            [&](uint64_t run_id,
                uint32_t copy_index)
            -> std::optional<OccurrenceId> {
            const auto copy_map_it =
                child_mapping_it->second.find(
                    run_id);
            const auto occurrence_it =
                occurrences_by_run.find(run_id);
            if (copy_map_it ==
                    child_mapping_it->second.end() ||
                occurrence_it ==
                    occurrences_by_run.end() ||
                copy_index >=
                    copy_map_it->second.size()) {
                return std::nullopt;
            }
            const uint32_t parent_copy_index =
                copy_map_it->second[copy_index];
            if (parent_copy_index >=
                occurrence_it->second.size()) {
                throw std::runtime_error(
                    "HAL child adjacency maps outside parent copy range");
            }
            const OccurrenceId occurrence_id =
                occurrence_it->second[
                    parent_copy_index];
            return candidate_index.contains(
                       occurrence_id)
                       ? std::optional{
                             occurrence_id}
                       : std::nullopt;
        };

        std::optional<ProjectedOccurrence> previous;
        auto reset_path = [&]() {
            previous.reset();
        };
        auto add_occurrence =
            [&](uint64_t run_id,
                uint32_t copy_index,
                bool forward,
                uint64_t start,
                uint32_t length) {
            const auto mapped =
                parent_occurrence(
                    run_id,
                    copy_index);
            if (!mapped) {
                return;
            }
            const ProjectedOccurrence current{
                *mapped,
                forward,
                start,
                length};
            if (previous) {
                if (previous->occurrence_id ==
                    current.occurrence_id) {
                    if (stats != nullptr) {
                        ++stats
                              ->paralogy_self_adjacency_count;
                    }
                } else {
                    const uint64_t previous_end =
                        previous->start +
                        previous->length;
                    const uint64_t gap_bases =
                        current.start > previous_end
                            ? current.start -
                                  previous_end
                            : 0;
                    const AdjacencyVote vote{
                        previous->occurrence_id,
                        previous->forward,
                        current.occurrence_id,
                        current.forward,
                        previous->length,
                        current.length,
                        gap_bases};
                    const RunPairKey key =
                        canonical_key(vote);
                    auto [edge_it, inserted] =
                        child_edges.try_emplace(
                            key);
                    auto& contribution =
                        edge_it->second;
                    if (inserted) {
                        contribution.from =
                            key.from;
                        contribution.to =
                            key.to;
                        contribution
                            .from_forward_to_canonical =
                            key.from_forward_to_canonical;
                        contribution
                            .to_forward_to_canonical =
                            key.to_forward_to_canonical;
                        contribution.child_id =
                            child_id;
                    }
                    const bool
                        canonical_preserves_order =
                            key.from ==
                            previous
                                ->occurrence_id;
                    contribution.left_length =
                        canonical_preserves_order
                            ? previous->length
                            : current.length;
                    contribution.right_length =
                        canonical_preserves_order
                            ? current.length
                            : previous->length;
                    contribution.minimum_gap =
                        std::min(
                            contribution
                                .minimum_gap,
                            gap_bases);
                    contribution
                        .weighted_support +=
                        calculateCactusZScore(
                            previous->length,
                            current.length,
                            gap_bases,
                            config
                                .adjacency_theta) *
                        child_weight;
                    ++contribution
                          .occurrence_support;
                }
            }
            previous = current;
        };

        const auto& child =
            tree.nodes[
                static_cast<size_t>(child_id)];
        if (child.is_leaf) {
            const auto path_it =
                leaf_paths.find(child.name);
            if (path_it ==
                leaf_paths.end()) {
                continue;
            }
            std::string_view sequence_name;
            for (const auto& occurrence :
                 path_it->second) {
                if (!sequence_name.empty() &&
                    sequence_name !=
                        occurrence.span
                            .hal_sequence_name) {
                    reset_path();
                }
                sequence_name =
                    occurrence.span
                        .hal_sequence_name;
                add_occurrence(
                    occurrence.run_id,
                    occurrence.copy_index,
                    occurrence
                        .forward_to_canonical,
                    occurrence.span.start,
                    occurrence.span.length);
            }
        } else {
            const auto model_it =
                prior_models.find(child_id);
            if (model_it ==
                prior_models.end()) {
                continue;
            }
            for (const auto& sequence :
                 model_it->second.sequences) {
                if (sequence.path.size() !=
                    sequence.joins.size() + 1) {
                    throw std::runtime_error(
                        "HAL internal child sequence has an invalid path/join layout");
                }
                reset_path();
                for (size_t occurrence_index = 0;
                     occurrence_index <
                     sequence.path.size();
                     ++occurrence_index) {
                    if (occurrence_index != 0 &&
                        sequence
                                .joins[
                                    occurrence_index -
                                    1]
                                .kind ==
                            ReferenceJoinKind::SCAFFOLD) {
                        reset_path();
                    }
                    const auto& occurrence =
                        sequence.path[
                            occurrence_index];
                    const auto& placement =
                        model_it->second
                            .placements.at(
                                occurrence
                                    .occurrence_id);
                    add_occurrence(
                        occurrence.run_id,
                        occurrence.copy_index,
                        occurrence.forward,
                        placement.start,
                        placement.length);
                }
            }
        }

        contributions.reserve(
            contributions.size() +
            child_edges.size());
        for (auto& [key, edge] :
             child_edges) {
            (void)key;
            contributions.push_back(
                std::move(edge));
        }
    }

    std::sort(
        contributions.begin(),
        contributions.end(),
        [](const auto& lhs, const auto& rhs) {
            return std::tie(
                       lhs.from,
                       lhs.to,
                       lhs.from_forward_to_canonical,
                       lhs.to_forward_to_canonical,
                       lhs.child_id) <
                   std::tie(
                       rhs.from,
                       rhs.to,
                       rhs.from_forward_to_canonical,
                       rhs.to_forward_to_canonical,
                       rhs.child_id);
        });

    std::vector<EdgeSupport> edges;
    edges.reserve(contributions.size());
    size_t index = 0;
    while (index < contributions.size()) {
        EdgeSupport edge;
        edge.from = contributions[index].from;
        edge.to = contributions[index].to;
        edge.from_forward_to_canonical =
            contributions[index]
                .from_forward_to_canonical;
        edge.to_forward_to_canonical =
            contributions[index]
                .to_forward_to_canonical;
        while (index < contributions.size() &&
               contributions[index].from ==
                   edge.from &&
               contributions[index].to ==
                   edge.to &&
               contributions[index]
                       .from_forward_to_canonical ==
                   edge
                       .from_forward_to_canonical &&
               contributions[index]
                       .to_forward_to_canonical ==
                   edge
                       .to_forward_to_canonical) {
            edge.occurrence_support +=
                contributions[index]
                    .occurrence_support;
            edge.supporting_children.push_back(
                contributions[index].child_id);
            edge.weighted_support +=
                contributions[index]
                    .weighted_support;
            edge.minimum_gap =
                std::min(
                    edge.minimum_gap,
                    contributions[index]
                        .minimum_gap);
            ++index;
        }
        edges.push_back(std::move(edge));
    }
    return edges;
}

std::unordered_map<std::string, std::vector<LeafOccurrence>> buildLeafPaths(
    const std::vector<ColumnRun>& runs) {

    std::unordered_map<std::string, std::vector<LeafOccurrence>> paths;
    for (const auto& run : runs) {
        for (const auto& span : run.leaf_spans) {
            paths[span.leaf_name].push_back(LeafOccurrence{run.run_id, span, !span.reversed});
        }
    }
    for (auto& [leaf_name, ordered] : paths) {
        std::sort(ordered.begin(), ordered.end(), [](const auto& a, const auto& b) {
            if (a.span.hal_sequence_name != b.span.hal_sequence_name) {
                return a.span.hal_sequence_name < b.span.hal_sequence_name;
            }
            if (a.span.start != b.span.start) {
                return a.span.start < b.span.start;
            }
            return a.run_id < b.run_id;
        });
        std::unordered_map<uint64_t, uint32_t> next_copy_index;
        for (auto& occurrence : ordered) {
            occurrence.copy_index =
                next_copy_index[occurrence.run_id]++;
        }
    }
    return paths;
}



std::vector<TerminalEndSupport> buildTerminalEndSupport(
    int node_id,
    const TreeMeta& tree,
    const std::unordered_map<uint64_t, std::vector<OccurrenceId>>&
        occurrences_by_run,
    const ChildRunCopyMap& parent_copy_by_child_run,
    const std::unordered_map<OccurrenceId, uint32_t>& candidate_index,
    const std::unordered_map<int, NodeModel>& prior_models,
    const std::unordered_map<std::string, std::vector<LeafOccurrence>>&
        leaf_paths,
    const ExportConfig& config) {
    using EndKey = std::pair<OccurrenceId, uint8_t>;
    struct LineageContribution {
        uint32_t occurrence_support = 0;
        long double weighted_support = 0.0L;
    };
    std::map<EndKey, std::map<int, LineageContribution>>
        contributions;

    auto parent_occurrence =
        [&](int child_id,
            uint64_t run_id,
            uint32_t copy_index)
        -> std::optional<OccurrenceId> {
        const auto child_it =
            parent_copy_by_child_run.find(child_id);
        if (child_it ==
            parent_copy_by_child_run.end()) {
            return std::nullopt;
        }
        const auto copy_map_it =
            child_it->second.find(run_id);
        if (copy_map_it ==
                child_it->second.end() ||
            copy_index >=
                copy_map_it->second.size()) {
            return std::nullopt;
        }
        const auto run_it =
            occurrences_by_run.find(run_id);
        if (run_it ==
                occurrences_by_run.end() ||
            copy_map_it->second[copy_index] >=
                run_it->second.size()) {
            return std::nullopt;
        }
        const OccurrenceId occurrence_id =
            run_it->second[
                copy_map_it->second[copy_index]];
        return candidate_index.contains(
                   occurrence_id)
                   ? std::optional{occurrence_id}
                   : std::nullopt;
    };
    auto add_terminal =
        [&](int child_id,
            int lineage_id,
            uint64_t run_id,
            uint32_t copy_index,
            OccurrenceEndSide side) {
        const auto mapped =
            parent_occurrence(
                child_id,
                run_id,
                copy_index);
        if (!mapped) {
            return;
        }
        auto& contribution =
            contributions[EndKey{
                *mapped,
                static_cast<uint8_t>(side)}]
                         [lineage_id];
        ++contribution.occurrence_support;
        contribution.weighted_support +=
            std::exp(
                -static_cast<long double>(
                    config.phylogenetic_phi) *
                static_cast<long double>(
                    tree.nodes[
                        static_cast<size_t>(
                            child_id)]
                        .branch_length_to_parent));
    };
    auto add_leaf_terminal =
        [&](int child_id,
            uint64_t run_id,
            uint32_t copy_index,
            bool forward,
            bool sequence_start) {
        const OccurrenceEndSide side =
            sequence_start
                ? (forward
                       ? OccurrenceEndSide::LEFT
                       : OccurrenceEndSide::RIGHT)
                : (forward
                       ? OccurrenceEndSide::RIGHT
                       : OccurrenceEndSide::LEFT);
        add_terminal(
            child_id,
            child_id,
            run_id,
            copy_index,
            side);
    };

    for (int child_id :
         tree.nodes[
             static_cast<size_t>(node_id)]
             .children) {
        const auto& child =
            tree.nodes[
                static_cast<size_t>(child_id)];
        if (child.is_leaf) {
            const auto path_it =
                leaf_paths.find(child.name);
            if (path_it ==
                leaf_paths.end()) {
                continue;
            }
            const LeafOccurrence* first =
                nullptr;
            const LeafOccurrence* last =
                nullptr;
            std::string_view sequence_name;
            auto flush_sequence = [&]() {
                if (first == nullptr ||
                    last == nullptr) {
                    return;
                }
                add_leaf_terminal(
                    child_id,
                    first->run_id,
                    first->copy_index,
                    first
                        ->forward_to_canonical,
                    true);
                add_leaf_terminal(
                    child_id,
                    last->run_id,
                    last->copy_index,
                    last
                        ->forward_to_canonical,
                    false);
            };
            for (const auto& occurrence :
                 path_it->second) {
                if (sequence_name !=
                    occurrence.span
                        .hal_sequence_name) {
                    flush_sequence();
                    sequence_name =
                        occurrence.span
                            .hal_sequence_name;
                    first = nullptr;
                    last = nullptr;
                }
                if (!parent_occurrence(
                        child_id,
                        occurrence.run_id,
                        occurrence.copy_index)) {
                    continue;
                }
                if (first == nullptr) {
                    first = &occurrence;
                }
                last = &occurrence;
            }
            flush_sequence();
            continue;
        }

        const auto model_it =
            prior_models.find(child_id);
        if (model_it ==
            prior_models.end()) {
            continue;
        }
        const NodeModel& child_model =
            model_it->second;
        for (const auto& terminal :
             child_model.terminal_ends) {
            const auto run_it =
                child_model.run_by_occurrence.find(
                    terminal.end.occurrence_id);
            const auto copy_it =
                child_model
                    .copy_index_by_occurrence.find(
                        terminal.end
                            .occurrence_id);
            if (run_it ==
                    child_model
                        .run_by_occurrence.end() ||
                copy_it ==
                    child_model
                        .copy_index_by_occurrence
                        .end()) {
                throw std::runtime_error(
                    "HAL terminal provenance references an unknown child occurrence");
            }
            if (terminal
                    .supporting_lineages.size() >=
                2) {
                add_terminal(
                    child_id,
                    child_id,
                    run_it->second,
                    copy_it->second,
                    terminal.end.side);
            }
        }
    }

    std::vector<TerminalEndSupport>
        terminal_ends;
    terminal_ends.reserve(
        contributions.size());
    for (const auto& [key, by_lineage] :
         contributions) {
        TerminalEndSupport support{};
        support.end = OccurrenceEnd{
            key.first,
            static_cast<
                OccurrenceEndSide>(
                key.second)};
        for (const auto& [lineage_id,
                          contribution] :
             by_lineage) {
            support
                .supporting_lineages
                .push_back(lineage_id);
            support.occurrence_support +=
                contribution
                    .occurrence_support;
            support.weighted_support +=
                contribution
                    .weighted_support;
        }
        terminal_ends.push_back(
            std::move(support));
    }
    return terminal_ends;
}


NodeModel buildNodeModel(
    int node_id,
    const TreeMeta& tree,
    const std::vector<ColumnRun>& runs,
    const std::unordered_map<uint64_t, const ColumnRun*>& run_by_id,
    const std::unordered_map<int, NodeModel>& prior_models,
    const std::unordered_map<std::string, std::vector<LeafOccurrence>>& leaf_paths,
    const ExportConfig& config,
    OccurrenceId* next_occurrence_id,
    ExportStats* stats,
    NodeModelBuildTimings* timings_out) {

    NodeModel model;
    model.genome_name = tree.nodes[node_id].name;
    using Clock = std::chrono::steady_clock;

    if (next_occurrence_id == nullptr || *next_occurrence_id == 0) {
        throw std::invalid_argument(
            "HAL ancestor occurrence allocator is not initialized");
    }
    auto candidate_filter_begin = Clock::now();
    std::vector<uint64_t> candidate_runs;
    candidate_runs.reserve(runs.size());
    for (const auto& run : runs) {
        if (run.present_by_node[node_id]) {
            candidate_runs.push_back(run.run_id);
        }
    }
    if (candidate_runs.empty()) {
        return model;
    }
    std::unordered_map<uint64_t, uint32_t> candidate_index;
    candidate_index.reserve(candidate_runs.size());
    for (uint32_t index = 0; index < candidate_runs.size(); ++index) {
        candidate_index.emplace(candidate_runs[index], index);
    }

    using ContextToken =
        std::tuple<uint8_t, uint64_t, uint8_t>;
    using ContextSignature =
        std::pair<ContextToken, ContextToken>;
    struct ContextOccurrence {
        uint64_t run_id = 0;
        uint32_t copy_index = 0;
        bool forward = true;
    };
    std::map<
        int,
        std::map<
            uint64_t,
            std::map<uint32_t, ContextSignature>>>
        contexts_by_child;
    const ContextToken terminal_context{0, 0, 0};
    auto add_context_path =
        [&](int child_id,
            const std::vector<ContextOccurrence>& path) {
            for (size_t index = 0;
                 index < path.size();
                 ++index) {
                const auto& occurrence = path[index];
                auto context_for_side =
                    [&](OccurrenceEndSide side) {
                        const bool neighbor_before =
                            side == OccurrenceEndSide::LEFT
                                ? occurrence.forward
                                : !occurrence.forward;
                        const bool has_neighbor =
                            neighbor_before
                                ? index != 0
                                : index + 1 < path.size();
                        if (!has_neighbor) {
                            return terminal_context;
                        }
                        const auto& neighbor =
                            path[neighbor_before
                                     ? index - 1
                                     : index + 1];
                        const auto neighbor_side =
                            neighbor_before
                                ? (neighbor.forward
                                       ? OccurrenceEndSide::RIGHT
                                       : OccurrenceEndSide::LEFT)
                                : (neighbor.forward
                                       ? OccurrenceEndSide::LEFT
                                       : OccurrenceEndSide::RIGHT);
                        return ContextToken{
                            1,
                            neighbor.run_id,
                            static_cast<uint8_t>(
                                neighbor_side)};
                    };
                const ContextSignature signature{
                    context_for_side(
                        OccurrenceEndSide::LEFT),
                    context_for_side(
                        OccurrenceEndSide::RIGHT)};
                auto [it, inserted] =
                    contexts_by_child[child_id]
                        [occurrence.run_id]
                            .emplace(
                                occurrence.copy_index,
                                signature);
                if (!inserted && it->second != signature) {
                    throw std::runtime_error(
                        "HAL child occurrence has inconsistent bilateral context");
                }
            }
        };

    for (int child_id :
         tree.nodes[static_cast<size_t>(node_id)].children) {
        if (tree.nodes[static_cast<size_t>(child_id)].is_leaf) {
            const auto& child_name =
                tree.nodes[static_cast<size_t>(child_id)].name;
            const auto path_it = leaf_paths.find(child_name);
            if (path_it == leaf_paths.end()) {
                continue;
            }
            std::vector<ContextOccurrence> path;
            std::string sequence_name;
            for (const auto& occurrence : path_it->second) {
                if (!sequence_name.empty() &&
                    occurrence.span.hal_sequence_name !=
                        sequence_name) {
                    add_context_path(child_id, path);
                    path.clear();
                }
                sequence_name =
                    occurrence.span.hal_sequence_name;
                path.push_back(ContextOccurrence{
                    occurrence.run_id,
                    occurrence.copy_index,
                    occurrence.forward_to_canonical});
            }
            add_context_path(child_id, path);
        } else {
            const auto child_model_it =
                prior_models.find(child_id);
            if (child_model_it == prior_models.end()) {
                continue;
            }
            for (const auto& sequence :
                 child_model_it->second.sequences) {
                std::vector<ContextOccurrence> path;
                path.reserve(sequence.path.size());
                for (const auto& occurrence :
                     sequence.path) {
                    path.push_back(ContextOccurrence{
                        occurrence.run_id,
                        occurrence.copy_index,
                        occurrence.forward});
                }
                add_context_path(child_id, path);
            }
        }
    }

    std::vector<OccurrenceId> candidate_occurrences;
    for (uint64_t run_id : candidate_runs) {
        using CopiesByChild =
            std::map<int, std::vector<uint32_t>>;
        std::map<ContextSignature, CopiesByChild>
            copies_by_signature;
        for (int child_id :
             tree.nodes[static_cast<size_t>(node_id)].children) {
            const auto child_it =
                contexts_by_child.find(child_id);
            if (child_it == contexts_by_child.end()) {
                continue;
            }
            const auto run_it =
                child_it->second.find(run_id);
            if (run_it == child_it->second.end()) {
                continue;
            }
            auto& copy_map =
                model.parent_copy_by_child_run[child_id]
                    [run_id];
            if (!run_it->second.empty()) {
                copy_map.resize(
                    static_cast<size_t>(
                        run_it->second.rbegin()->first) +
                        1,
                    0);
            }
            for (const auto& [copy_index, signature] :
                 run_it->second) {
                copies_by_signature[signature][child_id]
                    .push_back(copy_index);
            }
        }

        uint32_t ancestor_copy_count = 0;
        for (auto& [signature, copies_by_child] :
             copies_by_signature) {
            (void)signature;
            std::vector<size_t> child_counts;
            child_counts.reserve(copies_by_child.size());
            for (auto& [child_id, copies] :
                 copies_by_child) {
                (void)child_id;
                std::sort(copies.begin(), copies.end());
                child_counts.push_back(copies.size());
            }
            if (child_counts.size() < 2) {
                continue;
            }
            std::sort(
                child_counts.begin(),
                child_counts.end(),
                std::greater<>());
            const size_t shared_copy_count =
                child_counts[1];
            if (shared_copy_count >
                std::numeric_limits<uint32_t>::max() -
                    ancestor_copy_count) {
                throw std::overflow_error(
                    "HAL reconciled ancestor copy count overflow");
            }
            for (size_t shared_index = 0;
                 shared_index < shared_copy_count;
                 ++shared_index) {
                const uint32_t parent_copy_index =
                    ancestor_copy_count++;
                for (const auto& [child_id, copies] :
                     copies_by_child) {
                    if (shared_index >= copies.size()) {
                        continue;
                    }
                    model.parent_copy_by_child_run
                        .at(child_id)
                        .at(run_id)
                        .at(copies[shared_index]) =
                        parent_copy_index;
                }
            }
        }
        ancestor_copy_count =
            std::max<uint32_t>(1, ancestor_copy_count);
        auto& run_occurrences =
            model.occurrences_by_run[run_id];
        run_occurrences.reserve(ancestor_copy_count);
        for (uint32_t copy_index = 0;
             copy_index < ancestor_copy_count;
             ++copy_index) {
            if (*next_occurrence_id ==
                std::numeric_limits<OccurrenceId>::max()) {
                throw std::overflow_error(
                    "HAL ancestor occurrence id overflow");
            }
            const OccurrenceId occurrence_id =
                (*next_occurrence_id)++;
            candidate_occurrences.push_back(occurrence_id);
            run_occurrences.push_back(occurrence_id);
            model.run_by_occurrence.emplace(
                occurrence_id,
                run_id);
            model.copy_index_by_occurrence.emplace(
                occurrence_id,
                copy_index);
        }
    }
    std::unordered_map<OccurrenceId, uint32_t> occurrence_index;
    occurrence_index.reserve(candidate_occurrences.size());
    for (uint32_t index = 0;
         index < candidate_occurrences.size();
         ++index) {
        occurrence_index.emplace(candidate_occurrences[index], index);
    }
    uint64_t candidate_filter_ms = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::milliseconds>(
            Clock::now() - candidate_filter_begin)
            .count());
    spdlog::info(
        "HAL export ancestor {} prepared {} runs and {} occurrences "
        "in {} ms",
        model.genome_name,
        candidate_runs.size(),
        candidate_occurrences.size(),
        candidate_filter_ms);

    model.run_dna.reserve(candidate_runs.size());
    model.placements.reserve(candidate_occurrences.size());
    if (stats != nullptr) {
        stats->ancestor_occurrence_count +=
            candidate_occurrences.size();
    }
    std::vector<BucketedDonor> donor_candidates;
    donor_candidates.reserve(tree.nodes[node_id].children.size());

    struct InternalChildLookup {
        std::vector<const std::string*> run_dna;
    };

    struct ChildBuildInput {
        int child_id = -1;
        double weight = 0.0;
        bool is_leaf = false;
        std::string leaf_name;
        const InternalChildLookup* internal_lookup = nullptr;
    };

    std::unordered_map<int, InternalChildLookup> internal_child_lookups;
    internal_child_lookups.reserve(tree.nodes[node_id].children.size());
    std::vector<ChildBuildInput> child_inputs;
    child_inputs.reserve(tree.nodes[node_id].children.size());

    for (int child_id : tree.nodes[node_id].children) {
        double branch = std::max(1e-6, tree.nodes[child_id].branch_length_to_parent);
        double weight = 1.0 / branch;
        if (tree.nodes[child_id].is_leaf) {
            child_inputs.push_back(ChildBuildInput{child_id, weight, true, tree.nodes[child_id].name, nullptr});
            continue;
        }
        auto child_model_it = prior_models.find(child_id);
        if (child_model_it == prior_models.end()) {
            continue;
        }
        auto [lookup_it, inserted] = internal_child_lookups.emplace(child_id, InternalChildLookup{});
        (void)inserted;
        auto& lookup = lookup_it->second;
        lookup.run_dna.assign(candidate_runs.size(), nullptr);
        for (const auto& [run_id, dna] : child_model_it->second.run_dna) {
            auto idx_it = candidate_index.find(run_id);
            if (idx_it != candidate_index.end()) {
                lookup.run_dna[idx_it->second] = &dna;
            }
        }
        child_inputs.push_back(ChildBuildInput{child_id, weight, false, {}, &lookup});
    }

    auto run_dna_begin = Clock::now();
    for (uint32_t candidate_idx = 0; candidate_idx < candidate_runs.size(); ++candidate_idx) {
        uint64_t run_id = candidate_runs[candidate_idx];
        const auto* run = run_by_id.at(run_id);
        donor_candidates.clear();
        for (const auto& child : child_inputs) {
            if (!run->present_by_node[child.child_id]) {
                continue;
            }
            if (child.is_leaf) {
                for (const auto& span : run->leaf_spans) {
                    if (span.leaf_name != child.leaf_name) {
                        continue;
                    }
                    donor_candidates.push_back(
                        BucketedDonor{child.child_id, false, span.dna, child.weight});
                }
                continue;
            }
            const auto* dna = child.internal_lookup->run_dna[candidate_idx];
            if (dna != nullptr) {
                donor_candidates.push_back(
                    BucketedDonor{child.child_id, true, *dna, child.weight});
            }
        }
        auto donors = selectBestDonorsByBucket(donor_candidates);
        if (donors.empty()) {
            for (const auto& span : run->leaf_spans) {
                donors.emplace_back(span.dna, 1.0);
            }
        }
        model.run_dna.emplace(
            run_id,
            buildConsensusDNA(donors, run->col_end - run->col_beg, config.consensus_threshold));
    }
    uint64_t run_dna_ms = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - run_dna_begin).count());
    spdlog::info(
        "HAL export ancestor {} reconstructed {} run sequences "
        "in {} ms",
        model.genome_name,
        model.run_dna.size(),
        run_dna_ms);

    auto edge_collect_begin = Clock::now();
    auto edges = collectAdjacencySupport(
        node_id,
        tree,
        model.occurrences_by_run,
        model.parent_copy_by_child_run,
        occurrence_index,
        prior_models,
        leaf_paths,
        config,
        stats);
    uint64_t edge_collect_ms = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::milliseconds>(
            Clock::now() - edge_collect_begin)
            .count());
    spdlog::info(
        "HAL export ancestor {} collected {} adjacency candidates "
        "in {} ms",
        model.genome_name,
        edges.size(),
        edge_collect_ms);
    std::unordered_map<OccurrenceId, RunOrderKey>
        occurrence_order_keys;
    occurrence_order_keys.reserve(candidate_occurrences.size());
    for (OccurrenceId occurrence_id : candidate_occurrences) {
        const uint64_t run_id =
            model.run_by_occurrence.at(occurrence_id);
        const auto* run = run_by_id.at(run_id);
        occurrence_order_keys.emplace(
            occurrence_id,
            RunOrderKey{run->block_id, run->col_beg});
    }
    auto path_decompose_begin = Clock::now();
    model.terminal_ends =
        buildTerminalEndSupport(
            node_id,
            tree,
            model.occurrences_by_run,
            model.parent_copy_by_child_run,
            occurrence_index,
            prior_models,
            leaf_paths,
            config);
    spdlog::info(
        "HAL export ancestor {} collected {} terminal-end candidates",
        model.genome_name,
        model.terminal_ends.size());
    AncestralSequenceAssembly assembly =
        buildAncestralSequenceAssembly(
            candidate_occurrences,
            edges,
            occurrence_order_keys,
            model.terminal_ends,
            config.scaffold_gap_length,
            stats);
    spdlog::info(
        "HAL export ancestor {} decomposed {} occurrences into {} "
        "reference intervals",
        model.genome_name,
        candidate_occurrences.size(),
        assembly.sequences.size());
    size_t supported_fragment_count = 0;
    for (const auto& sequence : assembly.sequences) {
        supported_fragment_count += sequence.supported_fragments.size();
    }
    if (assembly.sequences.size() < supported_fragment_count) {
        spdlog::info(
            "HAL export assembled {} supported chains into {} reference intervals for {}",
            supported_fragment_count,
            assembly.sequences.size(),
            model.genome_name);
    }
    uint64_t path_decompose_ms = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - path_decompose_begin).count());
    if (stats != nullptr && node_id == tree.root_id) {
        stats->root_sequence_count = assembly.sequences.size();
        stats->root_singleton_sequence_count = std::count_if(
            assembly.sequences.begin(),
            assembly.sequences.end(),
            [](const auto& sequence) {
                size_t run_count = 0;
                for (const auto& fragment : sequence.supported_fragments) {
                    run_count += fragment.size();
                }
                return run_count == 1;
            });
    }

    auto sequence_materialize_begin = Clock::now();
    model.sequences.reserve(assembly.sequences.size());
    size_t scaffold_index = 0;
    for (const auto& assembled_sequence : assembly.sequences) {
        SequenceModel sequence;
        sequence.seq_name = formatScaffoldName(
            model.genome_name, scaffold_index++);
        if (assembled_sequence.supported_fragments.empty() ||
            assembled_sequence.joins.size() + 1 !=
                assembled_sequence.supported_fragments.size()) {
            throw std::runtime_error(
                "Reference interval has an invalid chain/join layout");
        }
        size_t occurrence_count = 0;
        size_t dna_size = 0;
        for (const auto& fragment :
             assembled_sequence.supported_fragments) {
            occurrence_count += fragment.size();
            for (OccurrenceId occurrence_id : fragment) {
                const uint64_t run_id =
                    model.run_by_occurrence.at(occurrence_id);
                dna_size += model.run_dna.at(run_id).size();
            }
        }
        if (occurrence_count !=
            assembled_sequence
                .supported_fragments.size()) {
            throw std::runtime_error(
                "Reference interval contains a non-elementary supported fragment");
        }
        sequence.joins =
            assembled_sequence.joins;
        for (const auto& join : assembled_sequence.joins) {
            dna_size += join.gap_length;
        }
        sequence.path.reserve(occurrence_count);
        sequence.gaps.reserve(assembled_sequence.joins.size());
        sequence.dna.reserve(dna_size);
        uint64_t cursor = 0;
        for (size_t fragment_index = 0;
             fragment_index <
                 assembled_sequence.supported_fragments.size();
             ++fragment_index) {
            if (fragment_index != 0) {
                const uint32_t gap_length =
                    assembled_sequence
                        .joins[fragment_index - 1]
                        .gap_length;
                if (gap_length != 0) {
                    sequence.gaps.push_back(
                        SequenceGap{cursor, gap_length});
                    sequence.dna.append(gap_length, 'N');
                    cursor += gap_length;
                }
            }
            for (OccurrenceId occurrence_id :
                 assembled_sequence
                     .supported_fragments[fragment_index]) {
                const uint64_t run_id =
                    model.run_by_occurrence.at(occurrence_id);
                const bool forward =
                    assembly.forward_by_occurrence.at(
                        occurrence_id);
                std::string dna = orientRunDNAForPlacement(
                    model.run_dna.at(run_id),
                    forward);
                OrientedOccurrence occurrence{
                    occurrence_id,
                    run_id,
                    model.copy_index_by_occurrence.at(
                        occurrence_id),
                    forward};
                sequence.path.push_back(occurrence);
                model.placements.emplace(
                    occurrence_id,
                    NodePlacement{
                        sequence.seq_name,
                        cursor,
                        static_cast<uint32_t>(dna.size()),
                        forward});
                sequence.dna += dna;
                cursor += dna.size();
            }
        }
        model.sequences.push_back(std::move(sequence));
    }
    uint64_t sequence_materialize_ms = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - sequence_materialize_begin).count());

    if (timings_out != nullptr) {
        timings_out->candidate_occurrence_count =
            candidate_occurrences.size();
        timings_out->edge_count = edges.size();
        timings_out->path_count = model.sequences.size();
        timings_out->candidate_filter_ms = candidate_filter_ms;
        timings_out->run_dna_ms = run_dna_ms;
        timings_out->edge_collect_ms = edge_collect_ms;
        timings_out->path_decompose_ms = path_decompose_ms;
        timings_out->sequence_materialize_ms = sequence_materialize_ms;
    }

    return model;
}

std::unordered_map<int, NodeModel> buildNodeModels(
    const TreeMeta& tree,
    const std::vector<ColumnRun>& runs,
    const std::unordered_map<uint64_t, const ColumnRun*>& run_by_id,
    const std::unordered_map<std::string, std::vector<LeafOccurrence>>& leaf_paths,
    const ExportConfig& config,
    ExportStats* stats) {

    using Clock = std::chrono::steady_clock;
    std::unordered_map<int, NodeModel> models;
    OccurrenceId next_occurrence_id = 1;
    for (int node_id : tree.internal_postorder) {
        NodeModelBuildTimings node_timings;
        auto node_begin = Clock::now();
        NodeModel model = buildNodeModel(
            node_id,
            tree,
            runs,
            run_by_id,
            models,
            leaf_paths,
            config,
            &next_occurrence_id,
            stats,
            &node_timings);
        uint64_t node_ms = static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - node_begin).count());
        node_timings.total_ms = node_ms;
        if (stats != nullptr) {
            if (node_id == tree.root_id) {
                stats->root_build_models_ms += node_ms;
                stats->root_candidate_filter_ms += node_timings.candidate_filter_ms;
                stats->root_run_dna_ms += node_timings.run_dna_ms;
                stats->root_edge_collect_ms += node_timings.edge_collect_ms;
                stats->root_path_decompose_ms += node_timings.path_decompose_ms;
                stats->root_sequence_materialize_ms += node_timings.sequence_materialize_ms;
            } else {
                stats->non_root_build_models_ms += node_ms;
                stats->non_root_candidate_filter_ms += node_timings.candidate_filter_ms;
                stats->non_root_run_dna_ms += node_timings.run_dna_ms;
                stats->non_root_edge_collect_ms += node_timings.edge_collect_ms;
                stats->non_root_path_decompose_ms += node_timings.path_decompose_ms;
                stats->non_root_sequence_materialize_ms += node_timings.sequence_materialize_ms;
            }
        }
        models.emplace(node_id, std::move(model));
    }
    return models;
}

void writeHalFasta(
    const std::filesystem::path& hal_fa_path,
    const std::vector<std::string>& genome_order,
    const std::vector<SequenceEmission>& emissions,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
    const SoftMask::IndexMap& softmask_indexes) {

    std::ofstream out(hal_fa_path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("Cannot open HAL FASTA output: " + hal_fa_path.string());
    }

    std::vector<GenomeSequenceName> genome_sequences;
    genome_sequences.reserve(emissions.size());
    std::unordered_map<std::string, const SequenceEmission*> emission_by_key;
    emission_by_key.reserve(emissions.size());
    for (const auto& emission : emissions) {
        genome_sequences.emplace_back(emission.genome_name, emission.seq_name);
        emission_by_key.emplace(emission.genome_name + '\t' + emission.seq_name, &emission);
    }

    for (const auto& [genome_name, seq_name] : buildOutputSequenceOrder(genome_order, genome_sequences)) {
        auto key = genome_name + '\t' + seq_name;
        auto emission_it = emission_by_key.find(key);
        if (emission_it == emission_by_key.end()) {
            throw std::runtime_error("Missing emission for HAL FASTA output key: " + key);
        }
        const auto& emission = *emission_it->second;
        if (emission.dna) {
            out << '>' << emission.seq_name << '\n' << *emission.dna << '\n';
            continue;
        }
        if (!emission.leaf_source) {
            continue;
        }
        const auto& [species_name, chr_name] = *emission.leaf_source;
        auto mgr_it = seqpro_managers.find(species_name);
        if (mgr_it == seqpro_managers.end()) {
            throw std::runtime_error("Missing SeqPro manager for leaf genome: " + species_name);
        }
        uint64_t chr_length = fetchSequenceLength(mgr_it->second, chr_name);
        std::string dna = fetchSubSequence(mgr_it->second, chr_name, 0, chr_length);
        const auto softmask_it = softmask_indexes.find(species_name);
        if (softmask_it == softmask_indexes.end() || !softmask_it->second) {
            throw std::runtime_error(
                "Missing soft-mask index for leaf genome: " + species_name);
        }
        softmask_it->second->restore(chr_name, 0, dna);
        out << '>' << emission.seq_name << '\n' << dna << '\n';
    }
}

void writeC2H(const std::filesystem::path& c2h_path,
              const std::vector<std::string>& genome_order,
              const std::vector<SequenceEmission>& emissions) {

    std::ofstream out(c2h_path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("Cannot open c2h output: " + c2h_path.string());
    }

    std::vector<GenomeSequenceName> genome_sequences;
    genome_sequences.reserve(emissions.size());
    std::unordered_map<std::string, const SequenceEmission*> emission_by_key;
    emission_by_key.reserve(emissions.size());
    for (const auto& emission : emissions) {
        genome_sequences.emplace_back(emission.genome_name, emission.seq_name);
        emission_by_key.emplace(emission.genome_name + '\t' + emission.seq_name, &emission);
    }

    for (const auto& [genome_name, seq_name] : buildOutputSequenceOrder(genome_order, genome_sequences)) {
        auto key = genome_name + '\t' + seq_name;
        auto emission_it = emission_by_key.find(key);
        if (emission_it == emission_by_key.end()) {
            throw std::runtime_error("Missing emission for c2h output key: " + key);
        }
        const auto& emission = *emission_it->second;
        // The cactus .c2h sequence header stores a 0/1 "has bottoms" flag,
        // not the raw count of bottom segments for the sequence.
        out << "s\t'" << emission.genome_name << "'\t'" << emission.seq_name << "'\t"
            << c2hHasBottomFlag(emission.bottom_count) << "\n";
        for (const auto& bottom : emission.bottoms) {
            out << "a\t" << bottom.name << '\t' << bottom.start << '\t' << bottom.length << "\n";
        }
        for (const auto& top : emission.tops) {
            if (!top.parent_bottom_name) {
                out << "a\t" << top.start << '\t' << top.length << "\n";
            } else {
                out << "a\t" << top.start << '\t' << top.length << '\t'
                    << *top.parent_bottom_name << '\t' << (top.forward_to_parent ? 1 : 0) << "\n";
            }
        }
        out << "\n";
    }
}

std::unordered_map<OccurrenceId, uint64_t>
appendBottomEmissionsForNode(
    int node_id,
    const TreeMeta& tree,
    const std::unordered_map<int, NodeModel>& models,
    std::vector<SequenceEmission>& emissions) {

    auto model_it = models.find(node_id);
    if (model_it == models.end()) {
        throw std::runtime_error(
            "Missing node model for local subtree root");
    }

    std::unordered_map<OccurrenceId, uint64_t> bottom_names;
    uint64_t next_bottom_name = 1;
    const auto& model = model_it->second;
    for (const auto& sequence : model.sequences) {
        SequenceEmission emission;
        emission.genome_name = tree.nodes[node_id].name;
        emission.seq_name = sequence.seq_name;
        emission.bottom_count = sequence.path.size() + sequence.gaps.size();
        emission.dna = sequence.dna;
        for (const auto& occurrence : sequence.path) {
            const auto& placement =
                model.placements.at(occurrence.occurrence_id);
            emission.bottoms.push_back(
                BottomSegmentLine{
                    occurrence.occurrence_id,
                    next_bottom_name,
                    placement.start,
                    placement.length});
            bottom_names.emplace(
                occurrence.occurrence_id,
                next_bottom_name);
            ++next_bottom_name;
        }
        for (const auto& gap : sequence.gaps) {
            emission.bottoms.push_back(
                BottomSegmentLine{
                    0,
                    next_bottom_name,
                    gap.start,
                    gap.length});
            ++next_bottom_name;
        }
        std::sort(
            emission.bottoms.begin(),
            emission.bottoms.end(),
            [](const auto& lhs, const auto& rhs) {
                return lhs.start < rhs.start;
            });
        emissions.push_back(std::move(emission));
    }
    return bottom_names;
}
OccurrenceId mapChildCopyToParentOccurrence(
    const NodeModel& parent_model,
    int child_id,
    uint64_t run_id,
    uint32_t child_copy_index) {
    const auto child_it =
        parent_model.parent_copy_by_child_run.find(
            child_id);
    const auto occurrences_it =
        parent_model.occurrences_by_run.find(run_id);
    if (child_it ==
            parent_model.parent_copy_by_child_run.end() ||
        occurrences_it ==
            parent_model.occurrences_by_run.end()) {
        throw std::runtime_error(
            "Missing reconciled child-to-parent occurrence mapping");
    }
    const auto run_it = child_it->second.find(run_id);
    if (run_it == child_it->second.end() ||
        child_copy_index >= run_it->second.size() ||
        run_it->second[child_copy_index] >=
            occurrences_it->second.size()) {
        throw std::runtime_error(
            "Invalid reconciled child-to-parent occurrence mapping");
    }
    return occurrences_it->second[
        run_it->second[child_copy_index]];
}


struct ParentProjectedFragment {
    std::vector<OrientedOccurrence> path;
    std::vector<ReferenceJoin> joins;
    size_t parent_sequence_rank = 0;
    uint64_t parent_begin = 0;
    uint64_t parent_end = 0;
    size_t source_sequence_rank = 0;
    size_t source_fragment_rank = 0;
};


void projectInternalChildContainer(
    int parent_id,
    int child_id,
    std::unordered_map<int, NodeModel>& models,
    uint32_t scaffold_gap_length,
    ExportStats* stats) {
    const auto parent_it = models.find(parent_id);
    const auto child_it = models.find(child_id);
    if (parent_it == models.end() ||
        child_it == models.end()) {
        throw std::runtime_error(
            "Missing node model for parent-constrained reference projection");
    }
    const NodeModel& parent_model =
        parent_it->second;
    NodeModel& child_model =
        child_it->second;

    std::unordered_map<std::string, size_t>
        parent_sequence_rank;
    parent_sequence_rank.reserve(
        parent_model.sequences.size());
    for (size_t rank = 0;
         rank < parent_model.sequences.size();
         ++rank) {
        const auto& sequence =
            parent_model.sequences[rank];
        if (!parent_sequence_rank
                 .emplace(sequence.seq_name, rank)
                 .second) {
            throw std::runtime_error(
                "Duplicate parent reference sequence name");
        }
    }

    const size_t old_sequence_count =
        child_model.sequences.size();
    uint64_t old_scaffold_join_count = 0;
    uint64_t old_scaffold_gap_bases = 0;
    size_t projected_occurrence_count = 0;
    std::vector<ParentProjectedFragment>
        fragments;
    fragments.reserve(
        child_model.sequences.size());

    std::function<void(
        SequenceModel&,
        size_t,
        size_t,
        size_t,
        size_t)>
        add_fragment;
    add_fragment =
        [&](SequenceModel& source,
            size_t sequence_rank,
            size_t fragment_rank,
            size_t begin,
            size_t end) {
            if (begin >= end ||
                end > source.path.size()) {
                throw std::logic_error(
                    "Invalid child reference fragment bounds");
            }
            ParentProjectedFragment fragment;
            fragment.path.reserve(end - begin);
            fragment.joins.reserve(
                end - begin - 1);
            fragment.source_sequence_rank =
                sequence_rank;
            fragment.source_fragment_rank =
                fragment_rank;
            auto find_parent_placement =
                [&](const OrientedOccurrence& occurrence)
                    -> const NodePlacement& {
                    const OccurrenceId
                        parent_occurrence_id =
                            mapChildCopyToParentOccurrence(
                                parent_model,
                                child_id,
                                occurrence.run_id,
                                occurrence.copy_index);
                    const auto placement_it =
                        parent_model.placements.find(
                            parent_occurrence_id);
                    if (placement_it ==
                        parent_model.placements.end()) {
                        throw std::runtime_error(
                            "Projected child occurrence has no parent placement");
                    }
                    return placement_it->second;
                };
            auto parent_rank =
                [&](const OrientedOccurrence& occurrence) {
                    const auto& placement =
                        find_parent_placement(
                            occurrence);
                    const auto rank_it =
                        parent_sequence_rank.find(
                            placement.seq_name);
                    if (rank_it ==
                        parent_sequence_rank.end()) {
                        throw std::runtime_error(
                            "Projected child occurrence references an unknown parent sequence");
                    }
                    return rank_it->second;
                };
            const size_t first_parent_rank =
                parent_rank(
                    source.path[begin]);
            size_t bounded_end = begin + 1;
            while (bounded_end < end &&
                   parent_rank(
                       source.path[bounded_end]) ==
                       first_parent_rank) {
                ++bounded_end;
            }
            if (bounded_end < end) {
                add_fragment(
                    source,
                    sequence_rank,
                    fragment_rank,
                    begin,
                    bounded_end);
                add_fragment(
                    source,
                    sequence_rank,
                    fragment_rank + 1,
                    bounded_end,
                    end);
                return;
            }
            for (size_t index = begin;
                 index < end;
                 ++index) {
                fragment.path.push_back(
                    std::move(
                        source.path[index]));
            }
            for (size_t index = begin;
                 index + 1 < end;
                 ++index) {
                if (source.joins[index].kind ==
                    ReferenceJoinKind::SCAFFOLD) {
                    throw std::logic_error(
                        "Child reference fragment retained an artificial scaffold join");
                }
                fragment.joins.push_back(
                    std::move(
                        source.joins[index]));
            }

            bool parent_rank_initialized =
                false;
            bool saw_increasing = false;
            bool saw_decreasing = false;
            uint64_t previous_parent_start = 0;
            bool have_previous_parent = false;
            fragment.parent_begin =
                std::numeric_limits<uint64_t>::max();
            fragment.parent_end = 0;
            for (const auto& occurrence :
                 fragment.path) {
                const auto& parent_placement =
                    find_parent_placement(
                        occurrence);
                const size_t occurrence_parent_rank =
                    parent_rank(
                        occurrence);
                if (!parent_rank_initialized) {
                    fragment
                        .parent_sequence_rank =
                        occurrence_parent_rank;
                    parent_rank_initialized =
                        true;
                } else if (
                    fragment
                        .parent_sequence_rank !=
                    occurrence_parent_rank) {
                    throw std::logic_error(
                        "Parent-boundary fragment splitting failed");
                }
                const uint64_t parent_start =
                    parent_placement.start;
                const uint64_t parent_length =
                    parent_placement.length;
                if (parent_start >
                    std::numeric_limits<
                        uint64_t>::max() -
                        parent_length) {
                    throw std::overflow_error(
                        "Projected parent occurrence coordinate overflow");
                }
                fragment.parent_begin =
                    std::min(
                        fragment.parent_begin,
                        parent_start);
                fragment.parent_end =
                    std::max(
                        fragment.parent_end,
                        parent_start +
                            parent_length);
                if (have_previous_parent) {
                    saw_increasing =
                        saw_increasing ||
                        parent_start >
                            previous_parent_start;
                    saw_decreasing =
                        saw_decreasing ||
                        parent_start <
                            previous_parent_start;
                }
                previous_parent_start =
                    parent_start;
                have_previous_parent = true;
            }
            if (!parent_rank_initialized) {
                throw std::logic_error(
                    "Parent-constrained projection emitted an empty child fragment");
            }
            if (saw_decreasing &&
                !saw_increasing) {
                std::reverse(
                    fragment.path.begin(),
                    fragment.path.end());
                std::reverse(
                    fragment.joins.begin(),
                    fragment.joins.end());
                for (auto& occurrence :
                     fragment.path) {
                    occurrence.forward =
                        !occurrence.forward;
                }
            }
            projected_occurrence_count +=
                fragment.path.size();
            fragments.push_back(
                std::move(fragment));
        };

    for (size_t sequence_rank = 0;
         sequence_rank <
         child_model.sequences.size();
         ++sequence_rank) {
        auto& sequence =
            child_model.sequences[
                sequence_rank];
        if (sequence.path.empty() ||
            sequence.path.size() !=
                sequence.joins.size() + 1) {
            throw std::runtime_error(
                "Child reference sequence has an invalid path/join layout");
        }
        size_t fragment_begin = 0;
        size_t fragment_rank = 0;
        for (size_t join_index = 0;
             join_index <
             sequence.joins.size();
             ++join_index) {
            const auto& join =
                sequence.joins[join_index];
            if (join.kind !=
                ReferenceJoinKind::SCAFFOLD) {
                continue;
            }
            ++old_scaffold_join_count;
            old_scaffold_gap_bases +=
                join.gap_length;
            add_fragment(
                sequence,
                sequence_rank,
                fragment_rank++,
                fragment_begin,
                join_index + 1);
            fragment_begin =
                join_index + 1;
        }
        add_fragment(
            sequence,
            sequence_rank,
            fragment_rank,
            fragment_begin,
            sequence.path.size());
    }

    if (projected_occurrence_count !=
        child_model
            .run_by_occurrence.size()) {
        throw std::runtime_error(
            "Parent-constrained projection did not retain every child occurrence");
    }

    std::map<
        size_t,
        std::vector<ParentProjectedFragment>>
        fragments_by_parent_sequence;
    for (auto& fragment : fragments) {
        fragments_by_parent_sequence[
            fragment.parent_sequence_rank]
                .push_back(
                    std::move(fragment));
    }
    std::vector<ParentProjectedFragment>{}
        .swap(fragments);
    std::vector<SequenceModel>{}.swap(
        child_model.sequences);
    std::unordered_map<
        OccurrenceId,
        NodePlacement>{}
        .swap(child_model.placements);

    uint64_t new_scaffold_join_count = 0;
    uint64_t new_scaffold_gap_bases = 0;
    child_model.sequences.reserve(
        fragments_by_parent_sequence.size());
    size_t scaffold_index = 0;
    for (auto& [parent_rank,
                parent_fragments] :
         fragments_by_parent_sequence) {
        (void)parent_rank;
        std::sort(
            parent_fragments.begin(),
            parent_fragments.end(),
            [](const auto& lhs,
               const auto& rhs) {
                return std::tie(
                           lhs.parent_begin,
                           lhs.parent_end,
                           lhs.source_sequence_rank,
                           lhs.source_fragment_rank) <
                       std::tie(
                           rhs.parent_begin,
                           rhs.parent_end,
                           rhs.source_sequence_rank,
                           rhs.source_fragment_rank);
            });

        SequenceModel projected;
        projected.seq_name =
            formatScaffoldName(
                child_model.genome_name,
                scaffold_index++);
        size_t path_size = 0;
        size_t dna_size = 0;
        for (const auto& fragment :
             parent_fragments) {
            path_size +=
                fragment.path.size();
            for (const auto& occurrence :
                 fragment.path) {
                const size_t run_length =
                    child_model.run_dna
                        .at(occurrence.run_id)
                        .size();
                if (dna_size >
                    std::numeric_limits<
                        size_t>::max() -
                        run_length) {
                    throw std::overflow_error(
                        "Projected child reference DNA size overflow");
                }
                dna_size += run_length;
            }
            for (const auto& join :
                 fragment.joins) {
                if (dna_size >
                    std::numeric_limits<
                        size_t>::max() -
                        join.gap_length) {
                    throw std::overflow_error(
                        "Projected child reference gap size overflow");
                }
                dna_size +=
                    join.gap_length;
            }
        }
        if (parent_fragments.size() > 1) {
            const size_t bridge_count =
                parent_fragments.size() - 1;
            if (bridge_count >
                std::numeric_limits<
                    size_t>::max() /
                    scaffold_gap_length ||
                dna_size >
                    std::numeric_limits<
                        size_t>::max() -
                        bridge_count *
                            scaffold_gap_length) {
                throw std::overflow_error(
                    "Projected child reference scaffold size overflow");
            }
            dna_size +=
                bridge_count *
                scaffold_gap_length;
        }
        projected.path.reserve(path_size);
        if (path_size != 0) {
            projected.joins.reserve(
                path_size - 1);
        }
        projected.gaps.reserve(
            path_size);
        projected.dna.reserve(dna_size);
        uint64_t cursor = 0;

        auto append_join =
            [&](const ReferenceJoin& join) {
                if (projected.path.empty()) {
                    throw std::logic_error(
                        "Projected child reference begins with a join");
                }
                projected.joins.push_back(
                    join);
                if (join.gap_length != 0) {
                    projected.gaps.push_back(
                        SequenceGap{
                            cursor,
                            join.gap_length});
                    projected.dna.append(
                        join.gap_length,
                        'N');
                    cursor +=
                        join.gap_length;
                }
            };
        bool first_fragment = true;
        for (auto& fragment :
             parent_fragments) {
            if (!first_fragment) {
                append_join(
                    ReferenceJoin{
                        ReferenceJoinKind::
                            SCAFFOLD,
                        scaffold_gap_length,
                        0.0L});
                ++new_scaffold_join_count;
                new_scaffold_gap_bases +=
                    scaffold_gap_length;
            }
            first_fragment = false;
            for (size_t occurrence_index = 0;
                 occurrence_index <
                 fragment.path.size();
                 ++occurrence_index) {
                if (occurrence_index != 0) {
                    append_join(
                        fragment.joins[
                            occurrence_index -
                            1]);
                }
                auto occurrence =
                    std::move(
                        fragment.path[
                            occurrence_index]);
                const auto& canonical_dna =
                    child_model.run_dna.at(
                        occurrence.run_id);
                std::string dna =
                    orientRunDNAForPlacement(
                        canonical_dna,
                        occurrence.forward);
                auto [placement_it,
                      inserted] =
                    child_model
                        .placements.emplace(
                            occurrence
                                .occurrence_id,
                            NodePlacement{
                                projected.seq_name,
                                cursor,
                                static_cast<
                                    uint32_t>(
                                    dna.size()),
                                occurrence
                                    .forward});
                (void)placement_it;
                if (!inserted) {
                    throw std::runtime_error(
                        "Parent-constrained projection duplicated a child occurrence");
                }
                projected.path.push_back(
                    std::move(occurrence));
                projected.dna += dna;
                cursor += dna.size();
            }
        }
        if (projected.path.empty() ||
            projected.path.size() !=
                projected.joins.size() + 1 ||
            projected.dna.size() != cursor) {
            throw std::runtime_error(
                "Parent-constrained projection emitted an invalid child reference sequence");
        }
        child_model.sequences.push_back(
            std::move(projected));
    }

    if (child_model.placements.size() !=
        child_model
            .run_by_occurrence.size()) {
        throw std::runtime_error(
            "Parent-constrained projection lost a child occurrence placement");
    }
    if (stats != nullptr) {
        if (stats->reference_interval_count <
                old_sequence_count ||
            stats->scaffold_join_count <
                old_scaffold_join_count ||
            stats->scaffold_gap_bases <
                old_scaffold_gap_bases) {
            throw std::logic_error(
                "HAL export statistics underflow during parent-constrained projection");
        }
        stats->reference_interval_count -=
            old_sequence_count;
        stats->reference_interval_count +=
            child_model.sequences.size();
        stats->scaffold_join_count -=
            old_scaffold_join_count;
        stats->scaffold_join_count +=
            new_scaffold_join_count;
        stats->scaffold_gap_bases -=
            old_scaffold_gap_bases;
        stats->scaffold_gap_bases +=
            new_scaffold_gap_bases;
    }
    spdlog::info(
        "HAL export projected {} supported fragments into {} parent-constrained reference intervals for {}",
        old_scaffold_join_count +
            old_sequence_count,
        child_model.sequences.size(),
        child_model.genome_name);
}


void projectInternalReferenceContainersTopDown(
    const TreeMeta& tree,
    const ExportConfig& config,
    std::unordered_map<int, NodeModel>& models,
    ExportStats* stats) {
    std::function<void(int)> project_children =
        [&](int parent_id) {
            for (int child_id :
                 tree.nodes[
                     static_cast<size_t>(
                         parent_id)]
                     .children) {
                if (tree.nodes[
                        static_cast<size_t>(
                            child_id)]
                        .is_leaf) {
                    continue;
                }
                projectInternalChildContainer(
                    parent_id,
                    child_id,
                    models,
                    config.scaffold_gap_length,
                    stats);
                project_children(child_id);
            }
        };
    project_children(tree.root_id);
}


void appendInternalChildTopEmissions(
    int parent_id,
    int child_id,
    const TreeMeta& tree,
    const std::unordered_map<int, NodeModel>& models,
    const std::unordered_map<uint64_t, const ColumnRun*>& run_by_id,
    const std::unordered_map<OccurrenceId, uint64_t>& bottom_names,
    std::vector<SequenceEmission>& emissions) {

    auto model_it = models.find(child_id);
    if (model_it == models.end()) {
        throw std::runtime_error("Missing child node model for local subtree");
    }
    const auto& model = model_it->second;
    const auto& parent_model = models.at(parent_id);

    for (const auto& sequence : model.sequences) {
        SequenceEmission emission;
        emission.genome_name = model.genome_name;
        emission.seq_name = sequence.seq_name;
        emission.bottom_count = 0;
        emission.dna = sequence.dna;
        for (const auto& occurrence : sequence.path) {
            const auto& placement =
                model.placements.at(occurrence.occurrence_id);
            auto run_it = run_by_id.find(occurrence.run_id);
            if (run_it == run_by_id.end()) {
                throw std::runtime_error(
                    "Missing run metadata while building local internal child emissions");
            }
            const auto& run = *run_it->second;
            if (run.present_by_node[parent_id]) {
                const auto parent_occurrences_it =
                    parent_model.occurrences_by_run.find(
                        occurrence.run_id);
                if (parent_occurrences_it ==
                        parent_model.occurrences_by_run.end() ||
                    parent_occurrences_it->second.empty()) {
                    throw std::runtime_error(
                        "Missing parent occurrence for local aligned run");
                }
                const OccurrenceId parent_occurrence_id =
                    mapChildCopyToParentOccurrence(
                        parent_model,
                        child_id,
                        occurrence.run_id,
                        occurrence.copy_index);
                auto name_it =
                    bottom_names.find(parent_occurrence_id);
                if (name_it == bottom_names.end()) {
                    throw std::runtime_error(
                        "Missing parent bottom segment name for local aligned occurrence");
                }
                const auto& parent_placement =
                    parent_model.placements.at(
                        parent_occurrence_id);
                emission.tops.push_back(
                    TopSegmentLine{
                        occurrence.occurrence_id,
                        placement.start,
                        placement.length,
                        name_it->second,
                        computeForwardToParent(
                            occurrence.forward,
                            parent_placement.forward)});
            } else {
                emission.tops.push_back(
                    TopSegmentLine{
                        occurrence.occurrence_id,
                        placement.start,
                        placement.length,
                        std::nullopt,
                        true});
            }
        }
        for (const auto& gap : sequence.gaps) {
            emission.tops.push_back(
                TopSegmentLine{
                    0,
                    gap.start,
                    gap.length,
                    std::nullopt,
                    true});
        }
        std::sort(
            emission.tops.begin(),
            emission.tops.end(),
            [](const auto& lhs, const auto& rhs) {
                return lhs.start < rhs.start;
            });
        emissions.push_back(std::move(emission));
    }
}

void appendLeafChildTopEmissions(
    int parent_id,
    int leaf_id,
    const TreeMeta& tree,
    const std::unordered_map<int, NodeModel>& models,
    const std::unordered_map<std::string, std::vector<LeafOccurrence>>& leaf_paths,
    const std::unordered_map<uint64_t, const ColumnRun*>& run_by_id,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
    const std::unordered_map<OccurrenceId, uint64_t>& bottom_names,
    std::vector<SequenceEmission>& emissions) {

    const std::string& species_name = tree.nodes[leaf_id].name;
    auto mgr_it = seqpro_managers.find(species_name);
    if (mgr_it == seqpro_managers.end()) {
        throw std::runtime_error("Missing SeqPro manager for leaf child genome: " + species_name);
    }
    const auto& parent_model = models.at(parent_id);

    std::unordered_map<std::string, std::vector<LeafOccurrence>> windows_by_sequence;
    auto path_it = leaf_paths.find(species_name);
    if (path_it != leaf_paths.end()) {
        for (const auto& occurrence : path_it->second) {
            windows_by_sequence[occurrence.span.hal_sequence_name].push_back(occurrence);
        }
    }

    for (const auto& chr_name : fetchSequenceNames(mgr_it->second)) {
        const std::string qualified_sequence_name =
            makeQualifiedLeafSequenceName(species_name, chr_name);
        SequenceEmission emission;
        emission.genome_name = species_name;
        emission.seq_name = chr_name;
        emission.bottom_count = 0;
        emission.leaf_source = std::make_pair(species_name, chr_name);

        uint64_t chr_length = fetchSequenceLength(mgr_it->second, chr_name);
        auto& windows = windows_by_sequence[qualified_sequence_name];
        std::sort(windows.begin(), windows.end(), [](const auto& lhs, const auto& rhs) {
            if (lhs.span.start != rhs.span.start) {
                return lhs.span.start < rhs.span.start;
            }
            return lhs.run_id < rhs.run_id;
        });

        uint64_t cursor = 0;
        for (const auto& occurrence : windows) {
            const auto& span = occurrence.span;
            if (span.start > cursor) {
                emission.tops.push_back(
                    TopSegmentLine{0, cursor, static_cast<uint32_t>(span.start - cursor), std::nullopt, true});
            }
            const auto& run =
                *run_by_id.at(occurrence.run_id);
            if (run.present_by_node[parent_id]) {
                const auto parent_occurrences_it =
                    parent_model.occurrences_by_run.find(
                        occurrence.run_id);
                if (parent_occurrences_it ==
                        parent_model.occurrences_by_run.end() ||
                    parent_occurrences_it->second.empty()) {
                    throw std::runtime_error(
                        "Missing parent occurrence for local leaf-aligned run");
                }
                const OccurrenceId parent_occurrence_id =
                    mapChildCopyToParentOccurrence(
                        parent_model,
                        leaf_id,
                        occurrence.run_id,
                        occurrence.copy_index);
                auto name_it =
                    bottom_names.find(parent_occurrence_id);
                if (name_it == bottom_names.end()) {
                    throw std::runtime_error(
                        "Missing parent bottom segment name for local leaf-aligned occurrence");
                }
                const auto& parent_placement =
                    parent_model.placements.at(
                        parent_occurrence_id);
                emission.tops.push_back(
                    TopSegmentLine{
                        0,
                        span.start,
                        span.length,
                        name_it->second,
                        computeForwardToParent(
                            occurrence.forward_to_canonical,
                            parent_placement.forward)});
            } else {
                emission.tops.push_back(
                    TopSegmentLine{
                        0,
                        span.start,
                        span.length,
                        std::nullopt,
                        true});
            }
            cursor = span.start + span.length;
        }
        if (cursor < chr_length) {
            emission.tops.push_back(
                TopSegmentLine{0, cursor, static_cast<uint32_t>(chr_length - cursor), std::nullopt, true});
        }
        emissions.push_back(std::move(emission));
    }
}

std::vector<SequenceEmission> buildLocalSubtreeEmissions(
    int node_id,
    const TreeMeta& tree,
    const std::unordered_map<int, NodeModel>& models,
    const std::unordered_map<std::string, std::vector<LeafOccurrence>>& leaf_paths,
    const std::unordered_map<uint64_t, const ColumnRun*>& run_by_id,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {

    std::vector<SequenceEmission> emissions;
    auto bottom_names = appendBottomEmissionsForNode(node_id, tree, models, emissions);
    for (int child_id : tree.nodes[node_id].children) {
        if (tree.nodes[child_id].is_leaf) {
            appendLeafChildTopEmissions(
                node_id, child_id, tree, models, leaf_paths, run_by_id, seqpro_managers, bottom_names, emissions);
        } else {
            appendInternalChildTopEmissions(
                node_id, child_id, tree, models, run_by_id, bottom_names, emissions);
        }
    }
    return emissions;
}

std::vector<std::string> buildLocalGenomeOrder(const TreeMeta& tree, int node_id) {
    std::vector<std::string> order;
    order.reserve(tree.nodes[node_id].children.size() + 1);
    order.push_back(tree.nodes[node_id].name);
    for (int child_id : tree.nodes[node_id].children) {
        order.push_back(tree.nodes[child_id].name);
    }
    return order;
}

std::string buildLocalNewick(const TreeMeta& tree, int node_id) {
    const auto& node = tree.nodes[node_id];
    std::ostringstream oss;
    if (!node.children.empty()) {
        oss << '(';
        for (size_t i = 0; i < node.children.size(); ++i) {
            if (i != 0) {
                oss << ',';
            }
            const auto& child = tree.nodes[static_cast<size_t>(node.children[i])];
            oss << child.name << ':' << std::fixed << std::setprecision(6) << child.branch_length_to_parent;
        }
        oss << ')';
    }
    oss << node.name << ';';
    return oss.str();
}

void accumulateEmissionStats(
    const std::vector<SequenceEmission>& emissions,
    ExportStats* stats) {

    if (stats == nullptr) {
        return;
    }
    std::unordered_map<std::string, std::unordered_map<uint64_t, uint64_t>>
        top_count_by_child_and_parent_bottom;
    for (const auto& emission : emissions) {
        for (const auto& top : emission.tops) {
            if (top.parent_bottom_name) {
                ++stats->aligned_top_count;
                ++top_count_by_child_and_parent_bottom[emission.genome_name][*top.parent_bottom_name];
                if (!top.forward_to_parent) {
                    ++stats->reverse_top_count;
                }
            } else {
                ++stats->insertion_top_count;
            }
        }
    }
    for (const auto& [child_name, counts] : top_count_by_child_and_parent_bottom) {
        (void)child_name;
        for (const auto& [parent_bottom, count] : counts) {
            (void)parent_bottom;
            if (count > 1) {
                stats->paralogous_top_count += count - 1;
            }
        }
    }
}

void writeTreeFile(const std::filesystem::path& tree_path, const std::string& newick) {
    std::ofstream out(tree_path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("Cannot open tree output: " + tree_path.string());
    }
    out << newick << '\n';
}

void runHalAppend(const std::filesystem::path& c2h_path,
                  const std::filesystem::path& hal_fa_path,
                  const std::filesystem::path& tree_path,
                  const std::filesystem::path& hal_path) {
    const auto executable =
        RaMAxDependencies::locateHalAppendCactusSubtreeExecutable();
    if (executable.empty()) {
        throw std::runtime_error(
            "halAppendCactusSubtree is required but is no longer available");
    }

    std::filesystem::path stdout_path = hal_path;
    stdout_path += ".halAppendCactusSubtree.stdout.log";
    std::filesystem::path stderr_path = hal_path;
    stderr_path += ".halAppendCactusSubtree.stderr.log";
    const auto result = RaMAxExternalTool::run(
        executable,
        {
            c2h_path.string(),
            hal_fa_path.string(),
            tree_path.string(),
            hal_path.string(),
            "--hdf5InMemory",
        },
        stdout_path,
        stderr_path);
    if (result.exit_code != 0) {
        std::ostringstream error;
        error << "halAppendCactusSubtree failed with exit code "
              << result.exit_code << "; stderr: " << stderr_path;
        throw std::runtime_error(error.str());
    }
}


} // namespace

void exportToMaf(
    const std::vector<std::weak_ptr<Block>>& blocks,
    const std::filesystem::path& maf_path,
    const std::map<
        SpeciesName,
        SeqPro::SharedManagerVariant>&
        seqpro_managers,
    bool pairwise_mode) {
    exportCanonicalMafImpl(
        blocks,
        maf_path,
        seqpro_managers,
        pairwise_mode);
}

std::unordered_set<uint64_t>
findRejectedSecondaryHomologyBlocks(
    const std::vector<std::weak_ptr<Block>>& blocks,
    const std::map<
        SpeciesName,
        SeqPro::SharedManagerVariant>&
        seqpro_managers) {
    return findRejectedSecondaryHomologyBlocksImpl(
        blocks,
        seqpro_managers);
}

TreeMeta buildTreeMeta(const NewickParser& parser) {
    TreeMeta meta;
    const auto& nodes = parser.getNodes();
    meta.nodes.reserve(nodes.size());
    for (const auto& node : nodes) {
        TreeNodeMeta meta_node;
        meta_node.id = node.id;
        meta_node.name = node.name.empty() ? ("internal_" + std::to_string(node.id)) : node.name;
        meta_node.parent = node.father;
        meta_node.branch_length_to_parent = node.branchLength;
        meta_node.is_leaf = node.isLeaf;
        meta.nodes.push_back(std::move(meta_node));
    }

    for (auto& node : meta.nodes) {
        meta.name_to_id.emplace(node.name, node.id);
        if (node.parent == -1) {
            meta.root_id = node.id;
        }
    }
    for (auto& node : meta.nodes) {
        if (node.parent != -1) {
            meta.nodes[static_cast<size_t>(node.parent)].children.push_back(node.id);
        }
    }

    int next_leaf_index = 0;
    for (auto& node : meta.nodes) {
        if (node.is_leaf) {
            node.leaf_index = next_leaf_index++;
            meta.leaf_ids.push_back(node.id);
        }
    }
    if (meta.root_id == -1) {
        throw std::runtime_error("TreeMeta build failed: root node not found");
    }
    meta.internal_postorder = computePostorderInternal(meta, meta.root_id);
    return meta;
}


BinaryInferenceResult inferDescendantUnion(
    const TreeMeta& tree,
    const std::unordered_map<std::string, bool>& leaf_presence) {
    std::vector<uint8_t> leaf_states(tree.leaf_ids.size(), 0);
    for (int leaf_id : tree.leaf_ids) {
        const auto& leaf = tree.nodes[leaf_id];
        auto state_it = leaf_presence.find(leaf.name);
        if (state_it != leaf_presence.end() && state_it->second) {
            leaf_states[static_cast<size_t>(leaf.leaf_index)] = 1;
        }
    }

    auto fast = inferDescendantUnionFast(tree, leaf_states);
    BinaryInferenceResult result;
    result.present_by_node = std::move(fast.present_by_node);
    result.margin = std::move(fast.margin);
    result.score0.resize(tree.nodes.size(), 0.0);
    result.score1.resize(tree.nodes.size(), 0.0);
    for (const auto& node : tree.nodes) {
        const uint8_t present = result.present_by_node[static_cast<size_t>(node.id)];
        result.score0[static_cast<size_t>(node.id)] = present ? 1.0 : 0.0;
        result.score1[static_cast<size_t>(node.id)] = present ? 0.0 : 1.0;
        result.present_by_name.emplace(node.name, present);
    }
    return result;
}

LeafInterval projectLeafInterval(
    uint64_t segment_start,
    uint32_t segment_length,
    bool reversed,
    uint32_t non_gap_before,
    uint32_t run_length) {

    LeafInterval interval;
    interval.length = run_length;
    if (!reversed) {
        interval.start = segment_start + non_gap_before;
        interval.forward_to_parent = true;
    } else {
        interval.start = segment_start + (segment_length - non_gap_before - run_length);
        interval.forward_to_parent = false;
    }
    return interval;
}

std::vector<ElementaryRunProjection> projectElementaryRuns(
    const std::vector<AlignedOccurrence>& rows) {
    if (rows.empty()) {
        return {};
    }

    const size_t column_count = rows.front().aligned_dna.size();
    if (column_count == 0) {
        return {};
    }

    std::unordered_set<std::string> row_ids;
    row_ids.reserve(rows.size());
    std::vector<std::vector<uint32_t>> prefixes(rows.size());
    for (size_t row_index = 0; row_index < rows.size(); ++row_index) {
        const auto& row = rows[row_index];
        if (!row_ids.insert(row.row_id).second) {
            throw std::runtime_error("Duplicate aligned occurrence row id: " + row.row_id);
        }
        if (row.aligned_dna.size() != column_count) {
            throw std::runtime_error("Aligned occurrence rows have inconsistent column counts");
        }
        auto& prefix = prefixes[row_index];
        prefix.assign(column_count + 1, 0);
        for (size_t col = 0; col < column_count; ++col) {
            prefix[col + 1] = prefix[col] + (row.aligned_dna[col] != '-' ? 1u : 0u);
        }
        if (prefix.back() != row.segment_length) {
            throw std::runtime_error(
                "Aligned occurrence non-gap length differs from its source segment: " + row.row_id);
        }
    }

    auto same_mask = [&](size_t lhs, size_t rhs) {
        for (const auto& row : rows) {
            if ((row.aligned_dna[lhs] != '-') != (row.aligned_dna[rhs] != '-')) {
                return false;
            }
        }
        return true;
    };

    std::vector<ElementaryRunProjection> projections;
    size_t col = 0;
    while (col < column_count) {
        size_t next = col + 1;
        while (next < column_count && same_mask(col, next)) {
            ++next;
        }

        ElementaryRunProjection projection;
        projection.col_beg = static_cast<uint32_t>(col);
        projection.col_end = static_cast<uint32_t>(next);
        projection.occurrences.reserve(rows.size());
        for (size_t row_index = 0; row_index < rows.size(); ++row_index) {
            const auto& row = rows[row_index];
            if (row.aligned_dna[col] == '-') {
                continue;
            }
            uint32_t non_gap_before = prefixes[row_index][col];
            uint32_t run_length = prefixes[row_index][next] - non_gap_before;
            LeafInterval interval = projectLeafInterval(
                row.segment_start,
                row.segment_length,
                row.reversed,
                non_gap_before,
                run_length);

            ProjectedOccurrence occurrence;
            occurrence.row_id = row.row_id;
            occurrence.genome_name = row.genome_name;
            occurrence.sequence_name = row.sequence_name;
            occurrence.start = interval.start;
            occurrence.length = interval.length;
            occurrence.reversed = !interval.forward_to_parent;
            occurrence.dna = stripGaps(row.aligned_dna.substr(col, next - col));
            projection.occurrences.push_back(std::move(occurrence));
        }
        if (!projection.occurrences.empty()) {
            projections.push_back(std::move(projection));
        }
        col = next;
    }
    return projections;
}

int c2hHasBottomFlag(size_t bottom_count) {
    return bottom_count == 0 ? 0 : 1;
}


bool computeForwardToParent(bool child_forward_to_canonical, bool parent_forward_to_canonical) {
    return child_forward_to_canonical == parent_forward_to_canonical;
}

std::optional<AdjacencyVote> orientAdjacencyVote(
    uint64_t left_run_id,
    bool left_forward_to_canonical,
    uint64_t right_run_id,
    bool right_forward_to_canonical,
    uint32_t left_length,
    uint32_t right_length,
    uint64_t gap_bases) {

    if (left_run_id == right_run_id) {
        return std::nullopt;
    }
    return AdjacencyVote{
        left_run_id,
        left_forward_to_canonical,
        right_run_id,
        right_forward_to_canonical,
        left_length,
        right_length,
        gap_bases};
}
long double calculateCactusZScore(
    uint64_t left_length,
    uint64_t right_length,
    uint64_t gap_length,
    long double theta) {

    if (theta < 0.0L || theta > 1.0L) {
        throw std::invalid_argument("Cactus adjacency theta must be in [0, 1]");
    }
    if (theta == 0.0L) {
        return static_cast<long double>(left_length) *
               static_cast<long double>(right_length);
    }
    const long double beta = 1.0L - theta;
    return ((1.0L - std::pow(beta, static_cast<long double>(left_length))) / theta) *
           std::pow(beta, static_cast<long double>(gap_length)) *
           ((1.0L - std::pow(beta, static_cast<long double>(right_length))) / theta);
}


std::string orientRunDNAForPlacement(const std::string& canonical_dna, bool forward) {
    std::string dna = canonical_dna;
    if (!forward) {
        reverseComplement(dna);
    }
    return dna;
}

std::string buildConsensusDNA(
    const std::vector<std::pair<std::string, double>>& donors,
    size_t expected_length,
    double consensus_threshold) {

    return buildConsensusDNAImpl(donors, expected_length, consensus_threshold);
}

std::vector<std::pair<std::string, double>> selectBestDonorsByBucket(
    const std::vector<BucketedDonor>& donors) {

    std::map<int, BucketedDonor> best_by_bucket;
    for (const auto& donor : donors) {
        if (donor.dna.empty()) {
            continue;
        }
        auto it = best_by_bucket.find(donor.bucket_id);
        if (it == best_by_bucket.end()) {
            best_by_bucket.emplace(donor.bucket_id, donor);
            continue;
        }
        const auto& current = it->second;
        bool replace = false;
        if (donor.prefer_internal != current.prefer_internal) {
            replace = donor.prefer_internal;
        } else if (donor.weight != current.weight) {
            replace = donor.weight > current.weight;
        } else if (donor.dna != current.dna) {
            replace = donor.dna < current.dna;
        }
        if (replace) {
            it->second = donor;
        }
    }

    std::vector<std::pair<std::string, double>> selected;
    selected.reserve(best_by_bucket.size());
    for (const auto& [bucket_id, donor] : best_by_bucket) {
        (void)bucket_id;
        selected.emplace_back(donor.dna, donor.weight);
    }
    return selected;
}

bool exportBlockOrderLess(const ExportBlockOrderKey& lhs, const ExportBlockOrderKey& rhs) {
    return std::tie(lhs.ref_species, lhs.ref_chr, lhs.ref_start, lhs.block_id) <
           std::tie(rhs.ref_species, rhs.ref_chr, rhs.ref_start, rhs.block_id);
}


namespace {

struct SparseMatchingEdge {
    uint32_t first_vertex = 0;
    uint32_t second_vertex = 0;
    long double weight = 0.0L;
    uint32_t payload = 0;
};

struct SparseMatchingSelection {
    long double weight = 0.0L;
    std::vector<uint32_t> payloads;
};

bool sparseSelectionIsBetter(
    const SparseMatchingSelection& lhs,
    const SparseMatchingSelection& rhs) {
    if (lhs.payloads.size() !=
        rhs.payloads.size()) {
        return lhs.payloads.size() >
               rhs.payloads.size();
    }
    if (lhs.weight != rhs.weight) {
        return lhs.weight > rhs.weight;
    }
    return std::lexicographical_compare(
        lhs.payloads.begin(),
        lhs.payloads.end(),
        rhs.payloads.begin(),
        rhs.payloads.end());
}

struct BinaryLongDouble {
    unsigned __int128 significand = 0;
    int scale_exponent = 0;
};

BinaryLongDouble decomposeBinaryLongDouble(
    long double value) {
    if (value < 0.0L ||
        !std::isfinite(value)) {
        throw std::runtime_error(
            "Ancestral matching edge has invalid weight");
    }
    if (value == 0.0L) {
        return {};
    }
    static_assert(
        std::numeric_limits<long double>::radix == 2);
    static_assert(
        std::numeric_limits<long double>::digits <=
        128);
    int exponent = 0;
    const long double fraction =
        std::frexp(value, &exponent);
    const long double scaled_significand =
        std::ldexp(
            fraction,
            std::numeric_limits<long double>::
                digits);
    const auto significand =
        static_cast<unsigned __int128>(
            scaled_significand);
    if (significand == 0) {
        throw std::runtime_error(
            "Positive ancestral matching weight "
            "lost its binary significand");
    }
    return BinaryLongDouble{
        significand,
        exponent -
            std::numeric_limits<long double>::
                digits};
}

template<class ExactWeight>
ExactWeight exactScaledWeight(
    const BinaryLongDouble& binary,
    int minimum_scale_exponent) {
    if (binary.significand == 0) {
        return 0;
    }
    const int shift =
        binary.scale_exponent -
        minimum_scale_exponent;
    if (shift < 0) {
        throw std::runtime_error(
            "Ancestral matching exact weight has "
            "a negative scale shift");
    }
    ExactWeight result =
        static_cast<uint64_t>(
            binary.significand);
    const uint64_t high =
        static_cast<uint64_t>(
            binary.significand >> 64U);
    if (high != 0) {
        if constexpr (
            std::numeric_limits<
                ExactWeight>::digits <= 64) {
            throw std::overflow_error(
                "Ancestral matching significand exceeds "
                "the selected exact weight type");
        } else {
            result +=
                static_cast<ExactWeight>(high) <<
                64U;
        }
    }
    result <<= static_cast<unsigned>(shift);
    return result;
}

template<class ExactWeight>
SparseMatchingSelection solveExactBlossomMatching(
    const std::vector<uint32_t>& vertices,
    const std::vector<
        std::pair<uint32_t, uint32_t>>&
        component_edges,
    size_t component_begin,
    size_t component_end,
    const std::vector<SparseMatchingEdge>& edges,
    int minimum_scale_exponent) {
    using MatchingGraph =
        boost::adjacency_list<
            boost::vecS,
            boost::vecS,
            boost::undirectedS,
            boost::no_property,
            boost::property<
                boost::edge_weight_t,
                ExactWeight>>;
    MatchingGraph graph(vertices.size());
    std::unordered_map<uint32_t, uint32_t>
        local_index;
    local_index.reserve(vertices.size());
    for (uint32_t index = 0;
         index < vertices.size();
         ++index) {
        local_index.emplace(vertices[index], index);
    }

    std::vector<ExactWeight> exact_scores;
    exact_scores.reserve(
        component_end - component_begin);
    ExactWeight score_sum = 0;
    for (size_t index = component_begin;
         index < component_end;
         ++index) {
        const auto binary =
            decomposeBinaryLongDouble(
                edges[
                    component_edges[index]
                        .second]
                    .weight);
        auto score =
            exactScaledWeight<ExactWeight>(
                binary,
                minimum_scale_exponent);
        score_sum += score;
        exact_scores.push_back(
            std::move(score));
    }
    const ExactWeight cardinality_bonus =
        score_sum + 1;

    std::map<
        std::pair<uint32_t, uint32_t>,
        uint32_t>
        payload_by_vertices;
    for (size_t index = component_begin;
         index < component_end;
         ++index) {
        const auto& edge =
            edges[component_edges[index].second];
        uint32_t first =
            local_index.at(edge.first_vertex);
        uint32_t second =
            local_index.at(edge.second_vertex);
        if (second < first) {
            std::swap(first, second);
        }
        payload_by_vertices.emplace(
            std::pair{first, second},
            edge.payload);
        boost::add_edge(
            first,
            second,
            cardinality_bonus +
                exact_scores[
                    index - component_begin],
            graph);
    }

    using Vertex =
        boost::graph_traits<
            MatchingGraph>::vertex_descriptor;
    std::vector<Vertex> mate(
        boost::num_vertices(graph));
    auto mate_map =
        boost::make_iterator_property_map(
            mate.begin(),
            boost::get(
                boost::vertex_index,
                graph));
    boost::maximum_weighted_matching(
        graph,
        mate_map);
    const Vertex null_vertex =
        boost::graph_traits<
            MatchingGraph>::null_vertex();

    SparseMatchingSelection result;
    for (uint32_t local = 0;
         local < mate.size();
         ++local) {
        if (mate[local] == null_vertex ||
            local >=
                static_cast<uint32_t>(
                    mate[local])) {
            continue;
        }
        result.payloads.push_back(
            payload_by_vertices.at(
                {local,
                 static_cast<uint32_t>(
                     mate[local])}));
    }
    std::sort(
        result.payloads.begin(),
        result.payloads.end());
    return result;
}

SparseMatchingSelection solveTreeMaximumWeightMatching(
    const std::vector<uint32_t>& vertices,
    const std::vector<
        std::pair<uint32_t, uint32_t>>&
        component_edges,
    size_t component_begin,
    size_t component_end,
    const std::vector<SparseMatchingEdge>& edges) {
    const uint32_t no_vertex =
        std::numeric_limits<uint32_t>::max();
    std::unordered_map<uint32_t, uint32_t>
        local_index;
    local_index.reserve(vertices.size());
    for (uint32_t index = 0;
         index < vertices.size();
         ++index) {
        local_index.emplace(vertices[index], index);
    }

    std::vector<
        std::vector<
            std::pair<uint32_t, uint32_t>>>
        adjacency(vertices.size());
    for (size_t index = component_begin;
         index < component_end;
         ++index) {
        const uint32_t edge_index =
            component_edges[index].second;
        const auto& edge = edges[edge_index];
        const uint32_t first =
            local_index.at(edge.first_vertex);
        const uint32_t second =
            local_index.at(edge.second_vertex);
        adjacency[first].emplace_back(
            second,
            edge_index);
        adjacency[second].emplace_back(
            first,
            edge_index);
    }

    std::vector<uint32_t> parent(
        vertices.size(),
        no_vertex);
    std::vector<uint32_t> parent_edge(
        vertices.size(),
        no_vertex);
    std::vector<uint32_t> order;
    order.reserve(vertices.size());
    parent[0] = 0;
    order.push_back(0);
    for (size_t cursor = 0;
         cursor < order.size();
         ++cursor) {
        const uint32_t vertex = order[cursor];
        for (const auto& [neighbor, edge_index] :
             adjacency[vertex]) {
            if (parent[neighbor] != no_vertex) {
                continue;
            }
            parent[neighbor] = vertex;
            parent_edge[neighbor] = edge_index;
            order.push_back(neighbor);
        }
    }
    if (order.size() != vertices.size()) {
        throw std::runtime_error(
            "Sparse tree matching component is disconnected");
    }

    struct TreeObjective {
        uint64_t cardinality = 0;
        long double weight = 0.0L;
    };
    const auto objective_is_better =
        [](const TreeObjective& lhs,
           const TreeObjective& rhs) {
            return lhs.cardinality != rhs.cardinality
                ? lhs.cardinality > rhs.cardinality
                : lhs.weight > rhs.weight;
        };

    std::vector<TreeObjective>
        free_at_parent(vertices.size());
    std::vector<TreeObjective>
        matched_to_parent(vertices.size());
    std::vector<uint32_t> chosen_child(
        vertices.size(),
        no_vertex);
    for (auto order_it = order.rbegin();
         order_it != order.rend();
         ++order_it) {
        const uint32_t vertex = *order_it;
        TreeObjective baseline;
        for (const auto& [neighbor, edge_index] :
             adjacency[vertex]) {
            (void)edge_index;
            if (parent[neighbor] != vertex) {
                continue;
            }
            baseline.cardinality +=
                free_at_parent[neighbor]
                    .cardinality;
            baseline.weight +=
                free_at_parent[neighbor]
                    .weight;
        }
        matched_to_parent[vertex] = baseline;
        free_at_parent[vertex] = baseline;
        for (const auto& [child, edge_index] :
             adjacency[vertex]) {
            if (parent[child] != vertex) {
                continue;
            }
            if (baseline.cardinality <
                free_at_parent[child]
                    .cardinality) {
                throw std::runtime_error(
                    "Sparse tree matching cardinality underflow");
            }
            TreeObjective candidate{
                baseline.cardinality -
                    free_at_parent[child]
                        .cardinality +
                    matched_to_parent[child]
                        .cardinality +
                    1,
                baseline.weight -
                    free_at_parent[child]
                        .weight +
                    matched_to_parent[child]
                        .weight +
                    edges[edge_index].weight};
            if (objective_is_better(
                    candidate,
                    free_at_parent[vertex])) {
                free_at_parent[vertex] =
                    candidate;
                chosen_child[vertex] =
                    child;
            }
        }
    }

    SparseMatchingSelection result;
    result.weight =
        free_at_parent.front().weight;
    std::vector<std::pair<uint32_t, bool>>
        reconstruction{{0, false}};
    while (!reconstruction.empty()) {
        const auto [vertex, matched_by_parent] =
            reconstruction.back();
        reconstruction.pop_back();
        const uint32_t selected_child =
            matched_by_parent
            ? no_vertex
            : chosen_child[vertex];
        for (const auto& [child, edge_index] :
             adjacency[vertex]) {
            if (parent[child] != vertex) {
                continue;
            }
            const bool child_matched =
                child == selected_child;
            if (child_matched) {
                result.payloads.push_back(
                    edges[edge_index].payload);
            }
            reconstruction.emplace_back(
                child,
                child_matched);
        }
    }
    if (result.payloads.size() !=
        free_at_parent.front().cardinality) {
        throw std::runtime_error(
            "Sparse tree matching reconstruction "
            "cardinality mismatch");
    }
    std::sort(
        result.payloads.begin(),
        result.payloads.end());
    return result;
}

std::vector<uint32_t> solveSparseMaximumWeightMatching(
    uint32_t vertex_count,
    const std::vector<SparseMatchingEdge>& edges) {
    if (edges.empty()) {
        return {};
    }
    constexpr uint32_t no_edge =
        std::numeric_limits<uint32_t>::max();
    std::vector<uint32_t> parent(vertex_count);
    std::iota(parent.begin(), parent.end(), 0);
    auto find_root = [&](uint32_t value) {
        uint32_t root = value;
        while (parent[root] != root) {
            root = parent[root];
        }
        while (parent[value] != value) {
            const uint32_t next = parent[value];
            parent[value] = root;
            value = next;
        }
        return root;
    };
    auto unite = [&](uint32_t lhs, uint32_t rhs) {
        lhs = find_root(lhs);
        rhs = find_root(rhs);
        if (lhs == rhs) {
            return;
        }
        if (lhs > rhs) {
            std::swap(lhs, rhs);
        }
        parent[rhs] = lhs;
    };

    std::vector<uint32_t> degree(vertex_count, 0);
    std::vector<std::array<uint32_t, 2>> incident(
        vertex_count,
        std::array<uint32_t, 2>{no_edge, no_edge});
    for (uint32_t edge_index = 0;
         edge_index < edges.size();
         ++edge_index) {
        const auto& edge = edges[edge_index];
        if (edge.first_vertex >= vertex_count ||
            edge.second_vertex >= vertex_count ||
            edge.first_vertex == edge.second_vertex) {
            throw std::invalid_argument(
                "Sparse matching edge has invalid vertices");
        }
        unite(edge.first_vertex, edge.second_vertex);
        for (uint32_t vertex :
             {edge.first_vertex, edge.second_vertex}) {
            if (degree[vertex] < 2) {
                incident[vertex][degree[vertex]] =
                    edge_index;
            }
            ++degree[vertex];
        }
    }

    std::vector<std::pair<uint32_t, uint32_t>>
        component_edges;
    component_edges.reserve(edges.size());
    for (uint32_t edge_index = 0;
         edge_index < edges.size();
         ++edge_index) {
        component_edges.emplace_back(
            find_root(edges[edge_index].first_vertex),
            edge_index);
    }
    std::sort(
        component_edges.begin(),
        component_edges.end());

    auto solve_linear =
        [&](const std::vector<uint32_t>& ordered_edges,
            size_t begin,
            size_t end) {
            const size_t count = end - begin;
            std::vector<uint32_t> best_cardinality(
                count + 1,
                0);
            std::vector<long double> best_weight(
                count + 1,
                0.0L);
            std::vector<uint8_t> take(
                count + 1,
                0);
            for (size_t prefix = 1;
                 prefix <= count;
                 ++prefix) {
                const uint32_t include_cardinality =
                    1U +
                    (prefix > 1
                         ? best_cardinality[prefix - 2]
                         : 0U);
                const long double include_weight =
                    edges[ordered_edges[
                              begin + prefix - 1]]
                        .weight +
                    (prefix > 1
                         ? best_weight[prefix - 2]
                         : 0.0L);
                const uint32_t exclude_cardinality =
                    best_cardinality[prefix - 1];
                const long double exclude_weight =
                    best_weight[prefix - 1];
                if (include_cardinality >
                        exclude_cardinality ||
                    (include_cardinality ==
                         exclude_cardinality &&
                     include_weight >
                         exclude_weight)) {
                    best_cardinality[prefix] =
                        include_cardinality;
                    best_weight[prefix] =
                        include_weight;
                    take[prefix] = 1;
                } else {
                    best_cardinality[prefix] =
                        exclude_cardinality;
                    best_weight[prefix] =
                        exclude_weight;
                }
            }
            SparseMatchingSelection selection;
            selection.weight = best_weight[count];
            size_t prefix = count;
            while (prefix != 0) {
                if (!take[prefix]) {
                    --prefix;
                    continue;
                }
                selection.payloads.push_back(
                    edges[ordered_edges[
                              begin + prefix - 1]]
                        .payload);
                prefix = prefix > 1
                             ? prefix - 2
                             : 0;
            }
            std::sort(
                selection.payloads.begin(),
                selection.payloads.end());
            return selection;
        };

    std::vector<uint32_t> selected_payloads;
    size_t singleton_components = 0;
    size_t degree_two_components = 0;
    size_t tree_components = 0;
    size_t blossom_components = 0;
    size_t maximum_blossom_vertices = 0;
    size_t component_begin = 0;
    while (component_begin <
           component_edges.size()) {
        size_t component_end =
            component_begin + 1;
        while (component_end <
                   component_edges.size() &&
               component_edges[component_end].first ==
                   component_edges[component_begin].first) {
            ++component_end;
        }
        if (component_end - component_begin == 1) {
            ++singleton_components;
            selected_payloads.push_back(
                edges[component_edges[
                          component_begin]
                          .second]
                    .payload);
            component_begin = component_end;
            continue;
        }

        bool degree_at_most_two = true;
        uint32_t start_vertex =
            edges[component_edges[
                      component_begin]
                      .second]
                .first_vertex;
        for (size_t index = component_begin;
             index < component_end;
             ++index) {
            const auto& edge =
                edges[component_edges[index].second];
            for (uint32_t vertex :
                 {edge.first_vertex,
                  edge.second_vertex}) {
                if (degree[vertex] > 2) {
                    degree_at_most_two = false;
                } else if (
                    degree[vertex] == 1 &&
                    (degree[start_vertex] != 1 ||
                     vertex < start_vertex)) {
                    start_vertex = vertex;
                }
            }
        }

        SparseMatchingSelection selection;
        if (degree_at_most_two) {
            ++degree_two_components;
            const bool is_cycle =
                degree[start_vertex] == 2;
            std::vector<uint32_t> ordered_edges;
            ordered_edges.reserve(
                component_end - component_begin);
            uint32_t previous_edge = no_edge;
            uint32_t current_vertex = start_vertex;
            while (ordered_edges.size() <
                   component_end - component_begin) {
                uint32_t next_edge =
                    incident[current_vertex][0] !=
                            previous_edge
                        ? incident[current_vertex][0]
                        : incident[current_vertex][1];
                if (next_edge == no_edge) {
                    break;
                }
                ordered_edges.push_back(next_edge);
                const auto& edge = edges[next_edge];
                current_vertex =
                    edge.first_vertex == current_vertex
                        ? edge.second_vertex
                        : edge.first_vertex;
                previous_edge = next_edge;
            }
            if (ordered_edges.size() !=
                component_end - component_begin) {
                throw std::runtime_error(
                    "Sparse degree-two matching component traversal failed");
            }
            if (!is_cycle) {
                selection = solve_linear(
                    ordered_edges,
                    0,
                    ordered_edges.size());
            } else {
                SparseMatchingSelection exclude_last =
                    solve_linear(
                        ordered_edges,
                        0,
                        ordered_edges.size() - 1);
                SparseMatchingSelection include_last =
                    solve_linear(
                        ordered_edges,
                        1,
                        ordered_edges.size() - 2);
                include_last.weight +=
                    edges[ordered_edges.back()].weight;
                include_last.payloads.push_back(
                    edges[ordered_edges.back()].payload);
                std::sort(
                    include_last.payloads.begin(),
                    include_last.payloads.end());
                selection = sparseSelectionIsBetter(
                                include_last,
                                exclude_last)
                                ? std::move(include_last)
                                : std::move(exclude_last);
            }
        } else {
            std::vector<uint32_t> vertices;
            vertices.reserve(
                (component_end - component_begin) *
                2);
            for (size_t index = component_begin;
                 index < component_end;
                 ++index) {
                const auto& edge =
                    edges[component_edges[index].second];
                vertices.push_back(edge.first_vertex);
                vertices.push_back(edge.second_vertex);
            }
            std::sort(vertices.begin(), vertices.end());
            vertices.erase(
                std::unique(
                    vertices.begin(),
                    vertices.end()),
                vertices.end());
            if (component_end -
                        component_begin +
                    1 ==
                vertices.size()) {
                ++tree_components;
                const auto tree_begin =
                    std::chrono::steady_clock::now();
                if (vertices.size() >= 64) {
                    spdlog::info(
                        "HAL ancestral sparse matching started tree "
                        "component with {} vertices and {} edges",
                        vertices.size(),
                        component_end - component_begin);
                }
                selection =
                    solveTreeMaximumWeightMatching(
                        vertices,
                        component_edges,
                        component_begin,
                        component_end,
                        edges);
                if (vertices.size() >= 64) {
                    spdlog::info(
                        "HAL ancestral sparse matching finished tree "
                        "component with {} vertices in {} ms",
                        vertices.size(),
                        std::chrono::duration_cast<
                            std::chrono::milliseconds>(
                            std::chrono::steady_clock::now() -
                            tree_begin)
                            .count());
                }
            } else {
                const auto blossom_begin =
                    std::chrono::steady_clock::now();
                if (vertices.size() >= 64) {
                    spdlog::info(
                        "HAL ancestral sparse matching started exact "
                        "blossom component with {} vertices and {} edges",
                        vertices.size(),
                        component_end - component_begin);
                }
                ++blossom_components;
                maximum_blossom_vertices =
                    std::max(
                        maximum_blossom_vertices,
                        vertices.size());
                constexpr size_t
                    maximum_blossom_component_vertices =
                        2048;
                if (vertices.size() >
                    maximum_blossom_component_vertices) {
                    throw std::runtime_error(
                        "Ancestral matching component exceeds the exact sparse-topology limit");
                }

                int minimum_scale_exponent =
                    std::numeric_limits<int>::max();
                int maximum_scale_exponent =
                    std::numeric_limits<int>::min();
                bool have_positive_weight = false;
                for (size_t index = component_begin;
                     index < component_end;
                     ++index) {
                    const auto binary =
                        decomposeBinaryLongDouble(
                            edges[
                                component_edges[index]
                                    .second]
                                .weight);
                    if (binary.significand == 0) {
                        continue;
                    }
                    have_positive_weight = true;
                    minimum_scale_exponent =
                        std::min(
                            minimum_scale_exponent,
                            binary.scale_exponent);
                    maximum_scale_exponent =
                        std::max(
                            maximum_scale_exponent,
                            binary.scale_exponent);
                }
                if (!have_positive_weight) {
                    minimum_scale_exponent = 0;
                    maximum_scale_exponent = 0;
                }
                size_t edge_count_bits = 0;
                for (size_t count =
                         component_end -
                         component_begin;
                     count != 0;
                     count >>= 1U) {
                    ++edge_count_bits;
                }
                const size_t required_bits =
                    static_cast<size_t>(
                        std::numeric_limits<
                            long double>::digits) +
                    static_cast<size_t>(
                        maximum_scale_exponent -
                        minimum_scale_exponent) +
                    edge_count_bits +
                    8;
                if (required_bits <= 62) {
                    selection =
                        solveExactBlossomMatching<
                            int64_t>(
                            vertices,
                            component_edges,
                            component_begin,
                            component_end,
                            edges,
                            minimum_scale_exponent);
                } else if (required_bits <= 126) {
                    selection =
                        solveExactBlossomMatching<
                            __int128_t>(
                            vertices,
                            component_edges,
                            component_begin,
                            component_end,
                            edges,
                            minimum_scale_exponent);
                } else {
                    constexpr size_t
                        wide_weight_bits = 32768;
                    if (required_bits >
                        wide_weight_bits - 8) {
                        throw std::overflow_error(
                            "Ancestral matching exact weight exceeds "
                            "the supported binary exponent range");
                    }
                    using WideExactWeight =
                        boost::multiprecision::number<
                            boost::multiprecision::
                                cpp_int_backend<
                                    wide_weight_bits,
                                    wide_weight_bits,
                                    boost::multiprecision::
                                        signed_magnitude,
                                    boost::multiprecision::
                                        unchecked,
                                    void>>;
                    selection =
                        solveExactBlossomMatching<
                            WideExactWeight>(
                            vertices,
                            component_edges,
                            component_begin,
                            component_end,
                            edges,
                            minimum_scale_exponent);
                }
                if (vertices.size() >= 64) {
                    spdlog::info(
                        "HAL ancestral sparse matching finished exact "
                        "blossom component with {} vertices in {} ms",
                        vertices.size(),
                        std::chrono::duration_cast<
                            std::chrono::milliseconds>(
                            std::chrono::steady_clock::now() -
                            blossom_begin)
                            .count());
                }
            }
        }
        selected_payloads.insert(
            selected_payloads.end(),
            selection.payloads.begin(),
            selection.payloads.end());
        component_begin = component_end;
    }
    std::sort(
        selected_payloads.begin(),
        selected_payloads.end());
    spdlog::info(
        "HAL ancestral sparse matching solved {} edges in {} "
        "singleton, {} degree-two, {} tree, and {} blossom "
        "components (largest blossom {} vertices)",
        edges.size(),
        singleton_components,
        degree_two_components,
        tree_components,
        blossom_components,
        maximum_blossom_vertices);
    return selected_payloads;
}

}  // namespace

PathCoverResult buildMaximumCardinalityWeightPathCoverDetailedImpl(
    const std::vector<uint64_t>& run_ids,
    const std::vector<EdgeSupport>& edges,
    const std::unordered_map<uint64_t, RunOrderKey>& run_order_keys,
    ExportStats* stats) {

    PathCoverResult result;
    if (run_ids.empty()) {
        return result;
    }

    auto run_less = [&](uint64_t lhs, uint64_t rhs) {
        const auto lhs_it = run_order_keys.find(lhs);
        const auto rhs_it = run_order_keys.find(rhs);
        if (lhs_it == run_order_keys.end() || rhs_it == run_order_keys.end()) {
            return lhs < rhs;
        }
        return std::tie(
                   lhs_it->second.block_id,
                   lhs_it->second.col_beg,
                   lhs) <
               std::tie(
                   rhs_it->second.block_id,
                   rhs_it->second.col_beg,
                   rhs);
    };

    std::vector<uint64_t> sorted_runs = run_ids;
    std::sort(sorted_runs.begin(), sorted_runs.end(), run_less);
    if (std::adjacent_find(sorted_runs.begin(), sorted_runs.end()) != sorted_runs.end()) {
        throw std::runtime_error("Duplicate run id in maximum-weight path cover");
    }
    if (sorted_runs.size() >
        static_cast<size_t>(std::numeric_limits<uint32_t>::max() / 2U)) {
        throw std::overflow_error("Too many runs for ancestral extremity matching");
    }

    const uint32_t run_count = static_cast<uint32_t>(sorted_runs.size());
    const uint32_t endpoint_count = run_count * 2U;
    std::unordered_map<uint64_t, uint32_t> run_index;
    run_index.reserve(sorted_runs.size());
    for (uint32_t index = 0; index < run_count; ++index) {
        run_index.emplace(sorted_runs[index], index);
    }

    struct CandidateEdge {
        uint32_t first_endpoint = 0;
        uint32_t second_endpoint = 0;
        uint32_t child_support = 0;
        uint64_t occurrence_support = 0;
        long double weighted_support = 0.0L;
    };

    std::vector<CandidateEdge> candidates;
    candidates.reserve(edges.size());
    for (const auto& edge : edges) {
        const auto from_it = run_index.find(edge.from);
        const auto to_it = run_index.find(edge.to);
        if (from_it == run_index.end() || to_it == run_index.end()) {
            continue;
        }
        if (edge.from == edge.to) {
            throw std::runtime_error(
                "Self adjacency in ancestral thread evidence for run " +
                std::to_string(edge.from));
        }
        if (edge.supporting_children.empty() && edge.occurrence_support == 0) {
            continue;
        }

        uint32_t first_endpoint =
            from_it->second * 2U +
            (edge.from_forward_to_canonical ? 1U : 0U);
        uint32_t second_endpoint =
            to_it->second * 2U +
            (edge.to_forward_to_canonical ? 0U : 1U);
        if (second_endpoint < first_endpoint) {
            std::swap(first_endpoint, second_endpoint);
        }
        candidates.push_back(CandidateEdge{
            first_endpoint,
            second_endpoint,
            static_cast<uint32_t>(edge.supporting_children.size()),
            edge.occurrence_support,
            edge.weighted_support});
    }
    std::sort(candidates.begin(), candidates.end(), [](const auto& lhs, const auto& rhs) {
        return std::tie(lhs.first_endpoint, lhs.second_endpoint) <
               std::tie(rhs.first_endpoint, rhs.second_endpoint);
    });

    std::vector<CandidateEdge> unique_candidates;
    unique_candidates.reserve(candidates.size());
    for (const auto& edge : candidates) {
        if (!unique_candidates.empty() &&
            unique_candidates.back().first_endpoint == edge.first_endpoint &&
            unique_candidates.back().second_endpoint == edge.second_endpoint) {
            auto& aggregate = unique_candidates.back();
            aggregate.child_support =
                std::max(aggregate.child_support, edge.child_support);
            aggregate.occurrence_support += edge.occurrence_support;
            aggregate.weighted_support += edge.weighted_support;
        } else {
            unique_candidates.push_back(edge);
        }
    }
    auto candidate_score = [](const CandidateEdge& edge) {
        return edge.weighted_support > 0.0L
                   ? edge.weighted_support
                   : static_cast<long double>(edge.child_support) *
                             1'000'000.0L +
                         static_cast<long double>(
                             edge.occurrence_support);
    };


    std::vector<SparseMatchingEdge> matching_edges;
    matching_edges.reserve(unique_candidates.size());
    for (uint32_t edge_index = 0;
         edge_index < unique_candidates.size();
         ++edge_index) {
        const auto& edge =
            unique_candidates[edge_index];
        matching_edges.push_back(
            SparseMatchingEdge{
                edge.first_endpoint,
                edge.second_endpoint,
                candidate_score(edge),
                edge_index});
    }
    std::vector<uint32_t> matched_edges =
        solveSparseMaximumWeightMatching(
            endpoint_count,
            matching_edges);
    std::vector<int64_t> edge_by_endpoint(
        endpoint_count,
        -1);
    std::sort(
        matched_edges.begin(),
        matched_edges.end(),
        [&](uint32_t lhs, uint32_t rhs) {
            const long double lhs_score =
                candidate_score(
                    unique_candidates[lhs]);
            const long double rhs_score =
                candidate_score(
                    unique_candidates[rhs]);
            if (lhs_score != rhs_score) {
                return lhs_score > rhs_score;
            }
            return std::tie(
                       unique_candidates[lhs]
                           .first_endpoint,
                       unique_candidates[lhs]
                           .second_endpoint) <
                   std::tie(
                       unique_candidates[rhs]
                           .first_endpoint,
                       unique_candidates[rhs]
                           .second_endpoint);
        });

    std::vector<uint32_t> chain_parent(run_count);
    std::iota(
        chain_parent.begin(),
        chain_parent.end(),
        0);
    auto find_chain_root = [&](uint32_t value) {
        uint32_t root = value;
        while (chain_parent[root] != root) {
            root = chain_parent[root];
        }
        while (chain_parent[value] != value) {
            const uint32_t next =
                chain_parent[value];
            chain_parent[value] = root;
            value = next;
        }
        return root;
    };
    for (uint32_t edge_index : matched_edges) {
        const auto& edge =
            unique_candidates[edge_index];
        uint32_t first_root = find_chain_root(
            edge.first_endpoint / 2U);
        uint32_t second_root = find_chain_root(
            edge.second_endpoint / 2U);
        if (first_root == second_root) {
            continue;
        }
        if (edge_by_endpoint[
                edge.first_endpoint] >= 0 ||
            edge_by_endpoint[
                edge.second_endpoint] >= 0) {
            throw std::runtime_error(
                "Maximum-weight adjacency matching emitted duplicate extremities");
        }
        if (second_root < first_root) {
            std::swap(first_root, second_root);
        }
        chain_parent[second_root] = first_root;
        edge_by_endpoint[edge.first_endpoint] =
            static_cast<int64_t>(edge_index);
        edge_by_endpoint[edge.second_endpoint] =
            static_cast<int64_t>(edge_index);
    }

    std::vector<uint32_t> start_endpoints;
    start_endpoints.reserve(run_count * 2U);
    for (uint32_t endpoint = 0; endpoint < endpoint_count; ++endpoint) {
        if (edge_by_endpoint[endpoint] < 0) {
            start_endpoints.push_back(endpoint);
        }
    }

    std::vector<uint8_t> visited(run_count, 0);
    result.forward_by_occurrence.reserve(run_count);
    for (uint32_t start_endpoint : start_endpoints) {
        if (visited[start_endpoint / 2U]) {
            continue;
        }
        std::vector<uint64_t> path;
        uint32_t entry_endpoint = start_endpoint;
        while (true) {
            const uint32_t run = entry_endpoint / 2U;
            if (visited[run]) {
                throw std::runtime_error(
                    "Ancestral extremity path cover contains a residual cycle");
            }
            visited[run] = 1;
            const uint64_t run_id = sorted_runs[run];
            path.push_back(run_id);
            result.forward_by_occurrence.emplace(
                run_id,
                (entry_endpoint & 1U) == 0);

            const uint32_t exit_endpoint = entry_endpoint ^ 1U;
            const int64_t edge_index = edge_by_endpoint[exit_endpoint];
            if (edge_index < 0) {
                break;
            }
            const auto& edge =
                unique_candidates[static_cast<size_t>(edge_index)];
            entry_endpoint =
                edge.first_endpoint == exit_endpoint
                    ? edge.second_endpoint
                    : edge.first_endpoint;
        }
        result.paths.push_back(std::move(path));
    }

    if (std::find(visited.begin(), visited.end(), 0) != visited.end()) {
        throw std::runtime_error(
            "Ancestral extremity path cover left an unreachable run");
    }
    std::sort(
        result.paths.begin(),
        result.paths.end(),
        [&](const auto& lhs, const auto& rhs) {
            return run_less(lhs.front(), rhs.front());
        });
    return result;
}


std::vector<std::vector<uint64_t>> buildMaximumCardinalityWeightPathCover(
    const std::vector<uint64_t>& run_ids,
    const std::vector<EdgeSupport>& edges,
    const std::unordered_map<uint64_t, RunOrderKey>& run_order_keys,
    ExportStats* stats) {
    return buildMaximumCardinalityWeightPathCoverDetailedImpl(
        run_ids, edges, run_order_keys, stats).paths;
}

AncestralSequenceAssembly buildAncestralSequenceAssembly(
    const std::vector<uint64_t>& occurrence_ids,
    const std::vector<EdgeSupport>& edges,
    const std::unordered_map<uint64_t, RunOrderKey>&
        occurrence_order_keys,
    const std::vector<TerminalEndSupport>& terminal_ends,
    uint32_t scaffold_gap_length,
    ExportStats* stats) {
    if (scaffold_gap_length == 0) {
        throw std::invalid_argument(
            "Reference indirect gaps must contain at least one base");
    }

    using EndKey = std::pair<OccurrenceId, uint8_t>;
    using EndPairKey = std::pair<EndKey, EndKey>;
    auto canonical_end_pair = [](EndKey first,
                                 EndKey second) {
        if (second < first) {
            std::swap(first, second);
        }
        return EndPairKey{first, second};
    };
    auto edge_ends = [&](const EdgeSupport& edge) {
        return canonical_end_pair(
            EndKey{
                edge.from,
                static_cast<uint8_t>(
                    edge.from_forward_to_canonical
                        ? OccurrenceEndSide::RIGHT
                        : OccurrenceEndSide::LEFT)},
            EndKey{
                edge.to,
                static_cast<uint8_t>(
                    edge.to_forward_to_canonical
                        ? OccurrenceEndSide::LEFT
                        : OccurrenceEndSide::RIGHT)});
    };

    std::unordered_set<OccurrenceId> occurrence_set(
        occurrence_ids.begin(), occurrence_ids.end());
    std::map<EndKey, std::set<int>>
        terminal_children;
    for (const auto& terminal : terminal_ends) {
        if (!occurrence_set.contains(
                terminal.end.occurrence_id)) {
            throw std::invalid_argument(
                "Terminal evidence references an unknown ancestor occurrence");
        }
        auto& children = terminal_children[EndKey{
            terminal.end.occurrence_id,
            static_cast<uint8_t>(terminal.end.side)}];
        children.insert(
            terminal.supporting_lineages.begin(),
            terminal.supporting_lineages.end());
    }
    std::set<EndKey> confirmed_terminal_ends;
    for (const auto& [end, children] :
         terminal_children) {
        if (children.size() >= 2) {
            confirmed_terminal_ends.insert(end);
        }
    }

    std::vector<EdgeSupport> admissible_edges;
    admissible_edges.reserve(edges.size());
    std::map<EndPairKey, const EdgeSupport*>
        evidence_by_ends;
    auto evidence_is_better =
        [](const EdgeSupport& lhs,
           const EdgeSupport& rhs) {
            if (lhs.weighted_support !=
                rhs.weighted_support) {
                return lhs.weighted_support >
                       rhs.weighted_support;
            }
            if (lhs.supporting_children.size() !=
                rhs.supporting_children.size()) {
                return lhs.supporting_children.size() >
                       rhs.supporting_children.size();
            }
            if (lhs.occurrence_support !=
                rhs.occurrence_support) {
                return lhs.occurrence_support >
                       rhs.occurrence_support;
            }
            return lhs.minimum_gap <
                   rhs.minimum_gap;
        };
    for (const auto& edge : edges) {
        const auto ends = edge_ends(edge);
        if (confirmed_terminal_ends.contains(
                ends.first) ||
            confirmed_terminal_ends.contains(
                ends.second)) {
            continue;
        }
        admissible_edges.push_back(edge);
        const EdgeSupport* stored =
            &admissible_edges.back();
        auto [it, inserted] =
            evidence_by_ends.emplace(ends, stored);
        if (!inserted &&
            evidence_is_better(*stored, *it->second)) {
            it->second = stored;
        }
    }

    PathCoverResult paths =
        buildMaximumCardinalityWeightPathCoverDetailedImpl(
            occurrence_ids,
            admissible_edges,
            occurrence_order_keys,
            nullptr);
    AncestralSequenceAssembly assembly;
    assembly.forward_by_occurrence =
        std::move(paths.forward_by_occurrence);

    uint64_t direct_join_count = 0;
    uint64_t indirect_join_count = 0;
    uint64_t scaffold_join_count = 0;
    uint64_t gap_bases = 0;
    for (const auto& path : paths.paths) {
        if (path.empty()) {
            throw std::runtime_error(
                "Evidence-constrained path cover emitted an empty path");
        }
        AncestralSequencePath sequence;
        sequence.supported_fragments.reserve(
            path.size());
        sequence.joins.reserve(path.size() - 1);
        for (size_t index = 0;
             index < path.size();
             ++index) {
            const OccurrenceId occurrence_id =
                path[index];
            sequence.supported_fragments.push_back(
                {occurrence_id});
            if (index == 0) {
                continue;
            }
            const OccurrenceId previous =
                path[index - 1];
            const bool previous_forward =
                assembly.forward_by_occurrence.at(
                    previous);
            const bool current_forward =
                assembly.forward_by_occurrence.at(
                    occurrence_id);
            const EndPairKey ends =
                canonical_end_pair(
                    EndKey{
                        previous,
                        static_cast<uint8_t>(
                            previous_forward
                                ? OccurrenceEndSide::RIGHT
                                : OccurrenceEndSide::LEFT)},
                    EndKey{
                        occurrence_id,
                        static_cast<uint8_t>(
                            current_forward
                                ? OccurrenceEndSide::LEFT
                                : OccurrenceEndSide::RIGHT)});
            const auto edge_it =
                evidence_by_ends.find(ends);
            if (edge_it == evidence_by_ends.end()) {
                throw std::runtime_error(
                    "Evidence-constrained path cover selected an unsupported join");
            }
            const EdgeSupport& edge =
                *edge_it->second;
            if (edge.minimum_gap == 0) {
                sequence.joins.push_back(
                    ReferenceJoin{
                        ReferenceJoinKind::DIRECT,
                        0,
                        edge.weighted_support});
                ++direct_join_count;
            } else {
                sequence.joins.push_back(
                    ReferenceJoin{
                        ReferenceJoinKind::INDIRECT,
                        scaffold_gap_length,
                        edge.weighted_support});
                ++indirect_join_count;
                gap_bases += scaffold_gap_length;
            }
        }
        assembly.sequences.push_back(
            std::move(sequence));
    }

    std::vector<AncestralSequencePath> reference_paths =
        std::move(assembly.sequences);
    for (const auto& sequence : reference_paths) {
        if (sequence.supported_fragments.empty() ||
            sequence.supported_fragments.front().empty()) {
            throw std::runtime_error(
                "Evidence-constrained path cover emitted an empty reference interval");
        }
    }
    assembly.sequences.clear();
    assembly.sequences.reserve(paths.paths.size());

    auto endpoint =
        [&](const AncestralSequencePath& sequence,
            bool left) {
            const auto& fragment =
                left
                    ? sequence.supported_fragments
                          .front()
                    : sequence.supported_fragments
                          .back();
            const OccurrenceId occurrence_id =
                left ? fragment.front()
                     : fragment.back();
            const bool forward =
                assembly.forward_by_occurrence.at(
                    occurrence_id);
            return EndKey{
                occurrence_id,
                static_cast<uint8_t>(
                    left
                        ? (forward
                               ? OccurrenceEndSide::LEFT
                               : OccurrenceEndSide::RIGHT)
                        : (forward
                               ? OccurrenceEndSide::RIGHT
                               : OccurrenceEndSide::LEFT))};
        };
    auto reverse_path =
        [&](AncestralSequencePath& sequence) {
            std::reverse(
                sequence.supported_fragments.begin(),
                sequence.supported_fragments.end());
            std::reverse(
                sequence.joins.begin(),
                sequence.joins.end());
            for (const auto& fragment :
                 sequence.supported_fragments) {
                for (OccurrenceId occurrence_id :
                     fragment) {
                    auto& forward =
                        assembly
                            .forward_by_occurrence
                            .at(occurrence_id);
                    forward = !forward;
                }
            }
        };
    auto orient_canonically =
        [&](AncestralSequencePath& sequence) {
            if (endpoint(sequence, false) <
                endpoint(sequence, true)) {
                reverse_path(sequence);
            }
        };
    auto append_scaffold =
        [&](AncestralSequencePath& merged,
            AncestralSequencePath&& next) {
            if (confirmed_terminal_ends.contains(
                    endpoint(merged, false)) ||
                confirmed_terminal_ends.contains(
                    endpoint(next, true))) {
                throw std::logic_error(
                    "Reference scaffold attempted to consume a confirmed terminal end");
            }
            merged.joins.push_back(
                ReferenceJoin{
                    ReferenceJoinKind::SCAFFOLD,
                    scaffold_gap_length,
                    0.0L});
            merged.supported_fragments.insert(
                merged.supported_fragments.end(),
                std::make_move_iterator(
                    next.supported_fragments
                        .begin()),
                std::make_move_iterator(
                    next.supported_fragments
                        .end()));
            merged.joins.insert(
                merged.joins.end(),
                std::make_move_iterator(
                    next.joins.begin()),
                std::make_move_iterator(
                    next.joins.end()));
            ++scaffold_join_count;
            gap_bases += scaffold_gap_length;
        };
    auto stable_path_less =
        [&](const AncestralSequencePath& lhs,
            const AncestralSequencePath& rhs) {
            return std::tuple{
                       endpoint(lhs, true),
                       endpoint(lhs, false)} <
                   std::tuple{
                       endpoint(rhs, true),
                       endpoint(rhs, false)};
        };

    {
        auto& sequences = reference_paths;
        std::vector<AncestralSequencePath>
            closed_paths;
        std::vector<AncestralSequencePath>
            single_terminal_paths;
        std::vector<AncestralSequencePath>
            free_paths;
        closed_paths.reserve(sequences.size());
        single_terminal_paths.reserve(
            sequences.size());
        free_paths.reserve(sequences.size());
        for (auto& sequence : sequences) {
            const bool left_terminal =
                confirmed_terminal_ends.contains(
                    endpoint(sequence, true));
            const bool right_terminal =
                confirmed_terminal_ends.contains(
                    endpoint(sequence, false));
            if (left_terminal &&
                right_terminal) {
                orient_canonically(sequence);
                closed_paths.push_back(
                    std::move(sequence));
            } else if (left_terminal ||
                       right_terminal) {
                single_terminal_paths.push_back(
                    std::move(sequence));
            } else {
                orient_canonically(sequence);
                free_paths.push_back(
                    std::move(sequence));
            }
        }
        std::sort(
            closed_paths.begin(),
            closed_paths.end(),
            stable_path_less);
        std::sort(
            free_paths.begin(),
            free_paths.end(),
            stable_path_less);
        std::sort(
            single_terminal_paths.begin(),
            single_terminal_paths.end(),
            [&](const auto& lhs,
                const auto& rhs) {
                const EndKey lhs_left =
                    endpoint(lhs, true);
                const EndKey lhs_right =
                    endpoint(lhs, false);
                const EndKey rhs_left =
                    endpoint(rhs, true);
                const EndKey rhs_right =
                    endpoint(rhs, false);
                const EndKey lhs_terminal =
                    confirmed_terminal_ends
                            .contains(lhs_left)
                        ? lhs_left
                        : lhs_right;
                const EndKey rhs_terminal =
                    confirmed_terminal_ends
                            .contains(rhs_left)
                        ? rhs_left
                        : rhs_right;
                const EndKey lhs_free =
                    lhs_terminal == lhs_left
                        ? lhs_right
                        : lhs_left;
                const EndKey rhs_free =
                    rhs_terminal == rhs_left
                        ? rhs_right
                        : rhs_left;
                return std::tie(
                           lhs_terminal,
                           lhs_free) <
                       std::tie(
                           rhs_terminal,
                           rhs_free);
            });

        std::vector<AncestralSequencePath>
            component_sequences;
        component_sequences.reserve(
            closed_paths.size() +
            (single_terminal_paths.size() + 1) /
                2 +
            static_cast<size_t>(
                single_terminal_paths.empty() &&
                !free_paths.empty()));
        for (auto& sequence : closed_paths) {
            component_sequences.push_back(
                std::move(sequence));
        }
        auto orient_terminal_left =
            [&](AncestralSequencePath& sequence) {
                if (!confirmed_terminal_ends
                         .contains(
                             endpoint(sequence,
                                      true))) {
                    reverse_path(sequence);
                }
                if (!confirmed_terminal_ends
                         .contains(
                             endpoint(sequence,
                                      true)) ||
                    confirmed_terminal_ends
                        .contains(
                            endpoint(sequence,
                                     false))) {
                    throw std::logic_error(
                        "Expected exactly one confirmed terminal path end");
                }
            };
        auto orient_terminal_right =
            [&](AncestralSequencePath& sequence) {
                if (!confirmed_terminal_ends
                         .contains(
                             endpoint(sequence,
                                      false))) {
                    reverse_path(sequence);
                }
                if (!confirmed_terminal_ends
                         .contains(
                             endpoint(sequence,
                                      false)) ||
                    confirmed_terminal_ends
                        .contains(
                            endpoint(sequence,
                                     true))) {
                    throw std::logic_error(
                        "Expected exactly one confirmed terminal path end");
                }
            };
        auto append_all_free_paths =
            [&](AncestralSequencePath& merged) {
                for (auto& free_path :
                     free_paths) {
                    append_scaffold(
                        merged,
                        std::move(free_path));
                }
                free_paths.clear();
            };

        if (single_terminal_paths.empty()) {
            if (!free_paths.empty()) {
                AncestralSequencePath merged =
                    std::move(free_paths.front());
                for (size_t index = 1;
                     index < free_paths.size();
                     ++index) {
                    append_scaffold(
                        merged,
                        std::move(
                            free_paths[index]));
                }
                component_sequences.push_back(
                    std::move(merged));
            }
        } else {
            size_t index = 0;
            for (;
                 index + 1 <
                 single_terminal_paths.size();
                 index += 2) {
                AncestralSequencePath merged =
                    std::move(
                        single_terminal_paths[
                            index]);
                orient_terminal_left(merged);
                if (index == 0) {
                    append_all_free_paths(
                        merged);
                }
                auto right =
                    std::move(
                        single_terminal_paths[
                            index + 1]);
                orient_terminal_right(right);
                append_scaffold(
                    merged,
                    std::move(right));
                component_sequences.push_back(
                    std::move(merged));
            }
            if (index <
                single_terminal_paths.size()) {
                AncestralSequencePath merged =
                    std::move(
                        single_terminal_paths[
                            index]);
                orient_terminal_left(merged);
                if (index == 0) {
                    append_all_free_paths(
                        merged);
                }
                component_sequences.push_back(
                    std::move(merged));
            }
        }
        std::sort(
            component_sequences.begin(),
            component_sequences.end(),
            stable_path_less);
        for (auto& sequence :
             component_sequences) {
            assembly.sequences.push_back(
                std::move(sequence));
        }
    }

    if (stats != nullptr) {
        stats->path_vertex_count +=
            occurrence_ids.size();
        stats->terminal_end_candidate_count +=
            terminal_children.size();
        stats->confirmed_terminal_end_count +=
            confirmed_terminal_ends.size();
        stats->candidate_adjacency_count +=
            edges.size();
        stats->reference_interval_count +=
            assembly.sequences.size();
        stats->supported_join_count +=
            direct_join_count +
            indirect_join_count;
        stats->direct_join_count +=
            direct_join_count;
        stats->indirect_join_count +=
            indirect_join_count;
        stats->scaffold_join_count +=
            scaffold_join_count;
        stats->scaffold_gap_bases +=
            gap_bases;
    }
    return assembly;
}


std::vector<GenomeSequenceName> buildOutputSequenceOrder(
    const std::vector<std::string>& genome_order,
    const std::vector<GenomeSequenceName>& genome_sequences) {

    std::unordered_map<std::string, size_t> genome_rank;
    genome_rank.reserve(genome_order.size());
    for (size_t i = 0; i < genome_order.size(); ++i) {
        genome_rank.emplace(genome_order[i], i);
    }

    std::vector<GenomeSequenceName> ordered = genome_sequences;
    std::sort(ordered.begin(), ordered.end(), [&](const auto& lhs, const auto& rhs) {
        auto lhs_rank_it = genome_rank.find(lhs.first);
        auto rhs_rank_it = genome_rank.find(rhs.first);
        if (lhs_rank_it == genome_rank.end()) {
            throw std::runtime_error("Unknown genome in output ordering: " + lhs.first);
        }
        if (rhs_rank_it == genome_rank.end()) {
            throw std::runtime_error("Unknown genome in output ordering: " + rhs.first);
        }
        if (lhs_rank_it->second != rhs_rank_it->second) {
            return lhs_rank_it->second < rhs_rank_it->second;
        }
        return lhs.second < rhs.second;
    });

    return ordered;
}

void exportToHal(
    const std::vector<std::weak_ptr<Block>>& blocks,
    const std::filesystem::path& hal_path,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
    NewickParser parser,
    const std::string& root_name,
    const SoftMask::IndexMap& softmask_indexes,
    const ExportConfig& config,
    ExportStats* stats_out) {

    ExportStats local_stats;
    using Clock = std::chrono::steady_clock;

    if (parser.getNodes().empty()) {
        throw std::runtime_error("HAL export requires a non-empty Newick tree");
    }
    if (seqpro_managers.empty()) {
        throw std::runtime_error("HAL export requires non-empty sequence managers");
    }
    if (softmask_indexes.size() != seqpro_managers.size()) {
        throw std::runtime_error(
            "HAL export requires one soft-mask index per leaf genome");
    }
    for (const auto& [species, unused_manager] : seqpro_managers) {
        (void)unused_manager;
        const auto softmask_it = softmask_indexes.find(species);
        if (softmask_it == softmask_indexes.end() || !softmask_it->second) {
            throw std::runtime_error(
                "Missing soft-mask index for leaf genome: " + species);
        }
    }
    if (config.adjacency_theta < 0.0 ||
        config.adjacency_theta > 1.0) {
        throw std::invalid_argument(
            "HAL adjacency theta must be in [0, 1]");
    }
    if (config.phylogenetic_phi < 0.0) {
        throw std::invalid_argument(
            "HAL phylogenetic phi must be non-negative");
    }
    if (config.scaffold_gap_length == 0) {
        throw std::invalid_argument(
            "HAL scaffold gap length must be positive");
    }
    if (normalizeRootNode(parser, root_name)) {
        spdlog::info("HAL export inserted or named the phylogeny root");
    }
    validateLeafNamesExact(parser, seqpro_managers);
    TreeMeta tree = buildTreeMeta(parser);
    local_stats.internal_node_count = tree.internal_postorder.size();

    auto build_msa_begin = Clock::now();
    std::vector<BlockPtr> live_blocks;
    live_blocks.reserve(blocks.size());
    for (const auto& weak_block : blocks) {
        if (auto block = weak_block.lock()) {
            live_blocks.push_back(std::move(block));
        }
    }
    local_stats.block_count = live_blocks.size();

    std::vector<std::optional<BlockMSA>> block_results(live_blocks.size());
    std::exception_ptr build_failure;
    std::mutex build_failure_mutex;
    const int build_threads = std::max(1, config.parallel_threads);
#pragma omp parallel for schedule(dynamic) num_threads(build_threads)
    for (int64_t i = 0; i < static_cast<int64_t>(live_blocks.size()); ++i) {
        try {
            BlockMSA msa = buildBlockMSA(
                live_blocks[static_cast<size_t>(i)],
                tree,
                seqpro_managers,
                &softmask_indexes);
            if (msa.block && !msa.leaf_rows.empty()) {
                block_results[static_cast<size_t>(i)] = std::move(msa);
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(build_failure_mutex);
            if (!build_failure) {
                build_failure = std::current_exception();
            }
        }
    }
    if (build_failure) {
        std::rethrow_exception(build_failure);
    }

    std::vector<BlockMSA> block_msas;
    block_msas.reserve(live_blocks.size());
    for (auto& result : block_results) {
        if (result) {
            block_msas.push_back(std::move(*result));
        }
    }
    std::sort(block_msas.begin(), block_msas.end(), [](const BlockMSA& a, const BlockMSA& b) {
        return exportBlockOrderLess(a.order_key, b.order_key);
    });
    local_stats.msa_count = block_msas.size();
    local_stats.build_msa_ms = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - build_msa_begin).count());
    spdlog::info("HAL export built {} leaf-only block MSAs", block_msas.size());

    auto build_runs_begin = Clock::now();
    auto runs = buildColumnRuns(block_msas, tree);
    const size_t projected_run_count = runs.size();
    const auto secondary_selection =
        selectCoordinateConsistentSecondaryRuns(runs);
    const auto normalization =
        normalizeOverlappingColumnRuns(runs, tree);
    local_stats.run_count = runs.size();
    for (const auto& run : runs) {
        local_stats.observed_occurrence_count +=
            run.leaf_spans.size();
    }
    spdlog::info(
        "HAL export split {} elementary column runs; selected {}/{} "
        "secondary runs for the coordinate-consistent maximum-support "
        "forest ({} conflicting and {} redundant runs rejected, {} bp), "
        "then normalized to {} runs across {} overlapping homology "
        "components ({} overlap constraints, {} duplicate leaf "
        "occurrences removed)",
        projected_run_count,
        secondary_selection.accepted_runs,
        secondary_selection.candidate_runs,
        secondary_selection.conflict_rejected_runs,
        secondary_selection.redundant_runs,
        secondary_selection.rejected_bases,
        normalization.output_runs,
        normalization.connected_components,
        normalization.overlap_constraints,
        normalization.duplicate_leaf_occurrences);
    sanitizeLeafCoverage(runs, seqpro_managers);
    local_stats.build_runs_ms = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - build_runs_begin).count());

    auto infer_presence_begin = Clock::now();
    for (size_t i = 0; i < runs.size(); ++i) {
        auto& run = runs[i];
        auto inference = inferDescendantUnionFast(tree, run.leaf_present);
        run.present_by_node = std::move(inference.present_by_node);
        run.presence_margin = std::move(inference.margin);
        if ((i + 1) % 100000 == 0 || i + 1 == runs.size()) {
            spdlog::info("HAL export inferred binary presence for {}/{} runs", i + 1, runs.size());
        }
    }
    local_stats.infer_presence_ms = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - infer_presence_begin).count());

    std::unordered_map<uint64_t, const ColumnRun*> run_by_id;
    run_by_id.reserve(runs.size());
    for (const auto& run : runs) {
        run_by_id.emplace(run.run_id, &run);
    }

    auto leaf_paths = buildLeafPaths(runs);
    auto build_models_begin = Clock::now();
    auto models = buildNodeModels(
        tree, runs, run_by_id, leaf_paths, config, &local_stats);
    projectInternalReferenceContainersTopDown(
        tree, config, models, &local_stats);
    local_stats.build_models_ms = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - build_models_begin).count());
    for (const auto& [node_id, model] : models) {
        (void)node_id;
        local_stats.scaffold_count += model.sequences.size();
    }

    std::filesystem::path abs_hal_path = std::filesystem::absolute(hal_path);
    if (abs_hal_path.has_parent_path() && !abs_hal_path.parent_path().empty()) {
        std::filesystem::create_directories(abs_hal_path.parent_path());
    }
    std::filesystem::path tmp_hal_path = abs_hal_path;
    tmp_hal_path += ".tmp";
    if (std::filesystem::exists(tmp_hal_path)) {
        std::filesystem::remove(tmp_hal_path);
    }
    struct TemporaryHalGuard {
        std::filesystem::path path;
        bool committed = false;
        ~TemporaryHalGuard() {
            if (!committed) {
                std::error_code error;
                std::filesystem::remove(path, error);
            }
        }
    } temporary_hal_guard{tmp_hal_path};

    std::filesystem::path base_prefix = abs_hal_path;
    base_prefix.replace_extension();

    std::function<void(int)> append_subtree = [&](int node_id) {
        auto emit_begin = Clock::now();
        auto emissions = buildLocalSubtreeEmissions(node_id, tree, models, leaf_paths, run_by_id, seqpro_managers);
        local_stats.emit_subtrees_ms += static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - emit_begin).count());
        accumulateEmissionStats(emissions, &local_stats);

        std::filesystem::path node_prefix = base_prefix;
        node_prefix += "." + tree.nodes[node_id].name;
        std::filesystem::path c2h_path = node_prefix;
        c2h_path += ".c2h";
        std::filesystem::path hal_fa_path = node_prefix;
        hal_fa_path += ".hal.fa";
        std::filesystem::path tree_path = node_prefix;
        tree_path += ".newick";
        struct GeneratedFileGuard {
            std::array<std::filesystem::path, 3> paths;
            ~GeneratedFileGuard() {
                for (const auto& path : paths) {
                    std::error_code error;
                    std::filesystem::remove(path, error);
                }
            }
        } generated_file_guard{{c2h_path, hal_fa_path, tree_path}};

        auto genome_order = buildLocalGenomeOrder(tree, node_id);
        writeC2H(c2h_path, genome_order, emissions);
        writeHalFasta(
            hal_fa_path,
            genome_order,
            emissions,
            seqpro_managers,
            softmask_indexes);
        writeTreeFile(tree_path, buildLocalNewick(tree, node_id));
        spdlog::info("HAL export wrote {}", c2h_path.string());
        spdlog::info("HAL export wrote {}", hal_fa_path.string());

        auto hal_append_begin = Clock::now();
        runHalAppend(c2h_path, hal_fa_path, tree_path, tmp_hal_path);
        local_stats.hal_append_ms += static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - hal_append_begin).count());
        for (int child_id : tree.nodes[node_id].children) {
            if (!tree.nodes[child_id].is_leaf) {
                append_subtree(child_id);
            }
        }
    };
    append_subtree(tree.root_id);


    if (std::filesystem::exists(abs_hal_path)) {
        std::filesystem::remove(abs_hal_path);
    }
    std::filesystem::rename(tmp_hal_path, abs_hal_path);
    temporary_hal_guard.committed = true;
    spdlog::info("HAL export finished {}", abs_hal_path.string());
    spdlog::info(
        "HAL export stats: {} sequences from {} path vertices "
        "({} supported joins, {} terminal-end candidates, "
        "{} confirmed terminal ends, {} candidate adjacencies); "
        "{} observed and {} ancestral occurrences; "
        "joins direct={}, indirect={}, scaffold={}, gap bases={}; "
        "{} aligned tops, {} paralogous tops, {} paralogy self-adjacencies",
        local_stats.scaffold_count,
        local_stats.path_vertex_count,
        local_stats.supported_join_count,
        local_stats.terminal_end_candidate_count,
        local_stats.confirmed_terminal_end_count,
        local_stats.candidate_adjacency_count,
        local_stats.observed_occurrence_count,
        local_stats.ancestor_occurrence_count,
        local_stats.direct_join_count,
        local_stats.indirect_join_count,
        local_stats.scaffold_join_count,
        local_stats.scaffold_gap_bases,
        local_stats.aligned_top_count,
        local_stats.paralogous_top_count,
        local_stats.paralogy_self_adjacency_count);
    if (stats_out != nullptr) {
        *stats_out = local_stats;
    }
}

} // namespace RaMesh::hal_export
