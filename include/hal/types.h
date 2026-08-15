#pragma once

#include "config.hpp"

#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace RaMesh::hal_export {
using OccurrenceId = uint64_t;


struct TreeNodeMeta {
    int id = -1;
    std::string name;
    int parent = -1;
    std::vector<int> children;
    double branch_length_to_parent = 0.0;
    bool is_leaf = false;
    int leaf_index = -1;
};

struct TreeMeta {
    std::vector<TreeNodeMeta> nodes;
    std::unordered_map<std::string, int> name_to_id;
    std::vector<int> leaf_ids;
    std::vector<int> internal_postorder;
    int root_id = -1;
};

struct ExportConfig {
    double consensus_threshold = 0.6;
    int parallel_threads = 1;
    double adjacency_theta = 0.000001;
    double phylogenetic_phi = 1.0;
    uint32_t scaffold_gap_length = 10;
};

struct ExportStats {
    uint64_t block_count = 0;
    uint64_t msa_count = 0;
    uint64_t run_count = 0;
    uint64_t internal_node_count = 0;
    uint64_t scaffold_count = 0;
    uint64_t aligned_top_count = 0;
    uint64_t insertion_top_count = 0;
    uint64_t reverse_top_count = 0;
    uint64_t paralogy_self_adjacency_count = 0;
    uint64_t terminal_end_candidate_count = 0;
    uint64_t supported_join_count = 0;
    uint64_t paralogous_top_count = 0;
    uint64_t observed_occurrence_count = 0;
    uint64_t ancestor_occurrence_count = 0;
    uint64_t path_vertex_count = 0;
    uint64_t confirmed_terminal_end_count = 0;
    uint64_t candidate_adjacency_count = 0;
    uint64_t reference_interval_count = 0;
    uint64_t direct_join_count = 0;
    uint64_t indirect_join_count = 0;
    uint64_t scaffold_join_count = 0;
    uint64_t scaffold_gap_bases = 0;
    uint64_t build_msa_ms = 0;
    uint64_t build_runs_ms = 0;
    uint64_t infer_presence_ms = 0;
    uint64_t build_models_ms = 0;
    uint64_t emit_subtrees_ms = 0;
    uint64_t hal_append_ms = 0;
    uint64_t root_sequence_count = 0;
    uint64_t root_singleton_sequence_count = 0;
    uint64_t root_build_models_ms = 0;
    uint64_t non_root_build_models_ms = 0;
    uint64_t root_candidate_filter_ms = 0;
    uint64_t non_root_candidate_filter_ms = 0;
    uint64_t root_run_dna_ms = 0;
    uint64_t non_root_run_dna_ms = 0;
    uint64_t root_edge_collect_ms = 0;
    uint64_t non_root_edge_collect_ms = 0;
    uint64_t root_path_decompose_ms = 0;
    uint64_t non_root_path_decompose_ms = 0;
    uint64_t root_sequence_materialize_ms = 0;
    uint64_t non_root_sequence_materialize_ms = 0;
};

struct BinaryInferenceResult {
    std::vector<uint8_t> present_by_node;
    std::unordered_map<std::string, uint8_t> present_by_name;
    std::vector<double> score0;
    std::vector<double> score1;
    std::vector<double> margin;
};

struct LeafInterval {
    uint64_t start = 0;
    uint32_t length = 0;
    bool forward_to_parent = true;
};

struct AlignedOccurrence {
    std::string row_id;
    std::string genome_name;
    std::string sequence_name;
    uint64_t segment_start = 0;
    uint32_t segment_length = 0;
    bool reversed = false;
    std::string aligned_dna;
};

struct ProjectedOccurrence {
    std::string row_id;
    std::string genome_name;
    std::string sequence_name;
    uint64_t start = 0;
    uint32_t length = 0;
    bool reversed = false;
    std::string dna;
};

struct ElementaryRunProjection {
    uint32_t col_beg = 0;
    uint32_t col_end = 0;
    std::vector<ProjectedOccurrence> occurrences;
};

using GenomeSequenceName = std::pair<std::string, std::string>;

struct RunOrderKey {
    uint64_t block_id = 0;
    uint32_t col_beg = 0;
};

struct AdjacencyVote {
    uint64_t left_run_id = 0;
    bool left_forward_to_canonical = true;
    uint64_t right_run_id = 0;
    bool right_forward_to_canonical = true;
    uint32_t left_length = 0;
    uint32_t right_length = 0;
    uint64_t gap_bases = 0;
    uint32_t left_copy_index = 0;
    uint32_t right_copy_index = 0;
};

struct EdgeSupport {
    OccurrenceId from = 0;
    OccurrenceId to = 0;
    uint32_t occurrence_support = 0;
    std::vector<int> supporting_children;
    bool from_forward_to_canonical = true;
    bool to_forward_to_canonical = true;
    long double weighted_support = 0.0;
    uint64_t minimum_gap = std::numeric_limits<uint64_t>::max();
};
enum class OccurrenceEndSide : uint8_t {
    LEFT = 0,
    RIGHT = 1,
};

struct OccurrenceEnd {
    OccurrenceId occurrence_id = 0;
    OccurrenceEndSide side = OccurrenceEndSide::LEFT;

    bool operator==(const OccurrenceEnd&) const = default;
};

struct TerminalEndSupport {
    OccurrenceEnd end;
    uint32_t occurrence_support = 0;
    std::vector<int> supporting_lineages;
    long double weighted_support = 0.0;
};


enum class ReferenceJoinKind : uint8_t {
    DIRECT = 0,
    INDIRECT = 1,
    SCAFFOLD = 2,
};

struct ReferenceJoin {
    ReferenceJoinKind kind = ReferenceJoinKind::DIRECT;
    uint32_t gap_length = 0;
    long double score = 0.0;
};

struct AncestralSequencePath {
    std::vector<std::vector<OccurrenceId>> supported_fragments;
    std::vector<ReferenceJoin> joins;
};

struct AncestralSequenceAssembly {
    std::vector<AncestralSequencePath> sequences;
    std::unordered_map<OccurrenceId, bool> forward_by_occurrence;
};

struct BucketedDonor {
    int bucket_id = -1;
    bool prefer_internal = false;
    std::string dna;
    double weight = 0.0;
};

struct ExportBlockOrderKey {
    std::string ref_species;
    std::string ref_chr;
    uint64_t ref_start = 0;
    uint64_t block_id = 0;
};

} // namespace RaMesh::hal_export
