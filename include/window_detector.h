#ifndef RAMAX_WINDOW_DETECTOR_H
#define RAMAX_WINDOW_DETECTOR_H

#include <cstdint>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include "ramesh.h"

namespace RaMesh::WindowDetection {

enum class DetectionMode : uint8_t {
    EACH_ROUND = 0,
    FINAL_ONLY = 1
};

struct Options {
    bool enabled = false;
    DetectionMode mode = DetectionMode::EACH_ROUND;
    std::filesystem::path report_dir;
    std::string threshold_profile = "alignathon-v1";

    uint64_t micro_block_max_bp = 10;
    uint64_t short_block_max_bp = 100;
    uint64_t weak_block_upper_bp = 500;
    uint64_t primary_gap_max_bp = 100;
    uint64_t extended_gap_max_bp = 500;
    uint64_t hard_boundary_gap_bp = 1000;
    uint64_t anchor_min_segment_bp = 100;
    uint64_t strong_anchor_bp = 500;
    uint64_t max_window_span_bp = 100000;
    size_t subset_search_budget = 100000;
};

struct SegmentSnapshot {
    std::string segment_id;
    std::string block_id;
    std::string species;
    std::string chromosome;
    uint64_t graph_start = 0;
    uint64_t graph_end = 0;
    uint64_t original_start = 0;
    uint64_t original_end = 0;
    bool original_coordinates_resolved = false;
    bool reverse = false;
    bool primary = true;
    bool left_extended = false;
    bool right_extended = false;
    size_t path_index = 0;
    std::string cigar_summary;
};

struct BlockSnapshot {
    std::string block_id;
    std::string canonical_signature;
    std::string diagnostic_ref_chromosome;
    std::vector<SegmentSnapshot> segments;
    std::vector<std::string> participating_genomes;
    uint64_t min_segment_length = 0;
    uint64_t max_segment_length = 0;
    bool copy_ambiguity = false;
    bool has_secondary = false;
};

struct PathSnapshot {
    std::string species;
    std::string chromosome;
    std::vector<SegmentSnapshot> segments;
};

struct AuditRecord {
    std::string severity;
    std::string code;
    std::string species;
    std::string chromosome;
    uint64_t graph_position = 0;
    std::string message;
};

struct GraphSnapshot {
    uint64_t round_id = 0;
    std::string current_reference;
    std::vector<std::string> input_genomes;
    std::vector<BlockSnapshot> blocks;
    std::vector<PathSnapshot> paths;
    std::vector<AuditRecord> audit_records;
    uint64_t active_segment_count = 0;
};

struct BoundaryGenomeEvidence {
    std::string target_species;
    bool left_present = false;
    bool right_present = false;
    bool unique_target = false;
    std::string left_chromosome;
    uint64_t left_graph_start = 0;
    uint64_t left_graph_end = 0;
    uint64_t left_original_start = 0;
    uint64_t left_original_end = 0;
    bool left_original_resolved = false;
    bool left_reverse = false;
    std::string right_chromosome;
    uint64_t right_graph_start = 0;
    uint64_t right_graph_end = 0;
    uint64_t right_original_start = 0;
    uint64_t right_original_end = 0;
    bool right_original_resolved = false;
    bool right_reverse = false;
    int64_t target_gap = 0;
    bool target_gap_resolved = false;
    bool order_consistent = false;
    bool strand_consistent = false;
    bool target_consistent = false;
};

struct BoundaryEvidence {
    std::string boundary_id;
    std::string reference_species;
    std::string reference_chromosome;
    size_t path_boundary_index = 0;
    std::string left_block_id;
    std::string right_block_id;
    SegmentSnapshot left_reference_segment;
    SegmentSnapshot right_reference_segment;
    int64_t reference_gap = 0;
    uint64_t reference_overlap = 0;
    uint64_t min_adjacent_block_length = 0;
    uint64_t max_adjacent_block_length = 0;
    bool hard_boundary = false;
    bool copy_ambiguity = false;
    bool detector_recommends_bridge = false;
    std::vector<std::string> signals;
    std::vector<std::string> dropout_genomes;
    std::vector<BoundaryGenomeEvidence> genomes;
};

struct Seed {
    std::string seed_id;
    std::string boundary_id;
    std::string species;
    std::string chromosome;
    size_t boundary_index = 0;
    std::vector<std::string> signals;
    std::string initial_priority;
    std::string reason;
};

struct WindowGenomeEvidence {
    std::string species;
    std::string chromosome;
    uint64_t graph_start = 0;
    uint64_t graph_end = 0;
    uint64_t original_start = 0;
    uint64_t original_end = 0;
    bool original_coordinates_resolved = false;
    std::string strand;
    bool left_anchor_present = false;
    bool right_anchor_present = false;
    bool currently_aligned = false;
    bool included_in_max_k = false;
    std::string excluded_reason;
};

struct CandidateWindow {
    std::string window_id;
    uint64_t round_id = 0;
    std::string current_reference;
    std::string report_species;
    std::string report_chromosome;
    size_t path_left_index = 0;
    size_t path_right_index = 0;
    uint64_t graph_start = 0;
    uint64_t graph_end = 0;
    uint64_t original_start = 0;
    uint64_t original_end = 0;
    bool original_coordinates_resolved = false;
    std::string left_anchor_id;
    std::string right_anchor_id;
    std::string left_anchor_level;
    std::string right_anchor_level;
    bool has_two_sided_anchors = false;
    std::vector<std::string> block_ids;
    std::vector<std::string> boundary_ids;
    std::vector<std::string> seed_boundary_ids;
    std::vector<std::string> signals;
    std::vector<std::string> complex_flags;
    std::string priority_tier;
    size_t input_genome_count = 0;
    size_t currently_aligned_genome_count = 0;
    size_t max_possible_k = 0;
    std::vector<std::string> max_k_genomes;
    std::vector<std::string> excluded_genomes;
    uint64_t min_block_length = 0;
    uint64_t max_block_length = 0;
    int64_t max_reference_gap = 0;
    std::vector<WindowGenomeEvidence> genomes;
    std::string detector_reason;
};

struct DetectionResult {
    GraphSnapshot snapshot;
    std::vector<BoundaryEvidence> boundaries;
    std::vector<Seed> seeds;
    std::vector<CandidateWindow> windows;
    double snapshot_seconds = 0.0;
    double detection_seconds = 0.0;
};

GraphSnapshot snapshotGraph(
    const RaMeshMultiGenomeGraph& graph,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    uint64_t round_id,
    const SpeciesName& current_reference);

DetectionResult detectProblemWindows(GraphSnapshot snapshot, const Options& options);

void writeDetectionReport(const DetectionResult& result, const Options& options);

DetectionResult detectAndWriteProblemWindows(
    const RaMeshMultiGenomeGraph& graph,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    uint64_t round_id,
    const SpeciesName& current_reference,
    const Options& options);

std::string detectionModeToString(DetectionMode mode);
DetectionMode detectionModeFromString(const std::string& value);

}  // namespace RaMesh::WindowDetection

#endif  // RAMAX_WINDOW_DETECTOR_H
