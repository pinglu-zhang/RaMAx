#pragma once

#include "SeqPro.h"
#include "data_process.h"
#include "hal/types.h"
#include "ramesh.h"

#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace RaMesh::hal_export {

TreeMeta buildTreeMeta(const NewickParser& parser);

BinaryInferenceResult inferDescendantUnion(
    const TreeMeta& tree,
    const std::unordered_map<std::string, bool>& leaf_presence);

LeafInterval projectLeafInterval(
    uint64_t segment_start,
    uint32_t segment_length,
    bool reversed,
    uint32_t non_gap_before,
    uint32_t run_length);

std::vector<ElementaryRunProjection> projectElementaryRuns(
    const std::vector<AlignedOccurrence>& rows);

int c2hHasBottomFlag(size_t bottom_count);

bool computeForwardToParent(bool child_forward_to_canonical, bool parent_forward_to_canonical);

std::optional<AdjacencyVote> orientAdjacencyVote(
    uint64_t left_run_id,
    bool left_forward_to_canonical,
    uint64_t right_run_id,
    bool right_forward_to_canonical,
    uint32_t left_length = 0,
    uint32_t right_length = 0,
    uint64_t gap_bases = 0);
long double calculateCactusZScore(
    uint64_t left_length,
    uint64_t right_length,
    uint64_t gap_length,
    long double theta);


std::string orientRunDNAForPlacement(const std::string& canonical_dna, bool forward);

std::string buildConsensusDNA(
    const std::vector<std::pair<std::string, double>>& donors,
    size_t expected_length,
    double consensus_threshold = 0.6);

std::vector<std::pair<std::string, double>> selectBestDonorsByBucket(
    const std::vector<BucketedDonor>& donors);

bool exportBlockOrderLess(const ExportBlockOrderKey& lhs, const ExportBlockOrderKey& rhs);

std::vector<std::vector<uint64_t>> buildMaximumCardinalityWeightPathCover(
    const std::vector<uint64_t>& run_ids,
    const std::vector<EdgeSupport>& edges,
    const std::unordered_map<uint64_t, RunOrderKey>& run_order_keys,
    ExportStats* stats = nullptr);

AncestralSequenceAssembly buildAncestralSequenceAssembly(
    const std::vector<uint64_t>& occurrence_ids,
    const std::vector<EdgeSupport>& edges,
    const std::unordered_map<uint64_t, RunOrderKey>& occurrence_order_keys,
    const std::vector<TerminalEndSupport>& terminal_ends,
    uint32_t scaffold_gap_length,
    ExportStats* stats = nullptr);

std::vector<GenomeSequenceName> buildOutputSequenceOrder(
    const std::vector<std::string>& genome_order,
    const std::vector<GenomeSequenceName>& genome_sequences);
void exportToMaf(
    const std::vector<std::weak_ptr<Block>>& blocks,
    const std::filesystem::path& maf_path,
    const std::map<
        SpeciesName,
        SeqPro::SharedManagerVariant>&
        seqpro_managers,
    bool pairwise_mode = false);

std::unordered_set<uint64_t>
findRejectedSecondaryHomologyBlocks(
    const std::vector<std::weak_ptr<Block>>& blocks,
    const std::map<
        SpeciesName,
        SeqPro::SharedManagerVariant>&
        seqpro_managers);


void exportToHal(
    const std::vector<std::weak_ptr<Block>>& blocks,
    const std::filesystem::path& hal_path,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
    NewickParser parser,
    const std::string& root_name,
    const SoftMask::IndexMap& softmask_indexes,
    const ExportConfig& config = {},
    ExportStats* stats_out = nullptr);

} // namespace RaMesh::hal_export
