#ifndef RARE_ALIGNER_H
#define RARE_ALIGNER_H

#include <optional>
#include <atomic>
#include <exception>
#include <limits>
#include <memory>
#include <mutex>
#include <vector>

#include "config.hpp"
#include "index.h"
#include "suffix_array_index.h"
#include "SeqPro.h"
#include "ramesh.h"
#include "structural_break_repair.h"
#include "short_block_repair.h"
#include "cache_manifest.h"
#include "mash_distance_estimator.h"
#include "wfmash_router.h"

void addAlignedRegionsAsMask(
    const RaMesh::RaMeshMultiGenomeGraph& graph,
    std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
    const SpeciesName& ref_name);

struct IndexCacheCounters {
    std::atomic_size_t memory_only_built{0};
};

// 多基因组比对核心调度类
class MultipleRareAligner {
public:
    FilePath work_dir;
    FilePath index_dir;

    SpeciesPathMap species_path_map;

    uint_t chunk_size;
    uint_t overlap_size;
    uint_t min_anchor_length;
    uint_t max_anchor_frequency;
    uint_t accurate_skip_threshold;
    uint_t sa_sampling_rate{1};
    bool allow_mem = false;
    SpeciesClusterMap secondary_cluster_map;
    SpeciesMatchVec3DPtrMap secondary_match_map;

    uint_t thread_num;
    uint_t group_id;
    uint_t round_id;
    bool trust_legacy_cache{false};
    std::shared_ptr<IndexCacheCounters> index_cache_counters{
        std::make_shared<IndexCacheCounters>()};

    bool enable_mask_export = false;
    bool mask_export_done = false;
    FilePath mask_export_dir;

    bool merge_exact_contiguous_blocks_enabled = true;
    uint_t merge_query_gap_max = 100;
    bool realign_single_missing_species_enabled = true;
    uint_t species_mismatch_realign_max_span = 3000;
    uint_t species_mismatch_zero_gap_max_span = 200;
    std::string species_mismatch_msa_executable;
    RaMesh::StructuralBreakRepair::Options structural_break_repair_options;
    RaMesh::ShortBlockRepair::Options short_block_repair_options;
    std::vector<MashDistanceRecord> first_reference_distances;
    double near_distance_threshold{0.01};
    double far_distance_threshold{0.02};

    // 构造函数声明：注意名称必须与类名完全一致
    MultipleRareAligner(
        const FilePath& work_dir,
        SpeciesPathMap& species_path_map,
        uint_t thread_num,
        uint_t chunk_size,
        uint_t overlap_size,
        uint_t min_anchor_length,
        uint_t max_anchor_frequency,
        uint_t accurate_skip_threshold,
        bool trust_legacy_cache = false
    );

    std::unique_ptr<RaMesh::RaMeshMultiGenomeGraph> starAlignment(
    std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers,
    std::string ref_name,
    bool only_one_round,
    bool                       fast_build,
    SeqPro::Length sampling_interval=32,
    uint_t min_span=65);

    SpeciesMatchVec3DPtrMapPtr alignMultipleGenome(SpeciesName ref_name, std::unordered_map<SpeciesName, SeqPro::SharedManagerVariant>& species_fasta_manager_map, SearchMode search_mode, bool fast_build, bool allow_MEM, bool allow_short_mum, sdsl::int_vector<0>& ref_global_cache, SeqPro::Length sampling_interval);

    void alignClusterConstructBounded(
        const SpeciesName& ref_name,
        std::unordered_map<SpeciesName, SeqPro::SharedManagerVariant>& species_fasta_manager_map,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
        SearchMode search_mode,
        bool fast_build,
        bool allow_MEM,
        bool allow_short_mum,
        sdsl::int_vector<0>& ref_global_cache,
        SeqPro::Length sampling_interval,
        uint_t min_span,
        bool is_first,
        RaMesh::RaMeshMultiGenomeGraph& graph);

    SpeciesClusterMapPtr filterMultipeSpeciesAnchors(SpeciesName ref_name, std::unordered_map<SpeciesName, SeqPro::SharedManagerVariant>& species_fm_map, SpeciesMatchVec3DPtrMapPtr species_match_map, uint_t min_span);

    void constructMultipleGraphsByGreedy(std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers, SpeciesName ref_name, const SpeciesClusterMap& species_cluster_map, RaMesh::RaMeshMultiGenomeGraph& graph, uint_t min_span);

    void constructMultipleGraphsByDpByRef(std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers, SpeciesName ref_name, const SpeciesClusterMap& species_cluster_map, RaMesh::RaMeshMultiGenomeGraph& graph, uint_t min_span);

    void constructMultipleGraphsByGreedyByRef(std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers, SpeciesName ref_name, const SpeciesClusterMap& species_cluster_map, RaMesh::RaMeshMultiGenomeGraph& graph, uint_t min_span);

    void constructMultipleGraphsByDp(const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers, SpeciesName ref_name, const SpeciesClusterMap& species_cluster_map, RaMesh::RaMeshMultiGenomeGraph& graph, uint_t min_span, bool is_first);

private:
    // 辅助函数：处理基本区间
    static RaMesh::BlockPtr processElementaryInterval(
        uint_t interval_start, uint_t interval_end,
        const std::vector<RaMesh::BlockPtr>& valid_blocks,
        const SpeciesName& ref_name
    );

};

struct PreparedAnchorSearch {
    struct Task {
        Region chunk;
        Strand strand =
            Strand::FORWARD;
    };

    SeqPro::ManagerVariant* query_manager =
        nullptr;
    SearchMode search_mode =
        ACCURATE_SEARCH;
    bool allow_mem = false;
    bool allow_short_mum = false;
    sdsl::int_vector<0>* ref_global_cache =
        nullptr;
    SeqPro::Length sampling_interval = 32;
    std::vector<Task> tasks;
    // Each task owns exactly one value slot.  This avoids one shared_ptr
    // control block and one exception_ptr slot per chunk/strand task while
    // retaining the original task order for collection.
    std::vector<std::optional<MatchVec2D>> task_results;
    std::mutex failure_mutex;
    std::exception_ptr first_failure;
    size_t first_failure_index = std::numeric_limits<size_t>::max();
};



class PairRareAligner {
public:
    FilePath work_dir;
    FilePath index_dir;

    SpeciesName ref_name;

    SeqPro::ManagerVariant* ref_seqpro_manager;
    std::optional<Suffix_Array_Index> ref_index;

    uint_t chunk_size;
    uint_t overlap_size;
    uint_t min_anchor_length;
    uint_t max_anchor_frequency;
    uint_t accurate_skip_threshold = 0;
    uint_t sa_sampling_rate = 1;

    uint_t group_id;
    uint_t round_id;

    uint_t thread_num;
    PairRareAligner(const FilePath work_dir, const uint_t thread_num, uint_t chunk_size,
        uint_t overlap_size, uint_t min_anchor_length, uint_t max_anchor_frequency,
        uint_t accurate_skip_threshold = 0,
        bool trust_legacy_cache = false);

    // 新增构造函数：从 MultipleRareAligner 初始化
    PairRareAligner(const MultipleRareAligner& mra)
        : work_dir(mra.work_dir),
        index_dir(mra.index_dir),
        chunk_size(mra.chunk_size),
        overlap_size(mra.overlap_size),
        min_anchor_length(mra.min_anchor_length),
        max_anchor_frequency(mra.max_anchor_frequency),
        accurate_skip_threshold(mra.accurate_skip_threshold),
        sa_sampling_rate(mra.sa_sampling_rate),
        thread_num(mra.thread_num),
        trust_legacy_cache(mra.trust_legacy_cache),
        index_cache_counters(mra.index_cache_counters)
    {
        this->group_id = mra.group_id;
        this->round_id = mra.round_id;
    }

    bool trust_legacy_cache{false};
    std::shared_ptr<IndexCacheCounters> index_cache_counters{
        std::make_shared<IndexCacheCounters>()};

    MatchVec3DPtr alignPairGenome(SpeciesName query_name, SeqPro::ManagerVariant& query_fasta_manager, SearchMode search_mode, bool allow_MEM, bool allow_short_mum, sdsl::int_vector<0>& ref_global_cache, SeqPro::Length sampling_interval);
    FilePath buildIndex(const std::string prefix, SeqPro::ManagerVariant& ref_fasta_manager_, bool fast_build);


    MatchVec3DPtr findQueryFileAnchor(const std::string prefix, SeqPro::ManagerVariant& query_fasta_manager, SearchMode search_mode, bool allow_MEM, bool allow_short_mum, sdsl::int_vector<0>& ref_global_cache, SeqPro::Length sampling_interval, bool isMultiple=false, bool include_masked_regions=false);
    std::shared_ptr<PreparedAnchorSearch>
    prepareQueryFileAnchor(
        SeqPro::ManagerVariant& query_fasta_manager,
        SearchMode search_mode,
        bool allow_mem,
        bool allow_short_mum,
        sdsl::int_vector<0>& ref_global_cache,
        SeqPro::Length sampling_interval,
        bool is_multiple = false,
        bool include_masked_regions = false);

    void executePreparedAnchorTask(
        PreparedAnchorSearch& plan,
        size_t task_index);

    MatchVec3DPtr collectPreparedAnchorSearch(
        PreparedAnchorSearch& plan);



    void constructGraphByGreedy(SpeciesName query_name, SeqPro::ManagerVariant& query_seqpro_manager, ClusterVecPtrByStrandByQueryRefPtr cluster_ptr, RaMesh::RaMeshMultiGenomeGraph& graph, uint_t min_span);

    void constructGraphByGreedyByRef(SpeciesName query_name, SeqPro::ManagerVariant& query_seqpro_manager, MatchClusterVecPtr cluster_vec_ptr, RaMesh::RaMeshMultiGenomeGraph& graph, 
        uint_t min_span, bool isMultiple=false);
    ClusterVecPtrByStrandByQueryRefPtr filterPairSpeciesAnchors(SpeciesName query_name, MatchVec3DPtr& anchors, SeqPro::ManagerVariant& query_fasta_manager, RaMesh::RaMeshMultiGenomeGraph& graph, uint_t min_span);

    AnchorBySQR_SparsePtr extendClusterToAnchorByChr(SpeciesName query_name, SeqPro::ManagerVariant& query_seqpro_manager, ClusterBySQR_SparsePtr cluster, bool is_first);


    void filterAnchorByDP(AnchorBySQR_SparsePtr anchor_map,uint_t ref_chr_cnt, uint_t qry_chr_cnt);
    AnchorPtrVec extendClusterGroupToAnchors(
        SeqPro::ManagerVariant& query_seqpro_manager,
        MatchClusterVec& cluster_group,
        bool is_first);

    void filterAnchorByDPDimension(
        AnchorBySQR_SparsePtr anchor_map,
        uint_t chromosome_id,
        bool filter_ref);

    void filterAnchorVectorByDP(
        AnchorPtrVec anchors,
        bool filter_ref);

    static uint64_t dpTreapFallbackCount();


    void constructGraphByDP(const SpeciesName& query_name,
                            SeqPro::ManagerVariant& query_seqpro_manager,
                            AnchorBySQR_SparsePtr anchor_ptr,
                            RaMesh::RaMeshMultiGenomeGraph& graph);
    void registerSecondaryAnchors(SpeciesName query_name, SeqPro::ManagerVariant& query_seqpro_manager, AnchorBySQR_SparsePtr anchor_ptr, RaMesh::RaMeshMultiGenomeGraph& graph, bool initial_round);

};

#ifdef RAMAX_PERFORMANCE_TEST_HOOKS
void ramaxFilterAnchorsByDpOptimizedForTesting(
    AnchorPtrVec anchors, bool filter_ref);
void ramaxFilterAnchorsByDpLegacyForTesting(
    AnchorPtrVec anchors, bool filter_ref);
#endif

#endif
