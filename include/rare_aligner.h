#ifndef RARE_ALIGNER_H
#define RARE_ALIGNER_H

#include <optional>

#include "config.hpp"
#include "index.h"
#include "SeqPro.h"
#include "threadpool.h"
#include "ramesh.h"
#include "window_detector.h"

// 多基因组比对核心调度类
class MultipleRareAligner {
public:
    FilePath work_dir;
    FilePath index_dir;

    SpeciesPathMap species_path_map;
    NewickParser newick_tree;

    uint_t chunk_size;
    uint_t overlap_size;
    uint_t min_anchor_length;
    uint_t max_anchor_frequency;

    uint_t thread_num;
    uint_t group_id;
    uint_t round_id;

    bool enable_mask_export = false;
    bool mask_export_done = false;
    FilePath mask_export_dir;

    RaMesh::WindowDetection::Options window_detection_options;

    // 构造函数声明：注意名称必须与类名完全一致
    MultipleRareAligner(
        const FilePath& work_dir,
        SpeciesPathMap& species_path_map,
        NewickParser& newick_tree,
        uint_t thread_num,
        uint_t chunk_size,
        uint_t overlap_size,
        uint_t min_anchor_length,
        uint_t max_anchor_frequency
    );

    std::unique_ptr<RaMesh::RaMeshMultiGenomeGraph> starAlignment(
    std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers,
    std::string ref_name,
    bool only_one_round,
    bool                       fast_build,
    SeqPro::Length sampling_interval=32,
    uint_t min_span=65);

    SpeciesMatchVec3DPtrMapPtr alignMultipleGenome(SpeciesName ref_name, std::unordered_map<SpeciesName, SeqPro::SharedManagerVariant>& species_fasta_manager_map, SearchMode search_mode, bool fast_build, bool allow_MEM, bool allow_short_mum, sdsl::int_vector<0>& ref_global_cache, SeqPro::Length sampling_interval);

    SpeciesClusterMapPtr filterMultipeSpeciesAnchors(SpeciesName ref_name, std::unordered_map<SpeciesName, SeqPro::SharedManagerVariant>& species_fm_map, SpeciesMatchVec3DPtrMapPtr species_match_map, uint_t min_span);

    void constructMultipleGraphsByGreedy(std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers, SpeciesName ref_name, const SpeciesClusterMap& species_cluster_map, RaMesh::RaMeshMultiGenomeGraph& graph, uint_t min_span);

    void constructMultipleGraphsByDpByRef(std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers, SpeciesName ref_name, const SpeciesClusterMap& species_cluster_map, RaMesh::RaMeshMultiGenomeGraph& graph, uint_t min_span);

    void constructMultipleGraphsByGreedyByRef(std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers, SpeciesName ref_name, const SpeciesClusterMap& species_cluster_map, RaMesh::RaMeshMultiGenomeGraph& graph, uint_t min_span);

    void constructMultipleGraphsByDp(std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers, SpeciesName ref_name, const SpeciesClusterMap& species_cluster_map, RaMesh::RaMeshMultiGenomeGraph& graph, uint_t min_span, bool is_first);

private:
    // 辅助函数：处理基本区间
    static RaMesh::BlockPtr processElementaryInterval(
        uint_t interval_start, uint_t interval_end,
        const std::vector<RaMesh::BlockPtr>& valid_blocks,
        const SpeciesName& ref_name
    );

};


class PairRareAligner {
public:
    FilePath work_dir;
    FilePath index_dir;

    SpeciesName ref_name;

    SeqPro::ManagerVariant* ref_seqpro_manager;
    std::optional<FM_Index> ref_index;

    uint_t chunk_size;
    uint_t overlap_size;
    uint_t min_anchor_length;
    uint_t max_anchor_frequency;

    uint_t group_id;
    uint_t round_id;

    uint_t thread_num;
    PairRareAligner(const FilePath work_dir, const uint_t thread_num, uint_t chunk_size, uint_t overlap_size, uint_t min_anchor_length, uint_t max_anchor_frequency);

    // 新增构造函数：从 MultipleRareAligner 初始化
    PairRareAligner(const MultipleRareAligner& mra)
        : work_dir(mra.work_dir),
        index_dir(mra.index_dir),
        chunk_size(mra.chunk_size),
        overlap_size(mra.overlap_size),
        min_anchor_length(mra.min_anchor_length),
        max_anchor_frequency(mra.max_anchor_frequency),
        thread_num(mra.thread_num)
    {
        this->group_id = mra.group_id;
        this->round_id = mra.round_id;
    }

    MatchVec3DPtr alignPairGenome(SpeciesName query_name, SeqPro::ManagerVariant& query_fasta_manager, SearchMode search_mode, bool allow_MEM, bool allow_short_mum, sdsl::int_vector<0>& ref_global_cache, SeqPro::Length sampling_interval);
    FilePath buildIndex(const std::string prefix, SeqPro::ManagerVariant& ref_fasta_manager_, bool fast_build);


    MatchVec3DPtr findQueryFileAnchor(const std::string prefix, SeqPro::ManagerVariant& query_fasta_manager, SearchMode search_mode, bool allow_MEM, bool allow_short_mum, sdsl::int_vector<0>& ref_global_cache, SeqPro::Length sampling_interval, bool isMultiple=false);

    void constructGraphByGreedy(SpeciesName query_name, SeqPro::ManagerVariant& query_seqpro_manager, ClusterVecPtrByStrandByQueryRefPtr cluster_ptr, RaMesh::RaMeshMultiGenomeGraph& graph, uint_t min_span);

    void constructGraphByGreedyByRef(SpeciesName query_name, SeqPro::ManagerVariant& query_seqpro_manager, MatchClusterVecPtr cluster_vec_ptr, RaMesh::RaMeshMultiGenomeGraph& graph, 
        uint_t min_span, bool isMultiple=false);
    void constructGraphByDpByRef(SpeciesName query_name, SeqPro::ManagerVariant& query_seqpro_manager, MatchClusterVecPtr cluster_vec_ptr, RaMesh::RaMeshMultiGenomeGraph& graph, ThreadPool& pool, uint_t thread_num, uint_t min_span, bool isMultiple);
    
    
    ClusterVecPtrByStrandByQueryRefPtr filterPairSpeciesAnchors(SpeciesName query_name, MatchVec3DPtr& anchors, SeqPro::ManagerVariant& query_fasta_manager, RaMesh::RaMeshMultiGenomeGraph& graph, uint_t min_span);

    AnchorBySQR_SparsePtr extendClusterToAnchorByChr(SpeciesName query_name, SeqPro::ManagerVariant& query_seqpro_manager, ClusterBySQR_SparsePtr cluster, bool is_first);

    void filterAnchorByDP(AnchorBySQR_SparsePtr anchor_map,uint_t ref_chr_cnt, uint_t qry_chr_cnt);

    void constructGraphByDP(SpeciesName query_name, SeqPro::ManagerVariant& query_seqpro_manager, AnchorBySQR_SparsePtr anchor_ptr, RaMesh::RaMeshMultiGenomeGraph& graph);

};


#endif


