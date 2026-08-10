#pragma once

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <functional>
#include <cctype>
#include <optional>
#include "halDefs.h"
#include <unordered_map>
#include <set>

#include <spdlog/spdlog.h>
// forward declare hal types where possible to reduce direct includes in header
// keep heavy HAL headers in .cpp to avoid header bloat

#include "ramesh.h"
#include "SeqPro.h"
#include "softmask_index.h"
#include "data_process.h"  // 使用现有的NewickParser

// 前向声明
namespace RaMesh {
    class RaMeshMultiGenomeGraph;
    class Block;
    struct Segment;
    using BlockPtr = std::shared_ptr<Block>;
    using SegPtr = std::shared_ptr<Segment>;
    using SequenceId = SeqPro::SequenceId;
}

namespace RaMesh {

namespace hal_converter {

    // ========================================
    // 祖先节点信息结构
    // ========================================

    struct AncestorNode {
        std::string node_name;                              // 祖先节点名称
        std::vector<std::string> descendant_leaves;         // 该祖先下的所有叶节点
        std::map<std::string, std::string> reconstructed_sequences;  // chr_name -> sequence
        std::map<std::string, hal_size_t> sequence_lengths; // chr_name -> length

        // 树结构信息
        int tree_depth;                                     // 在树中的深度
        std::string parent_name;                            // 父节点名称
        std::vector<std::string> children_names;            // 子节点名称（所有后代）
        std::vector<std::string> direct_children_names;     // 直接子节点名称
        bool is_generated_root;                             // 是否为自动生成的根节点
        double branch_length;                               // 到父节点的分支长度

        AncestorNode() : tree_depth(0), is_generated_root(false), branch_length(0.0) {}
    };

    // ========================================
    // 系统发育树解析和处理
    // ========================================

    /**
     * 解析Newick树并提取祖先节点信息
     * @param newick_tree Newick格式的系统发育树字符串
     * @param seqpro_managers 序列管理器映射，用于验证叶节点名称
     * @return 祖先节点列表
     */
    std::vector<AncestorNode> parsePhylogeneticTree(
        const std::string& newick_tree,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
        const std::string& root_name);

    /**
     * 检查树是否有显式根节点，如果没有则自动生成
     * @param parser 已解析的Newick树
     * @param seqpro_managers 序列管理器映射
     * @return 是否添加了新的根节点
     */
    bool ensureRootNode(NewickParser& parser,
                       const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
                       const std::string& root_name);

    // ========================================
    // 祖先序列重建规划
    // ========================================

    /**
     * 规划祖先序列重建的顺序和参考叶子选择
     * @param ancestor_nodes 祖先节点列表
     * @param parser 系统发育树解析器
     * @return 重建计划：ancestor_name -> reference_leaf_name 的映射
     */
    std::vector<std::pair<std::string, std::string>> planAncestorReconstruction(
        const std::vector<AncestorNode>& ancestor_nodes,
        const NewickParser& parser);

    /**
     * 为指定祖先节点找到系统发育距离最近的叶子节点
     * @param ancestor 祖先节点
     * @param parser 系统发育树解析器
     * @return 最近叶子节点的名称，如果找不到则返回空字符串
     */
    std::string findClosestLeafForAncestor(
        const AncestorNode& ancestor,
        const NewickParser& parser);

    // ========================================
    // 祖先序列重建数据结构
    // ========================================

    /**
     * 祖先segment信息，用于重建序列
     */
    struct AncestorSegmentInfo {
        uint_t start;
        uint_t length;
        SequenceId chr_id;              // 使用ID而不是字符串，节省内存
        BlockPtr source_block;          // 来源block
        bool is_from_ref;               // 是否直接来自参考叶子

        // 间隙信息（用于序列重建时添加N）
        bool need_gap_before;           // 前面是否需要添加N
        bool need_gap_after;            // 后面是否需要添加N
    };

    /**
     * 祖先重建数据
     */
    struct AncestorReconstructionData {
        std::vector<AncestorSegmentInfo> segments;
        std::set<BlockPtr> processed_blocks;        // 避免重复添加
        std::string reference_leaf;                 // 参考叶子名称

        // 染色体ID映射函数
        std::function<SequenceId(const std::string&)> getSequenceId;
    };

    /**
     * 执行祖先序列重建的第二阶段
     * @param reconstruction_plan 第一阶段生成的重建计划
     * @param ancestor_nodes 祖先节点列表
     * @param seqpro_managers 序列管理器映射
     * @param graph 多基因组图
     * @return 祖先重建数据映射
     */
    std::map<std::string, AncestorReconstructionData> reconstructAncestorSequences(
        const std::vector<std::pair<std::string, std::string>>& reconstruction_plan,
        const std::vector<AncestorNode>& ancestor_nodes,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
        RaMeshMultiGenomeGraph& graph);

    // ========================================
    // 祖先序列构建
    // ========================================

    /**
     * 从segment信息中提取DNA序列
     * @param segment 祖先segment信息
     * @param seqpro_managers 序列管理器映射
     * @param reference_leaf 参考叶子名称
     * @return 提取的DNA序列
     */
    std::string extractSegmentDNA(const AncestorSegmentInfo& segment,
                                 const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
                                 const std::string& reference_leaf);

    /**
     * 根据重建数据构建完整的祖先序列
     * @param data 祖先重建数据
     * @param seqpro_managers 序列管理器映射
     * @return 构建的完整序列
     */
    std::string buildAncestorSequence(const AncestorReconstructionData& data,
                                     const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers);

    /**
     * 为所有祖先构建序列（使用投票法）
     * @param ancestor_reconstruction_data 所有祖先的重建数据
     * @param ancestor_nodes 祖先节点信息
     * @param seqpro_managers 序列管理器映射
     * @param reconstruction_plan 重建计划，确保正确的重建顺序
     * @param alignment 可选的HAL alignment对象，如果提供则直接存储到HAL中
     * @return 祖先序列映射 (ancestor_name -> chr_name -> chr_sequence)
     */
    std::map<std::string, std::map<std::string, std::string>> buildAllAncestorSequencesByVoting(
        const std::map<std::string, AncestorReconstructionData>& ancestor_reconstruction_data,
        const std::vector<AncestorNode>& ancestor_nodes,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
        const std::vector<std::pair<std::string, std::string>>& reconstruction_plan,
        const SoftMask::IndexMap& softmask_indexes,
        hal::AlignmentPtr alignment = nullptr);

    // ========================================
    // 投票法祖先序列重建
    // ========================================

    /**
     * 从block中提取所有物种的序列和CIGAR信息
     * @param block 比对块
     * @param seqpro_managers 序列管理器映射
     * @return {序列映射, CIGAR映射}
     */
    std::pair<std::unordered_map<std::string, std::string>, std::unordered_map<ChrName, Cigar_t>>
    extractSequencesAndCigarsFromBlock(
        BlockPtr block,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
        const SoftMask::IndexMap* softmask_indexes = nullptr);

    /**
     * 按列投票重建祖先序列
     * @param aligned_sequences 对齐后的序列映射
     * @param ancestor 祖先节点信息
     * @return 重建的祖先序列（不含gap）
     */
    std::string voteForAncestorSequence(
        const std::unordered_map<std::string, std::string>& aligned_sequences,
        const AncestorNode& ancestor);

    /**
     * 使用投票法重建单个segment的序列
     * @param segment 祖先segment信息
     * @param ancestor 祖先节点信息
     * @param seqpro_managers 序列管理器映射
     * @return 重建的segment序列
     */
    std::string reconstructSegmentByVoting(
        const AncestorSegmentInfo& segment,
        const AncestorNode& ancestor,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
        const SoftMask::IndexMap& softmask_indexes);

    /**
     * 使用投票法构建完整的祖先序列，同时创建祖先segment并加入到对应的block中
     * @param data 祖先重建数据
     * @param ancestor 祖先节点信息
     * @param seqpro_managers 序列管理器映射
     * @param chr_name 染色体名称
     * @return 构建的完整序列
     */
    std::string buildAncestorSequenceByVoting(
        const AncestorReconstructionData& data,
        const AncestorNode& ancestor,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
        const SoftMask::IndexMap& softmask_indexes,
        const std::string& chr_name);

    /**
     * 从NewickParser提取祖先节点信息
     * @param parser 已解析的Newick树
     * @param seqpro_managers 序列管理器映射
     * @return 祖先节点列表
     */
    std::vector<AncestorNode> extractAncestorNodes(
        const NewickParser& parser,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
        const std::string& root_name);

    /**
     * 验证叶节点名称是否与提供的物种名称匹配
     * @param parser 已解析的Newick树
     * @param seqpro_managers 序列管理器映射
     * @return 验证结果和错误信息
     */
    std::pair<bool, std::string> validateLeafNames(
        const NewickParser& parser,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers);

    // ========================================
    // HAL基础结构创建
    // ========================================

    /**
     * 创建基本的基因组结构（仅叶节点）
     * @param alignment HAL alignment对象
     * @param seqpro_managers 序列管理器映射
     */
    void createBasicGenomeStructure(
        hal::AlignmentPtr alignment,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers);

    /**
     * 设置基因组的序列维度和DNA数据
     * @param alignment HAL alignment对象
     * @param seqpro_managers 序列管理器映射
     */
    void setupGenomeSequences(
        hal::AlignmentPtr alignment,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers);

    /**
     * 为所有叶基因组设置真实的序列维度与DNA
     * @param alignment HAL alignment对象
     * @param seqpro_managers 序列管理器映射（每个叶物种）
     */
    void setupLeafGenomesWithRealDNA(
        hal::AlignmentPtr alignment,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
        const SoftMask::IndexMap& softmask_indexes);





    /**
     * 基于系统发育树创建基因组骨架（根/中间祖先/叶）
     * 仅创建拓扑结构，不设置维度或DNA；后续阶段再 setDimensions/setString
     * @param alignment HAL alignment对象
     * @param ancestor_nodes 祖先节点列表（内部祖先信息）
     * @param parser NewickParser对象（用于解析叶节点父子与分支长度）
     */
    void createGenomesFromPhylogeny(
        hal::AlignmentPtr alignment,
        const std::vector<AncestorNode>& ancestor_nodes,
        const NewickParser& parser,
        const std::string& root_name);

    /**
     * 用正确的维度创建祖先基因组
     * @param alignment HAL alignment对象
     * @param ancestor_sequences 重建的祖先序列
     * @param ancestor_nodes 祖先节点列表
     * @param seqpro_managers 序列管理器映射
     */
    void createAncestorGenomesWithCorrectDimensions(
        hal::AlignmentPtr alignment,
        const std::map<std::string, std::string>& ancestor_sequences,
        const std::vector<AncestorNode>& ancestor_nodes,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers);

    /**
     * 从NewickParser重建Newick字符串
     * @param parser 已解析的Newick树
     * @return 重建的Newick字符串
     */
    std::string reconstructNewickFromParser(const NewickParser& parser);

    /**
     * 应用系统发育树结构到HAL alignment
     * @param alignment HAL alignment对象
     * @param parser 已解析的Newick树
     */
    void applyPhylogeneticTree(
        hal::AlignmentPtr alignment,
        const NewickParser& parser);

    // ========================================
    // 合并的blocks分析和HAL结构构建
    // ========================================

    /**
     * 当前block的segment映射信息（以父段为聚合单位，一个bottom可对应多个top）
     */
    struct CurrentBlockMapping {
        struct ChildInfo {
            std::string child_genome;
            std::string child_chr_name;
            hal_size_t child_start;
            hal_size_t child_length;
            bool is_reversed;
        };

        std::string parent_genome;
        std::string parent_chr_name;  // parent基因组的染色体名称
        hal_size_t parent_start;
        hal_size_t parent_length;

        std::vector<ChildInfo> children;  // 与该父段相连的所有直系子top段
    };

    /**
     * HAL segment索引管理器
     */
    struct SegmentIndexManager {
        // genome_name -> chr_name -> current_index
        std::map<std::string, std::map<std::string, hal_index_t>> top_segment_indices;
        std::map<std::string, std::map<std::string, hal_index_t>> bottom_segment_indices;

        // genome_name -> chr_name -> total_count
        std::map<std::string, std::map<std::string, hal_size_t>> top_segment_counts;
        std::map<std::string, std::map<std::string, hal_size_t>> bottom_segment_counts;

        // 已更新维度的基因组记录
        std::set<std::string> dimensions_updated_genomes;
    };

    /**
     * 流式分析blocks并直接构建HAL结构（合并第3、4、5阶段）
     * @param blocks 所有的alignment blocks
     * @param ancestor_nodes 祖先节点信息
     * @param alignment HAL alignment对象
     * @param ancestor_data 祖先重建数据（用于gap-aware映射）
     * @param ancestor_sequences 祖先序列映射 (ancestor_name -> chr_name -> chr_sequence)
     * @param seqpro_managers 序列管理器（用于提取序列）
     */
    void analyzeBlocksAndBuildHalStructure(
        const std::vector<std::weak_ptr<Block>>& blocks,
        const std::vector<AncestorNode>& ancestor_nodes,
        hal::AlignmentPtr alignment,
        const std::map<std::string, AncestorReconstructionData>& ancestor_data,
        const std::map<std::string, std::map<std::string, std::string>>& ancestor_sequences,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers);

    /**
     * 分析当前block中的segment关系（按你的算法思路）
     * @param block 当前block
     * @param ancestor_nodes 祖先节点信息
     * @return 当前block的映射关系
     */
    std::vector<CurrentBlockMapping> analyzeCurrentBlock(
        BlockPtr block,
        const std::vector<AncestorNode>& ancestor_nodes);

    /**
     * 使用gap-aware方法分析当前block的映射关系
     * @param block 当前block
     * @param ancestor_nodes 祖先节点信息
     * @param ancestor_data 祖先重建数据
     * @param ancestor_sequences 祖先序列映射 (ancestor_name -> chr_name -> chr_sequence)
     * @param seqpro_managers 序列管理器
     * @return 拆分后的映射关系列表
     */
    std::vector<CurrentBlockMapping> analyzeBlockWithGapHandling(
        BlockPtr block,
        const std::vector<AncestorNode>& ancestor_nodes,
        const std::map<std::string, AncestorReconstructionData>& ancestor_data,
        const std::map<std::string, std::map<std::string, std::string>>& ancestor_sequences,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers);

    /**
     * 连续段信息结构体
     */
    struct ContinuousSegmentInfo {
        hal_index_t start;
        hal_size_t length;
        bool is_gap;                    // 是否是填充的间隙段
        hal_index_t array_index = -1;   // HAL 数组索引（写入后记录）

        // 原始映射信息（仅对非间隙段有效）
        std::optional<std::pair<BlockPtr, size_t>> mapping_info;

        bool operator<(const ContinuousSegmentInfo& other) const {
            return start < other.start;
        }
    };

    // 连续段映射类型
    using ContinuousSegmentMap = std::map<std::pair<std::string, std::string>, std::vector<ContinuousSegmentInfo>>;

    /**
     * 创建连续段并更新HAL维度（整合方案）
     * @param alignment HAL alignment对象
     * @param refined_mappings 精炼的映射关系
     * @param seqpro_managers 序列管理器
     * @return 连续的bottom和top段映射
     */
    std::pair<ContinuousSegmentMap, ContinuousSegmentMap> checkAndUpdateHalDimensions(
        hal::AlignmentPtr alignment,
        const std::map<BlockPtr, std::vector<CurrentBlockMapping>>& refined_mappings,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers);

    /**
     */
    inline void createSegmentsForCurrentBlock(
        hal::AlignmentPtr,
        const std::vector<CurrentBlockMapping>&,
        SegmentIndexManager&) {}

    /**
     * 从缓存的映射关系创建段并建立链接
     * @param alignment HAL alignment对象
     * @param refined_mappings 缓存的映射关系
     * @param index_manager segment索引管理器
     */
    void createSegmentsAndLinksFromMappings(
        hal::AlignmentPtr alignment,
        const std::map<BlockPtr, std::vector<CurrentBlockMapping>>& refined_mappings,
        SegmentIndexManager& index_manager);

    /**
     * 为单个映射关系创建段并建立链接
     * @param alignment HAL alignment对象
     * @param mapping 单个映射关系
     * @param current_top_indices 当前top段索引
     * @param current_bottom_indices 当前bottom段索引
     */
    void createSegmentAndLinkForMapping(
        hal::AlignmentPtr alignment,
        const CurrentBlockMapping& mapping,
        std::map<std::string, std::map<std::string, hal_index_t>>& current_top_indices,
        std::map<std::string, std::map<std::string, hal_index_t>>& current_bottom_indices);

    /**
     * 辅助函数：检查是否为叶节点
     */
    bool isLeafGenome(const std::string& genome_name, const std::vector<AncestorNode>& ancestor_nodes);

    /**
     * 辅助函数：检查是否为根节点
     */
    bool isRootGenome(const std::string& genome_name, const std::vector<AncestorNode>& ancestor_nodes);

    /**
     * 辅助函数：从block中获取指定物种的segment
     */
    SegPtr getSegmentFromBlock(BlockPtr block, const std::string& species_name);



    // ========================================
    // 验证和工具函数
    // ========================================

    /**
     * 验证HAL文件的完整性
     * @param alignment HAL alignment对象
     * @param hal_path HAL文件路径
     */
    void validateHalFile(hal::AlignmentPtr alignment, const FilePath& hal_path);

    /**
     * 生成默认的根节点名称
     * @param leaf_names 叶节点名称列表
     * @return 生成的根节点名称
     */
    std::string generateRootName(const std::vector<std::string>& leaf_names);

} // namespace hal_converter

} // namespace RaMesh
