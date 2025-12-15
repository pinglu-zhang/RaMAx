#include "hal_converter.h"
#include "align.h"
#include <spdlog/spdlog.h>
#include "../submodule/hal/api/inc/hal.h"
#include "../submodule/hal/api/inc/halCommon.h"
#include "../../include/threadpool.h"
#include <mutex>
#include <chrono>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>

namespace std {
    template<>
    struct hash<std::tuple<std::string, std::string, hal_index_t, hal_size_t>> {
        size_t operator()(const std::tuple<std::string, std::string, hal_index_t, hal_size_t>& t) const {
            auto hash_combine = [](size_t& seed, const auto& v) {
                seed ^= std::hash<std::decay_t<decltype(v)>>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            };
            size_t result = 0;
            std::apply([&](const auto&... args) {
                (hash_combine(result, args), ...);
            }, t);
            return result;
        }
    };

    template<>
    struct hash<std::pair<std::string, std::string>> {
        size_t operator()(const std::pair<std::string, std::string>& p) const {
            auto h1 = std::hash<std::string>{}(p.first);
            auto h2 = std::hash<std::string>{}(p.second);
            // A simple way to combine hashes
            return h1 ^ (h2 << 1);
        }
    };
}

namespace RaMesh {
namespace hal_converter {

    // 全局 HAL 写锁：保护所有对 hal::Alignment/Genome 的写操作（HDF5 非线程安全）
    static std::mutex g_hal_write_mutex;

    // ========================================
    // 系统发育树解析和处理
    // ========================================

    std::vector<AncestorNode> parsePhylogeneticTree(
        const std::string& newick_tree,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
        const std::string& root_name) {

        spdlog::info("Parsing phylogenetic tree to identify ancestor nodes...");

        std::vector<AncestorNode> ancestor_nodes;

        if (newick_tree.empty()) {
            throw std::runtime_error("Newick tree is required but empty; simplify assumptions: please provide a valid tree");
        }

        // 解析Newick树
        try {
            NewickParser parser(newick_tree);

            // 验证叶节点名称
            auto [is_valid, error_msg] = validateLeafNames(parser, seqpro_managers);
            if (!is_valid) {
                spdlog::warn("Leaf name validation failed: {}", error_msg);
                spdlog::warn("Proceeding with available species...");
            }

            // 检查并确保有根节点
            NewickParser mutable_parser = parser;  // 创建可修改的副本
            bool added_root = ensureRootNode(mutable_parser, seqpro_managers, root_name);
            if (added_root) {
                spdlog::info("Added artificial root node '{}'", root_name);
            }

            // 提取祖先节点信息
            ancestor_nodes = extractAncestorNodes(mutable_parser, seqpro_managers, root_name);

            spdlog::info("Found {} ancestor nodes from tree", ancestor_nodes.size());
            for (const auto& ancestor : ancestor_nodes) {
                spdlog::debug("  Ancestor '{}' with {} descendants (depth: {})",
                             ancestor.node_name, ancestor.descendant_leaves.size(), ancestor.tree_depth);
            }

        } catch (const std::exception& e) {
            spdlog::error("Failed to parse Newick tree: {}", e.what());
            throw;
        }

        return ancestor_nodes;
    }

    bool ensureRootNode(NewickParser& parser,
                       const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
                       const std::string& root_name) {

        const auto& nodes = parser.getNodes();
        if (nodes.empty()) {
            return false;
        }

        // 找到根节点（father == -1的节点）
        int root_id = -1;
        for (const auto& node : nodes) {
            if (node.father == -1) {
                root_id = node.id;
                break;
            }
        }

        if (root_id == -1) {
            spdlog::error("No root node found in parsed tree");
            return false;
        }

        const auto& root_node = nodes[root_id];

        // 检查根节点是否是叶节点
        if (root_node.isLeaf) {
            spdlog::info("Root is a leaf node, adding artificial root '{}'", root_name);

            // 创建新的根节点
            NewickTreeNode new_root;
            new_root.id = parser.currentIndex_++;
            new_root.name = root_name.empty() ? std::string("ancestor") : root_name;
            new_root.father = -1;
            new_root.isLeaf = false;
            new_root.branchLength = 0.0;
            new_root.leftChild = root_id;
            new_root.rightChild = -1;  // 只有一个子节点

            // 更新原根节点
            parser.nodes_[root_id].father = new_root.id;
            parser.nodes_[root_id].branchLength = 1.0;  // 设置默认分支长度

            // 添加新根节点
            parser.nodes_.push_back(new_root);

            return true;
        }

        // 检查根节点是否有名称，如果没有则设置默认名称
        if (root_node.name.empty()) {
            parser.nodes_[root_id].name = root_name.empty() ? std::string("ancestor") : root_name;
            spdlog::info("Set root node name to '{}'", parser.nodes_[root_id].name);
            return true;
        }

        return false;
    }

    std::vector<AncestorNode> extractAncestorNodes(
        const NewickParser& parser,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
        const std::string& root_name) {

        std::vector<AncestorNode> ancestor_nodes;
        const auto& nodes = parser.getNodes();

        if (nodes.empty()) {
            return ancestor_nodes;
        }

        // 为每个内部节点创建AncestorNode
        for (const auto& node : nodes) {
            if (!node.isLeaf) {  // 只处理内部节点
                AncestorNode ancestor;
                ancestor.node_name = node.name.empty() ? ("internal_" + std::to_string(node.id)) : node.name;
                ancestor.branch_length = node.branchLength;
                ancestor.is_generated_root = ((node.name == root_name || (root_name.empty() && node.name == "ancestor")) && node.father == -1);

                // 计算树深度
                ancestor.tree_depth = 0;
                int current_id = node.id;
                while (current_id != -1) {
                    bool found = false;
                    for (const auto& n : nodes) {
                        if (n.id == current_id && n.father != -1) {
                            ancestor.tree_depth++;
                            current_id = n.father;
                            found = true;
                            break;
                        }
                    }
                    if (!found) break;
                }

                // 设置父节点名称
                if (node.father != -1) {
                    for (const auto& n : nodes) {
                        if (n.id == node.father) {
                            ancestor.parent_name = n.name.empty() ? ("internal_" + std::to_string(n.id)) : n.name;
                            break;
                        }
                    }
                }

                // 收集所有后代叶节点
                std::function<void(int)> collectLeaves = [&](int node_id) {
                    for (const auto& n : nodes) {
                        if (n.father == node_id) {
                            if (n.isLeaf) {
                                // 只添加在seqpro_managers中存在的叶节点
                                if (seqpro_managers.find(n.name) != seqpro_managers.end()) {
                                    ancestor.descendant_leaves.push_back(n.name);
                                }
                            } else {
                                collectLeaves(n.id);
                            }

                            // 添加到直接子节点列表
                            std::string child_name = n.name.empty() ? ("internal_" + std::to_string(n.id)) : n.name;
                            ancestor.children_names.push_back(child_name);
                        }
                    }
                };

                // 正确填充"直接子节点"：仅 father == node.id 的一级孩子
                ancestor.direct_children_names.clear();
                for (const auto& n : nodes) {
                    if (n.father == node.id) {
                        std::string child_name = n.name.empty() ? ("internal_" + std::to_string(n.id)) : n.name;
                        ancestor.direct_children_names.push_back(child_name);
                    }
                }

                collectLeaves(node.id);

                // 只有当祖先节点有后代叶节点时才添加
                if (!ancestor.descendant_leaves.empty()) {
                    ancestor_nodes.push_back(ancestor);
                }
            }
        }

        return ancestor_nodes;
    }

    // ========================================
    // 祖先序列重建规划
    // ========================================

    std::vector<std::pair<std::string, std::string>> planAncestorReconstruction(
        const std::vector<AncestorNode>& ancestor_nodes,
        const NewickParser& parser) {

        spdlog::debug("Planning ancestor sequence reconstruction...");

        // 第一步：按树深度排序祖先节点（深度大的先处理，即叶子优先）
        std::vector<AncestorNode> sorted_ancestors = ancestor_nodes;
        std::sort(sorted_ancestors.begin(), sorted_ancestors.end(),
            [](const AncestorNode& a, const AncestorNode& b) {
                return a.tree_depth > b.tree_depth; // 深度大的先处理
            });

        spdlog::debug("Sorted {} ancestors by tree depth (deepest first)", sorted_ancestors.size());
        for (const auto& ancestor : sorted_ancestors) {
            spdlog::debug("  Ancestor '{}' at depth {}", ancestor.node_name, ancestor.tree_depth);
        }

        // 第二步：为每个祖先确定参考叶子
        std::vector<std::pair<std::string, std::string>> reconstruction_plan;

        for (const auto& ancestor : sorted_ancestors) {
            std::string reference_leaf = findClosestLeafForAncestor(ancestor, parser);

            if (!reference_leaf.empty()) {
                reconstruction_plan.emplace_back(ancestor.node_name, reference_leaf);
                spdlog::debug("  Ancestor '{}' -> Reference leaf '{}'",
                             ancestor.node_name, reference_leaf);
            } else {
                spdlog::warn("  Could not find reference leaf for ancestor '{}'", ancestor.node_name);
            }
        }

        spdlog::info("Reconstruction plan created for {} ancestors", reconstruction_plan.size());
        return reconstruction_plan;
    }

    std::string findClosestLeafForAncestor(
        const AncestorNode& ancestor,
        const NewickParser& parser) {

        spdlog::debug("Finding closest leaf for ancestor '{}'", ancestor.node_name);

        // 如果祖先没有后代叶子，返回空
        if (ancestor.descendant_leaves.empty()) {
            spdlog::debug("  No descendant leaves found for ancestor '{}'", ancestor.node_name);
            return "";
        }

        // 策略1：简单版本 - 选择第一个后代叶子作为参考
        // 在更复杂的实现中，可以基于分支长度计算真正的系统发育距离
        std::string closest_leaf = ancestor.descendant_leaves[0];

        spdlog::debug("  Selected '{}' as reference leaf for ancestor '{}' (simple strategy)",
                     closest_leaf, ancestor.node_name);

        // 策略2：基于分支长度的选择（如果需要更精确的选择）
        // 这里可以扩展为计算从祖先到每个后代叶子的累积分支长度
        // 并选择距离最短的叶子

        // TODO: 实现基于分支长度的精确距离计算
        // double min_distance = std::numeric_limits<double>::max();
        // std::string best_leaf;
        // for (const auto& leaf : ancestor.descendant_leaves) {
        //     double distance = calculatePhylogeneticDistance(ancestor.node_name, leaf, parser);
        //     if (distance < min_distance) {
        //         min_distance = distance;
        //         best_leaf = leaf;
        //     }
        // }

        return closest_leaf;
    }

    // ========================================
    // 祖先序列重建实现
    // ========================================

    /**
     * 获取染色体ID的辅助函数
     */
    SequenceId getChrId(const std::string& chr_name, const std::string& species_name,
                       const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {
        auto it = seqpro_managers.find(species_name);
        if (it == seqpro_managers.end()) {
            return SeqPro::SequenceIndex::INVALID_ID;
        }

        return std::visit([&chr_name](const auto& mgr) -> SequenceId {
            return mgr->getSequenceId(chr_name);
        }, *it->second);
    }

    /**
     * 根据名称查找祖先节点
     */
    const AncestorNode* findAncestorByName(const std::string& ancestor_name,
                                          const std::vector<AncestorNode>& ancestor_nodes) {
        for (const auto& ancestor : ancestor_nodes) {
            if (ancestor.node_name == ancestor_name) {
                return &ancestor;
            }
        }
        return nullptr;
    }

    /**
     * 检查block是否包含指定物种
     */
    bool blockContainsSpecies(BlockPtr block, const std::string& species_name) {
        if (!block) return false;

        std::shared_lock lock(block->rw);
        for (const auto& [species_chr, segment] : block->anchors) {
            if (species_chr.first == species_name) {
                return true;
            }
        }
        return false;
    }

    /**
     * 在block中查找指定物种的segment
     */
    SegPtr findSegmentInBlock(BlockPtr block, const std::string& species_name) {
        if (!block) return nullptr;

        std::shared_lock lock(block->rw);
        for (const auto& [species_chr, segment] : block->anchors) {
            if (species_chr.first == species_name) {
                return segment;
            }
        }
        return nullptr;
    }

    /**
     * 获取当前block中属于该祖先的其他物种（直接子节点）
     */
    std::set<std::string> getAncestorSpeciesInBlock(BlockPtr block, const AncestorNode& ancestor) {
        std::set<std::string> result;
        if (!block) return result;

        std::shared_lock lock(block->rw);

        for (const auto& [species_chr, segment] : block->anchors) {
            // 检查是否是祖先的直接子节点
            if (std::find(ancestor.children_names.begin(), ancestor.children_names.end(),
                         species_chr.first) != ancestor.children_names.end()) {
                result.insert(species_chr.first);
            }
        }
        return result;
    }

    /**
     * 创建来自参考segment的祖先segment信息
     */
    AncestorSegmentInfo createFromRefSegment(SegPtr segment, SequenceId chr_id, bool is_from_ref) {
        return {
            .start = segment->start,
            .length = segment->length,
            .chr_id = chr_id,
            .source_block = segment->parent_block,
            .is_from_ref = is_from_ref,
            .need_gap_before = false,
            .need_gap_after = false
        };
    }

    /**
     * 检查并设置间隙信息
     */
    void checkAndSetGapInfo(AncestorSegmentInfo& segment_info,
                           SegPtr prev_segment, SegPtr current_segment,
                           bool is_gap_before) {

        uint_t prev_end = prev_segment->start + prev_segment->length;
        uint_t current_start = current_segment->start;

        if (current_start > prev_end) {
            // 有间隙，需要添加N
            if (is_gap_before) {
                segment_info.need_gap_before = true;
            } else {
                segment_info.need_gap_after = true;
            }
        }
    }

    /**
     * 为祖先填补缺失区域
     */
    void fillGapsForAncestor(SegPtr ref_segment, const AncestorNode& ancestor,
                            AncestorReconstructionData& data, SequenceId chr_id) {

        // 获取当前block中属于该祖先的其他物种
        auto other_species = getAncestorSpeciesInBlock(ref_segment->parent_block, ancestor);

        SegPtr ref_next = ref_segment->primary_path.next.load(std::memory_order_acquire);
        if (!ref_next || ref_next->isTail()) return;

        // 检查其他物种的next segment路径
        for (const auto& species : other_species) {
            SegPtr species_segment = findSegmentInBlock(ref_segment->parent_block, species);
            if (!species_segment) continue;

            SegPtr species_next = species_segment->primary_path.next.load(std::memory_order_acquire);

            // 沿着species的next链遍历，直到找到包含ref的block
            while (species_next && !species_next->isTail()) {
                if (!blockContainsSpecies(species_next->parent_block, data.reference_leaf)) {
                    // 这是一个gap block，需要添加
                    if (data.processed_blocks.find(species_next->parent_block) == data.processed_blocks.end()) {

                        AncestorSegmentInfo gap_segment = createFromRefSegment(species_next, chr_id, false);

                        // 检查gap segment的间隙信息
                        SegPtr gap_prev = species_next->primary_path.prev.load(std::memory_order_acquire);
                        SegPtr gap_next = species_next->primary_path.next.load(std::memory_order_acquire);

                        if (gap_prev && !gap_prev->isHead()) {
                            checkAndSetGapInfo(gap_segment, gap_prev, species_next, true);
                        }
                        if (gap_next && !gap_next->isTail()) {
                            checkAndSetGapInfo(gap_segment, species_next, gap_next, false);
                        }

                        data.segments.push_back(gap_segment);
                        data.processed_blocks.insert(species_next->parent_block);

                        // spdlog::debug("    Added gap segment from species '{}' at {}:{}",
                        //              species, gap_segment.start, gap_segment.start + gap_segment.length);
                    }
                    species_next = species_next->primary_path.next.load(std::memory_order_acquire);
                } else {
                    break; // 找到了包含ref的block，停止
                }
            }
        }
    }

    /**
     * 重建单个染色体的祖先序列
     */
    void reconstructAncestorChromosome(
        const std::string& ancestor_name,
        const std::string& ref_leaf,
        const std::string& chr_name,
        const GenomeEnd& ref_genome_end,
        const AncestorNode& ancestor,
        AncestorReconstructionData& data) {

        SequenceId chr_id = data.getSequenceId(chr_name);
        if (chr_id == SeqPro::SequenceIndex::INVALID_ID) {
            spdlog::warn("Cannot get chromosome ID for {} in ancestor {}", chr_name, ancestor_name);
            return;
        }

        // spdlog::debug("  Reconstructing chromosome '{}' for ancestor '{}'", chr_name, ancestor_name);

        std::shared_lock end_lock(ref_genome_end.rw);

        // 遍历参考叶子的segment链表
        SegPtr current = ref_genome_end.head->primary_path.next.load(std::memory_order_acquire);
        SegPtr prev_segment = nullptr;
        size_t segment_count = 0;

        while (current && !current->isTail()) {
            if (current->isSegment()) {
                // 1. 添加当前ref segment
                AncestorSegmentInfo segment_info = createFromRefSegment(current, chr_id, true);

                // 2. 检查与前一个segment的间隙
                if (prev_segment && !prev_segment->isHead()) {
                    checkAndSetGapInfo(segment_info, prev_segment, current, true); // gap_before
                }

                data.segments.push_back(segment_info);
                data.processed_blocks.insert(current->parent_block);
                segment_count++;

                // spdlog::debug("    Added ref segment at {}:{} (gap_before: {})",
                //              segment_info.start, segment_info.start + segment_info.length,
                //              segment_info.need_gap_before);

                // 3. 检查并填补缺失区域
                fillGapsForAncestor(current, ancestor, data, chr_id);

                prev_segment = current;
            }
            current = current->primary_path.next.load(std::memory_order_acquire);
        }

        // 检查最后一个segment的gap_after
        if (!data.segments.empty() && prev_segment) {
            auto& last_segment = data.segments.back();
            SegPtr next_segment = prev_segment->primary_path.next.load(std::memory_order_acquire);
            if (next_segment && !next_segment->isTail()) {
                checkAndSetGapInfo(last_segment, prev_segment, next_segment, false); // gap_after
                spdlog::debug("    Last segment gap_after: {}", last_segment.need_gap_after);
            }
        }

        // spdlog::debug("  Completed chromosome '{}': {} segments", chr_name, segment_count);
    }

    /**
     * 重建单个祖先的序列
     */
    void reconstructSingleAncestor(
        const std::string& ancestor_name,
        const std::string& ref_leaf,
        const std::vector<AncestorNode>& ancestor_nodes,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
        RaMeshMultiGenomeGraph& graph,
        AncestorReconstructionData& data) {

        spdlog::info("Reconstructing ancestor '{}' using reference leaf '{}'", ancestor_name, ref_leaf);

        const AncestorNode* ancestor = findAncestorByName(ancestor_name, ancestor_nodes);
        if (!ancestor) {
            spdlog::error("Cannot find ancestor node: {}", ancestor_name);
            return;
        }

        // 遍历参考叶子的所有染色体
        auto& ref_genome = graph.species_graphs[ref_leaf];
        std::shared_lock ref_lock(ref_genome.rw);

        size_t total_segments = 0;
        for (const auto& [chr_name, genome_end] : ref_genome.chr2end) {
            size_t segments_before = data.segments.size();

            reconstructAncestorChromosome(ancestor_name, ref_leaf, chr_name,
                                        genome_end, *ancestor, data);

            size_t segments_added = data.segments.size() - segments_before;
            total_segments += segments_added;

            // spdlog::debug("  Chromosome '{}': added {} segments", chr_name, segments_added);
        }

        spdlog::info("Completed ancestor '{}': {} total segments across {} chromosomes",
                    ancestor_name, total_segments, ref_genome.chr2end.size());
    }

    /**
     * 执行祖先序列重建的第二阶段
     */
    std::map<std::string, AncestorReconstructionData> reconstructAncestorSequences(
        const std::vector<std::pair<std::string, std::string>>& reconstruction_plan,
        const std::vector<AncestorNode>& ancestor_nodes,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
        RaMeshMultiGenomeGraph& graph) {

        spdlog::info("Starting ancestor sequence reconstruction for {} ancestors", reconstruction_plan.size());

        std::map<std::string, AncestorReconstructionData> ancestor_reconstruction_data;

        // 初始化每个祖先的重建数据
        for (const auto& [ancestor_name, ref_leaf] : reconstruction_plan) {
            auto& data = ancestor_reconstruction_data[ancestor_name];
            data.reference_leaf = ref_leaf;

            // 设置染色体ID获取函数
            data.getSequenceId = [&ref_leaf, &seqpro_managers](const std::string& chr_name) {
                return getChrId(chr_name, ref_leaf, seqpro_managers);
            };

            spdlog::debug("Initialized reconstruction data for ancestor '{}' with reference '{}'",
                         ancestor_name, ref_leaf);
        }

        // 并行处理每个祖先重建（计算密集，安全并行）
        ThreadPool pool(std::max(1u, std::thread::hardware_concurrency()));
        for (const auto& [ancestor_name, ref_leaf] : reconstruction_plan) {
            auto data_ptr = &ancestor_reconstruction_data[ancestor_name];
            pool.enqueue([&, ancestor_name, ref_leaf, data_ptr]() {
            reconstructSingleAncestor(ancestor_name, ref_leaf, ancestor_nodes,
                                          seqpro_managers, graph, *data_ptr);
            });
        }
        pool.waitAllTasksDone();

        // 输出统计信息
        spdlog::info("Ancestor sequence reconstruction completed:");
        for (const auto& [ancestor_name, data] : ancestor_reconstruction_data) {
            spdlog::info("  Ancestor '{}': {} segments, {} processed blocks",
                        ancestor_name, data.segments.size(), data.processed_blocks.size());
        }

        return ancestor_reconstruction_data;
    }

    // ========================================
    // 祖先序列构建实现
    // ========================================

    std::string extractSegmentDNA(const AncestorSegmentInfo& segment,
                                 const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
                                 const std::string& reference_leaf) {

        auto it = seqpro_managers.find(reference_leaf);
        if (it == seqpro_managers.end()) {
            throw std::runtime_error("Reference leaf not found: " + reference_leaf);
        }

        return std::visit([&](const auto& mgr) -> std::string {
            using PtrType = std::decay_t<decltype(mgr)>;
            if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
                return mgr->getSubSequence(segment.chr_id, segment.start, segment.length);
            } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                return mgr->getSubSequence(segment.chr_id, segment.start, segment.length);
            } else {
                throw std::runtime_error("Unhandled manager type in variant.");
            }
        }, *it->second);
    }

    std::string buildAncestorSequence(const AncestorReconstructionData& data,
                                     const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {

        // spdlog::debug("Building ancestor sequence from {} segments", data.segments.size());

        std::string full_sequence;
        size_t total_segments = data.segments.size();
        size_t gaps_added = 0;

        for (size_t i = 0; i < total_segments; ++i) {
            const auto& segment = data.segments[i];

            // 处理gap_before（避免重复N）
            if (segment.need_gap_before) {
                if (full_sequence.back() != 'N') {
                    full_sequence += 'N';
                    gaps_added++;
                }
            }

            // 提取并添加segment的DNA序列
            try {
                std::string segment_dna = extractSegmentDNA(segment, seqpro_managers, data.reference_leaf);
                full_sequence += segment_dna;

                // spdlog::debug("  Added segment {} ({}:{}) length={}, is_from_ref={}",
                //              i, segment.start, segment.start + segment.length,
                //              segment.length, segment.is_from_ref);
            } catch (const std::exception& e) {
                spdlog::error("Failed to extract DNA for segment {}: {}", i, e.what());
                throw;
            }

            // 处理gap_after（避免重复N）
            if (segment.need_gap_after) {
                full_sequence += 'N';
                gaps_added++;
            }
        }

        spdlog::info("Built ancestor sequence: {} segments, {} gaps, {} total length",
                    total_segments, gaps_added, full_sequence.length());

        return full_sequence;
    }

    std::map<std::string, std::string> buildAllAncestorSequences(
        const std::map<std::string, AncestorReconstructionData>& ancestor_reconstruction_data,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {

        spdlog::info("Building sequences for {} ancestors", ancestor_reconstruction_data.size());

        std::map<std::string, std::string> ancestor_sequences;

        for (const auto& [ancestor_name, data] : ancestor_reconstruction_data) {
            spdlog::info("Building sequence for ancestor '{}'", ancestor_name);

            try {
                std::string sequence = buildAncestorSequence(data, seqpro_managers);
                ancestor_sequences[ancestor_name] = std::move(sequence);

                spdlog::info("Successfully built sequence for ancestor '{}': {} bp",
                            ancestor_name, ancestor_sequences[ancestor_name].length());
            } catch (const std::exception& e) {
                spdlog::error("Failed to build sequence for ancestor '{}': {}", ancestor_name, e.what());
                throw;
            }
        }

        // 输出统计信息
        spdlog::info("Ancestor sequence construction completed:");
        size_t total_length = 0;
        for (const auto& [ancestor_name, sequence] : ancestor_sequences) {
            spdlog::info("  Ancestor '{}': {} bp", ancestor_name, sequence.length());
            total_length += sequence.length();
        }
        spdlog::info("  Total ancestor sequence length: {} bp", total_length);

        return ancestor_sequences;
    }

    // ========================================
    // 投票法祖先序列重建实现
    // ========================================

    std::pair<std::unordered_map<std::string, std::string>, std::unordered_map<ChrName, Cigar_t>>
    extractSequencesAndCigarsFromBlock(
        BlockPtr block,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {

        std::unordered_map<std::string, std::string> sequences;
        std::unordered_map<ChrName, Cigar_t> cigars;

        if (!block) {
            return {sequences, cigars};
        }

        std::shared_lock lock(block->rw);

        // lambda函数：提取序列（复用ramesh_export.cpp中的逻辑）
        auto fetchSeq = [](const SeqPro::SharedManagerVariant& shared_mv,
            const ChrName& chr, Coord_t start, Coord_t length) -> std::string {
                return std::visit([&](const auto& mgr) -> std::string {
                    using PtrType = std::decay_t<decltype(mgr)>;
                    if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
                        auto chr_id = mgr->getSequenceId(chr);
                        return mgr->getSubSequence(chr_id, start, length);
                    } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                        auto chr_id = mgr->getSequenceId(chr);
                        return mgr->getOriginalManager().getSubSequence(chr_id, start, length);
                    } else {
                        throw std::runtime_error("Unhandled manager type in variant.");
                    }
                }, *shared_mv);
            };

        // 遍历block中的所有物种segment
        for (const auto& [species_chr, segment] : block->anchors) {
            const std::string& species_name = species_chr.first;
            const std::string& chr_name = species_chr.second;

            auto it = seqpro_managers.find(species_name);
            if (it == seqpro_managers.end()) {
                // spdlog::warn("Species '{}' not found in seqpro_managers", species_name);
                continue;
            }

            try {
                // 提取DNA序列
                std::string sequence = fetchSeq(it->second, chr_name, segment->start, segment->length);

                // 如果是反向链，进行反向互补转换
                if (segment->strand == Strand::REVERSE) {
                    hal::reverseComplement(sequence);
                }

                sequences[species_name] = sequence;

                // 直接使用segment中已有的CIGAR数据
                if (!segment->cigar.empty()) {
                    cigars[species_name] = segment->cigar;
                }

                // spdlog::debug("Extracted sequence for species '{}': {} bp, CIGAR: {} ops, strand: {}",
                //              species_name, sequence.length(), segment->cigar.size(),
                //              (segment->strand == Strand::REVERSE ? "REVERSE" : "FORWARD"));

            } catch (const std::exception& e) {
                spdlog::error("Failed to extract sequence for species '{}': {}", species_name, e.what());
            }
        }

        // spdlog::debug("Extracted {} sequences and {} CIGARs from block",
        //              sequences.size(), cigars.size());

        return {sequences, cigars};
    }

    std::string voteForAncestorSequence(
        const std::unordered_map<std::string, std::string>& aligned_sequences,
        const AncestorNode& ancestor) {

        if (aligned_sequences.empty()) {
            return "";
        }

        // 获取比对长度
        size_t alignment_length = aligned_sequences.begin()->second.length();
        std::string ancestor_sequence;
        ancestor_sequence.reserve(alignment_length);

        // spdlog::debug("Voting for ancestor sequence from {} aligned sequences, length: {}",
        //              aligned_sequences.size(), alignment_length);

        // 按列进行投票
        for (size_t pos = 0; pos < alignment_length; ++pos) {
            std::map<char, int> base_counts;

            // 统计该位置所有后代叶子的碱基
            for (const std::string& leaf : ancestor.descendant_leaves) {
                auto it = aligned_sequences.find(leaf);
                if (it != aligned_sequences.end() && pos < it->second.length()) {
                    char base = std::toupper(it->second[pos]);

                    // 只统计有效碱基，忽略gap和N
                    if (base != '-' && base != 'N' &&
                        (base == 'A' || base == 'C' || base == 'G' || base == 'T')) {
                        base_counts[base]++;
                    }
                }
            }

            // 选择出现次数最多的碱基
            char best_base = 'N';
            int max_count = 0;

            for (const auto& [base, count] : base_counts) {
                if (count > max_count) {
                    max_count = count;
                    best_base = base;
                }
            }

            // 只添加非gap字符到最终序列
            if (best_base != '-') {
                ancestor_sequence += best_base;
            }
        }

        // spdlog::debug("Voting completed: {} -> {} bp", alignment_length, ancestor_sequence.length());
        return ancestor_sequence;
    }

    std::string reconstructSegmentByVoting(
        const AncestorSegmentInfo& segment,
        const AncestorNode& ancestor,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {

        // 1. 从block中提取所有物种的序列和CIGAR
        auto [sequences, cigars] = extractSequencesAndCigarsFromBlock(segment.source_block, seqpro_managers);

        if (sequences.empty()) {
            spdlog::warn("No sequences extracted from block for segment");
            return "";
        }

        // 如果只有一个序列，直接返回（优先选择祖先后代叶子）
        if (sequences.size() == 1) {
            return sequences.begin()->second;
        }

        // 2. 选择参考序列：直接使用 source_block 的参考染色体对应的物种
        std::string ref_key;
        if (segment.source_block) {
            std::shared_lock blk_lock(segment.source_block->rw);
            const auto& ref_chr = segment.source_block->ref_chr;
            for (const auto& [species_chr, _head] : segment.source_block->anchors) {
                if (species_chr.second == ref_chr) {
                    // 以物种名作为 key（与 sequences/cigars 的 key 一致）
                    const std::string& species_name = species_chr.first;
                    if (sequences.find(species_name) != sequences.end()) {
                        ref_key = species_name;
                        break;
                    }
                }
            }
        }
        // 回退：若异常未找到，退到任意一个已提取的序列，避免崩溃
        if (ref_key.empty()) {
            ref_key = sequences.begin()->first;
            spdlog::debug("Selected reference '{}' (fallback to any)", ref_key);
        }

        // 3. 准备CIGAR数据
        std::unordered_map<ChrName, Cigar_t> final_cigars;
        for (const auto& [species, cigar] : cigars) {
            final_cigars[species] = cigar;
        }

        // spdlog::debug("Using {} sequences for alignment, {} CIGARs for non-reference",
        //              sequences.size(), final_cigars.size());

        // 4. 使用mergeAlignmentByRef创建多序列比对
        try {
            mergeAlignmentByRef(ref_key, sequences, final_cigars);
        } catch (const std::exception& e) {
            spdlog::warn("mergeAlignmentByRef failed for segment: {}", e.what());
            // 回退到参考序列
            auto it = sequences.find(ref_key);
            return it != sequences.end() ? it->second : "";
        }

        // 5. 按列投票重建祖先序列
        return voteForAncestorSequence(sequences, ancestor);
    }

    std::string buildAncestorSequenceByVoting(
        const AncestorReconstructionData& data,
        const AncestorNode& ancestor,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
        const std::string& chr_name) {

        // spdlog::debug("Building ancestor sequence using voting method from {} segments for chr '{}'",
        //              data.segments.size(), chr_name);

        std::string full_sequence;
        size_t total_segments = data.segments.size();
        size_t gaps_added = 0;
        size_t segments_with_blocks_added = 0;
        hal_index_t current_pos = 0;  // 在祖先序列中的当前位置

        for (size_t i = 0; i < total_segments; ++i) {
            const auto& segment_info = data.segments[i];

            // 处理gap_before
            if (segment_info.need_gap_before) {
                if (full_sequence.empty() || full_sequence.back() != 'N') {
                    full_sequence += 'N';
                    current_pos++;
                    gaps_added++;
                }
            }

            // 使用投票法重建该segment的序列
            try {
                std::string segment_sequence = reconstructSegmentByVoting(segment_info, ancestor, seqpro_managers);
                // 创建祖先segment并加入到block中
                if (segment_info.source_block && !segment_sequence.empty()) {
                    SegPtr ancestor_seg = Segment::create(
                        current_pos,                    // 在祖先序列中的位置
                        segment_sequence.length(),      // 实际重建的长度
                        Strand::FORWARD,               // 祖先序列总是正向
                        Cigar_t{},                     // 祖先segment没有CIGAR
                        AlignRole::PRIMARY,
                        SegmentRole::SEGMENT,
                        segment_info.source_block
                    );

                    // 将祖先segment注册到block中（使用统一命名 ancestorName.chrN）
                    {
                        std::unique_lock lk(segment_info.source_block->rw);
                        SpeciesChrPair ancestor_key{ancestor.node_name, chr_name};
                        segment_info.source_block->anchors[ancestor_key] = ancestor_seg;
                    }

                    segments_with_blocks_added++;
                    // spdlog::debug("  Added ancestor segment to block: pos={}, len={}",
                    //              current_pos, segment_sequence.length());
                }

                full_sequence += segment_sequence;
                current_pos += segment_sequence.length();

                // spdlog::debug("  Added segment {} using voting: {} bp", i, segment_sequence.length());
            } catch (const std::exception& e) {
                spdlog::error("Failed to reconstruct segment {} using voting: {}", i, e.what());
                throw;
            }

            // 处理gap_after
            if (segment_info.need_gap_after) {
                full_sequence += 'N';
                current_pos++;
                gaps_added++;
            }
        }

        // spdlog::info("Built ancestor sequence using voting: {} segments, {} gaps, {} total length, {} segments added to blocks",
        //             total_segments, gaps_added, full_sequence.length(), segments_with_blocks_added);

        return full_sequence;
    }



    std::map<std::string, std::map<std::string, std::string>> buildAllAncestorSequencesByVoting(
        const std::map<std::string, AncestorReconstructionData>& ancestor_reconstruction_data,
        const std::vector<AncestorNode>& ancestor_nodes,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
        const std::vector<std::pair<std::string, std::string>>& reconstruction_plan,
        hal::AlignmentPtr alignment) {

        spdlog::info("Building sequences for {} ancestors using voting method", ancestor_reconstruction_data.size());

        std::map<std::string, std::map<std::string, std::string>> ancestor_sequences;

        // 并行为每个祖先构建序列（计算并行，HAL 写入加锁）
        ThreadPool pool(std::max(1u, std::thread::hardware_concurrency()));
        std::mutex seq_write_mutex; // 保护 ancestor_sequences 的并发写入

        for (const auto& [ancestor_name, ref_leaf] : reconstruction_plan) {
            pool.enqueue([&, ancestor_name, ref_leaf]() {
            auto it = ancestor_reconstruction_data.find(ancestor_name);
            if (it == ancestor_reconstruction_data.end()) {
                spdlog::warn("Ancestor '{}' not found in reconstruction data", ancestor_name);
                    return;
            }
            const auto& data = it->second;
            spdlog::info("Building sequence for ancestor '{}' using voting", ancestor_name);

            // 找到对应的祖先节点
            const AncestorNode* ancestor = nullptr;
            for (const auto& node : ancestor_nodes) {
                    if (node.node_name == ancestor_name) { ancestor = &node; break; }
                }
            if (!ancestor) {
                spdlog::error("Cannot find ancestor node: {}", ancestor_name);
                    return;
                }

                // 读取参考叶子的所有染色体名称
                std::vector<std::string> chr_names;
                if (auto ref_it = seqpro_managers.find(data.reference_leaf); ref_it != seqpro_managers.end()) {
                    std::visit([&chr_names](const auto& mgr) { chr_names = mgr->getSequenceNames(); }, *ref_it->second);
                }
                if (chr_names.empty()) {
                    spdlog::error("No chromosomes found for reference leaf: {}", data.reference_leaf);
                    return;
                }

                struct ChrSeq { std::string name; std::string seq; size_t segs; };
                std::vector<ChrSeq> built;
                built.reserve(chr_names.size());

                std::map<std::string, std::string> chr_sequences;
                size_t total_length = 0;

                for (size_t chr_idx = 0; chr_idx < chr_names.size(); ++chr_idx) {
                    const auto& chr_name = chr_names[chr_idx];
                    std::string ancestor_chr_name = ancestor_name + ".chr" + std::to_string(chr_idx + 1);

                    AncestorReconstructionData chr_data;
                    chr_data.reference_leaf = data.reference_leaf;
                    chr_data.getSequenceId = data.getSequenceId;
                    SequenceId chr_id = data.getSequenceId(chr_name);
                    for (const auto& segment : data.segments) if (segment.chr_id == chr_id) chr_data.segments.push_back(segment);

                    if (!chr_data.segments.empty()) {
                        std::string chr_sequence = buildAncestorSequenceByVoting(chr_data, *ancestor, seqpro_managers, ancestor_chr_name);
                        built.push_back({ancestor_chr_name, chr_sequence, chr_data.segments.size()});
                        chr_sequences[ancestor_chr_name] = chr_sequence;
                        total_length += chr_sequence.length();
                    }
                }

                // 写入 HAL（需要加全局锁）
                if (alignment && !built.empty()) {
                    std::lock_guard<std::mutex> lk(g_hal_write_mutex);
                    hal::Genome* genome = alignment->openGenome(ancestor_name);
                    if (!genome) {
                        std::string parent = ancestor->parent_name.empty() ? alignment->getRootName() : ancestor->parent_name;
                        if (parent.empty()) parent = ancestor_name;
                        genome = alignment->addLeafGenome(ancestor_name, parent, ancestor->branch_length);
                    }
                    if (genome) {
                    std::vector<hal::Sequence::Info> dims;
                    dims.reserve(built.size());
                        for (const auto& cs : built) dims.emplace_back(cs.name, static_cast<hal_size_t>(cs.seq.size()), 0, 0);
                    genome->setDimensions(dims);
                        for (const auto& cs : built) if (auto* hal_seq = genome->getSequence(cs.name)) hal_seq->setString(cs.seq);
                        alignment->closeGenome(genome);
                    }
                }

                {
                    std::lock_guard<std::mutex> lk(seq_write_mutex);
                ancestor_sequences[ancestor_name] = std::move(chr_sequences);
                }

                spdlog::info("Successfully built sequence for ancestor '{}' using voting: {} chromosomes, {} bp total",
                            ancestor_name, built.size(), total_length);
            });
            }
        pool.waitAllTasksDone();

        // 输出统计信息
        spdlog::info("Ancestor sequence construction using voting completed:");
        size_t total_length = 0;
        for (const auto& [ancestor_name, chr_sequences] : ancestor_sequences) {
            size_t ancestor_total = 0;
            for (const auto& [chr_name, chr_seq] : chr_sequences) {
                ancestor_total += chr_seq.length();
            }
            spdlog::info("  Ancestor '{}': {} chromosomes, {} bp total", ancestor_name, chr_sequences.size(), ancestor_total);
            total_length += ancestor_total;
        }
        spdlog::info("  Total ancestor sequence length: {} bp", total_length);

        return ancestor_sequences;
    }



    std::pair<bool, std::string> validateLeafNames(
        const NewickParser& parser,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {

        std::vector<std::string> tree_leaves = parser.getLeafNames();
        std::set<std::string> manager_species;

        for (const auto& [species_name, _] : seqpro_managers) {
            manager_species.insert(species_name);
        }

        std::vector<std::string> missing_in_tree;
        std::vector<std::string> missing_in_managers;

        // 检查树中的叶节点是否都在managers中
        for (const auto& leaf : tree_leaves) {
            if (manager_species.find(leaf) == manager_species.end()) {
                missing_in_managers.push_back(leaf);
            }
        }

        // 检查managers中的物种是否都在树中
        std::set<std::string> tree_leaf_set(tree_leaves.begin(), tree_leaves.end());
        for (const auto& species : manager_species) {
            if (tree_leaf_set.find(species) == tree_leaf_set.end()) {
                missing_in_tree.push_back(species);
            }
        }

        if (missing_in_tree.empty() && missing_in_managers.empty()) {
            return {true, "All leaf names match"};
        }

        std::string error_msg = "Leaf name mismatches: ";
        if (!missing_in_tree.empty()) {
            error_msg += "Missing in tree: ";
            for (const auto& name : missing_in_tree) {
                error_msg += name + " ";
            }
        }
        if (!missing_in_managers.empty()) {
            error_msg += "Missing in managers: ";
            for (const auto& name : missing_in_managers) {
                error_msg += name + " ";
            }
        }

        return {false, error_msg};
    }

    // ========================================
    // HAL基础结构创建
    // ========================================
    void setupGenomeSequences(
        hal::AlignmentPtr alignment,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {

        spdlog::info("Setting up sequence data for all genomes...");

        for (const auto& [species_name, seq_mgr] : seqpro_managers) {
            hal::Genome* genome = alignment->openGenome(species_name);
            if (!genome) {
                spdlog::warn("Cannot open genome: {}", species_name);
                continue;
            }

            // 获取序列信息并创建维度
            std::vector<hal::Sequence::Info> sequence_dimensions;

            std::visit([&](const auto& mgr) {
                auto chr_names = mgr->getSequenceNames();
                for (const auto& chr_name : chr_names) {
                    hal_size_t length = mgr->getSequenceLength(chr_name);
                    sequence_dimensions.emplace_back(chr_name, length, 0, 0);
                }
            }, *seq_mgr);

            // 设置基因组维度和DNA数据
            if (!sequence_dimensions.empty()) {
                genome->setDimensions(sequence_dimensions);

                // 添加DNA数据（暂时用N填充）
                std::string dna_data(genome->getSequenceLength(), 'N');
                genome->setString(dna_data);

                spdlog::info("  Setup completed for genome: {} ({} sequences, {} bp)",
                           species_name, sequence_dimensions.size(), genome->getSequenceLength());
            } else {
                spdlog::warn("  No sequences found for genome: {}", species_name);
            }

            alignment->closeGenome(genome);
        }
    }

    void setupLeafGenomesWithRealDNA(
        hal::AlignmentPtr alignment,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {

        spdlog::info("Setting up leaf genomes with real chromosome dimensions and DNA (parallel read, locked write)...");

        ThreadPool pool(std::max(1u, std::thread::hardware_concurrency()));
        for (const auto& [species_name, seq_mgr] : seqpro_managers) {
            pool.enqueue([alignment, species_name, seq_mgr]() {
                // 1) 读取该叶物种的所有染色体名称与长度（无 HAL 访问，可并行）
            std::vector<std::string> chr_names;
            std::vector<hal::Sequence::Info> dims;
            std::visit([&](const auto& mgr) {
                chr_names = mgr->getSequenceNames();
                    dims.reserve(chr_names.size());
                for (const auto& chr : chr_names) {
                    hal_size_t len = mgr->getSequenceLength(chr);
                    dims.emplace_back(chr, len, 0, 0);
                }
            }, *seq_mgr);

            if (dims.empty()) {
                spdlog::warn("  No chromosomes found for leaf genome: {}", species_name);
                    return;
                }

                // 2) 打开基因组并设置维度（HAL 写：需加锁）
                hal::Genome* genome = nullptr;
                {
                    std::lock_guard<std::mutex> lk(g_hal_write_mutex);
                    genome = alignment->openGenome(species_name);
                    if (!genome) {
                        spdlog::warn("Cannot open leaf genome: {}", species_name);
                        return;
                    }
            genome->setDimensions(dims);
                }

                // 3) 逐条染色体读取 DNA 并写入（读取无锁，写 HAL 加锁）
            std::visit([&](const auto& mgr) {
                    using PtrType = std::decay_t<decltype(mgr)>;
                for (const auto& chr : chr_names) {
                    auto chr_id = mgr->getSequenceId(chr);
                    hal_size_t len = mgr->getSequenceLength(chr);
                    std::string dna;
                    if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
                        dna = mgr->getSubSequence(chr_id, 0, len);
                    } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                        dna = mgr->getOriginalManager().getSubSequence(chr_id, 0, len);
                    }

                        if (!dna.empty()) {
                            std::lock_guard<std::mutex> lk(g_hal_write_mutex);
                            if (auto* hal_seq = genome->getSequence(chr)) {
                        hal_seq->setString(dna);
                            }
                    }
                }
            }, *seq_mgr);

                {
                    std::lock_guard<std::mutex> lk(g_hal_write_mutex);
            spdlog::info("  Leaf genome '{}' set: {} sequences, {} bp",
                         species_name, genome->getNumSequences(), genome->getSequenceLength());
            alignment->closeGenome(genome);
        }
            });
        }

        pool.waitAllTasksDone();
    }





    void createGenomesFromPhylogeny(
        hal::AlignmentPtr alignment,
        const std::vector<AncestorNode>& ancestor_nodes,
        const NewickParser& parser,
        const std::string& preferred_root_name) {

        spdlog::info("Creating genomes from phylogeny (root / internal ancestors / leaves)...");

        // 1) 确定真实根节点名称（parent_name 为空者）
        std::string root_name;
        for (const auto& anc : ancestor_nodes) {
            if (anc.parent_name.empty()) {
                root_name = anc.node_name;
                break;
            }
        }
        if (root_name.empty()) {
            root_name = preferred_root_name.empty() ? std::string("ancestor") : preferred_root_name; // 兜底
            spdlog::warn("No explicit root found from ancestor_nodes, fallback to '{}'", root_name);
        }

        // 2) 创建或复用真实根节点
        hal::Genome* root = alignment->openGenome(root_name);
        if (!root) {
            root = alignment->addRootGenome(root_name);
            if (!root) {
                throw std::runtime_error("Failed to create real root genome: " + root_name);
            }
            spdlog::info("Created real root genome: {}", root_name);
        } else {
            spdlog::info("Reusing existing root genome: {}", root_name);
        }

        // 3) 按 Newick 左->右顺序递归创建内部祖先与叶，保证拓扑顺序与输入一致
        const auto& nodes = parser.getNodes();

        auto getNameById = [&](int id) -> std::string {
            if (id < 0) return std::string();
            for (const auto& n : nodes) if (n.id == id) return n.name.empty() ? (std::string("internal_") + std::to_string(n.id)) : n.name;
            return std::string();
        };

        // 找到根 id
        int root_id = -1;
        for (const auto& n : nodes) if (n.father == -1) { root_id = n.id; break; }
        if (root_id == -1) {
            spdlog::error("Failed to locate root id from parser when creating genomes");
            return;
        }

        // 递归创建函数，严格按 leftChild -> rightChild 顺序
        std::function<void(int)> createSubtree = [&](int node_id) {
            // 当前节点信息
            const NewickTreeNode* cur = nullptr;
            for (const auto& n : nodes) { if (n.id == node_id) { cur = &n; break; } }
            if (cur == nullptr) return;

            std::string cur_name = getNameById(cur->id);
            std::string parent_name = getNameById(cur->father);
            double branch_len = cur->branchLength;

            // 创建当前节点（除根外均作为父的子节点）
            if (cur->father == -1) {
                // 已在上面创建/复用 root
            } else {
                hal::Genome* g = alignment->openGenome(cur_name);
                if (!g) {
                    hal::Genome* parent_g = alignment->openGenome(parent_name);
                    if (!parent_g) {
                        spdlog::error("Parent genome not found when creating '{}': {}", cur_name, parent_name);
                    } else {
                        g = alignment->addLeafGenome(cur_name, parent_name, branch_len);
                        if (!g) {
                            spdlog::error("Failed to create genome: {} (parent: {})", cur_name, parent_name);
                        } else {
                            spdlog::info("Created genome: {} (parent: {}, bl={})", cur_name, parent_name, branch_len);
                        }
                    }
                }
            }

            // 递归创建子节点，按左->右，保证输出顺序
            if (cur->leftChild != -1) {
                createSubtree(cur->leftChild);
            } else {
                // 兼容非严格二叉：遍历所有以 cur 为父的孩子，按 id 升序
                for (const auto& n : nodes) if (n.father == cur->id && n.id != cur->rightChild) {
                    createSubtree(n.id);
                }
            }
            if (cur->rightChild != -1) {
                createSubtree(cur->rightChild);
            }
        };

        // 从根开始递归，严格保序
        createSubtree(root_id);
    }

    void createAncestorGenomesWithCorrectDimensions(
        hal::AlignmentPtr alignment,
        const std::map<std::string, std::string>& ancestor_sequences,
        const std::vector<AncestorNode>& ancestor_nodes,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {

        spdlog::info("Creating ancestor genomes with correct dimensions...");

        // 获取参考基因组的染色体信息
        std::vector<std::string> chr_names;
        if (!seqpro_managers.empty()) {
            const auto& [ref_species, ref_mgr] = *seqpro_managers.begin();
            std::visit([&chr_names](const auto& mgr) {
                chr_names = mgr->getSequenceNames();
            }, *ref_mgr);
        }

        // 首先创建根节点
        for (const auto& ancestor : ancestor_nodes) {
            if (ancestor.parent_name.empty()) {
                // 这是根节点，替换临时根节点
                std::string temp_root_name = "temp_root";

                // 计算根节点的正确维度
                auto seq_it = ancestor_sequences.find(ancestor.node_name);
                if (seq_it == ancestor_sequences.end()) {
                    spdlog::error("No sequence found for root ancestor: {}", ancestor.node_name);
                    continue;
                }

                std::vector<hal::Sequence::Info> root_dimensions;
                hal_size_t total_length = seq_it->second.length();

                // 按比例分配到各个染色体
                hal_size_t offset = 0;
                for (size_t i = 0; i < chr_names.size(); ++i) {
                    hal_size_t chr_length;
                    if (i == chr_names.size() - 1) {
                        chr_length = total_length - offset;
                    } else {
                        chr_length = total_length / chr_names.size(); // 简单平均分配
                    }
                    root_dimensions.emplace_back(chr_names[i], chr_length, 0, 0);
                    offset += chr_length;
                }

                // 创建真正的根节点
                hal::Genome* root_genome = alignment->addRootGenome(ancestor.node_name);
                if (!root_genome) {
                    throw std::runtime_error("Failed to create root genome: " + ancestor.node_name);
                }

                // 设置正确的维度
                root_genome->setDimensions(root_dimensions);

                // 设置DNA序列
                offset = 0;
                for (size_t i = 0; i < chr_names.size(); ++i) {
                    hal::Sequence* sequence = root_genome->getSequence(chr_names[i]);
                    if (sequence) {
                        hal_size_t chr_length = root_dimensions[i]._length;
                        std::string chr_sequence = seq_it->second.substr(offset, chr_length);
                        sequence->setString(chr_sequence);
                        offset += chr_length;
                    }
                }

                spdlog::info("Created root genome '{}' with correct dimensions: {} bp",
                           ancestor.node_name, total_length);
                break;
            }
        }

        // 然后创建内部节点（祖先）
        for (const auto& ancestor : ancestor_nodes) {
            if (!ancestor.parent_name.empty()) {
                auto seq_it = ancestor_sequences.find(ancestor.node_name);
                if (seq_it == ancestor_sequences.end()) {
                    spdlog::error("No sequence found for ancestor: {}", ancestor.node_name);
                    continue;
                }

                // 计算正确的维度
                std::vector<hal::Sequence::Info> ancestor_dimensions;
                hal_size_t total_length = seq_it->second.length();

                hal_size_t offset = 0;
                for (size_t i = 0; i < chr_names.size(); ++i) {
                    hal_size_t chr_length;
                    if (i == chr_names.size() - 1) {
                        chr_length = total_length - offset;
                    } else {
                        chr_length = total_length / chr_names.size();
                    }
                    ancestor_dimensions.emplace_back(chr_names[i], chr_length, 0, 0);
                    offset += chr_length;
                }

                // 创建祖先基因组
                hal::Genome* ancestor_genome = alignment->addLeafGenome(ancestor.node_name, ancestor.parent_name, ancestor.branch_length);
                if (!ancestor_genome) {
                    spdlog::error("Failed to create ancestor genome: {}", ancestor.node_name);
                    continue;
                }

                // 设置正确的维度
                ancestor_genome->setDimensions(ancestor_dimensions);

                // 设置DNA序列
                offset = 0;
                for (size_t i = 0; i < chr_names.size(); ++i) {
                    hal::Sequence* sequence = ancestor_genome->getSequence(chr_names[i]);
                    if (sequence) {
                        hal_size_t chr_length = ancestor_dimensions[i]._length;
                        std::string chr_sequence = seq_it->second.substr(offset, chr_length);
                        sequence->setString(chr_sequence);
                        offset += chr_length;
                    }
                }

                spdlog::info("Created ancestor genome '{}' with correct dimensions: {} bp",
                           ancestor.node_name, total_length);
            }
        }

        spdlog::info("All ancestor genomes created with correct dimensions");
    }

    std::string reconstructNewickFromParser(const NewickParser& parser) {
        const auto& nodes = parser.getNodes();
        if (nodes.empty()) {
            return "";
        }

        // 找到根节点
        int root_id = -1;
        for (const auto& node : nodes) {
            if (node.father == -1) {
                root_id = node.id;
                break;
            }
        }

        if (root_id == -1) {
            return "";
        }

        // 递归构建Newick字符串
        std::function<std::string(int)> buildNewick = [&](int node_id) -> std::string {
            // 找到对应的节点
            const NewickTreeNode* current_node = nullptr;
            for (const auto& node : nodes) {
                if (node.id == node_id) {
                    current_node = &node;
                    break;
                }
            }

            if (!current_node) {
                return "";
            }

            std::string result;

            // 如果不是叶节点，需要处理子节点
            if (!current_node->isLeaf) {
                result += "(";
                std::vector<std::string> children_strs;

                // 收集所有子节点
                for (const auto& node : nodes) {
                    if (node.father == node_id) {
                        std::string child_str = buildNewick(node.id);
                        if (!child_str.empty()) {
                            children_strs.push_back(child_str);
                        }
                    }
                }

                // 连接子节点字符串
                for (size_t i = 0; i < children_strs.size(); ++i) {
                    if (i > 0) result += ",";
                    result += children_strs[i];
                }

                result += ")";
            }

            // 节点名称（为避免 HAL 内部对空名处理异常，这里为所有节点提供兜底名）
            std::string node_name = current_node->name.empty() ? (std::string("internal_") + std::to_string(current_node->id)) : current_node->name;
            result += node_name;

            // 添加分支长度（除了根节点）
            if (current_node->father != -1) {
                result += ":" + std::to_string(current_node->branchLength);
            }

            return result;
        };

        std::string newick = buildNewick(root_id) + ";";
        return newick;
    }

    void applyPhylogeneticTree(
        hal::AlignmentPtr alignment,
        const NewickParser& parser) {

        spdlog::info("Applying phylogenetic tree structure to HAL alignment...");

        try {
            // 重建Newick字符串从解析的树结构
            std::string newick_tree = reconstructNewickFromParser(parser);

            if (!newick_tree.empty()) {
                // 使用HAL API设置树结构
                alignment->replaceNewickTree(newick_tree);
                spdlog::info("Successfully applied phylogenetic tree: {}", newick_tree);
            } else {
                spdlog::warn("Failed to reconstruct Newick tree from parser");
            }
        } catch (const std::exception& e) {
            spdlog::error("Failed to apply phylogenetic tree: {}", e.what());
        }
    }

    // ========================================
    // 验证和工具函数
    // ========================================

    void validateHalFile(hal::AlignmentPtr alignment, const FilePath& hal_path) {
        spdlog::info("Validating HAL file structure...");
        spdlog::info("Total genomes created: {}", alignment->getNumGenomes());

        try {
            std::string root_name = alignment->getRootName();
            if (!root_name.empty()) {
                spdlog::info("Root genome: {}", root_name);
            }

            std::string tree = alignment->getNewickTree();
            if (!tree.empty()) {
                spdlog::info("Phylogenetic tree: {}", tree);
            }
        } catch (const std::exception& e) {
            spdlog::warn("Error accessing HAL structure: {}", e.what());
        }

        spdlog::info("HAL file validation completed: {}", hal_path.string());
    }

    std::string generateRootName(const std::vector<std::string>& leaf_names) {
        if (leaf_names.empty()) {
            return "ancestor";
        }

        // 生成一个不与任何叶节点冲突的根节点名称
        std::set<std::string> leaf_set(leaf_names.begin(), leaf_names.end());

        std::string base_name = "ancestor";
        std::string root_name = base_name;
        int counter = 1;

        while (leaf_set.find(root_name) != leaf_set.end()) {
            root_name = base_name + "_" + std::to_string(counter);
            counter++;
        }

        return root_name;
    }

    // ========================================
    // 第三阶段（第一遍）：统计段数并更新HAL维度
    // ========================================

    std::vector<CurrentBlockMapping> analyzeCurrentBlock(
        BlockPtr block,
        const std::vector<AncestorNode>& ancestor_nodes) {

        std::vector<CurrentBlockMapping> mappings;
        if (!block) return mappings;

        // 收集当前 block 中的 (species -> {chr, segment})
        std::unordered_map<std::string, std::pair<std::string, SegPtr>> species_to_entry;
        {
            std::shared_lock lk(block->rw);
            for (const auto& [species_chr, segment] : block->anchors) {
                const std::string& species = species_chr.first;
                const std::string& chr = species_chr.second;
                if (segment) {
                    species_to_entry[species] = {chr, segment};
                }
            }
        }

        if (species_to_entry.empty()) return mappings;

        // 对每个在本 block 出现的祖先，聚合其直系子到同一个父段上
        for (const auto& anc : ancestor_nodes) {
            auto itParent = species_to_entry.find(anc.node_name);
            if (itParent == species_to_entry.end()) continue; // 祖先不在当前块

            const auto& [pChr, pSeg] = itParent->second;
            if (!pSeg) continue;

            CurrentBlockMapping map{};
            map.parent_genome = anc.node_name;
            map.parent_chr_name = pChr;
            map.parent_start = static_cast<hal_size_t>(pSeg->start);
            map.parent_length = static_cast<hal_size_t>(pSeg->length);

            for (const auto& child_name : anc.direct_children_names) {
                auto itChild = species_to_entry.find(child_name);
                if (itChild == species_to_entry.end()) continue; // 子不在当前块

                const auto& [cChr, cSeg] = itChild->second;
                if (!cSeg) continue;

                CurrentBlockMapping::ChildInfo ci{};
                ci.child_genome = child_name;
                ci.child_chr_name = cChr;
                ci.child_start = static_cast<hal_size_t>(cSeg->start);
                ci.child_length = static_cast<hal_size_t>(cSeg->length);
                ci.is_reversed = (cSeg->strand != pSeg->strand);
                map.children.push_back(std::move(ci));
            }

            if (!map.children.empty()) {
                mappings.push_back(std::move(map));
            }
        }

        return mappings;
    }

    /**
     * 从祖先染色体序列中提取特定block的序列片段
     */
    std::string extractAncestorSequenceForBlock(
        BlockPtr block,
        const AncestorNode& ancestor,
        const std::map<std::string, std::map<std::string, std::string>>& ancestor_sequences) {

        if (!block) return "";

        // 1. 查找祖先的染色体序列映射
        auto ancestor_it = ancestor_sequences.find(ancestor.node_name);
        if (ancestor_it == ancestor_sequences.end()) {
            spdlog::debug("Ancestor '{}' sequences not found", ancestor.node_name);
            return "";
        }

        // 2. 在block中查找该祖先的segment信息
        SegPtr ancestor_segment = nullptr;
        std::string chr_name;
        {
            std::shared_lock blk_lock(block->rw);
            for (const auto& [species_chr, segment] : block->anchors) {
                if (species_chr.first == ancestor.node_name) {
                    ancestor_segment = segment;
                    chr_name = species_chr.second;  // 获取染色体名称
                    break;
                }
            }
        }

        if (!ancestor_segment) {
            spdlog::debug("Ancestor '{}' segment not found in block", ancestor.node_name);
            return "";
        }

        // 3. 查找对应染色体的序列
        const auto& chr_sequences = ancestor_it->second;
        auto chr_it = chr_sequences.find(chr_name);
        if (chr_it == chr_sequences.end()) {
            spdlog::debug("Ancestor '{}' chromosome '{}' sequence not found", ancestor.node_name, chr_name);
            return "";
        }

        // 4. 从染色体序列中提取片段
        const std::string& chr_sequence = chr_it->second;
        size_t start = ancestor_segment->start;
        size_t length = ancestor_segment->length;

        if (start + length > chr_sequence.length()) {
            spdlog::warn("Ancestor '{}' segment range [{}:{}] exceeds chromosome '{}' length {}",
                        ancestor.node_name, start, start + length, chr_name, chr_sequence.length());
            return "";
        }

        std::string fragment = chr_sequence.substr(start, length);

        // 5. 处理反向链（如果需要）
        if (ancestor_segment->strand == Strand::REVERSE) {
            hal::reverseComplement(fragment);
        }

        // spdlog::debug("Extracted ancestor '{}' fragment from chr '{}': {}:{} ({} bp, strand: {})",
        //              ancestor.node_name, chr_name, start, start + length, fragment.length(),
        //              (ancestor_segment->strand == Strand::REVERSE ? "REVERSE" : "FORWARD"));

        return fragment;
    }

    /**
     * 添加祖先序列到已对齐的叶子序列中
     */
    void addAncestorSequencesToAlignment(
        BlockPtr block,
        const std::vector<AncestorNode>& ancestor_nodes,
        const std::map<std::string, AncestorReconstructionData>& ancestor_data,
        const std::map<std::string, std::map<std::string, std::string>>& ancestor_sequences,
        std::unordered_map<std::string, std::string>& aligned_sequences) {

        if (!block) return;

        // 遍历所有祖先节点
        for (const auto& ancestor : ancestor_nodes) {
            // 检查该祖先是否在当前block中有segment
            bool ancestor_in_block = false;
            {
                std::shared_lock blk_lock(block->rw);
                for (const auto& [species_chr, segment] : block->anchors) {
                    if (species_chr.first == ancestor.node_name) {
                        ancestor_in_block = true;
                        break;
                    }
                }
            }

            if (!ancestor_in_block) continue;

            // 从祖先完整序列中提取该祖先在此block的序列片段
            std::string ancestor_seq = extractAncestorSequenceForBlock(block, ancestor, ancestor_sequences);
            if (!ancestor_seq.empty()) {
                aligned_sequences[ancestor.node_name] = ancestor_seq;
                // spdlog::debug("Added ancestor '{}' sequence to alignment: {} bp",
                //              ancestor.node_name, ancestor_seq.length());
            }
        }
    }

    /**
     * 区域信息结构
     */
    struct Region {
        size_t start_col;
        size_t end_col;
        std::set<std::string> participants;

        size_t length() const { return end_col - start_col; }
    };

    /**
     * 分析参与者集合变化，拆分成连续区域
     */
    std::vector<Region> analyzeRegionsByParticipants(
        const std::unordered_map<std::string, std::string>& aligned_sequences,
        const std::unordered_set<std::string>& leaf_species,
        size_t alignment_length) {

        std::vector<Region> regions;
        if (alignment_length == 0) return regions;

        std::set<std::string> current_participants;
        size_t region_start = 0;

        for (size_t col = 0; col < alignment_length; ++col) {
            // 分析当前列的参与者
            std::set<std::string> col_participants;
            for (const auto& [species, sequence] : aligned_sequences) {
                if (col >= sequence.length()) continue;

                const char c = sequence[col];
                // 重要：对叶子物种而言，'N' 代表未知碱基，仍然占用基因组坐标；不能当作 gap。
                //      对祖先序列而言，'N' 常被用作“缺失/空位”占位符，需要当作 gap 才能正确拆分映射区域。
                const bool is_gap = (c == '-') || (!leaf_species.count(species) && (c == 'N' || c == 'n'));
                if (!is_gap) {
                    col_participants.insert(species);
                }
            }

            // 检查参与者集合是否发生变化
            if (col_participants != current_participants) {
                // 如果不是第一列，先保存前一个区域
                if (col > 0) {
                    Region region;
                    region.start_col = region_start;
                    region.end_col = col;
                    region.participants = current_participants;
                    if (!region.participants.empty()) {
                        regions.push_back(region);
                    }
                }

                // 开始新区域
                region_start = col;
                current_participants = col_participants;
            }
        }

        // 保存最后一个区域
        if (!current_participants.empty()) {
            Region region;
            region.start_col = region_start;
            region.end_col = alignment_length;
            region.participants = current_participants;
            regions.push_back(region);
        }

        return regions;
    }

    /**
     * 从block中获取指定物种的segment信息
     */
    SegPtr getSegmentFromBlock(BlockPtr block, const std::string& species_name) {
        if (!block) return nullptr;

        std::shared_lock blk_lock(block->rw);
        for (const auto& [species_chr, segment] : block->anchors) {
            if (species_chr.first == species_name) {
                return segment;
            }
        }
        return nullptr;
    }

    /**
     * 从block中获取指定物种的染色体名称
     */
    std::string getChrNameFromBlock(BlockPtr block, const std::string& species_name) {
        if (!block) return "";

        std::shared_lock blk_lock(block->rw);
        for (const auto& [species_chr, segment] : block->anchors) {
            if (species_chr.first == species_name) {
                return species_chr.second;
            }
        }
        return "";
    }

    /**
     * 计算由于gap导致的坐标偏移
     */
    hal_size_t calculateGapOffset(
        const std::string& species_name,
        const Region& region,
        const std::unordered_map<std::string, std::string>& aligned_sequences,
        const std::unordered_set<std::string>& leaf_species) {

        auto it = aligned_sequences.find(species_name);
        if (it == aligned_sequences.end()) return 0;

        const std::string& sequence = it->second;

        // 计算区域开始位置之前有多少个“占用基因组坐标”的字符
        // - 叶子：'N' 仍占位（未知碱基），不能当作 gap
        // - 祖先：'N' 作为占位符时不占位，需要当作 gap
        hal_size_t offset = 0;
        for (size_t i = 0; i < region.start_col && i < sequence.length(); ++i) {
            char c = sequence[i];
            const bool is_gap = (c == '-') || (!leaf_species.count(species_name) && (c == 'N' || c == 'n'));
            if (!is_gap) {
                offset++;
            }
        }

        return offset;
    }

    /**
     * 为单个区域创建映射块
     */
    std::vector<CurrentBlockMapping> createMappingsForRegion(
        const Region& region,
        const std::unordered_map<std::string, std::string>& aligned_sequences,
        const std::unordered_set<std::string>& leaf_species,
        const std::vector<AncestorNode>& ancestor_nodes,
        BlockPtr block) {

        std::vector<CurrentBlockMapping> mappings;

        // 在该区域的参与者中找到所有祖先
        for (const auto& ancestor : ancestor_nodes) {
            if (region.participants.count(ancestor.node_name) == 0) continue;

            // 从block中获取该祖先的原始segment信息
            SegPtr parent_segment = getSegmentFromBlock(block, ancestor.node_name);
            if (!parent_segment) continue;

            CurrentBlockMapping mapping;
            mapping.parent_genome = ancestor.node_name;
            mapping.parent_chr_name = getChrNameFromBlock(block, ancestor.node_name);

            // 计算由于gap导致的坐标偏移（注意反向链需从右端计算）
            hal_size_t parent_offset = calculateGapOffset(ancestor.node_name, region, aligned_sequences, leaf_species);
            {
                const hal_size_t rlen = region.length();
                if (parent_segment->strand == Strand::REVERSE) {
                    mapping.parent_start = parent_segment->start + (parent_segment->length - (parent_offset + rlen));
                } else {
            mapping.parent_start = parent_segment->start + parent_offset;
                }
            }

            // 区域长度就是区域的列数（因为参与者在该区域内没有gap）
            mapping.parent_length = region.length();

            // 处理直系子
            for (const auto& child_name : ancestor.direct_children_names) {
                if (region.participants.count(child_name) == 0) continue;

                SegPtr child_segment = getSegmentFromBlock(block, child_name);
                if (!child_segment) continue;

                CurrentBlockMapping::ChildInfo child_info;
                child_info.child_genome = child_name;
                child_info.child_chr_name = getChrNameFromBlock(block, child_name);

                hal_size_t child_offset = calculateGapOffset(child_name, region, aligned_sequences, leaf_species);
                {
                    const hal_size_t rlen = region.length();
                    if (child_segment->strand == Strand::REVERSE) {
                        child_info.child_start = child_segment->start + (child_segment->length - (child_offset + rlen));
                    } else {
                child_info.child_start = child_segment->start + child_offset;
                    }
                }
                child_info.child_length = region.length();  // 与parent长度相同
                child_info.is_reversed = (child_segment->strand != parent_segment->strand);

                mapping.children.push_back(child_info);
            }

            // 只有当有直系子时才添加映射
            if (!mapping.children.empty()) {
                mappings.push_back(mapping);
            }
        }

        return mappings;
    }

    /**
     * 按列分析并拆分映射关系
     */
    std::vector<CurrentBlockMapping> splitMappingsByColumns(
        const std::unordered_map<std::string, std::string>& aligned_sequences,
        const std::unordered_set<std::string>& leaf_species,
        const std::vector<AncestorNode>& ancestor_nodes,
        const std::string& ref_key,
        BlockPtr block) {

        std::vector<CurrentBlockMapping> result;
        if (aligned_sequences.empty()) return result;

        // aligned_sequences 是 unordered_map，begin() 的元素顺序不稳定；这里取最大长度作为对齐列数。
        size_t alignment_length = 0;
        for (const auto& [_, seq] : aligned_sequences) {
            alignment_length = std::max(alignment_length, seq.length());
        }
        if (alignment_length == 0) return result;

        // 1. 按列分析参与者集合变化，拆分成连续区域
        std::vector<Region> regions = analyzeRegionsByParticipants(aligned_sequences, leaf_species, alignment_length);

        // 2. 为每个区域生成映射块
        for (const auto& region : regions) {
            auto mappings = createMappingsForRegion(region, aligned_sequences, leaf_species, ancestor_nodes, block);
            result.insert(result.end(), mappings.begin(), mappings.end());
        }

        // spdlog::debug("Split alignment into {} regions, generated {} mappings", regions.size(), result.size());
        return result;
    }

    std::vector<CurrentBlockMapping> analyzeBlockWithGapHandling(
        BlockPtr block,
        const std::vector<AncestorNode>& ancestor_nodes,
        const std::map<std::string, AncestorReconstructionData>& ancestor_data,
        const std::map<std::string, std::map<std::string, std::string>>& ancestor_sequences,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {

        std::vector<CurrentBlockMapping> result;
        if (!block) return result;

        // 1. 提取叶子序列和CIGAR（只包含叶子物种）
        auto [leaf_sequences, leaf_cigars] = extractSequencesAndCigarsFromBlock(block, seqpro_managers);
        if (leaf_sequences.empty()) {
            spdlog::warn("No leaf sequences extracted from block");
            return result;
        }

        // 2. 选择参考序列（优先选择block的ref_chr对应的物种）
        std::string ref_key;
        std::string block_ref_chr;
        {
            std::shared_lock blk_lock(block->rw);
            block_ref_chr = block->ref_chr;
            const auto& ref_chr = block_ref_chr;
            for (const auto& [species_chr, segment] : block->anchors) {
                if (species_chr.second == ref_chr) {
                    const std::string& species_name = species_chr.first;
                    if (leaf_sequences.find(species_name) != leaf_sequences.end()) {
                        ref_key = species_name;
                        break;
                    }
                }
            }
        }
        // 回退：如果未找到，使用第一个叶子序列
        if (ref_key.empty()) {
            ref_key = leaf_sequences.begin()->first;
            spdlog::debug("Using fallback reference sequence: {}", ref_key);
        }

        // 3. 对叶子序列进行多序列比对
        //    注意：mergeAlignmentByRef 依赖 cigar 与 ref/qry 序列长度一致，否则可能出现越界异常。
        //    科研软件必须严格：一旦发现不一致，直接报错，避免生成错误 HAL。
        const size_t ref_len = leaf_sequences.at(ref_key).size();
        for (const auto& [species, cigar] : leaf_cigars) {
            if (species == ref_key) {
                continue;
            }
            auto seq_it = leaf_sequences.find(species);
            if (seq_it == leaf_sequences.end()) {
                throw std::runtime_error(
                    "HAL构建失败：Block(ref_chr=" + block_ref_chr +
                    ") 中 cigar 物种 '" + species + "' 未找到对应序列，无法进行多序列对齐");
            }
            AlignCount cnt = countAlignedBases(cigar);
            const size_t qry_len = seq_it->second.size();
            if (cnt.ref_bases != ref_len || cnt.query_bases != qry_len) {
                throw std::runtime_error(
                    "HAL构建失败：Block(ref_chr=" + block_ref_chr + ") 中 cigar 长度与序列长度不一致：species='" +
                    species + "' cigar(ref=" + std::to_string(cnt.ref_bases) + ",qry=" +
                    std::to_string(cnt.query_bases) + ") seq(ref=" + std::to_string(ref_len) + ",qry=" +
                    std::to_string(qry_len) + ")");
            }
        }

        try {
            mergeAlignmentByRef(ref_key, leaf_sequences, leaf_cigars);
            // spdlog::debug("Multi-sequence alignment completed for {} leaf sequences", leaf_sequences.size());
        } catch (const std::exception& e) {
            throw std::runtime_error(
                "HAL构建失败：mergeAlignmentByRef 在 Block(ref_chr=" + block_ref_chr + ", ref_species=" + ref_key +
                ") 上失败：" + std::string(e.what()));
        }

        // 4. 添加祖先序列到已对齐的序列中
        addAncestorSequencesToAlignment(block, ancestor_nodes, ancestor_data, ancestor_sequences, leaf_sequences);

        // 5. 按列分析并拆分映射
        std::unordered_set<std::string> leaf_species;
        leaf_species.reserve(leaf_sequences.size());
        for (const auto& [sp, _] : leaf_sequences) leaf_species.insert(sp);
        result = splitMappingsByColumns(leaf_sequences, leaf_species, ancestor_nodes, ref_key, block);

        return result;
    }


    void checkAndUpdateHalDimensions(
        hal::AlignmentPtr alignment,
        SegmentIndexManager& index_manager) {

        // 聚合所有需要更新的基因组名称
        std::set<std::string> genomes_to_update;
        for (const auto& [g, _] : index_manager.top_segment_counts) genomes_to_update.insert(g);
        for (const auto& [g, _] : index_manager.bottom_segment_counts) genomes_to_update.insert(g);

        for (const auto& genome_name : genomes_to_update) {
            hal::Genome* genome = alignment->openGenome(genome_name);
            if (!genome) {
                spdlog::warn("Cannot open genome for dimension update: {}", genome_name);
                continue;
            }

            std::vector<hal::Sequence::UpdateInfo> topUpdates;
            std::vector<hal::Sequence::UpdateInfo> bottomUpdates;

            // 仅对已统计到的序列进行更新，避免无关修改
            if (auto itG = index_manager.top_segment_counts.find(genome_name); itG != index_manager.top_segment_counts.end()) {
                for (const auto& [chr, cnt] : itG->second) {
                    topUpdates.emplace_back(chr, cnt);
                }
            }
            if (auto itG = index_manager.bottom_segment_counts.find(genome_name); itG != index_manager.bottom_segment_counts.end()) {
                for (const auto& [chr, cnt] : itG->second) {
                    bottomUpdates.emplace_back(chr, cnt);
                }
            }

            if (!topUpdates.empty()) {
                genome->updateTopDimensions(topUpdates);
                spdlog::debug("  Updated top dims for '{}': {} sequences", genome_name, topUpdates.size());
            }
            if (!bottomUpdates.empty()) {
                genome->updateBottomDimensions(bottomUpdates);
                spdlog::debug("  Updated bottom dims for '{}': {} sequences", genome_name, bottomUpdates.size());
            }

            index_manager.dimensions_updated_genomes.insert(genome_name);
            alignment->closeGenome(genome);
        }
    }

    void analyzeBlocksAndBuildHalStructure(
        const std::vector<std::weak_ptr<Block>>& blocks,
        const std::vector<AncestorNode>& ancestor_nodes,
        hal::AlignmentPtr alignment,
        const std::map<std::string, AncestorReconstructionData>& ancestor_data,
        const std::map<std::string, std::map<std::string, std::string>>& ancestor_sequences,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {

        spdlog::info("Phase 3 (pass 1): Counting segments with gap-aware mapping...");

        auto __now = [](){ return std::chrono::steady_clock::now(); };
        auto __ms = [](auto a, auto b){ return std::chrono::duration_cast<std::chrono::milliseconds>(b - a).count(); };

        auto t_total_start = __now();
        auto t1_start = __now();

	        size_t totalBlocks = 0;
	        SegmentIndexManager idxMgr;
	        std::map<BlockPtr, std::vector<CurrentBlockMapping>> refined_mappings; // 缓存拆分结果
	        ThreadPool pool(std::max(1u, std::thread::hardware_concurrency()));
	        std::mutex mapping_mutex;
	        std::atomic<bool> pass1_failed{false};
	        std::vector<std::future<void>> pass1_futures;
	        pass1_futures.reserve(blocks.size());

	        // Pass 1: 并行分析每个 block 的 gap-aware 映射
	        for (const auto& wb : blocks) {
	            if (auto block = wb.lock()) {
	                totalBlocks++;
	                pass1_futures.emplace_back(pool.enqueue([&, block]() {
	                    if (pass1_failed.load(std::memory_order_relaxed)) {
	                        return;
	                    }
	                    try {
	                        auto split_mappings = analyzeBlockWithGapHandling(
	                            block, ancestor_nodes, ancestor_data, ancestor_sequences, seqpro_managers);
	                        if (pass1_failed.load(std::memory_order_relaxed)) {
	                            return;
	                        }
	                        std::lock_guard<std::mutex> lk(mapping_mutex);
	                        refined_mappings[block] = std::move(split_mappings);
	                    } catch (...) {
	                        pass1_failed.store(true, std::memory_order_relaxed);
	                        throw;
	                    }
	                }));
	            }
	        }
	        pool.waitAllTasksDone();
	        for (auto& f : pass1_futures) {
	            f.get(); // 传播 worker 中的异常，保证科研软件 fail-fast
	        }

	        auto t1_end = __now();
	        spdlog::info("  Pass1: analyzed {} blocks in {} ms", totalBlocks, __ms(t1_start, t1_end));

        // Pass 2 (并行聚合): 收集、排序、计算分裂点、生成完整段列表
        auto t2_start = __now();
        spdlog::info("Preparing raw segments (grouped by genome/chr, sorted by start)...");
        struct SimpleSegmentInfo {
            hal_index_t start;
            hal_size_t length;
            bool operator<(const SimpleSegmentInfo& other) const { return start < other.start; }
        };
        using SegMap = std::map<std::pair<std::string, std::string>, std::vector<SimpleSegmentInfo>>;
        using BreakpointMap = std::map<std::pair<std::string, std::string>, std::set<hal_index_t>>;

        // 2.1 并行聚合原始段和分裂点
        auto t21_start = __now();
        std::vector<std::future<std::tuple<SegMap, SegMap, BreakpointMap, BreakpointMap>>> futures;
        std::vector<BlockPtr> block_vec;
        for (const auto& [b, _] : refined_mappings) block_vec.push_back(b);

        for (auto& block : block_vec) {
            futures.emplace_back(pool.enqueue([&, block]() {
                SegMap local_raw_bottom, local_raw_top;
                BreakpointMap local_bp_bottom, local_bp_top;
                const auto& mappings = refined_mappings.at(block);

                for (const auto& m : mappings) {
                    auto pKey = std::make_pair(m.parent_genome, m.parent_chr_name);
                    local_raw_bottom[pKey].push_back({(hal_index_t)m.parent_start, m.parent_length});
                    local_bp_bottom[pKey].insert((hal_index_t)m.parent_start);
                    local_bp_bottom[pKey].insert((hal_index_t)m.parent_start + m.parent_length);

                    for (const auto& c : m.children) {
                        auto cKey = std::make_pair(c.child_genome, c.child_chr_name);
                        local_raw_top[cKey].push_back({(hal_index_t)c.child_start, c.child_length});
                        local_bp_top[cKey].insert((hal_index_t)c.child_start);
                        local_bp_top[cKey].insert((hal_index_t)c.child_start + c.child_length);
                    }
                }
                return std::make_tuple(local_raw_bottom, local_raw_top, local_bp_bottom, local_bp_top);
            }));
        }

        SegMap rawBottomSegments, rawTopSegments;
        BreakpointMap bpBottom, bpTop;
        for (auto& fut : futures) {
            auto [lrb, lrt, lbb, lbt] = fut.get();
            for (auto& [k, v] : lrb) rawBottomSegments[k].insert(rawBottomSegments[k].end(), v.begin(), v.end());
            for (auto& [k, v] : lrt) rawTopSegments[k].insert(rawTopSegments[k].end(), v.begin(), v.end());
            for (auto& [k, v] : lbb) bpBottom[k].insert(v.begin(), v.end());
            for (auto& [k, v] : lbt) bpTop[k].insert(v.begin(), v.end());
        }
        auto t21_end = __now();
        spdlog::info("  Pass2.1: aggregated segments from {} blocks in {} ms (rawBottom groups: {}, rawTop groups: {})",
                      block_vec.size(), __ms(t21_start, t21_end), rawBottomSegments.size(), rawTopSegments.size());

        // 2.2 并行排序
        auto t22_start = __now();
        // 注意：不能捕获结构化绑定变量 `v` 的引用（会形成悬垂引用，导致未定义行为/内存破坏）
        for (auto& kv : rawBottomSegments) {
            auto* vec = &kv.second;
            pool.enqueue([vec] { std::sort(vec->begin(), vec->end()); });
        }
        for (auto& kv : rawTopSegments) {
            auto* vec = &kv.second;
            pool.enqueue([vec] { std::sort(vec->begin(), vec->end()); });
        }
        pool.waitAllTasksDone();
        auto t22_end = __now();
        spdlog::info("  Pass2.2: sorted raw segments in {} ms", __ms(t22_start, t22_end));

        std::map<std::pair<std::string, std::string>, std::vector<SimpleSegmentInfo>> bottomSegmentsFull;
        std::map<std::pair<std::string, std::string>, std::vector<SimpleSegmentInfo>> topSegmentsFull;

        // 2.3 并行生成完整段列表
        auto t23_start = __now();
        // 重要：HAL 要求（有 parent 的 genome）每条 sequence 的 top segments 总长度必须等于 sequence 长度；
        //      （有 children 的 genome）每条 sequence 的 bottom segments 总长度必须等于 sequence 长度。
        //      因此即使某条染色体完全没有参与任何 block（bpTop/bpBottom 不含该 key），也必须生成“全长 gap 段”进行覆盖，
        //      否则会因为 HDF5 段坐标的边界存储方式导致前一条染色体的末段被错误拉长，进而触发 parent/child length mismatch。
        std::map<std::pair<std::string, std::string>, hal_size_t> seqLengths; // (genome, chr) -> chr length
        std::set<std::pair<std::string, std::string>> keysNeedTop;
        std::set<std::pair<std::string, std::string>> keysNeedBottom;
        std::set<std::pair<std::string, std::string>> allKeys;
        {
            std::set<std::string> allGenomeNames;
            for (const auto& [leafName, _] : seqpro_managers) allGenomeNames.insert(leafName);
            for (const auto& anc : ancestor_nodes) allGenomeNames.insert(anc.node_name);

            for (const auto& genomeName : allGenomeNames) {
                hal::Genome* g = alignment->openGenome(genomeName);
                if (!g) continue;

                const bool needTop = (g->getParent() != nullptr);
                const bool needBottom = (g->getNumChildren() > 0);

                for (auto seqIt = g->getSequenceIterator(0); !seqIt->atEnd(); seqIt->toNext()) {
                    hal::Sequence* s = seqIt->getSequence();
                    if (!s) continue;
                    auto key = std::make_pair(genomeName, s->getName());
                    seqLengths[key] = s->getSequenceLength();
                    if (needTop) keysNeedTop.insert(key);
                    if (needBottom) keysNeedBottom.insert(key);
                    if (needTop || needBottom) allKeys.insert(key);
                }

                alignment->closeGenome(g);
            }
        }

        std::mutex bottom_full_mutex, top_full_mutex;
        std::vector<std::future<void>> pass23_futures;
        pass23_futures.reserve(allKeys.size());
        for (const auto& key : allKeys) {
            pass23_futures.emplace_back(pool.enqueue([&, key]() {
                std::set<hal_index_t> uni_bp;
                if (auto it = bpBottom.find(key); it != bpBottom.end()) {
                    uni_bp.insert(it->second.begin(), it->second.end());
                }
                if (auto it = bpTop.find(key); it != bpTop.end()) {
                    uni_bp.insert(it->second.begin(), it->second.end());
                }

                const auto lenIt = seqLengths.find(key);
                if (lenIt == seqLengths.end()) {
                    throw std::runtime_error("HAL构建失败：无法获取 sequence 长度（内部错误）：genome=" + key.first + " seq=" + key.second);
                }
                const hal_index_t seqLen = static_cast<hal_index_t>(lenIt->second);
                uni_bp.insert(0);
                uni_bp.insert(seqLen);

                if (uni_bp.size() < 2) return;

                std::vector<SimpleSegmentInfo> segs;
                segs.reserve(uni_bp.size());
                hal_index_t prev = -1;
                for (hal_index_t x : uni_bp) {
                    if (prev != -1 && x > prev) segs.push_back({prev, static_cast<hal_size_t>(x - prev)});
                    prev = x;
                }

                if (keysNeedBottom.count(key)) {
                    std::lock_guard<std::mutex> lk(bottom_full_mutex);
                    bottomSegmentsFull[key] = segs;
                }
                if (keysNeedTop.count(key)) {
                    std::lock_guard<std::mutex> lk(top_full_mutex);
                    topSegmentsFull[key] = segs;
                }
            }));
        }
        pool.waitAllTasksDone();
        for (auto& f : pass23_futures) {
            f.get(); // 传播 worker 中的异常，保证科研软件 fail-fast
        }
        auto t23_end = __now();
        spdlog::info("  Pass2.3: built full segments in {} ms (bottom groups: {}, top groups: {})",
                      __ms(t23_start, t23_end), bottomSegmentsFull.size(), topSegmentsFull.size());

        size_t bottomFullGroups = 0, topFullGroups = 0;
        size_t bottomFullCount = 0, topFullCount = 0;
        for (const auto& [key, segs] : bottomSegmentsFull) { bottomFullGroups++; bottomFullCount += segs.size(); }
        for (const auto& [key, segs] : topSegmentsFull) { topFullGroups++; topFullCount += segs.size(); }

        spdlog::info("Full segments (with gaps) ready: bottom groups = {}, top groups = {}, total bottom segs = {}, total top segs = {}",
                      bottomFullGroups, topFullGroups, bottomFullCount, topFullCount);

        // 使用 setDimensions 一次性为每个基因组设置所有序列的 top/bottom 维度，避免区间重叠
        auto t_dim_start = __now();
        spdlog::info("Phase 3 (pass 1): Updating segment dimensions (no DNA writes)...");

        // 需要更新的基因组集合
        std::set<std::string> genomesToUpdate;
        for (const auto& [gc, _] : bottomSegmentsFull) genomesToUpdate.insert(gc.first);
        for (const auto& [gc, _] : topSegmentsFull) genomesToUpdate.insert(gc.first);

        size_t genomesUpdated = 0;
        for (const auto& genomeName : genomesToUpdate) {
            hal::Genome* genome = alignment->openGenome(genomeName);
            if (!genome) {
                spdlog::warn("Cannot open genome '{}' for updateDimensions", genomeName);
                continue;
            }

            std::vector<hal::Sequence::UpdateInfo> topUpdates;
            std::vector<hal::Sequence::UpdateInfo> bottomUpdates;

            for (auto seqIt = genome->getSequenceIterator(0); !seqIt->atEnd(); seqIt->toNext()) {
                hal::Sequence* seq = seqIt->getSequence();
                const std::string chrName = seq->getName();

                if (auto it = topSegmentsFull.find(std::make_pair(genomeName, chrName)); it != topSegmentsFull.end()) {
                    topUpdates.emplace_back(chrName, static_cast<hal_size_t>(it->second.size()));
                }
                if (auto it = bottomSegmentsFull.find(std::make_pair(genomeName, chrName)); it != bottomSegmentsFull.end()) {
                    bottomUpdates.emplace_back(chrName, static_cast<hal_size_t>(it->second.size()));
                }
            }

            if (!topUpdates.empty()) genome->updateTopDimensions(topUpdates);
            if (!bottomUpdates.empty()) genome->updateBottomDimensions(bottomUpdates);

            alignment->closeGenome(genome);
            genomesUpdated++;
        }

        auto t_dim_end = __now();
        spdlog::info("  Pass3.dim.update: updated top/bottom dims for {} genomes in {} ms", genomesUpdated, __ms(t_dim_start, t_dim_end));

        // 第二遍：从包含缝隙段的完整集合写入 HAL
        auto t_write_start = __now();
        spdlog::info("Phase 3 (pass 2): Writing segments from full (gap-filled) sets...");

        using SegmentInfo = SimpleSegmentInfo; // 直接复用已排序的 SimpleSegmentInfo
        std::map<std::pair<std::string, std::string>, std::vector<SegmentInfo>> bottomSegments = bottomSegmentsFull;
        std::map<std::pair<std::string, std::string>, std::vector<SegmentInfo>> topSegments = topSegmentsFull;

        // 2. 预先打开所有需要的基因组，避免频繁开启/关闭
        std::set<std::string> allGenomes;
        for (const auto& [genomeChr, segments] : bottomSegments) {
            allGenomes.insert(genomeChr.first);
        }
        for (const auto& [genomeChr, segments] : topSegments) {
            allGenomes.insert(genomeChr.first);
        }

        std::map<std::string, hal::Genome*> openGenomes;
        for (const std::string& genomeName : allGenomes) {
            hal::Genome* genome = alignment->openGenome(genomeName);
            if (genome) {
                openGenomes[genomeName] = genome;
            } else {
                spdlog::warn("Failed to open genome '{}'", genomeName);
            }
        }

        // 3. 建立段索引映射表（分开存储 top 与 bottom，避免键冲突）
        using SegmentKey = std::tuple<std::string, std::string, hal_index_t, hal_size_t>;
        std::unordered_map<SegmentKey, hal_index_t> bottomIndexMap;
        std::unordered_map<SegmentKey, hal_index_t> topIndexMap;

        // 4. 批量写入所有 BottomSegments
        auto t_write_bottom_start = __now();
        spdlog::info("Writing {} parent genome bottom segments...", bottomSegments.size());
        for (auto& [genomeChr, segments] : bottomSegments) {
            const std::string& genomeName = genomeChr.first;
            const std::string& chrName = genomeChr.second;

            std::sort(segments.begin(), segments.end());

            auto genomeIt = openGenomes.find(genomeName);
            if (genomeIt == openGenomes.end()) {
                spdlog::warn("Genome '{}' not found in opened genomes", genomeName);
                continue;
            }
            hal::Genome* genome = genomeIt->second;

            hal::Sequence* seq = genome->getSequence(chrName);
            if (!seq) {
                spdlog::warn("Genome '{}' has no sequence '{}' for bottom segments", genomeName, chrName);
                continue;
            }

            const hal_index_t seqBottomStart = seq->getBottomSegmentArrayIndex();
            auto botIt = genome->getBottomSegmentIterator(seqBottomStart);
            const hal_index_t seqGenomeStart = seq->getStartPosition();
            for (const auto& seg : segments) {
                auto* bs = botIt->getBottomSegment();
                bs->setCoordinates(seqGenomeStart + seg.start, seg.length);
                hal_size_t numChildren = bs->getNumChildren();
                for (hal_size_t i = 0; i < numChildren; ++i) bs->setChildIndex(i, hal::NULL_INDEX);
                    bs->setTopParseIndex(hal::NULL_INDEX);
                SegmentKey key = std::make_tuple(genomeName, chrName, seqGenomeStart + seg.start, seg.length);
                bottomIndexMap[key] = botIt->getArrayIndex();
                botIt->toRight();
            }
            // spdlog::debug("Written {} bottom segments for {}:{}", segments.size(), genomeName, chrName);
        }
        auto t_write_bottom_end = __now();
        spdlog::info("  Pass3.writeBottom: {} ms", __ms(t_write_bottom_start, t_write_bottom_end));

        // 5. 批量写入所有 TopSegments
        auto t_write_top_start = __now();
        spdlog::info("Writing {} child genome top segments...", topSegments.size());
        for (auto& [genomeChr, segments] : topSegments) {
            const std::string& genomeName = genomeChr.first;
            const std::string& chrName = genomeChr.second;

            std::sort(segments.begin(), segments.end());

            auto genomeIt = openGenomes.find(genomeName);
            if (genomeIt == openGenomes.end()) {
                spdlog::warn("Genome '{}' not found in opened genomes", genomeName);
                continue;
            }
            hal::Genome* genome = genomeIt->second;

            hal::Sequence* seq = genome->getSequence(chrName);
            if (!seq) {
                spdlog::warn("Genome '{}' has no sequence '{}' for top segments", genomeName, chrName);
                continue;
            }

            const hal_index_t seqTopStart = seq->getTopSegmentArrayIndex();
            auto topIt = genome->getTopSegmentIterator(seqTopStart);
            const hal_index_t seqGenomeStart = seq->getStartPosition();
            for (const auto& seg : segments) {
                auto* ts = topIt->getTopSegment();
                ts->setCoordinates(seqGenomeStart + seg.start, seg.length);
                ts->setParentIndex(hal::NULL_INDEX);
                ts->setParentReversed(false);
                ts->setNextParalogyIndex(hal::NULL_INDEX);
                ts->setBottomParseIndex(hal::NULL_INDEX);
                SegmentKey key = std::make_tuple(genomeName, chrName, seqGenomeStart + seg.start, seg.length);
                topIndexMap[key] = topIt->getArrayIndex();
                topIt->toRight();
            }

            // fail-fast: 写入后立即校验该 sequence 的 top segments 覆盖长度必须等于 sequence 长度（HAL 强制）
            // 注意：HDF5TopSegment 的坐标采用“边界数组”存储，若存在断档/越界，会导致前一段长度被错误拉长，从而在后续建链阶段表现为 length mismatch。
            {
                const hal_size_t expected_len = seq->getSequenceLength();
                hal_size_t total_len = 0;
                auto verifyIt = seq->getTopSegmentIterator();
                for (hal_size_t i = 0; i < seq->getNumTopSegments(); ++i) {
                    total_len += verifyIt->getTopSegment()->getLength();
                    verifyIt->toRight();
                }
                if (total_len != expected_len) {
                    throw std::runtime_error(
                        "HAL构建失败：写入 top segments 后 sequence 覆盖长度不等于 sequence 长度（会导致后续 parent/child 长度异常）：genome=" +
                        genomeName + " seq=" + chrName +
                        " seq_len=" + std::to_string(expected_len) +
                        " top_total_len=" + std::to_string(total_len) +
                        " numTopSegments=" + std::to_string(seq->getNumTopSegments()));
                }
            }
            // spdlog::debug("Written {} top segments for {}:{}", segments.size(), genomeName, chrName);
        }
        auto t_write_top_end = __now();
        spdlog::info("  Pass3.writeTop: {} ms", __ms(t_write_top_start, t_write_top_end));

        // 6.5 重建各基因组的 parse 信息（库内批量实现）
        auto t_fixparse_start = __now();
        spdlog::info("Phase 3 (pass 2.5): Fixing parse info per genome...");
        for (auto& [genomeName, genome] : openGenomes) {
            genome->fixParseInfo();
        }
        auto t_fixparse_end = __now();
        spdlog::info("  Pass3.fixParseInfo: {} ms", __ms(t_fixparse_start, t_fixparse_end));

        // 6. 建立父子映射关系 (re-organized for data locality)
        auto t_map_start = __now();
        spdlog::info("Phase 3 (pass 3): Establishing parent-child segment mappings...");

        // 压缩版链接结构，避免复制字符串，显著降低内存
        struct LinkIdx {
            uint32_t parent_gid;
            uint32_t child_gid;
            hal_index_t parent_bidx;  // bottom segment array index
            hal_index_t child_tidx;   // top segment array index
            hal_index_t child_idx_in_parent; // parentGenome->getChildIndex(childGenome)
            bool is_reversed;
        };

	
	        std::vector<LinkIdx> all_links;
        size_t totalMappings = 0;
        size_t totalChildren = 0;
        for (const auto& [block, mappings] : refined_mappings) {
            for (const auto& m : mappings) {
                totalMappings++;
                totalChildren += m.children.size();
            }
        }
        all_links.reserve(totalChildren);

        // 为 openGenomes 分配整数 id
        std::unordered_map<std::string, uint32_t> genomeId;
        genomeId.reserve(openGenomes.size());
        uint32_t gidCounter = 0;
        for (const auto& [name, _gen] : openGenomes) genomeId.emplace(name, gidCounter++);

        // 预缓存 parent->child 的 childIndex
        std::unordered_map<uint64_t, hal_index_t> parentChildIdx; // key = (parent_gid<<32) | child_gid
        parentChildIdx.reserve(openGenomes.size() * 4);
        auto keyPC = [](uint32_t pg, uint32_t cg) -> uint64_t { return (uint64_t(pg) << 32) | uint64_t(cg); };
        for (const auto& [pname, pgen] : openGenomes) {
            uint32_t pg = genomeId[pname];
            for (const auto& [cname, cgen] : openGenomes) {
                hal_index_t idx = pgen->getChildIndex(cgen);
                if (idx != hal::NULL_INDEX) parentChildIdx.emplace(keyPC(pg, genomeId[cname]), idx);
            }
        }

        // 生成压缩链接，直接解析成段数组索引，避免保存字符串键
        for (const auto& [block, mappings] : refined_mappings) {
            std::string blk_ref_chr;
            if (block) {
                std::shared_lock blk_lock(block->rw);
                blk_ref_chr = block->ref_chr;
            }
            for (const auto& mapping : mappings) {
                auto pgIt = openGenomes.find(mapping.parent_genome);
                if (pgIt == openGenomes.end()) {
                    throw std::runtime_error(
                        "HAL构建失败(ref_chr=" + blk_ref_chr + ")：parent genome 未打开/不存在：parent=" + mapping.parent_genome);
                }
                hal::Genome* parentGenome = pgIt->second;
                hal::Sequence* parentSeq = parentGenome->getSequence(mapping.parent_chr_name);
                if (!parentSeq) {
                    throw std::runtime_error(
                        "HAL构建失败(ref_chr=" + blk_ref_chr + ")：parent sequence 不存在：parent=" + mapping.parent_genome + " chr=" +
                        mapping.parent_chr_name);
                }
                hal_index_t pSeqStart = parentSeq->getStartPosition();
                SegmentKey pKey = std::make_tuple(mapping.parent_genome, mapping.parent_chr_name, pSeqStart + (hal_index_t)mapping.parent_start, mapping.parent_length);
                auto pFound = bottomIndexMap.find(pKey);
                if (pFound == bottomIndexMap.end()) {
                    throw std::runtime_error(
                        "HAL构建失败(ref_chr=" + blk_ref_chr + ")：未找到 parent bottom segment 索引：parent=" + mapping.parent_genome + ":" +
                        mapping.parent_chr_name + " start=" + std::to_string(mapping.parent_start) +
                        " len=" + std::to_string(mapping.parent_length));
                }
                hal_index_t pIdx = pFound->second;
                uint32_t pg = genomeId[mapping.parent_genome];

                // 强校验：bottomIndexMap 返回的索引必须与 key 对应的段一致（否则说明段写入/索引映射已损坏）
                {
                    auto pIt = parentGenome->getBottomSegmentIterator(pIdx);
                    auto* pSeg = pIt->getBottomSegment();
                    const hal_index_t expect_start = pSeqStart + (hal_index_t)mapping.parent_start;
                    if (pSeg->getStartPosition() != expect_start || pSeg->getLength() != mapping.parent_length ||
                        pSeg->getSequence()->getName() != mapping.parent_chr_name) {
                        throw std::runtime_error(
                            "HAL构建失败(ref_chr=" + blk_ref_chr +
                            ")：parent bottom segment 索引映射不一致（常见原因：某些sequence未生成bottom segments导致边界坐标被覆盖，从而把上一条sequence末段错误拉长；也可能是写入越界/并发导致的内存破坏）：parent=" +
                            mapping.parent_genome + " key_seq=" + mapping.parent_chr_name +
                            " parent_seqStart=" + std::to_string(pSeqStart) +
                            " parent_seqLen=" + std::to_string(parentSeq->getSequenceLength()) +
                            " key_start_in_seq=" + std::to_string(mapping.parent_start) +
                            " key_start=" + std::to_string(expect_start) +
                            " key_len=" + std::to_string(mapping.parent_length) +
                            " -> actual_seq=" + pSeg->getSequence()->getName() +
                            " actual_start=" + std::to_string(pSeg->getStartPosition()) +
                            " actual_len=" + std::to_string(pSeg->getLength()) +
                            " parent_bidx=" + std::to_string(pIdx));
                    }
                }

                for (const auto& child : mapping.children) {
                    if (child.child_length != mapping.parent_length) {
                        throw std::runtime_error(
                            "HAL构建失败(ref_chr=" + blk_ref_chr + ")：parent/child 长度不一致（HAL要求必须等长）：parent=" +
                            mapping.parent_genome + ":" + mapping.parent_chr_name +
                            " start=" + std::to_string(mapping.parent_start) +
                            " len=" + std::to_string(mapping.parent_length) +
                            " -> child=" + child.child_genome + ":" + child.child_chr_name +
                            " start=" + std::to_string(child.child_start) +
                            " len=" + std::to_string(child.child_length));
                    }
                    auto cgIt = openGenomes.find(child.child_genome);
                    if (cgIt == openGenomes.end()) {
                        throw std::runtime_error(
                            "HAL构建失败(ref_chr=" + blk_ref_chr + ")：child genome 未打开/不存在：child=" + child.child_genome);
                    }
                    hal::Genome* childGenome = cgIt->second;
                    hal::Sequence* childSeq = childGenome->getSequence(child.child_chr_name);
                    if (!childSeq) {
                        throw std::runtime_error(
                            "HAL构建失败(ref_chr=" + blk_ref_chr + ")：child sequence 不存在：child=" + child.child_genome + " chr=" +
                            child.child_chr_name);
                    }
                    hal_index_t cSeqStart = childSeq->getStartPosition();
                    SegmentKey cKey = std::make_tuple(child.child_genome, child.child_chr_name, cSeqStart + (hal_index_t)child.child_start, child.child_length);
                    auto cFound = topIndexMap.find(cKey);
                    if (cFound == topIndexMap.end()) {
                        throw std::runtime_error(
                            "HAL构建失败(ref_chr=" + blk_ref_chr + ")：未找到 child top segment 索引：child=" + child.child_genome + ":" +
                            child.child_chr_name + " start=" + std::to_string(child.child_start) +
                            " len=" + std::to_string(child.child_length));
                    }
                    hal_index_t cIdx = cFound->second;

                    // 强校验：topIndexMap 返回的索引必须与 key 对应的段一致（否则说明段写入/索引映射已损坏）
                    {
                        auto cIt = childGenome->getTopSegmentIterator(cIdx);
                        auto* cSeg = cIt->getTopSegment();
                        const hal_index_t expect_start = cSeqStart + (hal_index_t)child.child_start;
                        if (cSeg->getStartPosition() != expect_start || cSeg->getLength() != child.child_length ||
                            cSeg->getSequence()->getName() != child.child_chr_name) {
                            const hal_index_t expect_end = expect_start + static_cast<hal_index_t>(child.child_length);
                            const hal_index_t actual_end = cSeg->getStartPosition() + static_cast<hal_index_t>(cSeg->getLength());
                            throw std::runtime_error(
                                "HAL构建失败(ref_chr=" + blk_ref_chr +
                                ")：child top segment 索引映射不一致（常见原因：某些sequence未生成top segments导致边界坐标被覆盖，从而把上一条sequence末段错误拉长；也可能是写入越界/并发导致的内存破坏）：child=" +
                                child.child_genome + " key_seq=" + child.child_chr_name +
                                " child_seqStart=" + std::to_string(cSeqStart) +
                                " child_seqLen=" + std::to_string(childSeq->getSequenceLength()) +
                                " key_start_in_seq=" + std::to_string(child.child_start) +
                                " key_start=" + std::to_string(expect_start) +
                                " key_len=" + std::to_string(child.child_length) +
                                " key_end=" + std::to_string(expect_end) +
                                " -> actual_seq=" + cSeg->getSequence()->getName() +
                                " actual_start=" + std::to_string(cSeg->getStartPosition()) +
                                " actual_len=" + std::to_string(cSeg->getLength()) +
                                " actual_end=" + std::to_string(actual_end) +
                                " child_tidx=" + std::to_string(cIdx));
                        }
                    }

                    uint32_t cg = genomeId[child.child_genome];
                    auto pci = parentChildIdx.find(keyPC(pg, cg));
                    if (pci == parentChildIdx.end()) {
                        throw std::runtime_error(
                            "HAL构建失败(ref_chr=" + blk_ref_chr + ")：parent genome 不包含该 child genome（系统发育树/基因组关系不一致）：parent=" +
                            mapping.parent_genome + " child=" + child.child_genome);
                    }

                    all_links.push_back(LinkIdx{pg, cg, pIdx, cIdx, pci->second, child.is_reversed});
                }
            }
        }

	        // 反向映射 gid -> genome*
	        std::vector<hal::Genome*> gid2genome(genomeId.size(), nullptr);
	        for (const auto& [name, g] : openGenomes) gid2genome[genomeId[name]] = g;
	
	        // Pass 3: Link Parent <-> Child，并处理 duplication（nextParalogy）
	        // 说明：HAL 要求 parent bottom 与 child top 必须等长；并且若同一 child genome 中有多个 top segments 指向同一 parentIndex，
	        //      必须设置 nextParalogyIndex（halValidate::validateDuplications 会强制检查）。
	        auto t_link_start = __now();
	        auto group_sorter = [](const LinkIdx& a, const LinkIdx& b) {
	            if (a.parent_gid != b.parent_gid) return a.parent_gid < b.parent_gid;
	            if (a.child_gid != b.child_gid) return a.child_gid < b.child_gid;
	            if (a.parent_bidx != b.parent_bidx) return a.parent_bidx < b.parent_bidx;
	            return a.child_tidx < b.child_tidx;
	        };
	        std::sort(all_links.begin(), all_links.end(), group_sorter);
	
	        size_t successful_p2c = 0;
	        size_t successful_c2p = 0;
	
	        for (size_t i = 0; i < all_links.size();) {
	            const LinkIdx& first = all_links[i];
	            size_t j = i + 1;
	
	            const bool is_reversed = first.is_reversed;
	            std::vector<hal_index_t> child_tidxs;
	            child_tidxs.push_back(first.child_tidx);
	
	            while (j < all_links.size()) {
	                const auto& cur = all_links[j];
	                if (cur.parent_gid != first.parent_gid ||
	                    cur.child_gid != first.child_gid ||
	                    cur.parent_bidx != first.parent_bidx ||
	                    cur.child_idx_in_parent != first.child_idx_in_parent) {
	                    break;
	                }
	                if (cur.is_reversed != is_reversed) {
	                    throw std::runtime_error(
	                        "HAL构建失败：同一 parent/child/parentSeg 上出现不同方向的 mapping（无法安全表示）：parent_gid=" +
	                        std::to_string(first.parent_gid) + " child_gid=" + std::to_string(first.child_gid) +
	                        " parent_bidx=" + std::to_string(first.parent_bidx) +
	                        " child_slot=" + std::to_string(first.child_idx_in_parent));
	                }
	                child_tidxs.push_back(cur.child_tidx);
	                ++j;
	            }
	
	            std::sort(child_tidxs.begin(), child_tidxs.end());
	            child_tidxs.erase(std::unique(child_tidxs.begin(), child_tidxs.end()), child_tidxs.end());
	
	            hal::Genome* parentGenome = gid2genome[first.parent_gid];
	            hal::Genome* childGenome = gid2genome[first.child_gid];
	            if (!parentGenome || !childGenome) {
	                throw std::runtime_error("HAL构建失败：gid->genome 映射为空（内部错误）");
	            }
	
	            auto parentBottomIt = parentGenome->getBottomSegmentIterator(first.parent_bidx);
	            auto* parentBottomSeg = parentBottomIt->getBottomSegment();
	            if (first.child_idx_in_parent >= parentBottomSeg->getNumChildren()) {
	                throw std::runtime_error(
	                    "HAL构建失败：child slot 超出 parent bottom segment children 数量：parent=" + parentGenome->getName() +
	                    " parent_bidx=" + std::to_string(first.parent_bidx) +
	                    " child_slot=" + std::to_string(first.child_idx_in_parent) +
	                    " numChildren=" + std::to_string(parentBottomSeg->getNumChildren()));
	            }
	
	            // bottom: childIndex 只保存一个入口；若存在 duplication，依赖 nextParalogy 环遍历
	            const hal_index_t rep_tidx = child_tidxs.front();
	            const hal_index_t existing_child = parentBottomSeg->getChildIndex(first.child_idx_in_parent);
	            if (existing_child == hal::NULL_INDEX) {
	                parentBottomSeg->setChildIndex(first.child_idx_in_parent, rep_tidx);
	                parentBottomSeg->setChildReversed(first.child_idx_in_parent, is_reversed);
	                successful_p2c++;
	            } else {
	                if (!std::binary_search(child_tidxs.begin(), child_tidxs.end(), existing_child)) {
	                    throw std::runtime_error(
	                        "HAL构建失败：同一 parent bottom segment 的 childIndex 被不同 top segment 冲突占用（会导致输出错误）：parent=" +
	                        parentGenome->getName() + " child=" + childGenome->getName() +
	                        " parent_bidx=" + std::to_string(parentBottomSeg->getArrayIndex()) +
	                        " child_slot=" + std::to_string(first.child_idx_in_parent) +
	                        " existing_child_tidx=" + std::to_string(existing_child) +
	                        " new_rep_child_tidx=" + std::to_string(rep_tidx));
	                }
	                if (parentBottomSeg->getChildReversed(first.child_idx_in_parent) != is_reversed) {
	                    throw std::runtime_error(
	                        "HAL构建失败：同一 parent bottom segment 上 childReversed 冲突（会导致输出错误）：parent=" +
	                        parentGenome->getName() + " child=" + childGenome->getName() +
	                        " parent_bidx=" + std::to_string(parentBottomSeg->getArrayIndex()) +
	                        " child_slot=" + std::to_string(first.child_idx_in_parent));
	                }
	            }
	
	            // top: 设置 parentIndex / parentReversed，并校验等长（HAL 强制要求）
	            for (const hal_index_t tidx : child_tidxs) {
	                auto childTopIt = childGenome->getTopSegmentIterator(tidx);
	                auto* childTopSeg = childTopIt->getTopSegment();
	
	                if (childTopSeg->getLength() != parentBottomSeg->getLength()) {
	                    throw std::runtime_error(
	                        "HAL构建失败：parent/child segment 长度不一致（HAL 不允许）：parent=" + parentGenome->getName() +
	                        " parent_bidx=" + std::to_string(parentBottomSeg->getArrayIndex()) +
	                        " parent_start=" + std::to_string(parentBottomSeg->getStartPosition()) +
	                        " parent_len=" + std::to_string(parentBottomSeg->getLength()) +
	                        " parent_seq=" + parentBottomSeg->getSequence()->getName() +
	                        " child=" + childGenome->getName() +
	                        " child_tidx=" + std::to_string(childTopSeg->getArrayIndex()) +
	                        " child_start=" + std::to_string(childTopSeg->getStartPosition()) +
	                        " child_len=" + std::to_string(childTopSeg->getLength()) +
	                        " child_seq=" + childTopSeg->getSequence()->getName());
	                }
	
	                const hal_index_t existing_parent = childTopSeg->getParentIndex();
	                if (existing_parent != hal::NULL_INDEX && existing_parent != parentBottomSeg->getArrayIndex()) {
	                    throw std::runtime_error(
	                        "HAL构建失败：同一 child top segment 被多个不同 parent bottom segment 复用（当前实现不支持；需先解决数据冲突或实现更复杂的表示）：child=" +
	                        childGenome->getName() + " child_tidx=" + std::to_string(childTopSeg->getArrayIndex()) +
	                        " child_seq=" + childTopSeg->getSequence()->getName() +
	                        " child_start=" + std::to_string(childTopSeg->getStartPosition()) +
	                        " existing_parent_bidx=" + std::to_string(existing_parent) +
	                        " new_parent_bidx=" + std::to_string(parentBottomSeg->getArrayIndex()));
	                }
	
	                childTopSeg->setParentIndex(parentBottomSeg->getArrayIndex());
	                childTopSeg->setParentReversed(is_reversed);
	                successful_c2p++;
	            }
	
	            // duplication: nextParalogyIndex 需要对共享同一 parentIndex 的 top segments 设为非 NULL（这里用环形链表）
	            if (child_tidxs.size() == 1) {
	                auto tIt = childGenome->getTopSegmentIterator(child_tidxs[0]);
	                tIt->getTopSegment()->setNextParalogyIndex(hal::NULL_INDEX);
	            } else {
	                for (size_t k = 0; k < child_tidxs.size(); ++k) {
	                    const hal_index_t cur = child_tidxs[k];
	                    const hal_index_t nxt = child_tidxs[(k + 1) % child_tidxs.size()];
	                    auto tIt = childGenome->getTopSegmentIterator(cur);
	                    tIt->getTopSegment()->setNextParalogyIndex(nxt);
	                }
	            }
	
	            i = j;
	        }
	
	        auto t_link_end = __now();
	        spdlog::info("  Pass3.parentChild.p2c: linked {} parent->child in {} ms", successful_p2c, __ms(t_link_start, t_link_end));
	        spdlog::info("  Pass3.parentChild.c2p: linked {} child->parent in {} ms", successful_c2p, __ms(t_link_start, t_link_end));
	
	        size_t successfulMappings = std::min(successful_p2c, successful_c2p);
	
	        auto t_map_end = __now();
	        spdlog::info("  Pass3.parentChild: mapped {}/{} links in {} ms (total)", successfulMappings, totalChildren, __ms(t_map_start, t_map_end));

        // 关闭所有打开的基因组
        auto t_close_start = __now();
        for (auto& [genomeName, genome] : openGenomes) {
            alignment->closeGenome(genome);
            // spdlog::debug("Closed genome '{}'", genomeName);
        }
        openGenomes.clear();
        auto t_close_end = __now();
        spdlog::info("  Pass3.closeGenomes: {} ms", __ms(t_close_start, t_close_end));

        auto t_total_end = __now();
        spdlog::info("Phase 3 completed in {} ms. Bottom segs: {}, Top segs: {}",
                      __ms(t_total_start, t_total_end), bottomIndexMap.size(), topIndexMap.size());

    }

} // namespace hal_converter
} // namespace RaMesh
