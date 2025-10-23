#include <sequence_utils.h>
#include <algorithm>

#include "rare_aligner.h"
#include "anchor.h"  // 包含 UnionFind 定义
#include "SeqPro.h"  // 包含 SeqPro 相关定义
#include "ramesh.h"  // 包含 RaMesh 图结构定义

// 辅助函数：根据CIGAR字符串计算query区间对应关系
namespace {

    /**
     * @brief 获取query segment在ref上的映射区间
     * @param query_segment query segment对象
     * @param ref_block_ptr segment所属的block，用于获取ref anchor
     * @param ref_name 参考物种名称
     * @return pair<ref_start, ref_end> 映射到ref上的区间，如果无法映射返回{0,0}
     */
    std::pair<uint_t, uint_t> getRefMappedInterval(
        const RaMesh::SegPtr& query_segment,
        const RaMesh::BlockPtr& ref_block_ptr,
        const SpeciesName& ref_name) {

        if (!query_segment || !ref_block_ptr) {
            return {0, 0};
        }

        // 查找对应的ref segment
        std::shared_lock block_lock(ref_block_ptr->rw);

        // 在block的anchors中查找ref segment
        for (const auto& [species_chr_pair, segment] : ref_block_ptr->anchors) {
            const auto& [species_name, chr_name] = species_chr_pair;
            if (species_name == ref_name && segment && segment->isSegment()) {
                // 找到了ref segment，返回其区间
                return {segment->start, segment->start + segment->length};
            }
        }

        spdlog::warn("[getRefMappedInterval] Cannot find corresponding ref segment");
        return {0, 0};
    }

    /**
     * @brief 根据CIGAR字符串计算ref区间对应的query区间
     * @param cigar 原始segment的CIGAR字符串
     * @param target_ref_start 目标ref区间的起始位置
     * @param target_ref_end 目标ref区间的结束位置
     * @param original_ref_start 原始segment在ref上的起始位置
     * @param original_query_start 原始segment在query上的起始位置
     * @return pair<query_start, query_length> 对应的query区间
     */
    std::pair<uint_t, uint_t> calculateQueryInterval(
        const Cigar_t& cigar,
        uint_t target_ref_start, uint_t target_ref_end,
        uint_t original_ref_start, uint_t original_query_start) {

        // 如果CIGAR为空，使用简单的线性映射
        if (cigar.empty()) {
            uint_t ref_offset = target_ref_start - original_ref_start;
            uint_t length = target_ref_end - target_ref_start;
            return {original_query_start + ref_offset, length};
        }

        uint_t current_ref_pos = original_ref_start;
        uint_t current_query_pos = original_query_start;
        uint_t target_query_start = 0;
        uint_t target_query_end = 0;
        bool found_start = false;
        bool found_end = false;

        // 遍历CIGAR操作
        for (const auto& cigar_unit : cigar) {
            uint_t length = cigar_unit >> 4;  // 高28位是长度
            uint_t op_code = cigar_unit & 0xF; // 低4位是操作码

            switch (op_code) {
                case 0x0: // M - match/mismatch
                case 0x7: // = - exact match
                case 0x8: // X - mismatch
                {
                    // 这些操作同时消耗ref和query位置
                    uint_t ref_end_this_op = current_ref_pos + length;

                    // 检查target_ref_start是否在这个操作范围内
                    if (!found_start && current_ref_pos <= target_ref_start && target_ref_start < ref_end_this_op) {
                        uint_t offset = target_ref_start - current_ref_pos;
                        target_query_start = current_query_pos + offset;
                        found_start = true;
                    }

                    // 检查target_ref_end是否在这个操作范围内
                    if (!found_end && current_ref_pos <= target_ref_end && target_ref_end <= ref_end_this_op) {
                        uint_t offset = target_ref_end - current_ref_pos;
                        target_query_end = current_query_pos + offset;
                        found_end = true;
                        break;
                    }

                    current_ref_pos += length;
                    current_query_pos += length;
                    break;
                }
                case 0x1: // I - insertion (query only)
                {
                    // 插入操作只消耗query位置，不影响ref坐标映射
                    current_query_pos += length;
                    break;
                }
                case 0x2: // D - deletion (ref only)
                {
                    // 删除操作只消耗ref位置
                    uint_t ref_end_this_op = current_ref_pos + length;

                    // 检查删除区间是否包含我们的目标区间
                    if (!found_start && current_ref_pos <= target_ref_start && target_ref_start < ref_end_this_op) {
                        target_query_start = current_query_pos;
                        found_start = true;
                    }

                    if (!found_end && current_ref_pos <= target_ref_end && target_ref_end <= ref_end_this_op) {
                        target_query_end = current_query_pos;
                        found_end = true;
                        break;
                    }

                    current_ref_pos += length;
                    // query位置不变
                    break;
                }
                default:
                    // 其他CIGAR操作暂时按match处理
                    current_ref_pos += length;
                    current_query_pos += length;
                    break;
            }

            if (found_start && found_end) break;
        }

        // 如果没有找到精确的映射，使用线性近似
        if (!found_start || !found_end) {
            spdlog::warn("[calculateQueryInterval] Cannot precisely map CIGAR interval, using linear approximation");
            uint_t ref_offset = target_ref_start - original_ref_start;
            uint_t length = target_ref_end - target_ref_start;
            return {original_query_start + ref_offset, length};
        }

        // 如果出现 start > end（可能来源于反向链坐标系），进行交换校正
        if (target_query_start > target_query_end) {
            std::swap(target_query_start, target_query_end);
        }

        uint_t query_length = target_query_end > target_query_start ? (target_query_end - target_query_start) : 0;
        return {target_query_start, query_length};
    }

    /**
     * @brief 创建分割后的segment CIGAR字符串
     * @param original_cigar 原始CIGAR
     * @param ref_start_offset ref区间在原始segment中的偏移
     * @param ref_length ref区间长度
     * @return 分割后的CIGAR字符串
     */
    Cigar_t extractCigarSubsequence(
        const Cigar_t& original_cigar,
        uint_t ref_start_offset, uint_t ref_length) {

        if (original_cigar.empty()) {
            // 如果原始CIGAR为空，创建一个简单的match操作
            return Cigar_t{cigarToInt('M', ref_length)};
        }

        Cigar_t result_cigar;
        uint_t current_ref_pos = 0;
        uint_t target_ref_start = ref_start_offset;
        uint_t target_ref_end = ref_start_offset + ref_length;

        for (const auto& cigar_unit : original_cigar) {
            uint_t length = cigar_unit >> 4;
            uint_t op_code = cigar_unit & 0xF;

            switch (op_code) {
                case 0x0: case 0x7: case 0x8: case 0x2: // M, =, X, D - 消耗ref位置
                {
                    uint_t ref_end_this_op = current_ref_pos + length;

                    // 检查这个操作是否与目标区间重叠
                    if (current_ref_pos < target_ref_end && ref_end_this_op > target_ref_start) {
                        uint_t overlap_start = std::max(current_ref_pos, target_ref_start);
                        uint_t overlap_end = std::min(ref_end_this_op, target_ref_end);
                        uint_t overlap_length = overlap_end - overlap_start;

                        if (overlap_length > 0) {
                            result_cigar.push_back(cigarToInt(
                                (op_code == 0x0) ? 'M' :
                                (op_code == 0x7) ? '=' :
                                (op_code == 0x8) ? 'X' : 'D',
                                overlap_length));
                        }
                    }

                    current_ref_pos += length;
                    break;
                }
                case 0x1: // I - 插入操作，不消耗ref位置
                {
                    // 插入操作在ref坐标范围内时保留
                    if (current_ref_pos >= target_ref_start && current_ref_pos < target_ref_end) {
                        result_cigar.push_back(cigarToInt('I', length));
                    }
                    break;
                }
                default:
                    // 其他操作暂时跳过
                    break;
            }
        }

        // 如果结果为空，创建一个默认的match操作
        if (result_cigar.empty()) {
            result_cigar.push_back(cigarToInt('M', ref_length));
        }

        return result_cigar;
    }

    /**
     * @brief 修剪同一染色体上一组 segment 之间的重叠。
     *        输入必须已按 start 升序。算法：线性扫描，若发现 seg.start < prevEnd，
     *        则将 seg.start 调整为 prevEnd 并相应减少 length；若 length<=0 则丢弃。
     * @param segs  Segment 向量（已排序）。
     * @return true 表示处理后仍然存在重叠（异常）；false 表示已无重叠。
     */
    bool trimSegments(std::vector<RaMesh::SegPtr>& segs) {
        if (segs.empty()) return false;

        std::vector<RaMesh::SegPtr> trimmed;
        trimmed.reserve(segs.size());

        RaMesh::SegPtr prev = nullptr;
        for (auto seg : segs) {
            if (!prev) {
                trimmed.push_back(seg);
                prev = seg;
                continue;
            }

            uint_t prev_end = prev->start + prev->length;
            if (seg->start < prev_end) {
                uint_t overlap_len = prev_end - seg->start;
                if (seg->length <= overlap_len) {
                    // 完全被覆盖，跳过
                    continue;
                }
                seg->start = prev_end;
                seg->length -= overlap_len;
            }

            trimmed.push_back(seg);
            prev = seg;
        }

        segs.swap(trimmed);

        // 再次检查是否还有重叠
        for (size_t i = 1; i < segs.size(); ++i) {
            if (segs[i]->start < segs[i-1]->start + segs[i-1]->length) {
                return true; // 仍有重叠
            }
        }
        return false;
    }
} // anonymous namespace

/**
 * @brief 确保 SeqPro manager 是 MaskedSequenceManager 类型
 * @param manager_variant 当前的 manager variant
 * @return 指向 MaskedSequenceManager 的指针
 */
SeqPro::MaskedSequenceManager* ensureMaskedManager(SeqPro::SharedManagerVariant& manager_variant) {
    auto& variant = *manager_variant;
    
    // 检查当前类型
    if (std::holds_alternative<std::unique_ptr<SeqPro::MaskedSequenceManager>>(variant)) {
        // 已经是 MaskedSequenceManager
        return std::get<std::unique_ptr<SeqPro::MaskedSequenceManager>>(variant).get();
    } 
    else if (std::holds_alternative<std::unique_ptr<SeqPro::SequenceManager>>(variant)) {
        // 需要转换为 MaskedSequenceManager
        auto seq_manager = std::move(std::get<std::unique_ptr<SeqPro::SequenceManager>>(variant));
        auto masked_manager = std::make_unique<SeqPro::MaskedSequenceManager>(std::move(seq_manager));
        auto* result_ptr = masked_manager.get();
        
        // 替换 variant 中的内容
        variant = std::move(masked_manager);
        
        return result_ptr;
    }
    
    throw std::runtime_error("Invalid SeqPro manager variant type");
}

/**
 * @brief 从比对结果图中提取已比对的区间，并作为遮蔽区间添加到对应的 SeqPro manager 中
 * @param graph 比对结果图
 * @param seqpro_managers SeqPro manager 映射
 * @param ref_name 参考物种名称
 */
void addAlignedRegionsAsMask(
    const RaMesh::RaMeshMultiGenomeGraph& graph,
    std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
    const SpeciesName& ref_name) {
    
    if (graph.blocks.empty()) {
        spdlog::info("[addAlignedRegionsAsMask] No blocks to process for masking");
        return;
    }
    
    #ifdef _DEBUG_
    spdlog::info("[addAlignedRegionsAsMask] Extracting aligned regions as mask intervals from {} blocks", 
                 graph.blocks.size());
    #endif
    // 按物种和染色体分组收集区间
    std::unordered_map<SpeciesName, std::unordered_map<ChrName, std::vector<SeqPro::MaskInterval>>> 
        species_chr_intervals;
    
    size_t total_intervals = 0;
    size_t valid_blocks = 0;
    
    // 遍历所有 blocks，提取 segment 区间
    for (const auto& weak_block : graph.blocks) {
        auto block_ptr = weak_block.lock();
        if (!block_ptr) continue;
        
        valid_blocks++;
        std::shared_lock block_lock(block_ptr->rw);
        
        // 处理该 block 中的所有 anchors
        for (const auto& [species_chr_pair, segment] : block_ptr->anchors) {
            const auto& [species_name, chr_name] = species_chr_pair;
            
            // 只处理有效的 segment
            if (!segment || !segment->isSegment() || segment->length == 0) {
                continue;
            }
            
            // 检查该物种是否在 seqpro_managers 中
            if (seqpro_managers.find(species_name) == seqpro_managers.end()) {
                continue;
            }
            
            // 创建遮蔽区间（使用原始坐标）
            SeqPro::MaskInterval interval(segment->start, segment->start + segment->length);
            if (segment->start > 1000000) {
                std::cout << "";
            }
            species_chr_intervals[species_name][chr_name].push_back(interval);
            total_intervals++;
        }
    }
    #ifdef _DEBUG_
    spdlog::info("[addAlignedRegionsAsMask] Collected {} intervals from {} valid blocks across {} species", 
                 total_intervals, valid_blocks, species_chr_intervals.size());
    #endif
    // 为每个物种批量添加遮蔽区间
    for (auto& [species_name, chr_intervals] : species_chr_intervals) {
        try {
            // 确保该物种的 manager 是 MaskedSequenceManager
            auto* masked_manager = ensureMaskedManager(seqpro_managers[species_name]);
            
            size_t species_total_intervals = 0;
            
            // 按染色体处理区间
            for (auto& [chr_name, intervals] : chr_intervals) {
                if (intervals.empty()) continue;
                
                // 构造序列名（假设格式为染色体名）
                std::string seq_name = chr_name;
                
                                // 检查序列是否存在
                if (masked_manager->getSequenceId(seq_name) == SeqPro::SequenceIndex::INVALID_ID) {
                    spdlog::warn("[addAlignedRegionsAsMask] Sequence not found: {}:{}, skipping", 
                                species_name, seq_name);
                    continue;
                }
                
                // 批量添加区间（segment中的坐标是遮蔽后的坐标，需要转换为原始坐标）
                masked_manager->addMaskIntervals(seq_name, intervals);
                species_total_intervals += intervals.size();
                #ifdef _DEBUG_
                spdlog::debug("[addAlignedRegionsAsMask] Added {} intervals for {}:{}", 
                             intervals.size(), species_name, seq_name);
                #endif
            }
            
            // 定案该物种的所有遮蔽区间
            masked_manager->finalizeMaskIntervals();
            
            spdlog::info("[addAlignedRegionsAsMask] Successfully added {} mask intervals for species {}", 
                        species_total_intervals, species_name);
        }
        catch (const std::exception& e) {
            spdlog::error("[addAlignedRegionsAsMask] Error processing species {}: {}", 
                         species_name, e.what());
        }
    }
    
    spdlog::info("[addAlignedRegionsAsMask] Mask interval addition completed for all species");
}

MultipleRareAligner::MultipleRareAligner(
    const FilePath& work_dir_,       // 与声明中的类型、顺序一致
    SpeciesPathMap& species_path_map_,
    NewickParser& newick_tree_,
    uint_t thread_num_,              // 同理
    uint_t chunk_size_,
    uint_t overlap_size_,
    uint_t min_anchor_length_,
    uint_t max_anchor_frequency_
)
    : work_dir(work_dir_),                                  // 初始化成员
    index_dir(work_dir_ / INDEX_DIR),
    species_path_map(species_path_map_),
    newick_tree(newick_tree_),
    chunk_size(chunk_size_),
    overlap_size(overlap_size_),
    min_anchor_length(min_anchor_length_),
    max_anchor_frequency(max_anchor_frequency_),
    thread_num(thread_num_)
{
    // 确保工作目录存在
    if (!std::filesystem::exists(work_dir)) {
        std::filesystem::create_directories(work_dir);
        spdlog::info("Created work directory: {}", work_dir.string());
    }

    // 确保索引目录存在（默认放在 work_dir/index）
    if (!std::filesystem::exists(index_dir)) {
        std::filesystem::create_directories(index_dir);
        spdlog::info("Created index directory: {}", index_dir.string());
    }

    this->group_id = 0;
    this->round_id = 0;

}

void compareMatchedSequences(
    const SpeciesMatchVec3DPtrMapPtr& match_ptr,
    const std::unordered_map<SpeciesName, std::shared_ptr<SeqPro::ManagerVariant>>& seqpro_managers,
    const SpeciesName& ref_name)
{
    // 取片段函数
    auto subSeq = [&](const SeqPro::ManagerVariant& mv,
        const ChrIndex& chr, Coord_t b, Coord_t l) -> std::string {
            return std::visit([&](auto& p) {
                using T = std::decay_t<decltype(p)>;
                if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::SequenceManager>>) {
                    return p->getSubSequence(chr, b, l);
                }
                else if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                    return p->getOriginalManager().getSubSequence(chr, b, l);
                }
                }, mv);
        };

    // 遍历所有匹配
    for (auto& kv : *match_ptr) {
        for (auto& mv2 : *kv.second) {
            for (auto& mv1 : mv2) {
                for (auto& m : mv1) {
                    // 判断参考是否为Masked
                    if (std::holds_alternative<std::unique_ptr<SeqPro::MaskedSequenceManager>>(*seqpro_managers.at(ref_name))) {
                        std::string ref_seq = subSeq(*seqpro_managers.at(ref_name), m.ref_chr_index, m.ref_start, m.match_len());
                        std::string query_seq = subSeq(*seqpro_managers.at(kv.first), m.qry_chr_index, m.qry_start, m.match_len());
                        if (m.strand() == Strand::REVERSE) reverseComplement(query_seq);
                        if (ref_seq != query_seq) {
                            spdlog::error("Ref and query sequences do not match for {}: {} vs {} (ref_start: {}, query_start: {})",
                                m.qry_chr_index, ref_seq, query_seq, m.ref_start, m.qry_start);
                        }
                    }
                    else {
                        std::string ref_seq = subSeq(*seqpro_managers.at(ref_name), m.ref_chr_index, m.ref_start, m.match_len());
                        std::string query_seq = subSeq(*seqpro_managers.at(kv.first), m.qry_chr_index, m.qry_start, m.match_len());
                        if (m.strand() == Strand::REVERSE) reverseComplement(query_seq);
                        if (ref_seq != query_seq) {
                            spdlog::error("Ref and query sequences do not match for {}: {} vs {} (ref_start: {}, query_start: {})",
                                m.qry_chr_index, ref_seq, query_seq, m.ref_start, m.qry_start);
                        }
                    }
                }
            }
        }
    }
}

std::unique_ptr<RaMesh::RaMeshMultiGenomeGraph> MultipleRareAligner::
starAlignment(
    std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers,
    uint_t tree_root,
    SearchMode                 search_mode,
    bool                       fast_build,
    bool                       allow_MEM,
    bool                       mask_mode,
    SeqPro::Length sampling_interval,
    uint_t min_span)
{
    std::vector<int> leaf_vec = newick_tree.orderLeavesGreedyMinSum(tree_root);
    //std::swap(leaf_vec.front(), leaf_vec.back());
	uint_t leaf_num = leaf_vec.size();
    // 初始化Ref缓存
    sdsl::int_vector<0> ref_global_cache;
    // 创建共享线程池，供比对和过滤过程共同使用
    // uint_t count = 0;
    // 创建当前迭代的多基因组图
    auto multi_graph = std::make_unique<RaMesh::RaMeshMultiGenomeGraph>(seqpro_managers);
    //for (uint_t i = 0; i < 1; i++) { 
    for (uint_t i = 0; i < leaf_num; i++) {
        //auto multi_graph = std::make_unique<RaMesh::RaMeshMultiGenomeGraph>(seqpro_managers);
        // 使用工具函数构建缓存
        spdlog::info("build ref global cache for {}", newick_tree.getNodes()[leaf_vec[i]].name);
        SequenceUtils::buildRefGlobalCache(seqpro_managers[newick_tree.getNodes()[leaf_vec[i]].name], sampling_interval, ref_global_cache);
		uint_t ref_id = leaf_vec[i];
		SpeciesName ref_name = newick_tree.getNodes()[ref_id].name;
        std::unordered_map<SpeciesName, SeqPro::SharedManagerVariant> species_fasta_manager_map;
        for (uint_t j = i; j < leaf_num; j++) {
            uint_t query_id = leaf_vec[j];
			SpeciesName query_name = newick_tree.getNodes()[query_id].name;
            auto query_fasta_manager = seqpro_managers.at(query_name);
            species_fasta_manager_map.emplace(query_name, query_fasta_manager);
        }
		bool allow_short_mum = (i == 0) ? true : false; // 允许短MUMs从第二轮开始 
        spdlog::info("align multiple genome for {}", ref_name);
        // TODO 不同模式下最小长度要不同
        SpeciesMatchVec3DPtrMapPtr match_ptr = alignMultipleGenome(
            ref_name, species_fasta_manager_map,
            ACCURATE_SEARCH, fast_build, allow_MEM, allow_short_mum, ref_global_cache, sampling_interval
        );
//#ifdef _DEBUG_
//        compareMatchedSequences(match_ptr, species_fasta_manager_map, ref_name);
//#endif
        spdlog::info("align multiple genome for {} done", ref_name);


        // 使用同一个线程池进行过滤比对结果，获取cluster数据
        spdlog::info("filter multiple species anchors for {}", ref_name);
		SpeciesClusterMapPtr cluster_map = filterMultipeSpeciesAnchors(
			ref_name, species_fasta_manager_map, match_ptr, min_span);
 		spdlog::info("filter multiple species anchors for {} done", ref_name);

        // 并行构建多个比对结果图，共用线程池
        spdlog::info("construct multiple genome graphs for {}", ref_name);


        //constructMultipleGraphsByGreedyByRef(
        //    seqpro_managers, ref_name, *cluster_map, *multi_graph, min_span);

        constructMultipleGraphsByDp(
            seqpro_managers, ref_name, *cluster_map, *multi_graph, min_span, i==0);


        multi_graph->extendRefNodes(ref_name, seqpro_managers, thread_num);


        multi_graph->optimizeGraphStructure();
// #ifdef _DEBUG_
//         multi_graph->verifyGraphCorrectness(ref_name, true);
// #endif // _DEBUG_
        multi_graph->verifyGraphCorrectness(ref_name, true);
        spdlog::info("construct multiple genome graphs for {} done", ref_name);

        // SpeciesName target_species = "simOrang";
        // auto it = multi_graph->species_graphs.find(target_species);
        // if (it != multi_graph->species_graphs.end()) {
        //     const auto& genome_graph = it->second;
        //     spdlog::info("Checking unaligned regions for species {}", target_species);
        //
        //     for (const auto& [chr_name, genome_end] : genome_graph.chr2end) {
        //         spdlog::info("Chromosome {}", chr_name);
        //
        //         // 直接调用前面定义的函数
        //         RaMesh::reportUnalignedRegions(
        //             genome_end,
        //             seqpro_managers.at(target_species),
        //             chr_name
        //         );
        //     }
        // }
        // else {
        //     spdlog::warn("Species {} not found in multi_graph", target_species);
        // }

        //target_species = "simGorilla";
        //it = multi_graph->species_graphs.find(target_species);
        //if (it != multi_graph->species_graphs.end()) {
        //    const auto& genome_graph = it->second;
        //    spdlog::info("Checking unaligned regions for species {}", target_species);

        //    for (const auto& [chr_name, genome_end] : genome_graph.chr2end) {
        //        spdlog::info("Chromosome {}", chr_name);

        //        // 直接调用前面定义的函数
        //        RaMesh::reportUnalignedRegions(
        //            genome_end,
        //            seqpro_managers.at(target_species),
        //            chr_name
        //        );
        //    }
        //}
        //else {
        //    spdlog::warn("Species {} not found in multi_graph", target_species);
        //}

        spdlog::info("merge multiple genome graphs for {}", ref_name);
        multi_graph->mergeMultipleGraphs(ref_name, thread_num);
        spdlog::info("merge multiple genome graphs for {} done", ref_name);
        multi_graph->verifyGraphCorrectness(true);
        multi_graph->optimizeGraphStructure();
        spdlog::info("optimize graph genome graphs for {} done", ref_name);
		multi_graph->markAllExtended();

#ifdef _DEBUG_
        multi_graph->verifyGraphCorrectness(ref_name, true, false, false, true, false);
#endif // _DEBUG_


        // 将当前轮次的比对结果作为遮蔽区间添加到 SeqPro managers 中
        // 这样后续轮次就不会重复比对已经成功比对的区间
        spdlog::info("Adding aligned regions as mask intervals for {}", ref_name);
        try {
            addAlignedRegionsAsMask(*multi_graph, seqpro_managers, ref_name);
            spdlog::info("Successfully added mask intervals for round with reference {}", ref_name);
        }
        catch (const std::exception& e) {
            spdlog::error("Failed to add mask intervals for {}: {}", ref_name, e.what());
        }
        // std::string s = std::to_string(count);
        // multi_graph->exportToMaf("/mnt/d/Result/RaMAx/Alignathon/result/primate-small"+ s + ".maf", seqpro_managers, true, false);
        
    }

    // 在所有迭代完成后，进行最终的图正确性验证
    spdlog::default_logger()->flush();


    return std::move(multi_graph);
    //return std::move(std::make_unique<RaMesh::RaMeshMultiGenomeGraph>(seqpro_managers));

}

SpeciesMatchVec3DPtrMapPtr MultipleRareAligner::alignMultipleGenome(
    SpeciesName                ref_name,
    std::unordered_map<SpeciesName, SeqPro::SharedManagerVariant>& species_fasta_manager_map,
    SearchMode                 search_mode,
    bool                       fast_build,
    bool                       allow_MEM,
    bool                       allow_short_mum,
    sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval)
{
    /* ---------- 0. 合法性检查 ---------- */
    if (!species_fasta_manager_map.count(ref_name))
        throw std::runtime_error("[alignMultipleQuerys] reference species not found: " + ref_name);

    if (species_fasta_manager_map.size() <= 1) {
        spdlog::warn("[alignMultipleQuerys] only reference genome present, nothing to align.");
        return std::make_shared<SpeciesMatchVec3DPtrMap>();
    }

    /* ---------- 1. 结果目录与缓存文件 ---------- */
    FilePath result_dir = work_dir / RESULT_DIR
        / ("group_" + std::to_string(group_id))
        / ("round_" + std::to_string(round_id));
    std::filesystem::create_directories(result_dir);
    round_id++;
    FilePath anchor_file = result_dir / (ref_name + "_"
        + SearchModeToString(search_mode) + "." + ANCHOR_EXTENSION);

    /* ---------- 2. 如果已存在结果文件直接读取 ---------- */
    //if (std::filesystem::exists(anchor_file)) {
    //    spdlog::info("[alignMultipleQuerys] Load from {}", anchor_file.string());
    //    auto mp = std::make_shared<SpeciesMatchVec3DPtrMap>();
    //    if (loadSpeciesMatchMap(anchor_file, mp))
    //        return mp;
    //    // 如果读取失败则继续重新计算
    //}

    /* ---------- 3. 准备参考基因组索引 ---------- */
    FilePath ref_index_path = index_dir / ref_name;
    std::filesystem::create_directories(ref_index_path);

    PairRareAligner pra(*this);
	pra.buildIndex(ref_name, *species_fasta_manager_map[ref_name], fast_build);
	spdlog::info("[alignMultipleQuerys] reference index built for {}.", ref_name);

    /* ---------- 4. 创建共享线程池 ---------- */
    ThreadPool shared_pool(thread_num);

    /* ---------- 5. 为每个 query 物种启动异步任务 ---------- */
    std::unordered_map<SpeciesName, std::future<MatchVec3DPtr>> fut_map;

    for (auto& kv : species_fasta_manager_map) {
        SpeciesName  sp = kv.first;
        if (sp == ref_name) continue;           // 跳过参考

        std::string   prefix = ref_name + "_vs_" + sp;
        auto& fm = kv.second;

        fut_map.emplace(
            sp,
            std::async(std::launch::async,
                [&pra, prefix, &fm, search_mode, allow_MEM,allow_short_mum,  &shared_pool, &ref_global_cache, sampling_interval]() -> MatchVec3DPtr {
                    return pra.findQueryFileAnchor(prefix, *fm, search_mode, allow_MEM, allow_short_mum, shared_pool, ref_global_cache, sampling_interval, true);
                })
        );
    }

    shared_pool.waitAllTasksDone();

    /* ---------- 6. 收集所有结果 ---------- */
    auto result_map = std::make_shared<SpeciesMatchVec3DPtrMap>();

    size_t total = fut_map.size();
    size_t count = 0;
    size_t next_progress = 1;  // 下一次打印进度的“阶段”（1 到 20）

    for (auto& kv : fut_map) {
        const SpeciesName& sp = kv.first;
        try {
            MatchVec3DPtr mv3 = kv.second.get();  // 等待并获取
            (*result_map)[sp] = std::move(mv3);
            spdlog::info("[alignMultipleQuerys] {} aligned.", sp);
        }
        catch (const std::exception& e) {
            spdlog::error("[alignMultipleQuerys] {} failed: {}", sp, e.what());
        }

        ++count;

        // 计算是否达到下一个进度段（总共 20 段）
        size_t progress_stage = (count * 20) / total;
        if (progress_stage >= next_progress) {
            int percent = static_cast<int>((progress_stage * 100) / 20);
            spdlog::info("[alignMultipleQuerys] Progress: {}%", percent);
            next_progress = progress_stage + 1;
        }
    }

    /* ---------- 7. 保存到文件 ---------- */
    // saveSpeciesMatchMap(anchor_file, result_map);
    spdlog::info("[alignMultipleQuerys] all species done. Saved to {}", anchor_file.string());

    //{
    //    SpeciesName target_species = "simOrang";
    //    //int qry_chr_index = 1;
    //    //int_t region_start = 54620705;
    //    //int_t region_end = 54636058;
  
    //    int qry_chr_index = 1;
    //    int_t region_start = 54320705;
    //    int_t region_end = 54336058;

    //    std::vector<Match> matched_vec; // 先收集

    //    auto it = result_map->find(target_species);
    //    if (it != result_map->end()) {
    //        auto match3d_ptr = it->second;
    //        if (match3d_ptr) {
    //            for (const auto& by_ref : *match3d_ptr) {
    //                for (const auto& by_qry : by_ref) {
    //                    for (const auto& m : by_qry) {
    //                        if (m.qry_chr_index == qry_chr_index) {
    //                            Coord_t m_start = m.ref_start;
    //                            Coord_t m_end = m.ref_start + m.match_len();
    //                            if (m_end > region_start && m_start < region_end) {
    //                                matched_vec.push_back(m);
    //                            }
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //    }
    //    else {
    //        spdlog::warn("Species {} not found in result_map", target_species);
    //    }

    //    if (!matched_vec.empty()) {
    //        std::sort(matched_vec.begin(), matched_vec.end(),
    //            [](const Match& a, const Match& b) {
    //                return a.ref_start < b.ref_start;
    //            });

    //        for (const auto& m : matched_vec) {
    //            spdlog::info(
    //                "[{}] chr={} 区间 [{}, {}) match: "
    //                "ref_chr={} ref_start={} len={} strand={}",
    //                target_species,
    //                m.qry_chr_index,
    //                m.qry_start,
    //                m.qry_start + m.match_len(),
    //                m.ref_chr_index,
    //                m.ref_start,
    //                m.match_len(),
    //                (m.strand() == Strand::FORWARD ? "FORWARD" : "REVERSE")
    //            );
    //        }
    //    }
    //    else {
    //        spdlog::info("[{}] chr={} 区间 [{}, {}) 没有任何 match",
    //            target_species, qry_chr_index, region_start, region_end);
    //    }
    //}

    return result_map;
}


SpeciesClusterMapPtr MultipleRareAligner::filterMultipeSpeciesAnchors(
    SpeciesName                       ref_name,
    std::unordered_map<SpeciesName, SeqPro::SharedManagerVariant>& species_fm_map,
    SpeciesMatchVec3DPtrMapPtr        species_match_map,
    uint_t min_span)
{
    if (!species_match_map || species_match_map->empty()) {
        return std::make_shared<SpeciesClusterMap>();
    }
    ThreadPool shared_pool(thread_num);
    /*-------------- 预分配表 ----------------------------------------*/
    std::unordered_map<SpeciesName, MatchByStrandByQueryRefPtr> unique_map;
    std::unordered_map<SpeciesName, MatchByStrandByQueryRefPtr> repeat_map;
    SpeciesClusterMap cluster_map;                    // 最终结果
    spdlog::info("Group Match By Query and Ref...");

    /*========================= Phase-1  : group =====================*/
    for (auto it = species_match_map->begin();
        it != species_match_map->end(); ++it)
    {
        const SpeciesName  species = it->first;
        if (species == ref_name) continue;

        MatchVec3DPtr mv3_ptr = it->second;
        auto& qfm = species_fm_map.at(species);
        auto& rfm = species_fm_map.at(ref_name);

        auto u_ptr = std::make_shared<MatchByStrandByQueryRef>();
        auto r_ptr = std::make_shared<MatchByStrandByQueryRef>();
        unique_map[species] = u_ptr;
        repeat_map[species] = r_ptr;

        // 把所有需要的变量显式捕获
        shared_pool.enqueue(
            [mv3_ptr,
            u_ptr,
            r_ptr,
            qfm,
            rfm,
            &shared_pool]() mutable
            {
                groupMatchByQueryRef(mv3_ptr,
                    u_ptr,
                    r_ptr,
                    *rfm,
                    *qfm,
                    shared_pool);          // 仍用同一池
            });
    }
    shared_pool.waitAllTasksDone();                          // —— Phase-1 完

    // 一步到位，最快释放
    for (auto it = species_match_map->begin();
        it != species_match_map->end(); ++it)
    {
        it->second.reset();   // 释放 MatchVec3D 对象占用的全部内存
    }
    species_match_map->clear();         // 清掉 map 自身节点
    spdlog::info("Group Match By Query and Ref Done");


 //   /*========================= Phase-2  : sort ======================*/
	//spdlog::info("Sort Match By Query Start...");   
 //   for (auto it = unique_map.begin(); it != unique_map.end(); ++it) {
 //       MatchByStrandByQueryRefPtr u_ptr = it->second;
 //       shared_pool.enqueue([u_ptr, &shared_pool]() mutable {
 //           sortMatchByQueryStart(u_ptr, shared_pool);
 //           });
 //   }
 //   for (auto it = repeat_map.begin(); it != repeat_map.end(); ++it) {
 //       MatchByStrandByQueryRefPtr r_ptr = it->second;
 //       shared_pool.enqueue([r_ptr, &shared_pool]() mutable {
 //           sortMatchByQueryStart(r_ptr, shared_pool);
 //           });
 //   }
 //   shared_pool.waitAllTasksDone();                          // —— Phase-2 完
	//spdlog::info("Sort Match By Query Done");

    /*========================= Phase-3  :  ===================*/
	spdlog::info("Cluster All Chr Match Start...");
    using Fut = std::future<ClusterVecPtrByStrandByQueryRefPtr>;
    std::unordered_map<SpeciesName, Fut> fut_map;

    for (auto it = unique_map.begin(); it != unique_map.end(); ++it) {
        const SpeciesName  species = it->first;
        auto               u_ptr = it->second;
        auto               r_ptr = repeat_map.at(species);

        fut_map.emplace(
            species,
            shared_pool.enqueue(
                [u_ptr, r_ptr, &shared_pool, min_span]() mutable -> ClusterVecPtrByStrandByQueryRefPtr {
                    return clusterAllChrMatch(u_ptr, r_ptr, shared_pool, min_span);
                })
        );
    }

    // shared_pool.waitAllTasksDone();
    // 收集 future
    size_t total = fut_map.size();
    size_t count = 0;
    size_t next_progress = 1;  // 下一个打印阶段（1~20）

    for (auto it = fut_map.begin(); it != fut_map.end(); ++it) {
        const SpeciesName& species = it->first;
        ClusterVecPtrByStrandByQueryRefPtr clus_ptr = it->second.get();
        cluster_map.emplace(species, std::move(clus_ptr));

        ++count;
        size_t progress_stage = (count * 20) / total;
        if (progress_stage >= next_progress || count == total) {
            int percent = static_cast<int>((progress_stage * 100) / 20);
            spdlog::info("[cluster] Progress: {}% ({} of {})", percent, count, total);
            next_progress = progress_stage + 1;
        }
    }


    shared_pool.waitAllTasksDone();                          // —— Phase-3 完

    // 返回cluster数据用于后续处理
    auto cluster_map_ptr = std::make_shared<SpeciesClusterMap>(std::move(cluster_map));
	spdlog::info("Cluster All Chr Match Done, species num: {}", cluster_map_ptr->size());
    return cluster_map_ptr;

}

/* ============================================================= *
 *  多基因组比对类的成员函数：并行构建多个比对结果图，共用一个线程池
 * ============================================================= */
/**
 * @brief  并行构建多个比对结果图，基于贪婪算法处理多个物种的cluster数据
 *
 * @param ref_name             参考物种名称
 * @param species_cluster_map  所有物种的cluster数据映射
 * @param graph                多基因组图对象
 * @param shared_pool          共享线程池
 * @param min_span             最小跨度阈值
 */
//void MultipleRareAligner::constructMultipleGraphsByGreedy(
//std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers,
//    SpeciesName ref_name,
//    const SpeciesClusterMap& species_cluster_map,
//    RaMesh::RaMeshMultiGenomeGraph& graph,
//    uint_t min_span)
//{
//    if (species_cluster_map.empty()) {
//        spdlog::warn("[constructMultipleGraphsByGreedy] Empty cluster map, nothing to process.");
//        return;
//    }
//
//    spdlog::info("[constructMultipleGraphsByGreedy] Processing {} species clusters",
//                species_cluster_map.size());
//
//    ThreadPool shared_pool(thread_num);
//
//    // 【修复】：添加互斥锁保护graph的并发访问
//    std::mutex graph_mutex;
//
//    /* ---------- 1. 为每个物种并行处理cluster数据 ---------- */
//    std::vector<std::future<void>> species_futures;
//    species_futures.reserve(species_cluster_map.size());
//
//    for (const auto& [species_name, cluster_ptr] : species_cluster_map) {
//        if (species_name == ref_name) continue;  // 跳过参考物种
//
//        // 为每个物种启动异步任务
//        auto species_future = std::async(std::launch::async,
//            [this, ref_name, species_name, cluster_ptr, &graph, &shared_pool, &graph_mutex, min_span,&seqpro_managers]() {
//                try {
//                    if (!cluster_ptr || cluster_ptr->empty()) {
//                        spdlog::warn("[constructMultipleGraphsByGreedy] Empty cluster data for species: {}",
//                                   species_name);
//                        return;
//                    }
//
//                    // 使用PairRareAligner的贪婪算法构建图
//                    PairRareAligner pra(*this);
//                    pra.ref_name = ref_name;
//                    // 【修复】：设置ref_seqpro_manager，避免空指针
//                    pra.ref_seqpro_manager = &(*seqpro_managers.at(ref_name));
//
//                    //// 【修复】：避免过度并行化，改为串行处理chromosome数据
//                    //// 每个物种内部串行处理，避免graph的深度并发访问
//                    //for (const auto& strand_data : *cluster_ptr) {
//                    //    for (const auto& query_ref_data : strand_data) {
//                    //        for (const auto& cluster_vec : query_ref_data) {
//                    //            if (cluster_vec && !cluster_vec->empty()) {
//                    //                // 将集合转换为向量以便处理
//                    //                auto cluster_vec_ptr = std::make_shared<MatchClusterVec>(
//                    //                    cluster_vec->begin(), cluster_vec->end());
//
//                    //                // 【修复】：使用互斥锁保护graph访问
//                    //                {
//                    //                    std::lock_guard<std::mutex> lock(graph_mutex);
//                    //                    pra.constructGraphByGreedy(species_name, *seqpro_managers[species_name],cluster_vec_ptr,
//                    //                                             graph, min_span);
//                    //                }
//                    //            }
//                    //        }
//                    //    }
//                    //}
//                    {
//                        std::lock_guard<std::mutex> lock(graph_mutex);
//                        pra.constructGraphByGreedy(species_name, *seqpro_managers[species_name], cluster_ptr,
//                                                                      graph, min_span);
//                    }
//
//                    spdlog::info("[constructMultipleGraphsByGreedy] Species {} processed successfully",
//                               species_name);
//                }
//                catch (const std::exception& e) {
//                    spdlog::error("[constructMultipleGraphsByGreedy] Error processing species {}: {}",
//                                species_name, e.what());
//                }
//            });
//
//        species_futures.emplace_back(std::move(species_future));
//    }
//
//    /* ---------- 2. 等待所有物种处理完成 ---------- */
//    for (auto& future : species_futures) {
//        try {
//            future.wait();
//        }
//        catch (const std::exception& e) {
//            spdlog::error("[constructMultipleGraphsByGreedy] Error waiting for species processing: {}",
//                        e.what());
//        }
//    }
//
//    // 确保共享线程池中的所有任务都完成
//    shared_pool.waitAllTasksDone();
//
//    spdlog::info("[constructMultipleGraphsByGreedy] All species graphs constructed successfully");
//}


/* ============================================================= *
 *  多基因组比对类的成员函数：并行构建多个比对结果图，共用一个线程池
 * ============================================================= */
 /**
  * @brief  并行构建多个比对结果图，基于贪婪算法处理多个物种的cluster数据
  *
  * @param ref_name             参考物种名称
  * @param species_cluster_map  所有物种的cluster数据映射
  * @param graph                多基因组图对象
  * @param shared_pool          共享线程池
  * @param min_span             最小跨度阈值
  */
void MultipleRareAligner::constructMultipleGraphsByDpByRef(
    std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers,
    SpeciesName ref_name,
    const SpeciesClusterMap& species_cluster_map,
    RaMesh::RaMeshMultiGenomeGraph& graph,
    uint_t min_span)
{
    if (species_cluster_map.empty()) {
        spdlog::warn("[constructMultipleGraphsByGreedy] Empty cluster map, nothing to process.");
        return;
    }

    spdlog::info("[constructMultipleGraphsByGreedy] Processing {} species clusters",
        species_cluster_map.size());

    ThreadPool pool(thread_num);
    std::map<SpeciesName, ClusterVecPtrByRefPtr> result_map;
    std::vector<std::future<void>> futures;

    for (const auto& [species_name, cluster_ptr_3d] : species_cluster_map) {
        // 为空跳过
        if (!cluster_ptr_3d) continue;

        futures.emplace_back(
            pool.enqueue([&, species_name, cluster_ptr_3d]() {
                try {
                    // 生成该物种的按 Ref 分组聚簇
                    ClusterVecPtrByRefPtr grouped_ref_clusters =
                        groupClustersToRefVec(cluster_ptr_3d, pool, thread_num);

                    result_map[species_name] = std::move(grouped_ref_clusters);

                    spdlog::info("[constructMultipleGraphsByGreedy] Finished clustering for species: {}", species_name);
                }
                catch (const std::exception& e) {
                    spdlog::error("[constructMultipleGraphsByGreedy] Error processing species {}: {}", species_name, e.what());
                }
                })
        );
    }

    // 等待所有任务完成
    for (auto& fut : futures) fut.get();
    pool.waitAllTasksDone();

    PairRareAligner pra(*this);
    pra.ref_name = ref_name;
    // 【修复】：设置ref_seqpro_manager，避免空指针
    pra.ref_seqpro_manager = &(*seqpro_managers.at(ref_name));

    for (auto& [species_name, cluster_ref_ptr] : result_map) {
        for (auto& cluster_ptr : *cluster_ref_ptr) {

            //pool.enqueue([&, species_name, cluster_ptr]() {
            //    pra.constructGraphByGreedyByRef(species_name, *seqpro_managers[species_name], cluster_ptr,
            //        graph, pool, min_span);
            //    });
            for (auto& cluster : *cluster_ptr) {
                pra.constructGraphByGreedyByRef(species_name, *seqpro_managers[species_name], cluster_ptr,
                    graph, min_span, false);
            }
        }
    }
    pool.waitAllTasksDone();

        for (auto& [species_name, genome_graph] : graph.species_graphs) {
            // if (species_name == ref_name) continue;
            for (auto& [chr_name, end] : genome_graph.chr2end) {
                //pool.enqueue([&]() {
                //    end.removeOverlap();
                //    });
                end.removeOverlap(species_name == ref_name);
                
            }

        }
        pool.waitAllTasksDone();

        spdlog::info("[constructMultipleGraphsByGreedy] All species graphs constructed successfully");
    }

void MultipleRareAligner::constructMultipleGraphsByGreedyByRef(
    std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers,
    SpeciesName ref_name,
    const SpeciesClusterMap& species_cluster_map,
    RaMesh::RaMeshMultiGenomeGraph& graph,
    uint_t min_span)
{
    if (species_cluster_map.empty()) {
        spdlog::warn("[constructMultipleGraphsByGreedy] Empty cluster map, nothing to process.");
        return;
    }

    spdlog::info("[constructMultipleGraphsByGreedy] Processing {} species clusters",
        species_cluster_map.size());


    ThreadPool pool(thread_num);
    std::map<SpeciesName, ClusterVecPtrByRefPtr> result_map;
    std::vector<std::future<void>> futures;

    for (const auto& [species_name, cluster_ptr_3d] : species_cluster_map) {
        // 为空跳过
        if (!cluster_ptr_3d) continue;

        futures.emplace_back(
            pool.enqueue([&, species_name, cluster_ptr_3d]() {
                try {
                    // 生成该物种的按 Ref 分组聚簇
                    ClusterVecPtrByRefPtr grouped_ref_clusters =
                        groupClustersToRefVec(cluster_ptr_3d, pool, thread_num);

                    result_map[species_name] = std::move(grouped_ref_clusters);

                    spdlog::info("[constructMultipleGraphsByGreedy] Finished clustering for species: {}", species_name);
                }
                catch (const std::exception& e) {
                    spdlog::error("[constructMultipleGraphsByGreedy] Error processing species {}: {}", species_name, e.what());
                }
                })
        );
    }

    // 等待所有任务完成
    for (auto& fut : futures) fut.get();
    pool.waitAllTasksDone();

    PairRareAligner pra(*this);
    pra.ref_name = ref_name;
    // 【修复】：设置ref_seqpro_manager，避免空指针
    pra.ref_seqpro_manager = &(*seqpro_managers.at(ref_name));

    for (auto& [species_name, cluster_ref_ptr] : result_map) {
        for (auto& cluster_ptr : *cluster_ref_ptr) {

            //pool.enqueue([&, species_name, cluster_ptr]() {
            //    pra.constructGraphByGreedyByRef(species_name, *seqpro_managers[species_name], cluster_ptr,
            //        graph, pool, min_span);
            //    });
            
            pra.constructGraphByGreedyByRef(species_name, *seqpro_managers[species_name], cluster_ptr,
                    graph, min_span, false);

                //pra.constructGraphByDpByRef(species_name, *seqpro_managers[species_name], cluster_ptr,
                //    graph, pool, thread_num, min_span, false);
            
        }
    }
    pool.waitAllTasksDone();

    for (auto& [species_name, genome_graph] : graph.species_graphs) {
        // if (species_name == ref_name) continue;
        for (auto& [chr_name, end] : genome_graph.chr2end) {
            //pool.enqueue([&]() {
            //    end.removeOverlap();
            //    });
            end.removeOverlap(species_name == ref_name);

        }

    }
    pool.waitAllTasksDone();

    spdlog::info("[constructMultipleGraphsByGreedy] All species graphs constructed successfully");
}


void MultipleRareAligner::constructMultipleGraphsByDp(
    std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers,
    SpeciesName ref_name,
    const SpeciesClusterMap& species_cluster_map,
    RaMesh::RaMeshMultiGenomeGraph& graph,
    uint_t min_span, bool is_first) {
	if (species_cluster_map.empty()) {
		spdlog::warn("[constructMultipleGraphsByGreedy] Empty cluster map, nothing to process.");
		return;
	}

    spdlog::info("[constructMultipleGraphsByGreedy] Processing {} species clusters",
        species_cluster_map.size());

    PairRareAligner pra(*this);
    pra.ref_name = ref_name;
    // 【修复】：设置ref_seqpro_manager，避免空指针
    pra.ref_seqpro_manager = &(*seqpro_managers.at(ref_name));

    std::map<SpeciesName, AnchorPtrVecByStrandByQueryByRefPtr> anchor_map;

    ThreadPool shared_pool(thread_num);
    
    for (const auto& [species_name, cluster_ptr_3d] : species_cluster_map) {
        // 为空跳过
        if (!cluster_ptr_3d) continue;
        anchor_map[species_name] = pra.extendClusterToAnchorByChr(species_name, *seqpro_managers[species_name], cluster_ptr_3d, shared_pool, is_first);
    }
    shared_pool.waitAllTasksDone();
    std::cout << "extend successfully with " << is_first << std::endl;
    
    for (auto& [species_name, anchor_ptr] : anchor_map) {
        pra.filterAnchorByDP(anchor_ptr);
		std::cout << "filter successfully for " << species_name << " with " << anchor_ptr->size() << std::endl;
		pra.constructGraphByDP(species_name, *seqpro_managers[species_name], anchor_ptr, graph);
    }
}
