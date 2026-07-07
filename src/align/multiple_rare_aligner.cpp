#include <sequence_utils.h>
#include <algorithm>
#include <stdexcept>
#include <unordered_set>

#include "rare_aligner.h"
#include "anchor.h"  // 包含 UnionFind 定义
#include "SeqPro.h"  // 包含 SeqPro 相关定义
#include "ramesh.h"  // 包含 RaMesh 图结构定义

// 辅助函数：根据CIGAR字符串计算query区间对应关系
namespace {

    constexpr size_t kMaxReferenceSequenceCount = 100000;

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

void exportMaskIntervalsToDirectory(
    const std::filesystem::path& export_dir,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {

    if (export_dir.empty()) {
        spdlog::warn("[mask-export] Export directory is empty, skip exporting.");
        return;
    }

    std::error_code ec;
    std::filesystem::create_directories(export_dir, ec);
    if (ec) {
        spdlog::error("[mask-export] Failed to create export directory {}: {}", export_dir.string(), ec.message());
        return;
    }

    for (const auto& [species_name, manager_variant] : seqpro_managers) {
        if (!manager_variant) {
            spdlog::warn("[mask-export] Skip species {} due to null manager.", species_name);
            continue;
        }

        auto* masked_manager = std::visit(
            [](auto& ptr) -> SeqPro::MaskedSequenceManager* {
                using T = std::decay_t<decltype(ptr)>;
                if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                    return ptr.get();
                }
                return nullptr;
            },
            *manager_variant);

        if (!masked_manager) {
            spdlog::warn("[mask-export] Species {} does not use MaskedSequenceManager, skip exporting.", species_name);
            continue;
        }

        std::ostringstream buffer;
        size_t species_total_intervals = 0;
        auto seq_names = masked_manager->getSequenceNames();

        for (const auto& seq_name : seq_names) {
            const auto& intervals = masked_manager->getMaskIntervals(seq_name);
            if (intervals.empty()) {
                continue;
            }

            buffer << '>' << seq_name << '\n';
            for (size_t i = 0; i < intervals.size(); ++i) {
                const auto& interval = intervals[i];
                buffer << (interval.start + 1) << '-' << interval.end;
                if (i + 1 < intervals.size()) {
                    buffer << ' ';
                }
            }
            buffer << '\n';
            species_total_intervals += intervals.size();
        }

        if (species_total_intervals == 0) {
            spdlog::info("[mask-export] No intervals found for species {}, skip writing file.", species_name);
            continue;
        }

        auto output_path = export_dir / (species_name + ".intervals");
        std::ofstream ofs(output_path);
        if (!ofs.is_open()) {
            spdlog::error("[mask-export] Failed to open file {} for species {}", output_path.string(), species_name);
            continue;
        }

        ofs << buffer.str();
        ofs.close();
        spdlog::info("[mask-export] Exported {} intervals for species {} to {}", species_total_intervals, species_name, output_path.string());
    }
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

    struct SpeciesAssemblyQuality {
        SpeciesName name;
        SeqPro::Length n50 = 0;
        SeqPro::Length total_length = 0;
        size_t sequence_count = 0;
        bool reference_eligible = false;
    };

    bool isReferenceEligibleSequenceCount(size_t sequence_count) {
        return sequence_count > 0 && sequence_count <= kMaxReferenceSequenceCount;
    }

    size_t getManagerSequenceCount(const SeqPro::SharedManagerVariant& manager_variant) {
        if (!manager_variant) {
            return 0;
        }

        return std::visit(
            [](auto const& manager_ptr) -> size_t {
                return manager_ptr ? manager_ptr->getSequenceCount() : 0;
            },
            *manager_variant
        );
    }

    void validateReferenceSequenceCount(
        const SpeciesName& ref_name,
        const SeqPro::SharedManagerVariant& manager_variant) {

        const size_t sequence_count = getManagerSequenceCount(manager_variant);
        if (isReferenceEligibleSequenceCount(sequence_count)) {
            return;
        }

        throw std::runtime_error(
            "[reference-selection] Genome " + ref_name + " cannot be used as reference: sequence_count=" +
            std::to_string(sequence_count) + ", max_allowed=" +
            std::to_string(kMaxReferenceSequenceCount));
    }

    SeqPro::Length calculateN50(std::vector<SeqPro::Length> lengths, SeqPro::Length total_length) {
        if (lengths.empty() || total_length == 0) {
            return 0;
        }

        std::sort(lengths.begin(), lengths.end(), [](SeqPro::Length a, SeqPro::Length b) {
            return a > b;
        });

        const SeqPro::Length half_total = (total_length + 1) / 2;
        SeqPro::Length cumulative = 0;

        for (SeqPro::Length length : lengths) {
            cumulative += length;
            if (cumulative >= half_total) {
                return length;
            }
        }

        return 0;
    }

    SpeciesAssemblyQuality getAssemblyQuality(
        const SpeciesName& species_name,
        const SeqPro::SequenceManager& manager) {

        SpeciesAssemblyQuality quality;
        quality.name = species_name;

        auto seq_names = manager.getSequenceNames();
        quality.sequence_count = seq_names.size();
        quality.reference_eligible = isReferenceEligibleSequenceCount(quality.sequence_count);

        std::vector<SeqPro::Length> lengths;
        lengths.reserve(seq_names.size());

        for (const auto& seq_name : seq_names) {
            SeqPro::Length length = manager.getSequenceLength(seq_name);
            if (length == 0) {
                continue;
            }
            lengths.push_back(length);
            quality.total_length += length;
        }

        quality.n50 = calculateN50(std::move(lengths), quality.total_length);
        return quality;
    }

    SpeciesAssemblyQuality getAssemblyQuality(
        const SpeciesName& species_name,
        const SeqPro::MaskedSequenceManager& manager) {
        return getAssemblyQuality(species_name, manager.getOriginalManager());
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
// ------------------------------------------------------------
// MultipleRareAligner::starAlignment
// 说明：
// - 输入 seqpro_managers：物种名 -> SeqPro::SharedManagerVariant（shared_ptr<variant<unique_ptr<...>>>）
// - ref_name：指定参考物种名（会被放到处理顺序最前）
// - only_one_round：若为 true 只跑一轮；否则按物种数跑多轮（每轮换一个 ref）
// - fast_build：索引构建是否走快速模式
// - sampling_interval：用于 ref_global_cache 的采样间隔（加速 global->local 映射）
// - min_span：用于聚簇/过滤时的最小跨度阈值
// 返回：最终构建完成的多基因组图（unique_ptr）
// ------------------------------------------------------------
std::unique_ptr<RaMesh::RaMeshMultiGenomeGraph> MultipleRareAligner::
starAlignment(
    std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers,
    std::string ref_name,
    bool only_one_round,
    bool                       fast_build,
    SeqPro::Length sampling_interval,
    uint_t min_span)
{
    // ------------------------------------------------------------
    // （注释代码）可选：从磁盘导入 mask intervals 到 MaskedSequenceManager
    // 目的：如果之前导出过遮蔽区间，则在后续轮次继续使用遮蔽信息，避免重复比对
    // 目前整段被注释，不参与逻辑执行
    // ------------------------------------------------------------
    // auto import_mask_if_needed = [&](SeqPro::MaskedSequenceManager& manager,
    //                                              const SpeciesName& species_name) {
    //
    //     std::filesystem::path mask_file = work_dir / "mask_interval" / std::to_string(0) / (species_name + ".intervals");
    //     if (!std::filesystem::exists(mask_file)) {
    //         spdlog::debug("[mask-import] Mask file not found for {} at {}", species_name, mask_file.string());
    //         return;
    //     }
    //
    //     try {
    //         if (manager.loadMaskIntervalsFromFile(mask_file, true)) {
    //             manager.finalizeMaskIntervals();
    //             spdlog::info("[mask-import] Loaded mask intervals for {} from {}", species_name, mask_file.string());
    //         } else {
    //             spdlog::warn("[mask-import] Failed to load mask intervals for {} from {}", species_name, mask_file.string());
    //         }
    //     } catch (const std::exception& e) {
    //         spdlog::error("[mask-import] Exception while loading mask for {}: {}", species_name, e.what());
    //     }
    // };
    // for (auto [sp, mgr_variant] : seqpro_managers) {
    //     if (!mgr_variant) continue;
    //
    //     auto* masked_mgr = std::visit(
    //         [](auto& ptr) -> SeqPro::MaskedSequenceManager* {
    //             using T = std::decay_t<decltype(ptr)>;
    //             if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
    //                 return ptr.get();
    //             }
    //             return nullptr;
    //         },
    //         *mgr_variant);
    //
    //     if (masked_mgr) {
    //         import_mask_if_needed(*masked_mgr, sp);
    //     }
    // }

    // ------------------------------------------------------------
    // 1) 统计每个物种的 assembly N50，用于决定处理顺序
    //    - seqpro_managers 里 value 是 SharedManagerVariant（shared_ptr<variant<...>>）
    //    - N50 根据原始序列长度计算；MaskedSequenceManager 使用原始 manager，避免 mask 影响组装质量
    // ------------------------------------------------------------
    std::vector<SpeciesAssemblyQuality> species_qualities;
    species_qualities.reserve(seqpro_managers.size());

    for (const auto& entry : seqpro_managers) {
        const SpeciesName& species_name = entry.first;
        const SeqPro::SharedManagerVariant& shared_mgr_variant = entry.second;

        if (!shared_mgr_variant) {
            // 空指针则 N50 和长度记为 0（保持原逻辑）
            species_qualities.push_back({species_name, 0, 0, 0});
            continue;
        }

        // 用 std::visit 计算该物种 assembly N50
        SpeciesAssemblyQuality quality = std::visit(
            [&species_name](auto const& mgrPtr) -> SpeciesAssemblyQuality {
                using ManagerPtrT = std::decay_t<decltype(mgrPtr)>;

                if constexpr (std::is_same_v<ManagerPtrT, std::unique_ptr<SeqPro::SequenceManager>>) {
                    return mgrPtr ? getAssemblyQuality(species_name, *mgrPtr)
                                  : SpeciesAssemblyQuality{species_name, 0, 0, 0};
                }
                else if constexpr (std::is_same_v<ManagerPtrT, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                    return mgrPtr ? getAssemblyQuality(species_name, *mgrPtr)
                                  : SpeciesAssemblyQuality{species_name, 0, 0, 0};
                }
                else {
                    // 理论上不会到这里（保持原逻辑）
                    return SpeciesAssemblyQuality{species_name, 0, 0, 0};
                }
            },
            *shared_mgr_variant // shared_ptr<variant<...>> 解引用得到 variant
        );

        species_qualities.push_back(std::move(quality));
    }

    // ------------------------------------------------------------
    // 2) 按 assembly N50 从大到小排序（降序）
    //    N50 相同时，使用总长度和物种名做确定性 tie-break
    // ------------------------------------------------------------
    auto quality_order_less = [](const auto& a, const auto& b) {
        if (a.n50 != b.n50) {
            return a.n50 > b.n50;
        }
        if (a.total_length != b.total_length) {
            return a.total_length > b.total_length;
        }
        return a.name < b.name;
    };
    std::sort(species_qualities.begin(), species_qualities.end(), quality_order_less);

    // ------------------------------------------------------------
    // 3) 将指定的 ref_name 移动到最前面（如果存在于列表）
    //    若显式指定的 reference 超过序列数阈值，则直接报错
    // ------------------------------------------------------------
    auto ref_it = std::find_if(
        species_qualities.begin(),
        species_qualities.end(),
        [&ref_name](const SpeciesAssemblyQuality& quality) {
            return quality.name == ref_name;
        }
    );
    if (ref_it != species_qualities.end()) {
        if (!ref_it->reference_eligible) {
            throw std::runtime_error(
                "[reference-selection] Explicit reference " + ref_it->name +
                " cannot be used as reference: sequence_count=" +
                std::to_string(ref_it->sequence_count) + ", max_allowed=" +
                std::to_string(kMaxReferenceSequenceCount));
        }

        SpeciesAssemblyQuality ref = std::move(*ref_it);
        species_qualities.erase(ref_it);
        species_qualities.insert(species_qualities.begin(), std::move(ref));
    }

    // ------------------------------------------------------------
    // 4) 提取排序后的物种名列表 species_order，并单独构建合法 reference 列表
    //    超过序列数阈值的物种保留为 query，但不会进入 reference_order
    // ------------------------------------------------------------
    std::vector<SpeciesName> species_order;
    species_order.reserve(species_qualities.size());
    std::vector<SpeciesName> reference_order;
    reference_order.reserve(species_qualities.size());

    for (const auto& quality : species_qualities) {
        species_order.push_back(quality.name);
        if (quality.reference_eligible) {
            reference_order.push_back(quality.name);
        }
    }

    uint_t leaf_num = species_order.size();
    uint_t reference_num = reference_order.size();
    if (reference_num == 0) {
        throw std::runtime_error(
            "[reference-selection] No eligible reference genome found. Genomes with sequence_count > " +
            std::to_string(kMaxReferenceSequenceCount) + " cannot be used as reference.");
    }

    // 打印物种处理顺序
    spdlog::info("Species processing order ({} total, sorted by assembly N50):", leaf_num);
    for (size_t i = 0; i < species_qualities.size(); ++i) {
        const auto& quality = species_qualities[i];
        spdlog::info("  [{}] {} (N50={}, total_length={}, sequences={}, reference_eligible={})",
            i, quality.name, quality.n50, quality.total_length, quality.sequence_count,
            quality.reference_eligible ? "true" : "false");
        if (!quality.reference_eligible) {
            spdlog::warn(
                "[reference-selection] Excluding {} as reference: sequence_count={}, max_allowed={}",
                quality.name, quality.sequence_count, kMaxReferenceSequenceCount);
        }
    }

    // ------------------------------------------------------------
    // 5) 初始化参考缓存 ref_global_cache
    //    - 用于加速 global 坐标定位到序列（避免频繁二分/查询）
    // ------------------------------------------------------------
    sdsl::int_vector<0> ref_global_cache;

    // ------------------------------------------------------------
    // 6) 创建多基因组图（整个 starAlignment 过程使用同一个 multi_graph）
    // ------------------------------------------------------------
    auto multi_graph = std::make_unique<RaMesh::RaMeshMultiGenomeGraph>(seqpro_managers);

    // ------------------------------------------------------------
    // 7) 决定轮数：only_one_round 只跑 1 轮，否则只遍历合法 reference
    //    超过序列数阈值的物种保留为 query，但不会作为 reference
    // ------------------------------------------------------------
    uint_t round = only_one_round ? 1 : reference_num;
    std::unordered_set<SpeciesName> processed_reference_species;

    for (uint_t i = 0; i < round; i++) {

        // --------------------------------------------------------
        // 7-1) 当前轮参考物种
        // --------------------------------------------------------
        SpeciesName current_ref_name = reference_order[i];
        spdlog::info("build ref global cache for {}", current_ref_name);

        // 构建 ref_global_cache（用 sampling_interval 采样）
        SequenceUtils::buildRefGlobalCache(seqpro_managers[current_ref_name], sampling_interval, ref_global_cache);

        // --------------------------------------------------------
        // 7-2) 构造本轮参与比对的物种集合：
        //      当前 reference + 尚未作为合法 reference 处理过的物种都作为 query
        //      超过序列数阈值的物种永不进入 processed_reference_species，因此仍作为 query 参与
        // --------------------------------------------------------
        std::unordered_map<SpeciesName, SeqPro::SharedManagerVariant> species_fasta_manager_map;
        for (const auto& query_name : species_order) {
            if (processed_reference_species.count(query_name) > 0) {
                continue;
            }
            auto query_fasta_manager = seqpro_managers.at(query_name);
            species_fasta_manager_map.emplace(query_name, query_fasta_manager);
        }

        // 从第二轮开始允许短 MUM（保持原逻辑：i==0 -> true，否则 false）
        bool allow_short_mum = (i == 0) ? true : false;

        spdlog::info("align multiple genome for {}", current_ref_name);

        // --------------------------------------------------------
        // 7-3) 多物种比对：对当前 ref 与各 query 进行 anchor 搜索
        //      search_mode 固定 ACCURATE_SEARCH（保持原逻辑）
        //      allow_MEM 固定 false
        // --------------------------------------------------------
        SpeciesMatchVec3DPtrMapPtr match_ptr = alignMultipleGenome(
            current_ref_name,
            species_fasta_manager_map,
            ACCURATE_SEARCH,
            fast_build,
            false,
            allow_short_mum,
            ref_global_cache,
            sampling_interval
        );

        spdlog::info("align multiple genome for {} done", current_ref_name);

        // --------------------------------------------------------
        // 7-4) 过滤 anchors：对多个物种的 anchors 聚簇/筛选，得到 cluster_map
        // --------------------------------------------------------
        spdlog::info("filter multiple species anchors for {}", current_ref_name);
        SpeciesClusterMapPtr cluster_map = filterMultipeSpeciesAnchors(
            current_ref_name,
            species_fasta_manager_map,
            match_ptr,
            min_span
        );
        spdlog::info("filter multiple species anchors for {} done", current_ref_name);

        // --------------------------------------------------------
        // 7-5) 构建多基因组图：DP 方式构图（i==0 作为 is_first）
        // --------------------------------------------------------
        spdlog::info("construct multiple genome graphs for {}", current_ref_name);

        constructMultipleGraphsByDp(
            seqpro_managers,
            current_ref_name,
            *cluster_map,
            *multi_graph,
            min_span,
            i == 0
        );

        // --------------------------------------------------------
        // 7-6) 扩展/优化/验证图结构
        // --------------------------------------------------------
        spdlog::info("begin to extend nodes for {}", current_ref_name);
        multi_graph->extendRefNodes(current_ref_name, seqpro_managers, 200);

        multi_graph->optimizeGraphStructure();

#ifdef _DEBUG_
        multi_graph->verifyGraphCorrectness(current_ref_name, true);
#endif // _DEBUG_

        spdlog::info("construct multiple genome graphs for {} done", current_ref_name);

        // --------------------------------------------------------
        // 7-7) 合并本轮生成的多个子图到 multi_graph
        // --------------------------------------------------------
        spdlog::info("merge multiple genome graphs for {}", current_ref_name);
        multi_graph->mergeMultipleGraphs(current_ref_name, thread_num);
        spdlog::info("merge multiple genome graphs for {} done", current_ref_name);

#ifdef _DEBUG_
        multi_graph->verifyGraphCorrectness(true);
#endif

        // 合并后再优化一次
        multi_graph->optimizeGraphStructure();
        spdlog::info("optimize graph genome graphs for {} done", current_ref_name);

        // 标记所有节点已扩展
        multi_graph->markAllExtended();

#ifdef _DEBUG_
        multi_graph->verifyGraphCorrectness(current_ref_name, true, false, false, true, false);
#endif // _DEBUG_

        // --------------------------------------------------------
        // 7-8) 将本轮已对齐区域加入遮蔽区间（mask intervals）
        //      目的：后续轮次避免重复比对已成功对齐的区间
        // --------------------------------------------------------
        spdlog::info("Adding aligned regions as mask intervals for {}", current_ref_name);
        try {
            addAlignedRegionsAsMask(*multi_graph, seqpro_managers, current_ref_name);
            spdlog::info("Successfully added mask intervals for round with reference {}", current_ref_name);
        }
        catch (const std::exception& e) {
            spdlog::error("Failed to add mask intervals for {}: {}", current_ref_name, e.what());
        }

        // --------------------------------------------------------
        // 7-9) 导出 mask intervals 到目录 work_dir/mask_interval/i
        // --------------------------------------------------------
        spdlog::info("[mask-export] Exporting mask intervals captured...");
        FilePath mask_export_dir = work_dir / "mask_interval" / std::to_string(i);
        exportMaskIntervalsToDirectory(mask_export_dir, seqpro_managers);
        processed_reference_species.insert(current_ref_name);
    }

    // 所有轮次完成后，flush logger
    spdlog::default_logger()->flush();

    // 返回最终 multi_graph
    return std::move(multi_graph);
}


// ------------------------------------------------------------
// MultipleRareAligner::alignMultipleGenome
// 功能：以 ref_name 为参考，对 species_fasta_manager_map 中其它物种并行寻找 anchors
// 输入：
// - ref_name：参考物种
// - species_fasta_manager_map：本轮需要参与的物种（包含 ref 与多个 query）
// - search_mode / fast_build / allow_MEM / allow_short_mum：比对参数
// - ref_global_cache / sampling_interval：用于 global 坐标快速定位（加速）
// 输出：
// - SpeciesName -> MatchVec3DPtr（按物种收集的 anchors）
// ------------------------------------------------------------
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

    validateReferenceSequenceCount(ref_name, species_fasta_manager_map.at(ref_name));

    if (species_fasta_manager_map.size() <= 1) {
        spdlog::warn("[alignMultipleQuerys] only reference genome present, nothing to align.");
        return std::make_shared<SpeciesMatchVec3DPtrMap>();
    }

    /* ---------- 1. 结果目录与缓存文件 ---------- */
    // result_dir：按 group_id/round_id 组织输出目录
    FilePath result_dir = work_dir / RESULT_DIR
        / ("group_" + std::to_string(group_id))
        / ("round_" + std::to_string(round_id));
    std::filesystem::create_directories(result_dir);

    // round_id 自增：用于下次调用创建新的 round 目录
    round_id++;

    // anchor_file：本轮 ref_name + search_mode 对应的结果文件（当前读取逻辑已注释）
    FilePath anchor_file = result_dir / (ref_name + "_"
        + SearchModeToString(search_mode) + "." + ANCHOR_EXTENSION);

    /* ---------- 2. 如果已存在结果文件直接读取 ---------- */
    // 当前逻辑被注释：即使文件存在也会重新计算
    //if (std::filesystem::exists(anchor_file)) {
    //    spdlog::info("[alignMultipleQuerys] Load from {}", anchor_file.string());
    //    auto mp = std::make_shared<SpeciesMatchVec3DPtrMap>();
    //    if (loadSpeciesMatchMap(anchor_file, mp))
    //        return mp;
    //}

    /* ---------- 3. 准备参考基因组索引 ---------- */
    FilePath ref_index_path = index_dir / ref_name;
    std::filesystem::create_directories(ref_index_path);

    // PairRareAligner：用于 ref vs query 的 pairwise anchor 查找
    PairRareAligner pra(*this);

    // 构建参考索引（fast_build 控制构建策略）
    pra.buildIndex(ref_name, *species_fasta_manager_map[ref_name], fast_build);
    spdlog::info("[alignMultipleQuerys] reference index built for {}.", ref_name);

    /* ---------- 4. 创建共享线程池 ---------- */
    ThreadPool shared_pool(thread_num);

    /* ---------- 5. 为每个 query 物种启动异步任务 ---------- */
    std::unordered_map<SpeciesName, std::future<MatchVec3DPtr>> fut_map;

    for (auto& kv : species_fasta_manager_map) {
        SpeciesName sp = kv.first;
        if (sp == ref_name) continue; // 跳过参考自身

        std::string prefix = ref_name + "_vs_" + sp;
        auto& fm = kv.second;

        // std::async 并行查找 anchors
        fut_map.emplace(
            sp,
            std::async(std::launch::async,
                [&pra, prefix, &fm, search_mode, allow_MEM, allow_short_mum,
                 &shared_pool, &ref_global_cache, sampling_interval]() -> MatchVec3DPtr {
                    return pra.findQueryFileAnchor(
                        prefix,
                        *fm,
                        search_mode,
                        allow_MEM,
                        allow_short_mum,
                        shared_pool,
                        ref_global_cache,
                        sampling_interval,
                        true
                    );
                })
        );
    }

    // 等待线程池内任务完成（注意：async 的 future 仍需 get 等待结果）
    shared_pool.waitAllTasksDone();

    /* ---------- 6. 收集所有结果 ---------- */
    auto result_map = std::make_shared<SpeciesMatchVec3DPtrMap>();

    size_t total = fut_map.size();
    size_t count = 0;
    size_t next_progress = 1; // 下一次打印进度的阶段（1..20）

    for (auto& kv : fut_map) {
        const SpeciesName& sp = kv.first;
        try {
            // get()：等待该物种异步任务结束并取出 MatchVec3DPtr
            MatchVec3DPtr mv3 = kv.second.get();
            (*result_map)[sp] = std::move(mv3);
            spdlog::info("[alignMultipleQuerys] {} aligned.", sp);
        }
        catch (const std::exception& e) {
            spdlog::error("[alignMultipleQuerys] {} failed: {}", sp, e.what());
        }

        ++count;

        // 进度分 20 段打印（0..100%）
        size_t progress_stage = (count * 20) / total;
        if (progress_stage >= next_progress) {
            int percent = static_cast<int>((progress_stage * 100) / 20);
            spdlog::info("[alignMultipleQuerys] Progress: {}%", percent);
            next_progress = progress_stage + 1;
        }
    }

    return result_map;
}


// ------------------------------------------------------------
// MultipleRareAligner::filterMultipeSpeciesAnchors
// 功能：对 alignMultipleGenome 的结果进行分组与聚簇，得到每个物种的 cluster_map
// 过程：
// 1) 从 species_match_map 收集物种列表（去除 ref）
// 2) 为每个物种预创建 unique/repeat sparse 容器
// 3) 并行 groupMatchByQueryRefSparse：把 anchors 按 (ref,query,chr) 分桶
// 4) 释放 MatchVec3D（节省内存）
// 5) 对每个物种调用 clusterAllChrMatchSparse 聚簇
// ------------------------------------------------------------
SpeciesClusterMapPtr MultipleRareAligner::filterMultipeSpeciesAnchors(
    SpeciesName                       ref_name,
    std::unordered_map<SpeciesName, SeqPro::SharedManagerVariant>& species_fm_map,
    SpeciesMatchVec3DPtrMapPtr        species_match_map,
    uint_t min_span)
{
    if (!species_match_map || species_match_map->empty()) {
        return std::make_shared<SpeciesClusterMap>();
    }

    // unique_map / repeat_map：每个物种分别存 unique/repeat anchors（稀疏 key→MatchVec）
    std::unordered_map<SpeciesName, MatchBySQR_SparsePtr> unique_map;
    std::unordered_map<SpeciesName, MatchBySQR_SparsePtr> repeat_map;

    SpeciesClusterMap cluster_map;

    // 1) 收集 species 列表（不含 ref）
    std::vector<SpeciesName> species_list;
    species_list.reserve(species_match_map->size());
    for (auto& kv : *species_match_map) {
        const SpeciesName& species = kv.first;
        if (species == ref_name) continue;
        species_list.push_back(species);
    }

    // 2) 串行创建容器（避免并发写 unordered_map）
    unique_map.reserve(species_list.size());
    repeat_map.reserve(species_list.size());
    for (const auto& species : species_list) {
        unique_map[species] = std::make_shared<MatchBySQR_Sparse>();
        repeat_map[species] = std::make_shared<MatchBySQR_Sparse>();
    }

    // rfm：参考物种的 SeqPro manager
    auto& rfm = species_fm_map.at(ref_name);

    spdlog::info("Group Match By Query and Ref (Sparse)...");

    // 3) 并行：对每个物种进行 anchors 分组（稀疏结构）
#pragma omp parallel for schedule(dynamic) num_threads(thread_num)
    for (long long idx = 0; idx < (long long)species_list.size(); ++idx) {
        const auto& species = species_list[(size_t)idx];

        MatchVec3DPtr mv3_ptr = species_match_map->at(species);
        auto u_ptr = unique_map.at(species);
        auto r_ptr = repeat_map.at(species);
        auto& qfm = species_fm_map.at(species);

        groupMatchByQueryRefSparse(mv3_ptr, u_ptr, r_ptr, *rfm, *qfm);
    }

    spdlog::info("Group Match By Query and Ref (Sparse) Done");

    // 4) 释放 MatchVec3D（节省内存）
    for (auto& kv : *species_match_map) {
        kv.second.reset();
    }
    species_match_map->clear();

    spdlog::info("Cluster All Chr Match (Sparse) Start...");

    // 5) 串行对每个物种做聚簇（clusterAllChrMatchSparse 内部可并行）
    for (size_t i = 0; i < species_list.size(); ++i) {
        const auto& species = species_list[i];
        auto u_ptr = unique_map.at(species);
        auto r_ptr = repeat_map.at(species);

        auto clusters = clusterAllChrMatchSparse(u_ptr, r_ptr, min_span, thread_num);
        cluster_map.emplace(species, std::move(clusters));
    }

    auto cluster_map_ptr = std::make_shared<SpeciesClusterMap>(std::move(cluster_map));
    spdlog::info("Cluster All Chr Match (Sparse) Done, species num: {}", cluster_map_ptr->size());

    return cluster_map_ptr;
}


// ------------------------------------------------------------
// MultipleRareAligner::constructMultipleGraphsByDp
// 功能：基于过滤后的 cluster_map，为每个 query 物种扩展 anchors 并 DP 过滤，最终构图
// 步骤：
// 1) 收集 species_list（cluster_map 非空的物种）
// 2) Phase-A：对每个物种 extendClusterToAnchorByChr（将 cluster 扩展为 anchor 结构）
// 3) 对每个物种：DP 过滤 anchors -> constructGraphByDP 构图
// ------------------------------------------------------------
void MultipleRareAligner::constructMultipleGraphsByDp(
    std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers,
    SpeciesName ref_name,
    const SpeciesClusterMap& species_cluster_map,
    RaMesh::RaMeshMultiGenomeGraph& graph,
    uint_t min_span, bool is_first)
{
    if (species_cluster_map.empty()) {
        spdlog::warn("[constructMultipleGraphsByDP] Empty species cluster map.");
        return;
    }

    // PairRareAligner：这里复用其 DP 过滤与构图能力
    PairRareAligner pra(*this);
    pra.ref_name = ref_name;

    // 记录参考 manager 的指针（注意：shared_ptr<variant> 解引用得到 variant）
    pra.ref_seqpro_manager = &(*seqpro_managers.at(ref_name));

    /*======================= 1) 收集 species 列表 =======================*/

    std::vector<SpeciesName> species_list;
    species_list.reserve(species_cluster_map.size());

    for (auto& kv : species_cluster_map) {
        if (kv.second) // 只处理非空 cluster map
            species_list.push_back(kv.first);
    }

    if (species_list.empty()) {
        spdlog::warn("[constructMultipleGraphsByDP] No valid species to process.");
        return;
    }

    /*============================================================
      2) Phase-A：extendClusterToAnchorByChr（扩展 cluster 到 anchor）
         - 输入：稀疏结构 ClusterBySQR_SparsePtr（key->clusters）
         - 输出：AnchorBySQR_SparsePtr（同样为稀疏扩展后的 anchor 结构）
         说明：原代码未使用并行（保持原逻辑）
    ============================================================*/

    std::vector<AnchorBySQR_SparsePtr> anchor_results(species_list.size());

    for (long long idx = 0; idx < (long long)species_list.size(); ++idx) {
        const auto& species = species_list[(size_t)idx];
        auto cluster_ptr_sparse = species_cluster_map.at(species);

        // 将 cluster 扩展为 anchor（内部会 decode_sqr_key 并构建对应维度结构）
        AnchorBySQR_SparsePtr result_ptr = pra.extendClusterToAnchorByChr(
            species,
            *seqpro_managers[species],
            cluster_ptr_sparse,
            is_first
        );

        anchor_results[(size_t)idx] = std::move(result_ptr);
    }

    spdlog::info("[constructMultipleGraphsByDP] Phase-A extend done.");

    // ------------------------------------------------------------
    // 3) 对每个物种：DP 过滤 + 构图
    // ------------------------------------------------------------
    for (size_t i = 0; i < species_list.size(); ++i) {
        const auto& species = species_list[i];
        auto& anchor_ptr = anchor_results[i];

        // ref/qry 染色体数量，用于 DP 过滤
        const uint_t ref_chr_cnt = std::visit(
            [](auto& m) { return m->getSequenceCount(); },
            *seqpro_managers[ref_name]
        );
        const uint_t qry_chr_cnt = std::visit(
            [](auto& m) { return m->getSequenceCount(); },
            *seqpro_managers[species]
        );

        if (!anchor_ptr)
            continue;

        // DP 过滤 anchors
        pra.filterAnchorByDP(anchor_ptr, ref_chr_cnt, qry_chr_cnt);
        spdlog::info("DP filter success for {}", species);

        // 用过滤后的 anchors 构建图结构
        pra.constructGraphByDP(
            species,
            *seqpro_managers[species],
            anchor_ptr,
            graph
        );
    }

    spdlog::info("[constructMultipleGraphsByDP] Completed all species.");
}
