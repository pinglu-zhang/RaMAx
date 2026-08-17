#include <sequence_utils.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <stdexcept>
#include <unordered_set>

#include "rare_aligner.h"
#include "anchor.h"  // 包含 UnionFind 定义
#include "SeqPro.h"  // 包含 SeqPro 相关定义
#include "ramesh.h"  // 包含 RaMesh 图结构定义

// 辅助函数：根据CIGAR字符串计算query区间对应关系
namespace {

    constexpr size_t kMaxReferenceSequenceCount = 100000;

    struct OpenMPStageActivity {
        std::atomic<size_t> active{0};
        std::atomic<size_t> max_active{0};

        void enter() {
            const size_t current = active.fetch_add(1,
                std::memory_order_relaxed) + 1;
            size_t observed = max_active.load(std::memory_order_relaxed);
            while (observed < current &&
                   !max_active.compare_exchange_weak(observed, current,
                       std::memory_order_relaxed)) {}
        }

        void leave() {
            active.fetch_sub(1, std::memory_order_relaxed);
        }
    };

    size_t stageWorkerCount(uint_t requested, size_t work_items) {
        if (work_items == 0) return 1;
        return std::max<size_t>(1,
            std::min<size_t>(requested, work_items));
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
    uint_t thread_num_,              // 同理
    uint_t chunk_size_,
    uint_t overlap_size_,
    uint_t min_anchor_length_,
    uint_t max_anchor_frequency_,
    uint_t accurate_skip_threshold_,
    bool trust_legacy_cache_
)
    : work_dir(work_dir_),                                  // 初始化成员
    index_dir(work_dir_ / INDEX_DIR),
    species_path_map(species_path_map_),
    chunk_size(chunk_size_),
    overlap_size(overlap_size_),
    min_anchor_length(min_anchor_length_),
    max_anchor_frequency(max_anchor_frequency_),
    accurate_skip_threshold(accurate_skip_threshold_),
    thread_num(thread_num_),
    trust_legacy_cache(trust_legacy_cache_)
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

    // Estimate whole-genome Mash distance once for the first selected
    // reference. The legacy FM-index aligner remains the only alignment
    // backend in this phase; these records are retained for the future router.
    MashDistanceEstimator mash_estimator(
        locateMashExecutable(), work_dir / "similarity", thread_num);
    first_reference_distances = mash_estimator.estimateFirstReference(
        reference_order.front(), seqpro_managers);

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
        if (merge_exact_contiguous_blocks_enabled) {
            const size_t eliminated_boundaries =
                multi_graph->mergeExactContiguousBlocks(
                    current_ref_name, 1000000, merge_query_gap_max);
            spdlog::debug(
                "[exact-block-merge] round={} reference={} optimization "
                "completed: eliminated_boundaries={} max_query_gap={}",
                i + 1, current_ref_name, eliminated_boundaries,
                merge_query_gap_max);
        }
        if (realign_single_missing_species_enabled) {
            const size_t replaced_windows =
                multi_graph->realignSingleMissingSpeciesWindows(
                    current_ref_name, seqpro_managers,
                    species_mismatch_msa_executable,
                    species_mismatch_realign_max_span,
                    thread_num,
                    species_mismatch_zero_gap_max_span,
                    merge_query_gap_max);
            size_t eliminated_after_realign = 0;
            if (replaced_windows > 0) {
                eliminated_after_realign =
                    multi_graph->mergeExactContiguousBlocks(
                        current_ref_name, 1000000,
                        merge_query_gap_max);
            }
            spdlog::debug(
                "[species-mismatch-realign] round={} reference={} "
                "replaced_windows={} post_realign_eliminated_boundaries={}",
                i + 1, current_ref_name, replaced_windows,
                eliminated_after_realign);
        }
        if (structural_break_repair_options.enabled) {
            auto repair_options = structural_break_repair_options;
            repair_options.parallel_threads = thread_num;
            const auto repair =
                RaMesh::StructuralBreakRepair::repairAnchorBoundedStructuralBreaks(
                    *multi_graph, current_ref_name, seqpro_managers,
                    repair_options);
            size_t eliminated_after_repair = 0;
            if (repair.committed > 0 &&
                merge_exact_contiguous_blocks_enabled) {
                eliminated_after_repair =
                    multi_graph->mergeExactContiguousBlocks(
                        current_ref_name, 1000000,
                        merge_query_gap_max);
            }
            spdlog::debug(
                "[structural-break-repair] round={} reference={} "
                "committed={} post_repair_eliminated_boundaries={}",
                i + 1, current_ref_name, repair.committed,
                eliminated_after_repair);
        }

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

    if (merge_exact_contiguous_blocks_enabled && !reference_order.empty() &&
        spdlog::should_log(spdlog::level::debug)) {
        for (uint_t reference_index = 0; reference_index < round;
             ++reference_index) {
            multi_graph->inspectExactContiguousBlockBoundaries(
                reference_order[reference_index], "final-graph-post-alignment");
        }
    }

    if (short_block_repair_options.enabled) {
        auto options = short_block_repair_options;
        options.parallel_threads = thread_num;
        options.maximum_query_gap = merge_query_gap_max;
        RaMesh::ShortBlockRepair::repairFinalShortBlocks(
            *multi_graph, reference_order, seqpro_managers, options);
    }

    spdlog::info(
        "[cache-summary] FM-index reused={} rebuilt={}",
        index_cache_counters->reused.load(),
        index_cache_counters->rebuilt.load());

    // 所有轮次完成后，flush logger
    spdlog::default_logger()->flush();

    // 返回最终 multi_graph
    return multi_graph;
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

    /* ---------- 4. 收集 query 物种 ---------- */
    std::vector<SpeciesName> query_species;
    query_species.reserve(species_fasta_manager_map.size() - 1);

    for (const auto& kv : species_fasta_manager_map) {
        const SpeciesName& sp = kv.first;
        if (sp == ref_name) continue; // 跳过参考自身
        query_species.push_back(sp);
    }

    /* ---------- 5. 展平全部 (species, chunk, strand) 任务 ---------- */
    std::vector<std::shared_ptr<PreparedAnchorSearch>> plans;
    plans.reserve(query_species.size());
    struct AnchorWorkItem { size_t species_index; size_t task_index; };
    std::vector<AnchorWorkItem> work_items;

    for (const auto& species : query_species) {
        auto& manager = species_fasta_manager_map.at(species);
        auto plan = pra.prepareQueryFileAnchor(
            ref_name + "_vs_" + species, *manager, search_mode,
            allow_MEM, allow_short_mum, ref_global_cache,
            sampling_interval, true);
        const size_t species_index = plans.size();
        for (size_t task_index = 0; task_index < plan->tasks.size(); ++task_index) {
            work_items.push_back({species_index, task_index});
        }
        plans.emplace_back(std::move(plan));
    }

    const size_t workers = stageWorkerCount(thread_num, work_items.size());
    OpenMPStageActivity activity;
    const auto stage_start = std::chrono::steady_clock::now();

    if (workers == 1) {
        for (const auto& item : work_items) {
            activity.enter();
            pra.executePreparedAnchorTask(
                *plans[item.species_index], item.task_index);
            activity.leave();
        }
    } else {
#pragma omp parallel for schedule(dynamic, 1) num_threads(workers)
        for (long long i = 0;
             i < static_cast<long long>(work_items.size()); ++i) {
            const auto& item = work_items[static_cast<size_t>(i)];
            activity.enter();
            pra.executePreparedAnchorTask(
                *plans[item.species_index], item.task_index);
            activity.leave();
        }
    }

    const double stage_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - stage_start).count();
    spdlog::info(
        "[parallel-stage] anchor-search: tasks={}, threads={}, "
        "max_active={}, elapsed_ms={:.3f}",
        work_items.size(), workers,
        activity.max_active.load(std::memory_order_relaxed), stage_ms);

    // Preserve the original species insertion order in the unordered map.
    auto result_map = std::make_shared<SpeciesMatchVec3DPtrMap>();
    for (size_t i = 0; i < query_species.size(); ++i) {
        const auto& species = query_species[i];
        try {
            (*result_map)[species] =
                pra.collectPreparedAnchorSearch(*plans[i]);
            spdlog::info("[alignMultipleQuerys] {} aligned.", species);
        } catch (const std::exception& e) {
            spdlog::error("[alignMultipleQuerys] {} failed: {}",
                species, e.what());
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

    // 5) Flatten all (species, SQR key) clustering work into one team.
    struct ClusterWorkItem { size_t species_index; size_t key_index; };
    std::vector<std::vector<uint64_t>> species_keys(species_list.size());
    std::vector<std::vector<MatchClusterVecPtr>> cluster_results(
        species_list.size());
    std::vector<ClusterWorkItem> cluster_work;

    for (size_t species_index = 0;
         species_index < species_list.size(); ++species_index) {
        const auto& species = species_list[species_index];
        auto& keys = species_keys[species_index];
        const auto& anchors = *unique_map.at(species);
        keys.reserve(anchors.size());
        for (const auto& [key, _] : anchors) keys.push_back(key);
        cluster_results[species_index].resize(keys.size());
        for (size_t key_index = 0; key_index < keys.size(); ++key_index) {
            cluster_work.push_back({species_index, key_index});
        }
    }

    const size_t cluster_workers =
        stageWorkerCount(thread_num, cluster_work.size());
    OpenMPStageActivity cluster_activity;
    const auto cluster_start = std::chrono::steady_clock::now();
    auto run_cluster_item = [&](const ClusterWorkItem& item) {
        const auto& species = species_list[item.species_index];
        const uint64_t key = species_keys[item.species_index][item.key_index];
        MatchVec matches = unique_map.at(species)->at(key);
        cluster_results[item.species_index][item.key_index] =
            clusterChrMatch(matches, min_span);
    };

    if (cluster_workers == 1) {
        for (const auto& item : cluster_work) {
            cluster_activity.enter();
            run_cluster_item(item);
            cluster_activity.leave();
        }
    } else {
#pragma omp parallel for schedule(dynamic, 1) num_threads(cluster_workers)
        for (long long i = 0;
             i < static_cast<long long>(cluster_work.size()); ++i) {
            cluster_activity.enter();
            run_cluster_item(cluster_work[static_cast<size_t>(i)]);
            cluster_activity.leave();
        }
    }

    const double cluster_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - cluster_start).count();
    spdlog::info(
        "[parallel-stage] clustering: tasks={}, threads={}, "
        "max_active={}, elapsed_ms={:.3f}",
        cluster_work.size(), cluster_workers,
        cluster_activity.max_active.load(std::memory_order_relaxed),
        cluster_ms);

    // Rebuild every sparse map in its original key order, then insert species
    // in the original order so downstream iteration remains unchanged.
    for (size_t species_index = 0;
         species_index < species_list.size(); ++species_index) {
        auto clusters = std::make_shared<ClusterBySQR_Sparse>();
        clusters->reserve(species_keys[species_index].size());
        for (size_t key_index = 0;
             key_index < species_keys[species_index].size(); ++key_index) {
            auto& result = cluster_results[species_index][key_index];
            if (result && !result->empty()) {
                (*clusters)[species_keys[species_index][key_index]] =
                    std::move(result);
            }
        }
        cluster_map.emplace(species_list[species_index], std::move(clusters));
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

    struct ExtensionWorkItem { size_t species_index; size_t group_index; };
    std::vector<std::vector<MatchClusterVecPtr>> cluster_groups(
        species_list.size());
    std::vector<std::vector<AnchorPtrVec>> extension_results(
        species_list.size());
    std::vector<ExtensionWorkItem> extension_work;

    for (size_t species_index = 0;
         species_index < species_list.size(); ++species_index) {
        const auto& sparse = *species_cluster_map.at(
            species_list[species_index]);
        auto& groups = cluster_groups[species_index];
        groups.reserve(sparse.size());
        for (const auto& [_, group] : sparse) groups.push_back(group);
        extension_results[species_index].resize(groups.size());
        for (size_t group_index = 0;
             group_index < groups.size(); ++group_index) {
            extension_work.push_back({species_index, group_index});
        }
    }

    const size_t extension_workers =
        stageWorkerCount(thread_num, extension_work.size());
    OpenMPStageActivity extension_activity;
    const auto extension_start = std::chrono::steady_clock::now();
    auto run_extension_item = [&](const ExtensionWorkItem& item) {
        const auto& species = species_list[item.species_index];
        extension_results[item.species_index][item.group_index] =
            pra.extendClusterGroupToAnchors(
                species, *seqpro_managers.at(species),
                cluster_groups[item.species_index][item.group_index],
                is_first);
    };

    if (extension_workers == 1) {
        for (const auto& item : extension_work) {
            extension_activity.enter();
            run_extension_item(item);
            extension_activity.leave();
        }
    } else {
#pragma omp parallel for schedule(dynamic, 1) num_threads(extension_workers)
        for (long long i = 0;
             i < static_cast<long long>(extension_work.size()); ++i) {
            extension_activity.enter();
            run_extension_item(extension_work[static_cast<size_t>(i)]);
            extension_activity.leave();
        }
    }

    std::vector<AnchorBySQR_SparsePtr> anchor_results(species_list.size());
    for (size_t species_index = 0;
         species_index < species_list.size(); ++species_index) {
        auto anchors = std::make_shared<AnchorBySQR_Sparse>();
        anchors->reserve(extension_results[species_index].size());
        for (auto& group_result : extension_results[species_index]) {
            if (!group_result.empty()) {
                anchors->emplace_back(std::move(group_result));
            }
        }
        anchor_results[species_index] = std::move(anchors);
    }

    const double extension_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - extension_start).count();
    spdlog::info(
        "[parallel-stage] cluster-extension: tasks={}, threads={}, "
        "max_active={}, elapsed_ms={:.3f}",
        extension_work.size(), extension_workers,
        extension_activity.max_active.load(std::memory_order_relaxed),
        extension_ms);

    // ------------------------------------------------------------
    // 3) 对每个物种：DP 过滤 + 构图
    // ------------------------------------------------------------
    struct DPWorkItem {
        size_t species_index;
        uint_t chromosome_id;
        bool filter_ref;
    };
    const uint_t ref_chr_cnt = std::visit(
        [](auto& manager) { return manager->getSequenceCount(); },
        *seqpro_managers.at(ref_name));
    std::vector<DPWorkItem> dp_work;
    for (size_t species_index = 0;
         species_index < species_list.size(); ++species_index) {
        if (!anchor_results[species_index]) continue;
        for (uint_t chr = 0; chr < ref_chr_cnt; ++chr) {
            dp_work.push_back({species_index, chr, true});
        }
        const uint_t query_chr_cnt = std::visit(
            [](auto& manager) { return manager->getSequenceCount(); },
            *seqpro_managers.at(species_list[species_index]));
        for (uint_t chr = 0; chr < query_chr_cnt; ++chr) {
            dp_work.push_back({species_index, chr, false});
        }
    }

    const size_t dp_workers = stageWorkerCount(thread_num, dp_work.size());
    OpenMPStageActivity dp_activity;
    const auto dp_start = std::chrono::steady_clock::now();
    auto run_dp_item = [&](const DPWorkItem& item) {
        pra.filterAnchorByDPDimension(
            anchor_results[item.species_index], item.chromosome_id,
            item.filter_ref);
    };

    if (dp_workers == 1) {
        for (const auto& item : dp_work) {
            dp_activity.enter();
            run_dp_item(item);
            dp_activity.leave();
        }
    } else {
#pragma omp parallel for schedule(dynamic, 1) num_threads(dp_workers)
        for (long long i = 0;
             i < static_cast<long long>(dp_work.size()); ++i) {
            dp_activity.enter();
            run_dp_item(dp_work[static_cast<size_t>(i)]);
            dp_activity.leave();
        }
    }

    const double dp_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - dp_start).count();
    spdlog::info(
        "[parallel-stage] anchor-dp: tasks={}, threads={}, "
        "max_active={}, elapsed_ms={:.3f}",
        dp_work.size(), dp_workers,
        dp_activity.max_active.load(std::memory_order_relaxed), dp_ms);

    // Graph insertion remains serial and follows the original species order.
    for (size_t i = 0; i < species_list.size(); ++i) {
        const auto& species = species_list[i];
        auto& anchor_ptr = anchor_results[i];
        if (!anchor_ptr) continue;
        spdlog::info("DP filter success for {}", species);
        pra.constructGraphByDP(
            species,
            *seqpro_managers[species],
            anchor_ptr,
            graph
        );
    }

    spdlog::info("[constructMultipleGraphsByDP] Completed all species.");
}
