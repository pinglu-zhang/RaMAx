#include <sequence_utils.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <exception>
#include <fstream>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <string_view>

#include <omp.h>

#include "rare_aligner.h"
#include "anchor.h"  // 包含 UnionFind 定义
#include "../anchor/anchor_link_internal.h"
#include "process_memory.h"
#include "reference_selection_internal.h"
#include "SeqPro.h"  // 包含 SeqPro 相关定义
#include "ramesh.h"  // 包含 RaMesh 图结构定义

// 辅助函数：根据CIGAR字符串计算query区间对应关系
namespace {

    constexpr size_t kMaxReferenceSequenceCount = 10000;

    void logStageMemory(std::string_view stage,
                        std::string_view event,
                        size_t items = 0,
                        size_t auxiliary = 0) {
        const RaMAxMemory::ProcessMemorySnapshot memory =
            RaMAxMemory::readProcessMemorySnapshot();
        if (!memory.available) {
            spdlog::info(
                "[stage-memory] stage={} event={} available=false items={} auxiliary={}",
                stage, event, items, auxiliary);
            return;
        }
        spdlog::info(
            "[stage-memory] stage={} event={} rss_kib={} peak_rss_kib={} "
            "virtual_kib={} cgroup_limit_bytes={} items={} auxiliary={}",
            stage, event, memory.rss_kib, memory.peak_rss_kib,
            memory.virtual_kib, memory.cgroup_limit_bytes, items, auxiliary);
    }
    struct OpenMPStageActivity {
        std::atomic<size_t> active{0};
        std::atomic<size_t> maximum{0};

        void enter() {
            const size_t current =
                active.fetch_add(
                    1,
                    std::memory_order_relaxed) +
                1;
            size_t observed =
                maximum.load(
                    std::memory_order_relaxed);
            while (observed < current &&
                   !maximum.compare_exchange_weak(
                       observed,
                       current,
                       std::memory_order_relaxed)) {
            }
        }

        void leave() {
            active.fetch_sub(
                1,
                std::memory_order_relaxed);
        }
    };

    size_t stageWorkerCount(
        uint_t requested,
        size_t work_items) {
        if (work_items == 0) {
            return 1;
        }
        return std::max<size_t>(
            1,
            std::min<size_t>(
                requested,
                work_items));
    }

    template <typename Function>
    void executeParallelStage(
        size_t work_items,
        uint_t requested_threads,
        OpenMPStageActivity& activity,
        Function&& function,
        int dynamic_chunk = 1) {
        std::mutex failure_mutex;
        std::exception_ptr first_failure;
        size_t first_failure_index = work_items;
        const size_t workers =
            stageWorkerCount(
                requested_threads,
                work_items);
        const auto execute =
            [&](size_t index) {
                activity.enter();
                try {
                    function(index);
                } catch (...) {
                    std::lock_guard<std::mutex> lock(
                        failure_mutex);
                    if (index < first_failure_index) {
                        first_failure_index = index;
                        first_failure =
                            std::current_exception();
                    }
                }
                activity.leave();
            };
        if (workers == 1 || omp_in_parallel()) {
            for (size_t index = 0;
                 index < work_items;
                 ++index) {
                execute(index);
            }
        } else {
#pragma omp parallel for schedule(dynamic, dynamic_chunk) num_threads(workers)
            for (long long index = 0;
                 index <
                     static_cast<long long>(
                         work_items);
                 ++index) {
                execute(
                    static_cast<size_t>(
                        index));
            }
        }
        if (first_failure) {
            std::rethrow_exception(first_failure);
        }
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
        const SeqPro::SharedManagerVariant& manager_variant,
        bool allow_reference_count_fallback) {

        const size_t sequence_count = getManagerSequenceCount(manager_variant);
        if (isReferenceEligibleSequenceCount(sequence_count) ||
            (allow_reference_count_fallback && sequence_count > 0)) {
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

    logStageMemory("aligned-mask", "start", graph.blocks.size());
    const auto collect_started = std::chrono::steady_clock::now();
    using ChromosomeJournal =
        std::unordered_map<ChrName, std::vector<RaMesh::Segment*>>;
    std::unordered_map<SpeciesName, ChromosomeJournal>
        species_chr_journal;

    size_t valid_blocks = 0;
    size_t scanned_segments = 0;
    size_t changed_segments = 0;
    size_t unchanged_segments = 0;

    for (const auto& weak_block : graph.blocks) {
        auto block_ptr = weak_block.lock();
        if (!block_ptr) continue;

        ++valid_blocks;
        std::shared_lock block_lock(block_ptr->rw);

        for (const auto& [species_chr_pair, segment] : block_ptr->anchors) {
            const auto& [species_name, chr_name] = species_chr_pair;
            ++scanned_segments;
            if (!segment || !segment->isSegment() || segment->length == 0) {
                continue;
            }

            if (!seqpro_managers.contains(species_name)) {
                continue;
            }

            if (segment->mask_journal_length != 0 &&
                segment->mask_journal_start == segment->start &&
                segment->mask_journal_length == segment->length) {
                ++unchanged_segments;
                continue;
            }

            species_chr_journal[species_name][chr_name].push_back(
                segment.get());
            ++changed_segments;
        }
    }

    const double collect_seconds =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - collect_started)
            .count();
    size_t committed_segments = 0;
    size_t missing_sequences = 0;
    size_t touched_chromosomes = 0;

    for (auto& [species_name, chr_journals] : species_chr_journal) {
        const auto species_started = std::chrono::steady_clock::now();
        try {
            auto* masked_manager =
                ensureMaskedManager(seqpro_managers.at(species_name));
            std::vector<SeqPro::MaskIntervalDelta> deltas;
            deltas.reserve(chr_journals.size());
            std::vector<RaMesh::Segment*> publish_journal;

            size_t species_candidates = 0;
            for (const auto& [unused_chr_name, journal] : chr_journals) {
                (void)unused_chr_name;
                species_candidates += journal.size();
            }
            publish_journal.reserve(species_candidates);

            for (auto& [chr_name, journal] : chr_journals) {
                if (journal.empty()) {
                    continue;
                }

                const SeqPro::SequenceId seq_id =
                    masked_manager->getSequenceId(chr_name);
                if (seq_id == SeqPro::SequenceIndex::INVALID_ID) {
                    spdlog::warn(
                        "[addAlignedRegionsAsMask] Sequence not found: {}:{}, skipping",
                        species_name, chr_name);
                    ++missing_sequences;
                    continue;
                }

                SeqPro::MaskIntervalDelta delta;
                delta.sequence_id = seq_id;
                delta.intervals.reserve(journal.size());
                for (RaMesh::Segment* segment : journal) {
                    delta.intervals.emplace_back(
                        segment->start,
                        segment->start + segment->length);
                }
                publish_journal.insert(
                    publish_journal.end(), journal.begin(), journal.end());
                deltas.push_back(std::move(delta));
            }

            const SeqPro::MaskBatchMergeStats stats =
                masked_manager->applyFinalizedMaskDeltas(
                    std::move(deltas));
            for (RaMesh::Segment* segment : publish_journal) {
                segment->mask_journal_start = segment->start;
                segment->mask_journal_length = segment->length;
            }

            committed_segments += publish_journal.size();
            touched_chromosomes += stats.touched_sequences;

            const double species_seconds =
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - species_started)
                    .count();
            const RaMAxMemory::ProcessMemorySnapshot memory =
                RaMAxMemory::readProcessMemorySnapshot();
            if (memory.available) {
                spdlog::info(
                    "[mask-finalize] species={} incoming_intervals={} "
                    "normalized_delta_intervals={} previous_intervals={} "
                    "final_intervals={} touched_chromosomes={} "
                    "sort_seconds={:.6f} merge_seconds={:.6f} "
                    "metadata_seconds={:.6f} wall_seconds={:.6f} "
                    "rss_kib={} peak_rss_kib={}",
                    species_name, stats.incoming_intervals,
                    stats.normalized_delta_intervals,
                    stats.previous_intervals, stats.final_intervals,
                    stats.touched_sequences, stats.sort_seconds,
                    stats.merge_seconds, stats.metadata_seconds,
                    species_seconds, memory.rss_kib, memory.peak_rss_kib);
            } else {
                spdlog::info(
                    "[mask-finalize] species={} incoming_intervals={} "
                    "normalized_delta_intervals={} previous_intervals={} "
                    "final_intervals={} touched_chromosomes={} "
                    "sort_seconds={:.6f} merge_seconds={:.6f} "
                    "metadata_seconds={:.6f} wall_seconds={:.6f} "
                    "memory_available=false",
                    species_name, stats.incoming_intervals,
                    stats.normalized_delta_intervals,
                    stats.previous_intervals, stats.final_intervals,
                    stats.touched_sequences, stats.sort_seconds,
                    stats.merge_seconds, stats.metadata_seconds,
                    species_seconds);
            }

            spdlog::info(
                "[addAlignedRegionsAsMask] Successfully added {} mask intervals for species {}",
                publish_journal.size(), species_name);
        } catch (const std::exception& e) {
            spdlog::error(
                "[addAlignedRegionsAsMask] Error processing species {}: {}",
                species_name, e.what());
        }

        ChromosomeJournal released;
        chr_journals.swap(released);
    }

    const RaMAxMemory::ProcessMemorySnapshot memory =
        RaMAxMemory::readProcessMemorySnapshot();
    if (memory.available) {
        spdlog::info(
            "[mask-journal] reference={} scanned_blocks={} scanned_segments={} "
            "changed_candidates={} appended={} unchanged_skipped={} "
            "missing_sequences={} touched_species={} touched_chromosomes={} "
            "collect_seconds={:.6f} rss_kib={} peak_rss_kib={}",
            ref_name, valid_blocks, scanned_segments, changed_segments,
            committed_segments, unchanged_segments, missing_sequences,
            species_chr_journal.size(), touched_chromosomes,
            collect_seconds, memory.rss_kib, memory.peak_rss_kib);
    } else {
        spdlog::info(
            "[mask-journal] reference={} scanned_blocks={} scanned_segments={} "
            "changed_candidates={} appended={} unchanged_skipped={} "
            "missing_sequences={} touched_species={} touched_chromosomes={} "
            "collect_seconds={:.6f} memory_available=false",
            ref_name, valid_blocks, scanned_segments, changed_segments,
            committed_segments, unchanged_segments, missing_sequences,
            species_chr_journal.size(), touched_chromosomes,
            collect_seconds);
    }
    logStageMemory(
        "aligned-mask", "complete",
        committed_segments, unchanged_segments);
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

    const size_t fallback_index =
        RaMAxReferenceSelection::enableMinimumSequenceFallback(species_qualities);
    const bool using_fallback_reference = fallback_index != species_qualities.size();
    if (using_fallback_reference) {
        const auto& fallback = species_qualities[fallback_index];
        spdlog::warn(
            "[reference-selection] No genome satisfies max_allowed={}. "
            "Using minimum-sequence genome {} (sequence_count={}) as the sole "
            "fallback reference for one alignment round.",
            kMaxReferenceSequenceCount, fallback.name, fallback.sequence_count);
        if (!ref_name.empty() && ref_name != fallback.name) {
            spdlog::warn(
                "[reference-selection] Requested reference {} is overridden by "
                "minimum-sequence fallback {}.", ref_name, fallback.name);
        }
    }

    // ------------------------------------------------------------
    // 3) 将指定的 ref_name 移动到最前面（如果存在于列表）
    //    有常规候选时，显式 reference 超限仍报错；兜底模式固定选择最少序列者
    // ------------------------------------------------------------
    auto ref_it = std::find_if(
        species_qualities.begin(),
        species_qualities.end(),
        [&ref_name](const SpeciesAssemblyQuality& quality) {
            return quality.name == ref_name;
        }
    );
    if (!using_fallback_reference && ref_it != species_qualities.end()) {
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
    //    超限物种保留为 query；全部超限时仅兜底物种进入 reference_order
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
            "[reference-selection] No nonempty reference genome found; "
            "minimum-sequence fallback is unavailable.");
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
    // reference. The native suffix-array aligner remains the only legacy
    // backend in this phase; these records are retained for the future router.
    MashDistanceEstimator mash_estimator(
        locateMashExecutable(), work_dir / "similarity", thread_num);
    first_reference_distances = mash_estimator.estimateFirstReference(
        reference_order.front(), seqpro_managers);

    // Route only first-round near genomes through wfmash. Tool/version
    // validation and the public-reference FAI are intentionally completed
    // before any native suffix-array index can be built.
    FirstRoundWfmashRouter wfmash_router(
        locateSamtoolsExecutable(), locateWfmashExecutable(),
        work_dir / "wfmash" / "round_0", thread_num);
    const FirstRoundWfmashResult first_round_wfmash = wfmash_router.run(
        reference_order.front(), first_reference_distances,
        near_distance_threshold, far_distance_threshold, seqpro_managers);
    spdlog::info(
        "[wfmash-router] tools: samtools={} wfmash={} successful_pairs={}",
        wfmash_router.samtoolsVersion(), wfmash_router.wfmashVersion(),
        first_round_wfmash.successful_species.size());

    // ------------------------------------------------------------
    // 5) 初始化参考缓存 ref_global_cache
    //    - 用于加速 global 坐标定位到序列（避免频繁二分/查询）
    // ------------------------------------------------------------
    sdsl::int_vector<0> ref_global_cache;

    // ------------------------------------------------------------
    // 6) 创建多基因组图（整个 starAlignment 过程使用同一个 multi_graph）
    // ------------------------------------------------------------
    auto multi_graph = std::make_unique<RaMesh::RaMeshMultiGenomeGraph>(seqpro_managers);
    multi_graph->reference_order = reference_order;

    // ------------------------------------------------------------
    // 7) 决定轮数：only_one_round 只跑 1 轮，否则只遍历合法 reference
    //    兜底模式 reference_num 为 1，因此也只运行一轮
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
            if (i == 0 && query_name != current_ref_name &&
                first_round_wfmash.successful_species.contains(query_name)) {
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
        const size_t legacy_query_count =
            species_fasta_manager_map.size() -
            (species_fasta_manager_map.contains(current_ref_name) ? 1 : 0);
        if (legacy_query_count > 0) {
            logStageMemory(
                "anchor-search", "start", legacy_query_count);
            SpeciesMatchVec3DPtrMapPtr match_ptr = alignMultipleGenome(
                current_ref_name,
                species_fasta_manager_map,
                ACCURATE_SEARCH,
                fast_build,
                allow_mem,
                allow_short_mum,
                ref_global_cache,
                sampling_interval,
                using_fallback_reference
            );
            logStageMemory(
                "anchor-search", "complete", match_ptr ? match_ptr->size() : 0);

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
        } else {
            if (i == 0) {
                spdlog::info(
                    "[wfmash-router] all first-round queries succeeded; skipping suffix-array index, cluster, and DP for {}",
                    current_ref_name);
            } else {
                spdlog::info(
                    "No remaining legacy queries for {}; skipping suffix-array index, cluster, and DP",
                    current_ref_name);
            }
        }

        if (i == 0) {
            auto& reference_manager = *seqpro_managers.at(current_ref_name);
            for (const auto& [query_name, anchors] :
                 first_round_wfmash.anchors_by_species) {
                auto& query_manager = *seqpro_managers.at(query_name);
                for (const auto& anchor : anchors) {
                    multi_graph->insertAnchorIntoGraph(
                        reference_manager, query_manager, current_ref_name,
                        query_name, anchor, true);
                }
                spdlog::info(
                    "[wfmash-router] imported {} anchors for {} before graph extension",
                    anchors.size(), query_name);
            }
        }

        // --------------------------------------------------------
        // 7-6) 扩展/优化/验证图结构
        // --------------------------------------------------------
        spdlog::info("begin to extend nodes for {}", current_ref_name);
        logStageMemory("extend-ref-nodes", "start");
        multi_graph->extendRefNodes(current_ref_name, seqpro_managers, 200);
        logStageMemory("extend-ref-nodes", "complete");

        multi_graph->optimizeGraphStructure();

#ifdef _DEBUG_
        multi_graph->verifyGraphCorrectness(current_ref_name, true);
#endif // _DEBUG_

        spdlog::info("construct multiple genome graphs for {} done", current_ref_name);

        // --------------------------------------------------------
        // 7-7) 合并本轮生成的多个子图到 multi_graph
        // --------------------------------------------------------
        spdlog::info("merge multiple genome graphs for {}", current_ref_name);
        logStageMemory("graph-merge", "start", multi_graph->blocks.size());
        multi_graph->mergeMultipleGraphs(current_ref_name, thread_num);
        logStageMemory("graph-merge", "complete", multi_graph->blocks.size());
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
        "[cache-summary] suffix-array storage=memory-only built={} "
        "disk-reused=0 disk-bytes=0",
        index_cache_counters->memory_only_built.load());

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
    SeqPro::Length sampling_interval,
    bool allow_reference_count_fallback)
{
    secondary_match_map.clear();
    /* ---------- 0. 合法性检查 ---------- */
    if (!species_fasta_manager_map.count(ref_name))
        throw std::runtime_error("[alignMultipleQuerys] reference species not found: " + ref_name);

    validateReferenceSequenceCount(ref_name, species_fasta_manager_map.at(ref_name),
                                   allow_reference_count_fallback);

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

    struct SpeciesSearchPlans {
        std::shared_ptr<
            PreparedAnchorSearch>
            primary;
        std::shared_ptr<
            PreparedAnchorSearch>
            repeat_masked;
        std::shared_ptr<
            PreparedAnchorSearch>
            repeat_full;
    };
    struct AnchorSearchWorkItem {
        PreparedAnchorSearch* plan =
            nullptr;
        size_t task_index = 0;
    };
    std::vector<SpeciesSearchPlans>
        plans(query_species.size());
    std::vector<AnchorSearchWorkItem>
        work_items;
    const auto append_plan =
        [&](const std::shared_ptr<
                PreparedAnchorSearch>&
                plan) {
            for (size_t task_index = 0;
                 task_index <
                     plan->tasks.size();
                 ++task_index) {
                work_items.push_back(
                    {plan.get(),
                     task_index});
            }
        };
    for (size_t species_index = 0;
         species_index <
             query_species.size();
         ++species_index) {
        const auto& species =
            query_species[species_index];
        auto& manager =
            species_fasta_manager_map.at(
                species);
        auto& species_plans =
            plans[species_index];
        species_plans.primary =
            pra.prepareQueryFileAnchor(
                *manager,
                search_mode,
                false,
                allow_short_mum,
                ref_global_cache,
                sampling_interval,
                true,
                false);
        append_plan(
            species_plans.primary);
        if (allow_MEM) {
            species_plans.repeat_masked =
                pra.prepareQueryFileAnchor(
                    *manager,
                    search_mode,
                    true,
                    allow_short_mum,
                    ref_global_cache,
                    sampling_interval,
                    true,
                    false);
            species_plans.repeat_full =
                pra.prepareQueryFileAnchor(
                    *manager,
                    search_mode,
                    true,
                    allow_short_mum,
                    ref_global_cache,
                    sampling_interval,
                    true,
                    true);
            append_plan(
                species_plans
                    .repeat_masked);
            append_plan(
                species_plans
                    .repeat_full);
        }
    }

    OpenMPStageActivity activity;
    const auto stage_start =
        std::chrono::steady_clock::now();
    executeParallelStage(
        work_items.size(),
        thread_num,
        activity,
        [&](size_t work_index) {
            const auto& item =
                work_items[work_index];
            pra.executePreparedAnchorTask(
                *item.plan,
                item.task_index);
        });
    const double elapsed_ms =
        std::chrono::duration<
            double,
            std::milli>(
            std::chrono::steady_clock::now() -
            stage_start)
            .count();
    spdlog::info(
        "[parallel-stage] anchor-search: "
        "tasks={}, threads={}, max_active={}, "
        "elapsed_ms={:.3f}",
        work_items.size(),
        stageWorkerCount(
            thread_num,
            work_items.size()),
        activity.maximum.load(
            std::memory_order_relaxed),
        elapsed_ms);

    auto result_map =
        std::make_shared<
            SpeciesMatchVec3DPtrMap>();
    for (size_t species_index = 0;
         species_index <
             query_species.size();
         ++species_index) {
        const auto& species =
            query_species[species_index];
        auto& species_plans =
            plans[species_index];
        try {
            (*result_map)[species] =
                pra.collectPreparedAnchorSearch(
                    *species_plans.primary);
            if (allow_MEM) {
                auto repeat_anchors =
                    pra.collectPreparedAnchorSearch(
                        *species_plans
                             .repeat_masked);
                auto full_repeat_anchors =
                    pra.collectPreparedAnchorSearch(
                        *species_plans
                             .repeat_full);
                repeat_anchors->insert(
                    repeat_anchors->end(),
                    std::make_move_iterator(
                        full_repeat_anchors
                            ->begin()),
                    std::make_move_iterator(
                        full_repeat_anchors
                            ->end()));
                secondary_match_map[species] =
                    std::move(
                        repeat_anchors);
            }
            spdlog::info(
                "[alignMultipleQuerys] {} aligned.",
                species);
        } catch (const std::exception& e) {
            spdlog::error(
                "[alignMultipleQuerys] {} failed: {}",
                species,
                e.what());
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
    logStageMemory(
        "sparse-clustering", "start", species_match_map->size());

    // unique_map / repeat_map：每个物种分别存 unique/repeat anchors（稀疏 key→MatchVec）
    std::unordered_map<SpeciesName, MatchBySQR_SparsePtr> unique_map;
    std::unordered_map<SpeciesName, MatchBySQR_SparsePtr> repeat_map;

    SpeciesClusterMap cluster_map;
    secondary_cluster_map.clear();

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
        if (const auto secondary_it =
                secondary_match_map.find(species);
            secondary_it != secondary_match_map.end()) {
            auto secondary_unique =
                std::make_shared<MatchBySQR_Sparse>();
            MatchVec3DPtr secondary_matches =
                secondary_it->second;
            groupMatchByQueryRefSparse(
                secondary_matches, secondary_unique, r_ptr,
                *rfm, *qfm);
        }
    }


    spdlog::info("Group Match By Query and Ref (Sparse) Done");

    // 4) 释放 MatchVec3D（节省内存）
    for (auto& kv : *species_match_map) {
        kv.second.reset();
    }
    species_match_map->clear();
    secondary_match_map.clear();

    spdlog::info("Cluster All Chr Match (Sparse) Start...");

    const auto cluster_sparse_maps =
        [&](const auto& source_maps,
            SpeciesClusterMap& destination,
            std::string_view stage_name) {
            struct ClusterWorkItem {
                size_t species_index = 0;
                size_t key_index = 0;
            };
            std::vector<std::vector<uint64_t>>
                keys(species_list.size());
            std::vector<std::vector<
                MatchClusterVecPtr>>
                results(species_list.size());
            std::vector<ClusterWorkItem>
                work_items;
            for (size_t species_index = 0;
                 species_index <
                     species_list.size();
                 ++species_index) {
                const auto& source =
                    *source_maps.at(
                        species_list[
                            species_index]);
                auto& species_keys =
                    keys[species_index];
                species_keys.reserve(
                    source.size());
                for (const auto& [key, unused] :
                     source) {
                    (void)unused;
                    species_keys.push_back(
                        key);
                }
                results[species_index].resize(
                    species_keys.size());
                for (size_t key_index = 0;
                     key_index <
                         species_keys.size();
                     ++key_index) {
                    work_items.push_back(
                        {species_index,
                         key_index});
                }
            }

            OpenMPStageActivity activity;
            const auto stage_start =
                std::chrono::steady_clock::now();
            executeParallelStage(
                work_items.size(),
                thread_num,
                activity,
                [&](size_t work_index) {
                    const auto& item =
                        work_items[
                            work_index];
                    const auto& species =
                        species_list[
                            item.species_index];
                    const uint64_t key =
                        keys[item.species_index]
                            [item.key_index];
                    auto& matches =
                        source_maps.at(
                            species)
                            ->at(key);
                    results[item.species_index]
                           [item.key_index] =
                        clusterChrMatch(
                            matches,
                            min_span);
                });
            const double elapsed_ms =
                std::chrono::duration<
                    double,
                    std::milli>(
                    std::chrono::steady_clock::now() -
                    stage_start)
                    .count();
            spdlog::info(
                "[parallel-stage] {}: tasks={}, "
                "threads={}, max_active={}, "
                "elapsed_ms={:.3f}",
                stage_name,
                work_items.size(),
                stageWorkerCount(
                    thread_num,
                    work_items.size()),
                activity.maximum.load(
                    std::memory_order_relaxed),
                elapsed_ms);

            for (size_t species_index = 0;
                 species_index <
                     species_list.size();
                 ++species_index) {
                auto clusters =
                    std::make_shared<
                        ClusterBySQR_Sparse>();
                clusters->reserve(
                    keys[species_index].size());
                for (size_t key_index = 0;
                     key_index <
                         keys[species_index]
                             .size();
                     ++key_index) {
                    auto& result =
                        results[species_index]
                               [key_index];
                    if (result &&
                        !result->empty()) {
                        (*clusters)
                            [keys[species_index]
                                 [key_index]] =
                            std::move(result);
                    }
                }
                destination.emplace(
                    species_list[
                        species_index],
                    std::move(clusters));
            }
        };

    cluster_sparse_maps(
        unique_map,
        cluster_map,
        "clustering-primary");
    if (allow_mem) {
        cluster_sparse_maps(
            repeat_map,
            secondary_cluster_map,
            "clustering-secondary");
    }
    auto cluster_map_ptr = std::make_shared<SpeciesClusterMap>(std::move(cluster_map));
    spdlog::info("Cluster All Chr Match (Sparse) Done, species num: {}", cluster_map_ptr->size());
    logStageMemory(
        "sparse-clustering", "complete", cluster_map_ptr->size());

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

    struct ExtensionWorkItem {
        uint32_t species_index = 0;
        uint32_t group_index = 0;
        MatchClusterVec* group = nullptr;
        uint64_t cost = 0;
        bool secondary = false;
        bool heavy = false;
    };
    std::vector<SeqPro::ManagerVariant*> species_managers;
    species_managers.reserve(species_list.size());
    for (const auto& species : species_list) {
        species_managers.push_back(seqpro_managers.at(species).get());
    }
    std::vector<std::vector<AnchorPtrVec>>
        primary_extensions(
            species_list.size());
    std::vector<std::vector<AnchorPtrVec>>
        secondary_extensions(
            species_list.size());
    std::vector<std::vector<AnchorVec>>
        primary_materialized(species_list.size());
    std::vector<std::vector<AnchorVec>>
        secondary_materialized(species_list.size());

    std::vector<ExtensionWorkItem>
        extension_work;

    const auto saturatingAdd = [](uint64_t left, uint64_t right) {
        return right > std::numeric_limits<uint64_t>::max() - left
            ? std::numeric_limits<uint64_t>::max()
            : left + right;
    };
    const auto estimateExtensionCost = [&](const MatchClusterVec& group) {
        uint64_t cost = 0;
        bool heavy = false;
        uint64_t total_matches = 0;
        for (const auto& cluster : group) {
            total_matches = saturatingAdd(
                total_matches, cluster.size());
        }
        cost = saturatingAdd(cost, total_matches);
        if (!is_first) {
            return std::pair<uint64_t, bool>{cost, false};
        }
        const uint64_t candidate_width =
            std::min<size_t>(group.size(), 2000);
        const uint64_t candidate_cost = candidate_width != 0 &&
                group.size() >
                    std::numeric_limits<uint64_t>::max() / candidate_width
            ? std::numeric_limits<uint64_t>::max()
            : static_cast<uint64_t>(group.size()) * candidate_width;
        cost = saturatingAdd(cost, candidate_cost);
        heavy = group.size() > 1;
        for (const auto& cluster : group) {
            for (size_t index = 1; index < cluster.size(); ++index) {
                const Match& previous = cluster[index - 1];
                const Match& current = cluster[index];
                const int64_t ref_gap = static_cast<int64_t>(current.ref_start) -
                    static_cast<int64_t>(previous.ref_start + previous.match_len());
                int64_t query_gap = 0;
                if (previous.strand() == FORWARD) {
                    query_gap = static_cast<int64_t>(current.qry_start) -
                        static_cast<int64_t>(previous.qry_start + previous.match_len());
                } else {
                    query_gap = static_cast<int64_t>(previous.qry_start) -
                        static_cast<int64_t>(current.qry_start + current.match_len());
                }
                if (ref_gap < 0 || query_gap < 0 ||
                    ref_gap > 10000 || query_gap > 10000) {
                    continue;
                }
                heavy = true;
                const uint64_t maximum_gap = static_cast<uint64_t>(
                    std::max(ref_gap, query_gap));
                const uint64_t band = static_cast<uint64_t>(auto_band(
                    static_cast<int>(ref_gap),
                    static_cast<int>(query_gap)));
                const uint64_t width = saturatingAdd(band, band) + 1;
                const uint64_t alignment_cost = maximum_gap != 0 &&
                        width > std::numeric_limits<uint64_t>::max() / maximum_gap
                    ? std::numeric_limits<uint64_t>::max()
                    : maximum_gap * width;
                cost = saturatingAdd(cost, alignment_cost);
            }
        }
        return std::pair<uint64_t, bool>{cost, heavy};
    };

    for (size_t species_index = 0;
         species_index <
             species_list.size();
         ++species_index) {
        const auto& species =
            species_list[species_index];
        const auto& primary_map =
            *species_cluster_map.at(
                species);
        primary_extensions[species_index]
            .resize(primary_map.size());
        size_t primary_index = 0;
        primary_materialized[species_index]
            .resize(primary_map.size());
        for (const auto& [key, group] :
             primary_map) {
            (void)key;
            const size_t output_index = primary_index++;
            if (!group || group->empty()) continue;
            const auto [cost, heavy] = estimateExtensionCost(*group);
            extension_work.push_back({
                static_cast<uint32_t>(species_index),
                static_cast<uint32_t>(output_index),
                group.get(), cost, false, heavy});
        }

        if (allow_mem) {
            const auto secondary_it =
                secondary_cluster_map.find(
                    species);
            if (secondary_it !=
                    secondary_cluster_map.end() &&
                secondary_it->second) {
                secondary_extensions[
                    species_index]
                    .resize(
                        secondary_it->second->size());
                size_t secondary_index = 0;
                secondary_materialized[
                    species_index]
                    .resize(
                        secondary_it->second->size());
                for (const auto& [key, group] :
                      *secondary_it->second) {
                    (void)key;
                    const size_t output_index = secondary_index++;
                    if (!group || group->empty()) continue;
                    const auto [cost, heavy] = estimateExtensionCost(*group);
                    extension_work.push_back({
                        static_cast<uint32_t>(species_index),
                        static_cast<uint32_t>(output_index),
                        group.get(), cost, true, heavy});
                }
            }
        }
    }

    std::sort(
        extension_work.begin(),
        extension_work.end(),
        [](const ExtensionWorkItem& left,
           const ExtensionWorkItem& right) {
            if (left.cost != right.cost) {
                return left.cost > right.cost;
            }
            if (left.species_index != right.species_index) {
                return left.species_index < right.species_index;
            }
            if (left.secondary != right.secondary) {
                return left.secondary < right.secondary;
            }
            return left.group_index < right.group_index;
        });
    const size_t heavy_work_items = static_cast<size_t>(
        std::count_if(extension_work.begin(), extension_work.end(),
            [](const ExtensionWorkItem& item) { return item.heavy; }));

    OpenMPStageActivity
        extension_activity;
    logStageMemory(
        "cluster-extension", "start", extension_work.size());
    const auto extension_start =
        std::chrono::steady_clock::now();
    const auto materialization_start = extension_start;
    std::atomic<size_t> materialized_groups{0};
    std::mutex materialization_progress_mutex;
    auto last_materialization_progress = materialization_start;
    executeParallelStage(
        extension_work.size(),
        thread_num,
        extension_activity,
        [&](size_t work_index) {
            const auto& item =
                extension_work[work_index];
            auto& manager = *species_managers[item.species_index];
            if (is_first) {
                AnchorVec& destination = item.secondary
                    ? secondary_materialized[item.species_index]
                          [item.group_index]
                    : primary_materialized[item.species_index]
                          [item.group_index];
                destination = AnchorLinkDetail::materializeClusterAnchors(
                    *item.group, *pra.ref_seqpro_manager, manager);
            } else if (item.secondary) {
                secondary_extensions[
                    item.species_index]
                    [item.group_index] =
                    pra.extendClusterGroupToAnchors(
                        manager,
                        *item.group,
                        false);
            } else {
                primary_extensions[
                    item.species_index]
                    [item.group_index] =
                    pra.extendClusterGroupToAnchors(
                        manager,
                        *item.group,
                        false);
            }
            const size_t completed = materialized_groups.fetch_add(
                1, std::memory_order_relaxed) + 1;
            if ((completed & 255U) == 0) {
                std::unique_lock<std::mutex> lock(
                    materialization_progress_mutex, std::try_to_lock);
                const auto now = std::chrono::steady_clock::now();
                if (lock.owns_lock() &&
                    now - last_materialization_progress >=
                        std::chrono::minutes(5)) {
                    last_materialization_progress = now;
                    spdlog::info(
                        "[cluster-extension-progress] phase=materialize "
                        "groups={}/{} active={} max_active={}",
                        completed, extension_work.size(),
                        extension_activity.active.load(
                            std::memory_order_relaxed),
                        extension_activity.maximum.load(
                            std::memory_order_relaxed));
                    logStageMemory(
                        "cluster-extension-materialize", "progress",
                        completed, extension_work.size());
                }
            }
        }, 1);
    const double materialization_seconds =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - materialization_start).count();
    logStageMemory(
        "cluster-extension-materialize", "complete",
        materialized_groups.load(std::memory_order_relaxed),
        extension_work.size());

    struct ComponentWorkItem {
        uint32_t species_index{0};
        uint32_t group_index{0};
        uint32_t component_index{0};
        size_t begin{0};
        size_t end{0};
        uint64_t cost{0};
        size_t output_slot{0};
        bool secondary{false};
    };
    std::vector<ComponentWorkItem> component_work;
    std::vector<AnchorPtrVec> component_outputs;
    std::vector<std::vector<std::vector<size_t>>> primary_component_slots(
        species_list.size());
    std::vector<std::vector<std::vector<size_t>>> secondary_component_slots(
        species_list.size());
    size_t component_count = 0;
    AnchorLinkDetail::Statistics split_statistics;

    if (is_first) {
        for (size_t species_index = 0;
             species_index < species_list.size(); ++species_index) {
            primary_component_slots[species_index].resize(
                primary_materialized[species_index].size());
            secondary_component_slots[species_index].resize(
                secondary_materialized[species_index].size());
            const auto append_components = [&](AnchorVec& anchors,
                                               size_t group_index,
                                               bool secondary) {
                const auto ranges =
                    AnchorLinkDetail::splitAnchorComponents(
                        anchors, &split_statistics);
                auto& slots = secondary
                    ? secondary_component_slots[species_index][group_index]
                    : primary_component_slots[species_index][group_index];
                slots.reserve(ranges.size());
                for (size_t component_index = 0;
                     component_index < ranges.size(); ++component_index) {
                    const auto& range = ranges[component_index];
                    const size_t output_slot = component_outputs.size();
                    component_outputs.emplace_back();
                    slots.push_back(output_slot);
                    component_work.push_back({
                        static_cast<uint32_t>(species_index),
                        static_cast<uint32_t>(group_index),
                        static_cast<uint32_t>(component_index),
                        range.begin, range.end, range.estimated_cost,
                        output_slot, secondary});
                }
            };
            for (size_t group = 0;
                 group < primary_materialized[species_index].size(); ++group) {
                append_components(
                    primary_materialized[species_index][group], group, false);
            }
            for (size_t group = 0;
                 group < secondary_materialized[species_index].size(); ++group) {
                append_components(
                    secondary_materialized[species_index][group], group, true);
            }
        }
        component_count = component_work.size();
        std::sort(component_work.begin(), component_work.end(),
            [](const ComponentWorkItem& left,
               const ComponentWorkItem& right) {
                if (left.cost != right.cost) return left.cost > right.cost;
                if (left.species_index != right.species_index) {
                    return left.species_index < right.species_index;
                }
                if (left.secondary != right.secondary) {
                    return left.secondary < right.secondary;
                }
                if (left.group_index != right.group_index) {
                    return left.group_index < right.group_index;
                }
                return left.component_index < right.component_index;
            });

        std::atomic<uint64_t> candidate_checks{0};
        std::atomic<uint64_t> sequence_extractions{0};
        std::atomic<uint64_t> direct_ksw_calls{0};
        std::atomic<uint64_t> fallback_ksw_calls{0};
        std::atomic<uint64_t> long_gap_rejections{
            split_statistics.long_gap_rejections};
        std::atomic<uint64_t> maximum_seen_gap{
            split_statistics.maximum_seen_gap};
        std::atomic<uint64_t> estimated_ksw_cells{0};
        std::atomic<size_t> completed_components{0};
        std::atomic<uint64_t> slowest_component_microseconds{0};
        std::mutex progress_mutex;
        auto last_progress = std::chrono::steady_clock::now();
        const auto linking_start = last_progress;
        executeParallelStage(
            component_work.size(), thread_num, extension_activity,
            [&](size_t work_index) {
                const auto& item = component_work[work_index];
                auto& manager = *species_managers[item.species_index];
                AnchorVec& anchors = item.secondary
                    ? secondary_materialized[item.species_index]
                          [item.group_index]
                    : primary_materialized[item.species_index]
                          [item.group_index];
                const auto component_started =
                    std::chrono::steady_clock::now();
                AnchorLinkDetail::Statistics statistics;
                component_outputs[item.output_slot] =
                    AnchorLinkDetail::linkAnchorRange(
                        anchors, item.begin, item.end,
                        *pra.ref_seqpro_manager, manager, &statistics);
                candidate_checks.fetch_add(
                    statistics.candidate_checks, std::memory_order_relaxed);
                sequence_extractions.fetch_add(
                    statistics.sequence_extractions, std::memory_order_relaxed);
                direct_ksw_calls.fetch_add(
                    statistics.direct_ksw_calls, std::memory_order_relaxed);
                fallback_ksw_calls.fetch_add(
                    statistics.fallback_ksw_calls, std::memory_order_relaxed);
                long_gap_rejections.fetch_add(
                    statistics.long_gap_rejections, std::memory_order_relaxed);
                estimated_ksw_cells.fetch_add(
                    statistics.estimated_ksw_cells, std::memory_order_relaxed);
                uint64_t observed = maximum_seen_gap.load(
                    std::memory_order_relaxed);
                while (observed < statistics.maximum_seen_gap &&
                       !maximum_seen_gap.compare_exchange_weak(
                           observed, statistics.maximum_seen_gap,
                           std::memory_order_relaxed)) {
                }
                const uint64_t component_microseconds =
                    static_cast<uint64_t>(std::chrono::duration_cast<
                        std::chrono::microseconds>(
                        std::chrono::steady_clock::now() - component_started)
                        .count());
                observed = slowest_component_microseconds.load(
                    std::memory_order_relaxed);
                while (observed < component_microseconds &&
                       !slowest_component_microseconds.compare_exchange_weak(
                           observed, component_microseconds,
                           std::memory_order_relaxed)) {
                }
                const size_t completed = completed_components.fetch_add(
                    1, std::memory_order_relaxed) + 1;
                if ((completed & 255U) == 0) {
                    std::unique_lock<std::mutex> lock(
                        progress_mutex, std::try_to_lock);
                    const auto now = std::chrono::steady_clock::now();
                    if (lock.owns_lock() &&
                        now - last_progress >= std::chrono::minutes(5)) {
                        last_progress = now;
                        spdlog::info(
                            "[cluster-extension-progress] phase=link "
                            "groups={}/{} components={}/{} "
                            "active={} max_active={} candidate_checks={} "
                            "sequence_extractions={} direct_ksw={} "
                            "fallback_ksw={} long_gap_rejections={} "
                            "maximum_seen_gap={} "
                            "completed_ksw_cells={}",
                            materialized_groups.load(
                                std::memory_order_relaxed),
                            extension_work.size(), completed,
                            component_work.size(),
                            extension_activity.active.load(
                                std::memory_order_relaxed),
                            extension_activity.maximum.load(
                                std::memory_order_relaxed),
                            candidate_checks.load(std::memory_order_relaxed),
                            sequence_extractions.load(
                                std::memory_order_relaxed),
                            direct_ksw_calls.load(std::memory_order_relaxed),
                            fallback_ksw_calls.load(std::memory_order_relaxed),
                            long_gap_rejections.load(
                                std::memory_order_relaxed),
                            maximum_seen_gap.load(std::memory_order_relaxed),
                            estimated_ksw_cells.load(
                                std::memory_order_relaxed));
                        logStageMemory(
                            "cluster-extension", "progress", completed,
                            component_work.size());
                    }
                }
            }, 1);
        const double linking_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - linking_start).count();

        for (size_t species_index = 0;
             species_index < species_list.size(); ++species_index) {
            const auto collect = [&](auto& slots_by_group,
                                     auto& destinations) {
                for (size_t group = 0;
                     group < slots_by_group.size(); ++group) {
                    auto& destination = destinations[group];
                    for (const size_t slot : slots_by_group[group]) {
                        auto& output = component_outputs[slot];
                        destination.insert(destination.end(),
                            std::make_move_iterator(output.begin()),
                            std::make_move_iterator(output.end()));
                    }
                }
            };
            collect(primary_component_slots[species_index],
                    primary_extensions[species_index]);
            collect(secondary_component_slots[species_index],
                    secondary_extensions[species_index]);
        }
        spdlog::info(
            "[cluster-extension-summary] groups={}/{} components={}/{} "
            "active_workers={} max_active_workers={} candidate_checks={} "
            "sequence_extractions={} direct_ksw={} fallback_ksw={} "
            "long_gap_rejections={} maximum_seen_gap={} "
            "estimated_ksw_cells={} completed_ksw_cells={} "
            "materialization_seconds={:.6f} linking_seconds={:.6f} "
            "slowest_component_seconds={:.6f}",
            materialized_groups.load(std::memory_order_relaxed),
            extension_work.size(),
            completed_components.load(std::memory_order_relaxed),
            component_work.size(),
            extension_activity.active.load(std::memory_order_relaxed),
            extension_activity.maximum.load(std::memory_order_relaxed),
            candidate_checks.load(std::memory_order_relaxed),
            sequence_extractions.load(std::memory_order_relaxed),
            direct_ksw_calls.load(std::memory_order_relaxed),
            fallback_ksw_calls.load(std::memory_order_relaxed),
            long_gap_rejections.load(std::memory_order_relaxed),
            maximum_seen_gap.load(std::memory_order_relaxed),
            estimated_ksw_cells.load(std::memory_order_relaxed),
            estimated_ksw_cells.load(std::memory_order_relaxed),
            materialization_seconds, linking_seconds,
            static_cast<double>(slowest_component_microseconds.load(
                std::memory_order_relaxed)) / 1000000.0);
    }
    const double extension_ms =
        std::chrono::duration<
            double,
            std::milli>(
            std::chrono::steady_clock::now() -
            extension_start)
            .count();
    std::vector<std::vector<AnchorVec>>().swap(primary_materialized);
    std::vector<std::vector<AnchorVec>>().swap(secondary_materialized);
    std::vector<ComponentWorkItem>().swap(component_work);
    std::vector<AnchorPtrVec>().swap(component_outputs);
    spdlog::info(
        "[parallel-stage] cluster-extension: "
        "groups={}, components={}, heavy={}, threads={}, max_active={}, "
        "elapsed_ms={:.3f}",
        extension_work.size(),
        component_count,
        heavy_work_items,
        stageWorkerCount(
            thread_num,
            extension_work.size()),
        extension_activity.maximum.load(
            std::memory_order_relaxed),
        extension_ms);
    logStageMemory(
        "cluster-extension", "complete", extension_work.size());
    std::vector<ExtensionWorkItem>().swap(extension_work);

    std::vector<AnchorBySQR_SparsePtr>
        anchor_results(
            species_list.size());
    std::vector<AnchorBySQR_SparsePtr>
        secondary_anchor_results(
            species_list.size());
    for (size_t species_index = 0;
         species_index <
             species_list.size();
         ++species_index) {
        auto anchors =
            std::make_shared<
                AnchorBySQR_Sparse>();
        anchors->reserve(
            primary_extensions[
                species_index]
                .size());
        for (auto& group :
             primary_extensions[
                 species_index]) {
            if (!group.empty()) {
                anchors->emplace_back(
                    std::move(group));
            }
        }
        anchor_results[species_index] =
            std::move(anchors);

        auto secondary_anchors =
            std::make_shared<
                AnchorBySQR_Sparse>();
        secondary_anchors->reserve(
            secondary_extensions[
                species_index]
                .size());
        for (auto& group :
             secondary_extensions[
                 species_index]) {
            if (!group.empty()) {
                secondary_anchors
                    ->emplace_back(
                        std::move(group));
            }
        }
        secondary_anchor_results[
            species_index] =
            std::move(secondary_anchors);
    }
    std::vector<std::vector<AnchorPtrVec>>().swap(primary_extensions);
    std::vector<std::vector<AnchorPtrVec>>().swap(secondary_extensions);
    logStageMemory(
        "cluster-extension", "released", species_list.size());

    if (allow_mem) {
        size_t secondary_anchor_count = 0;
        for (size_t species_index = 0;
             species_index <
                 species_list.size();
             ++species_index) {
            const auto& secondary_anchors =
                secondary_anchor_results[
                    species_index];
            for (const auto& group :
                 *secondary_anchors) {
                secondary_anchor_count +=
                    group.size();
            }
            pra.registerSecondaryAnchors(
                species_list[species_index],
                *seqpro_managers.at(
                    species_list[
                        species_index]),
                secondary_anchors,
                graph,
                is_first);
            for (auto& group : *secondary_anchors) {
                AnchorPtrVec().swap(group);
            }
            secondary_anchor_results[species_index].reset();
        }
        spdlog::info(
            "[secondary-mem] extended_anchors={}",
            secondary_anchor_count);
    }

    const uint_t reference_count =
        std::visit(
            [](auto& manager) {
                return manager
                    ->getSequenceCount();
            },
            *seqpro_managers.at(
                ref_name));
    const auto run_dp_stage =
        [&](bool filter_reference,
            std::string_view stage_name) {
            struct DPGroupRef {
                size_t species_index = 0;
                size_t group_index = 0;
                uint_t chromosome_id = 0;
            };
            struct DPWorkItem {
                size_t species_index = 0;
                uint_t chromosome_id = 0;
                size_t begin = 0;
                size_t end = 0;
                size_t anchor_count = 0;
                size_t original_task_index = 0;
            };
            size_t possible_chromosomes = 0;
            size_t total_groups = 0;
            for (size_t species_index = 0;
                 species_index < species_list.size();
                 ++species_index) {
                const auto& anchors = anchor_results[species_index];
                if (!anchors) continue;
                total_groups += anchors->size();
                if (filter_reference) {
                    possible_chromosomes += reference_count;
                } else {
                    possible_chromosomes += std::visit(
                        [](auto& manager) {
                            return static_cast<size_t>(
                                manager->getSequenceCount());
                        },
                        *seqpro_managers.at(
                            species_list[species_index]));
                }
            }

            std::vector<DPGroupRef> group_refs;
            group_refs.reserve(total_groups);
            for (size_t species_index = 0;
                 species_index < species_list.size();
                 ++species_index) {
                const auto& anchors = anchor_results[species_index];
                if (!anchors) continue;
                for (size_t group_index = 0;
                     group_index < anchors->size();
                     ++group_index) {
                    const auto& group = (*anchors)[group_index];
                    if (group.empty()) continue;
                    group_refs.push_back({
                        species_index,
                        group_index,
                        filter_reference
                            ? group.front()->ref_chr_index
                            : group.front()->qry_chr_index});
                }
            }
            std::sort(
                group_refs.begin(),
                group_refs.end(),
                [](const DPGroupRef& left,
                   const DPGroupRef& right) {
                    if (left.species_index != right.species_index) {
                        return left.species_index < right.species_index;
                    }
                    if (left.chromosome_id != right.chromosome_id) {
                        return left.chromosome_id < right.chromosome_id;
                    }
                    return left.group_index < right.group_index;
                });

            std::vector<DPWorkItem> work_items;
            work_items.reserve(group_refs.size());
            for (size_t index = 0; index < group_refs.size(); ++index) {
                const auto& ref = group_refs[index];
                if (work_items.empty() ||
                    work_items.back().species_index != ref.species_index ||
                    work_items.back().chromosome_id != ref.chromosome_id) {
                    work_items.push_back({
                        ref.species_index,
                        ref.chromosome_id,
                        index,
                        index,
                        0,
                        work_items.size()});
                }
                auto& item = work_items.back();
                item.end = index + 1;
                item.anchor_count +=
                    (*anchor_results[ref.species_index])
                        [ref.group_index]
                            .size();
            }
            std::stable_sort(
                work_items.begin(), work_items.end(),
                [](const DPWorkItem& left, const DPWorkItem& right) {
                    if (left.anchor_count != right.anchor_count) {
                        return left.anchor_count > right.anchor_count;
                    }
                    return left.original_task_index < right.original_task_index;
                });

            OpenMPStageActivity activity;
            logStageMemory(
                stage_name, "start", work_items.size(), group_refs.size());
            const auto stage_start =
                std::chrono::steady_clock::now();
            const uint64_t fallback_before =
                PairRareAligner::dpTreapFallbackCount();
            executeParallelStage(
                work_items.size(),
                thread_num,
                activity,
                [&](size_t work_index) {
                    const auto& item =
                        work_items[
                            work_index];
                    AnchorPtrVec anchors;
                    anchors.reserve(item.anchor_count);
                    for (size_t index = item.begin;
                         index < item.end;
                         ++index) {
                        const auto& ref = group_refs[index];
                        const auto& group =
                            (*anchor_results[ref.species_index])
                                [ref.group_index];
                        anchors.insert(
                            anchors.end(),
                            group.begin(),
                            group.end());
                    }
                    pra.filterAnchorVectorByDP(
                        std::move(anchors),
                        filter_reference);
                });
            const uint64_t fallback_count =
                PairRareAligner::dpTreapFallbackCount() - fallback_before;
            const double elapsed_ms =
                std::chrono::duration<
                    double,
                    std::milli>(
                    std::chrono::steady_clock::now() -
                    stage_start)
                    .count();
            spdlog::info(
                "[parallel-stage] {}: tasks={}, "
                "threads={}, max_active={}, groups={}, "
                "possible_chromosomes={}, skipped_empty={}, "
                "treap_fallbacks={}, "
                "elapsed_ms={:.3f}",
                stage_name,
                work_items.size(),
                stageWorkerCount(
                    thread_num,
                    work_items.size()),
                activity.maximum.load(
                    std::memory_order_relaxed),
                group_refs.size(),
                possible_chromosomes,
                possible_chromosomes >= work_items.size()
                    ? possible_chromosomes - work_items.size()
                    : 0,
                fallback_count,
                elapsed_ms);
            logStageMemory(
                stage_name, "complete", work_items.size(), group_refs.size());
        };
    run_dp_stage(
        true,
        "anchor-dp-reference");
    run_dp_stage(
        false,
        "anchor-dp-query");

    logStageMemory(
        "serial-graph-insertion", "start", species_list.size());
    for (size_t species_index = 0;
         species_index <
             species_list.size();
         ++species_index) {
        const auto& species =
            species_list[species_index];
        auto& anchors =
            anchor_results[species_index];
        if (!anchors) {
            continue;
        }

        spdlog::info(
            "DP filter success for {}",
            species);
        pra.constructGraphByDP(
            species,
            *seqpro_managers.at(species),
            anchors,
            graph);
        anchor_results[species_index].reset();
    }
    logStageMemory(
        "serial-graph-insertion", "complete", species_list.size());

    spdlog::info("[constructMultipleGraphsByDP] Completed all species.");
}
