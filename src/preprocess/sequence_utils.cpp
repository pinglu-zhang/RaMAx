#include "sequence_utils.h"
#include <spdlog/spdlog.h>
#include <algorithm>
#include <limits>

namespace SequenceUtils {

// ------------------------------------------------------------------
// 记录参考序列统计信息（模板版本）
// 说明：
// - 对 reference 物种（species_name == "reference"）额外统计：最短/最长序列长度
// - 对非 reference 物种仅记录序列条数
// - reference_min_seq_length 用于后续限制 sampling_interval（避免采样间隔超过最短序列长度）
// ------------------------------------------------------------------
template<typename ManagerType>
void recordReferenceSequenceStats(const std::string& species_name,
                                 const std::unique_ptr<ManagerType>& manager,
                                 SeqPro::Length& reference_min_seq_length) {
    // 获取序列条数
    auto seq_count = manager->getSequenceCount();

    // 仅对 reference 计算长度统计信息
    if (species_name == "reference") {
        // 获取所有序列名
        auto seq_names = manager->getSequenceNames();

        // 初始化最短/最长长度
        SeqPro::Length min_seq_length = std::numeric_limits<SeqPro::Length>::max();
        SeqPro::Length max_seq_length = 0;

        // 遍历每条序列，统计最短/最长
        for (const auto& seq_name : seq_names) {
            auto seq_length = manager->getSequenceLength(seq_name);
            min_seq_length = std::min(min_seq_length, seq_length);
            max_seq_length = std::max(max_seq_length, seq_length);
        }

        // 打印统计信息
        spdlog::info("[{}] Loaded {} sequences, min length: {}, max length: {}",
                     species_name, seq_count, min_seq_length, max_seq_length);

        // 更新 reference 最短序列长度（用于后续采样间隔截断）
        reference_min_seq_length = min_seq_length;
    } else {
        // 非 reference：只记录序列条数
        spdlog::info("[{}] Loaded {} sequences", species_name, seq_count);
    }
}

// ------------------------------------------------------------------
// 构建 ref_global_cache（基于 ManagerVariant）
// 目的：
// - 为 reference 的“全局坐标 -> 序列ID”查询提供快速近似索引（按 sampling_interval 采样）
// - 避免每次查询都做二分搜索，通过顺序扫描序列区间实现更快的填充
//
// 参数：
// - manager_variant: 参考序列管理器（可能是 SequenceManager 或 MaskedSequenceManager）
// - sampling_interval: 采样间隔（越小越精细但 cache 越大）
// - ref_global_cache: 输出缓存，长度约为 total_length / sampling_interval
//
// 返回值：
// - 构建耗时（秒）
// ------------------------------------------------------------------
double buildRefGlobalCache(const SeqPro::ManagerVariant& manager_variant,
                          SeqPro::Length sampling_interval,
                          sdsl::int_vector<0>& ref_global_cache) {
    // 目前没有清空：理论上上层应保证不会重复构建导致覆盖问题
    auto t_start_cache = std::chrono::steady_clock::now();

    // 获取参考序列总长度，并据此计算 cache 大小
    // - SequenceManager：getTotalLength()
    // - MaskedSequenceManager：getTotalLengthWithSeparators()（包含分隔符等额外长度）
    auto total_length = std::visit([](auto&& manager_ptr) {
        using T = std::decay_t<decltype(manager_ptr)>;
        if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::SequenceManager>>) {
            return manager_ptr->getTotalLength();
        }
        else if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            return manager_ptr->getTotalLengthWithSeparators();
        }
    }, manager_variant);

    // cache_size = total_length / sampling_interval + 1
    auto cache_size = (total_length / sampling_interval) + 1;

    spdlog::info("Building ref_global_cache, sampling_interval={}, cache_size={}",
                 sampling_interval, cache_size);

    // 分配 cache 空间（int_vector 适合后续 bit_compress）
    ref_global_cache.resize(cache_size);

    // 构建 cache：对每个采样点 i，计算 sample_global_pos = i * sampling_interval
    // 并找到其落在哪条序列区间内，记录该序列的 id
    std::visit([&](auto&& manager_ptr) {
        using T = std::decay_t<decltype(manager_ptr)>;

        // 获取所有序列名，并收集对应的 SequenceInfo 指针
        // 目的：按照全局起点排序后进行顺序扫描填充 cache
        auto seq_names = manager_ptr->getSequenceNames();
        std::vector<const SeqPro::SequenceInfo*> seq_infos;
        seq_infos.reserve(seq_names.size());

        for (const auto& name : seq_names) {
            if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::SequenceManager>>) {
                // 普通管理器：直接从索引中取 SequenceInfo
                const auto* info = manager_ptr->getIndex().getSequenceInfo(name);
                if (info) seq_infos.push_back(info);
            } else if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                // 遮蔽管理器：从 original manager 的索引里取 SequenceInfo
                const auto* info = manager_ptr->getOriginalManager().getIndex().getSequenceInfo(name);
                if (info) seq_infos.push_back(info);
            }
        }

        // 按 masked_global_start_pos 排序，保证顺序扫描可行
        std::sort(seq_infos.begin(), seq_infos.end(),
                  [](const SeqPro::SequenceInfo* a, const SeqPro::SequenceInfo* b) {
                      return a->masked_global_start_pos < b->masked_global_start_pos;
                  });

        // 顺序填充：
        // - current_seq_idx 指向当前可能覆盖 sample_global_pos 的序列
        // - 随着采样点递增，current_seq_idx 只会单调递增（避免重复二分搜索）
        size_t current_seq_idx = 0;
        for (SeqPro::Position i = 0; i < cache_size; ++i) {
            SeqPro::Position sample_global_pos = i * sampling_interval;

            // 超过总长度：标记为 INVALID
            if (sample_global_pos >= total_length) {
                ref_global_cache[i] = SeqPro::SequenceIndex::INVALID_ID;
                continue;
            }

            // 查找包含该采样点的序列区间
            while (current_seq_idx < seq_infos.size()) {
                const auto* current_seq = seq_infos[current_seq_idx];

                // 计算当前序列的区间终点
                // 注意：此处对 Masked / 非 Masked 使用同一字段 masked_*（保持原逻辑）
                SeqPro::Position seq_end;
                if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                    seq_end = current_seq->masked_global_start_pos + current_seq->masked_length + 1;
                }
                else {
                    seq_end = current_seq->masked_global_start_pos + current_seq->masked_length + 1;
                }

                // sample_global_pos 落在 [start, end) 内：记录序列 id
                if (sample_global_pos >= current_seq->masked_global_start_pos &&
                    sample_global_pos < seq_end) {
                    ref_global_cache[i] = current_seq->id;
                    break;
                }
                // sample_global_pos 已在该序列之后：推进到下一条序列
                else if (sample_global_pos >= seq_end) {
                    current_seq_idx++;
                }
                // sample_global_pos 在 start 之前：理论上不应发生（序列区间应有序且不回退）
                else {
                    spdlog::warn("Unexpected coordinate order: sample_pos={}, seq_start={}",
                                 sample_global_pos, current_seq->masked_global_start_pos);
                    ref_global_cache[i] = SeqPro::SequenceIndex::INVALID_ID;
                    break;
                }
            }

            // 若所有序列都无法覆盖该采样点：标记为 INVALID
            if (current_seq_idx >= seq_infos.size()) {
                ref_global_cache[i] = SeqPro::SequenceIndex::INVALID_ID;
            }
        }
    }, manager_variant);

    // 对 cache 进行 bit 压缩，节省内存
    sdsl::util::bit_compress(ref_global_cache);

    auto t_end_cache = std::chrono::steady_clock::now();
    std::chrono::duration<double> cache_time = t_end_cache - t_start_cache;

    spdlog::info("ref_global_cache building completed in {:.3f} seconds", cache_time.count());

    // 返回构建耗时（秒）
    return cache_time.count();
}

// ------------------------------------------------------------------
// 记录参考序列统计信息（SharedManagerVariant 版本）
// 说明：
// - SharedManagerVariant 本质上是 shared_ptr<ManagerVariant>
// - 这里通过 std::visit 解包后复用模板版本 recordReferenceSequenceStats
// ------------------------------------------------------------------
void recordReferenceSequenceStats(const std::string& species_name,
                                 const SeqPro::SharedManagerVariant& shared_manager,
                                 SeqPro::Length& reference_min_seq_length) {
    std::visit([&](auto&& manager_ptr) {
        recordReferenceSequenceStats(species_name, manager_ptr, reference_min_seq_length);
    }, *shared_manager);
}

// ------------------------------------------------------------------
// 构建 ref_global_cache（SharedManagerVariant 版本）
// 说明：
// - 直接解引用 shared_ptr，复用 ManagerVariant 版本 buildRefGlobalCache
// ------------------------------------------------------------------
double buildRefGlobalCache(const SeqPro::SharedManagerVariant& shared_manager_variant,
                          SeqPro::Length sampling_interval,
                          sdsl::int_vector<0>& ref_global_cache) {
    // Delegate to the ManagerVariant version by dereferencing the shared_ptr
    return buildRefGlobalCache(*shared_manager_variant, sampling_interval, ref_global_cache);
}

// ------------------------------------------------------------------
// 显式模板实例化（Explicit template instantiations）
// 目的：
// - 在 cpp 中显式实例化常用的模板组合，避免链接阶段找不到符号
// ------------------------------------------------------------------
template void recordReferenceSequenceStats<SeqPro::SequenceManager>(
    const std::string&, const std::unique_ptr<SeqPro::SequenceManager>&, SeqPro::Length&);

template void recordReferenceSequenceStats<SeqPro::MaskedSequenceManager>(
    const std::string&, const std::unique_ptr<SeqPro::MaskedSequenceManager>&, SeqPro::Length&);

} // namespace SequenceUtils
