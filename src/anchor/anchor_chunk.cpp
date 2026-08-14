#include "anchor.h"
#include "SeqPro.h"
#include <algorithm>
#include <spdlog/spdlog.h>

// ------------------------------------------------------------------
// 智能序列分块函数实现
// ------------------------------------------------------------------

RegionVec preAllocateChunks(const SeqPro::ManagerVariant& seq_manager,
                           uint_t chunk_size,
                           uint_t overlap_size,
                           size_t max_sequences_threshold,
                           uint_t min_chunk_size) {
    
    size_t seq_count = std::visit([](auto&& manager_ptr) -> size_t {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
            return manager_ptr->getSequenceCount();
        } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            return manager_ptr->getSequenceCount();
        } else {
            throw std::runtime_error("Unhandled manager type in variant.");
        }
    }, seq_manager);

    uint64_t total_length = std::visit([](auto&& manager_ptr) -> uint64_t {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
            return manager_ptr->getTotalLength();
        } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            return manager_ptr->getTotalLength();
        } else {
            throw std::runtime_error("Unhandled manager type in variant.");
        }
    }, seq_manager);
    
    // 智能策略选择
    bool use_sequence_based_strategy = false;
    
    // 策略1: 如果序列数量很多，优先按序列划分
    if (seq_count >= max_sequences_threshold) {
        use_sequence_based_strategy = true;
        spdlog::info("Using sequence-based chunking strategy: {} sequences (>= threshold {})", 
                    seq_count, max_sequences_threshold);
    }
    // 策略2: 如果序列数量较少但每个序列都很短，也按序列划分
    else if (seq_count > 0 && total_length / seq_count < chunk_size) {
        use_sequence_based_strategy = true;
        spdlog::info("Using sequence-based chunking strategy: average sequence length ({}) < chunk_size ({})", 
                    total_length / seq_count, chunk_size);
    }
    // 策略3: 默认按大小划分
    else {
        spdlog::info("Using size-based chunking strategy: {} sequences, total {} bases", 
                    seq_count, total_length);
    }
    
    if (use_sequence_based_strategy) {
        return preAllocateChunksBySequence(seq_manager);
    } else {
        return preAllocateChunksBySize(seq_manager, chunk_size, overlap_size, min_chunk_size, false);
    }
}

RegionVec preAllocateChunksBySequence(const SeqPro::ManagerVariant& seq_manager) {
    RegionVec chunks;
    auto seq_names = std::visit([](auto&& manager_ptr) -> std::vector<std::string> {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
            return manager_ptr->getSequenceNames();
        }
        else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            return manager_ptr->getSequenceNames();
        }
        else {
            throw std::runtime_error("Unhandled manager type in variant.");
        }
        }, seq_manager);

    chunks.reserve(seq_names.size());

    for (const auto& seq_name : seq_names) {
        uint64_t seq_length = std::visit([&seq_name](auto&& manager_ptr) -> uint64_t {
            using PtrType = std::decay_t<decltype(manager_ptr)>;
            if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
                return manager_ptr->getSequenceLength(seq_name);
            }
            else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                return manager_ptr->getSequenceLength(seq_name);
            }
            else {
                throw std::runtime_error("Unhandled manager type in variant.");
            }
            }, seq_manager);

        auto seq_id = std::visit([&seq_name](auto&& manager_ptr) -> SeqPro::SequenceId {
            return manager_ptr->getSequenceId(seq_name);
            }, seq_manager);

        chunks.push_back({ seq_id, 0, static_cast<Coord_t>(seq_length) });
    }

    spdlog::info("Generated {} chunks by sequence (one chunk per sequence)", chunks.size());
    return chunks;
}

RegionVec preAllocateChunksBySize(const SeqPro::ManagerVariant& seq_manager,
                                 uint_t chunk_size,
                                 uint_t overlap_size,
                                 uint_t min_chunk_size,
                                 bool isMultiple) {
    RegionVec chunks;
    
    // 预估总chunk数以减少realloc
    uint64_t total_length = std::visit([](auto&& manager_ptr) -> uint64_t {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
            return manager_ptr->getTotalLength();
        } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            return manager_ptr->getTotalLengthWithSeparators();
        } else {
            throw std::runtime_error("Unhandled manager type in variant.");
        }
    }, seq_manager);
    size_t est_chunks = (total_length + chunk_size - 1) / chunk_size;
    chunks.reserve(est_chunks);
    
    auto seq_names = std::visit([](auto&& manager_ptr) -> std::vector<std::string> {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
            return manager_ptr->getSequenceNames();
        } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            return manager_ptr->getSequenceNames();
        } else {
            throw std::runtime_error("Unhandled manager type in variant.");
        }
    }, seq_manager);
    
    // 对每条序列分别切分
    for (const auto& seq_name : seq_names) {
        uint64_t seq_length = std::visit([&seq_name](auto&& manager_ptr) -> uint64_t {
            using PtrType = std::decay_t<decltype(manager_ptr)>;
            if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
                return manager_ptr->getSequenceLength(seq_name);
            } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                return manager_ptr->getSequenceLengthWithSeparators(seq_name);
            } else {
                throw std::runtime_error("Unhandled manager type in variant.");
            }
        }, seq_manager);

        // 如果是多基因组对齐且是MaskedSequenceManager，使用遮蔽区间预分割
        if (isMultiple) {
            std::visit([&](auto&& manager_ptr) {
                using PtrType = std::decay_t<decltype(manager_ptr)>;
                if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                    splitByMaskedRegions(chunks, seq_name, seq_length, *manager_ptr, chunk_size, overlap_size, seq_manager);
                } else {
                    // 对于非MaskedSequenceManager，使用常规分割
                    normalSizeBasedChunking(chunks, seq_name, seq_length, chunk_size, overlap_size, min_chunk_size, seq_manager);
                }
            }, seq_manager);
        } else {
            // 常规分割逻辑
            normalSizeBasedChunking(chunks, seq_name, seq_length, chunk_size, overlap_size, min_chunk_size, seq_manager);
        }
    }
    
    spdlog::info("Generated {} chunks by size (chunk_size: {}, overlap_size: {}, isMultiple: {})", 
                chunks.size(), chunk_size, overlap_size, isMultiple);
    return chunks;
}

// 辅助函数：常规的基于大小的分块
void normalSizeBasedChunking(RegionVec& chunks, const std::string& seq_name,
    uint64_t seq_length, uint_t chunk_size,
    uint_t overlap_size, uint_t min_chunk_size, const SeqPro::ManagerVariant& seq_manager) {

    auto seq_id = std::visit([&seq_name](auto&& manager_ptr) -> SeqPro::SequenceId {
        return manager_ptr->getSequenceId(seq_name);
        }, seq_manager);


    // 如果序列长度小于等于chunk_size或min_chunk_size，则只生成一个不重叠chunk
    if (seq_length <= chunk_size || seq_length <= min_chunk_size) {
        chunks.push_back({ seq_id, 0, static_cast<Coord_t>(seq_length) });
        return;
    }

    // 多个chunk，需要在它们之间保留overlap_size
    uint64_t start = 0;
    while (start < seq_length) {
        // 本chunk的实际长度（最后一块可能不足chunk_size）
        uint64_t this_len = std::min<uint64_t>(chunk_size, seq_length - start);
        chunks.push_back({ seq_id, static_cast<Coord_t>(start), static_cast<Coord_t>(this_len) });

        // 计算下一个chunk的起始位置：前进chunk_size，然后回退overlap_size
        if (start + this_len >= seq_length) {
            break;
        }
        start += chunk_size;
        // 保证不回退超过已经走过的距离
        start = (start >= overlap_size ? start - overlap_size : 0);
    }
}

// 辅助函数：基于遮蔽区间的预分割
void splitByMaskedRegions(RegionVec& chunks, const std::string& seq_name, 
                        uint64_t seq_length, 
                        const SeqPro::MaskedSequenceManager& masked_manager,
                        uint_t chunk_size, uint_t overlap_size, const SeqPro::ManagerVariant& seq_manager) {
    try {
        // 获取该序列的遮蔽区间
        auto masked_intervals = masked_manager.getMaskIntervals(seq_name);
        
        if (masked_intervals.empty()) {
            //spdlog::debug("No masked intervals for sequence {}, using normal chunking", seq_name);
            normalSizeBasedChunking(chunks, seq_name, seq_length, chunk_size, overlap_size, 0, seq_manager);
            return;
        }
        
        // spdlog::info("Processing {} masked intervals for sequence {}", masked_intervals.size(), seq_name);
        
        // // 按起始位置排序遮蔽区间
        // std::sort(masked_intervals.begin(), masked_intervals.end(),
        //          [](const SeqPro::MaskInterval& a, const SeqPro::MaskInterval& b) {
        //              return a.start < b.start;
        //          });
        
        // // 合并重叠的遮蔽区间
        // std::vector<SeqPro::MaskInterval> merged_intervals;
        // for (const auto& interval : masked_intervals) {
        //     if (merged_intervals.empty() || merged_intervals.back().end < interval.start) {
        //         merged_intervals.push_back(interval);
        //     } else {
        //         merged_intervals.back().end = std::max(merged_intervals.back().end, interval.end);
        //     }
        // }
        
        // 在非遮蔽区间中进行分块
        uint64_t current_pos = 0;
        for (const auto& masked_interval : masked_intervals) {
            uint64_t unmasked_start = current_pos;
            uint64_t unmasked_end = masked_interval.start;
            
            // 处理遮蔽区间前的非遮蔽区域
            if (unmasked_start < unmasked_end) {
                chunkUnmaskedRegion(chunks, seq_name, unmasked_start, unmasked_end, chunk_size, overlap_size, seq_manager);
            }
            
            // 跳过遮蔽区间
            current_pos = masked_interval.end;
        }
        
        // 处理最后一个遮蔽区间之后的非遮蔽区域
        if (current_pos < seq_length) {
            chunkUnmaskedRegion(chunks, seq_name, current_pos, seq_length, chunk_size, overlap_size, seq_manager);
        }
        
        // spdlog::info("Generated {} chunks for sequence {} with masked region pre-splitting",
        //              chunks.size(), seq_name);
                     
    } catch (const std::exception& e) {
        spdlog::warn("Failed to get masked intervals for sequence {}: {}, using normal chunking", 
                    seq_name, e.what());
        normalSizeBasedChunking(chunks, seq_name, seq_length, chunk_size, overlap_size, 0, seq_manager);
    }
}

// 辅助函数：对非遮蔽区域进行分块
void chunkUnmaskedRegion(RegionVec& chunks, const std::string& seq_name,
    uint64_t region_start, uint64_t region_end,
    uint_t chunk_size, uint_t overlap_size, const SeqPro::ManagerVariant& seq_manager) {
    uint64_t region_length = region_end - region_start;

    auto seq_id = std::visit([&seq_name](auto&& manager_ptr) -> SeqPro::SequenceId {
        return manager_ptr->getSequenceId(seq_name);
        }, seq_manager);

    // 如果区域太小，直接作为一个chunk
    if (region_length <= chunk_size) {
        chunks.push_back({ seq_id, static_cast<Coord_t>(region_start), static_cast<Coord_t>(region_length) });
        return;
    }

    // 在该区域内进行分块
    uint64_t start = region_start;
    while (start < region_end) {
        uint64_t this_len = std::min<uint64_t>(chunk_size, region_end - start);
        chunks.push_back({ seq_id, static_cast<Coord_t>(start), static_cast<Coord_t>(this_len) });

        // 计算下一个chunk的起始位置
        if (start + this_len >= region_end) {
            break;
        }
        start += chunk_size;
        // 保证不回退超过已经走过的距离，同时不超出区域边界
        start = (start >= overlap_size ? start - overlap_size : region_start);
    }
}
