#include "anchor.h"

#include <SeqPro.h>

#include "data_process.h"

////--------------------------------------------------------------------
//// 适配 4‑D 结构的聚类分发
////--------------------------------------------------------------------
//void groupMatchByQueryRef(MatchVec3DPtr& anchors,
//    MatchByStrandByQueryRefPtr unique_anchors,
//    MatchByStrandByQueryRefPtr repeat_anchors,
//    SeqPro::ManagerVariant& ref_fasta_manager,
//    SeqPro::ManagerVariant& query_fasta_manager,
//	ThreadPool& pool)
//{
//    //------------------------------------------------------------
//    // 0) 初始化输出矩阵 [strand][query][ref]
//    //------------------------------------------------------------
//    constexpr uint_t STRAND_CNT = 2;                         // 0=FWD,1=REV
//    const uint_t ref_chr_cnt = std::visit([](auto&& manager) {
//        using PtrType = std::decay_t<decltype(manager)>;
//        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
//            return manager->getSequenceCount();
//        } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
//            return manager->getSequenceCount();
//        } else {
//            throw std::runtime_error("Unhandled manager type in variant.");
//        }
//    }, ref_fasta_manager);
//
//    const uint_t query_chr_cnt = std::visit([](auto&& manager) {
//        using PtrType = std::decay_t<decltype(manager)>;
//        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
//            return manager->getSequenceCount();
//        } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
//            return manager->getSequenceCount();
//        } else {
//            throw std::runtime_error("Unhandled manager type in variant.");
//        }
//    }, query_fasta_manager);
//    auto initTarget = [&](MatchByStrandByQueryRefPtr& tgt) {
//        tgt->resize(STRAND_CNT);
//        for (uint_t s = 0; s < STRAND_CNT; ++s) {
//            (*tgt)[s].resize(query_chr_cnt);
//            for (uint_t q = 0; q < query_chr_cnt; ++q)
//                (*tgt)[s][q].resize(ref_chr_cnt);
//        }
//        };
//    initTarget(unique_anchors);
//    initTarget(repeat_anchors);
//
//    //------------------------------------------------------------
//    // 1) 行级互斥锁：放到堆上，用 shared_ptr 延长生命周期
//    //------------------------------------------------------------
//    auto rowLocks = std::make_shared<std::vector<std::mutex>>(STRAND_CNT * query_chr_cnt);
//
//    //------------------------------------------------------------
//    // 2) 并行遍历 3-D 数据，slice 级别开任务
//    //------------------------------------------------------------
//    for (auto& slice : *anchors) {
//        // 按值捕获 slice、rowLocks、unique_anchors 等
//        pool.enqueue(
//            [slice,                                               // 值
//            rowLocks,                                            // 值 (shared_ptr)
//            unique_anchors, repeat_anchors,                      // 值 (shared_ptr)
//            &ref_fasta_manager, &query_fasta_manager,            // 引用
//            query_chr_cnt] () mutable
//            {
//                if (slice.empty()) return;
//
//                const Match& first = slice.front().front();
//                uint_t sIdx = (first.strand == REVERSE ? 1u : 0u);
//
//            uint_t qIdx = std::visit([&first](auto&& manager) {
//                using PtrType = std::decay_t<decltype(manager)>;
//                if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
//                    return manager->getSequenceId(first.query_region.chr_name);
//                } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
//                    return manager->getSequenceId(first.query_region.chr_name);
//                } else {
//                    throw std::runtime_error("Unhandled manager type in variant.");
//                }
//            }, query_fasta_manager);
//
//            if(qIdx == SeqPro::SequenceIndex::INVALID_ID) return;
//
//            for (auto& vec : slice)
//            {
//                if (vec.empty()) continue;
//
//                uint_t rIdx = std::visit([&vec](auto&& manager) {
//                    using PtrType = std::decay_t<decltype(manager)>;
//                    if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
//                        return manager->getSequenceId(vec.front().ref_region.chr_name);
//                    } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
//                        return manager->getSequenceId(vec.front().ref_region.chr_name);
//                    } else {
//                        throw std::runtime_error("Unhandled manager type in variant.");
//                    }
//                }, ref_fasta_manager);
//
//                if(rIdx == SeqPro::SequenceIndex::INVALID_ID) continue;
//
//                // 计算锁下标并加锁
//                uint_t lockIdx = sIdx * query_chr_cnt + qIdx;
//                std::lock_guard<std::mutex> lk((*rowLocks)[lockIdx]);
//
//                MatchVec& tgt = (vec.size() == 1)
//                    ? (*unique_anchors)[sIdx][qIdx][rIdx]
//                    : (*repeat_anchors)[sIdx][qIdx][rIdx];
//
//                    tgt.insert(tgt.end(),
//                        std::make_move_iterator(vec.begin()),
//                        std::make_move_iterator(vec.end()));
//                    vec.clear();
//                    vec.shrink_to_fit();
//                }
//            });
//    }
//    
//}
// ---------------------------------------------------------------------------
// groupMatchByQueryRef – SERIAL, LOW‑MEMORY VERSION
// ---------------------------------------------------------------------------
// - 完全取消并行，不依赖 ThreadPool，也不需要锁。
// - 按需扩展 3‑D 目标结构，避免一次性为空槽分配 vector。
// - 每次把 slice 内部 vec 的元素 move 到目标后，立即 clear()+shrink_to_fit()
//   以释放多余 capacity，降低峰值 RSS。
// ---------------------------------------------------------------------------
// 依赖类型说明（保持与原工程一致）：
//   using Match               = ...;
//   using MatchVec            = std::vector<Match>;
//   using MatchVec2D          = std::vector<MatchVec>;                // [ref]
//   using MatchVec3D          = std::vector<MatchVec2D>;              // [query][ref]
//   using MatchVec3DPtr       = std::shared_ptr<MatchVec3D>;          // slice 别名
//   using MatchByStrandByQueryRef = std::vector<MatchVec3D>;          // [strand][query][ref]
//   using MatchByStrandByQueryRefPtr = std::shared_ptr<MatchByStrandByQueryRef>;
//   enum Strand { FORWARD, REVERSE };
// ---------------------------------------------------------------------------

void groupMatchByQueryRef(
    MatchVec3DPtr& anchors,
    MatchByStrandByQueryRefPtr unique_anchors,
    MatchByStrandByQueryRefPtr repeat_anchors,
    SeqPro::ManagerVariant& ref_fasta_manager,
    SeqPro::ManagerVariant& query_fasta_manager,
    ThreadPool & pool)
{
    constexpr uint_t STRAND_CNT = 2; // 0 = FWD, 1 = REV

    // ----------------- 0) 获取染色体数 -----------------
    const uint_t ref_chr_cnt = std::visit([](auto& m) { return m->getSequenceCount(); }, ref_fasta_manager);
    const uint_t qry_chr_cnt = std::visit([](auto& m) { return m->getSequenceCount(); }, query_fasta_manager);

    // ----------------- 1) 按需扩展工具 lambda -----------------
    auto ensure_slot = [&](MatchByStrandByQueryRefPtr& tgt,
        uint_t s, uint_t q, uint_t r) -> MatchVec&
        {
            //if (s >= tgt->size())            tgt->resize(STRAND_CNT);
            //if (q >= (*tgt)[s].size())       (*tgt)[s].resize(qry_chr_cnt);
            //if (r >= (*tgt)[s][q].size())    (*tgt)[s][q].resize(ref_chr_cnt);
            return (*tgt)[s][q][r];
        };

    auto initTarget = [&](MatchByStrandByQueryRefPtr& tgt) {
    tgt->resize(STRAND_CNT);
    for (uint_t s = 0; s < STRAND_CNT; ++s) {
        (*tgt)[s].resize(qry_chr_cnt);
        for (uint_t q = 0; q < qry_chr_cnt; ++q)
            (*tgt)[s][q].resize(ref_chr_cnt);
    }
    };
    initTarget(unique_anchors);
    initTarget(repeat_anchors);

    // ----------------- 2) 顺序遍历所有 slice -----------------
    for (auto& slice : *anchors) {
        if (slice.empty()) continue;

        // sIdx & qIdx 只需要从 slice 第一个元素即可得出
        const Match& first = slice.front().front();
        const uint_t sIdx = (first.strand() == REVERSE ? 1u : 0u);

        ChrIndex qIdx = first.qry_chr_index;
        if (qIdx == SeqPro::SequenceIndex::INVALID_ID) continue;

        // --- 遍历 slice 内的每个 MatchVec（同一 query & strand, 不同 ref）
        for (auto& vec : slice) {
            if (vec.empty()) continue;

            ChrIndex rIdx = vec.front().ref_chr_index;
            if (rIdx == SeqPro::SequenceIndex::INVALID_ID) continue;

            MatchVec& dest = (vec.size() == 1)
                ? ensure_slot(unique_anchors, sIdx, qIdx, rIdx)
                : ensure_slot(repeat_anchors, sIdx, qIdx, rIdx);

            if (dest.empty()) dest.reserve(vec.size()); // 减少后续扩容

            dest.insert(dest.end(),
                std::make_move_iterator(vec.begin()),
                std::make_move_iterator(vec.end()));

            // --- 立即回收 vec 占用容量 ---
            vec.clear();
            vec.shrink_to_fit();
        }
    }

    // ----------------- 3) anchors 自身可以释放 -----------------
    anchors->clear();
    anchors->shrink_to_fit();

    // 遍历unique_anchors，使用shrink_to_fit()
    auto shrink_matrix = [](MatchByStrandByQueryRefPtr& matrix) {
        for (auto& strand_v : *matrix) {
            for (auto& qry_v : strand_v) {
                for (auto& ref_v : qry_v) {
                    if (ref_v.capacity() > ref_v.size())
                        ref_v.shrink_to_fit();
                }
            }
        }
        };
    shrink_matrix(unique_anchors);
    shrink_matrix(repeat_anchors);
}



// 重载版本：支持 SharedManagerVariant，通过解引用调用原始版本
void groupMatchByQueryRef(MatchVec3DPtr& anchors,
    MatchByStrandByQueryRefPtr unique_anchors,
    MatchByStrandByQueryRefPtr repeat_anchors,
    SeqPro::SharedManagerVariant& ref_fasta_manager,
    SeqPro::SharedManagerVariant& query_fasta_manager,
    ThreadPool& pool)
{
    // 解引用 SharedManagerVariant 并调用原始函数
    groupMatchByQueryRef(anchors, unique_anchors, repeat_anchors,
                        *ref_fasta_manager, *query_fasta_manager, pool);
}

//void sortMatchByRefStart(MatchByQueryRefPtr& anchors, ThreadPool& pool) {
//
//    for (auto& row : *anchors) {
//        for (auto& vec : row) {
//            pool.enqueue([&vec]() {
//                if (vec.size() == 0) return;
//                std::sort(vec.begin(), vec.end(), [](const Match& a, const Match& b) {
//                    return (a.ref_region.start < b.ref_region.start) || (a.ref_region.start == b.ref_region.start && a.query_region.start < b.query_region.start);
//                    });
//                });
//        }
//    }
//
//    pool.waitAllTasksDone();
//}

// anchors[strand][query][ref] -- MatchByStrandByQueryRef
void sortMatchByQueryStart(MatchByStrandByQueryRefPtr& anchors, ThreadPool& pool)
{
    constexpr size_t MIN_SIZE_FOR_PARALLEL = 100;  // 阈值：小于100个元素直接串行排序
    constexpr size_t BATCH_SIZE = 50;              // 批处理大小
    
    std::vector<std::future<void>> futures;
    std::vector<std::vector<MatchVec*>> batches;
    std::vector<MatchVec*> current_batch;
    
    // 收集所有需要排序的向量，并按大小分类
    for (auto& strand_layer : *anchors) {
        for (auto& query_row : strand_layer) {
            for (auto& vec : query_row) {
                if (vec.empty()) continue;
                
                if (vec.size() >= MIN_SIZE_FOR_PARALLEL) {
                    // 大向量：单独提交任务
                    futures.emplace_back(pool.enqueue([v = &vec]() {
                        std::sort(v->begin(), v->end(),
                            [](const Match& a, const Match& b) {
                                return a.qry_start < b.qry_start;
                            });
                    }));
                } else {
                    // 小向量：加入批处理
                    current_batch.push_back(&vec);
                    if (current_batch.size() >= BATCH_SIZE) {
                        batches.push_back(std::move(current_batch));
                        current_batch.clear();
                    }
                }
            }
        }
    }
    
    // 处理剩余的小向量
    if (!current_batch.empty()) {
        batches.push_back(std::move(current_batch));
    }
    
    // 批处理小向量
    for (auto& batch : batches) {
        futures.emplace_back(pool.enqueue([batch = std::move(batch)]() {
            for (auto* vec : batch) {
                std::sort(vec->begin(), vec->end(),
                    [](const Match& a, const Match& b) {
                        return a.qry_start < b.qry_start;
                    });
            }
        }));
    }
    
    // 等待所有任务完成
    for (auto& future : futures) {
        future.get();
    }
}

//
//AnchorPtrVec findNonOverlapAnchors(const AnchorVec& anchors)
//{
//    AnchorPtrVec uniques;
//    const size_t n = anchors.size();
//    if (n == 0) return uniques;
//
//    Coord_t max_end_so_far = 0;  // 记录到当前位置前的最大 end
//
//    for (size_t i = 0; i < n; ++i) {
//        const Anchor& curr = anchors[i];
//        Coord_t curr_start = curr.match.ref_region.start;
//        Coord_t curr_end = curr_start + curr.match.ref_region.length;
//
//        bool overlap_left = (i > 0 && max_end_so_far > curr_start);
//        bool overlap_right = (i + 1 < n &&
//            anchors[i + 1].match.ref_region.start < curr_end);
//
//        if (!overlap_left && !overlap_right)
//            uniques.emplace_back(std::make_shared<Anchor>(curr));
//
//        // 更新左侧最大 end
//        if (curr_end > max_end_so_far) max_end_so_far = curr_end;
//    }
//    return std::move(uniques);
//}

MatchClusterVecPtr
groupClustersToVec(const ClusterVecPtrByStrandByQueryRefPtr& src,
    ThreadPool& pool, uint_t thread_num)
{
    auto result = std::make_shared<MatchClusterVec>();
    
    if (!src || src->empty()) {
        return result;
    }
    
    // 预计算总的簇数量，用于预分配空间
    size_t total_clusters = 0;
    for (const auto& strand_layer : *src) {
        for (const auto& query_row : strand_layer) {
            for (const auto& cluster_vec_ptr : query_row) {
                if (cluster_vec_ptr && !cluster_vec_ptr->empty()) {
                    total_clusters += cluster_vec_ptr->size();
                }
            }
        }
    }
    
    result->reserve(total_clusters);
    
    // 使用线程安全的方式收集所有簇
    std::mutex result_mutex;
    std::vector<std::future<std::vector<MatchCluster>>> futures;
    
    // 将收集任务分配给多个线程
    size_t tasks_per_thread = std::max(size_t(1), total_clusters / thread_num);
    size_t current_task_count = 0;
    std::vector<MatchCluster> current_batch;
    
    for (const auto& strand_layer : *src) {
        for (const auto& query_row : strand_layer) {
            for (const auto& cluster_vec_ptr : query_row) {
                if (!cluster_vec_ptr || cluster_vec_ptr->empty()) continue;
                
                for (const auto& cluster : *cluster_vec_ptr) {
                    if (!cluster.empty()) {
                        current_batch.push_back(cluster);
                        current_task_count++;
                        
                        // 当批次足够大时，提交到线程池
                        if (current_task_count >= tasks_per_thread) {
                            auto batch_copy = current_batch;  // 复制用于线程安全
                            futures.emplace_back(pool.enqueue([batch_copy]() -> std::vector<MatchCluster> {
                                return batch_copy;
                            }));
                            
                            current_batch.clear();
                            current_task_count = 0;
                        }
                    }
                }
            }
        }
    }
    
    // 处理剩余的批次
    if (!current_batch.empty()) {
        futures.emplace_back(pool.enqueue([current_batch]() -> std::vector<MatchCluster> {
            return current_batch;
        }));
    }
    
    // 收集所有结果
    for (auto& future : futures) {
        auto batch_result = future.get();
        for (auto& cluster : batch_result) {
            result->emplace_back(std::move(cluster));
        }
    }
    
    return result;
}

// ------------------------------------------------------------------
// 把 3D: [strand][queryRef][ref] 的聚簇重新组织成 1D: [ref]
// ------------------------------------------------------------------
ClusterVecPtrByRefPtr groupClustersToRefVec(
    const ClusterVecPtrByStrandByQueryRefPtr& src,
    ThreadPool& pool, uint_t thread_num)
{
    // ---------- 1. 空输入 ----------
    if (!src || src->empty())
        return std::make_shared<ClusterVecPtrByRef>();

    // ---------- 2. 探测 ref 数 ----------
    size_t n_ref = 0;
    for (const auto& strandVec : *src) {
        for (const auto& queryVec : strandVec) {
            if (!queryVec.empty()) {
                n_ref = queryVec.size();
                break;
            }
        }
        if (n_ref) break;
    }
    if (n_ref == 0)
        return std::make_shared<ClusterVecPtrByRef>();

    // ---------- 3. 目标 ----------
    auto dst = std::make_shared<ClusterVecPtrByRef>(n_ref, nullptr);

    /* ******************************************************************
     * 4. 决定是否并行
     *    只能使用 pool 已有的信息，不能改 ThreadPool 接口。
     *    - 大部分线程池都会暴露  thread_count()/size()/parallelism() 之类的读接口；
     *      如果你的实现没有，可以在 ThreadPool 里加一个 const 方法，
     *      但这属于“实现细节”而不是“函数接口”变更，外部签名不动。
     ******************************************************************/
    const bool can_parallel =
        n_ref > 1 &&   // 工作量足够大才值得并行
        thread_num > 1;   // 剩余线程数至少 1 条

    /* =================== 顺序分支 =================== */
    if (!can_parallel) {
        for (size_t ref_id = 0; ref_id < n_ref; ++ref_id) {
            auto combined = std::make_shared<MatchClusterVec>();
            for (const auto& strandVec : *src)
                for (const auto& queryVec : strandVec)
                    if (ref_id < queryVec.size() && queryVec[ref_id])
                        combined->insert(combined->end(),
                            queryVec[ref_id]->begin(),
                            queryVec[ref_id]->end());
            (*dst)[ref_id] = std::move(combined);
        }
        return dst;                   // ✅ 单线程直接结束
    }

    /* =================== 并行分支（原逻辑） =================== */
    std::vector<std::future<void>> futures;
    futures.reserve(n_ref);

    for (size_t ref_id = 0; ref_id < n_ref; ++ref_id) {
        futures.emplace_back(
            pool.enqueue([&, ref_id] {
                auto combined = std::make_shared<MatchClusterVec>();
                for (const auto& strandVec : *src)
                    for (const auto& queryVec : strandVec)
                        if (ref_id < queryVec.size() && queryVec[ref_id])
                            combined->insert(combined->end(),
                                queryVec[ref_id]->begin(),
                                queryVec[ref_id]->end());
                (*dst)[ref_id] = std::move(combined);
                }));
    }
    for (auto& fut : futures) fut.get();
    return dst;
}

// ------------------------------------------------------------------
// 把 3D: [strand][queryRef][ref] 的聚簇重新组织成 1D: [ref]
// ------------------------------------------------------------------
ClusterVecPtrByRefPtr
groupClustersToRefVec(const ClusterVecPtrByStrandByQueryRefPtr& src,
    ThreadPool& pool)
{
    // ---------- 1. 空输入快速返回 ----------
    if (!src || src->empty()) {
        return std::make_shared<ClusterVecPtrByRef>();
    }

    // ---------- 2. 探测参考染色体 (Ref) 数 ----------
    size_t n_ref = 0;
    for (const auto& strandVec : *src) {
        for (const auto& queryVec : strandVec) {
            if (!queryVec.empty()) {
                n_ref = queryVec.size();   // 第一次遇到非空即足够
                break;
            }
        }
        if (n_ref) break;
    }
    if (n_ref == 0) {
        return std::make_shared<ClusterVecPtrByRef>();   // 没有 Ref
    }

    // ---------- 3. 准备目标结构 ----------
    auto dst = std::make_shared<ClusterVecPtrByRef>(n_ref, nullptr);

    // ---------- 4. 为 “每个 Ref ⇨ 一条任务” ----------
    std::vector<std::future<void>> futures;
    futures.reserve(n_ref);

    for (size_t ref_id = 0; ref_id < n_ref; ++ref_id) {
        futures.emplace_back(
            pool.enqueue([&, ref_id] {
                auto combined = std::make_shared<MatchClusterVec>();

                // 聚合：遍历所有 [strand][query]，把 ref_id 处的聚簇并入
                for (const auto& strandVec : *src) {
                    for (const auto& queryVec : strandVec) {
                        if (ref_id < queryVec.size() && queryVec[ref_id]) {
                            combined->insert(combined->end(),
                                queryVec[ref_id]->begin(),
                                queryVec[ref_id]->end());
                        }
                    }
                }
                // 无竞态：独占写入各自槽位
                (*dst)[ref_id] = std::move(combined);
                })
        );
    }

    // ---------- 5. 等待全部任务完成 ----------
    for (auto& fut : futures) fut.get();

    return dst;
}




//--------------------------------------------------------------------
// Anchor 验证功能实现
//--------------------------------------------------------------------

ValidationResult validateAnchorsCorrectness(
    const MatchVec3DPtr& anchors,
    const SeqPro::ManagerVariant& ref_manager,
    const SeqPro::ManagerVariant& query_manager
) {
    spdlog::info("开始验证 anchors 结果的正确性…");
    
    ValidationResult result;
    
    // 反向互补函数
    auto reverseComplement = [](const std::string &seq) -> std::string {
        std::string result;
        result.reserve(seq.length());
        for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
            result.push_back(BASE_COMPLEMENT[static_cast<unsigned char>(*it)]);
        }
        return result;
    };

    // 计算总工作项数
    uint64_t total_work_items = 0;
    for (const auto &level1: *anchors) {
        for (const auto &level2: level1) {
            total_work_items += level2.size();
        }
    }

    if (total_work_items == 0) {
        spdlog::info("没有需要验证的匹配项。");
        return result;
    }

    spdlog::info("总计需要验证 {} 个匹配项。", total_work_items);

    // 进度报告相关变量
    std::atomic<uint64_t> processed_items = 0;
    const uint64_t report_interval = std::max(1ULL, total_work_items / 100ULL);
    std::atomic<bool> diagnostic_dump_done = false;

    // 用于OpenMP reduction的临时变量
    uint64_t total_matches = 0;
    uint64_t correct_matches = 0;
    uint64_t incorrect_matches = 0;

#pragma omp parallel reduction(+: total_matches, correct_matches, incorrect_matches)
    {
#pragma omp for schedule(dynamic) nowait
        for (size_t i = 0; i < anchors->size(); ++i) {
            for (size_t j = 0; j < (*anchors)[i].size(); ++j) {
                for (size_t k = 0; k < (*anchors)[i][j].size(); ++k) {
                    uint64_t current_processed = processed_items.fetch_add(1, std::memory_order_relaxed) + 1;

                    if (current_processed % report_interval == 0) {
                        spdlog::info("验证进度: {} / {} ({:.2f}%)",
                                     current_processed,
                                     total_work_items,
                                     (100.0 * current_processed / total_work_items));
                    }

                    const Match &match = (*anchors)[i][j][k];
                    ++total_matches;

                    try {
                        // 提取参考序列
                        std::string ref_seq = std::visit([&match](auto &&manager_ptr) -> std::string {
                            using PtrType = std::decay_t<decltype(manager_ptr)>;
                            if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
                                return manager_ptr->getSubSequence(match.ref_chr_index, match.ref_start,
                                                                   match.match_len());
                            } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<
                                SeqPro::MaskedSequenceManager> >) {
                                return manager_ptr->getSubSequence(match.ref_chr_index, match.ref_start,
                                    match.match_len());
                            } else {
                                throw std::runtime_error("Unhandled manager type in variant.");
                            }
                        }, ref_manager);
                        
                        // 提取查询序列
                        std::string query_seq = std::visit([&match](auto &&manager_ptr) -> std::string {
                            using PtrType = std::decay_t<decltype(manager_ptr)>;
                            if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
                                return manager_ptr->getSubSequence(match.qry_chr_index,
                                                                   match.qry_start, match.match_len());
                            } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<
                                SeqPro::MaskedSequenceManager> >) {
                                return manager_ptr->getSubSequence(match.qry_chr_index,
                                    match.qry_start, match.match_len());
                            } else {
                                throw std::runtime_error("Unhandled manager type in variant.");
                            }
                        }, query_manager);

                        // 如果是反向链，进行反向互补
                        if (match.strand() == REVERSE) {
                            query_seq = reverseComplement(query_seq);
                        }

                        // 比较序列
                        if (ref_seq == query_seq) {
                            ++correct_matches;
                        } else {
                            ++incorrect_matches;

                            // 线程安全的首错捕获逻辑
                            bool expected = false;
                            if (!diagnostic_dump_done.load(std::memory_order_relaxed) &&
                                diagnostic_dump_done.compare_exchange_strong(expected, true)) {
                                spdlog::warn(
                                    "序列不匹配: ref_chr={}, ref_start={}, query_chr={}, "
                                    "query_start={}, length={}, strand={}\n"
                                    "  Ref Seq:    {}\n"
                                    "  Query Seq{}: {}",
                                    match.ref_chr_index, match.ref_start, match.qry_chr_index,
                                    match.qry_start, match.match_len(),
                                    (match.strand() == FORWARD ? "FORWARD" : "REVERSE"), ref_seq,
                                    (match.strand() == REVERSE ? " (RC)" : ""), query_seq
                                );
                                spdlog::error("--- [CAPTURED FIRST MISMATCH] INITIATING DIAGNOSTIC DUMP ---");

                                // 打印错误匹配的详细信息
                                spdlog::error("Failing Match Details:");
                                spdlog::error("  - Reference: {}:{} (len:{})", match.ref_chr_index,
                                              match.ref_start, match.match_len());
                                spdlog::error("  - Query:     {}:{} (len:{})", match.qry_chr_index,
                                              match.qry_start, match.match_len());
                                spdlog::error("  - Strand:    {}", (match.strand() == FORWARD ? "FORWARD" : "REVERSE"));
                            }
                        }
                    } catch (const std::exception &e) {
                        ++incorrect_matches;
                        spdlog::warn("处理匹配项时发生异常: {}", e.what());
                    }
                }
            }
        } // omp for
    } // omp parallel

    // 将OpenMP reduction结果赋值给结果结构体
    result.total_matches = total_matches;
    result.correct_matches = correct_matches;
    result.incorrect_matches = incorrect_matches;

    // 确保最终进度是100%
    spdlog::info("验证进度: {} / {} (100.00%)", total_work_items, total_work_items);

    spdlog::info("验证完成: 总匹配数={}, 正确匹配数={}, 错误匹配数={}, 正确率={:.2f}%",
                 result.total_matches,
                 result.correct_matches,
                 result.incorrect_matches,
                 result.accuracy());

    return result;
}

// 重载版本：支持 SharedManagerVariant
ValidationResult validateAnchorsCorrectness(
    const MatchVec3DPtr& anchors,
    const SeqPro::SharedManagerVariant& ref_manager,
    const SeqPro::SharedManagerVariant& query_manager
) {
    // 解引用 SharedManagerVariant 并调用原始函数
    return validateAnchorsCorrectness(anchors, *ref_manager, *query_manager);
}


AnchorVec extendClusterToAnchorVec(const MatchCluster& cluster,
    const SeqPro::ManagerVariant& ref_mgr,
    const SeqPro::ManagerVariant& query_mgr)
{
    AnchorVec anchors;
    if (cluster.empty()) return anchors;

    // todo 为了修复Ramax的BUG作了修改，需要加RamaG的判定
    // -- 快速 slice 提取：visit 一次，避免重复 λ 创建 --
    auto subSeq = [&](const SeqPro::ManagerVariant& mv,
        const ChrIndex& chr, Coord_t b, Coord_t l) -> std::string {
    return std::visit([&](auto& p) {
        using T = std::decay_t<decltype(p)>;
        if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::SequenceManager>>) {
            return p->getSubSequence(chr, b, l);
        } else if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            return p->getOriginalManager().getSubSequence(chr, b, l);
        }
    }, mv);
};

    /* ===== 初始 anchor 状态 ===== */
    const Match& first = cluster.front();
    Strand strand = first.strand();
    bool   fwd         = (strand == FORWARD); 
    ChrIndex ref_chr = first.ref_chr_index;
    ChrIndex qry_chr = first.qry_chr_index;

    Coord_t ref_beg = start1(first);
    Coord_t ref_end = ref_beg;

    Coord_t qry_beg = 0;
    Coord_t qry_end = 0;
	if (fwd) {
		qry_beg = start2(first);
	}
	else {
		qry_beg = start2(first) + len2(first);
	}
    qry_end = qry_beg;

    Cigar_t cigar; cigar.reserve(cluster.size() * 2);  // 预估
    Coord_t aln_len = 0;
    Coord_t match_len = 0;
    auto pushEq = [&](uint32_t len) {
        appendCigarOp(cigar, 'M', len);
        aln_len += len;  
        match_len += len;
        ref_end += len;  
        if (fwd) {
            qry_end += len;
        }
        else {
            qry_end -= len;
        }

        };
    auto flush = [&] {
        if (fwd) {
            anchors.emplace_back(ref_chr, ref_beg, ref_end - ref_beg, qry_chr, qry_beg, qry_end - qry_beg, strand, aln_len, 0, std::move(cigar));
        }
        else {
            anchors.emplace_back(ref_chr, ref_beg, ref_end - ref_beg, qry_chr, qry_end, qry_beg - qry_end, strand, aln_len, 0, std::move(cigar));
        }
        cigar.clear(); cigar.shrink_to_fit(); cigar.reserve(16);
        aln_len = 0;
        match_len = 0;
        ref_beg = ref_end;          // 推进到下一段起点
        qry_beg = qry_end;          // 同理（fwd 递增，rev 递减）
        };


    //wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    //attributes.distance_metric = gap_affine;
    //attributes.affine_penalties.mismatch = 2;      // X > 0
    //attributes.affine_penalties.gap_opening = 3;   // O >= 0
    //attributes.affine_penalties.gap_extension = 1; // E > 0
    //attributes.memory_mode = wavefront_memory_ultralow;
    //attributes.heuristic.strategy = wf_heuristic_wfadaptive;
    //attributes.heuristic.min_wavefront_length = 10;
    //attributes.heuristic.max_distance_threshold = 50;
    //attributes.heuristic.steps_between_cutoffs = 1;

    //wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);

    /* ==================== 遍历 cluster ==================== */
    for (size_t i = 0;i < cluster.size();++i) {
        const Match& m = cluster[i];
        pushEq(len1(m));                                  // 精确 match

        if (i + 1 == cluster.size()) break;

        const Match& nxt = cluster[i + 1];
        Coord_t rgBeg = start1(m) + len1(m), rgEnd = start1(nxt);
        Coord_t qgBeg = 0;
        Coord_t qgEnd = 0;
        if (fwd) {
            qgBeg = start2(m) + len2(m);
            qgEnd = start2(nxt);
        } else {
            qgBeg = start2(nxt) + len2(nxt);
            qgEnd = start2(m);
        } 
        
        uint32_t rgLen = rgEnd > rgBeg ? rgEnd - rgBeg : 0;
        uint32_t qgLen = qgEnd > qgBeg ? qgEnd - qgBeg : 0;
        // 无 gap
        if (rgLen == 0 && qgLen == 0) {
            continue;                       
        }
        Cigar_t buf;
        Cigar_t gap = {};
        if (rgLen == 0 || qgLen == 0) {
            // 纯 I / 纯 D，不跑比对
            char op = (rgLen == 0 ? 'I' : 'D');
            uint32_t len = (rgLen == 0 ? qgLen : rgLen);
            gap.push_back(cigarToInt(op, len));
        }
        else {
            // 2) 获取 gap 片段
            std::string ref_gap = subSeq(ref_mgr, ref_chr, rgBeg, rgEnd - rgBeg);
            std::string qry_gap = subSeq(query_mgr, qry_chr, qgBeg, qgEnd - qgBeg);
            if (strand == REVERSE) reverseComplement(qry_gap);

            uint_t Lt = ref_gap.size();
            uint_t Lq = qry_gap.size();
            int64_t d = std::abs(static_cast<int64_t>(Lt) - static_cast<int64_t>(Lq));

            //double rho = double(d) / std::min(Lt, Lq);

            //
            //if (rho <= 0.3 && Lt > 10 && Lq > 10) {
            //    int cigar_len;
            //    uint32_t* cigar_tmp;
            //    wavefront_align(wf_aligner, ref_gap.c_str(), ref_gap.length(), qry_gap.c_str(), qry_gap.length());
            //    cigar_get_CIGAR(wf_aligner->cigar, false, &cigar_tmp, &cigar_len);
            //    for (uint_t j = 0; j < cigar_len; ++j) {
            //        gap.push_back(cigar_tmp[j]);
            //    }

            //}
            //else {
            //    gap = globalAlignKSW2(ref_gap, qry_gap);
            //}
            gap = globalAlignKSW2(ref_gap, qry_gap);

            /* ---- 扫描 gap-cigar，遇 >50bp I/D 即分段 ---- */
            buf.reserve(gap.size());
        }

        for (auto unit : gap) {
            uint32_t len = unit >> 4;
            uint8_t  op = unit & 0xf;          // 0=M,1=I,2=D,7='=',8='X'
            bool big = ((op == 1 || op == 2) && len > 50);

            if (big) {
                // 先把已有片段 merge
                if (!buf.empty()) { appendCigar(cigar, buf); buf.clear(); }
                flush();                        // 输出 anchor
                // 移动起点：I 影响 query，D 影响 ref
                if (op == 1) {
                    if (fwd) {
                        qry_beg += len;
                    }
                    else {
                        qry_beg -= len;
                    }
                }
                else {
                    ref_beg += len;
                }       
                ref_end = ref_beg; qry_end = qry_beg;
            }
            else {
                buf.push_back(unit);
                // 更新末端坐标
                if (op == 1) {
                    if (fwd) {
                        qry_end += len;
                    }
                    else {
                        qry_end -= len;
                    }
                }            
                else if (op == 2)       ref_end += len;
                else { 
                    ref_end += len; 
                    if (fwd) {
                        qry_end += len;
                    }
                    else {
                        qry_end -= len;
                    }
                }
                if (op != 3) aln_len += len;     // 3(N)不会出现
            }
        }
        if (!buf.empty()) appendCigar(cigar, buf);
    }
    flush();   
    //wavefront_aligner_delete(wf_aligner);// 收尾
    return anchors;
}

std::pair<std::string, std::string> renderAlignment(
    const std::string& ref, const std::string& qry, const Cigar_t& cigar)
{
    std::ostringstream ref_line, qry_line;
    size_t rpos = 0, qpos = 0;
    for (auto c : cigar) {
        char op; uint32_t len;
        intToCigar(c, op, len);

        switch (op) {
        case 'M': case '=': case 'X':
            for (uint32_t i = 0; i < len; ++i) {
                char r = (rpos < ref.size()) ? ref[rpos++] : '-';
                char q = (qpos < qry.size()) ? qry[qpos++] : '-';
                ref_line << r;
                qry_line << q;
            }
            break;
        case 'I':
            for (uint32_t i = 0; i < len; ++i) {
                ref_line << '-';
                qry_line << ((qpos < qry.size()) ? qry[qpos++] : '-');
            }
            break;
        case 'D':
            for (uint32_t i = 0; i < len; ++i) {
                ref_line << ((rpos < ref.size()) ? ref[rpos++] : '-');
                qry_line << '-';
            }
            break;
        default:
            break;
        }
    }
    return { ref_line.str(), qry_line.str() };
}

Anchor extendClusterToAnchor(MatchCluster& cluster,
    const SeqPro::ManagerVariant& ref_mgr,
    const SeqPro::ManagerVariant& query_mgr) {
    if (cluster.empty()) return Anchor();
    Anchor anchor;
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
    const Match& first = cluster.front();
    if (first.qry_start == 212285) {
        std::cout << "";
    }
    Strand strand = first.strand();
    bool   fwd = (strand == FORWARD);

    if (!fwd) {
        std::sort(cluster.begin(), cluster.end(),
            [](auto& a, auto& b) { return a.ref_start < b.ref_start; });
    }

    ChrIndex ref_chr = first.ref_chr_index;
    ChrIndex qry_chr = first.qry_chr_index;

    Cigar_t cigar; cigar.reserve(cluster.size() * 2);  // 预估
    Coord_t aln_len = 0;
    Coord_t match_len = 0;
    for (size_t i = 0;i < cluster.size();++i) {
        const Match& m = cluster[i];

		uint_t len = len1(m); 
        uint_t ref_start = m.ref_start + len;
        uint_t qry_start = m.qry_start + len;

        appendCigarOp(cigar, 'M', len);
        aln_len += len;
        match_len += len;

        if (i + 1 == cluster.size()) break;

        const Match& nxt = cluster[i + 1];
        uint_t len2 = len1(nxt);
		uint_t ref_end = nxt.ref_start; 
        uint_t qry_end = nxt.qry_start;
        std::string ref_gap;
        std::string qry_gap;
        if (fwd) {
            ref_gap = subSeq(ref_mgr, ref_chr, ref_start, ref_end - ref_start);
            qry_gap = subSeq(query_mgr, qry_chr, qry_start, qry_end - qry_start);
        }
        else {
			ref_gap = subSeq(ref_mgr, ref_chr, ref_start, ref_end - ref_start);
			qry_gap = subSeq(query_mgr, qry_chr, qry_end + len2, m.qry_start - len2 - qry_end);
			reverseComplement(qry_gap);
        }


        Cigar_t gap = globalAlignKSW2(ref_gap, qry_gap);
        appendCigar(cigar, gap);
        match_len += countMatchOperations(gap);
        aln_len += countAlignmentLength(gap);

  //      //如果identity小于70% ，打印出具体比对结果
  //      uint_t min_len = std::min(ref_end - ref_start, qry_end - qry_start);
		//double identity = (min_len > 0) ? static_cast<double>(countMatchOperations(gap)) / min_len : 0.0;
  //      if (first.qry_start == 212285) {
  //          spdlog::warn("Low identity gap detected: Ref[{}:{}-{}], Qry[{}:{}-{}], Identity: {:.2f}%",
  //              ref_chr, ref_start, ref_end,
  //              qry_chr, qry_start, qry_end,
  //              identity * 100.0);
  //          auto [ref_aln, qry_aln] = renderAlignment(ref_gap, qry_gap, gap);

  //          spdlog::info("Ref: {}", ref_aln);
  //          spdlog::info("Qry: {}", qry_aln);
  //          std::cout << "";
  //      }
    }
	const Match last = cluster.back();
    if (fwd) {
        anchor = Anchor(ref_chr, first.ref_start, last.ref_start + last.match_len() - first.ref_start, qry_chr, first.qry_start, last.qry_start + last.match_len() - first.qry_start, strand, aln_len, match_len, std::move(cigar));
    }
    else {
        anchor = Anchor(ref_chr, first.ref_start, last.ref_start + last.match_len() - first.ref_start, qry_chr, last.qry_start, first.qry_start + first.match_len() - last.qry_start, strand, aln_len, match_len, std::move(cigar));
    }

    //if (first.qry_start == 212285) {
    //    // 提取 anchor 对应的序列
    //    std::string ref_seq = subSeq(ref_mgr, ref_chr, anchor.ref_start, anchor.ref_len);
    //    std::string qry_seq = subSeq(query_mgr, qry_chr, anchor.qry_start, anchor.qry_len);

    //    // 渲染对齐结果
    //    auto [ref_aln,  qry_aln] = renderAlignment(ref_seq, qry_seq, anchor.cigar);

    //    spdlog::info("Ref: {}", ref_aln);
    //    spdlog::info("Qry: {}", qry_aln);
    //    spdlog::info("CIGAR: {}", cigarToString(anchor.cigar));
    //    std::cout << "";
    //}

    

    return anchor;
}
//
//Anchor extendClusterToAnchor(MatchCluster& cluster,
//    const SeqPro::ManagerVariant& ref_mgr,
//    const SeqPro::ManagerVariant& query_mgr)
//{
//    if (cluster.empty()) return Anchor();
//    AnchorVec anchors;
//    // todo 为了修复Ramax的BUG作了修改，需要加RamaG的判定
//    // -- 快速 slice 提取：visit 一次，避免重复 λ 创建 --
//    auto subSeq = [&](const SeqPro::ManagerVariant& mv,
//        const ChrIndex& chr, Coord_t b, Coord_t l) -> std::string {
//            return std::visit([&](auto& p) {
//                using T = std::decay_t<decltype(p)>;
//                if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::SequenceManager>>) {
//                    return p->getSubSequence(chr, b, l);
//                }
//                else if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
//                    return p->getOriginalManager().getSubSequence(chr, b, l);
//                }
//                }, mv);
//        };
//
//    /* ===== 初始 anchor 状态 ===== */
//    const Match& first = cluster.front();
//    Strand strand = first.strand();
//    bool   fwd = (strand == FORWARD);
//
//    if (!fwd) {
//        std::sort(cluster.begin(), cluster.end(),
//            [](auto& a, auto& b) { return a.ref_start < b.ref_start; });
//    }
//    ChrIndex ref_chr = first.ref_chr_index;
//    ChrIndex qry_chr = first.qry_chr_index;
//
//    Coord_t ref_beg = start1(first);
//    Coord_t ref_end = ref_beg;
//
//    Coord_t qry_beg = 0;
//    Coord_t qry_end = 0;
//    if (fwd) {
//        qry_beg = start2(first);
//    }
//    else {
//        qry_beg = start2(first) + len2(first);
//    }
//    qry_end = qry_beg;
//
//    Cigar_t cigar; cigar.reserve(cluster.size() * 2);  // 预估
//    Coord_t aln_len = 0;
//    Coord_t match_len = 0;
//
//    auto pushEq = [&](uint32_t len) {
//        appendCigarOp(cigar, 'M', len);
//        aln_len += len;
//        match_len += len;
//        ref_end += len;
//        if (fwd) {
//            qry_end += len;
//        }
//        else {
//            qry_end -= len;
//        }
//
//        };
//    auto flush = [&] {
//        if (fwd) {
//            anchors.emplace_back(ref_chr, ref_beg, ref_end - ref_beg, qry_chr, qry_beg, qry_end - qry_beg, strand, aln_len, match_len, std::move(cigar));
//        }
//        else {
//            anchors.emplace_back(ref_chr, ref_beg, ref_end - ref_beg, qry_chr, qry_end, qry_beg - qry_end, strand, aln_len, match_len, std::move(cigar));
//        }
//        cigar.clear(); cigar.shrink_to_fit(); cigar.reserve(16);
//        aln_len = 0;
//        match_len = 0;
//        ref_beg = ref_end;          // 推进到下一段起点
//        qry_beg = qry_end;          // 同理（fwd 递增，rev 递减）
//        };
//
//
//    //wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
//    //attributes.distance_metric = gap_affine;
//    //attributes.affine_penalties.mismatch = 2;      // X > 0
//    //attributes.affine_penalties.gap_opening = 3;   // O >= 0
//    //attributes.affine_penalties.gap_extension = 1; // E > 0
//    //attributes.memory_mode = wavefront_memory_ultralow;
//    //attributes.heuristic.strategy = wf_heuristic_wfadaptive;
//    //attributes.heuristic.min_wavefront_length = 10;
//    //attributes.heuristic.max_distance_threshold = 50;
//    //attributes.heuristic.steps_between_cutoffs = 1;
//
//    //wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
//
//    /* ==================== 遍历 cluster ==================== */
//    for (size_t i = 0;i < cluster.size();++i) {
//        const Match& m = cluster[i];
//        pushEq(len1(m));                                  // 精确 match
//
//        if (i + 1 == cluster.size()) break;
//
//        const Match& nxt = cluster[i + 1];
//        Coord_t rgBeg = start1(m) + len1(m), rgEnd = start1(nxt);
//        Coord_t qgBeg = 0;
//        Coord_t qgEnd = 0;
//        if (fwd) {
//            qgBeg = start2(m) + len2(m);
//            qgEnd = start2(nxt);
//        }
//        else {
//            qgBeg = start2(nxt) + len2(nxt);
//            qgEnd = start2(m);
//        }
//
//        uint32_t rgLen = rgEnd > rgBeg ? rgEnd - rgBeg : 0;
//        uint32_t qgLen = qgEnd > qgBeg ? qgEnd - qgBeg : 0;
//        // 无 gap
//        if (rgLen == 0 && qgLen == 0) {
//            continue;
//        }
//        Cigar_t buf;
//        Cigar_t gap = {};
//        if (rgLen == 0 || qgLen == 0) {
//            // 纯 I / 纯 D，不跑比对
//            char op = (rgLen == 0 ? 'I' : 'D');
//            uint32_t len = (rgLen == 0 ? qgLen : rgLen);
//            gap.push_back(cigarToInt(op, len));
//        }
//        else {
//            // 2) 获取 gap 片段
//            std::string ref_gap = subSeq(ref_mgr, ref_chr, rgBeg, rgEnd - rgBeg);
//            std::string qry_gap = subSeq(query_mgr, qry_chr, qgBeg, qgEnd - qgBeg);
//            if (strand == REVERSE) reverseComplement(qry_gap);
//
//            uint_t Lt = ref_gap.size();
//            uint_t Lq = qry_gap.size();
//            int64_t d = std::abs(static_cast<int64_t>(Lt) - static_cast<int64_t>(Lq));
//
//            double rho = double(d) / std::min(Lt, Lq);
//
//
//            //if (rho <= 0.3 && Lt > 10 && Lq > 10) {
//            //    int cigar_len;
//            //    uint32_t* cigar_tmp;
//            //    wavefront_align(wf_aligner, ref_gap.c_str(), ref_gap.length(), qry_gap.c_str(), qry_gap.length());
//            //    cigar_get_CIGAR(wf_aligner->cigar, false, &cigar_tmp, &cigar_len);
//            //    for (uint_t j = 0; j < cigar_len; ++j) {
//            //        gap.push_back(cigar_tmp[j]);
//            //    }
//
//            //}
//            //else {
//            //    gap = globalAlignKSW2_2(ref_gap, qry_gap);
//            //}
//
//            gap = globalAlignKSW2(ref_gap, qry_gap);
//
//            /* ---- 扫描 gap-cigar，遇 >50bp I/D 即分段 ---- */
//            buf.reserve(gap.size());
//        }
//
//        for (auto unit : gap) {
//            uint32_t len = unit >> 4;
//            uint8_t  op = unit & 0xf;          // 0=M,1=I,2=D,7='=',8='X'
//            //bool big = ((op == 1 || op == 2) && len > 50);
//
//            if (false) {
//                // 先把已有片段 merge
//                if (!buf.empty()) { appendCigar(cigar, buf); buf.clear(); }
//                flush();                        // 输出 anchor
//                // 移动起点：I 影响 query，D 影响 ref
//                if (op == 1) {
//                    if (fwd) {
//                        qry_beg += len;
//                    }
//                    else {
//                        qry_beg -= len;
//                    }
//                }
//                else {
//                    ref_beg += len;
//                }
//                ref_end = ref_beg; qry_end = qry_beg;
//            }
//            else {
//                buf.push_back(unit);
//                // 更新末端坐标
//                if (op == 1) {
//                    if (fwd) {
//                        qry_end += len;
//                    }
//                    else {
//                        qry_end -= len;
//                    }
//                }
//                else if (op == 2)       ref_end += len;
//                else {
//                    ref_end += len;
//                    if (fwd) {
//                        qry_end += len;
//                    }
//                    else {
//                        qry_end -= len;
//                    }
//                }
//                if (op != 3) aln_len += len;     // 3(N)不会出现
//                if (op == 0 || op == 7) match_len += len;
//            }
//        }
//        if (!buf.empty()) appendCigar(cigar, buf);
//    }
//    flush();
//    //wavefront_aligner_delete(wf_aligner);// 收尾
//    return anchors[0];
//}

/// 同时验证 ref/query，两者都通过才算成功
void validateClusters(const ClusterVecPtrByStrandByQueryRefPtr& cluster_vec_ptr)
{
    if (!cluster_vec_ptr) {
        spdlog::debug("validateClusters: cluster_vec_ptr is null, nothing to check.");
        return;
    }

    std::size_t total_clusters = 0;
    std::size_t failed_clusters = 0;

    bool reverse_cluster = false;

    for (std::size_t strand_i = 0; strand_i < cluster_vec_ptr->size(); ++strand_i) {
        const auto& by_query = (*cluster_vec_ptr)[strand_i];

        for (std::size_t q_i = 0; q_i < by_query.size(); ++q_i) {
            const auto& by_ref = by_query[q_i];

            for (std::size_t r_i = 0; r_i < by_ref.size(); ++r_i) {
                const auto& clusters_ptr = by_ref[r_i];
                if (!clusters_ptr) continue;

                for (std::size_t c_i = 0; c_i < clusters_ptr->size(); ++c_i) {
                    ++total_clusters;
                    const MatchCluster& cluster = (*clusters_ptr)[c_i];
 
                    Strand strand = cluster[0].strand();

					if (strand == REVERSE && cluster.size() > 5) {
						reverse_cluster = true;
					}   

                    for (std::size_t m_i = 0; m_i < cluster.size(); ++m_i) {
                        if (cluster[m_i].strand() != strand) {
                            spdlog::debug(
                                "❌ Cluster FAILED (strand={},query={},ref={},cluster={}): "
                                "strand mismatch (expected {}, got {})",
                                strand_i, q_i, r_i, c_i, static_cast<int>(strand), static_cast<int>(cluster[m_i].strand()));
                        }
                    }


                    Coord_t ref_last_end = std::numeric_limits<Coord_t>::min();
                    Coord_t qry_last_pos;                 // 启动值依赖方向
					if (strand == FORWARD) {
						qry_last_pos = std::numeric_limits<Coord_t>::min();
					}
					else { // REVERSE
						qry_last_pos = std::numeric_limits<Coord_t>::max();
					}
                    

                    for (std::size_t m_i = 0; m_i < cluster.size(); ++m_i) {
                        const Match& m = cluster[m_i];

                        //-------------------//
                        // 1) 参考坐标检查
                        //-------------------//
                        Coord_t ref_start = m.ref_start;
                        Coord_t ref_end = ref_start + m.match_len(); // 右开

                        if (ref_start < ref_last_end) {
                            spdlog::debug(
                                "❌ Cluster FAILED (strand={},query={},ref={},cluster={}): "
                                "ref_start < ref_last_end ({} < {})",
                                strand_i, q_i, r_i, c_i, ref_start, ref_last_end);
                        }
                            
                        ref_last_end =  ref_end;

                        //-------------------//
                        // 2) 查询坐标检查
                        //-------------------//
                        Coord_t qry_start = m.qry_start;
                        Coord_t qry_end = qry_start + m.match_len();


                        
                        if (m.strand() == FORWARD) {
                            // 正向：要求升序且不重叠
                            if (qry_start < qry_last_pos) {
                                spdlog::debug(
                                    "❌ Cluster FAILED (strand={},query={},ref={},cluster={}): "
                                    "qry_start < qry_last_pos ({} < {})",
                                    strand_i, q_i, r_i, c_i, qry_start, qry_last_pos);
                            }
                            qry_last_pos = qry_end;
                        }
                        else { // REVERSE
                            // 反向：查询坐标应当递减，区间不重叠
                            if (qry_end > qry_last_pos) {
                                spdlog::debug(
                                    "❌ Cluster FAILED (strand={},query={},ref={},cluster={}): "
                                    "qry_start > qry_last_pos ({} > {})",
                                    strand_i, q_i, r_i, c_i, qry_start, qry_last_pos);
                            }
                            qry_last_pos = qry_start;
                        }
                        
                    } // end matches loop


                }
            }
        }
    }
    if (reverse_cluster == false) {
		spdlog::debug("reverse chain may fail");
		return;
    }

    spdlog::debug("validateClusters finished: {} clusters checked, {} failed.",
        total_clusters, failed_clusters);
}

void linkClusters(AnchorPtrVec& anchors,
    const SeqPro::ManagerVariant& ref_mgr,
    const SeqPro::ManagerVariant& qry_mgr)
{
    AnchorPtrVec linked;
    if (anchors.empty()) return;
	Strand strand = anchors[0]->strand;
    if (strand == FORWARD) {
        std::sort(anchors.begin(), anchors.end(),
            [](auto& a, auto& b) { return a->qry_start < b->qry_start; });
    }
    else {
        std::sort(anchors.begin(), anchors.end(),
            [](auto& a, auto& b) { return a->qry_start > b->qry_start; });
    }
    //std::sort(anchors.begin(), anchors.end(),
    //    [](auto& a, auto& b) { return a->ref_start < b->ref_start; });


    const int K = 2000; // 最多向前看 50 个 anchor

    auto curr = anchors.begin();

    while (true) {
        if (curr == anchors.end()) {
            break;
        }
        if ((*curr)->is_linked) {
            ++curr;
            /*it = curr + 1;*/
            continue;
        }


        //if ((*curr)->qry_chr_index == 4 && (*curr)->qry_start < 357618 + 24610 && (*curr)->qry_start + (*curr)->qry_len > 357618) {
        //    std::cout << "";
        //}


        auto best = anchors.end();

        int_t break_len = 200;

		int_t best_score = std::numeric_limits<int_t>::max();

        int looked = 0;

        for (auto it = curr + 1; it != anchors.end() && looked < K; ++it) {
            if ((*it)->is_linked) continue;
            ++looked;

            // ---- gap 计算 ----
            int_t ref_gap = static_cast<int_t>((*it)->ref_start)
                - static_cast<int_t>((*curr)->ref_start + (*curr)->ref_len);

            int_t qry_gap = 0;
            if ((*curr)->strand == FORWARD) {
                qry_gap = static_cast<int_t>((*it)->qry_start)
                    - static_cast<int_t>((*curr)->qry_start + (*curr)->qry_len);
            }
            else { // 反向链
                qry_gap = static_cast<int_t>((*curr)->qry_start)
                    - static_cast<int_t>((*it)->qry_start + (*it)->qry_len);
            }


            if (ref_gap < 0 || qry_gap < 0) {
                //if (ref_gap < 0 && qry_gap < 0) {
                //    (*it)->is_linked = true;
                //}
                continue; // overlap
            }
            long greater = std::max(ref_gap, qry_gap);
            long lesser = std::min(ref_gap, qry_gap);

            if (greater < break_len || greater - lesser <= break_len) {
                best = it;
                break;
            }
            int_t this_score = (greater << 1) - lesser;
            if (best_score > this_score) {
				best_score = this_score;
                best = it;
            }
        }

        if (best == anchors.end()) {
            ++curr;
            continue;
        }

        // ========== 提取 gap 序列 ==========
        Coord_t ref_gap_beg = (*curr)->ref_start + (*curr)->ref_len;
        Coord_t ref_gap_len = (*best)->ref_start - ref_gap_beg;
        //ref_gap_len = std::min(ref_gap_len, (Coord_t)10000);

        Coord_t qry_gap_beg, qry_gap_len;
        if ((*curr)->strand == FORWARD) {
            qry_gap_beg = (*curr)->qry_start + (*curr)->qry_len;
            qry_gap_len = (*best)->qry_start - qry_gap_beg;
        }
        else { // REVERSE
            qry_gap_beg = (*best)->qry_start + (*best)->qry_len;
            qry_gap_len = (*curr)->qry_start - qry_gap_beg;
        }
        //qry_gap_len = std::min(qry_gap_len, (Coord_t)10000);

   //     if (ref_gap_len > 30000 || qry_gap_len > 30000) {
			//(*curr)->is_linked = true;
   //         linked.push_back(*curr);
   //         curr++;
   //         continue;
   //     }

        std::string ref_gap_seq = std::visit([&](auto& p) {
            using T = std::decay_t<decltype(p)>;
            if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::SequenceManager>>) {
                return p->getSubSequence((*curr)->ref_chr_index, ref_gap_beg, ref_gap_len);
            }
            else {
                return p->getOriginalManager().getSubSequence((*curr)->ref_chr_index, ref_gap_beg, ref_gap_len);
            }
            }, ref_mgr);

        std::string qry_gap_seq = std::visit([&](auto& p) {
            using T = std::decay_t<decltype(p)>;
            if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::SequenceManager>>) {
                return p->getSubSequence((*curr)->qry_chr_index, qry_gap_beg, qry_gap_len);
            }
            else {
                return p->getOriginalManager().getSubSequence((*curr)->qry_chr_index, qry_gap_beg, qry_gap_len);
            }
            }, qry_mgr);

        if ((*curr)->strand == REVERSE) {
            reverseComplement(qry_gap_seq); // 保证方向一致
        }
        bool reach = false;
        uint_t ref_len = 0;
        uint_t qry_len = 0;
        Cigar_t gap_cigar;
        if (qry_gap_len > 10000 || ref_gap_len > 10000) {
            reach = false;
        }
        else {
            //gap_cigar = globalAlignKSW2_2(ref_gap_seq, qry_gap_seq);
            gap_cigar = extendAlignKSW2(ref_gap_seq, qry_gap_seq, 2 * break_len);
            //gap_cigar = extendAlignWFA2(ref_gap_seq, qry_gap_seq, break_len);
            ref_len = countRefLength(gap_cigar);
            qry_len = countQryLength(gap_cigar);
        }
        // ========== 调用 WFA 进行 gap 比对 ==========

        //if (checkGapCigarQuality(gap_cigar, ref_gap_len, qry_gap_len, 0.6)){
        if (ref_len == ref_gap_len && qry_len == qry_gap_len) {
			reach = true;
            // ---- 更新 curr 坐标 ----
            (*curr)->ref_len = ((*best)->ref_start + (*best)->ref_len) - (*curr)->ref_start;

            if ((*curr)->strand == FORWARD) {
                (*curr)->qry_len = ((*best)->qry_start + (*best)->qry_len) - (*curr)->qry_start;
                appendCigar((*curr)->cigar, gap_cigar);
                appendCigar((*curr)->cigar, (*best)->cigar);
            }
            else {
                (*curr)->qry_len = (*curr)->qry_start + (*curr)->qry_len - (*best)->qry_start;
				(*curr)->qry_start = (*best)->qry_start; // 反向链，更新起点
                Cigar_t c1 = gap_cigar;
                //std::reverse(c1.begin(), c1.end());
                Cigar_t c2 = (*best)->cigar;
                Cigar_t c3 = (*curr)->cigar;
                //std::reverse(c2.begin(), c2.end());
                (*curr)->cigar = c3;

                appendCigar((*curr)->cigar, c1);
                appendCigar((*curr)->cigar, c2);
            }

            (*curr)->alignment_length += (*best)->alignment_length + countAlignmentLength(gap_cigar);
            (*curr)->aligned_base += countMatchOperations(gap_cigar) + (*best)->aligned_base;

            (*best)->is_linked = true;
		}
		else {
			reach = false;
		}

        if(reach == false) {
            // ========== 在 push_back 前，尝试与 linked.back() 的 gap 比对 ==========
            if (!linked.empty()) {
                const int LOOK_BACK = 2000;
                const int_t MAX_GAP = 100000;

                auto best_it = linked.rend();  // 初始化为无效
                int_t best_score = std::numeric_limits<int_t>::max();
                int checked = 0;
                bool found_best = false;

                // 🔁 倒序遍历最近的 linked anchors
                for (auto it = linked.rbegin(); it != linked.rend() && checked < LOOK_BACK; ++it, ++checked) {
                    auto& prev = *it;

                    // strand、染色体一致才考虑
                    if (prev->strand != (*curr)->strand ||
                        prev->ref_chr_index != (*curr)->ref_chr_index ||
                        prev->qry_chr_index != (*curr)->qry_chr_index)
                        continue;

                    // ---- 计算 gap ----
                    int_t ref_gap = static_cast<int_t>((*curr)->ref_start)
                        - static_cast<int_t>(prev->ref_start + prev->ref_len);

                    int_t qry_gap;
                    if ((*curr)->strand == FORWARD) {
                        qry_gap = static_cast<int_t>((*curr)->qry_start)
                            - static_cast<int_t>(prev->qry_start + prev->qry_len);
                    }
                    else { // REVERSE
                        qry_gap = static_cast<int_t>(prev->qry_start)
                            - static_cast<int_t>((*curr)->qry_start + (*curr)->qry_len);
                    }

                    if (ref_gap < 0 || qry_gap < 0 || ref_gap > MAX_GAP || qry_gap > MAX_GAP)
                        continue;

                    long greater = std::max(ref_gap, qry_gap);
                    long lesser = std::min(ref_gap, qry_gap);
                    int_t this_score = (greater << 1) - lesser;

                    if (this_score < best_score) {
                        best_score = this_score;
                        best_it = it;
                        found_best = true;
                    }
                }

// 🚀 找到了最合适的 prev
                if (found_best) {
                    auto& prev = *best_it;

                    // ---- 提取 gap ----
                    Coord_t ref_gap_beg = prev->ref_start + prev->ref_len;
                    Coord_t ref_gap_len = (*curr)->ref_start - ref_gap_beg;

                    Coord_t qry_gap_beg, qry_gap_len;
                    if ((*curr)->strand == FORWARD) {
                        qry_gap_beg = prev->qry_start + prev->qry_len;
                        qry_gap_len = (*curr)->qry_start - qry_gap_beg;
                    }
                    else {
                        qry_gap_beg = (*curr)->qry_start + (*curr)->qry_len;
                        qry_gap_len = prev->qry_start - qry_gap_beg;
                    }

                    std::string ref_gap_seq = std::visit([&](auto& p) {
                        using T = std::decay_t<decltype(p)>;
                        if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::SequenceManager>>) {
                            return p->getSubSequence((*curr)->ref_chr_index, ref_gap_beg, ref_gap_len);
                        }
                        else {
                            return p->getOriginalManager().getSubSequence((*curr)->ref_chr_index, ref_gap_beg, ref_gap_len);
                        }
                        }, ref_mgr);

                    std::string qry_gap_seq = std::visit([&](auto& p) {
                        using T = std::decay_t<decltype(p)>;
                        if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::SequenceManager>>) {
                            return p->getSubSequence((*curr)->qry_chr_index, qry_gap_beg, qry_gap_len);
                        }
                        else {
                            return p->getOriginalManager().getSubSequence((*curr)->qry_chr_index, qry_gap_beg, qry_gap_len);
                        }
                        }, qry_mgr);

                    if ((*curr)->strand == REVERSE)
                        reverseComplement(qry_gap_seq);

                    // ---- gap 比对 ----
                    Cigar_t gap_cigar = extendAlignKSW2(ref_gap_seq, qry_gap_seq, 2 * break_len);
                    //Cigar_t gap_cigar = globalAlignKSW2_2(ref_gap_seq, qry_gap_seq);
                    uint_t ref_len = countRefLength(gap_cigar);
                    uint_t qry_len = countQryLength(gap_cigar);

                    // ✅ gap 完全比对上：把 gap_cigar 融入 curr
                    //if (checkGapCigarQuality(gap_cigar, ref_gap_len, qry_gap_len, 0.6))
                    if (ref_len == ref_gap_len && qry_len == qry_gap_len) {
                        (*curr)->ref_len += ref_len;
                        (*curr)->qry_len += qry_len;
                        (*curr)->alignment_length += countAlignmentLength(gap_cigar);
                        (*curr)->aligned_base += countMatchOperations(gap_cigar);

                        if ((*curr)->strand == FORWARD) {
                            (*curr)->ref_start -= ref_len;
                            (*curr)->qry_start -= qry_len;
                            prependCigar((*curr)->cigar, gap_cigar);
                        }
                        else {
                            (*curr)->ref_start -= ref_len;
                            //std::reverse(gap_cigar.begin(), gap_cigar.end());
                            prependCigar((*curr)->cigar, gap_cigar);
                        }
                    }
       //             if (!linked.empty()) {
       //                 // 打印curr和best
       //                 if (strand == FORWARD) {
       //                     std::cout << (int_t)(*curr)->qry_start - (int_t)((*best_it)->qry_start + (*best_it)->qry_len) << std::endl;
       //                 }
       //                 else {
							//std::cout << (int_t)((*best_it)->qry_start) - (int_t)((*curr)->qry_start + (*curr)->qry_len) << std::endl;
       //                 }
       //
       //             }
                }
            }

            linked.push_back(*curr);


            Coord_t curr_ref_end = (*curr)->ref_start + (*curr)->ref_len;
            Coord_t curr_qry_end =
                ((*curr)->strand == FORWARD)
                ? (*curr)->qry_start + (*curr)->qry_len
                : (*curr)->qry_start;

            for (auto it2 = curr + 1; it2 != best; ++it2) {
                if ((*it2)->is_linked) continue;
                Coord_t it_ref_end = (*it2)->ref_start + (*it2)->ref_len;
                Coord_t it_qry_end =
                    ((*curr)->strand == FORWARD)
                    ? (*it2)->qry_start + (*it2)->qry_len
                    : (*it2)->qry_start;

                if (it_ref_end <= curr_ref_end &&
                    (((*curr)->strand == FORWARD && it_qry_end <= curr_qry_end) ||
                        ((*curr)->strand == REVERSE && it_qry_end >= curr_qry_end))) {
                    (*it2)->is_linked = true;
                }


            }
            // 将新的cur设置为和老cur不冲突的anchor
            // TODO考虑反向链

            curr_qry_end = (*curr)->qry_start + (*curr)->qry_len;
            curr_ref_end = (*curr)->ref_start + (*curr)->ref_len;
            int_t cur_qry_start = (*curr)->qry_start;
            //curr = best;
            while (true) {
                if (strand == FORWARD) {
                    //if (!(*curr)->is_linked && (*curr)->qry_start >= curr_qry_end && (*curr)->ref_start >= curr_ref_end) {
                    //    break;
                    //}
                    if (!(*curr)->is_linked && (*curr)->qry_start >= curr_qry_end) {
                        break;
                    }
                }
                else {
                    //if (!(*curr)->is_linked && (*curr)->qry_start + (*curr)->qry_len <= cur_qry_start && (*curr)->ref_start >= curr_ref_end) {
                    //    break;
                    //}
                    if (!(*curr)->is_linked && (*curr)->qry_start + (*curr)->qry_len <= cur_qry_start) {
                        break;
                    }
                }

                ++curr;
            }

            std::cout << "";
            // it = curr + 1;
        }

        //++curr;
    }

    anchors.swap(linked);
    return;



}




