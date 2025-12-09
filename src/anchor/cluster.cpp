#include "anchor.h"
#include "data_process.h"
#include <omp.h>

// ────────────────────────────────────────────
//  Constructors / reset
// ────────────────────────────────────────────
UnionFind::UnionFind(std::size_t n)
{
    reset(n);
}

void UnionFind::reset(std::size_t n)
{
    parent_.assign(n, -1);                   // -1 表示单元素集合
    component_cnt_ = static_cast<int_t>(n);    // 初始有 n 个独立集合
}

// ────────────────────────────────────────────
//  Query
// ────────────────────────────────────────────
int_t UnionFind::find(int_t x)
{
    // 路径折半（迭代版）
    while (parent_[x] >= 0 && parent_[parent_[x]] >= 0)
    {
        parent_[x] = parent_[parent_[x]];
        x = parent_[x];
    }
    return (parent_[x] < 0) ? x : parent_[x];
}

int_t UnionFind::set_size(int_t x)
{
    return -parent_[find(x)];
}

bool UnionFind::same(int_t a, int_t b)
{
    return find(a) == find(b);
}

// ────────────────────────────────────────────
//  Modification
// ────────────────────────────────────────────
bool UnionFind::unite(int_t a, int_t b)
{
    a = find(a);
    b = find(b);
    if (a == b) return false;

    // 按大小合并：把小树挂到大树
    if (parent_[a] > parent_[b]) std::swap(a, b);

    parent_[a] += parent_[b];  // 更新根 a 的大小（负值）
    parent_[b] = a;           // b 挂到 a
    --component_cnt_;
    return true;
}

void filterAndMergeMatches(MatchVec& matches) {
    if (matches.empty()) return;

    // 假设 matches 已经按 qry_start 排好序
    const size_t N = matches.size();
    std::vector<bool> good(N, true);
    std::vector<bool> tentative(N, false);

    for (size_t i = 0; i < N; ++i) {
        if (!good[i]) continue;

        const Match& mi = matches[i];
        int_t i_diag = mi.qry_start - mi.ref_start;
        Coord_t i_end = mi.qry_start + mi.match_len();

        for (size_t j = i + 1; j < N && matches[j].qry_start <= i_end; ++j) {
            if (!good[j]) continue;

            const Match& mj = matches[j];
            int_t j_diag = mj.qry_start - mj.ref_start;

            // --- Case 1: 同一 diagonal ---
            if (i_diag == j_diag &&
                mi.strand() == mj.strand() &&
                mi.ref_chr_index == mj.ref_chr_index &&
                mi.qry_chr_index == mj.qry_chr_index)
            {
                // 合并为更长的
                Coord_t j_extent = mj.match_len() + mj.qry_start - mi.qry_start;
                if (j_extent > matches[i].match_len()) {
                    matches[i].set_match_len(j_extent);
                    i_end = mi.qry_start + j_extent;
                }
                good[j] = false;
            }

            // --- Case 2: 同一 ref 起点 ---
            else if (mi.ref_start == mj.ref_start &&
                mi.ref_chr_index == mj.ref_chr_index)
            {
                int_t overlap = mi.qry_start + mi.match_len() - mj.qry_start;
                if (mi.match_len() < mj.match_len()) {
                    if (overlap >= mi.match_len() / 2) {
                        good[i] = false;
                        break;
                    }
                }
                else if (mj.match_len() < mi.match_len()) {
                    if (overlap >= mj.match_len() / 2)
                        good[j] = false;
                }
                else { // 长度相等
                    if (overlap >= mi.match_len() / 2) {
                        tentative[j] = true;
                        if (tentative[i]) {
                            good[i] = false;
                            break;
                        }
                    }
                }
            }

            // --- Case 3: 同一 qry 起点 ---
            else if (mi.qry_start == mj.qry_start &&
                mi.qry_chr_index == mj.qry_chr_index)
            {
                int overlap = mi.ref_start + mi.match_len() - mj.ref_start;
                if (mi.match_len() < mj.match_len()) {
                    if (overlap >= mi.match_len() / 2) {
                        good[i] = false;
                        break;
                    }
                }
                else if (mj.match_len() < mi.match_len()) {
                    if (overlap >= mj.match_len() / 2)
                        good[j] = false;
                }
                else { // 长度相等
                    if (overlap >= mi.match_len() / 2) {
                        tentative[j] = true;
                        if (tentative[i]) {
                            good[i] = false;
                            break;
                        }
                    }
                }
            }
        }
    }

    // 最后收集 Good 的 matches
    MatchVec filtered;
    filtered.reserve(matches.size());
    for (size_t i = 0; i < N; ++i) {
        if (good[i]) filtered.push_back(matches[i]);
    }
    matches.swap(filtered);
}


// 1. 根据 max_gap / diagdiff / diagfactor 把 unique_match 聚成若干簇 - 优化版本
MatchClusterVec buildClusters(MatchVec& unique_match,
    int_t      max_gap,
    int_t      diagdiff,
    double     diagfactor)
{
    MatchClusterVec clusters;

    // 特殊情况：0 或 1 个元素，直接返回
    if (unique_match.size() < 2) {
        if (unique_match.size() == 1) {
            clusters.emplace_back();
            clusters.back().push_back(std::move(unique_match[0]));
        }
        return clusters;
    }

    const bool is_forward = (unique_match.front().strand() == FORWARD);

    // 先排序
    std::sort(unique_match.begin(), unique_match.end(),
        [](const Match& a, const Match& b) {
            // 这里我也顺手修个小笔误，后面解释
            if (start2(a) < start2(b)) return true;
            if (start2(a) > start2(b)) return false;
            return start1(a) < start1(b);
        });

    // 再合并压缩（可能会删元素、缩短 vector）
    filterAndMergeMatches(unique_match);

    // 现在重新拿 N
    const uint_t N = static_cast<uint_t>(unique_match.size());
    if (N < 2) {
        if (N == 1) {
            clusters.emplace_back();
            clusters.back().push_back(std::move(unique_match[0]));
        }
        return clusters;
    }

    UnionFind uf(N);

 
    // 优化的聚类算法：利用排序后的局部性
    for (uint_t i = 0; i < N; ++i) {
        uint_t i_end = start2(unique_match[i]) + len2(unique_match[i]);
        int_t  i_diag = 0;
        if (is_forward) {
            i_diag = diag(unique_match[i]);
        }
        else {
			i_diag = diag_reverse(unique_match[i]);
        }
        //i_diag = diag(unique_match[i]);
        

        // 只检查后续可能的匹配，利用排序优化
        for (uint_t j = i + 1; j < N; ++j) {
            int_t sep = start2(unique_match[j]) - i_end;

            // 早期退出：如果gap太大，后续的j也不可能匹配
            if (sep > static_cast<int_t>(max_gap)) break;
            int_t diag_diff = 0;
            //diag_diff = std::abs(diag(unique_match[j]) - i_diag);
            if (is_forward) {
				diag_diff = std::abs(diag(unique_match[j]) - i_diag);
			}
            else {
                diag_diff = std::abs(diag_reverse(unique_match[j]) - i_diag);
            }

            // int_t diag_diff = std::abs((is_forward ? diag(unique_match[j]):diag_reverse(unique_match[j])) - i_diag);
            int_t th = std::max(diagdiff, static_cast<int_t>(diagfactor * sep));

            if (diag_diff <= th) {
                uf.unite(i, j);
            }
        }
    }
    
    // 高效的簇构建：使用哈希表而不是线性搜索
    std::unordered_map<int_t, int_t> root_to_cluster_id;
    root_to_cluster_id.reserve(N / 4);  // 预估簇数量
    
    for (uint_t idx = 0; idx < unique_match.size(); idx++) {
        int_t root = uf.find(idx);
        auto it = root_to_cluster_id.find(root);
        
        int_t cid;
        if (it == root_to_cluster_id.end()) {
            cid = static_cast<int_t>(clusters.size());
            clusters.emplace_back();
            clusters.back().reserve(1);  // 预分配空间
            root_to_cluster_id[root] = cid;
        } else {
            cid = it->second;
        }
        if (unique_match[idx].ref_chr_index > 0) {
            std::cout << "";
        }
        clusters[cid].push_back(std::move(unique_match[idx]));
    }
    
    return clusters;
}

// 稀疏版：只聚类真实存在的 key
ClusterBySQR_SparsePtr clusterAllChrMatchSparse(
    MatchBySQR_SparsePtr& unique_anchors,
    MatchBySQR_SparsePtr& repeat_anchors,  // 按你的要求：占位不用
    uint_t min_span,
    uint_t thread_num)
{
    // 输出：每个 key -> 一组 MatchCluster
    auto out = std::make_shared<ClusterBySQR_Sparse>();

    // 如果 unique_anchors 为空指针或本身没有任何 key，直接返回空结果
    if (!unique_anchors || unique_anchors->empty()) {
        return out;
    }

    // 预留空间，减少 unordered_map rehash
    out->reserve(unique_anchors->size());
    std::vector<uint64_t> keys;
    keys.reserve(unique_anchors->size());

    for (const auto& kv : *unique_anchors) {
        keys.push_back(kv.first);
    }


    std::vector<std::pair<uint64_t, MatchClusterVecPtr>> local_results(keys.size());

#pragma omp parallel for schedule(dynamic) num_threads(thread_num)
    for (long long i = 0; i < (long long)keys.size(); ++i) {
        uint64_t key = keys[i];

        MatchVec mv = unique_anchors->at(key);
        auto clusters = clusterChrMatch(mv, min_span);

        local_results[i] = {key, std::move(clusters)};
    }

    // 并行结束后，串行写回 out
    out->reserve(keys.size());
    for (auto& kv : local_results) {
        if (!kv.second->empty()) {
            (*out)[kv.first] = std::move(kv.second);
        }
    }
    return out;
}

//
// ClusterVecPtrByStrandByQueryRefPtr
// clusterAllChrMatch(const MatchByStrandByQueryRefPtr& unique_anchors,
//                    const MatchByStrandByQueryRefPtr& repeat_anchors,
//                    uint_t min_cluster_length)
// {
//     // ---------- 0. 判空 ----------
//     if ((!unique_anchors || unique_anchors->empty()) &&
//         (!repeat_anchors || repeat_anchors->empty())) {
//         spdlog::warn("[clusterAllChrAnchors] empty anchor ptr");
//         return std::make_shared<ClusterVecPtrByStrandByQueryRef>();
//     }
//
//     // ---------- 1. 选择维度基准（避免 unique 为空时 front 越界） ----------
//     const auto& base = (unique_anchors && !unique_anchors->empty())
//                         ? unique_anchors
//                         : repeat_anchors;
//
//     // 防御式检查
//     if (!base || base->empty() || (*base)[0].empty() || (*base)[0][0].empty()) {
//         return std::make_shared<ClusterVecPtrByStrandByQueryRef>();
//     }
//
//     const size_t query_chr_n = (*base)[0].size();
//     const size_t ref_chr_n   = (*base)[0][0].size();
//
//     // ---------- 2. 构建结果桶 ----------
//     auto cluster_ptr = std::make_shared<ClusterVecPtrByStrandByQueryRef>();
//     cluster_ptr->resize(2);  // strand 维度
//
//     for (auto& query_layer : *cluster_ptr) {
//         query_layer.resize(query_chr_n);
//         for (auto& ref_row : query_layer) {
//             ref_row.resize(ref_chr_n);
//         }
//     }
//
//     // ---------- 3. 收集非空任务 ----------
//     struct ClusterTask {
//         uint_t k, i, j;
//         MatchVec* uniq_ptr;
//         MatchVec* rept_ptr;
//         size_t total_size;
//
//         ClusterTask(uint_t _k, uint_t _i, uint_t _j, MatchVec* _uniq, MatchVec* _rept)
//             : k(_k), i(_i), j(_j), uniq_ptr(_uniq), rept_ptr(_rept),
//               total_size((_uniq ? _uniq->size() : 0) + (_rept ? _rept->size() : 0)) {}
//     };
//
//     std::vector<ClusterTask> tasks;
//     tasks.reserve(1000);
//
//     for (uint_t k = 0; k < 2; ++k) {
//         // 假设 unique/repeat 的形状一致；若不一致可加更严格检查
//         for (uint_t i = 0; i < query_chr_n; ++i) {
//             for (uint_t j = 0; j < ref_chr_n; ++j) {
//
//                 MatchVec& uniq = (*unique_anchors)[k][i][j];
//                 MatchVec& rept = (*repeat_anchors)[k][i][j];
//
//                 if (!uniq.empty() || !rept.empty()) {
//                     tasks.emplace_back(k, i, j, &uniq, &rept);
//                 }
//             }
//         }
//     }
//
//     if (tasks.empty()) {
//         return cluster_ptr;
//     }
//
//     // ---------- 4. 按任务大小排序（可保留，也可去掉） ----------
//     std::sort(tasks.begin(), tasks.end(),
//         [](const ClusterTask& a, const ClusterTask& b) {
//             return a.total_size > b.total_size;
//         });
//
//     // ---------- 5. OpenMP 并行处理 tasks ----------
//     // 说明：
//     // - schedule(dynamic) 适合你这种任务大小不均
//     // - 是否显式 num_threads 取决于你的上层是否已经 omp_set_num_threads(thread_num)
//     //   如果你希望在这里强制线程数，可以取消注释 num_threads(...)
//     //
//     // #pragma omp parallel for schedule(dynamic) num_threads(thread_num)
//     #pragma omp parallel for schedule(dynamic)
//     for (long long t = 0; t < static_cast<long long>(tasks.size()); ++t) {
//         const auto& task = tasks[static_cast<size_t>(t)];
//
//         auto cluster_result = clusterChrMatch(
//             *task.uniq_ptr,
//             *task.rept_ptr,
//             min_cluster_length
//         );
//
//         // 每个 task 写不同的 [k][i][j]，不冲突
//         (*cluster_ptr)[task.k][task.i][task.j] = std::move(cluster_result);
//     }
//
//     return cluster_ptr;
// }


/*!
 * \brief  把 src（三维：strand × query × ref）按 ref 合并
 *         每个 ref 得到一个合并后的 MatchClusterVecPtr
 *
 * \param  src  输入三维结构（shared_ptr 外围保持不变）
 * \return      shared_ptr< vector< MatchClusterVecPtr > >
 */
ClusterVecPtrByRefPtr
groupClustersByRef(const ClusterVecPtrByStrandByQueryRefPtr& src)
{
    /* ---------- 0. 早退 ---------- */
    auto by_ref = std::make_shared<ClusterVecPtrByRef>();   // 最终返回对象
    if (!src || src->empty() || (*src)[0].empty()) return by_ref;

    /* ---------- 1. 基础信息 ---------- */
    const size_t ref_n = (*src)[0][0].size();               // 所有行列数一致

    /* ---------- 2. 统计每个 ref 的簇总量 & 记录首矢量指针 ---------- */
    std::vector<size_t> cluster_counts(ref_n, 0);
    std::vector<MatchClusterVec*> first_src_vec(ref_n, nullptr);

    for (auto& query_layer : *src)          // strand
        for (auto& ref_row : query_layer)   // query
            for (size_t r = 0; r < ref_n; ++r) {
                auto& vec = *ref_row[r];
                cluster_counts[r] += vec.size();
                if (!first_src_vec[r] && !vec.empty())
                    first_src_vec[r] = &vec;             // 只记录第一次出现
            }

    /* ---------- 3. 构造目标数组，首矢量直接 swap ---------- */
    by_ref->resize(ref_n);
    std::vector<bool> swapped(ref_n, false);

    for (size_t r = 0; r < ref_n; ++r) {
        auto dst = std::make_shared<MatchClusterVec>();
        if (first_src_vec[r]) {                          // 有首矢量就直接搬进来
            dst->swap(*first_src_vec[r]);                // 零拷贝
            swapped[r] = true;
        }
        dst->reserve(cluster_counts[r]);                 // 再保证容量充足
        (*by_ref)[r] = std::move(dst);
    }

    /* ---------- 4. 把其余簇 append 进去 ---------- */
    for (auto& query_layer : *src)
        for (auto& ref_row : query_layer)
            for (size_t r = 0; r < ref_n; ++r) {
                auto& src_vec = *ref_row[r];
                if (src_vec.empty()) continue;

                auto& dst = *(*by_ref)[r];
                if (swapped[r]) {
                    dst.insert(dst.end(),
                        std::make_move_iterator(src_vec.begin()),
                        std::make_move_iterator(src_vec.end()));
                }
                else {                                 // 之前全空，直接 swap
                    dst.swap(src_vec);
                    swapped[r] = true;
                }
            }

    return by_ref;                                       // shared_ptr 返回
}



// 检查两个锚点是否overlap
inline bool isOverlap(const Match& a, const Match& b) {
    // 检查ref维度overlap
    bool ref_overlap = !(a.ref_start + a.match_len() <= b.ref_start ||
        b.ref_start + b.match_len() <= a.ref_start);
    
    // 检查query维度overlap
    bool query_overlap = !(a.qry_start + a.match_len() <= b.qry_start ||
                          b.qry_start + b.match_len() <= a.qry_start);
    
    return ref_overlap && query_overlap;
}

// 2. 给定一个 cluster，返回其最佳非交叉链（DP O(N^2)）- 支持overlap合并
MatchVec bestChainDP(MatchVec& cluster, double diagfactor)
{
    if (cluster.empty()) return {};
    if (cluster.size() == 1) return cluster;

    Strand strand = cluster.front().strand();

    std::sort(cluster.begin(), cluster.end(),
        [](const Match& a, const Match& b) { return start2(a) < start2(b); });

    const uint_t N = static_cast<uint_t>(cluster.size());
    std::vector<int_t> score(N), pred(N, -1);
    uint_t best_idx = 0;

    for (uint_t i = 0; i < N; ++i) {
        score[i] = len2(cluster[i]);
        for (uint_t j = 0; j < i; ++j) {
            if (start2(cluster[i]) <= start2(cluster[j]) + len2(cluster[j])) continue;

            int_t sep = 0;
			int_t d = 0;
            if (strand == FORWARD) {
                int_t prev_endj = start1(cluster[j]) + len1(cluster[j]);
                if (start1(cluster[i]) <= prev_endj) continue;
                sep = start1(cluster[i]) - prev_endj;
                d = std::abs(diag(cluster[i]) - diag(cluster[j]));
            } else {
                int_t prev_endi = start1(cluster[i]) + len1(cluster[i]);
                if(prev_endi >= start1(cluster[j])) continue;
				sep = start1(cluster[j]) - prev_endi;
                d = std::abs(diag_reverse(cluster[i]) - diag_reverse(cluster[j]));
            }

             
            int_t cand = score[j] + len2(cluster[i]) - d;
            if (cand > score[i]) { score[i] = cand; pred[i] = static_cast<int_t>(j); }
        }
        if (score[i] > score[best_idx]) best_idx = i;
    }

    MatchVec chain;
    for (int_t k = static_cast<int_t>(best_idx); k != -1; k = pred[k])
        chain.emplace_back(cluster[k]);
    std::reverse(chain.begin(), chain.end());

    return chain;
}

/* ───────────────────────────────────────────────────────── *
 * 对外主函数                                              *
 * ───────────────────────────────────────────────────────── */
MatchClusterVecPtr clusterChrMatch(MatchVec& unique_match, uint_t min_cluster_length, int_t max_gap, int_t diagdiff, double diagfactor)
{

    auto best_chain_clusters = std::make_shared<MatchClusterVec>();
    if (unique_match.size() == 0) {
        return best_chain_clusters;
    }
    // 1. 聚簇
    MatchClusterVec clusters = buildClusters(unique_match, max_gap, diagdiff, diagfactor);

   // 释放unique_match的内存
	unique_match.clear();
	unique_match.shrink_to_fit();  // 释放内存

    
    best_chain_clusters->reserve(clusters.size());
    for (auto& cluster : clusters) {
        if (cluster.empty()) continue;

        MatchVec best_chain = bestChainDP(cluster, diagfactor);

        if (best_chain.empty()) { releaseCluster(cluster); continue; }
        int_t span = 0;
		// 遍历 best_chain 计算 span
        for (auto& m : best_chain) {
			span += m.match_len();
        }
        //int_t span = start1(best_chain.back()) + len1(best_chain.back()) - start1(best_chain.front());
        if (span >= min_cluster_length) {

            best_chain_clusters->emplace_back(std::move(best_chain));
        }
        /*else {
            std::move(best_chain.begin(), best_chain.end(), std::back_inserter(repeat_match));
        }*/
        releaseCluster(cluster);     // 回收 cluster 剩余元素
    }
    best_chain_clusters->shrink_to_fit();

    return best_chain_clusters;
}

 /*!
  * \brief  把 cl 依照冲突矩形切成 ≤2 段
  *         · 仅当 ref_hit 为 true 时才检查 [bad_r_beg,bad_r_end)
  *         · 仅当 query_hit 为 true 时才检查 [bad_q_beg,bad_q_end)
  *
  * \return 0 段：整簇都冲突；1 段：完全无冲突；2 段：被矩形切成前/后两段
  */
MatchClusterVec
splitCluster(const MatchCluster& cl,
    bool  ref_hit,
    int_t bad_r_beg, int_t bad_r_end,
    bool  query_hit,
    int_t bad_q_beg, int_t bad_q_end)
{
    MatchClusterVec parts;
    if (cl.empty()) return parts;

    /* ---------- 1. λ 判断“命中” ---------- */
    auto hit = [&](const Match& m) -> bool
        {
            bool r_overlap = false;   
            bool q_overlap = false;

            if (ref_hit) {
                int_t rb = start1(m), re = rb + len1(m);
                r_overlap = !(re <= bad_r_beg || rb >= bad_r_end);
            }
            if (query_hit) {
                int_t qb = start2(m), qe = qb + len2(m);
                q_overlap = !(qe <= bad_q_beg || qb >= bad_q_end);
            }

            /* 两个维度都要满足要检查的那部分 */
            return r_overlap || q_overlap;
        };

    /* ---------- 2. 找第一次 / 最后一次命中 ---------- */
    size_t first_hit = cl.size();   // 未命中
    size_t last_hit = 0;

    for (size_t i = 0; i < cl.size(); ++i)
        if (hit(cl[i])) {
            first_hit = std::min(first_hit, i);
            last_hit = i;
        }

    /* ---------- 3. 根据命中情况切分 ---------- */
    if (first_hit == cl.size()) {              // 完全不冲突
        parts.emplace_back(cl);
        return parts;
    }

    if (first_hit > 0)                         // 前段
        parts.emplace_back(cl.begin(), cl.begin() + first_hit);

    if (last_hit + 1 < cl.size())              // 后段
        parts.emplace_back(cl.begin() + last_hit + 1, cl.end());

    return parts;                              // 0、1 或 2 段
}



/*!
 * \brief  在原 vector 上直接执行“两级 map + 贪婪拆分”过滤
 * \param  clusters_ptr   ref 上的所有簇
 * \param  MIN_SPAN       最小跨度阈值
 * \param  MIN_MATCH      最小匹配数（若不想改，可写死 50）
 */
void keepWithSplitGreedy(MatchClusterVecPtr clusters_ptr,
    int_t              MIN_SPAN)
{
    if (!clusters_ptr || clusters_ptr->empty()) return;

    /* ---------- 1. 最大堆 ---------- */
    struct Node { MatchCluster cl; int_t sc; };
    auto cmp = [](const Node& a, const Node& b) { return a.sc < b.sc; };

    std::vector<Node> heap;  heap.reserve(clusters_ptr->size());
    for (auto& cl : *clusters_ptr) {
        if (cl.empty()) continue;

        int_t sc = clusterSpan(cl);                // ① 先安全读取
        heap.push_back({ std::move(cl), sc });     // ② 再移动进堆
    }
    clusters_ptr->clear();                              // **清空原容器**
    std::make_heap(heap.begin(), heap.end(), cmp);

    /* ---------- 2. 两级 interval map ---------- */
    IntervalMap refMap;
    std::vector<IntervalMap> qMaps;

    /* ---------- 3. 贪婪循环 ---------- */
    while (!heap.empty())
    {
        std::pop_heap(heap.begin(), heap.end(), cmp);
        Node cur = std::move(heap.back()); heap.pop_back();
        if (cur.sc < MIN_SPAN) continue;

        uint_t rb = start1(cur.cl.front());
        uint_t re = start1(cur.cl.back()) + len1(cur.cl.back());
        uint_t qb = start2(cur.cl.front());
        uint_t qe = start2(cur.cl.back()) + len2(cur.cl.back());
        const ChrIndex& qChr = cur.cl.front().qry_chr_index;

        int_t RL = 0, RR = 0, QL = 0, QR = 0;
        bool ref_hit = overlap1D(refMap, rb, re, RL, RR);
        bool query_hit = overlap1D(qMaps[qChr], qb, qe, QL, QR);

        if (!ref_hit && !query_hit) {
            insertInterval(refMap, rb, re);
            insertInterval(qMaps[qChr], qb, qe);
            clusters_ptr->emplace_back(std::move(cur.cl));   // **写回原 vector**
            continue;
        }


        for (auto& part : splitCluster(cur.cl, ref_hit, RL, RR, query_hit, QL, QR)) {
            int_t sc = clusterSpan(part);
            if (sc >= MIN_SPAN) {
                heap.push_back({ std::move(part), sc });
                std::push_heap(heap.begin(), heap.end(), cmp);
            }

        }
        /* clusters_ptr 现在就是过滤后的结果 */
    }
}

/*------------------------------------------------------------------*
 *  by_ref        : shared_ptr<vector<MatchClusterVecPtr>>  (ref 维)
 *  返回值 by_refQ : shared_ptr<vector< vector<MatchClusterVecPtr> >>
 *                   └──ref──┘└────query────┘
 *------------------------------------------------------------------*/
ClusterVecPtrByRefQueryPtr
groupClustersByRefQuery(const ClusterVecPtrByRefPtr& by_ref,
    SeqPro::ManagerVariant& query_fasta_manager,
    ThreadPool& pool)
{
    /* ---------- 0. 提前创建返回对象 ---------- */
    auto by_ref_query = std::make_shared<ClusterVecPtrByRefQuery>();
    if (!by_ref || by_ref->empty()) return by_ref_query;

    /* ---------- 1. 查询染色体总数 ---------- */
    const uint32_t query_cnt = std::visit([](auto&& mgr) {
        return static_cast<uint32_t>(mgr->getSequenceCount());
        }, query_fasta_manager);

    /* ---------- 2. 预分配 [ref][query] 矩阵 ---------- */
    const uint32_t ref_cnt = static_cast<uint32_t>(by_ref->size());
    by_ref_query->resize(ref_cnt);
    for (uint32_t r = 0; r < ref_cnt; ++r)
        (*by_ref_query)[r].resize(query_cnt);          // 初值 = 空 vector

    /* ---------- 3. 每个 ref 开一个线程任务 ---------- */
    for (uint32_t r = 0; r < ref_cnt; ++r) {
        auto src_ptr = (*by_ref)[r];                   // 可能为空
        if (!src_ptr || src_ptr->empty()) continue;

        /* 捕获：src_ptr 按值、by_ref_query 按 shared_ptr 值 */
        pool.enqueue([src_ptr,
            by_ref_query,
            &query_fasta_manager,
            r, query_cnt]
            {
                /* 指向目标行，避免 lambda 内多级索引 */
                ClusterVecPtrByRef* dst_row = &(*by_ref_query)[r];

                for (auto& cluster : *src_ptr) {
                    if (cluster.empty()) continue;

                    /* ---- 3-1 计算 query 染色体 ID ---- */
                    const auto& m = cluster.front();
                    uint32_t qIdx = m.qry_chr_index;

                    if (qIdx >= query_cnt) continue;      // 无效或未收录

                    /* ---- 3-2 若目标指针为空就创建 ---- */
                    if (!dst_row->at(qIdx))
                        dst_row->at(qIdx) = std::make_shared<MatchClusterVec>();

                    /* ---- 3-3 把簇移动到目标桶 ---- */
                    dst_row->at(qIdx)->emplace_back(std::move(cluster));
                }
                /* src_ptr 中元素已 move，src_ptr 本身随后被析构 */
            });
    }

    /* ---------- 4. 同步，保证数据就绪 ---------- */
    return by_ref_query;
}


/* ============================================================= *
 *  把三维 clusters  ->  按 ref 的一维 clusters
 *  并行对每个 ref 走 keepWithSplitGreedy，再写回三维结构
 * ============================================================= */
void filterClustersByGreedy(ClusterVecPtrByRefPtr by_ref,
    ThreadPool& pool,
    int_t                                min_span)
{

    if (!by_ref || by_ref->empty()) return;

    for (size_t r = 0; r < by_ref->size(); ++r) {
        pool.enqueue([&, r] {
            keepWithSplitGreedy((*by_ref)[r], min_span);
            });
    }

    return;

}

