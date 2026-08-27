#include "anchor.h"
#include "data_process.h"
#include <omp.h>
#include <unordered_map>
#include <iostream>

// ────────────────────────────────────────────
// UnionFind：并查集
// - parent_[x] < 0 表示 x 为根，且集合大小为 -parent_[x]
// - parent_[x] >=0 表示 parent_[x] 为父节点下标
// ────────────────────────────────────────────

// 构造：创建 n 个独立集合
UnionFind::UnionFind(std::size_t n) {
    reset(n);
}

// 重置并查集为 n 个独立集合
void UnionFind::reset(std::size_t n) {
    parent_.assign(n, -1);                      // -1 表示单元素集合（根节点，大小=1）
    component_cnt_ = static_cast<int_t>(n);     // 初始有 n 个连通分量
}

// ────────────────────────────────────────────
// Query（查询接口）
// ────────────────────────────────────────────

// 查找：返回元素 x 所在集合的根
int_t UnionFind::find(int_t x) {
    // 路径折半（迭代版）：不断把 x 挂到祖父节点上，加速后续查询
    while (parent_[x] >= 0 && parent_[parent_[x]] >= 0) {
        parent_[x] = parent_[parent_[x]];
        x = parent_[x];
    }
    return (parent_[x] < 0) ? x : parent_[x];
}

// 返回 x 所在集合大小
int_t UnionFind::set_size(int_t x) {
    return -parent_[find(x)];
}

// 判断 a 和 b 是否属于同一集合
bool UnionFind::same(int_t a, int_t b) {
    return find(a) == find(b);
}

// ────────────────────────────────────────────
// Modification（修改接口）
// ────────────────────────────────────────────

// 合并两个集合（按大小合并）
// 返回 true 表示确实发生了合并；false 表示原本就在同一集合
bool UnionFind::unite(int_t a, int_t b) {
    a = find(a);
    b = find(b);
    if (a == b) return false;

    // 按大小合并：parent_ 为负表示大小，数值越小(更负)集合越大
    // 若 a 的集合更小，则交换，让 a 做大根
    if (parent_[a] > parent_[b]) std::swap(a, b);

    parent_[a] += parent_[b];   // 更新根 a 的大小（负值累加）
    parent_[b] = a;             // b 挂到 a
    --component_cnt_;
    return true;
}

// ────────────────────────────────────────────
// filterAndMergeMatches：对 match 列表进行过滤与合并
// 假设：matches 已按 qry_start 排序
// 规则：
// 1) 同一 diagonal 且同向同染色体：合并成更长 match，删除后者
// 2) 同 ref 起点：根据重叠比例过滤较短者；同长时用 tentative 机制处理
// 3) 同 qry 起点：同上
// ────────────────────────────────────────────
void filterAndMergeMatches(MatchVec& matches) {
    if (matches.empty()) return;

    const size_t N = matches.size();
    std::vector<bool> good(N, true);        // good[i]=true 表示保留
    std::vector<bool> tentative(N, false);  // tentative 标记：用于“长度相等且重叠较大”的歧义处理

    for (size_t i = 0; i < N; ++i) {
        if (!good[i]) continue;

        const Match& mi = matches[i];
        int_t  i_diag = mi.qry_start - mi.ref_start;
        Coord_t i_end = mi.qry_start + mi.match_len();  // i 在 query 维度的结束位置（开区间右端）

        // 由于 matches 按 qry_start 排序，只需检查 qry_start <= i_end 的后续元素
        for (size_t j = i + 1; j < N && matches[j].qry_start <= i_end; ++j) {
            if (!good[j]) continue;

            const Match& mj = matches[j];
            int_t j_diag = mj.qry_start - mj.ref_start;

            // --- Case 1: 同一 diagonal（同一条对角线） ---
            // 条件：diag 相等 + strand 相等 + ref/qry 染色体一致
            if (i_diag == j_diag &&
                mi.strand() == mj.strand() &&
                mi.ref_chr_index == mj.ref_chr_index &&
                mi.qry_chr_index == mj.qry_chr_index) {

                // 合并为更长的：计算 mj 相对 mi 的延伸长度
                Coord_t j_extent = mj.match_len() + mj.qry_start - mi.qry_start;
                if (j_extent > matches[i].match_len()) {
                    matches[i].set_match_len(j_extent);
                    i_end = mi.qry_start + j_extent;
                }
                good[j] = false; // 删除 mj
            }

            // --- Case 2: 同一 ref 起点 ---
            else if (mi.ref_start == mj.ref_start &&
                     mi.ref_chr_index == mj.ref_chr_index) {

                int_t overlap = mi.qry_start + mi.match_len() - mj.qry_start;

                if (mi.match_len() < mj.match_len()) {
                    // i 更短：若重叠超过 i 一半，丢弃 i
                    if (overlap >= mi.match_len() / 2) {
                        good[i] = false;
                        break;
                    }
                }
                else if (mj.match_len() < mi.match_len()) {
                    // j 更短：若重叠超过 j 一半，丢弃 j
                    if (overlap >= mj.match_len() / 2)
                        good[j] = false;
                }
                else {
                    // 长度相等：若重叠超过一半，标记 tentative
                    if (overlap >= mi.match_len() / 2) {
                        tentative[j] = true;
                        if (tentative[i]) {
                            // 如果 i 也已被 tentative，则丢弃 i
                            good[i] = false;
                            break;
                        }
                    }
                }
            }

            // --- Case 3: 同一 qry 起点 ---
            else if (mi.qry_start == mj.qry_start &&
                     mi.qry_chr_index == mj.qry_chr_index) {

                int64_t overlap = static_cast<int64_t>(mi.ref_start) +
                                  static_cast<int64_t>(mi.match_len()) -
                                  static_cast<int64_t>(mj.ref_start);

                if (mi.match_len() < mj.match_len()) {
                    if (overlap >= static_cast<int64_t>(mi.match_len() / 2)) {
                        good[i] = false;
                        break;
                    }
                }
                else if (mj.match_len() < mi.match_len()) {
                    if (overlap >= static_cast<int64_t>(mj.match_len() / 2))
                        good[j] = false;
                }
                else {
                    if (overlap >= static_cast<int64_t>(mi.match_len() / 2)) {
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

    // 收集 good 的 matches，生成新数组
    MatchVec filtered;
    filtered.reserve(matches.size());
    for (size_t i = 0; i < N; ++i) {
        if (good[i]) filtered.push_back(matches[i]);
    }
    matches.swap(filtered);
}

// ────────────────────────────────────────────
// buildClusters：将 unique_match 根据 max_gap / diagdiff / diagfactor 聚类
// 逻辑：
// 1) 对 match 排序（按 start2，再按 start1）
// 2) filterAndMergeMatches 做合并压缩
// 3) 使用并查集：若两个 match 在 query 维度距离 sep <= max_gap
//    且 diagonal 差 <= max(diagdiff, diagfactor * sep)，则归为一类
// 4) 最后根据并查集根构建簇列表
// ────────────────────────────────────────────
MatchClusterVec buildClusters(MatchVec& unique_match,
                             int_t  max_gap,
                             int_t  diagdiff,
                             double diagfactor) {
    MatchClusterVec clusters;

    // 特殊情况：0 或 1 个元素，直接返回
    if (unique_match.size() < 2) {
        if (unique_match.size() == 1) {
            clusters.emplace_back();
            clusters.back().push_back(std::move(unique_match[0]));
        }
        return clusters;
    }

    // 判断链方向（假设同一批次一致）
    const bool is_forward = (unique_match.front().strand() == FORWARD);

    // 先排序：按 start2 再按 start1
    std::sort(unique_match.begin(), unique_match.end(),
        [](const Match& a, const Match& b) {
            if (start2(a) < start2(b)) return true;
            if (start2(a) > start2(b)) return false;
            return start1(a) < start1(b);
        });

    // 再合并压缩（可能会删除元素、缩短 vector）
    filterAndMergeMatches(unique_match);

    // 重新获取 N
    const uint_t N = static_cast<uint_t>(unique_match.size());
    if (N < 2) {
        if (N == 1) {
            clusters.emplace_back();
            clusters.back().push_back(std::move(unique_match[0]));
        }
        return clusters;
    }

    UnionFind uf(N);

    // 利用排序后的局部性进行聚类：对每个 i，只检查后续 sep 不超过 max_gap 的 j
    for (uint_t i = 0; i < N; ++i) {
        uint_t i_end = start2(unique_match[i]) + len2(unique_match[i]);

        int_t i_diag = 0;
        if (is_forward) {
            i_diag = diag(unique_match[i]);
        }
        else {
            i_diag = diag_reverse(unique_match[i]);
        }

        for (uint_t j = i + 1; j < N; ++j) {
            int_t sep = start2(unique_match[j]) - i_end;

            // 早停：gap 太大则后面的 j 也不可能满足
            if (sep > static_cast<int_t>(max_gap)) break;

            int_t diag_diff = 0;
            if (is_forward) {
                diag_diff = std::abs(diag(unique_match[j]) - i_diag);
            }
            else {
                diag_diff = std::abs(diag_reverse(unique_match[j]) - i_diag);
            }

            // 阈值：max(diagdiff, diagfactor * sep)
            int_t th = std::max(diagdiff, static_cast<int_t>(diagfactor * sep));

            // diagonal 差满足阈值：归并到同一簇
            if (diag_diff <= th) {
                uf.unite(i, j);
            }
        }
    }

    // 根据并查集根构建簇：root -> cluster_id
    std::unordered_map<int_t, int_t> root_to_cluster_id;
    root_to_cluster_id.reserve(N / 4);

    for (uint_t idx = 0; idx < unique_match.size(); idx++) {
        int_t root = uf.find(idx);
        auto it = root_to_cluster_id.find(root);

        int_t cid;
        if (it == root_to_cluster_id.end()) {
            cid = static_cast<int_t>(clusters.size());
            clusters.emplace_back();
            clusters.back().reserve(1);
            root_to_cluster_id[root] = cid;
        }
        else {
            cid = it->second;
        }

        clusters[cid].push_back(std::move(unique_match[idx]));
    }

    return clusters;
}

// ────────────────────────────────────────────
// clusterAllChrMatchSparse：稀疏版聚类（仅处理 unique_anchors 里实际存在的 key）
// 输入：
// - unique_anchors：key -> MatchVec
// - repeat_anchors：占位不用
// 输出：
// - key -> MatchClusterVecPtr（每个 key 的聚类结果）
// 说明：
// - 先收集所有 keys
// - OpenMP 并行处理每个 key
// - 并行结束后串行写回 out
// ────────────────────────────────────────────
ClusterBySQR_SparsePtr clusterAllChrMatchSparse(
    MatchBySQR_SparsePtr& unique_anchors,
    MatchBySQR_SparsePtr& repeat_anchors,
    uint_t min_span,
    uint_t thread_num,
    bool include_repeats) {
    auto out = std::make_shared<ClusterBySQR_Sparse>();

    const bool unique_empty = !unique_anchors || unique_anchors->empty();
    const bool repeat_empty = !repeat_anchors || repeat_anchors->empty();
    if (unique_empty && (!include_repeats || repeat_empty)) {
        return out;
    }

    const size_t key_capacity =
        (unique_empty ? 0 : unique_anchors->size()) +
        (include_repeats && !repeat_empty ? repeat_anchors->size() : 0);
    out->reserve(key_capacity);

    std::vector<uint64_t> keys;
    keys.reserve(key_capacity);
    if (!unique_empty) {
        for (const auto& kv : *unique_anchors) {
            keys.push_back(kv.first);
        }
    }
    if (include_repeats && !repeat_empty) {
        for (const auto& kv : *repeat_anchors) {
            if (unique_empty || !unique_anchors->contains(kv.first)) {
                keys.push_back(kv.first);
            }
        }
    }

    // 每个 key 对应一个结果（key, clusters）
    std::vector<std::pair<uint64_t, MatchClusterVecPtr>> local_results(keys.size());

#pragma omp parallel for schedule(dynamic) num_threads(thread_num)
    for (long long i = 0; i < (long long)keys.size(); ++i) {
        uint64_t key = keys[i];

        MatchVec mv;
        if (!unique_empty) {
            if (const auto it = unique_anchors->find(key);
                it != unique_anchors->end()) {
                mv = it->second;
            }
        }
        if (include_repeats && !repeat_empty) {
            if (const auto it = repeat_anchors->find(key);
                it != repeat_anchors->end()) {
                mv.insert(mv.end(), it->second.begin(), it->second.end());
            }
        }

        auto clusters = clusterChrMatch(mv, min_span);

        local_results[i] = { key, std::move(clusters) };
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

// ────────────────────────────────────────────
// isOverlap：检查两个锚点是否在 ref 与 query 两个维度上都发生重叠
// ────────────────────────────────────────────
inline bool isOverlap(const Match& a, const Match& b) {
    // ref 维度重叠：区间 [start, start+len) 是否相交
    bool ref_overlap = !(a.ref_start + a.match_len() <= b.ref_start ||
                         b.ref_start + b.match_len() <= a.ref_start);

    // query 维度重叠
    bool query_overlap = !(a.qry_start + a.match_len() <= b.qry_start ||
                           b.qry_start + b.match_len() <= a.qry_start);

    return ref_overlap && query_overlap;
}

// ────────────────────────────────────────────
// bestChainDP：给定一个 cluster，返回其最佳非交叉链（DP O(N^2)）
// 说明：
// - 对 cluster 按 start2 排序
// - DP：score[i] 表示以 i 结尾的最佳链分数
// - pred[i] 记录前驱
// - 打分：cand = score[j] + len2(i) - d（d 为 diagonal 差惩罚）
// 注意：
// - 目前代码以“非交叉”为条件（通过 start1/start2 及 len 判断）
// - isOverlap() 在此函数中未使用（保持原逻辑）
// ────────────────────────────────────────────
MatchVec bestChainDP(MatchVec& cluster, double diagfactor) {
    if (cluster.empty()) return {};
    if (cluster.size() == 1) return cluster;

    Strand strand = cluster.front().strand();

    std::sort(cluster.begin(), cluster.end(),
        [](const Match& a, const Match& b) { return start2(a) < start2(b); });

    const uint_t N = static_cast<uint_t>(cluster.size());
    std::vector<int_t> score(N), pred(N, -1);
    uint_t best_idx = 0;

    for (uint_t i = 0; i < N; ++i) {
        // 初始分：自身长度
        score[i] = len2(cluster[i]);

        for (uint_t j = 0; j < i; ++j) {
            // query 维度必须不重叠（i 的 start2 需在 j 的末尾之后）
            if (start2(cluster[i]) <= start2(cluster[j]) + len2(cluster[j])) continue;

            int_t d = 0;

            if (strand == FORWARD) {
                // ref 维度也不重叠：i 的 start1 需在 j 的 ref 末尾之后
                int_t prev_endj = start1(cluster[j]) + len1(cluster[j]);
                if (start1(cluster[i]) <= prev_endj) continue;

                d = std::abs(diag(cluster[i]) - diag(cluster[j]));
            }
            else {
                // 反向链：按原逻辑计算不交叉条件与 sep
                int_t prev_endi = start1(cluster[i]) + len1(cluster[i]);
                if (prev_endi >= start1(cluster[j])) continue;

                d = std::abs(diag_reverse(cluster[i]) - diag_reverse(cluster[j]));
            }

            // 候选分数：前链分数 + 当前长度 - diagonal 差惩罚
            int_t cand = score[j] + len2(cluster[i]) - d;
            if (cand > score[i]) {
                score[i] = cand;
                pred[i] = static_cast<int_t>(j);
            }
        }

        if (score[i] > score[best_idx]) best_idx = i;
    }

    // 回溯构建最优链
    MatchVec chain;
    for (int_t k = static_cast<int_t>(best_idx); k != -1; k = pred[k])
        chain.emplace_back(cluster[k]);
    std::reverse(chain.begin(), chain.end());

    return chain;
}

/* ───────────────────────────────────────────────────────── *
 * clusterChrMatch：对外主函数
 * - 对输入 unique_match 聚簇
 * - 对每个簇做 bestChainDP
 * - 计算 span（这里按 match_len 求和）
 * - span >= min_cluster_length 则保留该 best_chain
 * ───────────────────────────────────────────────────────── */
MatchClusterVecPtr clusterChrMatch(MatchVec& unique_match,
                                  uint_t min_cluster_length,
                                  int_t  max_gap,
                                  int_t  diagdiff,
                                  double diagfactor) {

    auto best_chain_clusters = std::make_shared<MatchClusterVec>();

    if (unique_match.size() == 0) {
        return best_chain_clusters;
    }

    // 1) 聚簇
    MatchClusterVec clusters = buildClusters(unique_match, max_gap, diagdiff, diagfactor);

    // 释放 unique_match 的内存（保持原逻辑）
    unique_match.clear();
    unique_match.shrink_to_fit();

    best_chain_clusters->reserve(clusters.size());

    // 2) 每个簇选最佳链
    for (auto& cluster : clusters) {
        if (cluster.empty()) continue;

        MatchVec best_chain = bestChainDP(cluster, diagfactor);

        if (best_chain.empty()) {
            releaseCluster(cluster);
            continue;
        }

        uint_t span = 0;
        // 遍历 best_chain 累加 span（保持原逻辑）
        for (auto& m : best_chain) {
            span += m.match_len();
        }

        // 满足最小簇长度才保留
        if (span >= min_cluster_length) {
            best_chain_clusters->emplace_back(std::move(best_chain));
        }

        // 回收 cluster 剩余元素
        releaseCluster(cluster);
    }

    best_chain_clusters->shrink_to_fit();
    return best_chain_clusters;
}

/*------------------------------------------------------------------*
 * groupClustersByRefQuery：
 * 输入：
 *   by_ref  : shared_ptr<vector<MatchClusterVecPtr>>  (ref 维)
 * 输出：
 *   by_refQ : shared_ptr<vector< vector<MatchClusterVecPtr> > >
 *             └──ref──┘└────────query────────┘
 * 说明：
 * - 先确定 query 染色体数 query_cnt
 * - 预分配 [ref][query] 的矩阵结构
 * - 每个 ref 维度提交一个线程任务，把 cluster 依 qry_chr_index 分桶
 *------------------------------------------------------------------*/
// The former dense ref-by-query regrouping helper was unused. The active
// sparse clustering path lives in groupMatchByQueryRefSparse().
