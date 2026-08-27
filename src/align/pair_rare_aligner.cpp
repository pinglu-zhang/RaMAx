#include "rare_aligner.h"
#include "anchor.h"
#include "align.h"
#include <omp.h>
#include <atomic>
#include <exception>
#include <mutex>
#include <cstdlib>
#include <fstream>

PairRareAligner::PairRareAligner(const FilePath work_dir,
	const uint_t thread_num,
	uint_t chunk_size,
	uint_t overlap_size,
	uint_t min_anchor_length,
	uint_t max_anchor_frequency,
	uint_t accurate_skip_threshold,
	bool trust_legacy_cache_)
	: work_dir(work_dir)
	, index_dir(work_dir / INDEX_DIR)
	, chunk_size(chunk_size)
	, overlap_size(overlap_size)
	, min_anchor_length(min_anchor_length)
	, max_anchor_frequency(max_anchor_frequency)
	, accurate_skip_threshold(accurate_skip_threshold)
	, thread_num(thread_num)
	, trust_legacy_cache(trust_legacy_cache_)
{
	this->group_id = 0;
	this->round_id = 0;

}

// MatchVec3DPtr PairRareAligner::alignPairGenome(
// 	SpeciesName query_name,
// 	SeqPro::ManagerVariant& query_fasta_manager,
// 	SearchMode         search_mode,
// 	bool allow_MEM,
// 	bool allow_short_mum,
// 	sdsl::int_vector<0>& ref_global_cache,
// 	SeqPro::Length sampling_interval) {
//
//
// 	MatchVec3DPtr anchors = findQueryFileAnchor(
// 		query_name, query_fasta_manager, search_mode, allow_MEM, allow_short_mum,  ref_global_cache, sampling_interval);
//
// 	return anchors;
//
// }

FilePath PairRareAligner::buildIndex(const std::string prefix,
    SeqPro::ManagerVariant& ref_fasta_manager_, bool fast_build) {
    ref_name = prefix;
    ref_seqpro_manager = &ref_fasta_manager_;
    const FilePath logical_index_path = index_dir / prefix;

    ref_index.emplace(prefix, ref_fasta_manager_, sa_sampling_rate);
    spdlog::info(
        "Suffix-array indexing with prefix: {}, storage: memory-only, "
        "disk-cache: disabled, sampling rate: {}",
        prefix, sa_sampling_rate);
    if (!ref_index->buildIndex({}, fast_build, thread_num)) {
        throw std::runtime_error(
            "Failed to build in-memory suffix-array index for " + prefix);
    }
    ++index_cache_counters->memory_only_built;
    spdlog::info(
        "Suffix-array indexing finished in memory for {} (K={}, rows={}, "
        "disk-bytes=0)",
        prefix, ref_index->samplingRate(), ref_index->storedSuffixCount());
    return logical_index_path;
}

MatchVec3DPtr PairRareAligner::findQueryFileAnchor(
	const std::string prefix,
	SeqPro::ManagerVariant& query_fasta_manager,
	SearchMode         search_mode,
	bool allow_MEM,
	bool allow_short_mum,
	sdsl::int_vector<0>& ref_global_cache,
	SeqPro::Length sampling_interval,
	bool isMultiple,
	bool include_masked_regions)
{
	/* ---------- 1. 结果文件路径，与多基因组保持同一目录 ---------- */
	FilePath result_dir = work_dir / RESULT_DIR
		/ ("group_" + std::to_string(group_id))
		/ ("round_" + std::to_string(round_id));
	std::filesystem::create_directories(result_dir);

	spdlog::info("[findQueryFileAnchor] begin to algin {}", prefix);

	FilePath anchor_file = result_dir /
		(prefix + "_" + SearchModeToString(search_mode) + "." + ANCHOR_EXTENSION);

	/* ---------- 若已存在结果文件，直接加载 ---------- */
	//if (std::filesystem::exists(anchor_file)) {
	//	MatchVec3DPtr result = std::make_shared<MatchVec3D>();
	//	loadMatchVec3D(anchor_file, result);
	//	return result;
	//}
	// TODO 对于二轮之后，低于20的chunk可以不用输入
	/* ---------- 读取 FASTA 并分片 ---------- */
	// 修改：使用新的预分割逻辑，支持多基因组模式
	RegionVec chunks;
	if (include_masked_regions) {
		chunks = preAllocateOriginalChunksBySize(
			query_fasta_manager, chunk_size, overlap_size, 10000);
	} else if (isMultiple) {
		// 多基因组 primary pass：使用遮蔽区间预分割
		chunks = preAllocateChunksBySize(
			query_fasta_manager, chunk_size, overlap_size, 10000, true);
	} else {
		// 双基因组模式：使用普通分割
		chunks = preAllocateChunks(
			query_fasta_manager, chunk_size, overlap_size, 1000, 10000);
	}
	// 智能分块策略：自动根据序列数量和长度选择最优的分块方式
	/* ---------- ① 计时：搜索 Anchor ---------- */
	auto t_search0 = std::chrono::steady_clock::now();

	struct AnchorSearchTask {
		Region chunk;
		Strand strand;
	};

	std::vector<AnchorSearchTask> tasks;
	tasks.reserve(chunks.size() * 2);
	for (const auto& ck : chunks) {
		tasks.push_back({ck, Strand::FORWARD});
		tasks.push_back({ck, Strand::REVERSE});
	}

	std::vector<MatchVec2DPtr> task_results(tasks.size());
	std::atomic<size_t> completed_tasks{0};
	size_t next_progress = 1; // 1~20
	std::exception_ptr task_exception = nullptr;
	std::mutex task_exception_mutex;

	MatchVec3DPtr result = std::make_shared<MatchVec3D>();

	const size_t total = tasks.size();

#pragma omp parallel for schedule(dynamic) num_threads(thread_num)
	for (long long task_idx = 0; task_idx < static_cast<long long>(tasks.size()); ++task_idx) {
		try {
			const auto& task = tasks[static_cast<size_t>(task_idx)];
			const auto& ck = task.chunk;

			std::string seq = std::visit([&ck](auto&& manager_ptr) -> std::string {
				using PtrType = std::decay_t<decltype(manager_ptr)>;
				if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
					return manager_ptr->getSubSequence(ck.chr_index, ck.start, ck.length);
				} else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
					// 不再使用分隔符，因为chunks已经预分割了
					return manager_ptr->getOriginalManager().getSubSequence(ck.chr_index, ck.start, ck.length);
				} else {
					throw std::runtime_error("Unhandled manager type in variant.");
				}
			}, query_fasta_manager);

			if (seq.length() < ck.length) {
				task_results[static_cast<size_t>(task_idx)] = std::make_shared<MatchVec2D>();
			} else {
				const SearchMode task_search_mode =
					(task.strand == Strand::FORWARD) ? search_mode : ACCURATE_SEARCH;
				task_results[static_cast<size_t>(task_idx)] = ref_index->findAnchors(
					ck.chr_index, seq, task_search_mode,
					task.strand,
					allow_MEM,
					ck.start,
					min_anchor_length,
					allow_short_mum,
					max_anchor_frequency,
					ref_global_cache,
					sampling_interval,
					accurate_skip_threshold);
			}
		} catch (...) {
			std::lock_guard<std::mutex> lock(task_exception_mutex);
			if (!task_exception) {
				task_exception = std::current_exception();
			}
		}

		const size_t count = ++completed_tasks;
		size_t progress_stage = (count * 20) / total;
#pragma omp critical(find_query_file_anchor_progress)
		{
			if (progress_stage >= next_progress || count == total) {
				int percent = static_cast<int>((progress_stage * 100) / 20);
				spdlog::info("[{}] Progress: {}% ({} of {})", prefix, percent, count, total);
				next_progress = progress_stage + 1;
			}
		}
	}

	if (task_exception) {
		std::rethrow_exception(task_exception);
	}

	result->reserve(task_results.size());
	for (auto& part : task_results) {
		if (!part) {
			part = std::make_shared<MatchVec2D>();
		}
		result->emplace_back(std::move(*part));
	}


	auto t_search1 = std::chrono::steady_clock::now();
	double search_ms = std::chrono::duration<double, std::milli>(t_search1 - t_search0).count();

	/* ---------- ③ 计时：保存 ---------- */
	auto t_save0 = std::chrono::steady_clock::now();
	// saveMatchVec3D(anchor_file, result);
	auto t_save1 = std::chrono::steady_clock::now();
	double save_ms = std::chrono::duration<double, std::milli>(t_save1 - t_save0).count();

	/* ---------- Performance Statistics ---------- */
	spdlog::info("");
	spdlog::info("┌─────────────────────────────────────────────────────────┐");
	spdlog::info("│               findQueryFileAnchor Performance           │");
	spdlog::info("├─────────────────────────────────────────────────────────┤");
	spdlog::info("│  Search phase    : {:>8.3f} ms                          │", search_ms);
	spdlog::info("│  Save phase      : {:>8.3f} ms                          │", save_ms);
	spdlog::info("│  Total time      : {:>8.3f} ms                          │", search_ms + save_ms);
	spdlog::info("└─────────────────────────────────────────────────────────┘");
	spdlog::info("");

	return result;          // NRVO / move-elided
}

std::shared_ptr<PreparedAnchorSearch>
PairRareAligner::prepareQueryFileAnchor(
    SeqPro::ManagerVariant& query_fasta_manager,
    SearchMode search_mode,
    bool allow_mem,
    bool allow_short_mum,
    sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval,
    bool is_multiple,
    bool include_masked_regions) {
    auto plan =
        std::make_shared<PreparedAnchorSearch>();
    plan->query_manager =
        &query_fasta_manager;
    plan->search_mode = search_mode;
    plan->allow_mem = allow_mem;
    plan->allow_short_mum =
        allow_short_mum;
    plan->ref_global_cache =
        &ref_global_cache;
    plan->sampling_interval =
        sampling_interval;

    RegionVec chunks;
    if (include_masked_regions) {
        chunks =
            preAllocateOriginalChunksBySize(
                query_fasta_manager,
                chunk_size,
                overlap_size,
                10000);
    } else if (is_multiple) {
        chunks =
            preAllocateChunksBySize(
                query_fasta_manager,
                chunk_size,
                overlap_size,
                10000,
                true);
    } else {
        chunks =
            preAllocateChunks(
                query_fasta_manager,
                chunk_size,
                overlap_size,
                1000,
                10000);
    }

    plan->tasks.reserve(
        chunks.size() * 2);
    for (const auto& chunk : chunks) {
        plan->tasks.push_back(
            {chunk, Strand::FORWARD});
        plan->tasks.push_back(
            {chunk, Strand::REVERSE});
    }
    plan->task_results.resize(
        plan->tasks.size());
    plan->task_errors.resize(
        plan->tasks.size());
    return plan;
}

void PairRareAligner::
executePreparedAnchorTask(
    PreparedAnchorSearch& plan,
    size_t task_index) {
    if (task_index >=
        plan.tasks.size()) {
        throw std::out_of_range(
            "Prepared anchor task index is out of range");
    }
    try {
        if (plan.query_manager == nullptr ||
            plan.ref_global_cache == nullptr ||
            !ref_index) {
            throw std::logic_error(
                "Prepared anchor search is not fully initialized");
        }
        const auto& task =
            plan.tasks[task_index];
        const auto& chunk =
            task.chunk;
        std::string sequence =
            std::visit(
                [&](auto&& manager_ptr)
                    -> std::string {
                    using PtrType =
                        std::decay_t<
                            decltype(manager_ptr)>;
                    if constexpr (
                        std::is_same_v<
                            PtrType,
                            std::unique_ptr<
                                SeqPro::SequenceManager>>) {
                        return manager_ptr
                            ->getSubSequence(
                                chunk.chr_index,
                                chunk.start,
                                chunk.length);
                    } else {
                        return manager_ptr
                            ->getOriginalManager()
                            .getSubSequence(
                                chunk.chr_index,
                                chunk.start,
                                chunk.length);
                    }
                },
                *plan.query_manager);

        if (sequence.length() <
            chunk.length) {
            plan.task_results[task_index] =
                std::make_shared<
                    MatchVec2D>();
            return;
        }
        const SearchMode task_search_mode =
            task.strand ==
                    Strand::FORWARD
                ? plan.search_mode
                : ACCURATE_SEARCH;
        plan.task_results[task_index] =
            ref_index->findAnchors(
                chunk.chr_index,
                sequence,
                task_search_mode,
                task.strand,
                plan.allow_mem,
                chunk.start,
                min_anchor_length,
                plan.allow_short_mum,
                max_anchor_frequency,
                *plan.ref_global_cache,
                plan.sampling_interval);
    } catch (...) {
        plan.task_errors[task_index] =
            std::current_exception();
    }
}

MatchVec3DPtr PairRareAligner::
collectPreparedAnchorSearch(
    PreparedAnchorSearch& plan) {
    for (const auto& error :
         plan.task_errors) {
        if (error) {
            std::rethrow_exception(
                error);
        }
    }
    auto result =
        std::make_shared<MatchVec3D>();
    result->reserve(
        plan.task_results.size());
    for (auto& task_result :
         plan.task_results) {
        if (!task_result) {
            task_result =
                std::make_shared<
                    MatchVec2D>();
        }
        result->emplace_back(
            std::move(*task_result));
    }
    return result;
}

// ========== 工具：快速取子串 ==========
static std::string subSeq(const SeqPro::ManagerVariant& mv,
                                 ChrIndex chr, int_t b, int_t l) {
    return std::visit([&](auto& p) -> std::string {
        using T = std::decay_t<decltype(p)>;
        if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::SequenceManager>>) {
            return p->getSubSequence(chr, b, l);
        } else {
            return p->getOriginalManager().getSubSequence(chr, b, l);
        }
    }, mv);
}

// ========== 工具：判断 cluster 是否被某个 anchor 完全“罩住”（shadow） ==========
static bool shadowedBy(const Anchor& a, const MatchCluster& cl) {
    if (cl.empty()) return true;
    const auto& first = cl.front();
    const auto& last  = cl.back();

    if (a.ref_chr_index != first.ref_chr_index
        || a.qry_chr_index != first.qry_chr_index
        || a.strand != first.strand()) return false;

    int_t c_ref_beg = start1(first);
    int_t c_ref_end = start1(last) + len1(last);

    int_t c_qry_beg_fwd = start2(first);
    int_t c_qry_end_fwd = start2(last) + len2(last);

    if (a.strand == FORWARD) {
        return a.ref_start <= c_ref_beg && (a.ref_start + a.ref_len) >= c_ref_end
            && a.qry_start <= c_qry_beg_fwd && (a.qry_start + a.qry_len) >= c_qry_end_fwd;
    } else {
        // 反链：anchor.qry_start 是“较小坐标”的起点
        int_t c_qry_low  = start2(last);
        int_t c_qry_high = start2(first) + len2(first);
        return a.ref_start <= c_ref_beg && (a.ref_start + a.ref_len) >= c_ref_end
            && a.qry_start <= c_qry_low && (a.qry_start + a.qry_len) >= c_qry_high;
    }
}

// 两端 exact match 已固定边界；这里必须优化完整区间，而不是做 ends-free extension。
static bool extend_gap_must_reach(const std::string& ref_gap,
                                  const std::string& qry_gap,
                                  Cigar_t& out) {
    out = globalAlignKSW2(ref_gap, qry_gap);
    if (countRefLength(out) != ref_gap.size() ||
        countQryLength(out) != qry_gap.size()) {
        return false;
    }
    return alignmentCigarPreferredToUnaligned(
        ref_gap, qry_gap, out);
}

// 判定 Match 与现有 Anchor 是否冲突
// 返回 true 表示“冲突/不能直接作为下一个目标 match”。
inline bool isMatchConflictingWithAnchor(
	const Match& m,
	const Anchor& a
) {

	// 2) 计算 ref gap（新 match 的 ref 起点应在 anchor 末端之后或相等）
	const Coord_t ref_end = a.ref_start + a.ref_len; // anchor 的 ref 前向端
	const int_t   ref_gap = static_cast<int_t>(m.ref_start) - static_cast<int_t>(ref_end);
	if (ref_gap < 0) return true;                    // 倒退/重叠 → 冲突


	// 3) 计算 qry gap（按链向不同，前向端定义不同）
	int_t qry_gap = 0;
	if (a.strand == FORWARD) {
		const Coord_t qry_end = a.qry_start + a.qry_len; // anchor 的 qry 前向端（大坐标端）
		qry_gap = static_cast<int_t>(m.qry_start) - static_cast<int_t>(qry_end);
		if (qry_gap < 0) return true;

	} else { // REVERSE
		// 反向链上，anchor 的“前向端”在小坐标端 = a.qry_start
		// 新 match 的“前向端”取 (sB + len)
		const Coord_t m_q_front = m.qry_start + len2(m);
		qry_gap = static_cast<int_t>(a.qry_start) - static_cast<int_t>(m_q_front);
		if (qry_gap < 0) return true;
	}

	// 通过全部检查：不冲突
	return false;
}


AnchorPtrVec extendClustersToAnchors(
    const MatchClusterVecPtr& cluster_vec_ptr,
    const SeqPro::ManagerVariant& ref_mgr,
    const SeqPro::ManagerVariant& qry_mgr) {
    AnchorPtrVec anchors;
    if (!cluster_vec_ptr || cluster_vec_ptr->empty()) {
        return anchors;
    }

    auto clusters = *cluster_vec_ptr;
    for (auto& cluster : clusters) {
        std::sort(
            cluster.begin(), cluster.end(),
            [](const Match& left, const Match& right) {
                return left.ref_start < right.ref_start;
            });
    }
    std::sort(
        clusters.begin(), clusters.end(),
        [](const MatchCluster& left, const MatchCluster& right) {
            if (left.empty() || right.empty()) {
                return left.size() < right.size();
            }
            return start1(left.front()) < start1(right.front());
        });

    const auto make_anchor = [](const Match& match) {
        auto anchor = std::make_shared<Anchor>();
        anchor->ref_chr_index = match.ref_chr_index;
        anchor->qry_chr_index = match.qry_chr_index;
        anchor->strand = match.strand();
        anchor->ref_start = start1(match);
        anchor->qry_start = start2(match);
        anchor->ref_len = len1(match);
        anchor->qry_len = len2(match);
        anchor->aligned_base = len1(match);
        anchor->alignment_length = len1(match);
        appendCigarOp(anchor->cigar, 'M', len1(match));
        return anchor;
    };
    const auto flush_anchor = [&anchors](AnchorPtr& anchor) {
        if (anchor) {
            anchors.push_back(std::move(anchor));
        }
    };

    constexpr int_t kMaximumInternalGap = 10000;
    for (auto& cluster : clusters) {
        if (cluster.empty()) {
            continue;
        }
        const bool shadowed = std::any_of(
            anchors.begin(), anchors.end(),
            [&cluster](const AnchorPtr& anchor) {
                return anchor && shadowedBy(*anchor, cluster);
            });
        if (shadowed) {
            continue;
        }

        AnchorPtr current;
        for (const Match& match : cluster) {
            if (len1(match) == 0 || len1(match) != len2(match)) {
                continue;
            }
            if (!current) {
                current = make_anchor(match);
                continue;
            }
            if (match.ref_chr_index != current->ref_chr_index ||
                match.qry_chr_index != current->qry_chr_index ||
                match.strand() != current->strand) {
                flush_anchor(current);
                current = make_anchor(match);
                continue;
            }
            if (isMatchConflictingWithAnchor(match, *current)) {
                continue;
            }

            const int_t ref_gap =
                static_cast<int_t>(start1(match)) -
                static_cast<int_t>(
                    current->ref_start + current->ref_len);
            const bool forward = current->strand == FORWARD;
            const int_t qry_gap =
                forward
                    ? static_cast<int_t>(start2(match)) -
                          static_cast<int_t>(
                              current->qry_start + current->qry_len)
                    : static_cast<int_t>(current->qry_start) -
                          static_cast<int_t>(
                              start2(match) + len2(match));

            Cigar_t gap_cigar;
            bool bridged = false;
            if (ref_gap >= 0 && qry_gap >= 0 &&
                ref_gap < kMaximumInternalGap &&
                qry_gap < kMaximumInternalGap) {
                if (ref_gap == 0 && qry_gap == 0) {
                    bridged = true;
                } else if (ref_gap == 0) {
                    appendCigarOp(
                        gap_cigar, 'I',
                        static_cast<uint32_t>(qry_gap));
                    bridged = true;
                } else if (qry_gap == 0) {
                    appendCigarOp(
                        gap_cigar, 'D',
                        static_cast<uint32_t>(ref_gap));
                    bridged = true;
                } else {
                    const int_t ref_gap_start =
                        current->ref_start + current->ref_len;
                    const int_t qry_gap_start =
                        forward
                            ? current->qry_start + current->qry_len
                            : start2(match) + len2(match);
                    std::string ref_sequence = subSeq(
                        ref_mgr, current->ref_chr_index,
                        ref_gap_start, ref_gap);
                    std::string qry_sequence = subSeq(
                        qry_mgr, current->qry_chr_index,
                        qry_gap_start, qry_gap);
                    if (!forward) {
                        reverseComplement(qry_sequence);
                    }
                    bridged = extend_gap_must_reach(
                        ref_sequence, qry_sequence, gap_cigar);
                }
            }

            if (!bridged) {
                flush_anchor(current);
                current = make_anchor(match);
                continue;
            }

            appendCigar(current->cigar, gap_cigar);
            current->alignment_length +=
                countAlignmentLength(gap_cigar);
            current->aligned_base +=
                countMatchOperations(gap_cigar);
            current->ref_len += ref_gap;
            current->qry_len += qry_gap;
            if (!forward) {
                current->qry_start -= qry_gap;
            }

            appendCigarOp(current->cigar, 'M', len1(match));
            current->aligned_base += len1(match);
            current->alignment_length += len1(match);
            current->ref_len += len1(match);
            current->qry_len += len2(match);
            if (!forward) {
                current->qry_start -= len2(match);
            }
        }
        flush_anchor(current);
    }
    return anchors;
}


AnchorPtrVec PairRareAligner::
extendClusterGroupToAnchors(
    SeqPro::ManagerVariant&
        query_seqpro_manager,
    MatchClusterVec& cluster_group,
    bool is_first) {
    AnchorPtrVec anchors;
    if (cluster_group.empty()) {
        return anchors;
    }
    if (!is_first) {
        for (auto& cluster :
             cluster_group) {
            for (auto& match : cluster) {
                MatchVec single_match{
                    match};
                Anchor anchor =
                    extendClusterToAnchor(
                        single_match,
                        *ref_seqpro_manager,
                        query_seqpro_manager);
                anchors.push_back(
                    std::make_shared<Anchor>(
                        std::move(anchor)));
            }
            MatchCluster().swap(cluster);
        }
    } else {
        anchors = linkClusters(
            cluster_group,
            *ref_seqpro_manager,
            query_seqpro_manager);
    }
    MatchClusterVec().swap(cluster_group);
    return anchors;
}


AnchorBySQR_SparsePtr PairRareAligner::extendClusterToAnchorByChr(SpeciesName query_name, SeqPro::ManagerVariant& query_seqpro_manager, ClusterBySQR_SparsePtr cluster, bool is_first)
{
    // 输出结构
    AnchorBySQR_SparsePtr result = std::make_shared<AnchorBySQR_Sparse>();


	const uint_t cluster_num = static_cast<uint_t>(cluster->size());
	std::vector<MatchClusterVecPtr> tmp_cluster_vec;
	AnchorBySQR_Sparse tmp_res(cluster_num);

	for (auto& [key, c_p] : *cluster)
	{
		tmp_cluster_vec.emplace_back(c_p);
	}


    std::atomic<size_t> completed{0};
    std::atomic<int> next_milestone{10}; // 10%,20%...

    #pragma omp parallel for schedule(dynamic) num_threads(thread_num)
    for (long long t = 0; t < static_cast<long long>(cluster_num); ++t) {
	MatchClusterVecPtr tmp_p = tmp_cluster_vec[t];
	if (!tmp_p || tmp_p->empty()) continue;
        AnchorPtrVec anchors =
            extendClusterGroupToAnchors(
                query_seqpro_manager,
                *tmp_p,
                is_first);

        if (!anchors.empty()) {
            tmp_res[t] = std::move(anchors);
        }

        // ===== 进度条2：10% 里程碑打印 =====
        size_t done = completed.fetch_add(1, std::memory_order_relaxed) + 1;

        // 避免 total=0 的极端情况
        if (cluster_num > 0) {
            int pct = static_cast<int>((done * 100) / cluster_num);
            int milestone = next_milestone.load(std::memory_order_relaxed);

            if (pct >= milestone && milestone <= 100) {
                #pragma omp critical
                {
                    milestone = next_milestone.load(std::memory_order_relaxed);
                    if (pct >= milestone && milestone <= 100) {
                        spdlog::info(
                            "extend cluster progress for {}: {}% ({}/{})",
                            query_name,
                            milestone,
                            done,
                            cluster_num
                        );
                        next_milestone.store(milestone + 10, std::memory_order_relaxed);
                    }
                }
            }
        }
    }

	for (uint_t i = 0; i < cluster_num; ++i)
	{
		if (tmp_res[i].size() > 0)
		{
			result->emplace_back(tmp_res[i]);
		}
	}

    spdlog::info("extend cluster to anchor successfully for {}", query_name);
    return result;
}




static void filterAnchorsByDP(AnchorPtrVec result, bool filter_ref)
{
	if (result.empty()) return;

	// ====== DP 部分 ======

	// 1. 排序
	std::sort(result.begin(), result.end(),
		[filter_ref](const AnchorPtr& a, const AnchorPtr& b) {
			return filter_ref ? (a->ref_start < b->ref_start)
				: (a->qry_start < b->qry_start);
		});

	//size_t n = result.size();
	std::vector<double> dp(result.size(), 0);
	std::vector<int_t> pre(result.size(), -1);

	constexpr double OVERLAP_MAX_FRAC = 0; // allow overlap up to 10% of the SHORTER interval
	constexpr double OVERLAP_PENALTY_W = 1.0; // per-base penalty weight to subtract when overlapping

	size_t K = 5000; // 可调参数
	auto interval = [&](size_t idx) -> std::pair<long long, long long> {
		if (filter_ref) {
			return { static_cast<long long>(result[idx]->ref_start), static_cast<long long>(result[idx]->ref_len) };
		} else {
			return { static_cast<long long>(result[idx]->qry_start), static_cast<long long>(result[idx]->qry_len) };
		}
	};
	///
	for (size_t i = 0; i < result.size(); ++i) {
		double idy = static_cast<float>(result[i]->aligned_base) / result[i]->alignment_length;
		double score_i = result[i]->alignment_length * pow(idy, 2);
		// if (score_i > 1000000 || result[i]->cigar.size() == 0) {
		// 	continue;
		// }
		// if (filter_ref == true && id != result[i]->ref_chr_index)
		// {
		// 	continue;
		// }
		// if (filter_ref == false && id != result[i]->qry_chr_index)
		// {
		// 	continue;
		// }
		dp[i] = score_i;


		// 只回看最近 K 个 j
		size_t j_start = (i > K ? i - K : 0);
		//size_t j_start = 0;
		for (size_t j = j_start; j < i; ++j) {
			auto [sj, lj] = interval(j);
			auto [si, li] = interval(i);
			long long ej = sj + lj;
			long long ei = si + li;


			// compute overlap length in [start, end)
			long long overlap = std::max(0LL, std::min(ej, ei) - std::max(sj, si));
			long long short_len = std::min(lj, li);
			double overlap_ratio = (short_len > 0) ? static_cast<double>(overlap) / static_cast<double>(short_len) : 0.0;


			// Accept if no overlap OR small (<=10%) overlap, otherwise skip this predecessor j
			if (overlap_ratio <= OVERLAP_MAX_FRAC) {
				// penalty proportional to the absolute overlap length
				double penalty = OVERLAP_PENALTY_W * static_cast<double>(overlap);
				double candidate = dp[j] + score_i - penalty;
				if (candidate > dp[i]) {
					dp[i] = candidate;
					pre[i] = j;
				}
			}
		}
		
	}
	// 3. 找最大值
	uint_t best = 0;
	size_t best_idx = 0;
	for (size_t i = 0; i < result.size(); ++i) {
		if (dp[i] > best) {
			best = dp[i];
			best_idx = i;
		}
	}

	// 4. 回溯并标记
	for (int i = static_cast<int>(best_idx); i >= 0; i = pre[i]) {
		if (filter_ref)
			result[i]->ref_selected = true;
		else
			result[i]->qry_selected = true;
		if (pre[i] == -1) break;
	}

	// ----------------------------------------------------
// 检查选中的比对是否存在重叠
// ----------------------------------------------------
#ifdef _DEBUG_
	bool hasOverlap = false;

	// 收集所有选中的 index
	std::vector<int> selected;
	for (int i = 0; i < (int)result.size(); ++i) {
		if ((filter_ref && result[i]->ref_selected) ||
			(!filter_ref && result[i]->qry_selected))
			selected.push_back(i);
	}

	// 按 ref_start 或 qry_start 排序
	std::sort(selected.begin(), selected.end(),
		[&](int a, int b) {
			return filter_ref
				? (result[a]->ref_start < result[b]->ref_start)
				: (result[a]->qry_start < result[b]->qry_start);
		});


	// 检查相邻区间是否重叠
	for (size_t i = 1; i < selected.size(); ++i) {
		auto& prev = result[selected[i - 1]];
		auto& curr = result[selected[i]];

		if (filter_ref) {
			if (curr->ref_start < prev->ref_start + prev->ref_len) {
				hasOverlap = true;
				std::cerr << "[Overlap] Ref overlap between "
					<< selected[i - 1] << " and " << selected[i]
					<< " (" << prev->ref_start << "-" << prev->ref_start + prev->ref_len
						<< " vs " << curr->ref_start << "-" << curr->ref_start + curr->ref_len
						<< ")\n";
			}
		}
		else {
			if (curr->qry_start < prev->qry_start + prev->qry_len) {
				hasOverlap = true;
				std::cerr << "[Overlap] Qry overlap between "
					<< selected[i - 1] << " and " << selected[i]
					<< " (" << prev->qry_start << "-" << prev->qry_start + prev->qry_len
						<< " vs " << curr->qry_start << "-" << curr->qry_start + curr->qry_len
						<< ")\n";
			}
		}
	}

	if (!hasOverlap)
		spdlog::debug("No overlaps detected in final selection.");
	else
		spdlog::debug("Overlaps found in final selection!");
#endif
}

static void filterChrByDP(
	AnchorBySQR_SparsePtr anchor_map,
	uint_t id,
	bool filter_ref)
{
	AnchorPtrVec result;
	if (!anchor_map) return;
	for (auto& a_vec : *anchor_map) {
		if (a_vec.empty()) continue;
		const uint_t chromosome_id = filter_ref
			? a_vec.front()->ref_chr_index
			: a_vec.front()->qry_chr_index;
		if (chromosome_id == id) {
			result.insert(result.end(), a_vec.begin(), a_vec.end());
		}
	}
	filterAnchorsByDP(std::move(result), filter_ref);
}

void PairRareAligner::
filterAnchorByDPDimension(
    AnchorBySQR_SparsePtr anchor_map,
    uint_t chromosome_id,
    bool filter_ref) {
    filterChrByDP(
        std::move(anchor_map),
        chromosome_id,
        filter_ref);
}

void PairRareAligner::filterAnchorVectorByDP(
    AnchorPtrVec anchors, bool filter_ref) {
    filterAnchorsByDP(std::move(anchors), filter_ref);
}

void PairRareAligner::filterAnchorByDP(AnchorBySQR_SparsePtr anchor_map, uint_t ref_chr_cnt, uint_t qry_chr_cnt) {

#pragma omp parallel for schedule(dynamic) num_threads(thread_num)
	for (uint_t i = 0; i < ref_chr_cnt; i++) {
		filterAnchorByDPDimension(anchor_map, i, true);
	}

#pragma omp parallel for schedule(dynamic) num_threads(thread_num)
	for (uint_t i = 0; i < qry_chr_cnt; i++) {
		filterAnchorByDPDimension(anchor_map, i, false);

	}

}

void PairRareAligner::constructGraphByDP(
    SpeciesName query_name,
    SeqPro::ManagerVariant& query_seqpro_manager,
    AnchorBySQR_SparsePtr anchor_ptr,
    RaMesh::RaMeshMultiGenomeGraph& graph) {
    for (auto& anchor_group : *anchor_ptr) {
        for (const auto& anchor : anchor_group) {
            if (!anchor->ref_selected ||
                !anchor->qry_selected) {
                continue;
            }
            try {
                graph.insertAnchorIntoGraph(
                    *ref_seqpro_manager,
                    query_seqpro_manager,
                    ref_name,
                    query_name,
                    *anchor,
                    false);
            } catch (const std::exception& error) {
                spdlog::error(
                    "Error inserting anchor into graph: {}",
                    error.what());
            }
        }
        AnchorPtrVec().swap(anchor_group);
    }
}

void PairRareAligner::registerSecondaryAnchors(
    SpeciesName query_name,
    SeqPro::ManagerVariant& query_seqpro_manager,
    AnchorBySQR_SparsePtr anchor_ptr,
    RaMesh::RaMeshMultiGenomeGraph& graph,
    bool initial_round) {
    if (!anchor_ptr) {
        return;
    }
    const auto chromosome_name = [](
        SeqPro::ManagerVariant& manager,
        const ChrIndex chromosome_index) -> ChrName {
        return std::visit(
            [&](auto& manager_ptr) -> ChrName {
                using ManagerPtr = std::decay_t<decltype(manager_ptr)>;
                if constexpr (std::is_same_v<
                                  ManagerPtr,
                                  std::unique_ptr<
                                      SeqPro::SequenceManager>>) {
                    return manager_ptr->getSequenceName(
                        chromosome_index);
                } else {
                    return manager_ptr->getOriginalManager()
                        .getSequenceName(chromosome_index);
                }
            },
            manager);
    };
    for (const auto& anchor_group : *anchor_ptr) {
        for (const auto& anchor : anchor_group) {
            const ChrName ref_chromosome = chromosome_name(
                *ref_seqpro_manager, anchor->ref_chr_index);
            const ChrName query_chromosome = chromosome_name(
                query_seqpro_manager, anchor->qry_chr_index);
            graph.registerSecondaryAnchorCandidate(
                ref_name, ref_chromosome, query_name,
                query_chromosome, *anchor, initial_round, true);
            graph.registerSecondaryAnchorCandidate(
                ref_name, ref_chromosome, query_name,
                query_chromosome, *anchor, initial_round, false);
        }
    }
}
