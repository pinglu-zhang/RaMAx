#include "rare_aligner.h"
#include "anchor.h"
#include "align.h"
#include <omp.h>
#include <atomic>

PairRareAligner::PairRareAligner(const FilePath work_dir,
	const uint_t thread_num,
	uint_t chunk_size,
	uint_t overlap_size,
	uint_t min_anchor_length,
	uint_t max_anchor_frequency)
	: work_dir(work_dir)
	, index_dir(work_dir / INDEX_DIR)
	, thread_num(thread_num)
	, chunk_size(chunk_size)
	, overlap_size(overlap_size)
	, min_anchor_length(min_anchor_length)
	, max_anchor_frequency(max_anchor_frequency)
{
	if (!std::filesystem::exists(index_dir)) {
		std::filesystem::create_directories(index_dir);
	}


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

FilePath PairRareAligner::buildIndex(const std::string prefix, SeqPro::ManagerVariant& ref_fasta_manager_, bool fast_build) {

	ref_name = prefix;
	FilePath ref_index_path = index_dir / prefix;

	if (!std::filesystem::exists(ref_index_path)) {
		std::filesystem::create_directories(ref_index_path);
	}

	ref_seqpro_manager = &ref_fasta_manager_;
	// 使用 std::visit 来获取 fasta_path
	std::string fasta_path_str;
	std::visit([&fasta_path_str](auto&& manager_ptr) {
		using PtrType = std::decay_t<decltype(manager_ptr)>;
		if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
			// 如果是 SequenceManager，直接调用 getFastaPath()
			fasta_path_str = manager_ptr->getFastaPath();
		} else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
			// 如果是 MaskedSequenceManager，先 getOriginalManager() 再 getFastaPath()
			fasta_path_str = manager_ptr->getOriginalManager().getFastaPath();
		} else {
			throw std::runtime_error("Unhandled manager type in variant.");
		}
	}, ref_fasta_manager_);

	// index_dir的路径加上prefix的前缀加上fasta_manager.fasta_path_的扩展名
	FilePath output_path = ref_index_path / (prefix + std::filesystem::path(fasta_path_str).extension().string());

	ref_index.emplace(prefix, ref_fasta_manager_);

	FilePath idx_file_path = ref_index_path / (prefix + "." + FMINDEX_EXTESION);

	spdlog::info("Indexing with prefix: {}, index path: {}", prefix, ref_index_path.string());

	if (!std::filesystem::exists(idx_file_path)) {
		ref_index->buildIndex(output_path, fast_build, thread_num);
		ref_index->saveToFile(idx_file_path.string());
	}
	else {
		ref_index->loadFromFile(idx_file_path.string());
	}
	spdlog::info("Indexing finished, index path: {}", ref_index_path.string());

	return ref_index_path;
}

MatchVec3DPtr PairRareAligner::findQueryFileAnchor(
	const std::string prefix,
	SeqPro::ManagerVariant& query_fasta_manager,
	SearchMode         search_mode,
	bool allow_MEM,
	bool allow_short_mum,
	ThreadPool& pool,
	sdsl::int_vector<0>& ref_global_cache,
	SeqPro::Length sampling_interval,
	bool isMultiple)
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
	if (isMultiple) {
		// 多基因组模式：使用遮蔽区间预分割
		chunks = preAllocateChunksBySize(query_fasta_manager, chunk_size, overlap_size, 10000, true);
	} else {
		// 双基因组模式：使用普通分割
		chunks = preAllocateChunks(query_fasta_manager, chunk_size, overlap_size, 1000, 10000);
	}
	// 智能分块策略：自动根据序列数量和长度选择最优的分块方式
	/* ---------- ① 计时：搜索 Anchor ---------- */
	auto t_search0 = std::chrono::steady_clock::now();

	// ThreadPool pool(thread_num);
	std::vector<std::future<MatchVec2DPtr>> futures;
	// 根据线程数量和chunk数量决定每个线程处理的chunk数量
	// 确保每个线程至少处理一个chunk，同时避免线程过多
	size_t num_chunks_per_thread = chunks.size() / thread_num;
	if (chunks.size() % thread_num != 0) {
		num_chunks_per_thread++; // 如果不能整除，则向上取整，确保所有chunk都被处理
	}
	if (num_chunks_per_thread == 0 && !chunks.empty()) { // 至少处理一个chunk
	    num_chunks_per_thread = 1;
	}

	futures.reserve(chunks.size() / num_chunks_per_thread + (chunks.size() % num_chunks_per_thread != 0 ? 1 : 0));

	for (size_t i = 0; i < chunks.size(); i += num_chunks_per_thread) {
		std::vector<Region> chunk_group;
		for (size_t j = i; j < std::min(i + num_chunks_per_thread, chunks.size()); ++j) {
			chunk_group.push_back(chunks[j]);
		}

		futures.emplace_back(
			pool.enqueue(
				[this, chunk_group, &query_fasta_manager, search_mode, allow_MEM, allow_short_mum, &ref_global_cache, sampling_interval, isMultiple]() -> MatchVec2DPtr {
					MatchVec2DPtr group_matches = std::make_shared<MatchVec2D>();
					for (const auto& ck : chunk_group) {
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
						if (seq.length() <ck.length) continue;
						MatchVec2DPtr forwoard_matches = ref_index->findAnchors(
							ck.chr_index, seq, search_mode,
							Strand::FORWARD,
							allow_MEM,
							ck.start,
							min_anchor_length,
							allow_short_mum,
							max_anchor_frequency,
							ref_global_cache,
							sampling_interval);

						// 合并当前chunk的matches到group_matches
						for (const auto& match_list : *forwoard_matches) {
							group_matches->push_back(match_list);
						}
					}
					return group_matches;
				}));
		futures.emplace_back(
			pool.enqueue(
				[this, chunk_group, &query_fasta_manager, search_mode, allow_MEM,allow_short_mum, &ref_global_cache, sampling_interval, isMultiple]() -> MatchVec2DPtr {
					MatchVec2DPtr group_matches = std::make_shared<MatchVec2D>();
					for (const auto& ck : chunk_group) {
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
						if (seq.length() <ck.length) continue;
						MatchVec2DPtr reverse_matches = ref_index->findAnchors(
							ck.chr_index, seq, ACCURATE_SEARCH,
							Strand::REVERSE,
							allow_MEM,
							ck.start,
							min_anchor_length,
							allow_short_mum,
							max_anchor_frequency,
							ref_global_cache,
							sampling_interval);
						// 合并当前chunk的matches到group_matches
						for (const auto& match_list : *reverse_matches) {
							group_matches->push_back(match_list);
						}
					}
					return group_matches;
				}));

	}

	MatchVec3DPtr result = std::make_shared<MatchVec3D>();

	result->reserve(futures.size());
	size_t total = futures.size();
	size_t count = 0;
	size_t next_progress = 1; // 1~20

	for (auto& fut : futures) {
		MatchVec2DPtr part = fut.get();
		result->emplace_back(std::move(*part));
		++count;
		size_t progress_stage = (count * 20) / total;
		if (progress_stage >= next_progress || count == total) {
			int percent = static_cast<int>((progress_stage * 100) / 20);
			spdlog::info("[{}] Progress: {}% ({} of {})", prefix, percent, count, total);
			next_progress = progress_stage + 1;
		}
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


/* ============================================================= *
 *  把三维 clusters  ->  按 ref 的一维 clusters
 *  并行对每个 ref 走 keepWithSplitGreedy，再写回三维结构
 * ============================================================= */
 /**
  * @brief  全局 MatchClusterVec 上执行“两级 map + 最大堆贪婪拆分”过滤
  *
  * @param cluster_vec_ptr  所有 clusters 的 shared_ptr
  * @param pool             ThreadPool（本实现单线程，参数仅留作占位）
  * @param min_span         最小跨度阈值
  */
void PairRareAligner::constructGraphByGreedy(SpeciesName query_name, SeqPro::ManagerVariant& query_seqpro_manager, ClusterVecPtrByStrandByQueryRefPtr cluster_ptr, RaMesh::RaMeshMultiGenomeGraph& graph, uint_t min_span)
{
#ifdef _DEBUG_
	using Clock = std::chrono::steady_clock;
	using nsec_t = uint64_t;
	static std::atomic<nsec_t> ns_producer{ 0 }, ns_extend{ 0 }, ns_insert{ 0 };
	static std::atomic<uint64_t> cnt_submit{ 0 }, cnt_extend{ 0 }, cnt_insert{ 0 }, cnt_finish{ 0 };

	static std::atomic<nsec_t> ns_check_overlap{ 0 }, ns_insert_interval{ 0 },
		ns_clone_cluster{ 0 }, ns_enqueue_task{ 0 }, ns_store_kept{ 0 };

	auto ns2ms = [](nsec_t ns) { return static_cast<double>(ns) / 1'000'000.0; };
#endif

	ThreadPool pool(thread_num);

	MatchClusterVecPtr cluster_vec_ptr = groupClustersToVec(cluster_ptr, pool, thread_num);

	pool.waitAllTasksDone();

	if (!cluster_vec_ptr || cluster_vec_ptr->empty()) return;

	struct Node {
		MatchCluster cl;
		int_t        span;
	};
	auto cmp = [](const Node& a, const Node& b) { return a.span < b.span; };

	std::vector<Node> heap;
	heap.reserve(cluster_vec_ptr->size());

	for (auto& cl : *cluster_vec_ptr) {
		if (cl.empty()) continue;
		int_t sc = clusterSpan(cl);
		if (sc >= min_span)
			heap.push_back({ std::move(cl), sc });
	}
	cluster_vec_ptr->clear();
	std::make_heap(heap.begin(), heap.end(), cmp);

	std::vector<IntervalMap> rMaps;
	std::vector<IntervalMap> qMaps;

	MatchClusterVec kept;
	kept.reserve(heap.size());

	while (!heap.empty()) {
#ifdef _DEBUG_
		auto t0 = Clock::now();
#endif

		std::pop_heap(heap.begin(), heap.end(), cmp);
		Node cur = std::move(heap.back());
		heap.pop_back();

		if (cur.span < min_span || cur.cl.empty()) continue;

		const ChrIndex refChr = cur.cl.front().ref_chr_index;
		const ChrIndex qChr = cur.cl.front().qry_chr_index;

		Strand strand = cur.cl.front().strand();
		uint_t rb = start1(cur.cl.front());
		uint_t re = start1(cur.cl.back()) + len1(cur.cl.back());
		uint_t qb = 0;
		uint_t qe = 0;
		if (strand == FORWARD) {
			qb = start2(cur.cl.front());
			qe = start2(cur.cl.back()) + len2(cur.cl.back());
		}
		else {
			qb = start2(cur.cl.back());
			qe = start2(cur.cl.front()) + len2(cur.cl.front());
		}

		int_t RL = 0, RR = 0, QL = 0, QR = 0;

#ifdef _DEBUG_
		auto t_chk = Clock::now();
#endif
		bool ref_hit = overlap1D(rMaps[refChr], rb, re, RL, RR);
		bool query_hit = overlap1D(qMaps[qChr], qb, qe, QL, QR);
#ifdef _DEBUG_
		ns_check_overlap += std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - t_chk).count();
#endif

		if (!ref_hit && !query_hit) {
#ifdef _DEBUG_
			auto t_ins = Clock::now();
#endif
			insertInterval(rMaps[refChr], rb, re);
			insertInterval(qMaps[qChr], qb, qe);
#ifdef _DEBUG_
			ns_insert_interval += std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - t_ins).count();
#endif

#ifdef _DEBUG_
			auto t_clone = Clock::now();
#endif
			auto task_cl = std::make_shared<MatchCluster>(cur.cl);
#ifdef _DEBUG_
			ns_clone_cluster += std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - t_clone).count();
#endif

#ifdef _DEBUG_
			auto t_enq = Clock::now();
#endif
			pool.enqueue([this, &graph, &query_name, task_cl, &query_seqpro_manager] {
				AnchorVec anchor_vec = extendClusterToAnchorVec(*task_cl, *ref_seqpro_manager, query_seqpro_manager);
				for (auto& anchor : anchor_vec) {
					graph.insertAnchorIntoGraph(*ref_seqpro_manager, query_seqpro_manager, ref_name, query_name, anchor, true);
				}	
				
				// graph.insertClusterIntoGraph(ref_name, query_name, *task_cl);
				});
#ifdef _DEBUG_
			ns_enqueue_task += std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - t_enq).count();
#endif

#ifdef _DEBUG_
			auto t_store = Clock::now();
#endif
			kept.emplace_back(cur.cl);
#ifdef _DEBUG_
			ns_store_kept += std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - t_store).count();
#endif

#ifdef _DEBUG_
			ns_producer += std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - t0).count();
			cnt_submit.fetch_add(1, std::memory_order_relaxed);
#endif

			continue;
		}

		for (auto& part : splitCluster(cur.cl, ref_hit, RL, RR, query_hit, QL, QR)) {
			int_t sc = clusterSpan(part);
			if (sc >= min_span) {
				heap.push_back({ part, sc });
				std::push_heap(heap.begin(), heap.end(), cmp);
			}
		}
	}

	pool.waitAllTasksDone();

#ifdef _DEBUG_
	spdlog::debug("");
	spdlog::debug("================================================================");
	spdlog::debug("                  constructGraphByGreedy Performance            ");
	spdlog::debug("================================================================");
	spdlog::debug("Main Operations:");
	spdlog::debug("  Producer (split/enqueue)  : {:>10.3f} ms", ns2ms(ns_producer.load()));
	spdlog::debug("  extendClusterToAnchor     : {:>10.3f} ms ({} calls, avg {:.3f} ms)",
	              ns2ms(ns_extend.load()), cnt_extend.load(),
	              ns2ms(ns_extend.load()) / std::max<uint64_t>(1, cnt_extend));
	spdlog::debug("  insertAnchorIntoGraph     : {:>10.3f} ms ({} calls, avg {:.3f} ms)",
	              ns2ms(ns_insert.load()), cnt_insert.load(),
	              ns2ms(ns_insert.load()) / std::max<uint64_t>(1, cnt_insert));
	spdlog::debug("");
	spdlog::debug("Task Statistics:");
	spdlog::debug("  Tasks submitted           : {}", cnt_submit.load());
	spdlog::debug("  Tasks finished            : {}", cnt_finish.load());
	spdlog::debug("");
	spdlog::debug("Producer Step Breakdown:");
	spdlog::debug("  check overlap             : {:>10.3f} ms", ns2ms(ns_check_overlap));
	spdlog::debug("  insertInterval            : {:>10.3f} ms", ns2ms(ns_insert_interval));
	spdlog::debug("  clone cluster (shared)    : {:>10.3f} ms", ns2ms(ns_clone_cluster));
	spdlog::debug("  enqueue task              : {:>10.3f} ms", ns2ms(ns_enqueue_task));
	spdlog::debug("  store kept                : {:>10.3f} ms", ns2ms(ns_store_kept));
	spdlog::debug("================================================================");
	spdlog::debug("");
#endif

	// *cluster_vec_ptr = std::move(kept);
}

/* ============================================================= *
 *  把三维 clusters  ->  按 ref 的一维 clusters
 *  并行对每个 ref 走 keepWithSplitGreedy，再写回三维结构
 * ============================================================= */
 /**
  * @brief  全局 MatchClusterVec 上执行“两级 map + 最大堆贪婪拆分”过滤
  *
  * @param cluster_vec_ptr  所有 clusters 的 shared_ptr
  * @param pool             ThreadPool（本实现单线程，参数仅留作占位）
  * @param min_span         最小跨度阈值
  */
void PairRareAligner::constructGraphByGreedyByRef(SpeciesName query_name, SeqPro::ManagerVariant& query_seqpro_manager, MatchClusterVecPtr cluster_vec_ptr, RaMesh::RaMeshMultiGenomeGraph& graph, uint_t min_span, bool isMultiple)
{
	if (!cluster_vec_ptr || cluster_vec_ptr->empty()) return;

	struct Node {
		MatchCluster cl;
		int_t        span;
	};
	auto cmp = [](const Node& a, const Node& b) { return a.span < b.span; };

	std::vector<Node> heap;
	heap.reserve(cluster_vec_ptr->size());

	for (auto& cl : *cluster_vec_ptr) {
		if (cl.empty()) continue;
		int_t sc = clusterSpan(cl);
		if (sc >= min_span)
			heap.push_back({ std::move(cl), sc });
	}
	cluster_vec_ptr->clear();
	std::make_heap(heap.begin(), heap.end(), cmp);

	std::vector<IntervalMap> rMaps;
	std::vector<IntervalMap> qMaps;
	

	// 初始化rMaps和qMaps

	MatchClusterVec kept;
	kept.reserve(heap.size());

	while (!heap.empty()) {

		std::pop_heap(heap.begin(), heap.end(), cmp);
		Node cur = std::move(heap.back());
		heap.pop_back();

		if (cur.span < min_span || cur.cl.empty()) continue;

		const ChrIndex refChr = cur.cl.front().ref_chr_index;
		const ChrIndex qChr = cur.cl.front().qry_chr_index;

		if (refChr >= rMaps.size()) {
			rMaps.resize(refChr + 1);
		}
		if (qChr >= qMaps.size()) {
			qMaps.resize(qChr + 1);
		}

		Strand strand = cur.cl.front().strand();
		uint_t rb = start1(cur.cl.front());
		uint_t re = start1(cur.cl.back()) + len1(cur.cl.back());
		uint_t qb = 0;
		uint_t qe = 0;
		if (strand == FORWARD) {
			qb = start2(cur.cl.front());
			qe = start2(cur.cl.back()) + len2(cur.cl.back());
		}
		else {
			qb = start2(cur.cl.back());
			qe = start2(cur.cl.front()) + len2(cur.cl.front());
		}

		int_t RL = 0, RR = 0, QL = 0, QR = 0;

		bool ref_hit = overlap1D(rMaps[refChr], rb, re, RL, RR);
		bool query_hit = overlap1D(qMaps[qChr], qb, qe, QL, QR);

		if (!ref_hit && !query_hit) {
			insertInterval(rMaps[refChr], rb, re);
			insertInterval(qMaps[qChr], qb, qe);
			auto task_cl = std::make_shared<MatchCluster>(cur.cl);
			
			AnchorVec anchor_vec = extendClusterToAnchorVec(*task_cl, *ref_seqpro_manager, query_seqpro_manager);
			for (auto& anchor : anchor_vec) {

				graph.insertAnchorIntoGraph(*ref_seqpro_manager,query_seqpro_manager, ref_name, query_name, anchor, isMultiple);

			}
			
			
			//pool.enqueue([this, &graph, query_name, task_cl, &query_seqpro_manager] {
			//	AnchorVec anchor_vec = extendClusterToAnchor(*task_cl, *ref_seqpro_manager, query_seqpro_manager);
			//	graph.insertAnchorIntoGraph(ref_name, query_name, anchor_vec);
			//	// graph.insertClusterIntoGraph(ref_name, query_name, *task_cl);
			//	});

			kept.emplace_back(cur.cl);

			continue;
		}

		for (auto& part : splitCluster(cur.cl, ref_hit, RL, RR, query_hit, QL, QR)) {
			int_t sc = clusterSpan(part);
			if (sc >= min_span) {
				heap.push_back({ part, sc });
				std::push_heap(heap.begin(), heap.end(), cmp);
			}
		}
	}
	return;
}

void PairRareAligner::constructGraphByDpByRef(SpeciesName query_name, SeqPro::ManagerVariant& query_seqpro_manager, MatchClusterVecPtr cluster_vec_ptr, RaMesh::RaMeshMultiGenomeGraph& graph, ThreadPool& pool,uint_t thread_num, uint_t min_span, bool isMultiple)
{
	if (!cluster_vec_ptr || cluster_vec_ptr->empty()) return;

	std::vector<Anchor> anchors;
	anchors.resize(cluster_vec_ptr->size());

	size_t cluster_num = cluster_vec_ptr->size();
	size_t chunk = (cluster_vec_ptr->size() + thread_num - 1) / thread_num;   // 向上取整

	std::vector<std::future<void>> futs;
	for (size_t t = 0; t < thread_num && t * chunk < cluster_num; ++t) {
		size_t beg = t * chunk;
		size_t end = std::min(cluster_num, beg + chunk);

		futs.emplace_back(
			pool.enqueue([&, beg, end] {
				for (size_t i = beg; i < end; ++i) {
					anchors[i] =
						extendClusterToAnchor((*cluster_vec_ptr)[i],
							*ref_seqpro_manager,
							query_seqpro_manager);
				}
				})
		);
	}
	for (auto& f : futs) f.get();

	// 释放cluster_vec_ptr
	if (cluster_vec_ptr) {
		cluster_vec_ptr->clear();   // 可选：先把元素析构掉
		cluster_vec_ptr->shrink_to_fit(); // 可选：把 capacity 也回收给 STL 实现
		cluster_vec_ptr.reset();    // 关键：降低引用计数 / 转移所有权
	}

	/* --------------------------------------------- *
	 *  1. 直接按参考起点升序排序 anchors             *
	 * --------------------------------------------- */
	std::sort(anchors.begin(), anchors.end(), [](const Anchor& a, const Anchor& b) {
		return a.ref_start < b.ref_start;
		});

	const size_t n = anchors.size();
	if (n == 0) return;

	constexpr size_t K = 50;                   // 只看最近 K 个前驱
	std::vector<int_t> dp(n, 0);           // dp[i] = 以 i 结尾的最大得分
	std::vector<int_t>       prev(n, -1);        // 链前驱索引

	int_t best_score = 0;
	int_t    best_pos = 0;

	for (size_t i = 0; i < n; ++i) {
		const auto ai_start = anchors[i].ref_start;
		const auto ai_end = ai_start + anchors[i].ref_len - 1;

		dp[i] = anchors[i].alignment_length;

		/* ---------- 只回溯最近 K 个元素 ---------- */
		size_t j_beg = (i > K) ? i - K : 0;
		for (size_t j = i; j-- > j_beg; ) {
			const auto aj_end =
				anchors[j].ref_start +
				anchors[j].ref_len - 1;

			if (aj_end < ai_start) {                  // 无重叠
				long long cand = dp[j] + anchors[i].alignment_length;
				if (cand > dp[i]) {
					dp[i] = cand;
					prev[i] = static_cast<int>(j);
				}
			}
		}

		if (dp[i] > best_score) {
			best_score = dp[i];
			best_pos = i;
		}
	}

	/* --------------------------------------------- *
	 *  2. 回溯得到最长链并替换 anchors               *
	 * --------------------------------------------- */
	std::vector<Anchor> lis_chain;
	for (int p = static_cast<int>(best_pos); p != -1; p = prev[p])
		lis_chain.emplace_back(std::move(anchors[p]));
	std::reverse(lis_chain.begin(), lis_chain.end());

	anchors.swap(lis_chain);   // anchors 只保留 LIS
	std::vector<Anchor>().swap(lis_chain);

	for (auto& anchor : anchors) {
		graph.insertAnchorIntoGraph(*ref_seqpro_manager, query_seqpro_manager, ref_name, query_name, anchor, isMultiple);
	}

	return;
}

// ClusterVecPtrByStrandByQueryRefPtr PairRareAligner::filterPairSpeciesAnchors(SpeciesName query_name, MatchVec3DPtr& anchors, SeqPro::ManagerVariant& query_fasta_manager, RaMesh::RaMeshMultiGenomeGraph& graph, uint_t min_span)
// {
//
// 	MatchByStrandByQueryRefPtr unique_anchors = std::make_shared<MatchByStrandByQueryRef>();;
// 	MatchByStrandByQueryRefPtr repeat_anchors = std::make_shared<MatchByStrandByQueryRef>();;
//
// 	groupMatchByQueryRef(anchors, unique_anchors, repeat_anchors,
// 		*ref_seqpro_manager, query_fasta_manager);
//
// 	anchors.reset();
// 	spdlog::info("groupMatchByQueryRef done");
//
// 	ClusterVecPtrByStrandByQueryRefPtr cluster_ptr = clusterAllChrMatch(
// 		unique_anchors,
// 		repeat_anchors,
// 		min_span);
//
//
// 	spdlog::info("clusterAllChrMatch done");
//
// 	return cluster_ptr;
//
// }
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
            && a.qry_start <= c_qry_low && (a.qry_start + a.qry_len) >= (c_qry_high - c_qry_low);
    }
}

// ========== 工具：把 gap 用 extendAlignKSW2 ends-free 对齐并“必须到达目标” ==========
static bool extend_gap_must_reach(const std::string& ref_gap,
                                         const std::string& qry_gap,
                                         Cigar_t& out) {
    // 你项目里的 ends-free 对齐实现：extendAlignKSW2(ref, query, zdrop)
    // （已在 align.h / align.cpp 定义）:contentReference[oaicite:0]{index=0} :contentReference[oaicite:1]{index=1}
    out = extendAlignKSW2(ref_gap, qry_gap, /*zdrop=*/200);
    return countRefLength(out) == ref_gap.size() && countQryLength(out) == qry_gap.size() && checkGapCigarQuality(out, ref_gap.size(), qry_gap.size(), 0.6);
}

Anchor* selectBestPrevAnchor(
	const AnchorPtrVec& anchors,
	Strand strand,
	Coord_t curr_ref_start,
	Coord_t curr_qry_for_gap,
	int look_back = 500
) {
	Anchor* best = nullptr;
	int_t   best_score = std::numeric_limits<int_t>::max();

	int checked = 0;
	for (auto it = anchors.rbegin(); it != anchors.rend() && checked < look_back; ++it, ++checked) {
		const Anchor* prev = it->get();
		if (!prev) continue;


		// 计算 gap（必须非负；负数代表重叠/倒退）
		const int_t ref_gap = static_cast<int_t>(curr_ref_start)
							- static_cast<int_t>(prev->ref_start + prev->ref_len);

		int_t qry_gap;
		if (strand == FORWARD) {
			// FORWARD: 当前 sB - prev.q_end
			qry_gap = static_cast<int_t>(curr_qry_for_gap)
					- static_cast<int_t>(prev->qry_start + prev->qry_len);
		} else {
			// REVERSE: prev.q_start - (当前 sB+len) —— 右端对右端
			qry_gap = static_cast<int_t>(prev->qry_start)
					- static_cast<int_t>(curr_qry_for_gap);
		}

		if (ref_gap < 0 || qry_gap < 0) continue;

		const long greater = std::max(ref_gap, qry_gap);
		const long lesser  = std::min(ref_gap, qry_gap);
		const int_t score  = static_cast<int_t>((greater << 1) - lesser);

		// 分数更小即更优；分数相等时保留当前 best（更靠近尾部者自然优先）
		if (score < best_score) {
			best_score = score;
			best = const_cast<Anchor*>(prev);
		}
	}

	return best;
}

// 依赖：AnchorPtr, MatchClusterVec, uint_t, int_t, len2(const Match&)
// Anchor 字段：ref_start, ref_len, qry_start, qry_len, strand, ref_chr_index, qry_chr_index
// Match  字段：ref_start, qry_start, strand(), ref_chr_index, qry_chr_index

uint_t selectBestForwardCluster(AnchorPtr ap, MatchClusterVec cluster_vec, uint_t cur_idx)
{
    const size_t N = cluster_vec.size();
    if (!ap || cur_idx + 1 >= N) return static_cast<uint_t>(N);

    const Strand s   = ap->strand;
    const auto  rchr = ap->ref_chr_index;
    const auto  qchr = ap->qry_chr_index;

    // 计算单个 match 的分数；若冲突（负 gap / 染色体或链向不一致）返回 false
    auto score_match = [&](const Match& m, int_t& out_score) -> bool {
        if (m.strand() != s || m.ref_chr_index != rchr || m.qry_chr_index != qchr) return false;

        // ref 方向：anchor 的前向端在 ref_end（半开区间端点），候选必须在其之后
        const int_t ref_end = static_cast<int_t>(ap->ref_start + ap->ref_len);
        const int_t ref_gap = static_cast<int_t>(m.ref_start) - ref_end;
        if (ref_gap < 0) return false;

        // query 方向：正向用 q_end；反向用 anchor 的小坐标端与 match 的“右端”比较
        int_t qry_gap = 0;
        if (s == FORWARD) {
            const int_t qry_end = static_cast<int_t>(ap->qry_start + ap->qry_len);
            qry_gap = static_cast<int_t>(m.qry_start) - qry_end;
        } else {
            const int_t m_q_front = static_cast<int_t>(m.qry_start + len2(m)); // match 的“右端”
            qry_gap = static_cast<int_t>(ap->qry_start) - m_q_front;           // anchor 的“前向端”（小坐标端）
        }
        if (qry_gap < 0) return false;

        const long greater = std::max(ref_gap, qry_gap);
        const long lesser  = std::min(ref_gap, qry_gap);
        out_score = static_cast<int_t>((greater << 1) - lesser);
        return true;
    };

    const size_t LOOK_AHEAD = 500;
    const size_t j_end = std::min(N, static_cast<size_t>(cur_idx + 1 + LOOK_AHEAD));

    bool found = false;
    uint_t best_idx = static_cast<uint_t>(N);
    int_t  best_score = std::numeric_limits<int_t>::max();

    for (size_t j = cur_idx + 1; j < j_end; ++j) {
        const auto& cl = cluster_vec[j];
        if (cl.empty()) continue;

        // 只考虑同链向/同 chr 的簇（用 front 做快速过滤）
        if (cl.front().strand() != s ||
            cl.front().ref_chr_index != rchr ||
            cl.front().qry_chr_index != qchr)
            continue;

        int_t sc = 0;
        bool ok  = false;

        // 先试 front match
        ok = score_match(cl.front(), sc);

        // 若 front 冲突，则遍历整簇找一个可行的最小分
        if (!ok) {
            int_t local_best = std::numeric_limits<int_t>::max();
            bool  any = false;
            for (const auto& m : cl) {
                int_t s_tmp;
                if (score_match(m, s_tmp)) {
                    any = true;
                    if (s_tmp < local_best) local_best = s_tmp;
                	break;
                }
            }
            if (any) { sc = local_best; ok = true; }
        }

        if (!ok) continue; // 整簇均冲突，跳过

        if (sc < best_score) {
            best_score = sc;
            best_idx   = static_cast<uint_t>(j);
            found      = true;
        }
    }

    return found ? best_idx : static_cast<uint_t>(N);
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


AnchorPtrVec extendClustersToAnchors(const MatchClusterVecPtr& cluster_vec_ptr,
                                            const SeqPro::ManagerVariant& ref_mgr,
                                            const SeqPro::ManagerVariant& qry_mgr) {
    AnchorPtrVec anchors;
    if (!cluster_vec_ptr || cluster_vec_ptr->empty()) return anchors;

    // 1) 参照 A（ref）起点升序排序（MUMmer 同步骤）:contentReference[oaicite:2]{index=2}
    auto clusters = *cluster_vec_ptr; // 本地拷贝可重排
    std::sort(clusters.begin(), clusters.end(), [](const MatchCluster& a, const MatchCluster& b){
        if (a.empty() || b.empty()) return a.size() < b.size();
        return start1(a.front()) < start1(b.front());
    });

    const size_t N = clusters.size();
    std::vector<char> fused(N, 0);

    // 2) 外层 while：按 MUMmer 逻辑在“当前簇”上做内部延伸 & 选取“前向目标簇”并尝试融合:contentReference[oaicite:3]{index=3}
    size_t curr_idx = 0;                  // CurrCp
    size_t prev_idx = 0;                  // PrevCp
    bool   target_reached = false;        // 是否已“接上目标”
    size_t target_idx = N;                // TargetCp（默认为 end）
	/// 80-90    100-110 120-130 (100-130)
	/// 120-130   100-110 80-90 (80-110)

	AnchorPtr cur_anchor = nullptr;


    while (curr_idx < N)
    {
	    auto& cl = clusters[curr_idx];
    	if (cl.empty())
    	{
    		++curr_idx;
    		continue;
    	}

    	// 按 ref 坐标排序该簇内部匹配；反链也统一按 ref 升序（你已有的实现也这么做）:contentReference[oaicite:4]{index=4}
    	std::sort(cl.begin(), cl.end(), [](const Match& x, const Match& y){ return x.ref_start < y.ref_start; });

    	// 跳过已经融合过或被“shadow”的 cluster（MUMmer：wasFused / isShadowedCluster）:contentReference[oaicite:5]{index=5}
    	bool skip = false;
    	if (fused[curr_idx]) skip = true;
    	if (!skip) {
    		for (const auto& ap : anchors) {
    			if (shadowedBy(*ap, cl)) { skip = true; break; }
    		}
    	}
    	if (skip) { fused[curr_idx] = 1; curr_idx++; target_reached = false; continue; }

    	// 3) 选择方向；对反链我们在提取 query gap 时做 reverseComplement，保证传给 KSW2 的方向一致:contentReference[oaicite:6]{index=6}
    	Strand strand = cl.front().strand();
    	bool   fwd    = (strand == FORWARD);
    	bool reach_target = false;

    	ChrIndex ref_chr = cl.front().ref_chr_index;
    	ChrIndex qry_chr = cl.front().qry_chr_index;
    	// 当前 anchor 的两端坐标
    	int_t ref_beg = start1(cl.front());
    	int_t ref_end = start1(cl.back()) + len1(cl.back());
    	int_t qry_beg = fwd ? start2(cl.front()) : start2(cl.back());

    	int_t qry_end = fwd ? (start2(cl.back()) + len2(cl.back()))
						 : (start2(cl.front()) + len2(cl.front()));


    	if (!anchors.empty() && reach_target == false) {

    		Anchor* best = nullptr;
    		// 在anchor中倒着遍历500个，计算
    		// long greater = std::max(ref_gap, qry_gap);
    		// long lesser = std::min(ref_gap, qry_gap);
    		// int_t this_score = (greater << 1) - lesser;
    		// 选择分数最小的anchor作为best
    		best = selectBestPrevAnchor(anchors, strand, ref_beg, qry_beg, 500);

    		if (best) {
    			// 取 gap：best 的末端 → 本簇首 match 起点
    			int_t ref_gap_beg = best->ref_start + best->ref_len;
    			int_t ref_gap_len = ref_beg - ref_gap_beg;

    			int_t qry_gap_beg, qry_gap_len;
    			if (fwd) {
    				qry_gap_beg = best->qry_start + best->qry_len;
    				qry_gap_len = qry_beg - qry_gap_beg;
    			} else {
    				// 反链：上一 anchor 在 query 上的“右端”为（qry_start），向本簇首的“左端”（qry_beg 是大坐标端）
    				qry_gap_beg = qry_end;                  // 先以低坐标取片段，再统一 RC
    				qry_gap_len = qry_end - best->qry_start ;
    			}

    			if (ref_gap_len > 0 && qry_gap_len > 0 && ref_gap_len < 10000 && qry_gap_len < 10000) {
    				std::string ref_gap = subSeq(ref_mgr, ref_chr, ref_gap_beg, ref_gap_len);
    				std::string qry_gap = fwd
						? subSeq(qry_mgr, qry_chr, qry_gap_beg, qry_gap_len)
						: subSeq(qry_mgr, qry_chr, qry_gap_beg, qry_gap_len);
    				if (!fwd) reverseComplement(qry_gap); // 反链统一方向:contentReference[oaicite:8]{index=8}

    				Cigar_t gap_cigar;
    				if (extend_gap_must_reach(ref_gap, qry_gap, gap_cigar)) {
    					// 合并到 best（等价 CurrAp = TargetAp）
    					appendCigar(best->cigar, gap_cigar);
    					best->alignment_length += countAlignmentLength(gap_cigar);
    					best->aligned_base     += countMatchOperations(gap_cigar);
    					best->ref_len          += ref_gap_len;
    					best->qry_len          += qry_gap_len;
    					if (!fwd)
    					{
    						best->qry_start -= qry_gap_len;
    					}
    				}
    			}
    		}
    	}

    	// // 6) 遍历本簇的每个 match：先追加 M，再把中间 gap 用 ends-free 延伸到“下一个 match 起点”（若成功即 target_reached=true）

    	if (reach_target == false)
    	{
    		cur_anchor = std::make_shared<Anchor>();
    		cur_anchor->ref_chr_index = ref_chr;
			cur_anchor->qry_chr_index = qry_chr;
    		cur_anchor->strand = strand;
    		cur_anchor->ref_start = start1(cl.front());
    		cur_anchor->qry_start = start2(cl.front());
			cur_anchor->ref_len = len1(cl.front());
    		cur_anchor->qry_len = len2(cl.front());
    		cur_anchor->aligned_base = len1(cl.front());
    		cur_anchor->alignment_length = len1(cl.front());;
    		cur_anchor->cigar = {};
    		appendCigarOp(cur_anchor->cigar, 'M',  len1(cl.front()));
    		prev_idx = curr_idx;
    	}

  //   	for (auto m : cl) {
  //   		// 判断和现有anchor是否冲突
		// 	if (isMatchConflictingWithAnchor(m, *cur_anchor)) continue;
		// 	// 追加 match
		// 	appendCigarOp(cur_anchor->cigar, 'M', len1(m));
		// 	cur_anchor->aligned_base += len1(m);
  //   		cur_anchor->alignment_length += len1(m);
	 //
		// 	// 计算anchor和m之间的gap，进行比对
		// }
    	// 建议把这段 for 改成：从第二个 match 开始（第一个已作为种子）
		for (size_t i = 0; i < cl.size(); ++i) {
		    const Match& m = cl[i];

		    // 1) 先做几何冲突判定（负 gap/倒退/不同 chr 或 strand）
		    if (isMatchConflictingWithAnchor(m, *cur_anchor)) {
		        continue; // 冲突直接跳过该 match
		    }

		    // 2) 计算 ref / qry 的 gap 长度与抽取范围
		    const bool fwd = (cur_anchor->strand == FORWARD);

		    // anchor 的“前向端”（ref）
		    Coord_t ca_ref_end = cur_anchor->ref_start + cur_anchor->ref_len;
		    // gap on ref: [ca_ref_end, m.ref_start)
		    int_t ref_gap_len = static_cast<int_t>(m.ref_start) - static_cast<int_t>(ca_ref_end);

		    // gap on qry
		    int_t qry_gap_len = 0;
		    Coord_t q_low = 0, q_high = 0, q_beg = 0; // 取片时使用
		    if (fwd) {
		        // FORWARD: 前向端在 (qry_start + qry_len)
		        Coord_t ca_qry_end = cur_anchor->qry_start + cur_anchor->qry_len;
		        qry_gap_len = static_cast<int_t>(m.qry_start) - static_cast<int_t>(ca_qry_end);
		        q_beg = ca_qry_end; // [ca_qry_end, m.qry_start)
		    } else {
		        // REVERSE: 前向端在小坐标端 = anchor.qry_start
		        // 目标 match 的“前向端”是 (m.qry_start + len(m))
		        Coord_t ca_qry_front = cur_anchor->qry_start;         // 小坐标端
		        Coord_t m_q_front    = m.qry_start + len2(m);          // 目标“右端”
		        q_low  = m_q_front;
		        q_high = ca_qry_front;
		        qry_gap_len = (q_high > q_low) ? static_cast<int_t>(q_high - q_low) : 0;
		    }

		    // 3) 若有 gap，需要 ends-free 吃满才能接上目标；失败就跳过该 match

		    Cigar_t gap_cigar;

		    if (ref_gap_len > 0 &&  qry_gap_len > 0 && ref_gap_len < 10000 && qry_gap_len < 10000) {
		        // 提取序列
		        std::string ref_gap_seq = (ref_gap_len > 0)
		            ? subSeq(ref_mgr, cur_anchor->ref_chr_index, ca_ref_end, static_cast<Coord_t>(ref_gap_len))
		            : std::string();

		        std::string qry_gap_seq;
		        if (fwd) {
		            if (qry_gap_len > 0)
		                qry_gap_seq = subSeq(qry_mgr, cur_anchor->qry_chr_index, q_beg, static_cast<Coord_t>(qry_gap_len));
		        } else {
		            if (qry_gap_len > 0) {
		                // 反向：以低→高取片，再做 RC
		                qry_gap_seq = subSeq(qry_mgr, cur_anchor->qry_chr_index, q_low, static_cast<Coord_t>(qry_gap_len));
		                reverseComplement(qry_gap_seq);
		            }
		        }
				gap_cigar = globalAlignKSW2(ref_gap_seq, qry_gap_seq);
		        //reached = extend_gap_must_reach(ref_gap_seq, qry_gap_seq, gap_cigar);


		        // 把 gap 的对齐并入 cur_anchor
		        appendCigar(cur_anchor->cigar, gap_cigar);
		        cur_anchor->alignment_length += countAlignmentLength(gap_cigar);
		        cur_anchor->aligned_base     += countMatchOperations(gap_cigar);
		        cur_anchor->ref_len          += countRefLength(gap_cigar);

		        if (fwd) {
		            cur_anchor->qry_len      += countQryLength(gap_cigar);
		        } else {
		            // 反向：向“更小坐标”扩展，需左移起点并增大长度
		            Coord_t qadd = countQryLength(gap_cigar);
		            cur_anchor->qry_start -= qadd;
		            cur_anchor->qry_len   += qadd;
		        }
		    }else
		    {
			    break;
		    }

		    // 4) 已到达目标 match 的起点：按 MUMmer 规则只追加 len-1
		    uint32_t L = len1(m);
		    if (L > 0) {
		        uint32_t addL = L;
		        if (addL > 0) {
		            appendCigarOp(cur_anchor->cigar, 'M', addL);
		            cur_anchor->alignment_length += addL;
		            cur_anchor->aligned_base     += addL;
		            cur_anchor->ref_len          += addL;
		            if (fwd) {
		                cur_anchor->qry_len      += addL;
		            } else {
		                cur_anchor->qry_start    -= addL;
		                cur_anchor->qry_len      += addL;
		            }
		        }
		    }
		}
    	fused[curr_idx] = 1;

    	// 找下一个最合适的cluster，找到之后尝试扩展
		uint_t target_idx = selectBestForwardCluster(cur_anchor, clusters, static_cast<uint_t>(curr_idx));
    	// 找到之后，找到不冲突的match，然后尝试将anchor进行扩展
    	if (target_idx < N)
    	{
    		MatchVec target_match = clusters[target_idx];

    		for (auto &m : target_match)
    		{
    			// 1) 先做几何冲突判定（负 gap/倒退/不同 chr 或 strand）
    			if (isMatchConflictingWithAnchor(m, *cur_anchor)) {
    				continue; // 冲突直接跳过该 match
    			}
    			// 2) 计算 ref / qry 的 gap 长度与抽取范围
			    const bool fwd = (cur_anchor->strand == FORWARD);

			    // anchor 的“前向端”（ref）
			    Coord_t ca_ref_end = cur_anchor->ref_start + cur_anchor->ref_len;
			    // gap on ref: [ca_ref_end, m.ref_start)
			    int_t ref_gap_len = static_cast<int_t>(m.ref_start) - static_cast<int_t>(ca_ref_end);

			    // gap on qry
			    int_t qry_gap_len = 0;
			    Coord_t q_low = 0, q_high = 0, q_beg = 0; // 取片时使用
			    if (fwd) {
			        // FORWARD: 前向端在 (qry_start + qry_len)
			        Coord_t ca_qry_end = cur_anchor->qry_start + cur_anchor->qry_len;
			        qry_gap_len = static_cast<int_t>(m.qry_start) - static_cast<int_t>(ca_qry_end);
			        q_beg = ca_qry_end; // [ca_qry_end, m.qry_start)
			    } else {
			        // REVERSE: 前向端在小坐标端 = anchor.qry_start
			        // 目标 match 的“前向端”是 (m.qry_start + len(m))
			        Coord_t ca_qry_front = cur_anchor->qry_start;         // 小坐标端
			        Coord_t m_q_front    = m.qry_start + len2(m);          // 目标“右端”
			        q_low  = m_q_front;
			        q_high = ca_qry_front;
			        qry_gap_len = (q_high > q_low) ? static_cast<int_t>(q_high - q_low) : 0;
			    }

			    // 3) 若有 gap，需要 ends-free 吃满才能接上目标；失败就跳过该 match
			    bool reached = true;
			    Cigar_t gap_cigar;

			    if (ref_gap_len > 0 && qry_gap_len > 0 && ref_gap_len < 10000 && qry_gap_len < 10000) {
			        // 提取序列
			        std::string ref_gap_seq = (ref_gap_len > 0)
			            ? subSeq(ref_mgr, cur_anchor->ref_chr_index, ca_ref_end, static_cast<Coord_t>(ref_gap_len))
			            : std::string();

			        std::string qry_gap_seq;
			        if (fwd) {
			            if (qry_gap_len > 0)
			                qry_gap_seq = subSeq(qry_mgr, cur_anchor->qry_chr_index, q_beg, static_cast<Coord_t>(qry_gap_len));
			        } else {
			            if (qry_gap_len > 0) {
			                // 反向：以低→高取片，再做 RC
			                qry_gap_seq = subSeq(qry_mgr, cur_anchor->qry_chr_index, q_low, static_cast<Coord_t>(qry_gap_len));
			                reverseComplement(qry_gap_seq);
			            }
			        }

			        reached = extend_gap_must_reach(ref_gap_seq, qry_gap_seq, gap_cigar);
			        if (!reached) {
		        		target_reached = false;
		        		break;
			        }

			        // 把 gap 的对齐并入 cur_anchor
			        appendCigar(cur_anchor->cigar, gap_cigar);
			        cur_anchor->alignment_length += countAlignmentLength(gap_cigar);
			        cur_anchor->aligned_base     += countMatchOperations(gap_cigar);
			        cur_anchor->ref_len          += countRefLength(gap_cigar);

			        if (fwd) {
			            cur_anchor->qry_len      += countQryLength(gap_cigar);
			        } else {
			            // 反向：向“更小坐标”扩展，需左移起点并增大长度
			            Coord_t qadd = countQryLength(gap_cigar);
			            cur_anchor->qry_start -= qadd;
			            cur_anchor->qry_len   += qadd;
			        }
		    		target_reached = true;
		    		break;
			    }else
			    {
			    	target_reached = false;
			    	break;
			    }

    		}
    	} else
    	{
    		target_reached = false;
    	}

    	if (target_reached == true)
    	{
    		curr_idx = target_idx;

    	}else
    	{
    		anchors.push_back(cur_anchor);
    		curr_idx = prev_idx+1;
    	}

    }
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
        AnchorPtrVec anchors;
        anchors.reserve(1);

        if (!is_first) {
            for (auto& c : (*tmp_p)) {
                for (auto& sub_c : c) {
                    MatchVec mc;
                    mc.push_back(sub_c);

                    Anchor anchor = extendClusterToAnchor(
                        mc, *ref_seqpro_manager, query_seqpro_manager
                    );
                    anchors.push_back(std::make_shared<Anchor>(std::move(anchor)));
                }
            }
        } else {
            anchors = linkClusters(
                *tmp_p, *ref_seqpro_manager, query_seqpro_manager
            );
        }

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



static void filterChrByDP(
	AnchorBySQR_SparsePtr anchor_map,
	uint_t id,
	bool filter_ref)
{
	AnchorPtrVec result;

	if (!anchor_map) return;

	for (auto & a_vec : *anchor_map)
	{
		if (a_vec.size() > 0)
		{
			uint_t ref_idx = a_vec[0]->ref_chr_index;
			uint_t qry_idx = a_vec[0]->qry_chr_index;
			if (filter_ref)
			{
				if (ref_idx == id)
				{
					result.insert(result.end(), a_vec.begin(), a_vec.end());
				}
			}else
			{
				if (qry_idx == id)
				{
					result.insert(result.end(), a_vec.begin(), a_vec.end());
				}
			}
		}
	}

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
		std::cout << "✅ No overlaps detected in final selection.\n";
	else
		std::cerr << "❌ Overlaps found in final selection!\n";
#endif
}

void PairRareAligner::filterAnchorByDP(AnchorBySQR_SparsePtr anchor_map, uint_t ref_chr_cnt, uint_t qry_chr_cnt) {

	const uint_t n = anchor_map->size();
#pragma omp parallel for schedule(dynamic) num_threads(thread_num)
	for (uint_t i = 0; i < ref_chr_cnt; i++) {
		filterChrByDP(anchor_map, i, true);
	}

#pragma omp parallel for schedule(dynamic) num_threads(thread_num)
	for (uint_t i = 0; i < qry_chr_cnt; i++) {
		filterChrByDP(anchor_map, i, false);
	}

	return;

}

void PairRareAligner::constructGraphByDP(SpeciesName query_name, SeqPro::ManagerVariant& query_seqpro_manager, AnchorBySQR_SparsePtr anchor_ptr, RaMesh::RaMeshMultiGenomeGraph& graph) {
	for (auto & a_vec : *anchor_ptr)
	{
		for (size_t a_idx = 0; a_idx < a_vec.size(); ++a_idx)
		{
			AnchorPtr anchor = a_vec.at(a_idx);
			if (anchor->ref_selected && anchor->qry_selected) {
				try {
					graph.insertAnchorIntoGraph(*ref_seqpro_manager, query_seqpro_manager, ref_name, query_name, *anchor, false);
				}
				catch (const std::exception& e) {
					spdlog::error("Error inserting anchor into graph: {}", e.what());
				}
			}
		}
	}

}






