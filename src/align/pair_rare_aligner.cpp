#include "rare_aligner.h"
#include "anchor.h"

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

MatchVec3DPtr PairRareAligner::alignPairGenome(
	SpeciesName query_name,
	SeqPro::ManagerVariant& query_fasta_manager,
	SearchMode         search_mode,
	bool allow_MEM,
	bool allow_short_mum,
	sdsl::int_vector<0>& ref_global_cache,
	SeqPro::Length sampling_interval) {

	ThreadPool shared_pool(thread_num);

	MatchVec3DPtr anchors = findQueryFileAnchor(
		query_name, query_fasta_manager, search_mode, allow_MEM, allow_short_mum, shared_pool, ref_global_cache, sampling_interval);

	return anchors;

}

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

ClusterVecPtrByStrandByQueryRefPtr PairRareAligner::filterPairSpeciesAnchors(SpeciesName query_name, MatchVec3DPtr& anchors, SeqPro::ManagerVariant& query_fasta_manager, RaMesh::RaMeshMultiGenomeGraph& graph, uint_t min_span)
{
	ThreadPool shared_pool(thread_num);
	MatchByStrandByQueryRefPtr unique_anchors = std::make_shared<MatchByStrandByQueryRef>();;
	MatchByStrandByQueryRefPtr repeat_anchors = std::make_shared<MatchByStrandByQueryRef>();;

	groupMatchByQueryRef(anchors, unique_anchors, repeat_anchors,
		*ref_seqpro_manager, query_fasta_manager, shared_pool);

	shared_pool.waitAllTasksDone();
	//anchors->clear();
	//anchors->shrink_to_fit();
	anchors.reset();
	spdlog::info("groupMatchByQueryRef done");
	sortMatchByQueryStart(unique_anchors, shared_pool);
	sortMatchByQueryStart(repeat_anchors, shared_pool);

	shared_pool.waitAllTasksDone();
	spdlog::info("sortMatchByQueryStart done");

	ClusterVecPtrByStrandByQueryRefPtr cluster_ptr = clusterAllChrMatch(
		unique_anchors,
		repeat_anchors,
		shared_pool, min_span);

	shared_pool.waitAllTasksDone();

	spdlog::info("clusterAllChrMatch done");

	return cluster_ptr;

}

AnchorPtrVecByStrandByQueryByRefPtr PairRareAligner::extendClusterToAnchorByChr(SpeciesName query_name, SeqPro::ManagerVariant& query_seqpro_manager, ClusterVecPtrByStrandByQueryRefPtr cluster, ThreadPool& pool, bool is_first) {
	// 输出结构
	auto result = std::make_shared<AnchorPtrVecByStrandByQueryByRef>();

	//ThreadPool pool(thread_num);
	std::vector<std::future<void>> futures;

	for (size_t strand = 0; strand < 2; ++strand) {
		const auto& byQueryRef = (*cluster)[strand];
		result->emplace_back(); // 对应 strand

		for (size_t q = 0; q < byQueryRef.size(); ++q) {
			const auto& byRef = byQueryRef[q];
			(*result)[strand].emplace_back(); // 对应 query_chr

			for (size_t r = 0; r < byRef.size(); ++r) {
				MatchClusterVecPtr cluster_vec_ptr = byRef[r];
				(*result)[strand][q].emplace_back(); // 对应 ref_chr

				if (!cluster_vec_ptr) continue;
				//if (strand == 1) continue;

				// 提交任务：处理一个 ClusterVec
				futures.emplace_back(pool.enqueue(
					[&, strand, q, r, cluster_vec_ptr]() {
						AnchorPtrVec anchors;
						anchors.reserve(1);
						// TODO 二轮之后的处理方式有点粗糙
						if (!is_first)
						{
							for (auto & c : (*cluster_vec_ptr)) {
								for (auto & sub_c : c)
								{	MatchVec mc;
									mc.push_back(sub_c);
									Anchor anchor = extendClusterToAnchor(mc, *ref_seqpro_manager, query_seqpro_manager);
									AnchorPtr p = std::make_shared<Anchor>(std::move(anchor));
									anchors.push_back(p);
								}
							}

						}else {
							for (auto & c : (*cluster_vec_ptr)) {
								Anchor anchor = extendClusterToAnchor(c, *ref_seqpro_manager, query_seqpro_manager);
								AnchorPtr p = std::make_shared<Anchor>(std::move(anchor));
								anchors.push_back(p);
							}
							linkClusters(anchors, *ref_seqpro_manager, query_seqpro_manager);
						}

						
						if (anchors.empty()) return;
						(*result)[strand][q][r] = anchors;

					}
				));
			}
		}
	}

	// 等待所有任务完成
	for (auto& f : futures) f.get();



    spdlog::info("extend cluster to anchor successfully for {}", query_name);
	return result;
}

static void filterChrByDP(
	AnchorPtrVecByStrandByQueryByRefPtr anchor_map,
	uint_t id,
	bool filter_ref)
{
	AnchorPtrVec result;

	if (!anchor_map) return;

	if (filter_ref) {
		// 按 ref 来过滤：只看 ref == id
		for (size_t strand_idx = 0; strand_idx < anchor_map->size(); ++strand_idx) {
			const auto& queries = anchor_map->at(strand_idx);
			for (size_t qry_idx = 0; qry_idx < queries.size(); ++qry_idx) {
				const auto& refs = queries.at(qry_idx);
				if (id < refs.size()) {
					const auto& anchors = refs[id];
					// AnchorPtrVec temp;
					// for (const auto& a : anchors) {
					// 	if (a->qry_selected) {
					// 		temp.push_back(a);
					// 	}
					// }
					// result.insert(result.end(), temp.begin(), temp.end());
					result.insert(result.end(), anchors.begin(), anchors.end());
				}
			}
		}
	}
	else {
		// 按 qry 来过滤：只看 qry == id
		for (size_t strand_idx = 0; strand_idx < anchor_map->size(); ++strand_idx) {
			const auto& queries = anchor_map->at(strand_idx);
			if (id < queries.size()) {
				const auto& refs = queries.at(id);
				for (const auto& anchors : refs) {
					/*AnchorPtrVec temp;
					for (const auto& a : anchors) {
						if (a->ref_selected) {
							temp.push_back(a);
						}
					}
					result.insert(result.end(), temp.begin(), temp.end());*/

					result.insert(result.end(), anchors.begin(), anchors.end());
				}
			}
		}
	}
	//if (filter_ref) {
	//	// 按 ref 来过滤：只看 ref == id
	//	for (size_t strand_idx = 0; strand_idx < anchor_map->size(); ++strand_idx) {
	//		const auto& queries = anchor_map->at(strand_idx);
	//		for (size_t qry_idx = 0; qry_idx < queries.size(); ++qry_idx) {
	//			const auto& refs = queries.at(qry_idx);
	//			if (id < refs.size()) {
	//				const auto& anchors = refs[id];
	//				for (const auto& a : anchors) {
	//					if (!a->is_linked) {              // ✅ 只保留 is_linked == false
	//						result.push_back(a);
	//					}
	//				}
	//			}
	//		}
	//	}
	//}
	//else {
	//	// 按 qry 来过滤：只看 qry == id
	//	for (size_t strand_idx = 0; strand_idx < anchor_map->size(); ++strand_idx) {
	//		const auto& queries = anchor_map->at(strand_idx);
	//		if (id < queries.size()) {
	//			const auto& refs = queries.at(id);
	//			for (const auto& anchors : refs) {
	//				for (const auto& a : anchors) {
	//					if (!a->is_linked) {              // ✅ 只保留 is_linked == false
	//						result.push_back(a);
	//					}
	//				}
	//			}
	//		}
	//	}
	//}


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

	size_t K = 5000; // 可调参数
	///
	for (size_t i = 0; i < result.size(); ++i) {
		double idy = static_cast<float>(result[i]->aligned_base) / result[i]->alignment_length;
		double score_i = result[i]->alignment_length * pow(idy, 2);
		if (score_i > 1000000 || result[i]->cigar.size() == 0) {
			continue;
		}
		if (filter_ref == true && id != result[i]->ref_chr_index)
		{
			continue;
		}
		if (filter_ref == false && id != result[i]->qry_chr_index)
		{
			continue;
		}
		dp[i] = score_i;


		// 只回看最近 K 个 j
		size_t j_start = (i > K ? i - K : 0);
		//size_t j_start = 0;
		for (size_t j = j_start; j < i; ++j) {
			bool non_overlap;
			if (filter_ref) {
				non_overlap = (result[j]->ref_start + result[j]->ref_len <= result[i]->ref_start);
			}
			else {
				non_overlap = (result[j]->qry_start + result[j]->qry_len <= result[i]->qry_start);
			}
			
			if (non_overlap && dp[j] + score_i > dp[i]) {
				dp[i] = dp[j] + score_i;
				pre[i] = j;
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
}

void PairRareAligner::filterAnchorByDP(AnchorPtrVecByStrandByQueryByRefPtr anchor_map) {

	ThreadPool pool(thread_num);
	// 获得ref和query的染色体个数，根据AnchorVecByStrandByQueryByRefPtr的维度
	uint_t qry_num = anchor_map->at(0).size();
	uint_t ref_num = anchor_map->at(0)[0].size();

	for (uint_t i = 0; i < qry_num; i++) {
		pool.enqueue([&anchor_map, i]() {
				filterChrByDP(anchor_map, i, false);
				}
		);
	}
	pool.waitAllTasksDone();
	for (uint_t i = 0; i < ref_num; i++) {
		pool.enqueue([&anchor_map, i]() {
			filterChrByDP(anchor_map, i, true);
			}
		);
	}

	pool.waitAllTasksDone();
	return;

}

void PairRareAligner::constructGraphByDP(SpeciesName query_name, SeqPro::ManagerVariant& query_seqpro_manager, AnchorPtrVecByStrandByQueryByRefPtr anchor_ptr, RaMesh::RaMeshMultiGenomeGraph& graph) {
	for (size_t strand_idx = 0; strand_idx < anchor_ptr->size(); ++strand_idx) {
		const auto& queries = anchor_ptr->at(strand_idx);
		for (size_t qry_idx = 0; qry_idx < queries.size(); ++qry_idx) {
			const auto& refs = queries.at(qry_idx);
			for (size_t ref_idx = 0; ref_idx < refs.size(); ++ref_idx) {
				const auto& anchors = refs.at(ref_idx);
				for (size_t a_idx = 0; a_idx < anchors.size(); ++a_idx) {
					AnchorPtr anchor = anchors.at(a_idx);
					// if (anchor->ref_selected && anchor->qry_selected)
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
	}
}






