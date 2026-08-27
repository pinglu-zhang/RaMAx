#include "anchor.h"
#include <SeqPro.h>
#include "data_process.h"

namespace {
constexpr Coord_t kMaximumLinkedGapAlignmentLength = 10000;

constexpr bool linkedGapCanBeAligned(Coord_t ref_gap, Coord_t query_gap) {
    return ref_gap <= kMaximumLinkedGapAlignmentLength &&
           query_gap <= kMaximumLinkedGapAlignmentLength;
}

static_assert(linkedGapCanBeAligned(9999, 9999));
static_assert(linkedGapCanBeAligned(10000, 10000));
static_assert(!linkedGapCanBeAligned(10001, 10000));
static_assert(!linkedGapCanBeAligned(10000, 10001));

struct ManagedSequenceSlice {
    std::string storage;
    KswSequenceView view;
};

const SeqPro::SequenceManager& originalSequenceManager(
    const SeqPro::ManagerVariant& manager) {
    return std::visit([](const auto& pointer) -> const SeqPro::SequenceManager& {
        using Pointer = std::decay_t<decltype(pointer)>;
        if constexpr (std::is_same_v<
                          Pointer,
                          std::unique_ptr<SeqPro::SequenceManager>>) {
            return *pointer;
        } else {
            return pointer->getOriginalManager();
        }
    }, manager);
}

void loadSequenceSlice(const SeqPro::ManagerVariant& manager,
                       ChrIndex chromosome, Coord_t start, Coord_t length,
                       bool reverse_complement,
                       ManagedSequenceSlice& output) {
    output.storage.clear();
    if (length == 0) {
        output.view = {{}, reverse_complement};
        return;
    }
    const auto& original = originalSequenceManager(manager);
    std::string_view contiguous;
    if (original.tryGetContiguousSubSequence(
            chromosome, start, length, contiguous)) {
        output.view = {contiguous, reverse_complement};
        return;
    }
    original.getSubSequenceInto(
        chromosome, start, length, output.storage);
    output.view = {output.storage, reverse_complement};
}

void appendCigarMove(Cigar_t& destination, Cigar_t& source) {
    if (source.empty()) return;
    if (destination.empty()) {
        destination = std::move(source);
        Cigar_t().swap(source);
        return;
    }
    destination.reserve(destination.size() + source.size());
    appendCigar(destination, source);
    Cigar_t().swap(source);
}
}  // namespace

// 稀疏版：不再 resize 出 [strand][qry][ref] 的满矩阵
void groupMatchByQueryRefSparse(
    MatchVec3DPtr& anchors,
    MatchBySQR_SparsePtr unique_anchors,
    MatchBySQR_SparsePtr repeat_anchors,
    SeqPro::ManagerVariant& /*ref_fasta_manager*/,
    SeqPro::ManagerVariant& /*query_fasta_manager*/)
{
    if (!anchors || anchors->empty()) return;

    // 经验性 reserve，避免大量 rehash
    // 你也可以按 anchors 的总 slice/vec 数粗略估算
    unique_anchors->reserve(1024);
    repeat_anchors->reserve(1024);

    for (auto& slice : *anchors) {
        if (slice.empty()) continue;

        for (auto& vec : slice) {
            if (vec.empty()) continue;

            ChrIndex rIdx = vec.front().ref_chr_index;
            ChrIndex qIdx = vec.front().qry_chr_index;
            const uint32_t sIdx = (vec.front().strand() == REVERSE ? 1u : 0u);

            if (rIdx == SeqPro::SequenceIndex::INVALID_ID) continue;
            if (qIdx == SeqPro::SequenceIndex::INVALID_ID) continue;

            const uint64_t key = make_sqr_key(
                sIdx,
                static_cast<uint32_t>(qIdx),
                static_cast<uint32_t>(rIdx)
            );

            // 如果你仍要区分 unique/repeat，可用 vec.size()==1 的规则
            MatchBySQR_SparsePtr& tgt = (vec.size() == 1) ? unique_anchors : repeat_anchors;
            auto& dest = (*tgt)[key];

            if (dest.empty()) dest.reserve(vec.size());
            dest.insert(dest.end(),
                        std::make_move_iterator(vec.begin()),
                        std::make_move_iterator(vec.end()));

            vec.clear();
            vec.shrink_to_fit();
        }
    }

    anchors->clear();
    anchors->shrink_to_fit();
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
    const Match& first = cluster.front();

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
    ManagedSequenceSlice reference_slice;
    ManagedSequenceSlice query_slice;
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
        Coord_t query_gap_begin = 0;
        Coord_t query_gap_length = 0;
        if (fwd) {
            query_gap_begin = qry_start;
            query_gap_length = qry_end - qry_start;
        }
        else {
            query_gap_begin = qry_end + len2;
            query_gap_length = m.qry_start - len2 - qry_end;
        }

        loadSequenceSlice(ref_mgr, ref_chr, ref_start,
            ref_end - ref_start, false, reference_slice);
        loadSequenceSlice(query_mgr, qry_chr, query_gap_begin,
            query_gap_length, !fwd, query_slice);
        AlignmentResult gap = globalAlignKSW2Result(
            reference_slice.view, query_slice.view);
        match_len += gap.summary.match_length;
        aln_len += gap.summary.alignment_length;
        appendCigarMove(cigar, gap.cigar);

    }
	const Match& last = cluster.back();
    if (fwd) {
        anchor = Anchor(ref_chr, first.ref_start, last.ref_start + last.match_len() - first.ref_start, qry_chr, first.qry_start, last.qry_start + last.match_len() - first.qry_start, strand, aln_len, match_len, std::move(cigar));
    }
    else {
        anchor = Anchor(ref_chr, first.ref_start, last.ref_start + last.match_len() - first.ref_start, qry_chr, last.qry_start, first.qry_start + first.match_len() - last.qry_start, strand, aln_len, match_len, std::move(cigar));
    }

    return anchor;
}

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
                                "Cluster FAILED (strand={},query={},ref={},cluster={}): "
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
                                "Cluster FAILED (strand={},query={},ref={},cluster={}): "
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
                                    "Cluster FAILED (strand={},query={},ref={},cluster={}): "
                                    "qry_start < qry_last_pos ({} < {})",
                                    strand_i, q_i, r_i, c_i, qry_start, qry_last_pos);
                            }
                            qry_last_pos = qry_end;
                        }
                        else { // REVERSE
                            // 反向：查询坐标应当递减，区间不重叠
                            if (qry_end > qry_last_pos) {
                                spdlog::debug(
                                    "Cluster FAILED (strand={},query={},ref={},cluster={}): "
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
    // if (strand == FORWARD) {
    //     std::sort(anchors.begin(), anchors.end(),
    //         [](auto& a, auto& b) { return a->qry_start < b->qry_start; });
    // }
    // else {
    //     std::sort(anchors.begin(), anchors.end(),
    //         [](auto& a, auto& b) { return a->qry_start > b->qry_start; });
    // }
    std::sort(anchors.begin(), anchors.end(),
        [](auto& a, auto& b) { return a->ref_start < b->ref_start; });


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

        bool reach = false;
        if (best != anchors.end()) {
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

            if (!linkedGapCanBeAligned(ref_gap_len, qry_gap_len)) {
                reach = false;
            } else {
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

            uint_t ref_len = 0;
            uint_t qry_len = 0;
            Cigar_t gap_cigar;
            gap_cigar = extendAlignKSW2(ref_gap_seq, qry_gap_seq, 2 * break_len);
            ref_len = countRefLength(gap_cigar);
            qry_len = countQryLength(gap_cigar);
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
            }
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

// Found the best preceding anchor.
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

                    // The gap aligns completely; append it to curr.
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

            curr_qry_end = (*curr)->qry_start + (*curr)->qry_len;
            curr_ref_end = (*curr)->ref_start + (*curr)->ref_len;
            //curr = best;
            while (curr != anchors.end()) {
                if (!(*curr)->is_linked && (*curr)->ref_start >= curr_ref_end) {
                    break;
                }

                ++curr;
            }
        }

        //++curr;
    }

    anchors.swap(linked);
    return;
}
inline bool intervalOverlap(Coord_t a_lo, Coord_t a_hi, Coord_t b_lo, Coord_t b_hi) {
    // 闭开区间 [lo,hi)，不重叠当且仅当 a_hi <= b_lo 或 b_hi <= a_lo
    return !(a_hi <= b_lo || b_hi <= a_lo);
}

namespace {

AnchorPtrVec linkClustersByValue(
    MatchClusterVec& clusters,
    const SeqPro::ManagerVariant& ref_mgr,
    const SeqPro::ManagerVariant& qry_mgr) {
    std::sort(clusters.begin(), clusters.end(),
        [](const MatchCluster& left, const MatchCluster& right) {
            return left.front().ref_start < right.front().ref_start;
        });

    MatchClusterVec cleaned;
    cleaned.reserve(clusters.size());
    bool have_previous = false;
    ChrIndex previous_ref_chromosome = 0;
    ChrIndex previous_query_chromosome = 0;
    Strand previous_strand = FORWARD;
    Coord_t previous_ref_end = 0;
    Coord_t previous_query_low = 0;
    Coord_t previous_query_high = 0;

    const auto queryBounds = [](const MatchCluster& cluster) {
        Coord_t low = std::numeric_limits<Coord_t>::max();
        Coord_t high = 0;
        for (const auto& match : cluster) {
            const Coord_t first = match.qry_start;
            const Coord_t second = match.qry_start + len2(match);
            low = std::min(low, std::min(first, second));
            high = std::max(high, std::max(first, second));
        }
        if (low == std::numeric_limits<Coord_t>::max()) low = 0;
        return std::pair<Coord_t, Coord_t>{low, high};
    };

    for (auto& cluster : clusters) {
        if (cluster.empty()) continue;
        if (!have_previous ||
            cluster.front().ref_chr_index != previous_ref_chromosome ||
            cluster.front().qry_chr_index != previous_query_chromosome ||
            cluster.front().strand() != previous_strand) {
            MatchCluster kept = std::move(cluster);
            previous_ref_chromosome = kept.front().ref_chr_index;
            previous_query_chromosome = kept.front().qry_chr_index;
            previous_strand = kept.front().strand();
            previous_ref_end = kept.back().ref_start + len1(kept.back());
            const auto [low, high] = queryBounds(kept);
            previous_query_low = low;
            previous_query_high = high;
            cleaned.push_back(std::move(kept));
            have_previous = true;
            continue;
        }

        MatchCluster pruned;
        pruned.reserve(cluster.size());
        for (auto& match : cluster) {
            const bool reference_ok = match.ref_start >= previous_ref_end;
            const Coord_t first = match.qry_start;
            const Coord_t second = match.qry_start + len2(match);
            const Coord_t low = std::min(first, second);
            const Coord_t high = std::max(first, second);
            const bool query_ok = !intervalOverlap(
                low, high, previous_query_low, previous_query_high);
            if (reference_ok && query_ok) {
                pruned.push_back(std::move(match));
            }
        }
        if (pruned.empty()) continue;
        previous_ref_end = pruned.back().ref_start + len1(pruned.back());
        const auto [low, high] = queryBounds(pruned);
        previous_query_low = low;
        previous_query_high = high;
        cleaned.push_back(std::move(pruned));
    }
    clusters.swap(cleaned);
    MatchClusterVec().swap(cleaned);

    AnchorVec anchors;
    anchors.reserve(clusters.size());
    for (auto& cluster : clusters) {
        if (cluster.empty()) continue;
        anchors.push_back(extendClusterToAnchor(cluster, ref_mgr, qry_mgr));
        anchors.back().is_linked = false;
        MatchCluster().swap(cluster);
    }
    MatchClusterVec().swap(clusters);

    AnchorPtrVec output;
    if (anchors.empty()) return output;
    std::vector<size_t> linked;
    linked.reserve(anchors.size());
    constexpr size_t kNoIndex = std::numeric_limits<size_t>::max();
    constexpr int kCandidateLimit = 2000;
    constexpr int kLookBack = 2000;
    constexpr int_t kMaximumFallbackGap = 100000;
    constexpr int_t kBreakLength = 200;
    ManagedSequenceSlice reference_slice;
    ManagedSequenceSlice query_slice;

    size_t current_index = 0;
    while (current_index < anchors.size()) {
        Anchor& current = anchors[current_index];
        if (current.is_linked) {
            ++current_index;
            continue;
        }

        size_t best_index = kNoIndex;
        int_t best_score = std::numeric_limits<int_t>::max();
        int looked = 0;
        for (size_t index = current_index + 1;
             index < anchors.size() && looked < kCandidateLimit; ++index) {
            Anchor& candidate = anchors[index];
            if (candidate.is_linked) continue;
            ++looked;
            const int_t ref_gap = static_cast<int_t>(candidate.ref_start) -
                static_cast<int_t>(current.ref_start + current.ref_len);
            const int_t query_gap = current.strand == FORWARD
                ? static_cast<int_t>(candidate.qry_start) -
                    static_cast<int_t>(current.qry_start + current.qry_len)
                : static_cast<int_t>(current.qry_start) -
                    static_cast<int_t>(candidate.qry_start + candidate.qry_len);
            if (ref_gap < 0 || query_gap < 0) continue;
            const long greater = std::max(ref_gap, query_gap);
            const long lesser = std::min(ref_gap, query_gap);
            if (greater < kBreakLength ||
                greater - lesser <= kBreakLength) {
                best_index = index;
                break;
            }
            const int_t score = (greater << 1) - lesser;
            if (best_score > score) {
                best_score = score;
                best_index = index;
            }
        }

        bool reached = false;
        if (best_index != kNoIndex) {
            Anchor& best = anchors[best_index];
            const Coord_t ref_gap_begin = current.ref_start + current.ref_len;
            const Coord_t ref_gap_length = best.ref_start - ref_gap_begin;
            Coord_t query_gap_begin = 0;
            Coord_t query_gap_length = 0;
            if (current.strand == FORWARD) {
                query_gap_begin = current.qry_start + current.qry_len;
                query_gap_length = best.qry_start - query_gap_begin;
            } else {
                query_gap_begin = best.qry_start + best.qry_len;
                query_gap_length = current.qry_start - query_gap_begin;
            }

            if (linkedGapCanBeAligned(ref_gap_length, query_gap_length)) {
                loadSequenceSlice(ref_mgr, current.ref_chr_index,
                    ref_gap_begin, ref_gap_length, false, reference_slice);
                loadSequenceSlice(qry_mgr, current.qry_chr_index,
                    query_gap_begin, query_gap_length,
                    current.strand == REVERSE, query_slice);
                AlignmentResult gap = extendAlignKSW2Result(
                    reference_slice.view, query_slice.view,
                    2 * kBreakLength);
                if (gap.summary.reference_length == ref_gap_length &&
                    gap.summary.query_length == query_gap_length) {
                    reached = true;
                    current.ref_len = best.ref_start + best.ref_len -
                        current.ref_start;
                    if (current.strand == FORWARD) {
                        current.qry_len = best.qry_start + best.qry_len -
                            current.qry_start;
                    } else {
                        current.qry_len = current.qry_start + current.qry_len -
                            best.qry_start;
                        current.qry_start = best.qry_start;
                    }
                    current.cigar.reserve(current.cigar.size() +
                        gap.cigar.size() + best.cigar.size());
                    appendCigarMove(current.cigar, gap.cigar);
                    appendCigarMove(current.cigar, best.cigar);
                    current.alignment_length += best.alignment_length +
                        gap.summary.alignment_length;
                    current.aligned_base += gap.summary.match_length +
                        best.aligned_base;
                    best.is_linked = true;
                }
            }
        }

        if (!reached) {
            if (!linked.empty()) {
                size_t best_linked_position = kNoIndex;
                int_t linked_best_score = std::numeric_limits<int_t>::max();
                int checked = 0;
                for (size_t position = linked.size();
                     position > 0 && checked < kLookBack; --position, ++checked) {
                    const Anchor& previous = anchors[linked[position - 1]];
                    if (previous.strand != current.strand ||
                        previous.ref_chr_index != current.ref_chr_index ||
                        previous.qry_chr_index != current.qry_chr_index) {
                        continue;
                    }
                    const int_t ref_gap = static_cast<int_t>(current.ref_start) -
                        static_cast<int_t>(previous.ref_start + previous.ref_len);
                    const int_t query_gap = current.strand == FORWARD
                        ? static_cast<int_t>(current.qry_start) -
                            static_cast<int_t>(previous.qry_start + previous.qry_len)
                        : static_cast<int_t>(previous.qry_start) -
                            static_cast<int_t>(current.qry_start + current.qry_len);
                    if (ref_gap < 0 || query_gap < 0 ||
                        ref_gap > kMaximumFallbackGap ||
                        query_gap > kMaximumFallbackGap) {
                        continue;
                    }
                    const long greater = std::max(ref_gap, query_gap);
                    const long lesser = std::min(ref_gap, query_gap);
                    const int_t score = (greater << 1) - lesser;
                    if (score < linked_best_score) {
                        linked_best_score = score;
                        best_linked_position = position - 1;
                    }
                }

                if (best_linked_position != kNoIndex) {
                    const Anchor& previous =
                        anchors[linked[best_linked_position]];
                    const Coord_t ref_gap_begin =
                        previous.ref_start + previous.ref_len;
                    const Coord_t ref_gap_length =
                        current.ref_start - ref_gap_begin;
                    Coord_t query_gap_begin = 0;
                    Coord_t query_gap_length = 0;
                    if (current.strand == FORWARD) {
                        query_gap_begin = previous.qry_start + previous.qry_len;
                        query_gap_length = current.qry_start - query_gap_begin;
                    } else {
                        query_gap_begin = current.qry_start + current.qry_len;
                        query_gap_length = previous.qry_start - query_gap_begin;
                    }
                    loadSequenceSlice(ref_mgr, current.ref_chr_index,
                        ref_gap_begin, ref_gap_length, false, reference_slice);
                    loadSequenceSlice(qry_mgr, current.qry_chr_index,
                        query_gap_begin, query_gap_length,
                        current.strand == REVERSE, query_slice);
                    AlignmentResult gap = extendAlignKSW2Result(
                        reference_slice.view, query_slice.view,
                        2 * kBreakLength);
                    if (gap.summary.reference_length == ref_gap_length &&
                        gap.summary.query_length == query_gap_length) {
                        current.ref_len += gap.summary.reference_length;
                        current.qry_len += gap.summary.query_length;
                        current.alignment_length +=
                            gap.summary.alignment_length;
                        current.aligned_base += gap.summary.match_length;
                        current.ref_start -= gap.summary.reference_length;
                        if (current.strand == FORWARD) {
                            current.qry_start -= gap.summary.query_length;
                        }
                        prependCigar(current.cigar, gap.cigar);
                    }
                }
            }

            linked.push_back(current_index);
            const Coord_t current_ref_end = current.ref_start + current.ref_len;
            const Coord_t current_query_end = current.strand == FORWARD
                ? current.qry_start + current.qry_len
                : current.qry_start;
            const size_t stop = best_index == kNoIndex
                ? anchors.size() : best_index;
            for (size_t index = current_index + 1; index < stop; ++index) {
                Anchor& candidate = anchors[index];
                if (candidate.is_linked) continue;
                const Coord_t candidate_ref_end =
                    candidate.ref_start + candidate.ref_len;
                const Coord_t candidate_query_end = current.strand == FORWARD
                    ? candidate.qry_start + candidate.qry_len
                    : candidate.qry_start;
                if (candidate_ref_end <= current_ref_end &&
                    ((current.strand == FORWARD &&
                      candidate_query_end <= current_query_end) ||
                     (current.strand == REVERSE &&
                      candidate_query_end >= current_query_end))) {
                    candidate.is_linked = true;
                    Cigar_t().swap(candidate.cigar);
                }
            }
            while (current_index < anchors.size()) {
                const Anchor& candidate = anchors[current_index];
                if (!candidate.is_linked &&
                    candidate.ref_start >= current_ref_end) {
                    break;
                }
                ++current_index;
            }
        }
    }

    output.reserve(linked.size());
    for (const size_t index : linked) {
        anchors[index].is_linked = false;
        output.push_back(std::make_shared<Anchor>(std::move(anchors[index])));
    }
    return output;
}

}  // namespace

// 以 cluster 为输入：内部先 extend 成 anchor，再执行你已有的 linking 流程
AnchorPtrVec linkClusters(MatchClusterVec& clusters,
                          const SeqPro::ManagerVariant& ref_mgr,
                          const SeqPro::ManagerVariant& qry_mgr)
{
    return linkClustersByValue(clusters, ref_mgr, qry_mgr);


}
