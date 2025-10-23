#include "align.h"

KSW2AlignConfig makeDefaultKSW2Config() {
    static int8_t simple_dna_mat[25];
    static bool initialized = false;
    if (!initialized) {
        // A C G T N -> 0 1 2 3 4
        const int match = 2;
        const int mismatch = -3;
        const int ambiguous = -1;  // 对N的惩罚较小

        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                if (i == 4 || j == 4) {
                    simple_dna_mat[i * 5 + j] = ambiguous;  // N匹配
                }
                else if (i == j) {
                    simple_dna_mat[i * 5 + j] = match;
                }
                else {
                    simple_dna_mat[i * 5 + j] = mismatch;
                }
            }
        }
        initialized = true;
    }

    return {
        .mat = simple_dna_mat,
        .alphabet_size = 5,
        .gap_open = 8,
        .gap_extend = 1,
        .end_bonus = 0,
        .zdrop = 100,           // 较远匹配终止
        .band_width = -1,      // 合理的band可提高性能
        .flag = KSW_EZ_GENERIC_SC | KSW_EZ_RIGHT
    };
}

Cigar_t globalAlignKSW2(const std::string& ref,
    const std::string& query)
{
    /* ---------- 1. 编码序列 ---------- */
    std::vector<uint8_t> ref_enc(ref.size());
    std::vector<uint8_t> qry_enc(query.size());

    for (size_t i = 0; i < ref.size(); ++i)
        ref_enc[i] = ScoreChar2Idx[static_cast<uint8_t>(ref[i])];
    for (size_t i = 0; i < query.size(); ++i)
        qry_enc[i] = ScoreChar2Idx[static_cast<uint8_t>(query[i])];

    /* ---------- 2. 复制 cfg 并修正常见坑 ---------- */
    // KSW2AlignConfig cfg = cfg_in;                   // 本地副本可调整
    init_simd_mat();

    KSW2AlignConfig cfg;
    cfg.mat = dna5_simd_mat;
    cfg.alphabet_size = 5;
    cfg.gap_open = 5;          // gap open penalty
    cfg.gap_extend = 2;        // gap extension penalty
    cfg.end_bonus = 0;         // ❌ 不需要 ends-free 奖励
    cfg.zdrop = -1;            // ❌ 禁用 z-drop（全局比对必须完整比完）
    cfg.band_width = -1;       // 启用全矩阵（也可设 auto_band）
    cfg.flag = KSW_EZ_RIGHT; // ✅ 通用矩阵 + gap右对齐

    /* ---------- 3. 调用 KSW2 ---------- */
    ksw_extz_t ez{};

    ksw_extz2_sse(0,
        static_cast<int>(qry_enc.size()), qry_enc.data(),
        static_cast<int>(ref_enc.size()), ref_enc.data(),
        cfg.alphabet_size, cfg.mat,
        cfg.gap_open, cfg.gap_extend,
        cfg.band_width, cfg.zdrop, cfg.end_bonus,
        cfg.flag, &ez);


    /* ---------- 4. 拷贝 / 释放 CIGAR ---------- */
    Cigar_t cigar;
    cigar.reserve(ez.n_cigar);
    for (int i = 0; i < ez.n_cigar; ++i)
        cigar.push_back(ez.cigar[i]);

    free(ez.cigar);           // KSW2 用 malloc()
    return cigar;
}

Cigar_t globalAlignKSW2_2(const std::string& ref,
    const std::string& query)
{
    /* ---------- 1. 编码序列 ---------- */
    std::vector<uint8_t> ref_enc(ref.size());
    std::vector<uint8_t> qry_enc(query.size());

    for (size_t i = 0; i < ref.size(); ++i)
        ref_enc[i] = ScoreChar2Idx[static_cast<uint8_t>(ref[i])];
    for (size_t i = 0; i < query.size(); ++i)
        qry_enc[i] = ScoreChar2Idx[static_cast<uint8_t>(query[i])];

    /* ---------- 2. 复制 cfg 并修正常见坑 ---------- */
    // KSW2AlignConfig cfg = cfg_in;                   // 本地副本可调整
    KSW2AlignConfig cfg = makeTurboKSW2Config2(query.size(), ref.size());

    /* ---------- 3. 调用 KSW2 ---------- */
    ksw_extz_t ez{};

    ksw_extz2_sse(0,
        static_cast<int>(qry_enc.size()), qry_enc.data(),
        static_cast<int>(ref_enc.size()), ref_enc.data(),
        cfg.alphabet_size, cfg.mat,
        cfg.gap_open, cfg.gap_extend,
        cfg.band_width, cfg.zdrop, cfg.end_bonus,
        cfg.flag, &ez);


    /* ---------- 4. 拷贝 / 释放 CIGAR ---------- */
    Cigar_t cigar;
    cigar.reserve(ez.n_cigar);
    for (int i = 0; i < ez.n_cigar; ++i)
        cigar.push_back(ez.cigar[i]);

    free(ez.cigar);           // KSW2 用 malloc()
    return cigar;
}

//Cigar_t globalAlignWFA2(const std::string& ref,
//    const std::string& query)
//{
//    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
//    attributes.distance_metric = gap_affine;
//    attributes.affine_penalties.mismatch = 2;      // X > 0
//    attributes.affine_penalties.gap_opening = 3;   // O >= 0
//    attributes.affine_penalties.gap_extension = 1; // E > 0
//    attributes.memory_mode = wavefront_memory_high;
//    attributes.heuristic.strategy = wf_heuristic_banded_adaptive;
//    attributes.heuristic.min_k = -50;
//    attributes.heuristic.max_k = +50;
//    attributes.heuristic.steps_between_cutoffs = 1;
//    //// Create a WFAligner
//    //
//    wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
//
//    wavefront_align(wf_aligner, ref.c_str(), ref.length(), query.c_str(), query.length());
//    /*wfa::WFAlignerGapAffine aligner(2, 3, 1, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryUltralow);
//
//    aligner.alignEnd2End(ref, query);*/
//    uint32_t* cigar_buffer; // Buffer to hold the resulting CIGAR operations.
//    int cigar_length = 0; // Length of the CIGAR string.
//    // Retrieve the CIGAR string from the wavefront aligner.
//    cigar_get_CIGAR(wf_aligner->cigar, true, &cigar_buffer, &cigar_length);
//
//    /* ---------- 4. 拷贝 / 释放 CIGAR ---------- */
//    Cigar_t cigar;
//
//    for (int i = 0; i < cigar_length; ++i)
//        cigar.push_back(cigar_buffer[i]);
//
//    wavefront_aligner_delete(wf_aligner);
//
//    return cigar;
//}
//
//Cigar_t extendAlignWFA2(const std::string& ref,
//    const std::string& query, int zdrop)
//{
//    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
//    attributes.distance_metric = gap_affine;
//    attributes.affine_penalties.mismatch = 2;      // X > 0
//    attributes.affine_penalties.gap_opening = 3;   // O >= 0
//    attributes.affine_penalties.gap_extension = 1; // E > 0
//    attributes.memory_mode = wavefront_memory_high;
//    attributes.heuristic.strategy = wf_heuristic_zdrop;
//    attributes.heuristic.zdrop = zdrop;
//    attributes.heuristic.steps_between_cutoffs = 1;
//    //// Create a WFAligner
//    //
//    wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
//
//    wavefront_align(wf_aligner, ref.c_str(), ref.length(), query.c_str(), query.length());
//    /*wfa::WFAlignerGapAffine aligner(2, 3, 1, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryUltralow);
//
//    aligner.alignEnd2End(ref, query);*/
//
//    uint32_t* cigar_buffer; // Buffer to hold the resulting CIGAR operations.
//    int cigar_length = 0; // Length of the CIGAR string.
//    // Retrieve the CIGAR string from the wavefront aligner.
//    cigar_get_CIGAR(wf_aligner->cigar, true, &cigar_buffer, &cigar_length);
//
//    /* ---------- 4. 拷贝 / 释放 CIGAR ---------- */
//    Cigar_t cigar;
//
//    for (int i = 0; i < cigar_length; ++i)
//        cigar.push_back(cigar_buffer[i]);
//
//    wavefront_aligner_delete(wf_aligner);
//
//    return cigar;
//}

/**********************************************************************
*  extendAlignKSW2  ——  ends-free（seed-and-extend）比对
*    @param ref        参考片段（目标方向）
*    @param query      查询片段（同方向；若反链请先反向互补）
*    @param zdrop      Z-drop 剪枝阈值（默认 200）
*    @param band       带宽限制；<0 表示不限制
*    @return           Cigar_t（BAM 编码）
**********************************************************************/
Cigar_t extendAlignKSW2(const std::string& ref,
    const std::string& query,
    int zdrop)
{
    /* ---------- 1. 序列编码 ---------- */
    std::vector<uint8_t> ref_enc(ref.size());
    std::vector<uint8_t> qry_enc(query.size());
    for (size_t i = 0; i < ref.size(); ++i) ref_enc[i] = ScoreChar2Idx[(uint8_t)ref[i]];
    for (size_t i = 0; i < query.size(); ++i) qry_enc[i] = ScoreChar2Idx[(uint8_t)query[i]];

    ///* ---------- 2. 配置 ---------- */
    //KSW2AlignConfig cfg = makeTurboKSW2Config(query.size(), ref.size());
    ////KSW2AlignConfig cfg;
    //cfg.zdrop = zdrop;       // 用于提前终止
    //cfg.flag = KSW_EZ_EXTZ_ONLY     // ends-free extension
    //    | KSW_EZ_APPROX_MAX    // 跟踪 ez.max_q/max_t
    //    | KSW_EZ_APPROX_DROP   // 在 approximate 模式下触发 z-drop 就中断
    //    | KSW_EZ_RIGHT;        // （可选）gap 右对齐     // **关键**：启用 extension/ends-free
    //// 若需要右对齐 gaps 建议保留 KSW_EZ_RIGHT
    //cfg.end_bonus = 100;
    //cfg.band_width = -1;
    init_simd_mat();
    KSW2AlignConfig cfg;
	cfg.mat = dna5_simd_mat;
    cfg.zdrop = zdrop;
    cfg.flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT | KSW_EZ_APPROX_DROP;
    cfg.end_bonus = 50;
    cfg.band_width = -1;
    cfg.alphabet_size = 5;
    cfg.gap_open = 5;
    cfg.gap_extend = 2;


    /* ---------- 3. 调用 KSW2 ---------- */
    ksw_extz_t ez{};
    ksw_extz2_sse(nullptr,
        static_cast<int>(qry_enc.size()), qry_enc.data(),
        static_cast<int>(ref_enc.size()), ref_enc.data(),
        cfg.alphabet_size, cfg.mat,
        cfg.gap_open, cfg.gap_extend,
        cfg.band_width, cfg.zdrop, cfg.end_bonus,
        cfg.flag, &ez);

    // 赋值bool& if_zdrop,int& ref_end,int& qry_end
    /* ---------- 4. 拷贝 & 释放 ---------- */
    Cigar_t cigar;
    cigar.reserve(ez.n_cigar);
    for (int i = 0; i < ez.n_cigar; ++i)
        cigar.push_back(ez.cigar[i]);

    free(ez.cigar);                    // ksw2 使用 malloc
    return cigar;                      // 返回的 CIGAR 即延伸片段
}

/* ──────────── 合并成 MSA (就地修改 seqs) ──────────── */
uint_t mergeAlignmentByRef(
    ChrName ref_name,
    std::unordered_map<ChrName, std::string>& seqs,
    const std::unordered_map<ChrName, Cigar_t>& cigars)
{
    auto ref_it = seqs.find(ref_name);
    if (ref_it == seqs.end())
        throw std::invalid_argument("mergeAlignmentByRef: ref not found");

    std::string& ref_raw = ref_it->second;
    uint_t total_aligned_length = ref_raw.size();
	RefAlignInfo insert_info;

    for (const auto& [key, cigar] : cigars) {
        if (key == ref_name) continue;
        auto q_it = seqs.find(key);
        if (q_it == seqs.end()) {
            throw std::invalid_argument("mergeAlignmentByRef: seq missing");
        }

		std::string& qry_raw = q_it->second;

		uint_t ref_pos = 0;
		uint_t qry_pos = 0;
        for (auto& unit : cigar) {
            uint32_t len;
            char op;
			intToCigar(unit, op, len);

            if (op == 'D') {
				qry_raw.insert(qry_pos, len, '-');
                ref_pos += len;
                qry_pos += len;
            }
            else if (op == 'I') {
                std::string ins = qry_raw.substr(qry_pos, len);

                // 2) 在 insert_info 里插入或更新
                auto it = insert_info.find(ref_pos);
                if (it != insert_info.end()) {
                    it->second.seqs[key] = ins;
                }
                else {
                    InsertInfo info;
                    info.seqs[key] = ins;
                    insert_info[ref_pos] = std::move(info);
                }

                // 3) 从原始 query 序列里移除这段已“消费”的子串
                qry_raw.erase(qry_pos, len);
            }
            else {
                ref_pos += len;
                qry_pos += len;
            }
        }
    }

    uint_t offset = 0;
	for (auto& [ref_pos, info] : insert_info) {
		info.alignSeqs(); // 对齐所有插入序列
		if (info.ref_name.empty()) continue; // 没有参考序列，跳过
		for (auto& [sp_name, seq] : seqs) {
			auto it = info.seqs.find(sp_name);
            if (it != info.seqs.end()) {
				seq.insert(ref_pos + offset, it->second); // 在 ref_pos 位置插入
            }
            else {
                seq.insert(ref_pos + offset, info.total_length, '-');   // 直接用 string::insert 重载
            }
		}
        offset += info.total_length; // 更新总长度
        total_aligned_length += info.total_length;
	}

    for (auto& [chr, seq] : seqs) {
        if (seq.size() != total_aligned_length) {
            std::cout << "";
        }
    }

    return total_aligned_length;

}

AlignCount countAlignedBases(const Cigar_t& cigar) {
    AlignCount cnt;
    for (auto op : cigar) {
        uint32_t len;
        char type;
        intToCigar(op, type, len);
        switch (type) {
        case 'M': // match or mismatch
        case '=': // match
        case 'X': // mismatch
            cnt.ref_bases += len;
            cnt.query_bases += len;
            break;
        case 'I': // insertion wrt ref
            cnt.query_bases += len;
            break;
        case 'D': // deletion wrt ref
            cnt.ref_bases += len;
            break;
            // 视情况处理 clip/skip
        case 'S': // soft clip
            cnt.query_bases += len;
            break;
        case 'H': // hard clip
            // 不计入
            break;
        case 'N': // skipped region in ref
            cnt.ref_bases += len;
            break;
        default:
            break;
        }
    }
    return cnt;
}





