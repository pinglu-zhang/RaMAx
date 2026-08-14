#include "index.h"
#include <sdsl/io.hpp>

// ------------------------------------------------------------------
// FM_Index：FM-index（基于 BWT + wavelet tree + 采样 SA）实现
// 说明：
// - fasta_manager 是一个 ManagerVariant（可能为 SequenceManager 或 MaskedSequenceManager）
// - buildIndex() 会根据序列规模/fast_mode 选择不同构建策略
// - findAnchors* 系列基于 backward search 寻找锚点（MUM/MEM）
// - ref_global_cache 用于加速 global 坐标到染色体/局部坐标的定位
// ------------------------------------------------------------------

// 构造函数：接收物种名、FastaManager 变体引用、采样率
FM_Index::FM_Index(SpeciesName species_name,
                   SeqPro::ManagerVariant& fasta_manager,
                   uint_t sample_rate)
    : species_name(species_name), fasta_manager(fasta_manager) {

    // 设置 SA 采样率（sample_rate 越大，sampled_sa 越小，定位 SA 需要更多 LF 步）
    this->sample_rate = sample_rate;

    // 若总序列长度小于采样率，则降采样率为 1，保证后续索引合法
    // 注意：total_size 在构造时可能尚未更新，但这里保留原逻辑不改变
    if (total_size < sample_rate) {
        this->sample_rate = 1;
    }
}

// ------------------------------------------------------------------
// 使用 CaPS 算法构建索引（适用于较大的序列）
// - 构建后缀数组 SA
// - 构建 BWT
// - 构建采样 SA（sampled_sa）
// - 构建 wavelet tree（wt_bwt）
// ------------------------------------------------------------------
bool FM_Index::buildIndexUsingCaPS(uint_t thread_count) {
    // 将所有序列拼接成一个大字符串 T（使用分隔符 '\1'）
    // - SequenceManager：concatAllSequences
    // - MaskedSequenceManager：concatAllSequencesSeparated（带分隔符/遮蔽相关处理）
    std::string T = std::visit([](auto&& manager_ptr) -> std::string {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if (!manager_ptr) {
            throw std::runtime_error("Manager pointer is null inside variant.");
        }
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
            return manager_ptr->concatAllSequences('\1');
        }
        else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager> >) {
            return manager_ptr->concatAllSequencesSeparated('\1');
        }
        else {
            throw std::runtime_error("Unhandled manager type in variant.");
        }
    }, fasta_manager);

    size_t n = T.size();
    this->total_size = n;

    // 若总序列长度小于采样率，则降采样率为 1，保证后续索引合法
    if (total_size < sample_rate) {
        this->sample_rate = 1;
    }

    if (n == 0)
        return false;

    // 根据序列长度选择 32 位/64 位索引类型，避免溢出并提升性能
    if (n <= std::numeric_limits<uint32_t>::max()) {
        return buildIndexUsingCaPSImpl<uint32_t>(T, thread_count);
    }
    else {
        return buildIndexUsingCaPSImpl<uint64_t>(T, thread_count);
    }
}

// ------------------------------------------------------------------
// 模板实现：CaPS 构建 SA + BWT + 采样 SA + wavelet tree
// index_t 通常为 uint32_t 或 uint64_t
// ------------------------------------------------------------------
template<typename index_t>
bool FM_Index::buildIndexUsingCaPSImpl(const std::string& T,
                                      uint_t thread_count) {
    using SA_t = index_t;

    // 这里 n = |T| + 1，通常用于后缀数组构建包含终止符的处理（保持原逻辑）
    size_t n = T.size() + 1;

    // 设置 parlay 线程数（CaPS 并行构建使用）
    std::string value = std::to_string(thread_count);
    setenv("PARLAY_NUM_THREADS", value.c_str(), 1);

    // 构造后缀数组（CaPS）
    CaPS_SA::Suffix_Array<SA_t> suf_arr(T.c_str(), static_cast<SA_t>(n), 0, 0);
    suf_arr.construct();

    // 拷贝 SA 到 std::vector
    std::vector<SA_t> SA(n);
    std::copy(suf_arr.SA(), suf_arr.SA() + n, SA.begin());

    // ------------------------------
    // 构建 BWT
    // BWT[i] = T[SA[i] - 1]（循环处理）
    // ------------------------------
    sdsl::int_vector<8> BWT(n);

    ThreadPool pool(thread_count);
    size_t chunk_size = (n + thread_count - 1) / thread_count;

    for (uint_t t = 0; t < thread_count; ++t) {
        size_t start = t * chunk_size;
        size_t end = std::min(start + chunk_size, n);

        pool.enqueue([&T, &SA, &BWT, n, start, end]() {
            for (size_t i = start; i < end; ++i) {
                size_t si = static_cast<size_t>(SA[i]);
                // (si + n - 1) % n：实现 SA[i]==0 时回绕
                BWT[i] = static_cast<uint8_t>(T[(si + n - 1) % n]);
            }
        });
    }
    pool.waitAllTasksDone();

    // ------------------------------
    // 构建采样 SA：每 sample_rate 个 SA 采一次
    // sampled_sa[idx] = SA[i]
    // ------------------------------
    size_t sampled_count = (n + sample_rate - 1) / sample_rate;
    sampled_sa.resize(sampled_count);

    for (size_t i = 0, idx = 0; i < n; i += sample_rate, ++idx) {
        sampled_sa[idx] = static_cast<uint64_t>(SA[i]);
    }
    sdsl::util::bit_compress(sampled_sa); // 压缩以节省空间

    // ------------------------------
    // 构建 wavelet tree：支持 rank/select 查询（FM-index 核心）
    // ------------------------------
    sdsl::construct_im(this->wt_bwt, BWT);

    return true;
}

// ------------------------------------------------------------------
// 使用 divsufsort 构建索引（适用于中小序列）
// - divsufsort：经典后缀数组构建
// - 后续同样构建 BWT / sampled SA / wavelet tree
// ------------------------------------------------------------------
bool FM_Index::buildIndexUsingDivsufsort(uint_t thread_count) {
    // 将所有序列拼接成一个大字符串 T（使用分隔符 '\1'）
    std::string T = std::visit([](auto&& manager_ptr) -> std::string {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if (!manager_ptr) {
            throw std::runtime_error("Manager pointer is null inside variant.");
        }
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
            return manager_ptr->concatAllSequences('\1');
        }
        else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager> >) {
            return manager_ptr->concatAllSequencesSeparated('\1');
        }
        else {
            throw std::runtime_error("Unhandled manager type in variant.");
        }
    }, fasta_manager);

    size_t n = T.size();
    this->total_size = n;

    // 若总序列长度小于采样率，则降采样率为 1，保证后续索引合法
    if (total_size < sample_rate) {
        this->sample_rate = 1;
    }

    if (n == 0)
        return false;

    // divsufsort 需要 sauchar_t 输入
    const auto* Tptr = reinterpret_cast<const sauchar_t*>(T.data());
    std::vector<saidx_t> SA(n);

    // 构建后缀数组
    if (divsufsort(Tptr, SA.data(), static_cast<saidx_t>(n)) < 0) {
        spdlog::error("divsufsort() failed");
        return false;
    }

    // ------------------------------
    // 构建 BWT（并行）
    // ------------------------------
    sdsl::int_vector<8> BWT(n);

    ThreadPool pool(thread_count);
    size_t chunk_size = (n + thread_count - 1) / thread_count;

    for (uint_t t = 0; t < thread_count; ++t) {
        size_t start = t * chunk_size;
        size_t end = std::min(start + chunk_size, n);

        pool.enqueue([&T, &SA, &BWT, n, start, end]() {
            for (size_t i = start; i < end; ++i) {
                saidx_t si = SA[i];
                BWT[i] = static_cast<uint8_t>(T[(static_cast<size_t>(si) + n - 1) % n]);
            }
        });
    }
    pool.waitAllTasksDone();

    // ------------------------------
    // 构建采样 SA
    // ------------------------------
    size_t sampled_count = (n + sample_rate - 1) / sample_rate;
    sampled_sa.resize(sampled_count);

    for (size_t i = 0, idx = 0; i < n; i += sample_rate, ++idx) {
        sampled_sa[idx] = SA[i];
    }
    sdsl::util::bit_compress(sampled_sa);

    // 构建 wavelet tree
    sdsl::construct_im(this->wt_bwt, BWT);

    return true;
}

// ------------------------------------------------------------------
// 总调度函数：根据 fast_mode 与输入规模选择索引构建方式
// - 很小文件：divsufsort
// - 其它：fast_mode -> CaPS；否则 -> BigBWT（未实现）
// - 构建完成后：初始化字母表映射（count_array / char2idx）
// ------------------------------------------------------------------
bool FM_Index::buildIndex(FilePath output_path, bool fast_mode, uint_t thread) {
    // 提取 fasta 路径字符串（SequenceManager / MaskedSequenceManager 的获取方式不同）
    std::string fasta_path_str;
    std::visit([&fasta_path_str](auto&& manager_ptr) {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
            fasta_path_str = manager_ptr->getFastaPath();
        }
        else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager> >) {
            fasta_path_str = manager_ptr->getOriginalManager().getFastaPath();
        }
        else {
            throw std::runtime_error("Unhandled manager type in variant.");
        }
    }, fasta_manager);

    // 根据文件大小选择构建策略
    if (isFileSmallerThan(fasta_path_str, 1024)) {
        // 小文件：使用 divsufsort
        buildIndexUsingDivsufsort(thread);
    }
    else {
        if (fast_mode) {
            // fast_mode：使用 CaPS
            // TODO 对齐 divsufsort
            buildIndexUsingCaPS(thread);
        }
        else {
            // 大文件慢构建模式：BigBWT（未实现）
            buildIndexUsingBigBWT(output_path, thread);
        }
    }

    // ---------- 1. 判断是否需要把 'N' 纳入字母表 ----------
    // 若参考序列中存在 ambiguous bases，则采用包含 'N' 的字母表
    bool has_ambiguous_bases = std::visit([](auto&& manager_ptr) -> bool {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>)
            return manager_ptr->hasAmbiguousBasesAll();
        else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>)
            return manager_ptr->getOriginalManager().hasAmbiguousBasesAll();
        else
            throw std::runtime_error("Unhandled manager type in variant.");
    }, fasta_manager);

    // ---------- 2. 选定实际使用的字母集合 ----------
    const std::span<const char> alph =
        has_ambiguous_bases
            ? std::span<const char>(alpha_set.data(), alpha_set.size())
            : std::span<const char>(alpha_set_without_N.data(), alpha_set_without_N.size());

    // ---------- 3. 初始化数据结构 ----------
    size_t cumulative = 0;   // C[c] 的前缀计数游标（小于 c 的字符总数）
    char2idx.fill(0xFF);     // 0xFF 表示“未映射”
    count_array.fill(0);     // C 表初始化清零

    // ---------- 4. 依序遍历 alph，填充 C 表和 char→rank 映射 ----------
    uint_t rank = 0;         // 连续字符 rank：0,1,2,...
    for (char c : alph) {
        // occ：BWT 中字符 c 出现次数（rank(size,c)）
        size_t occ = wt_bwt.rank(wt_bwt.size(), c);

        // C[c]：小于 c 的所有字符出现总数
        count_array[static_cast<uint8_t>(c)] = cumulative;
        cumulative += occ;

        // char2idx：字符 → rank（用于 wavelet tree 查询/编码）
        char2idx[static_cast<uint8_t>(c)] = rank++;
    }

    return true;
}

// ------------------------------------------------------------------
// LF 映射：给定 BWT 位置 pos，返回上一列对应位置
// LF(pos) = C[c] + Occ(c, pos)
// ------------------------------------------------------------------
uint_t FM_Index::LF(uint_t pos) const {
    uint8_t c = wt_bwt[pos];
    return count_array[c] + wt_bwt.rank(pos, c);
}

// ------------------------------------------------------------------
// 根据位置 pos 逆向跳转恢复 SA 值（利用 sampled_sa）
// - 若 pos 不是采样点，则反复 LF 直到采样点，并累计 steps
// - SA = sampled_sa[pos/sample_rate] + steps（超出 total_size 需回绕）
// ------------------------------------------------------------------
uint_t FM_Index::getSA(size_t pos) const {
    uint_t steps = 0;
    size_t cur = pos;

    // 跳到最近的采样点（cur % sample_rate == 0）
    while ((cur) % sample_rate != 0) {
        cur = LF(cur);
        ++steps;
    }

    size_t sa_sample = sampled_sa[cur / sample_rate];
    size_t sa_val = sa_sample + steps;

    // SA 回绕修正
    if (sa_val >= total_size)
        sa_val -= total_size;

    return sa_val;
}

// ------------------------------------------------------------------
// Backward Extend：Backward Search 核心
// 给定 SA 区间 I=[l,r) 与字符 c，扩展得到新的区间
// new_l = C[c] + Occ(c, l)
// new_r = C[c] + Occ(c, r)
// ------------------------------------------------------------------
SAInterval FM_Index::backwardExtend(const SAInterval& I, char c) {
    uint8_t ch = static_cast<uint8_t>(c);
    uint_t new_l = count_array[ch] + wt_bwt.rank(I.l, ch);
    uint_t new_r = count_array[ch] + wt_bwt.rank(I.r, ch);
    return { new_l, new_r };
}

// ------------------------------------------------------------------
// 根据查询模式调度 anchor 搜索函数
// - FAST_SEARCH / MIDDLE_SEARCH / ACCURATE_SEARCH
// ------------------------------------------------------------------
MatchVec2DPtr FM_Index::findAnchors(ChrIndex query_chr_index, std::string query,
    SearchMode search_mode, Strand strand,
    bool allow_MEM, uint_t query_offset,
    uint_t min_anchor_length,
    bool allow_short_mum,
    uint_t max_anchor_frequency,
    sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval) {

    if (search_mode == FAST_SEARCH) {
        return findAnchorsFast(query_chr_index, query, strand, allow_MEM, query_offset,
            min_anchor_length, allow_short_mum, max_anchor_frequency, ref_global_cache, sampling_interval);
    }
    else if (search_mode == ACCURATE_SEARCH) {
        return findAnchorsAccurate(query_chr_index, query, strand, allow_MEM,
            query_offset, min_anchor_length, allow_short_mum,
            max_anchor_frequency, ref_global_cache, sampling_interval);
    }
    else if (search_mode == MIDDLE_SEARCH) {
        return findAnchorsMiddle(query_chr_index, query, strand, allow_MEM, query_offset,
            min_anchor_length, allow_short_mum, max_anchor_frequency, ref_global_cache, sampling_interval);
    }
    else {
        throw std::invalid_argument("Invalid search mode");
    }
}

// ------------------------------------------------------------------
// Middle 模式：逐段匹配 query 子串
// 说明：
// - FORWARD：reverse(query) 以适配 backward search
// - 其它链：对碱基做互补（保持原逻辑）
// - 循环中每次调用 findSubSeqAnchors 获取当前段最长匹配
// ------------------------------------------------------------------
MatchVec2DPtr FM_Index::findAnchorsMiddle(ChrIndex query_chr_index, std::string query, Strand strand, bool allow_MEM,
    uint_t query_offset, uint_t min_anchor_length, bool allow_short_mum, uint_t max_anchor_frequency,
    sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval) {

    // query 反向或互补（FM 索引默认 backward search）
    if (strand == Strand::FORWARD) {
        std::reverse(query.begin(), query.end());
    }
    else {
        //std::reverse(query.begin(), query.end());
        for (char& ch : query) {
            ch = BASE_COMPLEMENT[static_cast<unsigned char>(ch)];
        }
    }

    MatchVec2DPtr anchor_ptr_list_vec = std::make_shared<MatchVec2D>();

    uint_t total_length = 0;
    uint_t query_length = query.length();

    // 主循环：从前往后逐段查找 anchor
    while (total_length < query_length) {
        RegionVec region_vec;

        uint_t match_length = findSubSeqAnchors(
            query.c_str() + total_length, query_length - total_length, allow_MEM,
            region_vec, min_anchor_length, allow_short_mum, max_anchor_frequency,
            ref_global_cache, sampling_interval);

        if (!region_vec.empty()) {
            Coord_t qry_start;
            if (strand == Strand::FORWARD) {
                // FORWARD 模式：由于 query 被 reverse，需换算回原坐标
                qry_start = query_length - match_length - total_length + query_offset;
            }
            else {
                qry_start = total_length + query_offset;
            }

            MatchVec anchor_ptr_list;
            anchor_ptr_list.reserve(region_vec.size());

            for (uint_t i = 0; i < region_vec.size(); i++) {
                Match match(region_vec[i].chr_index, region_vec[i].start,
                            query_chr_index, qry_start, match_length, strand);
                anchor_ptr_list.push_back(match);
            }

            anchor_ptr_list_vec->push_back(anchor_ptr_list);
        }

        // 按 match_length 推进（保持原逻辑）
        total_length += match_length;
    }

    anchor_ptr_list_vec->shrink_to_fit();
    return anchor_ptr_list_vec;
}

// ------------------------------------------------------------------
// Fast 模式：跳跃推进，减少重复匹配
// 说明：
// - 与 Middle 类似，但通过 last_pos 做简单去重
// - 推进策略：total_length += min(min_anchor_length, match_length==1?min_anchor_length:match_length)
// ------------------------------------------------------------------
MatchVec2DPtr FM_Index::findAnchorsFast(ChrIndex query_chr_index, std::string query,
    Strand strand, bool allow_MEM,
    uint_t query_offset,
    uint_t min_anchor_length,
    bool allow_short_mum,
    uint_t max_anchor_frequency,
    sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval) {

    if (strand == Strand::FORWARD) {
        std::reverse(query.begin(), query.end());
    }
    else {
        for (char& ch : query) {
            ch = BASE_COMPLEMENT[static_cast<unsigned char>(ch)];
        }
    }

    MatchVec2DPtr anchor_ptr_list_vec = std::make_shared<MatchVec2D>();

    uint_t total_length = 0;
    uint_t query_length = query.length();
    uint_t last_pos = 0;

    while (total_length < query_length) {
        RegionVec region_vec;

        uint_t match_length = findSubSeqAnchors(
            query.c_str() + total_length, query_length - total_length, allow_MEM,
            region_vec, min_anchor_length, allow_short_mum, max_anchor_frequency,
            ref_global_cache, sampling_interval);

        Coord_t qry_start;
        if (strand == Strand::FORWARD) {
            qry_start = query_length - match_length - total_length + query_offset;
        }
        else {
            qry_start = total_length + query_offset;
        }

        if (!region_vec.empty()) {
            MatchVec anchor_ptr_list;
            anchor_ptr_list.reserve(region_vec.size());

            uint_t ref_end_pos = region_vec[0].start + match_length;

            // 仅当与上一段区域不重复时才添加
            if (ref_end_pos != last_pos) {
                for (uint_t i = 0; i < region_vec.size(); i++) {
                    Match match(region_vec[i].chr_index, region_vec[i].start,
                                query_chr_index, qry_start, match_length, strand);
                    anchor_ptr_list.push_back(match);
                }
                if (!anchor_ptr_list.empty())
                    anchor_ptr_list_vec->push_back(anchor_ptr_list);
            }
            last_pos = ref_end_pos;
        }

        // 快速模式推进（保持原逻辑）
        total_length += std::min(min_anchor_length,
                                 match_length == 1 ? min_anchor_length : match_length);
    }

    anchor_ptr_list_vec->shrink_to_fit();
    return anchor_ptr_list_vec;
}

// ------------------------------------------------------------------
// Accurate 模式：更细粒度推进（每次推进 1）
// 说明：
// - 同样使用 last_pos 做简单去重
// - 推进策略：total_length += 1
// ------------------------------------------------------------------
MatchVec2DPtr FM_Index::findAnchorsAccurate(ChrIndex query_chr_index, std::string query,
    Strand strand, bool allow_MEM,
    uint_t query_offset,
    uint_t min_anchor_length,
    bool allow_short_mum,
    uint_t max_anchor_frequency,
    sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval) {

    if (strand == Strand::FORWARD) {
        std::reverse(query.begin(), query.end());
    }
    else {
        for (char& ch : query) {
            ch = BASE_COMPLEMENT[static_cast<unsigned char>(ch)];
        }
    }

    MatchVec2DPtr anchor_ptr_list_vec = std::make_shared<MatchVec2D>();

    uint_t total_length = 0;
    uint_t query_length = query.length();
    uint_t last_pos = 0;

    while (total_length < query_length) {
        RegionVec region_vec;

        uint_t match_length = findSubSeqAnchors(
            query.c_str() + total_length, query_length - total_length, allow_MEM,
            region_vec, min_anchor_length, allow_short_mum, max_anchor_frequency,
            ref_global_cache, sampling_interval);

        Coord_t qry_start;
        if (strand == Strand::FORWARD) {
            qry_start = query_length - match_length - total_length + query_offset;
        }
        else {
            qry_start = total_length + query_offset;
        }

        if (!region_vec.empty()) {
            MatchVec anchor_ptr_list;
            anchor_ptr_list.reserve(region_vec.size());

            uint_t ref_end_pos = region_vec[0].start;

            // 仅当与上一段区域不重复时才添加
            if (ref_end_pos != last_pos) {
                for (uint_t i = 0; i < region_vec.size(); i++) {
                    Match match(region_vec[i].chr_index, region_vec[i].start,
                                query_chr_index, qry_start, match_length, strand);
                    anchor_ptr_list.push_back(match);
                }
                if (!anchor_ptr_list.empty())
                    anchor_ptr_list_vec->push_back(anchor_ptr_list);
            }
            last_pos = ref_end_pos;
        }

        // 精确模式推进（保持原逻辑）
        total_length += 1;
    }

    anchor_ptr_list_vec->shrink_to_fit();
    return anchor_ptr_list_vec;
}

// -------------------------------------------------------------
// 递归二分：把 (left.pos, right.pos) 区间彻底搜索干净
// -------------------------------------------------------------
void FM_Index::bisectAnchors(const std::string& query, ChrIndex query_chr_index,
    Strand strand, bool allow_MEM, uint_t query_offset,
    uint_t query_length, uint_t min_len, bool allow_short_mum,
    uint_t max_freq, const MUMInfo& left,
    const MUMInfo& right, MatchVec2D& out,
    sdsl::int_vector<0>& ref_global_cache,
    uint_t sampling_interval) {

    if (right.pos <= left.pos + 1)
        return; // 区间不足 1bp

    uint_t mid = left.pos + (right.pos - left.pos) / 2;

    RegionVec regs;
    uint_t mid_len = findSubSeqAnchors(query.c_str() + mid, query_length - mid,
        allow_MEM, regs, min_len, allow_short_mum, max_freq, ref_global_cache, sampling_interval);

    bool mid_is_mum = (regs.size() == 1);
    bool same_as_left = (left.pos + left.len == mid + mid_len);
    bool same_as_right = (right.pos + right.len == mid + mid_len);

    // 如果 mid_is_mum 为 false 写入；如果为 true，则要求与左右都不同才写入（保持原逻辑）
    if (((regs.size() == 1 || allow_MEM) && !mid_is_mum) ||
        (!same_as_left && !same_as_right)) {

        Coord_t qry_start;
        if (strand == Strand::FORWARD) {
            qry_start = query_length - mid_len - mid + query_offset;
        }
        else {
            qry_start = query_offset + mid;
        }

        MatchVec lst;
        lst.reserve(regs.size());
        for (auto const& rg : regs) {
            Match match(rg.chr_index, rg.start, query_chr_index, qry_start, mid_len, strand);
            lst.emplace_back(match);
        }
        if (!lst.empty())
            out.emplace_back(std::move(lst));
    }

    // 递归左侧
    if (!(left.is_mum && same_as_left)) {
        bisectAnchors(query, query_chr_index, strand, allow_MEM, query_offset,
            query_length, min_len, allow_short_mum, max_freq, left,
            { mid, mid_len, mid_is_mum }, out, ref_global_cache, sampling_interval);
    }

    // 递归右侧
    if (!(right.is_mum && same_as_right)) {
        bisectAnchors(query, query_chr_index, strand, allow_MEM, query_offset,
            query_length, min_len, allow_short_mum, max_freq, { mid, mid_len, mid_is_mum },
            right, out, ref_global_cache, sampling_interval);
    }
}

uint_t FM_Index::findSubSeqAnchors(const char* query, uint_t query_length,
    bool allow_MEM, RegionVec& region_vec,
    uint_t min_anchor_length,
    bool allow_short_mum,
    uint_t max_anchor_frequency,
    sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval) {

    uint_t match_length = 0;
    SAInterval I = { 0, total_size - 1 };
    SAInterval next_I = { 0, total_size - 1 };

    // backward search：不断向前扩展区间
    while (match_length < query_length) {
        if (query[match_length] == 'N' || query[match_length] == 'n') break;
        next_I = backwardExtend(I, query[match_length]);
        if (next_I.l == next_I.r) break;
        match_length++;
        I = next_I;
    }

    uint_t frequency = I.r - I.l;

    // 频率过滤：过高则直接返回 1（保持原逻辑）
    if (frequency > max_anchor_frequency) {
        return 1;
    }
    // 若不允许 MEM，则只接受唯一匹配（frequency==1）
    if (allow_MEM == false && frequency > 1) {
        return 1;
    }

    // 长度过滤：小于最短锚长度则返回 1
    if (match_length < min_anchor_length) {
        return 1;
    }

    if (true) {
        region_vec.reserve(region_vec.size() + frequency);

        // 使用传入 cache 的采样间隔参数定位染色体 id
        auto cache_size = ref_global_cache.size();

        for (uint_t i = I.l; i < I.r; i++) {
            uint_t ref_global_pos = getSA(i);

            std::visit([&](auto&& manager_ptr) {
                SeqPro::SequenceId seq_id = SeqPro::SequenceIndex::INVALID_ID;
                SeqPro::Position local_pos = 0;

                // 先尝试用 ref_global_cache 快速命中候选序列
                if (cache_size > 0) {
                    auto cache_index = ref_global_pos / sampling_interval + 1;

                    if (cache_index < cache_size) {
                        auto candidate_seq_id = ref_global_cache[cache_index];

                        if (candidate_seq_id != SeqPro::SequenceIndex::INVALID_ID) {
                            // 获取候选序列信息并验证 global 坐标落在该序列范围内
                            using PtrType = std::decay_t<decltype(manager_ptr)>;
                            const SeqPro::SequenceInfo* candidate_info = nullptr;

                            if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
                                candidate_info = manager_ptr->getIndex().getSequenceInfo(candidate_seq_id);
                            }
                            else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager> >) {
                                candidate_info = manager_ptr->getOriginalManager().getIndex().getSequenceInfo(candidate_seq_id);
                            }

                            // 快速验证：global 坐标是否落在候选序列 masked 区间内
                            if (candidate_info &&
                                ref_global_pos >= candidate_info->masked_global_start_pos &&
                                ref_global_pos < candidate_info->masked_global_start_pos + candidate_info->masked_length) {

                                // 缓存命中：设置 seq_id
                                seq_id = candidate_seq_id;

                                if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
                                    // SequenceManager：直接 globalToLocal
                                    auto [id, pos] = manager_ptr->globalToLocal(ref_global_pos);
                                    local_pos = pos;
                                    region_vec.emplace_back(id, local_pos, match_length);
                                }
                                else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager> >) {
                                    // MaskedSequenceManager：先转 masked local，再映射回 original position
                                    local_pos = ref_global_pos - candidate_info->masked_global_start_pos;
                                    local_pos = manager_ptr->toOriginalPositionSeparated(seq_id, local_pos);
                                    region_vec.emplace_back(seq_id, local_pos, match_length);
                                }
                            }
                        }
                    }
                }

                // 若缓存未命中：回退到原始定位方法（可能包含二分/索引定位）
                if (seq_id == SeqPro::SequenceIndex::INVALID_ID) {
                    using PtrType = std::decay_t<decltype(manager_ptr)>;

                    if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
                        auto [fallback_seq_id, fallback_local_pos] = manager_ptr->globalToLocal(ref_global_pos);
                        region_vec.emplace_back(fallback_seq_id, fallback_local_pos, match_length);
                    }
                    else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager> >) {
                        auto [fallback_seq_id, fallback_local_pos] = manager_ptr->globalToLocalSeparated(ref_global_pos);
                        region_vec.emplace_back(fallback_seq_id, fallback_local_pos, match_length);
                    }
                }

                // 若依旧无效则返回（保持原逻辑结构）
                if (seq_id == SeqPro::SequenceIndex::INVALID_ID) {
                    return;
                }
            }, fasta_manager);
        }
    }

    return match_length;
}

// ------------------------------------------------------------------
// 序列化保存 FM-index 到文件
// - 基本字段写入一个主文件（binary）
// - sampled_sa / wt_bwt 使用 SDSL 的 store_to_file 分别存成 .sa / .wt
// ------------------------------------------------------------------
bool FM_Index::saveToFile(const std::string& filename) const {
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) return false;

    // 保存基本信息
    ofs.write(reinterpret_cast<const char*>(&sample_rate), sizeof(sample_rate));
    ofs.write(reinterpret_cast<const char*>(&total_size), sizeof(total_size));
    ofs.write(reinterpret_cast<const char*>(&alpha_set), sizeof(alpha_set));
    ofs.write(reinterpret_cast<const char*>(&alpha_set_without_N), sizeof(alpha_set_without_N));
    ofs.write(reinterpret_cast<const char*>(&count_array), sizeof(count_array));
    ofs.write(reinterpret_cast<const char*>(&char2idx), sizeof(char2idx));

    // 使用 SDSL 的 store_to_file 分别写 SA / WT
    std::string sa_file = filename + ".sa";
    std::string wt_file = filename + ".wt";

    if (!sdsl::store_to_file(sampled_sa, sa_file)) {
        return false;
    }
    if (!sdsl::store_to_file(wt_bwt, wt_file)) {
        return false;
    }

    return static_cast<bool>(ofs);
}

// ------------------------------------------------------------------
// 从文件加载 FM-index
// - 读取主文件中的基本信息
// - sampled_sa / wt_bwt 使用 SDSL 的 load_from_file 读取 .sa / .wt
// ------------------------------------------------------------------
bool FM_Index::loadFromFile(const std::string& filename) {
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) return false;

    // 加载基本信息
    ifs.read(reinterpret_cast<char*>(&sample_rate), sizeof(sample_rate));
    ifs.read(reinterpret_cast<char*>(&total_size), sizeof(total_size));
    ifs.read(reinterpret_cast<char*>(&alpha_set), sizeof(alpha_set));
    ifs.read(reinterpret_cast<char*>(&alpha_set_without_N), sizeof(alpha_set_without_N));
    ifs.read(reinterpret_cast<char*>(&count_array), sizeof(count_array));
    ifs.read(reinterpret_cast<char*>(&char2idx), sizeof(char2idx));

    // 使用 SDSL 的 load_from_file 读取 SA / WT
    std::string sa_file = filename + ".sa";
    std::string wt_file = filename + ".wt";

    if (!sdsl::load_from_file(sampled_sa, sa_file)) {
        return false;
    }
    if (!sdsl::load_from_file(wt_bwt, wt_file)) {
        return false;
    }

    return static_cast<bool>(ifs);
}
