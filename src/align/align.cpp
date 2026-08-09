#include "align.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>

#include <spdlog/spdlog.h>

namespace {

std::string external_insertion_msa_executable;
std::atomic<uint64_t> external_msa_file_counter{0};
std::atomic<uint64_t> external_msa_completed{0};

std::string shellQuote(const std::string& value) {
    std::string quoted = "'";
    for (const char c : value) {
        if (c == '\'') {
            quoted += "'\\''";
        } else {
            quoted += c;
        }
    }
    quoted += "'";
    return quoted;
}

std::string ungappedUpper(std::string sequence) {
    sequence.erase(
        std::remove_if(
            sequence.begin(), sequence.end(),
            [](unsigned char c) {
                return c == '-' || std::isspace(c);
            }),
        sequence.end());
    std::transform(
        sequence.begin(), sequence.end(), sequence.begin(),
        [](unsigned char c) {
            return static_cast<char>(std::toupper(c));
        });
    return sequence;
}

class TemporaryMsaFiles {
public:
    std::filesystem::path input;
    std::filesystem::path output;

    ~TemporaryMsaFiles() {
        std::error_code error;
        if (!input.empty()) {
            std::filesystem::remove(input, error);
        }
        error.clear();
        if (!output.empty()) {
            std::filesystem::remove(output, error);
        }
    }
};

bool runExternalMsa(
    const std::string& executable,
    std::unordered_map<ChrName, std::string>& sequences) {
    if (executable.empty() ||
        sequences.size() < 2) {
        return false;
    }

    std::vector<ChrName> keys;
    keys.reserve(sequences.size());
    for (const auto& [key, unused] : sequences) {
        (void)unused;
        keys.push_back(key);
    }
    std::sort(keys.begin(), keys.end());

    const uint64_t serial =
        external_msa_file_counter.fetch_add(1, std::memory_order_relaxed);
    const uint64_t stamp = static_cast<uint64_t>(
        std::chrono::steady_clock::now().time_since_epoch().count());
    const auto prefix =
        std::filesystem::temp_directory_path() /
        ("ramax-insertion-msa-" + std::to_string(stamp) + "-" +
         std::to_string(serial));

    TemporaryMsaFiles temporary{
        prefix.string() + ".input.fa",
        prefix.string() + ".output.fa"};

    {
        std::ofstream input(temporary.input, std::ios::binary);
        if (!input) {
            throw std::runtime_error(
                "Cannot create external MSA input: " +
                temporary.input.string());
        }
        for (size_t index = 0; index < keys.size(); ++index) {
            input << ">s" << index << "\n"
                  << sequences.at(keys[index]) << "\n";
        }
    }

    const std::string command =
        shellQuote(executable) +
        " -r 1 -t 1 " + shellQuote(temporary.input.string()) +
        " > " + shellQuote(temporary.output.string());
    if (std::system(command.c_str()) != 0) {
        spdlog::warn(
            "[external-msa] command failed: {}",
            executable);
        return false;
    }

    std::ifstream output(temporary.output, std::ios::binary);
    if (!output) {
        spdlog::warn(
            "[external-msa] output is unavailable");
        return false;
    }

    std::unordered_map<std::string, std::string> aligned_by_id;
    std::string current_id;
    std::string line;
    while (std::getline(output, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (!line.empty() && line.front() == '>') {
            std::istringstream header(line.substr(1));
            header >> current_id;
            if (current_id.empty() ||
                !aligned_by_id.emplace(current_id, std::string{}).second) {
                spdlog::warn(
                    "[external-msa] invalid FASTA header");
                return false;
            }
            continue;
        }
        if (current_id.empty()) {
            if (line.empty()) {
                continue;
            }
            spdlog::warn(
                "[external-msa] sequence before FASTA header");
            return false;
        }
        for (const unsigned char c : line) {
            if (!std::isspace(c)) {
                aligned_by_id[current_id].push_back(
                    static_cast<char>(c));
            }
        }
    }

    if (aligned_by_id.size() != keys.size()) {
        spdlog::warn("[external-msa] output row count mismatch");
        return false;
    }

    size_t aligned_length = 0;
    std::vector<std::string> aligned(keys.size());
    for (size_t index = 0; index < keys.size(); ++index) {
        const std::string id = "s" + std::to_string(index);
        const auto aligned_it = aligned_by_id.find(id);
        if (aligned_it == aligned_by_id.end() ||
            ungappedUpper(aligned_it->second) !=
                ungappedUpper(sequences.at(keys[index]))) {
            spdlog::warn(
                "[external-msa] output sequence validation failed");
            return false;
        }
        if (index == 0) {
            aligned_length = aligned_it->second.size();
        } else if (aligned_it->second.size() != aligned_length) {
            spdlog::warn(
                "[external-msa] output rows have unequal lengths");
            return false;
        }
        aligned[index] = aligned_it->second;
    }

    for (size_t index = 0; index < keys.size(); ++index) {
        sequences[keys[index]] = std::move(aligned[index]);
    }

    const uint64_t completed =
        external_msa_completed.fetch_add(1, std::memory_order_relaxed) + 1;
    if (completed % 1000 == 0) {
        spdlog::info(
            "[external-msa] completed {} MSAs with {}",
            completed, executable);
    }
    return true;
}

}  // namespace

bool alignSequencesWithExternalMsa(
    const std::string& executable,
    std::unordered_map<ChrName, std::string>& sequences) {
    return runExternalMsa(executable, sequences);
}

void configureExternalInsertionMsa(const std::string& executable) {
    external_insertion_msa_executable = executable;
    if (!executable.empty()) {
        spdlog::info(
            "[external-msa] insertion aligner configured: {}",
            executable);
    }
}

void InsertInfo::alignSeqs() {
    if (aligned || seqs.empty()) {
        return;
    }

    if (seqs.size() == 1) {
        ref_name = seqs.begin()->first;
        total_length = seqs.begin()->second.size();
        aligned = true;
        return;
    }

    if (runExternalMsa(
            external_insertion_msa_executable, seqs)) {
        ref_name = seqs.begin()->first;
        total_length = seqs.begin()->second.size();
        aligned = true;
        return;
    }

    size_t max_len = 0;
    ChrName longest_key;
    for (const auto& [key, sequence] : seqs) {
        if (sequence.size() > max_len) {
            max_len = sequence.size();
            longest_key = key;
        }
    }
    ref_name = longest_key;

    std::unordered_map<ChrName, Cigar_t> cigars;
    for (const auto& [key, sequence] : seqs) {
        (void)sequence;
        if (key == ref_name) {
            continue;
        }
        cigars[key] = globalAlignKSW2(seqs[ref_name], seqs.at(key));
    }
    total_length = mergeAlignmentByRef(ref_name, seqs, cigars);
    aligned = true;
}

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
    cfg.alphabet_size = 5;
    cfg.gap_open = 5;
    cfg.gap_extend = 2;
    cfg.band_width = auto_band(ref.size(), query.size());


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





