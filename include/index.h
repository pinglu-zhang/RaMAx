#ifndef INDEX_H
#define INDEX_H

// -------------------------
// 依赖头文件
// -------------------------
#include "anchor.h"
#include "bwtparse.hpp"
#include "data_process.h"
#include "newscan.hpp"
#include "pfbwt.hpp"

extern "C" {
#include "gsa/gsacak_wrapper.h" // C语言实现的后缀数组构造器
}

#include "CaPS-SA/Suffix_Array.hpp" // 并行后缀数组构造器（CaPS）
#include "sdsl/int_vector.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/util.hpp"
#include "sdsl/wt_huff.hpp"

#include "SeqPro.h"
#include "divsufsort.h"
#include "divsufsort64.h"
#include <cstdlib>
#include <sstream>
#include <variant>

// -------------------------
// 常量定义
// -------------------------
#define WINDOW_SIZE 10u          // 默认滑动窗口大小
#define STOP_MODULUS 100u        // 停止匹配阈值（用于 MEM/MUM）
#define FMINDEX_EXTESION "fmidx" // 索引文件后缀名

// -------------------------
// 类型定义
// -------------------------
using SampledSA_t = sdsl::int_vector<0>; // 压缩采样 SA
using WtHuff_t =
sdsl::wt_huff<sdsl::bit_vector>; // BWT 的波形树表示（Huffman 编码）

// -------------------------
// 后缀数组区间结构（闭区间）
// -------------------------
struct SAInterval {
    uint_t l{ 0 }; // 左端点
    uint_t r{ 0 }; // 右端点（不包含）

    bool empty() const noexcept { return l == r; }
    bool isUnique() const noexcept { return (r - l) == 1; }
};

// -------------------------
// 枚举：锚点查找模式
// -------------------------
enum SearchMode {
    FAST_SEARCH,    // 快速模式：一次一段，跳跃推进
    MIDDLE_SEARCH,  // 中速模式：允许一定重复
    ACCURATE_SEARCH // 精确匹配：全局最优
};

// 将枚举转为字符串表示
inline const char* SearchModeToString(SearchMode mode) {
    switch (mode) {
    case FAST_SEARCH:
        return "fast";
    case MIDDLE_SEARCH:
        return "middle";
    case ACCURATE_SEARCH:
        return "accurate";
    default:
        return "unknown";
    }
}

// -------------------------
// FM-Index 主类
// -------------------------
class FM_Index {
public:
    // -------------------------
    // 构造与初始化
    // -------------------------
    FM_Index() = default;
    FM_Index(SpeciesName species_name, SeqPro::ManagerVariant& fasta_manager,
        uint_t sample_rate = 32);

    // -------------------------
    // 索引构建函数
    // -------------------------
    bool buildIndex(FilePath output_path, bool fast_mode, uint_t thread);

    // 使用外部工具构建大型 BWT
    bool buildIndexUsingBigBWT(const FilePath& output_path, uint_t thread);

    // 使用 divsufsort 构建索引（适合小数据）
    bool buildIndexUsingDivsufsort(uint_t thread);

    // 使用 CaPS 构建索引（并行后缀数组）
    bool buildIndexUsingCaPS(uint_t thread);

    // 模板化的 CaPS 实现（适配 32/64 位）
    template <typename index_t>
    bool buildIndexUsingCaPSImpl(const std::string& T, uint_t thread_count);

    // -------------------------
    // BWT/SA 从文件加载（快速恢复索引）
    // -------------------------
    bool read_and_build_sampled_sa(const FilePath& sa_file_path);
    bool read_and_build_bwt(const FilePath& bwt_file_path);

    // -------------------------
    // 高级构建方案（性能优化）- 实验性
    // -------------------------
    bool newScan(const FilePath& fasta_path, const FilePath& output_path,
        uint_t thread);
    bool bwtParse(const FilePath& fasta_path, const FilePath& output_path,
        uint_t thread);
    bool pfBWT(const FilePath& fasta_path, const FilePath& output_path,
        uint_t thread);

    // -------------------------
    // SA / BWT 操作接口
    // -------------------------
    BWTParse::sa_index_t* compute_SA(uint32_t* Text, long n,
        long k); // 计算 SA（实验性）
    uint_t getSA(uint_t pos) const;           // 从 BWT 位置反向求原始位置
    uint_t LF(uint_t pos) const;              // LF 映射函数
    SAInterval backwardExtend(const SAInterval& I, char c); // 扩展后缀区间

    MatchVec2DPtr findAnchors(ChrIndex query_chr_index, std::string query, SearchMode search_mode, Strand strand, bool allow_MEM, uint_t query_offset, uint_t min_anchor_length, bool allow_short_mum, uint_t max_anchor_frequency, sdsl::int_vector<0>& ref_global_cache, SeqPro::Length sampling_interval);

    MatchVec2DPtr findAnchorsFast(ChrIndex query_chr_index, std::string query,
        Strand strand, bool allow_MEM,
        uint_t query_offset, uint_t min_anchor_length,
        bool allow_short_mum,
        uint_t max_anchor_frequency,
        sdsl::int_vector<0>& ref_global_cache,
        SeqPro::Length sampling_interval);

    MatchVec2DPtr findAnchorsAccurate(ChrIndex query_chr_index, std::string query,
        Strand strand, bool allow_MEM,
        uint_t query_offset,
        uint_t min_anchor_length,
        bool allow_short_mum,
        uint_t max_anchor_frequency,
        sdsl::int_vector<0>& ref_global_cache,
        SeqPro::Length sampling_interval);

    MatchVec2DPtr findAnchorsMiddle(ChrIndex query_chr_index, std::string query,
        Strand strand, bool allow_MEM,
        uint_t query_offset, uint_t min_anchor_length,
        bool allow_short_mum,
        uint_t max_anchor_frequency,
        sdsl::int_vector<0>& ref_global_cache,
        SeqPro::Length sampling_interval);

    // -------------------------
    // MEM / MUM 区间拆分辅助结构与算法
    // -------------------------
    struct MUMInfo {
        uint_t pos;  // query 上的起点
        uint_t len;  // 匹配长度
        bool is_mum; // true=MUM, false=MEM
    };

    void bisectAnchors(const std::string& query, ChrIndex query_chr_index, Strand strand,
        bool allow_MEM, uint_t query_offset, uint_t query_length,
        uint_t min_anchor_length, bool allow_short_mum, uint_t max_anchor_frequency,
        const MUMInfo& left, const MUMInfo& right,
        MatchVec2D& out,
        sdsl::int_vector<0>& ref_global_cache,
        SeqPro::Length sampling_interval);

    // -------------------------
    // 子串查找（返回匹配长度与区域）
    // -------------------------
    // uint_t findSubSeqAnchorsFast(
    //    const char* query,
    //    uint_t query_length,
    //    bool allow_MEM,
    //    RegionVec& region_vec,
    //    uint_t min_anchor_length,
    //    uint_t max_anchor_frequency);

    uint_t findSubSeqAnchors(const char* query, uint_t query_length,
        bool allow_MEM, RegionVec& region_vec,
        uint_t min_anchor_length,
        bool allow_short_mum,
        uint_t max_anchor_frequency,
        sdsl::int_vector<0>& ref_global_cache,
        SeqPro::Length sampling_interval);

    // -------------------------
    // 索引序列化 / 反序列化
    // -------------------------
    bool saveToFile(const std::string& filename) const;
    bool loadFromFile(const std::string& filename);

    /* bool encode_kmer(const char* kmer, uint_t length, uint64_t code);

     void build_kmer_table(uint_t k, uint_t thread_count);*/

     // -------- 工具函数 --------
    template <size_t N>
    static std::array<char, N> repositionNullAfter(const std::array<char, N>& arr,
        char c) {
        std::vector<char> temp;
        temp.reserve(N);
        for (char ch : arr) {
            if (ch != '\0') {
                temp.push_back(ch);
            }
        }
        auto it = std::find(temp.begin(), temp.end(), c);
        if (it == temp.end()) {
            throw std::runtime_error("Character not found in array");
        }
        size_t pos = std::distance(temp.begin(), it);
        temp.insert(temp.begin() + pos + 1, '\0');

        if (temp.size() != N) {
            throw std::runtime_error("Unexpected size mismatch after reordering");
        }
        std::array<char, N> result;
        std::copy(temp.begin(), temp.end(), result.begin());
        return result;
    }

    // 完全使用SDSL原生序列化，避免与cereal混合使用

    SpeciesName species_name;
    SeqPro::ManagerVariant& fasta_manager;
    std::array<char, 7> alpha_set{ '\0', '\1', 'A', 'C', 'G', 'N', 'T' };
    std::array<char, 6> alpha_set_without_N{ '\0', '\1', 'A', 'C', 'G', 'T' };
    SampledSA_t sampled_sa;
    WtHuff_t wt_bwt;
    std::array<uint_t, 256> count_array{};
    uint_t sample_rate{ 32 };
    uint_t total_size{ 0 };
    std::array<uint8_t, 256> char2idx{};
};

#endif
