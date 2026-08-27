#ifndef ANCHOR_H
#define ANCHOR_H

// ------------------------------------------------------------------
// 引入必要的头文件
// ------------------------------------------------------------------
#include "config.hpp"              // 包含通用配置（如类型定义）
#include "align.h"
#include <memory>                 // 用于智能指针（如 std::shared_ptr）
#include <vector>                 // 用于 std::vector 容器
#include <SeqPro.h>
#include <xxhash.h>


#define ANCHOR_EXTENSION "anchor"  // Anchor 文件保存使用的默认扩展名

// ------------------------------------------------------------------
// 类型别名定义
// ------------------------------------------------------------------
// 定义坐标类型（可在 config.hpp 中设为 uint32_t 或 uint64_t）
using Coord_t = uint32_t;
using Length_t = uint32_t;

/* 区间与 map 类型 */
using IntervalMap = std::map<int_t, int_t>;          // key = beg, value = end   (右开)

/* 插入非重叠区间到 map（调用处已保证不重叠，只需 emplace 即可） */
inline void insertInterval(IntervalMap& m, int_t l, int_t r) { m.emplace(l, r); }

/* 判断 [l,r) 是否与 map 中已有区间重叠；若重叠，把那条区间的 [L,R) 返到 outBox */
inline bool overlap1D(const IntervalMap& m,
    int_t l, int_t r,
    int_t& L, int_t& R)
{
    if (m.empty()) return false;

    auto it = m.lower_bound(l);          // 第一个 beg ≥ l
    if (it != m.end() && it->first < r) {  // 与右邻区间重叠
        L = it->first; R = it->second;
        return true;
    }
    if (it != m.begin()) {                 // 与左邻区间重叠？
        --it;
        if (it->second > l) {
            L = it->first; R = it->second;
            return true;
        }
    }
    return false;
}


// ------------------------------------------------------------------
// 区域 Region 与 比对 Match 的结构定义
// ------------------------------------------------------------------

// 表示基因组上的一个区域
struct Region {
    ChrIndex chr_index{};     // 染色体名称
    Coord_t start{ 0 };     // 起始坐标
    Coord_t length{ 0 };    // 区域长度

    Region() = default;
    Region(ChrIndex chr_index, Coord_t s, Coord_t len)
        : chr_index{ chr_index }, start{ s }, length{ len } {
    }
};

using RegionVec = std::vector<Region>;  // 区域集合

// 链接方向：正向或反向
enum Strand { FORWARD, REVERSE };

inline uint64_t make_sqr_key(uint32_t s, uint32_t q, uint32_t r) noexcept {
    return (uint64_t(s & 1u) << 62) | (uint64_t(q & 0x7FFFFFFFu) << 31) | uint64_t(r & 0x7FFFFFFFu);
}
inline uint32_t key_s(uint64_t k) noexcept { return uint32_t((k >> 62) & 1u); }
inline uint32_t key_q(uint64_t k) noexcept { return uint32_t((k >> 31) & 0x7FFFFFFFu); }
inline uint32_t key_r(uint64_t k) noexcept { return uint32_t(k & 0x7FFFFFFFu); }

// -------------------- 2) xxHash hasher --------------------
struct XXH64Hasher {
    size_t operator()(uint64_t k) const noexcept {
        return static_cast<size_t>(XXH64(&k, sizeof(k), 0));
    }
};

// 表示参考与查询序列之间的一段匹配区域
struct Match {
    ChrIndex ref_chr_index;
    Coord_t  ref_start;
    ChrIndex qry_chr_index;
    Coord_t  qry_start;
    uint32_t len_and_strand; // 高位存 strand，低 31 位存 length

    Match() = default;

    Match(ChrIndex r_chr, Coord_t r_start,
        ChrIndex q_chr, Coord_t q_start,
        uint32_t len, Strand strand)
        : ref_chr_index(r_chr),
        ref_start(r_start),
        qry_chr_index(q_chr),
        qry_start(q_start),
        len_and_strand((len & 0x7FFFFFFF) | (static_cast<uint32_t>(strand) << 31)) {
    }

    uint32_t match_len() const {
        return len_and_strand & 0x7FFFFFFF;
    }

    Strand strand() const {
        return (len_and_strand >> 31) ? Strand::REVERSE : Strand::FORWARD;
    }

    void set_match_len(uint32_t len) {
        len_and_strand = (len & 0x7FFFFFFF) | (static_cast<uint32_t>(strand()) << 31);
    }

    void set_strand(Strand s) {
        len_and_strand = (match_len() & 0x7FFFFFFF) | (static_cast<uint32_t>(s) << 31);
    }
};

using MatchVec = std::vector<Match>;
using MatchPtr = std::shared_ptr<Match>; // 智能指针类型
using MatchPtrVec = std::vector<MatchPtr>; // 智能指针集合
using MatchBySQR_Sparse = std::unordered_map<uint64_t, MatchVec, XXH64Hasher>;
using MatchBySQR_SparsePtr = std::shared_ptr<MatchBySQR_Sparse>;

using MatchVec2D = std::vector<MatchVec>;
using MatchVec2DPtr = std::shared_ptr<MatchVec2D>;
using MatchVec3D = std::vector<std::vector<MatchVec>>;
using MatchVec3DPtr = std::shared_ptr<MatchVec3D>;

using SpeciesMatchVec3DPtrMap = std::unordered_map<std::string, MatchVec3DPtr>;
using SpeciesMatchVec3DPtrMapPtr = std::shared_ptr<SpeciesMatchVec3DPtrMap>;

using MatchByRef = std::vector<MatchVec>;
using MatchByQueryRef = std::vector<MatchByRef>;
using MatchByStrandByQueryRef = std::vector<MatchByQueryRef>;
using MatchByStrandByQueryRefPtr = std::shared_ptr<MatchByStrandByQueryRef>;

using SpeciesMatchByStrandByQueryRefPtrMap = std::unordered_map<SpeciesName, MatchByStrandByQueryRefPtr>;

using MatchCluster = MatchVec; // 匹配簇，包含多个匹配向量
using MatchClusterVec = std::vector<MatchCluster>;
using MatchClusterVecPtr = std::shared_ptr<MatchClusterVec>;

using ClusterSynteny = std::vector<MatchCluster>;
using ClusterSyntenyVec = std::vector<ClusterSynteny>;
using ClusterSyntenyVecPtr = std::shared_ptr<ClusterSyntenyVec>;
using ClusterSyntenyVecPtrByRef = std::vector<ClusterSyntenyVecPtr>;

using ClusterBySQR_Sparse = std::unordered_map<uint64_t, MatchClusterVecPtr, XXH64Hasher>;
using ClusterBySQR_SparsePtr = std::shared_ptr<ClusterBySQR_Sparse>;

using ClusterVecPtrByRef = std::vector<MatchClusterVecPtr>;
using ClusterVecPtrByQueryRef = std::vector<ClusterVecPtrByRef>;
using ClusterVecPtrByStrandByQueryRef = std::vector<ClusterVecPtrByQueryRef>;

using ClusterVecPtrByRefPtr = std::shared_ptr<ClusterVecPtrByRef>;
using ClusterVecPtrByRefQuery = std::vector<ClusterVecPtrByRef>;
using ClusterVecPtrByRefQueryPtr = std::shared_ptr<ClusterVecPtrByRefQuery>;
using ClusterVecPtrByStrandByQueryRefPtr = std::shared_ptr<ClusterVecPtrByStrandByQueryRef>;

using SpeciesClusterMap =
std::unordered_map<SpeciesName, ClusterBySQR_SparsePtr>;
using SpeciesClusterMapPtr = std::shared_ptr<SpeciesClusterMap>;

inline Coord_t start1(const Match& m) { return (m.ref_start); }
inline Coord_t start2(const Match& m) { return (m.qry_start); }
inline Length_t len1(const Match& m) { return (m.match_len()); }
inline Length_t len2(const Match& m) { return (m.match_len()); }
inline int_t diag(const Match& m) {
    return start2(m) - start1(m);
}

inline int_t diag_reverse(const Match& m) {
    return -start2(m) - start1(m);
}

// 3. 彻底释放 cluster 内剩余元素的内存（swap 技巧）
inline void releaseCluster(MatchVec& cluster) {
    MatchVec tmp; tmp.swap(cluster);
}

/* ---------- 辅助：计算簇跨度 ---------- */
inline int_t clusterSpan(const MatchCluster& cl)
{
    return start1(cl.back()) + len1(cl.back())
        - start1(cl.front());
}

void groupMatchByQueryRef(MatchVec3DPtr& anchors,
    MatchByStrandByQueryRefPtr unique_anchors,
    MatchByStrandByQueryRefPtr repeat_anchors,
    SeqPro::ManagerVariant& ref_fasta_manager,
    SeqPro::ManagerVariant& query_fasta_manager);

void groupMatchByQueryRefSparse(
    MatchVec3DPtr& anchors,
    MatchBySQR_SparsePtr unique_anchors,
    MatchBySQR_SparsePtr repeat_anchors,
    SeqPro::ManagerVariant& ref_fasta_manager,
    SeqPro::ManagerVariant& query_fasta_manager);

ClusterVecPtrByRefPtr
groupClustersByRef(const ClusterVecPtrByStrandByQueryRefPtr& src);

MatchClusterVecPtr clusterChrMatch(MatchVec& unique_match, uint_t min_cluster_length, int_t max_gap = 90, int_t diagdiff = 5, double diagfactor = 0.12);

MatchVec bestChainDP(MatchVec& cluster, double diagfactor);

MatchClusterVec buildClusters(MatchVec& unique_match,
    int_t      max_gap,
    int_t      diagdiff,
    double     diagfactor);

ClusterVecPtrByStrandByQueryRefPtr
clusterAllChrMatch(const MatchByStrandByQueryRefPtr& unique_anchors,
    const MatchByStrandByQueryRefPtr& repeat_anchors, uint_t min_span);

ClusterBySQR_SparsePtr clusterAllChrMatchSparse(
    MatchBySQR_SparsePtr& unique_anchors,
    MatchBySQR_SparsePtr& repeat_anchors,  // 按你的要求：占位不用
    uint_t min_span,
    uint_t thread_num,
    bool include_repeats);


void validateClusters(const ClusterVecPtrByStrandByQueryRefPtr& cluster_vec_ptr);

struct Anchor {
    ChrIndex ref_chr_index;
    Coord_t  ref_start;
    Length_t ref_len;
    ChrIndex qry_chr_index;
    Coord_t  qry_start;
    Length_t qry_len;

    Strand strand;
    uint_t alignment_length{}; // 对齐长度
    uint_t aligned_base{};
    Cigar_t cigar{};                 // 对齐的 CIGAR 字符串

    bool ref_selected = false;
    bool qry_selected = false;
    bool is_linked = false;

    Anchor() = default;

    Anchor(ChrIndex ref_chr, Coord_t ref_start, Length_t ref_len,
        ChrIndex qry_chr, Coord_t qry_start, Length_t qry_len,
        Strand strand, uint_t align_len, uint_t aligned_base, Cigar_t cigar_str)
        : ref_chr_index(ref_chr),
        ref_start(ref_start),
        ref_len(ref_len),
        qry_chr_index(qry_chr),
        qry_start(qry_start),
        qry_len(qry_len),
        strand(strand),
        alignment_length(align_len),
        aligned_base(aligned_base),
        cigar(cigar_str) {
    }


};

using AnchorPtr = std::shared_ptr<Anchor>;
using AnchorPtrVec = std::vector<AnchorPtr>;
using AnchorVec = std::vector<Anchor>;
using AnchorPtrVecByStrandByQueryByRef = std::vector<std::vector<std::vector<AnchorPtrVec>>>;
using AnchorPtrVecByStrandByQueryByRefPtr = std::shared_ptr<AnchorPtrVecByStrandByQueryByRef>;

using AnchorBySQR_Sparse = std::vector<AnchorPtrVec>;
using AnchorBySQR_SparsePtr = std::shared_ptr<AnchorBySQR_Sparse>;

using AnchorPtrVecByQueryByRef = std::vector<std::vector<AnchorPtrVec>>;
using AnchorPtrVecByQueryByRefPtr = std::shared_ptr<AnchorPtrVecByQueryByRef>;
AnchorVec extendClusterToAnchorVec(const MatchCluster& cluster,
    const SeqPro::ManagerVariant& ref_mgr,
    const SeqPro::ManagerVariant& query_mgr);

Anchor extendClusterToAnchor(MatchCluster& cluster,
    const SeqPro::ManagerVariant& ref_mgr,
    const SeqPro::ManagerVariant& query_mgr);

// ------------------------------------------------------------------
// 智能序列分块函数
// ------------------------------------------------------------------

/**
 * @brief 智能预分配序列分块，支持自动分块策略选择
 * 
 * 该函数会根据数据量大小智能选择分块策略：
 * - 当序列数量较多时（>= max_sequences_threshold），按序列划分而不是按长度切分
 * - 当序列数量较少但总长度很大时，按指定大小切分序列并保留重叠区域
 * 
 * @param seq_manager SeqPro::Manager对象，用于获取序列信息
 * @param chunk_size 每个chunk的目标大小（碱基数）
 * @param overlap_size chunk之间的重叠大小（碱基数）
 * @param max_sequences_threshold 序列数量阈值，超过此值时按序列划分
 * @param min_chunk_size 最小chunk大小，小于此值的序列不会被进一步切分
 * @return RegionVec 分块结果，每个Region代表一个chunk
 */
RegionVec preAllocateChunks(const SeqPro::ManagerVariant& seq_manager,
                           uint_t chunk_size,
                           uint_t overlap_size = 0,
                           size_t max_sequences_threshold = 1000,
                           uint_t min_chunk_size = 10000);

/**
 * @brief 按序列划分的简单分块策略
 * 
 * 当序列数量较多时使用，每个序列作为一个独立的chunk
 * 
 * @param seq_manager SeqPro::Manager对象
 * @return RegionVec 分块结果，每个序列一个chunk
 */
RegionVec preAllocateChunksBySequence(const SeqPro::ManagerVariant& seq_manager);

/**
 * @brief 按大小划分的分块策略（原始FastaManager逻辑）
 * 
 * 当序列数量较少但总长度很大时使用，会对长序列进行切分
 * 
 * @param seq_manager SeqPro::Manager对象
 * @param chunk_size 每个chunk的目标大小
 * @param overlap_size chunk之间的重叠大小
 * @param min_chunk_size 最小chunk大小，小于此值的序列不会被切分
 * @return RegionVec 分块结果
 */
RegionVec preAllocateChunksBySize(const SeqPro::ManagerVariant& seq_manager,
                                 uint_t chunk_size,
                                 uint_t overlap_size = 0,
                                 uint_t min_chunk_size = 10000,
                                 bool isMultiple = false);

/**
 * @brief 按原始坐标切分全部序列，包括 MaskedSequenceManager 的遮蔽区间
 */
RegionVec preAllocateOriginalChunksBySize(
    const SeqPro::ManagerVariant& seq_manager,
    uint_t chunk_size,
    uint_t overlap_size = 0,
    uint_t min_chunk_size = 10000);

/**
 * @brief 常规的基于大小的分块辅助函数
 * @param chunks 输出的分块结果
 * @param seq_name 序列名称
 * @param seq_length 序列长度
 * @param chunk_size 分块大小
 * @param overlap_size 重叠大小
 * @param min_chunk_size 最小分块大小
 */
void normalSizeBasedChunking(RegionVec& chunks, const std::string& seq_name, 
                           uint64_t seq_length, uint_t chunk_size, 
                           uint_t overlap_size, uint_t min_chunk_size, const SeqPro::ManagerVariant& seq_manager);

/**
 * @brief 基于遮蔽区间的预分割函数
 * @param chunks 输出的分块结果
 * @param seq_name 序列名称
 * @param seq_length 序列长度
 * @param masked_manager 遮蔽序列管理器
 * @param chunk_size 分块大小
 * @param overlap_size 重叠大小
 */
void splitByMaskedRegions(RegionVec& chunks, const std::string& seq_name, 
                        uint64_t seq_length, 
                        const SeqPro::MaskedSequenceManager& masked_manager,
                        uint_t chunk_size, uint_t overlap_size, const SeqPro::ManagerVariant& seq_manager);

/**
 * @brief 对非遮蔽区域进行分块
 * @param chunks 输出的分块结果
 * @param seq_name 序列名称
 * @param region_start 区域起始位置
 * @param region_end 区域结束位置
 * @param chunk_size 分块大小
 * @param overlap_size 重叠大小
 */
void chunkUnmaskedRegion(RegionVec& chunks, const std::string& seq_name,
                       uint64_t region_start, uint64_t region_end,
                       uint_t chunk_size, uint_t overlap_size, const SeqPro::ManagerVariant& seq_manager);

// ------------------------------------------------------------------
// Cereal 序列化支持：用于将结构写入文件或从文件读取
// ------------------------------------------------------------------
namespace cereal {

    //// Region 的序列化
    //template <class Archive>
    //void serialize(Archive& ar, Region& r) {
    //    ar(r.chr_name, r.start, r.length);
    //}

    // Match 的序列化
    //template <class Archive>
    //void serialize(Archive& ar, Match& m) {
    //    ar(m.ref_region, m.query_region, m.strand);
    //}

    // Anchor 的序列化
    //template <class Archive>
    //void serialize(Archive& ar, Anchor& a) {
    //    ar(a.match,
    //        a.alignment_length,
    //        a.cigar);
    //}

} // namespace cereal

// ------------------------------------------------------------------
// Anchor 文件的读写接口
// ------------------------------------------------------------------

//bool saveMatchVec3D(const std::string& filename, const MatchVec3DPtr& data);
//bool loadMatchVec3D(const std::string& filename, MatchVec3DPtr& data);
//
//bool loadSpeciesMatchMap(const std::string& filename, SpeciesMatchVec3DPtrMapPtr& map_ptr);
//bool saveSpeciesMatchMap(const std::string& filename, const SpeciesMatchVec3DPtrMapPtr& map_ptr);

class UnionFind
{
public:
    /// 构造并初始化为 n 个独立元素
    explicit UnionFind(std::size_t n = 0);

    /// 清空并重建大小为 n 的并查集
    void reset(std::size_t n);

    /// 查找 x 所在集合的代表（根）
    int_t find(int_t x);

    /// 返回 x 所在集合的元素个数
    int_t set_size(int_t x);

    /// 判断 a 与 b 是否位于同一集合
    bool same(int_t a, int_t b);

    /// 当前独立集合个数
    int_t components() const { return component_cnt_; }

    /**
     * @brief 合并 a、b 所在集合
     * @return true  如果二者原本不连通并成功合并
     *         false 如果二者已在同一集合
     */
    bool unite(int_t a, int_t b);

private:
    std::vector<int_t> parent_;   // 负数：根 & -size；非负：父索引
    int_t component_cnt_{ 0 };      // 当前集合数
};

// ------------------------------------------------------------------
// Anchor 验证功能
// ------------------------------------------------------------------

/**
 * @brief 验证 anchors 结果的序列匹配正确性
 * 
 * 通过从原始序列中提取匹配区域并比较序列内容来验证anchor的正确性。
 * 支持反向互补序列的验证。
 * 
 * @param anchors 待验证的三维匹配向量
 * @param ref_manager 参考序列管理器
 * @param query_manager 查询序列管理器
 * @return 包含验证统计信息的结构体
 */
struct ValidationResult {
    uint64_t total_matches = 0;
    uint64_t correct_matches = 0;
    uint64_t incorrect_matches = 0;
    
    double accuracy() const {
        return total_matches ? (100.0 * correct_matches / total_matches) : 0.0;
    }
};

ValidationResult validateAnchorsCorrectness(
    const MatchVec3DPtr& anchors,
    const SeqPro::ManagerVariant& ref_manager,
    const SeqPro::ManagerVariant& query_manager
);

// 重载版本：支持 SharedManagerVariant
ValidationResult validateAnchorsCorrectness(
    const MatchVec3DPtr& anchors,
    const SeqPro::SharedManagerVariant& ref_manager,
    const SeqPro::SharedManagerVariant& query_manager
);


void linkClusters(AnchorPtrVec& anchors,
    const SeqPro::ManagerVariant& ref_mgr,
    const SeqPro::ManagerVariant& qry_mgr);

AnchorPtrVec linkClusters(MatchClusterVec& clusters,
                          const SeqPro::ManagerVariant& ref_mgr,
                          const SeqPro::ManagerVariant& qry_mgr);
#endif // ANCHOR_H
