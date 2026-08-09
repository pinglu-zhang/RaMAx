// =============================================================
//  File: ramesh.h   –  Public interface (v0.6‑alpha, 2025‑06‑25)
//  ✧  Header for the lock‑free multi‑genome graph  ✧
// =============================================================
#ifndef RAMESH_H
#define RAMESH_H

#include <atomic>
#include <memory>
#include <unordered_map>
#include <vector>
#include <array>
#include <map>
#include <shared_mutex>
#include <chrono>
#include <algorithm>

#include "config.hpp"
#include "anchor.h"
#include "softmask_index.h"

// 前向声明：NewickParser 位于全局命名空间（见 data_process.h）
class NewickParser;

namespace RaMesh {

    // ────────────────────────────────────────────────
    // Forward declarations
    // ────────────────────────────────────────────────
    class  RaMeshNode;
    class  Block;
    struct Segment;

    using BlockPtr = std::shared_ptr<Block>;
    using WeakBlock = std::weak_ptr<Block>;
    using SegPtr = std::shared_ptr<Segment>;

    // ────────────────────────────────────────────────
    // Helper types / hashing
    // ────────────────────────────────────────────────
    using SpeciesChrPair = std::pair<SpeciesName, ChrName>;

    struct SpeciesChrPairHash {
        std::size_t operator()(const SpeciesChrPair& k) const noexcept {
            std::size_t h1 = std::hash<std::string>{}(k.first);
            std::size_t h2 = std::hash<std::string>{}(k.second);
            return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
        }
    };

    using ChrHeadMap = std::unordered_map<SpeciesChrPair, SegPtr, SpeciesChrPairHash>;

    // ────────────────────────────────────────────────
    // Intrusive concurrent list node
    // ────────────────────────────────────────────────
    struct RaMeshPath {
        std::atomic<SegPtr> next{ nullptr };
        std::atomic<SegPtr> prev{ nullptr };
    };

    // ────────────────────────────────────────────────
    // Domain enums
    // ────────────────────────────────────────────────
    enum class AlignRole : uint8_t { PRIMARY = 0, SECONDARY = 1 };
    enum class SegmentRole : uint8_t { SEGMENT = 0, HEAD = 1, TAIL = 2 };

    // ────────────────────────────────────────────────
    // Segment / Sentinel
    // ────────────────────────────────────────────────
    class Segment {
    public:
        uint_t      start{ 0 };
        uint_t      length{ 0 };
        Cigar_t     cigar;
        Strand      strand{ Strand::FORWARD };
        AlignRole   align_role{ AlignRole::PRIMARY };
        SegmentRole seg_role{ SegmentRole::SEGMENT };

        RaMeshPath  primary_path;
        BlockPtr   parent_block;

        bool left_extend{ false };
		bool right_extend{ false };
        mutable std::shared_mutex rw;        // guards non‑list fields

        // ――― predicates ―――
        [[nodiscard]] bool isHead()      const noexcept { return seg_role == SegmentRole::HEAD; }
        [[nodiscard]] bool isTail()      const noexcept { return seg_role == SegmentRole::TAIL; }
        [[nodiscard]] bool isSegment()   const noexcept { return seg_role == SegmentRole::SEGMENT; }
        [[nodiscard]] bool isPrimary()   const noexcept { return align_role == AlignRole::PRIMARY; }
        [[nodiscard]] bool isSecondary() const noexcept { return align_role == AlignRole::SECONDARY; }

        // ――― factories ―――
        static SegPtr create(uint_t start, uint_t len, Strand sd,
            Cigar_t cg, AlignRole rl, SegmentRole sl,
            const BlockPtr& bp);

        static SegPtr createFromRegion(Region& region, Strand sd,
            Cigar_t cg, AlignRole rl, SegmentRole sl,
            const BlockPtr& bp);

        static SegPtr createHead();
        static SegPtr createTail();

        // ――― utilities ―――
        static void  linkChain(const std::vector<SegPtr>& segments);
        
        // ――― deletion utilities ―――
        static void unlinkSegment(SegPtr segment);
        static void deleteSegment(SegPtr segment);
        static void deleteBatch(const std::vector<SegPtr>& segments);
    };

    // ────────────────────────────────────────────────
    // Block – a bi‑species alignment block
    // ────────────────────────────────────────────────
    class Block : public std::enable_shared_from_this<Block> {
    public:
        ChrName   ref_chr;      // guard: rw
        ChrHeadMap anchors;     // guard: rw (head sentinel of every chr)
        mutable std::shared_mutex rw;

        static BlockPtr create(std::size_t hint = 1);
        static BlockPtr createEmpty(const ChrName& chr, std::size_t hint = 1);

        // Convenience helper – create both ref&qry segments, register anchors
        static std::pair<SegPtr, SegPtr> createSegmentPair(const Match& match,
            const SpeciesName& ref_name,
            const SpeciesName& qry_name,
            const ChrName& ref_chr,
            const ChrName& qry_chr,
            const BlockPtr& blk);

        static std::pair<SegPtr, SegPtr> createSegmentPair(const Anchor& anchor,
            const SpeciesName& ref_name,
            const SpeciesName& qry_name,
            const ChrName& ref_chr,
            const ChrName& qry_chr,
            const BlockPtr& blk);
            
        // ――― deletion methods ―――
        void removeAllSegments();
        void removeSegmentsBySpecies(const SpeciesName& species);
        bool removeSegment(const SpeciesName& species, const ChrName& chr);

        Block() = default;
        ~Block() = default;
    };

    // ────────────────────────────────────────────────
    // GenomeEnd  – head/tail sentinels per chromosome
    // ────────────────────────────────────────────────
    class GenomeEnd {
    public:
        static constexpr uint_t kSampleStep = 100000;
        GenomeEnd();
        GenomeEnd(GenomeEnd&&)            noexcept = default;
        GenomeEnd& operator=(GenomeEnd&&) noexcept = default;
        GenomeEnd(const GenomeEnd&) = delete;
        GenomeEnd& operator=(const GenomeEnd&) = delete;

        SegPtr findSurrounding(uint_t range_start);

        void insertSegment(const SegPtr seg);

        void clearAllSegments();
        
        // ――― deletion methods ―――
        bool removeSegment(SegPtr segment);
        bool removeSegmentRange(uint_t start, uint_t end);
        void removeBatch(const std::vector<SegPtr>& segments);
        void invalidateSampling(uint_t start, uint_t end);

        void resortSegments();

        void alignInterval(const SpeciesName ref_name, const SpeciesName query_name, const ChrName query_chr_name, SegPtr cur_node, std::map<SpeciesName, SeqPro::SharedManagerVariant> managers, bool is_left_extend, int_t zdrop);

        void removeOverlap(bool if_ref);

        SegPtr head{ nullptr };
        SegPtr tail{ nullptr };

        mutable std::shared_mutex rw;


        std::unique_ptr<Segment> head_holder;
        std::unique_ptr<Segment> tail_holder;

        std::vector<SegPtr> sample_vec;

        static bool spliceRange(SegPtr prev, SegPtr next,
            SegPtr first, SegPtr last);

        void ensureSampleSize(uint_t pos);          // 自动扩容


        
    
    public:
        void updateSampling(const std::vector<SegPtr>& segs); // 插入后修补

        void setToSampling(SegPtr cur);

    };

    // ────────────────────────────────────────────────
    // RaMeshGenomeGraph / RaMeshMultiGenomeGraph
    // ────────────────────────────────────────────────
    class RaMeshGenomeGraph {
    public:
        RaMeshGenomeGraph() = default;
        explicit RaMeshGenomeGraph(const SpeciesName& sp);
        RaMeshGenomeGraph(const SpeciesName& sp, const std::vector<ChrName>& chrs);

        size_t debugPrint(bool show_detail) const;

        SpeciesName                             species_name;
        std::unordered_map<ChrName, GenomeEnd>  chr2end;   // guard: rw
        mutable std::shared_mutex               rw;        // multi‑reader / single‑writer
    };

    class RaMeshMultiGenomeGraph {
    public:
        // RaMeshMultiGenomeGraph() = default;
        explicit RaMeshMultiGenomeGraph(std::map<SpeciesName, SeqPro::ManagerVariant>& seqpro_map);

        explicit RaMeshMultiGenomeGraph(std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_map);

   //     void insertClusterIntoGraph(SpeciesName ref_name, SpeciesName qry_name,
   //         const MatchCluster& cluster);

   //     void insertAnchorVecIntoGraph(SpeciesName ref_name, SpeciesName qry_name,
			//const AnchorVec& anchor_vec);

        void insertAnchorIntoGraph(SeqPro::ManagerVariant& ref_mgr, SeqPro::ManagerVariant& qry_mgr, SpeciesName ref_name, SpeciesName qry_name, const Anchor& anchor, bool isMultiple=false);

        void markAllExtended() {
            std::unique_lock lock(rw); // 锁保护整个 species_graphs
            for (auto& [species, genome_graph] : species_graphs) {
                std::shared_lock g_lock(genome_graph.rw);
                for (auto& [chr, end] : genome_graph.chr2end) {
                    std::shared_lock e_lock(end.rw);

                    SegPtr cur = end.head;
                    while (cur) {
                        if (!cur->isHead() && !cur->isTail()) {
                            std::unique_lock s_lock(cur->rw);
                            cur->left_extend = true;
                            cur->right_extend = true;
                        }
                        cur = cur->primary_path.next.load(std::memory_order_acquire);
                    }
                }
            }
        };


		void extendRefNodes(const SpeciesName& ref_name, std::map<SpeciesName, SeqPro::SharedManagerVariant> managers, int_t zdrop);
        void debugPrint(bool show_detail) const;
        
        // 图正确性验证函数
        bool verifyGraphCorrectness(bool verbose = false, bool show_detailed_segments = false) const;
        bool verifyGraphCorrectness(const SpeciesName& reference_species,
                                    bool verbose = false,
                                    bool show_detailed_segments = false,
                                    bool require_reference_overlap = true,
                                    bool forbid_non_reference_overlap = true,
                                    bool allow_reference_overlap = true) const;

        // ――― Enhanced verification system ―――
        enum class VerificationType : uint32_t {
            POINTER_VALIDITY = 1 << 0,        // 指针有效性
            LINKED_LIST_INTEGRITY = 1 << 1,   // 链表完整性
            COORDINATE_OVERLAP = 1 << 2,      // 坐标重叠检测
            COORDINATE_ORDERING = 1 << 3,     // 坐标排序检测
            BLOCK_CONSISTENCY = 1 << 4,       // Block一致性
            MEMORY_INTEGRITY = 1 << 5,        // 内存完整性
            THREAD_SAFETY = 1 << 6,           // 线程安全
            PERFORMANCE_ISSUES = 1 << 7,      // 性能问题
            BLOCK_REFERENCE_CHR = 1 << 8      // Block 缺失参考染色体或未注册参考锚点
        };

        enum class ErrorSeverity {
            INFO,     // 信息性问题
            WARNING,  // 警告
            ERROR,    // 错误
            CRITICAL  // 严重错误
        };

        struct VerificationError {
            VerificationType type;
            ErrorSeverity severity;
            std::string species;
            std::string chromosome;
            size_t segment_index;
            uint_t position;
            std::string message;
            std::string details;

            VerificationError(VerificationType t, ErrorSeverity s, const std::string& sp,
                            const std::string& chr, size_t idx, uint_t pos,
                            const std::string& msg, const std::string& det)
                : type(t), severity(s), species(sp), chromosome(chr),
                  segment_index(idx), position(pos), message(msg), details(det) {}
        };

        struct VerificationOptions {
            uint32_t enabled_checks = 0;          // 默认按需启用检查
            bool verbose = false;
            size_t max_errors_per_type = 100000;   // 完整统计所有错误
            size_t max_total_errors = 500000;      // 完整统计所有错误
            size_t max_verbose_errors_per_type = 5; // 每种类型只显示前5条详细信息
            bool stop_on_critical = false;
            bool include_performance_checks = false;
            bool show_detailed_segments = false;   // 是否显示详细的段信息调试日志

            struct ReferenceOverlapPolicy {
                bool enabled = false;
                SpeciesName reference_species;
                bool allow_reference_overlap = true;
                bool require_reference_overlap = false;
                bool forbid_non_reference_overlap = true;

                void reset() {
                    enabled = false;
                    reference_species.clear();
                    allow_reference_overlap = true;
                    require_reference_overlap = false;
                    forbid_non_reference_overlap = true;
                }
            };

            ReferenceOverlapPolicy reference_policy;

            VerificationOptions() {
                // 默认启用关键的结构与坐标检查
                enable(VerificationType::POINTER_VALIDITY);
                enable(VerificationType::LINKED_LIST_INTEGRITY);
                enable(VerificationType::COORDINATE_OVERLAP);
                enable(VerificationType::COORDINATE_ORDERING);
            }

            void enable(VerificationType type) {
                enabled_checks |= static_cast<uint32_t>(type);
            }

            void disable(VerificationType type) {
                enabled_checks &= ~static_cast<uint32_t>(type);
            }

            bool isEnabled(VerificationType type) const {
                return (enabled_checks & static_cast<uint32_t>(type)) != 0;
            }

            void enableReferenceOverlapPolicy(const SpeciesName& ref,
                                              bool allow_reference_overlap = true,
                                              bool require_reference_overlap = false,
                                              bool forbid_non_reference_overlap = true) {
                reference_policy.enabled = true;
                reference_policy.reference_species = ref;
                reference_policy.allow_reference_overlap = allow_reference_overlap;
                reference_policy.require_reference_overlap = require_reference_overlap;
                reference_policy.forbid_non_reference_overlap = forbid_non_reference_overlap;
            }

            void disableReferenceOverlapPolicy() {
                reference_policy.reset();
            }
        };

        struct VerificationResult {
            std::vector<VerificationError> errors;
            bool is_valid = true;
            std::chrono::microseconds verification_time{0};

            // 优化的错误计数器：避免每次都遍历整个错误向量
            mutable std::unordered_map<VerificationType, size_t> error_counts;
            mutable std::unordered_map<std::string,
                   std::unordered_map<VerificationType, size_t>> verbose_counts_by_species;
            mutable std::unordered_map<VerificationType, std::array<size_t, 4>>
                severity_verbose_counts;

            size_t getErrorCount(VerificationType type) const {
                return std::count_if(errors.begin(), errors.end(),
                    [type](const VerificationError& error) {
                        return error.type == type;
                    });
            }

            // 快速错误计数方法
            size_t getErrorCountFast(VerificationType type) const {
                auto it = error_counts.find(type);
                return it != error_counts.end() ? it->second : 0;
            }

            // 增加错误计数
            void incrementErrorCount(VerificationType type) {
                error_counts[type]++;
            }

            size_t getSpeciesVerboseCount(const std::string& species,
                                           VerificationType type) const {
                auto sit = verbose_counts_by_species.find(species);
                if (sit == verbose_counts_by_species.end()) return 0;
                auto tit = sit->second.find(type);
                return tit != sit->second.end() ? tit->second : 0;
            }

            void incrementSpeciesVerboseCount(const std::string& species,
                                              VerificationType type) {
                verbose_counts_by_species[species][type]++;
            }

            size_t getSeverityVerboseCount(VerificationType type,
                                           ErrorSeverity severity) const {
                auto it = severity_verbose_counts.find(type);
                if (it == severity_verbose_counts.end()) {
                    return 0;
                }
                auto idx = static_cast<size_t>(severity);
                if (idx >= it->second.size()) {
                    return 0;
                }
                return it->second[idx];
            }

            void incrementSeverityVerboseCount(VerificationType type,
                                               ErrorSeverity severity) {
                auto idx = static_cast<size_t>(severity);
                severity_verbose_counts[type][idx]++;
            }
        };

        VerificationResult verifyGraphCorrectness(const VerificationOptions& options) const;

        void safeLink(SegPtr prev, SegPtr next);

        void mergeMultipleGraphs(const SpeciesName& ref_name, uint_t thread_num);
        size_t mergeExactContiguousBlocks(
            const SpeciesName& reference_species,
            uint_t maximum_reference_span = 1000000,
            uint_t maximum_query_gap = 0);
        size_t realignSingleMissingSpeciesWindows(
            const SpeciesName& reference_species,
            const std::map<SpeciesName, SeqPro::SharedManagerVariant>&
                seqpro_managers,
            const std::string& msa_executable,
            uint_t maximum_span = 10000,
            uint_t parallel_threads = 1);
        void inspectExactContiguousBlockBoundaries(
            const SpeciesName& reference_species,
            const std::string& stage,
            uint_t maximum_reference_span = 1000000) const;


        std::unordered_map<SpeciesName, RaMeshGenomeGraph> species_graphs; // guard: rw
        std::vector<WeakBlock>                             blocks;         // guard: rw
        mutable std::shared_mutex                          rw;             // multi‑reader / single‑writer

        void exportToMaf(const FilePath& maf_path, const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers, bool only_primary, bool is_pairwise) const;

        void exportToMafWithoutReverse(const FilePath& maf_path, const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seq_mgrs, bool only_primary, bool pairwise_mode) const;

        void exportToMultipleMaf(const std::vector<std::pair<SpeciesName, FilePath>>& outs, const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seq_mgrs, bool only_primary, bool pairwise_mode) const;

        void exportToHal(const FilePath& hal_path,
                        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
                        const std::string& newick_tree = "",
                        bool only_primary = true,
                        const std::string& root_name = "root",
                        const SoftMask::PathMap& softmask_paths = {}) const;

        // 重载：直接复用已解析的 NewickParser，避免重复读取导致子树选择失效
        void exportToHal(const FilePath& hal_path,
                        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
                        const NewickParser& parser,
                        bool only_primary = true,
                        const std::string& root_name = "root",
                        const SoftMask::PathMap& softmask_paths = {}) const;

        
        // ――― high-performance deletion methods ―――
        bool removeBlock(const BlockPtr& block);
        bool removeBlockById(size_t block_index);
        void removeBlocksBatch(const std::vector<BlockPtr>& blocks);
        size_t removeExpiredBlocks();
        
        void removeSpecies(const SpeciesName& species);
        void removeChromosome(const SpeciesName& species, const ChrName& chr);
        void clearAllGraphs();
        
        // ――― garbage collection and maintenance ―――
        size_t compactBlockPool();
        void optimizeGraphStructure();
        
        // ――― deletion statistics and diagnostics ―――
        struct DeletionStats {
            size_t segments_removed = 0;
            size_t blocks_removed = 0;
            size_t expired_blocks_cleaned = 0;
            size_t sampling_updates = 0;
            std::chrono::microseconds total_time{0};
        };
        
        DeletionStats performMaintenance(bool full_gc = false);

    private:
        // ――― Enhanced verification helper methods ―――
        void addVerificationError(VerificationResult& result, const VerificationOptions& options,
                                VerificationType type, ErrorSeverity severity,
                                const std::string& species, const std::string& chr,
                                size_t segment_index, uint_t position,
                                const std::string& message, const std::string& details) const;

        bool shouldStopVerification(const VerificationResult& result, const VerificationOptions& options) const;
        void logVerificationSummary(const VerificationResult& result, const VerificationOptions& options) const;

        void verifyPointerValidity(VerificationResult& result, const VerificationOptions& options) const;
        void verifyLinkedListIntegrity(VerificationResult& result, const VerificationOptions& options) const;
        void verifyCoordinateOverlap(VerificationResult& result, const VerificationOptions& options) const;
        void verifyCoordinateOrdering(VerificationResult& result, const VerificationOptions& options) const;
        void verifyBlockConsistency(VerificationResult& result, const VerificationOptions& options) const;
        void verifyBlockReferenceChromosome(VerificationResult& result, const VerificationOptions& options) const;
        void verifyMemoryIntegrity(VerificationResult& result, const VerificationOptions& options) const;
        void verifyThreadSafety(VerificationResult& result, const VerificationOptions& options) const;
        void verifyPerformanceIssues(VerificationResult& result, const VerificationOptions& options) const;

        // 优化的统一遍历函数
        void verifyWithUnifiedTraversal(VerificationResult& result, const VerificationOptions& options) const;
    };

    void reportUnalignedRegions(const GenomeEnd& end,
        const SeqPro::SharedManagerVariant& mgr,
		const ChrName& chr_name);

} // namespace RaMesh

#endif /* RAMESH_H */
