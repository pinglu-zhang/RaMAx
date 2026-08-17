// RaMeshExport.cpp
// ------------------------------------------------------------
// 统一导出实现：emitMafBlock + 三个外层接口
// ------------------------------------------------------------
#include "ramesh.h"
#include "hal_converter.h"
#include "data_process.h" // for NewickParser definition
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <unordered_map>
#include <exception>
#include <map>
#include <sstream>
#include <vector>
#include <spdlog/spdlog.h>
#include <chrono>
#include <functional>
#include <thread>
#include <variant>
#include "halAlignmentInstance.h"
#include "halGenome.h"

// ============================================================
// emitMafBlock —— 所有导出函数共享的“写一个 MAF 块”实现
// ============================================================
struct MafExportStats {
    size_t invalid_blocks{ 0 };
    std::string first_error;

    void recordInvalid(const std::string& error) {
        ++invalid_blocks;
        if (first_error.empty()) first_error = error;
    }
};

static void logMafExportSummary(const FilePath& maf_path,
    const MafExportStats& stats)
{
    if (stats.invalid_blocks == 0) return;
    spdlog::warn(
        "MAF export skipped {} invalid block(s): output={}, first_error={}",
        stats.invalid_blocks, maf_path.string(), stats.first_error);
}

static std::filesystem::path mafTemporaryPath(
    const std::filesystem::path& output_path) {
    const auto stamp = std::chrono::steady_clock::now()
                           .time_since_epoch().count();
    const auto thread_hash =
        std::hash<std::thread::id>{}(std::this_thread::get_id());
    for (std::uint64_t attempt = 0; ; ++attempt) {
        auto candidate = output_path;
        candidate += ".tmp." + std::to_string(stamp) + "." +
                     std::to_string(thread_hash) + "." +
                     std::to_string(attempt);
        if (!std::filesystem::exists(candidate)) return candidate;
    }
}

static bool emitMafBlock(std::ostream& os,
    const RaMesh::BlockPtr& blk,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seq_mgrs,
    bool only_primary,
    bool pairwise_mode,
    bool allow_reverse,
    MafExportStats& stats,
    const SpeciesName* first_sp = nullptr,
    bool ensure_forward = true)
{
    if (!blk) return false;

    // ---------- 1. 收集 segment ----------
    struct Rec { SpeciesName sp; ChrName chr; RaMesh::SegPtr seg; };
    std::vector<Rec> recs;
    recs.reserve(blk->anchors.size());
    bool have_reverse = false;
    for (auto& [sp_chr, seg] : blk->anchors) {
        if (only_primary && !seg->isPrimary()) continue;
        if (!allow_reverse && seg->strand == Strand::REVERSE) have_reverse = true;
        recs.push_back({ sp_chr.first, sp_chr.second, seg });
    }
    if (recs.size() < 2 || (!allow_reverse && have_reverse)) return false;

    const auto ref_rec_it = std::find_if(recs.begin(), recs.end(),
        [&](const auto& r) {
            return r.sp == blk->ref_species && r.chr == blk->ref_chr;
        });
    if (ref_rec_it == recs.end()) {
        std::ostringstream oss;
        oss << "Block reference segment missing: ref_species="
            << blk->ref_species << ", ref_chr=" << blk->ref_chr
            << ", records=" << recs.size();
        stats.recordInvalid(oss.str());
        return false;
    }

    // ---------- 2. 决定首行 ----------
    if (first_sp) {
        auto it_first = std::find_if(recs.begin(), recs.end(),
            [&](auto& r) { return r.sp == *first_sp; });
        if (it_first == recs.end()) {
            auto it_ref = std::find_if(recs.begin(), recs.end(),
                [&](auto& r) {
                    return r.sp == blk->ref_species && r.chr == blk->ref_chr;
                });
            if (it_ref != recs.end()) std::swap(*recs.begin(), *it_ref);
        }
        else {
            if (it_first != recs.begin()) std::swap(*recs.begin(), *it_first);
        }
        
    }
    else {
        auto it_ref = std::find_if(recs.begin(), recs.end(),
            [&](auto& r) {
                return r.sp == blk->ref_species && r.chr == blk->ref_chr;
            });
        if (it_ref != recs.end()) std::swap(*recs.begin(), *it_ref);
    }

    // ---------- 3. lambda: fetchSeq / fetchLen ----------
    auto fetchSeq = [](const SeqPro::ManagerVariant& mv,
        const ChrName& chr, Coord_t b, Coord_t l) {
            return std::visit([&](auto& p) {
                using T = std::decay_t<decltype(p)>;
                if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::SequenceManager>>)
                    return p->getSubSequence(chr, b, l);
                else
                    return p->getOriginalManager().getSubSequence(chr, b, l);
                }, mv);
        };
    auto fetchLen = [](const SeqPro::ManagerVariant& mv,
        const ChrName& chr) {
            return std::visit([&](auto& p) {
                using T = std::decay_t<decltype(p)>;
                if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::SequenceManager>>)
                    return p->getSequenceLength(chr);
                else
                    return p->getOriginalManager().getSequenceLength(chr);
                }, mv);
        };

    // ---------- 4. 准备原始序列 ----------
    std::unordered_map<ChrName, std::string> seqs;
    std::unordered_map<ChrName, Cigar_t>     cigars;
    for (auto& r : recs) {
        auto mit = seq_mgrs.find(r.sp);
        if (mit == seq_mgrs.end()) {
            stats.recordInvalid(
                "Sequence manager missing for species=" + r.sp);
            return false;
        }
        std::string raw = fetchSeq(*mit->second, r.chr, r.seg->start, r.seg->length);

        if (raw.size() != r.seg->length) {
            std::ostringstream oss;
            oss << "Segment sequence length mismatch: ref_species="
                << blk->ref_species << ", ref_chr=" << blk->ref_chr
                << ", key=" << r.sp << '.' << r.chr
                << ", segment_length=" << r.seg->length
                << ", fetched_length=" << raw.size();
            stats.recordInvalid(oss.str());
            return false;
        }

        if (r.seg->strand == Strand::REVERSE) reverseComplement(raw);
        ChrName key = pairwise_mode ? r.chr : r.sp + "." + r.chr;
        if (!seqs.emplace(key, std::move(raw)).second ||
            !cigars.emplace(key, r.seg->cigar).second) {
            stats.recordInvalid(
                "Duplicate MAF sequence key in block: key=" + key);
            return false;
        }

    }

    // ---------- 5. 归并对齐 ----------
    const ChrName ref_key = pairwise_mode
        ? blk->ref_chr
        : blk->ref_species + "." + blk->ref_chr;
    const auto ref_seq_it = seqs.find(ref_key);
    if (ref_seq_it == seqs.end()) {
        stats.recordInvalid(
            "Reference sequence missing after extraction: key=" + ref_key);
        return false;
    }

    for (const auto& [key, cigar] : cigars) {
        if (key == ref_key) continue;
        const auto seq_it = seqs.find(key);
        if (seq_it == seqs.end()) {
            stats.recordInvalid(
                "CIGAR sequence missing after extraction: key=" + key);
            return false;
        }

        const AlignCount count = countAlignedBases(cigar);
        if (count.ref_bases != ref_seq_it->second.size() ||
            count.query_bases != seq_it->second.size()) {
            std::ostringstream oss;
            oss << "CIGAR consumption mismatch: ref=" << ref_key
                << ", key=" << key
                << ", cigar_ref=" << count.ref_bases
                << ", ref_size=" << ref_seq_it->second.size()
                << ", cigar_query=" << count.query_bases
                << ", query_size=" << seq_it->second.size();
            stats.recordInvalid(oss.str());
            return false;
        }
    }
    try {
        mergeAlignmentByRef(ref_key, seqs, cigars);
    }
    catch (const std::exception& e) {
        std::ostringstream oss;
        oss << "mergeAlignmentByRef failed: ref=" << ref_key
            << ", records=" << recs.size()
            << ", pairwise_mode=" << (pairwise_mode ? "true" : "false")
            << ", error=" << e.what();
        stats.recordInvalid(oss.str());
        return false;
    }

    // ---------- 6. 判断是否整体翻转 ----------
    bool need_flip = false;
    if (first_sp && ensure_forward)
        need_flip = (recs.front().seg->strand == Strand::REVERSE);
    std::unordered_map<ChrName, std::string> flipped;
    
    if (need_flip) {
        //auto rc = [](auto& s) { std::string t(s.rbegin(), s.rend()); for (char& c : t) { switch (c) { case 'A':c = 'T';break;case 'a':c = 't';break;case 'T':c = 'A';break;case 't':c = 'a';break;case 'C':c = 'G';break;case 'c':c = 'g';break;case 'G':c = 'C';break;case 'g':c = 'c';break;default:; } } return t; };
        for (auto& kv : seqs) reverseComplement(kv.second);
    }
    const auto* view = &seqs;

    // ---------- 7. 写块 ----------
    os << "a score=0\n";
    for (size_t i = 0;i < recs.size();++i) {
        auto& r = recs[i]; auto& mgr = *seq_mgrs.at(r.sp);
        uint64_t chr_len = fetchLen(mgr, r.chr);
        bool orig_rev = (r.seg->strand == Strand::REVERSE);
        bool final_rev = need_flip ? !orig_rev : orig_rev;
        uint64_t start = final_rev ? (chr_len - r.seg->start - r.seg->length) : r.seg->start;
        ChrName key = pairwise_mode ? r.chr : r.sp + "." + r.chr;
        os << "s " << std::left << std::setw(20) << key << std::right << std::setw(12) << start << std::setw(12) << r.seg->length
            << ' ' << (final_rev ? '-' : '+') << std::setw(12) << chr_len << ' ' << view->at(key) << "\n";
    }
    os << "\n";
    return true;
}

// ============================================================
// RaMeshMultiGenomeGraph 导出接口
// ============================================================
namespace RaMesh {

    void RaMeshMultiGenomeGraph::exportToMaf(
        const FilePath& maf_path,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seq_mgrs,
        bool only_primary,
        bool pairwise_mode) const
    {
        namespace fs = std::filesystem;
        const FilePath output_path = fs::absolute(maf_path);
        if (!output_path.parent_path().empty()) {
            fs::create_directories(output_path.parent_path());
        }
        const FilePath temporary_path = mafTemporaryPath(output_path);
        struct TemporaryGuard {
            FilePath path;
            bool keep{false};
            ~TemporaryGuard() {
                if (!keep) {
                    std::error_code error;
                    fs::remove(path, error);
                }
            }
        } guard{temporary_path};

        std::ofstream ofs(temporary_path, std::ios::binary | std::ios::trunc);
        if (!ofs) {
            throw std::runtime_error(
                "Cannot open MAF temporary output: " +
                temporary_path.string());
        }
        ofs << "##maf version=1 scoring=none\n";
        MafExportStats stats;
        for (auto& wblk : blocks) emitMafBlock(ofs, wblk.lock(), seq_mgrs, only_primary, pairwise_mode, true, stats, nullptr, true);
        ofs.close();
        if (!ofs) {
            throw std::runtime_error(
                "Failed to finalize MAF output: " + output_path.string());
        }

        std::error_code rename_error;
        fs::rename(temporary_path, output_path, rename_error);
        if (rename_error) {
            throw std::runtime_error(
                "Failed to atomically replace MAF output: " +
                rename_error.message());
        }
        guard.keep = true;
        logMafExportSummary(output_path, stats);
    }

    void RaMeshMultiGenomeGraph::exportToMafWithoutReverse(
        const FilePath& maf_path,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seq_mgrs,
        bool only_primary,
        bool pairwise_mode) const
    {
        namespace fs = std::filesystem;
        if (!maf_path.parent_path().empty()) fs::create_directories(maf_path.parent_path());
        std::ofstream ofs(maf_path, std::ios::binary | std::ios::trunc);
        if (!ofs) throw std::runtime_error("Cannot open: " + maf_path.string());
        ofs << "##maf version=1 scoring=none\n";
        MafExportStats stats;
        for (auto& wblk : blocks) emitMafBlock(ofs, wblk.lock(), seq_mgrs, only_primary, pairwise_mode, false, stats, nullptr, true);
        logMafExportSummary(maf_path, stats);
    }

    void RaMeshMultiGenomeGraph::exportToMultipleMaf(
        const std::vector<std::pair<SpeciesName, FilePath>>& outs,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seq_mgrs,
        bool only_primary,
        bool pairwise_mode) const
    {
        namespace fs = std::filesystem;
        struct Out { SpeciesName sp; FilePath path; std::ofstream ofs; MafExportStats stats; }; std::vector<Out> dests;
        dests.reserve(outs.size());
        for (auto& pr : outs) {
            if (!pr.second.parent_path().empty()) fs::create_directories(pr.second.parent_path());
            std::ofstream f(pr.second, std::ios::binary | std::ios::trunc);
            if (!f) throw std::runtime_error("Cannot open: " + pr.second.string());
            f << "##maf version=1 scoring=none\n";
            dests.push_back({ pr.first, pr.second, std::move(f), {} });
        }
        for (auto& wblk : blocks) {
            auto blk = wblk.lock(); if (!blk) continue;
            for (auto& out : dests) {
                emitMafBlock(out.ofs, blk, seq_mgrs, only_primary, pairwise_mode, true, out.stats, &out.sp, true);
            }
        }
        for (const auto& out : dests) logMafExportSummary(out.path, out.stats);
    }

    void RaMeshMultiGenomeGraph::exportToHal(const FilePath& hal_path,
                                            const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
                                            const std::string& newick_tree,
                                            bool only_primary,
                                            const std::string& root_name,
                                            const SoftMask::PathMap& softmask_paths) const
    {
        auto __now = []() { return std::chrono::steady_clock::now(); };
        auto __ms = [](auto a, auto b) { return std::chrono::duration_cast<std::chrono::milliseconds>(b - a).count(); };
        auto __us = [](auto a, auto b) { return std::chrono::duration_cast<std::chrono::microseconds>(b - a).count(); };
        auto __start_total = __now();
        spdlog::info("Starting HAL export to: {}", hal_path.string());
        spdlog::info("Export configuration: {} species, tree={}, only_primary={}",
                     seqpro_managers.size(), newick_tree.empty() ? "none" : "provided", only_primary);

        // ========================================
        // 基础验证和准备
        // ========================================
        if (seqpro_managers.empty()) {
            throw std::runtime_error("No sequence managers provided for HAL export");
        }

        if (softmask_paths.size() != seqpro_managers.size()) {
            throw std::runtime_error("HAL export requires one soft-mask index per leaf genome");
        }
        for (const auto& [species, unused_manager] : seqpro_managers) {
            (void)unused_manager;
            if (!softmask_paths.contains(species)) {
                throw std::runtime_error("Missing HAL soft-mask index for species: " + species);
            }
        }
        // Sidecars are deliberately opened only after alignment has completed
        // and HAL export has begun.
        const SoftMask::IndexMap softmask_indexes = SoftMask::loadIndexes(softmask_paths);

        if (blocks.empty()) {
            spdlog::warn("No alignment blocks found - creating HAL with genomes only");
        }

            // ========================================
            // 第一阶段：建立正确拓扑并写入叶基因组维度与DNA
            // - 基于真实Newick创建 root/内部祖先/叶 的拓扑
            // - 立刻为所有叶基因组 setDimensions + setString（来自 seqpro_managers）
            // ========================================
            auto __t1_start = __now();
            spdlog::info("Phase 1: Establishing basic HAL structure...");

            // 1.1 创建HAL文件和alignment对象
            auto __t11_start = __now();
            std::filesystem::path abs_hal_path = std::filesystem::absolute(hal_path);
            if (abs_hal_path.has_parent_path() && !abs_hal_path.parent_path().empty()) {
                std::filesystem::create_directories(abs_hal_path.parent_path());
            }

            // 重要：为避免生成结构不合法的 HAL 文件后仍残留在用户指定路径，采用“临时文件 + 校验通过后原子替换”策略。
            // - 成功：validate 通过后再 rename 为用户指定输出
            // - 失败：自动删除临时文件，保证用户不会误用半成品
            std::filesystem::path tmp_hal_path = abs_hal_path;
            tmp_hal_path += ".tmp";
            if (std::filesystem::exists(tmp_hal_path)) {
                std::filesystem::remove(tmp_hal_path);
            }
            struct HalTmpGuard {
                std::filesystem::path tmp;
                bool keep = false;
                ~HalTmpGuard() {
                    if (!keep) {
                        std::error_code ec;
                        std::filesystem::remove(tmp, ec);
                    }
                }
            } tmp_guard{tmp_hal_path};

            hal::AlignmentPtr alignment = hal::openHalAlignment(tmp_hal_path.string(), nullptr, hal::CREATE_ACCESS);
            if (!alignment) {
                throw std::runtime_error("Failed to create HAL alignment instance");
            }
            spdlog::info("  HAL temporary file created: {}", tmp_hal_path.string());
            auto __t11_end = __now();
            spdlog::debug("  Phase 1.1 (open/create HAL) took {} ms", __ms(__t11_start, __t11_end));

            // 1.2 解析系统发育树并识别祖先节点（单次解析并复用parser）
            std::vector<hal_converter::AncestorNode> ancestor_nodes;
            ::NewickParser parser; // 创建NewickParser对象用于获取分支长度

            spdlog::info("  Parsing phylogenetic tree...");
            auto __t12_start = __now();
            if (newick_tree.empty()) {
                throw std::runtime_error("Newick tree is required for HAL export but was empty");
            }
            auto __t12_end = __now();
            spdlog::debug("  Phase 1.2 (parse/extract tree) took {} ms", __ms(__t12_start, __t12_end));
            try {
                parser = ::NewickParser(newick_tree);
                hal_converter::ensureRootNode(parser, seqpro_managers, root_name);
                auto [is_valid, error_msg] = hal_converter::validateLeafNames(parser, seqpro_managers);
                if (!is_valid) {
                    spdlog::warn("Leaf name validation failed: {}", error_msg);
                    spdlog::warn("Proceeding with available species...");
                }
                ancestor_nodes = hal_converter::extractAncestorNodes(parser, seqpro_managers, root_name);
                spdlog::info("  Found {} ancestor nodes from tree", ancestor_nodes.size());
            } catch (const std::exception& e) {
                spdlog::error("Failed to parse Newick tree: {}", e.what());
                throw;
            }

            // 1.3 创建拓扑并为叶落盘维度与DNA
            if (!ancestor_nodes.empty()) {
                auto __t13_start = __now();
                hal_converter::createGenomesFromPhylogeny(alignment, ancestor_nodes, parser, root_name);
                spdlog::info("  Created genomes from phylogeny (topology only)");
                // 为所有叶基因组设置真实维度与DNA
                auto __t13a_start = __now();
                hal_converter::setupLeafGenomesWithRealDNA(
                    alignment, seqpro_managers, softmask_indexes);
                spdlog::info("  Set up leaf genomes with real DNA");
                auto __t13_end = __now();
                spdlog::debug("  Phase 1.3 (create topology) took {} ms; leaves DNA took {} ms",
                             __ms(__t13_start, __t13a_start), __ms(__t13a_start, __t13_end));
            }
            auto __t1_end = __now();
            spdlog::debug("Phase 1 completed successfully ({} ms)", __ms(__t1_start, __t1_end));

            // ========================================
            // 第二阶段：祖先序列重建并写入祖先基因组（维度+DNA）
            // - 2.1 生成祖先->参考叶计划
            // - 2.2 收集祖先重建片段（沿参考叶主链 + 缺口补齐）
            // - 2.3 按 voting 构建每条染色体的祖先序列，统一命名 ancestorName.chrN
            //       同步对祖先基因组 setDimensions + setString
            // ========================================
            auto __t2_start = __now();
            spdlog::info("Phase 2: Reconstructing ancestor sequences...");

            // 声明祖先序列变量，供后续阶段使用 (ancestor_name -> chr_name -> chr_sequence)
            std::map<std::string, std::map<std::string, std::string>> ancestor_sequences;
            std::map<std::string, hal_converter::AncestorReconstructionData> ancestor_reconstruction_data;

            if (!ancestor_nodes.empty() && !blocks.empty()) {
                // ========================================
                // 子阶段 2.1：确定重建顺序与参考叶
                // ========================================
                auto __t21_start = __now();
                spdlog::info("  Phase 2.1: Determining reconstruction order and reference leaves...");

                auto reconstruction_plan = planAncestorReconstruction(ancestor_nodes, parser);
                spdlog::info("  Planned reconstruction for {} ancestors", reconstruction_plan.size());

                for (const auto& [ancestor_name, ref_leaf] : reconstruction_plan) {
                    spdlog::debug("    Ancestor '{}' -> Reference leaf '{}'", ancestor_name, ref_leaf);
                }

                auto __t21_end = __now();
                spdlog::debug("  Phase 2.1 completed successfully ({} ms)", __ms(__t21_start, __t21_end));

                // ========================================
                // 子阶段 2.2：执行祖先节点重建
                // ========================================
                auto __t22_start = __now();
                spdlog::info("  Phase 2.2: Executing ancestor sequence reconstruction...");

                ancestor_reconstruction_data = hal_converter::reconstructAncestorSequences(
                    reconstruction_plan, ancestor_nodes, seqpro_managers, const_cast<RaMeshMultiGenomeGraph&>(*this));

                auto __t22_end = __now();
                spdlog::debug("  Phase 2.2 completed successfully ({} ms)", __ms(__t22_start, __t22_end));

                // ========================================
                // 子阶段 2.3：构建祖先序列（投票法）
                // ========================================
                auto __t23_start = __now();
                spdlog::info("  Phase 2.3: Building ancestor sequences using voting method...");

                ancestor_sequences = hal_converter::buildAllAncestorSequencesByVoting(
                    ancestor_reconstruction_data, ancestor_nodes, seqpro_managers,
                    reconstruction_plan, softmask_indexes, alignment);
                auto __t23_end = __now();
                spdlog::debug("  Phase 2.3 completed successfully ({} ms)", __ms(__t23_start, __t23_end));

            } else {
                spdlog::info("  Skipping ancestor reconstruction (no ancestors or blocks)");
            }
            auto __t2_end = __now();
            spdlog::debug("Phase 2 completed successfully ({} ms)", __ms(__t2_start, __t2_end));

            // ========================================
            // 第三阶段：映射阶段（从 RaMAx 块构建 HAL 段结构）
            // - 统计并更新每个基因组/染色体的段维度（top/bottom counts）
            // - 为各基因组创建并填充 Top/Bottom 段坐标
            // - 建立 parent-child 关系与 parse 信息
            // ========================================
	            auto __t3_start = __now();
	            spdlog::info("Phase 3: Building HAL segments from blocks (gap-aware mapping phase)...");

	            // 重要：在建立 parent/child segment 链接之前，必须先确定最终系统发育树。
	            // HAL 的 BottomSegment::childIndex 槽位是按“当前树中 child 的顺序”解释的；
	            // 若建链后再 replaceNewickTree/rename，会改变 child 顺序，导致槽位语义变化，从而生成结构不一致的 HAL。
	            try {
	                auto __t_tree_start = __now();
	                ::NewickParser final_parser = parser; // 以原解析器为基础
	                hal_converter::ensureRootNode(final_parser, seqpro_managers, root_name);
	                std::string final_newick = hal_converter::reconstructNewickFromParser(final_parser);
	                if (final_newick.empty()) {
	                    throw std::runtime_error("Final Newick reconstruction returned empty string");
	                }
	                alignment->replaceNewickTree(final_newick);
	                spdlog::info("Applied final Newick tree with root '{}': {}", root_name, final_newick);
	                auto __t_tree_end = __now();
	                spdlog::debug("Final tree application took {} ms", __ms(__t_tree_start, __t_tree_end));
	            } catch (const std::exception& e) {
	                spdlog::warn("Failed to apply final Newick tree before segment linking: {}", e.what());
	                // 回退：基于当前 HAL 拓扑重建一个简化 Newick（不含真实分支长度，统一 :0）
	                try {
	                    auto buildNewickFromAlignment = [&](auto&& self, const std::string& name) -> std::string {
	                        auto children = alignment->getChildNames(name);
	                        std::string res;
	                        if (!children.empty()) {
	                            res += "(";
	                            for (size_t i = 0; i < children.size(); ++i) {
	                                if (i) res += ",";
	                                res += self(self, children[i]);
	                            }
	                            res += ")";
	                        }
	                        res += name + ":0"; // 无法可靠获取分支长度，给定占位 0
	                        return res;
	                    };
	                    std::string root = alignment->getRootName();
	                    if (root.empty()) {
	                        root = root_name.empty() ? std::string("ancestor") : root_name;
	                    }
	                    std::string fallback_newick = buildNewickFromAlignment(buildNewickFromAlignment, root) + ";";
	                    alignment->replaceNewickTree(fallback_newick);
	                    spdlog::info("Applied fallback Newick tree: {}", fallback_newick);
	                } catch (const std::exception& e2) {
	                    throw std::runtime_error(std::string("Failed to apply Newick tree (final + fallback): ") + e2.what());
	                }
	            }
	            if (!blocks.empty() && !ancestor_nodes.empty()) {
	                hal_converter::analyzeBlocksAndBuildHalStructure(blocks, ancestor_nodes, alignment,
	                                                               ancestor_reconstruction_data, ancestor_sequences, seqpro_managers);
	            } else {
	                spdlog::info("  Skipping mapping phase (no blocks or ancestors)");
            }
            auto __t3_end = __now();
            spdlog::debug("Phase 3 completed successfully ({} ms)", __ms(__t3_start, __t3_end));


		            spdlog::info("Skipping internal HAL validation; run external halValidate for structure checks.");

		            // Write the final HAL path via atomic replacement.
		            alignment.reset(); // Close HDF5 handles before rename.
		            if (std::filesystem::exists(abs_hal_path)) {
		                std::filesystem::remove(abs_hal_path);
		            }
		            std::filesystem::rename(tmp_hal_path, abs_hal_path);
		            tmp_guard.keep = true;
		            spdlog::info("HAL file written: {}", abs_hal_path.string());

	            auto __end_total = __now();
	            spdlog::debug("HAL export finished. Total time: {} ms ({} us)", __ms(__start_total, __end_total), __us(__start_total, __end_total));

	    }

    // 重载实现：直接使用已解析的 NewickParser，避免重复解析
    void RaMeshMultiGenomeGraph::exportToHal(const FilePath& hal_path,
                                            const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
                                            const NewickParser& parser,
                                            bool only_primary,
                                            const std::string& root_name,
                                            const SoftMask::PathMap& softmask_paths) const
    {
        // 通过重建 Newick 字符串委托给字符串版本，避免重复读取原始文件，保留子树裁剪
        std::string rebuilt_newick = hal_converter::reconstructNewickFromParser(parser);
        if (rebuilt_newick.empty()) {
            throw std::runtime_error("exportToHal: reconstructed Newick from parser is empty");
        }
        exportToHal(hal_path, seqpro_managers, rebuilt_newick, only_primary,
                    root_name, softmask_paths);
    }
} // namespace RaMesh
