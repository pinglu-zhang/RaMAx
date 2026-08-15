// RaMeshExport.cpp
// ------------------------------------------------------------
// 统一导出实现：emitMafBlock + 三个外层接口
// ------------------------------------------------------------
#include "ramesh.h"
#include "hal/export.h"
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <unordered_map>
#include <exception>
#include <map>
#include <vector>
#include <spdlog/spdlog.h>
#include <variant>

// ============================================================
// emitMafBlock —— 所有导出函数共享的“写一个 MAF 块”实现
// ============================================================
static bool isMafReferenceCompatibleCigar(
    const Cigar_t& cigar,
    size_t sequence_length) {
    size_t consumed = 0;
    for (const auto unit : cigar) {
        uint32_t length = 0;
        char operation = '\0';
        intToCigar(unit, operation, length);
        if (operation != 'M' && operation != '=' && operation != 'X') {
            return false;
        }
        consumed += length;
    }
    return consumed == sequence_length;
}

static bool emitMafBlock(std::ostream& os,
    const RaMesh::BlockPtr& blk,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seq_mgrs,
    bool only_primary,
    bool pairwise_mode,
    bool allow_reverse,
    const SpeciesName* first_sp = nullptr,
    bool ensure_forward = true)
{
    if (!blk) return false;

    // ---------- 1. 收集 segment ----------
    struct Rec {
        SpeciesName sp;
        ChrName chr;
        RaMesh::SegPtr seg;
        ChrName row_key;
    };
    std::vector<Rec> recs;
    recs.reserve(blk->anchors.size());
    bool have_reverse = false;
    for (auto& [sp_chr, seg] : blk->anchors) {
        if (only_primary && !seg->isPrimary()) continue;
        if (!allow_reverse && seg->strand == Strand::REVERSE) have_reverse = true;
        recs.push_back({ sp_chr.first, sp_chr.second, seg, {} });
    }
    if (recs.size() < 2 || (!allow_reverse && have_reverse)) return false;

    // ---------- 2. 选择对齐参考并决定首行 ----------
    // ref_chr does not identify a species.  When several genomes use the
    // same chromosome name, unordered anchor iteration must not choose a
    // deletion-only hybrid row as the reference backbone.
    const size_t declared_reference_occurrences = std::count_if(
        recs.begin(), recs.end(), [&](const Rec& record) {
            return record.sp == blk->ref_species &&
                record.chr == blk->ref_chr;
        });
    if (declared_reference_occurrences != 1) {
        throw std::runtime_error(
            "MAF export failed: block " + std::to_string(blk->block_id) +
            " must contain exactly one declared reference occurrence '" +
            blk->ref_species + "." + blk->ref_chr + "', found " +
            std::to_string(declared_reference_occurrences));
    }
    auto alignment_ref = std::find_if(
        recs.begin(), recs.end(), [&](const Rec& record) {
            return record.sp == blk->ref_species &&
                record.chr == blk->ref_chr &&
                isMafReferenceCompatibleCigar(
                    record.seg->cigar, record.seg->length);
        });
    if (alignment_ref == recs.end()) {
        alignment_ref = std::find_if(
            recs.begin(), recs.end(),
            [&](const Rec& record) {
                return record.sp == blk->ref_species &&
                    record.chr == blk->ref_chr;
            });
    }
    if (alignment_ref == recs.end()) {
        throw std::runtime_error(
            "MAF export failed: declared reference row disappeared");
    }
    const RaMesh::SegPtr alignment_ref_segment = alignment_ref->seg;

    if (first_sp) {
        auto it_first = std::find_if(recs.begin(), recs.end(),
            [&](auto& r) { return r.sp == *first_sp; });
        if (it_first == recs.end()) {
            auto it_ref = std::find_if(
                recs.begin(), recs.end(),
                [&](const Rec& record) {
                    return record.seg == alignment_ref_segment;
                });
            if (it_ref != recs.end()) std::swap(*recs.begin(), *it_ref);
        }
        else {
            if (it_first != recs.begin()) std::swap(*recs.begin(), *it_first);
        }
        
    }
    else {
        auto it_ref = std::find_if(
            recs.begin(), recs.end(),
            [&](const Rec& record) {
                return record.seg == alignment_ref_segment;
            });
        if (it_ref != recs.end()) std::swap(*recs.begin(), *it_ref);
    }

    std::unordered_map<ChrName, uint32_t> next_occurrence;
    for (auto& record : recs) {
        const ChrName source = pairwise_mode
            ? record.chr
            : record.sp + "." + record.chr;
        record.row_key =
            source + '\1' + std::to_string(next_occurrence[source]++);
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
        if (mit == seq_mgrs.end()) { seqs.clear(); break; }
        std::string raw = fetchSeq(*mit->second, r.chr, r.seg->start, r.seg->length);

        if (r.seg->strand == Strand::REVERSE) reverseComplement(raw);
        seqs.emplace(r.row_key, std::move(raw));
        cigars.emplace(r.row_key, r.seg->cigar);

    }
    if (seqs.empty()) return false;

    // ---------- 5. 归并对齐 ----------
    const auto reference_record = std::find_if(
        recs.begin(), recs.end(),
        [&](const Rec& record) {
            return record.seg == alignment_ref_segment;
        });
    if (reference_record == recs.end()) return false;
    const ChrName& ref_key = reference_record->row_key;
    try {
        mergeAlignmentByRef(ref_key, seqs, cigars);
    }
    catch (const std::exception& e) {
        spdlog::warn("mergeAlignmentByRef failed: {}", e.what());
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
        const ChrName source =
            pairwise_mode ? r.chr : r.sp + "." + r.chr;
        os << "s " << std::left << std::setw(20) << source
            << std::right << std::setw(12) << start
            << std::setw(12) << r.seg->length
            << ' ' << (final_rev ? '-' : '+') << std::setw(12) << chr_len
            << ' ' << view->at(r.row_key) << "\n";
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
        if (!maf_path.parent_path().empty()) fs::create_directories(maf_path.parent_path());
        std::ofstream ofs(maf_path, std::ios::binary | std::ios::trunc);
        if (!ofs) throw std::runtime_error("Cannot open: " + maf_path.string());
        ofs << "##maf version=1 scoring=none\n";
        for (auto& wblk : blocks) emitMafBlock(ofs, wblk.lock(), seq_mgrs, only_primary, pairwise_mode, true, nullptr, true);
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
        for (auto& wblk : blocks) emitMafBlock(ofs, wblk.lock(), seq_mgrs, only_primary, pairwise_mode, false, nullptr, true);
    }

    void RaMeshMultiGenomeGraph::exportToMultipleMaf(
        const std::vector<std::pair<SpeciesName, FilePath>>& outs,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seq_mgrs,
        bool only_primary,
        bool pairwise_mode) const
    {
        namespace fs = std::filesystem;
        struct Out { SpeciesName sp; std::ofstream ofs; }; std::vector<Out> dests;
        dests.reserve(outs.size());
        for (auto& pr : outs) {
            if (!pr.second.parent_path().empty()) fs::create_directories(pr.second.parent_path());
            std::ofstream f(pr.second, std::ios::binary | std::ios::trunc);
            if (!f) throw std::runtime_error("Cannot open: " + pr.second.string());
            f << "##maf version=1 scoring=none\n";
            dests.push_back({ pr.first,std::move(f) });
        }
        for (auto& wblk : blocks) {
            auto blk = wblk.lock(); if (!blk) continue;
            for (auto& out : dests) {
                emitMafBlock(out.ofs, blk, seq_mgrs, only_primary, pairwise_mode, true, &out.sp, true);
            }
        }
    }

    void RaMeshMultiGenomeGraph::exportToHal(
        const FilePath& hal_path,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
        const NewickParser& parser,
        const std::string& root_name,
        int parallel_threads,
        const SoftMask::PathMap& softmask_paths) const {
        if (softmask_paths.size() != seqpro_managers.size()) {
            throw std::runtime_error(
                "HAL export requires one soft-mask index per leaf genome");
        }
        for (const auto& [species, unused_manager] : seqpro_managers) {
            (void)unused_manager;
            if (!softmask_paths.contains(species)) {
                throw std::runtime_error(
                    "Missing HAL soft-mask index for species: " + species);
            }
        }
        const SoftMask::IndexMap softmask_indexes =
            SoftMask::loadIndexes(softmask_paths);
        hal_export::exportToHal(
            blocks,
            hal_path,
            seqpro_managers,
            parser,
            root_name,
            softmask_indexes,
            hal_export::ExportConfig{
                .parallel_threads = parallel_threads,
            });
    }

} // namespace RaMesh
