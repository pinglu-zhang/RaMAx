#include "ramesh.h"

#include "aligned_block_view.h"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <functional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <unordered_map>
#include <variant>

#include <spdlog/spdlog.h>

namespace {

using PafRow = RaMesh::Export::AlignedBlockRow;
using PreparedBlock = RaMesh::Export::AlignedBlockView;

std::string referenceName(const RaMesh::BlockPtr& block) {
    if (!block) return "<expired>";
    return RaMesh::Export::qualifiedName(
        block->ref_species, block->ref_chr);
}

void validatePafName(const std::string& name) {
    if (name.empty()) throw std::runtime_error("PAF sequence name is empty");
    if (std::any_of(name.begin(), name.end(), [](char value) {
            return std::isspace(static_cast<unsigned char>(value)) != 0;
        })) {
        throw std::runtime_error(
            "PAF sequence name contains whitespace: " + name);
    }
}

PreparedBlock prepareBlock(
    const RaMesh::BlockPtr& block,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    bool only_primary) {
    return RaMesh::Export::prepareAlignedBlock(
        block, managers, only_primary);
}

std::filesystem::path temporaryPath(
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

void recordInvalid(RaMesh::Paf::PafExportStats& stats,
                   const RaMesh::BlockPtr& block,
                   const std::string& reason) {
    ++stats.invalid_blocks_skipped;
    if (stats.first_invalid_reason.empty()) {
        stats.first_invalid_reference = referenceName(block);
        stats.first_invalid_reason = reason;
    }
}

}  // namespace

namespace RaMesh {

Paf::PafExportStats RaMeshMultiGenomeGraph::exportToPaf(
    const FilePath& paf_path,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    const Paf::PafExportOptions& options) const {
    namespace fs = std::filesystem;
    const auto started = std::chrono::steady_clock::now();
    if (managers.empty()) {
        throw std::runtime_error("No sequence managers provided for PAF export");
    }

    std::set<std::string> all_names;
    for (const auto& [species, manager] : managers) {
        for (const auto& chromosome : RaMesh::Export::fetchNames(*manager)) {
            const std::string name = RaMesh::Export::qualifiedName(
                species, chromosome);
            validatePafName(name);
            if (!all_names.insert(name).second) {
                throw std::runtime_error(
                    "Duplicate qualified PAF name: " + name);
            }
        }
    }

    const fs::path output_path = fs::absolute(paf_path);
    if (!output_path.parent_path().empty()) {
        fs::create_directories(output_path.parent_path());
    }
    const fs::path temporary_output = temporaryPath(output_path);
    struct TemporaryGuard {
        fs::path path;
        bool keep = false;
        ~TemporaryGuard() {
            if (!keep) {
                std::error_code error;
                fs::remove(path, error);
            }
        }
    } guard{temporary_output};

    std::ofstream output(
        temporary_output, std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error(
            "Cannot open PAF output: " + temporary_output.string());
    }

    Paf::PafExportStats stats;
    for (const auto& weak_block : blocks) {
        ++stats.blocks_seen;
        const auto block = weak_block.lock();
        if (!block) {
            ++stats.blocks_expired;
            continue;
        }

        try {
            PreparedBlock prepared = prepareBlock(
                block, managers, options.only_primary);
            if (prepared.rows.size() < 2) {
                ++stats.blocks_too_small;
                continue;
            }

            std::vector<std::string> aligned_rows;
            std::vector<std::string> names;
            aligned_rows.reserve(prepared.rows.size());
            names.reserve(prepared.rows.size());
            for (const auto& row : prepared.rows) {
                aligned_rows.push_back(row.aligned);
                names.push_back(row.name);
            }

            Paf::PairSelectionResult selection = Paf::selectPairs(
                aligned_rows, names, prepared.reference_index, options.mode);
            const std::size_t theoretical_pairs = selection.theoretical_pairs;
            const std::size_t eligible_pairs = selection.eligible_pairs;
            const std::size_t base_pairs = selection.base_pairs;
            const std::size_t supplemental_pairs =
                selection.supplemental_pairs;
            bool used_fallback = false;

            if (!selection.verified && options.mode == Paf::Mode::CONNECTED) {
                selection = Paf::selectPairs(
                    aligned_rows, names, prepared.reference_index,
                    Paf::Mode::ALL);
                used_fallback = true;
            }
            if (!selection.verified) {
                throw std::runtime_error(
                    "all-pairs fallback failed column-connectivity validation");
            }
            if (selection.pairs.empty()) {
                throw std::runtime_error(
                    "PAF Block has no pair with a shared non-gap column");
            }

            std::ostringstream block_output;
            std::size_t block_records = 0;
            std::size_t block_reverse_records = 0;
            std::uint64_t block_matching_bases = 0;
            std::uint64_t block_alignment_columns = 0;
            for (const auto& [left, right] : selection.pairs) {
                std::size_t target_index = left;
                std::size_t query_index = right;
                if (right == prepared.reference_index) {
                    target_index = right;
                    query_index = left;
                } else if (left != prepared.reference_index &&
                           prepared.rows[right].name <
                               prepared.rows[left].name) {
                    target_index = right;
                    query_index = left;
                }

                const PafRow& target = prepared.rows[target_index];
                const PafRow& query = prepared.rows[query_index];
                const Paf::PairProjection projection = Paf::projectPair(
                    target.aligned, query.aligned,
                    target.segment->strand == Strand::REVERSE);
                if (!projection.valid ||
                    projection.target_consumed != target.segment->length ||
                    projection.query_consumed != query.segment->length) {
                    throw std::runtime_error(
                        "Invalid PAF projection for " + query.name +
                        " vs " + target.name);
                }

                const char strand =
                    query.segment->strand == target.segment->strand
                        ? '+' : '-';
                block_output
                    << query.name << '\t'
                    << query.sequence_length << '\t'
                    << query.segment->start << '\t'
                    << query.segment->start + query.segment->length << '\t'
                    << strand << '\t'
                    << target.name << '\t'
                    << target.sequence_length << '\t'
                    << target.segment->start << '\t'
                    << target.segment->start + target.segment->length << '\t'
                    << projection.matching_bases << '\t'
                    << projection.block_length << '\t'
                    << 255 << '\t'
                    << "cg:Z:" << projection.cigar << '\t'
                    << "tp:A:P" << '\t'
                    << "NM:i:" << projection.edit_distance << '\n';
                if (!block_output) {
                    throw std::runtime_error(
                        "Failed while formatting PAF Block");
                }

                ++block_records;
                if (strand == '-') ++block_reverse_records;
                block_matching_bases += projection.matching_bases;
                block_alignment_columns += projection.block_length;
            }

            output << block_output.str();
            if (!output) {
                throw std::runtime_error("Failed while writing PAF output");
            }
            stats.theoretical_all_pairs += theoretical_pairs;
            stats.eligible_all_pairs += eligible_pairs;
            stats.base_pairs += base_pairs;
            stats.supplemental_pairs += supplemental_pairs;
            if (used_fallback) ++stats.fallback_blocks;
            stats.records_written += block_records;
            stats.reverse_records += block_reverse_records;
            stats.matching_bases += block_matching_bases;
            stats.alignment_columns += block_alignment_columns;
            ++stats.blocks_exported;
        } catch (const std::exception& error) {
            recordInvalid(stats, block, error.what());
        }
    }

    if (stats.records_written == 0) {
        throw std::runtime_error("PAF export produced no alignment records");
    }
    output.close();
    if (!output) throw std::runtime_error("Failed to finalize PAF output");

    std::error_code rename_error;
    fs::rename(temporary_output, output_path, rename_error);
    if (rename_error) {
        throw std::runtime_error(
            "Failed to atomically replace PAF output: " +
            rename_error.message());
    }
    guard.keep = true;

    stats.elapsed_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    const double compression_ratio = stats.eligible_all_pairs == 0
        ? 0.0
        : static_cast<double>(stats.records_written) /
              static_cast<double>(stats.eligible_all_pairs);
    spdlog::info(
        "PAF export complete: blocks(seen/exported/expired/too_small/invalid)="
        "{}/{}/{}/{}/{}, fallback_blocks={}, pairs(theoretical/eligible/base/"
        "supplemental/actual)={}/{}/{}/{}/{}, compression_ratio={:.6f}, "
        "reverse_records={}, matching_bases={}, elapsed_seconds={:.3f}",
        stats.blocks_seen, stats.blocks_exported, stats.blocks_expired,
        stats.blocks_too_small, stats.invalid_blocks_skipped,
        stats.fallback_blocks, stats.theoretical_all_pairs,
        stats.eligible_all_pairs, stats.base_pairs,
        stats.supplemental_pairs, stats.records_written,
        compression_ratio, stats.reverse_records, stats.matching_bases,
        stats.elapsed_seconds);
    if (stats.invalid_blocks_skipped != 0) {
        spdlog::warn(
            "PAF export skipped {} invalid Block(s): first_reference={}, "
            "first_error={}",
            stats.invalid_blocks_skipped, stats.first_invalid_reference,
            stats.first_invalid_reason);
    }
    return stats;
}

}  // namespace RaMesh
