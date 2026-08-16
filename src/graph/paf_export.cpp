#include "ramesh.h"

#include "align.h"

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

struct PafRow {
    SpeciesName species;
    ChrName chromosome;
    std::string name;
    RaMesh::SegPtr segment;
    std::uint64_t sequence_length = 0;
    std::string aligned;
};

struct PreparedBlock {
    std::vector<PafRow> rows;
    std::size_t reference_index = 0;
};

std::string qualifiedName(const SpeciesName& species,
                          const ChrName& chromosome) {
    return species + "." + chromosome;
}

std::string referenceName(const RaMesh::BlockPtr& block) {
    if (!block) return "<expired>";
    return qualifiedName(block->ref_species, block->ref_chr);
}

std::string fetchSequence(
    const SeqPro::ManagerVariant& manager,
    const ChrName& chromosome,
    Coord_t start,
    Coord_t length) {
    return std::visit([&](const auto& pointer) {
        using T = std::decay_t<decltype(pointer)>;
        if constexpr (std::is_same_v<
                          T, std::unique_ptr<SeqPro::SequenceManager>>) {
            return pointer->getSubSequence(chromosome, start, length);
        } else {
            return pointer->getOriginalManager().getSubSequence(
                chromosome, start, length);
        }
    }, manager);
}

std::uint64_t fetchLength(
    const SeqPro::ManagerVariant& manager,
    const ChrName& chromosome) {
    return std::visit([&](const auto& pointer) -> std::uint64_t {
        using T = std::decay_t<decltype(pointer)>;
        if constexpr (std::is_same_v<
                          T, std::unique_ptr<SeqPro::SequenceManager>>) {
            return pointer->getSequenceLength(chromosome);
        } else {
            return pointer->getOriginalManager().getSequenceLength(chromosome);
        }
    }, manager);
}

std::vector<std::string> fetchNames(
    const SeqPro::ManagerVariant& manager) {
    return std::visit([](const auto& pointer) {
        return pointer->getSequenceNames();
    }, manager);
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

std::size_t countUngapped(const std::string& sequence) {
    return static_cast<std::size_t>(std::count_if(
        sequence.begin(), sequence.end(),
        [](char base) { return base != '-'; }));
}

PreparedBlock prepareBlock(
    const RaMesh::BlockPtr& block,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    bool only_primary) {
    PreparedBlock prepared;
    if (!block) return prepared;

    for (const auto& [key, segment] : block->anchors) {
        if (!segment || (only_primary && !segment->isPrimary())) continue;
        const auto manager = managers.find(key.first);
        if (manager == managers.end()) {
            throw std::runtime_error(
                "Missing sequence manager for PAF species: " + key.first);
        }
        PafRow row;
        row.species = key.first;
        row.chromosome = key.second;
        row.name = qualifiedName(row.species, row.chromosome);
        row.segment = segment;
        row.sequence_length = fetchLength(*manager->second, row.chromosome);
        prepared.rows.push_back(std::move(row));
    }

    std::sort(prepared.rows.begin(), prepared.rows.end(),
              [](const PafRow& left, const PafRow& right) {
                  return left.name < right.name;
              });
    if (prepared.rows.size() < 2) return prepared;

    const auto reference = std::find_if(
        prepared.rows.begin(), prepared.rows.end(),
        [&](const PafRow& row) {
            return row.species == block->ref_species &&
                   row.chromosome == block->ref_chr;
        });
    if (reference == prepared.rows.end()) {
        throw std::runtime_error(
            "Block reference segment is missing from primary PAF rows");
    }
    prepared.reference_index = static_cast<std::size_t>(
        std::distance(prepared.rows.begin(), reference));

    std::unordered_map<ChrName, std::string> sequences;
    std::unordered_map<ChrName, Cigar_t> cigars;
    sequences.reserve(prepared.rows.size());
    cigars.reserve(prepared.rows.size());

    for (auto& row : prepared.rows) {
        if (row.segment->start > row.sequence_length ||
            row.segment->length > row.sequence_length - row.segment->start) {
            throw std::runtime_error(
                "PAF segment exceeds sequence bounds: " + row.name);
        }
        const auto manager = managers.find(row.species);
        std::string sequence = fetchSequence(
            *manager->second, row.chromosome,
            row.segment->start, row.segment->length);
        if (sequence.size() != row.segment->length) {
            throw std::runtime_error(
                "PAF extracted sequence length mismatch: " + row.name);
        }
        if (row.segment->strand == Strand::REVERSE) {
            reverseComplement(sequence);
        }
        if (!sequences.emplace(row.name, std::move(sequence)).second ||
            !cigars.emplace(row.name, row.segment->cigar).second) {
            throw std::runtime_error(
                "Duplicate PAF row name in Block: " + row.name);
        }
    }

    const std::string reference_name =
        prepared.rows[prepared.reference_index].name;
    const auto reference_sequence = sequences.find(reference_name);
    if (reference_sequence == sequences.end()) {
        throw std::runtime_error("PAF reference sequence is missing");
    }
    for (const auto& row : prepared.rows) {
        if (row.name == reference_name) continue;
        const AlignCount count = countAlignedBases(row.segment->cigar);
        if (count.ref_bases != reference_sequence->second.size() ||
            count.query_bases != row.segment->length) {
            std::ostringstream message;
            message << "PAF CIGAR consumption mismatch: key=" << row.name
                    << ", cigar_ref=" << count.ref_bases
                    << ", ref_size=" << reference_sequence->second.size()
                    << ", cigar_query=" << count.query_bases
                    << ", query_size=" << row.segment->length;
            throw std::runtime_error(message.str());
        }
    }

    mergeAlignmentByRef(reference_name, sequences, cigars);

    std::size_t width = 0;
    for (auto& row : prepared.rows) {
        const auto aligned = sequences.find(row.name);
        if (aligned == sequences.end()) {
            throw std::runtime_error(
                "PAF merged alignment lost row: " + row.name);
        }
        row.aligned = aligned->second;
        if (width == 0) width = row.aligned.size();
        if (row.aligned.size() != width) {
            throw std::runtime_error(
                "PAF merged alignment has unequal row widths");
        }
        if (countUngapped(row.aligned) != row.segment->length) {
            throw std::runtime_error(
                "PAF merged alignment does not preserve row length: " +
                row.name);
        }
    }
    return prepared;
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
        for (const auto& chromosome : fetchNames(*manager)) {
            const std::string name = qualifiedName(species, chromosome);
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
