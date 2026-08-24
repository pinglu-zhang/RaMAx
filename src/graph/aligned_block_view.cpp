#include "aligned_block_view.h"

#include "align.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <variant>

namespace RaMesh::Export {
namespace {

std::size_t countUngapped(const std::string& sequence) {
    return static_cast<std::size_t>(std::count_if(
        sequence.begin(), sequence.end(),
        [](char base) { return base != '-'; }));
}

}  // namespace

std::string qualifiedName(const SpeciesName& species,
                          const ChrName& chromosome) {
    return species + "." + chromosome;
}

std::string fetchSequence(const SeqPro::ManagerVariant& manager,
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

std::uint64_t fetchLength(const SeqPro::ManagerVariant& manager,
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

std::vector<std::string> fetchNames(const SeqPro::ManagerVariant& manager) {
    return std::visit([](const auto& pointer) {
        return pointer->getSequenceNames();
    }, manager);
}

AlignedBlockView prepareAlignedBlock(
    const BlockPtr& block,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    bool only_primary) {
    AlignedBlockView prepared;
    if (!block) return prepared;

    for (const auto& [key, segment] : block->anchors) {
        if (!segment || (only_primary && !segment->isPrimary())) continue;
        const auto manager = managers.find(key.first);
        if (manager == managers.end()) {
            throw std::runtime_error(
                "Missing sequence manager for Block species: " + key.first);
        }
        AlignedBlockRow row;
        row.species = key.first;
        row.chromosome = key.second;
        row.name = qualifiedName(row.species, row.chromosome);
        row.segment = segment;
        row.sequence_length = fetchLength(*manager->second, row.chromosome);
        prepared.rows.push_back(std::move(row));
    }

    std::sort(prepared.rows.begin(), prepared.rows.end(),
              [](const AlignedBlockRow& left, const AlignedBlockRow& right) {
                  return left.name < right.name;
              });
    if (prepared.rows.size() < 2) return prepared;

    const auto reference = std::find_if(
        prepared.rows.begin(), prepared.rows.end(),
        [&](const AlignedBlockRow& row) {
            return row.species == block->ref_species &&
                   row.chromosome == block->ref_chr;
        });
    if (reference == prepared.rows.end()) {
        throw std::runtime_error(
            "Block reference segment is missing from primary rows");
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
                "Block segment exceeds sequence bounds: " +
                row.name);
        }
        const auto manager = managers.find(row.species);
        std::string sequence = fetchSequence(
            *manager->second, row.chromosome,
            row.segment->start, row.segment->length);
        if (sequence.size() != row.segment->length) {
            throw std::runtime_error(
                "Block extracted sequence length mismatch: " +
                row.name);
        }
        if (row.segment->strand == Strand::REVERSE) {
            reverseComplement(sequence);
        }
        if (!sequences.emplace(row.name, std::move(sequence)).second ||
            !cigars.emplace(row.name, row.segment->cigar).second) {
            throw std::runtime_error(
                "Duplicate row name in Block: " + row.name);
        }
    }

    const std::string reference_name =
        prepared.rows[prepared.reference_index].name;
    const auto reference_sequence = sequences.find(reference_name);
    if (reference_sequence == sequences.end()) {
        throw std::runtime_error("Block reference sequence is missing");
    }
    for (const auto& row : prepared.rows) {
        if (row.name == reference_name) continue;
        const AlignCount count = countAlignedBases(row.segment->cigar);
        if (count.ref_bases != reference_sequence->second.size() ||
            count.query_bases != row.segment->length) {
            std::ostringstream message;
            message << "Block CIGAR consumption mismatch: key="
                    << row.name
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
                "Merged alignment lost row: " + row.name);
        }
        row.aligned = aligned->second;
        if (width == 0) width = row.aligned.size();
        if (row.aligned.size() != width) {
            throw std::runtime_error(
                "Merged alignment has unequal row widths");
        }
        if (countUngapped(row.aligned) != row.segment->length) {
            throw std::runtime_error(
                "Merged alignment does not preserve row length: " +
                row.name);
        }
    }
    return prepared;
}

}  // namespace RaMesh::Export
