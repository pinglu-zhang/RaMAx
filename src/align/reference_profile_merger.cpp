#include "reference_profile_merger.h"

#include "external_msa_runner.h"

#include <stdexcept>
#include <utility>

namespace RaMesh::Alignment {

uint_t mergeReferenceProfile(
    const ChrName& reference_name,
    std::unordered_map<ChrName, std::string>& sequences,
    const std::unordered_map<ChrName, Cigar_t>& cigars) {
    const auto reference = sequences.find(reference_name);
    if (reference == sequences.end()) {
        throw std::invalid_argument("mergeAlignmentByRef: ref not found");
    }

    std::string& reference_raw = reference->second;
    uint_t total_aligned_length = reference_raw.size();
    RefAlignInfo insertions;

    for (const auto& [key, cigar] : cigars) {
        if (key == reference_name) continue;
        const auto query = sequences.find(key);
        if (query == sequences.end()) {
            throw std::invalid_argument("mergeAlignmentByRef: seq missing");
        }

        std::string& query_raw = query->second;
        uint_t reference_position = 0;
        uint_t query_position = 0;
        for (const auto unit : cigar) {
            uint32_t length = 0;
            char operation = '\0';
            intToCigar(unit, operation, length);
            if (operation == 'D') {
                query_raw.insert(query_position, length, '-');
                reference_position += length;
                query_position += length;
            } else if (operation == 'I') {
                std::string insertion = query_raw.substr(
                    query_position, length);
                insertions[reference_position].seqs[key] =
                    std::move(insertion);
                query_raw.erase(query_position, length);
            } else {
                reference_position += length;
                query_position += length;
            }
        }
    }

    uint_t offset = 0;
    for (auto& [reference_position, insertion] : insertions) {
        insertion.alignSeqs();
        if (insertion.ref_name.empty()) continue;
        for (auto& [key, sequence] : sequences) {
            const auto row = insertion.seqs.find(key);
            if (row != insertion.seqs.end()) {
                sequence.insert(reference_position + offset, row->second);
            } else {
                sequence.insert(reference_position + offset,
                                insertion.total_length, '-');
            }
        }
        offset += insertion.total_length;
        total_aligned_length += insertion.total_length;
    }
    return total_aligned_length;
}

}  // namespace RaMesh::Alignment

bool alignSequencesWithExternalMsa(
    const std::string& executable,
    std::unordered_map<ChrName, std::string>& sequences) {
    return RaMesh::Alignment::ExternalMsaRunner::instance().align(
        executable, sequences);
}

void configureExternalInsertionMsa(const std::string& executable) {
    RaMesh::Alignment::ExternalMsaRunner::instance()
        .configureDefaultExecutable(executable);
}

void InsertInfo::alignSeqs() {
    if (aligned || seqs.empty()) return;
    if (seqs.size() == 1) {
        ref_name = seqs.begin()->first;
        total_length = seqs.begin()->second.size();
        aligned = true;
        return;
    }

    if (RaMesh::Alignment::ExternalMsaRunner::instance()
            .alignWithDefault(seqs)) {
        ref_name = seqs.begin()->first;
        total_length = seqs.begin()->second.size();
        aligned = true;
        return;
    }

    size_t maximum_length = 0;
    for (const auto& [key, sequence] : seqs) {
        if (sequence.size() > maximum_length) {
            maximum_length = sequence.size();
            ref_name = key;
        }
    }
    std::unordered_map<ChrName, Cigar_t> cigars;
    for (const auto& [key, sequence] : seqs) {
        if (key != ref_name) {
            cigars[key] = globalAlignKSW2(seqs[ref_name], sequence);
        }
    }
    total_length = RaMesh::Alignment::mergeReferenceProfile(
        ref_name, seqs, cigars);
    aligned = true;
}
