#pragma once

#include "hal/export.h"
#include "halAlignmentInstance.h"

#include <cstdint>
#include <filesystem>
#include <optional>
#include <map>
#include <string>
#include <utility>
#include <string_view>
#include <vector>

namespace RaMesh::hal_export {

struct TopSegmentLine {
    OccurrenceId occurrence_id = 0;
    uint64_t start = 0;
    uint32_t length = 0;
    std::optional<uint64_t> parent_bottom_name;
    bool forward_to_parent = true;
};

struct BottomSegmentLine {
    OccurrenceId occurrence_id = 0;
    uint64_t name = 0;
    uint64_t start = 0;
    uint32_t length = 0;
};

struct SequenceEmission {
    std::string genome_name;
    std::string seq_name;
    size_t bottom_count = 0;
    std::vector<BottomSegmentLine> bottoms;
    std::vector<TopSegmentLine> tops;
    std::optional<std::string_view> dna;
    std::optional<std::pair<std::string, std::string>> leaf_source;
};

// Emission DNA is borrowed for appendSubtree only. The writer does not retain
// emission records. It keeps the root arrays open while each child's arrays
// are populated and flushed in turn, rather than retaining all child arrays.
class NativeHalWriter {
public:
    NativeHalWriter(
        const std::filesystem::path& path,
        const TreeMeta& tree,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
        const SoftMask::IndexMap& softmask_indexes);
    ~NativeHalWriter();

    NativeHalWriter(const NativeHalWriter&) = delete;
    NativeHalWriter& operator=(const NativeHalWriter&) = delete;

    void appendSubtree(int node_id, const std::vector<SequenceEmission>& emissions);
    void close();

private:
    hal::AlignmentPtr alignment_;
    const TreeMeta& tree_;
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers_;
    const SoftMask::IndexMap& softmask_indexes_;
};

}  // namespace RaMesh::hal_export
