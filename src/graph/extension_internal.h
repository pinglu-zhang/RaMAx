#ifndef RAMAX_GRAPH_EXTENSION_INTERNAL_H
#define RAMAX_GRAPH_EXTENSION_INTERNAL_H

#include "ramesh.h"

#include <cstdint>
#include <string>

namespace RaMesh::detail {

constexpr int_t kMaximumExtensionGap = 10000;

struct ExtensionSequenceContext {
    const SeqPro::SequenceManager* manager{nullptr};
    SeqPro::SequenceId sequence_id{SeqPro::SequenceIndex::INVALID_ID};
    SeqPro::Length chromosome_length{0};
};

struct ExtensionSequenceScratch {
    std::string reference;
    std::string query;

    void reserve() {
        reference.reserve(static_cast<std::size_t>(kMaximumExtensionGap));
        query.reserve(static_cast<std::size_t>(kMaximumExtensionGap));
    }
};

enum class PreparedExtensionStatus : uint8_t {
    SKIPPED,
    NO_CANDIDATE,
    NONPOSITIVE_GAP,
    GAP_TOO_LONG,
    ALIGNMENT_REJECTED,
    ALIGNMENT_ACCEPTED,
};

struct PreparedExtensionResult {
    PreparedExtensionStatus status{PreparedExtensionStatus::SKIPPED};

    [[nodiscard]] bool alignmentCalled() const noexcept {
        return status == PreparedExtensionStatus::ALIGNMENT_REJECTED ||
               status == PreparedExtensionStatus::ALIGNMENT_ACCEPTED;
    }

    [[nodiscard]] bool accepted() const noexcept {
        return status == PreparedExtensionStatus::ALIGNMENT_ACCEPTED;
    }
};

const SeqPro::SequenceManager& originalSequenceManager(
    const SeqPro::SharedManagerVariant& manager);

PreparedExtensionResult alignIntervalPrepared(
    Segment* query_segment,
    Segment* reference_segment,
    Segment* reference_candidate,
    const ExtensionSequenceContext& query_context,
    const ExtensionSequenceContext& reference_context,
    bool is_left_extend,
    int_t zdrop,
    ExtensionSequenceScratch& scratch);

}  // namespace RaMesh::detail

#endif  // RAMAX_GRAPH_EXTENSION_INTERNAL_H
