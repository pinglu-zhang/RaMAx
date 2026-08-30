#ifndef RAMAX_ANCHOR_LINK_INTERNAL_H
#define RAMAX_ANCHOR_LINK_INTERNAL_H

#include "anchor.h"

namespace AnchorLinkDetail {

struct Statistics {
    uint64_t candidate_checks{0};
    uint64_t sequence_extractions{0};
    uint64_t direct_ksw_calls{0};
    uint64_t fallback_ksw_calls{0};
    uint64_t long_gap_rejections{0};
    uint64_t maximum_seen_gap{0};
    uint64_t estimated_ksw_cells{0};
};

struct ComponentRange {
    size_t begin{0};
    size_t end{0};
    uint64_t estimated_cost{0};
};

AnchorVec materializeClusterAnchors(
    MatchClusterVec& clusters,
    const SeqPro::ManagerVariant& ref_mgr,
    const SeqPro::ManagerVariant& qry_mgr);

std::vector<ComponentRange> splitAnchorComponents(
    const AnchorVec& anchors,
    Statistics* statistics = nullptr);

AnchorPtrVec linkAnchorRange(
    AnchorVec& anchors,
    size_t begin,
    size_t end,
    const SeqPro::ManagerVariant& ref_mgr,
    const SeqPro::ManagerVariant& qry_mgr,
    Statistics* statistics = nullptr);

}  // namespace AnchorLinkDetail

#endif
