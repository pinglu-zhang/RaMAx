#pragma once

#include <cstddef>

namespace RaMAxReferenceSelection {

// Qualities are already ordered by N50, total length, and name. Retaining
// the first minimum therefore provides a deterministic tie-break.
// Return size() when a regular reference exists or no nonempty input exists.
template <typename Qualities>
std::size_t enableMinimumSequenceFallback(Qualities& qualities) {
    for (const auto& quality : qualities) {
        if (quality.reference_eligible) {
            return qualities.size();
        }
    }
    std::size_t selected = qualities.size();
    for (std::size_t i = 0; i < qualities.size(); ++i) {
        if (qualities[i].sequence_count > 0 &&
            (selected == qualities.size() ||
             qualities[i].sequence_count < qualities[selected].sequence_count)) {
            selected = i;
        }
    }
    if (selected != qualities.size()) {
        qualities[selected].reference_eligible = true;
    }
    return selected;
}

}  // namespace RaMAxReferenceSelection
