#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace RaMAxReferenceSelection {

// Explicit references take priority over automatic eligibility and ordering.
// Only the requested genome is promoted; other candidates retain their order.
template <typename Qualities>
void prioritizeExplicitReference(Qualities& qualities, const std::string& name) {
    auto selected = std::find_if(qualities.begin(), qualities.end(),
        [&name](const auto& quality) { return quality.name == name; });
    if (selected == qualities.end()) {
        throw std::runtime_error("[reference-selection] Explicit reference " + name +
                                 " was not found in the input genomes.");
    }
    if (selected->sequence_count == 0 || selected->total_length == 0) {
        throw std::runtime_error("[reference-selection] Explicit reference " + name +
                                 " has no usable sequence bases.");
    }
    selected->reference_eligible = true;
    std::rotate(qualities.begin(), selected, selected + 1);
}

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
        if (qualities[i].sequence_count > 0 && qualities[i].total_length > 0 &&
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
