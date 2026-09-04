#include "../src/align/reference_selection_internal.h"

#include <algorithm>
#include <cassert>
#include <string>
#include <vector>

struct Quality {
    std::string name;
    std::size_t sequence_count;
    bool reference_eligible;
};

int main() {
    using RaMAxReferenceSelection::enableMinimumSequenceFallback;
    std::vector<Quality> all_large = {
        {"Pongo_abelii", 22794, false},
        {"Gorilla_gorilla_gorilla", 26146, false},
        {"Pan_troglodytes", 20847, false},
        {"Homo_sapiens", 23274, false}};
    assert(enableMinimumSequenceFallback(all_large) == 2);
    assert(all_large[2].reference_eligible);
    assert(std::count_if(all_large.begin(), all_large.end(),
                        [](const auto& q) { return q.reference_eligible; }) == 1);

    std::vector<Quality> mixed = {{"large", 10001, false}, {"valid", 10000, true}};
    assert(enableMinimumSequenceFallback(mixed) == mixed.size());
    assert(!mixed[0].reference_eligible && mixed[1].reference_eligible);
    std::vector<Quality> ties = {{"first_by_quality", 10001, false},
                                {"second_by_quality", 10001, false}};
    assert(enableMinimumSequenceFallback(ties) == 0);
    assert(!ties[1].reference_eligible);
    std::vector<Quality> with_empty = {{"empty", 0, false}, {"large", 20000, false}};
    assert(enableMinimumSequenceFallback(with_empty) == 1);
    std::vector<Quality> empty_only = {{"empty", 0, false}};
    assert(enableMinimumSequenceFallback(empty_only) == empty_only.size());
    std::vector<Quality> no_inputs;
    assert(enableMinimumSequenceFallback(no_inputs) == no_inputs.size());
}
