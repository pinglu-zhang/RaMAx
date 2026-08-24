#include "paf_export.h"

#include <algorithm>
#include <cctype>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

namespace RaMesh::Paf {
namespace {

class DisjointSet {
public:
    explicit DisjointSet(std::size_t size) : parent_(size), rank_(size, 0) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }

    std::size_t find(std::size_t value) {
        if (parent_[value] != value) parent_[value] = find(parent_[value]);
        return parent_[value];
    }

    void unite(std::size_t left, std::size_t right) {
        left = find(left);
        right = find(right);
        if (left == right) return;
        if (rank_[left] < rank_[right]) std::swap(left, right);
        parent_[right] = left;
        if (rank_[left] == rank_[right]) ++rank_[left];
    }

private:
    std::vector<std::size_t> parent_;
    std::vector<unsigned char> rank_;
};

char normalizedBase(char base) {
    if (base == '-') return '-';
    return static_cast<char>(
        std::toupper(static_cast<unsigned char>(base)));
}

bool isConnectableBase(char base) {
    const char normalized = normalizedBase(base);
    return normalized != '-' && normalized != 'X';
}

SequencePair orderedPair(std::size_t left, std::size_t right) {
    return left < right ? SequencePair{left, right}
                        : SequencePair{right, left};
}

void validateRows(const std::vector<std::string>& rows,
                  const std::vector<std::string>* names = nullptr) {
    if (rows.empty()) throw std::invalid_argument("PAF alignment has no rows");
    const std::size_t width = rows.front().size();
    if (width == 0) throw std::invalid_argument("PAF alignment has zero columns");
    for (const auto& row : rows) {
        if (row.size() != width) {
            throw std::invalid_argument("PAF alignment rows have unequal widths");
        }
    }
    if (names && names->size() != rows.size()) {
        throw std::invalid_argument("PAF alignment names do not match rows");
    }
}

class PairOverlapMatrix {
public:
    explicit PairOverlapMatrix(std::size_t size)
        : size_(size), bits_((pairCount(size) + 63) / 64, 0) {}

    void mark(std::size_t left, std::size_t right) {
        if (left == right) return;
        const std::size_t index = pairIndex(left, right);
        const std::uint64_t mask = std::uint64_t{1} << (index & 63);
        auto& word = bits_[index >> 6];
        if ((word & mask) == 0) {
            word |= mask;
            ++count_;
        }
    }

    bool contains(std::size_t left, std::size_t right) const {
        if (left == right) return false;
        const std::size_t index = pairIndex(left, right);
        return (bits_[index >> 6] &
                (std::uint64_t{1} << (index & 63))) != 0;
    }

    std::size_t count() const { return count_; }

private:
    static std::size_t pairCount(std::size_t size) {
        return size < 2 ? 0 : size * (size - 1) / 2;
    }

    std::size_t pairIndex(std::size_t left, std::size_t right) const {
        if (left > right) std::swap(left, right);
        return left * (2 * size_ - left - 1) / 2 + (right - left - 1);
    }

    std::size_t size_ = 0;
    std::vector<std::uint64_t> bits_;
    std::size_t count_ = 0;
};

PairOverlapMatrix buildOverlapMatrix(const std::vector<std::string>& rows) {
    PairOverlapMatrix overlap(rows.size());
    std::vector<std::size_t> previous;
    std::vector<std::size_t> participants;
    for (std::size_t column = 0; column < rows.front().size(); ++column) {
        participants.clear();
        for (std::size_t row = 0; row < rows.size(); ++row) {
            if (rows[row][column] != '-') participants.push_back(row);
        }
        if (participants == previous) continue;
        previous = participants;
        for (std::size_t left = 0; left < participants.size(); ++left) {
            for (std::size_t right = left + 1; right < participants.size(); ++right) {
                overlap.mark(participants[left], participants[right]);
            }
        }
    }
    return overlap;
}

std::vector<std::vector<std::size_t>> columnGroups(
    const std::vector<std::string>& rows,
    std::size_t column) {
    std::map<char, std::vector<std::size_t>> by_base;
    for (std::size_t row = 0; row < rows.size(); ++row) {
        const char base = normalizedBase(rows[row][column]);
        if (isConnectableBase(base)) by_base[base].push_back(row);
    }
    std::vector<std::vector<std::size_t>> groups;
    for (auto& [unused, members] : by_base) {
        (void)unused;
        if (members.size() > 1) groups.push_back(std::move(members));
    }
    std::sort(groups.begin(), groups.end());
    return groups;
}

std::vector<std::vector<std::size_t>> componentsForGroup(
    std::size_t row_count,
    const std::vector<std::size_t>& group,
    const std::set<SequencePair>& selected) {
    std::vector<std::size_t> local(row_count, row_count);
    for (std::size_t index = 0; index < group.size(); ++index) {
        local[group[index]] = index;
    }

    DisjointSet dsu(group.size());
    for (const auto& [left, right] : selected) {
        if (local[left] != row_count && local[right] != row_count) {
            dsu.unite(local[left], local[right]);
        }
    }

    std::map<std::size_t, std::vector<std::size_t>> by_root;
    for (std::size_t index = 0; index < group.size(); ++index) {
        by_root[dsu.find(index)].push_back(group[index]);
    }

    std::vector<std::vector<std::size_t>> components;
    components.reserve(by_root.size());
    for (auto& [unused, members] : by_root) {
        (void)unused;
        components.push_back(std::move(members));
    }
    return components;
}

std::size_t bestMember(const std::vector<std::size_t>& members,
                       const std::vector<std::size_t>& coverage,
                       const std::vector<std::string>& names) {
    return *std::min_element(
        members.begin(), members.end(),
        [&](std::size_t left, std::size_t right) {
            if (coverage[left] != coverage[right]) {
                return coverage[left] > coverage[right];
            }
            return names[left] < names[right];
        });
}

void appendOperation(std::string& cigar, char& last_operation,
                     std::uint64_t& last_length, char operation) {
    if (last_operation == operation) {
        ++last_length;
        return;
    }
    if (last_length != 0) {
        cigar += std::to_string(last_length);
        cigar.push_back(last_operation);
    }
    last_operation = operation;
    last_length = 1;
}

}  // namespace

bool verifyColumnConnectivity(
    const std::vector<std::string>& aligned_rows,
    const std::vector<SequencePair>& pairs) {
    validateRows(aligned_rows);
    std::set<SequencePair> selected;
    for (const auto& [left, right] : pairs) {
        if (left >= aligned_rows.size() || right >= aligned_rows.size() ||
            left == right) {
            return false;
        }
        selected.insert(orderedPair(left, right));
    }

    std::vector<std::vector<std::size_t>> previous;
    for (std::size_t column = 0; column < aligned_rows.front().size(); ++column) {
        auto groups = columnGroups(aligned_rows, column);
        if (groups == previous) continue;
        previous = groups;
        for (const auto& group : groups) {
            if (componentsForGroup(
                    aligned_rows.size(), group, selected).size() != 1) {
                return false;
            }
        }
    }
    return true;
}

PairSelectionResult selectPairs(
    const std::vector<std::string>& aligned_rows,
    const std::vector<std::string>& names,
    std::size_t reference_index,
    Mode mode) {
    validateRows(aligned_rows, &names);
    if (reference_index >= aligned_rows.size()) {
        throw std::invalid_argument("PAF reference row is out of range");
    }

    const PairOverlapMatrix overlap = buildOverlapMatrix(aligned_rows);
    PairSelectionResult result;
    result.theoretical_pairs =
        aligned_rows.size() * (aligned_rows.size() - 1) / 2;
    result.eligible_pairs = overlap.count();
    std::set<SequencePair> selected;

    if (mode == Mode::ALL) {
        for (std::size_t left = 0; left < aligned_rows.size(); ++left) {
            for (std::size_t right = left + 1; right < aligned_rows.size(); ++right) {
                if (overlap.contains(left, right)) selected.emplace(left, right);
            }
        }
        result.pairs.assign(selected.begin(), selected.end());
        result.verified = verifyColumnConnectivity(aligned_rows, result.pairs);
        return result;
    }

    for (std::size_t row = 0; row < aligned_rows.size(); ++row) {
        if (row != reference_index && overlap.contains(reference_index, row)) {
            selected.insert(orderedPair(reference_index, row));
        }
    }
    result.base_pairs = selected.size();

    std::vector<std::size_t> coverage(aligned_rows.size(), 0);
    for (std::size_t row = 0; row < aligned_rows.size(); ++row) {
        coverage[row] = static_cast<std::size_t>(std::count_if(
            aligned_rows[row].begin(), aligned_rows[row].end(),
            [](char base) { return base != '-'; }));
    }

    std::vector<std::vector<std::size_t>> previous;
    for (std::size_t column = 0; column < aligned_rows.front().size(); ++column) {
        auto groups = columnGroups(aligned_rows, column);
        if (groups == previous) continue;
        previous = groups;
        for (const auto& group : groups) {
            auto components = componentsForGroup(
                aligned_rows.size(), group, selected);
            if (components.size() < 2) continue;

            const std::size_t hub = bestMember(group, coverage, names);
            auto hub_component = std::find_if(
                components.begin(), components.end(),
                [&](const auto& component) {
                    return std::find(component.begin(), component.end(), hub) !=
                           component.end();
                });
            if (hub_component == components.end()) {
                throw std::runtime_error(
                    "PAF local hub was not found in its component");
            }

            const std::size_t hub_component_index =
                static_cast<std::size_t>(std::distance(
                    components.begin(), hub_component));
            for (std::size_t component_index = 0;
                 component_index < components.size(); ++component_index) {
                if (component_index == hub_component_index) continue;
                const auto& component = components[component_index];
                selected.insert(orderedPair(
                    hub, bestMember(component, coverage, names)));
            }
        }
    }

    result.pairs.assign(selected.begin(), selected.end());
    result.supplemental_pairs = result.pairs.size() - result.base_pairs;
    result.verified = verifyColumnConnectivity(aligned_rows, result.pairs);
    return result;
}

PairProjection projectPair(
    const std::string& target_aligned,
    const std::string& query_aligned,
    bool reverse_columns) {
    if (target_aligned.size() != query_aligned.size() ||
        target_aligned.empty()) {
        throw std::invalid_argument("Cannot project unequal or empty PAF rows");
    }

    PairProjection result;
    char last_operation = '\0';
    std::uint64_t last_length = 0;
    const std::size_t width = target_aligned.size();

    for (std::size_t step = 0; step < width; ++step) {
        const std::size_t column = reverse_columns ? width - step - 1 : step;
        const char target = target_aligned[column];
        const char query = query_aligned[column];
        if (target == '-' && query == '-') continue;

        char operation = '\0';
        if (target == '-') {
            operation = 'I';
            ++result.query_consumed;
            ++result.edit_distance;
        } else if (query == '-') {
            operation = 'D';
            ++result.target_consumed;
            ++result.edit_distance;
        } else {
            ++result.query_consumed;
            ++result.target_consumed;
            ++result.overlap_columns;
            const char normalized_target = normalizedBase(target);
            const char normalized_query = normalizedBase(query);
            if (normalized_target == normalized_query &&
                normalized_target != 'X') {
                operation = '=';
                ++result.matching_bases;
            } else {
                operation = 'X';
                ++result.edit_distance;
            }
        }
        ++result.block_length;
        appendOperation(result.cigar, last_operation, last_length, operation);
    }

    if (last_length != 0) {
        result.cigar += std::to_string(last_length);
        result.cigar.push_back(last_operation);
    }
    result.valid = !result.cigar.empty() && result.overlap_columns != 0 &&
                   result.query_consumed != 0 && result.target_consumed != 0;
    return result;
}

}  // namespace RaMesh::Paf
