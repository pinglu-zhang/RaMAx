#include "window_detector.h"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>

#include "align.h"
#include "spdlog/spdlog.h"

namespace RaMesh::WindowDetection {
namespace {

using Clock = std::chrono::steady_clock;

struct OriginalInterval {
    uint64_t start = 0;
    uint64_t end = 0;
    bool resolved = false;
};

struct RawWindow {
    size_t path_index = 0;
    size_t left_index = 0;
    size_t right_index = 0;
    std::vector<size_t> seed_boundary_indices;
};

uint64_t fnv1a64(const std::string& value) {
    uint64_t hash = 14695981039346656037ULL;
    for (unsigned char c : value) {
        hash ^= static_cast<uint64_t>(c);
        hash *= 1099511628211ULL;
    }
    return hash;
}

std::string hashId(const char prefix, const std::string& value) {
    std::ostringstream out;
    out << prefix << '_' << std::hex << std::setw(16) << std::setfill('0')
        << fnv1a64(value);
    return out.str();
}

std::string strandText(bool reverse) {
    return reverse ? "-" : "+";
}

std::string boolText(bool value) {
    return value ? "1" : "0";
}

std::string tsvSafe(std::string value) {
    for (char& c : value) {
        if (c == '\t' || c == '\n' || c == '\r') {
            c = ' ';
        }
    }
    return value;
}

std::string jsonEscape(const std::string& value) {
    std::ostringstream out;
    for (unsigned char c : value) {
        switch (c) {
        case '"': out << "\\\""; break;
        case '\\': out << "\\\\"; break;
        case '\b': out << "\\b"; break;
        case '\f': out << "\\f"; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default:
            if (c < 0x20) {
                out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                    << static_cast<int>(c) << std::dec;
            } else {
                out << static_cast<char>(c);
            }
        }
    }
    return out.str();
}

template <typename T>
void sortUnique(std::vector<T>& values) {
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
}

std::string join(const std::vector<std::string>& values,
                 const std::string& separator = ",") {
    std::ostringstream out;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i != 0) {
            out << separator;
        }
        out << values[i];
    }
    return out.str();
}

bool contains(const std::vector<std::string>& values, const std::string& value) {
    return std::find(values.begin(), values.end(), value) != values.end();
}

int64_t signedGap(uint64_t left_end, uint64_t right_start) {
    if (right_start >= left_end) {
        const uint64_t gap = right_start - left_end;
        return gap > static_cast<uint64_t>(std::numeric_limits<int64_t>::max())
                   ? std::numeric_limits<int64_t>::max()
                   : static_cast<int64_t>(gap);
    }
    const uint64_t overlap = left_end - right_start;
    return overlap > static_cast<uint64_t>(std::numeric_limits<int64_t>::max())
               ? std::numeric_limits<int64_t>::min()
               : -static_cast<int64_t>(overlap);
}

OriginalInterval toOriginalInterval(
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    const std::string& species,
    const std::string& chromosome,
    uint64_t graph_start,
    uint64_t graph_end) {
    OriginalInterval result;
    if (graph_end <= graph_start) {
        return result;
    }
    const auto manager_it = managers.find(species);
    if (manager_it == managers.end() || !manager_it->second) {
        return result;
    }

    try {
        return std::visit(
            [&](const auto& manager) -> OriginalInterval {
                using ManagerPtr = std::decay_t<decltype(manager)>;
                if (!manager) {
                    return {};
                }
                if constexpr (std::is_same_v<
                                  ManagerPtr,
                                  std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                    const auto seq_id = manager->getSequenceId(chromosome);
                    if (seq_id == SeqPro::SequenceIndex::INVALID_ID) {
                        return {};
                    }
                    const uint64_t start = manager->toOriginalPosition(
                        chromosome, static_cast<SeqPro::Position>(graph_start));
                    const uint64_t last = manager->toOriginalPosition(
                        chromosome, static_cast<SeqPro::Position>(graph_end - 1));
                    if (start == SeqPro::MaskedSequenceManager::INVALID_POSITION ||
                        last == SeqPro::MaskedSequenceManager::INVALID_POSITION ||
                        last < start) {
                        return {};
                    }
                    return {start, last + 1, true};
                } else {
                    if (manager->getSequenceId(chromosome) ==
                        SeqPro::SequenceIndex::INVALID_ID) {
                        return {};
                    }
                    if (graph_end > manager->getSequenceLength(chromosome)) {
                        return {};
                    }
                    return {graph_start, graph_end, true};
                }
            },
            *manager_it->second);
    } catch (const std::exception&) {
        return result;
    }
}

std::vector<const SegmentSnapshot*> segmentsForSpecies(
    const BlockSnapshot& block,
    const std::string& species) {
    std::vector<const SegmentSnapshot*> result;
    for (const auto& segment : block.segments) {
        if (segment.species == species) {
            result.push_back(&segment);
        }
    }
    return result;
}

const SegmentSnapshot* uniqueSegmentForSpecies(
    const BlockSnapshot& block,
    const std::string& species) {
    const SegmentSnapshot* result = nullptr;
    for (const auto& segment : block.segments) {
        if (segment.species != species) {
            continue;
        }
        if (result != nullptr) {
            return nullptr;
        }
        result = &segment;
    }
    return result;
}

std::unordered_map<std::string, const BlockSnapshot*> blockIndex(
    const GraphSnapshot& snapshot) {
    std::unordered_map<std::string, const BlockSnapshot*> result;
    result.reserve(snapshot.blocks.size());
    for (const auto& block : snapshot.blocks) {
        result.emplace(block.block_id, &block);
    }
    return result;
}

std::string anchorLevel(const BlockSnapshot& block, size_t genome_count) {
    const size_t support = block.participating_genomes.size();
    if (support == genome_count) {
        return "ANCHOR_N";
    }
    if (support + 1 >= genome_count && support >= 2) {
        return "ANCHOR_HIGH_K";
    }
    if (support > 2) {
        return "ANCHOR_SUBSET";
    }
    return "ANCHOR_PAIR";
}

bool isReliableAnchor(const BlockSnapshot& block,
                      const Options& options) {
    return !block.copy_ambiguity && block.participating_genomes.size() >= 2 &&
           block.min_segment_length >= options.anchor_min_segment_bp;
}

std::vector<std::string> participantIntersection(
    const BlockSnapshot& left,
    const BlockSnapshot& right) {
    std::vector<std::string> result;
    std::set_intersection(left.participating_genomes.begin(),
                          left.participating_genomes.end(),
                          right.participating_genomes.begin(),
                          right.participating_genomes.end(),
                          std::back_inserter(result));
    return result;
}

BoundaryEvidence buildBoundary(
    const PathSnapshot& path,
    size_t boundary_index,
    const std::unordered_map<std::string, const BlockSnapshot*>& blocks,
    const std::vector<std::string>& input_genomes,
    const Options& options) {
    BoundaryEvidence boundary;
    boundary.reference_species = path.species;
    boundary.reference_chromosome = path.chromosome;
    boundary.path_boundary_index = boundary_index;
    boundary.left_reference_segment = path.segments.at(boundary_index);
    boundary.right_reference_segment = path.segments.at(boundary_index + 1);
    boundary.left_block_id = boundary.left_reference_segment.block_id;
    boundary.right_block_id = boundary.right_reference_segment.block_id;
    boundary.reference_gap = signedGap(boundary.left_reference_segment.graph_end,
                                      boundary.right_reference_segment.graph_start);
    boundary.reference_overlap = boundary.reference_gap < 0
                                     ? static_cast<uint64_t>(-boundary.reference_gap)
                                     : 0;
    boundary.hard_boundary = boundary.reference_gap >
                             static_cast<int64_t>(options.hard_boundary_gap_bp);

    const auto* left = blocks.at(boundary.left_block_id);
    const auto* right = blocks.at(boundary.right_block_id);
    boundary.min_adjacent_block_length =
        std::min(left->min_segment_length, right->min_segment_length);
    boundary.max_adjacent_block_length =
        std::max(left->max_segment_length, right->max_segment_length);
    boundary.copy_ambiguity = left->copy_ambiguity || right->copy_ambiguity;

    bool any_compatible_target = false;
    bool any_order_break = false;
    bool any_target_switch = false;
    bool any_strand_switch = false;

    for (const auto& species : input_genomes) {
        if (species == path.species) {
            continue;
        }
        BoundaryGenomeEvidence evidence;
        evidence.target_species = species;
        const auto left_segments = segmentsForSpecies(*left, species);
        const auto right_segments = segmentsForSpecies(*right, species);
        evidence.left_present = !left_segments.empty();
        evidence.right_present = !right_segments.empty();
        evidence.unique_target = left_segments.size() == 1 &&
                                 right_segments.size() == 1;

        if (left_segments.size() == 1) {
            const auto& segment = *left_segments.front();
            evidence.left_chromosome = segment.chromosome;
            evidence.left_graph_start = segment.graph_start;
            evidence.left_graph_end = segment.graph_end;
            evidence.left_original_start = segment.original_start;
            evidence.left_original_end = segment.original_end;
            evidence.left_original_resolved =
                segment.original_coordinates_resolved;
            evidence.left_reverse =
                segment.reverse != boundary.left_reference_segment.reverse;
        }
        if (right_segments.size() == 1) {
            const auto& segment = *right_segments.front();
            evidence.right_chromosome = segment.chromosome;
            evidence.right_graph_start = segment.graph_start;
            evidence.right_graph_end = segment.graph_end;
            evidence.right_original_start = segment.original_start;
            evidence.right_original_end = segment.original_end;
            evidence.right_original_resolved =
                segment.original_coordinates_resolved;
            evidence.right_reverse =
                segment.reverse != boundary.right_reference_segment.reverse;
        }

        if (evidence.unique_target) {
            const auto& left_target = *left_segments.front();
            const auto& right_target = *right_segments.front();
            evidence.target_consistent =
                left_target.chromosome == right_target.chromosome;

            const bool left_relative_reverse =
                boundary.left_reference_segment.reverse != left_target.reverse;
            const bool right_relative_reverse =
                boundary.right_reference_segment.reverse != right_target.reverse;
            evidence.strand_consistent =
                left_relative_reverse == right_relative_reverse;

            if (evidence.target_consistent && evidence.strand_consistent) {
                if (!left_relative_reverse) {
                    evidence.order_consistent =
                        left_target.graph_start <= right_target.graph_start;
                    evidence.target_gap = signedGap(left_target.graph_end,
                                                    right_target.graph_start);
                } else {
                    evidence.order_consistent =
                        left_target.graph_start >= right_target.graph_start;
                    evidence.target_gap = signedGap(right_target.graph_end,
                                                    left_target.graph_start);
                }
                evidence.target_gap_resolved = true;
            }

            any_target_switch = any_target_switch || !evidence.target_consistent;
            any_strand_switch = any_strand_switch ||
                                (evidence.target_consistent &&
                                 !evidence.strand_consistent);
            any_order_break = any_order_break ||
                              (evidence.target_consistent &&
                               evidence.strand_consistent &&
                               !evidence.order_consistent);
            any_compatible_target = any_compatible_target ||
                                    (evidence.target_consistent &&
                                     evidence.strand_consistent &&
                                     evidence.order_consistent);
        } else if (left_segments.size() > 1 || right_segments.size() > 1) {
            boundary.copy_ambiguity = true;
        }
        boundary.genomes.push_back(std::move(evidence));
    }

    if (boundary.min_adjacent_block_length <= options.micro_block_max_bp) {
        boundary.signals.push_back("MICRO_BLOCK");
    }
    if (boundary.min_adjacent_block_length > options.micro_block_max_bp &&
        boundary.min_adjacent_block_length <= options.short_block_max_bp &&
        ((boundary.reference_gap >= 0 &&
          boundary.reference_gap <=
              static_cast<int64_t>(options.primary_gap_max_bp)) ||
         any_order_break || any_target_switch || any_strand_switch)) {
        boundary.signals.push_back("SHORT_BLOCK");
    }
    if (boundary.reference_gap >= 1 &&
        boundary.reference_gap <=
            static_cast<int64_t>(options.primary_gap_max_bp) &&
        any_compatible_target) {
        boundary.signals.push_back("SHORT_GAP");
    }
    if (boundary.reference_gap >
            static_cast<int64_t>(options.primary_gap_max_bp) &&
        boundary.reference_gap <=
            static_cast<int64_t>(options.extended_gap_max_bp)) {
        boundary.signals.push_back("EXTENDED_GAP");
    }
    if (any_order_break) {
        boundary.signals.push_back("ORDER_BREAK");
    }
    if (any_target_switch) {
        boundary.signals.push_back("TARGET_SWITCH");
    }
    if (any_strand_switch) {
        boundary.signals.push_back("STRAND_SWITCH");
    }
    if (boundary.copy_ambiguity) {
        boundary.signals.push_back("COPY_AMBIGUITY");
    }
    sortUnique(boundary.signals);

    std::ostringstream signature;
    signature << path.species << '|' << path.chromosome << '|'
              << boundary.left_block_id << '|'
              << boundary.left_reference_segment.graph_start << '|'
              << boundary.left_reference_segment.graph_end << '|'
              << boundary.right_block_id << '|'
              << boundary.right_reference_segment.graph_start << '|'
              << boundary.right_reference_segment.graph_end;
    boundary.boundary_id = hashId('Y', signature.str());
    return boundary;
}

void addDropoutAndTransientSignals(
    const PathSnapshot& path,
    const std::unordered_map<std::string, const BlockSnapshot*>& blocks,
    std::vector<BoundaryEvidence>& boundaries) {
    if (path.segments.size() < 3) {
        return;
    }
    for (size_t center_index = 1; center_index + 1 < path.segments.size();
         ++center_index) {
        const auto* left = blocks.at(path.segments[center_index - 1].block_id);
        const auto* center = blocks.at(path.segments[center_index].block_id);
        const auto* right = blocks.at(path.segments[center_index + 1].block_id);

        const auto outer_common = participantIntersection(*left, *right);
        for (const auto& species : outer_common) {
            if (!std::binary_search(center->participating_genomes.begin(),
                                    center->participating_genomes.end(), species)) {
                boundaries[center_index - 1].dropout_genomes.push_back(species);
                boundaries[center_index].dropout_genomes.push_back(species);
                boundaries[center_index - 1].signals.push_back("GENOME_DROPOUT");
                boundaries[center_index].signals.push_back("GENOME_DROPOUT");
                boundaries[center_index - 1].signals.push_back("SUPPORT_VOLATILITY");
                boundaries[center_index].signals.push_back("SUPPORT_VOLATILITY");
            }

            const auto* left_segment = uniqueSegmentForSpecies(*left, species);
            const auto* center_segment = uniqueSegmentForSpecies(*center, species);
            const auto* right_segment = uniqueSegmentForSpecies(*right, species);
            if (left_segment && center_segment && right_segment &&
                left_segment->chromosome == right_segment->chromosome &&
                center_segment->chromosome != left_segment->chromosome) {
                boundaries[center_index - 1].signals.push_back(
                    "TRANSIENT_TARGET_SWITCH");
                boundaries[center_index].signals.push_back(
                    "TRANSIENT_TARGET_SWITCH");
            }
        }
    }
    for (auto& boundary : boundaries) {
        sortUnique(boundary.dropout_genomes);
        sortUnique(boundary.signals);
    }
}

std::string initialPriority(const BoundaryEvidence& boundary,
                            const Options& options) {
    if (boundary.hard_boundary || boundary.copy_ambiguity) {
        return "PROTECTED";
    }
    if (contains(boundary.signals, "MICRO_BLOCK") ||
        (contains(boundary.signals, "ORDER_BREAK") &&
         boundary.min_adjacent_block_length <= options.short_block_max_bp) ||
        contains(boundary.signals, "TRANSIENT_TARGET_SWITCH") ||
        (contains(boundary.signals, "SHORT_GAP") &&
         contains(boundary.signals, "GENOME_DROPOUT"))) {
        return "A";
    }
    if (contains(boundary.signals, "SHORT_BLOCK") ||
        contains(boundary.signals, "EXTENDED_GAP") ||
        contains(boundary.signals, "TARGET_SWITCH") ||
        contains(boundary.signals, "GENOME_DROPOUT")) {
        return "B";
    }
    return "C";
}

std::vector<RawWindow> buildRawWindows(
    const GraphSnapshot& snapshot,
    const std::vector<std::vector<BoundaryEvidence>>& path_boundaries,
    const std::unordered_map<std::string, const BlockSnapshot*>& blocks,
    const Options& options) {
    std::vector<RawWindow> result;
    for (size_t path_index = 0; path_index < snapshot.paths.size(); ++path_index) {
        const auto& path = snapshot.paths[path_index];
        const auto& boundaries = path_boundaries[path_index];
        if (path.segments.size() < 2) {
            continue;
        }
        std::vector<RawWindow> path_windows;
        for (size_t boundary_index = 0; boundary_index < boundaries.size();
             ++boundary_index) {
            if (boundaries[boundary_index].signals.empty()) {
                continue;
            }

            size_t left_index = boundary_index;
            while (!isReliableAnchor(*blocks.at(path.segments[left_index].block_id),
                                     options)) {
                if (left_index == 0 || boundaries[left_index - 1].hard_boundary) {
                    break;
                }
                --left_index;
            }

            size_t right_index = boundary_index + 1;
            while (!isReliableAnchor(*blocks.at(path.segments[right_index].block_id),
                                     options)) {
                if (right_index + 1 >= path.segments.size() ||
                    boundaries[right_index].hard_boundary) {
                    break;
                }
                ++right_index;
            }

            path_windows.push_back(
                {path_index, left_index, right_index, {boundary_index}});
        }

        std::sort(path_windows.begin(), path_windows.end(),
                  [](const RawWindow& left, const RawWindow& right) {
                      if (left.left_index != right.left_index) {
                          return left.left_index < right.left_index;
                      }
                      return left.right_index < right.right_index;
                  });
        for (auto& window : path_windows) {
            if (!result.empty() && result.back().path_index == path_index &&
                window.left_index < result.back().right_index) {
                result.back().left_index =
                    std::min(result.back().left_index, window.left_index);
                result.back().right_index =
                    std::max(result.back().right_index, window.right_index);
                result.back().seed_boundary_indices.insert(
                    result.back().seed_boundary_indices.end(),
                    window.seed_boundary_indices.begin(),
                    window.seed_boundary_indices.end());
                std::sort(result.back().seed_boundary_indices.begin(),
                          result.back().seed_boundary_indices.end());
                result.back().seed_boundary_indices.erase(
                    std::unique(result.back().seed_boundary_indices.begin(),
                                result.back().seed_boundary_indices.end()),
                    result.back().seed_boundary_indices.end());
            } else {
                result.push_back(std::move(window));
            }
        }
    }
    return result;
}

struct GenomeAnchorPair {
    const SegmentSnapshot* left = nullptr;
    const SegmentSnapshot* right = nullptr;
    std::string invalid_reason;
};

bool anchorPairsCompatible(const GenomeAnchorPair& first,
                           const GenomeAnchorPair& second) {
    if (!first.left || !first.right || !second.left || !second.right ||
        !first.invalid_reason.empty() || !second.invalid_reason.empty()) {
        return false;
    }
    const bool first_forward =
        first.left->graph_start <= first.right->graph_start;
    const bool second_forward =
        second.left->graph_start <= second.right->graph_start;
    const bool left_relative_reverse =
        first.left->reverse != second.left->reverse;
    const bool right_relative_reverse =
        first.right->reverse != second.right->reverse;
    if (left_relative_reverse != right_relative_reverse) {
        return false;
    }
    return left_relative_reverse ? first_forward != second_forward
                                 : first_forward == second_forward;
}

struct CompatibleSubsetResult {
    std::vector<std::string> genomes;
    bool truncated = false;
};

CompatibleSubsetResult maximumCompatibleSubset(
    const std::vector<std::string>& input_genomes,
    const std::map<std::string, GenomeAnchorPair>& anchor_pairs,
    size_t search_budget) {
    std::vector<std::string> valid;
    for (const auto& species : input_genomes) {
        const auto found = anchor_pairs.find(species);
        if (found != anchor_pairs.end() && found->second.invalid_reason.empty()) {
            valid.push_back(species);
        }
    }

    std::map<std::pair<std::string, std::string>, bool> compatible;
    for (size_t i = 0; i < valid.size(); ++i) {
        for (size_t j = i + 1; j < valid.size(); ++j) {
            compatible[{valid[i], valid[j]}] = anchorPairsCompatible(
                anchor_pairs.at(valid[i]), anchor_pairs.at(valid[j]));
        }
    }
    auto pairCompatible = [&](const std::string& left,
                              const std::string& right) {
        if (left == right) {
            return true;
        }
        return compatible.at(std::minmax(left, right));
    };

    CompatibleSubsetResult result;
    // Deterministic greedy solutions provide a useful lower bound and a stable
    // fallback if branch-and-bound reaches its configured budget.
    for (size_t start = 0; start < valid.size(); ++start) {
        std::vector<std::string> candidate;
        candidate.push_back(valid[start]);
        for (size_t index = 0; index < valid.size(); ++index) {
            if (index == start) {
                continue;
            }
            if (std::all_of(candidate.begin(), candidate.end(),
                            [&](const std::string& member) {
                                return pairCompatible(member, valid[index]);
                            })) {
                candidate.push_back(valid[index]);
            }
        }
        std::sort(candidate.begin(), candidate.end());
        if (candidate.size() > result.genomes.size() ||
            (candidate.size() == result.genomes.size() &&
             candidate < result.genomes)) {
            result.genomes = std::move(candidate);
        }
    }

    size_t visited = 0;
    std::function<void(std::vector<std::string>&,
                       const std::vector<std::string>&)>
        search;
    search = [&](std::vector<std::string>& current,
                 const std::vector<std::string>& candidates) {
        if (result.truncated) {
            return;
        }
        if (++visited > search_budget) {
            result.truncated = true;
            return;
        }
        if (current.size() + candidates.size() < result.genomes.size()) {
            return;
        }
        if (candidates.empty()) {
            auto sorted = current;
            std::sort(sorted.begin(), sorted.end());
            if (sorted.size() > result.genomes.size() ||
                (sorted.size() == result.genomes.size() &&
                 sorted < result.genomes)) {
                result.genomes = std::move(sorted);
            }
            return;
        }

        const std::string selected = candidates.front();
        std::vector<std::string> remaining(candidates.begin() + 1,
                                           candidates.end());
        std::vector<std::string> compatible_remaining;
        compatible_remaining.reserve(remaining.size());
        for (const auto& candidate : remaining) {
            if (pairCompatible(selected, candidate)) {
                compatible_remaining.push_back(candidate);
            }
        }
        current.push_back(selected);
        search(current, compatible_remaining);
        current.pop_back();
        search(current, remaining);
    };

    if (!valid.empty() && search_budget > 0) {
        std::vector<std::string> current;
        search(current, valid);
    } else if (!valid.empty()) {
        result.truncated = true;
    }
    sortUnique(result.genomes);
    return result;
}

CandidateWindow materializeWindow(
    const RawWindow& raw,
    const GraphSnapshot& snapshot,
    const std::vector<std::vector<BoundaryEvidence>>& path_boundaries,
    const std::unordered_map<std::string, const BlockSnapshot*>& blocks,
    const Options& options) {
    const auto& path = snapshot.paths.at(raw.path_index);
    const auto& boundaries = path_boundaries.at(raw.path_index);
    CandidateWindow window;
    window.round_id = snapshot.round_id;
    window.current_reference = snapshot.current_reference;
    window.report_species = path.species;
    window.report_chromosome = path.chromosome;
    window.path_left_index = raw.left_index;
    window.path_right_index = raw.right_index;
    window.input_genome_count = snapshot.input_genomes.size();
    window.graph_start = path.segments.at(raw.left_index).graph_start;
    window.graph_end = path.segments.at(raw.right_index).graph_end;
    window.left_anchor_id = path.segments.at(raw.left_index).block_id;
    window.right_anchor_id = path.segments.at(raw.right_index).block_id;

    const auto* left_anchor = blocks.at(window.left_anchor_id);
    const auto* right_anchor = blocks.at(window.right_anchor_id);
    const bool left_reliable = isReliableAnchor(*left_anchor, options);
    const bool right_reliable = isReliableAnchor(*right_anchor, options);
    window.has_two_sided_anchors = left_reliable && right_reliable;
    window.left_anchor_level = left_reliable
                                   ? anchorLevel(*left_anchor,
                                                 snapshot.input_genomes.size())
                                   : "NONE";
    window.right_anchor_level = right_reliable
                                    ? anchorLevel(*right_anchor,
                                                  snapshot.input_genomes.size())
                                    : "NONE";

    window.original_coordinates_resolved = true;
    window.original_start = std::numeric_limits<uint64_t>::max();
    window.original_end = 0;
    window.min_block_length = std::numeric_limits<uint64_t>::max();
    window.max_reference_gap = std::numeric_limits<int64_t>::min();

    for (size_t index = raw.left_index; index <= raw.right_index; ++index) {
        const auto& segment = path.segments[index];
        const auto* block = blocks.at(segment.block_id);
        window.block_ids.push_back(segment.block_id);
        window.min_block_length =
            std::min(window.min_block_length, block->min_segment_length);
        window.max_block_length =
            std::max(window.max_block_length, block->max_segment_length);
        if (!segment.original_coordinates_resolved) {
            window.original_coordinates_resolved = false;
        } else {
            window.original_start =
                std::min(window.original_start, segment.original_start);
            window.original_end =
                std::max(window.original_end, segment.original_end);
        }
    }
    sortUnique(window.block_ids);
    if (!window.original_coordinates_resolved) {
        window.original_start = 0;
        window.original_end = 0;
        window.complex_flags.push_back("COORDINATE_UNRESOLVED");
    }

    for (size_t boundary_index = raw.left_index;
         boundary_index < raw.right_index; ++boundary_index) {
        const auto& boundary = boundaries.at(boundary_index);
        window.boundary_ids.push_back(boundary.boundary_id);
        window.max_reference_gap =
            std::max(window.max_reference_gap, boundary.reference_gap);
        if (boundary.hard_boundary) {
            window.complex_flags.push_back("HARD_BOUNDARY");
        }
        if (boundary.copy_ambiguity) {
            window.complex_flags.push_back("COPY_AMBIGUITY");
        }
    }
    for (size_t boundary_index : raw.seed_boundary_indices) {
        const auto& boundary = boundaries.at(boundary_index);
        window.seed_boundary_ids.push_back(boundary.boundary_id);
        window.signals.insert(window.signals.end(), boundary.signals.begin(),
                              boundary.signals.end());
    }
    sortUnique(window.boundary_ids);
    sortUnique(window.seed_boundary_ids);
    sortUnique(window.signals);

    const uint64_t span = window.graph_end >= window.graph_start
                              ? window.graph_end - window.graph_start
                              : 0;
    if (span > options.max_window_span_bp) {
        window.complex_flags.push_back("COMPLEX_LONG_WINDOW");
    }
    if (!window.has_two_sided_anchors) {
        window.complex_flags.push_back("ONE_SIDED_ANCHOR");
    }
    sortUnique(window.complex_flags);

    std::map<std::string, GenomeAnchorPair> anchor_pairs;
    for (const auto& species : snapshot.input_genomes) {
        WindowGenomeEvidence evidence;
        evidence.species = species;
        const auto left_segments = segmentsForSpecies(*left_anchor, species);
        const auto right_segments = segmentsForSpecies(*right_anchor, species);
        evidence.left_anchor_present = left_segments.size() == 1;
        evidence.right_anchor_present = right_segments.size() == 1;

        std::set<std::string> chromosomes;
        bool all_original_resolved = true;
        uint64_t graph_start = std::numeric_limits<uint64_t>::max();
        uint64_t graph_end = 0;
        uint64_t original_start = std::numeric_limits<uint64_t>::max();
        uint64_t original_end = 0;
        for (const auto& block_id : window.block_ids) {
            const auto* block = blocks.at(block_id);
            for (const auto& segment : block->segments) {
                if (segment.species != species) {
                    continue;
                }
                evidence.currently_aligned = true;
                chromosomes.insert(segment.chromosome);
                graph_start = std::min(graph_start, segment.graph_start);
                graph_end = std::max(graph_end, segment.graph_end);
                if (!segment.original_coordinates_resolved) {
                    all_original_resolved = false;
                } else {
                    original_start =
                        std::min(original_start, segment.original_start);
                    original_end = std::max(original_end, segment.original_end);
                }
            }
        }
        if (evidence.currently_aligned) {
            evidence.graph_start = graph_start;
            evidence.graph_end = graph_end;
            if (chromosomes.size() == 1) {
                evidence.chromosome = *chromosomes.begin();
            } else {
                evidence.chromosome = "MULTIPLE";
                all_original_resolved = false;
            }
            if (all_original_resolved) {
                evidence.original_start = original_start;
                evidence.original_end = original_end;
                evidence.original_coordinates_resolved = true;
            }
        }

        GenomeAnchorPair anchor_pair;
        if (left_segments.size() != 1) {
            anchor_pair.invalid_reason = left_segments.empty()
                                             ? "NO_LEFT_ANCHOR"
                                             : "COPY_AMBIGUITY";
        } else if (right_segments.size() != 1) {
            anchor_pair.invalid_reason = right_segments.empty()
                                             ? "NO_RIGHT_ANCHOR"
                                             : "COPY_AMBIGUITY";
        } else {
            anchor_pair.left = left_segments.front();
            anchor_pair.right = right_segments.front();
            if (anchor_pair.left->chromosome != anchor_pair.right->chromosome) {
                anchor_pair.invalid_reason = "TARGET_CONFLICT";
            }
        }
        evidence.excluded_reason = anchor_pair.invalid_reason;
        anchor_pairs.emplace(species, anchor_pair);

        if (evidence.currently_aligned) {
            ++window.currently_aligned_genome_count;
        }
        window.genomes.push_back(std::move(evidence));
    }

    const auto compatible_subset = maximumCompatibleSubset(
        snapshot.input_genomes, anchor_pairs, options.subset_search_budget);
    window.max_k_genomes = compatible_subset.genomes;
    window.max_possible_k = window.max_k_genomes.size();
    if (compatible_subset.truncated) {
        window.complex_flags.push_back("SUBSET_SEARCH_TRUNCATED");
    }

    const std::string baseline_species =
        window.max_k_genomes.empty() ? std::string{} : window.max_k_genomes.front();
    for (auto& evidence : window.genomes) {
        evidence.included_in_max_k = std::binary_search(
            window.max_k_genomes.begin(), window.max_k_genomes.end(),
            evidence.species);
        if (evidence.included_in_max_k) {
            const auto& pair = anchor_pairs.at(evidence.species);
            if (evidence.species == baseline_species) {
                evidence.strand = "+";
            } else {
                const auto& baseline = anchor_pairs.at(baseline_species);
                evidence.strand =
                    (baseline.left->reverse != pair.left->reverse) ? "-" : "+";
            }
        } else {
            window.excluded_genomes.push_back(evidence.species);
            if (evidence.excluded_reason.empty()) {
                evidence.excluded_reason = "COMPATIBILITY_CONFLICT";
                const auto& excluded_pair = anchor_pairs.at(evidence.species);
                for (const auto& included_species : window.max_k_genomes) {
                    const auto& included_pair = anchor_pairs.at(included_species);
                    if (anchorPairsCompatible(excluded_pair, included_pair)) {
                        continue;
                    }
                    const bool left_relative_reverse =
                        excluded_pair.left->reverse != included_pair.left->reverse;
                    const bool right_relative_reverse =
                        excluded_pair.right->reverse != included_pair.right->reverse;
                    evidence.excluded_reason =
                        left_relative_reverse != right_relative_reverse
                            ? "STRAND_CONFLICT"
                            : "ORDER_CONFLICT";
                    break;
                }
            }
        }
    }
    sortUnique(window.excluded_genomes);
    sortUnique(window.complex_flags);

    const bool has_dropout = contains(window.signals, "GENOME_DROPOUT");
    const bool high_k = window.max_possible_k == window.input_genome_count ||
                        (window.max_possible_k + 1 ==
                             window.input_genome_count &&
                         has_dropout);
    const bool tier_a_signal =
        contains(window.signals, "MICRO_BLOCK") ||
        (contains(window.signals, "ORDER_BREAK") &&
         window.min_block_length <= options.short_block_max_bp) ||
        (contains(window.signals, "TRANSIENT_TARGET_SWITCH") &&
         window.min_block_length <= options.short_block_max_bp) ||
        (contains(window.signals, "SHORT_GAP") && has_dropout);

    if (contains(window.complex_flags, "HARD_BOUNDARY") ||
        contains(window.complex_flags, "COPY_AMBIGUITY") ||
        contains(window.complex_flags, "COORDINATE_UNRESOLVED") ||
        contains(window.complex_flags, "COMPLEX_LONG_WINDOW") ||
        contains(window.complex_flags, "SUBSET_SEARCH_TRUNCATED")) {
        window.priority_tier = "PROTECTED";
    } else if (window.has_two_sided_anchors && high_k && tier_a_signal) {
        window.priority_tier = "A";
    } else if (window.has_two_sided_anchors &&
               (contains(window.signals, "SHORT_BLOCK") ||
                contains(window.signals, "EXTENDED_GAP") ||
                contains(window.signals, "TARGET_SWITCH") ||
                has_dropout || window.max_possible_k + 1 >=
                                   window.input_genome_count)) {
        window.priority_tier = "B";
    } else {
        window.priority_tier = "C";
    }

    std::ostringstream signature;
    signature << snapshot.round_id << '|' << window.report_species << '|'
              << window.report_chromosome << '|' << window.left_anchor_id << '|'
              << window.right_anchor_id << '|' << join(window.block_ids, ";")
              << '|' << join(window.seed_boundary_ids, ";");
    window.window_id = hashId('W', signature.str());
    window.detector_reason = join(window.signals, ";");
    return window;
}

std::string windowCanonicalKey(const CandidateWindow& window) {
    std::vector<std::string> anchors{window.left_anchor_id,
                                     window.right_anchor_id};
    sortUnique(anchors);
    return join(anchors, "|") + "|" + join(window.block_ids, ";");
}

bool bridgeSignal(const BoundaryEvidence& boundary) {
    return contains(boundary.signals, "MICRO_BLOCK") ||
           contains(boundary.signals, "SHORT_BLOCK") ||
           contains(boundary.signals, "SHORT_GAP") ||
           contains(boundary.signals, "ORDER_BREAK") ||
           contains(boundary.signals, "TRANSIENT_TARGET_SWITCH") ||
           contains(boundary.signals, "GENOME_DROPOUT");
}

void writeFile(const std::filesystem::path& path,
               const std::function<void(std::ostream&)>& writer) {
    std::filesystem::create_directories(path.parent_path());
    const auto partial = path.string() + ".partial";
    {
        std::ofstream output(partial, std::ios::binary | std::ios::trunc);
        if (!output) {
            throw std::runtime_error("Cannot open window report file: " + partial);
        }
        writer(output);
        output.flush();
        if (!output) {
            throw std::runtime_error("Failed writing window report file: " +
                                     partial);
        }
    }
    std::error_code ec;
    std::filesystem::remove(path, ec);
    ec.clear();
    std::filesystem::rename(partial, path, ec);
    if (ec) {
        throw std::runtime_error("Cannot finalize window report file " +
                                 path.string() + ": " + ec.message());
    }
}

std::string roundName(uint64_t round_id) {
    std::ostringstream out;
    out << "round_" << std::setw(3) << std::setfill('0') << round_id;
    return out.str();
}

std::map<std::string, size_t> tierCounts(
    const std::vector<CandidateWindow>& windows) {
    std::map<std::string, size_t> result;
    for (const auto& window : windows) {
        ++result[window.priority_tier];
    }
    return result;
}

}  // namespace

std::string detectionModeToString(DetectionMode mode) {
    return mode == DetectionMode::FINAL_ONLY ? "final-only" : "each-round";
}

DetectionMode detectionModeFromString(const std::string& value) {
    if (value == "each-round") {
        return DetectionMode::EACH_ROUND;
    }
    if (value == "final-only") {
        return DetectionMode::FINAL_ONLY;
    }
    throw std::invalid_argument(
        "window detection mode must be each-round or final-only");
}

GraphSnapshot snapshotGraph(
    const RaMeshMultiGenomeGraph& graph,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    uint64_t round_id,
    const SpeciesName& current_reference) {
    GraphSnapshot snapshot;
    snapshot.round_id = round_id;
    snapshot.current_reference = current_reference;

    struct TemporaryBlock {
        BlockPtr owner;
        const Block* raw = nullptr;
        BlockSnapshot snapshot;
        std::vector<const Segment*> raw_segments;
    };

    std::shared_lock graph_lock(graph.rw);
    for (const auto& [species, unused] : graph.species_graphs) {
        (void)unused;
        snapshot.input_genomes.push_back(species);
    }
    std::sort(snapshot.input_genomes.begin(), snapshot.input_genomes.end());

    std::vector<TemporaryBlock> temporary_blocks;
    temporary_blocks.reserve(graph.blocks.size());
    for (const auto& weak_block : graph.blocks) {
        auto block = weak_block.lock();
        if (!block) {
            continue;
        }
        TemporaryBlock temporary;
        temporary.owner = block;
        temporary.raw = block.get();
        {
            std::shared_lock block_lock(block->rw);
            temporary.snapshot.diagnostic_ref_chromosome = block->ref_chr;
            std::vector<std::pair<SpeciesChrPair, SegPtr>> anchors;
            anchors.reserve(block->anchors.size());
            for (const auto& anchor : block->anchors) {
                anchors.push_back(anchor);
            }
            std::sort(anchors.begin(), anchors.end(),
                      [](const auto& left, const auto& right) {
                          return left.first < right.first;
                      });
            for (const auto& [species_chr, segment] : anchors) {
                if (!segment || !segment->isSegment() || segment->length == 0) {
                    snapshot.audit_records.push_back(
                        {"ERROR", "INVALID_BLOCK_ANCHOR", species_chr.first,
                         species_chr.second, segment ? segment->start : 0,
                         "Block contains a null, sentinel, or zero-length anchor"});
                    continue;
                }
                SegmentSnapshot copied;
                {
                    std::shared_lock segment_lock(segment->rw);
                    copied.species = species_chr.first;
                    copied.chromosome = species_chr.second;
                    copied.graph_start = segment->start;
                    copied.graph_end = static_cast<uint64_t>(segment->start) +
                                       static_cast<uint64_t>(segment->length);
                    copied.reverse = segment->strand == Strand::REVERSE;
                    copied.primary = segment->isPrimary();
                    copied.left_extended = segment->left_extend;
                    copied.right_extended = segment->right_extend;
                    copied.cigar_summary = cigarToString(segment->cigar);
                    if (segment->parent_block.get() != block.get()) {
                        snapshot.audit_records.push_back(
                            {"ERROR", "PARENT_BLOCK_MISMATCH",
                             species_chr.first, species_chr.second,
                             segment->start,
                             "Block anchor Segment points to a different parent Block"});
                    }
                }
                const auto original = toOriginalInterval(
                    managers, copied.species, copied.chromosome,
                    copied.graph_start, copied.graph_end);
                copied.original_start = original.start;
                copied.original_end = original.end;
                copied.original_coordinates_resolved = original.resolved;
                temporary.snapshot.segments.push_back(std::move(copied));
                temporary.raw_segments.push_back(segment.get());
            }
        }

        std::vector<size_t> order(temporary.snapshot.segments.size());
        for (size_t i = 0; i < order.size(); ++i) {
            order[i] = i;
        }
        std::sort(order.begin(), order.end(), [&](size_t left, size_t right) {
            const auto& a = temporary.snapshot.segments[left];
            const auto& b = temporary.snapshot.segments[right];
            return std::tie(a.species, a.chromosome, a.graph_start, a.graph_end,
                            a.reverse, a.primary) <
                   std::tie(b.species, b.chromosome, b.graph_start, b.graph_end,
                            b.reverse, b.primary);
        });

        std::vector<SegmentSnapshot> sorted_segments;
        std::vector<const Segment*> sorted_raw;
        sorted_segments.reserve(order.size());
        sorted_raw.reserve(order.size());
        std::ostringstream signature;
        for (size_t index : order) {
            auto segment = std::move(temporary.snapshot.segments[index]);
            signature << segment.species << ':' << segment.chromosome << ':'
                      << segment.graph_start << '-' << segment.graph_end << ':'
                      << strandText(segment.reverse) << ':'
                      << (segment.primary ? 'P' : 'S') << '|';
            sorted_segments.push_back(std::move(segment));
            sorted_raw.push_back(temporary.raw_segments[index]);
        }
        temporary.snapshot.segments = std::move(sorted_segments);
        temporary.raw_segments = std::move(sorted_raw);
        temporary.snapshot.canonical_signature = signature.str();
        temporary_blocks.push_back(std::move(temporary));
    }

    std::stable_sort(
        temporary_blocks.begin(), temporary_blocks.end(),
        [](const TemporaryBlock& left, const TemporaryBlock& right) {
            return left.snapshot.canonical_signature <
                   right.snapshot.canonical_signature;
        });

    std::unordered_map<const Segment*, std::pair<size_t, size_t>> segment_lookup;
    std::string previous_signature;
    size_t duplicate_ordinal = 0;
    bool has_previous_signature = false;
    for (auto& temporary : temporary_blocks) {
        if (has_previous_signature &&
            temporary.snapshot.canonical_signature == previous_signature) {
            ++duplicate_ordinal;
            snapshot.audit_records.push_back(
                {"WARNING", "DUPLICATE_BLOCK_SIGNATURE", "", "", 0,
                 "Multiple active Blocks have the same canonical signature"});
        } else {
            previous_signature = temporary.snapshot.canonical_signature;
            duplicate_ordinal = 0;
            has_previous_signature = true;
        }
        temporary.snapshot.block_id =
            hashId('B', temporary.snapshot.canonical_signature);
        if (duplicate_ordinal != 0) {
            temporary.snapshot.block_id += "_" +
                                           std::to_string(duplicate_ordinal);
        }

        std::map<std::string, size_t> species_counts;
        temporary.snapshot.min_segment_length =
            std::numeric_limits<uint64_t>::max();
        for (size_t segment_index = 0;
             segment_index < temporary.snapshot.segments.size();
             ++segment_index) {
            auto& segment = temporary.snapshot.segments[segment_index];
            segment.block_id = temporary.snapshot.block_id;
            segment.segment_id = hashId(
                'S', segment.block_id + '|' + segment.species + '|' +
                         segment.chromosome + '|' +
                         std::to_string(segment.graph_start) + '|' +
                         std::to_string(segment.graph_end) + '|' +
                         strandText(segment.reverse));
            ++species_counts[segment.species];
            temporary.snapshot.min_segment_length = std::min(
                temporary.snapshot.min_segment_length,
                segment.graph_end - segment.graph_start);
            temporary.snapshot.max_segment_length = std::max(
                temporary.snapshot.max_segment_length,
                segment.graph_end - segment.graph_start);
            temporary.snapshot.has_secondary =
                temporary.snapshot.has_secondary || !segment.primary;
        }
        if (temporary.snapshot.segments.empty()) {
            temporary.snapshot.min_segment_length = 0;
        }
        for (const auto& [species, count] : species_counts) {
            temporary.snapshot.participating_genomes.push_back(species);
            temporary.snapshot.copy_ambiguity =
                temporary.snapshot.copy_ambiguity || count > 1;
        }

        const size_t block_index = snapshot.blocks.size();
        snapshot.blocks.push_back(std::move(temporary.snapshot));
        for (size_t segment_index = 0;
             segment_index < temporary.raw_segments.size(); ++segment_index) {
            segment_lookup.emplace(temporary.raw_segments[segment_index],
                                   std::make_pair(block_index, segment_index));
        }
    }

    std::vector<std::pair<std::string, const RaMeshGenomeGraph*>> genomes;
    genomes.reserve(graph.species_graphs.size());
    for (const auto& [species, genome] : graph.species_graphs) {
        genomes.emplace_back(species, &genome);
    }
    std::sort(genomes.begin(), genomes.end(),
              [](const auto& left, const auto& right) {
                  return left.first < right.first;
              });

    for (const auto& [species, genome] : genomes) {
        std::shared_lock genome_lock(genome->rw);
        std::vector<std::pair<std::string, const GenomeEnd*>> chromosomes;
        chromosomes.reserve(genome->chr2end.size());
        for (const auto& [chromosome, end] : genome->chr2end) {
            chromosomes.emplace_back(chromosome, &end);
        }
        std::sort(chromosomes.begin(), chromosomes.end(),
                  [](const auto& left, const auto& right) {
                      return left.first < right.first;
                  });

        for (const auto& [chromosome, end] : chromosomes) {
            std::shared_lock end_lock(end->rw);
            PathSnapshot path;
            path.species = species;
            path.chromosome = chromosome;
            std::unordered_set<const Segment*> seen;
            SegPtr previous = end->head;
            SegPtr current = end->head->primary_path.next.load(
                std::memory_order_acquire);
            while (current && !current->isTail()) {
                if (!seen.insert(current.get()).second) {
                    snapshot.audit_records.push_back(
                        {"CRITICAL", "PATH_CYCLE", species, chromosome,
                         current->start, "Cycle detected in primary path"});
                    break;
                }
                const auto actual_previous = current->primary_path.prev.load(
                    std::memory_order_acquire);
                if (actual_previous != previous) {
                    snapshot.audit_records.push_back(
                        {"ERROR", "PREV_NEXT_MISMATCH", species, chromosome,
                         current->start,
                         "Segment previous pointer does not match traversal"});
                }
                const auto lookup = segment_lookup.find(current.get());
                if (lookup == segment_lookup.end()) {
                    snapshot.audit_records.push_back(
                        {"ERROR", "SEGMENT_NOT_IN_BLOCK_POOL", species,
                         chromosome, current->start,
                         "Path segment is not registered in an active Block"});
                } else {
                    const auto [block_index, segment_index] = lookup->second;
                    auto copied = snapshot.blocks[block_index].segments[segment_index];
                    copied.path_index = path.segments.size();
                    snapshot.blocks[block_index].segments[segment_index].path_index =
                        copied.path_index;
                    if (copied.species != species ||
                        copied.chromosome != chromosome) {
                        snapshot.audit_records.push_back(
                            {"ERROR", "PATH_ANCHOR_KEY_MISMATCH", species,
                             chromosome, copied.graph_start,
                             "Path location disagrees with Block anchor key"});
                    }
                    if (!path.segments.empty() &&
                        copied.graph_start < path.segments.back().graph_start) {
                        snapshot.audit_records.push_back(
                            {"ERROR", "PATH_COORDINATE_ORDER", species,
                             chromosome, copied.graph_start,
                             "Primary path Segment starts before its predecessor"});
                    }
                    path.segments.push_back(std::move(copied));
                    ++snapshot.active_segment_count;
                }
                previous = current;
                current = current->primary_path.next.load(
                    std::memory_order_acquire);
            }
            if (current && current->isTail()) {
                const auto tail_previous = current->primary_path.prev.load(
                    std::memory_order_acquire);
                if (tail_previous != previous) {
                    snapshot.audit_records.push_back(
                        {"ERROR", "TAIL_PREV_MISMATCH", species, chromosome, 0,
                         "Tail previous pointer does not match traversal"});
                }
            }
            if (!path.segments.empty()) {
                snapshot.paths.push_back(std::move(path));
            }
        }
    }
    return snapshot;
}

DetectionResult detectProblemWindows(GraphSnapshot snapshot,
                                     const Options& options) {
    const auto start = Clock::now();
    DetectionResult result;
    result.snapshot = std::move(snapshot);
    const auto blocks = blockIndex(result.snapshot);

    std::vector<std::vector<BoundaryEvidence>> path_boundaries;
    path_boundaries.resize(result.snapshot.paths.size());
    for (size_t path_index = 0; path_index < result.snapshot.paths.size();
         ++path_index) {
        const auto& path = result.snapshot.paths[path_index];
        auto& boundaries = path_boundaries[path_index];
        if (path.segments.size() < 2) {
            continue;
        }
        boundaries.reserve(path.segments.size() - 1);
        for (size_t boundary_index = 0;
             boundary_index + 1 < path.segments.size(); ++boundary_index) {
            boundaries.push_back(buildBoundary(
                path, boundary_index, blocks, result.snapshot.input_genomes,
                options));
        }
        addDropoutAndTransientSignals(path, blocks, boundaries);
        for (const auto& boundary : boundaries) {
            if (boundary.signals.empty()) {
                continue;
            }
            Seed seed;
            seed.boundary_id = boundary.boundary_id;
            seed.species = path.species;
            seed.chromosome = path.chromosome;
            seed.boundary_index = boundary.path_boundary_index;
            seed.signals = boundary.signals;
            seed.initial_priority = initialPriority(boundary, options);
            seed.reason = join(boundary.signals, ";");
            seed.seed_id = hashId(
                'D', seed.boundary_id + '|' + join(seed.signals, ";"));
            result.seeds.push_back(std::move(seed));
        }
    }

    auto raw_windows =
        buildRawWindows(result.snapshot, path_boundaries, blocks, options);
    std::map<std::string, CandidateWindow> deduplicated;
    for (const auto& raw : raw_windows) {
        auto window = materializeWindow(raw, result.snapshot, path_boundaries,
                                        blocks, options);
        const auto key = windowCanonicalKey(window);
        const auto found = deduplicated.find(key);
        if (found == deduplicated.end() ||
            std::tie(window.report_species, window.report_chromosome,
                     window.graph_start, window.graph_end) <
                std::tie(found->second.report_species,
                         found->second.report_chromosome,
                         found->second.graph_start, found->second.graph_end)) {
            deduplicated[key] = std::move(window);
        }
    }
    for (auto& [key, window] : deduplicated) {
        (void)key;
        result.windows.push_back(std::move(window));
    }
    std::sort(result.windows.begin(), result.windows.end(),
              [](const CandidateWindow& left, const CandidateWindow& right) {
                  return std::tie(left.report_species, left.report_chromosome,
                                  left.graph_start, left.graph_end,
                                  left.window_id) <
                         std::tie(right.report_species, right.report_chromosome,
                                  right.graph_start, right.graph_end,
                                  right.window_id);
              });

    std::unordered_set<std::string> bridge_boundaries;
    std::unordered_set<std::string> reported_boundaries;
    for (const auto& window : result.windows) {
        reported_boundaries.insert(window.boundary_ids.begin(),
                                   window.boundary_ids.end());
        if (window.priority_tier != "A" && window.priority_tier != "B") {
            continue;
        }
        for (const auto& seed_id : window.seed_boundary_ids) {
            bridge_boundaries.insert(seed_id);
        }
    }

    for (auto& boundaries : path_boundaries) {
        for (auto& boundary : boundaries) {
            if (!reported_boundaries.contains(boundary.boundary_id)) {
                continue;
            }
            boundary.detector_recommends_bridge =
                bridge_boundaries.contains(boundary.boundary_id) &&
                bridgeSignal(boundary) && !boundary.hard_boundary &&
                !boundary.copy_ambiguity;
            result.boundaries.push_back(std::move(boundary));
        }
    }
    std::sort(result.boundaries.begin(), result.boundaries.end(),
              [](const BoundaryEvidence& left, const BoundaryEvidence& right) {
                  return std::tie(left.reference_species,
                                  left.reference_chromosome,
                                  left.left_reference_segment.graph_start,
                                  left.right_reference_segment.graph_start,
                                  left.boundary_id) <
                         std::tie(right.reference_species,
                                  right.reference_chromosome,
                                  right.left_reference_segment.graph_start,
                                  right.right_reference_segment.graph_start,
                                  right.boundary_id);
              });
    std::sort(result.seeds.begin(), result.seeds.end(),
              [](const Seed& left, const Seed& right) {
                  return std::tie(left.species, left.chromosome,
                                  left.boundary_index, left.seed_id) <
                         std::tie(right.species, right.chromosome,
                                  right.boundary_index, right.seed_id);
              });
    result.detection_seconds =
        std::chrono::duration<double>(Clock::now() - start).count();
    return result;
}

void writeDetectionReport(const DetectionResult& result,
                          const Options& options) {
    if (options.report_dir.empty()) {
        throw std::invalid_argument(
            "window report directory cannot be empty when detection is enabled");
    }
    const auto round_dir = options.report_dir / roundName(result.snapshot.round_id);
    std::filesystem::create_directories(round_dir);

    std::unordered_map<std::string, const BoundaryEvidence*> boundaries;
    for (const auto& boundary : result.boundaries) {
        boundaries.emplace(boundary.boundary_id, &boundary);
    }
    const auto blocks = blockIndex(result.snapshot);

    writeFile(round_dir / "windows.tsv", [&](std::ostream& out) {
        out << "window_id\tdetector_version\tthreshold_profile\tround_id\t"
               "current_reference\tpriority_tier\treport_species\t"
               "report_chromosome\tgraph_start\tgraph_end\toriginal_start\t"
               "original_end\tcoordinate_status\twindow_span\tleft_anchor_id\t"
               "right_anchor_id\tleft_anchor_level\tright_anchor_level\t"
               "has_two_sided_anchors\tinternal_block_count\t"
               "internal_boundary_count\tmin_block_length\tmax_block_length\t"
               "max_reference_gap\tsignal_list\tcomplex_flag_list\t"
               "input_genome_count\tcurrently_aligned_genome_count\t"
               "max_possible_k\tcandidate_subset_count\tcandidate_genome_sets\t"
               "excluded_genomes\tcopy_ambiguity\tsoftmask_fraction\t"
               "identity_summary\tdetector_reason\n";
        for (const auto& window : result.windows) {
            const uint64_t span = window.original_coordinates_resolved
                                      ? window.original_end - window.original_start
                                      : window.graph_end - window.graph_start;
            out << window.window_id << "\tphase1-v1\t"
                << tsvSafe(options.threshold_profile) << '\t' << window.round_id
                << '\t' << tsvSafe(window.current_reference) << '\t'
                << window.priority_tier << '\t'
                << tsvSafe(window.report_species) << '\t'
                << tsvSafe(window.report_chromosome) << '\t'
                << window.graph_start << '\t' << window.graph_end << '\t'
                << window.original_start << '\t' << window.original_end << '\t'
                << (window.original_coordinates_resolved ? "resolved"
                                                         : "unresolved")
                << '\t' << span << '\t' << window.left_anchor_id << '\t'
                << window.right_anchor_id << '\t' << window.left_anchor_level
                << '\t' << window.right_anchor_level << '\t'
                << boolText(window.has_two_sided_anchors) << '\t'
                << window.block_ids.size() << '\t' << window.boundary_ids.size()
                << '\t' << window.min_block_length << '\t'
                << window.max_block_length << '\t' << window.max_reference_gap
                << '\t' << join(window.signals, ";") << '\t'
                << join(window.complex_flags, ";") << '\t'
                << window.input_genome_count << '\t'
                << window.currently_aligned_genome_count << '\t'
                << window.max_possible_k << '\t'
                << (window.max_k_genomes.empty() ? 0 : 1) << '\t'
                << join(window.max_k_genomes, ",") << '\t'
                << join(window.excluded_genomes, ",") << '\t'
                << boolText(contains(window.complex_flags, "COPY_AMBIGUITY"))
                << "\tunknown\tunknown\t"
                << tsvSafe(window.detector_reason) << '\n';
        }
    });

    writeFile(round_dir / "window_boundaries.tsv", [&](std::ostream& out) {
        out << "window_id\tboundary_id\tleft_block_id\tright_block_id\t"
               "signal_type\tdetector_priority\tdetector_recommends_bridge\t"
               "reference_species\treference_chromosome\tleft_ref_graph_start\t"
               "left_ref_graph_end\tright_ref_graph_start\tright_ref_graph_end\t"
               "left_ref_start\tleft_ref_end\tright_ref_start\tright_ref_end\t"
               "coordinate_status\treference_gap\treference_overlap\t"
               "min_adjacent_block_length\tmax_adjacent_block_length\t"
               "support_genome_count\tdropout_genomes\thard_boundary\t"
               "copy_ambiguity\n";
        for (const auto& window : result.windows) {
            for (const auto& boundary_id : window.boundary_ids) {
                const auto& boundary = *boundaries.at(boundary_id);
                const bool resolved =
                    boundary.left_reference_segment.original_coordinates_resolved &&
                    boundary.right_reference_segment.original_coordinates_resolved;
                size_t support = 1;
                for (const auto& genome : boundary.genomes) {
                    if (genome.unique_target) {
                        ++support;
                    }
                }
                out << window.window_id << '\t' << boundary.boundary_id << '\t'
                    << boundary.left_block_id << '\t'
                    << boundary.right_block_id << '\t'
                    << join(boundary.signals, ";") << '\t'
                    << window.priority_tier << '\t'
                    << boolText(boundary.detector_recommends_bridge) << '\t'
                    << tsvSafe(boundary.reference_species) << '\t'
                    << tsvSafe(boundary.reference_chromosome) << '\t'
                    << boundary.left_reference_segment.graph_start << '\t'
                    << boundary.left_reference_segment.graph_end << '\t'
                    << boundary.right_reference_segment.graph_start << '\t'
                    << boundary.right_reference_segment.graph_end << '\t'
                    << boundary.left_reference_segment.original_start << '\t'
                    << boundary.left_reference_segment.original_end << '\t'
                    << boundary.right_reference_segment.original_start << '\t'
                    << boundary.right_reference_segment.original_end << '\t'
                    << (resolved ? "resolved" : "unresolved") << '\t'
                    << boundary.reference_gap << '\t'
                    << boundary.reference_overlap << '\t'
                    << boundary.min_adjacent_block_length << '\t'
                    << boundary.max_adjacent_block_length << '\t' << support
                    << '\t' << join(boundary.dropout_genomes, ",") << '\t'
                    << boolText(boundary.hard_boundary) << '\t'
                    << boolText(boundary.copy_ambiguity) << '\n';
            }
        }
    });

    writeFile(round_dir / "boundary_genomes.tsv", [&](std::ostream& out) {
        out << "window_id\tboundary_id\treference_species\treference_chromosome\t"
               "target_species\tleft_present\tright_present\tunique_target\t"
               "left_target_chromosome\tleft_target_graph_start\t"
               "left_target_graph_end\tleft_target_start\tleft_target_end\t"
               "left_strand\tright_target_chromosome\tright_target_graph_start\t"
               "right_target_graph_end\tright_target_start\tright_target_end\t"
               "right_strand\tcoordinate_status\ttarget_gap\t"
               "target_gap_resolved\torder_consistent\tstrand_consistent\t"
               "target_consistent\n";
        for (const auto& window : result.windows) {
            for (const auto& boundary_id : window.boundary_ids) {
                const auto& boundary = *boundaries.at(boundary_id);
                for (const auto& genome : boundary.genomes) {
                    const bool resolved = genome.left_original_resolved &&
                                          genome.right_original_resolved;
                    out << window.window_id << '\t' << boundary.boundary_id << '\t'
                        << tsvSafe(boundary.reference_species) << '\t'
                        << tsvSafe(boundary.reference_chromosome) << '\t'
                        << tsvSafe(genome.target_species) << '\t'
                        << boolText(genome.left_present) << '\t'
                        << boolText(genome.right_present) << '\t'
                        << boolText(genome.unique_target) << '\t'
                        << tsvSafe(genome.left_chromosome) << '\t'
                        << genome.left_graph_start << '\t'
                        << genome.left_graph_end << '\t'
                        << genome.left_original_start << '\t'
                        << genome.left_original_end << '\t'
                        << strandText(genome.left_reverse) << '\t'
                        << tsvSafe(genome.right_chromosome) << '\t'
                        << genome.right_graph_start << '\t'
                        << genome.right_graph_end << '\t'
                        << genome.right_original_start << '\t'
                        << genome.right_original_end << '\t'
                        << strandText(genome.right_reverse) << '\t'
                        << (resolved ? "resolved" : "unresolved") << '\t'
                        << genome.target_gap << '\t'
                        << boolText(genome.target_gap_resolved) << '\t'
                        << boolText(genome.order_consistent) << '\t'
                        << boolText(genome.strand_consistent) << '\t'
                        << boolText(genome.target_consistent) << '\n';
                }
            }
        }
    });

    writeFile(round_dir / "window_blocks.tsv", [&](std::ostream& out) {
        out << "window_id\tblock_id\tcanonical_block_signature\twindow_role\t"
               "min_segment_length\tmax_segment_length\tparticipating_genomes\t"
               "participating_genome_count\tsegment_count\tprimary_or_secondary\t"
               "local_identity\tgap_rate\tmismatch_rate\tsoftmask_fraction\t"
               "suspected_wrong\tsuspected_wrong_reason\n";
        for (const auto& window : result.windows) {
            for (const auto& block_id : window.block_ids) {
                const auto& block = *blocks.at(block_id);
                std::string role = "INTERIOR";
                if (block_id == window.left_anchor_id) {
                    role = "LEFT_ANCHOR";
                }
                if (block_id == window.right_anchor_id) {
                    role = role == "LEFT_ANCHOR" ? "BOTH_ANCHORS"
                                                   : "RIGHT_ANCHOR";
                }
                out << window.window_id << '\t' << block.block_id << '\t'
                    << tsvSafe(block.canonical_signature) << '\t' << role << '\t'
                    << block.min_segment_length << '\t'
                    << block.max_segment_length << '\t'
                    << join(block.participating_genomes, ",") << '\t'
                    << block.participating_genomes.size() << '\t'
                    << block.segments.size() << '\t'
                    << (block.has_secondary ? "mixed" : "primary")
                    << "\tunknown\tunknown\tunknown\tunknown\t0\t\n";
            }
        }
    });

    writeFile(round_dir / "window_block_segments.tsv", [&](std::ostream& out) {
        out << "window_id\tblock_id\tsegment_id\tspecies\tchromosome\t"
               "graph_start\tgraph_end\toriginal_start\toriginal_end\t"
               "coordinate_status\tstrand\talign_role\tpath_index\tcigar\n";
        for (const auto& window : result.windows) {
            for (const auto& block_id : window.block_ids) {
                const auto& block = *blocks.at(block_id);
                for (const auto& segment : block.segments) {
                    out << window.window_id << '\t' << block.block_id << '\t'
                        << segment.segment_id << '\t' << tsvSafe(segment.species)
                        << '\t' << tsvSafe(segment.chromosome) << '\t'
                        << segment.graph_start << '\t' << segment.graph_end << '\t'
                        << segment.original_start << '\t' << segment.original_end
                        << '\t'
                        << (segment.original_coordinates_resolved ? "resolved"
                                                                  : "unresolved")
                        << '\t' << strandText(segment.reverse) << '\t'
                        << (segment.primary ? "PRIMARY" : "SECONDARY") << '\t'
                        << segment.path_index << '\t'
                        << tsvSafe(segment.cigar_summary) << '\n';
                }
            }
        }
    });

    writeFile(round_dir / "window_genomes.tsv", [&](std::ostream& out) {
        out << "window_id\tspecies\tchromosome\tgraph_start\tgraph_end\t"
               "original_start\toriginal_end\tcoordinate_status\tstrand\t"
               "left_anchor_present\tright_anchor_present\tcurrently_aligned\t"
               "included_in_max_k\texcluded_reason\traw_sequence_available\t"
               "softmask_fraction\n";
        for (const auto& window : result.windows) {
            for (const auto& genome : window.genomes) {
                out << window.window_id << '\t' << tsvSafe(genome.species) << '\t'
                    << tsvSafe(genome.chromosome) << '\t' << genome.graph_start
                    << '\t' << genome.graph_end << '\t' << genome.original_start
                    << '\t' << genome.original_end << '\t'
                    << (genome.original_coordinates_resolved ? "resolved"
                                                             : "unresolved")
                    << '\t' << genome.strand << '\t'
                    << boolText(genome.left_anchor_present) << '\t'
                    << boolText(genome.right_anchor_present) << '\t'
                    << boolText(genome.currently_aligned) << '\t'
                    << boolText(genome.included_in_max_k) << '\t'
                    << genome.excluded_reason << "\t1\tunknown\n";
            }
        }
    });

    writeFile(round_dir / "seeds.tsv", [&](std::ostream& out) {
        out << "seed_id\tboundary_id\tspecies\tchromosome\tboundary_index\t"
               "signals\tinitial_priority\treason\n";
        for (const auto& seed : result.seeds) {
            out << seed.seed_id << '\t' << seed.boundary_id << '\t'
                << tsvSafe(seed.species) << '\t' << tsvSafe(seed.chromosome)
                << '\t' << seed.boundary_index << '\t'
                << join(seed.signals, ";") << '\t' << seed.initial_priority
                << '\t' << tsvSafe(seed.reason) << '\n';
        }
    });

    writeFile(round_dir / "graph_audit.tsv", [&](std::ostream& out) {
        out << "severity\tcode\tspecies\tchromosome\tgraph_position\tmessage\n";
        for (const auto& record : result.snapshot.audit_records) {
            out << record.severity << '\t' << record.code << '\t'
                << tsvSafe(record.species) << '\t'
                << tsvSafe(record.chromosome) << '\t' << record.graph_position
                << '\t' << tsvSafe(record.message) << '\n';
        }
    });

    const auto counts = tierCounts(result.windows);
    auto writeSummary = [&](std::ostream& out) {
        out << "{\n"
            << "  \"detector_version\": \"phase1-v1\",\n"
            << "  \"round_id\": " << result.snapshot.round_id << ",\n"
            << "  \"current_reference\": \""
            << jsonEscape(result.snapshot.current_reference) << "\",\n"
            << "  \"active_block_count\": " << result.snapshot.blocks.size()
            << ",\n"
            << "  \"active_segment_count\": "
            << result.snapshot.active_segment_count << ",\n"
            << "  \"path_count\": " << result.snapshot.paths.size() << ",\n"
            << "  \"audit_record_count\": "
            << result.snapshot.audit_records.size() << ",\n"
            << "  \"seed_count\": " << result.seeds.size() << ",\n"
            << "  \"window_count\": " << result.windows.size() << ",\n"
            << "  \"reported_boundary_count\": " << result.boundaries.size()
            << ",\n"
            << "  \"tier_counts\": {";
        bool first = true;
        for (const auto& [tier, count] : counts) {
            if (!first) out << ',';
            out << "\n    \"" << jsonEscape(tier) << "\": " << count;
            first = false;
        }
        if (!counts.empty()) out << '\n';
        out << "  },\n"
            << "  \"snapshot_seconds\": " << result.snapshot_seconds << ",\n"
            << "  \"detection_seconds\": " << result.detection_seconds << '\n'
            << "}\n";
    };
    writeFile(round_dir / "summary.json", writeSummary);
    writeFile(options.report_dir / "summary.json", writeSummary);

    writeFile(options.report_dir / "run_manifest.json", [&](std::ostream& out) {
        out << "{\n"
            << "  \"status\": \"complete\",\n"
            << "  \"detector_version\": \"phase1-v1\",\n"
            << "  \"threshold_profile\": \""
            << jsonEscape(options.threshold_profile) << "\",\n"
            << "  \"mode\": \""
            << detectionModeToString(options.mode) << "\",\n"
            << "  \"latest_round\": " << result.snapshot.round_id << ",\n"
            << "  \"current_reference\": \""
            << jsonEscape(result.snapshot.current_reference) << "\",\n"
            << "  \"coordinate_system\": {\"graph\": \"masked-local\", "
               "\"original\": \"0-based-half-open-forward\"},\n"
            << "  \"parameters\": {\n"
            << "    \"micro_block_max_bp\": "
            << options.micro_block_max_bp << ",\n"
            << "    \"short_block_max_bp\": "
            << options.short_block_max_bp << ",\n"
            << "    \"weak_block_upper_bp\": "
            << options.weak_block_upper_bp << ",\n"
            << "    \"primary_gap_max_bp\": "
            << options.primary_gap_max_bp << ",\n"
            << "    \"extended_gap_max_bp\": "
            << options.extended_gap_max_bp << ",\n"
            << "    \"hard_boundary_gap_bp\": "
            << options.hard_boundary_gap_bp << ",\n"
            << "    \"anchor_min_segment_bp\": "
            << options.anchor_min_segment_bp << ",\n"
            << "    \"strong_anchor_bp\": " << options.strong_anchor_bp
            << ",\n"
            << "    \"max_window_span_bp\": "
            << options.max_window_span_bp << ",\n"
            << "    \"subset_search_budget\": "
            << options.subset_search_budget << "\n"
            << "  },\n"
            << "  \"input_genomes\": [";
        for (size_t i = 0; i < result.snapshot.input_genomes.size(); ++i) {
            if (i != 0) out << ", ";
            out << '"' << jsonEscape(result.snapshot.input_genomes[i]) << '"';
        }
        out << "],\n"
            << "  \"round_directory\": \""
            << jsonEscape(roundName(result.snapshot.round_id)) << "\",\n"
            << "  \"files\": [\"windows.tsv\", "
               "\"window_boundaries.tsv\", \"boundary_genomes.tsv\", "
               "\"window_blocks.tsv\", \"window_block_segments.tsv\", "
               "\"window_genomes.tsv\", \"seeds.tsv\", "
               "\"graph_audit.tsv\", \"summary.json\"]\n"
            << "}\n";
    });
}

DetectionResult detectAndWriteProblemWindows(
    const RaMeshMultiGenomeGraph& graph,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    uint64_t round_id,
    const SpeciesName& current_reference,
    const Options& options) {
    const auto snapshot_start = Clock::now();
    auto snapshot =
        snapshotGraph(graph, managers, round_id, current_reference);
    const double snapshot_seconds =
        std::chrono::duration<double>(Clock::now() - snapshot_start).count();
    auto result = detectProblemWindows(std::move(snapshot), options);
    result.snapshot_seconds = snapshot_seconds;
    writeDetectionReport(result, options);
    spdlog::info(
        "[window-detection] round={} reference={} blocks={} segments={} "
        "seeds={} windows={} report={}",
        round_id, current_reference, result.snapshot.blocks.size(),
        result.snapshot.active_segment_count, result.seeds.size(),
        result.windows.size(),
        (options.report_dir / roundName(round_id)).string());
    return result;
}

}  // namespace RaMesh::WindowDetection
