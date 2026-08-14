#include "structural_break_repair.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <future>
#include <limits>
#include <mutex>
#include <optional>
#include <set>
#include <shared_mutex>
#include <sstream>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "align.h"
#include "ramesh.h"
#include "spdlog/spdlog.h"
#include "threadpool.h"

namespace RaMesh::StructuralBreakRepair {
namespace {

using Clock = std::chrono::steady_clock;

enum Signal : uint8_t {
    TARGET_SWITCH = 1U,
    STRAND_SWITCH = 2U,
    ORDER_BREAK = 4U
};

struct PairQuality {
    uint64_t common = 0;
    uint64_t matches = 0;
    double coverage = 0.0;
    double identity = 0.0;
    int64_t score = std::numeric_limits<int64_t>::min();
};

struct CandidatePath {
    SpeciesChrPair key;
    SegPtr left;
    SegPtr right;
    uint_t start = 0;
    uint_t length = 0;
    Strand strand = Strand::FORWARD;
    bool anomalous = false;
};

struct Candidate {
    SpeciesChrPair reference_key;
    SegPtr left_reference;
    SegPtr right_reference;
    std::vector<CandidatePath> paths;
    size_t anomaly_path_index = 0;
    std::vector<SegPtr> reference_interior;
    std::vector<std::pair<SpeciesChrPair, SegPtr>> anomalous_segments;
    std::vector<BlockPtr> interior_blocks;
    uint_t reference_start = 0;
    uint_t reference_length = 0;
    uint8_t signals = 0;
    uint8_t strong_anchors = 0;
    uint64_t serial = 0;
    std::string signature;
};

struct PathEdit {
    SpeciesChrPair key;
    SegPtr previous;
    SegPtr next;
    std::vector<SegPtr> old_segments;
    SegPtr replacement;
};

struct Prepared {
    Candidate candidate;
    BlockPtr replacement_block;
    std::vector<SegPtr> replacement_segments;
    std::vector<PathEdit> edits;
    std::vector<Cigar_t> cigars;
    std::vector<PairQuality> new_qualities;
    std::vector<double> anchor_identity_floors;
    int64_t old_score = std::numeric_limits<int64_t>::min();
    bool msa_ok = false;
    bool quality_ok = false;
    std::string rejection;
    double sequence_seconds = 0.0;
    double msa_seconds = 0.0;
    double quality_seconds = 0.0;
};

double secondsSince(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

uint64_t segmentEnd(const SegPtr& segment) {
    return static_cast<uint64_t>(segment->start) + segment->length;
}

bool fetchSequence(
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    const SpeciesChrPair& key,
    uint_t start,
    uint_t length,
    Strand strand,
    std::string& sequence) {
    const auto found = managers.find(key.first);
    if (found == managers.end() || !found->second || length == 0) return false;
    try {
        sequence = std::visit(
            [&](auto& manager) -> std::string {
                using Manager = std::decay_t<decltype(manager)>;
                if (!manager) return {};
                if constexpr (std::is_same_v<
                                  Manager,
                                  std::unique_ptr<SeqPro::SequenceManager>>) {
                    return manager->getSubSequence(key.second, start, length);
                } else {
                    return manager->getOriginalManager().getSubSequence(
                        key.second, start, length);
                }
            },
            *found->second);
    } catch (const std::exception&) {
        return false;
    }
    if (sequence.size() != length) return false;
    if (strand == Strand::REVERSE) reverseComplement(sequence);
    return true;
}

bool cigarFromRows(const std::string& reference,
                   const std::string& query,
                   Cigar_t& cigar) {
    cigar.clear();
    if (reference.empty() || reference.size() != query.size()) return false;
    for (size_t column = 0; column < reference.size(); ++column) {
        const bool reference_gap = reference[column] == '-';
        const bool query_gap = query[column] == '-';
        if (reference_gap && query_gap) continue;
        appendCigarOp(cigar,
                      reference_gap ? 'I' : (query_gap ? 'D' : 'M'), 1);
    }
    return !cigar.empty();
}

PairQuality evaluateRows(const std::string& reference,
                         const std::string& query,
                         size_t reference_length,
                         size_t query_length) {
    PairQuality quality;
    if (reference.size() != query.size() || reference.empty()) return quality;
    bool in_reference_gap = false;
    bool in_query_gap = false;
    int64_t score = 0;
    for (size_t column = 0; column < reference.size(); ++column) {
        const char r = reference[column];
        const char q = query[column];
        if (r == '-' && q == '-') continue;
        if (r == '-') {
            score -= in_reference_gap ? GAP_EXTEND_PENALTY : GAP_OPEN_PENALTY;
            in_reference_gap = true;
            in_query_gap = false;
        } else if (q == '-') {
            score -= in_query_gap ? GAP_EXTEND_PENALTY : GAP_OPEN_PENALTY;
            in_query_gap = true;
            in_reference_gap = false;
        } else {
            ++quality.common;
            if (std::toupper(static_cast<unsigned char>(r)) ==
                std::toupper(static_cast<unsigned char>(q))) {
                ++quality.matches;
            }
            score += subsScore(r, q);
            in_reference_gap = false;
            in_query_gap = false;
        }
    }
    const size_t shorter = std::min(reference_length, query_length);
    quality.coverage = shorter == 0
                           ? 0.0
                           : static_cast<double>(quality.common) / shorter;
    quality.identity = quality.common == 0
                           ? 0.0
                           : static_cast<double>(quality.matches) /
                                 quality.common;
    quality.score = score;
    return quality;
}

std::optional<PairQuality> qualityFromSegment(
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    const SpeciesChrPair& reference_key,
    const SegPtr& reference,
    const SpeciesChrPair& query_key,
    const SegPtr& query) {
    if (!reference || !query) return std::nullopt;
    std::string reference_sequence;
    std::string query_sequence;
    if (!fetchSequence(managers, reference_key, reference->start,
                       reference->length, Strand::FORWARD,
                       reference_sequence) ||
        !fetchSequence(managers, query_key, query->start, query->length,
                       query->strand, query_sequence)) {
        return std::nullopt;
    }
    Cigar_t cigar = query->cigar;
    if (cigar.empty()) {
        if (reference->length != query->length) return std::nullopt;
        appendCigarOp(cigar, 'M', reference->length);
    }
    if (countRefLength(cigar) != reference->length ||
        countQryLength(cigar) != query->length) {
        return std::nullopt;
    }
    const auto rows = buildAlignment(reference_sequence, query_sequence, cigar);
    return evaluateRows(rows.first, rows.second, reference_sequence.size(),
                        query_sequence.size());
}

using OuterAnchors =
    std::map<SpeciesName, std::pair<SpeciesChrPair, SegPtr>>;

bool uniqueOuterBlock(const BlockPtr& block,
                      const SpeciesName& reference_species,
                      const ChrName& reference_chromosome,
                      OuterAnchors& anchors,
                      SegPtr& reference) {
    if (!block || block->anchors.size() < 2 ||
        block->ref_chr != reference_chromosome) {
        return false;
    }
    anchors.clear();
    reference.reset();
    for (const auto& [key, segment] : block->anchors) {
        if (!segment || !segment->isSegment() || !segment->isPrimary() ||
            segment->parent_block.get() != block.get()) {
            return false;
        }
        if (!anchors.emplace(key.first, std::pair{key, segment}).second) {
            return false;
        }
        if (key.first == reference_species && key.second == reference_chromosome) {
            if (reference || segment->strand != Strand::FORWARD) return false;
            reference = segment;
        } else if (key.first == reference_species) {
            return false;
        }
    }
    return reference && anchors.size() >= 2;
}

bool collectBetween(const SegPtr& previous,
                    const SegPtr& next,
                    size_t maximum,
                    std::vector<SegPtr>& segments) {
    segments.clear();
    if (!previous || !next || previous == next) return false;
    auto current = previous->primary_path.next.load(std::memory_order_acquire);
    while (current && current != next && !current->isTail()) {
        if (!current->isSegment() || segments.size() >= maximum) return false;
        segments.push_back(current);
        current = current->primary_path.next.load(std::memory_order_acquire);
    }
    return current == next;
}

std::string candidateSignature(const Candidate& candidate) {
    std::ostringstream out;
    out << candidate.reference_key.first << '|' << candidate.reference_key.second
        << '|' << candidate.reference_start << '|' << candidate.reference_length
        << '|';
    for (const auto& path : candidate.paths) {
        out << path.key.first << ':' << path.key.second << ':' << path.start
            << ':' << path.length << ':' << static_cast<int>(path.strand)
            << ':' << path.anomalous << ';';
    }
    for (const auto& block : candidate.interior_blocks) out << block.get() << ',';
    return out.str();
}

std::vector<Candidate> scanCandidates(
    RaMeshMultiGenomeGraph& graph,
    const SpeciesName& reference_species,
    const Options& options,
    Result& result) {
    std::vector<Candidate> candidates;
    const auto reference_graph = graph.species_graphs.find(reference_species);
    if (reference_graph == graph.species_graphs.end()) return candidates;
    uint64_t serial = 0;
    for (const auto& [chromosome, end] : reference_graph->second.chr2end) {
        std::vector<SegPtr> path;
        for (auto current = end.head->primary_path.next.load();
             current && !current->isTail();
             current = current->primary_path.next.load()) {
            if (current->isSegment()) path.push_back(current);
        }
        for (size_t left_index = 0; left_index < path.size(); ++left_index) {
            const size_t maximum_right = std::min(
                path.size(), left_index + options.maximum_interior_blocks + 2);
            for (size_t right_index = left_index + 2;
                 right_index < maximum_right; ++right_index) {
                ++result.scanned_windows;
                const auto& left_reference = path[left_index];
                const auto& right_reference = path[right_index];
                OuterAnchors left_anchors;
                OuterAnchors right_anchors;
                SegPtr checked_left_reference;
                SegPtr checked_right_reference;
                if (!uniqueOuterBlock(left_reference->parent_block,
                                      reference_species, chromosome,
                                      left_anchors, checked_left_reference) ||
                    !uniqueOuterBlock(right_reference->parent_block,
                                      reference_species, chromosome,
                                      right_anchors, checked_right_reference) ||
                    checked_left_reference != left_reference ||
                    checked_right_reference != right_reference ||
                    left_anchors.size() != right_anchors.size() ||
                    left_reference->length < options.minimum_outer_anchor ||
                    right_reference->length < options.minimum_outer_anchor) {
                    ++result.outer_anchor_invalid;
                    continue;
                }
                bool same_participants = true;
                for (const auto& [species, left_anchor] : left_anchors) {
                    const auto right_anchor = right_anchors.find(species);
                    if (right_anchor == right_anchors.end() ||
                        left_anchor.first != right_anchor->second.first) {
                        same_participants = false;
                        break;
                    }
                }
                if (!same_participants) {
                    ++result.outer_anchor_invalid;
                    continue;
                }
                const uint64_t reference_start = segmentEnd(left_reference);
                if (right_reference->start <= reference_start) {
                    ++result.reference_empty;
                    continue;
                }
                const uint64_t reference_length =
                    right_reference->start - reference_start;
                if (reference_length == 0) {
                    ++result.reference_empty;
                    continue;
                }
                if (reference_length > options.maximum_span) {
                    ++result.span_exceeded;
                    continue;
                }

                Candidate candidate;
                candidate.reference_key = {reference_species, chromosome};
                candidate.left_reference = left_reference;
                candidate.right_reference = right_reference;
                candidate.reference_start = static_cast<uint_t>(reference_start);
                candidate.reference_length = static_cast<uint_t>(reference_length);
                candidate.serial = serial++;
                candidate.strong_anchors = static_cast<uint8_t>(
                    left_reference->length >= options.strong_outer_anchor) +
                    static_cast<uint8_t>(
                        right_reference->length >= options.strong_outer_anchor);

                bool invalid = false;
                for (const auto& [species, left_anchor] : left_anchors) {
                    if (species == reference_species) continue;
                    const auto& right_anchor = right_anchors.at(species);
                    const auto& left_segment = left_anchor.second;
                    const auto& right_segment = right_anchor.second;
                    if (left_segment->strand != right_segment->strand ||
                        left_segment->length < options.minimum_outer_anchor ||
                        right_segment->length < options.minimum_outer_anchor) {
                        invalid = true;
                        break;
                    }
                    uint64_t interval_start = 0;
                    uint64_t interval_end = 0;
                    SegPtr physical_left;
                    SegPtr physical_right;
                    if (left_segment->strand == Strand::FORWARD) {
                        interval_start = segmentEnd(left_segment);
                        interval_end = right_segment->start;
                        physical_left = left_segment;
                        physical_right = right_segment;
                    } else {
                        interval_start = segmentEnd(right_segment);
                        interval_end = left_segment->start;
                        physical_left = right_segment;
                        physical_right = left_segment;
                    }
                    if (interval_end <= interval_start) {
                        invalid = true;
                        break;
                    }
                    if (interval_end - interval_start > options.maximum_span) {
                        ++result.span_exceeded;
                        invalid = true;
                        break;
                    }
                    std::vector<SegPtr> existing;
                    if (!collectBetween(physical_left, physical_right,
                                        options.maximum_interior_blocks + 1,
                                        existing) ||
                        !existing.empty()) {
                        ++result.too_many_interior_blocks;
                        invalid = true;
                        break;
                    }
                    candidate.paths.push_back(
                        {left_anchor.first, left_segment, right_segment,
                         static_cast<uint_t>(interval_start),
                         static_cast<uint_t>(interval_end - interval_start),
                         left_segment->strand, false});
                    candidate.strong_anchors += static_cast<uint8_t>(
                        left_segment->length >= options.strong_outer_anchor) +
                        static_cast<uint8_t>(right_segment->length >=
                                             options.strong_outer_anchor);
                }
                if (invalid || candidate.paths.empty()) {
                    ++result.outer_anchor_invalid;
                    continue;
                }

                std::optional<size_t> anomaly_path;
                for (size_t index = left_index + 1; index < right_index; ++index) {
                    const auto& reference_segment = path[index];
                    const auto& block = reference_segment->parent_block;
                    if (!block || block->anchors.size() != 2 ||
                        block->ref_chr != chromosome) {
                        invalid = true;
                        break;
                    }
                    candidate.reference_interior.push_back(reference_segment);
                    candidate.interior_blocks.push_back(block);
                    bool found_anomaly = false;
                    for (const auto& [key, segment] : block->anchors) {
                        if (key == candidate.reference_key) continue;
                        const auto path_it = std::find_if(
                            candidate.paths.begin(), candidate.paths.end(),
                            [&](const CandidatePath& candidate_path) {
                                return candidate_path.key.first == key.first;
                            });
                        if (!segment || path_it == candidate.paths.end() ||
                            !segment->isPrimary() ||
                            segment->parent_block.get() != block.get()) {
                            invalid = true;
                            break;
                        }
                        const size_t path_index = static_cast<size_t>(
                            std::distance(candidate.paths.begin(), path_it));
                        uint8_t signals = 0;
                        if (key.second != path_it->key.second) {
                            signals |= TARGET_SWITCH;
                        }
                        if (segment->strand != path_it->strand) {
                            signals |= STRAND_SWITCH;
                        }
                        if (key.second == path_it->key.second &&
                            (segment->start < path_it->start ||
                             segmentEnd(segment) >
                                 static_cast<uint64_t>(path_it->start) +
                                     path_it->length)) {
                            signals |= ORDER_BREAK;
                        }
                        if (signals == 0) {
                            invalid = true;
                            break;
                        }
                        if (anomaly_path && *anomaly_path != path_index) {
                            invalid = true;
                            break;
                        }
                        anomaly_path = path_index;
                        found_anomaly = true;
                        candidate.signals |= signals;
                        candidate.anomalous_segments.emplace_back(key, segment);
                    }
                    if (invalid || !found_anomaly) {
                        invalid = true;
                        break;
                    }
                }
                if (invalid || candidate.signals == 0 ||
                    candidate.anomalous_segments.empty() || !anomaly_path) {
                    ++result.large_gap_only;
                    continue;
                }
                candidate.anomaly_path_index = *anomaly_path;
                candidate.paths[*anomaly_path].anomalous = true;
                candidate.signature = candidateSignature(candidate);
                candidates.push_back(std::move(candidate));
            }
        }
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate& left, const Candidate& right) {
                  return std::tuple{-left.strong_anchors,
                                    left.reference_length,
                                    left.interior_blocks.size(),
                                    left.reference_key.second,
                                    left.reference_start,
                                    left.serial} <
                         std::tuple{-right.strong_anchors,
                                    right.reference_length,
                                    right.interior_blocks.size(),
                                    right.reference_key.second,
                                    right.reference_start,
                                    right.serial};
              });
    return candidates;
}

std::vector<Candidate> selectConflictFree(std::vector<Candidate> candidates,
                                          Result& result) {
    std::unordered_set<const Block*> reserved;
    std::vector<Candidate> selected;
    for (auto& candidate : candidates) {
        bool conflict = reserved.count(candidate.left_reference->parent_block.get()) ||
                        reserved.count(candidate.right_reference->parent_block.get());
        for (const auto& block : candidate.interior_blocks) {
            conflict = conflict || reserved.count(block.get()) != 0;
        }
        if (conflict) {
            ++result.conflict_deferred;
            continue;
        }
        reserved.insert(candidate.left_reference->parent_block.get());
        reserved.insert(candidate.right_reference->parent_block.get());
        for (const auto& block : candidate.interior_blocks) {
            reserved.insert(block.get());
        }
        selected.push_back(std::move(candidate));
    }
    return selected;
}

void detachPrepared(Prepared& prepared) {
    for (auto& segment : prepared.replacement_segments) {
        if (!segment) continue;
        segment->parent_block.reset();
        segment->primary_path.prev.store(nullptr);
        segment->primary_path.next.store(nullptr);
    }
    if (prepared.replacement_block) prepared.replacement_block->anchors.clear();
}

bool buildRemovalEdits(const Candidate& candidate,
                       const std::vector<SegPtr>& replacements,
                       std::vector<PathEdit>& edits) {
    if (replacements.size() != candidate.paths.size() + 1) return false;
    edits.clear();
    PathEdit reference_edit;
    reference_edit.key = candidate.reference_key;
    reference_edit.previous = candidate.left_reference;
    reference_edit.next = candidate.right_reference;
    reference_edit.old_segments = candidate.reference_interior;
    reference_edit.replacement = replacements.front();
    edits.push_back(std::move(reference_edit));

    for (size_t index = 0; index < candidate.paths.size(); ++index) {
        const auto& path = candidate.paths[index];
        PathEdit insertion;
        insertion.key = path.key;
        if (path.strand == Strand::FORWARD) {
            insertion.previous = path.left;
            insertion.next = path.right;
        } else {
            insertion.previous = path.right;
            insertion.next = path.left;
        }
        insertion.replacement = replacements[index + 1];
        edits.push_back(std::move(insertion));
    }
    return true;
}

Prepared prepareCandidate(
    Candidate candidate,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    const Options& options) {
    Prepared prepared;
    prepared.candidate = std::move(candidate);
    const auto sequence_start = Clock::now();
    std::string reference_sequence;
    if (!fetchSequence(managers, prepared.candidate.reference_key,
                       prepared.candidate.reference_start,
                       prepared.candidate.reference_length, Strand::FORWARD,
                       reference_sequence)) {
        prepared.rejection = "sequence_fetch_failed";
        prepared.sequence_seconds = secondsSince(sequence_start);
        return prepared;
    }

    std::vector<std::string> raw_sequences;
    raw_sequences.reserve(prepared.candidate.paths.size());
    prepared.anchor_identity_floors.reserve(prepared.candidate.paths.size());
    for (const auto& path : prepared.candidate.paths) {
        std::string sequence;
        if (!fetchSequence(managers, path.key, path.start, path.length,
                           path.strand, sequence)) {
            prepared.rejection = "sequence_fetch_failed";
            prepared.sequence_seconds = secondsSince(sequence_start);
            return prepared;
        }
        const auto left_quality = qualityFromSegment(
            managers, prepared.candidate.reference_key,
            prepared.candidate.left_reference, path.key, path.left);
        const auto right_quality = qualityFromSegment(
            managers, prepared.candidate.reference_key,
            prepared.candidate.right_reference, path.key, path.right);
        if (!left_quality || !right_quality) {
            prepared.rejection = "outer_anchor_quality_unavailable";
            prepared.sequence_seconds = secondsSince(sequence_start);
            return prepared;
        }
        prepared.anchor_identity_floors.push_back(
            std::min(left_quality->identity, right_quality->identity));
        raw_sequences.push_back(std::move(sequence));
    }

    int64_t old_score = 0;
    bool have_old_score = false;
    for (size_t i = 0; i < prepared.candidate.reference_interior.size(); ++i) {
        const auto& reference = prepared.candidate.reference_interior[i];
        const auto& block = prepared.candidate.interior_blocks[i];
        for (const auto& [key, segment] : block->anchors) {
            if (key.first != prepared.candidate.paths[
                    prepared.candidate.anomaly_path_index].key.first) continue;
            const auto quality = qualityFromSegment(
                managers, prepared.candidate.reference_key, reference, key,
                segment);
            if (!quality) {
                prepared.rejection = "old_projection_invalid";
                prepared.sequence_seconds = secondsSince(sequence_start);
                return prepared;
            }
            old_score += quality->score;
            have_old_score = true;
        }
    }
    if (!have_old_score) {
        prepared.rejection = "old_projection_missing";
        prepared.sequence_seconds = secondsSince(sequence_start);
        return prepared;
    }
    prepared.old_score = old_score;
    prepared.sequence_seconds = secondsSince(sequence_start);

    std::unordered_map<ChrName, std::string> rows;
    rows.emplace("0_reference", reference_sequence);
    std::vector<ChrName> row_ids;
    row_ids.reserve(prepared.candidate.paths.size());
    for (size_t index = 0; index < prepared.candidate.paths.size(); ++index) {
        const ChrName id = "q" + std::to_string(index);
        row_ids.push_back(id);
        rows.emplace(id, raw_sequences[index]);
    }
    const auto msa_start = Clock::now();
    prepared.msa_ok = alignSequencesWithExternalMsa(
        options.msa_executable, rows);
    prepared.msa_seconds = secondsSince(msa_start);
    if (!prepared.msa_ok) {
        prepared.rejection = "msa_failed";
        return prepared;
    }
    const auto reference_row = rows.find("0_reference");
    if (reference_row == rows.end()) {
        prepared.rejection = "msa_invalid";
        return prepared;
    }
    prepared.cigars.resize(prepared.candidate.paths.size());
    for (size_t index = 0; index < row_ids.size(); ++index) {
        const auto query_row = rows.find(row_ids[index]);
        if (query_row == rows.end() ||
            reference_row->second.size() != query_row->second.size() ||
            !cigarFromRows(reference_row->second, query_row->second,
                           prepared.cigars[index]) ||
            countRefLength(prepared.cigars[index]) !=
                prepared.candidate.reference_length ||
            countQryLength(prepared.cigars[index]) !=
                prepared.candidate.paths[index].length) {
            prepared.rejection = "msa_invalid";
            return prepared;
        }
    }

    const auto quality_start = Clock::now();
    prepared.new_qualities.reserve(prepared.candidate.paths.size());
    for (size_t index = 0; index < row_ids.size(); ++index) {
        const auto& query_row = rows.at(row_ids[index]);
        const auto quality = evaluateRows(
            reference_row->second, query_row,
            prepared.candidate.reference_length,
            prepared.candidate.paths[index].length);
        prepared.new_qualities.push_back(quality);
        const uint64_t minimum_common = std::min<uint64_t>(
            10, std::min(prepared.candidate.reference_length,
                         prepared.candidate.paths[index].length));
        const double identity_floor = std::max(
            options.minimum_identity,
            prepared.anchor_identity_floors[index] -
                options.maximum_anchor_identity_drop);
        if (quality.common < minimum_common) {
            prepared.rejection = "common_columns_below_threshold";
            break;
        }
        if (quality.coverage + 1e-12 < options.minimum_coverage) {
            prepared.rejection = "coverage_below_threshold";
            break;
        }
        if (quality.identity + 1e-12 < identity_floor) {
            prepared.rejection = prepared.candidate.paths[index].anomalous
                                     ? "identity_below_threshold"
                                     : "protected_species_regressed";
            break;
        }
        if (prepared.candidate.paths[index].anomalous &&
            quality.score <= prepared.old_score) {
            prepared.rejection = "score_not_improved";
            break;
        }
    }
    prepared.quality_ok = prepared.rejection.empty();
    prepared.quality_seconds = secondsSince(quality_start);
    if (!prepared.quality_ok) return prepared;

    prepared.replacement_block = Block::createEmpty(
        prepared.candidate.reference_key.second,
        prepared.candidate.paths.size() + 1);
    Cigar_t reference_cigar;
    appendCigarOp(reference_cigar, 'M',
                  prepared.candidate.reference_length);
    prepared.replacement_segments.push_back(Segment::create(
        prepared.candidate.reference_start,
        prepared.candidate.reference_length, Strand::FORWARD,
        std::move(reference_cigar), AlignRole::PRIMARY,
        SegmentRole::SEGMENT, prepared.replacement_block));
    prepared.replacement_block->anchors.emplace(
        prepared.candidate.reference_key, prepared.replacement_segments.front());
    for (size_t index = 0; index < prepared.candidate.paths.size(); ++index) {
        const auto& path = prepared.candidate.paths[index];
        auto replacement = Segment::create(
            path.start, path.length, path.strand,
            std::move(prepared.cigars[index]), AlignRole::PRIMARY,
            SegmentRole::SEGMENT, prepared.replacement_block);
        prepared.replacement_block->anchors.emplace(path.key, replacement);
        prepared.replacement_segments.push_back(std::move(replacement));
    }
    if (!buildRemovalEdits(prepared.candidate,
                           prepared.replacement_segments, prepared.edits)) {
        prepared.quality_ok = false;
        prepared.rejection = "transaction_footprint_invalid";
        detachPrepared(prepared);
    }
    return prepared;
}

bool editStillValid(const PathEdit& edit) {
    if (!edit.previous || !edit.next) return false;
    if (edit.old_segments.empty()) {
        return edit.previous->primary_path.next.load() == edit.next &&
               edit.next->primary_path.prev.load() == edit.previous;
    }
    if (edit.previous->primary_path.next.load() != edit.old_segments.front() ||
        edit.old_segments.front()->primary_path.prev.load() != edit.previous ||
        edit.old_segments.back()->primary_path.next.load() != edit.next ||
        edit.next->primary_path.prev.load() != edit.old_segments.back()) {
        return false;
    }
    for (size_t i = 1; i < edit.old_segments.size(); ++i) {
        if (edit.old_segments[i - 1]->primary_path.next.load() !=
                edit.old_segments[i] ||
            edit.old_segments[i]->primary_path.prev.load() !=
                edit.old_segments[i - 1]) {
            return false;
        }
    }
    return true;
}

void applyEdit(PathEdit& edit) {
    if (edit.replacement) {
        edit.previous->primary_path.next.store(edit.replacement);
        edit.replacement->primary_path.prev.store(edit.previous);
        edit.replacement->primary_path.next.store(edit.next);
        edit.next->primary_path.prev.store(edit.replacement);
    } else {
        edit.previous->primary_path.next.store(edit.next);
        edit.next->primary_path.prev.store(edit.previous);
    }
}

void rollbackEdit(PathEdit& edit) {
    if (edit.old_segments.empty()) {
        edit.previous->primary_path.next.store(edit.next);
        edit.next->primary_path.prev.store(edit.previous);
    } else {
        edit.previous->primary_path.next.store(edit.old_segments.front());
        edit.old_segments.front()->primary_path.prev.store(edit.previous);
        edit.old_segments.back()->primary_path.next.store(edit.next);
        edit.next->primary_path.prev.store(edit.old_segments.back());
    }
    if (edit.replacement) {
        edit.replacement->primary_path.prev.store(nullptr);
        edit.replacement->primary_path.next.store(nullptr);
    }
}

GenomeEnd* findEnd(RaMeshMultiGenomeGraph& graph,
                   const SpeciesChrPair& key) {
    const auto species = graph.species_graphs.find(key.first);
    if (species == graph.species_graphs.end()) return nullptr;
    const auto chromosome = species->second.chr2end.find(key.second);
    return chromosome == species->second.chr2end.end()
               ? nullptr
               : &chromosome->second;
}

void rebuildSampling(GenomeEnd& end) {
    end.sample_vec.clear();
    end.sample_vec.resize(1, end.head);
    for (auto current = end.head->primary_path.next.load();
         current && !current->isTail();
         current = current->primary_path.next.load()) {
        end.setToSampling(current);
    }
}

bool auditPathsAndBlocks(RaMeshMultiGenomeGraph& graph,
                         const std::set<SpeciesChrPair>& paths) {
    std::unordered_set<const Block*> pool;
    for (const auto& weak : graph.blocks) {
        const auto block = weak.lock();
        if (!block || !pool.insert(block.get()).second) return false;
        for (const auto& [key, segment] : block->anchors) {
            if (!segment || segment->parent_block.get() != block.get()) return false;
            (void)key;
        }
    }
    for (const auto& key : paths) {
        auto* end = findEnd(graph, key);
        if (!end) return false;
        auto previous = end->head;
        uint64_t previous_start = 0;
        bool have_previous = false;
        std::unordered_set<const Segment*> seen;
        for (auto current = previous->primary_path.next.load();
             current && !current->isTail();
             current = current->primary_path.next.load()) {
            if (!current->isSegment() || !seen.insert(current.get()).second ||
                current->primary_path.prev.load() != previous ||
                !current->parent_block ||
                pool.count(current->parent_block.get()) == 0) {
                return false;
            }
            if (have_previous && current->start < previous_start) return false;
            const auto anchor = current->parent_block->anchors.find(key);
            if (anchor == current->parent_block->anchors.end() ||
                anchor->second != current) return false;
            previous_start = current->start;
            have_previous = true;
            previous = current;
        }
        if (previous->primary_path.next.load() != end->tail ||
            end->tail->primary_path.prev.load() != previous) return false;
    }
    return true;
}

bool commitPrepared(RaMeshMultiGenomeGraph& graph, Prepared& prepared) {
    if (!prepared.quality_ok || !prepared.replacement_block) return false;
    for (const auto& edit : prepared.edits) {
        if (!editStillValid(edit)) return false;
    }
    std::unordered_set<const Block*> old_set;
    for (const auto& block : prepared.candidate.interior_blocks) {
        old_set.insert(block.get());
    }
    std::vector<BlockPtr> old_blocks;
    old_blocks.reserve(old_set.size());
    struct Residual {
        BlockPtr block;
        SpeciesChrPair key;
        SegPtr segment;
        BlockPtr old_parent;
        Cigar_t old_cigar;
        Strand old_strand = Strand::FORWARD;
    };
    std::vector<Residual> residuals;
    residuals.reserve(old_set.size());
    for (const auto& [key, segment] : prepared.candidate.anomalous_segments) {
        auto residual = Block::createEmpty(key.second, 1);
        residual->anchors.emplace(key, segment);
        residuals.push_back({residual, key, segment, segment->parent_block,
                             segment->cigar, segment->strand});
    }
    if (residuals.size() != old_set.size()) return false;

    std::vector<WeakBlock> new_pool;
    new_pool.reserve(graph.blocks.size() + 1);
    bool inserted = false;
    size_t residual_index = 0;
    for (const auto& weak : graph.blocks) {
        const auto block = weak.lock();
        if (!block) return false;
        if (old_set.count(block.get()) != 0) {
            old_blocks.push_back(block);
            if (!inserted) {
                new_pool.emplace_back(prepared.replacement_block);
                inserted = true;
            }
            if (residual_index >= residuals.size()) return false;
            new_pool.emplace_back(residuals[residual_index++].block);
        } else {
            new_pool.emplace_back(block);
        }
    }
    if (!inserted || old_blocks.size() != old_set.size() ||
        residual_index != residuals.size()) return false;

    std::set<SpeciesChrPair> affected;
    std::map<SpeciesChrPair, std::vector<SegPtr>> sampling_snapshots;
    for (const auto& edit : prepared.edits) affected.insert(edit.key);
    for (const auto& residual : residuals) affected.insert(residual.key);
    for (const auto& key : affected) {
        auto* end = findEnd(graph, key);
        if (!end) return false;
        sampling_snapshots.emplace(key, end->sample_vec);
    }
    const auto old_pool = graph.blocks;
    for (auto& edit : prepared.edits) applyEdit(edit);
    for (auto& residual : residuals) {
        residual.segment->parent_block = residual.block;
        residual.segment->cigar.clear();
        residual.segment->strand = Strand::FORWARD;
    }
    graph.blocks.swap(new_pool);
    for (const auto& key : affected) rebuildSampling(*findEnd(graph, key));
    if (!auditPathsAndBlocks(graph, affected)) {
        graph.blocks = old_pool;
        for (auto edit = prepared.edits.rbegin(); edit != prepared.edits.rend(); ++edit) {
            rollbackEdit(*edit);
        }
        for (auto& residual : residuals) {
            residual.segment->parent_block = residual.old_parent;
            residual.segment->cigar = std::move(residual.old_cigar);
            residual.segment->strand = residual.old_strand;
            residual.block->anchors.clear();
        }
        for (auto& [key, sampling] : sampling_snapshots) {
            findEnd(graph, key)->sample_vec = std::move(sampling);
        }
        return false;
    }

    for (auto& edit : prepared.edits) {
        for (auto& segment : edit.old_segments) {
            segment->primary_path.prev.store(nullptr);
            segment->primary_path.next.store(nullptr);
            segment->parent_block.reset();
        }
    }
    for (auto& block : old_blocks) block->anchors.clear();
    return true;
}

}  // namespace

bool FailureCache::contains(const std::string& signature) const {
    std::lock_guard lock(mutex_);
    return failures_.contains(signature);
}

void FailureCache::remember(std::string signature, std::string reason) {
    std::lock_guard lock(mutex_);
    failures_.insert_or_assign(std::move(signature), std::move(reason));
}

void FailureCache::clear() {
    std::lock_guard lock(mutex_);
    failures_.clear();
}

Result repairAnchorBoundedStructuralBreaks(
    RaMeshMultiGenomeGraph& graph,
    const SpeciesName& reference_species,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    const Options& options) {
    Result result;
    if (!options.enabled) return result;
    if (options.maximum_span == 0 || options.minimum_outer_anchor == 0 ||
        options.maximum_interior_blocks == 0 ||
        options.parallel_threads == 0 || options.msa_executable.empty()) {
        throw std::invalid_argument("invalid structural-break repair options");
    }

    std::unique_lock graph_lock(graph.rw);
    const auto scan_start = Clock::now();
    auto candidates = scanCandidates(graph, reference_species, options, result);
    if (options.failure_cache) {
        candidates.erase(
            std::remove_if(
                candidates.begin(), candidates.end(),
                [&](const Candidate& candidate) {
                    const bool cached =
                        options.failure_cache->contains(candidate.signature);
                    result.failure_cache_hits += cached;
                    return cached;
                }),
            candidates.end());
    }
    result.scan_seconds = secondsSince(scan_start);
    result.structural_candidates = candidates.size();
    for (const auto& candidate : candidates) {
        ++result.candidates_by_k[candidate.paths.size() + 1];
        result.target_switch += (candidate.signals & TARGET_SWITCH) != 0;
        result.strand_switch += (candidate.signals & STRAND_SWITCH) != 0;
        result.order_break += (candidate.signals & ORDER_BREAK) != 0;
    }
    auto selected = selectConflictFree(std::move(candidates), result);

    std::vector<Prepared> prepared(selected.size());
    if (selected.size() == 1 || options.parallel_threads == 1) {
        for (size_t index = 0; index < selected.size(); ++index) {
            prepared[index] = prepareCandidate(
                std::move(selected[index]), managers, options);
        }
    } else if (!selected.empty()) {
        ThreadPool pool(std::min<size_t>(options.parallel_threads,
                                         selected.size()));
        std::vector<std::future<Prepared>> futures;
        futures.reserve(selected.size());
        for (auto& candidate : selected) {
            futures.push_back(pool.enqueue(
                [&managers, &options](Candidate value) {
                    return prepareCandidate(std::move(value), managers, options);
                },
                std::move(candidate)));
        }
        for (size_t index = 0; index < futures.size(); ++index) {
            prepared[index] = futures[index].get();
        }
    }

    for (auto& item : prepared) {
        result.sequence_seconds += item.sequence_seconds;
        result.msa_seconds += item.msa_seconds;
        result.quality_seconds += item.quality_seconds;
        if (item.sequence_seconds > 0.0) ++result.sequence_prepared;
        if (item.sequence_seconds > 0.0) {
            ++result.prepared_by_k[item.candidate.paths.size() + 1];
        }
        if (item.msa_seconds > 0.0) ++result.msa_calls;
        if (!item.msa_ok && item.rejection == "msa_failed") ++result.msa_failed;
        if (!item.quality_ok) {
            ++result.quality_rejected;
            ++result.rejection_reasons[item.rejection.empty()
                                           ? "quality_rejected"
                                           : item.rejection];
            if (options.failure_cache) {
                options.failure_cache->remember(
                    item.candidate.signature,
                    item.rejection.empty() ? "quality_rejected"
                                           : item.rejection);
            }
            detachPrepared(item);
            continue;
        }
        const auto transaction_start = Clock::now();
        if (commitPrepared(graph, item)) {
            ++result.committed;
            ++result.committed_by_k[item.candidate.paths.size() + 1];
        } else {
            ++result.transaction_rejected;
            ++result.rejection_reasons["transaction_audit_failed"];
            if (options.failure_cache) {
                options.failure_cache->remember(
                    item.candidate.signature,
                    "transaction_audit_failed");
            }
            detachPrepared(item);
        }
        result.transaction_seconds += secondsSince(transaction_start);
    }

    spdlog::info(
        "[structural-break-repair] scanned={} candidates={} selected={} "
        "deferred={} cache_hits={} prepared={} msa_calls={} msa_failed={} quality_rejected={} "
        "transaction_rejected={} committed={} signals(target/strand/order)="
        "{}/{}/{} rejects(anchor/span/ref_empty/interior/gap_only)="
        "{}/{}/{}/{}/{} time(scan/sequence/msa/quality/transaction)="
        "{:.3f}/{:.3f}/{:.3f}/{:.3f}/{:.3f}s",
        result.scanned_windows, result.structural_candidates,
        selected.size(), result.conflict_deferred, result.failure_cache_hits,
        result.sequence_prepared,
        result.msa_calls, result.msa_failed, result.quality_rejected,
        result.transaction_rejected, result.committed, result.target_switch,
        result.strand_switch, result.order_break, result.outer_anchor_invalid,
        result.span_exceeded, result.reference_empty,
        result.too_many_interior_blocks, result.large_gap_only,
        result.scan_seconds,
        result.sequence_seconds, result.msa_seconds, result.quality_seconds,
        result.transaction_seconds);
    for (const auto& [reason, count] : result.rejection_reasons) {
        spdlog::debug("[structural-break-repair] rejection {}={}",
                      reason, count);
    }
    for (const auto& [participants, count] : result.candidates_by_k) {
        spdlog::debug(
            "[structural-break-repair] K={} candidates={} prepared={} "
            "committed={}",
            participants, count, result.prepared_by_k[participants],
            result.committed_by_k[participants]);
    }
    return result;
}

}  // namespace RaMesh::StructuralBreakRepair
