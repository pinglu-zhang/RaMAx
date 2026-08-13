#include "short_block_repair.h"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <deque>
#include <limits>
#include <optional>
#include <set>
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

namespace RaMesh::ShortBlockRepair {
namespace {

using Clock = std::chrono::steady_clock;
enum class Side : uint8_t { Left, Right };

struct Quality {
    uint64_t common = 0;
    uint64_t matches = 0;
    double coverage = 0.0;
    double identity = 0.0;
};

struct PathReplacement {
    SpeciesChrPair key;
    std::vector<SegPtr> old_segments;
    SegPtr replacement;
    SegPtr previous;
    SegPtr next;
};

struct PreparedMerge {
    BlockPtr short_block;
    BlockPtr neighbor_block;
    BlockPtr merged_block;
    SpeciesName reference_species;
    ChrName reference_chromosome;
    uint_t reference_start = 0;
    Side side = Side::Left;
    size_t missing_count = 0;
    double minimum_identity = 1.0;
    double average_coverage = 1.0;
    uint64_t missing_span = 0;
    uint_t neighbor_reference_length = 0;
    size_t reference_rank = 0;
    std::vector<PathReplacement> paths;
    uint64_t ksw2_calls = 0;
    uint64_t ksw2_passed = 0;
};

struct WorkItem {
    SpeciesName reference_species;
    SegPtr reference_segment;
};

void detachPrepared(PreparedMerge& prepared) {
    for (auto& path : prepared.paths) {
        if (!path.replacement) continue;
        path.replacement->primary_path.prev.store(nullptr);
        path.replacement->primary_path.next.store(nullptr);
        path.replacement->parent_block.reset();
    }
    if (prepared.merged_block) prepared.merged_block->anchors.clear();
}

double secondsSince(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

uint64_t segmentEnd(const SegPtr& segment) {
    return static_cast<uint64_t>(segment->start) + segment->length;
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
    end.sample_vec.assign(1, end.head);
    for (auto current = end.head->primary_path.next.load();
         current && !current->isTail();
         current = current->primary_path.next.load()) {
        end.setToSampling(current);
    }
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

bool appendValidated(Cigar_t& destination, const Cigar_t& source) {
    for (const auto unit : source) {
        char operation = '?';
        uint32_t length = 0;
        intToCigar(unit, operation, length);
        if (length == 0 || (operation != 'M' && operation != '=' &&
                            operation != 'X' && operation != 'I' &&
                            operation != 'D')) {
            return false;
        }
        appendCigarOp(destination, operation, length);
    }
    return true;
}

bool normalizedCigar(const SegPtr& segment,
                     uint_t reference_length,
                     Cigar_t& result) {
    result.clear();
    if (!segment || segment->length == 0 || reference_length == 0) return false;
    if (segment->cigar.empty()) {
        if (segment->length != reference_length) return false;
        appendCigarOp(result, 'M', reference_length);
    } else if (!appendValidated(result, segment->cigar)) {
        return false;
    }
    return countRefLength(result) == reference_length &&
           countQryLength(result) == segment->length;
}

bool pathAdjacent(const SegPtr& left, const SegPtr& right) {
    if (!left || !right || left->strand != right->strand) return false;
    if (left->strand == Strand::FORWARD) {
        return left->primary_path.next.load() == right &&
               right->primary_path.prev.load() == left;
    }
    return right->primary_path.next.load() == left &&
           left->primary_path.prev.load() == right;
}

int64_t gapInReferenceOrder(const SegPtr& left, const SegPtr& right) {
    if (!left || !right || left->strand != right->strand) {
        return std::numeric_limits<int64_t>::min();
    }
    return left->strand == Strand::FORWARD
               ? static_cast<int64_t>(right->start) -
                     static_cast<int64_t>(segmentEnd(left))
               : static_cast<int64_t>(left->start) -
                     static_cast<int64_t>(segmentEnd(right));
}

Quality evaluate(const std::string& reference,
                 const std::string& query,
                 const Cigar_t& cigar) {
    Quality quality;
    if (countRefLength(cigar) != reference.size() ||
        countQryLength(cigar) != query.size()) {
        return quality;
    }
    uint64_t reference_position = 0;
    uint64_t query_position = 0;
    for (const auto unit : cigar) {
        char operation = '?';
        uint32_t length = 0;
        intToCigar(unit, operation, length);
        if (operation == 'M' || operation == '=' || operation == 'X') {
            for (uint32_t offset = 0; offset < length; ++offset) {
                ++quality.common;
                if (std::toupper(static_cast<unsigned char>(
                        reference[reference_position + offset])) ==
                    std::toupper(static_cast<unsigned char>(
                        query[query_position + offset]))) {
                    ++quality.matches;
                }
            }
            reference_position += length;
            query_position += length;
        } else if (operation == 'D') {
            reference_position += length;
        } else if (operation == 'I') {
            query_position += length;
        }
    }
    const uint64_t shorter = std::min(reference.size(), query.size());
    quality.coverage = shorter == 0
                           ? 0.0
                           : static_cast<double>(quality.common) / shorter;
    quality.identity = quality.common == 0
                           ? 0.0
                           : static_cast<double>(quality.matches) /
                                 quality.common;
    return quality;
}

std::optional<std::pair<uint_t, uint_t>> missingInterval(
    const SegPtr& neighbor,
    Side side,
    uint_t maximum_span) {
    if (!neighbor) return std::nullopt;
    SegPtr boundary;
    uint64_t start = 0;
    uint64_t end = 0;
    if (side == Side::Left) {
        if (neighbor->strand == Strand::FORWARD) {
            boundary = neighbor->primary_path.next.load();
            if (!boundary || boundary->isTail()) return std::nullopt;
            start = segmentEnd(neighbor);
            end = boundary->start;
        } else {
            boundary = neighbor->primary_path.prev.load();
            if (!boundary || boundary->isHead()) return std::nullopt;
            start = segmentEnd(boundary);
            end = neighbor->start;
        }
    } else if (neighbor->strand == Strand::FORWARD) {
        boundary = neighbor->primary_path.prev.load();
        if (!boundary || boundary->isHead()) return std::nullopt;
        start = segmentEnd(boundary);
        end = neighbor->start;
    } else {
        boundary = neighbor->primary_path.next.load();
        if (!boundary || boundary->isTail()) return std::nullopt;
        start = segmentEnd(neighbor);
        end = boundary->start;
    }
    if (end <= start || end - start > maximum_span ||
        end > std::numeric_limits<uint_t>::max()) {
        return std::nullopt;
    }
    return std::pair{static_cast<uint_t>(start),
                     static_cast<uint_t>(end - start)};
}

bool subsetKeys(const BlockPtr& short_block,
                const BlockPtr& neighbor_block) {
    if (!short_block || !neighbor_block ||
        short_block->anchors.size() > neighbor_block->anchors.size()) {
        return false;
    }
    for (const auto& [key, unused] : short_block->anchors) {
        (void)unused;
        if (!neighbor_block->anchors.contains(key)) return false;
    }
    return true;
}

std::optional<PreparedMerge> prepareMerge(
    const BlockPtr& short_block,
    const BlockPtr& neighbor_block,
    const SpeciesName& reference_species,
    const ChrName& reference_chromosome,
    Side side,
    size_t reference_rank,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    const Options& options,
    Result& result) {
    if (!subsetKeys(short_block, neighbor_block)) {
        ++result.participant_rejected;
        return std::nullopt;
    }
    const SpeciesChrPair reference_key{reference_species,
                                       reference_chromosome};
    const auto short_reference_it = short_block->anchors.find(reference_key);
    const auto neighbor_reference_it = neighbor_block->anchors.find(reference_key);
    if (short_reference_it == short_block->anchors.end() ||
        neighbor_reference_it == neighbor_block->anchors.end()) {
        ++result.path_rejected;
        return std::nullopt;
    }
    const auto short_reference = short_reference_it->second;
    const auto neighbor_reference = neighbor_reference_it->second;
    const auto reference_left = side == Side::Left
                                    ? neighbor_reference
                                    : short_reference;
    const auto reference_right = side == Side::Left
                                     ? short_reference
                                     : neighbor_reference;
    if (short_reference->strand != Strand::FORWARD ||
        neighbor_reference->strand != Strand::FORWARD ||
        gapInReferenceOrder(reference_left, reference_right) != 0 ||
        !pathAdjacent(reference_left, reference_right)) {
        ++result.path_rejected;
        return std::nullopt;
    }

    std::string short_reference_sequence;
    if (!fetchSequence(managers, reference_key, short_reference->start,
                       short_reference->length, Strand::FORWARD,
                       short_reference_sequence)) {
        ++result.path_rejected;
        return std::nullopt;
    }

    PreparedMerge prepared;
    prepared.short_block = short_block;
    prepared.neighbor_block = neighbor_block;
    prepared.reference_species = reference_species;
    prepared.reference_chromosome = reference_chromosome;
    prepared.reference_start = reference_left->start;
    prepared.side = side;
    prepared.reference_rank = reference_rank;
    prepared.neighbor_reference_length = neighbor_reference->length;
    prepared.merged_block = Block::createEmpty(
        neighbor_block->ref_chr, neighbor_block->anchors.size());
    double coverage_sum = 0.0;
    const auto reject = [&]() -> std::optional<PreparedMerge> {
        detachPrepared(prepared);
        return std::nullopt;
    };

    for (const auto& [key, neighbor_segment] : neighbor_block->anchors) {
        if (!neighbor_segment || !neighbor_segment->isPrimary()) {
            ++result.path_rejected;
            return reject();
        }
        const auto short_it = short_block->anchors.find(key);
        SegPtr short_segment = short_it == short_block->anchors.end()
                                   ? nullptr
                                   : short_it->second;
        Cigar_t neighbor_cigar;
        if (!normalizedCigar(neighbor_segment,
                             neighbor_reference->length,
                             neighbor_cigar)) {
            ++result.cigar_rejected;
            return reject();
        }
        Cigar_t short_cigar;
        uint_t missing_start = 0;
        uint_t missing_length = 0;
        if (short_segment) {
            if (short_segment->strand != neighbor_segment->strand) {
                ++result.strand_rejected;
                return reject();
            }
            const auto query_left = side == Side::Left
                                        ? neighbor_segment
                                        : short_segment;
            const auto query_right = side == Side::Left
                                         ? short_segment
                                         : neighbor_segment;
            const int64_t gap = gapInReferenceOrder(query_left, query_right);
            if (gap < 0 || static_cast<uint64_t>(gap) >
                               options.maximum_query_gap ||
                !pathAdjacent(query_left, query_right) ||
                !normalizedCigar(short_segment, short_reference->length,
                                 short_cigar)) {
                ++result.path_rejected;
                return reject();
            }
            if (gap > 0) {
                Cigar_t with_gap;
                if (side == Side::Left) {
                    with_gap = neighbor_cigar;
                    appendCigarOp(with_gap, 'I', static_cast<uint32_t>(gap));
                    appendValidated(with_gap, short_cigar);
                } else {
                    with_gap = short_cigar;
                    appendCigarOp(with_gap, 'I', static_cast<uint32_t>(gap));
                    appendValidated(with_gap, neighbor_cigar);
                }
                short_cigar = std::move(with_gap);
            } else {
                Cigar_t combined;
                if (side == Side::Left) {
                    combined = neighbor_cigar;
                    appendValidated(combined, short_cigar);
                } else {
                    combined = short_cigar;
                    appendValidated(combined, neighbor_cigar);
                }
                short_cigar = std::move(combined);
            }
        } else {
            const auto interval = missingInterval(
                neighbor_segment, side, options.maximum_missing_span);
            if (!interval) {
                ++result.span_rejected;
                return reject();
            }
            missing_start = interval->first;
            missing_length = interval->second;
            std::string query_sequence;
            if (!fetchSequence(managers, key, missing_start, missing_length,
                               neighbor_segment->strand, query_sequence)) {
                ++result.span_rejected;
                return reject();
            }
            const auto ksw_start = Clock::now();
            Cigar_t aligned = globalAlignKSW2BandedPublic(
                short_reference_sequence, query_sequence);
            result.ksw2_seconds += secondsSince(ksw_start);
            ++result.ksw2_calls;
            ++prepared.ksw2_calls;
            const Quality quality = evaluate(
                short_reference_sequence, query_sequence, aligned);
            const uint64_t minimum_common = std::min<uint64_t>(
                10, std::min(short_reference_sequence.size(),
                             query_sequence.size()));
            if (quality.common < minimum_common ||
                quality.coverage + 1e-12 < options.minimum_coverage ||
                quality.identity + 1e-12 < options.minimum_identity) {
                ++result.ksw2_failed;
                return reject();
            }
            ++result.ksw2_passed;
            ++prepared.ksw2_passed;
            ++prepared.missing_count;
            prepared.missing_span += missing_length;
            prepared.minimum_identity = std::min(
                prepared.minimum_identity, quality.identity);
            coverage_sum += quality.coverage;
            if (side == Side::Left) {
                short_cigar = neighbor_cigar;
                appendValidated(short_cigar, aligned);
            } else {
                short_cigar = std::move(aligned);
                appendValidated(short_cigar, neighbor_cigar);
            }
        }

        uint64_t merged_start = neighbor_segment->start;
        uint64_t merged_end = segmentEnd(neighbor_segment);
        std::vector<SegPtr> old_segments{neighbor_segment};
        if (short_segment) {
            merged_start = std::min<uint64_t>(merged_start,
                                              short_segment->start);
            merged_end = std::max<uint64_t>(merged_end,
                                            segmentEnd(short_segment));
            old_segments.push_back(short_segment);
        } else {
            merged_start = std::min<uint64_t>(merged_start, missing_start);
            merged_end = std::max<uint64_t>(merged_end,
                                            static_cast<uint64_t>(missing_start) +
                                                missing_length);
        }
        if (merged_end <= merged_start ||
            merged_end > std::numeric_limits<uint_t>::max() ||
            countRefLength(short_cigar) !=
                neighbor_reference->length + short_reference->length ||
            countQryLength(short_cigar) != merged_end - merged_start) {
            ++result.cigar_rejected;
            return reject();
        }
        auto replacement = Segment::create(
            static_cast<uint_t>(merged_start),
            static_cast<uint_t>(merged_end - merged_start),
            neighbor_segment->strand, std::move(short_cigar),
            AlignRole::PRIMARY, SegmentRole::SEGMENT,
            prepared.merged_block);
        prepared.merged_block->anchors.emplace(key, replacement);
        prepared.paths.push_back(
            {key, std::move(old_segments), replacement, nullptr, nullptr});
    }
    if (prepared.missing_count == 0) {
        prepared.minimum_identity = 1.0;
        prepared.average_coverage = 1.0;
    } else {
        prepared.average_coverage = coverage_sum / prepared.missing_count;
    }
    return prepared;
}

bool betterCandidate(const PreparedMerge& left,
                     const PreparedMerge& right) {
    return std::tuple{left.missing_count,
                      -left.minimum_identity,
                      -left.average_coverage,
                      left.missing_span,
                      -static_cast<int64_t>(left.neighbor_reference_length),
                      left.reference_rank,
                      left.reference_start,
                      left.side == Side::Right} <
           std::tuple{right.missing_count,
                      -right.minimum_identity,
                      -right.average_coverage,
                      right.missing_span,
                      -static_cast<int64_t>(right.neighbor_reference_length),
                      right.reference_rank,
                      right.reference_start,
                      right.side == Side::Right};
}

bool auditGraph(RaMeshMultiGenomeGraph& graph,
                const std::set<SpeciesChrPair>& affected) {
    std::unordered_set<const Block*> pool;
    for (const auto& weak : graph.blocks) {
        const auto block = weak.lock();
        if (!block || !pool.insert(block.get()).second) return false;
        for (const auto& [key, segment] : block->anchors) {
            if (!segment || segment->parent_block.get() != block.get()) {
                return false;
            }
            (void)key;
        }
    }
    for (const auto& key : affected) {
        auto* end = findEnd(graph, key);
        if (!end) return false;
        auto previous = end->head;
        uint64_t previous_start = 0;
        bool have_previous = false;
        std::unordered_set<const Segment*> seen;
        for (auto current = previous->primary_path.next.load();
             current && !current->isTail();
             current = current->primary_path.next.load()) {
            if (!seen.insert(current.get()).second ||
                current->primary_path.prev.load() != previous ||
                !current->parent_block ||
                !pool.contains(current->parent_block.get()) ||
                (have_previous && current->start < previous_start)) {
                return false;
            }
            const auto anchor = current->parent_block->anchors.find(key);
            if (anchor == current->parent_block->anchors.end() ||
                anchor->second != current) {
                return false;
            }
            previous_start = current->start;
            have_previous = true;
            previous = current;
        }
        if (previous->primary_path.next.load() != end->tail ||
            end->tail->primary_path.prev.load() != previous) {
            return false;
        }
    }
    return true;
}

bool commitMerge(RaMeshMultiGenomeGraph& graph,
                 PreparedMerge& prepared,
                 const std::unordered_set<const Block*>& active_blocks) {
    if (!prepared.short_block || !prepared.neighbor_block ||
        !active_blocks.contains(prepared.short_block.get()) ||
        !active_blocks.contains(prepared.neighbor_block.get())) {
        return false;
    }
    std::set<SpeciesChrPair> affected;
    std::map<SpeciesChrPair, std::vector<SegPtr>> sampling;
    for (auto& path : prepared.paths) {
        std::sort(path.old_segments.begin(), path.old_segments.end(),
                  [](const SegPtr& left, const SegPtr& right) {
                      return left->start < right->start;
                  });
        path.previous = path.old_segments.front()->primary_path.prev.load();
        path.next = path.old_segments.back()->primary_path.next.load();
        if (!path.previous || !path.next) return false;
        auto current = path.previous->primary_path.next.load();
        for (const auto& segment : path.old_segments) {
            if (current != segment) return false;
            current = current->primary_path.next.load();
        }
        if (current != path.next) return false;
        affected.insert(path.key);
    }
    for (const auto& key : affected) {
        auto* end = findEnd(graph, key);
        if (!end) return false;
        sampling.emplace(key, end->sample_vec);
    }
    const auto old_pool = graph.blocks;
    std::vector<WeakBlock> new_pool;
    new_pool.reserve(graph.blocks.size() - 1);
    bool inserted = false;
    for (const auto& weak : graph.blocks) {
        const auto block = weak.lock();
        if (!block) return false;
        if (block == prepared.short_block || block == prepared.neighbor_block) {
            if (!inserted) {
                new_pool.emplace_back(prepared.merged_block);
                inserted = true;
            }
        } else {
            new_pool.emplace_back(block);
        }
    }
    if (!inserted) return false;
    for (auto& path : prepared.paths) {
        path.previous->primary_path.next.store(path.replacement);
        path.replacement->primary_path.prev.store(path.previous);
        path.replacement->primary_path.next.store(path.next);
        path.next->primary_path.prev.store(path.replacement);
    }
    graph.blocks.swap(new_pool);
    for (const auto& key : affected) rebuildSampling(*findEnd(graph, key));
    if (!auditGraph(graph, affected)) {
        graph.blocks = old_pool;
        for (auto& path : prepared.paths) {
            path.previous->primary_path.next.store(path.old_segments.front());
            path.old_segments.front()->primary_path.prev.store(path.previous);
            for (size_t index = 1; index < path.old_segments.size(); ++index) {
                path.old_segments[index - 1]->primary_path.next.store(
                    path.old_segments[index]);
                path.old_segments[index]->primary_path.prev.store(
                    path.old_segments[index - 1]);
            }
            path.old_segments.back()->primary_path.next.store(path.next);
            path.next->primary_path.prev.store(path.old_segments.back());
            path.replacement->primary_path.prev.store(nullptr);
            path.replacement->primary_path.next.store(nullptr);
        }
        for (auto& [key, snapshot] : sampling) {
            findEnd(graph, key)->sample_vec = std::move(snapshot);
        }
        return false;
    }
    for (auto& path : prepared.paths) {
        for (auto& segment : path.old_segments) {
            segment->primary_path.prev.store(nullptr);
            segment->primary_path.next.store(nullptr);
            segment->parent_block.reset();
        }
    }
    prepared.short_block->anchors.clear();
    prepared.neighbor_block->anchors.clear();
    return true;
}

bool commitDeleteBatch(RaMeshMultiGenomeGraph& graph,
                       const std::vector<BlockPtr>& blocks,
                       const std::unordered_set<const Block*>& active_blocks,
                       Result& result) {
    if (blocks.empty()) return true;
    std::unordered_set<const Block*> targets;
    std::set<SpeciesChrPair> affected;
    for (const auto& block : blocks) {
        if (!block || !active_blocks.contains(block.get()) ||
            !targets.insert(block.get()).second) {
            continue;
        }
        for (const auto& [key, segment] : block->anchors) {
            if (!segment || segment->parent_block != block) return false;
            affected.insert(key);
        }
    }
    if (targets.empty()) return true;
    std::map<SpeciesChrPair, std::vector<SegPtr>> sampling;
    std::map<SpeciesChrPair, std::vector<SegPtr>> old_paths;
    for (const auto& key : affected) {
        auto* end = findEnd(graph, key);
        if (!end) return false;
        sampling.emplace(key, end->sample_vec);
        auto& path = old_paths[key];
        for (auto current = end->head->primary_path.next.load();
             current && !current->isTail();
             current = current->primary_path.next.load()) {
            path.push_back(current);
        }
    }
    const auto old_pool = graph.blocks;
    std::vector<WeakBlock> new_pool;
    new_pool.reserve(graph.blocks.size() - targets.size());
    for (const auto& weak : graph.blocks) {
        const auto block = weak.lock();
        if (!block) return false;
        if (!targets.contains(block.get())) new_pool.emplace_back(block);
    }
    for (const auto& [key, old_path] : old_paths) {
        auto* end = findEnd(graph, key);
        auto previous = end->head;
        for (const auto& segment : old_path) {
            if (targets.contains(segment->parent_block.get())) continue;
            previous->primary_path.next.store(segment);
            segment->primary_path.prev.store(previous);
            previous = segment;
        }
        previous->primary_path.next.store(end->tail);
        end->tail->primary_path.prev.store(previous);
    }
    for (const auto& block : blocks) {
        if (!targets.contains(block.get())) continue;
        for (const auto& [key, segment] : block->anchors) {
            (void)key;
            segment->parent_block.reset();
        }
    }
    graph.blocks.swap(new_pool);
    for (const auto& key : affected) rebuildSampling(*findEnd(graph, key));
    if (!auditGraph(graph, affected)) {
        graph.blocks = old_pool;
        for (const auto& block : blocks) {
            if (!targets.contains(block.get())) continue;
            for (const auto& [key, segment] : block->anchors) {
                (void)key;
                segment->parent_block = block;
            }
        }
        for (const auto& [key, old_path] : old_paths) {
            auto* end = findEnd(graph, key);
            auto previous = end->head;
            for (const auto& segment : old_path) {
                previous->primary_path.next.store(segment);
                segment->primary_path.prev.store(previous);
                previous = segment;
            }
            previous->primary_path.next.store(end->tail);
            end->tail->primary_path.prev.store(previous);
        }
        for (auto& [key, saved] : sampling) {
            findEnd(graph, key)->sample_vec = std::move(saved);
        }
        return false;
    }
    for (const auto& block : blocks) {
        if (!targets.contains(block.get())) continue;
        for (const auto& [key, segment] : block->anchors) {
            ++result.deleted_segments;
            result.deleted_bases_by_species[key.first] += segment->length;
            segment->primary_path.prev.store(nullptr);
            segment->primary_path.next.store(nullptr);
        }
        block->anchors.clear();
        ++result.deleted_blocks;
    }
    return true;
}

void countBlocks(const RaMeshMultiGenomeGraph& graph,
                 const std::vector<SpeciesName>& reference_order,
                 uint64_t& total,
                 uint64_t& le10,
                 uint64_t& le50,
                 uint64_t& le100) {
    total = le10 = le50 = le100 = 0;
    std::unordered_set<SpeciesName> references(reference_order.begin(),
                                               reference_order.end());
    for (const auto& weak : graph.blocks) {
        const auto block = weak.lock();
        if (!block) continue;
        ++total;
        uint_t minimum = std::numeric_limits<uint_t>::max();
        for (const auto& [key, segment] : block->anchors) {
            if (references.contains(key.first)) {
                minimum = std::min(minimum, segment->length);
            }
        }
        if (minimum <= 10) ++le10;
        if (minimum <= 50) ++le50;
        if (minimum <= 100) ++le100;
    }
}

}  // namespace

Result repairFinalShortBlocks(
    RaMeshMultiGenomeGraph& graph,
    const std::vector<SpeciesName>& reference_order,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    const Options& options) {
    Result result;
    if (!options.enabled) return result;
    if (reference_order.empty() || options.maximum_short_length == 0 ||
        options.maximum_missing_span == 0 || options.parallel_threads == 0) {
        throw std::invalid_argument("invalid short-Block repair options");
    }
    const auto total_start = Clock::now();
    std::unique_lock graph_lock(graph.rw);
    countBlocks(graph, reference_order, result.blocks_before,
                result.blocks_le_10_before, result.blocks_le_50_before,
                result.blocks_le_100_before);

    std::unordered_map<SpeciesName, size_t> ranks;
    for (size_t index = 0; index < reference_order.size(); ++index) {
        ranks.emplace(reference_order[index], index);
    }
    std::unordered_set<const Block*> active_blocks;
    active_blocks.reserve(graph.blocks.size() * 2 + 1);
    for (const auto& weak : graph.blocks) {
        const auto block = weak.lock();
        if (block) active_blocks.insert(block.get());
    }
    std::deque<WorkItem> queue;
    std::unordered_set<const Block*> initial_short_blocks;
    const auto scan_start = Clock::now();
    for (const auto& reference_species : reference_order) {
        const auto species = graph.species_graphs.find(reference_species);
        if (species == graph.species_graphs.end()) continue;
        ++result.reference_paths;
        for (const auto& [chromosome, end] : species->second.chr2end) {
            (void)chromosome;
            for (auto current = end.head->primary_path.next.load();
                 current && !current->isTail();
                 current = current->primary_path.next.load()) {
                ++result.scanned_blocks;
                if (current->length <= options.maximum_short_length) {
                    queue.push_back({reference_species, current});
                    initial_short_blocks.insert(current->parent_block.get());
                }
            }
        }
    }
    result.unique_short_blocks = initial_short_blocks.size();
    result.scan_seconds = secondsSince(scan_start);

    uint64_t generation = 0;
    std::unordered_map<std::string, uint64_t> attempted;
    auto enqueueDirty = [&](const BlockPtr& block) {
        if (!block) return;
        for (const auto& reference_species : reference_order) {
            for (const auto& [key, segment] : block->anchors) {
                if (key.first == reference_species &&
                    segment->length <= options.maximum_short_length) {
                    queue.push_back({reference_species, segment});
                }
            }
        }
    };

    while (!queue.empty()) {
        WorkItem item = std::move(queue.front());
        queue.pop_front();
        const auto segment = item.reference_segment;
        const auto block = segment ? segment->parent_block : nullptr;
        if (!block || !active_blocks.contains(block.get()) ||
            segment->length > options.maximum_short_length) {
            continue;
        }
        SpeciesChrPair reference_key;
        bool found_reference = false;
        for (const auto& [key, current] : block->anchors) {
            if (key.first == item.reference_species && current == segment) {
                reference_key = key;
                found_reference = true;
                break;
            }
        }
        if (!found_reference) continue;
        std::ostringstream signature;
        signature << block.get() << '|' << item.reference_species << '|'
                  << generation;
        if (attempted.contains(signature.str())) continue;
        attempted.emplace(signature.str(), generation);

        std::optional<PreparedMerge> left;
        std::optional<PreparedMerge> right;
        const auto previous = segment->primary_path.prev.load();
        if (previous && !previous->isHead() && previous->parent_block != block) {
            ++result.left_candidates;
            left = prepareMerge(block, previous->parent_block,
                                item.reference_species, reference_key.second,
                                Side::Left, ranks.at(item.reference_species),
                                managers, options, result);
        }
        const auto next = segment->primary_path.next.load();
        if (next && !next->isTail() && next->parent_block != block) {
            ++result.right_candidates;
            right = prepareMerge(block, next->parent_block,
                                 item.reference_species, reference_key.second,
                                 Side::Right, ranks.at(item.reference_species),
                                 managers, options, result);
        }
        if (!left && !right) continue;
        PreparedMerge selected;
        if (left && right) {
            if (betterCandidate(*left, *right)) {
                selected = std::move(*left);
                detachPrepared(*right);
            } else {
                selected = std::move(*right);
                detachPrepared(*left);
            }
        } else if (left) {
            selected = std::move(*left);
        } else {
            selected = std::move(*right);
        }
        const auto transaction_start = Clock::now();
        if (!commitMerge(graph, selected, active_blocks)) {
            ++result.transaction_rollbacks;
            result.transaction_seconds += secondsSince(transaction_start);
            detachPrepared(selected);
            continue;
        }
        result.transaction_seconds += secondsSince(transaction_start);
        active_blocks.erase(selected.short_block.get());
        active_blocks.erase(selected.neighbor_block.get());
        active_blocks.insert(selected.merged_block.get());
        ++generation;
        ++result.fixed_point_generations;
        if (selected.side == Side::Left) ++result.left_merged;
        else ++result.right_merged;
        enqueueDirty(selected.merged_block);
        for (const auto& [key, merged] : selected.merged_block->anchors) {
            (void)key;
            const auto dirty_previous = merged->primary_path.prev.load();
            const auto dirty_next = merged->primary_path.next.load();
            if (dirty_previous && !dirty_previous->isHead()) {
                enqueueDirty(dirty_previous->parent_block);
            }
            if (dirty_next && !dirty_next->isTail()) {
                enqueueDirty(dirty_next->parent_block);
            }
        }
    }

    std::unordered_set<const Block*> delete_seen;
    std::vector<BlockPtr> delete_blocks;
    for (const auto& reference_species : reference_order) {
        const auto species = graph.species_graphs.find(reference_species);
        if (species == graph.species_graphs.end()) continue;
        for (const auto& [chromosome, end] : species->second.chr2end) {
            (void)chromosome;
            for (auto current = end.head->primary_path.next.load();
                 current && !current->isTail();
                 current = current->primary_path.next.load()) {
                if (current->length <= options.maximum_short_length &&
                    current->parent_block &&
                    delete_seen.insert(current->parent_block.get()).second) {
                    delete_blocks.push_back(current->parent_block);
                }
            }
        }
    }
    const auto delete_start = Clock::now();
    if (!commitDeleteBatch(graph, delete_blocks, active_blocks, result)) {
        ++result.transaction_rollbacks;
    }
    result.transaction_seconds += secondsSince(delete_start);
    countBlocks(graph, reference_order, result.blocks_after,
                result.blocks_le_10_after, result.blocks_le_50_after,
                result.blocks_le_100_after);
    result.total_seconds = secondsSince(total_start);
    spdlog::info(
        "[short-block-repair] references={} scanned={} unique_short={} "
        "candidates(left/right)={}/{} rejects(participant/path/strand/span/cigar)="
        "{}/{}/{}/{}/{} ksw2(calls/pass/fail)={}/{}/{} "
        "merged(left/right)={}/{} generations={} deleted(blocks/segments)="
        "{}/{} rollbacks={} blocks(before/after)={}/{} "
        "short<=10/50/100(before)={}/{}/{} short<=10/50/100(after)="
        "{}/{}/{} time(scan/ksw2/transaction/total)="
        "{:.3f}/{:.3f}/{:.3f}/{:.3f}s",
        result.reference_paths, result.scanned_blocks,
        result.unique_short_blocks, result.left_candidates,
        result.right_candidates, result.participant_rejected,
        result.path_rejected, result.strand_rejected, result.span_rejected,
        result.cigar_rejected, result.ksw2_calls, result.ksw2_passed,
        result.ksw2_failed, result.left_merged, result.right_merged,
        result.fixed_point_generations, result.deleted_blocks,
        result.deleted_segments, result.transaction_rollbacks,
        result.blocks_before, result.blocks_after,
        result.blocks_le_10_before, result.blocks_le_50_before,
        result.blocks_le_100_before, result.blocks_le_10_after,
        result.blocks_le_50_after, result.blocks_le_100_after,
        result.scan_seconds, result.ksw2_seconds,
        result.transaction_seconds, result.total_seconds);
    for (const auto& [species, bases] : result.deleted_bases_by_species) {
        spdlog::info("[short-block-repair] deleted_species={} bases={}",
                     species, bases);
    }
    return result;
}

}  // namespace RaMesh::ShortBlockRepair
