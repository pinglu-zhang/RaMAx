#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

#include "SeqPro.h"
#include "align.h"
#include "anchor.h"
#include "rare_aligner.h"
#include "ramesh.h"
#include "../src/anchor/anchor_link_internal.h"
#include "../src/graph/merge_internal.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <type_traits>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

namespace {

std::filesystem::path makeTestDirectory() {
    const auto stamp = std::chrono::steady_clock::now()
        .time_since_epoch().count();
    auto path = std::filesystem::temp_directory_path() /
        ("ramax-core-performance-" + std::to_string(stamp));
    std::filesystem::create_directories(path);
    return path;
}

void writeFile(const std::filesystem::path& path,
               const std::string& contents) {
    std::ofstream output(path, std::ios::binary);
    assert(output);
    output.write(contents.data(), static_cast<std::streamsize>(contents.size()));
    assert(output.good());
}

void testSequenceViewsAndKswSummary(const std::filesystem::path& directory) {
    const auto single_line = directory / "single.fa";
    const auto wrapped = directory / "wrapped.fa";
    writeFile(single_line, ">chr1\nACGTACGTACGT\n");
    writeFile(wrapped, ">chr1\r\nACGT\r\nACGT\r\nACGT\r\n");

    SeqPro::SequenceManager single_manager(single_line);
    SeqPro::SequenceManager wrapped_manager(wrapped);
    std::string_view view;
    assert(single_manager.tryGetContiguousSubSequence(0, 2, 8, view));
    assert(view == "GTACGTAC");
    assert(!wrapped_manager.tryGetContiguousSubSequence(0, 2, 8, view));
    std::string buffer = "retained-capacity";
    wrapped_manager.getSubSequenceInto(0, 2, 8, buffer);
    assert(buffer == wrapped_manager.getSubSequence(0, 2, 8));
    assert(buffer == "GTACGTAC");

    const std::vector<std::pair<std::string, std::string>> cases = {
        {"", ""}, {"ACGT", "ACGT"}, {"ACNT", "ACGT"},
        {"acgt", "acgt"}, {"ACGTT", "ACGT"}, {"ACGT", "ACGTT"}};
    for (const auto& [reference, query] : cases) {
        AlignmentResult result = extendAlignKSW2Result(
            {reference, false}, {query, false}, 400);
        const CigarSummary expected = summarizeCigar(result.cigar);
        assert(result.summary.reference_length == expected.reference_length);
        assert(result.summary.query_length == expected.query_length);
        assert(result.summary.alignment_length == expected.alignment_length);
        assert(result.summary.match_length == expected.match_length);
        assert(result.summary.reference_length == countRefLength(result.cigar));
        assert(result.summary.query_length == countQryLength(result.cigar));
        assert(result.summary.alignment_length == countAlignmentLength(result.cigar));
        assert(result.summary.match_length == countMatchOperations(result.cigar));
        if (!reference.empty() || !query.empty()) {
            AlignmentResult raw = ramaxExtendAlignKSW2RawForTesting(
                {reference, false}, {query, false}, 400);
            assert(result.cigar == raw.cigar);
        }
    }
    AlignmentResult reverse = extendAlignKSW2Result(
        {"AACG", false}, {"CGTT", true}, 400);
    assert(reverse.cigar.size() == 1);
    assert(cigarToString(reverse.cigar) == "4M");

    std::mt19937 generator(0xc1a01234u);
    constexpr std::string_view alphabet = "ACGTNacgt";
    for (size_t fixture = 0; fixture < 10000; ++fixture) {
        std::string reference(1 + generator() % 96, 'A');
        std::string query(1 + generator() % 96, 'A');
        for (char& base : reference) base = alphabet[generator() % alphabet.size()];
        for (char& base : query) base = alphabet[generator() % alphabet.size()];
        const bool reverse_query = generator() % 2 != 0;
        AlignmentResult optimized = extendAlignKSW2Result(
            {reference, false}, {query, reverse_query}, 400);
        AlignmentResult raw = ramaxExtendAlignKSW2RawForTesting(
            {reference, false}, {query, reverse_query}, 400);
        assert(optimized.cigar == raw.cigar);
        assert(optimized.summary.reference_length ==
               summarizeCigar(optimized.cigar).reference_length);
        if (fixture < 1000) {
            std::string oriented_query = query;
            if (reverse_query) reverseComplement(oriented_query);
            AlignmentResult global_view = globalAlignKSW2Result(
                {reference, false}, {query, reverse_query});
            assert(global_view.cigar ==
                   globalAlignKSW2(reference, oriented_query));
            assert(global_view.summary.match_length ==
                   countMatchOperations(global_view.cigar));
        }
    }

    const std::string exact(10000, 'A');
    AlignmentResult exact_optimized = extendAlignKSW2Result(
        {exact, false}, {exact, false}, 400);
    AlignmentResult exact_raw = ramaxExtendAlignKSW2RawForTesting(
        {exact, false}, {exact, false}, 400);
    assert(exact_optimized.cigar == exact_raw.cigar);
}

Anchor makeLinkingAnchor(uint32_t ref_start, uint32_t query_start,
                         Strand strand) {
    constexpr uint32_t length = 20;
    return Anchor(
        0, ref_start, length, 0, query_start, length, strand,
        length, length, Cigar_t{cigarToInt('M', length)});
}

AnchorVec makeLinkingPair(uint32_t gap, Strand strand) {
    constexpr uint32_t length = 20;
    constexpr uint32_t first_ref = 100;
    constexpr uint32_t first_forward_query = 100;
    constexpr uint32_t first_reverse_query = 220000;
    const uint32_t second_ref = first_ref + length + gap;
    const uint32_t first_query = strand == FORWARD
        ? first_forward_query : first_reverse_query;
    const uint32_t second_query = strand == FORWARD
        ? first_forward_query + length + gap
        : first_reverse_query - length - gap;
    AnchorVec anchors;
    anchors.push_back(makeLinkingAnchor(first_ref, first_query, strand));
    anchors.push_back(makeLinkingAnchor(second_ref, second_query, strand));
    return anchors;
}

void testAnchorLinkGapBoundaries(const std::filesystem::path& directory) {
    const auto reference_path = directory / "anchor-link-reference.fa";
    const auto query_path = directory / "anchor-link-query.fa";
    writeFile(reference_path,
              ">chr1\n" + std::string(250000, 'A') + "\n");
    writeFile(query_path,
              ">chr1\n" + std::string(250000, 'A') + "\n");
    SeqPro::ManagerVariant reference_manager{
        std::make_unique<SeqPro::SequenceManager>(reference_path)};
    SeqPro::ManagerVariant query_manager{
        std::make_unique<SeqPro::SequenceManager>(query_path)};

    for (const uint32_t gap :
         {0U, 1U, 9999U, 10000U, 10001U, 50000U, 100000U, 100001U}) {
        AnchorVec anchors = makeLinkingPair(gap, FORWARD);
        AnchorLinkDetail::Statistics split_statistics;
        const auto components =
            AnchorLinkDetail::splitAnchorComponents(anchors, &split_statistics);
        const size_t expected_components = gap <= 10000 ? 1 : 2;
        assert(components.size() == expected_components);
        assert(components.front().begin == 0);
        assert(components.back().end == 2);
        assert(split_statistics.long_gap_rejections ==
               (gap <= 10000 ? 0U : 1U));
        assert(split_statistics.maximum_seen_gap ==
               (gap <= 10000 ? 0U : gap));

        AnchorLinkDetail::Statistics statistics;
        auto output = AnchorLinkDetail::linkAnchorRange(
            anchors, 0, anchors.size(), reference_manager, query_manager,
            &statistics);
        if (gap <= 10000) {
            assert(output.size() == 1);
            assert(statistics.sequence_extractions == 2);
            assert(statistics.direct_ksw_calls == 1);
            assert(statistics.fallback_ksw_calls == 0);
            assert(statistics.long_gap_rejections == 0);
            assert(output.front()->ref_len == 40 + gap);
            assert(output.front()->qry_len == 40 + gap);
        } else {
            assert(output.size() == 2);
            assert(statistics.sequence_extractions == 0);
            assert(statistics.direct_ksw_calls == 0);
            assert(statistics.fallback_ksw_calls == 0);
            assert(statistics.long_gap_rejections >= 1);
            assert(statistics.maximum_seen_gap == gap);
        }
    }

    for (const uint32_t gap : {9999U, 10000U, 10001U, 100000U}) {
        const AnchorVec reverse = makeLinkingPair(gap, REVERSE);
        const auto components =
            AnchorLinkDetail::splitAnchorComponents(reverse);
        assert(components.size() == (gap <= 10000 ? 1U : 2U));
    }

    AnchorVec prefix_maximum;
    prefix_maximum.emplace_back(
        0, 100, 20000, 0, 100, 20000, FORWARD,
        20000, 20000, Cigar_t{cigarToInt('M', 20000)});
    prefix_maximum.push_back(makeLinkingAnchor(200, 200, FORWARD));
    prefix_maximum.push_back(makeLinkingAnchor(30100, 30100, FORWARD));
    AnchorLinkDetail::Statistics split_statistics;
    auto components =
        AnchorLinkDetail::splitAnchorComponents(
            prefix_maximum, &split_statistics);
    assert(components.size() == 1);
    assert(split_statistics.long_gap_rejections == 0);
    prefix_maximum.back().ref_start = 30101;
    split_statistics = {};
    components = AnchorLinkDetail::splitAnchorComponents(
        prefix_maximum, &split_statistics);
    assert(components.size() == 2);
    assert(split_statistics.long_gap_rejections == 1);
    assert(split_statistics.maximum_seen_gap == 10001);
    assert(components[0].begin == 0 && components[0].end == 2);
    assert(components[1].begin == 2 && components[1].end == 3);
}

AnchorPtr makeAnchor(uint32_t ref_start, uint32_t ref_length,
                     uint32_t query_start, uint32_t query_length,
                     uint32_t alignment_length,
                     uint32_t aligned_bases) {
    return std::make_shared<Anchor>(
        0, ref_start, ref_length, 0, query_start, query_length,
        FORWARD, alignment_length, aligned_bases,
        Cigar_t{cigarToInt('M', std::max<uint32_t>(1, aligned_bases))});
}

AnchorPtrVec cloneAnchors(const AnchorPtrVec& source) {
    AnchorPtrVec result;
    result.reserve(source.size());
    for (const auto& anchor : source) {
        result.push_back(std::make_shared<Anchor>(*anchor));
        result.back()->ref_selected = false;
        result.back()->qry_selected = false;
    }
    return result;
}

void compareDpSelection(const AnchorPtrVec& input, bool filter_reference) {
    AnchorPtrVec legacy = cloneAnchors(input);
    AnchorPtrVec optimized = cloneAnchors(input);
    ramaxFilterAnchorsByDpLegacyForTesting(legacy, filter_reference);
    ramaxFilterAnchorsByDpOptimizedForTesting(optimized, filter_reference);
    assert(legacy.size() == optimized.size());
    for (size_t index = 0; index < legacy.size(); ++index) {
        assert(legacy[index]->ref_selected == optimized[index]->ref_selected);
        assert(legacy[index]->qry_selected == optimized[index]->qry_selected);
    }
}

void testDpEquivalence() {
    std::mt19937 generator(0x5eed1234u);
    for (size_t fixture = 0; fixture < 10000; ++fixture) {
        const size_t count = 1 + generator() % 64;
        AnchorPtrVec anchors;
        anchors.reserve(count);
        for (size_t index = 0; index < count; ++index) {
            const uint32_t ref_start = generator() % 2000;
            const uint32_t query_start = generator() % 2000;
            const uint32_t ref_length = generator() % 80;
            const uint32_t query_length = generator() % 80;
            const uint32_t alignment_length = 1 + generator() % 120;
            const uint32_t aligned_bases = generator() %
                (alignment_length + 1);
            anchors.push_back(makeAnchor(
                ref_start, ref_length, query_start, query_length,
                alignment_length, aligned_bases));
        }
        compareDpSelection(anchors, true);
        compareDpSelection(anchors, false);
    }

    AnchorPtrVec boundary;
    boundary.reserve(5002);
    for (uint32_t index = 0; index < 5002; ++index) {
        boundary.push_back(makeAnchor(
            index * 3, index % 17 == 0 ? 0 : 2,
            index * 4, index % 19 == 0 ? 0 : 3,
            50 + index % 11, 40 + index % 7));
    }
    compareDpSelection(boundary, true);
    compareDpSelection(boundary, false);

    AnchorPtrVec anomalous = {
        makeAnchor(0, 0, 0, 0, 0, 0),
        makeAnchor(0, 10, 0, 10, 10, 10),
        makeAnchor(10, 10, 10, 10, 10, 10),
        makeAnchor(20, 0, 20, 0, 1, 0)};
    compareDpSelection(anomalous, true);
    compareDpSelection(anomalous, false);
}

void testDpPerformanceSmoke() {
    AnchorPtrVec input;
    input.reserve(50000);
    for (uint32_t index = 0; index < 50000; ++index) {
        input.push_back(makeAnchor(
            index * 3, 2, index * 4, 3, 100, 90));
    }
    AnchorPtrVec legacy = cloneAnchors(input);
    AnchorPtrVec optimized = cloneAnchors(input);
    const auto legacy_start = std::chrono::steady_clock::now();
    ramaxFilterAnchorsByDpLegacyForTesting(legacy, true);
    const auto optimized_start = std::chrono::steady_clock::now();
    ramaxFilterAnchorsByDpOptimizedForTesting(optimized, true);
    const auto finish = std::chrono::steady_clock::now();
    const double legacy_ms = std::chrono::duration<double, std::milli>(
        optimized_start - legacy_start).count();
    const double optimized_ms = std::chrono::duration<double, std::milli>(
        finish - optimized_start).count();
    for (size_t index = 0; index < legacy.size(); ++index) {
        assert(legacy[index]->ref_selected == optimized[index]->ref_selected);
    }
    const double speedup = legacy_ms / std::max(optimized_ms, 0.001);
    std::cout << "DP 50000-anchor speedup: " << speedup
              << "x (legacy=" << legacy_ms
              << " ms, optimized=" << optimized_ms << " ms)\n";
    assert(speedup >= 8.0);
}

using SegmentState = std::tuple<uint32_t, uint32_t, Strand,
                                std::string, size_t>;

std::vector<SegmentState> graphChromosomeState(
    RaMesh::RaMeshMultiGenomeGraph& graph,
    const SpeciesName& species, const ChrName& chromosome) {
    std::unordered_map<uint64_t, size_t> block_order;
    for (size_t index = 0; index < graph.blocks.size(); ++index) {
        if (auto block = graph.blocks[index].lock()) {
            block_order.emplace(block->block_id, index);
        }
    }
    std::vector<SegmentState> state;
    auto& end = graph.species_graphs.at(species).chr2end.at(chromosome);
    for (auto segment = end.head->primary_path.next.load();
         segment && !segment->isTail();
         segment = segment->primary_path.next.load()) {
        assert(segment->parent_block);
        state.emplace_back(segment->start, segment->length,
            segment->strand, cigarToString(segment->cigar),
            block_order.at(segment->parent_block->block_id));
    }
    return state;
}

using ExtendedSegmentState =
    std::tuple<uint_t, uint_t, Strand, std::string, bool, bool, size_t>;

std::vector<ExtendedSegmentState> extendedGraphChromosomeState(
    RaMesh::RaMeshMultiGenomeGraph& graph,
    const SpeciesName& species,
    const ChrName& chromosome) {
    std::unordered_map<uint64_t, size_t> block_order;
    for (size_t index = 0; index < graph.blocks.size(); ++index) {
        if (auto block = graph.blocks[index].lock()) {
            block_order.emplace(block->block_id, index);
        }
    }

    std::vector<ExtendedSegmentState> state;
    auto& end = graph.species_graphs.at(species).chr2end.at(chromosome);
    for (auto segment = end.head->primary_path.next.load();
         segment && !segment->isTail();
         segment = segment->primary_path.next.load()) {
        assert(segment->parent_block);
        state.emplace_back(
            segment->start,
            segment->length,
            segment->strand,
            cigarToString(segment->cigar),
            segment->left_extend,
            segment->right_extend,
            block_order.at(segment->parent_block->block_id));
    }
    return state;
}

using BlockOccurrenceState =
    std::tuple<SpeciesName, ChrName, uint_t, uint_t, Strand,
               std::string, bool, bool>;

std::vector<std::vector<BlockOccurrenceState>> graphBlockState(
    RaMesh::RaMeshMultiGenomeGraph& graph) {
    std::vector<std::vector<BlockOccurrenceState>> result;
    result.reserve(graph.blocks.size());
    for (const auto& weak_block : graph.blocks) {
        auto block = weak_block.lock();
        assert(block);
        std::vector<BlockOccurrenceState> occurrences;
        occurrences.reserve(block->anchors.size());
        for (const auto& [key, segment] : block->anchors) {
            assert(segment);
            occurrences.emplace_back(
                key.first, key.second, segment->start, segment->length,
                segment->strand, cigarToString(segment->cigar),
                segment->left_extend, segment->right_extend);
        }
        std::sort(occurrences.begin(), occurrences.end());
        result.push_back(std::move(occurrences));
    }
    return result;
}

void legacyAlignIntervalForTesting(
    RaMesh::GenomeEnd& end,
    const SpeciesName& reference_name,
    const SpeciesName& query_name,
    const ChrName& query_chromosome,
    RaMesh::SegPtr current,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    bool is_left_extend,
    int_t zdrop) {
    if (!current || current == end.head || current == end.tail ||
        (current->left_extend && current->right_extend)) {
        return;
    }

    const auto fetch_sequence = [](
        const SeqPro::ManagerVariant& manager,
        const ChrName& chromosome,
        Coord_t start,
        Coord_t length) {
        return std::visit(
            [&](const auto& pointer) {
                using Pointer = std::decay_t<decltype(pointer)>;
                if constexpr (std::is_same_v<
                                  Pointer,
                                  std::unique_ptr<SeqPro::SequenceManager>>) {
                    return pointer->getSubSequence(
                        chromosome, start, length);
                } else {
                    return pointer->getOriginalManager().getSubSequence(
                        chromosome, start, length);
                }
            },
            manager);
    };
    const auto chromosome_length = [](
        const SeqPro::ManagerVariant& manager,
        const ChrName& chromosome) {
        return std::visit(
            [&](const auto& pointer) {
                using Pointer = std::decay_t<decltype(pointer)>;
                if constexpr (std::is_same_v<
                                  Pointer,
                                  std::unique_ptr<SeqPro::SequenceManager>>) {
                    return pointer->getSequenceLength(chromosome);
                } else {
                    return pointer->getOriginalManager().getSequenceLength(
                        chromosome);
                }
            },
            manager);
    };

    if (is_left_extend && !current->left_extend) {
        int_t query_start = 0;
        int_t query_length = 0;
        const Strand strand = current->strand;
        if (strand == FORWARD) {
            const auto query_left = current->primary_path.prev.load();
            query_start = !query_left->isHead()
                ? query_left->start + query_left->length
                : 0;
            query_length =
                static_cast<int_t>(current->start) - query_start;
        } else {
            const auto query_right = current->primary_path.next.load();
            query_start = current->start + current->length;
            query_length = !query_right->isTail()
                ? static_cast<int_t>(query_right->start) - query_start
                : static_cast<int_t>(chromosome_length(
                      *managers.at(query_name), query_chromosome)) -
                      query_start;
        }

        const auto block = current->parent_block;
        const auto reference_occurrence = block->anchors.find(
            {reference_name, block->ref_chr});
        assert(reference_occurrence != block->anchors.end());
        const auto reference_current = reference_occurrence->second;
        auto reference_left =
            reference_current->primary_path.prev.load();
        bool found = false;
        while (true) {
            if (reference_left->isHead()) break;
            if (reference_left->left_extend &&
                reference_left->right_extend) {
                break;
            }
            for (const auto& [key, segment] :
                 reference_left->parent_block->anchors) {
                (void)segment;
                if (key.first == query_name) {
                    found = true;
                    break;
                }
            }
            if (found) break;
            reference_left = reference_left->primary_path.prev.load();
        }

        const int_t reference_start = !reference_left->isHead()
            ? reference_left->start + reference_left->length
            : 0;
        const int_t reference_length =
            static_cast<int_t>(reference_current->start) -
            reference_start;
        if (found && query_length > 0 && reference_length > 0) {
            if (query_length > 10000 || reference_length > 10000) {
                return;
            }
            std::string query_sequence = fetch_sequence(
                *managers.at(query_name), query_chromosome,
                query_start, query_length);
            std::string reference_sequence = fetch_sequence(
                *managers.at(reference_name), block->ref_chr,
                reference_start, reference_length);
            std::reverse(
                reference_sequence.begin(), reference_sequence.end());
            if (strand == FORWARD) {
                std::reverse(
                    query_sequence.begin(), query_sequence.end());
            } else {
                baseComplement(query_sequence);
            }

            Cigar_t result = extendAlignKSW2(
                reference_sequence, query_sequence, zdrop);
            if (alignmentCigarPreferredToUnaligned(
                    reference_sequence, query_sequence, result)) {
                std::reverse(result.begin(), result.end());
                const AlignCount count = countAlignedBases(result);
                if (strand == FORWARD) {
                    current->start -= count.query_bases;
                }
                current->length += count.query_bases;
                reference_current->start -= count.ref_bases;
                reference_current->length += count.ref_bases;
                prependCigar(current->cigar, result);
            }
        }
    }

    if (!is_left_extend && !current->right_extend) {
        int_t query_start = 0;
        int_t query_length = 0;
        const Strand strand = current->strand;
        if (strand == FORWARD) {
            const auto query_right = current->primary_path.next.load();
            query_start = current->start + current->length;
            query_length = !query_right->isTail()
                ? static_cast<int_t>(query_right->start) - query_start
                : static_cast<int_t>(chromosome_length(
                      *managers.at(query_name), query_chromosome)) -
                      query_start;
        } else {
            const auto query_left = current->primary_path.prev.load();
            query_start = !query_left->isHead()
                ? query_left->start + query_left->length
                : 0;
            query_length =
                static_cast<int_t>(current->start) - query_start;
        }

        const auto block = current->parent_block;
        const auto reference_occurrence = block->anchors.find(
            {reference_name, block->ref_chr});
        assert(reference_occurrence != block->anchors.end());
        const auto reference_current = reference_occurrence->second;
        auto reference_right =
            reference_current->primary_path.next.load();
        bool found = false;
        while (true) {
            if (reference_right->isTail()) break;
            if (reference_right->left_extend &&
                reference_right->right_extend) {
                break;
            }
            for (const auto& [key, segment] :
                 reference_right->parent_block->anchors) {
                (void)segment;
                if (key.first == query_name) {
                    found = true;
                    break;
                }
            }
            if (found) break;
            reference_right = reference_right->primary_path.next.load();
        }

        const int_t reference_start =
            reference_current->start + reference_current->length;
        const int_t reference_length = !reference_right->isTail()
            ? static_cast<int_t>(reference_right->start) -
                  reference_start
            : static_cast<int_t>(chromosome_length(
                  *managers.at(reference_name), block->ref_chr)) -
                  reference_start;
        if (found && query_length > 0 && reference_length > 0) {
            if (query_length > 10000 || reference_length > 10000) {
                return;
            }
            std::string query_sequence = fetch_sequence(
                *managers.at(query_name), query_chromosome,
                query_start, query_length);
            std::string reference_sequence = fetch_sequence(
                *managers.at(reference_name), block->ref_chr,
                reference_start, reference_length);
            if (strand == REVERSE) {
                reverseComplement(query_sequence);
            }

            Cigar_t result = extendAlignKSW2(
                reference_sequence, query_sequence, zdrop);
            if (alignmentCigarPreferredToUnaligned(
                    reference_sequence, query_sequence, result)) {
                const AlignCount count = countAlignedBases(result);
                if (strand == REVERSE) {
                    current->start -= count.query_bases;
                }
                current->length += count.query_bases;
                reference_current->length += count.ref_bases;
                appendCigar(current->cigar, result);
            }
        }
    }
}

void runLegacyExtendRefNodes(
    RaMesh::RaMeshMultiGenomeGraph& graph,
    const SpeciesName& reference_name,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    int_t zdrop) {
    for (auto& [query_name, query_graph] : graph.species_graphs) {
        if (query_name == reference_name) continue;
        for (auto& [chromosome, end] : query_graph.chr2end) {
            auto segment = end.head;
            while (segment) {
                if (segment != end.head && segment != end.tail &&
                    !(segment->left_extend && segment->right_extend)) {
                    legacyAlignIntervalForTesting(
                        end,
                        reference_name, query_name, chromosome, segment,
                        managers, false, zdrop);
                }
                segment = segment->primary_path.next.load();
            }
        }
        for (auto& [chromosome, end] : query_graph.chr2end) {
            auto segment = end.head;
            while (segment) {
                if (segment != end.head && segment != end.tail &&
                    !(segment->left_extend && segment->right_extend)) {
                    legacyAlignIntervalForTesting(
                        end,
                        reference_name, query_name, chromosome, segment,
                        managers, true, zdrop);
                }
                segment = segment->primary_path.next.load();
            }
        }
    }
    for (auto& [chromosome, end] :
         graph.species_graphs.at(reference_name).chr2end) {
        (void)chromosome;
        end.resortSegments();
    }
}

using SharedManagerMap =
    std::map<SpeciesName, SeqPro::SharedManagerVariant>;

SharedManagerMap makeExtensionManagers(
    const std::filesystem::path& reference_path,
    const std::filesystem::path& query_path,
    const std::filesystem::path& other_path) {
    SharedManagerMap managers;
    managers.emplace(
        "ref",
        std::make_shared<SeqPro::ManagerVariant>(
            std::make_unique<SeqPro::SequenceManager>(reference_path)));
    managers.emplace(
        "qry",
        std::make_shared<SeqPro::ManagerVariant>(
            std::make_unique<SeqPro::SequenceManager>(query_path)));
    managers.emplace(
        "other",
        std::make_shared<SeqPro::ManagerVariant>(
            std::make_unique<SeqPro::SequenceManager>(other_path)));
    return managers;
}

void insertExtensionFixture(
    RaMesh::RaMeshMultiGenomeGraph& graph,
    SharedManagerMap& managers) {
    struct FixtureAnchor {
        SpeciesName query;
        Anchor anchor;
    };

    const auto make_anchor = [](
        uint32_t chromosome_id,
        uint32_t reference_start,
        uint32_t query_start,
        Strand strand) {
        constexpr uint32_t length = 20;
        return Anchor(
            chromosome_id, reference_start, length,
            chromosome_id, query_start, length, strand,
            length, length, Cigar_t{cigarToInt('M', length)});
    };

    std::vector<FixtureAnchor> anchors;
    anchors.push_back({"qry", make_anchor(0, 100, 100, FORWARD)});
    anchors.push_back({"other", make_anchor(0, 115, 116, FORWARD)});
    anchors.push_back({"qry", make_anchor(0, 131, 132, FORWARD)});
    anchors.push_back({"other", make_anchor(0, 5000, 5001, REVERSE)});
    anchors.push_back({"qry", make_anchor(0, 10150, 10152, FORWARD)});
    anchors.push_back({"other", make_anchor(0, 15000, 15005, FORWARD)});
    anchors.push_back({"qry", make_anchor(0, 20170, 20172, REVERSE)});
    anchors.push_back({"qry", make_anchor(0, 30300, 30304, FORWARD)});
    anchors.push_back({"other", make_anchor(0, 35000, 35000, REVERSE)});
    anchors.push_back({"qry", make_anchor(1, 50, 60, FORWARD)});
    anchors.push_back({"other", make_anchor(1, 4000, 4010, FORWARD)});
    anchors.push_back({"qry", make_anchor(1, 10070, 10080, REVERSE)});

    for (const FixtureAnchor& fixture : anchors) {
        graph.insertAnchorIntoGraph(
            *managers.at("ref"), *managers.at(fixture.query),
            "ref", fixture.query, fixture.anchor, false);
    }

    auto& reference_end =
        graph.species_graphs.at("ref").chr2end.at("chr1");
    RaMesh::SegPtr barrier;
    for (auto segment = reference_end.head->primary_path.next.load();
         segment && !segment->isTail();
         segment = segment->primary_path.next.load()) {
        if (segment->start == 10150) {
            barrier = segment;
            break;
        }
    }
    assert(barrier);
    barrier->left_extend = true;
    barrier->right_extend = true;

    auto duplicate = RaMesh::Segment::create(
        15000, 20, REVERSE, Cigar_t{cigarToInt('M', 20)},
        RaMesh::AlignRole::PRIMARY, RaMesh::SegmentRole::SEGMENT,
        barrier->parent_block);
    barrier->parent_block->anchors.emplace(
        RaMesh::SpeciesChrPair{"qry", "chr1"}, duplicate);
    graph.species_graphs.at("qry")
        .chr2end.at("chr1")
        .insertSegment(duplicate);
}

void testExtendRefNodesEquivalence(
    const std::filesystem::path& directory) {
    const auto reference_path = directory / "extension-ref.fa";
    const auto query_path = directory / "extension-qry.fa";
    const auto other_path = directory / "extension-other.fa";
    writeFile(
        reference_path,
        ">chr1\n" + std::string(40000, 'A') +
        "\n>chr2\n" + std::string(20000, 'A') + "\n");
    writeFile(
        query_path,
        ">chr1\n" + std::string(40000, 'A') +
        "\n>chr2\n" + std::string(20000, 'C') + "\n");
    writeFile(
        other_path,
        ">chr1\n" + std::string(40000, 'A') +
        "\n>chr2\n" + std::string(20000, 'A') + "\n");

    SharedManagerMap managers = makeExtensionManagers(
        reference_path, query_path, other_path);
    RaMesh::RaMeshMultiGenomeGraph legacy_graph(managers);
    RaMesh::RaMeshMultiGenomeGraph optimized_graph(managers);
    insertExtensionFixture(legacy_graph, managers);
    insertExtensionFixture(optimized_graph, managers);

    runLegacyExtendRefNodes(legacy_graph, "ref", managers, 400);
    optimized_graph.extendRefNodes("ref", managers, 400);

    for (const SpeciesName& species :
         std::vector<SpeciesName>{"ref", "qry", "other"}) {
        for (const ChrName& chromosome :
             std::vector<ChrName>{"chr1", "chr2"}) {
            assert(extendedGraphChromosomeState(
                       legacy_graph, species, chromosome) ==
                   extendedGraphChromosomeState(
                       optimized_graph, species, chromosome));
        }
    }
    assert(graphBlockState(legacy_graph) ==
           graphBlockState(optimized_graph));

    for (const int32_t gap : {-1, 0, 1, 9999, 10000, 10001}) {
        RaMesh::RaMeshMultiGenomeGraph boundary_legacy(managers);
        RaMesh::RaMeshMultiGenomeGraph boundary_optimized(managers);
        constexpr uint32_t first_start = 100;
        constexpr uint32_t anchor_length = 20;
        const uint32_t second_start = static_cast<uint32_t>(
            static_cast<int64_t>(first_start) + anchor_length + gap);
        const Anchor first_anchor(
            0, first_start, anchor_length,
            0, first_start, anchor_length, FORWARD,
            anchor_length, anchor_length,
            Cigar_t{cigarToInt('M', anchor_length)});
        const Anchor second_anchor(
            0, second_start, anchor_length,
            0, second_start, anchor_length, FORWARD,
            anchor_length, anchor_length,
            Cigar_t{cigarToInt('M', anchor_length)});

        for (RaMesh::RaMeshMultiGenomeGraph* graph :
             {&boundary_legacy, &boundary_optimized}) {
            graph->insertAnchorIntoGraph(
                *managers.at("ref"), *managers.at("qry"),
                "ref", "qry", first_anchor, false);
            graph->insertAnchorIntoGraph(
                *managers.at("ref"), *managers.at("qry"),
                "ref", "qry", second_anchor, false);
        }

        auto& legacy_end =
            boundary_legacy.species_graphs.at("qry").chr2end.at("chr1");
        auto& optimized_end =
            boundary_optimized.species_graphs.at("qry").chr2end.at("chr1");
        auto legacy_first = legacy_end.head->primary_path.next.load();
        auto optimized_first =
            optimized_end.head->primary_path.next.load();
        legacyAlignIntervalForTesting(
            legacy_end, "ref", "qry", "chr1", legacy_first,
            managers, false, 400);
        optimized_end.alignInterval(
            "ref", "qry", "chr1", optimized_first,
            managers, false, 400);

        assert(extendedGraphChromosomeState(
                   boundary_legacy, "ref", "chr1") ==
               extendedGraphChromosomeState(
                   boundary_optimized, "ref", "chr1"));
        assert(extendedGraphChromosomeState(
                   boundary_legacy, "qry", "chr1") ==
               extendedGraphChromosomeState(
                   boundary_optimized, "qry", "chr1"));
        assert(graphBlockState(boundary_legacy) ==
               graphBlockState(boundary_optimized));

        const bool should_extend = gap >= 1 && gap <= 10000;
        assert((optimized_first->length > anchor_length) ==
               should_extend);
    }
}

void testResortSegmentsShortcut() {
    RaMesh::GenomeEnd increasing;
    auto first = RaMesh::Segment::create(
        10, 2, FORWARD, {}, RaMesh::AlignRole::PRIMARY,
        RaMesh::SegmentRole::SEGMENT, nullptr);
    auto second = RaMesh::Segment::create(
        20, 2, FORWARD, {}, RaMesh::AlignRole::PRIMARY,
        RaMesh::SegmentRole::SEGMENT, nullptr);
    auto third = RaMesh::Segment::create(
        30, 2, FORWARD, {}, RaMesh::AlignRole::PRIMARY,
        RaMesh::SegmentRole::SEGMENT, nullptr);
    RaMesh::Segment::linkChain(
        {increasing.head, first, second, third, increasing.tail});
    increasing.resortSegments();
    assert(increasing.head->primary_path.next.load() == first);
    assert(first->primary_path.next.load() == second);
    assert(second->primary_path.next.load() == third);

    RaMesh::GenomeEnd unordered;
    auto high = RaMesh::Segment::create(
        30, 2, FORWARD, {}, RaMesh::AlignRole::PRIMARY,
        RaMesh::SegmentRole::SEGMENT, nullptr);
    auto low = RaMesh::Segment::create(
        10, 2, FORWARD, {}, RaMesh::AlignRole::PRIMARY,
        RaMesh::SegmentRole::SEGMENT, nullptr);
    auto middle = RaMesh::Segment::create(
        20, 2, FORWARD, {}, RaMesh::AlignRole::PRIMARY,
        RaMesh::SegmentRole::SEGMENT, nullptr);
    RaMesh::Segment::linkChain(
        {unordered.head, high, low, middle, unordered.tail});
    unordered.resortSegments();
    assert(unordered.head->primary_path.next.load() == low);
    assert(low->primary_path.next.load() == middle);
    assert(middle->primary_path.next.load() == high);

    RaMesh::GenomeEnd equal_starts;
    auto equal_a = RaMesh::Segment::create(
        10, 1, FORWARD, {}, RaMesh::AlignRole::PRIMARY,
        RaMesh::SegmentRole::SEGMENT, nullptr);
    auto equal_b = RaMesh::Segment::create(
        10, 2, FORWARD, {}, RaMesh::AlignRole::PRIMARY,
        RaMesh::SegmentRole::SEGMENT, nullptr);
    auto later = RaMesh::Segment::create(
        20, 1, FORWARD, {}, RaMesh::AlignRole::PRIMARY,
        RaMesh::SegmentRole::SEGMENT, nullptr);
    std::vector<RaMesh::SegPtr> expected{equal_a, equal_b, later};
    std::sort(
        expected.begin(), expected.end(),
        [](const RaMesh::SegPtr& lhs, const RaMesh::SegPtr& rhs) {
            return lhs->start < rhs->start;
        });
    RaMesh::Segment::linkChain(
        {equal_starts.head, equal_a, equal_b, later, equal_starts.tail});
    equal_starts.resortSegments();
    auto current = equal_starts.head->primary_path.next.load();
    for (const auto& segment : expected) {
        assert(current == segment);
        current = current->primary_path.next.load();
    }
    assert(current == equal_starts.tail);
}

struct LegacyMergeCigarSplit {
    bool has_prefix = false;
    bool has_suffix = false;
    RaMesh::detail::MergeCigarPiece prefix;
    RaMesh::detail::MergeCigarPiece overlap;
    RaMesh::detail::MergeCigarPiece suffix;
};

LegacyMergeCigarSplit legacySplitCigarForMerge(
    const Cigar_t& source,
    bool has_prefix,
    uint32_t prefix_reference_length,
    uint32_t overlap_reference_length,
    bool has_suffix,
    uint32_t suffix_reference_length) {
    LegacyMergeCigarSplit result;
    result.has_prefix = has_prefix;
    result.has_suffix = has_suffix;
    std::string remainder = cigarToString(source);

    const auto take_piece = [&](uint32_t reference_length,
                                RaMesh::detail::MergeCigarPiece& piece) {
        auto [cigar, next_remainder] =
            splitCigarMixed(remainder, reference_length);
        remainder = std::move(next_remainder);
        piece.cigar = std::move(cigar);
        piece.reference_length = countRefLength(piece.cigar);
        piece.query_length = countNonDeletionOperations(piece.cigar);
        if (piece.reference_length != reference_length) {
            throw std::logic_error(
                "legacy merge CIGAR ends before the requested span");
        }
    };

    if (has_prefix) {
        take_piece(prefix_reference_length, result.prefix);
    }
    take_piece(overlap_reference_length, result.overlap);
    if (has_suffix) {
        take_piece(suffix_reference_length, result.suffix);
    }

    Cigar_t terminal;
    parseCigarString(remainder, terminal);
    if (countRefLength(terminal) != 0) {
        throw std::logic_error(
            "legacy merge left reference-consuming operations");
    }
    auto& terminal_piece = has_suffix ? result.suffix : result.overlap;
    terminal_piece.query_length += countQryLength(terminal);
    appendCigar(terminal_piece.cigar, terminal);
    return result;
}

void assertMergeCigarPieceEqual(
    const RaMesh::detail::MergeCigarPiece& expected,
    const RaMesh::detail::MergeCigarPiece& observed) {
    assert(expected.cigar == observed.cigar);
    assert(expected.reference_length == observed.reference_length);
    assert(expected.query_length == observed.query_length);
}

void compareMergeCigarSplit(
    const Cigar_t& source,
    uint32_t prefix_reference_length,
    uint32_t overlap_reference_length,
    uint32_t suffix_reference_length) {
    const bool has_prefix = prefix_reference_length != 0;
    const bool has_suffix = suffix_reference_length != 0;
    bool legacy_threw = false;
    bool optimized_threw = false;
    LegacyMergeCigarSplit legacy;
    RaMesh::detail::MergeCigarSplit optimized;
    try {
        legacy = legacySplitCigarForMerge(
            source, has_prefix, prefix_reference_length,
            overlap_reference_length, has_suffix,
            suffix_reference_length);
    } catch (const std::logic_error&) {
        legacy_threw = true;
    }
    try {
        optimized = RaMesh::detail::splitCigarForOverlapMerge(
            source, has_prefix, prefix_reference_length,
            overlap_reference_length, has_suffix,
            suffix_reference_length);
    } catch (const std::logic_error&) {
        optimized_threw = true;
    }
    assert(legacy_threw == optimized_threw);
    if (legacy_threw) return;
    assert(optimized.has_prefix == has_prefix);
    assert(optimized.has_suffix == has_suffix);
    assert(optimized.source_units_scanned == source.size());
    assertMergeCigarPieceEqual(legacy.prefix, optimized.prefix);
    assertMergeCigarPieceEqual(legacy.overlap, optimized.overlap);
    assertMergeCigarPieceEqual(legacy.suffix, optimized.suffix);
}

void testMergeCigarSplitter() {
    compareMergeCigarSplit({}, 0, 1, 0);
    compareMergeCigarSplit({cigarToInt('M', 100)}, 25, 50, 25);
    compareMergeCigarSplit(
        {cigarToInt('I', 7), cigarToInt('M', 100)}, 25, 50, 25);
    compareMergeCigarSplit(
        {cigarToInt('M', 25), cigarToInt('I', 7),
         cigarToInt('D', 10), cigarToInt('=', 30),
         cigarToInt('X', 35), cigarToInt('I', 9)},
        25, 40, 35);
    compareMergeCigarSplit(
        {cigarToInt('M', 100), cigarToInt('I', 10)}, 0, 50, 50);
    compareMergeCigarSplit(
        {cigarToInt('I', 3), cigarToInt('I', 4)}, 0, 0, 0);

    std::mt19937 generator(0x6d657267u);
    constexpr std::array<char, 5> operations{'M', 'I', 'D', '=', 'X'};
    for (size_t fixture = 0; fixture < 100000; ++fixture) {
        Cigar_t source;
        uint32_t reference_length = 0;
        const size_t operation_count = 1 + generator() % 30;
        source.reserve(operation_count + 1);
        for (size_t index = 0; index < operation_count; ++index) {
            const char operation = operations[generator() % operations.size()];
            const uint32_t length = 1 + generator() % 40;
            source.push_back(cigarToInt(operation, length));
            if (operation != 'I') reference_length += length;
        }
        if (reference_length == 0) {
            source.push_back(cigarToInt('M', 1));
            reference_length = 1;
        }
        if (generator() % 2 == 0) {
            source.push_back(cigarToInt('I', 1 + generator() % 20));
        }

        const uint32_t prefix = generator() % reference_length;
        const uint32_t remaining = reference_length - prefix;
        const uint32_t overlap = 1 + generator() % remaining;
        const uint32_t suffix = remaining - overlap;
        compareMergeCigarSplit(source, prefix, overlap, suffix);
    }
}

std::map<SpeciesName, SeqPro::ManagerVariant> makeMergeManagers(
    const std::filesystem::path& directory,
    const std::vector<SpeciesName>& species) {
    std::map<SpeciesName, SeqPro::ManagerVariant> managers;
    for (const SpeciesName& name : species) {
        const auto fasta = directory / ("merge-" + name + ".fa");
        writeFile(fasta, ">chr1\n" + std::string(10000, 'A') + "\n");
        managers.emplace(
            name, std::make_unique<SeqPro::SequenceManager>(fasta));
    }
    return managers;
}

RaMesh::SegPtr addMergeSegment(
    const RaMesh::BlockPtr& block,
    const RaMesh::SpeciesChrPair& key,
    uint_t start,
    uint_t length,
    Strand strand,
    Cigar_t cigar) {
    auto segment = RaMesh::Segment::create(
        start, length, strand, std::move(cigar),
        RaMesh::AlignRole::PRIMARY, RaMesh::SegmentRole::SEGMENT, block);
    block->anchors.emplace(key, segment);
    return segment;
}

void linkMergePath(
    RaMesh::GenomeEnd& end,
    std::initializer_list<RaMesh::SegPtr> segments) {
    std::vector<RaMesh::SegPtr> chain;
    chain.reserve(segments.size() + 2);
    chain.push_back(end.head);
    chain.insert(chain.end(), segments.begin(), segments.end());
    chain.push_back(end.tail);
    RaMesh::Segment::linkChain(chain);
    end.sample_vec.assign(1, end.head);
    for (const auto& segment : segments) end.setToSampling(segment);
}

void assertSamplingReferencesLivePath(RaMesh::GenomeEnd& end) {
    std::unordered_set<const RaMesh::Segment*> live;
    auto segment = end.head;
    while (segment) {
        assert(live.emplace(segment.get()).second);
        if (segment == end.tail) break;
        segment = segment->primary_path.next.load();
    }
    assert(segment == end.tail);
    for (const auto& sampled : end.sample_vec) {
        assert(!sampled || live.count(sampled.get()) == 1);
    }
}

void testMergeGraphRetirementAndOrdering(
    const std::filesystem::path& directory) {
    auto managers = makeMergeManagers(
        directory, {"ref", "left", "right"});
    RaMesh::RaMeshMultiGenomeGraph graph(managers);

    auto left_block = RaMesh::Block::createEmpty("ref", "chr1", 2);
    auto right_block = RaMesh::Block::createEmpty("ref", "chr1", 2);
    auto left_reference = addMergeSegment(
        left_block, {"ref", "chr1"}, 100, 100, FORWARD, {});
    auto right_reference = addMergeSegment(
        right_block, {"ref", "chr1"}, 150, 100, FORWARD, {});
    auto left_query = addMergeSegment(
        left_block, {"left", "chr1"}, 1000, 105, FORWARD,
        {cigarToInt('M', 50), cigarToInt('I', 5),
         cigarToInt('M', 50)});
    auto right_query = addMergeSegment(
        right_block, {"right", "chr1"}, 2000, 105, REVERSE,
        {cigarToInt('M', 50), cigarToInt('I', 5),
         cigarToInt('M', 50)});

    linkMergePath(
        graph.species_graphs.at("ref").chr2end.at("chr1"),
        {left_reference, right_reference});
    linkMergePath(
        graph.species_graphs.at("left").chr2end.at("chr1"),
        {left_query});
    linkMergePath(
        graph.species_graphs.at("right").chr2end.at("chr1"),
        {right_query});
    graph.blocks = {left_block, right_block};

    const std::weak_ptr<RaMesh::Block> weak_left_block = left_block;
    const std::weak_ptr<RaMesh::Block> weak_right_block = right_block;
    const std::weak_ptr<RaMesh::Segment> weak_left_reference = left_reference;
    const std::weak_ptr<RaMesh::Segment> weak_right_reference = right_reference;
    const std::weak_ptr<RaMesh::Segment> weak_left_query = left_query;
    const std::weak_ptr<RaMesh::Segment> weak_right_query = right_query;
    left_block.reset();
    right_block.reset();
    left_reference.reset();
    right_reference.reset();
    left_query.reset();
    right_query.reset();

    graph.mergeMultipleGraphs("ref", 32);

    assert(weak_left_block.expired());
    assert(weak_right_block.expired());
    assert(weak_left_reference.expired());
    assert(weak_right_reference.expired());
    assert(weak_left_query.expired());
    assert(weak_right_query.expired());
    assert(graph.blocks.size() == 3);
    for (const auto& weak_block : graph.blocks) assert(!weak_block.expired());

    const auto reference_state = graphChromosomeState(
        graph, "ref", "chr1");
    assert(reference_state.size() == 3);
    assert(std::get<0>(reference_state[0]) == 100);
    assert(std::get<1>(reference_state[0]) == 50);
    assert(std::get<0>(reference_state[1]) == 150);
    assert(std::get<1>(reference_state[1]) == 50);
    assert(std::get<0>(reference_state[2]) == 200);
    assert(std::get<1>(reference_state[2]) == 50);

    const auto left_state = graphChromosomeState(
        graph, "left", "chr1");
    assert(left_state.size() == 2);
    assert(std::get<0>(left_state[0]) == 1000);
    assert(std::get<3>(left_state[0]) == "50M");
    assert(std::get<0>(left_state[1]) == 1050);
    assert(std::get<3>(left_state[1]) == "5I50M");

    const auto right_state = graphChromosomeState(
        graph, "right", "chr1");
    assert(right_state.size() == 2);
    assert(std::get<0>(right_state[0]) == 2000);
    assert(std::get<3>(right_state[0]) == "5I50M");
    assert(std::get<0>(right_state[1]) == 2055);
    assert(std::get<3>(right_state[1]) == "50M");

    for (const char* species : {"ref", "left", "right"}) {
        assertSamplingReferencesLivePath(
            graph.species_graphs.at(species).chr2end.at("chr1"));
    }
    assert(graph.verifyGraphCorrectness(false));
    assert(graph.compactBlockPool() == 0);
}

void testMergeParticipantConflict(const std::filesystem::path& directory) {
    auto managers = makeMergeManagers(directory, {"ref", "query"});
    RaMesh::RaMeshMultiGenomeGraph graph(managers);
    auto left_block = RaMesh::Block::createEmpty("ref", "chr1", 2);
    auto right_block = RaMesh::Block::createEmpty("ref", "chr1", 2);
    auto left_reference = addMergeSegment(
        left_block, {"ref", "chr1"}, 0, 100, FORWARD, {});
    auto right_reference = addMergeSegment(
        right_block, {"ref", "chr1"}, 50, 100, FORWARD, {});
    auto left_query = addMergeSegment(
        left_block, {"query", "chr1"}, 0, 100, FORWARD,
        {cigarToInt('M', 100)});
    auto right_query = addMergeSegment(
        right_block, {"query", "chr1"}, 200, 100, FORWARD,
        {cigarToInt('M', 100)});
    linkMergePath(
        graph.species_graphs.at("ref").chr2end.at("chr1"),
        {left_reference, right_reference});
    linkMergePath(
        graph.species_graphs.at("query").chr2end.at("chr1"),
        {left_query, right_query});
    graph.blocks = {left_block, right_block};
    graph.mergeMultipleGraphs("ref", 1);
    assert(graph.blocks.size() == 2);
    assert(graph.blocks[0].lock() == left_block);
    assert(graph.blocks[1].lock() == right_block);
    assert(graphChromosomeState(graph, "ref", "chr1").size() == 2);
    assert(graphChromosomeState(graph, "query", "chr1").size() == 2);
}

void testMergeRetirementStress(const std::filesystem::path& directory) {
    constexpr size_t merge_count = 100000;
    auto managers = makeMergeManagers(directory, {"ref"});
    RaMesh::RaMeshMultiGenomeGraph graph(managers);
    auto& end = graph.species_graphs.at("ref").chr2end.at("chr1");
    std::vector<std::weak_ptr<RaMesh::Block>> source_blocks;
    std::vector<std::weak_ptr<RaMesh::Segment>> source_segments;
    source_blocks.reserve(merge_count + 1);
    source_segments.reserve(merge_count + 1);
    graph.blocks.reserve(merge_count + 1);
    RaMesh::SegPtr previous = end.head;
    for (size_t index = 0; index <= merge_count; ++index) {
        auto block = RaMesh::Block::createEmpty("ref", "chr1", 1);
        auto segment = addMergeSegment(
            block, {"ref", "chr1"}, 1000, 100, FORWARD, {});
        previous->primary_path.next.store(segment);
        segment->primary_path.prev.store(previous);
        previous = segment;
        graph.blocks.emplace_back(block);
        source_blocks.emplace_back(block);
        source_segments.emplace_back(segment);
    }
    previous->primary_path.next.store(end.tail);
    end.tail->primary_path.prev.store(previous);
    previous.reset();
    end.sample_vec.assign(1, end.head);

    graph.mergeMultipleGraphs("ref", 32);
    for (const auto& block : source_blocks) assert(block.expired());
    for (const auto& segment : source_segments) assert(segment.expired());
    assert(graph.blocks.size() == 1);
    assert(!graph.blocks.front().expired());
    assert(graph.verifyGraphCorrectness(false));
    assert(graph.compactBlockPool() == 0);
}

void testBatchInsertion(const std::filesystem::path& directory) {
    const auto ref_path = directory / "ref.fa";
    const auto qry_path = directory / "qry.fa";
    writeFile(ref_path, ">chr1\n" + std::string(120000, 'A') + "\n");
    writeFile(qry_path, ">chr1\n" + std::string(120000, 'A') + "\n");
    std::map<SpeciesName, SeqPro::ManagerVariant> managers;
    managers.emplace("ref", std::make_unique<SeqPro::SequenceManager>(ref_path));
    managers.emplace("qry", std::make_unique<SeqPro::SequenceManager>(qry_path));
    RaMesh::RaMeshMultiGenomeGraph legacy_graph(managers);
    RaMesh::RaMeshMultiGenomeGraph batch_graph(managers);

    std::mt19937 generator(0x12345678u);
    AnchorVec legacy_anchors;
    AnchorVec batch_anchors;
    size_t graph_anchor_count =
        std::getenv("RAMAX_RUN_MILLION_SEGMENT_TEST") ? 1000000 : 100000;
    if (const char* configured = std::getenv("RAMAX_GRAPH_TEST_ANCHORS")) {
        graph_anchor_count = std::stoull(configured);
    }
    legacy_anchors.reserve(graph_anchor_count);
    batch_anchors.reserve(graph_anchor_count);
    for (size_t index = 0; index < graph_anchor_count; ++index) {
        const uint32_t ref_start = (generator() % 10000) * 8;
        const uint32_t query_start = (generator() % 10000) * 8;
        const uint32_t length = 1 + generator() % 80;
        Anchor anchor(0, ref_start, length, 0, query_start, length,
            index % 7 == 0 ? REVERSE : FORWARD,
            length, length, Cigar_t{cigarToInt('M', length)});
        legacy_anchors.push_back(anchor);
        batch_anchors.push_back(std::move(anchor));
    }
    const auto legacy_start = std::chrono::steady_clock::now();
    for (const Anchor& anchor : legacy_anchors) {
        legacy_graph.insertAnchorIntoGraph(
            managers.at("ref"), managers.at("qry"),
            "ref", "qry", anchor, false);
    }
    std::vector<Anchor*> pointers;
    pointers.reserve(batch_anchors.size());
    for (Anchor& anchor : batch_anchors) pointers.push_back(&anchor);
    const auto batch_start = std::chrono::steady_clock::now();
    batch_graph.insertAnchorsIntoGraphBatch(
        managers.at("ref"), managers.at("qry"),
        "ref", "qry", pointers);
    const auto batch_finish = std::chrono::steady_clock::now();

    assert(graphChromosomeState(legacy_graph, "ref", "chr1") ==
           graphChromosomeState(batch_graph, "ref", "chr1"));
    assert(graphChromosomeState(legacy_graph, "qry", "chr1") ==
           graphChromosomeState(batch_graph, "qry", "chr1"));
    assert(legacy_graph.blocks.size() == batch_graph.blocks.size());
    for (const Anchor& anchor : batch_anchors) assert(anchor.cigar.empty());
    const double legacy_ms = std::chrono::duration<double, std::milli>(
        batch_start - legacy_start).count();
    const double batch_ms = std::chrono::duration<double, std::milli>(
        batch_finish - batch_start).count();
    std::cout << "Graph " << graph_anchor_count
              << "-anchor insertion speedup: "
              << legacy_ms / std::max(batch_ms, 0.001)
              << "x (legacy=" << legacy_ms
              << " ms, batch=" << batch_ms << " ms)\n";
}

}  // namespace

int main() {
    const auto directory = makeTestDirectory();
    if (std::getenv("RAMAX_RUN_MERGE_STRESS_ONLY")) {
        testMergeRetirementStress(directory);
        std::filesystem::remove_all(directory);
        return 0;
    }
    if (std::getenv("RAMAX_RUN_MERGE_ONLY")) {
        testMergeCigarSplitter();
        testMergeGraphRetirementAndOrdering(directory);
        testMergeParticipantConflict(directory);
        std::filesystem::remove_all(directory);
        return 0;
    }
    testSequenceViewsAndKswSummary(directory);
    testAnchorLinkGapBoundaries(directory);
    testDpEquivalence();
    testDpPerformanceSmoke();
    testExtendRefNodesEquivalence(directory);
    testResortSegmentsShortcut();
    testMergeCigarSplitter();
    testMergeGraphRetirementAndOrdering(directory);
    testMergeParticipantConflict(directory);
    testBatchInsertion(directory);
    std::filesystem::remove_all(directory);
    return 0;
}
