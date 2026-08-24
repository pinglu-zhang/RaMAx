#include "suffix_array_index.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <type_traits>

#include <omp.h>

#if defined(__SSE2__)
#include <emmintrin.h>
#endif

extern "C" {
#include "divsufsort.h"
#include "divsufsort64.h"
}

namespace {

// The SIMD comparison kernel, Kasai LCP construction, and generalized
// suffix-link interval derivation are adapted from sufkit (MIT).  They are
// kept local to RaMAx so the anchor and graph layers do not depend on sufkit
// types or on a second copy of libdivsufsort.
struct ByteComparison {
    int order{0};
    size_t lcp{0};
};

ByteComparison comparePatternBytes(const uint8_t* text, size_t text_size,
    const uint8_t* pattern, size_t pattern_size, size_t known_lcp = 0) {
    size_t index = std::min(known_lcp, pattern_size);
    const size_t scalar_end =
        std::min(pattern_size, index + static_cast<size_t>(8));
    while (index < scalar_end) {
        if (index >= text_size) return {-1, index};
        const uint8_t left = text[index];
        const uint8_t right = pattern[index];
        if (left != right) return {left < right ? -1 : 1, index};
        ++index;
    }

#if defined(__SSE2__)
    while (index + 16 <= pattern_size && index + 16 <= text_size) {
        const auto left = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(text + index));
        const auto right = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(pattern + index));
        const auto equal = _mm_cmpeq_epi8(left, right);
        const unsigned mask = static_cast<unsigned>(_mm_movemask_epi8(equal));
        if (mask != 0xFFFFU) {
            const size_t mismatch = static_cast<size_t>(
                __builtin_ctz(static_cast<unsigned>(~mask) & 0xFFFFU));
            const size_t position = index + mismatch;
            return {text[position] < pattern[position] ? -1 : 1, position};
        }
        index += 16;
    }
#endif

    while (index < pattern_size) {
        if (index >= text_size) return {-1, index};
        const uint8_t left = text[index];
        const uint8_t right = pattern[index];
        if (left != right) return {left < right ? -1 : 1, index};
        ++index;
    }
    return {0, index};
}

size_t longestCommonPrefixBytes(const uint8_t* left, const uint8_t* right,
    size_t length) {
    return comparePatternBytes(left, length, right, length).lcp;
}

constexpr uint_t chooseAccurateSearchAdvance(
    uint_t match_length, size_t accepted_region_count,
    uint_t threshold) noexcept {
    return threshold > 0 && accepted_region_count == 1 &&
            match_length > threshold
        ? match_length
        : 1;
}

static_assert(chooseAccurateSearchAdvance(10000, 1, 10000) == 1);
static_assert(chooseAccurateSearchAdvance(10001, 1, 10000) == 10001);
static_assert(chooseAccurateSearchAdvance(10001, 1, 0) == 1);

constexpr uint_t kSaPrefixLength = 8;
constexpr size_t kSaPrefixCount = 1ULL << (2 * kSaPrefixLength);

bool encodeSaPrefix(const char* sequence, uint_t length, size_t& code) {
    if (length < kSaPrefixLength) return false;
    code = 0;
    for (uint_t index = 0; index < kSaPrefixLength; ++index) {
        code <<= 2;
        switch (sequence[index]) {
        case 'A': break;
        case 'C': code |= 1; break;
        case 'G': code |= 2; break;
        case 'T': code |= 3; break;
        default: return false;
        }
    }
    return true;
}

bool isSaCanonicalBase(char base) noexcept {
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

constexpr std::array<char, 8> kIndexMagic{
    'R', 'M', 'S', 'A', 'I', 'D', 'X', '1'};
constexpr uint32_t kIndexVersion = 2;

template<class T>
void writePod(std::ostream& output, const T& value) {
    output.write(reinterpret_cast<const char*>(&value), sizeof(value));
    if (!output) throw std::runtime_error("Failed to write suffix-array index");
}

template<class T>
void readPod(std::istream& input, T& value) {
    input.read(reinterpret_cast<char*>(&value), sizeof(value));
    if (!input) throw std::runtime_error("Truncated suffix-array index");
}

template<class Coordinate>
void writeCoordinates(std::ostream& output,
    const std::vector<Coordinate>& values) {
    if (!values.empty()) {
        output.write(reinterpret_cast<const char*>(values.data()),
            static_cast<std::streamsize>(values.size() * sizeof(Coordinate)));
    }
    if (!output) throw std::runtime_error("Failed to write suffix-array coordinates");
}

template<class Coordinate>
void readCoordinates(std::istream& input, std::vector<Coordinate>& values,
    size_t count) {
    values.resize(count);
    if (!values.empty()) {
        input.read(reinterpret_cast<char*>(values.data()),
            static_cast<std::streamsize>(values.size() * sizeof(Coordinate)));
    }
    if (!input) throw std::runtime_error("Truncated suffix-array coordinates");
}

template<class Coordinate>
std::vector<Coordinate> buildLcp(const std::string& text,
    const std::vector<Coordinate>& suffix_array,
    const std::vector<Coordinate>& inverse_suffix_array) {
    std::vector<Coordinate> lcp(suffix_array.size(), 0);
    uint64_t common = 0;
    for (uint64_t suffix = 0; suffix < suffix_array.size(); ++suffix) {
        const uint64_t row = inverse_suffix_array[static_cast<size_t>(suffix)];
        if (row == 0) {
            common = 0;
            continue;
        }
        const uint64_t previous =
            suffix_array[static_cast<size_t>(row - 1)];
        if (suffix + common < text.size() && previous + common < text.size()) {
            const size_t remaining = static_cast<size_t>(std::min<uint64_t>(
                text.size() - suffix - common,
                text.size() - previous - common));
            common += longestCommonPrefixBytes(
                reinterpret_cast<const uint8_t*>(text.data()) + suffix + common,
                reinterpret_cast<const uint8_t*>(text.data()) + previous + common,
                remaining);
        }
        lcp[static_cast<size_t>(row)] = static_cast<Coordinate>(common);
        if (common > 0) --common;
    }
    return lcp;
}

template<class Coordinate>
void validateCoordinates(const std::vector<Coordinate>& suffix_array,
    const std::vector<Coordinate>& inverse_suffix_array,
    const std::vector<Coordinate>& lcp) {
    const uint64_t size = suffix_array.size();
    if (inverse_suffix_array.size() != size || lcp.size() != size) {
        throw std::runtime_error("Suffix-array component sizes do not match");
    }
    if (size != 0 && lcp.front() != 0) {
        throw std::runtime_error("Suffix-array LCP row zero is invalid");
    }
    for (uint64_t row = 0; row < size; ++row) {
        const uint64_t suffix = suffix_array[static_cast<size_t>(row)];
        if (suffix >= size ||
            inverse_suffix_array[static_cast<size_t>(suffix)] != row ||
            lcp[static_cast<size_t>(row)] > size) {
            throw std::runtime_error("Suffix-array index is inconsistent");
        }
    }
}

} // namespace

Suffix_Array_Index::Suffix_Array_Index(SpeciesName species_name,
    SeqPro::ManagerVariant& fasta_manager)
    : species_name(std::move(species_name)), fasta_manager(fasta_manager) {}

bool Suffix_Array_Index::buildIndex(FilePath output_path, bool fast_mode,
    uint_t thread_count) {
    (void)output_path;

    std::string forward_text = std::visit([](auto&& manager_ptr) -> std::string {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if (!manager_ptr) {
            throw std::runtime_error("Manager pointer is null inside variant");
        }
        if constexpr (std::is_same_v<PtrType,
                          std::unique_ptr<SeqPro::SequenceManager>>) {
            return manager_ptr->concatAllSequences('\1');
        } else if constexpr (std::is_same_v<PtrType,
                                 std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            return manager_ptr->concatAllSequencesSeparated('\1');
        } else {
            throw std::runtime_error("Unhandled manager type in variant");
        }
    }, fasta_manager);

    if (forward_text.empty()) return false;
    if (forward_text.size() > std::numeric_limits<uint_t>::max()) {
        throw std::runtime_error("Reference is too large for RaMAx coordinates");
    }
    const uint_t logical_size = static_cast<uint_t>(forward_text.size());
    std::string fasta_path;
    std::visit([&](auto&& manager_ptr) {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if constexpr (std::is_same_v<PtrType,
                          std::unique_ptr<SeqPro::SequenceManager>>) {
            fasta_path = manager_ptr->getFastaPath();
        } else if constexpr (std::is_same_v<PtrType,
                                 std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            fasta_path = manager_ptr->getOriginalManager().getFastaPath();
        }
    }, fasta_manager);
    const bool legacy_caps_layout =
        fast_mode && !isFileSmallerThan(fasta_path, 1024);
    if (legacy_caps_layout) {
        // FM_Index::buildIndexUsingCaPS indexes T.size()+1 symbols. Keep the
        // same double-sentinel text so the replacement preserves every
        // historical interval and occurrence-order boundary.
        forward_text.push_back('\0');
    }
    if (forward_text.size() > std::numeric_limits<uint_t>::max()) {
        throw std::runtime_error("Reference is too large for RaMAx coordinates");
    }

    total_size = logical_size;
    suffix_count = static_cast<uint_t>(forward_text.size());
    reverse_text.assign(forward_text.rbegin(), forward_text.rend());

    const auto build_begin = std::chrono::steady_clock::now();
    const uint_t requested_threads = std::max<uint_t>(1, thread_count);
    excluded_reverse_positions.clear();

    suffix_array_32.clear();
    inverse_suffix_array_32.clear();
    lcp_32.clear();
    suffix_array_64.clear();
    inverse_suffix_array_64.clear();
    lcp_64.clear();

    if (reverse_text.size() <=
        static_cast<size_t>(std::numeric_limits<saidx_t>::max())) {
        coordinates_are_64_bit = false;
        std::vector<uint32_t> forward_suffix_array(forward_text.size());
        int status = divsufsort(
            reinterpret_cast<const sauchar_t*>(forward_text.data()),
            reinterpret_cast<saidx_t*>(forward_suffix_array.data()),
            static_cast<saidx_t>(forward_text.size()));
        if (status != 0) {
            throw std::runtime_error("divsufsort32 failed to build legacy suffix order");
        }
        const uint64_t legacy_root_end = total_size - 1;
        for (uint64_t row = legacy_root_end; row < suffix_count; ++row) {
            const uint64_t suffix = forward_suffix_array[static_cast<size_t>(row)];
            const uint64_t predecessor = (suffix + suffix_count - 1) % suffix_count;
            excluded_reverse_positions.push_back(suffix_count - 1 - predecessor);
        }
        suffix_array_32.resize(reverse_text.size());
        status = divsufsort(
            reinterpret_cast<const sauchar_t*>(reverse_text.data()),
            reinterpret_cast<saidx_t*>(suffix_array_32.data()),
            static_cast<saidx_t>(reverse_text.size()));
        if (status != 0) {
            throw std::runtime_error("divsufsort32 failed to build suffix array");
        }
        inverse_suffix_array_32.resize(reverse_text.size());
#pragma omp parallel for schedule(static) num_threads(static_cast<int>(requested_threads))
        for (long long row = 0;
             row < static_cast<long long>(suffix_array_32.size()); ++row) {
            inverse_suffix_array_32[suffix_array_32[static_cast<size_t>(row)]] =
                static_cast<uint32_t>(row);
        }
        lcp_32 = buildLcp(reverse_text, suffix_array_32,
            inverse_suffix_array_32);
    } else {
        coordinates_are_64_bit = true;
        std::vector<uint64_t> forward_suffix_array(forward_text.size());
        int status = divsufsort64(
            reinterpret_cast<const sauchar_t*>(forward_text.data()),
            reinterpret_cast<saidx64_t*>(forward_suffix_array.data()),
            static_cast<saidx64_t>(forward_text.size()));
        if (status != 0) {
            throw std::runtime_error("divsufsort64 failed to build legacy suffix order");
        }
        const uint64_t legacy_root_end = total_size - 1;
        for (uint64_t row = legacy_root_end; row < suffix_count; ++row) {
            const uint64_t suffix = forward_suffix_array[static_cast<size_t>(row)];
            const uint64_t predecessor = (suffix + suffix_count - 1) % suffix_count;
            excluded_reverse_positions.push_back(suffix_count - 1 - predecessor);
        }
        suffix_array_64.resize(reverse_text.size());
        status = divsufsort64(
            reinterpret_cast<const sauchar_t*>(reverse_text.data()),
            reinterpret_cast<saidx64_t*>(suffix_array_64.data()),
            static_cast<saidx64_t>(reverse_text.size()));
        if (status != 0) {
            throw std::runtime_error("divsufsort64 failed to build suffix array");
        }
        inverse_suffix_array_64.resize(reverse_text.size());
#pragma omp parallel for schedule(static) num_threads(static_cast<int>(requested_threads))
        for (long long row = 0;
             row < static_cast<long long>(suffix_array_64.size()); ++row) {
            inverse_suffix_array_64[suffix_array_64[static_cast<size_t>(row)]] =
                static_cast<uint64_t>(row);
        }
        lcp_64 = buildLcp(reverse_text, suffix_array_64,
            inverse_suffix_array_64);
    }
    std::sort(excluded_reverse_positions.begin(),
        excluded_reverse_positions.end());
    excluded_reverse_positions.erase(std::unique(excluded_reverse_positions.begin(),
        excluded_reverse_positions.end()), excluded_reverse_positions.end());
    buildPrefixDirectory();

    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - build_begin).count();
    spdlog::info(
        "Suffix-array index built for {}: symbols={}, suffixes={}, "
        "coordinates={} bit, SA+ISA+LCP, suffix-links=enabled, "
        "prefix-k={}, SIMD-LCP={}, seconds={:.3f}",
        species_name, total_size, suffix_count,
        coordinates_are_64_bit ? 64 : 32, kSaPrefixLength,
#if defined(__SSE2__)
        "SSE2",
#else
        "scalar",
#endif
        seconds);
    return true;
}

uint64_t Suffix_Array_Index::suffixAt(uint64_t row) const noexcept {
    return coordinates_are_64_bit
        ? suffix_array_64[static_cast<size_t>(row)]
        : suffix_array_32[static_cast<size_t>(row)];
}

uint64_t Suffix_Array_Index::inverseAt(uint64_t position) const noexcept {
    return coordinates_are_64_bit
        ? inverse_suffix_array_64[static_cast<size_t>(position)]
        : inverse_suffix_array_32[static_cast<size_t>(position)];
}

uint64_t Suffix_Array_Index::lcpAt(uint64_t row) const noexcept {
    return coordinates_are_64_bit
        ? lcp_64[static_cast<size_t>(row)]
        : lcp_32[static_cast<size_t>(row)];
}

bool Suffix_Array_Index::isExcludedReversePosition(
    uint64_t position) const noexcept {
    return std::binary_search(excluded_reverse_positions.begin(),
        excluded_reverse_positions.end(), position);
}

bool Suffix_Array_Index::isExcludedRow(uint64_t row) const noexcept {
    return isExcludedReversePosition(suffixAt(row));
}

void Suffix_Array_Index::buildPrefixDirectory() {
    prefix_directory.assign(kSaPrefixCount, SAInterval{});
    for (uint64_t row = 0; row < suffix_count; ++row) {
        const uint64_t suffix = suffixAt(row);
        if (suffix + kSaPrefixLength > suffix_count) continue;
        size_t code = 0;
        if (!encodeSaPrefix(reverse_text.data() + suffix,
                kSaPrefixLength, code)) {
            continue;
        }
        SAInterval& interval = prefix_directory[code];
        if (interval.empty()) interval.l = static_cast<uint_t>(row);
        interval.r = static_cast<uint_t>(row + 1);
    }
}

SAInterval Suffix_Array_Index::suffixLinkInterval(SAInterval previous,
    uint_t depth, uint_t shift) const {
    if (previous.empty() || shift == 0 || depth <= shift) return {};

    uint64_t first = previous.l;
    while (first < previous.r && isExcludedRow(first)) ++first;
    uint64_t last = previous.r;
    while (last > first && isExcludedRow(last - 1)) --last;
    if (first == last) return {};

    const uint64_t left_suffix = suffixAt(first);
    const uint64_t right_suffix = suffixAt(last - 1);
    if (left_suffix + shift >= total_size ||
        right_suffix + shift >= total_size) {
        return {};
    }

    uint64_t left = std::min(inverseAt(left_suffix + shift),
        inverseAt(right_suffix + shift));
    uint64_t right = std::max(inverseAt(left_suffix + shift),
        inverseAt(right_suffix + shift)) + 1;
    const uint64_t target_depth = depth - shift;

    while (left > 0 && lcpAt(left) >= target_depth) --left;
    while (right < suffix_count && lcpAt(right) >= target_depth) ++right;
    return {static_cast<uint_t>(left), static_cast<uint_t>(right)};
}

Suffix_Array_Index::PrefixMatch Suffix_Array_Index::longestPrefix(
    const char* query, uint_t query_length, SAInterval search_range,
    uint_t known_depth, uint_t accepted_frequency_limit) const {
    PrefixMatch result;
    if (query == nullptr || query_length == 0 || search_range.empty()) {
        return result;
    }

    const uint_t searchable_length = query_length;
    if (searchable_length == 0 || known_depth > searchable_length) return result;

    const auto* pattern = reinterpret_cast<const uint8_t*>(query);
    auto compare_row = [&](uint64_t row, uint_t known) {
        const uint64_t suffix = suffixAt(row);
        return comparePatternBytes(
            reinterpret_cast<const uint8_t*>(reverse_text.data()) + suffix,
            reverse_text.size() - static_cast<size_t>(suffix), pattern,
            searchable_length, known);
    };

    uint64_t lower = search_range.l;
    uint64_t upper = search_range.r;
    uint_t left_lcp = known_depth;
    uint_t right_lcp = known_depth;
    while (lower < upper) {
        const uint64_t middle = lower + (upper - lower) / 2;
        const uint_t known = std::min(left_lcp, right_lcp);
        const ByteComparison comparison = compare_row(middle, known);
        if (comparison.order < 0) {
            lower = middle + 1;
            left_lcp = static_cast<uint_t>(comparison.lcp);
        } else {
            upper = middle;
            right_lcp = static_cast<uint_t>(comparison.lcp);
        }
    }

    bool have_best = false;
    uint64_t best_row = search_range.l;
    uint_t best_length = 0;
    const auto consider = [&](uint64_t row, uint_t& length,
                              uint64_t& selected, bool& present) {
        if (row < search_range.l || row >= search_range.r) return;
        const uint_t candidate = static_cast<uint_t>(
            compare_row(row, known_depth).lcp);
        if (!present || candidate > length ||
            (candidate == length && row < selected)) {
            present = true;
            length = candidate;
            selected = row;
        }
    };
    uint64_t candidate = lower;
    while (candidate < search_range.r && isExcludedRow(candidate)) ++candidate;
    consider(candidate, best_length, best_row, have_best);
    candidate = lower;
    while (candidate > search_range.l) {
        --candidate;
        if (!isExcludedRow(candidate)) {
            consider(candidate, best_length, best_row, have_best);
            break;
        }
    }
    if (!have_best || best_length == 0) return result;

    result.length = best_length;
    uint64_t left = best_row;
    uint64_t right = best_row + 1;
    uint64_t frequency = 1;
    if (frequency > accepted_frequency_limit) return result;

    while (left > search_range.l && lcpAt(left) >= best_length) {
        --left;
        if (!isExcludedRow(left) &&
            ++frequency > accepted_frequency_limit) return result;
    }
    while (right < search_range.r && lcpAt(right) >= best_length) {
        if (!isExcludedRow(right) &&
            ++frequency > accepted_frequency_limit) return result;
        ++right;
    }

    result.interval = {
        static_cast<uint_t>(left), static_cast<uint_t>(right)};
    result.frequency = static_cast<uint_t>(frequency);
    result.interval_complete = true;
    return result;
}

bool Suffix_Array_Index::originalSuffixLess(uint_t left, uint_t right,
    uint_t shared_prefix) const noexcept {
    if (left == right) return false;
    uint64_t left_pos = left + shared_prefix;
    uint64_t right_pos = right + shared_prefix;
    while (left_pos < suffix_count && right_pos < suffix_count) {
        const uint8_t left_byte = static_cast<uint8_t>(
            reverse_text[static_cast<size_t>(suffix_count - 1 - left_pos)]);
        const uint8_t right_byte = static_cast<uint8_t>(
            reverse_text[static_cast<size_t>(suffix_count - 1 - right_pos)]);
        if (left_byte != right_byte) return left_byte < right_byte;
        ++left_pos;
        ++right_pos;
    }
    if (left_pos == suffix_count && right_pos != suffix_count) return true;
    if (right_pos == suffix_count && left_pos != suffix_count) return false;
    return left < right;
}

void Suffix_Array_Index::appendRegion(uint_t ref_global_pos,
    uint_t match_length, RegionVec& region_vec,
    sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval) const {
    const auto cache_size = ref_global_cache.size();
    std::visit([&](auto&& manager_ptr) {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        SeqPro::SequenceId seq_id = SeqPro::SequenceIndex::INVALID_ID;

        if (cache_size > 0) {
            const auto cache_index = ref_global_pos / sampling_interval + 1;
            if (cache_index < cache_size) {
                const auto candidate_seq_id = ref_global_cache[cache_index];
                if (candidate_seq_id != SeqPro::SequenceIndex::INVALID_ID) {
                    const SeqPro::SequenceInfo* candidate_info = nullptr;
                    if constexpr (std::is_same_v<PtrType,
                                      std::unique_ptr<SeqPro::SequenceManager>>) {
                        candidate_info = manager_ptr->getIndex()
                            .getSequenceInfo(candidate_seq_id);
                    } else if constexpr (std::is_same_v<PtrType,
                                             std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                        candidate_info = manager_ptr->getOriginalManager()
                            .getIndex().getSequenceInfo(candidate_seq_id);
                    }
                    if (candidate_info &&
                        ref_global_pos >= candidate_info->masked_global_start_pos &&
                        ref_global_pos < candidate_info->masked_global_start_pos +
                            candidate_info->masked_length) {
                        seq_id = candidate_seq_id;
                        if constexpr (std::is_same_v<PtrType,
                                          std::unique_ptr<SeqPro::SequenceManager>>) {
                            const auto [id, local_pos] =
                                manager_ptr->globalToLocal(ref_global_pos);
                            region_vec.emplace_back(id, local_pos, match_length);
                        } else if constexpr (std::is_same_v<PtrType,
                                                 std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                            auto local_pos = ref_global_pos -
                                candidate_info->masked_global_start_pos;
                            local_pos = manager_ptr->toOriginalPositionSeparated(
                                seq_id, local_pos);
                            region_vec.emplace_back(seq_id, local_pos, match_length);
                        }
                    }
                }
            }
        }

        if (seq_id == SeqPro::SequenceIndex::INVALID_ID) {
            if constexpr (std::is_same_v<PtrType,
                              std::unique_ptr<SeqPro::SequenceManager>>) {
                const auto [fallback_seq_id, fallback_local_pos] =
                    manager_ptr->globalToLocal(ref_global_pos);
                region_vec.emplace_back(
                    fallback_seq_id, fallback_local_pos, match_length);
            } else if constexpr (std::is_same_v<PtrType,
                                     std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                const auto [fallback_seq_id, fallback_local_pos] =
                    manager_ptr->globalToLocalSeparated(ref_global_pos);
                region_vec.emplace_back(
                    fallback_seq_id, fallback_local_pos, match_length);
            }
        }
    }, fasta_manager);
}

uint_t Suffix_Array_Index::findSubSeqAnchorsWithCursor(const char* query,
    uint_t query_length, uint_t shift_from_previous, SearchCursor* cursor,
    bool allow_MEM, RegionVec& region_vec, uint_t min_anchor_length,
    bool allow_short_mum, uint_t max_anchor_frequency,
    sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval) const {
    (void)allow_short_mum;

    SAInterval search_range{0, suffix_count};
    uint_t known_depth = 0;
    if (cursor != nullptr && cursor->reusable && shift_from_previous > 0 &&
        cursor->raw_match_length > shift_from_previous) {
        const SAInterval linked = suffixLinkInterval(cursor->interval,
            cursor->raw_match_length, shift_from_previous);
        if (!linked.empty()) {
            search_range = linked;
            known_depth = cursor->raw_match_length - shift_from_previous;
        }
    }
    if (cursor != nullptr) *cursor = SearchCursor{};

    if (known_depth == 0 && min_anchor_length >= kSaPrefixLength) {
        size_t prefix_code = 0;
        if (!encodeSaPrefix(query, query_length, prefix_code)) return 1;
        const SAInterval prefix_range = prefix_directory[prefix_code];
        if (prefix_range.empty()) return 1;
        search_range = prefix_range;
        known_depth = kSaPrefixLength;
    }

    uint_t accepted_frequency_limit = max_anchor_frequency;
    if (!allow_MEM) {
        accepted_frequency_limit =
            std::min<uint_t>(accepted_frequency_limit, 1);
    }

    PrefixMatch prefix = longestPrefix(query, query_length, search_range,
        known_depth, accepted_frequency_limit);
    if (prefix.length < min_anchor_length || !prefix.interval_complete) {
        return 1;
    }

    const uint_t frequency = prefix.frequency;
    if (frequency > max_anchor_frequency || (!allow_MEM && frequency > 1)) {
        return 1;
    }

    if (cursor != nullptr) {
        cursor->interval = prefix.interval;
        cursor->raw_match_length = prefix.length;
        cursor->reusable = true;
    }

    std::vector<uint_t> global_positions;
    global_positions.reserve(frequency);
    for (uint_t row = prefix.interval.l; row < prefix.interval.r; ++row) {
        const uint64_t reverse_position = suffixAt(row);
        if (isExcludedReversePosition(reverse_position)) continue;
        if (reverse_position + prefix.length > suffix_count) {
            throw std::runtime_error("Suffix-array match exceeds reference text");
        }
        global_positions.push_back(static_cast<uint_t>(
            suffix_count - reverse_position - prefix.length));
    }
    if (global_positions.size() > 1) {
        std::sort(global_positions.begin(), global_positions.end(),
            [&](uint_t left, uint_t right) {
                return originalSuffixLess(left, right, prefix.length);
            });
    }

    region_vec.reserve(region_vec.size() + global_positions.size());
    for (const uint_t position : global_positions) {
        appendRegion(position, prefix.length, region_vec,
            ref_global_cache, sampling_interval);
    }
    return prefix.length;
}

uint_t Suffix_Array_Index::findSubSeqAnchors(const char* query,
    uint_t query_length, bool allow_MEM, RegionVec& region_vec,
    uint_t min_anchor_length, bool allow_short_mum,
    uint_t max_anchor_frequency, sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval) const {
    uint_t searchable_length = 0;
    while (searchable_length < query_length &&
           isSaCanonicalBase(query[searchable_length])) ++searchable_length;
    return findSubSeqAnchorsWithCursor(query, searchable_length, 0, nullptr,
        allow_MEM, region_vec, min_anchor_length, allow_short_mum,
        max_anchor_frequency, ref_global_cache, sampling_interval);
}

MatchVec2DPtr Suffix_Array_Index::findAnchorsImpl(
    ChrIndex query_chr_index, std::string query, SearchMode search_mode,
    Strand strand, bool allow_MEM, uint_t query_offset,
    uint_t min_anchor_length, bool allow_short_mum,
    uint_t max_anchor_frequency, sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval,
    uint_t accurate_skip_threshold) const {
    if (strand == Strand::FORWARD) {
        std::reverse(query.begin(), query.end());
    } else {
        for (char& ch : query) {
            ch = BASE_COMPLEMENT[static_cast<unsigned char>(ch)];
        }
    }

    auto anchors = std::make_shared<MatchVec2D>();
    const uint_t query_length = query.length();
    uint_t total_length = 0;
    uint_t last_pos = 0;
    uint_t shift_from_previous = 0;
    uint_t searchable_end = 0;
    SearchCursor cursor;

    while (total_length < query_length) {
        if (total_length >= searchable_end) {
            searchable_end = total_length;
            while (searchable_end < query_length &&
                   isSaCanonicalBase(query[searchable_end])) ++searchable_end;
        }
        RegionVec regions;
        const uint_t match_length = findSubSeqAnchorsWithCursor(
            query.c_str() + total_length, searchable_end - total_length,
            shift_from_previous, &cursor, allow_MEM, regions,
            min_anchor_length, allow_short_mum, max_anchor_frequency,
            ref_global_cache, sampling_interval);

        const Coord_t query_start = strand == Strand::FORWARD
            ? query_length - match_length - total_length + query_offset
            : total_length + query_offset;

        if (!regions.empty()) {
            bool append = true;
            if (search_mode == FAST_SEARCH) {
                const uint_t ref_end_pos = regions[0].start + match_length;
                append = ref_end_pos != last_pos;
                last_pos = ref_end_pos;
            } else if (search_mode == ACCURATE_SEARCH) {
                const uint_t ref_end_pos = regions[0].start;
                append = ref_end_pos != last_pos;
                last_pos = ref_end_pos;
            }

            if (append) {
                MatchVec matches;
                matches.reserve(regions.size());
                for (const auto& region : regions) {
                    matches.emplace_back(region.chr_index, region.start,
                        query_chr_index, query_start, match_length, strand);
                }
                if (!matches.empty()) anchors->emplace_back(std::move(matches));
            }
        }

        uint_t advance = match_length;
        if (search_mode == FAST_SEARCH) {
            advance = std::min(min_anchor_length,
                match_length == 1 ? min_anchor_length : match_length);
        } else if (search_mode == ACCURATE_SEARCH) {
            advance = chooseAccurateSearchAdvance(match_length,
                regions.size(), accurate_skip_threshold);
        }
        shift_from_previous = advance;
        total_length += advance;
    }

    anchors->shrink_to_fit();
    return anchors;
}

MatchVec2DPtr Suffix_Array_Index::findAnchors(ChrIndex query_chr_index,
    std::string query, SearchMode search_mode, Strand strand,
    bool allow_MEM, uint_t query_offset, uint_t min_anchor_length,
    bool allow_short_mum, uint_t max_anchor_frequency,
    sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval,
    uint_t accurate_skip_threshold) const {
    if (search_mode != FAST_SEARCH && search_mode != MIDDLE_SEARCH &&
        search_mode != ACCURATE_SEARCH) {
        throw std::invalid_argument("Invalid search mode");
    }
    return findAnchorsImpl(query_chr_index, std::move(query), search_mode,
        strand, allow_MEM, query_offset, min_anchor_length,
        allow_short_mum, max_anchor_frequency, ref_global_cache,
        sampling_interval, accurate_skip_threshold);
}

MatchVec2DPtr Suffix_Array_Index::findAnchorsFast(
    ChrIndex query_chr_index, std::string query, Strand strand,
    bool allow_MEM, uint_t query_offset, uint_t min_anchor_length,
    bool allow_short_mum, uint_t max_anchor_frequency,
    sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval) const {
    return findAnchorsImpl(query_chr_index, std::move(query), FAST_SEARCH,
        strand, allow_MEM, query_offset, min_anchor_length,
        allow_short_mum, max_anchor_frequency, ref_global_cache,
        sampling_interval, 0);
}

MatchVec2DPtr Suffix_Array_Index::findAnchorsMiddle(
    ChrIndex query_chr_index, std::string query, Strand strand,
    bool allow_MEM, uint_t query_offset, uint_t min_anchor_length,
    bool allow_short_mum, uint_t max_anchor_frequency,
    sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval) const {
    return findAnchorsImpl(query_chr_index, std::move(query), MIDDLE_SEARCH,
        strand, allow_MEM, query_offset, min_anchor_length,
        allow_short_mum, max_anchor_frequency, ref_global_cache,
        sampling_interval, 0);
}

MatchVec2DPtr Suffix_Array_Index::findAnchorsAccurate(
    ChrIndex query_chr_index, std::string query, Strand strand,
    bool allow_MEM, uint_t query_offset, uint_t min_anchor_length,
    bool allow_short_mum, uint_t max_anchor_frequency,
    sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval,
    uint_t accurate_skip_threshold) const {
    return findAnchorsImpl(query_chr_index, std::move(query), ACCURATE_SEARCH,
        strand, allow_MEM, query_offset, min_anchor_length,
        allow_short_mum, max_anchor_frequency, ref_global_cache,
        sampling_interval, accurate_skip_threshold);
}

void Suffix_Array_Index::bisectAnchors(const std::string& query,
    ChrIndex query_chr_index, Strand strand, bool allow_MEM,
    uint_t query_offset, uint_t query_length, uint_t min_anchor_length,
    bool allow_short_mum, uint_t max_anchor_frequency,
    const MUMInfo& left, const MUMInfo& right, MatchVec2D& out,
    sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval) const {
    if (right.pos <= left.pos + 1) return;
    const uint_t mid = left.pos + (right.pos - left.pos) / 2;
    RegionVec regions;
    const uint_t mid_length = findSubSeqAnchors(
        query.c_str() + mid, query_length - mid, allow_MEM, regions,
        min_anchor_length, allow_short_mum, max_anchor_frequency,
        ref_global_cache, sampling_interval);
    const bool mid_is_mum = regions.size() == 1;
    const bool same_as_left = left.pos + left.len == mid + mid_length;
    const bool same_as_right = right.pos + right.len == mid + mid_length;

    if (((regions.size() == 1 || allow_MEM) && !mid_is_mum) ||
        (!same_as_left && !same_as_right)) {
        const Coord_t query_start = strand == Strand::FORWARD
            ? query_length - mid_length - mid + query_offset
            : query_offset + mid;
        MatchVec matches;
        matches.reserve(regions.size());
        for (const auto& region : regions) {
            matches.emplace_back(region.chr_index, region.start,
                query_chr_index, query_start, mid_length, strand);
        }
        if (!matches.empty()) out.emplace_back(std::move(matches));
    }

    if (!(left.is_mum && same_as_left)) {
        bisectAnchors(query, query_chr_index, strand, allow_MEM, query_offset,
            query_length, min_anchor_length, allow_short_mum,
            max_anchor_frequency, left,
            {mid, mid_length, mid_is_mum}, out, ref_global_cache,
            sampling_interval);
    }
    if (!(right.is_mum && same_as_right)) {
        bisectAnchors(query, query_chr_index, strand, allow_MEM, query_offset,
            query_length, min_anchor_length, allow_short_mum,
            max_anchor_frequency, {mid, mid_length, mid_is_mum}, right,
            out, ref_global_cache, sampling_interval);
    }
}

bool Suffix_Array_Index::saveToFile(const std::string& filename) const {
    std::ofstream output(filename, std::ios::binary | std::ios::trunc);
    if (!output) return false;
    output.write(kIndexMagic.data(), kIndexMagic.size());
    writePod(output, kIndexVersion);
    const uint8_t width = coordinates_are_64_bit ? 64 : 32;
    writePod(output, width);
    const uint64_t size = total_size;
    writePod(output, size);
    const uint64_t stored_suffix_count = suffix_count;
    writePod(output, stored_suffix_count);
    output.write(reverse_text.data(), static_cast<std::streamsize>(stored_suffix_count));
    if (!output) throw std::runtime_error("Failed to write suffix-array text");
    if (coordinates_are_64_bit) {
        writeCoordinates(output, suffix_array_64);
        writeCoordinates(output, inverse_suffix_array_64);
        writeCoordinates(output, lcp_64);
    } else {
        writeCoordinates(output, suffix_array_32);
        writeCoordinates(output, inverse_suffix_array_32);
        writeCoordinates(output, lcp_32);
    }
    const uint64_t excluded_count = excluded_reverse_positions.size();
    writePod(output, excluded_count);
    for (const uint64_t position : excluded_reverse_positions) {
        writePod(output, position);
    }
    output.flush();
    return static_cast<bool>(output);
}

bool Suffix_Array_Index::loadFromFile(const std::string& filename) {
    std::ifstream input(filename, std::ios::binary);
    if (!input) return false;
    std::array<char, kIndexMagic.size()> magic{};
    input.read(magic.data(), magic.size());
    if (!input || magic != kIndexMagic) {
        throw std::runtime_error("Invalid suffix-array index magic");
    }
    uint32_t version = 0;
    uint8_t width = 0;
    uint64_t size = 0;
    uint64_t stored_suffix_count = 0;
    readPod(input, version);
    readPod(input, width);
    readPod(input, size);
    readPod(input, stored_suffix_count);
    if (version != kIndexVersion || (width != 32 && width != 64) ||
        size == 0 || stored_suffix_count < size ||
        stored_suffix_count > std::numeric_limits<uint_t>::max() ||
        stored_suffix_count > std::numeric_limits<size_t>::max()) {
        throw std::runtime_error("Unsupported suffix-array index header");
    }

    total_size = static_cast<uint_t>(size);
    suffix_count = static_cast<uint_t>(stored_suffix_count);
    reverse_text.resize(static_cast<size_t>(stored_suffix_count));
    input.read(reverse_text.data(), static_cast<std::streamsize>(stored_suffix_count));
    if (!input) throw std::runtime_error("Truncated suffix-array text");

    coordinates_are_64_bit = width == 64;
    suffix_array_32.clear();
    inverse_suffix_array_32.clear();
    lcp_32.clear();
    suffix_array_64.clear();
    inverse_suffix_array_64.clear();
    lcp_64.clear();
    if (coordinates_are_64_bit) {
        readCoordinates(input, suffix_array_64, static_cast<size_t>(stored_suffix_count));
        readCoordinates(input, inverse_suffix_array_64,
            static_cast<size_t>(stored_suffix_count));
        readCoordinates(input, lcp_64, static_cast<size_t>(stored_suffix_count));
        validateCoordinates(suffix_array_64, inverse_suffix_array_64, lcp_64);
    } else {
        readCoordinates(input, suffix_array_32, static_cast<size_t>(stored_suffix_count));
        readCoordinates(input, inverse_suffix_array_32,
            static_cast<size_t>(stored_suffix_count));
        readCoordinates(input, lcp_32, static_cast<size_t>(stored_suffix_count));
        validateCoordinates(suffix_array_32, inverse_suffix_array_32, lcp_32);
    }
    uint64_t excluded_count = 0;
    readPod(input, excluded_count);
    if (excluded_count > suffix_count) {
        throw std::runtime_error("Invalid excluded suffix count");
    }
    excluded_reverse_positions.resize(static_cast<size_t>(excluded_count));
    for (uint64_t& position : excluded_reverse_positions) {
        readPod(input, position);
    }
    if (!std::is_sorted(excluded_reverse_positions.begin(),
            excluded_reverse_positions.end()) ||
        (!excluded_reverse_positions.empty() &&
         excluded_reverse_positions.back() >= suffix_count)) {
        throw std::runtime_error("Invalid excluded suffix positions");
    }
    buildPrefixDirectory();
    return true;
}
