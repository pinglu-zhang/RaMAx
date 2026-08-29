#include "suffix_array_index.h"
#include "runtime_resources.h"
#include "simd_bytes.h"
#include "CaPS-SA/Suffix_Array.hpp"
#include "parlay/parallel.h"

#include <algorithm>
#include <atomic>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <sys/resource.h>

#include <omp.h>

extern "C" {
#include "divsufsort.h"
#include "divsufsort64.h"
}

namespace {

// The SIMD comparison kernel and generalized suffix-link interval derivation
// are adapted from sufkit (MIT). CaPS constructs the complete SA and LCP arrays
// together; RaMAx keeps only the query/search adapter local.
using ByteComparison = RaMAxSimd::ByteComparison;

ByteComparison comparePatternBytes(const uint8_t* text, size_t text_size,
    const uint8_t* pattern, size_t pattern_size, size_t known_lcp = 0) {
    return RaMAxSimd::compareBytes(
        text, text_size, pattern, pattern_size, known_lcp);
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

template<class Coordinate>
double buildCompleteLcp(const std::string& text,
    const RaMAxSuffixDetail::UninitializedBuffer<Coordinate>& suffix_array,
    const RaMAxSuffixDetail::UninitializedBuffer<Coordinate>& inverse_suffix_array,
    RaMAxSuffixDetail::UninitializedBuffer<Coordinate>& lcp) {
    const auto begin = std::chrono::steady_clock::now();
    uint64_t common = 0;
    for (uint64_t suffix = 0; suffix < suffix_array.size(); ++suffix) {
        const uint64_t row = inverse_suffix_array[static_cast<size_t>(suffix)];
        if (row == 0) {
            lcp[0] = 0;
            common = 0;
            continue;
        }
        const uint64_t previous = suffix_array[static_cast<size_t>(row - 1)];
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
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now() - begin).count();
}

template<class Coordinate>
double buildDivsufsortSa(const std::string& text,
    RaMAxSuffixDetail::UninitializedBuffer<Coordinate>& suffix_array,
    const char* orientation) {
    if (text.empty() || suffix_array.size() != text.size()) {
        throw std::invalid_argument("divsufsort requires a complete output buffer");
    }
    const auto begin = std::chrono::steady_clock::now();
    int status = 0;
    if constexpr (std::is_same_v<Coordinate, uint32_t>) {
        if (text.size() > static_cast<size_t>(
                std::numeric_limits<saidx_t>::max())) {
            throw std::runtime_error(
                "Reference text exceeds divsufsort32 coordinate width");
        }
        status = divsufsort(
            reinterpret_cast<const sauchar_t*>(text.data()),
            reinterpret_cast<saidx_t*>(suffix_array.data()),
            static_cast<saidx_t>(text.size()));
    } else {
        status = divsufsort64(
            reinterpret_cast<const sauchar_t*>(text.data()),
            reinterpret_cast<saidx64_t*>(suffix_array.data()),
            static_cast<saidx64_t>(text.size()));
    }
    if (status != 0) {
        throw std::runtime_error(
            std::string("divsufsort failed for ") + orientation);
    }
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - begin).count();
    spdlog::info(
        "divsufsort SA built: orientation={}, symbols={}, seconds={:.3f}",
        orientation, text.size(), seconds);
    return seconds;
}

template<class Coordinate>
double buildCaPsSaLcp(const std::string& text,
    RaMAxSuffixDetail::UninitializedBuffer<Coordinate>& suffix_array,
    RaMAxSuffixDetail::UninitializedBuffer<Coordinate>& lcp,
    uint_t thread_count, const char* orientation) {
    if (text.empty() || suffix_array.size() != text.size() ||
        lcp.size() != text.size()) {
        throw std::invalid_argument("CaPS requires complete output buffers");
    }
    if (text.size() > static_cast<uint64_t>(
            std::numeric_limits<Coordinate>::max())) {
        throw std::runtime_error("Reference text exceeds CaPS coordinate width");
    }

    const uint_t requested_threads = std::max<uint_t>(1, thread_count);
    const std::string workers = std::to_string(requested_threads);
    if (::setenv("PARLAY_NUM_THREADS", workers.c_str(), 1) != 0) {
        throw std::runtime_error("Failed to configure CaPS worker count");
    }
    const size_t actual_workers = parlay::num_workers();
    if (actual_workers != requested_threads) {
        throw std::runtime_error(
            "CaPS worker count was initialized to " +
            std::to_string(actual_workers) + ", requested " + workers);
    }

    const auto begin = std::chrono::steady_clock::now();
    CaPS_SA::Suffix_Array<Coordinate> builder(text.data(),
        static_cast<Coordinate>(text.size()), 0, 0,
        suffix_array.data(), lcp.data());
    builder.construct();
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - begin).count();
    spdlog::info(
        "CaPS SA+LCP built: orientation={}, symbols={}, workers={}, seconds={:.3f}",
        orientation, text.size(), actual_workers, seconds);
    return seconds;
}

struct SuffixTopTwo {
    static constexpr uint64_t invalid = std::numeric_limits<uint64_t>::max();
    uint64_t first{invalid};
    uint64_t second{invalid};
    uint64_t compared_bytes{0};
    bool budget_exceeded{false};
};

int compareForwardSuffixes(const std::string& text, uint64_t left,
    uint64_t right, uint64_t& compared_bytes) {
    if (left == right) return 0;
    const size_t left_size = text.size() - static_cast<size_t>(left);
    const size_t right_size = text.size() - static_cast<size_t>(right);
    const size_t common_limit = std::min(left_size, right_size);
    const ByteComparison comparison = comparePatternBytes(
        reinterpret_cast<const uint8_t*>(text.data()) + left, left_size,
        reinterpret_cast<const uint8_t*>(text.data()) + right, common_limit);
    compared_bytes += comparison.lcp +
        static_cast<uint64_t>(comparison.lcp < common_limit);
    if (comparison.order != 0) return comparison.order;
    if (left_size == right_size) return left < right ? -1 : 1;
    return left_size < right_size ? -1 : 1;
}

void considerSuffix(const std::string& text, uint64_t candidate,
    SuffixTopTwo& top) {
    if (top.first == SuffixTopTwo::invalid ||
        compareForwardSuffixes(text, candidate, top.first,
            top.compared_bytes) > 0) {
        if (candidate != top.first) top.second = top.first;
        top.first = candidate;
        return;
    }
    if (candidate != top.first &&
        (top.second == SuffixTopTwo::invalid ||
         compareForwardSuffixes(text, candidate, top.second,
             top.compared_bytes) > 0)) {
        top.second = candidate;
    }
}

SuffixTopTwo findLargestForwardSuffixes(const std::string& text,
    uint_t thread_count) {
    if (text.empty()) throw std::invalid_argument("Cannot rank empty suffix text");
    const size_t workers = static_cast<size_t>(std::max<uint_t>(1, thread_count));
    std::vector<std::array<uint64_t, 256>> byte_counts(workers);
#pragma omp parallel num_threads(static_cast<int>(workers))
    {
        const size_t worker = static_cast<size_t>(omp_get_thread_num());
        auto& counts = byte_counts[worker];
        counts.fill(0);
        const uint64_t begin = text.size() * worker / workers;
        const uint64_t end = text.size() * (worker + 1) / workers;
        for (uint64_t position = begin; position < end; ++position) {
            ++counts[static_cast<uint8_t>(text[static_cast<size_t>(position)])];
        }
    }
    std::array<uint64_t, 256> counts{};
    for (const auto& local : byte_counts) {
        for (size_t byte = 0; byte < counts.size(); ++byte) {
            counts[byte] += local[byte];
        }
    }
    int largest_byte = 255;
    while (largest_byte >= 0 && counts[static_cast<size_t>(largest_byte)] == 0) {
        --largest_byte;
    }
    int second_byte = largest_byte - 1;
    while (second_byte >= 0 && counts[static_cast<size_t>(second_byte)] == 0) {
        --second_byte;
    }
    const bool need_second_byte =
        counts[static_cast<size_t>(largest_byte)] < 2 && second_byte >= 0;
    const uint64_t per_worker_budget =
        (32ULL * static_cast<uint64_t>(text.size())) / workers + text.size();
    std::vector<SuffixTopTwo> local_top(workers);
#pragma omp parallel num_threads(static_cast<int>(workers))
    {
        const size_t worker = static_cast<size_t>(omp_get_thread_num());
        SuffixTopTwo& top = local_top[worker];
        const uint64_t begin = text.size() * worker / workers;
        const uint64_t end = text.size() * (worker + 1) / workers;
        for (uint64_t position = begin; position < end; ++position) {
            const int byte = static_cast<uint8_t>(text[static_cast<size_t>(position)]);
            if (byte != largest_byte && !(need_second_byte && byte == second_byte)) {
                continue;
            }
            considerSuffix(text, position, top);
            if (top.compared_bytes > per_worker_budget) {
                top.budget_exceeded = true;
                break;
            }
        }
    }
    SuffixTopTwo result;
    for (const SuffixTopTwo& local : local_top) {
        result.compared_bytes += local.compared_bytes;
        result.budget_exceeded = result.budget_exceeded || local.budget_exceeded;
        if (result.budget_exceeded) continue;
        if (local.first != SuffixTopTwo::invalid) {
            considerSuffix(text, local.first, result);
        }
        if (local.second != SuffixTopTwo::invalid) {
            considerSuffix(text, local.second, result);
        }
    }
    return result;
}

uint64_t peakResidentBytes() noexcept {
    rusage usage{};
    return getrusage(RUSAGE_SELF, &usage) == 0
        ? static_cast<uint64_t>(usage.ru_maxrss) * 1024ULL
        : 0;
}

constexpr bool requires64BitCoordinates(
    uint64_t text_symbols, bool use_caps_builder) noexcept {
    return use_caps_builder
        ? text_symbols > static_cast<uint64_t>(
              std::numeric_limits<uint32_t>::max())
        : text_symbols > static_cast<uint64_t>(
              std::numeric_limits<saidx_t>::max());
}

static_assert(!requires64BitCoordinates(
    std::numeric_limits<uint32_t>::max(), true));
static_assert(requires64BitCoordinates(
    static_cast<uint64_t>(std::numeric_limits<uint32_t>::max()) + 1, true));
static_assert(!requires64BitCoordinates(
    std::numeric_limits<saidx_t>::max(), false));
static_assert(requires64BitCoordinates(
    static_cast<uint64_t>(std::numeric_limits<saidx_t>::max()) + 1, false));

} // namespace

Suffix_Array_Index::Suffix_Array_Index(SpeciesName species_name,
    SeqPro::ManagerVariant& fasta_manager, uint_t sampling_rate_)
    : species_name(std::move(species_name)), fasta_manager(fasta_manager),
      sampling_rate(sampling_rate_) {
    if (sampling_rate != 1) {
        throw std::invalid_argument(
            "Suffix-array sampling rate must be exactly 1 in the hybrid backend");
    }
}

bool Suffix_Array_Index::buildIndex(FilePath output_path, bool fast_mode,
    uint_t thread_count) {
    (void)output_path;
    const auto total_begin = std::chrono::steady_clock::now();
    const uint_t requested_threads = std::max<uint_t>(1, thread_count);

    const auto concat_begin = std::chrono::steady_clock::now();
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
    const double concat_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - concat_begin).count();
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
    const bool use_caps_builder = !isFileSmallerThan(fasta_path, 1024);
    const bool legacy_caps_layout = fast_mode && use_caps_builder;
    if (legacy_caps_layout) forward_text.push_back('\0');
    if (forward_text.size() > std::numeric_limits<uint_t>::max()) {
        throw std::runtime_error("Reference is too large for RaMAx coordinates");
    }

    total_size = logical_size;
    text_size = static_cast<uint_t>(forward_text.size());
    const uint64_t excluded_count = text_size - (static_cast<uint64_t>(total_size) - 1);
    if (excluded_count == 0 || excluded_count > 2) {
        throw std::runtime_error("Unexpected FM-compatible suffix exclusion count");
    }

    excluded_reverse_positions.clear();
    suffix_array_32.clear();
    inverse_suffix_array_32.clear();
    lcp_32.clear();
    suffix_array_64.clear();
    inverse_suffix_array_64.clear();
    lcp_64.clear();

    const auto largest_begin = std::chrono::steady_clock::now();
    SuffixTopTwo largest = findLargestForwardSuffixes(
        forward_text, requested_threads);
    bool forward_fallback = largest.budget_exceeded ||
        largest.first == SuffixTopTwo::invalid ||
        (excluded_count == 2 && largest.second == SuffixTopTwo::invalid);
    double forward_fallback_seconds = 0.0;
    std::array<uint64_t, 2> excluded_suffixes{
        largest.first, largest.second};
    if (forward_fallback) {
        const bool forward_64 = requires64BitCoordinates(
            forward_text.size(), use_caps_builder);
        if (forward_64) {
            RaMAxSuffixDetail::UninitializedBuffer<uint64_t> forward_sa;
            RaMAxSuffixDetail::UninitializedBuffer<uint64_t> forward_lcp;
            forward_sa.allocate(
                forward_text.size(), "suffix-forward-fallback-sa64");
            if (use_caps_builder) {
                forward_lcp.allocate(
                    forward_text.size(), "suffix-forward-fallback-lcp64");
                forward_fallback_seconds = buildCaPsSaLcp(forward_text,
                    forward_sa, forward_lcp, requested_threads,
                    "forward-fallback");
            } else {
                forward_fallback_seconds = buildDivsufsortSa(forward_text,
                    forward_sa, "forward-fallback");
            }
            excluded_suffixes[0] = forward_sa[forward_sa.size() - 1];
            if (excluded_count == 2) {
                excluded_suffixes[1] = forward_sa[forward_sa.size() - 2];
            }
        } else {
            RaMAxSuffixDetail::UninitializedBuffer<uint32_t> forward_sa;
            RaMAxSuffixDetail::UninitializedBuffer<uint32_t> forward_lcp;
            forward_sa.allocate(
                forward_text.size(), "suffix-forward-fallback-sa32");
            if (use_caps_builder) {
                forward_lcp.allocate(
                    forward_text.size(), "suffix-forward-fallback-lcp32");
                forward_fallback_seconds = buildCaPsSaLcp(forward_text,
                    forward_sa, forward_lcp, requested_threads,
                    "forward-fallback");
            } else {
                forward_fallback_seconds = buildDivsufsortSa(forward_text,
                    forward_sa, "forward-fallback");
            }
            excluded_suffixes[0] = forward_sa[forward_sa.size() - 1];
            if (excluded_count == 2) {
                excluded_suffixes[1] = forward_sa[forward_sa.size() - 2];
            }
        }
    }
    for (uint64_t index = 0; index < excluded_count; ++index) {
        const uint64_t suffix = excluded_suffixes[static_cast<size_t>(index)];
        const uint64_t predecessor = (suffix + text_size - 1) % text_size;
        excluded_reverse_positions.push_back(text_size - 1 - predecessor);
    }
    std::sort(excluded_reverse_positions.begin(),
        excluded_reverse_positions.end());
    excluded_reverse_positions.erase(std::unique(excluded_reverse_positions.begin(),
        excluded_reverse_positions.end()), excluded_reverse_positions.end());
    const double largest_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - largest_begin).count();

    const auto reverse_begin = std::chrono::steady_clock::now();
    std::reverse(forward_text.begin(), forward_text.end());
    reverse_text = std::move(forward_text);
    const double reverse_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - reverse_begin).count();

    coordinates_are_64_bit = requires64BitCoordinates(
        reverse_text.size(), use_caps_builder);
    const auto allocation_begin = std::chrono::steady_clock::now();
    const uint64_t coordinate_width_bytes = coordinates_are_64_bit
        ? sizeof(uint64_t) : sizeof(uint32_t);
    const uint64_t suffix_storage_bytes =
        static_cast<uint64_t>(reverse_text.size()) *
        coordinate_width_bytes * 3ULL;
    auto& resources =
        RaMAxResources::RuntimeResourceManager::instance();
    if (resources.configured()) {
        resources.requireAllocation(
            suffix_storage_bytes, "suffix-array-sa-isa-lcp");
        resources.requireTempSpace(
            suffix_storage_bytes, "suffix-array-sa-isa-lcp");
    }
    if (coordinates_are_64_bit) {
        suffix_array_64.allocate(reverse_text.size(), "suffix-sa64");
        inverse_suffix_array_64.allocate(reverse_text.size(), "suffix-isa64");
        lcp_64.allocate(reverse_text.size(), "suffix-lcp64");
        suffix_array_64.adviseSequential();
        inverse_suffix_array_64.adviseSequential();
        lcp_64.adviseSequential();
    } else {
        suffix_array_32.allocate(reverse_text.size(), "suffix-sa32");
        inverse_suffix_array_32.allocate(reverse_text.size(), "suffix-isa32");
        lcp_32.allocate(reverse_text.size(), "suffix-lcp32");
        suffix_array_32.adviseSequential();
        inverse_suffix_array_32.adviseSequential();
        lcp_32.adviseSequential();
    }
    const double allocation_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - allocation_begin).count();

    double builder_seconds = 0.0;
    double lcp_seconds = 0.0;
    double isa_seconds = 0.0;
    std::atomic<bool> invalid_suffix{false};
    if (coordinates_are_64_bit) {
        if (use_caps_builder) {
            builder_seconds = buildCaPsSaLcp(reverse_text,
                suffix_array_64, lcp_64, requested_threads,
                "reverse-search");
        } else {
            builder_seconds = buildDivsufsortSa(reverse_text,
                suffix_array_64, "reverse-search");
        }
        stored_suffix_count = static_cast<uint_t>(suffix_array_64.size());
        const auto isa_begin = std::chrono::steady_clock::now();
#pragma omp parallel for schedule(static) num_threads(static_cast<int>(requested_threads))
        for (long long row = 0;
             row < static_cast<long long>(stored_suffix_count); ++row) {
            const uint64_t suffix = suffix_array_64[static_cast<size_t>(row)];
            if (suffix >= text_size) {
                invalid_suffix.store(true, std::memory_order_relaxed);
            } else {
                inverse_suffix_array_64[static_cast<size_t>(suffix)] =
                    static_cast<uint64_t>(row);
            }
        }
        isa_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - isa_begin).count();
        if (!use_caps_builder) {
            lcp_seconds = buildCompleteLcp(reverse_text, suffix_array_64,
                inverse_suffix_array_64, lcp_64);
        }
    } else {
        if (use_caps_builder) {
            builder_seconds = buildCaPsSaLcp(reverse_text,
                suffix_array_32, lcp_32, requested_threads,
                "reverse-search");
        } else {
            builder_seconds = buildDivsufsortSa(reverse_text,
                suffix_array_32, "reverse-search");
        }
        stored_suffix_count = static_cast<uint_t>(suffix_array_32.size());
        const auto isa_begin = std::chrono::steady_clock::now();
#pragma omp parallel for schedule(static) num_threads(static_cast<int>(requested_threads))
        for (long long row = 0;
             row < static_cast<long long>(stored_suffix_count); ++row) {
            const uint32_t suffix = suffix_array_32[static_cast<size_t>(row)];
            if (static_cast<uint64_t>(suffix) >= text_size) {
                invalid_suffix.store(true, std::memory_order_relaxed);
            } else {
                inverse_suffix_array_32[static_cast<size_t>(suffix)] =
                    static_cast<uint32_t>(row);
            }
        }
        isa_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - isa_begin).count();
        if (!use_caps_builder) {
            lcp_seconds = buildCompleteLcp(reverse_text, suffix_array_32,
                inverse_suffix_array_32, lcp_32);
        }
    }
    if (invalid_suffix.load(std::memory_order_relaxed) ||
        stored_suffix_count != text_size ||
        (coordinates_are_64_bit ? lcp_64[0] : lcp_32[0]) != 0) {
        throw std::runtime_error("Suffix-array builder returned invalid dimensions");
    }
    if (coordinates_are_64_bit) {
        suffix_array_64.adviseRandom();
        inverse_suffix_array_64.adviseRandom();
        lcp_64.adviseRandom();
    } else {
        suffix_array_32.adviseRandom();
        inverse_suffix_array_32.adviseRandom();
        lcp_32.adviseRandom();
    }

    const auto prefix_begin = std::chrono::steady_clock::now();
    buildPrefixDirectory(requested_threads);
    const double prefix_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - prefix_begin).count();
    const double total_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - total_begin).count();
    const uint64_t coordinate_bytes = coordinates_are_64_bit ? 8 : 4;
    const uint64_t component_bytes =
        static_cast<uint64_t>(stored_suffix_count) * coordinate_bytes;
    const uint64_t caps_workspace_estimate = use_caps_builder
        ? 2 * component_bytes : 0;
    spdlog::info(
        "Suffix-array index built for {}: storage={}, persistent-cache=disabled, "
        "persistent-bytes=0, mapped-bytes={}, builder={}, threshold=1024MiB, "
        "logical-symbols={}, "
        "text-symbols={}, stored-rows={}, sampling-rate=1, coordinates={} bit, "
        "SA-bytes={}, ISA-bytes={}, LCP-bytes={}, caps-workspace-estimate={}, "
        "LCP-source={}, suffix-links=enabled, prefix-k={}, SIMD-search={}, "
        "text-concat-seconds={:.3f}, largest-suffix-seconds={:.3f}, "
        "largest-suffix-compared-bytes={}, forward-fallback={}, "
        "forward-fallback-seconds={:.3f}, text-reverse-seconds={:.3f}, "
        "buffer-allocation-seconds={:.3f}, reverse-builder-seconds={:.3f}, "
        "lcp-seconds={:.3f}, isa-seconds={:.3f}, prefix-directory-seconds={:.3f}, "
        "total-index-build-seconds={:.3f}, peak-rss-bytes={}",
        species_name,
        ((coordinates_are_64_bit && suffix_array_64.fileBacked()) ||
         (!coordinates_are_64_bit && suffix_array_32.fileBacked()))
            ? "file-mapped" : "heap",
        ((coordinates_are_64_bit && suffix_array_64.fileBacked()) ||
         (!coordinates_are_64_bit && suffix_array_32.fileBacked()))
            ? suffix_storage_bytes : 0,
        use_caps_builder ? "CaPS" : "divsufsort",
        total_size, text_size, stored_suffix_count,
        coordinates_are_64_bit ? 64 : 32, component_bytes, component_bytes,
        component_bytes, caps_workspace_estimate,
        use_caps_builder ? "CaPS-joint" : "complete-Kasai",
        kSaPrefixLength,
        RaMAxSimd::byteKernelName(RaMAxSimd::selectedByteKernel()),
        concat_seconds, largest_seconds, largest.compared_bytes,
        forward_fallback ? "yes" : "no", forward_fallback_seconds,
        reverse_seconds, allocation_seconds, builder_seconds, lcp_seconds,
        isa_seconds, prefix_seconds, total_seconds, peakResidentBytes());
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

void Suffix_Array_Index::buildPrefixDirectory(uint_t thread_count) {
    prefix_directory.assign(kSaPrefixCount, SAInterval{});
    const auto decode = [](size_t code) {
        std::array<char, kSaPrefixLength> pattern{};
        static constexpr std::array<char, 4> bases{'A', 'C', 'G', 'T'};
        for (size_t index = kSaPrefixLength; index-- > 0;) {
            pattern[index] = bases[code & 3U];
            code >>= 2;
        }
        return pattern;
    };
    const auto compare = [&](uint64_t row,
                             const std::array<char, kSaPrefixLength>& pattern) {
        const uint64_t suffix = suffixAt(row);
        const uint64_t available = text_size - suffix;
        const size_t length = static_cast<size_t>(
            std::min<uint64_t>(available, kSaPrefixLength));
        const auto* text = reinterpret_cast<const uint8_t*>(
            reverse_text.data() + suffix);
        for (size_t index = 0; index < length; ++index) {
            const uint8_t left = text[index];
            const uint8_t right = static_cast<uint8_t>(pattern[index]);
            if (left != right) return left < right ? -1 : 1;
        }
        return length < kSaPrefixLength ? -1 : 0;
    };
#pragma omp parallel for schedule(static) num_threads(static_cast<int>(std::max<uint_t>(1, thread_count)))
    for (long long raw_code = 0;
         raw_code < static_cast<long long>(kSaPrefixCount); ++raw_code) {
        const size_t code = static_cast<size_t>(raw_code);
        const auto pattern = decode(code);
        uint64_t lower = 0;
        uint64_t upper = stored_suffix_count;
        while (lower < upper) {
            const uint64_t middle = lower + (upper - lower) / 2;
            if (compare(middle, pattern) < 0) lower = middle + 1;
            else upper = middle;
        }
        const uint64_t begin = lower;
        upper = stored_suffix_count;
        while (lower < upper) {
            const uint64_t middle = lower + (upper - lower) / 2;
            if (compare(middle, pattern) <= 0) lower = middle + 1;
            else upper = middle;
        }
        if (begin != lower) {
            prefix_directory[code] = {
                static_cast<uint_t>(begin), static_cast<uint_t>(lower)};
        }
    }
}

SAInterval Suffix_Array_Index::suffixLinkInterval(SAInterval previous,
    uint_t depth, uint_t shift) const {
    if (previous.empty() || shift == 0 || depth <= shift) {
        return {};
    }

    uint64_t first = previous.l;
    uint64_t last = previous.r;
    while (first < previous.r && isExcludedRow(first)) ++first;
    while (last > first && isExcludedRow(last - 1)) --last;
    if (first == last) return {};

    const uint64_t left_suffix = suffixAt(first);
    const uint64_t right_suffix = suffixAt(last - 1);
    const uint64_t link_limit = total_size;
    if (left_suffix + shift >= link_limit ||
        right_suffix + shift >= link_limit) {
        return {};
    }

    uint64_t left = std::min(inverseAt(left_suffix + shift),
        inverseAt(right_suffix + shift));
    uint64_t right = std::max(inverseAt(left_suffix + shift),
        inverseAt(right_suffix + shift)) + 1;
    const uint64_t target_depth = depth - shift;

    while (left > 0 && lcpAt(left) >= target_depth) --left;
    while (right < stored_suffix_count && lcpAt(right) >= target_depth) ++right;
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
    while (left_pos < text_size && right_pos < text_size) {
        const uint8_t left_byte = static_cast<uint8_t>(
            reverse_text[static_cast<size_t>(text_size - 1 - left_pos)]);
        const uint8_t right_byte = static_cast<uint8_t>(
            reverse_text[static_cast<size_t>(text_size - 1 - right_pos)]);
        if (left_byte != right_byte) return left_byte < right_byte;
        ++left_pos;
        ++right_pos;
    }
    if (left_pos == text_size && right_pos != text_size) return true;
    if (right_pos == text_size && left_pos != text_size) return false;
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
    uint_t query_length, uint_t shift_from_previous,
    SearchCursor* cursor, bool allow_MEM, RegionVec& region_vec,
    uint_t min_anchor_length, bool allow_short_mum,
    uint_t max_anchor_frequency, sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval) const {
    (void)allow_short_mum;

    uint_t accepted_frequency_limit = max_anchor_frequency;
    if (!allow_MEM) {
        accepted_frequency_limit =
            std::min<uint_t>(accepted_frequency_limit, 1);
    }

    SAInterval search_range{0, stored_suffix_count};
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
        if (reverse_position + prefix.length > text_size) {
            throw std::runtime_error("Suffix-array match exceeds reference text");
        }
        global_positions.push_back(static_cast<uint_t>(
            text_size - reverse_position - prefix.length));
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
