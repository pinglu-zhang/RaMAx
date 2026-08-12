#include "align.h"
#include "cross_anchor_repair.h"
#include "external_msa_runner.h"
#include "reference_profile_merger.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>

#include <spdlog/spdlog.h>

namespace {

constexpr uint32_t kCrossAnchorMinimumInsertion = 10;
constexpr uint32_t kCrossAnchorMaximumDistance = 5;
constexpr uint32_t kCrossAnchorFlankLength = 16;
constexpr double kCrossAnchorMinimumCoverage = 0.70;
constexpr double kCrossAnchorMinimumIdentity = 0.60;
constexpr size_t kCrossAnchorCacheLimit = 1024;

using RaMesh::Alignment::CrossAnchorRepairConfiguration;

auto& cross_anchor_session =
    RaMesh::Alignment::CrossAnchorInsertionRepairSession::instance();
auto& cross_anchor_cache_mutex = cross_anchor_session.cache_mutex;
auto& cross_anchor_msa_cache = cross_anchor_session.msa_cache;
auto& cross_anchor_counters = cross_anchor_session.counters;

uint64_t elapsedNanoseconds(
    std::chrono::steady_clock::time_point start,
    std::chrono::steady_clock::time_point finish) {
    return static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            finish - start)
            .count());
}

}  // namespace

KSW2AlignConfig makeDefaultKSW2Config() {
    static int8_t simple_dna_mat[25];
    static bool initialized = false;
    if (!initialized) {
        // A C G T N -> 0 1 2 3 4
        const int match = 2;
        const int mismatch = -3;
        const int ambiguous = -1;  // 对N的惩罚较小

        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                if (i == 4 || j == 4) {
                    simple_dna_mat[i * 5 + j] = ambiguous;  // N匹配
                }
                else if (i == j) {
                    simple_dna_mat[i * 5 + j] = match;
                }
                else {
                    simple_dna_mat[i * 5 + j] = mismatch;
                }
            }
        }
        initialized = true;
    }

    return {
        .mat = simple_dna_mat,
        .alphabet_size = 5,
        .gap_open = 8,
        .gap_extend = 1,
        .end_bonus = 0,
        .zdrop = 100,           // 较远匹配终止
        .band_width = -1,      // 合理的band可提高性能
        .flag = KSW_EZ_GENERIC_SC | KSW_EZ_RIGHT
    };
}

Cigar_t globalAlignKSW2(const std::string& ref,
    const std::string& query)
{
    /* ---------- 1. 编码序列 ---------- */
    std::vector<uint8_t> ref_enc(ref.size());
    std::vector<uint8_t> qry_enc(query.size());

    for (size_t i = 0; i < ref.size(); ++i)
        ref_enc[i] = ScoreChar2Idx[static_cast<uint8_t>(ref[i])];
    for (size_t i = 0; i < query.size(); ++i)
        qry_enc[i] = ScoreChar2Idx[static_cast<uint8_t>(query[i])];

    /* ---------- 2. 复制 cfg 并修正常见坑 ---------- */
    // KSW2AlignConfig cfg = cfg_in;                   // 本地副本可调整
    init_simd_mat();

    KSW2AlignConfig cfg;
    cfg.mat = dna5_simd_mat;
    cfg.alphabet_size = 5;
    cfg.gap_open = 5;          // gap open penalty
    cfg.gap_extend = 2;        // gap extension penalty
    cfg.end_bonus = 0;         // ❌ 不需要 ends-free 奖励
    cfg.zdrop = -1;            // ❌ 禁用 z-drop（全局比对必须完整比完）
    cfg.band_width = -1;       // 启用全矩阵（也可设 auto_band）
    cfg.flag = KSW_EZ_RIGHT; // ✅ 通用矩阵 + gap右对齐

    /* ---------- 3. 调用 KSW2 ---------- */
    ksw_extz_t ez{};

    ksw_extz2_sse(0,
        static_cast<int>(qry_enc.size()), qry_enc.data(),
        static_cast<int>(ref_enc.size()), ref_enc.data(),
        cfg.alphabet_size, cfg.mat,
        cfg.gap_open, cfg.gap_extend,
        cfg.band_width, cfg.zdrop, cfg.end_bonus,
        cfg.flag, &ez);


    /* ---------- 4. 拷贝 / 释放 CIGAR ---------- */
    Cigar_t cigar;
    cigar.reserve(ez.n_cigar);
    for (int i = 0; i < ez.n_cigar; ++i)
        cigar.push_back(ez.cigar[i]);

    free(ez.cigar);           // KSW2 用 malloc()
    return cigar;
}

Cigar_t globalAlignKSW2_2(const std::string& ref,
    const std::string& query)
{
    /* ---------- 1. 编码序列 ---------- */
    std::vector<uint8_t> ref_enc(ref.size());
    std::vector<uint8_t> qry_enc(query.size());

    for (size_t i = 0; i < ref.size(); ++i)
        ref_enc[i] = ScoreChar2Idx[static_cast<uint8_t>(ref[i])];
    for (size_t i = 0; i < query.size(); ++i)
        qry_enc[i] = ScoreChar2Idx[static_cast<uint8_t>(query[i])];

    /* ---------- 2. 复制 cfg 并修正常见坑 ---------- */
    // KSW2AlignConfig cfg = cfg_in;                   // 本地副本可调整
    KSW2AlignConfig cfg = makeTurboKSW2Config2(query.size(), ref.size());

    /* ---------- 3. 调用 KSW2 ---------- */
    ksw_extz_t ez{};

    ksw_extz2_sse(0,
        static_cast<int>(qry_enc.size()), qry_enc.data(),
        static_cast<int>(ref_enc.size()), ref_enc.data(),
        cfg.alphabet_size, cfg.mat,
        cfg.gap_open, cfg.gap_extend,
        cfg.band_width, cfg.zdrop, cfg.end_bonus,
        cfg.flag, &ez);


    /* ---------- 4. 拷贝 / 释放 CIGAR ---------- */
    Cigar_t cigar;
    cigar.reserve(ez.n_cigar);
    for (int i = 0; i < ez.n_cigar; ++i)
        cigar.push_back(ez.cigar[i]);

    free(ez.cigar);           // KSW2 用 malloc()
    return cigar;
}


/**********************************************************************
*  extendAlignKSW2  ——  ends-free（seed-and-extend）比对
*    @param ref        参考片段（目标方向）
*    @param query      查询片段（同方向；若反链请先反向互补）
*    @param zdrop      Z-drop 剪枝阈值（默认 200）
*    @param band       带宽限制；<0 表示不限制
*    @return           Cigar_t（BAM 编码）
**********************************************************************/
Cigar_t extendAlignKSW2(const std::string& ref,
    const std::string& query,
    int zdrop)
{
    /* ---------- 1. 序列编码 ---------- */
    std::vector<uint8_t> ref_enc(ref.size());
    std::vector<uint8_t> qry_enc(query.size());
    for (size_t i = 0; i < ref.size(); ++i) ref_enc[i] = ScoreChar2Idx[(uint8_t)ref[i]];
    for (size_t i = 0; i < query.size(); ++i) qry_enc[i] = ScoreChar2Idx[(uint8_t)query[i]];

    ///* ---------- 2. 配置 ---------- */
    //KSW2AlignConfig cfg = makeTurboKSW2Config(query.size(), ref.size());
    ////KSW2AlignConfig cfg;
    //cfg.zdrop = zdrop;       // 用于提前终止
    //cfg.flag = KSW_EZ_EXTZ_ONLY     // ends-free extension
    //    | KSW_EZ_APPROX_MAX    // 跟踪 ez.max_q/max_t
    //    | KSW_EZ_APPROX_DROP   // 在 approximate 模式下触发 z-drop 就中断
    //    | KSW_EZ_RIGHT;        // （可选）gap 右对齐     // **关键**：启用 extension/ends-free
    //// 若需要右对齐 gaps 建议保留 KSW_EZ_RIGHT
    //cfg.end_bonus = 100;
    //cfg.band_width = -1;
    init_simd_mat();
    KSW2AlignConfig cfg;
	cfg.mat = dna5_simd_mat;
    cfg.zdrop = zdrop;
    cfg.flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT | KSW_EZ_APPROX_DROP;
    cfg.end_bonus = 50;
    cfg.alphabet_size = 5;
    cfg.gap_open = 5;
    cfg.gap_extend = 2;
    cfg.band_width = auto_band(ref.size(), query.size());


    /* ---------- 3. 调用 KSW2 ---------- */
    ksw_extz_t ez{};
    ksw_extz2_sse(nullptr,
        static_cast<int>(qry_enc.size()), qry_enc.data(),
        static_cast<int>(ref_enc.size()), ref_enc.data(),
        cfg.alphabet_size, cfg.mat,
        cfg.gap_open, cfg.gap_extend,
        cfg.band_width, cfg.zdrop, cfg.end_bonus,
        cfg.flag, &ez);

    // 赋值bool& if_zdrop,int& ref_end,int& qry_end
    /* ---------- 4. 拷贝 & 释放 ---------- */
    Cigar_t cigar;
    cigar.reserve(ez.n_cigar);
    for (int i = 0; i < ez.n_cigar; ++i)
        cigar.push_back(ez.cigar[i]);

    free(ez.cigar);                    // ksw2 使用 malloc
    return cigar;                      // 返回的 CIGAR 即延伸片段
}

namespace {

struct CrossAnchorInsertionEvent {
    ChrName key;
    uint32_t reference_position = 0;
    uint32_t query_position = 0;
    uint32_t length = 0;
};

struct CrossAnchorSimilarity {
    size_t left = 0;
    size_t right = 0;
    size_t paired_bases = 0;
    size_t matching_bases = 0;
    double coverage = 0.0;
    double identity = 0.0;
};

struct CrossAnchorGroup {
    std::vector<size_t> event_indices;
    std::vector<CrossAnchorSimilarity> edges;
    uint32_t minimum_position = 0;
    uint32_t maximum_position = 0;
    uint32_t left_reference = 0;
    uint32_t right_reference = 0;
    size_t species_count = 0;
    size_t shared_bases = 0;
    double average_identity = 0.0;
};

struct PairAlignmentMetrics {
    size_t paired_bases = 0;
    size_t matching_bases = 0;
};

char normalizedBase(char base) {
    return static_cast<char>(
        std::toupper(static_cast<unsigned char>(base)));
}

std::string removeAlignmentGaps(const std::string& sequence) {
    std::string result;
    result.reserve(sequence.size());
    for (const char base : sequence) {
        if (base != '-') {
            result.push_back(base);
        }
    }
    return result;
}

Cigar_t globalAlignKSW2Banded(const std::string& reference,
                              const std::string& query) {
    if (reference.empty()) {
        Cigar_t result;
        if (!query.empty()) {
            appendCigarOp(
                result, 'I', static_cast<uint32_t>(query.size()));
        }
        return result;
    }
    if (query.empty()) {
        Cigar_t result;
        appendCigarOp(
            result, 'D', static_cast<uint32_t>(reference.size()));
        return result;
    }

    std::vector<uint8_t> reference_encoded(reference.size());
    std::vector<uint8_t> query_encoded(query.size());
    for (size_t index = 0; index < reference.size(); ++index) {
        reference_encoded[index] =
            ScoreChar2Idx[static_cast<uint8_t>(reference[index])];
    }
    for (size_t index = 0; index < query.size(); ++index) {
        query_encoded[index] =
            ScoreChar2Idx[static_cast<uint8_t>(query[index])];
    }

    init_simd_mat();
    const size_t difference = reference.size() > query.size()
        ? reference.size() - query.size()
        : query.size() - reference.size();
    const size_t requested_band = difference + 32;
    const int band = requested_band >
            static_cast<size_t>(std::numeric_limits<int>::max())
        ? std::numeric_limits<int>::max()
        : static_cast<int>(requested_band);

    ksw_extz_t alignment{};
    ksw_extz2_sse(
        nullptr,
        static_cast<int>(query_encoded.size()), query_encoded.data(),
        static_cast<int>(reference_encoded.size()), reference_encoded.data(),
        5, dna5_simd_mat, 5, 2, band, -1, 0,
        KSW_EZ_RIGHT, &alignment);

    Cigar_t result;
    result.reserve(alignment.n_cigar);
    for (int index = 0; index < alignment.n_cigar; ++index) {
        result.push_back(alignment.cigar[index]);
    }
    free(alignment.cigar);
    return result;
}

std::optional<CrossAnchorSimilarity> insertionSimilarity(
    size_t left_index,
    size_t right_index,
    const CrossAnchorInsertionEvent& left_event,
    const CrossAnchorInsertionEvent& right_event,
    const std::unordered_map<ChrName, std::string>& sequences) {
    const auto left_sequence_it = sequences.find(left_event.key);
    const auto right_sequence_it = sequences.find(right_event.key);
    if (left_sequence_it == sequences.end() ||
        right_sequence_it == sequences.end() ||
        left_event.query_position + left_event.length >
            left_sequence_it->second.size() ||
        right_event.query_position + right_event.length >
            right_sequence_it->second.size()) {
        return std::nullopt;
    }

    const std::string left_sequence = left_sequence_it->second.substr(
        left_event.query_position, left_event.length);
    const std::string right_sequence = right_sequence_it->second.substr(
        right_event.query_position, right_event.length);
    const Cigar_t cigar = globalAlignKSW2Banded(
        left_sequence, right_sequence);

    size_t left_position = 0;
    size_t right_position = 0;
    size_t paired = 0;
    size_t matching = 0;
    for (const auto unit : cigar) {
        uint32_t length = 0;
        char operation = '\0';
        intToCigar(unit, operation, length);
        if (operation == 'M' || operation == '=' || operation == 'X') {
            if (left_position + length > left_sequence.size() ||
                right_position + length > right_sequence.size()) {
                return std::nullopt;
            }
            for (uint32_t offset = 0; offset < length; ++offset) {
                matching += normalizedBase(
                    left_sequence[left_position + offset]) ==
                    normalizedBase(
                        right_sequence[right_position + offset]);
            }
            paired += length;
            left_position += length;
            right_position += length;
        } else if (operation == 'I') {
            right_position += length;
        } else if (operation == 'D') {
            left_position += length;
        } else {
            return std::nullopt;
        }
    }
    if (left_position != left_sequence.size() ||
        right_position != right_sequence.size() || paired == 0) {
        return std::nullopt;
    }

    const size_t shorter = std::min(
        left_sequence.size(), right_sequence.size());
    const double coverage = shorter == 0
        ? 0.0
        : static_cast<double>(paired) /
            static_cast<double>(shorter);
    const double identity = static_cast<double>(matching) /
        static_cast<double>(paired);
    if (paired < kCrossAnchorMinimumInsertion ||
        coverage < kCrossAnchorMinimumCoverage ||
        identity < kCrossAnchorMinimumIdentity) {
        return std::nullopt;
    }
    return CrossAnchorSimilarity{
        left_index, right_index, paired, matching, coverage, identity};
}

PairAlignmentMetrics alignedPairMetrics(const std::string& left,
                                        const std::string& right) {
    PairAlignmentMetrics result;
    if (left.size() != right.size()) {
        return result;
    }
    for (size_t column = 0; column < left.size(); ++column) {
        if (left[column] == '-' || right[column] == '-') {
            continue;
        }
        ++result.paired_bases;
        result.matching_bases +=
            normalizedBase(left[column]) == normalizedBase(right[column]);
    }
    return result;
}

std::vector<CrossAnchorInsertionEvent> collectInsertionEvents(
    const ChrName& reference_name,
    const std::unordered_map<ChrName, std::string>& sequences,
    const std::unordered_map<ChrName, Cigar_t>& cigars) {
    std::vector<CrossAnchorInsertionEvent> events;
    for (const auto& [key, cigar] : cigars) {
        if (key == reference_name || sequences.find(key) == sequences.end()) {
            continue;
        }
        uint64_t reference_position = 0;
        uint64_t query_position = 0;
        for (const auto unit : cigar) {
            uint32_t length = 0;
            char operation = '\0';
            intToCigar(unit, operation, length);
            if (operation == 'M' || operation == '=' || operation == 'X') {
                reference_position += length;
                query_position += length;
            } else if (operation == 'I') {
                if (length >= kCrossAnchorMinimumInsertion &&
                    reference_position <=
                        std::numeric_limits<uint32_t>::max() &&
                    query_position <=
                        std::numeric_limits<uint32_t>::max()) {
                    events.push_back(CrossAnchorInsertionEvent{
                        key,
                        static_cast<uint32_t>(reference_position),
                        static_cast<uint32_t>(query_position),
                        length});
                } else if (length < kCrossAnchorMinimumInsertion) {
                    cross_anchor_counters.short_insertions_skipped.fetch_add(
                        1, std::memory_order_relaxed);
                }
                query_position += length;
            } else if (operation == 'D' || operation == 'N') {
                reference_position += length;
            } else if (operation == 'S') {
                query_position += length;
            }
        }
    }
    return events;
}

class ConstrainedDisjointSet {
public:
    explicit ConstrainedDisjointSet(
        const std::vector<CrossAnchorInsertionEvent>& events)
        : parent_(events.size()), rank_(events.size(), 0),
          minimum_(events.size()), maximum_(events.size()) {
        std::iota(parent_.begin(), parent_.end(), 0);
        for (size_t index = 0; index < events.size(); ++index) {
            minimum_[index] = events[index].reference_position;
            maximum_[index] = events[index].reference_position;
        }
    }

    size_t find(size_t value) {
        if (parent_[value] != value) {
            parent_[value] = find(parent_[value]);
        }
        return parent_[value];
    }

    bool unite(size_t left, size_t right) {
        left = find(left);
        right = find(right);
        if (left == right) {
            return true;
        }
        const uint32_t merged_minimum =
            std::min(minimum_[left], minimum_[right]);
        const uint32_t merged_maximum =
            std::max(maximum_[left], maximum_[right]);
        if (merged_maximum - merged_minimum >
            kCrossAnchorMaximumDistance) {
            return false;
        }
        if (rank_[left] < rank_[right]) {
            std::swap(left, right);
        }
        parent_[right] = left;
        minimum_[left] = merged_minimum;
        maximum_[left] = merged_maximum;
        if (rank_[left] == rank_[right]) {
            ++rank_[left];
        }
        return true;
    }

private:
    std::vector<size_t> parent_;
    std::vector<uint8_t> rank_;
    std::vector<uint32_t> minimum_;
    std::vector<uint32_t> maximum_;
};

std::vector<CrossAnchorGroup> buildCrossAnchorGroups(
    const std::vector<CrossAnchorInsertionEvent>& events,
    const std::unordered_map<ChrName, std::string>& sequences) {
    if (events.size() < 2) {
        return {};
    }

    std::vector<CrossAnchorSimilarity> edges;
    bool structural_candidate = false;
    for (size_t left = 0; left < events.size(); ++left) {
        for (size_t right = left + 1; right < events.size(); ++right) {
            if (events[left].key == events[right].key) {
                continue;
            }
            const uint32_t left_position = events[left].reference_position;
            const uint32_t right_position = events[right].reference_position;
            if (left_position == right_position) {
                cross_anchor_counters.same_anchor_skipped.fetch_add(
                    1, std::memory_order_relaxed);
                continue;
            }
            const uint32_t distance = left_position > right_position
                ? left_position - right_position
                : right_position - left_position;
            if (distance > kCrossAnchorMaximumDistance) {
                cross_anchor_counters.distance_skipped.fetch_add(
                    1, std::memory_order_relaxed);
                continue;
            }
            structural_candidate = true;
            const auto similarity = insertionSimilarity(
                left, right, events[left], events[right], sequences);
            if (similarity) {
                edges.push_back(*similarity);
                cross_anchor_counters.similarity_pairs.fetch_add(
                    1, std::memory_order_relaxed);
            }
        }
    }
    if (structural_candidate) {
        cross_anchor_counters.cigar_candidates.fetch_add(
            1, std::memory_order_relaxed);
    }
    if (edges.empty()) {
        return {};
    }

    std::sort(edges.begin(), edges.end(),
        [&](const auto& left, const auto& right) {
            if (left.paired_bases != right.paired_bases) {
                return left.paired_bases > right.paired_bases;
            }
            if (left.identity != right.identity) {
                return left.identity > right.identity;
            }
            const uint32_t left_coordinate = std::min(
                events[left.left].reference_position,
                events[left.right].reference_position);
            const uint32_t right_coordinate = std::min(
                events[right.left].reference_position,
                events[right.right].reference_position);
            return left_coordinate < right_coordinate;
        });

    ConstrainedDisjointSet components(events);
    std::vector<CrossAnchorSimilarity> accepted_edges;
    for (const auto& edge : edges) {
        if (components.unite(edge.left, edge.right)) {
            accepted_edges.push_back(edge);
        }
    }

    std::map<size_t, CrossAnchorGroup> grouped;
    for (size_t index = 0; index < events.size(); ++index) {
        const size_t root = components.find(index);
        grouped[root].event_indices.push_back(index);
    }
    for (const auto& edge : accepted_edges) {
        const size_t root = components.find(edge.left);
        if (components.find(edge.right) == root) {
            grouped[root].edges.push_back(edge);
        }
    }

    std::vector<CrossAnchorGroup> result;
    for (auto& [root, group] : grouped) {
        (void)root;
        if (group.edges.empty()) {
            continue;
        }
        std::set<ChrName> species;
        group.minimum_position = std::numeric_limits<uint32_t>::max();
        group.maximum_position = 0;
        double identity_sum = 0.0;
        for (const size_t index : group.event_indices) {
            species.insert(events[index].key);
            group.minimum_position = std::min(
                group.minimum_position,
                events[index].reference_position);
            group.maximum_position = std::max(
                group.maximum_position,
                events[index].reference_position);
        }
        if (species.size() < 2 ||
            group.minimum_position == group.maximum_position ||
            group.maximum_position - group.minimum_position >
                kCrossAnchorMaximumDistance) {
            continue;
        }
        for (const auto& edge : group.edges) {
            group.shared_bases += edge.paired_bases;
            identity_sum += edge.identity;
        }
        group.species_count = species.size();
        group.average_identity = identity_sum /
            static_cast<double>(group.edges.size());
        result.push_back(std::move(group));
    }
    return result;
}

std::string buildCrossAnchorCacheKey(
    const std::string& executable,
    const std::map<ChrName, std::string>& alleles) {
    std::string key;
    key.reserve(executable.size() + 64);
    key.append(executable);
    key.push_back('\0');
    for (const auto& [name, sequence] : alleles) {
        key.append(std::to_string(name.size()));
        key.push_back(':');
        key.append(name);
        key.append(std::to_string(sequence.size()));
        key.push_back(':');
        key.append(sequence);
        key.push_back('\0');
    }
    return key;
}

bool alignLocalAlleles(
    const std::map<ChrName, std::string>& raw_alleles,
    size_t insertion_species_count,
    const CrossAnchorRepairConfiguration& configuration,
    std::unordered_map<ChrName, std::string>& aligned) {
    if (raw_alleles.size() < 2) {
        return false;
    }
    aligned.clear();
    for (const auto& [key, sequence] : raw_alleles) {
        aligned.emplace(key, sequence);
    }

    if (insertion_species_count >= 3) {
        if (configuration.executable.empty()) {
            return false;
        }
        const std::string cache_key = buildCrossAnchorCacheKey(
            configuration.executable, raw_alleles);
        {
            std::lock_guard lock(cross_anchor_cache_mutex);
            const auto cached = cross_anchor_msa_cache.find(cache_key);
            if (cached != cross_anchor_msa_cache.end()) {
                aligned.clear();
                for (const auto& [key, sequence] : cached->second) {
                    aligned.emplace(key, sequence);
                }
                cross_anchor_counters.cache_hits.fetch_add(
                    1, std::memory_order_relaxed);
                return true;
            }
        }
        cross_anchor_counters.minipoa_calls.fetch_add(
            1, std::memory_order_relaxed);
        if (!RaMesh::Alignment::ExternalMsaRunner::instance().align(
                configuration.executable, aligned)) {
            return false;
        }
        std::vector<std::pair<ChrName, std::string>> cached_rows;
        cached_rows.reserve(aligned.size());
        for (const auto& [key, sequence] : aligned) {
            cached_rows.emplace_back(key, sequence);
        }
        std::sort(cached_rows.begin(), cached_rows.end());
        {
            std::lock_guard lock(cross_anchor_cache_mutex);
            if (cross_anchor_msa_cache.size() >= kCrossAnchorCacheLimit) {
                cross_anchor_msa_cache.clear();
            }
            cross_anchor_msa_cache.emplace(
                cache_key, std::move(cached_rows));
        }
        return true;
    }

    auto longest = raw_alleles.begin();
    for (auto iterator = std::next(raw_alleles.begin());
         iterator != raw_alleles.end(); ++iterator) {
        if (iterator->second.size() > longest->second.size()) {
            longest = iterator;
        }
    }
    std::unordered_map<ChrName, Cigar_t> local_cigars;
    for (const auto& [key, sequence] : raw_alleles) {
        if (key != longest->first) {
            local_cigars.emplace(
                key,
                globalAlignKSW2Banded(longest->second, sequence));
        }
    }
    try {
        RaMesh::Alignment::mergeReferenceProfile(
            longest->first, aligned, local_cigars);
    } catch (const std::exception&) {
        return false;
    }
    cross_anchor_counters.ksw2_repairs.fetch_add(
        1, std::memory_order_relaxed);
    return true;
}

bool rangesOverlap(const CrossAnchorGroup& left,
                   const CrossAnchorGroup& right) {
    return left.left_reference < right.right_reference &&
        right.left_reference < left.right_reference;
}

struct ReferenceAlignmentLayout {
    std::vector<std::pair<size_t, size_t>> gap_ranges;
    std::vector<size_t> base_columns;
};

std::optional<ReferenceAlignmentLayout> buildReferenceLayout(
    const std::string& aligned_reference) {
    ReferenceAlignmentLayout layout;
    size_t gap_start = 0;
    for (size_t column = 0; column < aligned_reference.size(); ++column) {
        if (aligned_reference[column] == '-') {
            continue;
        }
        layout.gap_ranges.emplace_back(gap_start, column);
        layout.base_columns.push_back(column);
        gap_start = column + 1;
    }
    layout.gap_ranges.emplace_back(gap_start, aligned_reference.size());
    if (layout.gap_ranges.size() != layout.base_columns.size() + 1) {
        return std::nullopt;
    }
    return layout;
}

bool mergeRepairedAndLegacyProfiles(
    const ChrName& reference_name,
    const std::set<ChrName>& participant_keys,
    const std::unordered_map<ChrName, std::string>& repaired,
    const std::unordered_map<ChrName, std::string>& legacy,
    std::unordered_map<ChrName, std::string>& merged) {
    const auto repaired_reference_it = repaired.find(reference_name);
    const auto legacy_reference_it = legacy.find(reference_name);
    if (repaired_reference_it == repaired.end() ||
        legacy_reference_it == legacy.end()) {
        return false;
    }
    const auto repaired_layout = buildReferenceLayout(
        repaired_reference_it->second);
    const auto legacy_layout = buildReferenceLayout(
        legacy_reference_it->second);
    if (!repaired_layout || !legacy_layout ||
        repaired_layout->base_columns.size() !=
            legacy_layout->base_columns.size()) {
        return false;
    }

    for (size_t index = 0;
         index < repaired_layout->base_columns.size(); ++index) {
        const char repaired_base = repaired_reference_it->second[
            repaired_layout->base_columns[index]];
        const char legacy_base = legacy_reference_it->second[
            legacy_layout->base_columns[index]];
        if (normalizedBase(repaired_base) != normalizedBase(legacy_base)) {
            return false;
        }
    }

    merged.clear();
    for (const auto& [key, sequence] : legacy) {
        (void)sequence;
        merged.emplace(key, std::string{});
    }
    for (const auto& [key, sequence] : repaired) {
        (void)sequence;
        if (merged.find(key) == merged.end()) {
            return false;
        }
    }

    const size_t reference_bases = repaired_layout->base_columns.size();
    for (size_t boundary = 0; boundary <= reference_bases; ++boundary) {
        std::vector<size_t> repaired_gap_columns;
        const auto [repaired_begin, repaired_end] =
            repaired_layout->gap_ranges[boundary];
        for (size_t column = repaired_begin;
             column < repaired_end; ++column) {
            const bool occupied = std::any_of(
                participant_keys.begin(), participant_keys.end(),
                [&](const ChrName& key) {
                    const auto row = repaired.find(key);
                    return row != repaired.end() &&
                        column < row->second.size() &&
                        row->second[column] != '-';
                });
            if (occupied) {
                repaired_gap_columns.push_back(column);
            }
        }

        std::vector<size_t> legacy_gap_columns;
        const auto [legacy_begin, legacy_end] =
            legacy_layout->gap_ranges[boundary];
        for (size_t column = legacy_begin;
             column < legacy_end; ++column) {
            bool occupied = false;
            for (const auto& [key, row] : legacy) {
                if (key != reference_name &&
                    participant_keys.find(key) == participant_keys.end() &&
                    column < row.size() && row[column] != '-') {
                    occupied = true;
                    break;
                }
            }
            if (occupied) {
                legacy_gap_columns.push_back(column);
            }
        }

        for (const size_t column : repaired_gap_columns) {
            for (auto& [key, row] : merged) {
                if (key == reference_name ||
                    participant_keys.find(key) == participant_keys.end()) {
                    row.push_back('-');
                } else {
                    row.push_back(repaired.at(key)[column]);
                }
            }
        }
        for (const size_t column : legacy_gap_columns) {
            for (auto& [key, row] : merged) {
                if (key == reference_name ||
                    participant_keys.find(key) != participant_keys.end()) {
                    row.push_back('-');
                } else {
                    row.push_back(legacy.at(key)[column]);
                }
            }
        }

        if (boundary == reference_bases) {
            continue;
        }
        const size_t repaired_column =
            repaired_layout->base_columns[boundary];
        const size_t legacy_column = legacy_layout->base_columns[boundary];
        for (auto& [key, row] : merged) {
            if (key == reference_name ||
                participant_keys.find(key) != participant_keys.end()) {
                row.push_back(repaired.at(key)[repaired_column]);
            } else {
                row.push_back(legacy.at(key)[legacy_column]);
            }
        }
    }
    return true;
}

bool repairCrossAnchorGroups(
    const ChrName& reference_name,
    size_t raw_reference_length,
    const std::vector<CrossAnchorInsertionEvent>& events,
    std::vector<CrossAnchorGroup> groups,
    const CrossAnchorRepairConfiguration& configuration,
    std::unordered_map<ChrName, std::string>& sequences) {
    const auto reference_it = sequences.find(reference_name);
    if (reference_it == sequences.end()) {
        return false;
    }
    const std::string& aligned_reference = reference_it->second;
    std::vector<size_t> reference_boundaries(
        raw_reference_length + 1, aligned_reference.size());
    reference_boundaries[0] = 0;
    size_t reference_position = 0;
    for (size_t column = 0; column < aligned_reference.size(); ++column) {
        if (aligned_reference[column] != '-') {
            if (++reference_position > raw_reference_length) {
                return false;
            }
            reference_boundaries[reference_position] = column + 1;
        }
    }
    if (reference_position != raw_reference_length) {
        return false;
    }

    for (auto& group : groups) {
        group.left_reference = group.minimum_position >
                kCrossAnchorFlankLength
            ? group.minimum_position - kCrossAnchorFlankLength
            : 0;
        const uint64_t requested_right =
            static_cast<uint64_t>(group.maximum_position) +
            kCrossAnchorFlankLength + 1;
        group.right_reference = static_cast<uint32_t>(std::min<uint64_t>(
            raw_reference_length, requested_right));
    }

    std::sort(groups.begin(), groups.end(),
        [](const auto& left, const auto& right) {
            if (left.species_count != right.species_count) {
                return left.species_count > right.species_count;
            }
            if (left.edges.size() != right.edges.size()) {
                return left.edges.size() > right.edges.size();
            }
            if (left.shared_bases != right.shared_bases) {
                return left.shared_bases > right.shared_bases;
            }
            if (left.average_identity != right.average_identity) {
                return left.average_identity > right.average_identity;
            }
            return left.minimum_position < right.minimum_position;
        });

    std::vector<CrossAnchorGroup> selected;
    for (auto& group : groups) {
        const bool overlaps = std::any_of(
            selected.begin(), selected.end(),
            [&](const auto& existing) {
                return rangesOverlap(group, existing);
            });
        if (!overlaps) {
            selected.push_back(std::move(group));
        }
    }
    std::sort(selected.begin(), selected.end(),
        [](const auto& left, const auto& right) {
            return left.left_reference > right.left_reference;
        });

    bool changed = false;
    for (const auto& group : selected) {
        const size_t start_column =
            reference_boundaries[group.left_reference];
        const size_t end_column =
            group.right_reference == raw_reference_length
            ? aligned_reference.size()
            : reference_boundaries[group.right_reference];
        if (start_column >= end_column) {
            cross_anchor_counters.fallback.fetch_add(
                1, std::memory_order_relaxed);
            continue;
        }

        std::unordered_map<ChrName, std::string> old_local;
        std::set<ChrName> insertion_species;
        for (const size_t event_index : group.event_indices) {
            if (event_index < events.size()) {
                insertion_species.insert(events[event_index].key);
            }
        }
        std::set<ChrName> repaired_keys = insertion_species;
        repaired_keys.insert(reference_name);
        std::map<ChrName, std::string> raw_alleles;
        bool invalid = false;
        for (const auto& [key, sequence] : sequences) {
            if (end_column > sequence.size()) {
                invalid = true;
                break;
            }
            const std::string local = sequence.substr(
                start_column, end_column - start_column);
            old_local.emplace(key, local);
            if (repaired_keys.find(key) == repaired_keys.end()) {
                continue;
            }
            std::string allele = removeAlignmentGaps(local);
            if (allele.size() > configuration.maximum_window_span) {
                invalid = true;
                break;
            }
            raw_alleles.emplace(key, std::move(allele));
        }
        if (invalid) {
            cross_anchor_counters.fallback.fetch_add(
                1, std::memory_order_relaxed);
            continue;
        }

        std::unordered_map<ChrName, std::string> repaired_local;
        if (!alignLocalAlleles(
                raw_alleles, insertion_species.size(),
                configuration, repaired_local) ||
            repaired_local.size() != raw_alleles.size()) {
            cross_anchor_counters.fallback.fetch_add(
                1, std::memory_order_relaxed);
            continue;
        }

        std::unordered_map<ChrName, std::string> new_local;
        if (!mergeRepairedAndLegacyProfiles(
                reference_name, insertion_species,
                repaired_local, old_local, new_local)) {
            cross_anchor_counters.fallback.fetch_add(
                1, std::memory_order_relaxed);
            continue;
        }

        size_t new_length = 0;
        for (const auto& [key, old_row] : old_local) {
            const auto aligned_it = new_local.find(key);
            if (aligned_it == new_local.end() ||
                removeAlignmentGaps(aligned_it->second) !=
                    removeAlignmentGaps(old_row) ||
                (new_length != 0 &&
                 aligned_it->second.size() != new_length)) {
                invalid = true;
                break;
            }
            new_length = aligned_it->second.size();
        }
        if (invalid || new_length == 0) {
            cross_anchor_counters.fallback.fetch_add(
                1, std::memory_order_relaxed);
            continue;
        }

        bool strict_improvement = false;
        bool degraded = false;
        for (const auto& edge : group.edges) {
            const auto& left_key = events[edge.left].key;
            const auto& right_key = events[edge.right].key;
            const auto old_metrics = alignedPairMetrics(
                old_local.at(left_key), old_local.at(right_key));
            const auto new_metrics = alignedPairMetrics(
                new_local.at(left_key), new_local.at(right_key));
            const double new_identity = new_metrics.paired_bases == 0
                ? 0.0
                : static_cast<double>(new_metrics.matching_bases) /
                    static_cast<double>(new_metrics.paired_bases);
            const uint32_t left_position =
                events[edge.left].reference_position;
            const uint32_t right_position =
                events[edge.right].reference_position;
            const size_t anchor_distance = left_position > right_position
                ? left_position - right_position
                : right_position - left_position;
            // A d-base reference bridge can force up to d homologous
            // insertion residues onto opposite sides of the bridge while
            // preserving sequence order.  Require the full 10-column gain
            // when possible, but do not make an exactly 10 bp insertion
            // impossible to repair merely because its anchors differ.
            const size_t recoverable_gain = edge.paired_bases >
                    anchor_distance
                ? edge.paired_bases - anchor_distance
                : 1;
            const size_t required_gain = std::min<size_t>(
                kCrossAnchorMinimumInsertion, recoverable_gain);
            degraded = degraded ||
                new_metrics.matching_bases < old_metrics.matching_bases;
            strict_improvement = strict_improvement ||
                (new_metrics.paired_bases >=
                     old_metrics.paired_bases +
                         required_gain &&
                 new_identity >= kCrossAnchorMinimumIdentity);
        }
        if (degraded || !strict_improvement) {
            cross_anchor_counters.fallback.fetch_add(
                1, std::memory_order_relaxed);
            continue;
        }

        for (auto& [key, sequence] : sequences) {
            sequence.replace(
                start_column, end_column - start_column,
                new_local.at(key));
        }
        cross_anchor_counters.accepted.fetch_add(
            1, std::memory_order_relaxed);
        changed = true;
    }
    return changed;
}

}  // namespace

void configureCrossAnchorInsertionRepair(
    const std::string& executable,
    uint_t maximum_window_span) {
    cross_anchor_session.configure(executable, maximum_window_span);
}

void logCrossAnchorInsertionRepairStats() {
    cross_anchor_session.logStats();
}

uint_t mergeAlignmentByRef(
    ChrName ref_name,
    std::unordered_map<ChrName, std::string>& seqs,
    const std::unordered_map<ChrName, Cigar_t>& cigars) {
    const auto started = std::chrono::steady_clock::now();
    cross_anchor_counters.blocks_scanned.fetch_add(
        1, std::memory_order_relaxed);
    const auto reference_it = seqs.find(ref_name);
    if (reference_it == seqs.end()) {
        throw std::invalid_argument("mergeAlignmentByRef: ref not found");
    }
    const size_t raw_reference_length = reference_it->second.size();
    const auto events = collectInsertionEvents(ref_name, seqs, cigars);
    auto groups = buildCrossAnchorGroups(events, seqs);

    const uint_t legacy_length = RaMesh::Alignment::mergeReferenceProfile(
        ref_name, seqs, cigars);
    if (!groups.empty()) {
        const CrossAnchorRepairConfiguration configuration =
            cross_anchor_session.configurationSnapshot();
        repairCrossAnchorGroups(
            ref_name, raw_reference_length, events,
            std::move(groups), configuration, seqs);
    }
    const auto finished = std::chrono::steady_clock::now();
    cross_anchor_counters.nanoseconds.fetch_add(
        elapsedNanoseconds(started, finished),
        std::memory_order_relaxed);
    const auto repaired_reference = seqs.find(ref_name);
    return repaired_reference == seqs.end()
        ? legacy_length
        : static_cast<uint_t>(repaired_reference->second.size());
}

AlignCount countAlignedBases(const Cigar_t& cigar) {
    AlignCount cnt;
    for (auto op : cigar) {
        uint32_t len;
        char type;
        intToCigar(op, type, len);
        switch (type) {
        case 'M': // match or mismatch
        case '=': // match
        case 'X': // mismatch
            cnt.ref_bases += len;
            cnt.query_bases += len;
            break;
        case 'I': // insertion wrt ref
            cnt.query_bases += len;
            break;
        case 'D': // deletion wrt ref
            cnt.ref_bases += len;
            break;
            // 视情况处理 clip/skip
        case 'S': // soft clip
            cnt.query_bases += len;
            break;
        case 'H': // hard clip
            // 不计入
            break;
        case 'N': // skipped region in ref
            cnt.ref_bases += len;
            break;
        default:
            break;
        }
    }
    return cnt;
}
