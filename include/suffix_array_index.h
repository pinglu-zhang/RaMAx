#ifndef SUFFIX_ARRAY_INDEX_H
#define SUFFIX_ARRAY_INDEX_H

#include "index.h"

#include <cstdint>
#include <string>
#include <vector>

#define SAINDEX_EXTENSION "saidx"

// Enhanced suffix-array anchor index.  Its public search contract mirrors
// FM_Index so PairRareAligner and every downstream stage can remain unchanged.
class Suffix_Array_Index {
public:
    Suffix_Array_Index(SpeciesName species_name,
        SeqPro::ManagerVariant& fasta_manager,
        uint_t sampling_rate = 1);

    bool buildIndex(FilePath output_path, bool fast_mode, uint_t thread_count);

    MatchVec2DPtr findAnchors(ChrIndex query_chr_index, std::string query,
        SearchMode search_mode, Strand strand, bool allow_MEM,
        uint_t query_offset, uint_t min_anchor_length,
        bool allow_short_mum, uint_t max_anchor_frequency,
        sdsl::int_vector<0>& ref_global_cache,
        SeqPro::Length sampling_interval,
        uint_t accurate_skip_threshold = 0) const;

    MatchVec2DPtr findAnchorsFast(ChrIndex query_chr_index,
        std::string query, Strand strand, bool allow_MEM,
        uint_t query_offset, uint_t min_anchor_length,
        bool allow_short_mum, uint_t max_anchor_frequency,
        sdsl::int_vector<0>& ref_global_cache,
        SeqPro::Length sampling_interval) const;

    MatchVec2DPtr findAnchorsMiddle(ChrIndex query_chr_index,
        std::string query, Strand strand, bool allow_MEM,
        uint_t query_offset, uint_t min_anchor_length,
        bool allow_short_mum, uint_t max_anchor_frequency,
        sdsl::int_vector<0>& ref_global_cache,
        SeqPro::Length sampling_interval) const;

    MatchVec2DPtr findAnchorsAccurate(ChrIndex query_chr_index,
        std::string query, Strand strand, bool allow_MEM,
        uint_t query_offset, uint_t min_anchor_length,
        bool allow_short_mum, uint_t max_anchor_frequency,
        sdsl::int_vector<0>& ref_global_cache,
        SeqPro::Length sampling_interval,
        uint_t accurate_skip_threshold) const;

    struct MUMInfo {
        uint_t pos;
        uint_t len;
        bool is_mum;
    };

    void bisectAnchors(const std::string& query,
        ChrIndex query_chr_index, Strand strand, bool allow_MEM,
        uint_t query_offset, uint_t query_length,
        uint_t min_anchor_length, bool allow_short_mum,
        uint_t max_anchor_frequency, const MUMInfo& left,
        const MUMInfo& right, MatchVec2D& out,
        sdsl::int_vector<0>& ref_global_cache,
        SeqPro::Length sampling_interval) const;

    uint_t findSubSeqAnchors(const char* query, uint_t query_length,
        bool allow_MEM, RegionVec& region_vec,
        uint_t min_anchor_length, bool allow_short_mum,
        uint_t max_anchor_frequency,
        sdsl::int_vector<0>& ref_global_cache,
        SeqPro::Length sampling_interval) const;

    bool saveToFile(const std::string& filename) const;
    bool loadFromFile(const std::string& filename);

    uint_t samplingRate() const noexcept { return sampling_rate; }
    uint_t textSize() const noexcept { return text_size; }
    uint_t storedSuffixCount() const noexcept {
        return stored_suffix_count;
    }

    SpeciesName species_name;
    SeqPro::ManagerVariant& fasta_manager;

private:
    struct SearchCursor {
        SAInterval interval{};
        uint_t raw_match_length{0};
        bool reusable{false};
    };

    struct PrefixMatch {
        SAInterval interval{};
        uint_t length{0};
        uint_t frequency{0};
        bool interval_complete{false};
    };

    MatchVec2DPtr findAnchorsImpl(ChrIndex query_chr_index,
        std::string query, SearchMode search_mode, Strand strand,
        bool allow_MEM, uint_t query_offset,
        uint_t min_anchor_length, bool allow_short_mum,
        uint_t max_anchor_frequency,
        sdsl::int_vector<0>& ref_global_cache,
        SeqPro::Length sampling_interval,
        uint_t accurate_skip_threshold) const;

    uint_t findSubSeqAnchorsWithCursor(const char* query,
        uint_t query_length, uint_t shift_from_previous,
        SearchCursor* cursor, bool allow_MEM, RegionVec& region_vec,
        uint_t min_anchor_length, bool allow_short_mum,
        uint_t max_anchor_frequency,
        sdsl::int_vector<0>& ref_global_cache,
        SeqPro::Length sampling_interval) const;

    PrefixMatch longestPrefix(const char* query, uint_t query_length,
        SAInterval search_range, uint_t known_depth,
        uint_t accepted_frequency_limit) const;
    SAInterval suffixLinkInterval(SAInterval previous, uint_t depth,
        uint_t shift) const;
    void buildPrefixDirectory();

    uint64_t suffixAt(uint64_t row) const noexcept;
    uint64_t inverseAt(uint64_t position) const noexcept;
    uint64_t lcpAt(uint64_t row) const noexcept;
    bool isExcludedReversePosition(uint64_t position) const noexcept;
    bool isExcludedRow(uint64_t row) const noexcept;
    bool originalSuffixLess(uint_t left, uint_t right,
        uint_t shared_prefix) const noexcept;
    void appendRegion(uint_t ref_global_pos, uint_t match_length,
        RegionVec& region_vec,
        sdsl::int_vector<0>& ref_global_cache,
        SeqPro::Length sampling_interval) const;

    std::string reverse_text;
    bool coordinates_are_64_bit{false};
    uint_t text_size{0};
    uint_t stored_suffix_count{0};
    uint_t sampling_rate{1};
    std::vector<uint32_t> suffix_array_32;
    std::vector<uint32_t> inverse_suffix_array_32;
    std::vector<uint32_t> lcp_32;
    std::vector<uint64_t> suffix_array_64;
    std::vector<uint64_t> inverse_suffix_array_64;
    std::vector<uint64_t> lcp_64;
    uint_t total_size{0};
    std::vector<uint64_t> excluded_reverse_positions;
    std::vector<SAInterval> prefix_directory;
};

#endif
