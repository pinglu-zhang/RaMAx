#include "softmask_index.h"

#include "kseq.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <limits>
#include <optional>
#include <regex>
#include <stdexcept>
#include <string_view>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

namespace SoftMask {
namespace {

constexpr std::array<char, 8> kMagic{'R', 'A', 'M', 'A', 'X', 'S', 'M', '1'};
constexpr uint32_t kVersion = 1;
constexpr uint32_t kEndianMarker = 0x01020304U;

struct Interval {
    uint64_t start = 0;
    uint64_t end = 0;
};

template <typename T>
void writePod(std::ostream& output, const T& value) {
    output.write(reinterpret_cast<const char*>(&value), sizeof(T));
    if (!output) {
        throw std::runtime_error("Failed to write soft-mask index");
    }
}

template <typename T>
T readPod(const std::byte*& cursor, const std::byte* end,
          const std::filesystem::path& path) {
    if (static_cast<size_t>(end - cursor) < sizeof(T)) {
        throw std::runtime_error("Truncated soft-mask index: " + path.string());
    }
    T value{};
    std::memcpy(&value, cursor, sizeof(T));
    cursor += sizeof(T);
    return value;
}

std::pair<uint64_t, int64_t> sourceMetadata(const std::filesystem::path& path) {
    std::error_code ec;
    const uint64_t size = std::filesystem::file_size(path, ec);
    if (ec) {
        throw std::runtime_error("Cannot stat input FASTA size: " + path.string());
    }
    const auto mtime = std::filesystem::last_write_time(path, ec);
    if (ec) {
        throw std::runtime_error("Cannot stat input FASTA mtime: " + path.string());
    }
    return {size, static_cast<int64_t>(mtime.time_since_epoch().count())};
}

char uppercaseForAlignment(char value) noexcept {
    const auto byte = static_cast<unsigned char>(value);
    const char upper = static_cast<char>(std::toupper(byte));
    switch (upper) {
    case 'A':
    case 'C':
    case 'G':
    case 'T':
        return upper;
    default:
        return 'N';
    }
}

char lowercaseForHal(char value) noexcept {
    switch (value) {
    case 'A': return 'a';
    case 'C': return 'c';
    case 'G': return 'g';
    case 'T': return 't';
    case 'N': return 'n';
    default: return value;
    }
}

void removeIfPresent(const std::filesystem::path& path) noexcept {
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

void publishFile(const std::filesystem::path& partial,
                 const std::filesystem::path& final_path) {
    std::error_code ec;
    std::filesystem::remove(final_path, ec);
    ec.clear();
    std::filesystem::rename(partial, final_path, ec);
    if (ec) {
        throw std::runtime_error("Cannot publish " + final_path.string() + ": " + ec.message());
    }
}

void buildArtifacts(const std::filesystem::path& input_fasta,
                    const std::filesystem::path& fasta_partial,
                    const std::filesystem::path& index_partial) {
    auto [source_size, source_mtime] = sourceMetadata(input_fasta);

    std::ofstream fasta_output(fasta_partial, std::ios::binary | std::ios::trunc);
    std::ofstream index_output(index_partial, std::ios::binary | std::ios::trunc);
    if (!fasta_output) {
        throw std::runtime_error("Cannot create alignment FASTA: " + fasta_partial.string());
    }
    if (!index_output) {
        throw std::runtime_error("Cannot create soft-mask index: " + index_partial.string());
    }

    index_output.write(kMagic.data(), static_cast<std::streamsize>(kMagic.size()));
    writePod(index_output, kVersion);
    writePod(index_output, kEndianMarker);
    writePod(index_output, source_size);
    writePod(index_output, source_mtime);
    const std::streampos sequence_count_pos = index_output.tellp();
    writePod(index_output, uint64_t{0});

    gzFile input = gzopen(input_fasta.c_str(), "r");
    if (input == nullptr) {
        throw std::runtime_error("Cannot open input FASTA: " + input_fasta.string());
    }
    kseq_t* sequence = kseq_init(input);
    if (sequence == nullptr) {
        gzclose(input);
        throw std::runtime_error("Cannot initialize FASTA reader: " + input_fasta.string());
    }

    uint64_t sequence_count = 0;
    try {
        while (kseq_read(sequence) >= 0) {
            const std::string name(sequence->name.s, sequence->name.l);
            const uint64_t length = static_cast<uint64_t>(sequence->seq.l);
            if (name.empty()) {
                throw std::runtime_error("FASTA contains an empty sequence name");
            }
            if (name.size() > std::numeric_limits<uint32_t>::max()) {
                throw std::runtime_error("FASTA sequence name is too long: " + name);
            }

            std::string uppercase;
            uppercase.resize(static_cast<size_t>(length));
            std::vector<Interval> intervals;
            bool in_lowercase_run = false;
            uint64_t run_start = 0;

            for (uint64_t i = 0; i < length; ++i) {
                const char input_base = sequence->seq.s[i];
                uppercase[static_cast<size_t>(i)] = uppercaseForAlignment(input_base);
                const bool lowercase = std::islower(static_cast<unsigned char>(input_base)) != 0;
                if (lowercase && !in_lowercase_run) {
                    in_lowercase_run = true;
                    run_start = i;
                } else if (!lowercase && in_lowercase_run) {
                    intervals.push_back({run_start, i});
                    in_lowercase_run = false;
                }
            }
            if (in_lowercase_run) {
                intervals.push_back({run_start, length});
            }

            fasta_output << '>' << name << '\n' << uppercase << '\n';
            if (!fasta_output) {
                throw std::runtime_error("Failed to write alignment FASTA");
            }

            const uint32_t name_size = static_cast<uint32_t>(name.size());
            writePod(index_output, name_size);
            index_output.write(name.data(), static_cast<std::streamsize>(name.size()));
            writePod(index_output, length);
            writePod(index_output, static_cast<uint64_t>(intervals.size()));
            for (const Interval& interval : intervals) {
                writePod(index_output, interval.start);
                writePod(index_output, interval.end);
            }
            ++sequence_count;
        }
    } catch (...) {
        kseq_destroy(sequence);
        gzclose(input);
        throw;
    }

    kseq_destroy(sequence);
    gzclose(input);

    if (sequence_count == 0) {
        throw std::runtime_error("Input FASTA contains no sequences: " + input_fasta.string());
    }

    index_output.seekp(sequence_count_pos);
    writePod(index_output, sequence_count);
    index_output.seekp(0, std::ios::end);
    fasta_output.flush();
    index_output.flush();
    if (!fasta_output || !index_output) {
        throw std::runtime_error("Failed to finalize soft-mask preprocessing artifacts");
    }
}

void writeCompletionMarker(const std::filesystem::path& marker_partial,
                           const std::filesystem::path& input_fasta,
                           const std::filesystem::path& output_fasta,
                           const std::filesystem::path& output_index) {
    auto [source_size, source_mtime] = sourceMetadata(input_fasta);
    std::ofstream marker(marker_partial, std::ios::binary | std::ios::trunc);
    if (!marker) {
        throw std::runtime_error("Cannot create soft-mask completion marker: " + marker_partial.string());
    }
    const auto fasta_metadata = sourceMetadata(output_fasta);
    const auto index_metadata = sourceMetadata(output_index);
    marker << "{\n"
           << "  \"version\": 2,\n"
           << "  \"source_size\": " << source_size << ",\n"
           << "  \"source_mtime\": " << source_mtime << ",\n"
           << "  \"alignment_fasta_size\": " << std::filesystem::file_size(output_fasta) << ",\n"
           << "  \"alignment_fasta_mtime\": " << fasta_metadata.second << ",\n"
           << "  \"softmask_index_size\": " << std::filesystem::file_size(output_index) << ",\n"
           << "  \"softmask_index_mtime\": " << index_metadata.second << "\n"
           << "}\n";
    marker.flush();
    if (!marker) {
        throw std::runtime_error("Failed to finalize soft-mask completion marker");
    }
}

std::optional<int64_t> markerNumber(const std::filesystem::path& marker,
                                    const std::string& key) {
    std::ifstream input(marker, std::ios::binary);
    if (!input) return std::nullopt;
    const std::string text((std::istreambuf_iterator<char>(input)),
                           std::istreambuf_iterator<char>());
    const std::regex pattern("\\\"" + key +
                             "\\\"\\s*:\\s*(-?[0-9]+)");
    std::smatch match;
    if (!std::regex_search(text, match, pattern)) return std::nullopt;
    try {
        return std::stoll(match[1].str());
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

bool currentMarkerMatches(const std::filesystem::path& marker,
                          const std::filesystem::path& input_fasta,
                          const std::filesystem::path& output_fasta,
                          const std::filesystem::path& output_index) {
    const auto version = markerNumber(marker, "version");
    const auto source_size = markerNumber(marker, "source_size");
    const auto source_mtime = markerNumber(marker, "source_mtime");
    const auto fasta_size = markerNumber(marker, "alignment_fasta_size");
    const auto fasta_mtime = markerNumber(marker, "alignment_fasta_mtime");
    const auto index_size = markerNumber(marker, "softmask_index_size");
    const auto index_mtime = markerNumber(marker, "softmask_index_mtime");
    if (!version || *version != 2 || !source_size || !source_mtime ||
        !fasta_size || !fasta_mtime || !index_size || !index_mtime) {
        return false;
    }
    const auto source = sourceMetadata(input_fasta);
    const auto fasta = sourceMetadata(output_fasta);
    const auto index = sourceMetadata(output_index);
    return *source_size == static_cast<int64_t>(source.first) &&
           *source_mtime == source.second &&
           *fasta_size == static_cast<int64_t>(fasta.first) &&
           *fasta_mtime == fasta.second &&
           *index_size == static_cast<int64_t>(index.first) &&
           *index_mtime == index.second;
}

}  // namespace

Index::Index(const std::filesystem::path& path) : path_(path) {
    fd_ = open(path.c_str(), O_RDONLY);
    if (fd_ < 0) {
        throw std::runtime_error("Cannot open soft-mask index: " + path.string());
    }

    struct stat stat_buffer {};
    if (fstat(fd_, &stat_buffer) != 0 || stat_buffer.st_size <= 0) {
        close(fd_);
        fd_ = -1;
        throw std::runtime_error("Cannot stat soft-mask index: " + path.string());
    }
    mapped_size_ = static_cast<size_t>(stat_buffer.st_size);
    void* mapping = mmap(nullptr, mapped_size_, PROT_READ, MAP_SHARED, fd_, 0);
    if (mapping == MAP_FAILED) {
        close(fd_);
        fd_ = -1;
        throw std::runtime_error("Cannot mmap soft-mask index: " + path.string());
    }
    mapped_data_ = static_cast<const std::byte*>(mapping);

    try {
        const std::byte* cursor = mapped_data_;
        const std::byte* end = mapped_data_ + mapped_size_;
        if (static_cast<size_t>(end - cursor) < kMagic.size() ||
            std::memcmp(cursor, kMagic.data(), kMagic.size()) != 0) {
            throw std::runtime_error("Invalid soft-mask index magic: " + path.string());
        }
        cursor += kMagic.size();
        const uint32_t version = readPod<uint32_t>(cursor, end, path);
        const uint32_t endian = readPod<uint32_t>(cursor, end, path);
        if (version != kVersion || endian != kEndianMarker) {
            throw std::runtime_error("Unsupported soft-mask index format: " + path.string());
        }
        source_size_ = readPod<uint64_t>(cursor, end, path);
        source_mtime_ = readPod<int64_t>(cursor, end, path);
        const uint64_t sequence_count = readPod<uint64_t>(cursor, end, path);

        for (uint64_t sequence_index = 0; sequence_index < sequence_count; ++sequence_index) {
            const uint32_t name_size = readPod<uint32_t>(cursor, end, path);
            if (static_cast<size_t>(end - cursor) < name_size) {
                throw std::runtime_error("Truncated sequence name in soft-mask index: " + path.string());
            }
            std::string name(reinterpret_cast<const char*>(cursor), name_size);
            cursor += name_size;
            SequenceEntry entry;
            entry.length = readPod<uint64_t>(cursor, end, path);
            entry.interval_count = readPod<uint64_t>(cursor, end, path);
            if (entry.interval_count > static_cast<uint64_t>(end - cursor) / (2 * sizeof(uint64_t))) {
                throw std::runtime_error("Truncated intervals in soft-mask index: " + path.string());
            }
            entry.interval_data = cursor;

            uint64_t previous_end = 0;
            for (uint64_t i = 0; i < entry.interval_count; ++i) {
                auto [start, finish] = intervalAt(entry, i);
                if (start >= finish || finish > entry.length || (i > 0 && start < previous_end)) {
                    throw std::runtime_error("Invalid interval ordering in soft-mask index: " + path.string());
                }
                previous_end = finish;
            }
            cursor += entry.interval_count * 2 * sizeof(uint64_t);
            if (!sequences_.emplace(std::move(name), entry).second) {
                throw std::runtime_error("Duplicate sequence in soft-mask index: " + path.string());
            }
        }
        if (cursor != end) {
            throw std::runtime_error("Trailing bytes in soft-mask index: " + path.string());
        }
    } catch (...) {
        munmap(const_cast<std::byte*>(mapped_data_), mapped_size_);
        close(fd_);
        mapped_data_ = nullptr;
        mapped_size_ = 0;
        fd_ = -1;
        throw;
    }
}

Index::~Index() {
    if (mapped_data_ != nullptr) {
        munmap(const_cast<std::byte*>(mapped_data_), mapped_size_);
    }
    if (fd_ >= 0) {
        close(fd_);
    }
}

std::pair<uint64_t, uint64_t>
Index::intervalAt(const SequenceEntry& entry, uint64_t index) const {
    if (index >= entry.interval_count) {
        throw std::out_of_range("Soft-mask interval index out of range");
    }
    const std::byte* position = entry.interval_data + index * 2 * sizeof(uint64_t);
    uint64_t start = 0;
    uint64_t end = 0;
    std::memcpy(&start, position, sizeof(uint64_t));
    std::memcpy(&end, position + sizeof(uint64_t), sizeof(uint64_t));
    return {start, end};
}

void Index::restore(const std::string& sequence_name,
                    uint64_t forward_start,
                    std::string& dna) const {
    const auto found = sequences_.find(sequence_name);
    if (found == sequences_.end()) {
        throw std::runtime_error("Sequence missing from soft-mask index: " + sequence_name);
    }
    const SequenceEntry& entry = found->second;
    if (forward_start > entry.length || dna.size() > entry.length - forward_start) {
        throw std::runtime_error("Soft-mask restore range is outside sequence: " + sequence_name);
    }
    if (dna.empty() || entry.interval_count == 0) {
        return;
    }

    uint64_t low = 0;
    uint64_t high = entry.interval_count;
    while (low < high) {
        const uint64_t middle = low + (high - low) / 2;
        const auto [unused_start, interval_end] = intervalAt(entry, middle);
        (void)unused_start;
        if (interval_end <= forward_start) {
            low = middle + 1;
        } else {
            high = middle;
        }
    }

    const uint64_t requested_end = forward_start + static_cast<uint64_t>(dna.size());
    for (uint64_t i = low; i < entry.interval_count; ++i) {
        const auto [interval_start, interval_end] = intervalAt(entry, i);
        if (interval_start >= requested_end) {
            break;
        }
        const uint64_t overlap_start = std::max(interval_start, forward_start);
        const uint64_t overlap_end = std::min(interval_end, requested_end);
        for (uint64_t position = overlap_start; position < overlap_end; ++position) {
            const size_t local = static_cast<size_t>(position - forward_start);
            dna[local] = lowercaseForHal(dna[local]);
        }
    }
}

bool Index::hasSequence(const std::string& sequence_name) const {
    return sequences_.contains(sequence_name);
}

uint64_t Index::sequenceLength(const std::string& sequence_name) const {
    const auto found = sequences_.find(sequence_name);
    if (found == sequences_.end()) {
        throw std::runtime_error("Sequence missing from soft-mask index: " + sequence_name);
    }
    return found->second.length;
}

uint64_t Index::intervalCount(const std::string& sequence_name) const {
    const auto found = sequences_.find(sequence_name);
    if (found == sequences_.end()) {
        throw std::runtime_error("Sequence missing from soft-mask index: " + sequence_name);
    }
    return found->second.interval_count;
}

void AncestorBaseVote::add(char observation) noexcept {
    static constexpr std::array<char, 4> bases{'A', 'C', 'G', 'T'};
    const auto raw = static_cast<unsigned char>(observation);
    const char base = static_cast<char>(std::toupper(raw));
    const bool lowercase = std::islower(raw) != 0;
    for (size_t index = 0; index < bases.size(); ++index) {
        if (base == bases[index]) {
            ++base_counts_[index];
            if (lowercase) ++lowercase_counts_[index];
            return;
        }
    }
    if (base == 'N') {
        ++n_count_;
        if (lowercase) ++lowercase_n_count_;
    }
}

char AncestorBaseVote::result() const noexcept {
    static constexpr std::array<char, 4> bases{'A', 'C', 'G', 'T'};
    char best_base = 'N';
    uint32_t best_count = 0;
    uint32_t best_lowercase_count = 0;
    for (size_t index = 0; index < bases.size(); ++index) {
        if (base_counts_[index] > best_count) {
            best_count = base_counts_[index];
            best_lowercase_count = lowercase_counts_[index];
            best_base = bases[index];
        }
    }
    if (best_count > 0) {
        if (best_lowercase_count * 2 >= best_count) {
            return static_cast<char>(std::tolower(static_cast<unsigned char>(best_base)));
        }
        return best_base;
    }
    if (n_count_ > 0 && lowercase_n_count_ * 2 >= n_count_) {
        return 'n';
    }
    return 'N';
}

bool ensureUppercaseFastaAndIndex(const std::filesystem::path& input_fasta,
                                  const std::filesystem::path& output_fasta,
                                  const std::filesystem::path& output_index,
                                  const std::filesystem::path& completion_marker,
                                  bool trust_legacy_cache) {
    std::filesystem::create_directories(output_fasta.parent_path());
    const auto [expected_size, expected_mtime] = sourceMetadata(input_fasta);

    if (std::filesystem::is_regular_file(output_fasta) &&
        std::filesystem::is_regular_file(output_index) &&
        std::filesystem::is_regular_file(completion_marker)) {
        try {
            Index existing(output_index);
            if (existing.sourceSize() == expected_size &&
                existing.sourceMtime() == expected_mtime &&
                currentMarkerMatches(completion_marker, input_fasta,
                                     output_fasta, output_index)) {
                return true;
            }
        } catch (const std::exception&) {
            // Rebuild exact versioned intermediates below.
        }
    }

    if (trust_legacy_cache &&
        std::filesystem::is_regular_file(output_fasta) &&
        std::filesystem::is_regular_file(output_index)) {
        try {
            Index existing(output_index);
            if (existing.sourceSize() == expected_size &&
                existing.sourceMtime() == expected_mtime) {
                std::filesystem::path marker_partial = completion_marker;
                marker_partial += ".partial";
                removeIfPresent(marker_partial);
                writeCompletionMarker(marker_partial, input_fasta,
                                      output_fasta, output_index);
                publishFile(marker_partial, completion_marker);
                return true;
            }
        } catch (const std::exception&) {
            // A legacy artifact that cannot be loaded is rebuilt below.
        }
    }

    std::filesystem::path fasta_partial = output_fasta;
    fasta_partial += ".partial";
    std::filesystem::path index_partial = output_index;
    index_partial += ".partial";
    std::filesystem::path marker_partial = completion_marker;
    marker_partial += ".partial";
    removeIfPresent(fasta_partial);
    removeIfPresent(index_partial);
    removeIfPresent(marker_partial);
    removeIfPresent(completion_marker);

    try {
        buildArtifacts(input_fasta, fasta_partial, index_partial);
        // Parse and validate the complete sidecar before publishing it.
        {
            Index validated(index_partial);
            if (validated.sourceSize() != expected_size || validated.sourceMtime() != expected_mtime) {
                throw std::runtime_error("Soft-mask source metadata changed during preprocessing");
            }
        }
        // Close the validation mmap before rename. This is required by some
        // mounted filesystems and is harmless on native Linux filesystems.
        publishFile(fasta_partial, output_fasta);
        publishFile(index_partial, output_index);
        writeCompletionMarker(marker_partial, input_fasta, output_fasta, output_index);
        publishFile(marker_partial, completion_marker);
    } catch (...) {
        removeIfPresent(fasta_partial);
        removeIfPresent(index_partial);
        removeIfPresent(marker_partial);
        throw;
    }
    return false;
}

IndexMap loadIndexes(const PathMap& paths) {
    IndexMap indexes;
    for (const auto& [species, path] : paths) {
        indexes.emplace(species, std::make_shared<const Index>(path));
    }
    return indexes;
}

}  // namespace SoftMask
