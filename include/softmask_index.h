#pragma once

#include <cstddef>
#include <cstdint>
#include <array>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace SoftMask {

using PathMap = std::map<std::string, std::filesystem::path>;

class Index {
public:
    explicit Index(const std::filesystem::path& path);
    ~Index();

    Index(const Index&) = delete;
    Index& operator=(const Index&) = delete;

    void restore(const std::string& sequence_name,
                 uint64_t forward_start,
                 std::string& dna) const;

    [[nodiscard]] bool hasSequence(const std::string& sequence_name) const;
    [[nodiscard]] uint64_t sequenceLength(const std::string& sequence_name) const;
    [[nodiscard]] uint64_t intervalCount(const std::string& sequence_name) const;
    [[nodiscard]] std::vector<std::pair<uint64_t, uint64_t>>
    intervals(const std::string& sequence_name) const;
    [[nodiscard]] uint64_t sourceSize() const noexcept { return source_size_; }
    [[nodiscard]] int64_t sourceMtime() const noexcept { return source_mtime_; }
    [[nodiscard]] const std::filesystem::path& path() const noexcept { return path_; }

private:
    struct SequenceEntry {
        uint64_t length = 0;
        uint64_t interval_count = 0;
        const std::byte* interval_data = nullptr;
    };

    std::filesystem::path path_;
    int fd_ = -1;
    const std::byte* mapped_data_ = nullptr;
    size_t mapped_size_ = 0;
    uint64_t source_size_ = 0;
    int64_t source_mtime_ = 0;
    std::unordered_map<std::string, SequenceEntry> sequences_;

    [[nodiscard]] std::pair<uint64_t, uint64_t>
    intervalAt(const SequenceEntry& entry, uint64_t index) const;
};

using IndexMap = std::map<std::string, std::shared_ptr<const Index>>;

/** Fixed-size, allocation-free vote used only while reconstructing HAL DNA. */
class AncestorBaseVote {
public:
    void add(char observation) noexcept;
    [[nodiscard]] char result() const noexcept;

private:
    std::array<uint32_t, 4> base_counts_{};
    std::array<uint32_t, 4> lowercase_counts_{};
    uint32_t n_count_ = 0;
    uint32_t lowercase_n_count_ = 0;
};

/**
 * Create the exact all-uppercase FASTA used by the existing alignment path,
 * plus a compact sidecar containing lowercase runs from the input FASTA.
 * Final artifacts are published atomically and are reused only when their
 * completion marker and source metadata match.
 */
void ensureUppercaseFastaAndIndex(const std::filesystem::path& input_fasta,
                                  const std::filesystem::path& output_fasta,
                                  const std::filesystem::path& output_index,
                                  const std::filesystem::path& completion_marker);

/** Load all species indexes read-only. This is intended to run at HAL export. */
IndexMap loadIndexes(const PathMap& paths);

}  // namespace SoftMask
