#ifndef SEQPRO_H
#define SEQPRO_H

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <execution>
#include <filesystem>
#include <fstream>
#include <future>
#include <memory>
#include <shared_mutex>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <variant>

// Cereal序列化库支持
#include "cereal/archives/binary.hpp"
#include "cereal/archives/json.hpp"
#include "cereal/types/memory.hpp"
#include "cereal/types/string.hpp"
#include "cereal/types/unordered_map.hpp"
#include "cereal/types/vector.hpp"

// C++20 span兼容性处理
#if __cplusplus >= 202002L
#include <span>
#else
namespace std {
template <typename T> class span {
public:
  constexpr span(T *data, size_t size) : data_(data), size_(size) {}
  constexpr T *data() const { return data_; }
  constexpr size_t size() const { return size_; }
  constexpr T &operator[](size_t idx) const { return data_[idx]; }
  constexpr T *begin() const { return data_; }
  constexpr T *end() const { return data_ + size_; }

private:
  T *data_;
  size_t size_;
};
} // namespace std
#endif

namespace SeqPro {

// 基础类型定义
using Position = uint64_t;
using Length = uint64_t;
using SequenceId = uint32_t;

// === 异常类 ===
class SeqProException : public std::runtime_error {
public:
  explicit SeqProException(const std::string &message) : std::runtime_error("SeqPro: " + message) {}
};

class FileException : public SeqProException {
public:
  explicit FileException(const std::string &message) : SeqProException("File error: " + message) {}
};

class SequenceException : public SeqProException {
public:
  explicit SequenceException(const std::string &message) : SeqProException("Sequence error: " + message) {}
};

class MaskException : public SeqProException {
public:
  explicit MaskException(const std::string &message) : SeqProException("Mask error: " + message) {}
};

// === 重复区间结构 ===
struct MaskInterval {
  Position start; // 0-based start position
  Position end;   // 0-based end position (exclusive)
  MaskInterval() = default;
  MaskInterval(Position s, Position e) : start(s), end(e) {}

  bool contains(Position pos) const { return pos >= start && pos < end; }
  bool overlaps(const MaskInterval &other) const {
    return start < other.end && end > other.start;
  }
  Length length() const { return end - start; }

  template <class Archive> void serialize(Archive &ar) { ar(start, end); }
};

struct MaskIntervalDelta {
  SequenceId sequence_id{UINT32_MAX};
  std::vector<MaskInterval> intervals;
};

struct MaskBatchMergeStats {
  size_t touched_sequences{0};
  size_t incoming_intervals{0};
  size_t normalized_delta_intervals{0};
  size_t previous_intervals{0};
  size_t final_intervals{0};
  double sort_seconds{0.0};
  double merge_seconds{0.0};
  double metadata_seconds{0.0};
};

struct FinalizedMaskSummary {
  Length masked_bases{0};
  Length masked_sequence_length{0};
  Length separator_count{0};
};

// === 序列信息结构 ===
struct SequenceInfo {
  SequenceId id;
  std::string name;
  Position length;
  Position masked_length;
  Position global_start_pos; // 全局起始位置 (0-based)
  Position masked_global_start_pos; // 遮蔽全局起始位置 (0-based)
  Position file_offset;      // 在文件中的偏移量
  uint32_t line_width;       // FASTA行宽
  uint32_t line_bytes;       // 每行字节数(包括换行符)

  SequenceInfo() = default;
  SequenceInfo(SequenceId seq_id, std::string seq_name, Position len,
               Position g_start_pos, Position m_start_pos, Position offset, uint32_t lw, uint32_t lb)
      : id(seq_id), name(std::move(seq_name)), length(len), masked_length(len),
        global_start_pos(g_start_pos), masked_global_start_pos(m_start_pos + seq_id), file_offset(offset), line_width(lw),
        line_bytes(lb) {}

  template <class Archive> void serialize(Archive &ar) {
    ar(id, name, length, global_start_pos, file_offset, line_width, line_bytes);
  }
};

// === 内存映射文件管理器 ===
class MemoryMapper {
public:
  explicit MemoryMapper(const std::filesystem::path &file_path);
  ~MemoryMapper();

  MemoryMapper(const MemoryMapper &) = delete;
  MemoryMapper &operator=(const MemoryMapper &) = delete;
  MemoryMapper(MemoryMapper &&) noexcept;
  MemoryMapper &operator=(MemoryMapper &&) noexcept;

  size_t getFileSize() const { return file_size_; }
  std::span<const char> getData(Position offset, Length length) const;
  bool isValid() const { return mapped_data_ != nullptr; }

private:
  void *mapped_data_;
  size_t file_size_;
  int fd_;
  void cleanup();
};

// === 序列索引（内部使用） ===
class SequenceIndex {
public:
  static constexpr SequenceId INVALID_ID = UINT32_MAX;

  SequenceIndex() = default;

  void addSequence(const SequenceInfo &info);
  const SequenceInfo *getSequenceInfo(const std::string &name) const;
  const SequenceInfo *getSequenceInfo(SequenceId id) const;
  const SequenceInfo *getSequenceInfoByGlobalPosition(Position global_pos) const;
  SequenceId getSequenceId(const std::string &name) const;
  std::vector<std::string> getSequenceNames() const;
  size_t getSequenceCount() const { return sequences_.size(); }
  Position getTotalGlobalLength() const;
  void clear();

  template <class Archive> void serialize(Archive &ar) {
    ar(sequences_, name_to_id_);
  }

private:
  std::vector<SequenceInfo> sequences_;
  std::unordered_map<std::string, SequenceId> name_to_id_;
};

// === 重复遮蔽管理器（内部使用） ===
class MaskManager {
public:
  MaskManager() = default;

  bool loadFromIntervalFile(const std::filesystem::path &interval_file, bool append = false);
  void addMaskInterval(const std::string &seq_name, const MaskInterval &interval);
  void addMaskInterval(SequenceId seq_id, const MaskInterval &interval);

  bool isMaskPosition(SequenceId seq_id, Position pos) const;
  const std::vector<MaskInterval> &getMaskIntervals(SequenceId seq_id) const;

  Position mapToMaskedPosition(SequenceId seq_id, Position original_pos) const;
  Position mapToOriginalPosition(SequenceId seq_id, Position masked_pos) const;
  Position mapToOriginalPositionSeparated(SequenceId seq_id, Position masked_pos) const;

  Length getMaskedSequenceLength(SequenceId seq_id, Length original_length) const;
  std::vector<std::pair<Position, Position>> getValidRanges(SequenceId seq_id,
                                                           Position start, Position end) const;

  void clear();
  bool hasData() const { return !mask_intervals_.empty(); }
  Length getTotalMaskedBases(SequenceId seq_id) const;
  Length getSeparatorCount(SequenceId seq_id, Length original_length) const;
  SequenceId getOrCreateSequenceId(const std::string &seq_name);


  template <class Archive> void serialize(Archive &ar) {
    ar(mask_intervals_, name_to_id_mapping_);
  }

  // === 批量操作支持 ===
  void addMaskIntervals(SequenceId seq_id, const std::vector<MaskInterval> &intervals);
  void finalizeMaskIntervals(SequenceId seq_id);
  MaskBatchMergeStats mergeFinalizedMaskIntervals(
      SequenceId seq_id, std::vector<MaskInterval> intervals);
  FinalizedMaskSummary summarizeFinalizedIntervals(
      SequenceId seq_id, Length original_length) const;
  void clearMaskIntervals(SequenceId seq_id);

private:
  std::unordered_map<SequenceId, std::vector<MaskInterval>> mask_intervals_;
  std::unordered_map<std::string, SequenceId> name_to_id_mapping_;
  static const std::vector<MaskInterval> empty_intervals_;

  void sortAndMergeIntervals(std::vector<MaskInterval> &intervals);
};

// === 原始序列管理器 ===
class SequenceManager {
public:
  // 构造函数
  explicit SequenceManager(const std::filesystem::path &fasta_path);
  ~SequenceManager() = default;

  // 禁止拷贝，允许移动
  SequenceManager(const SequenceManager &) = delete;
  SequenceManager &operator=(const SequenceManager &) = delete;
  SequenceManager(SequenceManager &&) = default;
  SequenceManager &operator=(SequenceManager &&) = default;

  // === 核心功能：基于原始坐标的序列获取 ===

  // 获取子序列（原始坐标）
  std::string getSubSequence(const std::string &seq_name, Position start, Length length) const;
  std::string getSubSequence(SequenceId seq_id, Position start, Length length) const;

  // Returns a non-owning view only when the requested range is physically
  // contiguous in the mapped FASTA. The view remains valid for this manager's
  // lifetime. Wrapped FASTA ranges that cross a physical line return false.
  bool tryGetContiguousSubSequence(SequenceId seq_id, Position start,
                                   Length length, std::string_view &view) const;

  // Owning fallback that reuses caller-managed storage. Existing string-return
  // APIs keep their original behavior.
  void getSubSequenceInto(SequenceId seq_id, Position start, Length length,
                          std::string &output) const;

  // 全局坐标获取（原始坐标）
  std::string getSubSequenceGlobal(Position global_start, Length length) const;

  // === 序列信息查询 ===

  std::filesystem::path getFastaPath() const { return fasta_path_; }
  SequenceId getSequenceId(const std::string &seq_name) const;
  std::vector<std::string> getSequenceNames() const;
  std::string getSequenceName(const uint32_t& seq_id) const;
  bool hasAmbiguousBasesAll() const { return has_ambiguous_bases; }

  // 序列长度（原始长度）
  Length getSequenceLength(const std::string &seq_name) const;
  Length getSequenceLength(SequenceId seq_id) const;
  size_t getSequenceCount() const;

  // 位置验证（原始坐标）
  bool isValidPosition(const std::string &seq_name, Position pos, Length len = 1) const;
  bool isValidPosition(SequenceId seq_id, Position pos, Length len = 1) const;

  // === 全局坐标系统 ===

  const SequenceInfo* getSequenceInfoFromGlobalPosition(Position global_pos) const;
  std::string getSequenceNameFromGlobalPosition(Position global_pos) const;
  std::pair<SequenceId, Position> globalToLocal(Position global_pos) const;
  Position localToGlobal(const std::string& seq_name, Position local_pos) const;
  Position localToGlobal(SequenceId seq_id, Position local_pos) const;
  Length getTotalLength() const;

  // === 序列拼接 ===

  std::string concatAllSequences(char separator = '\1') const;
  std::string concatSequences(const std::vector<std::string> &seq_names, char separator = '\1') const;
  void streamSequences(std::ostream &output, const std::vector<std::string> &seq_names,
                      char separator = '\1', size_t buffer_size = 1024 * 1024) const;

  // === 批量操作 ===

  struct Query {
    std::string seq_name;
    Position start;
    Length length;
    size_t id = 0;
  };

  struct Result {
    std::string sequence;
    size_t id = 0;
    bool success = true;
    std::string error;
  };

  std::vector<Result> batchQuery(const std::vector<Query> &queries, size_t num_threads = 0) const;

  // === 线程控制 ===

  void setMaxThreads(size_t max_threads);
  size_t getMaxThreads() const { return max_threads_; }

  // === 遮蔽区间管理 ===

  // 清理所有遮蔽区间（对于SequenceManager为空操作）
  void clearMaskedRegions();

  // === 内部访问（供MaskedSequenceManager使用） ===

  const SequenceIndex& getIndex() const { return sequence_index_; }
  const MemoryMapper& getMemoryMapper() const { return *memory_mapper_; }

private:
  std::filesystem::path fasta_path_;
  std::unique_ptr<MemoryMapper> memory_mapper_;
  SequenceIndex sequence_index_;
  bool has_ambiguous_bases = false;

  // 线程安全
  mutable std::shared_mutex mutex_;
  size_t max_threads_;

  // 内部方法
  void buildIndex();
  std::string extractSequence(const SequenceInfo &info, Position start, Length length) const;
  void extractSequenceInto(const SequenceInfo &info, Position start, Length length,
                           std::string &output) const;
};

// === 遮蔽序列管理器 ===
class MaskedSequenceManager {
public:
  // 构造函数
  MaskedSequenceManager(std::unique_ptr<SequenceManager> seq_manager,
                       const std::filesystem::path &mask_file);
  MaskedSequenceManager(std::unique_ptr<SequenceManager> seq_manager);
  ~MaskedSequenceManager() = default;

  // 禁止拷贝
  MaskedSequenceManager(const MaskedSequenceManager &) = delete;
  MaskedSequenceManager &operator=(const MaskedSequenceManager &) = delete;
  MaskedSequenceManager(MaskedSequenceManager &&) = default; // 允许移动构造
  MaskedSequenceManager &operator=(MaskedSequenceManager &&) = default; // 允许移动赋值

  // === 核心功能：基于遮蔽坐标的序列获取 ===

  // 获取子序列（遮蔽坐标）
  std::string getSubSequence(const std::string &seq_name, Position start, Length length) const;
  std::string getSubSequence(SequenceId seq_id, Position start, Length length) const;
  std::string getSubSequenceSeparated(const std::string &seq_name, Position start, Length length, char separator = '\1') const;
  std::string getSubSequenceSeparated(SequenceId seq_id, Position start, Length length, char separator = '\1') const;

  // 全局坐标获取（遮蔽坐标）
  std::string getSubSequenceGlobal(Position global_start, Length length) const;

  // === 序列信息查询（遮蔽后） ===

  // 序列长度（遮蔽后长度）
  Length getSequenceLength(const std::string &seq_name) const;
  Length getSequenceLength(SequenceId seq_id) const;
  Length getTotalLength() const;

  // 计算包含间隔符的总长度
  Length getTotalLengthWithSeparators() const;
  Length getSequenceLengthWithSeparators(const std::string &seq_name) const;
  Length getSequenceLengthWithSeparators(SequenceId seq_id) const;

  // 位置验证（遮蔽坐标）
  bool isValidPosition(const std::string &seq_name, Position pos, Length len = 1) const;
  bool isValidPosition(SequenceId seq_id, Position pos, Length len = 1) const;

  // === 坐标转换（核心功能） ===

  // 遮蔽坐标 -> 原始坐标
  Position toOriginalPosition(const std::string &seq_name, Position masked_pos) const;
  Position toOriginalPosition(SequenceId seq_id, Position masked_pos) const;
  Position toOriginalPositionSeparated(const std::string &seq_name, Position masked_pos) const;
  Position toOriginalPositionSeparated(SequenceId seq_id, Position masked_pos) const;

  // 原始坐标 -> 遮蔽坐标（如果在重复区域内返回INVALID_POSITION）
  static constexpr Position INVALID_POSITION = UINT64_MAX;
  Position toMaskedPosition(const std::string &seq_name, Position original_pos) const;
  Position toMaskedPosition(SequenceId seq_id, Position original_pos) const;

  // 批量坐标转换
  std::vector<Position> toOriginalPositions(const std::string &seq_name,
                                           const std::vector<Position> &masked_positions) const;
  std::vector<Position> toMaskedPositions(const std::string &seq_name,
                                         const std::vector<Position> &original_positions) const;

  // === 全局坐标系统（遮蔽坐标） ===

  std::pair<SequenceId, Position> globalToLocal(Position global_masked_pos) const;
  Position localToGlobal(const std::string& seq_name, Position local_masked_pos) const;
  Position localToGlobal(SequenceId seq_id, Position local_masked_pos) const;

  // 支持间隔符的全局坐标系统
  std::pair<SequenceId, Position> globalToLocalSeparated(Position global_pos_with_separators) const;
  Position localToGlobalSeparated(const std::string& seq_name, Position local_masked_pos) const;
  Position localToGlobalSeparated(SequenceId seq_id, Position local_masked_pos) const;

  // === 重复区域信息 ===

  bool isMaskPosition(const std::string &seq_name, Position original_pos) const;
  bool isMaskPosition(SequenceId seq_id, Position original_pos) const;

  const std::vector<MaskInterval>& getMaskIntervals(const std::string &seq_name) const;
  const std::vector<MaskInterval>& getMaskIntervals(SequenceId seq_id) const;

  Length getMaskedBases(const std::string &seq_name) const;
  Length getMaskedBases(SequenceId seq_id) const;
  Length getTotalMaskedBases() const;

  // 计算间隔符数量
  Length getSeparatorCount(const std::string &seq_name) const;
  Length getSeparatorCount(SequenceId seq_id) const;
  Length getTotalSeparatorCount() const;

  // === 序列拼接（遮蔽后） ===

  std::string concatAllSequences(char separator = '\1') const;
  std::string concatAllSequencesSeparated(char separator = '\1') const;
  std::string concatSequences(const std::vector<std::string> &seq_names, char separator = '\1') const;
  std::string concatSequencesSeparated(const std::vector<std::string> &seq_names, char separator = '\1') const;
  void streamSequences(std::ostream &output, const std::vector<std::string> &seq_names,
                      char separator = '\1', size_t buffer_size = 1024 * 1024) const;

  // === 批量操作 ===

  using Query = SequenceManager::Query;
  using Result = SequenceManager::Result;

  std::vector<Result> batchQuery(const std::vector<Query> &queries, size_t num_threads = 0) const;

  // === 访问原始管理器 ===

  const SequenceManager& getOriginalManager() const { return *original_manager_; }
  const MaskManager& getMaskManager() const { return mask_manager_; }

  // === 实用功能 ===

  // 获取序列名称列表（与原始管理器相同）
  std::vector<std::string> getSequenceNames() const;
  SequenceId getSequenceId(const std::string &seq_name) const;
  size_t getSequenceCount() const;

  // === 批量遮蔽区间管理 ===

  // 批量添加区间（基于遮蔽坐标，内部自动转换为原始坐标）
  void addMaskInterval(const std::string &seq_name, const MaskInterval &interval);
  void addMaskIntervals(const std::string &seq_name, const std::vector<MaskInterval> &intervals);
  void addMaskIntervals(SequenceId seq_id, const std::vector<MaskInterval> &intervals);
  MaskBatchMergeStats applyFinalizedMaskDeltas(
      std::vector<MaskIntervalDelta> deltas);

  // 从文件批量加载区间
  bool loadMaskIntervalsFromFile(const std::filesystem::path &file_path, bool append = false);

  // 定案/优化所有未处理的区间
  void finalizeMaskIntervals();
  void finalizeMaskIntervals(const std::string &seq_name);
  void finalizeMaskIntervals(SequenceId seq_id);

  // 清除特定序列的mask区间
  void clearMaskIntervals(const std::string &seq_name);
  void clearMaskIntervals(SequenceId seq_id);

  // 清除所有序列的遮蔽区间
  void clearMaskedRegions();

  // 状态查询
  bool hasUnfinalizedIntervals() const;
  size_t getUnfinalizedSequenceCount() const;

private:
  std::unique_ptr<SequenceManager> original_manager_;
  MaskManager mask_manager_;

  // 缓存全局坐标映射
  mutable std::unordered_map<SequenceId, Position> global_offset_cache_;
  mutable Position total_masked_length_ = 0;
  mutable bool cache_valid_ = false;

  // 未定案序列的标记
  mutable std::unordered_set<SequenceId> unfinalized_sequences_;

  // 构建缓存
  void buildGlobalOffsetCache() const;
  void ensureCacheValid() const;

  // 内部辅助方法
  const SequenceInfo* getSequenceInfoByMaskedGlobalPosition(Position global_masked_pos) const;

  // 坐标转换和定案辅助函数
  MaskInterval convertMaskedToOriginalInterval(SequenceId seq_id, const MaskInterval &masked_interval) const;
  void ensureFinalized(SequenceId seq_id) const;
  void refreshMaskMetadataAndCache();

  // 间隔符相关的坐标转换辅助函数
  Position convertSeparatedToMaskedPosition(SequenceId seq_id, Position separated_pos) const;
  Position convertMaskedToSeparatedPosition(SequenceId seq_id, Position masked_pos) const;
};


using ManagerVariant = std::variant<std::unique_ptr<SequenceManager>,std::unique_ptr<MaskedSequenceManager>>;
using SharedManagerVariant = std::shared_ptr<ManagerVariant>;
  
// === 工具函数 ===
namespace utils {
  std::filesystem::path cleanFastaFile(const std::filesystem::path &input_fasta,
                                                const std::filesystem::path &output_fasta,
                                                uint64_t line_width);
bool isValidFastaFile(const std::filesystem::path &file_path);
bool isValidMaskFile(const std::filesystem::path &mask_file);
std::string getReadableFileSize(const std::filesystem::path &file_path);
std::string cleanSequenceName(const std::string &name);
bool isValidDNASequence(const std::string &sequence);

} // namespace utils

} // namespace SeqPro

#endif // SEQPRO_H
