#include "SeqPro.h"
#include <algorithm>
#include <cctype>
#include <fcntl.h>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <sys/mman.h>
#include <unistd.h>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"

// 初始化 kseq，用于读取 gzFile 类型的文件
KSEQ_INIT(gzFile, gzread)

namespace SeqPro {

// ========================
// MaskManager 实现
// ========================

const std::vector<MaskInterval> MaskManager::empty_intervals_;

namespace {

size_t compactSortedMaskIntervals(std::vector<MaskInterval>& intervals) {
  if (intervals.empty()) {
    return 0;
  }

  size_t write = 0;
  for (size_t read = 1; read < intervals.size(); ++read) {
    if (intervals[write].end >= intervals[read].start) {
      intervals[write].end =
          std::max(intervals[write].end, intervals[read].end);
    } else {
      ++write;
      if (write != read) {
        intervals[write] = std::move(intervals[read]);
      }
    }
  }
  intervals.resize(write + 1);
  return intervals.size();
}

}  // namespace

bool MaskManager::loadFromIntervalFile(const std::filesystem::path &interval_file, bool append) {
  std::ifstream file(interval_file);
  if (!file.is_open()) {
    return false;
  }

  if (!append) {
    clear();
  }

  std::string line;
  std::string current_seq_name;
  SequenceId current_seq_id = 0;

  while (std::getline(file, line)) {
    if (line.empty()) continue;

    if (line[0] == '>') {
      // 序列头部
      current_seq_name = line.substr(1);
      // 清理序列名
      size_t space_pos = current_seq_name.find_first_of(" \t");
      if (space_pos != std::string::npos) {
        current_seq_name = current_seq_name.substr(0, space_pos);
      }
      current_seq_id = getOrCreateSequenceId(current_seq_name);
    } else {
      // 解析区间行
      std::istringstream iss(line);
      std::string token;
      while (iss >> token) {
        size_t dash_pos = token.find('-');
        if (dash_pos != std::string::npos) {
          try {
            // 转换为0-based坐标
            Position start = std::stoull(token.substr(0, dash_pos)) - 1;
            Position end = std::stoull(token.substr(dash_pos + 1));
            if (start < end) {
              addMaskInterval(current_seq_id, MaskInterval(start, end));
            }
          } catch (const std::exception &) {
            // 忽略解析错误
          }
        }
      }
    }
  }

  // 排序并合并重叠区间
  for (auto &pair : mask_intervals_) {
    sortAndMergeIntervals(pair.second);
  }

  return true;
}

void MaskManager::addMaskInterval(const std::string &seq_name, const MaskInterval &interval) {
  SequenceId seq_id = getOrCreateSequenceId(seq_name);
  addMaskInterval(seq_id, interval);
}

void MaskManager::addMaskInterval(SequenceId seq_id, const MaskInterval &interval) {
  mask_intervals_[seq_id].push_back(interval);
}

bool MaskManager::isMaskPosition(SequenceId seq_id, Position pos) const {
  auto it = mask_intervals_.find(seq_id);
  if (it == mask_intervals_.end()) {
    return false;
  }

  const auto &intervals = it->second;
  // 二分查找
  auto comp = [](const MaskInterval &interval, Position pos) {
    return interval.end <= pos;
  };
  auto lower = std::lower_bound(intervals.begin(), intervals.end(), pos, comp);
  return lower != intervals.end() && lower->contains(pos);
}

const std::vector<MaskInterval> &MaskManager::getMaskIntervals(SequenceId seq_id) const {
  auto it = mask_intervals_.find(seq_id);
  return (it != mask_intervals_.end()) ? it->second : empty_intervals_;
}

Position MaskManager::mapToMaskedPosition(SequenceId seq_id, Position original_pos) const {
  auto it = mask_intervals_.find(seq_id);
  if (it == mask_intervals_.end()) {
    return original_pos;
  }

  const auto &intervals = it->second;
  Position masked_pos = original_pos;
  
  for (const auto &interval : intervals) {
    if (interval.start > original_pos) {
      break;
    }
    if (interval.contains(original_pos)) {
      // 位置在重复区域内
      return MaskedSequenceManager::INVALID_POSITION;
    }
    if (interval.end <= original_pos) {
      masked_pos -= interval.length();
    }
  }

  return masked_pos;
}

Position MaskManager::mapToOriginalPosition(SequenceId seq_id, Position masked_pos) const {
  auto it = mask_intervals_.find(seq_id);
  if (it == mask_intervals_.end()) {
    return masked_pos;
  }

  const auto &intervals = it->second;
  Position original_pos = masked_pos;
  Position accumulated_masked = 0;

  for (const auto &interval : intervals) {
    Position interval_start_in_masked = interval.start - accumulated_masked;
    
    if (masked_pos < interval_start_in_masked) {
      // 目标位置在当前区间之前
      break;
    }
    
    // 累加跳过的遮蔽碱基
    original_pos += interval.length();
    accumulated_masked += interval.length();
  }

  return original_pos;
}

Position MaskManager::mapToOriginalPositionSeparated(SequenceId seq_id, Position masked_pos) const {
  auto it = mask_intervals_.find(seq_id);
  if (it == mask_intervals_.end()) {
    return masked_pos;
  }
  // TODO 性能可以优化
  // TODO 考虑intetval首个不是0的情况

  const auto &intervals = it->second;

  Position origin_pos = masked_pos;
  
  long long accumulated_unmasked = 0;
  if (intervals.size() > 0 && intervals[0].start == 0) {
      accumulated_unmasked = -1;
  }


  Position last_unmasked_interval_end = 0;

  Position cur_interval_len = 0;


  for (const auto &interval : intervals) {
      cur_interval_len = interval.start - last_unmasked_interval_end + 1;
      accumulated_unmasked += cur_interval_len;
      if (static_cast<long long>(origin_pos) < accumulated_unmasked) {
          accumulated_unmasked -= cur_interval_len;
          break;
      }
      //accumulated_masked += cur_interval_length;
	  last_unmasked_interval_end = interval.end;
  }

  // 返回时减去前面遮蔽区间的数量
  return origin_pos - accumulated_unmasked + last_unmasked_interval_end;
}



Length MaskManager::getMaskedSequenceLength(SequenceId seq_id, Length original_length) const {
  auto it = mask_intervals_.find(seq_id);
  if (it == mask_intervals_.end()) {
    return original_length;
  }

  const auto &intervals = it->second;
  Length masked_bases = 0;
  
  for (const auto &interval : intervals) {
    if (interval.start >= original_length) {
      break;
    }
    Position interval_end = std::min(interval.end, original_length);
    masked_bases += (interval_end - interval.start);
  }
  
  return original_length - masked_bases;
}

std::vector<std::pair<Position, Position>> MaskManager::getValidRanges(
    SequenceId seq_id, Position start, Position end) const {
  
  auto it = mask_intervals_.find(seq_id);
  if (it == mask_intervals_.end()) {
    return {{start, end}};
  }

  const auto &intervals = it->second;
  std::vector<std::pair<Position, Position>> valid_ranges;
  Position current_pos = start;

  for (const auto &interval : intervals) {
    if (interval.start >= end) {
      break;
    }
    if (interval.end <= start) {
      continue;
    }

    if (current_pos < interval.start && interval.start < end) {
      valid_ranges.emplace_back(current_pos, std::min(interval.start, end));
    }

    current_pos = std::max(current_pos, interval.end);
  }

  if (current_pos < end) {
    valid_ranges.emplace_back(current_pos, end);
  }

  return valid_ranges;
}

void MaskManager::clear() {
  mask_intervals_.clear();
  name_to_id_mapping_.clear();
}

Length MaskManager::getTotalMaskedBases(SequenceId seq_id) const {
  auto it = mask_intervals_.find(seq_id);
  if (it == mask_intervals_.end()) {
    return 0;
  }

  const auto &intervals = it->second;
  Length total = 0;
  for (const auto &interval : intervals) {
    total += interval.length();
  }
  return total;
}

  Length MaskManager::getSeparatorCount(SequenceId seq_id, Length original_length) const {
  // 如果序列本身没有长度，则不可能有任何区间或间隔符。
  if (original_length == 0) {
    return 0;
  }

  auto it = mask_intervals_.find(seq_id);

  // 如果没有遮蔽区间，则整个序列是一个大的有效区间。
  // 一个有效区间对应一个间隔符。
  if (it == mask_intervals_.end() || it->second.empty()) {
    return 0;
  }

  const auto &intervals = it->second;
  Length separator_count = 0;
  Position current_pos = 0; // 从序列的起始位置开始
  const Position end = original_length; // 序列的结束位置

  // 遍历所有已排序并合并的遮蔽区间
  for (const auto &interval : intervals) {
    // 如果当前位置和下一个遮蔽区间的起点之间存在间隙，
    // 这就是一个有效区间，需要计入一个间隔符。
    if (current_pos < interval.start) {
      separator_count++;
    }

    // 将当前位置“跳”到遮蔽区间的末尾之后，
    // 以便寻找下一个有效区间。
    current_pos = std::max(current_pos, interval.end);
  }

  // 循环结束后，如果当前位置还未到达序列末尾，
  // 说明最后一个遮蔽区间到序列末尾之间还有一个有效区间。
  if (current_pos < end) {
    separator_count++;
  }

  return separator_count;
}

SequenceId MaskManager::getOrCreateSequenceId(const std::string &seq_name) {
  auto it = name_to_id_mapping_.find(seq_name);
  if (it != name_to_id_mapping_.end()) {
    return it->second;
  }

  SequenceId new_id = static_cast<SequenceId>(name_to_id_mapping_.size());
  name_to_id_mapping_[seq_name] = new_id;
  return new_id;
}

void MaskManager::sortAndMergeIntervals(std::vector<MaskInterval> &intervals) {
  if (intervals.empty()) return;

  std::sort(intervals.begin(), intervals.end(),
            [](const MaskInterval &a, const MaskInterval &b) {
              return a.start < b.start;
            });
  compactSortedMaskIntervals(intervals);
}

// === MaskManager 批量操作支持实现 ===

void MaskManager::addMaskIntervals(SequenceId seq_id, const std::vector<MaskInterval> &intervals) {
  if (intervals.empty()) return;
  
  auto& seq_intervals = mask_intervals_[seq_id];
  seq_intervals.reserve(seq_intervals.size() + intervals.size());
  
  for (const auto& interval : intervals) {
    seq_intervals.push_back(interval);
  }
}

void MaskManager::finalizeMaskIntervals(SequenceId seq_id) {
  auto it = mask_intervals_.find(seq_id);
  if (it != mask_intervals_.end()) {
    sortAndMergeIntervals(it->second);
  }
}

MaskBatchMergeStats MaskManager::mergeFinalizedMaskIntervals(
    SequenceId seq_id, std::vector<MaskInterval> intervals) {
  MaskBatchMergeStats stats;
  stats.incoming_intervals = intervals.size();

  const auto existing_it = mask_intervals_.find(seq_id);
  stats.previous_intervals =
      existing_it == mask_intervals_.end() ? 0 : existing_it->second.size();
  stats.final_intervals = stats.previous_intervals;
  if (intervals.empty()) {
    return stats;
  }

  stats.touched_sequences = 1;
  const auto sort_started = std::chrono::steady_clock::now();
  std::sort(intervals.begin(), intervals.end(),
            [](const MaskInterval& left, const MaskInterval& right) {
              return left.start < right.start;
            });
  compactSortedMaskIntervals(intervals);
  stats.sort_seconds =
      std::chrono::duration<double>(
          std::chrono::steady_clock::now() - sort_started)
          .count();
  stats.normalized_delta_intervals = intervals.size();

  const auto merge_started = std::chrono::steady_clock::now();
  auto& existing = mask_intervals_[seq_id];
  if (existing.empty()) {
    existing = std::move(intervals);
  } else {
    const size_t existing_size = existing.size();
    const size_t delta_size = intervals.size();
    if (delta_size > std::numeric_limits<size_t>::max() - existing_size) {
      throw std::length_error("Mask interval merge exceeds vector capacity");
    }

    existing.resize(existing_size + delta_size);
    size_t existing_index = existing_size;
    size_t delta_index = delta_size;
    size_t write_index = existing_size + delta_size;

    while (existing_index > 0 && delta_index > 0) {
      if (existing[existing_index - 1].start >
          intervals[delta_index - 1].start) {
        existing[--write_index] =
            std::move(existing[--existing_index]);
      } else {
        existing[--write_index] =
            std::move(intervals[--delta_index]);
      }
    }
    while (delta_index > 0) {
      existing[--write_index] =
          std::move(intervals[--delta_index]);
    }
    while (existing_index > 0) {
      existing[--write_index] =
          std::move(existing[--existing_index]);
    }

    compactSortedMaskIntervals(existing);
  }

  stats.final_intervals = existing.size();
  stats.merge_seconds =
      std::chrono::duration<double>(
          std::chrono::steady_clock::now() - merge_started)
          .count();
  return stats;
}

FinalizedMaskSummary MaskManager::summarizeFinalizedIntervals(
    SequenceId seq_id, Length original_length) const {
  FinalizedMaskSummary summary;
  summary.masked_sequence_length = original_length;
  if (original_length == 0) {
    return summary;
  }

  const auto it = mask_intervals_.find(seq_id);
  if (it == mask_intervals_.end() || it->second.empty()) {
    return summary;
  }

  Position current_pos = 0;
  bool count_masked_bases = true;
  for (const auto& interval : it->second) {
    if (count_masked_bases) {
      if (interval.start >= original_length) {
        count_masked_bases = false;
      } else {
        const Position interval_end =
            std::min(interval.end, original_length);
        summary.masked_bases += interval_end - interval.start;
      }
    }

    if (current_pos < interval.start) {
      ++summary.separator_count;
    }
    current_pos = std::max(current_pos, interval.end);
  }
  if (current_pos < original_length) {
    ++summary.separator_count;
  }

  summary.masked_sequence_length =
      original_length - summary.masked_bases;
  return summary;
}

void MaskManager::clearMaskIntervals(SequenceId seq_id) {
  mask_intervals_.erase(seq_id);
}

// ========================
// SequenceIndex 实现
// ========================

void SequenceIndex::addSequence(const SequenceInfo &info) {
  if (name_to_id_.count(info.name)) {
    throw SequenceException("Sequence name already exists: " + info.name);
  }

  sequences_.push_back(info);
  name_to_id_[info.name] = info.id;
}

const SequenceInfo *SequenceIndex::getSequenceInfo(const std::string &name) const {
  auto it = name_to_id_.find(name);
  if (it == name_to_id_.end()) {
    return nullptr;
  }
  return getSequenceInfo(it->second);
}

const SequenceInfo *SequenceIndex::getSequenceInfo(SequenceId id) const {
  for (const auto &seq : sequences_) {
    if (seq.id == id) {
      return &seq;
    }
  }
  return nullptr;
}

SequenceId SequenceIndex::getSequenceId(const std::string &name) const {
  auto it = name_to_id_.find(name);
  return (it != name_to_id_.end()) ? it->second : INVALID_ID;
}

std::vector<std::string> SequenceIndex::getSequenceNames() const {
  std::vector<std::string> names;
  names.reserve(sequences_.size());
  for (const auto &seq : sequences_) {
    names.push_back(seq.name);
  }
  return names;
}

void SequenceIndex::clear() {
  sequences_.clear();
  name_to_id_.clear();
}

Position SequenceIndex::getTotalGlobalLength() const {
  if (sequences_.empty()) {
    return 0;
  }
  const auto &last_seq = sequences_.back();
  return last_seq.global_start_pos + last_seq.length;
}

const SequenceInfo *SequenceIndex::getSequenceInfoByGlobalPosition(Position global_pos) const {
  if (sequences_.empty()) {
    return nullptr;
  }

  // 二分查找
  auto it = std::upper_bound(sequences_.begin(), sequences_.end(), global_pos,
                             [](Position pos, const SequenceInfo &info) {
                               return pos < info.global_start_pos;
                             });

  if (it == sequences_.begin()) {
    return nullptr;
  }

  --it;
  if (global_pos < it->global_start_pos + it->length) {
    return &(*it);
  }

  return nullptr;
}

// ========================
// MemoryMapper 实现
// ========================

MemoryMapper::MemoryMapper(const std::filesystem::path &file_path)
    : mapped_data_(nullptr), file_size_(0), fd_(-1) {

  fd_ = open(file_path.c_str(), O_RDONLY);
  if (fd_ == -1) {
    throw FileException("Cannot open file: " + file_path.string());
  }

  file_size_ = lseek(fd_, 0, SEEK_END);
  if (file_size_ == static_cast<size_t>(-1)) {
    close(fd_);
    throw FileException("Cannot determine file size: " + file_path.string());
  }

  mapped_data_ = mmap(nullptr, file_size_, PROT_READ, MAP_SHARED, fd_, 0);
  if (mapped_data_ == MAP_FAILED) {
    close(fd_);
    throw FileException("Cannot map file to memory: " + file_path.string());
  }
}

MemoryMapper::~MemoryMapper() { 
  cleanup(); 
}

MemoryMapper::MemoryMapper(MemoryMapper &&other) noexcept
    : mapped_data_(other.mapped_data_), file_size_(other.file_size_), fd_(other.fd_) {
  other.mapped_data_ = nullptr;
  other.file_size_ = 0;
  other.fd_ = -1;
}

MemoryMapper &MemoryMapper::operator=(MemoryMapper &&other) noexcept {
  if (this != &other) {
    cleanup();
    mapped_data_ = other.mapped_data_;
    file_size_ = other.file_size_;
    fd_ = other.fd_;
    other.mapped_data_ = nullptr;
    other.file_size_ = 0;
    other.fd_ = -1;
  }
  return *this;
}

std::span<const char> MemoryMapper::getData(Position offset, Length length) const {
  if (!isValid()) {
    throw FileException("Memory mapper is not valid");
  }

  if (offset + length > file_size_) {
    throw FileException("Request exceeds file size");
  }

  const char *data_ptr = static_cast<const char *>(mapped_data_) + offset;
  return std::span<const char>(data_ptr, length);
}

void MemoryMapper::cleanup() {
  if (mapped_data_ && mapped_data_ != MAP_FAILED) {
    munmap(mapped_data_, file_size_);
  }
  if (fd_ != -1) {
    close(fd_);
  }
  mapped_data_ = nullptr;
  file_size_ = 0;
  fd_ = -1;
}

// ========================
// SequenceManager 实现
// ========================

SequenceManager::SequenceManager(const std::filesystem::path &fasta_path)
    : fasta_path_(fasta_path), max_threads_(std::thread::hardware_concurrency()) {

  if (!utils::isValidFastaFile(fasta_path)) {
    throw FileException("Invalid FASTA file: " + fasta_path.string());
  }

  memory_mapper_ = std::make_unique<MemoryMapper>(fasta_path);
  buildIndex();
}

void SequenceManager::buildIndex() {
  has_ambiguous_bases = false;
  auto file_data = memory_mapper_->getData(0, memory_mapper_->getFileSize());
  sequence_index_.clear();

  std::string current_name;
  Position current_length = 0;
  Position sequence_start_offset = 0;
  SequenceId seq_id = 0;
  uint32_t line_width = 0;
  uint32_t line_bytes = 0;
  bool first_seq_line = true;
  Position global_offset = 0;

  Position i = 0;
  while (i < file_data.size()) {
    Position line_start = i;

    // 找到行尾
    while (i < file_data.size() && file_data[i] != '\n') {
      ++i;
    }
        Position line_end_pos = i;
    if (i < file_data.size()) {
      ++i; // 跳过换行符
    }

     Position current_line_length = line_end_pos - line_start;
    if (current_line_length > 0 && file_data[line_start + current_line_length - 1] == '\r') {
      current_line_length--; // 如果行尾是\r，则长度减一
    }

    if (current_line_length > 0 && file_data[line_start] == '>') {
      // 保存前一个序列
      if (!current_name.empty()) {
        SequenceInfo info(seq_id++, current_name, current_length, global_offset,global_offset,
                          sequence_start_offset, line_width, line_bytes);
        sequence_index_.addSequence(info);
        global_offset += current_length;
      }

      // 解析新序列名
      std::string header(file_data.data() + line_start + 1, current_line_length - 1);
      current_name = utils::cleanSequenceName(header);
      current_length = 0;
      sequence_start_offset = i;
      first_seq_line = true;

    } else if (!current_name.empty() && current_line_length > 0) {

      if (!has_ambiguous_bases) {
        const char* line_ptr = file_data.data() + line_start;
        for (Position j = 0; j < current_line_length; ++j) {
          const char base = line_ptr[j];
          if (base == 'N' || base == 'n') {
            has_ambiguous_bases = true;
            break; // 找到后立即退出内层循环
          }
        }
      }

      if (first_seq_line) {
        line_width = current_line_length;
        line_bytes = i - line_start;
        first_seq_line = false;
      }
      current_length += current_line_length;
    }
  }

  // 保存最后一个序列
  if (!current_name.empty()) {
    SequenceInfo info(seq_id, current_name, current_length, global_offset,global_offset,
                      sequence_start_offset, line_width, line_bytes);
    sequence_index_.addSequence(info);
  }
}

std::string SequenceManager::getSubSequence(const std::string &seq_name, Position start, Length length) const {
  SequenceId seq_id = getSequenceId(seq_name);
  if (seq_id == SequenceIndex::INVALID_ID) {
    throw SequenceException("Sequence not found: " + seq_name);
  }
  return getSubSequence(seq_id, start, length);
}

std::string SequenceManager::getSubSequence(SequenceId seq_id, Position start, Length length) const {
  std::shared_lock<std::shared_mutex> lock(mutex_);
  
  const auto *info = sequence_index_.getSequenceInfo(seq_id);
  if (!info) {
    throw SequenceException("Invalid sequence ID: " + std::to_string(seq_id));
  }

  if (start >= info->length) {
    throw SequenceException("Start position out of bounds");
  }

  length = std::min(length, info->length - start);
  return extractSequence(*info, start, length);
}

bool SequenceManager::tryGetContiguousSubSequence(
    SequenceId seq_id, Position start, Length length,
    std::string_view &view) const {
  std::shared_lock<std::shared_mutex> lock(mutex_);
  const auto *info = sequence_index_.getSequenceInfo(seq_id);
  if (!info) {
    throw SequenceException("Invalid sequence ID: " + std::to_string(seq_id));
  }
  if (start >= info->length) {
    throw SequenceException("Start position out of bounds");
  }

  length = std::min(length, info->length - start);
  if (length == 0) {
    view = {};
    return true;
  }
  if (info->line_width == 0 || info->line_bytes < info->line_width) {
    view = {};
    return false;
  }

  const Position offset_in_line = start % info->line_width;
  if (length > info->line_width - offset_in_line) {
    view = {};
    return false;
  }
  const Position line_number = start / info->line_width;
  const Position file_position = info->file_offset +
      line_number * info->line_bytes + offset_in_line;
  const auto data = memory_mapper_->getData(file_position, length);
  view = std::string_view(data.data(), data.size());
  return true;
}

void SequenceManager::getSubSequenceInto(
    SequenceId seq_id, Position start, Length length,
    std::string &output) const {
  std::shared_lock<std::shared_mutex> lock(mutex_);
  const auto *info = sequence_index_.getSequenceInfo(seq_id);
  if (!info) {
    throw SequenceException("Invalid sequence ID: " + std::to_string(seq_id));
  }
  if (start >= info->length) {
    throw SequenceException("Start position out of bounds");
  }
  length = std::min(length, info->length - start);
  extractSequenceInto(*info, start, length, output);
}

std::string SequenceManager::extractSequence(const SequenceInfo &info, Position start, Length length) const {
  std::string result;
  extractSequenceInto(info, start, length, result);
  return result;
}

void SequenceManager::extractSequenceInto(
    const SequenceInfo &info, Position start, Length length,
    std::string &result) const {
  result.clear();
  if (length == 0) return;
  result.reserve(length);

  Position seq_pos = start;
  Length remaining = length;

  while (remaining > 0) {
    // 计算在文件中的位置
    Position line_num = seq_pos / info.line_width;
    Position char_offset_in_line = seq_pos % info.line_width;
    Position file_pos = info.file_offset + line_num * info.line_bytes + char_offset_in_line;
    
    // 计算本次读取的长度
    Position chars_remaining_in_line = info.line_width - char_offset_in_line;
    Length chars_to_read = std::min(remaining, chars_remaining_in_line);

    // 读取数据
    auto data = memory_mapper_->getData(file_pos, chars_to_read);
    result.append(data.data(), chars_to_read);

    seq_pos += chars_to_read;
    remaining -= chars_to_read;
  }

}

std::string SequenceManager::getSubSequenceGlobal(Position global_start, Length length) const {
  std::shared_lock<std::shared_mutex> lock(mutex_);

  if (length == 0) return "";

  std::string result;
  result.reserve(length);

  Position current_pos = global_start;
  Length remaining = length;

  while (remaining > 0) {
    const auto *seq_info = sequence_index_.getSequenceInfoByGlobalPosition(current_pos);
    if (!seq_info) {
      throw SequenceException("Invalid global position: " + std::to_string(current_pos));
    }

    Position local_start = current_pos - seq_info->global_start_pos;
    Length read_length = std::min(remaining, seq_info->length - local_start);

    std::string segment = extractSequence(*seq_info, local_start, read_length);
    result.append(segment);

    current_pos += read_length;
    remaining -= read_length;
  }

  return result;
}

SequenceId SequenceManager::getSequenceId(const std::string &seq_name) const {
  std::shared_lock<std::shared_mutex> lock(mutex_);
  return sequence_index_.getSequenceId(seq_name);
}

std::vector<std::string> SequenceManager::getSequenceNames() const {
  std::shared_lock<std::shared_mutex> lock(mutex_);
  return sequence_index_.getSequenceNames();
}

std::string SequenceManager::getSequenceName(const uint32_t& seq_id) const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    const auto* info = sequence_index_.getSequenceInfo(seq_id);
    if (info) {
        return info->name;
    }
    else {
		throw SequenceException("Invalid sequence ID: " + std::to_string(seq_id));
    }
}

Length SequenceManager::getSequenceLength(const std::string &seq_name) const {
  std::shared_lock<std::shared_mutex> lock(mutex_);
  const auto *info = sequence_index_.getSequenceInfo(seq_name);
  return info ? info->length : 0;
}

Length SequenceManager::getSequenceLength(SequenceId seq_id) const {
  std::shared_lock<std::shared_mutex> lock(mutex_);
  const auto *info = sequence_index_.getSequenceInfo(seq_id);
  return info ? info->length : 0;
}

size_t SequenceManager::getSequenceCount() const {
  std::shared_lock<std::shared_mutex> lock(mutex_);
  return sequence_index_.getSequenceCount();
}

bool SequenceManager::isValidPosition(const std::string &seq_name, Position pos, Length len) const {
  std::shared_lock<std::shared_mutex> lock(mutex_);
  const auto *info = sequence_index_.getSequenceInfo(seq_name);
  if (!info) return false;
  return pos < info->length && pos + len <= info->length;
}

bool SequenceManager::isValidPosition(SequenceId seq_id, Position pos, Length len) const {
  std::shared_lock<std::shared_mutex> lock(mutex_);
  const auto *info = sequence_index_.getSequenceInfo(seq_id);
  if (!info) return false;
  return pos < info->length && pos + len <= info->length;
}

const SequenceInfo* SequenceManager::getSequenceInfoFromGlobalPosition(Position global_pos) const {
  std::shared_lock<std::shared_mutex> lock(mutex_);
  return sequence_index_.getSequenceInfoByGlobalPosition(global_pos);
}

std::string SequenceManager::getSequenceNameFromGlobalPosition(Position global_pos) const {
  const auto *info = getSequenceInfoFromGlobalPosition(global_pos);
  return info ? info->name : "";
}

std::pair<SequenceId, Position> SequenceManager::globalToLocal(Position global_pos) const {
  const auto *info = getSequenceInfoFromGlobalPosition(global_pos);
  if (info) {
    Position local_pos = global_pos - info->global_start_pos;
    return {info->id, local_pos};
  }
  return {SequenceIndex::INVALID_ID, 0};
}

Position SequenceManager::localToGlobal(const std::string& seq_name, Position local_pos) const {
  std::shared_lock<std::shared_mutex> lock(mutex_);
  const auto *info = sequence_index_.getSequenceInfo(seq_name);
  if (!info) {
    throw SequenceException("Sequence not found: " + seq_name);
  }
  if (local_pos >= info->length) {
    throw SequenceException("Local position out of bounds");
  }
  return info->global_start_pos + local_pos;
}

Position SequenceManager::localToGlobal(SequenceId seq_id, Position local_pos) const {
  std::shared_lock<std::shared_mutex> lock(mutex_);
  const auto *info = sequence_index_.getSequenceInfo(seq_id);
  if (!info) {
    throw SequenceException("Invalid sequence ID: " + std::to_string(seq_id));
  }
  if (local_pos >= info->length) {
    throw SequenceException("Local position out of bounds");
  }
  return info->global_start_pos + local_pos;
}

Length SequenceManager::getTotalLength() const {
  std::shared_lock<std::shared_mutex> lock(mutex_);
  return sequence_index_.getTotalGlobalLength();
}

std::string SequenceManager::concatAllSequences(char separator) const {
  auto seq_names = getSequenceNames();
  auto result = concatSequences(seq_names, separator);
  return result;
}



std::string SequenceManager::concatSequences(const std::vector<std::string> &seq_names, char separator) const {
  std::string result;
  
  // 预估大小
  size_t estimated_size = 0;
  for (const auto& seq_name : seq_names) {
    estimated_size += getSequenceLength(seq_name);
    estimated_size += 1;
  }
  estimated_size += 1;
  result.reserve(estimated_size);

  for (size_t i = 0; i < seq_names.size(); ++i) {
    try {
      Length length = getSequenceLength(seq_names[i]);
      std::string sequence = getSubSequence(seq_names[i], 0, length);
      result.append(sequence);
     
      result.push_back('\1');
      
    } catch (const std::exception&) {
      // 忽略错误
    }
  }
  result.push_back('\0');
  return result;
}


void SequenceManager::streamSequences(std::ostream &output, const std::vector<std::string> &seq_names,
                                     char separator, size_t buffer_size) const {
  for (size_t i = 0; i < seq_names.size(); ++i) {
    try {
      Length length = getSequenceLength(seq_names[i]);
      Position pos = 0;
      
      while (pos < length) {
        Length read_length = std::min(buffer_size, static_cast<size_t>(length - pos));
        std::string chunk = getSubSequence(seq_names[i], pos, read_length);
        output.write(chunk.c_str(), chunk.size());
        pos += read_length;
      }

      if (separator != '\0' && i < seq_names.size() - 1) {
        output.put(separator);
      }
    } catch (const std::exception&) {
      // 忽略错误
    }
  }
}

std::vector<SequenceManager::Result> SequenceManager::batchQuery(const std::vector<Query> &queries, 
                                                               size_t num_threads) const {
  if (num_threads == 0) {
    num_threads = max_threads_;
  }

  std::vector<Result> results;
  results.reserve(queries.size());

  if (num_threads <= 1 || queries.size() < num_threads) {
    // 单线程处理
    for (const auto& query : queries) {
      try {
        std::string sequence = getSubSequence(query.seq_name, query.start, query.length);
        results.emplace_back(Result{std::move(sequence), query.id, true, ""});
      } catch (const std::exception& e) {
        results.emplace_back(Result{"", query.id, false, e.what()});
      }
    }
  } else {
    // 多线程处理
    std::vector<std::future<Result>> futures;
    futures.reserve(queries.size());

    for (const auto& query : queries) {
      futures.emplace_back(std::async(std::launch::async, [this, &query]() -> Result {
        try {
          std::string sequence = getSubSequence(query.seq_name, query.start, query.length);
          return Result{std::move(sequence), query.id, true, ""};
        } catch (const std::exception& e) {
          return Result{"", query.id, false, e.what()};
        }
      }));
    }

    for (auto& future : futures) {
      results.push_back(future.get());
    }
  }

  return results;
}

void SequenceManager::setMaxThreads(size_t max_threads) {
  max_threads_ = std::max(size_t(1), max_threads);
}

void SequenceManager::clearMaskedRegions() {
  // 对于SequenceManager，这是一个空操作
  // 因为它不处理遮蔽区间
}

// ========================
// MaskedSequenceManager 实现
// ========================

MaskedSequenceManager::MaskedSequenceManager(std::unique_ptr<SequenceManager> seq_manager,
                                             const std::filesystem::path &mask_file)
      : original_manager_(std::move(seq_manager)), cache_valid_(false) {

  // 检查 original_manager_ 是否有效
  if (!original_manager_) {
    throw SeqProException("MaskedSequenceManager received a null SequenceManager.");
  }

  if (!mask_manager_.loadFromIntervalFile(mask_file)) {
    throw FileException("Cannot load mask file: " + mask_file.string());
  }

  auto seq_names = original_manager_->getSequenceNames();
  for (const auto &name : seq_names) {
    mask_manager_.getOrCreateSequenceId(name);
  }
}

// 创建一个空的遮蔽manager
MaskedSequenceManager::MaskedSequenceManager(std::unique_ptr<SequenceManager> seq_manager)
      : original_manager_(std::move(seq_manager)), cache_valid_(false) {

  if (!original_manager_) {
    throw SeqProException("MaskedSequenceManager received a null SequenceManager.");
  }

  // 初始化空的 name -> id 映射，确保与原始序列管理器一致
  auto seq_names = original_manager_->getSequenceNames();
  for (const auto &name : seq_names) {
    mask_manager_.getOrCreateSequenceId(name);
  }
}

std::string MaskedSequenceManager::getSubSequence(const std::string &seq_name, 
                                                Position start, Length length) const {
  SequenceId seq_id = getSequenceId(seq_name);
  if (seq_id == SequenceIndex::INVALID_ID) {
    throw SequenceException("Sequence not found: " + seq_name);
  }
  return getSubSequence(seq_id, start, length);
}

std::string MaskedSequenceManager::getSubSequence(SequenceId seq_id, 
                                                Position start, Length length) const {
  // 验证遮蔽坐标
  if (!isValidPosition(seq_id, start, length)) {
    throw SequenceException("Invalid masked position");
  }

  // 转换为原始坐标并获取有效区间
  Position original_start = toOriginalPosition(seq_id, start);
  // 有分隔符的最后一位碱基局部坐标
  Position original_end = toOriginalPosition(seq_id, start + length - 1);
  
  auto valid_ranges = mask_manager_.getValidRanges(seq_id, original_start, original_end);
  
  std::string result;
  result.reserve(length);
  
  for (const auto &range : valid_ranges) {
    Position range_length = range.second - range.first;
    std::string segment = original_manager_->getSubSequence(seq_id, range.first, range_length);
    result.append(segment);
  }

  return result;
}

std::string MaskedSequenceManager::getSubSequenceSeparated(const std::string &seq_name, Position start, Length length, char separator) const {
  SequenceId seq_id = getSequenceId(seq_name);
  if (seq_id == SequenceIndex::INVALID_ID) {
    throw SequenceException("Sequence not found: " + seq_name);
  }
  return getSubSequenceSeparated(seq_id, start, length, separator);
}

std::string MaskedSequenceManager::getSubSequenceSeparated(SequenceId seq_id, Position start, Length length, char separator) const {

  // 转换为原始坐标并获取有效区间
  Position original_start = toOriginalPositionSeparated(seq_id, start);
  Position original_end = toOriginalPositionSeparated(seq_id, start + length - 1);
  
  auto valid_ranges = mask_manager_.getValidRanges(seq_id, original_start, original_end + 1);
  
  std::string result;
  result.reserve(length);
  
  for (const auto &range : valid_ranges) {
    Position range_length = range.second - range.first;
    std::string segment = original_manager_->getSubSequence(seq_id, range.first, range_length);
    result.append(segment);
    result.push_back(separator);
  }

  return result;
}

std::string MaskedSequenceManager::getSubSequenceGlobal(Position global_start, Length length) const {
  ensureCacheValid();

  if (length == 0) return "";

  std::string result;
  result.reserve(length);

  Position current_pos = global_start;
  Length remaining = length;

  while (remaining > 0) {
    auto [seq_id, local_pos] = globalToLocal(current_pos);
    if (seq_id == SequenceIndex::INVALID_ID) {
      throw SequenceException("Invalid global masked position: " + std::to_string(current_pos));
    }

    Length seq_remaining = getSequenceLength(seq_id) - local_pos;
    Length read_length = std::min(remaining, seq_remaining);

    std::string segment = getSubSequence(seq_id, local_pos, read_length);
    result.append(segment);

    current_pos += read_length;
    remaining -= read_length;
  }

  return result;
}

Length MaskedSequenceManager::getSequenceLength(const std::string &seq_name) const {
  SequenceId seq_id = getSequenceId(seq_name);
  if (seq_id == SequenceIndex::INVALID_ID) {
    return 0;
  }
  return getSequenceLength(seq_id);
}

Length MaskedSequenceManager::getSequenceLength(SequenceId seq_id) const {
  Length original_length = original_manager_->getSequenceLength(seq_id);
  return mask_manager_.getMaskedSequenceLength(seq_id, original_length);
}

Length MaskedSequenceManager::getTotalLength() const {
  ensureCacheValid();
  return total_masked_length_;
}

Length MaskedSequenceManager::getTotalLengthWithSeparators() const {
  // 总长度 = 遮蔽后的序列长度 + 间隔符数量 + 序列数量 - 1 （因为拼接时还加了间隔符）
  return getTotalLength() + getTotalSeparatorCount() + getSequenceCount() - 1;
}

Length MaskedSequenceManager::getSequenceLengthWithSeparators(const std::string &seq_name) const {
  SequenceId seq_id = getSequenceId(seq_name);
  if (seq_id == SequenceIndex::INVALID_ID) {
    return 0;
  }
  return getSequenceLengthWithSeparators(seq_id);
}

Length MaskedSequenceManager::getSequenceLengthWithSeparators(SequenceId seq_id) const {
  // 序列长度 = 遮蔽后的序列长度 + 该序列的间隔符数量
  return getSequenceLength(seq_id) + getSeparatorCount(seq_id);
}

bool MaskedSequenceManager::isValidPosition(const std::string &seq_name, Position pos, Length len) const {
  SequenceId seq_id = getSequenceId(seq_name);
  if (seq_id == SequenceIndex::INVALID_ID) {
    return false;
  }
  return isValidPosition(seq_id, pos, len);
}

bool MaskedSequenceManager::isValidPosition(SequenceId seq_id, Position pos, Length len) const {
  Length masked_length = getSequenceLength(seq_id);
  return pos < masked_length && pos + len <= masked_length;
}

Position MaskedSequenceManager::toOriginalPosition(const std::string &seq_name, Position masked_pos) const {
  SequenceId seq_id = getSequenceId(seq_name);
  if (seq_id == SequenceIndex::INVALID_ID) {
    throw SequenceException("Sequence not found: " + seq_name);
  }
  return toOriginalPosition(seq_id, masked_pos);
}

Position MaskedSequenceManager::toOriginalPosition(SequenceId seq_id, Position masked_pos) const {
  return mask_manager_.mapToOriginalPosition(seq_id, masked_pos);
}

Position MaskedSequenceManager::toOriginalPositionSeparated(const std::string &seq_name, Position masked_pos) const {
  SequenceId seq_id = getSequenceId(seq_name);
  if (seq_id == SequenceIndex::INVALID_ID) {
    throw SequenceException("Sequence not found: " + seq_name);
  }
  return toOriginalPositionSeparated(seq_id, masked_pos);
}
// 输入拼接序列的局部坐标（包括\01），返回原始的局部坐标(不包含01)
Position MaskedSequenceManager::toOriginalPositionSeparated(SequenceId seq_id, Position masked_pos) const {
  return mask_manager_.mapToOriginalPositionSeparated(seq_id, masked_pos);
}

Position MaskedSequenceManager::toMaskedPosition(const std::string &seq_name, Position original_pos) const {
  SequenceId seq_id = getSequenceId(seq_name);
  if (seq_id == SequenceIndex::INVALID_ID) {
    throw SequenceException("Sequence not found: " + seq_name);
  }
  return toMaskedPosition(seq_id, original_pos);
}

Position MaskedSequenceManager::toMaskedPosition(SequenceId seq_id, Position original_pos) const {
  return mask_manager_.mapToMaskedPosition(seq_id, original_pos);
}

std::vector<Position> MaskedSequenceManager::toOriginalPositions(const std::string &seq_name, 
                                                               const std::vector<Position> &masked_positions) const {
  std::vector<Position> results;
  results.reserve(masked_positions.size());
  for (Position pos : masked_positions) {
    results.push_back(toOriginalPosition(seq_name, pos));
  }
  return results;
}

std::vector<Position> MaskedSequenceManager::toMaskedPositions(const std::string &seq_name, 
                                                             const std::vector<Position> &original_positions) const {
  std::vector<Position> results;
  results.reserve(original_positions.size());
  for (Position pos : original_positions) {
    results.push_back(toMaskedPosition(seq_name, pos));
  }
  return results;
}

std::pair<SequenceId, Position> MaskedSequenceManager::globalToLocal(Position global_masked_pos) const {
  ensureCacheValid();
  
  const auto *seq_info = getSequenceInfoByMaskedGlobalPosition(global_masked_pos);
  if (seq_info) {
    auto it = global_offset_cache_.find(seq_info->id);
    if (it != global_offset_cache_.end()) {
      Position local_pos = global_masked_pos - it->second;
      return {seq_info->id, local_pos};
    }
  }
  
  return {SequenceIndex::INVALID_ID, 0};
}

Position MaskedSequenceManager::localToGlobal(const std::string& seq_name, Position local_masked_pos) const {
  SequenceId seq_id = getSequenceId(seq_name);
  if (seq_id == SequenceIndex::INVALID_ID) {
    throw SequenceException("Sequence not found: " + seq_name);
  }
  return localToGlobal(seq_id, local_masked_pos);
}

Position MaskedSequenceManager::localToGlobal(SequenceId seq_id, Position local_masked_pos) const {
  ensureCacheValid();

  if (!isValidPosition(seq_id, local_masked_pos)) {
    throw SequenceException("Local masked position out of bounds");
  }

  auto it = global_offset_cache_.find(seq_id);
  if (it != global_offset_cache_.end()) {
    return it->second + local_masked_pos;
  }

  throw SequenceException("Invalid sequence ID: " + std::to_string(seq_id));
}

std::pair<SequenceId, Position> MaskedSequenceManager::globalToLocalSeparated(Position global_pos_with_separators) const {
  ensureCacheValid();
  // 构建包含间隔符的全局坐标映射
  Position current_global_pos = 0;
  auto seq_names = getSequenceNames();

  for (const auto& seq_name : seq_names) {
    SequenceId seq_id = getSequenceId(seq_name);
    Length seq_length_with_separators = getSequenceLengthWithSeparators(seq_id);

    if (global_pos_with_separators >= current_global_pos &&
        global_pos_with_separators < current_global_pos + seq_length_with_separators) {
      // 位置在当前序列范围内
      Position local_pos_with_separators = global_pos_with_separators - current_global_pos;

      // 需要将包含间隔符的本地位置转换为不含间隔符的遮蔽位置
      //Position local_masked_pos = convertSeparatedToMaskedPosition(seq_id, local_pos_with_separators);
      Position local_original_pos = toOriginalPositionSeparated(seq_id, local_pos_with_separators);
      // local_original_pos -= seq_id;
      return {seq_id , local_original_pos};
    }

    current_global_pos += seq_length_with_separators;

    // 添加染色体间的间隔符（除了最后一个序列）
    if (seq_name != seq_names.back()) {
      if (global_pos_with_separators == current_global_pos) {
        // 位置正好在染色体间隔符上，返回下一个序列的开始位置
        auto next_seq_it = std::find(seq_names.begin(), seq_names.end(), seq_name);
        if (next_seq_it != seq_names.end() && ++next_seq_it != seq_names.end()) {
          return {getSequenceId(*next_seq_it), 0};
        }
      }
      current_global_pos += 1; // 染色体间隔符
    }
  }

  return {SequenceIndex::INVALID_ID, INVALID_POSITION};
}

Position MaskedSequenceManager::localToGlobalSeparated(const std::string& seq_name, Position local_masked_pos) const {
  SequenceId seq_id = getSequenceId(seq_name);
  if (seq_id == SequenceIndex::INVALID_ID) {
    throw SequenceException("Sequence not found: " + seq_name);
  }
  return localToGlobalSeparated(seq_id, local_masked_pos);
}

Position MaskedSequenceManager::localToGlobalSeparated(SequenceId seq_id, Position local_masked_pos) const {
  ensureCacheValid();

  if (!isValidPosition(seq_id, local_masked_pos)) {
    throw SequenceException("Local masked position out of bounds");
  }

  // 计算目标序列之前所有序列的累积长度（包含间隔符）
  Position global_offset = 0;
  auto seq_names = getSequenceNames();

  for (const auto& seq_name : seq_names) {
    SequenceId current_seq_id = getSequenceId(seq_name);

    if (current_seq_id == seq_id) {
      // 找到目标序列，将本地遮蔽位置转换为包含间隔符的位置
      Position local_pos_with_separators = convertMaskedToSeparatedPosition(seq_id, local_masked_pos);
      return global_offset + local_pos_with_separators;
    }

    // 累加当前序列的长度（包含间隔符）
    global_offset += getSequenceLengthWithSeparators(current_seq_id);

    // 添加染色体间的间隔符（除了最后一个序列）
    if (seq_name != seq_names.back()) {
      global_offset += 1;
    }
  }

  throw SequenceException("Invalid sequence ID: " + std::to_string(seq_id));
}

bool MaskedSequenceManager::isMaskPosition(const std::string &seq_name, Position original_pos) const {
  SequenceId seq_id = getSequenceId(seq_name);
  if (seq_id == SequenceIndex::INVALID_ID) {
    return false;
  }
  return isMaskPosition(seq_id, original_pos);
}

bool MaskedSequenceManager::isMaskPosition(SequenceId seq_id, Position original_pos) const {
  return mask_manager_.isMaskPosition(seq_id, original_pos);
}

const std::vector<MaskInterval>& MaskedSequenceManager::getMaskIntervals(const std::string &seq_name) const {
  SequenceId seq_id = getSequenceId(seq_name);
  return mask_manager_.getMaskIntervals(seq_id);
}

const std::vector<MaskInterval>& MaskedSequenceManager::getMaskIntervals(SequenceId seq_id) const {
  return mask_manager_.getMaskIntervals(seq_id);
}

Length MaskedSequenceManager::getMaskedBases(const std::string &seq_name) const {
  SequenceId seq_id = getSequenceId(seq_name);
  if (seq_id == SequenceIndex::INVALID_ID) {
    return 0;
  }
  return getMaskedBases(seq_id);
}

Length MaskedSequenceManager::getMaskedBases(SequenceId seq_id) const {
  return mask_manager_.getTotalMaskedBases(seq_id);
}

Length MaskedSequenceManager::getTotalMaskedBases() const {
  Length total = 0;
  auto seq_names = getSequenceNames();
  for (const auto& seq_name : seq_names) {
    total += getMaskedBases(seq_name);
  }
  return total;
}

Length MaskedSequenceManager::getSeparatorCount(const std::string &seq_name) const {
  SequenceId seq_id = getSequenceId(seq_name);
  if (seq_id == SequenceIndex::INVALID_ID) {
    return 0;
  }
  return getSeparatorCount(seq_id);
}

Length MaskedSequenceManager::getSeparatorCount(SequenceId seq_id) const {
  return mask_manager_.getSeparatorCount(seq_id, original_manager_->getSequenceLength(seq_id));
}

Length MaskedSequenceManager::getTotalSeparatorCount() const {
  Length total = 0;
  auto seq_names = getSequenceNames();
  for (const auto& seq_name : seq_names) {
    total += getSeparatorCount(seq_name);
  }
  // 每个染色体之间有1个间隔符
  total += getSequenceCount() - 1;
  return total;
}

std::string MaskedSequenceManager::concatAllSequences(char separator) const {
  auto seq_names = getSequenceNames();
  auto result = concatSequences(seq_names, separator);
  if(separator != '\0'){
    result.push_back('\0');
  }
  return result;
}

std::string MaskedSequenceManager::concatAllSequencesSeparated(char separator) const {
  auto seq_names = getSequenceNames();
  auto result = concatSequencesSeparated(seq_names, separator);
  return result;
}

std::string MaskedSequenceManager::concatSequences(const std::vector<std::string> &seq_names, char separator) const {
  std::string result;
  
  // 预估大小
  size_t estimated_size = 0;
  for (const auto& seq_name : seq_names) {
    estimated_size += getSequenceLength(seq_name);
    if (separator != '\1') {
      estimated_size += 1;
    }
  }
  result.reserve(estimated_size);

  for (size_t i = 0; i < seq_names.size(); ++i) {
    try {
      Length length = getSequenceLength(seq_names[i]);
      std::string sequence = getSubSequence (seq_names[i], 0, length);
      result.append(sequence);
      
      if (separator != '\1' && i < seq_names.size() - 1) {
        result.push_back(separator);
      }
    } catch (const std::exception&) {
      // 忽略错误
    }
  }

  return result;
}

  std::string MaskedSequenceManager::concatSequencesSeparated(const std::vector<std::string> &seq_names, char separator) const {
  std::string result;
  // TODO 性能优化
  // 预估大小
  //size_t estimated_size = 0;
  //for (const auto& seq_name : seq_names) {
  //  estimated_size += getSequenceLengthWithSeparators(seq_name);
  //  estimated_size += 1;
  //}
  //estimated_size += 1;
  //result.reserve(estimated_size);

  for (size_t i = 0; i < seq_names.size(); ++i)
  {
    Length length = getSequenceLengthWithSeparators(seq_names[i]);
    std::string sequence;
    if (length > 0)
    {
      sequence = getSubSequenceSeparated(seq_names[i], 0, length - 1, separator);
    }else
    {
      sequence = "";
    }

    result.append(sequence);

    result.push_back('\1');

  }
  result.push_back('\0');

  return result;
}

void MaskedSequenceManager::streamSequences(std::ostream &output, const std::vector<std::string> &seq_names,
                                          char separator, size_t buffer_size) const {
  for (size_t i = 0; i < seq_names.size(); ++i) {
    try {
      Length length = getSequenceLength(seq_names[i]);
      Position pos = 0;
      
      while (pos < length) {
        Length read_length = std::min(buffer_size, static_cast<size_t>(length - pos));
        std::string chunk = getSubSequence(seq_names[i], pos, read_length);
        output.write(chunk.c_str(), chunk.size());
        pos += read_length;
      }

      if (separator != '\0' && i < seq_names.size() - 1) {
        output.put(separator);
      }
    } catch (const std::exception&) {
      // 忽略错误
    }
  }
}

std::vector<MaskedSequenceManager::Result> MaskedSequenceManager::batchQuery(const std::vector<Query> &queries, 
                                                                           size_t num_threads) const {
  if (num_threads == 0) {
    num_threads = original_manager_->getMaxThreads();
  }

  std::vector<Result> results;
  results.reserve(queries.size());

  if (num_threads <= 1 || queries.size() < num_threads) {
    // 单线程处理
    for (const auto& query : queries) {
      try {
        std::string sequence = getSubSequence(query.seq_name, query.start, query.length);
        results.emplace_back(Result{std::move(sequence), query.id, true, ""});
      } catch (const std::exception& e) {
        results.emplace_back(Result{"", query.id, false, e.what()});
      }
    }
  } else {
    // 多线程处理
    std::vector<std::future<Result>> futures;
    futures.reserve(queries.size());

    for (const auto& query : queries) {
      futures.emplace_back(std::async(std::launch::async, [this, &query]() -> Result {
        try {
          std::string sequence = getSubSequence(query.seq_name, query.start, query.length);
          return Result{std::move(sequence), query.id, true, ""};
        } catch (const std::exception& e) {
          return Result{"", query.id, false, e.what()};
        }
      }));
    }

    for (auto& future : futures) {
      results.push_back(future.get());
    }
  }

  return results;
}

std::vector<std::string> MaskedSequenceManager::getSequenceNames() const {
  return original_manager_->getSequenceNames();
}

SequenceId MaskedSequenceManager::getSequenceId(const std::string &seq_name) const {
  return original_manager_->getSequenceId(seq_name);
}

size_t MaskedSequenceManager::getSequenceCount() const {
  return original_manager_->getSequenceCount();
}

void MaskedSequenceManager::buildGlobalOffsetCache() const {
  global_offset_cache_.clear();
  total_masked_length_ = 0;

  auto seq_names = getSequenceNames();
  for (const auto& seq_name : seq_names) {
    SequenceId seq_id = getSequenceId(seq_name);
    global_offset_cache_[seq_id] = total_masked_length_;
    total_masked_length_ += getSequenceLength(seq_id);
  }

  cache_valid_ = true;
}

void MaskedSequenceManager::addMaskInterval(const std::string &seq_name, const MaskInterval &interval) {
  // 转换为内部ID
  SequenceId seq_id = original_manager_->getSequenceId(seq_name);
  if (seq_id == SequenceIndex::INVALID_ID) {
    throw SequenceException("Sequence not found: " + seq_name);
  }
  
  // 直接调用内部 mask_manager_ 的方法
  mask_manager_.addMaskInterval(seq_id, interval);
  
  // 添加区间后，全局坐标缓存会失效
  cache_valid_ = false;
}

// === 批量遮蔽区间管理实现 ===

void MaskedSequenceManager::addMaskIntervals(const std::string &seq_name, 
                                           const std::vector<MaskInterval> &intervals) {
  SequenceId seq_id = getSequenceId(seq_name);
  if (seq_id == SequenceIndex::INVALID_ID) {
    throw SequenceException("Sequence not found: " + seq_name);
  }
  addMaskIntervals(seq_id, intervals);
}

void MaskedSequenceManager::addMaskIntervals(SequenceId seq_id, 
                                           const std::vector<MaskInterval> &intervals) {
  if (intervals.empty()) return;
  
  // 批量转换遮蔽坐标为原始坐标
  std::vector<MaskInterval> original_intervals;
  original_intervals.reserve(intervals.size());
  
  for (const auto& interval : intervals) {
    // MaskInterval original_interval = convertMaskedToOriginalInterval(seq_id, interval);
    //if (original_interval.start > 1000000) {
    //}
    original_intervals.push_back(interval);
  }
  
  // 使用 MaskManager 的公共接口添加
  mask_manager_.addMaskIntervals(seq_id, original_intervals);
  
  // 标记为未定案
  unfinalized_sequences_.insert(seq_id);
  
  // 缓存失效
  cache_valid_ = false;
}

MaskBatchMergeStats MaskedSequenceManager::applyFinalizedMaskDeltas(
    std::vector<MaskIntervalDelta> deltas) {
  MaskBatchMergeStats total;
  bool changed = false;

  for (auto& delta : deltas) {
    if (delta.intervals.empty()) {
      continue;
    }
    if (delta.sequence_id == SequenceIndex::INVALID_ID ||
        delta.sequence_id >= getSequenceCount()) {
      throw SequenceException(
          "Invalid sequence ID: " + std::to_string(delta.sequence_id));
    }

    changed = true;
    cache_valid_ = false;
    if (unfinalized_sequences_.count(delta.sequence_id)) {
      mask_manager_.finalizeMaskIntervals(delta.sequence_id);
      unfinalized_sequences_.erase(delta.sequence_id);
    }

    const MaskBatchMergeStats sequence_stats =
        mask_manager_.mergeFinalizedMaskIntervals(
            delta.sequence_id, std::move(delta.intervals));
    total.touched_sequences += sequence_stats.touched_sequences;
    total.incoming_intervals += sequence_stats.incoming_intervals;
    total.normalized_delta_intervals +=
        sequence_stats.normalized_delta_intervals;
    total.previous_intervals += sequence_stats.previous_intervals;
    total.final_intervals += sequence_stats.final_intervals;
    total.sort_seconds += sequence_stats.sort_seconds;
    total.merge_seconds += sequence_stats.merge_seconds;
  }

  if (!changed) {
    return total;
  }

  const auto metadata_started = std::chrono::steady_clock::now();
  refreshMaskMetadataAndCache();
  total.metadata_seconds =
      std::chrono::duration<double>(
          std::chrono::steady_clock::now() - metadata_started)
          .count();
  return total;
}

bool MaskedSequenceManager::loadMaskIntervalsFromFile(const std::filesystem::path &file_path, 
                                                    bool append) {
  if (!append) {
    mask_manager_.clear();
    unfinalized_sequences_.clear();
  }
  
  bool success = mask_manager_.loadFromIntervalFile(file_path, append);
  if (success) {
    // 标记所有序列为已定案（从文件加载的区间已经是原始坐标）
    unfinalized_sequences_.clear();
    cache_valid_ = false;
  }
  return success;
}

void MaskedSequenceManager::finalizeMaskIntervals() {
  // 复制集合以避免在迭代过程中修改
  auto sequences_to_finalize = unfinalized_sequences_;
  for (SequenceId seq_id : sequences_to_finalize) {
    finalizeMaskIntervals(seq_id);
  }
  refreshMaskMetadataAndCache();
}

void MaskedSequenceManager::finalizeMaskIntervals(const std::string &seq_name) {
  SequenceId seq_id = getSequenceId(seq_name);
  if (seq_id != SequenceIndex::INVALID_ID) {
    finalizeMaskIntervals(seq_id);
  }
}

void MaskedSequenceManager::finalizeMaskIntervals(SequenceId seq_id) {
  if (unfinalized_sequences_.count(seq_id)) {
    mask_manager_.finalizeMaskIntervals(seq_id);
    unfinalized_sequences_.erase(seq_id);
    cache_valid_ = false;
  }
}

void MaskedSequenceManager::refreshMaskMetadataAndCache() {
  const size_t sequence_count = getSequenceCount();
  global_offset_cache_.clear();
  global_offset_cache_.reserve(sequence_count);
  total_masked_length_ = 0;
  Position current_masked_pos_with_separators = 0;

  for (size_t index = 0; index < sequence_count; ++index) {
    const SequenceId seq_id = static_cast<SequenceId>(index);
    auto* seq_info = const_cast<SequenceInfo*>(
        original_manager_->getIndex().getSequenceInfo(seq_id));
    if (!seq_info) {
      throw SequenceException(
          "Invalid sequence ID while refreshing mask metadata: " +
          std::to_string(seq_id));
    }

    const FinalizedMaskSummary summary =
        mask_manager_.summarizeFinalizedIntervals(
            seq_id, seq_info->length);
    global_offset_cache_[seq_id] = total_masked_length_;
    total_masked_length_ += summary.masked_sequence_length;

    seq_info->masked_global_start_pos =
        current_masked_pos_with_separators;
    seq_info->masked_length =
        summary.masked_sequence_length + summary.separator_count;
    current_masked_pos_with_separators += seq_info->masked_length;
    if (index + 1 < sequence_count) {
      ++current_masked_pos_with_separators;
    }
  }

  cache_valid_ = true;
}

void MaskedSequenceManager::clearMaskIntervals(const std::string &seq_name) {
  SequenceId seq_id = getSequenceId(seq_name);
  if (seq_id != SequenceIndex::INVALID_ID) {
    clearMaskIntervals(seq_id);
  }
}

void MaskedSequenceManager::clearMaskIntervals(SequenceId seq_id) {
  mask_manager_.clearMaskIntervals(seq_id);
  unfinalized_sequences_.erase(seq_id);
  cache_valid_ = false;
}

void MaskedSequenceManager::clearMaskedRegions() {
  // 清除所有序列的遮蔽区间
  std::vector<std::string> seq_names = getSequenceNames();
  for (const auto& seq_name : seq_names) {
    clearMaskIntervals(seq_name);
  }
}

bool MaskedSequenceManager::hasUnfinalizedIntervals() const {
  return !unfinalized_sequences_.empty();
}

size_t MaskedSequenceManager::getUnfinalizedSequenceCount() const {
  return unfinalized_sequences_.size();
}

MaskInterval MaskedSequenceManager::convertMaskedToOriginalInterval(SequenceId seq_id, 
                                                                  const MaskInterval &masked_interval) const {
  Position orig_start = toOriginalPosition(seq_id, masked_interval.start);
  Position orig_end = toOriginalPosition(seq_id, masked_interval.end - 1) + 1;
  return MaskInterval(orig_start, orig_end);
}

void MaskedSequenceManager::ensureFinalized(SequenceId seq_id) const {
  if (unfinalized_sequences_.count(seq_id)) {
    // 由于这是const函数，需要cast掉const
    const_cast<MaskedSequenceManager*>(this)->finalizeMaskIntervals(seq_id);
  }
}

void MaskedSequenceManager::ensureCacheValid() const {
  if (!cache_valid_) {
    buildGlobalOffsetCache();
  }
}

const SequenceInfo* MaskedSequenceManager::getSequenceInfoByMaskedGlobalPosition(Position global_masked_pos) const {
  ensureCacheValid();

  // 查找包含该全局遮蔽位置的序列
  const SequenceInfo* found_info = nullptr;

  for (const auto& [seq_id, global_offset] : global_offset_cache_) {
    Length masked_length = getSequenceLength(seq_id);
    if (global_masked_pos >= global_offset &&
        global_masked_pos < global_offset + masked_length) {
      found_info = original_manager_->getIndex().getSequenceInfo(seq_id);
      break;
    }
  }

  return found_info;
}

Position MaskedSequenceManager::convertSeparatedToMaskedPosition(SequenceId seq_id, Position separated_pos) const {
  // 获取该序列的遮蔽区间
  const auto& mask_intervals = getMaskIntervals(seq_id);

  if (mask_intervals.empty()) {
    // 没有遮蔽区间，直接返回
    return separated_pos;
  }

  Position separator_count_before = 0;

  // 遍历遮蔽区间，计算在当前位置之前有多少个间隔符需要减去
  for (const auto& interval : mask_intervals) {
    // 将原始坐标的遮蔽区间转换为遮蔽坐标系统中的位置
    Position masked_interval_start = toMaskedPosition(seq_id, interval.start);

    if (masked_interval_start == INVALID_POSITION) {
      // 如果区间开始位置无效，跳过
      continue;
    }

    // 计算在包含间隔符的坐标系统中，这个遮蔽区间的间隔符位置
    Position interval_separator_pos = masked_interval_start + separator_count_before;

    if (separated_pos > interval_separator_pos) {
      // 当前位置在这个遮蔽区间的间隔符之后，需要减去1个间隔符
      separator_count_before += 1;
    } else {
      // 当前位置在这个遮蔽区间之前，不需要继续计算
      break;
    }
  }

  // 减去间隔符的数量得到遮蔽位置
  return separated_pos - separator_count_before;
}

Position MaskedSequenceManager::convertMaskedToSeparatedPosition(SequenceId seq_id, Position masked_pos) const {
  // 获取该序列的遮蔽区间
  const auto& mask_intervals = getMaskIntervals(seq_id);

  if (mask_intervals.empty()) {
    // 没有遮蔽区间，直接返回
    return masked_pos;
  }

  Position separator_count_before = 0;

  // 遍历遮蔽区间，计算在当前位置之前有多少个间隔符需要加上
  for (const auto& interval : mask_intervals) {
    // 将原始坐标的遮蔽区间转换为遮蔽坐标系统中的位置
    Position masked_interval_start = toMaskedPosition(seq_id, interval.start);

    if (masked_interval_start == INVALID_POSITION) {
      // 如果区间开始位置无效，跳过
      continue;
    }

    if (masked_pos >= masked_interval_start) {
      // 当前位置在这个遮蔽区间之后或之上，需要加上1个间隔符
      separator_count_before += 1;
    } else {
      // 当前位置在这个遮蔽区间之前，不需要继续计算
      break;
    }
  }

  // 加上间隔符的数量得到包含间隔符的位置
  return masked_pos + separator_count_before;
}

// ========================
// 工具函数实现
// ========================

namespace utils {
std::filesystem::path cleanFastaFile(const std::filesystem::path &input_fasta,
                                              const std::filesystem::path &output_fasta,
                                              uint64_t line_width) {
  gzFile fp; // 文件指针
  kseq_t *seq; // 序列结构体
  int l; // 用于存储读取序列的长度
  std::ofstream output(output_fasta);

  fp = gzopen(input_fasta.c_str(), "r");
  if (fp == NULL) {
    throw FileException("Cannot open input or output file");
  }

  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    std::string line;
    std::string current_sequence;
    std::string current_header;
    line = seq->seq.s;
    current_header = seq->name.s;
    for (char c : line)
    {
      char upper_c = std::toupper(c);

      if (upper_c == 'A' || upper_c == 'T' || upper_c == 'G' || upper_c == 'C') {
          current_sequence += upper_c;
        } else {
          current_sequence += 'N'; // 非法字符转为N
        }
    }
    output << '>' << current_header << "\n";
    output << current_sequence << "\n";

  }
  // 释放 kseq 结构体
  kseq_destroy(seq);

  // 关闭文件
  gzclose(fp);

  return output_fasta;
}


bool isValidFastaFile(const std::filesystem::path &file_path) {
  if (!std::filesystem::exists(file_path)) {
    return false;
  }

  std::ifstream file(file_path);
  if (!file) {
    return false;
  }

  std::string first_line;
  std::getline(file, first_line);
  return !first_line.empty() && first_line[0] == '>';
}

bool isValidMaskFile(const std::filesystem::path &mask_file) {
  if (!std::filesystem::exists(mask_file)) {
    return false;
  }

  std::ifstream file(mask_file);
  if (!file) {
    return false;
  }

  std::string first_line;
  std::getline(file, first_line);
  return !first_line.empty() && (first_line[0] == '>' || std::isdigit(first_line[0]));
}

std::string getReadableFileSize(const std::filesystem::path &file_path) {
  if (!std::filesystem::exists(file_path)) {
    return "File not found";
  }

  auto size = std::filesystem::file_size(file_path);
  const char *units[] = {"B", "KB", "MB", "GB", "TB"};
  int unit = 0;
  double readable_size = static_cast<double>(size);

  while (readable_size >= 1024.0 && unit < 4) {
    readable_size /= 1024.0;
    unit++;
  }

  std::ostringstream oss;
  oss.precision(2);
  oss << std::fixed << readable_size << " " << units[unit];
  return oss.str();
}

std::string cleanSequenceName(const std::string &name) {
  std::string clean_name = name;
  
  // 去除首尾空白
  clean_name.erase(0, clean_name.find_first_not_of(" \t\n\r\f\v"));
  clean_name.erase(clean_name.find_last_not_of(" \t\n\r\f\v") + 1);

  // 截取到第一个空格或制表符
  size_t space_pos = clean_name.find_first_of(" \t");
  if (space_pos != std::string::npos) {
    clean_name = clean_name.substr(0, space_pos);
  }

  return clean_name;
}

bool isValidDNASequence(const std::string &sequence) {
  return std::all_of(sequence.begin(), sequence.end(), [](char c) {
    char upper_c = std::toupper(c);
    return upper_c == 'A' || upper_c == 'T' || upper_c == 'G' ||
           upper_c == 'C' || upper_c == 'N';
  });
}

} // namespace utils

} // namespace SeqPro
