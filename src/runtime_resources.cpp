#include "runtime_resources.h"

#include "spdlog/spdlog.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <system_error>

namespace RaMAxResources {
namespace {

constexpr uint64_t kKiB = 1024ULL;
constexpr uint64_t kMiB = 1024ULL * kKiB;
constexpr uint64_t kGiB = 1024ULL * kMiB;
constexpr uint64_t kTiB = 1024ULL * kGiB;

uint64_t saturatedAdd(uint64_t left, uint64_t right) noexcept {
  return right > std::numeric_limits<uint64_t>::max() - left
      ? std::numeric_limits<uint64_t>::max()
      : left + right;
}

uint64_t multiplyChecked(uint64_t value, uint64_t scale) {
  if (value > std::numeric_limits<uint64_t>::max() / scale) {
    throw std::invalid_argument("Memory limit is too large");
  }
  return value * scale;
}

uint64_t effectiveHostLimit(
    const RaMAxMemory::ProcessMemorySnapshot& memory) noexcept {
  if (memory.physical_memory_bytes == 0) {
    return memory.cgroup_limit_bytes;
  }
  if (memory.cgroup_limit_bytes == 0) {
    return memory.physical_memory_bytes;
  }
  return std::min(memory.physical_memory_bytes,
                  memory.cgroup_limit_bytes);
}

const char* pressureName(MemoryPressure pressure) noexcept {
  switch (pressure) {
    case MemoryPressure::NORMAL: return "normal";
    case MemoryPressure::HOLD_BATCH: return "hold-batch";
    case MemoryPressure::SPILL: return "spill";
    case MemoryPressure::RECLAIM: return "reclaim";
    case MemoryPressure::CRITICAL: return "critical";
  }
  return "unknown";
}

}  // namespace

uint64_t parseMemorySize(std::string_view text) {
  if (text.empty()) {
    throw std::invalid_argument("--memory-limit must not be empty");
  }
  if (text == "auto" || text == "AUTO" || text == "Auto") return 0;

  size_t digit_end = 0;
  while (digit_end < text.size() &&
         std::isdigit(static_cast<unsigned char>(text[digit_end]))) {
    ++digit_end;
  }
  if (digit_end == 0) {
    throw std::invalid_argument(
        "--memory-limit must be auto or an integer followed by "
        "KiB, MiB, GiB, or TiB");
  }
  const std::string digits(text.substr(0, digit_end));
  size_t consumed = 0;
  const uint64_t value = std::stoull(digits, &consumed);
  if (consumed != digits.size() || value == 0) {
    throw std::invalid_argument("--memory-limit must be positive");
  }

  std::string suffix(text.substr(digit_end));
  for (char& character : suffix) {
    character = static_cast<char>(
        std::toupper(static_cast<unsigned char>(character)));
  }
  uint64_t scale = 0;
  if (suffix == "KIB") {
    scale = kKiB;
  } else if (suffix == "MIB") {
    scale = kMiB;
  } else if (suffix == "GIB") {
    scale = kGiB;
  } else if (suffix == "TIB") {
    scale = kTiB;
  } else {
    throw std::invalid_argument(
        "--memory-limit suffix must be KiB, MiB, GiB, or TiB");
  }
  return multiplyChecked(value, scale);
}

std::string formatBytes(uint64_t bytes) {
  static constexpr std::array<const char*, 5> units{
      "B", "KiB", "MiB", "GiB", "TiB"};
  double value = static_cast<double>(bytes);
  size_t unit = 0;
  while (value >= 1024.0 && unit + 1 < units.size()) {
    value /= 1024.0;
    ++unit;
  }
  char buffer[64];
  std::snprintf(buffer, sizeof(buffer), "%.3f %s", value, units[unit]);
  return buffer;
}

RuntimeResourceManager& RuntimeResourceManager::instance() {
  static RuntimeResourceManager manager;
  return manager;
}

void RuntimeResourceManager::configure(const RuntimeResourceConfig& config) {
  if (config.work_dir.empty()) {
    throw std::invalid_argument(
        "Runtime resource configuration requires a work directory");
  }
  const auto memory = RaMAxMemory::readProcessMemorySnapshot();
  const uint64_t host_limit = effectiveHostLimit(memory);
  const uint64_t requested = parseMemorySize(config.memory_limit);

  if (requested != 0) {
    if (host_limit != 0 && requested > host_limit) {
      throw std::invalid_argument(
          "--memory-limit " + formatBytes(requested) +
          " exceeds the effective physical/cgroup limit " +
          formatBytes(host_limit));
    }
    memory_budget_bytes_ = requested;
  } else {
    if (host_limit == 0) {
      throw std::runtime_error(
          "Cannot detect physical or cgroup memory; specify --memory-limit");
    }
    const uint64_t percentage_reserve = host_limit * 15ULL / 100ULL;
    const uint64_t bounded_fixed_reserve =
        std::min<uint64_t>(64ULL * kGiB, host_limit / 4ULL);
    const uint64_t reserve =
        std::max(percentage_reserve, bounded_fixed_reserve);
    memory_budget_bytes_ = host_limit > reserve
        ? host_limit - reserve
        : host_limit * 3ULL / 4ULL;
  }
  if (memory_budget_bytes_ < 256ULL * kMiB) {
    throw std::invalid_argument(
        "Effective --memory-limit is below the 256 MiB minimum");
  }

  temp_directory_ = config.temp_dir.empty()
      ? config.work_dir / "perf-spill"
      : config.temp_dir;
  std::error_code error;
  std::filesystem::create_directories(temp_directory_, error);
  if (error || !std::filesystem::is_directory(temp_directory_)) {
    throw std::runtime_error(
        "Cannot create runtime temporary directory " +
        temp_directory_.string() + ": " + error.message());
  }
  const auto space = std::filesystem::space(temp_directory_, error);
  if (error) {
    throw std::runtime_error(
        "Cannot inspect free space for " + temp_directory_.string() +
        ": " + error.message());
  }
  if (space.available < 64ULL * kMiB) {
    throw std::runtime_error(
        "Runtime temporary directory has less than 64 MiB free: " +
        temp_directory_.string());
  }

  requested_threads_ = std::max<size_t>(1, config.requested_threads);
  configured_ = true;
}

uint64_t RuntimeResourceManager::accountedMemoryBytes(
    const RaMAxMemory::ProcessMemorySnapshot& memory) const noexcept {
  const uint64_t process_bytes = memory.rss_kib * kKiB;
  return std::max(process_bytes, memory.cgroup_current_bytes);
}

ResourceSnapshot RuntimeResourceManager::snapshot() const noexcept {
  ResourceSnapshot result;
  result.memory = RaMAxMemory::readProcessMemorySnapshot();
  result.accounted_memory_bytes = accountedMemoryBytes(result.memory);
  try {
    std::error_code error;
    const auto space = std::filesystem::space(temp_directory_, error);
    if (!error) {
      result.temp_available_bytes = space.available;
      result.temp_capacity_bytes = space.capacity;
    }
  } catch (...) {
  }
  result.pressure = pressure();
  return result;
}

MemoryPressure RuntimeResourceManager::pressure(
    uint64_t projected_bytes) const noexcept {
  if (!configured_ || memory_budget_bytes_ == 0) {
    return MemoryPressure::NORMAL;
  }
  const auto memory = RaMAxMemory::readProcessMemorySnapshot();
  const uint64_t projected = saturatedAdd(
      accountedMemoryBytes(memory), projected_bytes);
  const long double fraction =
      static_cast<long double>(projected) /
      static_cast<long double>(memory_budget_bytes_);
  if (fraction >= 0.95L) return MemoryPressure::CRITICAL;
  if (fraction >= 0.92L) return MemoryPressure::RECLAIM;
  if (fraction >= 0.85L) return MemoryPressure::SPILL;
  if (fraction >= 0.75L) return MemoryPressure::HOLD_BATCH;
  return MemoryPressure::NORMAL;
}

bool RuntimeResourceManager::shouldSpill(
    uint64_t projected_bytes) const noexcept {
  return pressure(projected_bytes) >= MemoryPressure::SPILL;
}

bool RuntimeResourceManager::shouldReclaim(
    uint64_t projected_bytes) const noexcept {
  return pressure(projected_bytes) >= MemoryPressure::RECLAIM;
}

size_t RuntimeResourceManager::boundedConcurrency(
    size_t requested_items, uint64_t bytes_per_item,
    uint64_t fixed_bytes) const noexcept {
  if (requested_items == 0) return 0;
  if (!configured_ || bytes_per_item == 0) {
    return std::min(requested_items, requested_threads_);
  }
  const auto memory = RaMAxMemory::readProcessMemorySnapshot();
  const uint64_t current = accountedMemoryBytes(memory);
  const uint64_t planning_limit = memory_budget_bytes_ * 75ULL / 100ULL;
  if (current >= planning_limit ||
      fixed_bytes >= planning_limit - current) {
    return 1;
  }
  const uint64_t available = planning_limit - current - fixed_bytes;
  const uint64_t by_memory = std::max<uint64_t>(1, available / bytes_per_item);
  return std::max<size_t>(
      1, std::min<size_t>({requested_items, requested_threads_,
                           static_cast<size_t>(std::min<uint64_t>(
                               by_memory,
                               std::numeric_limits<size_t>::max()))}));
}

void RuntimeResourceManager::requireAllocation(
    uint64_t projected_bytes, std::string_view stage) const {
  if (pressure(projected_bytes) != MemoryPressure::CRITICAL) return;
  const auto current = snapshot();
  throw std::runtime_error(
      "Memory budget would be exceeded before stage " +
      std::string(stage) + ": current=" +
      formatBytes(current.accounted_memory_bytes) + ", requested=" +
      formatBytes(projected_bytes) + ", budget=" +
      formatBytes(memory_budget_bytes_));
}

void RuntimeResourceManager::requireTempSpace(
    uint64_t projected_bytes, std::string_view stage) const {
  std::error_code error;
  const auto space = std::filesystem::space(temp_directory_, error);
  if (error) {
    throw std::runtime_error(
        "Cannot inspect temporary space before stage " +
        std::string(stage) + ": " + error.message());
  }
  const uint64_t safety_margin =
      std::max<uint64_t>(space.capacity / 20ULL, 64ULL * kMiB);
  if (space.available <= safety_margin ||
      projected_bytes > space.available - safety_margin) {
    throw std::runtime_error(
        "Temporary space would be exceeded before stage " +
        std::string(stage) + ": available=" +
        formatBytes(space.available) + ", requested=" +
        formatBytes(projected_bytes) + ", safety_margin=" +
        formatBytes(safety_margin) + ", path=" +
        temp_directory_.string());
  }
}

void RuntimeResourceManager::logConfiguration() const {
  const auto current = snapshot();
  spdlog::info(
      "[resource-budget] memory_budget_bytes={} memory_budget={} "
      "physical_bytes={} cgroup_limit_bytes={} cgroup_current_bytes={} "
      "cgroup_peak_bytes={} cgroup_swap_current_bytes={} "
      "temp_dir={} temp_available_bytes={} requested_threads={}",
      memory_budget_bytes_, formatBytes(memory_budget_bytes_),
      current.memory.physical_memory_bytes,
      current.memory.cgroup_limit_bytes,
      current.memory.cgroup_current_bytes,
      current.memory.cgroup_peak_bytes,
      current.memory.cgroup_swap_current_bytes,
      temp_directory_.string(), current.temp_available_bytes,
      requested_threads_);
}

void RuntimeResourceManager::logSnapshot(
    std::string_view stage, std::string_view event,
    uint64_t projected_bytes) const {
  const auto current = snapshot();
  const MemoryPressure projected_pressure = pressure(projected_bytes);
  spdlog::info(
      "[resource-budget] stage={} event={} pressure={} "
      "rss_kib={} peak_rss_kib={} cgroup_current_bytes={} "
      "cgroup_peak_bytes={} cgroup_swap_current_bytes={} "
      "budget_bytes={} projected_bytes={} temp_available_bytes={}",
      stage, event, pressureName(projected_pressure),
      current.memory.rss_kib, current.memory.peak_rss_kib,
      current.memory.cgroup_current_bytes,
      current.memory.cgroup_peak_bytes,
      current.memory.cgroup_swap_current_bytes,
      memory_budget_bytes_, projected_bytes,
      current.temp_available_bytes);
}

}  // namespace RaMAxResources
