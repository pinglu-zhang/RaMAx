#ifndef RAMAX_RUNTIME_RESOURCES_H
#define RAMAX_RUNTIME_RESOURCES_H

#include "process_memory.h"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>

namespace RaMAxResources {

struct RuntimeResourceConfig {
  std::string memory_limit{"auto"};
  std::filesystem::path temp_dir;
  std::filesystem::path work_dir;
  size_t requested_threads{1};
};

enum class MemoryPressure : uint8_t {
  NORMAL = 0,
  HOLD_BATCH = 1,
  SPILL = 2,
  RECLAIM = 3,
  CRITICAL = 4
};

struct ResourceSnapshot {
  RaMAxMemory::ProcessMemorySnapshot memory;
  uint64_t temp_available_bytes{0};
  uint64_t temp_capacity_bytes{0};
  uint64_t accounted_memory_bytes{0};
  MemoryPressure pressure{MemoryPressure::NORMAL};
};

uint64_t parseMemorySize(std::string_view text);
std::string formatBytes(uint64_t bytes);

class RuntimeResourceManager {
 public:
  static RuntimeResourceManager& instance();

  void configure(const RuntimeResourceConfig& config);

  [[nodiscard]] bool configured() const noexcept { return configured_; }
  [[nodiscard]] uint64_t memoryBudgetBytes() const noexcept {
    return memory_budget_bytes_;
  }
  [[nodiscard]] const std::filesystem::path& tempDirectory() const noexcept {
    return temp_directory_;
  }
  [[nodiscard]] size_t requestedThreads() const noexcept {
    return requested_threads_;
  }

  [[nodiscard]] ResourceSnapshot snapshot() const noexcept;
  [[nodiscard]] MemoryPressure pressure(uint64_t projected_bytes = 0) const noexcept;
  [[nodiscard]] bool shouldSpill(uint64_t projected_bytes = 0) const noexcept;
  [[nodiscard]] bool shouldReclaim(uint64_t projected_bytes = 0) const noexcept;

  [[nodiscard]] size_t boundedConcurrency(size_t requested_items,
                                          uint64_t bytes_per_item,
                                          uint64_t fixed_bytes = 0) const noexcept;

  void requireAllocation(uint64_t projected_bytes,
                         std::string_view stage) const;
  void requireTempSpace(uint64_t projected_bytes,
                        std::string_view stage) const;
  void logConfiguration() const;
  void logSnapshot(std::string_view stage,
                   std::string_view event,
                   uint64_t projected_bytes = 0) const;

 private:
  RuntimeResourceManager() = default;

  uint64_t accountedMemoryBytes(
      const RaMAxMemory::ProcessMemorySnapshot& memory) const noexcept;

  bool configured_{false};
  uint64_t memory_budget_bytes_{0};
  std::filesystem::path temp_directory_;
  size_t requested_threads_{1};
};

}  // namespace RaMAxResources

#endif  // RAMAX_RUNTIME_RESOURCES_H
