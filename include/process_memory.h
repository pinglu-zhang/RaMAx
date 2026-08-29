#ifndef RAMAX_PROCESS_MEMORY_H
#define RAMAX_PROCESS_MEMORY_H

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

#if defined(__linux__)
#include <unistd.h>
#endif

namespace RaMAxMemory {

struct ProcessMemorySnapshot {
  uint64_t rss_kib = 0;
  uint64_t peak_rss_kib = 0;
  uint64_t virtual_kib = 0;
  uint64_t physical_memory_bytes = 0;
  uint64_t cgroup_current_bytes = 0;
  uint64_t cgroup_peak_bytes = 0;
  uint64_t cgroup_limit_bytes = 0;
  uint64_t cgroup_swap_current_bytes = 0;
  uint64_t cgroup_swap_limit_bytes = 0;
  bool available = false;
};

namespace detail {

inline bool readUnsignedFile(const std::filesystem::path& path,
                             uint64_t& value) noexcept {
  try {
    std::ifstream input(path);
    std::string text;
    if (!(input >> text) || text == "max") return false;
    size_t consumed = 0;
    const uint64_t parsed = std::stoull(text, &consumed);
    if (consumed != text.size()) return false;
    value = parsed;
    return true;
  } catch (...) {
    return false;
  }
}

inline std::filesystem::path cgroupV2Directory() noexcept {
#if defined(__linux__)
  try {
    std::ifstream cgroup("/proc/self/cgroup");
    std::string line;
    while (std::getline(cgroup, line)) {
      if (!line.starts_with("0::")) continue;
      std::string relative = line.substr(3);
      while (!relative.empty() && relative.front() == '/') {
        relative.erase(relative.begin());
      }
      return std::filesystem::path("/sys/fs/cgroup") / relative;
    }
  } catch (...) {
  }
#endif
  return "/sys/fs/cgroup";
}

inline std::filesystem::path cgroupV1MemoryDirectory() noexcept {
#if defined(__linux__)
  try {
    std::ifstream cgroup("/proc/self/cgroup");
    std::string line;
    while (std::getline(cgroup, line)) {
      const size_t first = line.find(':');
      const size_t second = first == std::string::npos
          ? std::string::npos
          : line.find(':', first + 1);
      if (second == std::string::npos) continue;
      const std::string controllers = line.substr(first + 1,
                                                   second - first - 1);
      if (controllers.find("memory") == std::string::npos) continue;
      std::string relative = line.substr(second + 1);
      while (!relative.empty() && relative.front() == '/') {
        relative.erase(relative.begin());
      }
      return std::filesystem::path("/sys/fs/cgroup/memory") / relative;
    }
  } catch (...) {
  }
#endif
  return "/sys/fs/cgroup/memory";
}

inline uint64_t usableCgroupLimit(uint64_t value) noexcept {
  // cgroup v1 commonly exposes an effectively-unlimited value near LONG_MAX.
  return value >= (1ULL << 60U) ? 0 : value;
}

}  // namespace detail

inline ProcessMemorySnapshot readProcessMemorySnapshot() noexcept {
  ProcessMemorySnapshot result;
#if defined(__linux__)
  try {
    std::ifstream status("/proc/self/status");
    result.available = status.is_open();
    std::string line;
    while (std::getline(status, line)) {
      std::istringstream fields(line);
      std::string key;
      uint64_t value = 0;
      std::string unit;
      if (!(fields >> key >> value >> unit)) {
        continue;
      }
      if (key == "VmRSS:") {
        result.rss_kib = value;
      } else if (key == "VmHWM:") {
        result.peak_rss_kib = value;
      } else if (key == "VmSize:") {
        result.virtual_kib = value;
      }
    }

    const long pages = ::sysconf(_SC_PHYS_PAGES);
    const long page_size = ::sysconf(_SC_PAGESIZE);
    if (pages > 0 && page_size > 0 &&
        static_cast<uint64_t>(pages) <=
            std::numeric_limits<uint64_t>::max() /
                static_cast<uint64_t>(page_size)) {
      result.physical_memory_bytes = static_cast<uint64_t>(pages) *
                                     static_cast<uint64_t>(page_size);
    }

    const auto unified = detail::cgroupV2Directory();
    uint64_t value = 0;
    if (detail::readUnsignedFile(unified / "memory.current", value)) {
      result.cgroup_current_bytes = value;
      if (detail::readUnsignedFile(unified / "memory.peak", value)) {
        result.cgroup_peak_bytes = value;
      }
      if (detail::readUnsignedFile(unified / "memory.max", value)) {
        result.cgroup_limit_bytes = detail::usableCgroupLimit(value);
      }
      if (detail::readUnsignedFile(unified / "memory.swap.current", value)) {
        result.cgroup_swap_current_bytes = value;
      }
      if (detail::readUnsignedFile(unified / "memory.swap.max", value)) {
        result.cgroup_swap_limit_bytes = detail::usableCgroupLimit(value);
      }
    } else {
      const auto legacy = detail::cgroupV1MemoryDirectory();
      if (detail::readUnsignedFile(legacy / "memory.usage_in_bytes", value)) {
        result.cgroup_current_bytes = value;
      }
      if (detail::readUnsignedFile(legacy / "memory.max_usage_in_bytes", value)) {
        result.cgroup_peak_bytes = value;
      }
      if (detail::readUnsignedFile(legacy / "memory.limit_in_bytes", value)) {
        result.cgroup_limit_bytes = detail::usableCgroupLimit(value);
      }
      if (detail::readUnsignedFile(legacy / "memory.memsw.usage_in_bytes", value) &&
          value > result.cgroup_current_bytes) {
        result.cgroup_swap_current_bytes = value - result.cgroup_current_bytes;
      }
    }
  } catch (...) {
  }
#endif
  return result;
}

}  // namespace RaMAxMemory

#endif  // RAMAX_PROCESS_MEMORY_H
