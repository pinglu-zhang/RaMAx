#ifndef RAMAX_PROCESS_MEMORY_H
#define RAMAX_PROCESS_MEMORY_H

#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>

namespace RaMAxMemory {

struct ProcessMemorySnapshot {
  uint64_t rss_kib = 0;
  uint64_t peak_rss_kib = 0;
  uint64_t virtual_kib = 0;
  uint64_t cgroup_limit_bytes = 0;
  bool available = false;
};

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

    for (const char* path : {
             "/sys/fs/cgroup/memory.max",
             "/sys/fs/cgroup/memory/memory.limit_in_bytes"}) {
      std::ifstream limit_file(path);
      std::string limit;
      if (limit_file >> limit && limit != "max") {
        result.cgroup_limit_bytes = std::stoull(limit);
        if (result.cgroup_limit_bytes >= (1ULL << 60U)) {
          result.cgroup_limit_bytes = 0;
        }
        break;
      }
    }
  } catch (...) {
  }
#endif
  return result;
}

}  // namespace RaMAxMemory

#endif  // RAMAX_PROCESS_MEMORY_H
