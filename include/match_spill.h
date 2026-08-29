#ifndef RAMAX_MATCH_SPILL_H
#define RAMAX_MATCH_SPILL_H

#include "anchor.h"

#include <cstdint>
#include <filesystem>
#include <string>

namespace RaMAxSpill {

enum class MatchSpillKind : uint32_t {
  PRIMARY = 1,
  REPEAT_MASKED = 2,
  REPEAT_FULL = 3,
  SECONDARY_COMBINED = 4
};

struct MatchSpillIdentity {
  uint64_t round{0};
  uint64_t query_ordinal{0};
  MatchSpillKind kind{MatchSpillKind::PRIMARY};
  std::string reference;
  std::string query;
};

struct MatchSpillStats {
  uint64_t outer_vectors{0};
  uint64_t inner_vectors{0};
  uint64_t matches{0};
  uint64_t bytes{0};
};

class MatchSpillStore {
 public:
  explicit MatchSpillStore(std::filesystem::path directory);

  [[nodiscard]] const std::filesystem::path& directory() const noexcept {
    return directory_;
  }

  [[nodiscard]] MatchSpillStats inspect(const MatchVec3D& matches) const;

  std::filesystem::path write(const MatchSpillIdentity& identity,
                              const MatchVec3D& matches) const;

  MatchVec3DPtr read(const std::filesystem::path& path,
                     const MatchSpillIdentity& expected) const;

  void consume(const std::filesystem::path& path) const noexcept;

 private:
  std::filesystem::path pathFor(const MatchSpillIdentity& identity) const;

  std::filesystem::path directory_;
};

}  // namespace RaMAxSpill

#endif  // RAMAX_MATCH_SPILL_H
