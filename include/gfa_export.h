#ifndef RAMAX_GFA_EXPORT_H
#define RAMAX_GFA_EXPORT_H

#include <cstddef>
#include <cstdint>
#include <string>

namespace RaMesh::Gfa {

enum class Version {
    V1_0,
    V1_1
};

Version parseVersion(const std::string& value);
const char* versionString(Version version) noexcept;

struct GfaExportOptions {
    Version version{Version::V1_1};
    bool only_primary{true};
    int threads{1};
};

struct GfaExportStats {
    std::size_t blocks_seen{0};
    std::size_t blocks_used{0};
    std::size_t exact_runs{0};
    std::size_t atomic_intervals{0};
    std::size_t nodes{0};
    std::size_t edges{0};
    std::size_t walks{0};
    std::uint64_t graph_sequence_bp{0};
    std::uint64_t shared_sequence_bp{0};
    std::uint64_t private_sequence_bp{0};
    double elapsed_seconds{0.0};
};

}  // namespace RaMesh::Gfa

#endif
