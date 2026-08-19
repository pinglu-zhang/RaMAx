#ifndef RAMAX_GFA_EXPORT_H
#define RAMAX_GFA_EXPORT_H

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>

namespace RaMesh::Gfa {

enum class Version {
    V1_0,
    V1_1
};

Version parseVersion(const std::string& value);
const char* versionString(Version version) noexcept;

enum class Profile {
    EXACT,
    COMPACT
};

Profile parseProfile(const std::string& value);
const char* profileString(Profile profile) noexcept;

// Internal compact-graph policy. These values are intentionally not exposed
// through the RaMAx command line or restart schema; they make the exporter
// testable and keep the selected compact profile reproducible.
struct GfaCompactConfig {
    std::uint64_t minimum_exact_run_bp{5};
    std::uint64_t short_node_bp{25};
    std::uint64_t stable_anchor_bp{200};
    std::uint64_t maximum_reference_span_bp{500};
    std::uint64_t maximum_path_span_bp{1500};
    std::size_t minimum_short_nodes{4};
    std::size_t maximum_distinct_alleles{16};
    double minimum_homology_retention{0.995};
    double maximum_graph_bp_ratio{1.02};
    bool enable_unitig_compaction{true};
    bool enable_compound_alleles{true};
};

struct GfaExportOptions {
    Version version{Version::V1_1};
    Profile profile{Profile::EXACT};
    bool only_primary{true};
    int threads{1};
    std::filesystem::path work_dir;
    GfaCompactConfig compact;
};

struct GfaExportStats {
    std::size_t blocks_seen{0};
    std::size_t blocks_used{0};
    std::size_t exact_runs{0};
    std::size_t atomic_intervals{0};
    std::size_t suppressed_exact_runs{0};
    std::size_t unitig_merges{0};
    std::size_t allele_islands{0};
    std::size_t rejected_allele_islands{0};
    std::size_t nodes{0};
    std::size_t edges{0};
    std::size_t walks{0};
    std::uint64_t graph_sequence_bp{0};
    std::uint64_t shared_sequence_bp{0};
    std::uint64_t private_sequence_bp{0};
    std::uint64_t homology_mass{0};
    std::uint64_t exact_homology_mass{0};
    double elapsed_seconds{0.0};
};

}  // namespace RaMesh::Gfa

#endif
