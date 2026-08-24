#ifndef RAMAX_MASH_DISTANCE_ESTIMATOR_H
#define RAMAX_MASH_DISTANCE_ESTIMATOR_H

#include "SeqPro.h"
#include "config.hpp"

#include <filesystem>
#include <map>
#include <string>
#include <string_view>
#include <vector>

struct MashDistanceRecord {
    SpeciesName reference;
    SpeciesName query;
    double distance{0.0};
    double p_value{0.0};
    uint64_t shared_hashes{0};
    uint64_t total_hashes{0};
    FilePath reference_fasta;
    FilePath query_fasta;
};

namespace MashDistanceDetail {

struct ParsedLine {
    std::string reference_id;
    std::string query_id;
    double distance{0.0};
    double p_value{0.0};
    uint64_t shared_hashes{0};
    uint64_t total_hashes{0};
};

ParsedLine parseLine(std::string_view line);
std::string validateVersion(std::string_view output);

}  // namespace MashDistanceDetail

class MashDistanceEstimator {
public:
    static constexpr uint_t kKmerSize = 31;
    static constexpr uint_t kSketchSize = 20000;

    MashDistanceEstimator(std::filesystem::path executable,
                          std::filesystem::path output_directory,
                          uint_t threads);

    std::vector<MashDistanceRecord> estimateFirstReference(
        const SpeciesName& reference,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers);

    const std::string& version() const { return version_; }
    const std::filesystem::path& executable() const { return executable_; }

private:
    std::filesystem::path executable_;
    std::filesystem::path output_directory_;
    uint_t threads_{1};
    std::string version_;
};

std::filesystem::path locateMashExecutable();

#endif
