#ifndef RAMAX_MM2PLUS_ROUTER_H
#define RAMAX_MM2PLUS_ROUTER_H

#include "SeqPro.h"
#include "anchor.h"
#include "mash_distance_estimator.h"
#include "wfmash_router.h"

#include <filesystem>
#include <map>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

struct FirstRoundMm2plusResult {
    std::unordered_set<SpeciesName> successful_species;
    std::map<SpeciesName, AnchorVec> anchors_by_species;
};

namespace Mm2plusRouterDetail {

struct NormalizationStats {
    size_t raw_records{0};
    size_t quality_discarded{0};
    size_t low_score_discarded{0};
    size_t low_identity_discarded{0};
    size_t exact_duplicates{0};
    size_t trimmed_records{0};
    size_t overlap_discarded{0};
    size_t short_discarded{0};
    size_t canonical_records{0};
};

bool passesGraphQualityFilters(
    const WfmashRouterDetail::ParsedPafRecord& record);

std::string validateMm2plusVersion(
    std::string_view version_output,
    std::string_view help_output);

std::vector<std::string> indexArguments(
    uint_t threads,
    uint64_t reference_bases,
    const std::filesystem::path& index,
    const std::filesystem::path& reference);

std::vector<std::string> alignmentArguments(
    uint_t threads,
    const std::filesystem::path& index,
    const std::filesystem::path& query);

uint_t threadsPerPair(uint_t total_threads);
size_t workerCount(size_t tasks, uint_t total_threads);

NormalizationStats normalizeForGraph(
    std::vector<WfmashRouterDetail::ParsedPafRecord>& records,
    uint_t min_span);

}  // namespace Mm2plusRouterDetail

class FirstRoundMm2plusRouter {
public:
    FirstRoundMm2plusRouter(
        std::filesystem::path executable,
        std::filesystem::path output_directory,
        uint_t threads);

    FirstRoundMm2plusResult run(
        const SpeciesName& reference,
        const std::vector<MashDistanceRecord>& distances,
        double far_distance,
        uint_t min_span,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers);

    const std::string& version() const { return version_; }

private:
    std::filesystem::path executable_;
    std::filesystem::path output_directory_;
    uint_t threads_{1};
    std::string version_;
};

std::filesystem::path locateMm2plusExecutable();

#endif
