#ifndef RAMAX_WFMASH_ROUTER_H
#define RAMAX_WFMASH_ROUTER_H

#include "SeqPro.h"
#include "anchor.h"
#include "mash_distance_estimator.h"

#include <chrono>
#include <filesystem>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

struct FirstRoundWfmashResult {
    std::unordered_set<SpeciesName> successful_species;
    std::map<SpeciesName, AnchorVec> anchors_by_species;
};

namespace WfmashRouterDetail {

struct SequenceRecord {
    std::string name;
    uint64_t length{0};

    bool operator==(const SequenceRecord&) const = default;
};

struct FaiRecord {
    std::string name;
    uint64_t length{0};
    uint64_t offset{0};
    uint64_t line_bases{0};
    uint64_t line_width{0};
};

struct ParsedPafRecord {
    std::string query_name;
    uint64_t query_length{0};
    uint64_t query_start{0};
    uint64_t query_end{0};
    Strand strand{Strand::FORWARD};
    std::string target_name;
    uint64_t target_length{0};
    uint64_t target_start{0};
    uint64_t target_end{0};
    uint64_t matches{0};
    uint64_t block_length{0};
    uint64_t cigar_columns{0};
    uint64_t mapq{0};
    Cigar_t cigar;
    std::string cigar_text;
};

struct PafNormalizationStats {
    size_t input_records{0};
    size_t trimmed_records{0};
    size_t skipped_records{0};
};

struct PairThreadSchedule {
    std::vector<uint_t> threads_per_worker;

    size_t workers() const { return threads_per_worker.size(); }
};

struct ExecutionPolicy {
    uint_t minimum_threads_per_process{4};
    std::chrono::milliseconds pair_timeout{std::chrono::hours(1)};
    std::chrono::milliseconds termination_grace{std::chrono::seconds(10)};
    std::chrono::milliseconds poll_interval{std::chrono::milliseconds(200)};
};

std::string validateSamtoolsVersion(std::string_view output);
std::string validateWfmashVersion(std::string_view output);
std::vector<FaiRecord> parseFai(std::istream& input);
ParsedPafRecord parsePafLine(std::string_view line, bool require_cigar);
PafNormalizationStats normalizePafForGraph(
    std::vector<ParsedPafRecord>& records);
std::vector<std::string> mappingArguments(
    uint_t threads, const std::filesystem::path& tmp_directory,
    const std::filesystem::path& reference,
    const std::filesystem::path& query);
std::vector<std::string> alignmentArguments(
    uint_t threads, const std::filesystem::path& tmp_directory,
    const std::filesystem::path& mapping,
    const std::filesystem::path& reference,
    const std::filesystem::path& query);
size_t workerCount(size_t tasks, uint_t total_threads);
uint_t threadsPerTask(size_t tasks, uint_t total_threads);
PairThreadSchedule pairThreadSchedule(
    size_t tasks, uint_t total_threads,
    uint_t minimum_threads_per_process = 4);

}  // namespace WfmashRouterDetail

class FirstRoundWfmashRouter {
public:
    FirstRoundWfmashRouter(
        std::filesystem::path samtools_executable,
        std::filesystem::path wfmash_executable,
        std::filesystem::path output_directory,
        uint_t threads,
        WfmashRouterDetail::ExecutionPolicy execution_policy = {});

    FirstRoundWfmashResult run(
        const SpeciesName& reference,
        const std::vector<MashDistanceRecord>& distances,
        double near_distance,
        double far_distance,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers);

    const std::string& samtoolsVersion() const { return samtools_version_; }
    const std::string& wfmashVersion() const { return wfmash_version_; }

private:
    std::filesystem::path samtools_executable_;
    std::filesystem::path wfmash_executable_;
    std::filesystem::path output_directory_;
    uint_t threads_{1};
    WfmashRouterDetail::ExecutionPolicy execution_policy_;
    std::string samtools_version_;
    std::string wfmash_version_;
};

std::filesystem::path locateSamtoolsExecutable();
std::filesystem::path locateWfmashExecutable();

#endif
