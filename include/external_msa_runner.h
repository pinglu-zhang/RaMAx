#ifndef RAMAX_EXTERNAL_MSA_RUNNER_H
#define RAMAX_EXTERNAL_MSA_RUNNER_H

#include "config.hpp"

#include <filesystem>
#include <string>
#include <unordered_map>

namespace RaMesh::Alignment {

class ExternalMsaRunner {
public:
    static ExternalMsaRunner& instance();

    void configureDefaultExecutable(std::string executable);
    void configureScratchDirectory(std::filesystem::path directory);
    bool align(const std::string& executable,
               std::unordered_map<ChrName, std::string>& sequences);
    bool alignWithDefault(
        std::unordered_map<ChrName, std::string>& sequences);
    void logSummary() const;

private:
    ExternalMsaRunner() = default;
};

}  // namespace RaMesh::Alignment

#endif
