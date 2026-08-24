#pragma once

#include <filesystem>

namespace RaMAxDependencies {

struct StartupDependencies {
    std::filesystem::path hal_append_cactus_subtree;
    std::filesystem::path minipoa;
    std::filesystem::path wfmash;
    std::filesystem::path mash;
};

std::filesystem::path locateHalAppendCactusSubtreeExecutable();
std::filesystem::path locateMinipoaExecutable();
std::filesystem::path locateWfmashExecutable();
std::filesystem::path locateMashExecutable();

StartupDependencies locateStartupDependencies();

// Throws one aggregated error for the three dependencies required by every
// normal run. halAppendCactusSubtree is output-dependent and checked
// separately after the effective output list is known.
void validateUnconditionalStartupDependencies(
    const StartupDependencies& dependencies);

void validateHalAppendCactusSubtree(
    const std::filesystem::path& executable,
    bool hal_output_requested);

StartupDependencies requireUnconditionalStartupDependencies();

}  // namespace RaMAxDependencies
