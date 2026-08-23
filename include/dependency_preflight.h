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

// Throws one aggregated error that names every missing dependency.
void validateStartupDependencies(
    const StartupDependencies& dependencies);

StartupDependencies requireStartupDependencies();

}  // namespace RaMAxDependencies
