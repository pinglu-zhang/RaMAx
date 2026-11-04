#pragma once

#include <filesystem>

namespace RaMAxDependencies {

struct StartupDependencies {
    std::filesystem::path minipoa;
    std::filesystem::path wfmash;
    std::filesystem::path mash;
};

std::filesystem::path locateMinipoaExecutable();
std::filesystem::path locateWfmashExecutable();
std::filesystem::path locateMashExecutable();

StartupDependencies locateStartupDependencies();

// Throws one aggregated error for the three dependencies required by every
// normal run.
void validateUnconditionalStartupDependencies(
    const StartupDependencies& dependencies);

StartupDependencies requireUnconditionalStartupDependencies();

}  // namespace RaMAxDependencies
