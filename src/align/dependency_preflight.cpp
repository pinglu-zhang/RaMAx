#include "dependency_preflight.h"

#include "external_tool.h"

#include <array>
#include <sstream>
#include <stdexcept>
#include <string_view>

#ifndef RAMAX_HAL_APPEND_CACTUS_SUBTREE_CONFIGURED_PATH
#define RAMAX_HAL_APPEND_CACTUS_SUBTREE_CONFIGURED_PATH ""
#endif

#ifndef RAMAX_MINIPOA_CONFIGURED_PATH
#define RAMAX_MINIPOA_CONFIGURED_PATH ""
#endif

#ifndef RAMAX_WFMASH_CONFIGURED_PATH
#define RAMAX_WFMASH_CONFIGURED_PATH ""
#endif

#ifndef RAMAX_MASH_CONFIGURED_PATH
#define RAMAX_MASH_CONFIGURED_PATH ""
#endif

namespace RaMAxDependencies {

std::filesystem::path locateHalAppendCactusSubtreeExecutable() {
    return RaMAxExternalTool::locateExecutable(
        "halAppendCactusSubtree",
        RAMAX_HAL_APPEND_CACTUS_SUBTREE_CONFIGURED_PATH);
}

std::filesystem::path locateMinipoaExecutable() {
    return RaMAxExternalTool::locateExecutable(
        "minipoa", RAMAX_MINIPOA_CONFIGURED_PATH);
}

std::filesystem::path locateWfmashExecutable() {
    return RaMAxExternalTool::locateExecutable(
        "wfmash", RAMAX_WFMASH_CONFIGURED_PATH);
}

std::filesystem::path locateMashExecutable() {
    return RaMAxExternalTool::locateExecutable(
        "mash", RAMAX_MASH_CONFIGURED_PATH);
}

StartupDependencies locateStartupDependencies() {
    return {
        .hal_append_cactus_subtree =
            locateHalAppendCactusSubtreeExecutable(),
        .minipoa = locateMinipoaExecutable(),
        .wfmash = locateWfmashExecutable(),
        .mash = locateMashExecutable(),
    };
}

void validateStartupDependencies(
    const StartupDependencies& dependencies) {
    struct Requirement {
        std::string_view executable;
        std::string_view cmake_variable;
        const std::filesystem::path* resolved;
    };
    const std::array requirements{
        Requirement{
            "halAppendCactusSubtree",
            "RAMAX_HAL_APPEND_CACTUS_SUBTREE_EXECUTABLE",
            &dependencies.hal_append_cactus_subtree},
        Requirement{
            "minipoa",
            "RAMAX_MINIPOA_EXECUTABLE",
            &dependencies.minipoa},
        Requirement{
            "wfmash",
            "RAMAX_WFMASH_EXECUTABLE",
            &dependencies.wfmash},
        Requirement{
            "mash",
            "RAMAX_MASH_EXECUTABLE",
            &dependencies.mash},
    };

    std::ostringstream error;
    size_t missing = 0;
    for (const auto& requirement : requirements) {
        if (!RaMAxExternalTool::isExecutable(
                *requirement.resolved)) {
            if (missing == 0) {
                error
                    << "RaMAx startup dependency check failed. "
                    << "Missing required executable(s):";
            }
            ++missing;
            error << "\n  - " << requirement.executable
                  << " (configure with -D"
                  << requirement.cmake_variable << "=<path>)";
        }
    }
    if (missing != 0) {
        error
            << "\nSearch order: CMake-configured path, directory "
            << "containing the ramax executable, then PATH.";
        throw std::runtime_error(error.str());
    }
}

StartupDependencies requireStartupDependencies() {
    StartupDependencies dependencies =
        locateStartupDependencies();
    validateStartupDependencies(dependencies);
    return dependencies;
}

}  // namespace RaMAxDependencies
