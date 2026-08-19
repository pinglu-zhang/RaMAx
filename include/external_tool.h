#ifndef RAMAX_EXTERNAL_TOOL_H
#define RAMAX_EXTERNAL_TOOL_H

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace RaMAxExternalTool {

struct CommandResult {
    int exit_code{0};
    double wall_seconds{0.0};
    double user_seconds{0.0};
    double system_seconds{0.0};
    long peak_rss_kb{0};
};

bool isExecutable(const std::filesystem::path& candidate);
std::filesystem::path executableDirectory();
std::filesystem::path searchPath(std::string_view name);
std::filesystem::path locateExecutable(
    std::string_view name,
    const std::filesystem::path& configured_path);

CommandResult run(
    const std::filesystem::path& executable,
    const std::vector<std::string>& arguments,
    const std::filesystem::path& stdout_path,
    const std::filesystem::path& stderr_path);

std::string readText(const std::filesystem::path& path);

}  // namespace RaMAxExternalTool

#endif
