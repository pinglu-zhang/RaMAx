#ifndef RAMAX_EXTERNAL_TOOL_H
#define RAMAX_EXTERNAL_TOOL_H

#include <chrono>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace RaMAxExternalTool {

struct CommandResult {
    int exit_code{0};
    bool timed_out{false};
    int termination_signal{0};
    std::chrono::milliseconds elapsed{0};
};

struct RunOptions {
    // A non-positive timeout preserves the historical blocking behaviour.
    std::chrono::milliseconds timeout{0};
    std::chrono::milliseconds termination_grace{std::chrono::seconds(10)};
    std::chrono::milliseconds poll_interval{std::chrono::milliseconds(200)};
    bool create_process_group{false};
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

CommandResult run(
    const std::filesystem::path& executable,
    const std::vector<std::string>& arguments,
    const std::filesystem::path& stdout_path,
    const std::filesystem::path& stderr_path,
    const RunOptions& options);

std::string readText(const std::filesystem::path& path);

}  // namespace RaMAxExternalTool

#endif
