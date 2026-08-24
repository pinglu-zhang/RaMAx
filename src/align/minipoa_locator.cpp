#include "minipoa_locator.h"

#include <array>
#include <cstdlib>
#include <limits.h>
#include <string_view>
#include <unistd.h>

#ifndef RAMAX_MINIPOA_CONFIGURED_PATH
#define RAMAX_MINIPOA_CONFIGURED_PATH ""
#endif
#ifndef RAMAX_TOOL_BIN_CONFIGURED_PATH
#define RAMAX_TOOL_BIN_CONFIGURED_PATH ""
#endif

namespace RaMesh::Alignment {
namespace {

bool isExecutable(const std::filesystem::path& candidate) {
    std::error_code error;
    return !candidate.empty() &&
           std::filesystem::is_regular_file(candidate, error) && !error &&
           ::access(candidate.c_str(), X_OK) == 0;
}

std::filesystem::path executableDirectory() {
    std::array<char, PATH_MAX + 1> buffer{};
    const ssize_t length = ::readlink("/proc/self/exe", buffer.data(), PATH_MAX);
    if (length <= 0) return {};
    buffer[static_cast<size_t>(length)] = '\0';
    return std::filesystem::path(buffer.data()).parent_path();
}

std::filesystem::path searchPath() {
    const char* raw_path = std::getenv("PATH");
    if (!raw_path) return {};
    std::string_view path(raw_path);
    size_t offset = 0;
    while (offset <= path.size()) {
        const size_t separator = path.find(':', offset);
        const size_t end = separator == std::string_view::npos
            ? path.size() : separator;
        if (end > offset) {
            const auto candidate =
                std::filesystem::path(path.substr(offset, end - offset)) /
                "minipoa";
            if (isExecutable(candidate)) return std::filesystem::absolute(candidate);
        }
        if (separator == std::string_view::npos) break;
        offset = separator + 1;
    }
    return {};
}

}  // namespace

std::filesystem::path locateMinipoaExecutable() {
    const std::filesystem::path configured(RAMAX_MINIPOA_CONFIGURED_PATH);
    if (isExecutable(configured)) return configured;
    const std::filesystem::path configured_directory(
        RAMAX_TOOL_BIN_CONFIGURED_PATH);
    if (!configured_directory.empty()) {
        const auto candidate = configured_directory / "minipoa";
        if (isExecutable(candidate)) return candidate;
    }
    const auto sibling = executableDirectory() / "minipoa";
    if (isExecutable(sibling)) return sibling;
    return searchPath();
}

}  // namespace RaMesh::Alignment
