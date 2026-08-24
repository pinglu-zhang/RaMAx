#ifndef RAMAX_OUTPUT_SPEC_HPP
#define RAMAX_OUTPUT_SPEC_HPP

#include "output_format.hpp"

#include <algorithm>
#include <filesystem>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace RaMAxOutput {

struct OutputSpec {
    std::filesystem::path path;
    MultipleGenomeOutputFormat format{MultipleGenomeOutputFormat::UNKNOWN};
};

inline const char* formatName(MultipleGenomeOutputFormat format) {
    switch (format) {
    case MultipleGenomeOutputFormat::MAF: return "MAF";
    case MultipleGenomeOutputFormat::PAF: return "PAF";
    case MultipleGenomeOutputFormat::GFA: return "GFA";
    case MultipleGenomeOutputFormat::HAL: return "HAL";
    case MultipleGenomeOutputFormat::UNKNOWN: return "UNKNOWN";
    }
    return "UNKNOWN";
}

inline int exportOrder(MultipleGenomeOutputFormat format) {
    switch (format) {
    case MultipleGenomeOutputFormat::MAF: return 0;
    case MultipleGenomeOutputFormat::PAF: return 1;
    case MultipleGenomeOutputFormat::GFA: return 2;
    case MultipleGenomeOutputFormat::HAL: return 3;
    case MultipleGenomeOutputFormat::UNKNOWN: return 4;
    }
    return 3;
}

inline std::vector<OutputSpec> validateOutputPaths(
    const std::vector<std::filesystem::path>& paths) {
    if (paths.empty()) {
        throw std::invalid_argument(
            "At least one output path is required: --output (-o)");
    }

    std::set<std::filesystem::path> normalized_paths;
    std::set<MultipleGenomeOutputFormat> formats;
    std::vector<OutputSpec> outputs;
    outputs.reserve(paths.size());

    for (const auto& path : paths) {
        if (path.empty()) {
            throw std::invalid_argument("Output path must not be empty");
        }
        const auto format = detectMultipleGenomeOutputFormat(path);
        if (format == MultipleGenomeOutputFormat::UNKNOWN) {
            throw std::invalid_argument(
                "Invalid output file extension for '" + path.string() +
                "'. Supported: .hal, .maf, .paf, .gfa");
        }

        const std::filesystem::path normalized =
            std::filesystem::absolute(path).lexically_normal();
        if (!normalized_paths.insert(normalized).second) {
            throw std::invalid_argument(
                "Duplicate output path: " + path.string());
        }
        if (!formats.insert(format).second) {
            throw std::invalid_argument(
                "Duplicate output format: " +
                std::string(formatName(format)));
        }
        outputs.push_back({path, format});
    }

    std::stable_sort(outputs.begin(), outputs.end(),
        [](const OutputSpec& left, const OutputSpec& right) {
            return exportOrder(left.format) < exportOrder(right.format);
        });
    return outputs;
}

inline bool hasFormat(const std::vector<OutputSpec>& outputs,
                      MultipleGenomeOutputFormat format) {
    return std::any_of(outputs.begin(), outputs.end(),
        [format](const OutputSpec& output) {
            return output.format == format;
        });
}

}  // namespace RaMAxOutput

#endif
