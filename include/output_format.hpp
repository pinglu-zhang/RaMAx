#ifndef RAMAX_OUTPUT_FORMAT_HPP
#define RAMAX_OUTPUT_FORMAT_HPP

#include <filesystem>
#include <string>

enum class MultipleGenomeOutputFormat {
    HAL = 0,
    MAF = 1,
    PAF = 2,
    GFA = 3,
    UNKNOWN = 255
};

inline MultipleGenomeOutputFormat detectMultipleGenomeOutputFormat(
    const std::filesystem::path& output_path) {
    const std::string ext = output_path.extension().string();
    if (ext == ".hal") return MultipleGenomeOutputFormat::HAL;
    if (ext == ".maf") return MultipleGenomeOutputFormat::MAF;
    if (ext == ".paf") return MultipleGenomeOutputFormat::PAF;
    if (ext == ".gfa") return MultipleGenomeOutputFormat::GFA;
    return MultipleGenomeOutputFormat::UNKNOWN;
}

#endif
