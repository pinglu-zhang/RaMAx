#ifndef RAMAX_OUTPUT_FORMAT_HPP
#define RAMAX_OUTPUT_FORMAT_HPP

#include <filesystem>
#include <string>

enum class MultipleGenomeOutputFormat {
    HAL,
    MAF,
    PAF,
    UNKNOWN
};

inline MultipleGenomeOutputFormat detectMultipleGenomeOutputFormat(
    const std::filesystem::path& output_path) {
    const std::string ext = output_path.extension().string();
    if (ext == ".hal") return MultipleGenomeOutputFormat::HAL;
    if (ext == ".maf") return MultipleGenomeOutputFormat::MAF;
    if (ext == ".paf") return MultipleGenomeOutputFormat::PAF;
    return MultipleGenomeOutputFormat::UNKNOWN;
}

#endif
