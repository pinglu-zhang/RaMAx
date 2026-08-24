#include "cache_manifest.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace RaMAxCache {
namespace {

constexpr uint32_t kMarkerSchema = 1;

struct MarkerData {
    uint32_t schema{0};
    std::string kind;
    uint32_t cache_version{0};
    std::string source_identity;
    bool source_is_url{false};
    FileMetadata source;
    FileMetadata output;
};

bool readMarker(const FilePath& marker, MarkerData& data) {
    std::ifstream input(marker, std::ios::binary);
    if (!input) return false;

    std::string magic;
    int source_is_url = 0;
    if (!(input >> magic >> data.schema) || magic != "RAMAX_CACHE" ||
        !(input >> std::quoted(data.kind)) ||
        !(input >> data.cache_version) ||
        !(input >> std::quoted(data.source_identity)) ||
        !(input >> source_is_url) ||
        !(input >> data.source.size >> data.source.mtime) ||
        !(input >> data.output.size >> data.output.mtime)) {
        return false;
    }
    input >> std::ws;
    data.source_is_url = source_is_url != 0;
    return input.eof();
}

}  // namespace

FileMetadata fileMetadata(const FilePath& path) {
    std::error_code error;
    if (!std::filesystem::is_regular_file(path, error) || error) {
        throw std::runtime_error("Cache artifact is not a regular file: " +
                                 path.string());
    }
    const auto size = std::filesystem::file_size(path, error);
    if (error) {
        throw std::runtime_error("Cannot read file size for " + path.string() +
                                 ": " + error.message());
    }
    const auto mtime = std::filesystem::last_write_time(path, error);
    if (error) {
        throw std::runtime_error("Cannot read file mtime for " + path.string() +
                                 ": " + error.message());
    }
    return {static_cast<uint64_t>(size),
            static_cast<int64_t>(mtime.time_since_epoch().count())};
}

FilePath completionMarkerPath(const FilePath& output) {
    FilePath marker = output;
    marker += ".ramax-complete";
    return marker;
}

bool markerMatches(const FilePath& marker,
                   const std::string& kind,
                   uint32_t cache_version,
                   const std::string& source_identity,
                   bool source_is_url,
                   const FilePath& source_path,
                   const FilePath& output_path) {
    try {
        MarkerData data;
        if (!readMarker(marker, data) || data.schema != kMarkerSchema ||
            data.kind != kind || data.cache_version != cache_version ||
            data.source_identity != source_identity ||
            data.source_is_url != source_is_url) {
            return false;
        }
        if (!source_is_url && data.source != fileMetadata(source_path)) {
            return false;
        }
        return data.output == fileMetadata(output_path);
    } catch (const std::exception&) {
        return false;
    }
}

void writeMarker(const FilePath& marker,
                 const std::string& kind,
                 uint32_t cache_version,
                 const std::string& source_identity,
                 bool source_is_url,
                 const FilePath& source_path,
                 const FilePath& output_path) {
    const FileMetadata source = source_is_url
        ? FileMetadata{}
        : fileMetadata(source_path);
    const FileMetadata output = fileMetadata(output_path);
    FilePath partial = marker;
    partial += ".partial";
    removeIfPresent(partial);

    std::ofstream stream(partial, std::ios::binary | std::ios::trunc);
    if (!stream) {
        throw std::runtime_error("Cannot create cache marker: " + partial.string());
    }
    stream << "RAMAX_CACHE " << kMarkerSchema << '\n'
           << std::quoted(kind) << '\n'
           << cache_version << '\n'
           << std::quoted(source_identity) << '\n'
           << (source_is_url ? 1 : 0) << '\n'
           << source.size << ' ' << source.mtime << '\n'
           << output.size << ' ' << output.mtime << '\n';
    stream.flush();
    if (!stream) {
        removeIfPresent(partial);
        throw std::runtime_error("Failed to finalize cache marker: " +
                                 partial.string());
    }
    stream.close();
    publishFile(partial, marker);
}

void publishFile(const FilePath& partial, const FilePath& final_path) {
    std::error_code error;
    std::filesystem::rename(partial, final_path, error);
    if (!error) return;

    // Windows-mounted filesystems may not replace an existing destination.
    std::filesystem::remove(final_path, error);
    error.clear();
    std::filesystem::rename(partial, final_path, error);
    if (error) {
        throw std::runtime_error("Cannot publish " + final_path.string() +
                                 ": " + error.message());
    }
}

void removeIfPresent(const FilePath& path) noexcept {
    std::error_code error;
    std::filesystem::remove(path, error);
}

}  // namespace RaMAxCache
