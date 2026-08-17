#ifndef RAMAX_CACHE_MANIFEST_H
#define RAMAX_CACHE_MANIFEST_H

#include "config.hpp"

#include <cstdint>
#include <string>

namespace RaMAxCache {

struct FileMetadata {
    uint64_t size{0};
    int64_t mtime{0};

    bool operator==(const FileMetadata&) const = default;
};

struct StageStats {
    size_t reused{0};
    size_t rebuilt{0};
};

FileMetadata fileMetadata(const FilePath& path);
FilePath completionMarkerPath(const FilePath& output);

bool markerMatches(const FilePath& marker,
                   const std::string& kind,
                   uint32_t cache_version,
                   const std::string& source_identity,
                   bool source_is_url,
                   const FilePath& source_path,
                   const FilePath& output_path);

void writeMarker(const FilePath& marker,
                 const std::string& kind,
                 uint32_t cache_version,
                 const std::string& source_identity,
                 bool source_is_url,
                 const FilePath& source_path,
                 const FilePath& output_path);

void publishFile(const FilePath& partial, const FilePath& final_path);
void removeIfPresent(const FilePath& path) noexcept;

}  // namespace RaMAxCache

#endif
