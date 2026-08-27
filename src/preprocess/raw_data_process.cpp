#include "SeqPro.h"
#include "data_process.h"

#include <algorithm>
#include <exception>
#include <mutex>
#include <vector>

#include <omp.h>

namespace {

template<class Function>
void parallelForIndexed(size_t count, int requested_threads, Function&& function) {
    if (count == 0) return;
    const int workers = omp_in_parallel()
        ? 1
        : std::max(1, std::min<int>(requested_threads,
              static_cast<int>(count)));
    std::mutex failure_mutex;
    std::exception_ptr first_failure;
    size_t first_failure_index = count;
#pragma omp parallel for schedule(dynamic, 1) num_threads(workers) if(workers > 1)
    for (long long raw_index = 0;
         raw_index < static_cast<long long>(count); ++raw_index) {
        const size_t index = static_cast<size_t>(raw_index);
        try {
            function(index);
        } catch (...) {
            std::lock_guard<std::mutex> lock(failure_mutex);
            if (index < first_failure_index) {
                first_failure_index = index;
                first_failure = std::current_exception();
            }
        }
    }
    if (first_failure) std::rethrow_exception(first_failure);
}

} // namespace

// -----------------------------
// 查找 windowmasker 可执行文件路径
// -----------------------------
std::string findWindowmaskerPath() {
    // 1. 首先尝试在系统 PATH 中查找
    if (std::system("which windowmasker > /dev/null 2>&1") == 0) {
        spdlog::info("Found windowmasker in system PATH");
        return "windowmasker";
    }

    // 2. 尝试相对于当前可执行文件的路径（安装后的情况）
    try {
        std::filesystem::path exe_path = std::filesystem::canonical("/proc/self/exe");
        std::filesystem::path exe_dir = exe_path.parent_path();
        std::filesystem::path windowmasker_path = exe_dir / "windowmasker";

        if (std::filesystem::exists(windowmasker_path) &&
            std::filesystem::is_regular_file(windowmasker_path)) {
            spdlog::info("Found windowmasker relative to executable: {}", windowmasker_path.string());
            return windowmasker_path.string();
        }
    }
    catch (const std::exception& e) {
        spdlog::warn("Failed to determine executable path: {}", e.what());
    }

    // 3. 如果都找不到，抛出异常
    throw std::runtime_error("windowmasker executable not found. Please ensure it is installed or available in PATH.");
}

// -----------------------------
// 工具函数：判断字符串是否为 URL
// -----------------------------
bool isUrl(const std::string& path_str) {
    static const std::regex url_pattern(R"(^(https?|ftp)://)");
    return std::regex_search(path_str, url_pattern);
}

// -----------------------------
// 验证 URL 是否可达（通过 CURL 发送 HEAD 请求）
// -----------------------------
bool verifyUrlReachable(const std::string& url) {
    if (url.empty()) {
        throw std::runtime_error("URL is empty and cannot be reached: " + url);
    }

    CURL* curl = curl_easy_init();
    if (!curl) {
        throw std::runtime_error("Failed to initialize CURL for URL verification" +
            url);
    }

    CURLcode res;
    long response_code = 0;

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_NOBODY, 1L);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, 10L);
    curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1L);

    res = curl_easy_perform(curl);
    if (res != CURLE_OK) {
        curl_easy_cleanup(curl);
        throw std::runtime_error("Failed to verify URL: " + url +
            std::string(curl_easy_strerror(res)));
    }

    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &response_code);
    curl_easy_cleanup(curl);

    if (response_code >= 200 && response_code < 400) {
        spdlog::info("Verified URL is reachable: {} (HTTP {})", url, response_code);
        return true;
    }
    else {
        throw std::runtime_error("URL verification failed with HTTP code: " +
            std::to_string(response_code));
    }
}

// -----------------------------
// 验证本地文件是否存在且为常规文件
// -----------------------------
void verifyLocalFile(const FilePath& file_path) {
    if (!std::filesystem::exists(file_path)) {
        throw std::runtime_error("Local file does not exist: " +
            file_path.string());
    }
    if (!std::filesystem::is_regular_file(file_path)) {
        throw std::runtime_error("Path is not a regular file: " +
            file_path.string());
    }
    spdlog::info("Verified local file exists: {}", file_path.string());
}

// ------------------------------
// CURL 写入回调函数，用于下载时写入文件流
// -----------------------------
size_t writeData(void* ptr, size_t size, size_t nmemb, void* stream) {
    std::ofstream* output_stream = static_cast<std::ofstream*>(stream);
    size_t written = size * nmemb;
    output_stream->write(static_cast<const char*>(ptr), written);
    return written;
}

// -----------------------------
// 使用 CURL 下载文件到指定路径
// -----------------------------
void downloadFile(const std::string& url, const FilePath& destination) {
    spdlog::info("Downloading {} to {}", url, destination.string());
    CURL* curl;
    CURLcode result;
    std::ofstream output_stream(destination, std::ios::binary);

    if (!output_stream) {
        throw std::runtime_error("Failed to create file at: " +
            destination.string());
    }

    curl = curl_easy_init();
    if (!curl) {
        throw std::runtime_error("Failed to initialize CURL");
    }

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writeData);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &output_stream);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1L);

    result = curl_easy_perform(curl);
    if (result != CURLE_OK) {
        curl_easy_cleanup(curl);
        throw std::runtime_error("CURL failed: " +
            std::string(curl_easy_strerror(result)));
    }

    curl_easy_cleanup(curl);
    output_stream.close();
    spdlog::info("Download completed: {}", destination.string());
}

// -----------------------------
// 拷贝本地文件到目标路径
// -----------------------------
void copyLocalFile(const FilePath& source, const FilePath& destination) {
    try {
        std::filesystem::copy_file(
            source, destination, std::filesystem::copy_options::overwrite_existing);
    }
    catch (const std::filesystem::filesystem_error& error) {
        throw std::runtime_error("Failed to copy file from " + source.string() +
            " to " + destination.string() + " : " +
            error.what());
    }
    spdlog::info("Copy local file from {} to {}", source.string(),
        destination.string());
}

// -----------------------------
// 获取文件扩展名字符串
// -----------------------------
std::string getFileExtension(const FilePath& file_path) {
    return file_path.extension().string();
}

// -----------------------------
// 主流程函数：复制或下载所有原始数据文件
// -----------------------------
bool copyRawData(const FilePath workdir_path, SpeciesPathMap& species_path_map,
    int thread_num, RaMAxCache::StageStats* cache_stats,
    bool trust_legacy_cache) {
    try {
        // 预验证文件路径或 URL
        for (const auto& [key, path] : species_path_map) {
            if (isUrl(path.string())) {
                verifyUrlReachable(path.string());
            }
            else {
                verifyLocalFile(path);
            }
        }

        // 创建原始数据文件夹
        FilePath raw_data_dir = workdir_path / DATA_DIR / RAW_DATA_DIR;
        std::filesystem::create_directories(raw_data_dir);
        spdlog::info("Created directory: {}", raw_data_dir.string());

        struct RawResult {
            SpeciesName species;
            FilePath path;
            bool reused{false};
        };

        struct RawInput {
            SpeciesName species;
            FilePath path;
        };
        std::vector<RawInput> inputs;
        inputs.reserve(species_path_map.size());
        for (const auto& [key, path] : species_path_map) {
            inputs.push_back({key, path});
        }
        std::vector<RawResult> results(inputs.size());
        parallelForIndexed(inputs.size(), thread_num, [&](size_t index) {
            const auto& key = inputs[index].species;
            const auto& path = inputs[index].path;
            std::string extension = getFileExtension(path);
            std::string final_name = key + extension;
            FilePath final_dest = raw_data_dir / final_name;
            const bool source_is_url = isUrl(path.string());
            const FilePath marker =
                RaMAxCache::completionMarkerPath(final_dest);
            if (RaMAxCache::markerMatches(
                    marker, "raw-fasta", 1, path.string(),
                    source_is_url, path, final_dest)) {
                spdlog::info("[cache] raw FASTA reused for {}: {}",
                             key, final_dest.string());
                results[index] = {key, final_dest, true};
                return;
            }
            if (trust_legacy_cache &&
                std::filesystem::is_regular_file(final_dest) &&
                !std::filesystem::exists(marker)) {
                RaMAxCache::writeMarker(
                    marker, "raw-fasta", 1, path.string(),
                    source_is_url, path, final_dest);
                spdlog::warn(
                    "[cache] trusting legacy raw FASTA for {}: {}",
                    key, final_dest.string());
                results[index] = {key, final_dest, true};
                return;
            }

            FilePath partial = final_dest;
            partial += ".partial";
            RaMAxCache::removeIfPresent(partial);
            try {
                if (source_is_url) {
                    downloadFile(path.string(), partial);
                } else {
                    copyLocalFile(path, partial);
                }
                RaMAxCache::publishFile(partial, final_dest);
                RaMAxCache::writeMarker(
                    marker, "raw-fasta", 1, path.string(),
                    source_is_url, path, final_dest);
            } catch (...) {
                RaMAxCache::removeIfPresent(partial);
                throw;
            }
            spdlog::info("Successfully processed species {} -> {}",
                         key, final_dest.string());
            results[index] = {key, final_dest, false};
        });

        RaMAxCache::StageStats local_stats;
        for (const auto& result : results) {
            species_path_map[result.species] = result.path;
            result.reused ? ++local_stats.reused : ++local_stats.rebuilt;
        }
        if (cache_stats) *cache_stats = local_stats;
        spdlog::info("All raw sequence data copy to work directory successfully.");
        return true;
    }
    catch (const std::exception& error) {
        spdlog::error("Error occurred: {}", error.what());
        throw std::runtime_error("Failed to copy raw data to work directory.");
    }
}

// -----------------------------
// 清洗原始数据：调用 SeqPro 清理并索引
// -----------------------------
bool cleanRawDataset(const FilePath workdir_path,
    SpeciesPathMap& species_path_map, int thread_num,
    RaMAxCache::StageStats* cache_stats, bool trust_legacy_cache) {
    try {
        // 确保清洗数据输出目录存在
        FilePath clean_data_dir = workdir_path / DATA_DIR / CLEAN_DATA_DIR;
        if (!std::filesystem::exists(clean_data_dir)) {
            std::filesystem::create_directories(clean_data_dir);
            spdlog::info("Created directory for cleaned data: {}", clean_data_dir.string());
        }
        else {
            spdlog::info("Cleaned data directory already exists: {}", clean_data_dir.string());
        }

        struct CleanResult {
            SpeciesName species;
            FilePath path;
            bool reused{false};
        };
        struct CleanInput {
            SpeciesName species;
            FilePath path;
        };
        std::vector<CleanInput> inputs;
        inputs.reserve(species_path_map.size());
        for (const auto& [species, raw_path] : species_path_map) {
            inputs.push_back({species, raw_path});
        }
        std::vector<CleanResult> results(inputs.size());
        parallelForIndexed(inputs.size(), thread_num, [&](size_t index) {
            const auto& species = inputs[index].species;
            const auto& raw_path = inputs[index].path;
            FilePath out_dir = workdir_path / DATA_DIR / CLEAN_DATA_DIR;
            FilePath out_fasta = out_dir / (species + ".fasta");
            const FilePath marker =
                RaMAxCache::completionMarkerPath(out_fasta);

            if (RaMAxCache::markerMatches(
                marker, "clean-fasta", 1, raw_path.string(), false,
                raw_path, out_fasta)) {
                spdlog::info("[cache] cleaned FASTA reused for {}: {}",
                             species, out_fasta.string());
                results[index] = {species, out_fasta, true};
                return;
            }
            if (trust_legacy_cache &&
                std::filesystem::is_regular_file(out_fasta) &&
                !std::filesystem::exists(marker)) {
                RaMAxCache::writeMarker(
                    marker, "clean-fasta", 1, raw_path.string(), false,
                    raw_path, out_fasta);
                spdlog::warn(
                    "[cache] trusting legacy cleaned FASTA for {}: {}",
                    species, out_fasta.string());
                results[index] = {species, out_fasta, true};
                return;
            }

            FilePath partial = out_fasta;
            partial += ".partial";
            RaMAxCache::removeIfPresent(partial);
            try {
                SeqPro::utils::cleanFastaFile(raw_path, partial, 60);
                RaMAxCache::publishFile(partial, out_fasta);
                RaMAxCache::writeMarker(
                    marker, "clean-fasta", 1, raw_path.string(), false,
                    raw_path, out_fasta);
            } catch (...) {
                RaMAxCache::removeIfPresent(partial);
                throw;
            }

            spdlog::info("Species {} cleaned file path updated to: {}", species,
                out_fasta.string());
            results[index] = {species, out_fasta, false};
        });

        RaMAxCache::StageStats local_stats;
        for (const auto& result : results) {
            species_path_map[result.species] = result.path;
            result.reused ? ++local_stats.reused : ++local_stats.rebuilt;
        }
        if (cache_stats) *cache_stats = local_stats;
        spdlog::info("All raw sequence data cleaned successfully.");
        return true;
    }
    catch (const std::exception& error) {
        spdlog::error("Error occurred during data cleaning: {}", error.what());
        throw std::runtime_error("Failed to clean raw dataset.");
    }
}

// -----------------------------
// HAL-only cleaning: alignment input remains uppercase, while original
// lowercase runs are stored in a compact sidecar that is not opened until
// HAL export.
// -----------------------------
bool cleanRawDatasetWithSoftMaskIndex(const FilePath workdir_path,
    SpeciesPathMap& species_path_map,
    SoftMask::PathMap& softmask_path_map,
    int thread_num, RaMAxCache::StageStats* cache_stats,
    bool trust_legacy_cache) {
    try {
        const FilePath clean_data_dir = workdir_path / DATA_DIR / CLEAN_DATA_DIR;
        std::filesystem::create_directories(clean_data_dir);
        spdlog::info("HAL soft-mask preprocessing directory: {}", clean_data_dir.string());

        struct CleanResult {
            SpeciesName species;
            FilePath fasta;
            FilePath index;
            bool reused{false};
        };

        struct CleanInput {
            SpeciesName species;
            FilePath path;
        };
        std::vector<CleanInput> inputs;
        inputs.reserve(species_path_map.size());
        for (const auto& [species, raw_path] : species_path_map) {
            inputs.push_back({species, raw_path});
        }
        std::vector<CleanResult> results(inputs.size());
        parallelForIndexed(inputs.size(), thread_num, [&](size_t index) {
            const auto& species = inputs[index].species;
            const auto& raw_path = inputs[index].path;
            const FilePath out_fasta =
                clean_data_dir / (species + ".align-v2.fasta");
            const FilePath out_index =
                clean_data_dir / (species + ".softmask-v1.bin");
            const FilePath marker =
                clean_data_dir / (species + ".softmask-v1.complete.json");

            const bool reused = SoftMask::ensureUppercaseFastaAndIndex(
                raw_path, out_fasta, out_index, marker,
                trust_legacy_cache);
            results[index] = {species, out_fasta, out_index, reused};
        });

        // Publish maps only after each task has completed successfully. This
        // avoids concurrent writes to unordered_map/map and propagates worker
        // exceptions instead of silently swallowing them.
        RaMAxCache::StageStats local_stats;
        for (const auto& result : results) {
            species_path_map[result.species] = result.fasta;
            softmask_path_map[result.species] = result.index;
            result.reused ? ++local_stats.reused : ++local_stats.rebuilt;
            spdlog::info("Species {} alignment FASTA: {}; HAL soft-mask index: {}",
                result.species, result.fasta.string(), result.index.string());
        }
        if (cache_stats) *cache_stats = local_stats;

        if (softmask_path_map.size() != species_path_map.size()) {
            throw std::runtime_error("Not every species received a HAL soft-mask index");
        }
        spdlog::info("All uppercase alignment FASTAs and HAL soft-mask indexes are ready.");
        return true;
    }
    catch (const std::exception& error) {
        spdlog::error("Error during HAL soft-mask preprocessing: {}", error.what());
        throw std::runtime_error("Failed to prepare uppercase FASTA and HAL soft-mask index: " +
            std::string(error.what()));
    }
}

// -----------------------------
// 运行 WindowMasker 对 FASTA 文件进行处理
// -----------------------------
std::map<SpeciesName, FilePath>
repeatSeqMasking(const FilePath workdir_path,
    const SpeciesPathMap& species_path_map_const, int thread_num) {
    SpeciesPathMap species_path_map =
        species_path_map_const; // Create a mutable copy
    std::map<SpeciesName, FilePath> interval_files_map;
    const std::string DATA_DIR_NAME = "data";
    const std::string MASKED_DATA_DIR_NAME = "masked_data";

    try {
        FilePath masked_data_dir =
            workdir_path / DATA_DIR_NAME / MASKED_DATA_DIR_NAME;
        if (!std::filesystem::exists(masked_data_dir)) {
            std::filesystem::create_directories(masked_data_dir);
            spdlog::info("Created directory for windowmasker output: {}",
                masked_data_dir.string());
        }
        else {
            spdlog::info("Windowmasker output directory already exists: {}",
                masked_data_dir.string());
        }

        // 获取 windowmasker 路径
        std::string windowmasker_path;
        try {
            windowmasker_path = findWindowmaskerPath();
        }
        catch (const std::exception& e) {
            spdlog::error("Failed to locate windowmasker: {}", e.what());
            throw std::runtime_error("Cannot proceed with repeat masking: " + std::string(e.what()));
        }

        struct MaskInput {
            SpeciesName species;
            FilePath fasta;
        };
        struct MaskResult {
            SpeciesName species;
            FilePath interval;
            bool success{false};
        };
        std::vector<MaskInput> inputs;
        inputs.reserve(species_path_map.size());
        for (const auto& [species, fasta] : species_path_map) {
            inputs.push_back({species, fasta});
        }
        std::vector<MaskResult> results(inputs.size());
        parallelForIndexed(inputs.size(), thread_num, [&](size_t index) {
                    const auto& species_key = inputs[index].species;
                    const auto& input_fasta_path = inputs[index].fasta;
                    FilePath counts_file = masked_data_dir / (species_key + ".counts");
                    FilePath interval_file = masked_data_dir / (species_key + ".interval");
                    results[index].species = species_key;
                    results[index].interval = interval_file;

                    // 检查最终的interval文件是否已存在
                    if (std::filesystem::exists(interval_file)) {
                        spdlog::warn("[{}] WindowMasker output interval file already exists, "
                            "skipping: {}",
                            species_key, interval_file.string());
                        results[index].success = true;
                        return;
                    }

                    // 命令1: windowmasker -mk_counts
                    std::string cmd1_str = "\"" + windowmasker_path + "\" -mk_counts -in \"" +
                        input_fasta_path.string() + "\" -out \"" +
                        counts_file.string() + "\" -sformat obinary";
                    spdlog::info("[{}] Executing: {}", species_key, cmd1_str);
                    int ret1 = std::system(cmd1_str.c_str());

                    if (ret1 != 0) {
                        spdlog::error("[{}] windowmasker -mk_counts failed with exit code "
                            "{}. Command: {}",
                            species_key, ret1, cmd1_str);
                        if (std::filesystem::exists(counts_file)) {
                            std::filesystem::remove(counts_file);
                        }
                        return; // 该物种任务失败
                    }
                    spdlog::info("[{}] windowmasker -mk_counts successful.", species_key);

                    // 命令2: windowmasker -ustat
                    std::string cmd2_str =
                        "\"" + windowmasker_path + "\" -ustat \"" + counts_file.string() + "\" -in \"" +
                        input_fasta_path.string() + "\" -out \"" + interval_file.string() +
                        "\" -outfmt interval -dust true";
                    spdlog::info("[{}] Executing: {}", species_key, cmd2_str);
                    int ret2 = std::system(cmd2_str.c_str());

                    if (ret2 != 0) {
                        spdlog::error(
                            "[{}] windowmasker -ustat failed with exit code {}. Command: {}",
                            species_key, ret2, cmd2_str);
                        if (std::filesystem::exists(interval_file)) {
                            std::filesystem::remove(interval_file);
                        }
                        return; // 该物种任务失败
                    }
                    spdlog::info("[{}] windowmasker -ustat successful. Output: {}",
                        species_key, interval_file.string());

                    spdlog::info("[{}] Interval file generated at: {}", species_key,
                        interval_file.string());
                    results[index].success = true;
        });

        for (const auto& result : results) {
            if (result.success) {
                interval_files_map[result.species] = result.interval;
            }
        }

        spdlog::info("All windowmasker tasks submitted and processing attempted. "
            "Check logs for individual species status.");
        return interval_files_map; // 返回包含interval文件路径的映射

    }
    catch (const std::exception& error) {
        spdlog::error("Critical error during windowmasker processing setup or "
            "thread pool execution: {}",
            error.what());
        throw std::runtime_error(
            std::string(
                "Failed to run windowmasker processing due to critical error: ") +
            error.what());
    }
    return {}; // 发生严重错误时返回空映射
}

// 辅助函数：解析区间字符串，格式如 "start - end"
std::pair<size_t, size_t> parseIntervalLine(const std::string& line) {
    std::stringstream ss(line);
    size_t start, end;
    char dash;
    ss >> start >> dash >> end;
    if (ss.fail() || dash != '-') {
        throw std::runtime_error("Invalid interval format: " + line);
    }
    return { start, end };
}

// -----------------------------
// 获取人类可读的文件大小字符串（自动转化为 KB / MB / GB）
// -----------------------------
std::string getReadableFileSize(const FilePath& filePath) {
    std::string pathStr = filePath.string();

    // 处理远程 URL 文件
    if (pathStr.rfind("http://", 0) == 0 || pathStr.rfind("https://", 0) == 0) {
        CURL* curl = curl_easy_init();
        if (!curl) {
            spdlog::error("curl_easy_init failed for URL: {}", pathStr);
            return "0 B";
        }
        curl_easy_setopt(curl, CURLOPT_URL, pathStr.c_str());
        curl_easy_setopt(curl, CURLOPT_NOBODY, 1L);
        curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
        curl_easy_setopt(curl, CURLOPT_TIMEOUT, 10L);

        CURLcode res = curl_easy_perform(curl);
        if (res != CURLE_OK) {
            spdlog::error("curl_easy_perform failed for URL {}: {}", pathStr,
                curl_easy_strerror(res));
            curl_easy_cleanup(curl);
            return "0 B";
        }

        curl_off_t content_length = 0;
        res = curl_easy_getinfo(curl, CURLINFO_CONTENT_LENGTH_DOWNLOAD_T,
            &content_length);
        curl_easy_cleanup(curl);

        if (res != CURLE_OK || content_length < 0) {
            spdlog::error("Failed to get Content-Length for URL: {}", pathStr);
            return "0 B";
        }

        // 转换为可读单位
        size_t size = static_cast<size_t>(content_length);
        const char* units[] = { "B", "KB", "MB", "GB", "TB" };
        int unit_index = 0;
        double display_size = static_cast<double>(size);
        while (display_size >= 1024 && unit_index < 4) {
            display_size /= 1024;
            ++unit_index;
        }
        char buf[64];
        snprintf(buf, sizeof(buf), "%.2f %s", display_size, units[unit_index]);
        return std::string(buf);

        // 处理本地文件
    }
    else {
        try {
            auto size = std::filesystem::file_size(filePath);
            const char* units[] = { "B", "KB", "MB", "GB", "TB" };
            int unit_index = 0;
            double display_size = static_cast<double>(size);
            while (display_size >= 1024 && unit_index < 4) {
                display_size /= 1024;
                ++unit_index;
            }
            char buf[64];
            snprintf(buf, sizeof(buf), "%.2f %s", display_size, units[unit_index]);
            return std::string(buf);
        }
        catch (const std::filesystem::filesystem_error& e) {
            spdlog::error("Failed to get file size for {}: {}", filePath.string(),
                e.what());
            return "0 B";
        }
    }
}

// -----------------------------
// 判断文件是否小于指定大小（单位：MB）
// -----------------------------
bool isFileSmallerThan(const FilePath& filePath, size_t maxSizeMB) {
    const std::string pathStr = filePath.string();
    const curl_off_t maxBytes = static_cast<curl_off_t>(maxSizeMB) * 1024 * 1024;

    // 对 URL 文件使用 CURL 获取 Content-Length
    if (pathStr.rfind("http://", 0) == 0 || pathStr.rfind("https://", 0) == 0) {
        CURL* curl = curl_easy_init();
        if (!curl) {
            spdlog::error("curl_easy_init failed for URL: {}", pathStr);
            return false;
        }
        curl_easy_setopt(curl, CURLOPT_URL, pathStr.c_str());
        curl_easy_setopt(curl, CURLOPT_NOBODY, 1L);
        curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
        curl_easy_setopt(curl, CURLOPT_TIMEOUT, 10L);

        if (curl_easy_perform(curl) != CURLE_OK) {
            spdlog::error("curl_easy_perform failed for URL: {}", pathStr);
            curl_easy_cleanup(curl);
            return false;
        }

        curl_off_t content_length = 0;
        if (curl_easy_getinfo(curl, CURLINFO_CONTENT_LENGTH_DOWNLOAD_T,
            &content_length) != CURLE_OK ||
            content_length < 0) {
            spdlog::error("Failed to get Content-Length for URL: {}", pathStr);
            curl_easy_cleanup(curl);
            return false;
        }
        curl_easy_cleanup(curl);

        return content_length < maxBytes;

        // 本地文件直接获取大小
    }
    else {
        try {
            auto size = static_cast<curl_off_t>(std::filesystem::file_size(filePath));
            return size < maxBytes;
        }
        catch (const std::filesystem::filesystem_error& e) {
            spdlog::error("Failed to get file size for {}: {}", pathStr, e.what());
            return false;
        }
    }
}

// -----------------------------
// 生成临时文件路径（用于处理中间文件）
// 例如 a.fasta -> a_in_process.fasta
// -----------------------------
FilePath getTempFilePath(const FilePath& input_path) {
    FilePath temp_file;
    if (input_path.has_extension()) {
        std::string stem = input_path.stem().string();
        std::string extension = input_path.extension().string();
        temp_file = input_path.parent_path() / (stem + "_in_process" + extension);
    }
    else {
        temp_file = input_path.parent_path() /
            (input_path.filename().string() + "_in_process");
    }
    return temp_file;
}

// -----------------------------
// 获取对应的 .fai 索引文件路径
// 例如 test.fasta -> test.fasta.fai
// -----------------------------
FilePath getFaiIndexPath(const FilePath& fasta_path) {
    FilePath fai_path = fasta_path;
    fai_path += ".fai";
    return fai_path;
}

// -----------------------------
// parseSeqfile：解析可选 Newick 树和物种路径映射
// -----------------------------
bool parseSeqfile(const FilePath& seqfile_path,
    NewickParser& newick_tree,
    SpeciesPathMap& species_map,
    const std::string& root /* = "" */) {
    if (!std::filesystem::exists(seqfile_path)) {
        throw std::runtime_error("parseSeqfile: cannot locate file: " +
            seqfile_path.string());
    }

    std::ifstream ifs(seqfile_path);
    if (!ifs.is_open()) {
        throw std::runtime_error("parseSeqfile: failed to open: " +
            seqfile_path.string());
    }

    std::string line;
    bool first_record = true;
    bool got_tree = false;

    species_map.clear();
    newick_tree.clear();

    // 用于子树过滤的叶子名集合
    std::unordered_set<std::string> allowed_leaves;

    while (std::getline(ifs, line)) {
        auto l = line.find_first_not_of(" \t\r\n");
        if (l == std::string::npos) continue;
        auto r = line.find_last_not_of(" \t\r\n");
        std::string trimmed = line.substr(l, r - l + 1);

        if (first_record) {
            first_record = false;

            // Multi-leaf trees start with '('; standard trees end with ';';
            // a single-field first record preserves single-leaf Newick input.
            const bool has_whitespace = std::any_of(
                trimmed.begin(), trimmed.end(), [](unsigned char c) {
                    return std::isspace(c) != 0;
                });
            const bool looks_like_newick =
                trimmed.front() == '(' || trimmed.back() == ';' || !has_whitespace;

            if (looks_like_newick) {
                try {
                    newick_tree.parse(trimmed);
                }
                catch (const std::runtime_error& e) {
                    throw std::runtime_error(
                        std::string("parseSeqfile: invalid Newick tree: ") + e.what());
                }
                got_tree = true;

                // 如果指定了 root，则裁剪到该子树
                if (!root.empty()) {
                    int rootId = newick_tree.findNodeIdByName(root);
                    if (rootId == -1) {
                        throw std::runtime_error("parseSeqfile: root name `" + root + "` not found in Newick tree");
                    }
                    newick_tree.restrictToSubtreeByRootId(rootId);
                }

                // 构建允许的叶子名集合（用于过滤 species）
                auto leafNames = newick_tree.getLeafNames();
                for (auto& nm : leafNames) {
                    std::string tmp = nm;
                    newick_tree.trimString(tmp);
                    if (!tmp.empty()) allowed_leaves.insert(std::move(tmp));
                }

                continue;
            }
        }

        // 物种映射：<speciesName> <path>
        std::istringstream iss(trimmed);
        std::string speciesName;
        std::string filePath;
        if (!(iss >> speciesName)) {
            continue;
        }
        if (!(iss >> std::ws) || !std::getline(iss, filePath)) {
            throw std::runtime_error("parseSeqfile: missing file path for species `" +
                speciesName + "`");
        }
        auto p_l = filePath.find_first_not_of(" \t\r\n");
        auto p_r = filePath.find_last_not_of(" \t\r\n");
        if (p_l == std::string::npos) {
            throw std::runtime_error("parseSeqfile: empty file path for species `" +
                speciesName + "`");
        }
        filePath = filePath.substr(p_l, p_r - p_l + 1);

        // 如果指定了 root，只接收子树中的叶子；否则全部接收
        if (!root.empty() && allowed_leaves.find(speciesName) == allowed_leaves.end()) {
            continue;
        }

        species_map.emplace(speciesName, FilePath(filePath));
    }

    if (first_record) {
        throw std::runtime_error("parseSeqfile: no non-empty records found");
    }
    if (species_map.empty()) {
        throw std::runtime_error("parseSeqfile: no species=>path mappings found"
            + std::string(root.empty() ? "" : " under root `" + root + "`"));
    }

    return got_tree;
}
