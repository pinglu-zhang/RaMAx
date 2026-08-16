#include "ramax_paf_fasta.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>

#include <curl/curl.h>
#include <zlib.h>

namespace RamaxPafFasta {
namespace {

namespace fs = std::filesystem;

std::string trim(const std::string& value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::string lowercase(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char character) {
                       return static_cast<char>(std::tolower(character));
                   });
    return value;
}

bool endsWith(const std::string& value, const std::string& suffix) {
    return value.size() >= suffix.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool isSupportedOutput(const fs::path& path, bool& gzip) {
    const std::string name = lowercase(path.filename().string());
    gzip = endsWith(name, ".gz");
    const std::string base = gzip ? name.substr(0, name.size() - 3) : name;
    return endsWith(base, ".fa") || endsWith(base, ".fasta") ||
           endsWith(base, ".fna");
}

bool isHttpUrl(const std::string& source) {
    const std::string normalized = lowercase(source);
    return normalized.rfind("http://", 0) == 0 ||
           normalized.rfind("https://", 0) == 0;
}

fs::path uniqueTemporaryPath(const fs::path& parent,
                             const std::string& stem) {
    const auto stamp = std::chrono::steady_clock::now()
                           .time_since_epoch().count();
    for (std::uint64_t attempt = 0; ; ++attempt) {
        const fs::path candidate = parent /
            ("." + stem + ".tmp." + std::to_string(stamp) + "." +
             std::to_string(attempt));
        if (!fs::exists(candidate)) return candidate;
    }
}

struct PathGuard {
    fs::path path;
    bool keep = false;
    ~PathGuard() {
        if (!keep && !path.empty()) {
            std::error_code error;
            fs::remove(path, error);
        }
    }
};

class OutputSink {
public:
    OutputSink(const fs::path& path, bool gzip) : gzip_(gzip) {
        if (gzip_) {
            gzip_file_ = gzopen(path.c_str(), "wb6");
            if (!gzip_file_) {
                throw std::runtime_error("cannot open gzip output: " +
                                         path.string());
            }
        } else {
            plain_.open(path, std::ios::binary | std::ios::trunc);
            if (!plain_) {
                throw std::runtime_error("cannot open output: " + path.string());
            }
        }
        buffer_.reserve(kBufferSize);
    }

    OutputSink(const OutputSink&) = delete;
    OutputSink& operator=(const OutputSink&) = delete;

    ~OutputSink() {
        if (!closed_) {
            try {
                close();
            } catch (...) {
            }
        }
    }

    void beginRecord(const std::string& name) {
        append(">");
        append(name);
        append("\n");
        column_ = 0;
    }

    void appendBase(char base) {
        buffer_.push_back(base);
        ++column_;
        if (column_ == 60) {
            buffer_.push_back('\n');
            column_ = 0;
        }
        flushIfNeeded();
    }

    void endRecord() {
        if (column_ != 0) {
            buffer_.push_back('\n');
            column_ = 0;
        }
        flushIfNeeded();
    }

    void close() {
        if (closed_) return;
        flush();
        if (gzip_) {
            const int status = gzclose(gzip_file_);
            gzip_file_ = nullptr;
            if (status != Z_OK) {
                throw std::runtime_error("failed to finalize gzip output");
            }
        } else {
            plain_.close();
            if (!plain_) throw std::runtime_error("failed to finalize output");
        }
        closed_ = true;
    }

private:
    static constexpr std::size_t kBufferSize = 1U << 20;

    void append(std::string_view value) {
        buffer_.append(value.data(), value.size());
        flushIfNeeded();
    }

    void flushIfNeeded() {
        if (buffer_.size() >= kBufferSize) flush();
    }

    void flush() {
        if (buffer_.empty()) return;
        if (gzip_) {
            std::size_t written = 0;
            while (written < buffer_.size()) {
                const unsigned int chunk = static_cast<unsigned int>(
                    std::min<std::size_t>(buffer_.size() - written,
                                          1U << 30));
                const int result = gzwrite(
                    gzip_file_, buffer_.data() + written, chunk);
                if (result <= 0 || static_cast<unsigned int>(result) != chunk) {
                    int error_number = Z_OK;
                    const char* message = gzerror(gzip_file_, &error_number);
                    throw std::runtime_error(
                        "failed to write gzip output: " +
                        std::string(message ? message : "unknown error"));
                }
                written += chunk;
            }
        } else {
            plain_.write(buffer_.data(),
                         static_cast<std::streamsize>(buffer_.size()));
            if (!plain_) throw std::runtime_error("failed to write output");
        }
        buffer_.clear();
    }

    bool gzip_ = false;
    bool closed_ = false;
    gzFile gzip_file_ = nullptr;
    std::ofstream plain_;
    std::string buffer_;
    std::size_t column_ = 0;
};

size_t curlWrite(void* data, size_t size, size_t count, void* context) {
    auto* output = static_cast<std::ofstream*>(context);
    const size_t bytes = size * count;
    output->write(static_cast<const char*>(data),
                  static_cast<std::streamsize>(bytes));
    return *output ? bytes : 0;
}

void download(const std::string& url, const fs::path& destination) {
    std::ofstream output(destination, std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error("cannot create URL temporary file: " +
                                 destination.string());
    }
    std::unique_ptr<CURL, decltype(&curl_easy_cleanup)> curl(
        curl_easy_init(), &curl_easy_cleanup);
    if (!curl) throw std::runtime_error("cannot initialize CURL");
    curl_easy_setopt(curl.get(), CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl.get(), CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl.get(), CURLOPT_FAILONERROR, 1L);
    curl_easy_setopt(curl.get(), CURLOPT_WRITEFUNCTION, curlWrite);
    curl_easy_setopt(curl.get(), CURLOPT_WRITEDATA, &output);
    const CURLcode result = curl_easy_perform(curl.get());
    output.close();
    if (result != CURLE_OK || !output) {
        throw std::runtime_error(
            "failed to download " + url + ": " +
            std::string(curl_easy_strerror(result)));
    }
}

char normalizeBase(char value, std::uint64_t& converted_to_n) {
    const unsigned char byte = static_cast<unsigned char>(value);
    const char upper = static_cast<char>(std::toupper(byte));
    if (upper == 'A' || upper == 'T' || upper == 'G' || upper == 'C') {
        return upper;
    }
    if (upper != 'N') ++converted_to_n;
    return 'N';
}

std::string headerId(const std::string& header) {
    const std::string cleaned = trim(header);
    if (cleaned.empty()) throw std::runtime_error("FASTA header is empty");
    const auto separator = cleaned.find_first_of(" \t");
    return cleaned.substr(0, separator);
}

void processFasta(const fs::path& input,
                  const std::string& species,
                  OutputSink& output,
                  std::set<std::string>& qualified_names,
                  Stats& stats) {
    gzFile file = gzopen(input.c_str(), "rb");
    if (!file) throw std::runtime_error("cannot open FASTA: " + input.string());
    struct GzGuard {
        gzFile file;
        ~GzGuard() { if (file) gzclose(file); }
    } guard{file};

    constexpr std::size_t kReadSize = 1U << 20;
    constexpr std::size_t kMaximumHeader = 1U << 20;
    std::array<char, kReadSize> buffer{};
    std::string header;
    bool line_start = true;
    bool reading_header = false;
    bool have_record = false;
    std::uint64_t record_bases = 0;
    std::uint64_t input_records = 0;

    auto finishRecord = [&]() {
        if (!have_record) return;
        if (record_bases == 0) {
            throw std::runtime_error(
                "FASTA record has no sequence bases: " + input.string());
        }
        output.endRecord();
        record_bases = 0;
    };

    auto finishHeader = [&]() {
        finishRecord();
        const std::string contig = headerId(header);
        const std::string qualified = species + "." + contig;
        if (std::any_of(qualified.begin(), qualified.end(),
                        [](unsigned char character) {
                            return std::isspace(character) != 0;
                        })) {
            throw std::runtime_error(
                "qualified FASTA name contains whitespace: " + qualified);
        }
        if (!qualified_names.insert(qualified).second) {
            throw std::runtime_error(
                "duplicate qualified FASTA name: " + qualified);
        }
        output.beginRecord(qualified);
        have_record = true;
        ++input_records;
        ++stats.contigs;
        header.clear();
    };

    for (;;) {
        const int count = gzread(file, buffer.data(),
                                 static_cast<unsigned int>(buffer.size()));
        if (count < 0) {
            int error_number = Z_OK;
            const char* message = gzerror(file, &error_number);
            throw std::runtime_error(
                "failed to read FASTA " + input.string() + ": " +
                std::string(message ? message : "unknown gzip error"));
        }
        if (count == 0) break;

        for (int index = 0; index < count; ++index) {
            const char character = buffer[static_cast<std::size_t>(index)];
            if (reading_header) {
                if (character == '\n') {
                    if (!header.empty() && header.back() == '\r') {
                        header.pop_back();
                    }
                    finishHeader();
                    reading_header = false;
                    line_start = true;
                } else {
                    if (header.size() == kMaximumHeader) {
                        throw std::runtime_error(
                            "FASTA header exceeds 1 MiB: " + input.string());
                    }
                    header.push_back(character);
                }
                continue;
            }

            if (line_start && character == '>') {
                reading_header = true;
                line_start = false;
                header.clear();
                continue;
            }
            if (character == '\n') {
                line_start = true;
                continue;
            }
            if (character == '\r') continue;
            if (std::isspace(static_cast<unsigned char>(character)) != 0) {
                continue;
            }
            if (!have_record) {
                throw std::runtime_error(
                    "sequence data appears before the first FASTA header: " +
                    input.string());
            }
            line_start = false;
            output.appendBase(normalizeBase(character, stats.converted_to_n));
            ++record_bases;
            ++stats.bases;
        }
    }

    if (reading_header) finishHeader();
    finishRecord();
    if (input_records == 0) {
        throw std::runtime_error("FASTA contains no records: " + input.string());
    }
}

fs::path normalizedAbsolute(const fs::path& path) {
    std::error_code error;
    fs::path result = fs::weakly_canonical(path, error);
    if (error) result = fs::absolute(path).lexically_normal();
    return result;
}

void publish(const fs::path& temporary, const fs::path& output, bool force) {
    if (!fs::exists(output)) {
        fs::rename(temporary, output);
        return;
    }
    if (!force) {
        throw std::runtime_error("output already exists; use --force: " +
                                 output.string());
    }
    if (!fs::is_regular_file(output)) {
        throw std::runtime_error("existing output is not a regular file: " +
                                 output.string());
    }

    const fs::path backup = uniqueTemporaryPath(
        output.parent_path(), output.filename().string() + ".backup");
    fs::rename(output, backup);
    try {
        fs::rename(temporary, output);
    } catch (...) {
        std::error_code rollback_error;
        fs::rename(backup, output, rollback_error);
        throw;
    }
    std::error_code remove_error;
    fs::remove(backup, remove_error);
    if (remove_error) {
        throw std::runtime_error(
            "output published but backup cleanup failed: " + backup.string());
    }
}

}  // namespace

std::vector<SeqfileEntry> parseSeqfile(const fs::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open seqfile: " + path.string());

    std::vector<SeqfileEntry> entries;
    std::set<std::string> species_seen;
    bool first_record = true;
    std::string line;
    while (std::getline(input, line)) {
        const std::string record = trim(line);
        if (record.empty()) continue;
        if (first_record) {
            first_record = false;
            const bool has_whitespace = std::any_of(
                record.begin(), record.end(), [](unsigned char character) {
                    return std::isspace(character) != 0;
                });
            if (record.front() == '(' || record.back() == ';' ||
                !has_whitespace) {
                continue;
            }
        }

        const auto separator = std::find_if(
            record.begin(), record.end(), [](unsigned char character) {
                return std::isspace(character) != 0;
            });
        if (separator == record.end()) {
            throw std::runtime_error(
                "seqfile mapping is missing a FASTA path: " + record);
        }
        const std::string species(record.begin(), separator);
        const std::string source = trim(
            std::string(separator, record.end()));
        if (species.empty() || source.empty()) {
            throw std::runtime_error("invalid seqfile mapping: " + record);
        }
        if (!species_seen.insert(species).second) {
            throw std::runtime_error(
                "duplicate species in seqfile: " + species);
        }
        entries.push_back({species, source});
    }
    if (first_record) throw std::runtime_error("seqfile contains no records");
    if (entries.empty()) throw std::runtime_error("seqfile contains no FASTA mappings");
    return entries;
}

Stats generate(const Options& options) {
    const auto started = std::chrono::steady_clock::now();
    bool gzip_output = false;
    if (!isSupportedOutput(options.output, gzip_output)) {
        throw std::runtime_error(
            "output must end in .fa, .fasta, or .fna, optionally followed by .gz");
    }
    if (!fs::exists(options.seqfile) || !fs::is_regular_file(options.seqfile)) {
        throw std::runtime_error("seqfile does not exist: " +
                                 options.seqfile.string());
    }
    if (fs::exists(options.output) && !options.force) {
        throw std::runtime_error("output already exists; use --force: " +
                                 options.output.string());
    }

    const auto entries = parseSeqfile(options.seqfile);
    const fs::path output = fs::absolute(options.output).lexically_normal();
    const fs::path parent = output.parent_path();
    fs::create_directories(parent);
    for (const auto& entry : entries) {
        if (!isHttpUrl(entry.source) &&
            normalizedAbsolute(entry.source) == normalizedAbsolute(output)) {
            throw std::runtime_error(
                "output path is also an input FASTA: " + output.string());
        }
    }

    const fs::path temporary = uniqueTemporaryPath(
        parent, output.filename().string());
    PathGuard output_guard{temporary};
    OutputSink sink(temporary, gzip_output);
    Stats stats;
    stats.species = entries.size();
    stats.gzip_output = gzip_output;
    std::set<std::string> qualified_names;

    if (curl_global_init(CURL_GLOBAL_DEFAULT) != CURLE_OK) {
        throw std::runtime_error("cannot initialize CURL globally");
    }
    struct CurlGlobalGuard {
        ~CurlGlobalGuard() { curl_global_cleanup(); }
    } curl_guard;

    for (const auto& entry : entries) {
        fs::path source = entry.source;
        std::unique_ptr<PathGuard> download_guard;
        if (isHttpUrl(entry.source)) {
            source = uniqueTemporaryPath(
                parent, "ramax-paf-fasta-download");
            download_guard = std::make_unique<PathGuard>();
            download_guard->path = source;
            download(entry.source, source);
        } else if (!fs::exists(source) || !fs::is_regular_file(source)) {
            throw std::runtime_error(
                "input FASTA does not exist: " + source.string());
        }
        processFasta(source, entry.species, sink, qualified_names, stats);
    }

    sink.close();
    publish(temporary, output, options.force);
    output_guard.keep = true;
    stats.output_bytes = fs::file_size(output);
    stats.elapsed_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    return stats;
}

}  // namespace RamaxPafFasta
