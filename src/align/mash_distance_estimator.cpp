#include "mash_distance_estimator.h"

#include <spdlog/spdlog.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <limits.h>
#include <spawn.h>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <sys/wait.h>
#include <unistd.h>
#include <unordered_map>

extern char** environ;

#ifndef RAMAX_MASH_CONFIGURED_PATH
#define RAMAX_MASH_CONFIGURED_PATH ""
#endif
#ifndef RAMAX_TOOL_BIN_CONFIGURED_PATH
#define RAMAX_TOOL_BIN_CONFIGURED_PATH ""
#endif

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

std::filesystem::path searchPath(std::string_view name) {
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
                name;
            if (isExecutable(candidate)) {
                return std::filesystem::absolute(candidate);
            }
        }
        if (separator == std::string_view::npos) break;
        offset = separator + 1;
    }
    return {};
}

std::string readText(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Cannot read Mash output: " + path.string());
    }
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

void runCommand(const std::filesystem::path& executable,
                const std::vector<std::string>& arguments,
                const std::filesystem::path& stdout_path,
                const std::filesystem::path& stderr_path) {
    const int stdout_fd = ::open(stdout_path.c_str(),
        O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (stdout_fd < 0) {
        throw std::runtime_error("Cannot open Mash stdout file: " +
                                 stdout_path.string());
    }
    const int stderr_fd = ::open(stderr_path.c_str(),
        O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (stderr_fd < 0) {
        ::close(stdout_fd);
        throw std::runtime_error("Cannot open Mash stderr file: " +
                                 stderr_path.string());
    }

    posix_spawn_file_actions_t actions;
    posix_spawn_file_actions_init(&actions);
    posix_spawn_file_actions_adddup2(&actions, stdout_fd, STDOUT_FILENO);
    posix_spawn_file_actions_adddup2(&actions, stderr_fd, STDERR_FILENO);
    posix_spawn_file_actions_addclose(&actions, stdout_fd);
    posix_spawn_file_actions_addclose(&actions, stderr_fd);

    std::vector<std::string> storage;
    storage.reserve(arguments.size() + 1);
    storage.push_back(executable.string());
    storage.insert(storage.end(), arguments.begin(), arguments.end());
    std::vector<char*> argv;
    argv.reserve(storage.size() + 1);
    for (auto& item : storage) argv.push_back(item.data());
    argv.push_back(nullptr);

    pid_t pid = 0;
    const int spawn_error = ::posix_spawn(
        &pid, executable.c_str(), &actions, nullptr, argv.data(), environ);
    posix_spawn_file_actions_destroy(&actions);
    ::close(stdout_fd);
    ::close(stderr_fd);
    if (spawn_error != 0) {
        throw std::runtime_error("Cannot start Mash: " +
                                 std::string(std::strerror(spawn_error)));
    }

    int status = 0;
    while (::waitpid(pid, &status, 0) < 0) {
        if (errno == EINTR) continue;
        throw std::runtime_error("Cannot wait for Mash process");
    }
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        const std::string detail = readText(stderr_path);
        throw std::runtime_error(
            "Mash exited unsuccessfully" +
            (detail.empty() ? std::string() : ": " + detail));
    }
}

std::string trim(std::string_view value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string_view::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return std::string(value.substr(first, last - first + 1));
}

double parseFiniteUnit(std::string_view text, const char* field) {
    const std::string value(text);
    size_t consumed = 0;
    double parsed = 0.0;
    try {
        parsed = std::stod(value, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error(std::string("Invalid Mash ") + field +
                                 ": " + value);
    }
    if (consumed != value.size() || !std::isfinite(parsed) ||
        parsed < 0.0 || parsed > 1.0) {
        throw std::runtime_error(std::string("Invalid Mash ") + field +
                                 ": " + value);
    }
    return parsed;
}

uint64_t parseUnsigned(std::string_view text, const char* field) {
    const std::string value(text);
    size_t consumed = 0;
    uint64_t parsed = 0;
    try {
        parsed = std::stoull(value, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error(std::string("Invalid Mash ") + field +
                                 ": " + value);
    }
    if (consumed != value.size()) {
        throw std::runtime_error(std::string("Invalid Mash ") + field +
                                 ": " + value);
    }
    return parsed;
}

FilePath fastaPath(const SeqPro::SharedManagerVariant& manager) {
    if (!manager) throw std::runtime_error("Null sequence manager");
    return std::visit([](const auto& pointer) -> FilePath {
        if (!pointer) throw std::runtime_error("Null sequence manager pointer");
        using Pointer = std::decay_t<decltype(pointer)>;
        if constexpr (std::is_same_v<Pointer,
                      std::unique_ptr<SeqPro::SequenceManager>>) {
            return pointer->getFastaPath();
        } else {
            return pointer->getOriginalManager().getFastaPath();
        }
    }, *manager);
}

std::string normalizedPath(const FilePath& path) {
    std::error_code error;
    auto absolute = std::filesystem::absolute(path, error);
    if (error) absolute = path;
    return absolute.lexically_normal().string();
}

}  // namespace

namespace MashDistanceDetail {

ParsedLine parseLine(std::string_view line) {
    std::array<std::string_view, 5> fields{};
    size_t start = 0;
    for (size_t i = 0; i < fields.size(); ++i) {
        const size_t end = line.find('\t', start);
        if (i + 1 == fields.size()) {
            if (end != std::string_view::npos) {
                throw std::runtime_error("Mash output has more than five fields");
            }
            fields[i] = line.substr(start);
        } else {
            if (end == std::string_view::npos) {
                throw std::runtime_error("Mash output has fewer than five fields");
            }
            fields[i] = line.substr(start, end - start);
            start = end + 1;
        }
    }
    ParsedLine result;
    result.reference_id = std::string(fields[0]);
    result.query_id = std::string(fields[1]);
    if (result.reference_id.empty() || result.query_id.empty()) {
        throw std::runtime_error("Mash output contains an empty sequence ID");
    }
    result.distance = parseFiniteUnit(fields[2], "distance");
    result.p_value = parseFiniteUnit(fields[3], "p-value");
    const size_t slash = fields[4].find('/');
    if (slash == std::string_view::npos ||
        fields[4].find('/', slash + 1) != std::string_view::npos) {
        throw std::runtime_error("Invalid Mash shared-hashes field");
    }
    result.shared_hashes =
        parseUnsigned(fields[4].substr(0, slash), "shared hash count");
    result.total_hashes =
        parseUnsigned(fields[4].substr(slash + 1), "total hash count");
    if (result.total_hashes == 0 ||
        result.shared_hashes > result.total_hashes) {
        throw std::runtime_error("Invalid Mash shared-hashes range");
    }
    return result;
}

std::string validateVersion(std::string_view output) {
    const std::string value = trim(output);
    if (value != "2.3" && value != "Mash version 2.3") {
        throw std::runtime_error(
            "RaMAx requires Mash 2.3, found: " + value);
    }
    return "2.3";
}

}  // namespace MashDistanceDetail

std::filesystem::path locateMashExecutable() {
    const std::filesystem::path configured(RAMAX_MASH_CONFIGURED_PATH);
    if (isExecutable(configured)) return configured;
    const std::filesystem::path configured_directory(
        RAMAX_TOOL_BIN_CONFIGURED_PATH);
    if (!configured_directory.empty()) {
        const auto candidate = configured_directory / "mash";
        if (isExecutable(candidate)) return candidate;
    }
    const auto sibling = executableDirectory() / "mash";
    if (isExecutable(sibling)) return sibling;
    return searchPath("mash");
}

MashDistanceEstimator::MashDistanceEstimator(
    std::filesystem::path executable,
    std::filesystem::path output_directory,
    uint_t threads)
    : executable_(std::move(executable)),
      output_directory_(std::move(output_directory)),
      threads_(std::max<uint_t>(1, threads)) {
    if (!isExecutable(executable_)) {
        throw std::runtime_error(
            "Mash 2.3 is required but was not found next to ramax, in the "
            "configured path, or in PATH");
    }
    std::filesystem::create_directories(output_directory_);
    const auto version_stdout = output_directory_ / "mash.version.stdout";
    const auto version_stderr = output_directory_ / "mash.version.stderr";
    runCommand(executable_, {"--version"}, version_stdout, version_stderr);
    version_ = MashDistanceDetail::validateVersion(readText(version_stdout));
}

std::vector<MashDistanceRecord> MashDistanceEstimator::estimateFirstReference(
    const SpeciesName& reference,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers) {
    const auto ref_it = managers.find(reference);
    if (ref_it == managers.end()) {
        throw std::runtime_error("Mash reference is absent: " + reference);
    }
    if (managers.size() < 2) {
        throw std::runtime_error(
            "Mash first-reference estimation requires at least two genomes");
    }

    const FilePath reference_fasta = fastaPath(ref_it->second);
    const std::string reference_path = normalizedPath(reference_fasta);
    std::map<std::string, std::vector<SpeciesName>> species_by_path;
    for (const auto& [species, manager] : managers) {
        if (species == reference) continue;
        const std::string path = normalizedPath(fastaPath(manager));
        if (path.find('\n') != std::string::npos ||
            path.find('\r') != std::string::npos) {
            throw std::runtime_error("FASTA path contains a newline: " + path);
        }
        species_by_path[path].push_back(species);
    }

    const auto query_list = output_directory_ / "mash_query_paths.txt";
    {
        std::ofstream output(query_list, std::ios::trunc);
        if (!output) {
            throw std::runtime_error("Cannot write Mash query list: " +
                                     query_list.string());
        }
        for (const auto& [path, unused] : species_by_path) output << path << '\n';
        output.flush();
        if (!output) {
            throw std::runtime_error("Cannot finalize Mash query list: " +
                                     query_list.string());
        }
    }

    const auto raw_stdout = output_directory_ / "mash.dist.stdout";
    const auto raw_stderr = output_directory_ / "mash.stderr.log";
    runCommand(executable_, {
        "dist", "-p", std::to_string(threads_),
        "-k", std::to_string(kKmerSize),
        "-s", std::to_string(kSketchSize),
        "-l", reference_path, query_list.string()
    }, raw_stdout, raw_stderr);

    std::unordered_map<std::string, MashDistanceDetail::ParsedLine> parsed_by_path;
    std::ifstream input(raw_stdout);
    if (!input) {
        throw std::runtime_error("Cannot read Mash distance output: " +
                                 raw_stdout.string());
    }
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        auto parsed = MashDistanceDetail::parseLine(line);
        if (normalizedPath(parsed.reference_id) != reference_path) {
            throw std::runtime_error(
                "Mash returned an unexpected reference: " +
                parsed.reference_id);
        }
        const std::string query_path = normalizedPath(parsed.query_id);
        if (!species_by_path.contains(query_path)) {
            throw std::runtime_error(
                "Mash returned an unexpected query: " + parsed.query_id);
        }
        if (!parsed_by_path.emplace(query_path, std::move(parsed)).second) {
            throw std::runtime_error(
                "Mash returned a duplicate query: " + query_path);
        }
    }
    if (parsed_by_path.size() != species_by_path.size()) {
        throw std::runtime_error(
            "Mash result count does not match the requested query count");
    }

    std::vector<MashDistanceRecord> records;
    for (const auto& [path, species_names] : species_by_path) {
        const auto& parsed = parsed_by_path.at(path);
        for (const auto& species : species_names) {
            records.push_back({
                reference, species, parsed.distance, parsed.p_value,
                parsed.shared_hashes, parsed.total_hashes,
                reference_fasta, FilePath(path)});
        }
    }
    std::sort(records.begin(), records.end(),
        [](const auto& left, const auto& right) {
            return left.query < right.query;
        });

    const auto final_path = output_directory_ / "mash_first_reference.tsv";
    const auto part_path = output_directory_ / "mash_first_reference.tsv.part";
    {
        std::ofstream output(part_path, std::ios::trunc);
        if (!output) {
            throw std::runtime_error("Cannot write Mash distance table");
        }
        output << "reference\tquery\tdistance\tp_value\tshared_hashes\t"
                  "kmer_size\tsketch_size\treference_fasta\tquery_fasta\t"
                  "mash_version\n";
        output << std::setprecision(17);
        for (const auto& record : records) {
            output << record.reference << '\t' << record.query << '\t'
                   << record.distance << '\t' << record.p_value << '\t'
                   << record.shared_hashes << '/' << record.total_hashes << '\t'
                   << kKmerSize << '\t' << kSketchSize << '\t'
                   << normalizedPath(record.reference_fasta) << '\t'
                   << normalizedPath(record.query_fasta) << '\t'
                   << version_ << '\n';
        }
        output.flush();
        if (!output) throw std::runtime_error("Cannot finalize Mash distance table");
    }
    std::error_code rename_error;
    std::filesystem::rename(part_path, final_path, rename_error);
    if (rename_error) {
        throw std::runtime_error("Cannot publish Mash distance table: " +
                                 rename_error.message());
    }

    const auto [minimum, maximum] = std::minmax_element(
        records.begin(), records.end(), [](const auto& left, const auto& right) {
            return left.distance < right.distance;
        });
    spdlog::info(
        "[mash-distance] executable={} version={} reference={} pairs={} "
        "k={} sketch_size={} distance_min={} distance_max={} table={}",
        executable_.string(), version_, reference, records.size(), kKmerSize,
        kSketchSize, minimum->distance, maximum->distance, final_path.string());
    return records;
}
