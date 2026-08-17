#include "wfmash_router.h"

#include "external_tool.h"

#include <spdlog/spdlog.h>

#include <algorithm>
#include <atomic>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mutex>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <system_error>
#include <tuple>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef RAMAX_SAMTOOLS_CONFIGURED_PATH
#define RAMAX_SAMTOOLS_CONFIGURED_PATH ""
#endif

#ifndef RAMAX_WFMASH_CONFIGURED_PATH
#define RAMAX_WFMASH_CONFIGURED_PATH ""
#endif

namespace {

using WfmashRouterDetail::ParsedPafRecord;
using WfmashRouterDetail::SequenceRecord;

constexpr std::string_view kWfmashParameterSummary =
    "-w 5000 -l 25000 -p 90 -n 1 -k 19 -F 0.001 "
    "--hg-filter=1.0,30,99.9";

// Some wfmash builds intermittently fail to reopen a freshly produced mapping
// PAF when several alignment processes start together. Keep the normal
// OpenMP-parallel fast path, but serialize one retry for failed alignments.
std::mutex wfmash_alignment_retry_mutex;

uint64_t parseUnsigned(std::string_view value, const char* field) {
    if (value.empty()) {
        throw std::runtime_error(std::string("Empty ") + field);
    }
    uint64_t result = 0;
    const auto [end, error] = std::from_chars(
        value.data(), value.data() + value.size(), result);
    if (error != std::errc{} || end != value.data() + value.size()) {
        throw std::runtime_error(std::string("Invalid ") + field + ": " +
                                 std::string(value));
    }
    return result;
}

std::vector<std::string_view> splitTabs(std::string_view line) {
    if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
    std::vector<std::string_view> fields;
    size_t start = 0;
    while (true) {
        const size_t tab = line.find('\t', start);
        if (tab == std::string_view::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, tab - start));
        start = tab + 1;
    }
    return fields;
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

std::vector<SequenceRecord> managerRecords(
    const SeqPro::SharedManagerVariant& manager) {
    if (!manager) throw std::runtime_error("Null sequence manager");
    return std::visit([](const auto& pointer) {
        if (!pointer) throw std::runtime_error("Null sequence manager pointer");
        std::vector<SequenceRecord> records;
        const auto names = pointer->getSequenceNames();
        records.reserve(names.size());
        for (const auto& name : names) {
            if (name.empty() || name.find_first_of("\t\r\n ") !=
                                    std::string::npos) {
                throw std::runtime_error(
                    "FASTA record name is empty or contains whitespace: " +
                    name);
            }
            using Pointer = std::decay_t<decltype(pointer)>;
            uint64_t length = 0;
            if constexpr (std::is_same_v<Pointer,
                          std::unique_ptr<SeqPro::SequenceManager>>) {
                length = pointer->getSequenceLength(name);
            } else {
                length = pointer->getOriginalManager().getSequenceLength(name);
            }
            records.push_back(SequenceRecord{name, length});
        }
        return records;
    }, *manager);
}

std::string managerSequence(const SeqPro::SharedManagerVariant& manager,
                            const std::string& name, uint64_t length) {
    if (length > std::numeric_limits<SeqPro::Length>::max()) {
        throw std::runtime_error("Sequence is too long for RaMAx: " + name);
    }
    return std::visit([&](const auto& pointer) -> std::string {
        using Pointer = std::decay_t<decltype(pointer)>;
        if constexpr (std::is_same_v<Pointer,
                      std::unique_ptr<SeqPro::SequenceManager>>) {
            return pointer->getSubSequence(name, 0,
                static_cast<SeqPro::Length>(length));
        } else {
            return pointer->getOriginalManager().getSubSequence(name, 0,
                static_cast<SeqPro::Length>(length));
        }
    }, *manager);
}

void ensureUniqueNames(const SpeciesName& species,
                       const std::vector<SequenceRecord>& records) {
    if (records.empty()) {
        throw std::runtime_error("Genome has no FASTA records: " + species);
    }
    std::unordered_set<std::string> seen;
    for (const auto& record : records) {
        if (!seen.emplace(record.name).second) {
            throw std::runtime_error("Duplicate FASTA record name in " +
                                     species + ": " + record.name);
        }
    }
}

std::string safeName(std::string_view value) {
    std::string result;
    result.reserve(value.size());
    for (const unsigned char c : value) {
        if (std::isalnum(c) || c == '.' || c == '_' || c == '-') {
            result.push_back(static_cast<char>(c));
        } else {
            result.push_back('_');
        }
    }
    if (result.empty()) result = "genome";
    const size_t hash = std::hash<std::string_view>{}(value);
    std::ostringstream suffix;
    suffix << '_' << std::hex << hash;
    result += suffix.str();
    return result;
}

int64_t fileMtime(const std::filesystem::path& path) {
    return static_cast<int64_t>(
        std::filesystem::last_write_time(path).time_since_epoch().count());
}

void atomicPublish(const std::filesystem::path& partial,
                   const std::filesystem::path& final_path) {
    std::error_code error;
    std::filesystem::rename(partial, final_path, error);
    if (!error) return;
    std::filesystem::remove(final_path, error);
    error.clear();
    std::filesystem::rename(partial, final_path, error);
    if (error) {
        throw std::runtime_error("Cannot publish " + final_path.string() +
                                 ": " + error.message());
    }
}

std::vector<WfmashRouterDetail::FaiRecord> readFai(
    const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot read FAI: " + path.string());
    return WfmashRouterDetail::parseFai(input);
}

void validateFaiRecords(
    const std::vector<WfmashRouterDetail::FaiRecord>& actual,
    const std::vector<SequenceRecord>& expected,
    const std::filesystem::path& path) {
    if (actual.size() != expected.size()) {
        throw std::runtime_error("FAI record count mismatch for " +
                                 path.string());
    }
    for (size_t i = 0; i < expected.size(); ++i) {
        if (actual[i].name != expected[i].name ||
            actual[i].length != expected[i].length) {
            throw std::runtime_error("FAI record order/name/length mismatch at " +
                                     std::to_string(i) + " for " + path.string());
        }
    }
}

std::string faiMarkerText(const std::filesystem::path& fasta,
                          std::string_view version) {
    return "schema\t1\npath\t" +
        std::filesystem::absolute(fasta).lexically_normal().string() +
        "\nsize\t" + std::to_string(std::filesystem::file_size(fasta)) +
        "\nmtime\t" + std::to_string(fileMtime(fasta)) +
        "\nsamtools\t" + std::string(version) + "\n";
}

std::string wfmashIndexMarkerText(const std::filesystem::path& fasta,
                                  std::string_view version) {
    return "schema\t1\npath\t" +
        std::filesystem::absolute(fasta).lexically_normal().string() +
        "\nsize\t" + std::to_string(std::filesystem::file_size(fasta)) +
        "\nmtime\t" + std::to_string(fileMtime(fasta)) +
        "\nwfmash\t" + std::string(version) +
        "\nparameters\t" + std::string(kWfmashParameterSummary) + "\n";
}

void writeAtomicText(const std::filesystem::path& path,
                     std::string_view text) {
    auto partial = path;
    partial += ".part";
    std::error_code ignored;
    std::filesystem::remove(partial, ignored);
    {
        std::ofstream output(partial, std::ios::binary | std::ios::trunc);
        if (!output) throw std::runtime_error("Cannot write " + partial.string());
        output.write(text.data(), static_cast<std::streamsize>(text.size()));
        output.flush();
        if (!output) throw std::runtime_error("Cannot finalize " + partial.string());
    }
    atomicPublish(partial, path);
}

void ensureFai(const std::filesystem::path& samtools,
               std::string_view samtools_version,
               const std::filesystem::path& fasta,
               const std::vector<SequenceRecord>& records,
               const std::filesystem::path& stderr_path) {
    std::filesystem::path fai = fasta;
    fai += ".fai";
    std::filesystem::path marker = fai;
    marker += ".ramax.complete";
    const std::string expected_marker = faiMarkerText(fasta, samtools_version);
    try {
        if (std::filesystem::is_regular_file(fai) &&
            std::filesystem::is_regular_file(marker) &&
            RaMAxExternalTool::readText(marker) == expected_marker) {
            validateFaiRecords(readFai(fai), records, fai);
            return;
        }
    } catch (const std::exception&) {
        // A stale or malformed cache is rebuilt by samtools below.
    }

    std::filesystem::path partial = fai;
    partial += ".part";
    std::filesystem::path stdout_path = stderr_path;
    stdout_path += ".stdout";
    std::error_code ignored;
    std::filesystem::remove(partial, ignored);
    const auto result = RaMAxExternalTool::run(
        samtools,
        {"faidx", "--fai-idx", partial.string(), fasta.string()},
        stdout_path, stderr_path);
    if (result.exit_code != 0) {
        throw std::runtime_error("samtools faidx failed for " + fasta.string() +
                                 " with exit code " +
                                 std::to_string(result.exit_code));
    }
    if (!std::filesystem::is_regular_file(partial) ||
        std::filesystem::file_size(partial) == 0) {
        throw std::runtime_error("samtools faidx did not create a non-empty FAI for " +
                                 fasta.string());
    }
    validateFaiRecords(readFai(partial), records, partial);
    atomicPublish(partial, fai);
    writeAtomicText(marker, expected_marker);
}

struct PreparedQuery {
    SpeciesName species;
    double distance{0.0};
    SeqPro::SharedManagerVariant manager;
    std::filesystem::path fasta;
    std::vector<SequenceRecord> records;
    std::unordered_map<std::string, std::string> alias_to_original;
    std::filesystem::path directory;
    std::string error;
};

void createAliasedView(
    PreparedQuery& query,
    const std::unordered_set<std::string>& reference_names,
    const std::filesystem::path& views_directory) {
    bool overlaps = false;
    for (const auto& record : query.records) {
        if (reference_names.contains(record.name)) {
            overlaps = true;
            break;
        }
    }
    if (!overlaps) return;

    std::filesystem::create_directories(views_directory);
    const std::string stem = safeName(query.species);
    const auto view = views_directory / (stem + ".fasta");
    auto partial = view;
    partial += ".part";
    std::vector<SequenceRecord> aliased;
    aliased.reserve(query.records.size());
    query.alias_to_original.clear();
    {
        std::ofstream output(partial, std::ios::binary | std::ios::trunc);
        if (!output) throw std::runtime_error("Cannot create wfmash FASTA view: " +
                                              partial.string());
        for (size_t i = 0; i < query.records.size(); ++i) {
            const auto& record = query.records[i];
            const std::string alias = "RAMAXQ_" + stem + "_" +
                                      std::to_string(i);
            if (reference_names.contains(alias) ||
                !query.alias_to_original.emplace(alias, record.name).second) {
                throw std::runtime_error("Cannot create unique wfmash alias for " +
                                         query.species);
            }
            aliased.push_back({alias, record.length});
            output << '>' << alias << '\n';
            const std::string sequence = managerSequence(
                query.manager, record.name, record.length);
            for (size_t offset = 0; offset < sequence.size(); offset += 80) {
                output.write(sequence.data() + offset,
                    static_cast<std::streamsize>(
                        std::min<size_t>(80, sequence.size() - offset)));
                output.put('\n');
            }
        }
        output.flush();
        if (!output) throw std::runtime_error("Cannot finalize wfmash FASTA view: " +
                                              partial.string());
    }
    atomicPublish(partial, view);
    query.fasta = view;
    query.records = std::move(aliased);

    std::ostringstream mapping;
    mapping << "wfmash_name\toriginal_name\n";
    for (const auto& record : query.records) {
        mapping << record.name << '\t'
                << query.alias_to_original.at(record.name) << '\n';
    }
    writeAtomicText(query.directory / "name_map.tsv", mapping.str());
}

std::vector<std::string> mappingArguments(
    uint_t threads, const std::filesystem::path& tmp_directory,
    const std::filesystem::path& reference,
    const std::filesystem::path& query,
    const std::optional<std::pair<bool, std::filesystem::path>>& index) {
    std::vector<std::string> args{
        "-m", "-w", "5000", "-l", "25000", "-p", "90", "-n", "1",
        "-k", "19", "-F", "0.001", "--hg-filter=1.0,30,99.9",
        "-t", std::to_string(threads), "--tmp-base", tmp_directory.string()
    };
    if (index) {
        args.push_back(std::string(index->first ? "--write-index=" :
                                                "--read-index=") +
                       index->second.string());
    }
    args.push_back(reference.string());
    args.push_back(query.string());
    return args;
}

std::vector<std::string> alignmentArguments(
    uint_t threads, const std::filesystem::path& tmp_directory,
    const std::filesystem::path& mapping,
    const std::filesystem::path& reference,
    const std::filesystem::path& query) {
    return {"-i", mapping.string(), "-t", std::to_string(threads),
            "--tmp-base", tmp_directory.string(), reference.string(),
            query.string()};
}

std::vector<std::string> indexArguments(
    uint_t threads, const std::filesystem::path& tmp_directory,
    const std::filesystem::path& index,
    const std::filesystem::path& reference) {
    return {
        "-m", "-w", "5000", "-l", "25000", "-p", "90", "-n", "1",
        "-k", "19", "-F", "0.001", "--hg-filter=1.0,30,99.9",
        "-t", std::to_string(threads), "--tmp-base", tmp_directory.string(),
        "--write-index=" + index.string(), reference.string()
    };
}

std::unordered_map<std::string, SequenceRecord> recordMap(
    const std::vector<SequenceRecord>& records) {
    std::unordered_map<std::string, SequenceRecord> result;
    for (const auto& record : records) result.emplace(record.name, record);
    return result;
}

std::vector<ParsedPafRecord> parseAndValidatePaf(
    const std::filesystem::path& path,
    const std::vector<SequenceRecord>& reference_records,
    const std::vector<SequenceRecord>& query_records,
    bool require_cigar) {
    const auto references = recordMap(reference_records);
    const auto queries = recordMap(query_records);
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot read wfmash PAF: " + path.string());
    std::vector<ParsedPafRecord> parsed;
    std::string line;
    size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty()) continue;
        try {
            auto record = WfmashRouterDetail::parsePafLine(line, require_cigar);
            const auto target = references.find(record.target_name);
            const auto query = queries.find(record.query_name);
            if (target == references.end() || query == queries.end()) {
                throw std::runtime_error("PAF contains an unknown sequence name");
            }
            if (target->second.length != record.target_length ||
                query->second.length != record.query_length) {
                throw std::runtime_error("PAF sequence length disagrees with FASTA");
            }
            parsed.push_back(std::move(record));
        } catch (const std::exception& error) {
            throw std::runtime_error("Invalid PAF " + path.string() + ":" +
                                     std::to_string(line_number) + ": " +
                                     error.what());
        }
    }
    if (!input.eof()) throw std::runtime_error("Cannot finish reading PAF: " + path.string());
    if (parsed.empty()) {
        throw std::runtime_error("wfmash PAF has no alignments: " + path.string());
    }
    return parsed;
}

void sortAndDeduplicatePaf(std::vector<ParsedPafRecord>& parsed) {
    std::sort(parsed.begin(), parsed.end(), [](const auto& left, const auto& right) {
        return std::tie(left.target_name, left.target_start, left.target_end,
                        left.query_name, left.query_start, left.query_end,
                        left.strand, left.cigar_text) <
               std::tie(right.target_name, right.target_start, right.target_end,
                        right.query_name, right.query_start, right.query_end,
                        right.strand, right.cigar_text);
    });
    parsed.erase(std::unique(parsed.begin(), parsed.end(),
        [](const auto& left, const auto& right) {
            return left.target_name == right.target_name &&
                   left.target_start == right.target_start &&
                   left.target_end == right.target_end &&
                   left.query_name == right.query_name &&
                   left.query_start == right.query_start &&
                   left.query_end == right.query_end &&
                   left.strand == right.strand &&
                   left.cigar_text == right.cigar_text;
        }), parsed.end());
}

void writeCanonicalPaf(const std::filesystem::path& path,
                       const std::vector<ParsedPafRecord>& parsed) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("Cannot rewrite validated PAF: " + path.string());
    for (const auto& record : parsed) {
        output << record.query_name << '\t' << record.query_length << '\t'
               << record.query_start << '\t' << record.query_end << '\t'
               << (record.strand == Strand::FORWARD ? '+' : '-') << '\t'
               << record.target_name << '\t' << record.target_length << '\t'
               << record.target_start << '\t' << record.target_end << '\t'
               << record.matches << '\t' << record.block_length << '\t'
               << record.mapq << "\tcg:Z:" << record.cigar_text << '\n';
    }
    output.flush();
    if (!output) throw std::runtime_error("Cannot finalize validated PAF: " + path.string());
}

AnchorVec makeAnchors(
    std::vector<ParsedPafRecord> parsed,
    const SeqPro::SharedManagerVariant& reference_manager,
    const SeqPro::SharedManagerVariant& query_manager,
    const std::unordered_map<std::string, std::string>& alias_to_original) {
    sortAndDeduplicatePaf(parsed);

    AnchorVec anchors;
    anchors.reserve(parsed.size());
    for (auto& record : parsed) {
        std::string original_query = record.query_name;
        if (!alias_to_original.empty()) {
            const auto alias = alias_to_original.find(record.query_name);
            if (alias == alias_to_original.end()) {
                throw std::runtime_error("PAF query alias is absent from name map: " +
                                         record.query_name);
            }
            original_query = alias->second;
        }
        auto sequenceId = [](const SeqPro::SharedManagerVariant& manager,
                             const std::string& name) -> ChrIndex {
            return std::visit([&](const auto& pointer) -> ChrIndex {
                return pointer->getSequenceId(name);
            }, *manager);
        };
        const uint64_t ref_span = record.target_end - record.target_start;
        const uint64_t query_span = record.query_end - record.query_start;
        const uint64_t maximum = std::numeric_limits<uint32_t>::max();
        if (record.target_start > maximum || record.query_start > maximum ||
            ref_span > maximum || query_span > maximum ||
            record.cigar_columns > maximum || record.matches > maximum) {
            throw std::runtime_error("PAF coordinate exceeds RaMAx 32-bit Anchor range");
        }
        anchors.emplace_back(
            sequenceId(reference_manager, record.target_name),
            static_cast<Coord_t>(record.target_start),
            static_cast<Length_t>(ref_span),
            sequenceId(query_manager, original_query),
            static_cast<Coord_t>(record.query_start),
            static_cast<Length_t>(query_span), record.strand,
            static_cast<uint_t>(record.cigar_columns),
            static_cast<uint_t>(record.matches), std::move(record.cigar));
    }
    return anchors;
}

struct PairExecutionResult {
    bool success{false};
    AnchorVec anchors;
    std::string error;
};

PairExecutionResult executePair(
    const std::filesystem::path& wfmash,
    uint_t pair_threads,
    const std::filesystem::path& reference_fasta,
    const std::vector<SequenceRecord>& reference_records,
    const SeqPro::SharedManagerVariant& reference_manager,
    const PreparedQuery& query,
    const std::optional<std::pair<bool, std::filesystem::path>>& index) {
    try {
        std::filesystem::create_directories(query.directory / "tmp");
        const auto mappings = query.directory / "mappings.paf";
        auto mappings_partial = mappings;
        mappings_partial += ".part";
        const auto mapping_stderr = query.directory / "mapping.stderr.log";
        const auto mapping_result = RaMAxExternalTool::run(
            wfmash,
            mappingArguments(pair_threads, query.directory / "tmp",
                             reference_fasta, query.fasta, index),
            mappings_partial, mapping_stderr);
        if (mapping_result.exit_code != 0) {
            throw std::runtime_error("wfmash mapping exited " +
                                     std::to_string(mapping_result.exit_code));
        }
        parseAndValidatePaf(mappings_partial, reference_records,
                            query.records, false);

        const auto alignment = query.directory / "alignment.paf";
        auto alignment_partial = alignment;
        alignment_partial += ".part";
        const auto alignment_stderr = query.directory / "alignment.stderr.log";
        const auto first_alignment_stderr =
            query.directory / "alignment.attempt1.stderr.log";
        auto alignment_result = RaMAxExternalTool::run(
            wfmash,
            alignmentArguments(pair_threads, query.directory / "tmp",
                               mappings_partial,
                               reference_fasta, query.fasta),
            alignment_partial, first_alignment_stderr);
        if (alignment_result.exit_code != 0) {
            const int first_exit_code = alignment_result.exit_code;
            spdlog::warn(
                "[wfmash-router] {} parallel alignment exited {}; "
                "retrying once under the process-wide serial guard",
                query.species, first_exit_code);
            std::lock_guard<std::mutex> retry_guard(
                wfmash_alignment_retry_mutex);
            std::ifstream mapping_check(mappings_partial, std::ios::binary);
            if (!mapping_check || mapping_check.peek() == std::ifstream::traits_type::eof()) {
                throw std::runtime_error(
                    "wfmash alignment retry cannot read non-empty mapping PAF after "
                    "initial exit " + std::to_string(first_exit_code));
            }
            mapping_check.close();
            alignment_result = RaMAxExternalTool::run(
                wfmash,
                alignmentArguments(pair_threads, query.directory / "tmp",
                                   mappings_partial,
                                   reference_fasta, query.fasta),
                alignment_partial, alignment_stderr);
            if (alignment_result.exit_code != 0) {
                throw std::runtime_error(
                    "wfmash alignment exited " +
                    std::to_string(first_exit_code) +
                    " and serial retry exited " +
                    std::to_string(alignment_result.exit_code));
            }
            spdlog::info(
                "[wfmash-router] {} serial alignment retry succeeded",
                query.species);
        } else {
            atomicPublish(first_alignment_stderr, alignment_stderr);
        }
        auto parsed = parseAndValidatePaf(
            alignment_partial, reference_records, query.records, true);
        sortAndDeduplicatePaf(parsed);
        const auto normalization =
            WfmashRouterDetail::normalizePafForGraph(parsed);
        spdlog::info(
            "[wfmash-router] {} graph-safe PAF normalization: input={} "
            "trimmed={} skipped={} output={}",
            query.species, normalization.input_records,
            normalization.trimmed_records, normalization.skipped_records,
            parsed.size());
        writeCanonicalPaf(alignment_partial, parsed);
        auto anchors = makeAnchors(std::move(parsed), reference_manager,
                                   query.manager, query.alias_to_original);
        // On DrvFS, a freshly renamed mapping file can be temporarily invisible
        // to another process. Alignment therefore consumes the validated .part
        // path and only publishes the stable user-facing name afterward.
        atomicPublish(mappings_partial, mappings);
        atomicPublish(alignment_partial, alignment);
        return {true, std::move(anchors), {}};
    } catch (const std::exception& error) {
        return {false, {}, error.what()};
    }
}

}  // namespace

namespace WfmashRouterDetail {

std::string validateSamtoolsVersion(std::string_view output) {
    std::istringstream input{std::string(output)};
    std::string first;
    std::string second;
    std::getline(input, first);
    std::getline(input, second);
    if (!first.empty() && first.back() == '\r') first.pop_back();
    if (!second.empty() && second.back() == '\r') second.pop_back();
    if (first != "samtools 1.24" || second != "Using htslib 1.24") {
        throw std::runtime_error(
            "RaMAx requires samtools 1.24 with HTSlib 1.24; observed: " +
            first + (second.empty() ? "" : " / " + second));
    }
    return "samtools 1.24; HTSlib 1.24";
}

std::string validateWfmashVersion(std::string_view output) {
    std::string first(output.substr(0, output.find_first_of("\r\n")));
    if (first != "v0.24.2" && first.rfind("v0.24.2-", 0) != 0) {
        throw std::runtime_error(
            "RaMAx requires wfmash 0.24.2; observed: " + first);
    }
    return first;
}

std::vector<FaiRecord> parseFai(std::istream& input) {
    std::vector<FaiRecord> result;
    std::string line;
    size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty()) throw std::runtime_error("Empty FAI line " + std::to_string(line_number));
        const auto fields = splitTabs(line);
        if (fields.size() != 5 || fields[0].empty()) {
            throw std::runtime_error("FAI line must contain exactly five fields at line " +
                                     std::to_string(line_number));
        }
        FaiRecord record;
        record.name = std::string(fields[0]);
        record.length = parseUnsigned(fields[1], "FAI LENGTH");
        record.offset = parseUnsigned(fields[2], "FAI OFFSET");
        record.line_bases = parseUnsigned(fields[3], "FAI LINEBASES");
        record.line_width = parseUnsigned(fields[4], "FAI LINEWIDTH");
        if (record.length == 0 || record.line_bases == 0 ||
            record.line_width < record.line_bases) {
            throw std::runtime_error("Invalid FAI numeric relationship at line " +
                                     std::to_string(line_number));
        }
        result.push_back(std::move(record));
    }
    if (result.empty()) throw std::runtime_error("FAI is empty");
    return result;
}

ParsedPafRecord parsePafLine(std::string_view line, bool require_cigar) {
    const auto fields = splitTabs(line);
    if (fields.size() < 12) {
        throw std::runtime_error("PAF has fewer than 12 fields");
    }
    ParsedPafRecord record;
    record.query_name = std::string(fields[0]);
    record.query_length = parseUnsigned(fields[1], "PAF query length");
    record.query_start = parseUnsigned(fields[2], "PAF query start");
    record.query_end = parseUnsigned(fields[3], "PAF query end");
    if (fields[4] == "+") record.strand = Strand::FORWARD;
    else if (fields[4] == "-") record.strand = Strand::REVERSE;
    else throw std::runtime_error("PAF strand must be + or -");
    record.target_name = std::string(fields[5]);
    record.target_length = parseUnsigned(fields[6], "PAF target length");
    record.target_start = parseUnsigned(fields[7], "PAF target start");
    record.target_end = parseUnsigned(fields[8], "PAF target end");
    record.matches = parseUnsigned(fields[9], "PAF matches");
    record.block_length = parseUnsigned(fields[10], "PAF block length");
    record.mapq = parseUnsigned(fields[11], "PAF mapping quality");
    if (record.query_name.empty() || record.target_name.empty() ||
        record.query_start >= record.query_end ||
        record.target_start >= record.target_end ||
        record.matches > record.block_length || record.mapq > 255) {
        throw std::runtime_error("PAF has invalid names, bounds, or counts");
    }
    if (require_cigar &&
        (record.query_end > record.query_length ||
         record.target_end > record.target_length)) {
        throw std::runtime_error("Final PAF coordinates exceed sequence bounds");
    }
    for (size_t i = 12; i < fields.size(); ++i) {
        if (fields[i].rfind("cg:Z:", 0) == 0) {
            if (!record.cigar_text.empty()) {
                throw std::runtime_error("PAF contains duplicate cg tags");
            }
            record.cigar_text = std::string(fields[i].substr(5));
        }
    }
    if (!require_cigar) return record;
    if (record.cigar_text.empty()) throw std::runtime_error("PAF is missing cg:Z");

    uint64_t target_consumed = 0;
    uint64_t query_consumed = 0;
    uint64_t columns = 0;
    size_t offset = 0;
    while (offset < record.cigar_text.size()) {
        const size_t digit_start = offset;
        while (offset < record.cigar_text.size() &&
               std::isdigit(static_cast<unsigned char>(record.cigar_text[offset]))) {
            ++offset;
        }
        if (digit_start == offset || offset == record.cigar_text.size()) {
            throw std::runtime_error("Malformed cg CIGAR");
        }
        const uint64_t length = parseUnsigned(
            std::string_view(record.cigar_text).substr(digit_start,
                                                       offset - digit_start),
            "CIGAR length");
        if (length == 0 || length > 0x0fffffffU) {
            throw std::runtime_error("CIGAR operation length is invalid");
        }
        const char operation = record.cigar_text[offset++];
        if (operation != 'M' && operation != '=' && operation != 'X' &&
            operation != 'I' && operation != 'D') {
            throw std::runtime_error("Unsupported CIGAR operation");
        }
        appendCigarOp(record.cigar, operation, static_cast<uint32_t>(length));
        if (operation == 'M' || operation == '=' || operation == 'X') {
            target_consumed += length;
            query_consumed += length;
        } else if (operation == 'I') {
            query_consumed += length;
        } else {
            target_consumed += length;
        }
        columns += length;
    }
    record.cigar_columns = columns;
    if (target_consumed != record.target_end - record.target_start ||
        query_consumed != record.query_end - record.query_start) {
        throw std::runtime_error("CIGAR consumption disagrees with PAF spans");
    }
    return record;
}

namespace {

using IntervalSet = std::map<uint64_t, uint64_t>;

uint64_t forwardPrefixOverlap(const IntervalSet& intervals,
                              uint64_t start) {
    auto next = intervals.upper_bound(start);
    if (next == intervals.begin()) return 0;
    const auto previous = std::prev(next);
    return previous->second > start ? previous->second - start : 0;
}

uint64_t reversePrefixOverlap(const IntervalSet& intervals,
                              uint64_t end) {
    auto next = intervals.lower_bound(end);
    if (next == intervals.begin()) return 0;
    const auto previous = std::prev(next);
    return previous->second >= end && previous->first < end
        ? end - previous->first : 0;
}

bool overlaps(const IntervalSet& intervals, uint64_t start, uint64_t end) {
    auto next = intervals.lower_bound(start);
    if (next != intervals.end() && next->first < end) return true;
    if (next == intervals.begin()) return false;
    return std::prev(next)->second > start;
}

bool trimPafPrefix(ParsedPafRecord& record,
                   uint64_t target_required,
                   uint64_t query_required) {
    if (target_required == 0 && query_required == 0) return true;
    Cigar_t remaining;
    uint64_t target_removed = 0;
    uint64_t query_removed = 0;
    uint64_t matches_removed = 0;
    bool trimming = true;
    for (const CigarUnit unit : record.cigar) {
        char operation = 0;
        uint32_t length = 0;
        intToCigar(unit, operation, length);
        const bool consumes_target =
            operation == 'M' || operation == '=' || operation == 'X' ||
            operation == 'D';
        const bool consumes_query =
            operation == 'M' || operation == '=' || operation == 'X' ||
            operation == 'I';
        uint32_t removed = 0;
        if (trimming) {
            const uint64_t target_need = target_removed < target_required
                ? target_required - target_removed : 0;
            const uint64_t query_need = query_removed < query_required
                ? query_required - query_removed : 0;
            if (target_need == 0 && query_need == 0) {
                trimming = false;
            } else {
                uint64_t useful_need = 0;
                if (consumes_target) useful_need = std::max(useful_need, target_need);
                if (consumes_query) useful_need = std::max(useful_need, query_need);
                removed = useful_need == 0
                    ? length
                    : static_cast<uint32_t>(std::min<uint64_t>(length, useful_need));
                if (consumes_target) target_removed += removed;
                if (consumes_query) query_removed += removed;
                if (operation == 'M' || operation == '=') {
                    matches_removed += removed;
                }
                // If this operation cannot satisfy the other coordinate's
                // outstanding trim, its remainder is still before the future
                // cut point and must be removed as well.
                if (removed < length &&
                    (target_removed < target_required ||
                     query_removed < query_required)) {
                    const uint32_t extra = length - removed;
                    removed = length;
                    if (consumes_target) target_removed += extra;
                    if (consumes_query) query_removed += extra;
                    if (operation == 'M' || operation == '=') {
                        matches_removed += extra;
                    }
                }
                if (removed < length &&
                    target_removed >= target_required &&
                    query_removed >= query_required) {
                    trimming = false;
                }
            }
        }
        if (!trimming || removed < length) {
            appendCigarOp(remaining, operation, length - removed);
        }
    }
    if (target_removed < target_required || query_removed < query_required ||
        remaining.empty()) {
        return false;
    }
    record.target_start += target_removed;
    if (record.strand == Strand::FORWARD) record.query_start += query_removed;
    else record.query_end -= query_removed;
    if (record.target_start >= record.target_end ||
        record.query_start >= record.query_end) {
        return false;
    }
    record.matches = matches_removed >= record.matches
        ? 0 : record.matches - matches_removed;
    record.cigar = std::move(remaining);
    record.cigar_text = cigarToString(record.cigar);
    record.cigar_columns = countAlignmentLength(record.cigar);
    record.block_length = std::max(record.target_end - record.target_start,
                                   record.query_end - record.query_start);
    return true;
}

}  // namespace

PafNormalizationStats normalizePafForGraph(
    std::vector<ParsedPafRecord>& records) {
    PafNormalizationStats stats;
    stats.input_records = records.size();
    std::map<std::string, IntervalSet> target_intervals;
    std::map<std::string, IntervalSet> query_intervals;
    std::vector<ParsedPafRecord> normalized;
    normalized.reserve(records.size());
    for (auto& record : records) {
        auto& targets = target_intervals[record.target_name];
        auto& queries = query_intervals[record.query_name];
        const uint64_t target_trim =
            forwardPrefixOverlap(targets, record.target_start);
        const uint64_t query_trim = record.strand == Strand::FORWARD
            ? forwardPrefixOverlap(queries, record.query_start)
            : reversePrefixOverlap(queries, record.query_end);
        if ((target_trim != 0 || query_trim != 0) &&
            !trimPafPrefix(record, target_trim, query_trim)) {
            ++stats.skipped_records;
            continue;
        }
        if (target_trim != 0 || query_trim != 0) ++stats.trimmed_records;
        if (overlaps(targets, record.target_start, record.target_end) ||
            overlaps(queries, record.query_start, record.query_end)) {
            ++stats.skipped_records;
            continue;
        }
        targets.emplace(record.target_start, record.target_end);
        queries.emplace(record.query_start, record.query_end);
        normalized.push_back(std::move(record));
    }
    records = std::move(normalized);
    return stats;
}

size_t workerCount(size_t tasks, uint_t total_threads) {
    if (tasks == 0) return 0;
    return std::max<size_t>(1, std::min<size_t>(tasks, total_threads));
}

uint_t threadsPerTask(size_t tasks, uint_t total_threads) {
    if (tasks == 0) return 1;
    return std::max<uint_t>(1, total_threads / static_cast<uint_t>(tasks));
}

}  // namespace WfmashRouterDetail

FirstRoundWfmashRouter::FirstRoundWfmashRouter(
    std::filesystem::path samtools_executable,
    std::filesystem::path wfmash_executable,
    std::filesystem::path output_directory,
    uint_t threads)
    : samtools_executable_(std::move(samtools_executable)),
      wfmash_executable_(std::move(wfmash_executable)),
      output_directory_(std::move(output_directory)),
      threads_(std::max<uint_t>(1, threads)) {
    if (!RaMAxExternalTool::isExecutable(samtools_executable_)) {
        throw std::runtime_error("samtools is required but was not found");
    }
    if (!RaMAxExternalTool::isExecutable(wfmash_executable_)) {
        throw std::runtime_error("wfmash is required but was not found");
    }
    std::filesystem::create_directories(output_directory_);
    const auto samtools_stdout = output_directory_ / "samtools.version.stdout";
    const auto samtools_stderr = output_directory_ / "samtools.version.stderr";
    const auto samtools_result = RaMAxExternalTool::run(
        samtools_executable_, {"--version"}, samtools_stdout, samtools_stderr);
    if (samtools_result.exit_code != 0) {
        throw std::runtime_error("samtools --version failed");
    }
    samtools_version_ = WfmashRouterDetail::validateSamtoolsVersion(
        RaMAxExternalTool::readText(samtools_stdout));

    const auto wfmash_stdout = output_directory_ / "wfmash.version.stdout";
    const auto wfmash_stderr = output_directory_ / "wfmash.version.stderr";
    const auto wfmash_result = RaMAxExternalTool::run(
        wfmash_executable_, {"--version"}, wfmash_stdout, wfmash_stderr);
    if (wfmash_result.exit_code != 0) {
        throw std::runtime_error("wfmash --version failed");
    }
    const std::string wfmash_version_output =
        RaMAxExternalTool::readText(wfmash_stdout) +
        RaMAxExternalTool::readText(wfmash_stderr);
    wfmash_version_ = WfmashRouterDetail::validateWfmashVersion(
        wfmash_version_output);
}

FirstRoundWfmashResult FirstRoundWfmashRouter::run(
    const SpeciesName& reference,
    const std::vector<MashDistanceRecord>& distances,
    double near_distance,
    double far_distance,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers) {
    FirstRoundWfmashResult result;
    const auto reference_it = managers.find(reference);
    if (reference_it == managers.end()) {
        throw std::runtime_error("wfmash reference manager is absent: " + reference);
    }
    const auto reference_records = managerRecords(reference_it->second);
    ensureUniqueNames(reference, reference_records);
    std::unordered_set<std::string> reference_names;
    for (const auto& record : reference_records) reference_names.emplace(record.name);

    std::vector<PreparedQuery> queries;
    for (const auto& distance : distances) {
        if (distance.reference != reference || distance.distance >= near_distance) continue;
        const auto manager = managers.find(distance.query);
        if (manager == managers.end()) {
            throw std::runtime_error("wfmash query manager is absent: " + distance.query);
        }
        PreparedQuery query;
        query.species = distance.query;
        query.distance = distance.distance;
        query.manager = manager->second;
        query.fasta = fastaPath(manager->second);
        query.records = managerRecords(manager->second);
        ensureUniqueNames(query.species, query.records);
        query.directory = output_directory_ / safeName(query.species);
        std::filesystem::create_directories(query.directory);
        createAliasedView(query, reference_names, output_directory_ / "views");
        queries.push_back(std::move(query));
    }
    std::sort(queries.begin(), queries.end(), [](const auto& left, const auto& right) {
        return left.species < right.species;
    });

    const auto reference_fasta = fastaPath(reference_it->second);
    ensureFai(samtools_executable_, samtools_version_, reference_fasta,
              reference_records, output_directory_ / "reference.faidx.stderr.log");

    std::map<std::string, std::vector<size_t>> query_indices_by_fasta;
    for (size_t i = 0; i < queries.size(); ++i) {
        query_indices_by_fasta[
            std::filesystem::absolute(queries[i].fasta).lexically_normal().string()]
            .push_back(i);
    }
    std::vector<std::vector<size_t>> fai_groups;
    for (auto& [unused, indices] : query_indices_by_fasta) {
        fai_groups.push_back(std::move(indices));
    }
    if (!fai_groups.empty()) {
        const int workers = static_cast<int>(
            WfmashRouterDetail::workerCount(fai_groups.size(), threads_));
#pragma omp parallel for schedule(dynamic) num_threads(workers)
        for (size_t group_index = 0; group_index < fai_groups.size(); ++group_index) {
            const auto& indices = fai_groups[group_index];
            const size_t representative = indices.front();
            try {
                ensureFai(samtools_executable_, samtools_version_,
                          queries[representative].fasta,
                          queries[representative].records,
                          queries[representative].directory / "faidx.stderr.log");
                for (const size_t index : indices) {
                    if (queries[index].records != queries[representative].records) {
                        throw std::runtime_error(
                            "Species sharing a FASTA path have inconsistent sequence metadata");
                    }
                }
            } catch (const std::exception& error) {
                for (const size_t index : indices) queries[index].error = error.what();
            }
        }
    }

    std::vector<size_t> runnable;
    for (size_t i = 0; i < queries.size(); ++i) {
        if (queries[i].error.empty()) runnable.push_back(i);
        else spdlog::warn("[wfmash-router] {} falls back to legacy: {}",
                          queries[i].species, queries[i].error);
    }
    const uint_t pair_threads = WfmashRouterDetail::threadsPerTask(
        std::max<size_t>(1, queries.size()), threads_);
    const size_t workers = WfmashRouterDetail::workerCount(runnable.size(), threads_);
    spdlog::info(
        "[wfmash-router] reference={} near={} far={} candidates={} workers={} threads_per_pair={} params={}",
        reference, near_distance, far_distance, queries.size(), workers,
        pair_threads, kWfmashParameterSummary);

    std::vector<PairExecutionResult> pair_results(queries.size());
    std::optional<std::filesystem::path> index_path;
    if (runnable.size() > 1) {
        const auto final_index = output_directory_ / "reference.wfmash.index";
        const auto index_marker = output_directory_ /
                                  "reference.wfmash.index.ramax.complete";
        auto partial_index = final_index;
        partial_index += ".part";
        const std::string expected_marker = wfmashIndexMarkerText(
            reference_fasta, wfmash_version_);
        bool reusable = false;
        try {
            reusable = std::filesystem::is_regular_file(final_index) &&
                std::filesystem::file_size(final_index) > 0 &&
                std::filesystem::is_regular_file(index_marker) &&
                RaMAxExternalTool::readText(index_marker) == expected_marker;
        } catch (const std::exception&) {
            reusable = false;
        }
        if (reusable) {
            index_path = final_index;
        } else {
            std::error_code ignored;
            std::filesystem::remove(partial_index, ignored);
            std::filesystem::create_directories(output_directory_ / "index_tmp");
            const auto index_result = RaMAxExternalTool::run(
                wfmash_executable_,
                indexArguments(pair_threads, output_directory_ / "index_tmp",
                               partial_index, reference_fasta),
                output_directory_ / "index.stdout",
                output_directory_ / "index.stderr.log");
            if (index_result.exit_code == 0 &&
                std::filesystem::is_regular_file(partial_index) &&
                std::filesystem::file_size(partial_index) > 0) {
                atomicPublish(partial_index, final_index);
                writeAtomicText(index_marker, expected_marker);
                index_path = final_index;
            } else {
                for (const size_t index : runnable) {
                    pair_results[index].error =
                        "shared reference index could not be constructed";
                }
                runnable.clear();
            }
        }
    }

    if (!runnable.empty()) {
        const int omp_workers = static_cast<int>(
            WfmashRouterDetail::workerCount(runnable.size(), threads_));
#pragma omp parallel for schedule(dynamic) num_threads(omp_workers)
        for (size_t position = 0; position < runnable.size(); ++position) {
            const size_t index = runnable[position];
            std::optional<std::pair<bool, std::filesystem::path>> index_option;
            if (index_path) index_option = std::pair{false, *index_path};
            pair_results[index] = executePair(
                wfmash_executable_, pair_threads, reference_fasta,
                reference_records, reference_it->second, queries[index],
                index_option);
        }
    }

    std::ostringstream routing;
    routing << "species\tdistance\tbackend\tstatus\tdetail\n";
    std::unordered_set<SpeciesName> routed_candidates;
    for (size_t i = 0; i < queries.size(); ++i) {
        routed_candidates.emplace(queries[i].species);
        if (pair_results[i].success) {
            result.successful_species.emplace(queries[i].species);
            result.anchors_by_species.emplace(
                queries[i].species, std::move(pair_results[i].anchors));
            routing << queries[i].species << '\t' << std::setprecision(17)
                    << queries[i].distance << "\twfmash\tsuccess\t"
                    << result.anchors_by_species.at(queries[i].species).size()
                    << " anchors\n";
        } else {
            const std::string detail = !queries[i].error.empty()
                ? queries[i].error : pair_results[i].error;
            routing << queries[i].species << '\t' << std::setprecision(17)
                    << queries[i].distance << "\tlegacy\tfallback\t"
                    << detail << '\n';
            spdlog::warn("[wfmash-router] {} falls back to legacy: {}",
                         queries[i].species, detail);
        }
    }
    for (const auto& distance : distances) {
        if (distance.reference != reference ||
            routed_candidates.contains(distance.query)) {
            continue;
        }
        routing << distance.query << '\t' << std::setprecision(17)
                << distance.distance
                << "\tlegacy\tnot_routed\tdistance_not_below_near\n";
    }
    writeAtomicText(output_directory_ / "routing.tsv", routing.str());
    std::ostringstream versions;
    versions << "samtools\t" << samtools_executable_.string() << '\t'
             << samtools_version_ << '\n'
             << "wfmash\t" << wfmash_executable_.string() << '\t'
             << wfmash_version_ << '\n'
             << "parameters\t" << kWfmashParameterSummary << '\n';
    writeAtomicText(output_directory_ / "tools.tsv", versions.str());
    return result;
}

std::filesystem::path locateSamtoolsExecutable() {
    const auto executable = RaMAxExternalTool::locateExecutable(
        "samtools", RAMAX_SAMTOOLS_CONFIGURED_PATH);
    if (executable.empty()) {
        throw std::runtime_error(
            "samtools 1.24 was not found at the configured path, next to ramax, or in PATH");
    }
    return executable;
}

std::filesystem::path locateWfmashExecutable() {
    const auto executable = RaMAxExternalTool::locateExecutable(
        "wfmash", RAMAX_WFMASH_CONFIGURED_PATH);
    if (executable.empty()) {
        throw std::runtime_error(
            "wfmash 0.24.2 was not found at the configured path, next to ramax, or in PATH");
    }
    return executable;
}
