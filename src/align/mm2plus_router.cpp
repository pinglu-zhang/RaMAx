#include "mm2plus_router.h"

#include "external_tool.h"
#include "xxhash.h"

#include <spdlog/spdlog.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef RAMAX_MM2PLUS_CONFIGURED_PATH
#define RAMAX_MM2PLUS_CONFIGURED_PATH ""
#endif

namespace {

using ParsedPafRecord = WfmashRouterDetail::ParsedPafRecord;
using SequenceRecord = WfmashRouterDetail::SequenceRecord;

constexpr uint64_t kMm2plusMaximumSequenceLength = 2147483647ULL;
constexpr std::string_view kMm2plusParameters =
    "-x asm20 -c --eqx --secondary=no";

std::string firstLine(std::string_view text) {
    const size_t end = text.find_first_of("\r\n");
    return std::string(text.substr(0, end));
}

std::string safeName(std::string_view value) {
    std::string result;
    result.reserve(value.size());
    for (const unsigned char c : value) {
        result.push_back(std::isalnum(c) || c == '.' || c == '_' || c == '-'
            ? static_cast<char>(c) : '_');
    }
    if (result.empty()) result = "genome";
    const uint64_t hash = XXH3_64bits(value.data(), value.size());
    std::ostringstream suffix;
    suffix << '_' << std::hex << hash;
    result += suffix.str();
    return result;
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

std::filesystem::path fastaPath(const SeqPro::SharedManagerVariant& manager) {
    if (!manager) throw std::runtime_error("Null sequence manager");
    return std::visit([](const auto& pointer) -> std::filesystem::path {
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
    const SeqPro::SharedManagerVariant& manager,
    const SpeciesName& species) {
    if (!manager) throw std::runtime_error("Null sequence manager");
    return std::visit([&](const auto& pointer) {
        if (!pointer) throw std::runtime_error("Null sequence manager pointer");
        std::vector<SequenceRecord> records;
        std::set<std::string> names;
        for (const auto& name : pointer->getSequenceNames()) {
            if (name.empty() || name.find_first_of("\t\r\n ") != std::string::npos) {
                throw std::runtime_error(
                    "FASTA record name is empty or contains whitespace in " +
                    species + ": " + name);
            }
            if (!names.emplace(name).second) {
                throw std::runtime_error("Duplicate FASTA record name in " +
                                         species + ": " + name);
            }
            using Pointer = std::decay_t<decltype(pointer)>;
            uint64_t length = 0;
            if constexpr (std::is_same_v<Pointer,
                          std::unique_ptr<SeqPro::SequenceManager>>) {
                length = pointer->getSequenceLength(name);
            } else {
                length = pointer->getOriginalManager().getSequenceLength(name);
            }
            if (length == 0 || length > kMm2plusMaximumSequenceLength) {
                throw std::runtime_error(
                    "mm2-plus requires each FASTA record to be 1..2147483647 bp: " +
                    species + "." + name);
            }
            records.push_back({name, length});
        }
        if (records.empty()) {
            throw std::runtime_error("Genome has no FASTA records: " + species);
        }
        return records;
    }, *manager);
}

uint64_t totalBases(const std::vector<SequenceRecord>& records) {
    uint64_t total = 0;
    for (const auto& record : records) {
        if (record.length > std::numeric_limits<uint64_t>::max() - total) {
            throw std::runtime_error("Genome length overflow");
        }
        total += record.length;
    }
    return total;
}

std::string fileFingerprint(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot fingerprint " + path.string());
    XXH3_state_t* state = XXH3_createState();
    if (!state || XXH3_128bits_reset(state) == XXH_ERROR) {
        if (state) XXH3_freeState(state);
        throw std::runtime_error("Cannot initialize XXH3 fingerprint");
    }
    std::array<char, 1 << 20> buffer{};
    while (input) {
        input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const auto count = input.gcount();
        if (count > 0 && XXH3_128bits_update(
                state, buffer.data(), static_cast<size_t>(count)) == XXH_ERROR) {
            XXH3_freeState(state);
            throw std::runtime_error("Cannot update XXH3 fingerprint");
        }
    }
    if (!input.eof()) {
        XXH3_freeState(state);
        throw std::runtime_error("Cannot finish fingerprinting " + path.string());
    }
    XXH128_canonical_t canonical{};
    XXH128_canonicalFromHash(&canonical, XXH3_128bits_digest(state));
    XXH3_freeState(state);
    std::ostringstream output;
    output << std::hex << std::setfill('0');
    for (const unsigned char byte : canonical.digest) {
        output << std::setw(2) << static_cast<unsigned>(byte);
    }
    return output.str();
}

std::string qualifiedName(const SpeciesName& species,
                          const std::string& contig) {
    std::string result = species + "." + contig;
    if (result.empty() || result.find('#') != std::string::npos) {
        throw std::runtime_error("Invalid qualified PAF name: " + result);
    }
    for (const unsigned char c : result) {
        if (std::iscntrl(c) || std::isspace(c)) {
            throw std::runtime_error(
                "Qualified PAF name contains whitespace/control characters: " +
                result);
        }
    }
    return result;
}

std::unordered_map<std::string, SequenceRecord> recordMap(
    const std::vector<SequenceRecord>& records) {
    std::unordered_map<std::string, SequenceRecord> result;
    for (const auto& record : records) result.emplace(record.name, record);
    return result;
}

bool containsMatchOperator(const ParsedPafRecord& record) {
    for (const CigarUnit unit : record.cigar) {
        char operation = 0;
        uint32_t length = 0;
        intToCigar(unit, operation, length);
        if (operation == 'M') return true;
    }
    return false;
}

std::vector<ParsedPafRecord> parsePaf(
    const std::filesystem::path& path,
    const std::vector<SequenceRecord>& references,
    const std::vector<SequenceRecord>& queries,
    uint_t min_span,
    Mm2plusRouterDetail::NormalizationStats& stats) {
    const auto reference_map = recordMap(references);
    const auto query_map = recordMap(queries);
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot read mm2-plus PAF: " + path.string());
    std::vector<ParsedPafRecord> parsed;
    std::string line;
    size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty()) continue;
        ++stats.raw_records;
        try {
            auto record = WfmashRouterDetail::parsePafLine(line, true);
            const auto target = reference_map.find(record.target_name);
            const auto query = query_map.find(record.query_name);
            if (target == reference_map.end() || query == query_map.end()) {
                throw std::runtime_error("unknown query or target name");
            }
            if (target->second.length != record.target_length ||
                query->second.length != record.query_length) {
                throw std::runtime_error("sequence length disagrees with FASTA");
            }
            if (!record.alignment_type ||
                (*record.alignment_type != 'P' && *record.alignment_type != 'I')) {
                throw std::runtime_error("tp:A must be primary P or primary inversion I");
            }
            if (!record.alignment_score) {
                throw std::runtime_error("AS:i is required");
            }
            if (containsMatchOperator(record)) {
                throw std::runtime_error("--eqx PAF unexpectedly contains M CIGAR operators");
            }
            if (record.matches != record.equal_bases) {
                throw std::runtime_error("PAF matches disagree with = CIGAR consumption");
            }
            if (record.block_length != record.cigar_columns) {
                throw std::runtime_error("PAF block length disagrees with CIGAR columns");
            }
            const bool low_score = *record.alignment_score <= 0;
            const bool low_identity = record.block_length == 0 ||
                static_cast<double>(record.matches) /
                    static_cast<double>(record.block_length) < 0.5;
            if (low_score || low_identity) {
                ++stats.quality_discarded;
                stats.low_score_discarded += low_score;
                stats.low_identity_discarded += low_identity;
                continue;
            }
            if (record.target_end - record.target_start < min_span ||
                record.query_end - record.query_start < min_span) {
                continue;
            }
            parsed.push_back(std::move(record));
        } catch (const std::exception& error) {
            throw std::runtime_error("Invalid mm2-plus PAF " + path.string() +
                ":" + std::to_string(line_number) + ": " + error.what());
        }
    }
    if (!input.eof()) throw std::runtime_error("Cannot finish reading " + path.string());
    if (parsed.empty()) throw std::runtime_error("mm2-plus PAF has no usable alignments");
    return parsed;
}

void validateNoOverlaps(const std::vector<ParsedPafRecord>& records) {
    using Interval = std::pair<uint64_t, uint64_t>;
    std::map<std::string, std::vector<Interval>> targets;
    std::map<std::string, std::vector<Interval>> queries;
    for (const auto& record : records) {
        targets[record.target_name].emplace_back(
            record.target_start, record.target_end);
        queries[record.query_name].emplace_back(
            record.query_start, record.query_end);
    }
    auto check = [](auto& by_name, const char* axis) {
        for (auto& [name, intervals] : by_name) {
            std::sort(intervals.begin(), intervals.end());
            for (size_t i = 1; i < intervals.size(); ++i) {
                if (intervals[i - 1].second > intervals[i].first) {
                    throw std::runtime_error(std::string("Canonical mm2-plus ") +
                        axis + " intervals overlap for " + name);
                }
            }
        }
    };
    check(targets, "target");
    check(queries, "query");
}

void sortCanonical(std::vector<ParsedPafRecord>& records) {
    std::sort(records.begin(), records.end(), [](const auto& left, const auto& right) {
        return std::tie(left.target_name, left.target_start, left.target_end,
                        left.query_name, left.query_start, left.query_end,
                        left.strand, left.cigar_text) <
               std::tie(right.target_name, right.target_start, right.target_end,
                        right.query_name, right.query_start, right.query_end,
                        right.strand, right.cigar_text);
    });
}

void writeCanonicalPaf(const std::filesystem::path& path,
                       const SpeciesName& reference_species,
                       const SpeciesName& query_species,
                       const std::vector<ParsedPafRecord>& records) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("Cannot write " + path.string());
    for (const auto& record : records) {
        output << qualifiedName(query_species, record.query_name) << '\t'
               << record.query_length << '\t' << record.query_start << '\t'
               << record.query_end << '\t'
               << (record.strand == Strand::FORWARD ? '+' : '-') << '\t'
               << qualifiedName(reference_species, record.target_name) << '\t'
               << record.target_length << '\t' << record.target_start << '\t'
               << record.target_end << '\t' << record.matches << '\t'
               << record.cigar_columns << '\t' << record.mapq
               << "\ttp:A:" << *record.alignment_type
               << "\tAS:i:" << *record.alignment_score
               << "\tcg:Z:" << record.cigar_text << '\n';
    }
    output.flush();
    if (!output) throw std::runtime_error("Cannot finalize " + path.string());
}

ChrIndex sequenceId(const SeqPro::SharedManagerVariant& manager,
                    const std::string& name) {
    return std::visit([&](const auto& pointer) -> ChrIndex {
        return pointer->getSequenceId(name);
    }, *manager);
}

AnchorVec makeAnchors(std::vector<ParsedPafRecord> records,
                      const SeqPro::SharedManagerVariant& reference_manager,
                      const SeqPro::SharedManagerVariant& query_manager) {
    AnchorVec anchors;
    anchors.reserve(records.size());
    for (auto& record : records) {
        const uint64_t ref_span = record.target_end - record.target_start;
        const uint64_t query_span = record.query_end - record.query_start;
        const uint64_t maximum = std::numeric_limits<uint32_t>::max();
        if (record.target_start > maximum || record.query_start > maximum ||
            ref_span > maximum || query_span > maximum ||
            record.cigar_columns > maximum || record.matches > maximum) {
            throw std::runtime_error("mm2-plus coordinate exceeds RaMAx Anchor range");
        }
        anchors.emplace_back(
            sequenceId(reference_manager, record.target_name),
            static_cast<Coord_t>(record.target_start),
            static_cast<Length_t>(ref_span),
            sequenceId(query_manager, record.query_name),
            static_cast<Coord_t>(record.query_start),
            static_cast<Length_t>(query_span), record.strand,
            static_cast<uint_t>(record.cigar_columns),
            static_cast<uint_t>(record.matches), std::move(record.cigar));
    }
    return anchors;
}

std::string namesText(const SpeciesName& reference,
                      const SpeciesName& query,
                      const std::vector<SequenceRecord>& reference_records,
                      const std::vector<SequenceRecord>& query_records) {
    std::ostringstream output;
    output << "raw_name\tqualified_name\trole\tspecies\tlength\n";
    for (const auto& record : reference_records) {
        output << record.name << '\t' << qualifiedName(reference, record.name)
               << "\treference\t" << reference << '\t' << record.length << '\n';
    }
    for (const auto& record : query_records) {
        output << record.name << '\t' << qualifiedName(query, record.name)
               << "\tquery\t" << query << '\t' << record.length << '\n';
    }
    return output.str();
}

std::string indexManifest(
    const std::filesystem::path& reference,
    const std::vector<SequenceRecord>& records,
    std::string_view version,
    uint64_t bases) {
    std::ostringstream output;
    output << "schema\t1\nversion\t" << version
           << "\nparameters\t-x asm20 -I " << (bases + 1)
           << "\nreference\t"
           << std::filesystem::absolute(reference).lexically_normal().string()
           << "\nfingerprint_xxh3_128\t" << fileFingerprint(reference)
           << "\nrecord_count\t" << records.size()
           << "\ntotal_bases\t" << bases << '\n';
    for (size_t i = 0; i < records.size(); ++i) {
        output << "record\t" << i << '\t' << records[i].name << '\t'
               << records[i].length << '\n';
    }
    return output.str();
}

RaMAxExternalTool::CommandResult ensureIndex(
    const std::filesystem::path& executable,
    std::string_view version,
    uint_t threads,
    const std::filesystem::path& output_directory,
    const std::filesystem::path& reference,
    const std::vector<SequenceRecord>& records) {
    const uint64_t bases = totalBases(records);
    if (bases == std::numeric_limits<uint64_t>::max()) {
        throw std::runtime_error("Reference size cannot be represented for -I");
    }
    const auto index = output_directory / "reference.asm20.mmi";
    const auto manifest = output_directory / "reference.index.manifest.tsv";
    const std::string expected = indexManifest(reference, records, version, bases);
    if (std::filesystem::is_regular_file(index) &&
        std::filesystem::file_size(index) > 0 &&
        std::filesystem::is_regular_file(manifest) &&
        RaMAxExternalTool::readText(manifest) == expected) {
        spdlog::info("[mm2plus-router] reusing reference index {}", index.string());
        return {};
    }
    auto partial = index;
    partial += ".part";
    std::error_code ignored;
    std::filesystem::remove(partial, ignored);
    const auto result = RaMAxExternalTool::run(
        executable,
        Mm2plusRouterDetail::indexArguments(
            std::min<uint_t>(3, threads), bases, partial, reference),
        output_directory / "reference.index.stdout.log",
        output_directory / "reference.index.stderr.log");
    if (result.exit_code != 0 || !std::filesystem::is_regular_file(partial) ||
        std::filesystem::file_size(partial) == 0) {
        throw std::runtime_error("mm2-plus reference indexing failed with exit " +
                                 std::to_string(result.exit_code));
    }
    atomicPublish(partial, index);
    writeAtomicText(manifest, expected);
    return result;
}

struct PreparedQuery {
    SpeciesName species;
    double distance{0.0};
    SeqPro::SharedManagerVariant manager;
    std::filesystem::path fasta;
    std::vector<SequenceRecord> records;
    std::filesystem::path directory;
    uint64_t bases{0};
};

struct PairResult {
    bool success{false};
    AnchorVec anchors;
    Mm2plusRouterDetail::NormalizationStats normalization;
    RaMAxExternalTool::CommandResult command;
    double weighted_divergence{0.0};
    std::string error;
};

PairResult executePair(
    const std::filesystem::path& executable,
    uint_t pair_threads,
    uint_t min_span,
    const SpeciesName& reference_species,
    const std::filesystem::path& index,
    const std::vector<SequenceRecord>& reference_records,
    const SeqPro::SharedManagerVariant& reference_manager,
    const PreparedQuery& query) {
    PairResult result;
    try {
        std::filesystem::create_directories(query.directory);
        writeAtomicText(query.directory / "names.tsv",
            namesText(reference_species, query.species,
                      reference_records, query.records));
        const auto raw = query.directory / "raw.paf";
        auto raw_partial = raw;
        raw_partial += ".part";
        std::error_code ignored;
        std::filesystem::remove(raw_partial, ignored);
        result.command = RaMAxExternalTool::run(
            executable,
            Mm2plusRouterDetail::alignmentArguments(
                pair_threads, index, query.fasta),
            raw_partial, query.directory / "stderr.log");
        if (result.command.exit_code != 0) {
            throw std::runtime_error("mm2-plus exited " +
                                     std::to_string(result.command.exit_code));
        }
        Mm2plusRouterDetail::NormalizationStats quality_stats;
        auto parsed = parsePaf(raw_partial, reference_records,
                               query.records, min_span, quality_stats);
        result.normalization =
            Mm2plusRouterDetail::normalizeForGraph(parsed, min_span);
        result.normalization.raw_records = quality_stats.raw_records;
        result.normalization.quality_discarded = quality_stats.quality_discarded;
        result.normalization.low_score_discarded =
            quality_stats.low_score_discarded;
        result.normalization.low_identity_discarded =
            quality_stats.low_identity_discarded;
        if (parsed.empty()) {
            throw std::runtime_error("mm2-plus normalization removed all alignments");
        }
        validateNoOverlaps(parsed);
        sortCanonical(parsed);
        uint64_t total_columns = 0;
        long double divergent_columns = 0.0;
        for (const auto& record : parsed) {
            total_columns += record.cigar_columns;
            divergent_columns += static_cast<long double>(
                record.cigar_columns - record.matches);
        }
        result.weighted_divergence = total_columns == 0 ? 0.0 :
            static_cast<double>(divergent_columns /
                                static_cast<long double>(total_columns));

        const auto canonical = query.directory / "alignment.paf";
        auto canonical_partial = canonical;
        canonical_partial += ".part";
        writeCanonicalPaf(canonical_partial, reference_species,
                          query.species, parsed);
        result.anchors = makeAnchors(
            parsed, reference_manager, query.manager);

        std::ostringstream stats;
        stats << "metric\tvalue\n"
              << "raw_records\t" << result.normalization.raw_records << '\n'
              << "quality_discarded\t" << result.normalization.quality_discarded << '\n'
              << "low_score_discarded\t" << result.normalization.low_score_discarded << '\n'
              << "low_identity_discarded\t" << result.normalization.low_identity_discarded << '\n'
              << "exact_duplicates\t" << result.normalization.exact_duplicates << '\n'
              << "trimmed_records\t" << result.normalization.trimmed_records << '\n'
              << "overlap_discarded\t" << result.normalization.overlap_discarded << '\n'
              << "short_discarded\t" << result.normalization.short_discarded << '\n'
              << "canonical_records\t" << result.normalization.canonical_records << '\n'
              << "alignment_weighted_divergence\t" << std::setprecision(17)
              << result.weighted_divergence << '\n';
        writeAtomicText(query.directory / "normalization.tsv", stats.str());
        atomicPublish(raw_partial, raw);
        atomicPublish(canonical_partial, canonical);
        result.success = true;
    } catch (const std::exception& error) {
        result.error = error.what();
    }
    return result;
}

}  // namespace

namespace Mm2plusRouterDetail {

bool passesGraphQualityFilters(const ParsedPafRecord& record) {
    if (!record.alignment_score) return false;
    if (*record.alignment_score <= 0 || record.block_length == 0) return false;
    return static_cast<double>(record.matches) /
        static_cast<double>(record.block_length) >= 0.5;
}

std::string validateMm2plusVersion(
    std::string_view version_output,
    std::string_view help_output) {
    const std::string version = firstLine(version_output);
    if (version != "1.3") {
        throw std::runtime_error(
            "RaMAx requires Bioconda mm2plus 1.3 (upstream Minimap2 2.31-r1302); observed: " +
            version);
    }
    if (help_output.find("Usage: mm2plus") == std::string_view::npos) {
        throw std::runtime_error(
            "The configured executable reports the minimap2 version but is not mm2plus");
    }
    return "mm2plus " + version + "; upstream Minimap2 2.31-r1302";
}

std::vector<std::string> indexArguments(
    uint_t threads,
    uint64_t reference_bases,
    const std::filesystem::path& index,
    const std::filesystem::path& reference) {
    if (reference_bases == std::numeric_limits<uint64_t>::max()) {
        throw std::runtime_error("Reference length overflows mm2-plus -I");
    }
    return {"-x", "asm20", "-I", std::to_string(reference_bases + 1),
            "-t", std::to_string(std::max<uint_t>(1, threads)),
            "-d", index.string(), reference.string()};
}

std::vector<std::string> alignmentArguments(
    uint_t threads,
    const std::filesystem::path& index,
    const std::filesystem::path& query) {
    return {"-x", "asm20", "-c", "--eqx", "--secondary=no",
            "-t", std::to_string(std::max<uint_t>(1, threads)),
            index.string(), query.string()};
}

uint_t threadsPerPair(uint_t total_threads) {
    return std::max<uint_t>(1, std::min<uint_t>(16, total_threads));
}

size_t workerCount(size_t tasks, uint_t total_threads) {
    if (tasks == 0) return 0;
    const uint_t per_pair = threadsPerPair(total_threads);
    return std::min<size_t>(tasks,
        std::max<uint_t>(1, total_threads / per_pair));
}

NormalizationStats normalizeForGraph(
    std::vector<ParsedPafRecord>& records,
    uint_t min_span) {
    NormalizationStats stats;
    stats.raw_records = records.size();

    std::sort(records.begin(), records.end(), [](const auto& left, const auto& right) {
        return std::tie(left.query_name, left.query_start, left.query_end,
                        left.target_name, left.target_start, left.target_end,
                        left.strand, left.cigar_text) <
               std::tie(right.query_name, right.query_start, right.query_end,
                        right.target_name, right.target_start, right.target_end,
                        right.strand, right.cigar_text);
    });
    const auto unique_end = std::unique(records.begin(), records.end(),
        [](const auto& left, const auto& right) {
            return left.query_name == right.query_name &&
                   left.query_start == right.query_start &&
                   left.query_end == right.query_end &&
                   left.target_name == right.target_name &&
                   left.target_start == right.target_start &&
                   left.target_end == right.target_end &&
                   left.strand == right.strand &&
                   left.cigar_text == right.cigar_text;
        });
    stats.exact_duplicates = static_cast<size_t>(records.end() - unique_end);
    records.erase(unique_end, records.end());

    std::sort(records.begin(), records.end(), [](const auto& left, const auto& right) {
        const uint64_t left_span = left.query_end - left.query_start +
                                   left.target_end - left.target_start;
        const uint64_t right_span = right.query_end - right.query_start +
                                    right.target_end - right.target_start;
        if (*left.alignment_score != *right.alignment_score) {
            return *left.alignment_score > *right.alignment_score;
        }
        if (left.matches != right.matches) return left.matches > right.matches;
        if (left.mapq != right.mapq) return left.mapq > right.mapq;
        if (left_span != right_span) return left_span > right_span;
        return std::tie(left.query_name, left.query_start,
                        left.target_name, left.target_start) <
               std::tie(right.query_name, right.query_start,
                        right.target_name, right.target_start);
    });
    const auto overlap = WfmashRouterDetail::normalizePafForGraph(records);
    stats.trimmed_records = overlap.trimmed_records;
    stats.overlap_discarded = overlap.skipped_records;

    const auto short_end = std::remove_if(records.begin(), records.end(),
        [&](const auto& record) {
            return record.query_end - record.query_start < min_span ||
                   record.target_end - record.target_start < min_span;
        });
    stats.short_discarded = static_cast<size_t>(records.end() - short_end);
    records.erase(short_end, records.end());
    stats.canonical_records = records.size();
    return stats;
}

}  // namespace Mm2plusRouterDetail

FirstRoundMm2plusRouter::FirstRoundMm2plusRouter(
    std::filesystem::path executable,
    std::filesystem::path output_directory,
    uint_t threads)
    : executable_(std::move(executable)),
      output_directory_(std::move(output_directory)),
      threads_(std::max<uint_t>(1, threads)) {
    if (!RaMAxExternalTool::isExecutable(executable_)) {
        throw std::runtime_error("mm2plus is required but was not found");
    }
    std::filesystem::create_directories(output_directory_);
    const auto version_result = RaMAxExternalTool::run(
        executable_, {"--version"},
        output_directory_ / "mm2plus.version.stdout",
        output_directory_ / "mm2plus.version.stderr");
    if (version_result.exit_code != 0) {
        throw std::runtime_error("mm2plus --version failed");
    }
    const auto help_result = RaMAxExternalTool::run(
        executable_, {"--help"},
        output_directory_ / "mm2plus.help.stdout",
        output_directory_ / "mm2plus.help.stderr");
    (void)help_result;
    version_ = Mm2plusRouterDetail::validateMm2plusVersion(
        RaMAxExternalTool::readText(output_directory_ / "mm2plus.version.stdout") +
            RaMAxExternalTool::readText(output_directory_ / "mm2plus.version.stderr"),
        RaMAxExternalTool::readText(output_directory_ / "mm2plus.help.stdout") +
            RaMAxExternalTool::readText(output_directory_ / "mm2plus.help.stderr"));
}

FirstRoundMm2plusResult FirstRoundMm2plusRouter::run(
    const SpeciesName& reference,
    const std::vector<MashDistanceRecord>& distances,
    double far_distance,
    uint_t min_span,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers) {
    FirstRoundMm2plusResult result;
    const auto reference_it = managers.find(reference);
    if (reference_it == managers.end()) {
        throw std::runtime_error("mm2-plus reference manager is absent: " + reference);
    }
    const auto reference_records = managerRecords(reference_it->second, reference);
    const auto reference_fasta = fastaPath(reference_it->second);
    const auto index_result = ensureIndex(
        executable_, version_, threads_, output_directory_,
        reference_fasta, reference_records);
    const auto index = output_directory_ / "reference.asm20.mmi";

    std::vector<PreparedQuery> queries;
    for (const auto& distance : distances) {
        if (distance.reference != reference || distance.distance <= far_distance) continue;
        const auto manager = managers.find(distance.query);
        if (manager == managers.end()) {
            throw std::runtime_error("mm2-plus query manager is absent: " + distance.query);
        }
        PreparedQuery query;
        query.species = distance.query;
        query.distance = distance.distance;
        query.manager = manager->second;
        query.fasta = fastaPath(manager->second);
        query.records = managerRecords(manager->second, query.species);
        query.bases = totalBases(query.records);
        query.directory = output_directory_ / safeName(query.species);
        queries.push_back(std::move(query));
    }
    std::sort(queries.begin(), queries.end(), [](const auto& left, const auto& right) {
        return left.species < right.species;
    });

    const uint_t pair_threads = Mm2plusRouterDetail::threadsPerPair(threads_);
    const size_t workers = Mm2plusRouterDetail::workerCount(queries.size(), threads_);
    spdlog::info(
        "[mm2plus-router] reference={} far={} candidates={} workers={} "
        "threads_per_pair={} params={}",
        reference, far_distance, queries.size(), workers, pair_threads,
        kMm2plusParameters);

    std::vector<PairResult> pair_results(queries.size());
    if (!queries.empty()) {
#pragma omp parallel for schedule(dynamic) num_threads(static_cast<int>(workers))
        for (size_t i = 0; i < queries.size(); ++i) {
            pair_results[i] = executePair(
                executable_, pair_threads, min_span, reference, index,
                reference_records, reference_it->second, queries[i]);
        }
    }

    std::ostringstream routing;
    routing << "species\tdistance\tbackend\tstatus\tdetail\n";
    std::ostringstream performance;
    performance << "stage\tspecies\tthreads\tquery_bp\texit_code\twall_seconds\t"
                   "user_seconds\tsystem_seconds\tpeak_rss_kb\traw_records\t"
                   "canonical_records\n";
    performance << "index\t" << reference << '\t' << std::min<uint_t>(3, threads_)
                << "\t0\t" << index_result.exit_code << '\t'
                << index_result.wall_seconds << '\t' << index_result.user_seconds
                << '\t' << index_result.system_seconds << '\t'
                << index_result.peak_rss_kb << "\t0\t0\n";
    for (size_t i = 0; i < queries.size(); ++i) {
        const auto& query = queries[i];
        auto& pair = pair_results[i];
        performance << "align\t" << query.species << '\t' << pair_threads << '\t'
                    << query.bases << '\t' << pair.command.exit_code << '\t'
                    << pair.command.wall_seconds << '\t' << pair.command.user_seconds
                    << '\t' << pair.command.system_seconds << '\t'
                    << pair.command.peak_rss_kb << '\t'
                    << pair.normalization.raw_records << '\t'
                    << pair.normalization.canonical_records << '\n';
        if (pair.success) {
            result.successful_species.emplace(query.species);
            result.anchors_by_species.emplace(
                query.species, std::move(pair.anchors));
            routing << query.species << '\t' << std::setprecision(17)
                    << query.distance << "\tmm2plus\tsuccess\t"
                    << result.anchors_by_species.at(query.species).size()
                    << " anchors; divergence=" << pair.weighted_divergence << '\n';
            if (query.distance > 0.05 || pair.weighted_divergence > 0.05) {
                spdlog::warn(
                    "[mm2plus-router] {} is outside the approximately 5% "
                    "validated divergence range (Mash d={}, alignment divergence={})",
                    query.species, query.distance, pair.weighted_divergence);
            }
        } else {
            routing << query.species << '\t' << std::setprecision(17)
                    << query.distance << "\tlegacy\tfallback\t"
                    << pair.error << '\n';
            spdlog::warn("[mm2plus-router] {} falls back to legacy: {}",
                         query.species, pair.error);
        }
    }
    writeAtomicText(output_directory_ / "routing.tsv", routing.str());
    writeAtomicText(output_directory_ / "performance.tsv", performance.str());
    std::ostringstream tools;
    tools << "tool\tpath\tversion_or_parameters\n"
          << "mm2plus\t" << executable_.string() << '\t' << version_ << '\n'
          << "alignment\t.\t" << kMm2plusParameters << '\n';
    writeAtomicText(output_directory_ / "tools.tsv", tools.str());
    return result;
}

std::filesystem::path locateMm2plusExecutable() {
    const auto executable = RaMAxExternalTool::locateExecutable(
        "mm2plus", RAMAX_MM2PLUS_CONFIGURED_PATH);
    if (executable.empty()) {
        throw std::runtime_error(
            "mm2plus 1.3 was not found at the configured path, "
            "next to ramax, or in PATH");
    }
    return executable;
}
