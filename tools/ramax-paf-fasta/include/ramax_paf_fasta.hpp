#ifndef RAMAX_PAF_FASTA_HPP
#define RAMAX_PAF_FASTA_HPP

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace RamaxPafFasta {

struct SeqfileEntry {
    std::string species;
    std::string source;
};

struct Options {
    std::filesystem::path seqfile;
    std::filesystem::path output;
    bool force = false;
};

struct Stats {
    std::uint64_t species = 0;
    std::uint64_t contigs = 0;
    std::uint64_t bases = 0;
    std::uint64_t converted_to_n = 0;
    std::uint64_t output_bytes = 0;
    double elapsed_seconds = 0.0;
    bool gzip_output = false;
};

std::vector<SeqfileEntry> parseSeqfile(const std::filesystem::path& path);
Stats generate(const Options& options);

}  // namespace RamaxPafFasta

#endif
