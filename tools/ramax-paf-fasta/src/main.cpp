#include "ramax_paf_fasta.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#ifndef RAMAX_PAF_FASTA_VERSION
#define RAMAX_PAF_FASTA_VERSION "unknown"
#endif

namespace {

void printUsage(std::ostream& output) {
    output
        << "ramax-paf-fasta: build the qualified FASTA paired with RaMAx PAF\n\n"
        << "Usage:\n"
        << "  ramax-paf-fasta -i <seqfile> -o <output.fa[.gz]> [--force]\n\n"
        << "Options:\n"
        << "  -i, --input <path>   RaMAx seqfile (optional Newick first record)\n"
        << "  -o, --output <path>  .fa/.fasta/.fna, optionally gzip-compressed\n"
        << "      --force          Replace an existing output after success\n"
        << "  -h, --help           Show this help\n"
        << "  -v, --version        Show the tool version\n";
}

RamaxPafFasta::Options parseArguments(int argc, char** argv) {
    RamaxPafFasta::Options options;
    bool have_input = false;
    bool have_output = false;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "-h" || argument == "--help") {
            printUsage(std::cout);
            std::exit(EXIT_SUCCESS);
        }
        if (argument == "-v" || argument == "--version") {
            std::cout << "ramax-paf-fasta " << RAMAX_PAF_FASTA_VERSION << '\n';
            std::exit(EXIT_SUCCESS);
        }
        if (argument == "--force") {
            options.force = true;
            continue;
        }
        if (argument == "-i" || argument == "--input") {
            if (have_input || index + 1 >= argc) {
                throw std::runtime_error("--input requires exactly one path");
            }
            options.seqfile = argv[++index];
            have_input = true;
            continue;
        }
        if (argument == "-o" || argument == "--output") {
            if (have_output || index + 1 >= argc) {
                throw std::runtime_error("--output requires exactly one path");
            }
            options.output = argv[++index];
            have_output = true;
            continue;
        }
        throw std::runtime_error("unknown option: " + argument);
    }
    if (!have_input) throw std::runtime_error("missing required --input");
    if (!have_output) throw std::runtime_error("missing required --output");
    return options;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const auto options = parseArguments(argc, argv);
        const auto stats = RamaxPafFasta::generate(options);
        std::cout
            << "ramax-paf-fasta complete: species=" << stats.species
            << ", contigs=" << stats.contigs
            << ", bases=" << stats.bases
            << ", converted_to_n=" << stats.converted_to_n
            << ", format=" << (stats.gzip_output ? "gzip" : "plain")
            << ", output_bytes=" << stats.output_bytes
            << ", elapsed_seconds=" << stats.elapsed_seconds << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "ramax-paf-fasta: error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
