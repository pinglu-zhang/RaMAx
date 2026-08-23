#include "ramax_paf_fasta.hpp"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <zlib.h>

namespace {

namespace fs = std::filesystem;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void writeFile(const fs::path& path, const std::string& content) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("cannot create fixture: " + path.string());
    output.write(content.data(), static_cast<std::streamsize>(content.size()));
    if (!output) throw std::runtime_error("cannot write fixture: " + path.string());
}

std::string readFile(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot read fixture: " + path.string());
    return std::string(std::istreambuf_iterator<char>(input),
                       std::istreambuf_iterator<char>());
}

void writeGzip(const fs::path& path, const std::string& content) {
    gzFile output = gzopen(path.c_str(), "wb6");
    if (!output) throw std::runtime_error("cannot create gzip fixture");
    const int written = gzwrite(output, content.data(),
                                static_cast<unsigned int>(content.size()));
    const int closed = gzclose(output);
    if (written != static_cast<int>(content.size()) || closed != Z_OK) {
        throw std::runtime_error("cannot write gzip fixture");
    }
}

std::string readGzip(const fs::path& path) {
    gzFile input = gzopen(path.c_str(), "rb");
    if (!input) throw std::runtime_error("cannot open gzip result");
    std::string result;
    char buffer[4096];
    for (;;) {
        const int count = gzread(input, buffer, sizeof(buffer));
        if (count < 0) {
            gzclose(input);
            throw std::runtime_error("cannot read gzip result");
        }
        if (count == 0) break;
        result.append(buffer, static_cast<std::size_t>(count));
    }
    if (gzclose(input) != Z_OK) throw std::runtime_error("cannot close gzip result");
    return result;
}

template<class Function>
void requireThrows(Function&& function, const std::string& fragment) {
    try {
        function();
    } catch (const std::exception& error) {
        require(std::string(error.what()).find(fragment) != std::string::npos,
                "unexpected error: " + std::string(error.what()));
        return;
    }
    throw std::runtime_error("expected failure containing: " + fragment);
}

void testPlainAndGzip(const fs::path& directory) {
    const fs::path plain_input = directory / "zeta genome.fa";
    const fs::path gzip_input = directory / "alpha.fa.gz";
    writeFile(plain_input,
              ">chr1 description\r\nacgtRyn\n"
              ">zeta.chr2 extra\nTT\nGG\n");
    writeGzip(gzip_input, ">ctg note\nnNa-x\n");

    const fs::path seqfile = directory / "genomes.seqfile";
    writeFile(seqfile,
              "((zeta:0.1,alpha:0.1):0.1);\n"
              "zeta " + plain_input.string() + "\n"
              "alpha " + gzip_input.string() + "\n");

    const auto entries = RamaxPafFasta::parseSeqfile(seqfile);
    require(entries.size() == 2, "ordered seqfile entry count");
    require(entries[0].species == "zeta" && entries[1].species == "alpha",
            "seqfile species order was not preserved");
    require(entries[0].source == plain_input.string(),
            "path containing spaces was not preserved");

    const std::string expected =
        ">zeta.chr1\nACGTNNN\n"
        ">zeta.zeta.chr2\nTTGG\n"
        ">alpha.ctg\nNNANN\n";
    const fs::path plain_output = directory / "combined.fa";
    const auto plain_stats = RamaxPafFasta::generate(
        {seqfile, plain_output, false});
    require(readFile(plain_output) == expected, "plain normalized FASTA");
    require(plain_stats.species == 2 && plain_stats.contigs == 3,
            "plain species/contig stats");
    require(plain_stats.bases == 16, "plain base stats");
    require(plain_stats.converted_to_n == 4, "plain conversion stats");
    require(!plain_stats.gzip_output, "plain output mode");

    const fs::path gzip_output = directory / "combined.fasta.gz";
    const auto gzip_stats = RamaxPafFasta::generate(
        {seqfile, gzip_output, false});
    require(readGzip(gzip_output) == expected,
            "gzip and plain decompressed content differ");
    require(gzip_stats.gzip_output, "gzip output mode");

    writeFile(plain_output, "sentinel\n");
    requireThrows(
        [&]() {
            (void)RamaxPafFasta::generate({seqfile, plain_output, false});
        },
        "output already exists");
    require(readFile(plain_output) == "sentinel\n",
            "default overwrite changed existing output");
    (void)RamaxPafFasta::generate({seqfile, plain_output, true});
    require(readFile(plain_output) == expected,
            "--force did not publish complete output");
}

void testSeqfileAndFastaFailures(const fs::path& directory) {
    const fs::path input = directory / "valid.fa";
    writeFile(input, ">c\nACGT\n");

    const fs::path duplicate_species = directory / "duplicate.seqfile";
    writeFile(duplicate_species,
              "sample " + input.string() + "\n"
              "sample " + input.string() + "\n");
    requireThrows(
        [&]() { (void)RamaxPafFasta::parseSeqfile(duplicate_species); },
        "duplicate species");

    const fs::path empty_record = directory / "empty.fa";
    writeFile(empty_record, ">empty\n>nonempty\nACGT\n");
    const fs::path empty_seqfile = directory / "empty.seqfile";
    writeFile(empty_seqfile, "sample " + empty_record.string() + "\n");
    requireThrows(
        [&]() {
            (void)RamaxPafFasta::generate(
                {empty_seqfile, directory / "empty-output.fa", false});
        },
        "no sequence bases");

    requireThrows(
        [&]() {
            (void)RamaxPafFasta::generate(
                {empty_seqfile, directory / "bad-extension.txt", false});
        },
        "output must end");

    const fs::path same_path_seqfile = directory / "same-path.seqfile";
    writeFile(same_path_seqfile, "sample " + input.string() + "\n");
    requireThrows(
        [&]() {
            (void)RamaxPafFasta::generate({same_path_seqfile, input, true});
        },
        "also an input FASTA");
    require(readFile(input) == ">c\nACGT\n",
            "same input/output check modified the source FASTA");

    const fs::path corrupt_gzip = directory / "corrupt.fa.gz";
    writeFile(corrupt_gzip,
              std::string("\x1f\x8b\x08\x00", 4) + "broken-gzip");
    const fs::path corrupt_seqfile = directory / "corrupt.seqfile";
    writeFile(corrupt_seqfile, "sample " + corrupt_gzip.string() + "\n");
    requireThrows(
        [&]() {
            (void)RamaxPafFasta::generate(
                {corrupt_seqfile, directory / "corrupt-output.fa", false});
        },
        "failed to read FASTA");
}

void testLineWrapping(const fs::path& directory) {
    const fs::path input = directory / "long-line.fa";
    writeFile(input, ">long\n" + std::string(61, 'a') + "\n");
    const fs::path seqfile = directory / "long-line.seqfile";
    writeFile(seqfile, "sample " + input.string() + "\n");
    const fs::path output = directory / "long-line-output.fa";
    (void)RamaxPafFasta::generate({seqfile, output, false});
    require(readFile(output) ==
                ">sample.long\n" + std::string(60, 'A') + "\nA\n",
            "FASTA output is not wrapped at 60 bases");
}

void testQualifiedCollisionIsAtomic(const fs::path& directory) {
    const fs::path first = directory / "collision-first.fa";
    const fs::path second = directory / "collision-second.fa";
    writeFile(first, ">c\nACGT\n");
    writeFile(second, ">b.c\nACGT\n");
    const fs::path seqfile = directory / "collision.seqfile";
    writeFile(seqfile,
              "a.b " + first.string() + "\n"
              "a " + second.string() + "\n");
    const fs::path output = directory / "collision.fa";
    writeFile(output, "sentinel\n");
    requireThrows(
        [&]() {
            (void)RamaxPafFasta::generate({seqfile, output, true});
        },
        "duplicate qualified FASTA name");
    require(readFile(output) == "sentinel\n",
            "collision replaced existing output");
    for (const auto& entry : fs::directory_iterator(directory)) {
        require(entry.path().filename().string().find(".collision.fa.tmp.") != 0,
                "collision left a temporary output");
    }
}

}  // namespace

int main() {
    const auto stamp = std::chrono::steady_clock::now()
                           .time_since_epoch().count();
    const fs::path directory = fs::temp_directory_path() /
        ("ramax-paf-fasta-test-" + std::to_string(stamp));
    struct Cleanup {
        fs::path path;
        ~Cleanup() {
            std::error_code error;
            fs::remove_all(path, error);
        }
    } cleanup{directory};

    try {
        fs::create_directories(directory);
        testPlainAndGzip(directory);
        testSeqfileAndFastaFailures(directory);
        testLineWrapping(directory);
        testQualifiedCollisionIsAtomic(directory);
    } catch (const std::exception& error) {
        std::cerr << "ramax-paf-fasta-tests: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    std::cout << "ramax-paf-fasta-tests: all checks passed\n";
    return EXIT_SUCCESS;
}
