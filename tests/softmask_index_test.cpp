#include "softmask_index.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unistd.h>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

std::string readFile(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(input),
                       std::istreambuf_iterator<char>());
}

}  // namespace

int main() {
    std::filesystem::path temp = std::filesystem::path("/tmp") /
        ("ramax-softmask-index-test-" + std::to_string(getpid()));
    std::filesystem::remove_all(temp);
    std::filesystem::create_directories(temp);

    try {
        const auto input = temp / "input.fa";
        const auto uppercase = temp / "clean.align-v2.fasta";
        const auto index_path = temp / "clean.softmask-v1.bin";
        const auto marker = temp / "clean.softmask-v1.complete.json";
        {
            std::ofstream output(input);
            output << ">chr1\nAAAccccGGrrT\n"
                   << ">chr2\nnNACgt\n";
        }

        SoftMask::ensureUppercaseFastaAndIndex(input, uppercase, index_path, marker);
        require(readFile(uppercase) == ">chr1\nAAACCCCGGNNT\n>chr2\nNNACGT\n",
                "uppercase alignment FASTA changed semantics");
        require(std::filesystem::is_regular_file(marker), "completion marker missing");

        SoftMask::Index index(index_path);
        require(index.sequenceLength("chr1") == 12, "chr1 length mismatch");
        require(index.intervalCount("chr1") == 2, "chr1 interval count mismatch");
        require(index.intervalCount("chr2") == 2, "chr2 interval count mismatch");

        std::string chr1 = "AAACCCCGGNNT";
        index.restore("chr1", 0, chr1);
        require(chr1 == "AAAccccGGnnT", "whole-chromosome mask restore failed");

        std::string slice = "CCCGGNN";  // chr1 [4, 11)
        index.restore("chr1", 4, slice);
        require(slice == "cccGGnn", "partial-range mask restore failed");

        std::string chr2 = "NNACGT";
        index.restore("chr2", 0, chr2);
        require(chr2 == "nNACgt", "N/n or terminal interval restore failed");

        std::filesystem::remove_all(temp);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "softmask_index_test: " << error.what() << '\n';
        std::filesystem::remove_all(temp);
        return 1;
    }
}
