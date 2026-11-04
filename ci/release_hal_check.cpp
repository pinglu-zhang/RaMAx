// Independent reader for release checks; it is not installed in the package.
#include "halAlignmentInstance.h"
#include "halGenome.h"
#include "halSequence.h"
#include "halValidate.h"
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    try {
        if (argc != 2 && argc != 3)
            throw std::runtime_error("usage: halValidate ALIGNMENT.hal [EXPECTED.tsv]");
        auto alignment = hal::openHalAlignment(argv[1], nullptr, hal::READ_ACCESS);
        if (!alignment) throw std::runtime_error("cannot open HAL");
        hal::validateAlignment(alignment.get());
        if (argc == 3) {
            // Expected rows are sample<TAB>contig<TAB>complete original DNA.
            std::ifstream expected(argv[2]);
            if (!expected) throw std::runtime_error("cannot open expected DNA");
            std::map<std::string, size_t> counts;
            std::string line;
            while (std::getline(expected, line)) {
                auto first = line.find('\t'), second = line.find('\t', first + 1);
                if (first == std::string::npos || second == std::string::npos)
                    throw std::runtime_error("invalid expected DNA row");
                const auto sample = line.substr(0, first);
                const auto contig = line.substr(first + 1, second - first - 1);
                const auto dna = line.substr(second + 1);
                const auto* genome = alignment->openGenome(sample);
                if (!genome) throw std::runtime_error("HAL missing sample " + sample);
                const auto* sequence = genome->getSequence(contig);
                if (!sequence) throw std::runtime_error("HAL missing contig " + contig);
                std::string observed;
                sequence->getString(observed);
                if (observed != dna)
                    throw std::runtime_error("HAL DNA mismatch: " + sample + "/" + contig);
                ++counts[sample];
                alignment->closeGenome(genome);
            }
            if (expected.bad() || counts.empty())
                throw std::runtime_error("empty or unreadable expected DNA");
            for (const auto& [sample, count] : counts) {
                const auto* genome = alignment->openGenome(sample);
                if (genome->getNumSequences() != count)
                    throw std::runtime_error("HAL contig count mismatch: " + sample);
                alignment->closeGenome(genome);
            }
        }
        std::cout << "HAL structure and requested leaf DNA checks passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "release_hal_check: " << e.what() << '\n';
        return 1;
    }
}
