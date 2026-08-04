#include "softmask_index.h"

#include "halAlignmentInstance.h"
#include "halGenome.h"
#include "halSequence.h"

#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

bool hasUpperAndLower(const std::string& dna) {
    bool uppercase = false;
    bool lowercase = false;
    for (char base : dna) {
        const auto byte = static_cast<unsigned char>(base);
        uppercase = uppercase || std::isupper(byte) != 0;
        lowercase = lowercase || std::islower(byte) != 0;
    }
    return uppercase && lowercase;
}

std::string readSingleSequence(const std::filesystem::path& path) {
    std::ifstream input(path);
    std::string header;
    std::string dna;
    std::getline(input, header);
    std::getline(input, dna);
    return dna;
}

void writeInput(const std::filesystem::path& path, const std::string& dna) {
    std::ofstream output(path);
    output << ">chr1\n" << dna << '\n';
}

std::string readGenome(hal::AlignmentPtr alignment, const std::string& name) {
    hal::Genome* genome = alignment->openGenome(name);
    if (!genome) throw std::runtime_error("missing genome in HAL: " + name);
    std::string dna;
    genome->getString(dna);
    alignment->closeGenome(genome);
    return dna;
}

}  // namespace

int main() {
    const auto temp = std::filesystem::path("/tmp") /
        ("ramax-hal-softmask-roundtrip-" + std::to_string(getpid()));
    std::filesystem::remove_all(temp);
    std::filesystem::create_directories(temp);

    try {
        const std::string leaf_a_expected = "AAAcccGGttTT";
        const std::string leaf_b_expected = "AAAcccGGTTtt";
        const std::string ancestor_expected = "AAAcccGGtttt";

        SoftMask::PathMap paths;
        for (const auto& [species, dna] :
             std::vector<std::pair<std::string, std::string>>{
                 {"leafA", leaf_a_expected}, {"leafB", leaf_b_expected}}) {
            const auto input = temp / (species + ".fa");
            const auto uppercase = temp / (species + ".align-v2.fasta");
            const auto index = temp / (species + ".softmask-v1.bin");
            const auto marker = temp / (species + ".complete.json");
            writeInput(input, dna);
            SoftMask::ensureUppercaseFastaAndIndex(input, uppercase, index, marker);
            paths[species] = index;
        }

        const SoftMask::IndexMap indexes = SoftMask::loadIndexes(paths);
        std::string leaf_a = "AAACCCGGTTTT";
        std::string leaf_b = "AAACCCGGTTTT";
        indexes.at("leafA")->restore("chr1", 0, leaf_a);
        indexes.at("leafB")->restore("chr1", 0, leaf_b);
        require(leaf_a == leaf_a_expected, "leafA restoration differs from input mask");
        require(leaf_b == leaf_b_expected, "leafB restoration differs from input mask");

        std::string ancestor;
        ancestor.reserve(leaf_a.size());
        for (size_t position = 0; position < leaf_a.size(); ++position) {
            SoftMask::AncestorBaseVote vote;
            vote.add(leaf_a[position]);
            vote.add(leaf_b[position]);
            ancestor.push_back(vote.result());
        }
        require(ancestor == ancestor_expected, "ancestor lowercase-majority vote is incorrect");

        const auto hal_path = temp / "case-roundtrip.hal";
        {
            hal::AlignmentPtr alignment = hal::openHalAlignment(
                hal_path.string(), nullptr, hal::CREATE_ACCESS);
            require(static_cast<bool>(alignment), "cannot create HAL");
            hal::Genome* root = alignment->addRootGenome("anc0");
            hal::Genome* leaf_a_genome = alignment->addLeafGenome("leafA", "anc0", 0.1);
            hal::Genome* leaf_b_genome = alignment->addLeafGenome("leafB", "anc0", 0.1);

            root->setDimensions({hal::Sequence::Info("anc0.chr1", ancestor.size(), 0, 0)});
            leaf_a_genome->setDimensions({hal::Sequence::Info("chr1", leaf_a.size(), 0, 0)});
            leaf_b_genome->setDimensions({hal::Sequence::Info("chr1", leaf_b.size(), 0, 0)});
            root->getSequence("anc0.chr1")->setString(ancestor);
            leaf_a_genome->getSequence("chr1")->setString(leaf_a);
            leaf_b_genome->getSequence("chr1")->setString(leaf_b);

            alignment->closeGenome(leaf_b_genome);
            alignment->closeGenome(leaf_a_genome);
            alignment->closeGenome(root);
            alignment.reset();
        }

        hal::AlignmentPtr reopened = hal::openHalAlignment(
            hal_path.string(), nullptr, hal::READ_ACCESS);
        require(static_cast<bool>(reopened), "cannot reopen HAL");
        const std::string hal_leaf_a = readGenome(reopened, "leafA");
        const std::string hal_leaf_b = readGenome(reopened, "leafB");
        const std::string hal_ancestor = readGenome(reopened, "anc0");
        require(hal_leaf_a == leaf_a_expected, "HAL leafA did not preserve case");
        require(hal_leaf_b == leaf_b_expected, "HAL leafB did not preserve case");
        require(hal_ancestor == ancestor_expected, "HAL ancestor did not preserve reconstructed case");
        require(hasUpperAndLower(hal_leaf_a), "HAL leafA lacks mixed case");
        require(hasUpperAndLower(hal_leaf_b), "HAL leafB lacks mixed case");
        require(hasUpperAndLower(hal_ancestor), "HAL ancestor lacks mixed case");
        reopened.reset();

        std::cout << "genome\tdna\n"
                  << "leafA\t" << hal_leaf_a << '\n'
                  << "leafB\t" << hal_leaf_b << '\n'
                  << "anc0\t" << hal_ancestor << '\n';
        std::filesystem::remove_all(temp);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "hal_softmask_roundtrip_test: " << error.what() << '\n';
        std::filesystem::remove_all(temp);
        return 1;
    }
}
