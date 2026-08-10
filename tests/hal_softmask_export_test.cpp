#include "ramesh.h"

#include "SeqPro.h"
#include "align.h"
#include "anchor.h"
#include "halAlignmentInstance.h"
#include "halGenome.h"

#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <unistd.h>

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

void writeInput(const std::filesystem::path& path, const std::string& dna) {
    std::ofstream output(path);
    output << ">chr1\n" << dna << '\n';
    if (!output) throw std::runtime_error("cannot write test FASTA");
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
        ("ramax-hal-softmask-export-" + std::to_string(getpid()));
    std::filesystem::remove_all(temp);
    std::filesystem::create_directories(temp);

    try {
        const std::string leaf_a_expected = "AAAcccGGttTT";
        const std::string leaf_b_expected = "AAAcccGGTTtt";
        const std::string ancestor_expected = "AAAcccGGtttt";

        SoftMask::PathMap softmask_paths;
        std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
        for (const auto& [species, dna] :
             std::map<std::string, std::string>{{"leafA", leaf_a_expected},
                                                 {"leafB", leaf_b_expected}}) {
            const auto input = temp / (species + ".input.fa");
            const auto uppercase = temp / (species + ".align-v2.fasta");
            const auto index = temp / (species + ".softmask-v1.bin");
            const auto marker = temp / (species + ".complete.json");
            writeInput(input, dna);
            SoftMask::ensureUppercaseFastaAndIndex(input, uppercase, index, marker);
            softmask_paths[species] = index;

            SeqPro::ManagerVariant manager =
                std::make_unique<SeqPro::SequenceManager>(uppercase);
            managers[species] =
                std::make_shared<SeqPro::ManagerVariant>(std::move(manager));
        }

        RaMesh::RaMeshMultiGenomeGraph graph(managers);
        constexpr uint32_t length = 12;
        Anchor anchor(0, 0, length, 0, 0, length, Strand::FORWARD,
                      length, length, Cigar_t{cigarToInt('M', length)});
        graph.insertAnchorIntoGraph(*managers.at("leafA"), *managers.at("leafB"),
                                    "leafA", "leafB", anchor);

        const auto hal_path = temp / "exported.hal";
        graph.exportToHal(hal_path, managers,
                          "(leafA:0.1,leafB:0.1)anc0;", true, "anc0",
                          softmask_paths);

        hal::AlignmentPtr alignment = hal::openHalAlignment(
            hal_path.string(), nullptr, hal::READ_ACCESS);
        require(static_cast<bool>(alignment), "cannot reopen exported HAL");
        const std::string leaf_a = readGenome(alignment, "leafA");
        const std::string leaf_b = readGenome(alignment, "leafB");
        const std::string ancestor = readGenome(alignment, "anc0");

        require(leaf_a == leaf_a_expected, "exported leafA mask differs from input");
        require(leaf_b == leaf_b_expected, "exported leafB mask differs from input");
        require(ancestor == ancestor_expected,
                "exported ancestor mask differs from lowercase vote");
        require(hasUpperAndLower(leaf_a), "exported leafA lacks mixed case");
        require(hasUpperAndLower(leaf_b), "exported leafB lacks mixed case");
        require(hasUpperAndLower(ancestor), "exported ancestor lacks mixed case");
        alignment.reset();

        std::cout << "genome\tdna\n"
                  << "leafA\t" << leaf_a << '\n'
                  << "leafB\t" << leaf_b << '\n'
                  << "anc0\t" << ancestor << '\n';
        std::filesystem::remove_all(temp);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "hal_softmask_export_test: " << error.what() << '\n';
        std::filesystem::remove_all(temp);
        return 1;
    }
}
