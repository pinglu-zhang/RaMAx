#include "mash_distance_estimator.h"

#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void expectThrow(const std::function<void()>& action) {
    bool threw = false;
    try {
        action();
    } catch (const std::exception&) {
        threw = true;
    }
    require(threw, "expected exception was not thrown");
}

void writeFasta(const std::filesystem::path& path,
                const std::string& name,
                std::string_view motif) {
    std::ofstream output(path);
    output << '>' << name << '\n';
    for (size_t i = 0; i < 2000; ++i) output << motif;
    output << '\n';
}

SeqPro::SharedManagerVariant makeManager(
    const std::filesystem::path& path) {
    auto original = std::make_unique<SeqPro::SequenceManager>(path);
    auto masked =
        std::make_unique<SeqPro::MaskedSequenceManager>(std::move(original));
    return std::make_shared<SeqPro::ManagerVariant>(std::move(masked));
}

}  // namespace

int main() {
    const auto parsed = MashDistanceDetail::parseLine(
        "reference.fa\tquery.fa\t0.0125\t1e-20\t15000/20000");
    require(parsed.reference_id == "reference.fa", "reference ID mismatch");
    require(parsed.query_id == "query.fa", "query ID mismatch");
    require(parsed.distance == 0.0125, "distance mismatch");
    require(parsed.shared_hashes == 15000, "shared hash mismatch");
    require(parsed.total_hashes == 20000, "total hash mismatch");
    require(MashDistanceDetail::validateVersion("2.3\n") == "2.3",
            "version mismatch");

    expectThrow([] {
        MashDistanceDetail::parseLine("too\tfew\tfields");
    });
    expectThrow([] {
        MashDistanceDetail::parseLine("r\tq\tnan\t0\t1/2");
    });
    expectThrow([] {
        MashDistanceDetail::parseLine("r\tq\t0.1\t0\t0/0");
    });
    expectThrow([] {
        MashDistanceDetail::parseLine("r\tq\t0.1\t0\t1x/2");
    });
    expectThrow([] {
        MashDistanceDetail::validateVersion("2.2");
    });
    expectThrow([] {
        MashDistanceEstimator missing(
            "/definitely/missing/mash", "missing-output", 1);
    });

    const auto mash = locateMashExecutable();
    if (mash.empty()) {
        std::cerr << "Mash 2.3 is not available for integration test\n";
        return 1;
    }

    const auto root = std::filesystem::current_path() /
        "mash estimator test tmp";
    std::filesystem::remove_all(root);
    std::filesystem::create_directories(root);
    const auto reference_fasta = root / "reference genome.fa";
    const auto query_fasta = root / "query genome.fa";
    writeFasta(reference_fasta, "reference", "ACGTGCAATGCC");
    writeFasta(query_fasta, "query", "ACGTGCAATGCT");

    std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
    managers.emplace("reference", makeManager(reference_fasta));
    managers.emplace("query_a", makeManager(query_fasta));
    managers.emplace("query_duplicate", makeManager(query_fasta));

    MashDistanceEstimator estimator(mash, root / "similarity output", 2);
    const auto records =
        estimator.estimateFirstReference("reference", managers);
    require(records.size() == 2, "duplicate paths were not expanded");
    require(records[0].query == "query_a", "first query mismatch");
    require(records[1].query == "query_duplicate", "second query mismatch");
    require(records[0].distance >= 0.0 && records[0].distance <= 1.0,
            "distance outside unit interval");
    require(records[0].distance == records[1].distance,
            "duplicate paths produced different distances");
    require(std::filesystem::is_regular_file(
        root / "similarity output" / "mash_first_reference.tsv"),
        "distance table was not written");

    std::filesystem::remove_all(root);
    std::cout << "Mash distance estimator tests passed\n";
    return 0;
}
