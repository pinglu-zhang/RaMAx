#include "structural_break_repair.h"

#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

#include "align.h"
#include "ramesh.h"

namespace {

using namespace RaMesh;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void writeFasta(const std::filesystem::path& path,
                char chr1_base,
                char chr2_base) {
    std::ofstream output(path);
    output << ">chr1\n" << std::string(1000, chr1_base) << '\n'
           << ">chr2\n" << std::string(1000, chr2_base) << '\n';
    require(static_cast<bool>(output), "cannot write structural test FASTA");
}

void writePassthroughMsa(const std::filesystem::path& path) {
    std::ofstream output(path);
    output << "#!/bin/sh\nexec /bin/cat \"$5\"\n";
    require(static_cast<bool>(output), "cannot write structural test MSA");
    output.close();
    std::filesystem::permissions(
        path,
        std::filesystem::perms::owner_read |
            std::filesystem::perms::owner_write |
            std::filesystem::perms::owner_exec,
        std::filesystem::perm_options::replace);
}

SegPtr add(const BlockPtr& block,
           const SpeciesChrPair& key,
           uint_t start,
           uint_t length,
           Strand strand = Strand::FORWARD) {
    auto segment = Segment::create(
        start, length, strand, Cigar_t{cigarToInt('M', length)},
        AlignRole::PRIMARY, SegmentRole::SEGMENT, block);
    require(block->anchors.emplace(key, segment).second,
            "duplicate structural test anchor");
    return segment;
}

void link(GenomeEnd& end, const std::vector<SegPtr>& segments) {
    auto previous = end.head;
    for (const auto& segment : segments) {
        previous->primary_path.next.store(segment);
        segment->primary_path.prev.store(previous);
        previous = segment;
    }
    previous->primary_path.next.store(end.tail);
    end.tail->primary_path.prev.store(previous);
    end.sample_vec.assign(1, end.head);
    for (const auto& segment : segments) end.setToSampling(segment);
}

struct Fixture {
    BlockPtr left;
    BlockPtr anomaly;
    BlockPtr right;
    SegPtr ref_left;
    SegPtr ref_anomaly;
    SegPtr ref_right;
    SegPtr qry_left;
    SegPtr qry_anomaly;
    SegPtr qry_right;
};

Fixture populate(RaMeshMultiGenomeGraph& graph) {
    Fixture f;
    f.left = Block::createEmpty("chr1", 2);
    f.anomaly = Block::createEmpty("chr1", 2);
    f.right = Block::createEmpty("chr1", 2);
    f.ref_left = add(f.left, {"ref", "chr1"}, 0, 100);
    f.ref_anomaly = add(f.anomaly, {"ref", "chr1"}, 100, 20);
    f.ref_right = add(f.right, {"ref", "chr1"}, 120, 100);
    f.qry_left = add(f.left, {"qry", "chr1"}, 0, 100);
    f.qry_anomaly = add(f.anomaly, {"qry", "chr2"}, 50, 20);
    f.qry_right = add(f.right, {"qry", "chr1"}, 120, 100);
    link(graph.species_graphs.at("ref").chr2end.at("chr1"),
         {f.ref_left, f.ref_anomaly, f.ref_right});
    link(graph.species_graphs.at("qry").chr2end.at("chr1"),
         {f.qry_left, f.qry_right});
    link(graph.species_graphs.at("qry").chr2end.at("chr2"),
         {f.qry_anomaly});
    graph.blocks = {f.left, f.anomaly, f.right};
    return f;
}

Fixture populateK3(RaMeshMultiGenomeGraph& graph) {
    auto fixture = populate(graph);
    const auto protected_left = add(
        fixture.left, {"protected", "chr1"}, 0, 100);
    const auto protected_right = add(
        fixture.right, {"protected", "chr1"}, 120, 100);
    link(graph.species_graphs.at("protected").chr2end.at("chr1"),
         {protected_left, protected_right});
    return fixture;
}

Fixture populateK4(RaMeshMultiGenomeGraph& graph) {
    auto fixture = populateK3(graph);
    const auto protected_left = add(
        fixture.left, {"protected2", "chr1"}, 0, 100);
    const auto protected_right = add(
        fixture.right, {"protected2", "chr1"}, 120, 100);
    link(graph.species_graphs.at("protected2").chr2end.at("chr1"),
         {protected_left, protected_right});
    return fixture;
}

size_t activeBlocks(const RaMeshMultiGenomeGraph& graph) {
    size_t count = 0;
    for (const auto& weak : graph.blocks) count += !weak.expired();
    return count;
}

}  // namespace

int main() {
    const auto temp = std::filesystem::path("/tmp") /
        ("ramax-structural-break-" + std::to_string(getpid()));
    std::filesystem::remove_all(temp);
    std::filesystem::create_directories(temp);
    try {
        std::map<SpeciesName, SeqPro::SharedManagerVariant> managers;
        for (const auto& species : {std::string("ref"), std::string("qry"),
                                    std::string("protected"),
                                    std::string("protected2")}) {
            const auto fasta = temp / (species + ".fa");
            writeFasta(fasta, 'A', 'C');
            auto manager = std::make_unique<SeqPro::SequenceManager>(fasta);
            managers.emplace(
                species,
                std::make_shared<SeqPro::ManagerVariant>(std::move(manager)));
        }
        const auto msa = temp / "passthrough-minipoa";
        writePassthroughMsa(msa);

        {
            RaMeshMultiGenomeGraph graph(managers);
            const auto fixture = populate(graph);
            StructuralBreakRepair::Options disabled;
            const auto result =
                StructuralBreakRepair::repairAnchorBoundedStructuralBreaks(
                    graph, "ref", managers, disabled);
            require(result.committed == 0 && activeBlocks(graph) == 3,
                    "disabled repair changed the graph");
            require(fixture.qry_anomaly->parent_block == fixture.anomaly,
                    "disabled repair changed anomaly ownership");
        }

        {
            RaMeshMultiGenomeGraph graph(managers);
            const auto fixture = populateK3(graph);
            StructuralBreakRepair::Options options;
            options.enabled = true;
            options.msa_executable = msa;
            const auto result =
                StructuralBreakRepair::repairAnchorBoundedStructuralBreaks(
                    graph, "ref", managers, options);
            require(result.candidates_by_k.at(3) == 1 &&
                        result.committed_by_k.at(3) == 1,
                    "K=3 target-switch window was not committed");
            const auto repaired =
                fixture.ref_left->primary_path.next.load();
            require(repaired && repaired->parent_block &&
                        repaired->parent_block->anchors.size() == 3,
                    "K=3 replacement did not retain the protected species");
            const auto protected_anchor = repaired->parent_block->anchors.find(
                {"protected", "chr1"});
            require(protected_anchor != repaired->parent_block->anchors.end() &&
                        protected_anchor->second->start == 100 &&
                        protected_anchor->second->length == 20,
                    "protected species interval was not rebuilt");
        }

        {
            RaMeshMultiGenomeGraph graph(managers);
            const auto fixture = populateK4(graph);
            StructuralBreakRepair::Options options;
            options.enabled = true;
            options.msa_executable = msa;
            const auto result =
                StructuralBreakRepair::repairAnchorBoundedStructuralBreaks(
                    graph, "ref", managers, options);
            require(result.candidates_by_k.at(4) == 1 &&
                        result.committed_by_k.at(4) == 1,
                    "K=4 target-switch window was not committed");
            const auto repaired =
                fixture.ref_left->primary_path.next.load();
            require(repaired && repaired->parent_block &&
                        repaired->parent_block->anchors.size() == 4,
                    "K=4 replacement did not retain both protected species");
        }

        {
            RaMeshMultiGenomeGraph graph(managers);
            const auto fixture = populate(graph);
            StructuralBreakRepair::Options options;
            options.enabled = true;
            options.msa_executable = msa;
            options.parallel_threads = 2;
            const auto result =
                StructuralBreakRepair::repairAnchorBoundedStructuralBreaks(
                    graph, "ref", managers, options);
            require(result.structural_candidates == 1 && result.committed == 1,
                    "target-switch window was not committed");
            require(activeBlocks(graph) == 4,
                    "repair did not retain the residual anomaly Block");
            const auto repaired =
                fixture.ref_left->primary_path.next.load();
            require(repaired && repaired != fixture.ref_anomaly &&
                        repaired->start == 100 && repaired->length == 20,
                    "reference interval was not replaced");
            const auto query_repaired =
                fixture.qry_left->primary_path.next.load();
            require(query_repaired && query_repaired != fixture.qry_right &&
                        query_repaired->start == 100 &&
                        query_repaired->length == 20 &&
                        query_repaired->primary_path.next.load() ==
                            fixture.qry_right,
                    "query interval was not inserted on the correct path");
            require(fixture.qry_anomaly->parent_block &&
                        fixture.qry_anomaly->parent_block->anchors.size() == 1 &&
                        fixture.qry_anomaly->parent_block != fixture.anomaly,
                    "wrong-location sequence was not retained as residual");
        }

        {
            RaMeshMultiGenomeGraph graph(managers);
            const auto fixture = populate(graph);
            StructuralBreakRepair::Options options;
            options.enabled = true;
            options.msa_executable = "/bin/false";
            const auto result =
                StructuralBreakRepair::repairAnchorBoundedStructuralBreaks(
                    graph, "ref", managers, options);
            require(result.msa_failed == 1 && result.committed == 0,
                    "failed minipoa was not rejected");
            require(activeBlocks(graph) == 3 &&
                        fixture.qry_anomaly->parent_block == fixture.anomaly &&
                        fixture.ref_left->primary_path.next.load() ==
                            fixture.ref_anomaly,
                    "failed minipoa changed the graph");
            const auto repeated =
                StructuralBreakRepair::repairAnchorBoundedStructuralBreaks(
                    graph, "ref", managers, options);
            require(repeated.failure_cache_hits == 1 &&
                        repeated.msa_calls == 0,
                    "failed window was realigned instead of using cache");
        }

        {
            RaMeshMultiGenomeGraph graph(managers);
            const auto fixture = populate(graph);
            StructuralBreakRepair::Options options;
            options.enabled = true;
            options.msa_executable = msa;
            options.maximum_span = 10;
            const auto result =
                StructuralBreakRepair::repairAnchorBoundedStructuralBreaks(
                    graph, "ref", managers, options);
            require(result.structural_candidates == 0 &&
                        result.msa_calls == 0 && result.span_exceeded > 0,
                    "span filter did not reject before minipoa");
            require(fixture.qry_anomaly->parent_block == fixture.anomaly,
                    "span rejection changed the graph");
        }
    } catch (...) {
        std::filesystem::remove_all(temp);
        throw;
    }
    std::filesystem::remove_all(temp);
    return 0;
}
