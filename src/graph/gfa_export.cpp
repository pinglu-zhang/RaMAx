#include "ramesh.h"

#include "aligned_block_view.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <spdlog/spdlog.h>

namespace RaMesh::Gfa {

Version parseVersion(const std::string& value) {
    if (value == "1.0") return Version::V1_0;
    if (value == "1.1") return Version::V1_1;
    throw std::invalid_argument("GFA version must be exactly 1.0 or 1.1");
}

const char* versionString(Version version) noexcept {
    return version == Version::V1_0 ? "1.0" : "1.1";
}

namespace {

using SequenceKey = std::pair<SpeciesName, ChrName>;

struct SourceSequence {
    SpeciesName species;
    ChrName chromosome;
    const SeqPro::ManagerVariant* manager{nullptr};
    std::uint64_t length{0};
    std::set<std::uint64_t> breakpoints;
    std::vector<std::uint64_t> ordered_breakpoints;
    std::unordered_map<std::uint64_t, std::size_t> atom_by_start;
    std::vector<std::size_t> atoms;
};

struct ExactRelation {
    std::size_t left_sequence{0};
    std::uint64_t left_start{0};
    std::size_t right_sequence{0};
    std::uint64_t right_start{0};
    std::uint64_t length{0};
    bool reverse{false};

    auto orderingKey() const {
        return std::tie(left_sequence, left_start, right_sequence,
                        right_start, length, reverse);
    }
};

struct RunBuilder {
    std::size_t left_sequence{0};
    std::size_t right_sequence{0};
    std::int64_t left_first{0};
    std::int64_t right_first{0};
    std::int64_t left_last{0};
    std::int64_t right_last{0};
    int left_direction{1};
    int right_direction{1};
    std::uint64_t length{0};
};

struct Atom {
    std::size_t sequence{0};
    std::uint64_t start{0};
    std::uint64_t end{0};
};

struct Step {
    std::uint64_t node{0};
    bool reverse{false};
};

struct Edge {
    std::uint64_t from{0};
    bool from_reverse{false};
    std::uint64_t to{0};
    bool to_reverse{false};

    auto key() const {
        return std::tie(from, from_reverse, to, to_reverse);
    }

    bool operator<(const Edge& other) const { return key() < other.key(); }
};

struct Node {
    std::uint64_t id{0};
    std::size_t representative_atom{0};
    std::string sequence;
    std::size_t occurrence_count{0};
};

class OrientedDisjointSet {
public:
    explicit OrientedDisjointSet(std::size_t size)
        : parent_(size), rank_(size, 0), parity_(size, 0) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }

    std::pair<std::size_t, bool> find(std::size_t value) {
        if (parent_[value] == value) return {value, false};
        const auto [root, parent_parity] = find(parent_[value]);
        parity_[value] = static_cast<unsigned char>(
            parity_[value] ^ static_cast<unsigned char>(parent_parity));
        parent_[value] = root;
        return {root, parity_[value] != 0};
    }

    void unite(std::size_t left, std::size_t right, bool reverse) {
        auto [left_root, left_parity] = find(left);
        auto [right_root, right_parity] = find(right);
        if (left_root == right_root) {
            if ((left_parity ^ right_parity) != reverse) {
                throw std::runtime_error(
                    "Contradictory orientation in exact-run closure");
            }
            return;
        }
        if (rank_[left_root] < rank_[right_root]) {
            std::swap(left_root, right_root);
            std::swap(left_parity, right_parity);
        }
        parent_[right_root] = left_root;
        parity_[right_root] = static_cast<unsigned char>(
            left_parity ^ right_parity ^ reverse);
        if (rank_[left_root] == rank_[right_root]) ++rank_[left_root];
    }

private:
    std::vector<std::size_t> parent_;
    std::vector<unsigned char> rank_;
    std::vector<unsigned char> parity_;
};

char normalizedBase(char value) {
    return static_cast<char>(
        std::toupper(static_cast<unsigned char>(value)));
}

bool isExactBase(char value) {
    const char base = normalizedBase(value);
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

char complement(char value) {
    switch (normalizedBase(value)) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    case 'N': return 'N';
    default: throw std::runtime_error("Invalid graph sequence base");
    }
}

bool validSequence(const std::string& sequence) {
    return !sequence.empty() &&
           std::all_of(sequence.begin(), sequence.end(), [](char value) {
               const char base = normalizedBase(value);
               return base == 'A' || base == 'C' || base == 'G' ||
                      base == 'T' || base == 'N';
           });
}

bool validWalkField(const std::string& value) {
    if (value.empty() || value.front() == '*' || value.front() == '=') {
        return false;
    }
    return std::all_of(value.begin(), value.end(), [](unsigned char ch) {
        return ch >= 33 && ch <= 126;
    });
}

bool validPathName(const std::string& value) {
    if (!validWalkField(value)) return false;
    return value.find("+,") == std::string::npos &&
           value.find("-,") == std::string::npos;
}

std::filesystem::path temporaryPath(
    const std::filesystem::path& output_path) {
    const auto stamp = std::chrono::steady_clock::now()
                           .time_since_epoch().count();
    const auto thread_hash =
        std::hash<std::thread::id>{}(std::this_thread::get_id());
    for (std::uint64_t attempt = 0; ; ++attempt) {
        auto candidate = output_path;
        candidate += ".tmp." + std::to_string(stamp) + "." +
                     std::to_string(thread_hash) + "." +
                     std::to_string(attempt);
        if (!std::filesystem::exists(candidate)) return candidate;
    }
}

ExactRelation finishRun(const RunBuilder& run) {
    if (run.length == 0) {
        throw std::runtime_error("Cannot finish an empty exact run");
    }
    const auto left_min = std::min(run.left_first, run.left_last);
    const auto right_min = std::min(run.right_first, run.right_last);
    if (left_min < 0 || right_min < 0) {
        throw std::runtime_error("Exact run has a negative source coordinate");
    }
    return {
        run.left_sequence,
        static_cast<std::uint64_t>(left_min),
        run.right_sequence,
        static_cast<std::uint64_t>(right_min),
        run.length,
        run.left_direction != run.right_direction
    };
}

std::vector<ExactRelation> collectExactRelations(
    const Export::AlignedBlockView& block,
    const std::map<SequenceKey, std::size_t>& sequence_index) {
    std::vector<ExactRelation> relations;
    if (block.rows.size() < 2) return relations;
    const std::size_t width = block.rows.front().aligned.size();
    std::vector<std::uint64_t> consumed(block.rows.size(), 0);
    std::vector<std::optional<std::uint64_t>> positions(block.rows.size());
    std::map<std::pair<std::size_t, std::size_t>, RunBuilder> active;

    for (std::size_t column = 0; column < width; ++column) {
        std::array<std::vector<std::size_t>, 4> groups;
        positions.assign(block.rows.size(), std::nullopt);
        for (std::size_t row_index = 0; row_index < block.rows.size();
             ++row_index) {
            const auto& row = block.rows[row_index];
            const char base = row.aligned[column];
            if (base == '-') continue;
            if (consumed[row_index] >= row.segment->length) {
                throw std::runtime_error(
                    "Aligned Block consumes beyond Segment length");
            }
            const std::uint64_t position =
                row.segment->strand == Strand::FORWARD
                    ? static_cast<std::uint64_t>(row.segment->start) +
                          consumed[row_index]
                    : static_cast<std::uint64_t>(row.segment->start) +
                          row.segment->length - 1 - consumed[row_index];
            positions[row_index] = position;
            const char normalized = normalizedBase(base);
            if (isExactBase(normalized)) {
                const std::size_t group = normalized == 'A' ? 0
                    : normalized == 'C' ? 1
                    : normalized == 'G' ? 2 : 3;
                groups[group].push_back(row_index);
            }
            ++consumed[row_index];
        }

        std::set<std::pair<std::size_t, std::size_t>> current;
        for (const auto& group : groups) {
            if (group.size() < 2) continue;
            const std::size_t anchor = group.front();
            for (std::size_t index = 1; index < group.size(); ++index) {
                const std::size_t member = group[index];
                const auto key = std::pair{anchor, member};
                current.insert(key);
                const auto left_position =
                    static_cast<std::int64_t>(*positions[anchor]);
                const auto right_position =
                    static_cast<std::int64_t>(*positions[member]);
                const int left_direction =
                    block.rows[anchor].segment->strand == Strand::FORWARD
                        ? 1 : -1;
                const int right_direction =
                    block.rows[member].segment->strand == Strand::FORWARD
                        ? 1 : -1;
                auto active_it = active.find(key);
                const bool extends = active_it != active.end() &&
                    active_it->second.left_last + left_direction ==
                        left_position &&
                    active_it->second.right_last + right_direction ==
                        right_position;
                if (!extends && active_it != active.end()) {
                    relations.push_back(finishRun(active_it->second));
                    active.erase(active_it);
                    active_it = active.end();
                }
                if (active_it == active.end()) {
                    const auto left_sequence = sequence_index.at(
                        {block.rows[anchor].species,
                         block.rows[anchor].chromosome});
                    const auto right_sequence = sequence_index.at(
                        {block.rows[member].species,
                         block.rows[member].chromosome});
                    active.emplace(key, RunBuilder{
                        left_sequence, right_sequence,
                        left_position, right_position,
                        left_position, right_position,
                        left_direction, right_direction, 1});
                } else {
                    active_it->second.left_last = left_position;
                    active_it->second.right_last = right_position;
                    ++active_it->second.length;
                }
            }
        }

        for (auto it = active.begin(); it != active.end(); ) {
            if (current.contains(it->first)) {
                ++it;
            } else {
                relations.push_back(finishRun(it->second));
                it = active.erase(it);
            }
        }
    }

    for (const auto& [unused, run] : active) {
        (void)unused;
        relations.push_back(finishRun(run));
    }
    for (std::size_t row = 0; row < consumed.size(); ++row) {
        if (consumed[row] != block.rows[row].segment->length) {
            throw std::runtime_error(
                "Aligned Block does not consume complete Segment");
        }
    }
    return relations;
}

std::uint64_t projectBreakpoint(const ExactRelation& relation,
                                bool from_left,
                                std::uint64_t position) {
    const std::uint64_t source_start =
        from_left ? relation.left_start : relation.right_start;
    const std::uint64_t target_start =
        from_left ? relation.right_start : relation.left_start;
    if (position < source_start ||
        position > source_start + relation.length) {
        throw std::runtime_error("Breakpoint is outside exact relation");
    }
    const std::uint64_t offset = position - source_start;
    return relation.reverse
        ? target_start + relation.length - offset
        : target_start + offset;
}

void propagateBreakpoints(std::vector<SourceSequence>& sequences,
                          const std::vector<ExactRelation>& relations) {
    bool changed = true;
    while (changed) {
        changed = false;
        for (const auto& relation : relations) {
            auto& left = sequences[relation.left_sequence].breakpoints;
            auto& right = sequences[relation.right_sequence].breakpoints;
            std::vector<std::uint64_t> left_points;
            std::vector<std::uint64_t> right_points;
            for (auto it = left.lower_bound(relation.left_start);
                 it != left.end() &&
                 *it <= relation.left_start + relation.length; ++it) {
                left_points.push_back(*it);
            }
            for (auto it = right.lower_bound(relation.right_start);
                 it != right.end() &&
                 *it <= relation.right_start + relation.length; ++it) {
                right_points.push_back(*it);
            }
            for (const auto point : left_points) {
                changed |= right.insert(
                    projectBreakpoint(relation, true, point)).second;
            }
            for (const auto point : right_points) {
                changed |= left.insert(
                    projectBreakpoint(relation, false, point)).second;
            }
        }
    }
}

std::string atomSequence(const SourceSequence& source, const Atom& atom) {
    const std::uint64_t length = atom.end - atom.start;
    if (atom.start > std::numeric_limits<Coord_t>::max() ||
        length > std::numeric_limits<Coord_t>::max()) {
        throw std::runtime_error(
            "GFA atom exceeds SequenceManager coordinate range");
    }
    std::string sequence = Export::fetchSequence(
        *source.manager, source.chromosome,
        static_cast<Coord_t>(atom.start), static_cast<Coord_t>(length));
    std::transform(sequence.begin(), sequence.end(), sequence.begin(),
                   [](unsigned char value) {
                       return static_cast<char>(std::toupper(value));
                   });
    return sequence;
}

bool orientedSequencesEqual(std::string left, std::string right,
                            bool reverse) {
    if (reverse) reverseComplement(right);
    return left == right;
}

Edge canonicalEdge(Edge edge) {
    Edge reversed{edge.to, !edge.to_reverse,
                  edge.from, !edge.from_reverse};
    return reversed < edge ? reversed : edge;
}

void validateGraph(
    const std::vector<SourceSequence>& sequences,
    const std::vector<Atom>& atoms,
    const std::vector<Node>& nodes,
    const std::vector<std::uint64_t>& node_for_atom,
    const std::vector<bool>& reverse_for_atom,
    const std::vector<std::vector<Step>>& walks,
    const std::set<Edge>& edges) {
    if (walks.size() != sequences.size()) {
        throw std::runtime_error("GFA walk count does not match contig count");
    }
    for (const auto& node : nodes) {
        if (node.id == 0 || node.id > nodes.size() ||
            !validSequence(node.sequence)) {
            throw std::runtime_error("GFA contains an invalid Segment");
        }
    }
    for (std::size_t sequence_index = 0;
         sequence_index < sequences.size(); ++sequence_index) {
        const auto& source = sequences[sequence_index];
        const auto& path_atoms = source.atoms;
        const auto& walk = walks[sequence_index];
        if (path_atoms.size() != walk.size()) {
            throw std::runtime_error("GFA walk lost an atomic interval");
        }
        std::string expected = Export::fetchSequence(
            *source.manager, source.chromosome, 0,
            static_cast<Coord_t>(source.length));
        std::transform(expected.begin(), expected.end(), expected.begin(),
                       [](unsigned char value) {
                           return static_cast<char>(std::toupper(value));
                       });
        std::uint64_t offset = 0;
        for (std::size_t index = 0; index < walk.size(); ++index) {
            const Step step = walk[index];
            if (step.node == 0 || step.node > nodes.size()) {
                throw std::runtime_error("GFA walk references an unknown node");
            }
            const Atom& atom = atoms[path_atoms[index]];
            if (node_for_atom[path_atoms[index]] != step.node ||
                reverse_for_atom[path_atoms[index]] != step.reverse) {
                throw std::runtime_error("GFA walk/atom mapping is inconsistent");
            }
            const std::string& node_sequence = nodes[step.node - 1].sequence;
            for (std::size_t base = 0; base < node_sequence.size(); ++base) {
                const char observed = step.reverse
                    ? complement(node_sequence[node_sequence.size() - 1 - base])
                    : node_sequence[base];
                if (offset >= expected.size() || observed != expected[offset]) {
                    throw std::runtime_error(
                        "GFA walk does not spell cleaned FASTA: " +
                        source.species + "." + source.chromosome);
                }
                ++offset;
            }
            if (index > 0) {
                const Step previous = walk[index - 1];
                if (!edges.contains(canonicalEdge(
                        {previous.node, previous.reverse,
                         step.node, step.reverse}))) {
                    throw std::runtime_error(
                        "GFA walk adjacency has no corresponding Link");
                }
            }
            if (atom.end <= atom.start) {
                throw std::runtime_error("GFA contains an empty atom");
            }
        }
        if (offset != source.length || offset != expected.size()) {
            throw std::runtime_error(
                "GFA walk content or length does not match cleaned FASTA: " +
                source.species + "." + source.chromosome);
        }
    }
}

}  // namespace

GfaExportStats exportGraph(
    const RaMeshMultiGenomeGraph& graph,
    const FilePath& gfa_path,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    const GfaExportOptions& options) {
    namespace fs = std::filesystem;
    const auto started = std::chrono::steady_clock::now();
    if (managers.empty()) {
        throw std::runtime_error("No sequence managers provided for GFA export");
    }

    std::vector<SourceSequence> sequences;
    std::map<SequenceKey, std::size_t> sequence_index;
    for (const auto& [species, manager] : managers) {
        if (!validWalkField(species)) {
            throw std::runtime_error("Invalid GFA sample identifier: " + species);
        }
        for (const auto& chromosome : Export::fetchNames(*manager)) {
            if (!validWalkField(chromosome)) {
                throw std::runtime_error(
                    "Invalid GFA sequence identifier: " + chromosome);
            }
            const std::uint64_t length = Export::fetchLength(
                *manager, chromosome);
            if (length == 0 || length > std::numeric_limits<Coord_t>::max()) {
                throw std::runtime_error(
                    "GFA requires a non-empty contig within the coordinate range: " +
                    species + "." + chromosome);
            }
            const std::size_t index = sequences.size();
            if (!sequence_index.emplace(SequenceKey{species, chromosome},
                                        index).second) {
                throw std::runtime_error(
                    "Duplicate GFA sample/contig identifier: " + species +
                    "." + chromosome);
            }
            SourceSequence source;
            source.species = species;
            source.chromosome = chromosome;
            source.manager = manager.get();
            source.length = length;
            source.breakpoints.insert(0);
            source.breakpoints.insert(length);
            sequences.push_back(std::move(source));
        }
    }

    GfaExportStats stats;
    stats.blocks_seen = graph.blocks.size();
    std::vector<BlockPtr> strong_blocks;
    strong_blocks.reserve(graph.blocks.size());
    for (const auto& weak : graph.blocks) {
        if (auto block = weak.lock()) strong_blocks.push_back(std::move(block));
    }

    std::vector<std::vector<ExactRelation>> relations_by_block(
        strong_blocks.size());
    std::vector<std::string> block_errors(strong_blocks.size());
    const int export_threads = std::max(1, options.threads);
#pragma omp parallel for schedule(dynamic) num_threads(export_threads)
    for (std::int64_t index = 0;
         index < static_cast<std::int64_t>(strong_blocks.size()); ++index) {
        try {
            const auto view = Export::prepareAlignedBlock(
                strong_blocks[static_cast<std::size_t>(index)], managers,
                options.only_primary);
            if (view.rows.size() >= 2) {
                relations_by_block[static_cast<std::size_t>(index)] =
                    collectExactRelations(view, sequence_index);
            }
        } catch (const std::exception& error) {
            block_errors[static_cast<std::size_t>(index)] = error.what();
        }
    }

    std::vector<ExactRelation> relations;
    for (std::size_t index = 0; index < relations_by_block.size(); ++index) {
        if (!block_errors[index].empty()) {
            throw std::runtime_error(
                "GFA Block projection failed: " + block_errors[index]);
        }
        if (!relations_by_block[index].empty()) ++stats.blocks_used;
        relations.insert(relations.end(),
                         relations_by_block[index].begin(),
                         relations_by_block[index].end());
    }
    std::sort(relations.begin(), relations.end(),
              [](const ExactRelation& left, const ExactRelation& right) {
                  return left.orderingKey() < right.orderingKey();
              });
    relations.erase(std::unique(relations.begin(), relations.end(),
        [](const ExactRelation& left, const ExactRelation& right) {
            return left.orderingKey() == right.orderingKey();
        }), relations.end());
    stats.exact_runs = relations.size();

    for (const auto& relation : relations) {
        auto validate_endpoint = [&](std::size_t sequence,
                                     std::uint64_t start) {
            if (sequence >= sequences.size() ||
                start > sequences[sequence].length ||
                relation.length > sequences[sequence].length - start) {
                throw std::runtime_error(
                    "Exact run exceeds cleaned sequence bounds");
            }
        };
        validate_endpoint(relation.left_sequence, relation.left_start);
        validate_endpoint(relation.right_sequence, relation.right_start);
        sequences[relation.left_sequence].breakpoints.insert(
            relation.left_start);
        sequences[relation.left_sequence].breakpoints.insert(
            relation.left_start + relation.length);
        sequences[relation.right_sequence].breakpoints.insert(
            relation.right_start);
        sequences[relation.right_sequence].breakpoints.insert(
            relation.right_start + relation.length);
    }
    propagateBreakpoints(sequences, relations);

    std::vector<Atom> atoms;
    for (std::size_t sequence = 0; sequence < sequences.size(); ++sequence) {
        auto& source = sequences[sequence];
        source.ordered_breakpoints.assign(
            source.breakpoints.begin(), source.breakpoints.end());
        for (std::size_t index = 1;
             index < source.ordered_breakpoints.size(); ++index) {
            const std::uint64_t start = source.ordered_breakpoints[index - 1];
            const std::uint64_t end = source.ordered_breakpoints[index];
            if (start >= end) continue;
            const std::size_t atom = atoms.size();
            atoms.push_back({sequence, start, end});
            source.atom_by_start.emplace(start, atom);
            source.atoms.push_back(atom);
        }
    }
    stats.atomic_intervals = atoms.size();
    OrientedDisjointSet dsu(atoms.size());

    for (const auto& relation : relations) {
        const auto& left_points =
            sequences[relation.left_sequence].ordered_breakpoints;
        auto it = std::lower_bound(
            left_points.begin(), left_points.end(), relation.left_start);
        if (it == left_points.end() || *it != relation.left_start) {
            throw std::runtime_error("Exact-run left boundary was not atomized");
        }
        while (*it < relation.left_start + relation.length) {
            const auto next = std::next(it);
            if (next == left_points.end() ||
                *next > relation.left_start + relation.length) {
                throw std::runtime_error("Exact run was incompletely atomized");
            }
            const std::uint64_t left_start = *it;
            const std::uint64_t left_end = *next;
            const std::uint64_t left_offset =
                left_start - relation.left_start;
            const std::uint64_t right_start = relation.reverse
                ? relation.right_start + relation.length -
                      (left_end - relation.left_start)
                : relation.right_start + left_offset;
            const auto left_atom =
                sequences[relation.left_sequence].atom_by_start.at(left_start);
            const auto right_atom =
                sequences[relation.right_sequence].atom_by_start.at(right_start);
            if (atoms[left_atom].end - atoms[left_atom].start !=
                atoms[right_atom].end - atoms[right_atom].start) {
                throw std::runtime_error(
                    "Exact-run atom lengths do not agree");
            }
            const std::string left_sequence = atomSequence(
                sequences[relation.left_sequence], atoms[left_atom]);
            const std::string right_sequence = atomSequence(
                sequences[relation.right_sequence], atoms[right_atom]);
            if (!orientedSequencesEqual(
                    left_sequence, right_sequence, relation.reverse)) {
                throw std::runtime_error(
                    "Exact-run atom sequences do not agree");
            }
            dsu.unite(left_atom, right_atom, relation.reverse);
            it = next;
        }
    }

    std::unordered_map<std::size_t, std::size_t> component_by_root;
    std::vector<Node> nodes;
    std::vector<std::uint64_t> node_for_atom(atoms.size(), 0);
    std::vector<bool> reverse_for_atom(atoms.size(), false);
    std::vector<bool> representative_parity;
    for (std::size_t atom = 0; atom < atoms.size(); ++atom) {
        const auto [root, parity] = dsu.find(atom);
        auto component = component_by_root.find(root);
        if (component == component_by_root.end()) {
            const std::size_t node_index = nodes.size();
            component_by_root.emplace(root, node_index);
            const Atom& representative = atoms[atom];
            std::string sequence = atomSequence(
                sequences[representative.sequence], representative);
            if (!validSequence(sequence)) {
                throw std::runtime_error(
                    "Cleaned FASTA contains a base outside A/C/G/T/N");
            }
            nodes.push_back({node_index + 1, atom, std::move(sequence), 0});
            representative_parity.push_back(parity);
            component = component_by_root.find(root);
        }
        Node& node = nodes[component->second];
        ++node.occurrence_count;
        node_for_atom[atom] = node.id;
        reverse_for_atom[atom] =
            parity ^ representative_parity[component->second];
    }

    std::vector<std::vector<Step>> walks(sequences.size());
    std::set<Edge> edges;
    for (std::size_t sequence = 0; sequence < sequences.size(); ++sequence) {
        auto& walk = walks[sequence];
        walk.reserve(sequences[sequence].atoms.size());
        for (const auto atom : sequences[sequence].atoms) {
            walk.push_back({node_for_atom[atom], reverse_for_atom[atom]});
            if (walk.size() > 1) {
                const Step& left = walk[walk.size() - 2];
                const Step& right = walk.back();
                edges.insert(canonicalEdge(
                    {left.node, left.reverse, right.node, right.reverse}));
            }
        }
    }
    validateGraph(sequences, atoms, nodes, node_for_atom,
                  reverse_for_atom, walks, edges);

    stats.nodes = nodes.size();
    stats.edges = edges.size();
    stats.walks = walks.size();
    for (const auto& node : nodes) {
        stats.graph_sequence_bp += node.sequence.size();
        if (node.occurrence_count > 1) {
            stats.shared_sequence_bp += node.sequence.size();
        } else {
            stats.private_sequence_bp += node.sequence.size();
        }
    }

    std::set<std::string> path_names;
    if (options.version == Version::V1_0) {
        for (const auto& source : sequences) {
            const std::string name = Export::qualifiedName(
                source.species, source.chromosome);
            if (!validPathName(name) || !path_names.insert(name).second) {
                throw std::runtime_error(
                    "Invalid or duplicate GFA 1.0 path name: " + name);
            }
        }
    }

    const std::string reference_species = !graph.reference_order.empty()
        ? graph.reference_order.front() : managers.begin()->first;
    if (!validWalkField(reference_species)) {
        throw std::runtime_error(
            "Invalid GFA reference sample identifier: " + reference_species);
    }

    const fs::path output_path = fs::absolute(gfa_path);
    if (!output_path.parent_path().empty()) {
        fs::create_directories(output_path.parent_path());
    }
    const fs::path temporary_output = temporaryPath(output_path);
    struct TemporaryGuard {
        fs::path path;
        bool keep{false};
        ~TemporaryGuard() {
            if (!keep) {
                std::error_code error;
                fs::remove(path, error);
            }
        }
    } guard{temporary_output};

    std::ofstream output(
        temporary_output, std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error(
            "Cannot open GFA output: " + temporary_output.string());
    }
    if (options.version == Version::V1_0) {
        output << "H\tVN:Z:1.0\n";
    } else {
        output << "H\tVN:Z:1.1\tRS:Z:" << reference_species << '\n';
    }
    for (const auto& node : nodes) {
        output << "S\t" << node.id << '\t' << node.sequence << '\n';
    }
    for (const auto& edge : edges) {
        output << "L\t" << edge.from << '\t'
               << (edge.from_reverse ? '-' : '+') << '\t'
               << edge.to << '\t' << (edge.to_reverse ? '-' : '+')
               << "\t0M\n";
    }
    for (std::size_t index = 0; index < sequences.size(); ++index) {
        const auto& source = sequences[index];
        const auto& walk = walks[index];
        if (options.version == Version::V1_0) {
            output << "P\t"
                   << Export::qualifiedName(source.species, source.chromosome)
                   << '\t';
            for (std::size_t step = 0; step < walk.size(); ++step) {
                if (step != 0) output << ',';
                output << walk[step].node
                       << (walk[step].reverse ? '-' : '+');
            }
            output << "\t*\n";
        } else {
            output << "W\t" << source.species << "\t0\t"
                   << source.chromosome << "\t0\t" << source.length
                   << '\t';
            for (const auto& step : walk) {
                output << (step.reverse ? '<' : '>') << step.node;
            }
            output << '\n';
        }
    }
    output.close();
    if (!output) throw std::runtime_error("Failed to finalize GFA output");

    std::error_code rename_error;
    fs::rename(temporary_output, output_path, rename_error);
    if (rename_error) {
        throw std::runtime_error(
            "Failed to atomically replace GFA output: " +
            rename_error.message());
    }
    guard.keep = true;

    stats.elapsed_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    spdlog::info(
        "GFA {} export complete: blocks(seen/used)={}/{}, exact_runs={}, "
        "atoms={}, nodes={}, edges={}, paths={}, graph_bp={}, "
        "shared_bp={}, private_bp={}, elapsed_seconds={:.3f}",
        versionString(options.version), stats.blocks_seen, stats.blocks_used,
        stats.exact_runs, stats.atomic_intervals, stats.nodes, stats.edges,
        stats.walks, stats.graph_sequence_bp, stats.shared_sequence_bp,
        stats.private_sequence_bp, stats.elapsed_seconds);
    return stats;
}

}  // namespace RaMesh::Gfa

namespace RaMesh {

Gfa::GfaExportStats RaMeshMultiGenomeGraph::exportToGfa(
    const FilePath& gfa_path,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    const Gfa::GfaExportOptions& options) const {
    return Gfa::exportGraph(*this, gfa_path, managers, options);
}

}  // namespace RaMesh
