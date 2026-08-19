#include "ramesh.h"

#include "aligned_block_view.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <deque>
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
#include <unordered_set>
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

Profile parseProfile(const std::string& value) {
    if (value == "exact") return Profile::EXACT;
    if (value == "compact") return Profile::COMPACT;
    throw std::invalid_argument("GFA profile must be exactly exact or compact");
}

const char* profileString(Profile profile) noexcept {
    return profile == Profile::EXACT ? "exact" : "compact";
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
    std::uint64_t source_start{0};
    std::uint64_t source_end{0};
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
    std::size_t source_sequence{0};
    std::uint64_t source_start{0};
};

struct GraphMetrics {
    std::size_t nodes{0};
    std::size_t edges{0};
    std::size_t steps{0};
    std::size_t nodes_le_10{0};
    std::size_t nodes_lt_25{0};
    std::size_t nodes_lt_100{0};
    std::uint64_t graph_bp{0};
    std::uint64_t homology_mass{0};
    std::uint64_t median_node_length{0};
};

struct GraphIR {
    std::vector<Node> nodes;
    std::vector<std::vector<Step>> walks;
    std::set<Edge> edges;
    std::size_t atomic_intervals{0};
};

struct IslandReplacement {
    std::size_t walk{0};
    std::size_t begin{0};
    std::size_t end{0};
    std::uint64_t replacement_node{0};
};

struct IslandCandidate {
    std::size_t reference_walk{0};
    std::size_t reference_begin{0};
    std::size_t reference_end{0};
    std::vector<IslandReplacement> replacements;
    std::vector<std::pair<std::string, std::vector<std::size_t>>> alleles;
    std::set<std::uint64_t> interior_nodes;
    std::uint64_t old_homology_mass{0};
    std::uint64_t new_homology_mass{0};
    std::uint64_t old_graph_bp{0};
    std::uint64_t new_graph_bp{0};
    std::size_t old_short_nodes{0};
    std::size_t new_short_nodes{0};
};

constexpr std::uint64_t COMPACT_MIN_EXACT_RUN = 25;
constexpr std::uint64_t COMPACT_SHORT_NODE = 25;
constexpr std::uint64_t COMPACT_STABLE_ANCHOR = 100;
constexpr std::uint64_t COMPACT_MAX_REFERENCE_SPAN = 1000;
constexpr std::uint64_t COMPACT_MAX_PATH_SPAN = 3000;
constexpr std::size_t COMPACT_MIN_SHORT_NODES = 2;
constexpr std::size_t COMPACT_MAX_ALLELES = 64;
constexpr double COMPACT_MIN_HOMOLOGY_RETENTION = 0.995;
constexpr double COMPACT_MAX_GRAPH_BP_RATIO = 1.02;

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

using Handle = std::pair<std::uint64_t, bool>;
using HandlePair = std::pair<Handle, Handle>;

Handle reverseHandle(const Handle& handle) {
    return {handle.first, !handle.second};
}

HandlePair canonicalHandlePair(HandlePair pair) {
    HandlePair reversed{reverseHandle(pair.second),
                        reverseHandle(pair.first)};
    return reversed < pair ? reversed : pair;
}

HandlePair orientPairByFirstOccurrence(
    const HandlePair& canonical,
    const std::vector<std::vector<Step>>& walks) {
    const HandlePair reversed{reverseHandle(canonical.second),
                              reverseHandle(canonical.first)};
    std::optional<std::tuple<std::size_t, std::uint64_t, HandlePair>> best;
    for (std::size_t walk_index = 0; walk_index < walks.size(); ++walk_index) {
        const auto& walk = walks[walk_index];
        for (std::size_t position = 1; position < walk.size(); ++position) {
            const HandlePair observed{
                {walk[position - 1].node, walk[position - 1].reverse},
                {walk[position].node, walk[position].reverse}};
            if (observed != canonical && observed != reversed) continue;
            const auto key = std::tuple{
                walk_index, walk[position - 1].source_start, observed};
            if (!best || key < *best) best = key;
        }
    }
    if (!best) {
        throw std::runtime_error("Cannot orient a unitig candidate");
    }
    return std::get<2>(*best);
}

std::string orientedNodeSequence(const Node& node, bool reverse) {
    if (!reverse) return node.sequence;
    std::string sequence = node.sequence;
    reverseComplement(sequence);
    return sequence;
}

std::set<Edge> rebuildEdges(const std::vector<std::vector<Step>>& walks) {
    std::set<Edge> edges;
    for (const auto& walk : walks) {
        for (std::size_t index = 1; index < walk.size(); ++index) {
            const auto& left = walk[index - 1];
            const auto& right = walk[index];
            edges.insert(canonicalEdge(
                {left.node, left.reverse, right.node, right.reverse}));
        }
    }
    return edges;
}

void rebuildOccurrences(std::vector<Node>& nodes,
                        const std::vector<std::vector<Step>>& walks) {
    for (auto& node : nodes) {
        node.occurrence_count = 0;
        node.source_sequence = std::numeric_limits<std::size_t>::max();
        node.source_start = std::numeric_limits<std::uint64_t>::max();
    }
    for (std::size_t sequence = 0; sequence < walks.size(); ++sequence) {
        for (const auto& step : walks[sequence]) {
            if (step.node == 0 || step.node > nodes.size()) {
                throw std::runtime_error(
                    "GFA occurrence references an unknown node");
            }
            auto& node = nodes[step.node - 1];
            ++node.occurrence_count;
            const auto key = std::pair{sequence, step.source_start};
            const auto current =
                std::pair{node.source_sequence, node.source_start};
            if (key < current) {
                node.source_sequence = sequence;
                node.source_start = step.source_start;
            }
        }
    }
}

void normalizeGraph(std::vector<Node>& nodes,
                    std::vector<std::vector<Step>>& walks,
                    std::set<Edge>& edges) {
    rebuildOccurrences(nodes, walks);
    std::vector<std::size_t> order;
    order.reserve(nodes.size());
    for (std::size_t index = 0; index < nodes.size(); ++index) {
        if (nodes[index].occurrence_count != 0) order.push_back(index);
    }
    std::sort(order.begin(), order.end(), [&](std::size_t left,
                                               std::size_t right) {
        const auto& left_node = nodes[left];
        const auto& right_node = nodes[right];
        return std::tie(left_node.source_sequence, left_node.source_start,
                        left_node.id) <
               std::tie(right_node.source_sequence, right_node.source_start,
                        right_node.id);
    });
    std::vector<std::uint64_t> new_id(nodes.size() + 1, 0);
    std::vector<Node> normalized;
    normalized.reserve(order.size());
    for (const auto old_index : order) {
        Node node = std::move(nodes[old_index]);
        node.id = normalized.size() + 1;
        new_id[old_index + 1] = node.id;
        normalized.push_back(std::move(node));
    }
    for (auto& walk : walks) {
        for (auto& step : walk) {
            if (step.node >= new_id.size() || new_id[step.node] == 0) {
                throw std::runtime_error(
                    "GFA normalization lost a referenced node");
            }
            step.node = new_id[step.node];
        }
    }
    nodes = std::move(normalized);
    rebuildOccurrences(nodes, walks);
    edges = rebuildEdges(walks);
}

GraphMetrics collectMetrics(const std::vector<Node>& nodes,
                            const std::vector<std::vector<Step>>& walks,
                            const std::set<Edge>& edges) {
    GraphMetrics metrics;
    metrics.nodes = nodes.size();
    metrics.edges = edges.size();
    std::vector<std::uint64_t> lengths;
    lengths.reserve(nodes.size());
    for (const auto& node : nodes) {
        const auto length = static_cast<std::uint64_t>(node.sequence.size());
        lengths.push_back(length);
        metrics.graph_bp += length;
        if (length <= 10) ++metrics.nodes_le_10;
        if (length < COMPACT_SHORT_NODE) ++metrics.nodes_lt_25;
        if (length < 100) ++metrics.nodes_lt_100;
        if (node.occurrence_count > 1) {
            metrics.homology_mass +=
                length * static_cast<std::uint64_t>(
                    node.occurrence_count - 1);
        }
    }
    for (const auto& walk : walks) metrics.steps += walk.size();
    if (!lengths.empty()) {
        const std::size_t middle = lengths.size() / 2;
        std::nth_element(lengths.begin(), lengths.begin() + middle,
                         lengths.end());
        metrics.median_node_length = lengths[middle];
    }
    return metrics;
}

void validateGraph(
    const std::vector<SourceSequence>& sequences,
    const std::vector<Node>& nodes,
    const std::vector<std::vector<Step>>& walks,
    const std::set<Edge>& edges) {
    if (walks.size() != sequences.size()) {
        throw std::runtime_error("GFA walk count does not match contig count");
    }
    std::vector<std::size_t> observed_occurrences(nodes.size(), 0);
    for (const auto& node : nodes) {
        if (node.id == 0 || node.id > nodes.size() ||
            !validSequence(node.sequence)) {
            throw std::runtime_error("GFA contains an invalid Segment");
        }
    }
    for (std::size_t sequence_index = 0;
         sequence_index < sequences.size(); ++sequence_index) {
        const auto& source = sequences[sequence_index];
        const auto& walk = walks[sequence_index];
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
            ++observed_occurrences[step.node - 1];
            if (step.source_start != offset ||
                step.source_end <= step.source_start) {
                throw std::runtime_error(
                    "GFA provenance is not a complete ordered partition: " +
                    source.species + "." + source.chromosome);
            }
            const std::string& node_sequence = nodes[step.node - 1].sequence;
            if (step.source_end - step.source_start != node_sequence.size()) {
                throw std::runtime_error(
                    "GFA provenance length does not match node length");
            }
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
        }
        if (offset != source.length || offset != expected.size()) {
            throw std::runtime_error(
                "GFA walk content or length does not match cleaned FASTA: " +
                source.species + "." + source.chromosome);
        }
    }
    for (std::size_t index = 0; index < nodes.size(); ++index) {
        if (observed_occurrences[index] != nodes[index].occurrence_count) {
            throw std::runtime_error(
                "GFA node occurrence count is inconsistent");
        }
    }
}

using Occurrence = std::pair<std::size_t, std::size_t>;

bool handleEquals(const Step& step, const Handle& handle) {
    return step.node == handle.first && step.reverse == handle.second;
}

bool deterministicPair(
    const HandlePair& pair,
    const std::vector<std::vector<Occurrence>>& occurrences,
    const std::vector<std::vector<Step>>& walks) {
    const auto& left = pair.first;
    const auto& right = pair.second;
    if (left.first == right.first) return false;
    for (const auto& [walk_index, position] : occurrences[left.first]) {
        const auto& walk = walks[walk_index];
        const auto& step = walk[position];
        if (step.reverse == left.second) {
            if (position + 1 >= walk.size() ||
                !handleEquals(walk[position + 1], right)) return false;
        } else {
            if (position == 0 ||
                !handleEquals(walk[position - 1], reverseHandle(right))) {
                return false;
            }
        }
    }
    for (const auto& [walk_index, position] : occurrences[right.first]) {
        const auto& walk = walks[walk_index];
        const auto& step = walk[position];
        if (step.reverse == right.second) {
            if (position == 0 ||
                !handleEquals(walk[position - 1], left)) return false;
        } else {
            if (position + 1 >= walk.size() ||
                !handleEquals(walk[position + 1], reverseHandle(left))) {
                return false;
            }
        }
    }
    return !occurrences[left.first].empty() &&
           occurrences[left.first].size() == occurrences[right.first].size();
}

std::size_t compactUnitigs(std::vector<Node>& nodes,
                           std::vector<std::vector<Step>>& walks,
                           std::set<Edge>& edges) {
    std::size_t merges = 0;
    while (true) {
        std::vector<std::vector<Occurrence>> occurrences(nodes.size() + 1);
        std::set<HandlePair> possible;
        for (std::size_t walk_index = 0; walk_index < walks.size();
             ++walk_index) {
            const auto& walk = walks[walk_index];
            for (std::size_t position = 0; position < walk.size(); ++position) {
                occurrences[walk[position].node].push_back(
                    {walk_index, position});
                if (position + 1 < walk.size()) {
                    possible.insert(canonicalHandlePair({
                        {walk[position].node, walk[position].reverse},
                        {walk[position + 1].node,
                         walk[position + 1].reverse}}));
                }
            }
        }
        std::vector<HandlePair> candidates;
        for (const auto& pair : possible) {
            if (deterministicPair(pair, occurrences, walks)) {
                candidates.push_back(pair);
            }
        }
        std::sort(candidates.begin(), candidates.end(),
                  [&](const HandlePair& left, const HandlePair& right) {
            const auto left_key = std::tuple{
                nodes[left.first.first - 1].source_sequence,
                nodes[left.first.first - 1].source_start,
                nodes[left.second.first - 1].source_sequence,
                nodes[left.second.first - 1].source_start,
                left};
            const auto right_key = std::tuple{
                nodes[right.first.first - 1].source_sequence,
                nodes[right.first.first - 1].source_start,
                nodes[right.second.first - 1].source_sequence,
                nodes[right.second.first - 1].source_start,
                right};
            return left_key < right_key;
        });
        std::vector<HandlePair> selected;
        std::unordered_set<std::uint64_t> used;
        for (const auto& canonical_pair : candidates) {
            const auto pair = orientPairByFirstOccurrence(
                canonical_pair, walks);
            if (used.contains(pair.first.first) ||
                used.contains(pair.second.first)) continue;
            used.insert(pair.first.first);
            used.insert(pair.second.first);
            selected.push_back(pair);
        }
        if (selected.empty()) break;

        std::map<HandlePair, std::pair<std::uint64_t, bool>> replacements;
        for (const auto& pair : selected) {
            const auto new_id = static_cast<std::uint64_t>(nodes.size() + 1);
            std::string sequence = orientedNodeSequence(
                nodes[pair.first.first - 1], pair.first.second);
            sequence += orientedNodeSequence(
                nodes[pair.second.first - 1], pair.second.second);
            const auto& first_node = nodes[pair.first.first - 1];
            const auto& second_node = nodes[pair.second.first - 1];
            const auto source_key = std::min(
                std::pair{first_node.source_sequence, first_node.source_start},
                std::pair{second_node.source_sequence,
                          second_node.source_start});
            nodes.push_back({new_id, 0, std::move(sequence), 0,
                             source_key.first, source_key.second});
            replacements.emplace(pair, std::pair{new_id, false});
            replacements.emplace(
                HandlePair{reverseHandle(pair.second),
                           reverseHandle(pair.first)},
                std::pair{new_id, true});
        }

        for (auto& walk : walks) {
            std::vector<Step> rewritten;
            rewritten.reserve(walk.size());
            for (std::size_t position = 0; position < walk.size();) {
                if (position + 1 < walk.size()) {
                    const HandlePair pair{
                        {walk[position].node, walk[position].reverse},
                        {walk[position + 1].node,
                         walk[position + 1].reverse}};
                    const auto replacement = replacements.find(pair);
                    if (replacement != replacements.end()) {
                        rewritten.push_back({replacement->second.first,
                                             replacement->second.second,
                                             walk[position].source_start,
                                             walk[position + 1].source_end});
                        position += 2;
                        continue;
                    }
                }
                rewritten.push_back(walk[position]);
                ++position;
            }
            walk = std::move(rewritten);
        }
        merges += selected.size();
        normalizeGraph(nodes, walks, edges);
    }
    return merges;
}

GraphIR buildInitialGraph(
    const std::vector<SourceSequence>& base_sequences,
    const std::vector<ExactRelation>& relations) {
    std::vector<SourceSequence> sequences = base_sequences;
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
                throw std::runtime_error("Exact-run atom lengths do not agree");
            }
            const std::string left_sequence = atomSequence(
                sequences[relation.left_sequence], atoms[left_atom]);
            const std::string right_sequence = atomSequence(
                sequences[relation.right_sequence], atoms[right_atom]);
            if (!orientedSequencesEqual(
                    left_sequence, right_sequence, relation.reverse)) {
                throw std::runtime_error("Exact-run atom sequences do not agree");
            }
            dsu.unite(left_atom, right_atom, relation.reverse);
            it = next;
        }
    }

    GraphIR graph;
    graph.atomic_intervals = atoms.size();
    std::unordered_map<std::size_t, std::size_t> component_by_root;
    std::vector<std::uint64_t> node_for_atom(atoms.size(), 0);
    std::vector<bool> reverse_for_atom(atoms.size(), false);
    std::vector<bool> representative_parity;
    for (std::size_t atom = 0; atom < atoms.size(); ++atom) {
        const auto [root, parity] = dsu.find(atom);
        auto component = component_by_root.find(root);
        if (component == component_by_root.end()) {
            const std::size_t node_index = graph.nodes.size();
            component_by_root.emplace(root, node_index);
            const Atom& representative = atoms[atom];
            std::string sequence = atomSequence(
                sequences[representative.sequence], representative);
            if (!validSequence(sequence)) {
                throw std::runtime_error(
                    "Cleaned FASTA contains a base outside A/C/G/T/N");
            }
            graph.nodes.push_back({node_index + 1, atom, std::move(sequence),
                                   0, representative.sequence,
                                   representative.start});
            representative_parity.push_back(parity);
            component = component_by_root.find(root);
        }
        Node& node = graph.nodes[component->second];
        ++node.occurrence_count;
        node_for_atom[atom] = node.id;
        reverse_for_atom[atom] =
            parity ^ representative_parity[component->second];
    }

    graph.walks.resize(sequences.size());
    for (std::size_t sequence = 0; sequence < sequences.size(); ++sequence) {
        auto& walk = graph.walks[sequence];
        walk.reserve(sequences[sequence].atoms.size());
        for (const auto atom : sequences[sequence].atoms) {
            walk.push_back({node_for_atom[atom], reverse_for_atom[atom],
                            atoms[atom].start, atoms[atom].end});
        }
    }
    graph.edges = rebuildEdges(graph.walks);
    validateGraph(sequences, graph.nodes, graph.walks, graph.edges);
    return graph;
}

void validatePathNames(const std::vector<SourceSequence>& sequences,
                       Version version) {
    if (version != Version::V1_0) return;
    std::set<std::string> path_names;
    for (const auto& source : sequences) {
        const std::string name = Export::qualifiedName(
            source.species, source.chromosome);
        if (!validPathName(name) || !path_names.insert(name).second) {
            throw std::runtime_error(
                "Invalid or duplicate GFA 1.0 path name: " + name);
        }
    }
}

void atomicWrite(const std::filesystem::path& output_path,
                 const std::function<void(std::ostream&)>& writer) {
    namespace fs = std::filesystem;
    const fs::path absolute = fs::absolute(output_path);
    if (!absolute.parent_path().empty()) {
        fs::create_directories(absolute.parent_path());
    }
    const fs::path temporary = temporaryPath(absolute);
    struct Guard {
        fs::path path;
        bool keep{false};
        ~Guard() {
            if (!keep) {
                std::error_code error;
                fs::remove(path, error);
            }
        }
    } guard{temporary};
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error("Cannot open temporary output: " +
                                 temporary.string());
    }
    writer(output);
    output.close();
    if (!output) {
        throw std::runtime_error("Failed to finalize output: " +
                                 absolute.string());
    }
    std::error_code rename_error;
    fs::rename(temporary, absolute, rename_error);
    if (rename_error) {
        throw std::runtime_error("Failed to atomically replace output: " +
                                 rename_error.message());
    }
    guard.keep = true;
}

void writeGraph(const std::filesystem::path& path,
                Version version,
                const std::string& reference_species,
                const std::vector<SourceSequence>& sequences,
                const GraphIR& graph) {
    atomicWrite(path, [&](std::ostream& output) {
        if (version == Version::V1_0) {
            output << "H\tVN:Z:1.0\n";
        } else {
            output << "H\tVN:Z:1.1\tRS:Z:" << reference_species << '\n';
        }
        for (const auto& node : graph.nodes) {
            output << "S\t" << node.id << '\t' << node.sequence << '\n';
        }
        for (const auto& edge : graph.edges) {
            output << "L\t" << edge.from << '\t'
                   << (edge.from_reverse ? '-' : '+') << '\t'
                   << edge.to << '\t' << (edge.to_reverse ? '-' : '+')
                   << "\t0M\n";
        }
        for (std::size_t index = 0; index < sequences.size(); ++index) {
            const auto& source = sequences[index];
            const auto& walk = graph.walks[index];
            if (version == Version::V1_0) {
                output << "P\t"
                       << Export::qualifiedName(source.species,
                                                source.chromosome)
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
    });
}

std::string spellInterval(const std::vector<Node>& nodes,
                          const std::vector<Step>& walk,
                          std::size_t begin,
                          std::size_t end) {
    std::string sequence;
    for (std::size_t index = begin; index < end; ++index) {
        sequence += orientedNodeSequence(
            nodes[walk[index].node - 1], walk[index].reverse);
    }
    return sequence;
}

bool hasN(const std::string& sequence) {
    return sequence.find('N') != std::string::npos;
}

std::optional<IslandCandidate> analyzeIsland(
    const GraphIR& graph,
    const std::vector<std::vector<Occurrence>>& occurrences,
    std::size_t reference_walk,
    std::size_t left_position,
    std::size_t right_position,
    std::string& rejection) {
    const auto& reference = graph.walks[reference_walk];
    if (right_position <= left_position + 1) {
        rejection = "empty_reference_interior";
        return std::nullopt;
    }
    const Step& left_anchor = reference[left_position];
    const Step& right_anchor = reference[right_position];
    if (left_anchor.reverse || right_anchor.reverse) {
        rejection = "anchor_orientation";
        return std::nullopt;
    }
    const auto reference_span =
        right_anchor.source_start - left_anchor.source_end;
    if (reference_span > COMPACT_MAX_REFERENCE_SPAN) {
        rejection = "reference_span";
        return std::nullopt;
    }

    IslandCandidate candidate;
    candidate.reference_walk = reference_walk;
    candidate.reference_begin = left_position;
    candidate.reference_end = right_position;
    for (std::size_t walk_index = 0; walk_index < graph.walks.size();
         ++walk_index) {
        const auto& walk = graph.walks[walk_index];
        std::optional<std::size_t> left;
        std::optional<std::size_t> right;
        for (std::size_t position = 0; position < walk.size(); ++position) {
            if (handleEquals(walk[position],
                             {left_anchor.node, left_anchor.reverse})) {
                if (left) {
                    rejection = "repeated_left_anchor";
                    return std::nullopt;
                }
                left = position;
            }
            if (handleEquals(walk[position],
                             {right_anchor.node, right_anchor.reverse})) {
                if (right) {
                    rejection = "repeated_right_anchor";
                    return std::nullopt;
                }
                right = position;
            }
        }
        if (left.has_value() != right.has_value()) {
            rejection = "unpaired_anchor";
            return std::nullopt;
        }
        if (!left) continue;
        if (*right <= *left) {
            rejection = "anchor_order_or_inversion";
            return std::nullopt;
        }
        const auto source_span =
            walk[*right].source_start - walk[*left].source_end;
        if (source_span > COMPACT_MAX_PATH_SPAN) {
            rejection = "path_span";
            return std::nullopt;
        }
        std::set<std::uint64_t> local_nodes;
        for (std::size_t position = *left + 1; position < *right; ++position) {
            const auto& step = walk[position];
            if (step.reverse || !local_nodes.insert(step.node).second) {
                rejection = "orientation_or_repeated_node";
                return std::nullopt;
            }
            if (hasN(graph.nodes[step.node - 1].sequence)) {
                rejection = "contains_N";
                return std::nullopt;
            }
            candidate.interior_nodes.insert(step.node);
        }
        const std::string allele = spellInterval(
            graph.nodes, walk, *left + 1, *right);
        auto allele_it = std::find_if(
            candidate.alleles.begin(), candidate.alleles.end(),
            [&](const auto& value) { return value.first == allele; });
        std::size_t allele_index = 0;
        if (allele_it == candidate.alleles.end()) {
            allele_index = candidate.alleles.size();
            candidate.alleles.push_back({allele, {walk_index}});
        } else {
            allele_index = static_cast<std::size_t>(
                std::distance(candidate.alleles.begin(), allele_it));
            allele_it->second.push_back(walk_index);
        }
        candidate.replacements.push_back(
            {walk_index, *left + 1, *right, allele_index + 1});
    }
    if (candidate.alleles.size() > COMPACT_MAX_ALLELES) {
        rejection = "allele_limit";
        return std::nullopt;
    }
    if (candidate.replacements.size() < 2) {
        rejection = "single_path";
        return std::nullopt;
    }

    std::size_t short_nodes = 0;
    for (const auto node_id : candidate.interior_nodes) {
        const auto& node = graph.nodes[node_id - 1];
        if (node.sequence.size() < COMPACT_SHORT_NODE) ++short_nodes;
        candidate.old_graph_bp += node.sequence.size();
        if (node.occurrence_count > 1) {
            candidate.old_homology_mass +=
                node.sequence.size() * (node.occurrence_count - 1);
        }
        for (const auto& occurrence : occurrences[node_id]) {
            const auto replacement = std::find_if(
                candidate.replacements.begin(), candidate.replacements.end(),
                [&](const IslandReplacement& value) {
                    return value.walk == occurrence.first &&
                           occurrence.second >= value.begin &&
                           occurrence.second < value.end;
                });
            if (replacement == candidate.replacements.end()) {
                rejection = "external_interior_occurrence";
                return std::nullopt;
            }
        }
    }
    std::map<std::uint64_t, std::set<std::uint64_t>> adjacency;
    std::map<std::uint64_t, std::size_t> indegree;
    for (const auto node_id : candidate.interior_nodes) indegree[node_id] = 0;
    for (const auto& replacement : candidate.replacements) {
        const auto& walk = graph.walks[replacement.walk];
        for (std::size_t position = replacement.begin + 1;
             position < replacement.end; ++position) {
            const auto left = walk[position - 1].node;
            const auto right = walk[position].node;
            if (left == right) {
                rejection = "cycle";
                return std::nullopt;
            }
            if (adjacency[left].insert(right).second) ++indegree[right];
        }
    }
    std::deque<std::uint64_t> ready;
    for (const auto& [node, degree] : indegree) {
        if (degree == 0) ready.push_back(node);
    }
    std::size_t visited = 0;
    while (!ready.empty()) {
        const auto node = ready.front();
        ready.pop_front();
        ++visited;
        for (const auto next : adjacency[node]) {
            if (--indegree[next] == 0) ready.push_back(next);
        }
    }
    if (visited != candidate.interior_nodes.size()) {
        rejection = "cycle";
        return std::nullopt;
    }
    candidate.old_short_nodes = short_nodes;
    if (short_nodes < COMPACT_MIN_SHORT_NODES) {
        rejection = "too_few_short_nodes";
        return std::nullopt;
    }
    for (const auto& [allele, member_walks] : candidate.alleles) {
        if (allele.empty()) continue;
        candidate.new_graph_bp += allele.size();
        if (allele.size() < COMPACT_SHORT_NODE) ++candidate.new_short_nodes;
        if (member_walks.size() > 1) {
            candidate.new_homology_mass +=
                allele.size() * (member_walks.size() - 1);
        }
    }
    const std::size_t new_nodes = std::count_if(
        candidate.alleles.begin(), candidate.alleles.end(),
        [](const auto& allele) { return !allele.first.empty(); });
    if (new_nodes >= candidate.interior_nodes.size() ||
        candidate.new_short_nodes >= candidate.old_short_nodes) {
        rejection = "no_fragmentation_gain";
        return std::nullopt;
    }
    return candidate;
}

std::size_t rewriteAlleleIslands(
    GraphIR& graph,
    const std::vector<SourceSequence>& sequences,
    const std::string& reference_species,
    std::uint64_t exact_homology_mass,
    std::uint64_t exact_graph_bp,
    int threads,
    std::vector<std::string>& rejections) {
    struct Request {
        std::size_t walk{0};
        std::size_t left{0};
        std::size_t right{0};
    };
    std::vector<std::vector<std::pair<std::uint64_t, Occurrence>>>
        occurrences_by_walk(graph.walks.size());
#pragma omp parallel for schedule(static) num_threads(threads)
    for (std::int64_t walk_index = 0;
         walk_index < static_cast<std::int64_t>(graph.walks.size());
         ++walk_index) {
        const auto& walk = graph.walks[static_cast<std::size_t>(walk_index)];
        auto& indexed = occurrences_by_walk[
            static_cast<std::size_t>(walk_index)];
        indexed.reserve(walk.size());
        for (std::size_t position = 0; position < walk.size(); ++position) {
            indexed.push_back({walk[position].node,
                               {static_cast<std::size_t>(walk_index),
                                position}});
        }
    }
    std::vector<std::vector<Occurrence>> occurrences(graph.nodes.size() + 1);
    for (const auto& indexed : occurrences_by_walk) {
        for (const auto& [node, occurrence] : indexed) {
            occurrences[node].push_back(occurrence);
        }
    }

    std::vector<Request> requests;
    for (std::size_t walk_index = 0; walk_index < sequences.size();
         ++walk_index) {
        if (sequences[walk_index].species != reference_species) continue;
        const auto& walk = graph.walks[walk_index];
        std::unordered_map<std::uint64_t, std::size_t> counts;
        for (const auto& step : walk) ++counts[step.node];
        std::vector<std::size_t> anchors;
        for (std::size_t position = 0; position < walk.size(); ++position) {
            const auto& node = graph.nodes[walk[position].node - 1];
            if (!walk[position].reverse &&
                node.sequence.size() >= COMPACT_STABLE_ANCHOR &&
                node.occurrence_count > 1 && counts[node.id] == 1) {
                anchors.push_back(position);
            }
        }
        for (std::size_t index = 1; index < anchors.size(); ++index) {
            requests.push_back(
                {walk_index, anchors[index - 1], anchors[index]});
        }
    }
    std::vector<std::optional<IslandCandidate>> analyses(requests.size());
    std::vector<std::string> analysis_rejections(requests.size());
#pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (std::int64_t index = 0;
         index < static_cast<std::int64_t>(requests.size()); ++index) {
        const auto& request = requests[static_cast<std::size_t>(index)];
        analyses[static_cast<std::size_t>(index)] = analyzeIsland(
            graph, occurrences, request.walk, request.left, request.right,
            analysis_rejections[static_cast<std::size_t>(index)]);
    }

    std::vector<IslandCandidate> selected;
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> occupied(
        graph.walks.size());
    std::uint64_t projected_homology = collectMetrics(
        graph.nodes, graph.walks, graph.edges).homology_mass;
    std::uint64_t projected_graph_bp = collectMetrics(
        graph.nodes, graph.walks, graph.edges).graph_bp;

    for (std::size_t index = 0; index < requests.size(); ++index) {
            const auto& request = requests[index];
            auto& candidate = analyses[index];
            if (!candidate) {
                rejections.push_back(
                    std::to_string(request.walk) + "\t" +
                    std::to_string(request.left) + "\t" +
                    std::to_string(request.right) + "\t" +
                    analysis_rejections[index]);
                continue;
            }
            bool overlaps = false;
            for (const auto& replacement : candidate->replacements) {
                for (const auto& interval : occupied[replacement.walk]) {
                    if (replacement.begin < interval.second &&
                        interval.first < replacement.end) {
                        overlaps = true;
                    }
                }
            }
            if (overlaps) {
                rejections.push_back(
                    std::to_string(request.walk) + "\t" +
                    std::to_string(request.left) + "\t" +
                    std::to_string(request.right) +
                    "\toverlapping_candidate");
                continue;
            }
            const auto next_homology = projected_homology -
                candidate->old_homology_mass +
                candidate->new_homology_mass;
            const auto next_graph_bp = projected_graph_bp -
                candidate->old_graph_bp + candidate->new_graph_bp;
            const double retention = exact_homology_mass == 0 ? 1.0 :
                static_cast<double>(next_homology) /
                static_cast<double>(exact_homology_mass);
            const double bp_ratio = exact_graph_bp == 0 ? 1.0 :
                static_cast<double>(next_graph_bp) /
                static_cast<double>(exact_graph_bp);
            if (retention < COMPACT_MIN_HOMOLOGY_RETENTION) {
                rejections.push_back(
                    std::to_string(request.walk) + "\t" +
                    std::to_string(request.left) + "\t" +
                    std::to_string(request.right) +
                    "\thomology_budget");
                continue;
            }
            if (bp_ratio > COMPACT_MAX_GRAPH_BP_RATIO) {
                rejections.push_back(
                    std::to_string(request.walk) + "\t" +
                    std::to_string(request.left) + "\t" +
                    std::to_string(request.right) + "\tgraph_bp_budget");
                continue;
            }
            projected_homology = next_homology;
            projected_graph_bp = next_graph_bp;
            for (const auto& replacement : candidate->replacements) {
                occupied[replacement.walk].push_back(
                    {replacement.begin, replacement.end});
            }
            selected.push_back(std::move(*candidate));
    }
    if (selected.empty()) return 0;

    std::vector<std::vector<IslandReplacement>> replacements(graph.walks.size());
    for (auto& candidate : selected) {
        std::vector<std::uint64_t> allele_nodes(candidate.alleles.size(), 0);
        for (std::size_t allele = 0; allele < candidate.alleles.size();
             ++allele) {
            const auto& sequence = candidate.alleles[allele].first;
            if (sequence.empty()) continue;
            const auto first_walk = candidate.alleles[allele].second.front();
            const auto replacement = std::find_if(
                candidate.replacements.begin(), candidate.replacements.end(),
                [&](const IslandReplacement& value) {
                    return value.walk == first_walk &&
                           value.replacement_node == allele + 1;
                });
            const auto source_start =
                graph.walks[first_walk][replacement->begin].source_start;
            const std::uint64_t id = graph.nodes.size() + 1;
            graph.nodes.push_back({id, 0, sequence, 0,
                                   first_walk, source_start});
            allele_nodes[allele] = id;
        }
        for (auto replacement : candidate.replacements) {
            const auto allele = replacement.replacement_node - 1;
            replacement.replacement_node = allele_nodes[allele];
            replacements[replacement.walk].push_back(replacement);
        }
    }
    for (std::size_t walk_index = 0; walk_index < graph.walks.size();
         ++walk_index) {
        auto& walk_replacements = replacements[walk_index];
        std::sort(walk_replacements.begin(), walk_replacements.end(),
                  [](const auto& left, const auto& right) {
                      return left.begin < right.begin;
                  });
        const auto old_walk = graph.walks[walk_index];
        std::vector<Step> rewritten;
        std::size_t position = 0;
        for (const auto& replacement : walk_replacements) {
            rewritten.insert(rewritten.end(), old_walk.begin() + position,
                             old_walk.begin() + replacement.begin);
            if (replacement.replacement_node != 0) {
                rewritten.push_back({replacement.replacement_node, false,
                                     old_walk[replacement.begin].source_start,
                                     old_walk[replacement.end - 1].source_end});
            }
            position = replacement.end;
        }
        rewritten.insert(rewritten.end(), old_walk.begin() + position,
                         old_walk.end());
        graph.walks[walk_index] = std::move(rewritten);
    }
    normalizeGraph(graph.nodes, graph.walks, graph.edges);
    return selected.size();
}

void fillStatsFromGraph(GfaExportStats& stats,
                        const GraphIR& graph,
                        const GraphMetrics& metrics) {
    stats.nodes = metrics.nodes;
    stats.edges = metrics.edges;
    stats.walks = graph.walks.size();
    stats.graph_sequence_bp = metrics.graph_bp;
    stats.homology_mass = metrics.homology_mass;
    stats.shared_sequence_bp = 0;
    stats.private_sequence_bp = 0;
    for (const auto& node : graph.nodes) {
        if (node.occurrence_count > 1) {
            stats.shared_sequence_bp += node.sequence.size();
        } else {
            stats.private_sequence_bp += node.sequence.size();
        }
    }
}

void writeCompactReports(
    const std::filesystem::path& directory,
    const std::vector<SourceSequence>& sequences,
    const GraphIR& exact,
    const GraphIR& compact,
    const GraphMetrics& exact_metrics,
    const GraphMetrics& compact_metrics,
    const std::vector<std::string>& rejections,
    const GfaExportStats& stats) {
    std::filesystem::create_directories(directory);
    atomicWrite(directory / "compact_parameters.tsv", [&](std::ostream& out) {
        out << "parameter\tvalue\n"
            << "profile\tcompact-v1\n"
            << "minimum_exact_run_bp\t" << COMPACT_MIN_EXACT_RUN << '\n'
            << "short_node_bp\t" << COMPACT_SHORT_NODE << '\n'
            << "stable_anchor_bp\t" << COMPACT_STABLE_ANCHOR << '\n'
            << "maximum_reference_span_bp\t" << COMPACT_MAX_REFERENCE_SPAN << '\n'
            << "maximum_path_span_bp\t" << COMPACT_MAX_PATH_SPAN << '\n'
            << "minimum_short_nodes\t" << COMPACT_MIN_SHORT_NODES << '\n'
            << "maximum_distinct_alleles\t" << COMPACT_MAX_ALLELES << '\n'
            << "minimum_homology_retention\t" << COMPACT_MIN_HOMOLOGY_RETENTION << '\n'
            << "maximum_graph_bp_ratio\t" << COMPACT_MAX_GRAPH_BP_RATIO << '\n';
    });
    atomicWrite(directory / "compact_stats.tsv", [&](std::ostream& out) {
        out << "metric\texact\tcompact\n"
            << "nodes\t" << exact_metrics.nodes << '\t' << compact_metrics.nodes << '\n'
            << "edges\t" << exact_metrics.edges << '\t' << compact_metrics.edges << '\n'
            << "steps\t" << exact_metrics.steps << '\t' << compact_metrics.steps << '\n'
            << "nodes_le_10\t" << exact_metrics.nodes_le_10 << '\t' << compact_metrics.nodes_le_10 << '\n'
            << "nodes_lt_25\t" << exact_metrics.nodes_lt_25 << '\t' << compact_metrics.nodes_lt_25 << '\n'
            << "nodes_lt_100\t" << exact_metrics.nodes_lt_100 << '\t' << compact_metrics.nodes_lt_100 << '\n'
            << "median_node_length\t" << exact_metrics.median_node_length << '\t' << compact_metrics.median_node_length << '\n'
            << "graph_sequence_bp\t" << exact_metrics.graph_bp << '\t' << compact_metrics.graph_bp << '\n'
            << "homology_mass\t" << exact_metrics.homology_mass << '\t' << compact_metrics.homology_mass << '\n'
            << "suppressed_exact_runs\t0\t" << stats.suppressed_exact_runs << '\n'
            << "unitig_merges\t0\t" << stats.unitig_merges << '\n'
            << "allele_islands\t0\t" << stats.allele_islands << '\n';
    });
    atomicWrite(directory / "compact_rejections.tsv", [&](std::ostream& out) {
        out << "reference_walk\tleft_step\tright_step\treason\n";
        for (const auto& rejection : rejections) out << rejection << '\n';
    });
    atomicWrite(directory / "compact_transform.tsv", [&](std::ostream& out) {
        out << "species\tcontig\tstart\tend\texact_node\tcompact_node\t"
               "exact_reverse\tcompact_reverse\n";
        for (std::size_t sequence = 0; sequence < sequences.size(); ++sequence) {
            const auto& exact_walk = exact.walks[sequence];
            const auto& compact_walk = compact.walks[sequence];
            std::size_t left = 0;
            std::size_t right = 0;
            while (left < exact_walk.size() && right < compact_walk.size()) {
                const auto start = std::max(exact_walk[left].source_start,
                                            compact_walk[right].source_start);
                const auto end = std::min(exact_walk[left].source_end,
                                          compact_walk[right].source_end);
                if (start < end) {
                    out << sequences[sequence].species << '\t'
                        << sequences[sequence].chromosome << '\t'
                        << start << '\t' << end << '\t'
                        << exact_walk[left].node << '\t'
                        << compact_walk[right].node << '\t'
                        << (exact_walk[left].reverse ? 1 : 0) << '\t'
                        << (compact_walk[right].reverse ? 1 : 0) << '\n';
                }
                const auto exact_end = exact_walk[left].source_end;
                const auto compact_end = compact_walk[right].source_end;
                if (exact_end <= compact_end) ++left;
                if (compact_end <= exact_end) ++right;
            }
        }
    });
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

    const std::string reference_species = !graph.reference_order.empty()
        ? graph.reference_order.front() : managers.begin()->first;
    if (!validWalkField(reference_species)) {
        throw std::runtime_error(
            "Invalid GFA reference sample identifier: " + reference_species);
    }

    validatePathNames(sequences, options.version);
    GraphIR exact = buildInitialGraph(sequences, relations);
    const GraphMetrics exact_metrics = collectMetrics(
        exact.nodes, exact.walks, exact.edges);
    stats.atomic_intervals = exact.atomic_intervals;
    stats.exact_homology_mass = exact_metrics.homology_mass;

    GraphIR result = exact;
    GraphMetrics result_metrics = exact_metrics;
    std::vector<std::string> rejections;
    if (options.profile == Profile::COMPACT) {
        if (options.work_dir.empty()) {
            throw std::runtime_error(
                "Compact GFA export requires a non-empty work directory");
        }
        const fs::path report_dir = fs::absolute(options.work_dir) / "gfa";
        const fs::path exact_shadow = report_dir / "exact.gfa";
        if (fs::absolute(gfa_path) == exact_shadow) {
            throw std::runtime_error(
                "Compact GFA output must differ from work/gfa/exact.gfa");
        }
        writeGraph(exact_shadow, options.version, reference_species,
                   sequences, exact);

        std::vector<unsigned char> keep(relations.size(), 0);
#pragma omp parallel for schedule(static) num_threads(export_threads)
        for (std::int64_t index = 0;
             index < static_cast<std::int64_t>(relations.size()); ++index) {
            keep[static_cast<std::size_t>(index)] =
                relations[static_cast<std::size_t>(index)].length >=
                COMPACT_MIN_EXACT_RUN;
        }
        std::vector<ExactRelation> retained;
        retained.reserve(relations.size());
        for (std::size_t index = 0; index < relations.size(); ++index) {
            if (!keep[index]) {
                ++stats.suppressed_exact_runs;
            } else {
                retained.push_back(relations[index]);
            }
        }
        result = buildInitialGraph(sequences, retained);
        auto initial_metrics = collectMetrics(
            result.nodes, result.walks, result.edges);
        const double initial_homology_retention =
            exact_metrics.homology_mass == 0 ? 1.0 :
            static_cast<double>(initial_metrics.homology_mass) /
            static_cast<double>(exact_metrics.homology_mass);
        const double initial_graph_bp_ratio = exact_metrics.graph_bp == 0 ? 1.0 :
            static_cast<double>(initial_metrics.graph_bp) /
            static_cast<double>(exact_metrics.graph_bp);
        if (initial_homology_retention < COMPACT_MIN_HOMOLOGY_RETENTION) {
            throw std::runtime_error(
                "Compact short-relation suppression violates the 99.5% "
                "homology-mass budget; exact shadow was retained");
        }
        if (initial_graph_bp_ratio > COMPACT_MAX_GRAPH_BP_RATIO) {
            throw std::runtime_error(
                "Compact short-relation suppression violates the 2% graph "
                "sequence budget; exact shadow was retained");
        }
        stats.unitig_merges += compactUnitigs(
            result.nodes, result.walks, result.edges);
        validateGraph(sequences, result.nodes, result.walks, result.edges);
        stats.allele_islands = rewriteAlleleIslands(
            result, sequences, reference_species,
            exact_metrics.homology_mass, exact_metrics.graph_bp,
            export_threads, rejections);
        stats.rejected_allele_islands = rejections.size();
        validateGraph(sequences, result.nodes, result.walks, result.edges);
        stats.unitig_merges += compactUnitigs(
            result.nodes, result.walks, result.edges);
        validateGraph(sequences, result.nodes, result.walks, result.edges);
        result_metrics = collectMetrics(
            result.nodes, result.walks, result.edges);
        const double final_homology_retention =
            exact_metrics.homology_mass == 0 ? 1.0 :
            static_cast<double>(result_metrics.homology_mass) /
            static_cast<double>(exact_metrics.homology_mass);
        const double final_graph_bp_ratio = exact_metrics.graph_bp == 0 ? 1.0 :
            static_cast<double>(result_metrics.graph_bp) /
            static_cast<double>(exact_metrics.graph_bp);
        if (final_homology_retention < COMPACT_MIN_HOMOLOGY_RETENTION ||
            final_graph_bp_ratio > COMPACT_MAX_GRAPH_BP_RATIO) {
            throw std::runtime_error(
                "Compact transforms exceeded a global safety budget; exact "
                "shadow was retained");
        }
        fillStatsFromGraph(stats, result, result_metrics);
        writeCompactReports(report_dir, sequences, exact, result,
                            exact_metrics, result_metrics, rejections, stats);
    } else {
        fillStatsFromGraph(stats, result, result_metrics);
    }

    writeGraph(gfa_path, options.version, reference_species,
               sequences, result);

    stats.elapsed_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    spdlog::info(
        "GFA {} {} export complete: blocks(seen/used)={}/{}, exact_runs={}, "
        "suppressed_runs={}, atoms={}, nodes={}, edges={}, paths={}, graph_bp={}, "
        "shared_bp={}, private_bp={}, homology_mass={}/{}, unitig_merges={}, "
        "allele_islands={}, elapsed_seconds={:.3f}",
        versionString(options.version), profileString(options.profile),
        stats.blocks_seen, stats.blocks_used,
        stats.exact_runs, stats.suppressed_exact_runs,
        stats.atomic_intervals, stats.nodes, stats.edges,
        stats.walks, stats.graph_sequence_bp,
        stats.shared_sequence_bp, stats.private_sequence_bp,
        stats.homology_mass, stats.exact_homology_mass,
        stats.unitig_merges, stats.allele_islands, stats.elapsed_seconds);
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
