// Segment conversion adapted from CactusHalConverter.
// Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu).
// MIT license; see THIRD_PARTY_NOTICES.md for the retained notice.

#include "subtree_writer.h"
#include "ramaxHdf5BulkDna.h"

#include "hal.h"
#include <H5Cpp.h>

#include <algorithm>
#include <exception>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace RaMesh::hal_export {
namespace {

constexpr uint64_t kDnaChunkSize = 1024ULL * 1024ULL;
std::runtime_error hdf5Failure(std::string_view operation,
                               const H5::Exception& error) {
    const std::string detail = error.getDetailMsg();
    if (detail.empty()) {
        return std::runtime_error(std::string(operation));
    }
    return std::runtime_error(std::string(operation) + ": " + detail);
}

hal::AlignmentPtr openNativeAlignment(const std::filesystem::path& path) {
    try {
        return hal::AlignmentPtr(hal::hdf5AlignmentInstance(
            path.string(),
            hal::READ_ACCESS | hal::WRITE_ACCESS | hal::CREATE_ACCESS,
            hal::hdf5DefaultFileCreatPropList(),
            hal::hdf5DefaultFileAccPropList(),
            hal::hdf5DefaultDSetCreatPropList(), true));
    } catch (const H5::Exception& error) {
        throw hdf5Failure("Unable to initialize native HAL writer", error);
    }
}

void discardAlignmentNoThrow(hal::AlignmentPtr& alignment) noexcept {
    if (!alignment) {
        return;
    }
    try {
        alignment->close();
    } catch (...) {
    }
    alignment.reset();
}

struct PreparedSequence {
    const SequenceEmission* emission = nullptr;
    hal_size_t length = 0;
};

using GenomeOrder = std::vector<const std::string*>;
using PreparedByGenome = std::vector<std::vector<PreparedSequence>>;

struct BottomState {
    hal_index_t array_index = hal::NULL_INDEX;
    hal_index_t duplicate_tail = hal::NULL_INDEX;
};
static_assert(sizeof(BottomState) == 2 * sizeof(hal_index_t));

class GenomeCloser {
public:
    explicit GenomeCloser(hal::Alignment* alignment) : alignment_(alignment) {}

    ~GenomeCloser() {
        closeNoThrow();
    }

    void add(hal::Genome* genome) {
        if (genome != nullptr) {
            genomes_.push_back(genome);
        }
    }

    void closeGenome(hal::Genome* genome) {
        const auto it = std::find(genomes_.begin(), genomes_.end(), genome);
        if (it == genomes_.end() || *it == nullptr) {
            throw std::logic_error(
                "Attempted to close an untracked HAL genome");
        }
        alignment_->closeGenome(*it);
        *it = nullptr;
    }

    void close() {
        std::exception_ptr first_failure;
        for (auto it = genomes_.rbegin(); it != genomes_.rend(); ++it) {
            if (*it == nullptr) {
                continue;
            }
            try {
                alignment_->closeGenome(*it);
                *it = nullptr;
            } catch (...) {
                if (!first_failure) {
                    first_failure = std::current_exception();
                }
            }
        }
        if (first_failure) {
            std::rethrow_exception(first_failure);
        }
        genomes_.clear();
    }

private:
    void closeNoThrow() noexcept {
        for (auto it = genomes_.rbegin(); it != genomes_.rend(); ++it) {
            if (*it == nullptr) {
                continue;
            }
            try {
                alignment_->closeGenome(*it);
            } catch (...) {
            }
        }
        genomes_.clear();
    }

    hal::Alignment* alignment_ = nullptr;
    std::vector<hal::Genome*> genomes_;
};

hal_size_t checkedAdd(hal_size_t lhs, hal_size_t rhs,
                      const std::string& context) {
    if (rhs > std::numeric_limits<hal_size_t>::max() - lhs) {
        throw std::overflow_error("HAL dimension overflow for " + context);
    }
    return lhs + rhs;
}

hal_index_t absoluteStart(const hal::Sequence& sequence, uint64_t local_start,
                          const std::string& context) {
    const hal_index_t sequence_start = sequence.getStartPosition();
    if (sequence_start < 0 ||
        local_start > static_cast<uint64_t>(
                          std::numeric_limits<hal_index_t>::max() -
                          sequence_start)) {
        throw std::overflow_error("HAL coordinate overflow for " + context);
    }
    return sequence_start + static_cast<hal_index_t>(local_start);
}

double roundedBranchLength(double branch_length) {
    // buildLocalNewick historically serialized this value at fixed precision
    // before cactus2hal parsed it back to double.  Preserve that roundtrip.
    std::ostringstream out;
    out << std::fixed << std::setprecision(6) << branch_length;
    size_t consumed = 0;
    const std::string text = out.str();
    const double rounded = std::stod(text, &consumed);
    if (consumed != text.size()) {
        throw std::runtime_error("Invalid HAL branch length");
    }
    return rounded;
}


const SeqPro::SequenceManager& originalManager(
    const SeqPro::SharedManagerVariant& shared_manager) {
    if (!shared_manager) {
        throw std::runtime_error("Null SeqPro manager in HAL export");
    }
    return std::visit(
        [](const auto& manager) -> const SeqPro::SequenceManager& {
            using Ptr = std::decay_t<decltype(manager)>;
            if (!manager) {
                throw std::runtime_error("Null SeqPro manager in HAL export");
            }
            if constexpr (
                std::is_same_v<Ptr,
                               std::unique_ptr<SeqPro::SequenceManager>>) {
                return *manager;
            } else if constexpr (
                std::is_same_v<Ptr,
                               std::unique_ptr<
                                   SeqPro::MaskedSequenceManager>>) {
                return manager->getOriginalManager();
            } else {
                static_assert(sizeof(Ptr) == 0,
                              "Unsupported SeqPro manager variant");
            }
        },
        *shared_manager);
}

void writeBorrowedDna(hal::Sequence& sequence, std::string_view dna) {
    if (dna.size() != sequence.getSequenceLength()) {
        throw std::runtime_error(
            "HAL DNA length does not match segment dimensions for " +
            sequence.getFullName());
    }
    hal::ramaxWriteBulkDna(
        *sequence.getGenome(), sequence.getStartPosition(), dna.data(),
        static_cast<hal_size_t>(dna.size()));
}

void writeLeafDna(
    hal::Sequence& sequence,
    const std::pair<std::string, std::string>& source,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    const SoftMask::IndexMap& softmask_indexes) {
    const auto manager_it = managers.find(source.first);
    if (manager_it == managers.end()) {
        throw std::runtime_error(
            "Missing SeqPro manager for leaf genome: " + source.first);
    }
    const auto softmask_it = softmask_indexes.find(source.first);
    if (softmask_it == softmask_indexes.end() || !softmask_it->second) {
        throw std::runtime_error(
            "Missing soft-mask index for leaf genome: " + source.first);
    }

    const SeqPro::SequenceManager& manager =
        originalManager(manager_it->second);
    const auto sequence_id = manager.getSequenceId(source.second);
    const uint64_t source_length = manager.getSequenceLength(sequence_id);
    if (source_length != sequence.getSequenceLength()) {
        throw std::runtime_error(
            "HAL leaf DNA length does not match segment dimensions for " +
            sequence.getFullName());
    }
    std::string buffer;
    buffer.reserve(static_cast<size_t>(
        std::min<uint64_t>(kDnaChunkSize, source_length)));
    for (uint64_t offset = 0; offset < source_length;) {
        const uint64_t length =
            std::min<uint64_t>(kDnaChunkSize, source_length - offset);
        manager.getSubSequenceInto(sequence_id, offset, length, buffer);
        softmask_it->second->restore(source.second, offset, buffer);
        hal::ramaxWriteBulkDna(
            *sequence.getGenome(),
            absoluteStart(sequence, offset, sequence.getFullName()),
            buffer.data(), static_cast<hal_size_t>(buffer.size()));
        offset += length;
    }
    manager.releaseMappedPages();
}

PreparedByGenome prepareEmissions(
    const GenomeOrder& genome_order,
    const std::vector<SequenceEmission>& emissions) {
    PreparedByGenome prepared(genome_order.size());

    for (const SequenceEmission& emission : emissions) {
        const auto genome_it =
            std::find_if(genome_order.begin(), genome_order.end(),
                         [&emission](const std::string* genome_name) {
                             return *genome_name == emission.genome_name;
                         });
        if (genome_it == genome_order.end()) {
            throw std::runtime_error(
                "Emission genome is outside the local HAL subtree: " +
                emission.genome_name);
        }
        if (emission.bottom_count != emission.bottoms.size()) {
            throw std::runtime_error(
                "HAL bottom count does not match emitted bottom segments for " +
                emission.genome_name + "." + emission.seq_name);
        }
        if (!emission.bottoms.empty() && !emission.tops.empty()) {
            throw std::runtime_error(
                "A local HAL sequence cannot emit top and bottom segments "
                "simultaneously: " +
                emission.genome_name + "." + emission.seq_name);
        }
        if (emission.dna && emission.leaf_source) {
            throw std::runtime_error(
                "HAL sequence has both borrowed and leaf-source DNA: " +
                emission.genome_name + "." + emission.seq_name);
        }

        const std::string context =
            emission.genome_name + "." + emission.seq_name;
        hal_size_t length = 0;
        if (!emission.bottoms.empty()) {
            for (const BottomSegmentLine& bottom : emission.bottoms) {
                length = checkedAdd(length, bottom.length, context);
            }
        } else {
            for (const TopSegmentLine& top : emission.tops) {
                length = checkedAdd(length, top.length, context);
            }
        }
        prepared[static_cast<size_t>(genome_it - genome_order.begin())]
            .push_back(PreparedSequence{&emission, length});
    }

    for (auto& sequences : prepared) {
        std::sort(sequences.begin(), sequences.end(),
                  [](const PreparedSequence& lhs,
                     const PreparedSequence& rhs) {
                      return lhs.emission->seq_name <
                             rhs.emission->seq_name;
                  });
        for (size_t index = 1; index < sequences.size(); ++index) {
            const SequenceEmission& previous =
                *sequences[index - 1].emission;
            const SequenceEmission& current = *sequences[index].emission;
            if (previous.seq_name == current.seq_name) {
                throw std::runtime_error(
                    "Duplicate HAL sequence emission: " +
                    current.genome_name + "." + current.seq_name);
            }
        }
    }
    return prepared;
}

std::vector<hal::Sequence::Info> makeDimensions(
    const std::vector<PreparedSequence>& sequences, bool bottom_genome) {
    std::vector<hal::Sequence::Info> dimensions;
    dimensions.reserve(sequences.size());
    for (const PreparedSequence& prepared : sequences) {
        const SequenceEmission& emission = *prepared.emission;
        if (bottom_genome && !emission.tops.empty()) {
            throw std::runtime_error(
                "Local subtree root emitted top segments for " +
                emission.genome_name + "." + emission.seq_name);
        }
        if (!bottom_genome && !emission.bottoms.empty()) {
            throw std::runtime_error(
                "Local subtree child emitted bottom segments for " +
                emission.genome_name + "." + emission.seq_name);
        }
        dimensions.emplace_back(
            emission.seq_name, prepared.length,
            bottom_genome ? 0 : emission.tops.size(),
            bottom_genome ? emission.bottoms.size() : 0);
    }
    return dimensions;
}

void setNewGenomeDimensionsAndDna(
    hal::Genome& genome,
    const std::vector<PreparedSequence>& sequences,
    bool bottom_genome,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    const SoftMask::IndexMap& softmask_indexes) {
    {
        const auto dimensions = makeDimensions(sequences, bottom_genome);
        genome.setDimensions(dimensions);
    }

    for (const PreparedSequence& prepared : sequences) {
        if (prepared.length == 0) {
            continue;
        }
        const SequenceEmission& emission = *prepared.emission;
        hal::Sequence* sequence = genome.getSequence(emission.seq_name);
        if (sequence == nullptr) {
            throw std::runtime_error(
                "HAL did not create sequence " + emission.genome_name + "." +
                emission.seq_name);
        }
        if (emission.dna) {
            writeBorrowedDna(*sequence, *emission.dna);
        } else if (emission.leaf_source) {
            writeLeafDna(*sequence, *emission.leaf_source, managers,
                         softmask_indexes);
        } else {
            throw std::runtime_error(
                "Missing DNA source for non-empty HAL sequence " +
                emission.genome_name + "." + emission.seq_name);
        }
    }
}

void updateExistingRootDimensions(
    hal::Genome& genome,
    const std::vector<PreparedSequence>& sequences) {
    std::vector<hal::Sequence::UpdateInfo> updates;
    updates.reserve(sequences.size());
    for (const PreparedSequence& prepared : sequences) {
        const SequenceEmission& emission = *prepared.emission;
        if (!emission.tops.empty()) {
            throw std::runtime_error(
                "Local subtree root emitted top segments for " +
                emission.genome_name + "." + emission.seq_name);
        }
        // Match CactusHalConverter: an update cannot add a sequence, so an
        // absent empty sequence is deliberately ignored.
        if (!emission.bottoms.empty() ||
            genome.getSequence(emission.seq_name) != nullptr) {
            updates.emplace_back(emission.seq_name, emission.bottoms.size());
        }
    }
    genome.updateBottomDimensions(updates);
}

std::vector<BottomState> initializeBottoms(
    hal::Genome& root, const std::vector<PreparedSequence>& sequences) {
    size_t total_bottoms = 0;
    for (const PreparedSequence& prepared : sequences) {
        const size_t sequence_bottoms = prepared.emission->bottoms.size();
        if (sequence_bottoms >
            std::numeric_limits<size_t>::max() - total_bottoms) {
            throw std::overflow_error("HAL bottom row count overflow");
        }
        total_bottoms += sequence_bottoms;
    }
    if (total_bottoms == std::numeric_limits<size_t>::max()) {
        throw std::overflow_error("HAL bottom state size overflow");
    }
    std::vector<BottomState> bottom_state(total_bottoms + 1);
    hal::RamaxHdf5SegmentWriter root_rows(root);

    for (const PreparedSequence& prepared : sequences) {
        const SequenceEmission& emission = *prepared.emission;
        hal::Sequence* sequence = root.getSequence(emission.seq_name);
        if (sequence == nullptr) {
            if (prepared.length == 0 && emission.bottoms.empty()) {
                continue;
            }
            throw std::runtime_error(
                "Missing root HAL sequence while writing bottoms: " +
                emission.genome_name + "." + emission.seq_name);
        }
        if (emission.bottoms.empty()) {
            continue;
        }
        auto iterator = sequence->getBottomSegmentIterator();
        const std::string context =
            emission.genome_name + "." + emission.seq_name;
        for (const BottomSegmentLine& line : emission.bottoms) {
            if (line.name == 0 || line.name > total_bottoms ||
                bottom_state[line.name].array_index != hal::NULL_INDEX) {
                throw std::runtime_error(
                    "HAL bottom names are not a bounded unique 1-based set");
            }
            hal::BottomSegment* segment = iterator->getBottomSegment();
            const hal_index_t array_index = segment->getArrayIndex();
            root_rows.initializeBottom(
                array_index,
                absoluteStart(*sequence, line.start, context),
                line.length);
            bottom_state[line.name].array_index = array_index;
            iterator->toRight();
        }
    }

    for (size_t name = 1; name < bottom_state.size(); ++name) {
        if (bottom_state[name].array_index == hal::NULL_INDEX) {
            throw std::runtime_error(
                "HAL bottom names are not a dense producer-assigned range");
        }
    }
    return bottom_state;
}

void initializeChildTops(
    hal::Genome& root, hal::Genome& child,
    const std::vector<PreparedSequence>& sequences,
    std::vector<BottomState>& bottom_state) {
    // One tail slot per bottom replaces the former hash node per duplicate
    // (bottom, child) pair.  Clear and reuse it for each child, so this state
    // never accumulates across children.
    for (BottomState& state : bottom_state) {
        state.duplicate_tail = hal::NULL_INDEX;
    }
    const hal_index_t child_index = root.getChildIndex(&child);
    if (child_index == hal::NULL_INDEX) {
        throw std::runtime_error("HAL child is absent from its parent: " +
                                 child.getName());
    }
    hal::RamaxHdf5SegmentWriter root_rows(root);
    hal::RamaxHdf5SegmentWriter child_rows(child);

    auto parent_iterator = root.getBottomSegmentIterator();
    auto paralogy_iterator = child.getTopSegmentIterator();

    for (const PreparedSequence& prepared : sequences) {
        const SequenceEmission& emission = *prepared.emission;
        hal::Sequence* sequence = child.getSequence(emission.seq_name);
        if (sequence == nullptr) {
            if (prepared.length == 0 && emission.tops.empty()) {
                continue;
            }
            throw std::runtime_error(
                "Missing child HAL sequence while writing tops: " +
                emission.genome_name + "." + emission.seq_name);
        }
        if (emission.tops.empty()) {
            continue;
        }

        auto top_iterator = sequence->getTopSegmentIterator();
        const std::string context = emission.genome_name + "." + emission.seq_name;
        for (const TopSegmentLine& line : emission.tops) {
            hal::TopSegment* top = top_iterator->getTopSegment();
            const hal_index_t top_index = top->getArrayIndex();
            const hal_index_t top_start =
                absoluteStart(*sequence, line.start, context);
            const bool parent_reversed =
                line.parent_bottom_name ? !line.forward_to_parent : true;
            child_rows.initializeTop(
                top_index, top_start, line.length, parent_reversed);

            if (line.parent_bottom_name) {
                const uint64_t name = *line.parent_bottom_name;
                if (name == 0 || name >= bottom_state.size() ||
                    bottom_state[name].array_index == hal::NULL_INDEX) {
                    throw std::runtime_error(
                        "Top segment refers to an unknown HAL bottom name");
                }
                BottomState& state = bottom_state[name];
                const hal_index_t parent_index = state.array_index;
                parent_iterator->getBottomSegment()->setArrayIndex(&root,
                                                                   parent_index);
                hal::BottomSegment* parent =
                    parent_iterator->getBottomSegment();
                if (parent->getLength() != line.length) {
                    throw std::runtime_error(
                        "HAL parent and child segment lengths differ");
                }
                child_rows.setTopParentIndex(top_index, parent_index);

                const hal_index_t first_index =
                    parent->getChildIndex(child_index);
                if (first_index == hal::NULL_INDEX) {
                    root_rows.setBottomChild(
                        parent_index, child_index, top_index,
                        parent_reversed);
                } else {
                    paralogy_iterator->getTopSegment()->setArrayIndex(
                        &child, first_index);
                    hal::TopSegment* first =
                        paralogy_iterator->getTopSegment();
                    if (first->getStartPosition() >= top_start) {
                        throw std::runtime_error(
                            "Paralogous HAL top segments are out of scan order");
                    }
                    child_rows.setTopNextParalogyIndex(
                        top_index, first_index);

                    if (first->getNextParalogyIndex() == hal::NULL_INDEX) {
                        if (state.duplicate_tail != hal::NULL_INDEX) {
                            throw std::runtime_error(
                                "Inconsistent HAL paralogy cache state");
                        }
                        child_rows.setTopNextParalogyIndex(
                            first_index, top_index);
                        state.duplicate_tail = top_index;
                    } else {
                        if (state.duplicate_tail == hal::NULL_INDEX) {
                            throw std::runtime_error(
                                "Missing HAL paralogy cache state");
                        }
                        paralogy_iterator->getTopSegment()->setArrayIndex(
                            &child, state.duplicate_tail);
                        hal::TopSegment* previous =
                            paralogy_iterator->getTopSegment();
                        if (previous->getParentIndex() != parent_index ||
                            previous->getNextParalogyIndex() != first_index) {
                            throw std::runtime_error(
                                "Inconsistent HAL paralogy chain");
                        }
                        child_rows.setTopNextParalogyIndex(
                            state.duplicate_tail, top_index);
                        state.duplicate_tail = top_index;
                    }
                }
            }
            top_iterator->toRight();
        }
    }
}

void updateRootParseInfo(hal::Genome& root) {
    auto bottom_iterator = root.getBottomSegmentIterator();
    auto top_iterator = root.getTopSegmentIterator();
    hal::RamaxHdf5SegmentWriter root_rows(root);

    while (!bottom_iterator->atEnd() && !top_iterator->atEnd()) {
        hal::BottomSegment* bottom = bottom_iterator->getBottomSegment();
        hal::TopSegment* top = top_iterator->getTopSegment();
        const hal_index_t bottom_start = bottom->getStartPosition();
        const hal_index_t bottom_end =
            bottom_start + static_cast<hal_index_t>(bottom->getLength());
        const hal_index_t top_start = top->getStartPosition();
        const hal_index_t top_end =
            top_start + static_cast<hal_index_t>(top->getLength());

        if (bottom_start >= top_start && bottom_start < top_end) {
            root_rows.setBottomTopParseIndex(
                bottom->getArrayIndex(), top->getArrayIndex());
        }
        const bool advance_bottom =
            bottom_end <= top_end || bottom_start == bottom_end;

        if (top_start >= bottom_start && top_start < bottom_end) {
            root_rows.setTopBottomParseIndex(
                top->getArrayIndex(), bottom->getArrayIndex());
        }
        const bool advance_top =
            top_end <= bottom_end || top_start == top_end;

        if (!advance_bottom && !advance_top) {
            throw std::runtime_error(
                "HAL top and bottom segment parses do not overlap");
        }
        if (advance_bottom) {
            bottom_iterator->toRight();
        }
        if (advance_top) {
            top_iterator->toRight();
        }
    }
}

}  // namespace

NativeHalWriter::NativeHalWriter(
    const std::filesystem::path& path,
    const TreeMeta& tree,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    const SoftMask::IndexMap& softmask_indexes)
    : alignment_(openNativeAlignment(path)),
      tree_(tree),
      managers_(managers),
      softmask_indexes_(softmask_indexes) {}

NativeHalWriter::~NativeHalWriter() {
    discardAlignmentNoThrow(alignment_);
}

void NativeHalWriter::appendSubtree(
    int node_id, const std::vector<SequenceEmission>& emissions) try {
    if (!alignment_) {
        throw std::runtime_error("Cannot append to a closed HAL writer");
    }
    if (node_id < 0 || static_cast<size_t>(node_id) >= tree_.nodes.size()) {
        throw std::out_of_range("HAL subtree node ID is out of range");
    }

    const TreeNodeMeta& node = tree_.nodes[static_cast<size_t>(node_id)];
    GenomeOrder genome_order;
    genome_order.reserve(node.children.size() + 1);
    genome_order.push_back(&node.name);
    for (int child_id : node.children) {
        if (child_id < 0 ||
            static_cast<size_t>(child_id) >= tree_.nodes.size()) {
            throw std::out_of_range("HAL child node ID is out of range");
        }
        genome_order.push_back(
            &tree_.nodes[static_cast<size_t>(child_id)].name);
    }
    PreparedByGenome prepared = prepareEmissions(genome_order, emissions);

    GenomeCloser closer(alignment_.get());
    hal::Genome* root = nullptr;
    const bool new_alignment = alignment_->getNumGenomes() == 0;
    if (new_alignment) {
        if (node_id != tree_.root_id) {
            throw std::runtime_error(
                "The first native HAL subtree must be the phylogeny root");
        }
        root = alignment_->addRootGenome(node.name);
    } else {
        root = alignment_->openGenome(node.name);
        if (root == nullptr) {
            throw std::runtime_error(
                "Cannot locate existing local-subtree root genome: " +
                node.name);
        }
    }
    closer.add(root);

    std::vector<hal::Genome*> children;
    children.reserve(node.children.size());
    for (int child_id : node.children) {
        const TreeNodeMeta& child_node =
            tree_.nodes[static_cast<size_t>(child_id)];
        hal::Genome* existing = alignment_->openGenome(child_node.name);
        if (existing != nullptr) {
            closer.add(existing);
            throw std::runtime_error(
                "Local-subtree child genome already exists in HAL: " +
                child_node.name);
        }
        hal::Genome* child = alignment_->addLeafGenome(
            child_node.name, node.name,
            roundedBranchLength(child_node.branch_length_to_parent));
        closer.add(child);
        children.push_back(child);
    }

    const auto& root_sequences = prepared.front();
    if (!root_sequences.empty()) {
        if (new_alignment) {
            setNewGenomeDimensionsAndDna(
                *root, root_sequences, true, managers_, softmask_indexes_);
        } else {
            updateExistingRootDimensions(*root, root_sequences);
        }
    }

    std::vector<BottomState> bottom_state =
        initializeBottoms(*root, root_sequences);
    // For B root bottoms the complete writer-owned row-linking allocation is
    // one BottomState per dense name plus the unused zero slot:
    // (B + 1) * (2 * sizeof(hal_index_t)); there is no T-row-sized writer
    // scratch and this allocation is reused for every child.  HAL's native
    // arrays retain the root's B bottom rows and any existing T_root top rows,
    // plus only the current child's T_child top rows.  Closing each child
    // flushes its arrays before the next child's arrays are populated.
    for (size_t child_offset = 0; child_offset < children.size();
         ++child_offset) {
        const auto& child_sequences = prepared[child_offset + 1];
        if (!child_sequences.empty()) {
            setNewGenomeDimensionsAndDna(
                *children[child_offset], child_sequences, false,
                managers_, softmask_indexes_);
        }
        initializeChildTops(*root, *children[child_offset], child_sequences,
                            bottom_state);
        closer.closeGenome(children[child_offset]);
        children[child_offset] = nullptr;
    }
    updateRootParseInfo(*root);

    // closeGenome writes all buffered arrays.  Keeping this per-subtree makes
    // the native one-process lifecycle no less bounded than the former helper
    // process lifecycle, and propagates write failures to the caller.
    closer.close();
} catch (const H5::Exception& error) {
    discardAlignmentNoThrow(alignment_);
    throw hdf5Failure("Unable to append native HAL subtree", error);
}

void NativeHalWriter::close() {
    if (!alignment_) {
        return;
    }

    hal::AlignmentPtr closing = std::move(alignment_);
    try {
        closing->close();
        closing.reset();
    } catch (const H5::Exception& error) {
        closing.reset();
        throw hdf5Failure("Unable to close native HAL writer", error);
    }
}

}  // namespace RaMesh::hal_export
