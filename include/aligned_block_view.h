#ifndef RAMAX_ALIGNED_BLOCK_VIEW_H
#define RAMAX_ALIGNED_BLOCK_VIEW_H

#include "ramesh.h"

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace RaMesh::Export {

struct AlignedBlockRow {
    SpeciesName species;
    ChrName chromosome;
    std::string name;
    SegPtr segment;
    std::uint64_t sequence_length{0};
    std::string aligned;
};

struct AlignedBlockView {
    std::vector<AlignedBlockRow> rows;
    std::size_t reference_index{0};
};

std::string qualifiedName(const SpeciesName& species,
                          const ChrName& chromosome);

std::string fetchSequence(const SeqPro::ManagerVariant& manager,
                          const ChrName& chromosome,
                          Coord_t start,
                          Coord_t length);

std::uint64_t fetchLength(const SeqPro::ManagerVariant& manager,
                          const ChrName& chromosome);

std::vector<std::string> fetchNames(const SeqPro::ManagerVariant& manager);

AlignedBlockView prepareAlignedBlock(
    const BlockPtr& block,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
    bool only_primary);

}  // namespace RaMesh::Export

#endif
