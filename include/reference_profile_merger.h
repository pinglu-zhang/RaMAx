#ifndef RAMAX_REFERENCE_PROFILE_MERGER_H
#define RAMAX_REFERENCE_PROFILE_MERGER_H

#include "align.h"

namespace RaMesh::Alignment {

uint_t mergeReferenceProfile(
    const ChrName& reference_name,
    std::unordered_map<ChrName, std::string>& sequences,
    const std::unordered_map<ChrName, Cigar_t>& cigars);

}  // namespace RaMesh::Alignment

#endif
