#include "ramesh.h"
#include "hal/export.h"

#include <stdexcept>

// ============================================================
// RaMeshMultiGenomeGraph 导出接口
// ============================================================
namespace RaMesh {

    void RaMeshMultiGenomeGraph::exportToMaf(
        const FilePath& maf_path,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seq_mgrs,
        bool pairwise_mode) const {
        hal_export::exportToMaf(
            blocks,
            maf_path,
            seq_mgrs,
            pairwise_mode);
    }

    void RaMeshMultiGenomeGraph::exportToHal(
        const FilePath& hal_path,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
        const NewickParser& parser,
        const std::string& root_name,
        int parallel_threads,
        const SoftMask::PathMap& softmask_paths) const {
        if (softmask_paths.size() != seqpro_managers.size()) {
            throw std::runtime_error(
                "HAL export requires one soft-mask index per leaf genome");
        }
        for (const auto& [species, unused_manager] : seqpro_managers) {
            (void)unused_manager;
            if (!softmask_paths.contains(species)) {
                throw std::runtime_error(
                    "Missing HAL soft-mask index for species: " + species);
            }
        }
        const SoftMask::IndexMap softmask_indexes =
            SoftMask::loadIndexes(softmask_paths);
        hal_export::exportToHal(
            blocks,
            hal_path,
            seqpro_managers,
            parser,
            root_name,
            softmask_indexes,
            hal_export::ExportConfig{
                .parallel_threads = parallel_threads,
            });
    }

} // namespace RaMesh
