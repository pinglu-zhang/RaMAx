# Third-party notices

RaMAx uses or distributes the following third-party components. Their original
copyright notices and license texts remain authoritative in their source
directories or upstream distributions.

| Component | Use | License |
|---|---|---|
| HAL | Hierarchical alignment storage and export | MIT |
| sonLib | HAL support library | MIT |
| libdivsufsort | Suffix-array construction | MIT |
| ParlayLib | Parallel primitives | MIT |
| spdlog | Logging | MIT |
| cereal | Serialization | BSD-3-Clause |
| SDSL | Succinct data structures | BSD-3-Clause |
| CLI11 | Command-line parsing | BSD-3-Clause |
| KSW2 | Banded pairwise alignment | MIT |
| bigBWT and bundled suffix-array implementations | FM-index construction | See bundled source notices |
| minipoa v1.4.2 | External partial-order MSA | MIT |

Bundled full license texts are available under `third_party/` where supplied.
A normal source build treats minipoa and halAppendCactusSubtree as external
runtime dependencies. The official Conda and Docker packaging configurations
install minipoa 1.4.2 from the malab Conda channel and copy the validated,
precompiled cactus-bin-v2.9.9 halAppendCactusSubtree executable together with
the applicable HAL license.
