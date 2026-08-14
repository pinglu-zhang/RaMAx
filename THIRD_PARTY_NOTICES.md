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
minipoa is a separate runtime dependency and is not downloaded or installed by
the RaMAx build.
