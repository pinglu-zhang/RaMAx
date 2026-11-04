# Third-party notices

RaMAx uses or distributes the following third-party components. Their original
copyright notices and license texts remain authoritative in their source
directories or upstream distributions.

| Component | Use | License |
|---|---|---|
| HAL | Hierarchical alignment storage and export | MIT |
| Cactus2HAL | Native subtree-conversion semantics adapted from CactusHalConverter | MIT |
| sonLib | HAL support library | MIT |
| libdivsufsort | Suffix-array construction | MIT |
| sufkit-derived code | SIMD comparison and suffix-link techniques adapted into RaMAx; no sufkit library is linked | MIT |
| ParlayLib | Parallel primitives | MIT |
| spdlog | Logging | MIT |
| cereal | Serialization | BSD-3-Clause |
| SDSL | Succinct data structures | BSD-3-Clause |
| CLI11 | Command-line parsing | BSD-3-Clause |
| KSW2 | Banded pairwise alignment | MIT |
| bigBWT and bundled suffix-array implementations | FM-index construction | See bundled source notices |
| minipoa v1.4.2 | External partial-order MSA | MIT |

Bundled full license texts are available under `third_party/` where supplied.
A normal source build treats minipoa as an external runtime dependency. The
official Conda and Docker configurations install minipoa 1.4.2 from the malab
Conda channel. HAL is linked into RaMAx; no Cactus append executable is required.
Build-private safety corrections retain the original HAL source notices.

## Cactus2HAL conversion logic

The native subtree writer follows the segment, paralogy, and parse-index
conversion in
[`CactusHalConverter`](https://github.com/ComparativeGenomicsToolkit/cactus2hal/blob/197ffc43c34278d8474f1a7a4ed6c2b125cc1e86/src/cactusHalConverter.cpp).
Its source identifies the MIT license and the following copyright:

Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)

The parent [Cactus MIT license](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/LICENSE.txt)
also carries:

Copyright (C) 2011 by Benedict Paten (benedictpaten@gmail.com)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
