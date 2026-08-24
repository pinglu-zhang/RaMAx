# sufkit-derived search kernels

- Source working tree: `/mnt/d/code/sufkit`
- Base commit: `fd1abbdb4a486abcbca3be7d915f43d3638f8b16`
- Local development branch inspected: `codex/v0.2.0-low-level-performance`
  (working tree contained uncommitted low-level performance changes)
- Inspected snapshot SHA-256:
  - `src/sequence_compare.hpp`: `8240ca49f277896ba993f9ea1fb2905c5f6073a75f9c0d77f4aaef80daf1d5a1`
  - `src/suffix_array.cpp`: `634480b98a921548ce908460ca8c2326b2358c5ef5f8fe9185c52fe150a82151`
- License: MIT

RaMAx does not vendor or link the complete sufkit library. The compact
`Suffix_Array_Index` implementation adapts the SIMD byte-comparison kernel,
Kasai SA/ISA/LCP construction logic, and generalized suffix-link interval
derivation from sufkit while reusing RaMAx's existing libdivsufsort targets.
This keeps the RaMAx anchor, coordinate, filtering, and downstream graph
interfaces unchanged.
