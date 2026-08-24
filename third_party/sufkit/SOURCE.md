# sufkit-derived search kernels

- Source working tree inspected: `/mnt/d/code/sufkit`
- Full-SA integration base: `fd1abbdb4a486abcbca3be7d915f43d3638f8b16`
- Current local sufkit main inspected: `e9a430ce6f3deba269961a3b901c5292036190c6`
- License: MIT

RaMAx does not vendor or link the complete sufkit library. The native
`Suffix_Array_Index` keeps the SIMD byte-comparison kernel and generalized
suffix-link interval derivation adapted from sufkit so Anchor filtering,
occurrence order, coordinate conversion, search-mode advancement, cluster, DP,
graph, and export interfaces remain unchanged.

The suffix-array backend stores complete K=1 SA/ISA/LCP arrays. References
whose FASTA file is smaller than 1024 MiB use the bundled libdivsufsort to build
SA followed by the complete SIMD/Kasai LCP pass; references at or above the
threshold use the bundled CaPS-SA implementation to construct SA and LCP
together. RaMAx builds complete ISA arrays with OpenMP in both cases. No
sampled-SA construction or residue-recovery path is active in this version.
