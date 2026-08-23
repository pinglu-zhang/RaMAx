# Changelog

## 1.0.7 - 2026-08-17

- Redefined `--restart` as raw/clean FASTA and FM-index cache reuse; alignment,
  clustering, graph construction, and export always restart from the beginning.
- Added explicit restart parameter overrides, complete repeated-output
  replacement, and atomic persistence of the latest effective schema-2 config.
- Added size/mtime completion markers and atomic publication for preprocessing,
  softmask, and all three FM-index components.
- Added schema-1 workdir migration, input-identity validation, stale post-index
  cleanup, previous-log archival, and cache reuse/rebuild summaries.

## 1.0.5 - 2026-08-14

- Added compatible neighboring-Block merging and bounded missing-sequence
  realignment, including zero-reference-gap and three-Block windows.
- Added conservative structural-discontinuity repair and short-Block merging.
- Added export-time repair for homologous insertions anchored at nearby
  reference positions; direct MAF and HAL export share the same implementation.
- Reduced fixed-point, graph-transaction, and external-MSA overhead while
  preserving candidate order, CIGARs, and graph rollback checks.
- Fixed concurrent external-MSA cache coordination.
- Parallelized cross-species anchor search, sparse clustering, cluster
  extension, and both coordinate-DP stages while preserving deterministic
  result collection and serial graph mutation.
- Prevented overlapping reference Blocks from merging when repeated
  non-reference participant species would be discarded by the single-occurrence
  Block model.
- Rebuilt HAL export around occurrence-level common refinement and
  evidence-constrained ancestral paths, eliminating leaf-contig-driven
  ancestor-sequence fragmentation while preserving direct-MAF leaf coverage.
- Restricted default graph construction to anchors selected on both reference
  and query coordinates; repeat-secondary occurrences remain explicit
  `--allow-mem` behavior.
- Preserved terminal query-only CIGAR operations when overlapping Blocks are
  split and merged on both forward and reverse strands.
- Isolated minipoa scratch files under the RaMAx work directory and removed
  them after each invocation.
- Simplified the Block-optimization command-line interface and introduced
  `--optimize-blocks`.
- Consolidated restart settings into a versioned `config.json` schema.
