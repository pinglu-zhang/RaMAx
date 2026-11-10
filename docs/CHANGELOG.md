# Changelog

## 1.1.1 - 2026-11-11

- Balance wfmash chunks within thread and pair limits.
- Fix legacy glibc tool startup and minipoa scratch handling.
- Validate the Conda package and HAL export on glibc 2.17.

## 1.1.0 - 2026-11-04

- Export HAL directly through the bundled HAL library, without requiring a
  Cactus append executable. Preserve leaf contig names and softmask.
- Share bounded export preparation between MAF and HAL and release the source
  graph when it is no longer needed by remaining outputs.
- Honor an explicit first-round reference and use the least fragmented
  nonempty genome as a fallback when all inputs exceed the automatic limit.
- Recover non-overlapping PAF fragments during wfmash graph normalization.
- Prepare Linux x86-64 Conda and Docker Hub packaging from the same versioned
  source, with installed-package and four-format integration checks.
- Remove the internal `qa/` source directory and its build targets; keep
  release-only validation helpers under `ci/`.
- This release uses the committed dev implementation; experimental dev-pan
  modes and local collinear refinement are not included.

## Unreleased

- Write HAL in-process, removing the external `halAppendCactusSubtree`
  dependency.
- Reduce MAF/HAL export memory use through shared preparation and incremental
  processing; accelerate HAL compression and writing.
- Preserve existing outputs on export failure and reject directory destinations
  without moving their contents.
- Fix HAL cleanup after write failures and sequence-name compatibility with
  modern HDF5 readers.
- Allow source builds without a bundled `bin/` directory, retaining external
  tool overrides and runtime lookup.

## 1.0.9 - 2026-08-30

- Improved end-to-end performance across anchor processing, cluster extension,
  chromosome-level dynamic programming, graph extension, and deterministic
  multi-graph merging.
- Reduced long-tail latency by scheduling expensive alignment work by estimated
  cost and improving parallel cluster-extension coordination.
- Added workload-aware wfmash scheduling that shares the RaMAx thread budget
  without oversubscribing the host.
- Added a fixed per-query wfmash time budget. Timed-out process groups are
  terminated safely, recorded for diagnosis, and rerouted to the native RaMAx
  alignment path.
- Reduced peak and retained memory by releasing consumed anchors, clusters,
  CIGAR buffers, Blocks, and Segments earlier.
- Accelerated reference-node extension scans and graph merging while preserving
  alignment semantics, deterministic ordering, and existing output formats.

## 1.0.8 - 2026-08-26

- Replaced native anchor search with a complete in-memory SA/ISA/LCP backend.
  `--sa-sampling-rate` currently accepts only `1`; the index is rebuilt for
  every process and is not persisted in the work directory.
- Build suffix arrays with divsufsort below 1024 MiB and CaPS at or above that
  threshold, while preserving the existing RaMAx anchor-search contract.
- Upgraded restart configuration to schema 6 and retained explicit migration
  defaults for schema-1 through schema-5 work directories.
- Fixed missing-window realignment coverage regressions without widening
  ordinary or hybrid merge spans.
- Reassembled MAF runs by source block so normalization preserves the intended
  source-block boundaries.
- Added a shared external-tool directory for source, Conda, and Docker builds,
  while retaining explicit per-tool overrides and runtime sibling lookup.
- Fixed Conda builds against a glibc 2.17 sysroot by compiling the
  `posix_spawn_file_actions_addchdir_np` path only where the interface is
  available.
- Updated the Linux x86-64 Conda and Docker packaging configuration for RaMAx
  1.0.8 and fixed exact companion-tool builds. These packaging checks establish
  build and runtime dependency compatibility, not cross-platform biological
  equivalence or performance.

## 1.0.7 - 2026-08-17

- Packaged environments use Samtools/HTSlib 1.23.1, the newest release that
  coexists with the PGGB-compatible wfmash 0.14 build. RaMAx uses only the
  stable `samtools faidx --fai-idx` interface.
- Added an early aggregated dependency preflight for `minipoa`, `wfmash`, and
  `mash`. `halAppendCactusSubtree` is required only for HAL output; non-HAL
  runs warn and continue when it is unavailable.
- Replaced the HAL append shell command with a parameterized subprocess call.
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
