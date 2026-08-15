# Changelog

## 1.0.5 - 2026-08-14

- Added compatible neighboring-Block merging and bounded missing-sequence
  realignment, including zero-reference-gap and three-Block windows.
- Added conservative structural-discontinuity repair and short-Block merging.
- Added export-time repair for homologous insertions anchored at nearby
  reference positions; direct MAF and HAL export share the same implementation.
- Reduced fixed-point, graph-transaction, and external-MSA overhead while
  preserving candidate order, CIGARs, and graph rollback checks.
- Fixed concurrent external-MSA cache coordination.
- Isolated minipoa scratch files under the RaMAx work directory and removed
  them after each invocation.
- Simplified the Block-optimization command-line interface and introduced
  `--optimize-blocks`.
- Consolidated restart settings into a versioned `config.json` schema.
