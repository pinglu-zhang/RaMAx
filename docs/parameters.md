# Command-line parameters

This page documents the public RaMAx 1.0.7 interface. Values and ranges follow
the current CLI implementation. Run `ramax --help` to inspect the installed
binary.

## General

| Option | Type/default | Description |
|---|---|---|
| `-h`, `--help` | flag | Print help and exit. |
| `-v`, `--version` | flag | Print the RaMAx version and exit. |

## Input Files

| Option | Type/default | Description |
|---|---|---|
| `-i`, `--input` | path; required for a new run | Seqfile containing genome mappings. A Newick tree is required only for HAL output. |

## Output

| Option | Type/default | Description |
|---|---|---|
| `-o`, `--output` | repeatable path; at least one required | Output file. Repeat `-o` to export several formats from one alignment. Each suffix must be `.maf`, `.hal`, `.paf`, or `.gfa`, and each format may appear once. Export order is MAF, PAF, GFA, HAL. |
| `-w`, `--workdir` | path; required | Intermediate and routing-artifact directory. It must be empty for a new Release run and is preserved after success. |
| `--paf-mode` | `connected`; `connected` or `all` | PAF pair-selection policy. Valid when any `-o` is `.paf`; `connected` adds the minimum deterministic supplemental pairs needed for column-wise same-base connectivity, while `all` is the all-pairs baseline. |
| `--gfa-version` | `1.1`; exactly `1.0` or `1.1` | Native GFA path encoding. `1.0` writes P-lines; `1.1` writes structured W-lines. Valid only when any `-o` is `.gfa`. |
| `--gfa-profile` | `exact`; exactly `exact` or `compact` | Native GFA graph construction. `exact` preserves every maximal exact-run relation. `compact` enables compact-v2-balanced while preserving every cleaned path base and writes an exact audit shadow plus staged transform reports under `work/gfa/`. Valid only with `.gfa` output. |
| `--root` | string; automatic | Preferred HAL root name. Valid when any `-o` is `.hal`; HAL requires a Newick tree. If an artificial unnamed root is required and no name is supplied, RaMAx uses `ancestor`. |

## Software Parameters

| Option | Type/default/range | Description |
|---|---|---|
| `--chunk_size` | integer; `10000000`; `1000000..INT_MAX` bp | Reference/query chunk size used during parallel anchor search. |
| `--ref` | string; automatic | Prefer a named genome as reference. When omitted, RaMAx determines reference order from assembly statistics. |
| `--overlap_size` | integer; `0`; `0..INT_MAX` bp | Overlap between adjacent chunks. It must be smaller than `--chunk_size`. |
| `--min_anchor_length` | integer; `20`; `1..INT_MAX` bp | Minimum anchor length. |
| `--max_anchor_frequency` | integer; `50`; `0..INT_MAX` | Maximum accepted anchor occurrence frequency. |
| `--search-mode` | `accurate`; `fast`, `middle`, or `accurate` | Anchor-search strategy. Numeric values `0`, `1`, and `2` map to the same modes. |
| `--allow-mem` | flag; off | Allow MEM anchors instead of restricting discovery to MUM anchors. |
| `--one-round` | flag; off | Stop after one reference-guided alignment round. |
| `--slow-build` | flag; off | Use the slower index-building implementation. |
| `--sampling-interval` | integer; `32`; `1..INT_MAX` bp | Sampling interval for the reference index. |
| `--min-span` | integer; `65`; `1..INT_MAX` bp | Minimum span used during graph construction. |
| `--near-distance` | float; `0.01`; `0..1` | In the first round, route a query to wfmash only when its Mash distance is strictly smaller than this value. Equality stays on the RaMAx backend. |
| `--far-distance` | float; `0.02`; `0..1` | Route first-round queries with strictly larger Mash distance to mm2-plus. Equality stays on legacy RaMAx; use `1` to disable mm2-plus routing. |

The thresholds must satisfy
`0 <= --near-distance < --far-distance <= 1`. First-round near pairwise mapping
uses PGGB-compatible wfmash `v0.14.0-0-g517e1bc` with
`-s 5000 -l 25000 -p 95 -n 1 -k 19 -H 0.001 -Y '#'
--hg-filter-ani-diff 30 --approx-map`; precise alignment reuses the mapping
PAF with `--invert-filtering`. Later reference rounds always use the RaMAx
FM-index backend. Unlike PGGB, RaMAx launches independent reference/query
pairs and therefore does not use `--lower-triangular`.

First-round distant pairwise mapping uses mm2plus 1.3 (Minimap2 2.31/r1302)
with `-x asm20 -c --eqx --secondary=no`. RaMAx builds a single-part reference
index once, validates the required `cg`, `tp`, and `AS` PAF tags, removes exact
duplicates and two-axis interval conflicts, and then imports the normalized
primary alignments as graph Anchors. A failed or empty query falls back to the
legacy first-round path. All later rounds use legacy RaMAx.

Parameter names containing underscores are retained in the current public CLI
for the anchor/chunk settings shown above.

## Graph Optimization

All four graph-optimization modules are enabled by default in 1.0.5.

| Option | Type/default/range | Description |
|---|---|---|
| `--optimize-blocks` | flag; default behavior | Explicitly enable the complete default optimization set. |
| `--merge-blocks` | flag; on | Merge compatible neighboring Blocks. |
| `--merge-gap` | integer; `100`; `0..10000` bp | Maximum query-coordinate gap accepted when merging existing neighboring Blocks. `0` requires coordinate continuity. |
| `--realign-missing` | flag; on | Realign bounded windows in which species participation differs between Blocks. |
| `--realign-span` | integer; `3000`; `1..10000` bp | Maximum complete local span considered by missing-sequence and external-MSA realignment. |
| `--zero-gap-span` | integer; `200`; `1..3000` bp | Maximum span for zero-reference-gap missing windows. It does not replace `--merge-gap`. |
| `--repair-breaks` | flag; on | Repair high-confidence target, strand, or order discontinuities after missing-window cleanup. |
| `--break-span` | integer; `1000`; `1..1000` bp | Maximum structural-discontinuity repair window. |
| `--merge-short-blocks` | flag; on | Try to merge reference-order short Blocks using banded KSW2. The internal merge threshold is 500 bp; failed merges are not deleted because the deletion threshold is 0 bp. |

The span parameters serve different candidate classes:

- `--merge-gap` checks coordinate gaps between already participating query
  Segments during ordinary Block merging.
- `--zero-gap-span` limits missing-species windows whose reference interval is
  empty or nearly empty.
- `--realign-span` is the broader bound for missing-species local realignment.
- `--break-span` applies only to structural-discontinuity repair.

The CLI validates that `--merge-gap` belongs to Block merging,
`--realign-span` and `--zero-gap-span` belong to missing-window realignment,
and `--break-span` belongs to structural-break repair. Since the corresponding
modules are enabled by default, the normal numeric-only overrides are valid.

## Performance

| Option | Type/default/range | Description |
|---|---|---|
| `-t`, `--threads` | integer; hardware concurrency; `1..INT_MAX` | Worker-thread count. It may also be supplied by `RAMAx_THREADS`. |
| `--restart` | flag; off | Reuse validated raw/clean FASTA and FM-index caches, then rerun alignment from the beginning. |

In restart mode, `-w` identifies the interrupted run and `-i` is forbidden.
Explicit output, algorithm, optimization, thread, root, PAF-mode, GFA-version, GFA-profile, and logging
options override the latest saved configuration. If any `-o` is supplied, the
complete repeated `-o` list replaces the saved output list. Existing boolean
flags retain their current one-way behavior; no new negative forms are added.

## Output Control

| Option | Type/default | Description |
|---|---|---|
| `--log-level` | `info`; `debug`, `info`, `warn`, or `error` | Select logging detail. |
| `--verbose` | flag; off | Enable additional diagnostics. Mutually exclusive with `--quiet`. |
| `--quiet` | flag; off | Emit errors only. Mutually exclusive with `--verbose`. |

## PAF companion FASTA

`ramax-paf-fasta` is a separate installed command and does not change the
`ramax` alignment interface.

| Option | Type/default | Description |
|---|---|---|
| `-i`, `--input` | seqfile; required | Input RaMAx seqfile, optionally beginning with a Newick record. |
| `-o`, `--output` | path; required | Combined `.fa`, `.fasta`, or `.fna`, optionally followed by `.gz`. |
| `--force` | flag; off | Atomically replace an existing output after complete generation. |
