# Realignment architecture

This is an internal design reference. For user-facing commands and defaults,
see [Command-line parameters](parameters.md) and [Usage guide](usage.md).

## Execution order

Block optimization runs after each reference-guided alignment round in a
deterministic order:

1. merge compatible contiguous Blocks;
2. realign bounded missing-sequence windows;
3. merge newly compatible contiguous Blocks;
4. repair high-confidence structural discontinuities;
5. merge newly compatible contiguous Blocks;
6. after all rounds, merge short Blocks;
7. export MAF, HAL, or PAF, including nearby cross-anchor insertion repair.

Disabling an optimization module avoids its scan and cache allocation.

## Missing-window planning

`BlockViewBuilder` parses Blocks and caches their participant keys, paths,
strands, coordinates, and CIGAR properties for one graph generation.
`MissingWindowPlanner` scans, ranks, and resolves conflicts between ordinary,
zero-reference-gap, K-0-K, K-K-K, and hybrid windows. Candidate ordering and
fixed-point barriers are stable. Failed preparation releases reservations so
the next compatible candidate can be reconsidered.

Reference-empty windows use the same pairwise graph representation without an
external MSA call. Other accepted windows invoke minipoa once with:

```text
minipoa input.fa -S -f 0 -r1
```

## Ownership and graph transactions

Prepared Blocks and Segments remain provisional until path splices, Block-pool
replacement, sampling rebuild, and graph audits all succeed. A failed check
restores the original Block order, parent links, paths, and sampling. Phase-wide
short-Block processing keeps the same invariants while batching the final pool
and sampling rebuild.

Every committed replacement must preserve genomic sequence coordinates and
produce CIGARs whose reference and query consumption matches the corresponding
Segment intervals.

## External MSA execution

`ExternalMsaRunner` launches the configured minipoa executable directly,
validates its output, and returns aligned rows in memory. Input normally uses
`memfd`; the file fallback and minipoa's own temporary files are confined to:

```text
<work>/minipoa_tmp/
```

Per-process files are removed after success, failure, or parse rejection.
Successful validated results may be cached in memory, with byte-for-byte key
confirmation and single-flight coordination for concurrent duplicate requests.

At runtime RaMAx searches for minipoa in this order:

1. the path selected during CMake configuration;
2. the directory containing the running `ramax` executable;
3. `PATH`.

If an enabled module requires minipoa and no executable is available, RaMAx
fails before graph construction.

## MAF and HAL consistency

Both exporters call the shared `mergeAlignmentByRef()` implementation. It first
constructs the reference profile from pairwise CIGARs, then repairs only
strictly improving nearby cross-anchor insertion groups. Any validation or
alignment failure retains the original rows. HAL-specific graph structures are
not modified by this repair.

## Required invariants

- aligned rows have equal length and ungap to their original sequences;
- graph paths remain coordinate ordered and bidirectionally linked;
- each Segment has exactly one valid parent Block;
- failed candidates and failed transactions leave the graph unchanged;
- candidate and commit order is deterministic across thread counts;
- MAF and HAL use the same reference-projection and insertion-repair code.
