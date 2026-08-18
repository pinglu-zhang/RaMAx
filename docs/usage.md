# Usage guide

## Input seqfile

A seqfile contains one genome name and local FASTA path or supported URL per
non-empty line. MAF and other non-HAL outputs accept a mapping-only seqfile:

```text
human      /data/genomes/human.fa
chimp      /data/genomes/chimp.fa
orangutan  /data/genomes/orangutan.fa.gz
```

HAL output requires a Newick species tree as the first non-empty record:

```text
((human:0.005,chimp:0.005):0.02,orangutan:0.025);
human      /data/genomes/human.fa
chimp      /data/genomes/chimp.fa
orangutan  /data/genomes/orangutan.fa.gz
```

A tree may still be supplied for non-HAL output, but it does not control
reference ordering. RaMAx selects references from assembly N50, total length,
genome name, and an optional `--ref`. After selecting the first reference,
RaMAx requires Mash 2.3 and estimates whole-genome distances to every other
input with `k=31` and sketch size 20,000. The normalized table is written to
`<workdir>/similarity/mash_first_reference.tsv`. First-round queries with
`d < --near-distance` (default `0.01`) use wfmash; only a validated final PAF
removes a query from the first-round FM-index path. Failed wfmash queries fall
back to RaMAx, and all later rounds use the original RaMAx algorithm. The
reserved `--far-distance` default is `0.02` and does not yet affect routing.
RaMAx requires Mash 2.3, PGGB-compatible wfmash
`v0.14.0-0-g517e1bc`, and Samtools/HTSlib 1.24. All wfmash
FAI files are produced by `samtools faidx`, including indexes for multi-FASTA
inputs and collision-safe query views. wfmash 0.14 has no persistent reference
index interface, so independent first-round pairs rebuild the reference index.
`--root` is valid whenever the output list contains HAL.
When a tree is present, its leaf names should match the genome mappings. Use
absolute FASTA paths when jobs may be launched from different directories.

## Basic commands

Write MAF:

```bash
ramax \
  -i /data/project/seqfile.txt \
  -o /data/project/results/alignment.maf \
  -w /data/project/work/ramax-maf \
  -t 16
```

Write HAL:

```bash
ramax \
  -i /data/project/seqfile.txt \
  -o /data/project/results/alignment.hal \
  -w /data/project/work/ramax-hal \
  -t 16 \
  --root ancestor
```

Write sparse, information-complete PAF:

```bash
ramax \
  -i /data/project/seqfile.txt \
  -o /data/project/results/alignment.paf \
  -w /data/project/work/ramax-paf \
  -t 16
```

Write all supported formats after one alignment and graph construction:

```bash
ramax \
  -i /data/project/seqfile.txt \
  -w /data/project/work/ramax-all \
  -o /data/project/results/alignment.maf \
  -o /data/project/results/alignment.paf \
  -o /data/project/results/alignment.hal \
  -t 16
```

`-o/--output` is repeatable, but each format may appear only once. RaMAx
validates every suffix before alignment, constructs the alignment graph once,
and exports in the fixed order MAF, PAF, HAL. A mixed output list containing
HAL requires a Newick tree even when MAF or PAF is the first `-o`.

PAF defaults to `--paf-mode connected`. Use `--paf-mode all` to emit every
primary path pair that shares at least one non-gap alignment column. A PAF-only
run does not require a Newick tree. `--paf-mode` is accepted whenever the
output list contains PAF and rejected otherwise.

Build the seqwish sequence input from the same seqfile. The companion tool
preserves species and contig order, writes the exact `species.contig` names
used by PAF, and applies the same base normalization:

```bash
ramax-paf-fasta \
  -i /data/project/seqfile.txt \
  -o /data/project/results/sequences.fa.gz

seqwish \
  -s /data/project/results/sequences.fa.gz \
  -p /data/project/results/alignment.paf \
  -g /data/project/results/graph.gfa \
  -k 0
```

Each output suffix selects its format. `.maf`, `.hal`, and `.paf` are accepted.
Duplicate paths, duplicate formats, and unsupported suffixes are rejected.

## Default Block optimization

RaMAx 1.0.5 enables the complete Block-optimization pipeline by default. The
following command is therefore equivalent to the basic MAF example above:

```bash
ramax \
  -i /data/project/seqfile.txt \
  -o /data/project/results/alignment.maf \
  -w /data/project/work/ramax-maf \
  -t 16 \
  --optimize-blocks
```

The default set is:

```text
--merge-blocks --merge-gap 100
--realign-missing --realign-span 3000 --zero-gap-span 200
--repair-breaks --break-span 1000
--merge-short-blocks
```

`--optimize-blocks` explicitly requests this same default set. The individual
module flags also remain public, but the four modules already start enabled in
1.0.5.

Short-Block processing considers Blocks whose reference Segment is at most
500 bp. Missing species are tested with banded KSW2 before a merge is accepted.
The deletion threshold is 0 bp, so an unmerged Block is retained.

## Selecting threads

`-t/--threads` controls RaMAx worker threads. If omitted, the default is the
number of hardware threads reported by the system. The `RAMAx_THREADS`
environment variable can also provide the value:

```bash
RAMAx_THREADS=16 ramax \
  -i /data/project/seqfile.txt \
  -o /data/project/results/alignment.maf \
  -w /data/project/work/ramax-maf
```

An explicit `--threads` value takes precedence over the environment default.
Use fewer threads when memory is the limiting resource.

## Work-directory lifecycle

The work directory stores copied/cleaned sequences, indexes, serialized
configuration, logs, mask intervals, and minipoa scratch files while a job is
running. Anchor, cluster, Block, DP, and graph objects are not serialized.

For a new Release run, the directory must not contain existing files. RaMAx
preserves the work directory after a successful export so that
`similarity/mash_first_reference.tsv` and `wfmash/round_0/` remain available.
Failed or interrupted runs also leave it in place for diagnosis and restart.

On restart, RaMAx retains only reusable preprocessing and FM-index caches. It
archives the previous log, removes stale `result/`, `mask_interval/`, and
`minipoa_tmp/` directories, and starts alignment from the first reference
round.

Always place the final MAF, HAL, or PAF outside the work directory:

```text
/data/project/results/alignment.maf   # final output
/data/project/work/ramax-run          # temporary work directory
```

See [Restart and work directories](restart.md) for recovery behavior.

## Logging

The default log level is `info`. Use `--log-level debug` when investigating a
candidate or graph transaction, `--quiet` for errors only, and `--verbose` for
additional diagnostics. `--verbose` and `--quiet` are mutually exclusive.
