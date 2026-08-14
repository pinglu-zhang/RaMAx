# Usage guide

## Input seqfile

RaMAx accepts a Cactus-compatible seqfile. The first non-comment line is a
Newick species tree. Each following non-empty line maps one leaf name to a
local FASTA file or supported URL:

```text
((human:0.005,chimp:0.005):0.02,orangutan:0.025);
human      /data/genomes/human.fa
chimp      /data/genomes/chimp.fa
orangutan  /data/genomes/orangutan.fa.gz
```

Every tree leaf must have exactly one genome mapping, and mapping names must
match the tree leaves. Use absolute FASTA paths when jobs may be launched from
different working directories.

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

The output suffix selects the format. Only `.maf` and `.hal` are accepted.

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
configuration, logs, graph state, mask intervals, and minipoa scratch files
while a job is running.

For a new Release run, the directory must not contain existing files. RaMAx
removes the entire work directory after a successful export. Failed or
interrupted runs leave it in place for diagnosis and restart.

Always place the final MAF or HAL outside the work directory:

```text
/data/project/results/alignment.maf   # final output
/data/project/work/ramax-run          # temporary work directory
```

See [Restart and work directories](restart.md) for recovery behavior.

## Logging

The default log level is `info`. Use `--log-level debug` when investigating a
candidate or graph transaction, `--quiet` for errors only, and `--verbose` for
additional diagnostics. `--verbose` and `--quiet` are mutually exclusive.
