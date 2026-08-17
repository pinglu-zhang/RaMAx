# Troubleshooting

## Input file or FASTA does not exist

RaMAx validates the seqfile and every local genome path before preprocessing.
Use absolute paths. MAF accepts mapping-only input; HAL requires a Newick tree
as the first non-empty record. When a tree is present, verify that its leaves
match the genome mappings.

```bash
test -r /data/project/seqfile.txt
test -r /data/genomes/human.fa
```

For URLs, verify network access from the environment that launches RaMAx.

## HAL reports that the Newick tree is missing

Only HAL requires a species tree. Add a valid Newick record before the genome
mappings, or change the output suffix to a non-HAL format. `--root` is also
HAL-only and is rejected for MAF and other non-HAL outputs.

## Work directory is not empty

A new Release run requires an empty work directory. Choose a new path rather
than mixing state from different experiments:

```bash
ramax -i seqfile.txt -o results/run2.maf -w work/run2 -t 16
```

Use `--restart -w <existing-work>` only for the same interrupted run. Remember
that a successful run removes its work directory automatically.

## minipoa is missing

RaMAx requires a separately installed minipoa executable for its MSA-based
realignment and shared export path. RaMAx does not download or install it.

Runtime lookup order is:

1. `RAMAX_MINIPOA_EXECUTABLE` selected when configuring the RaMAx build;
2. an executable named `minipoa` beside the running `ramax` binary;
3. `minipoa` in `PATH`.

Check the installation with:

```bash
command -v minipoa
minipoa --help
```

For a source build with an explicit executable:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DRAMAX_MINIPOA_EXECUTABLE=/opt/minipoa/bin/minipoa
```

## Restart schema is incompatible

RaMAx 1.0.5 accepts `schema_version: 1`. A missing or different version means
the work directory came from an incompatible build. Start a new run in a fresh
directory; do not add fields manually to an old configuration because its
other intermediate files may also be incompatible.

## MAF validation fails

For each Block, verify:

- equal aligned-row lengths;
- `size` equals the ungapped row length;
- strand is `+` or `-`;
- coordinates remain inside the declared source length.

Confirm that the file was completely written and that downstream name
normalization did not edit sequence columns or coordinates.

## HAL validation fails

Run:

```bash
halValidate alignment.hal
halStats alignment.hal
```

Check the first reported invalid genome/segment, confirm that the output was
not truncated, and retain the failed work directory and `RaMAx.log`. When
comparing a HAL-derived MAF to Truth, also verify the full leaf, ancestor,
root, and chromosome name mapping.

## PAF export skips Blocks or seqwish cannot find paths

Check the final `PAF export complete` summary first. Formal acceptance requires
`invalid=0`; the following warning reports the first malformed Block and its
explicit `ref_species.ref_chr` key without writing a partial Block.

Every PAF sequence name is `species.original_contig_name`. Build the seqwish
input FASTA with exactly those headers, including repeated-looking prefixes,
and use `seqwish -k 0` when comparing `connected` with `all`. Nonzero `-k`
filters short exact matches and is outside RaMAx connected-mode completeness.

Prefer `ramax-paf-fasta -i <seqfile> -o <sequences.fa.gz>` over manual header
rewriting. If it reports a duplicate qualified name, change the conflicting
species or contig name at the source; the tool intentionally does not invent
suffixes because they would no longer match PAF.

## High memory use or poor scaling

- Reduce `--threads`; the default is all reported hardware threads.
- Keep `--chunk_size` at its default unless profiling shows a reason to change
  it. Smaller chunks add scheduling and boundary overhead.
- Place work data on a fast local filesystem with sufficient free space.
- Use a Release build. `RAMAX_NATIVE_ARCH=OFF` improves binary portability but
  may reduce peak performance on the build host.

## Diagnosing a realignment failure

Use debug logging only for focused diagnosis:

```bash
ramax ... --log-level debug
```

Default `info` summarizes candidate counts, commits, failures, cache behavior,
and timing. Repeated external-MSA failures are aggregated, so inspect the final
`[external-msa]` summary as well as the first warning.
