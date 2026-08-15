# Output formats

RaMAx selects its exporter from the output filename suffix:

- `.maf` writes a Multiple Alignment Format file;
- `.hal` writes a Hierarchical Alignment Format file;
- other suffixes are rejected before alignment begins.

## MAF

MAF output does not require a species tree. The seqfile may begin directly with
genome-to-FASTA mappings. If a tree is supplied for compatibility, it does not
control alignment reference ordering.

Each MAF alignment Block begins with an `a` line followed by one or more `s`
rows. RaMAx writes standard zero-based, half-open coordinates:

```text
s source start size strand sourceSize alignedText
```

- `size` equals the number of non-gap characters in `alignedText`;
- all rows in one Block have equal aligned-text length;
- `strand` is `+` or `-`;
- `start + size` must not exceed `sourceSize`;
- removing `-` from a row recovers the represented genomic sequence in its
  alignment orientation.

## HAL

HAL output stores the same optimized graph in a hierarchy defined by the
seqfile tree. A valid Newick tree is therefore mandatory for HAL. RaMAx rejects
a mapping-only seqfile before preprocessing begins. RaMAx reconstructs ancestor
sequences and writes leaf/ancestor segments without changing the source FASTA
coordinates. `--root` controls the preferred artificial root name and is valid
only for HAL output.

Validate a HAL file with the HAL utilities:

```bash
halValidate /data/project/results/alignment.hal
halStats /data/project/results/alignment.hal
```

## Shared MAF/HAL realignment

Direct MAF export and HAL construction both use the shared
`mergeAlignmentByRef()` reference projection. Nearby homologous insertions may
be repaired during export, but failure or lack of strict improvement retains
the original rows. HAL-specific code does not implement a separate insertion
repair algorithm.

## Comparing MAF results

Use normalized MAF files for reproducible comparisons. Apply the same
normalization to every candidate and Truth file, including a stable source-name
mapping, Block/row ordering policy, and whitespace policy. Do not compare a
normalized candidate against an unnormalized Truth file.

For sampled Truth evaluation, keep the comparator version, sample count, seed,
and near-distance setting fixed. RaMAx development comparisons use
`Overall (w/o self)` as the primary F-score rather than inferring accuracy from
Block count alone.

## Temporary files

External MSA input normally uses an in-memory file descriptor. Any fallback
input file and minipoa internal file is confined to:

```text
<work>/minipoa_tmp/
```

RaMAx removes per-process minipoa files after successful execution, non-zero
exit, or output-parse failure. A forced process termination can leave files in
this directory, but should not create them beside the output MAF/HAL or in the
launch directory.
