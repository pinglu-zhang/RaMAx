# Output formats

RaMAx selects exporters from the suffixes of one or more repeated `-o` paths:

- `.maf` writes a Multiple Alignment Format file;
- `.hal` writes a Hierarchical Alignment Format file;
- `.paf` writes pairwise projections of the primary Block alignments;
- `.gfa` writes a native lossless sequence graph;
- other suffixes are rejected before alignment begins.

Each format may appear once. RaMAx constructs the optimized graph once and
exports in the fixed order MAF, PAF, GFA, HAL. Each file is published independently
through a temporary file. If one exporter fails, later formats are still
attempted; successful files remain, the command exits nonzero, and the work
directory is retained.

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

Leaf sequence names preserve the sequence IDs read from the source FASTA. RaMAx
does not add a species prefix, remove an existing prefix, or otherwise rewrite
those IDs. Sequence names are scoped by their HAL genome, so different leaf
genomes may each contain a sequence such as `chr1`. Internally, RaMAx still uses
species-qualified keys to keep those sequences distinct. Reconstructed ancestor
sequences follow the Cactus-compatible, zero-based naming convention
`<genome>refChr<N>`, for example `rootrefChr0` and `anc1refChr0`.

This HAL naming policy does not change MAF source names. Existing HAL files are
also unchanged; the policy applies when a new HAL is exported.

Validate a HAL file with the HAL utilities:

```bash
halValidate /data/project/results/alignment.hal
halStats /data/project/results/alignment.hal
```

## PAF

PAF output does not require a species tree. Sequence names are always
`species.original_contig_name`; RaMAx does not remove an existing species
prefix from a contig. The FASTA passed to a downstream tool must use exactly
the same qualified headers. Empty names, whitespace, and qualified-name
collisions are rejected before export.

Generate that FASTA without manually rewriting headers:

```bash
ramax-paf-fasta -i genomes.seqfile -o genomes.fa.gz
```

The tool streams plain or gzip inputs, preserves seqfile/contig order, converts
sequences to uppercase `A/T/G/C/N` with every non-ATGC symbol represented as
`N`, and publishes the result atomically. It refuses to overwrite an existing
file unless `--force` is supplied.

`--paf-mode connected` is the default. It first selects the explicit Block
reference against every overlapping primary row, then deterministically adds
pairs until every non-gap, non-`X`, case-normalized same-base set in every MSA
column is connected. `--paf-mode all` emits all primary row pairs sharing at
least one non-gap column and is the correctness baseline. A pair selected in a
Block is emitted once as a complete Block projection.

Records contain the standard 12 PAF fields followed by `cg:Z:`, `tp:A:P`, and
`NM:i:` tags. MAPQ is 255 and CIGAR operations are limited to `=`, `X`, `I`,
and `D`. Coordinates always refer to the forward source sequences; the PAF
strand is the XOR of the two Segment orientations. Column 10 counts `=` bases,
column 11 counts non-double-gap columns, and `NM` counts `X+I+D` bases.

Malformed Blocks are skipped as complete units and summarized once after
export. Treat any nonzero `invalid` count as a failed scientific acceptance
check. Global naming or I/O errors fail the export without replacing an
existing output file.

The connected-mode completeness guarantee applies to the homologous
same-base relations already present in valid primary Blocks and to
`seqwish -k 0`. A nonzero seqwish `-k` intentionally removes short exact
matches and therefore changes that relation set. See the seqwish
[graph-induction algorithm](https://github.com/pangenome/seqwish#squish-graph-induction-algorithm)
and [`-k` usage](https://github.com/pangenome/seqwish#usage).

## GFA

Native GFA export uses the same cleaned FASTA and primary Block projection as
the other direct exporters. `--gfa-version 1.1` is the default and emits one
W-line per input contig. `--gfa-version 1.0` emits the same logical walks as
P-lines for compatibility with GFA 1.0 consumers. The option is rejected when
the output list has no `.gfa` path.

`--gfa-profile exact` is the current default and preserves every maximal exact
run as the audit graph. `--gfa-profile compact` enables
compact-v2-balanced: exact runs shorter than the selected balanced threshold
are left private, occurrence-equivalent chains are unitigged, and conservative
dense small-variant intervals may be rewritten as observed compound alleles.
The profile favors SV accuracy among configurations with PGGB-like graph
granularity. It never invents an allele or removes sequence from an input path.
Compact export keeps the same-version exact graph at `work/gfa/exact.gfa` and
writes `compact_transform.tsv`, `compact_stats.tsv`,
`compact_rejections.tsv`, and `compact_parameters.tsv` beside it.

`compact_stats.tsv` reports the exact graph, short-relation filtering, first
unitig pass, compound-allele pass, and final unitig pass separately. The
effective fixed policy and profile version are recorded in
`compact_parameters.tsv`; these implementation parameters are not runtime CLI
options.

Both versions use byte-identical S-lines, L-lines, deterministic numeric node
IDs, and `0M` link overlap. They differ only in path serialization:

- GFA 1.0 writes `H VN:Z:1.0` and path names `species.contig`;
- GFA 1.1 writes `H VN:Z:1.1 RS:Z:<first-reference>` and structured W fields
  `SampleId=species`, `HapIndex=0`, `SeqId=contig`, `SeqStart=0`, and
  `SeqEnd=contig-length`.

The exporter divides each cleaned contig at propagated exact-run boundaries.
Only homologous, column-supported, identical A/C/G/T runs are shared; N bases,
unaligned intervals, and private sequence remain private. It never merges
unrelated repeats merely because their strings are equal. Before publication,
every logical walk is reconstructed in memory and compared base-for-base and
length-for-length with its cleaned input contig. Thus every input base occurs
once in its walk, including contigs with no alignment Block.

Compact transforms additionally retain at least 99.5% of the exact graph's
direction-aware homology mass and limit graph sequence growth to 2%. A budget
or invariant failure leaves the exact shadow available, does not publish the
requested compact GFA, and makes RaMAx exit nonzero.

GFA 1.0 rejects a collision between qualified `species.contig` path names.
GFA 1.1 keeps species and contig in separate W fields and can represent that
case as long as each field is valid. A failed GFA export does not publish its
temporary file; other successful output formats remain available and RaMAx
exits nonzero.

## Shared Block projection

Direct MAF/PAF/GFA export and HAL construction use the shared
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
this directory, but should not create them beside the output MAF/HAL/PAF/GFA or in the
launch directory.
