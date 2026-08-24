# ramax-paf-fasta

`ramax-paf-fasta` converts a RaMAx seqfile into the single qualified FASTA
needed beside RaMAx PAF output:

```bash
ramax-paf-fasta -i genomes.seqfile -o genomes.fa.gz
ramax -i genomes.seqfile -o alignment.paf -w work
seqwish -s genomes.fa.gz -p alignment.paf -g graph.gfa -k 0
```

Headers are written as `species.contig`, where `contig` is the first field of
the source FASTA header. Sequences use the same normalization as RaMAx PAF:
uppercase `A/T/G/C` are preserved and every other base symbol becomes `N`.

The output suffix must be `.fa`, `.fasta`, or `.fna`, optionally followed by
`.gz`. Existing outputs are rejected unless `--force` is supplied.

## Standalone build

```bash
cmake -S tools/ramax-paf-fasta -B tools/ramax-paf-fasta/build \
  -DCMAKE_BUILD_TYPE=Release
cmake --build tools/ramax-paf-fasta/build -j
```

When this source-local binary already exists, the root RaMAx CMake build
detects it, skips rebuilding the tool, and installs the existing executable.
Use `-DRAMAX_PAF_FASTA_PREBUILT=/path/to/ramax-paf-fasta` for a binary built in
another directory, or `-DRAMAX_BUILD_TOOLS=OFF` to omit companion tools.
