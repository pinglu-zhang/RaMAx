# RaMAx

RaMAx aligns multiple genomes and writes whole-genome alignments in MAF, HAL,
or PAF format. A seqfile contains one genome name and FASTA path per line. HAL
output additionally requires a Newick species tree as the first record.

## Install

### Conda

```bash
conda install -c conda-forge -c malab ramax
ramax --version
```

### Docker

```bash
docker pull pingluzhang/ramax:latest
docker run --rm pingluzhang/ramax:latest --version
```

### Build from source

RaMAx requires a C++23 compiler, CMake 3.20 or newer, zlib, libcurl, TBB,
OpenMP, HDF5, and the bundled HAL/sonLib sources. On Ubuntu or Debian:

```bash
sudo apt update
sudo apt install -y build-essential cmake libcurl4-openssl-dev zlib1g-dev \
  libtbb-dev libhdf5-dev
```

Configure, build, and install:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"
cmake --install build --prefix "$HOME/.local"
```

Install `minipoa` separately and place it in `PATH`, next to the installed
`ramax` executable, or provide its path during configuration:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DRAMAX_MINIPOA_EXECUTABLE=/opt/minipoa/bin/minipoa
```

RaMAx also requires Mash 2.3 at runtime. Place `mash` in `PATH`, next to the
installed `ramax`, or configure its exact location with:

```bash
cmake -S . -B build \
  -DRAMAX_MASH_EXECUTABLE=/opt/mash/bin/mash
```

After selecting the first reference, RaMAx records whole-genome Mash distances
using `k=31` and sketch size 20,000 before starting the legacy aligner.

Set `-DRAMAX_NATIVE_ARCH=OFF` for a portable x86-64 RaMAx build. Configuration
and compilation do not write to `/usr/local`; installation is controlled only
by `cmake --install` and its prefix.

## Quick start

```bash
ramax \
  -i seqfile.txt \
  -o alignment.maf \
  -w work \
  -t 16 \
  --optimize-blocks
```

Use an output name ending in `.maf`, `.hal`, `.paf`, or `.gfa`. PAF defaults to the
information-complete sparse `connected` mode; use `--paf-mode all` for the
all-pairs baseline. Logs and reusable preprocessing/index caches are kept under
the work directory; in-memory graph state is not serialized for restart.

Repeat `-o` to export several formats from the same completed alignment:

```bash
ramax -i seqfile.txt -w work -t 16 \
  -o alignment.maf \
  -o alignment.paf \
  -o alignment.hal
```

Each format may appear once. If HAL is requested anywhere in the output list,
the seqfile must contain a Newick tree. Export order is MAF, PAF, GFA, then HAL.
Native GFA export defaults to GFA 1.1 W-lines; use `--gfa-version 1.0` for
GFA 1.0 P-line compatibility. The graph profile currently defaults to the
audit-preserving `--gfa-profile exact`; `--gfa-profile compact` enables the
lossless compact-v1 transforms and keeps an exact shadow under `work/gfa/`.

For seqwish, generate the matching qualified FASTA directly from the same
seqfile. Keep `-k 0` to preserve the aligned-base relationships guaranteed by
the default `connected` PAF mode:

```bash
ramax-paf-fasta -i seqfile.txt -o sequences.fa.gz
ramax -i seqfile.txt -w work -t 16 \
  -o alignment.paf \
  -o alignment.maf
seqwish -s sequences.fa.gz -p alignment.paf -g graph.gfa -k 0
```

For a normalized PGGB graph, write an uncompressed FASTA, index it, and reuse
the RaMAx PAF with `pggb -a` so PGGB skips wfmash:

```bash
ramax-paf-fasta -i seqfile.txt -o sequences.fa
samtools faidx sequences.fa
mkdir -p pggb-out pggb-tmp
pggb -i sequences.fa -o pggb-out -a alignment.paf -n <genome-count> \
  -t 16 -T 8 -k 0 -D pggb-tmp
```

The detailed operator guide covers directory layout, validation, failure
handling, the direct seqwish route, the PGGB route, and a completed
seven-genome Chr09 example: [PAF-to-GFA workflow](docs/paf-to-gfa-workflow.md).

The graph optimizations are enabled by default. `--optimize-blocks` is an
explicit, repeatable way to request the same default set. The defaults are:

```text
--merge-blocks --merge-gap 100
--realign-missing --realign-span 3000 --zero-gap-span 200
--repair-breaks --break-span 1000
--merge-short-blocks
```

Run `ramax --help` for the complete core-alignment options.

## Documentation

- [Usage guide](docs/usage.md)
- [Command-line parameters](docs/parameters.md)
- [Restart and work directories](docs/restart.md)
- [MAF, HAL, and PAF output](docs/output-formats.md)
- [PAF-to-GFA workflow](docs/paf-to-gfa-workflow.md)
- [Troubleshooting](docs/troubleshooting.md)
- [Realignment architecture](docs/realignment_architecture.md)

See the [documentation index](docs/README.md) for the complete set.

## Input format

MAF and PAF output accept a mapping-only seqfile:

```text
human      /data/genomes/human.fa
chimp      /data/genomes/chimp.fa
orangutan  /data/genomes/orangutan.fa.gz
```

HAL output requires a Newick tree as the first record:

```text
((human:0.005,chimp:0.005):0.02,orangutan:0.025);
human      /data/genomes/human.fa
chimp      /data/genomes/chimp.fa
orangutan  /data/genomes/orangutan.fa.gz
```

Tree input remains optional for MAF and other non-HAL formats. When a tree is
provided, its leaf names should match the genome mappings. `--root` is valid
only for HAL output.

## Restart compatibility

RaMAx 1.0.7 treats restart as cache reuse, not alignment checkpointing. It
reuses validated raw/clean FASTA and FM-index artifacts, then reruns anchor
search, clustering, graph construction, and export from the beginning:

```bash
ramax --restart -w work

# Explicit settings override the latest saved values.
ramax --restart -w work -t 24 --min_anchor_length 30 \
  -o results/retry.maf -o results/retry.paf
```

The seqfile cannot be replaced. Schema-1 through schema-4 work directories are
loaded with explicit compatibility defaults and migrated to schema 5; older
work directories use `--gfa-profile exact`, and schemas before 4 default to
GFA 1.1 when no GFA version was persisted.

## License and dependencies

RaMAx is released under the MIT License. See [LICENSE](LICENSE) and
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
