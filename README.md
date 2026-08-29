# RaMAx

RaMAx aligns multiple genomes and writes whole-genome alignments in MAF, HAL,
PAF, or GFA format. A seqfile contains one genome name and FASTA path per line.
HAL output additionally requires a Newick species tree as the first record.

## Install

### Conda

```bash
conda install -c malab -c conda-forge -c bioconda ramax=1.0.8 -y

ramax --version
```

### Docker

```bash
docker pull pingluzhang/ramax:1.0.8
docker run --rm pingluzhang/ramax:1.0.8 --version
```

## Quick start

```bash
ramax -i seqfile.txt -o alignment.maf -w work -t 16 

```

Use an output name ending in `.maf`, `.hal`, `.paf`, or `.gfa`. PAF defaults to the
information-complete sparse `connected` mode; use `--paf-mode all` for the
all-pairs baseline. Logs and reusable preprocessing/index caches are kept under
the work directory, together with intermediate graph state and minipoa scratch
files; in-memory graph state is not serialized for restart.
MAF and HAL use the same normalized multiway homology relation, so changing
only the output suffix does not change leaf-to-leaf alignment coverage. HAL
soft-mask restoration happens after alignment decisions; letter case cannot
change the homology graph.

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
lossless compact-v2-balanced transforms and keeps an exact shadow under
`work/gfa/`.


```bash
ramax-paf-fasta -i seqfile.txt -o sequences.fa.gz
ramax -i seqfile.txt -w work -t 16 \
  -o alignment.paf \
  -o alignment.maf
seqwish -s sequences.fa.gz -p alignment.paf -g graph.gfa
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

### Build from source

RaMAx requires a C++23 compiler, CMake 3.20 or newer, zlib, libcurl, TBB,
OpenMP, HDF5, and the bundled HAL/sonLib sources. On Ubuntu or Debian:

```bash
sudo apt update
sudo apt install -y build-essential cmake libboost-graph-dev \
  libcurl4-openssl-dev zlib1g-dev \
  libtbb-dev libhdf5-dev
```

Configure, build, and install:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"
cmake --install build --prefix "$HOME/.local"
```

Before creating a work directory or reading input genomes, RaMAx requires
`minipoa`, PGGB-compatible wfmash `v0.14.0-0-g517e1bc`, and Mash 2.3.
`halAppendCactusSubtree` is required only when the output list contains HAL;
otherwise a missing HAL helper produces a warning and the run continues.
The source build uses `RAMAX_TOOL_BIN_DIR` as a shared external-tool
directory. It defaults to the source tree's `bin/` and expects `mash`,
`minipoa`, `wfmash`, `samtools`, and `halAppendCactusSubtree` there:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DRAMAX_TOOL_BIN_DIR="$PWD/bin"
cmake --build build --parallel 16
```

The explicit `RAMAX_TOOL_BIN_DIR` argument can be omitted for the standard
source layout. Runtime lookup order is a per-tool CMake override, the shared
tool directory, the directory containing the running `ramax` executable,
then `PATH`. Per-tool overrides remain available:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DRAMAX_HAL_APPEND_CACTUS_SUBTREE_EXECUTABLE=/opt/cactus/bin/halAppendCactusSubtree \
  -DRAMAX_MINIPOA_EXECUTABLE=/opt/minipoa/bin/minipoa \
  -DRAMAX_WFMASH_EXECUTABLE=/opt/wfmash/bin/wfmash \
  -DRAMAX_MASH_EXECUTABLE=/opt/mash/bin/mash \
  -DRAMAX_SAMTOOLS_EXECUTABLE=/opt/samtools/bin/samtools
```

To move a source-built bundle to another machine, place `ramax` and the five
tools in the same destination `bin/`; sibling lookup remains available if the
original configured directory no longer exists. CMake records tool locations
but does not make dynamically linked executables self-contained. A Conda
executable must retain its complete prefix/library closure or be reinstalled
on the destination. The bundled `halAppendCactusSubtree` is static and can be
copied directly.

Missing unconditional startup dependencies are reported together and RaMAx
exits before creating or modifying the work directory. A missing HAL helper
also stops HAL runs before normal work-directory initialization. `ramax --help` and
`ramax --version` remain available without these tools. Mash and wfmash
retain their strict version checks before use. Samtools/HTSlib 1.23.1 remains
required by the wfmash routing stage and can be configured with
`RAMAX_SAMTOOLS_EXECUTABLE`.
The official Conda and Docker packages install minipoa 1.4.2 from the
`malab` channel and bundle the validated cactus-bin-v2.9.9
`halAppendCactusSubtree` helper.
Conda builds set `RAMAX_EMBED_TOOL_PATHS=OFF`; the installed executable finds
all dependencies beside itself in the active environment's `$PREFIX/bin`
without retaining conda-build's temporary path.

After selecting the first reference, RaMAx records whole-genome Mash distances
using `k=31` and sketch size 20,000 before starting the legacy aligner.

Set `-DRAMAX_NATIVE_ARCH=OFF` for a portable x86-64 RaMAx build. Configuration
and compilation do not write to `/usr/local`; installation is controlled only
by `cmake --install` and its prefix.


## Restart compatibility

RaMAx 1.0.8 treats restart as cache reuse, not alignment checkpointing. It
reuses validated raw/clean FASTA artifacts, rebuilds the complete suffix-array
index, then reruns anchor search, clustering, graph construction, and export
from the beginning. On Linux, the complete SA/ISA/LCP buffers use temporary
file-backed mappings so the resource manager can reclaim cold pages; they are
not a persistent index cache:

```bash
ramax --restart -w work

# Explicit settings override the latest saved values.
ramax --restart -w work -t 24 --min_anchor_length 30 \
  --sa-sampling-rate 1 \
  -o results/retry.maf -o results/retry.paf
```

The seqfile cannot be replaced. Schema-1 through schema-5 work directories are
loaded with explicit compatibility defaults and migrated to schema 6; all
older workdirs use `--sa-sampling-rate 1`; the suffix-array backend accepts
only this complete-array value. References below 1024 MiB use divsufsort,
while references at or above the threshold use CaPS. The suffix-array index is
never restored across runs and its temporary mapping is unlinked immediately;
restart always rebuilds it. Older work directories
use `--gfa-profile exact`, and schemas before 4 default to
GFA 1.1 when no GFA version was persisted.

## License and dependencies

RaMAx is released under the MIT License. See [LICENSE](LICENSE) and
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
