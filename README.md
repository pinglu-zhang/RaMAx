# RaMAx

RaMAx aligns multiple genomes and writes whole-genome alignments in MAF or HAL
format. It accepts the Cactus-style seqfile format: a Newick species tree
followed by one genome name and FASTA path per line.

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

Use an output name ending in `.maf` or `.hal`. Intermediate graph state, logs,
and minipoa scratch files are kept under the work directory.

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
- [MAF and HAL output](docs/output-formats.md)
- [Troubleshooting](docs/troubleshooting.md)
- [Realignment architecture](docs/realignment_architecture.md)

See the [documentation index](docs/README.md) for the complete set.

## Input format

```text
((human:0.005,chimp:0.005):0.02,orangutan:0.025);
human      /data/genomes/human.fa
chimp      /data/genomes/chimp.fa
orangutan  /data/genomes/orangutan.fa.gz
```

Leaf names in the tree must match the genome names below it.

## Restart compatibility

RaMAx 1.0.5 stores all restart settings in `<work>/config.json` with
`schema_version: 1`. Restart with:

```bash
ramax --restart -w work
```

Restart directories created by older experimental configurations are not
compatible with 1.0.5 and are rejected with an explicit schema error.

## License and dependencies

RaMAx is released under the MIT License. See [LICENSE](LICENSE) and
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
