# RaMAx

RaMAx aligns multiple genomes and produces whole‑genome alignments in standard formats. Its input and output conventions are **compatible with [Cactus]** (i.e., it accepts the same `seqfile` style input and can emit **MAF** and **HAL** alignments), so you can drop RaMAx into existing Cactus‑based workflows with minimal changes.

> If you already have a Cactus `seqfile`, you can pass it directly to RaMAx.

---

## Install with Conda

Conda is the recommended installation method for RaMAx. Install it from the `malab` channel:

```bash
conda install -c conda-forge -c malab ramax
ramax --help
```

The Linux conda package is built against the conda-forge sysroot and is intended to run on systems with `glibc >= 2.17`.

---

## Run with Docker

Docker is an alternative for older or heterogeneous Linux systems, especially when the host glibc is incompatible. The image includes its own userspace libraries.

Pull the published image and check that RaMAx starts:

```bash
docker pull pinglu-zhang/ramax:latest
docker run pinglu-zhang/ramax:latest --help
```

## Build from Source

### Prerequisites (Ubuntu/Debian)

Install the following packages first:

```bash 
sudo apt update
sudo apt install -y \
    build-essential \
    cmake \
    libcurl4-openssl-dev \
    zlib1g-dev \
    libhdf5-dev
```
> **Notes**
>
> * `libhdf5-dev` is required for HAL support.
> * On non‑Debian systems, install the equivalent development packages for **CMake**, **libcurl**, **zlib**, and **HDF5**.

### Build Steps

```bash
# Clone the repository
git clone https://github.com/pinglu-zhang/RaMAx
cd RaMAx

# Create a build directory
mkdir build && cd build

# Configure
cmake ..

# Build using all available CPU cores
make -j"$(nproc)"
```

After a successful build, the executable **`ramax`** will be located in the `build/` directory.

---

## Quick Start

### Basic usage

```bash
./ramax \
  -i <path/to/seqfile.txt> \
  -o <path/to/output.{maf|hal}> \
  -w <path/to/workdir> \
  -t <threads>
```

---

## Command‑line Options (common)

| Flag        | Description                                                                  |
| ----------- | ---------------------------------------------------------------------------- |
| `-i <file>` | Input **seqfile** (same format used by Cactus).                              |
| `-o <file>` | Output alignment path. RaMAx supports **MAF** and **HAL**.                   |
| `-w <dir>`  | Working directory to store intermediate files and logs (created if missing). |
| `-t <N>`    | Number of worker threads to use.                                             |

**Input format (seqfile)**

RaMAx uses the same `seqfile` format as Cactus. The file has two parts:

1. **A species tree in Newick** on the first non‑comment line (leaf names must match genomes below; branch lengths optional).
2. **Genome mappings**, one per line: `<genome_name> <path/to/genome.{fa|fasta|fa.gz}>`.

**Minimal example**

```text
# seqfile.txt
# 1) Species tree (Newick)
((human:0.005,chimp:0.005):0.02,orangutan:0.025);

# 2) Genome mappings (one per leaf)
human      /data/genomes/hg38.fa
chimp      /data/genomes/panTro6.fa
orangutan  /data/genomes/ponAbe3.fa.gz
```

> **Output formats**
> RaMAx can write **MAF** for broad downstream compatibility and **HAL** for hierarchical alignment workflows.
