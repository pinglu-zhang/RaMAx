# Packaging and releasing RaMAx

This document describes the current Linux x86-64 release workflow. The
commands below target RaMAx 1.0.8. Update the version-specific values together
for a later release rather than copying an older packaging guide.

## Release contract

The Conda package and Docker image contain RaMAx and the external tools needed
by the installed workflow:

```text
RaMAx      1.0.8
Mash       2.3=hb105d93_9
wfmash     0.14.0=h11f254b_0 (v0.14.0-0-g517e1bc)
Samtools   1.23.1=ha83d96e_0
HTSlib     1.23.1=h633afcb_0
minipoa    1.4.2=hd5d28ae_0
Cactus     cactus-bin-v2.9.9 halAppendCactusSubtree
```

Mash build `hb105d93_9` and wfmash 0.14 use the same GSL 2.7 series.
Samtools/HTSlib 1.23.1 is the exact pair used by the wfmash routing stage.
The Conda recipe uses HDF5 1.14, while the Ubuntu 22.04 Docker builder uses
the distribution HDF5 development package.

Packaging validation establishes that the binaries build, start, resolve their
runtime libraries, and report the expected versions. It is not evidence of
cross-platform biological equivalence, alignment quality, or performance.

## Source preflight

Run release commands from the intended release checkout:

```bash
cd /mnt/d/code/RaMAx

test "$(sed -n 's/^project(RaMAx VERSION \([^ ]*\).*/\1/p' CMakeLists.txt)" = "1.0.8"
test -x bin/halAppendCactusSubtree
git status --short
git rev-parse HEAD
```

Do not create or push the release tag until both package formats pass their
local validation.

## Build the Conda package

Use flexible channel priority so the Bioconda tools and current conda-forge
runtime dependency chain can be solved together:

```bash
source /home/zpl/miniconda3/etc/profile.d/conda.sh
conda activate conda

cd /mnt/d/code/RaMAx/recipe

export CONDA_CHANNEL_PRIORITY=flexible
export CPU_COUNT=2
export CONDA_BUILD_MAX_JOBS=2
export RAMAX_CONDA_SOURCE=path

conda build . \
  --override-channels \
  -c malab \
  -c bioconda \
  -c conda-forge \
  --no-anaconda-upload
```

`--no-anaconda-upload` keeps build and upload as separate release gates. The
expected package is build number 0:

```bash
RAMAX_CONDA_PACKAGE=$(find "$CONDA_PREFIX/conda-bld/linux-64" \
  -maxdepth 1 \
  -name 'ramax-1.0.8-*_0.conda' \
  -print \
  | sort \
  | tail -n 1)

test -n "$RAMAX_CONDA_PACKAGE"
printf '%s\n' "$RAMAX_CONDA_PACKAGE"
```

### Validate a clean Conda installation

Install the local package in a new environment outside the source tree:

```bash
CONDA_CHANNEL_PRIORITY=flexible \
conda create -y \
  -p /mnt/d/Result/RaMAx/package/ramax-1.0.8-test \
  --override-channels \
  -c "file://$CONDA_PREFIX/conda-bld" \
  -c malab \
  -c bioconda \
  -c conda-forge \
  ramax=1.0.8

conda run -p /mnt/d/Result/RaMAx/package/ramax-1.0.8-test ramax --version
conda run -p /mnt/d/Result/RaMAx/package/ramax-1.0.8-test mash --version
conda run -p /mnt/d/Result/RaMAx/package/ramax-1.0.8-test wfmash --version
conda run -p /mnt/d/Result/RaMAx/package/ramax-1.0.8-test samtools --version
conda run -p /mnt/d/Result/RaMAx/package/ramax-1.0.8-test minipoa -v
```

Expected key output:

```text
RaMAx version 1.0.8
2.3
v0.14.0-0-g517e1bc
samtools 1.23.1
Using htslib 1.23.1
minipoa version 1.4.2
```

### Upload the Conda package

Upload only after the clean installation passes:

```bash
anaconda upload --user malab "$RAMAX_CONDA_PACKAGE"
```

After upload, solve and install `ramax=1.0.8` from the remote channels in a
new environment before declaring the Conda release complete.

## Build the Docker image

The Dockerfile copies the bundled `bin/halAppendCactusSubtree` and installs the
exact companion-tool builds into `/opt/ramax`. It uses the local Docker
frontend and does not require named build contexts or user-supplied integrity
arguments.

```bash
cd /mnt/d/code/RaMAx

docker build \
  --progress=plain \
  --build-arg BUILD_JOBS=16 \
  -t ramax:1.0.8 \
  .
```

Reduce `BUILD_JOBS` when memory is limited. A normal rebuild may reuse Docker
layers; `--no-cache` is not part of the standard release command.

### Validate the local Docker image

```bash
docker run --rm ramax:1.0.8 --version

docker run --rm --entrypoint bash ramax:1.0.8 -lc '
set -euo pipefail

ramax --version
ramax-paf-fasta --version
mash --version
wfmash --version
samtools --version | head -n 2
minipoa -v
command -v halAppendCactusSubtree

for executable in ramax wfmash minipoa halAppendCactusSubtree; do
    report="$(ldd "$(command -v "$executable")" 2>&1 || true)"
    printf "%s\n" "$report"
    if printf "%s\n" "$report" | grep -q "not found"; then
        exit 1
    fi
done
'
```

## Upload to Docker Hub

Docker Hub is the canonical container registry for this release. Push the
fixed version first:

```bash
docker login -u pingluzhang

docker tag ramax:1.0.8 pingluzhang/ramax:1.0.8
docker push pingluzhang/ramax:1.0.8
```

Verify that the remote manifest exists:

```bash
docker buildx imagetools inspect pingluzhang/ramax:1.0.8
```

Only after the fixed tag passes remote pull and runtime validation should the
moving `latest` tag be updated:

```bash
docker tag ramax:1.0.8 pingluzhang/ramax:latest
docker push pingluzhang/ramax:latest
```

Users should prefer the fixed tag for reproducibility:

```bash
docker pull pingluzhang/ramax:1.0.8
docker run --rm pingluzhang/ramax:1.0.8 --version
```

## Create the release tag

Create a new annotated tag after the release commit is final, the branch is
pushed normally, and both remote package formats pass validation:

```bash
cd /mnt/d/code/RaMAx

git status --short
git tag -a v1.0.8 -m "RaMAx 1.0.8" HEAD
git push origin dev-suffix
git push origin v1.0.8
```

Do not move or force-push an older release tag. Record the final Git commit,
Conda package filename, Docker manifest digest, and validation results in the
release notes.
