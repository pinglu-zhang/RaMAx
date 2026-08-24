# Packaging RaMAx 1.0.7

The Conda package and Docker image are Linux x86-64 releases. They contain
RaMAx 1.0.7, minipoa 1.4.2 from the malab Conda channel, and the precompiled
`halAppendCactusSubtree` from cactus-bin-v2.9.9. Mash, wfmash, Samtools, and
HTSlib use exact Conda builds.

## Fixed tools

```text
Mash       2.3=hb105d93_10
wfmash     0.14.0=h11f254b_0 (v0.14.0-0-g517e1bc)
Samtools   1.23.1=ha83d96e_0
HTSlib     1.23.1=h633afcb_0
minipoa    1.4.2=hd5d28ae_0 (malab)
Cactus     cactus-bin-v2.9.9 helper
```

The packaged helper has SHA-256
`e36439892f7f84e4a8cd86e8e22d9650b0f68ecf2baf0d964f4e9ef03b8a0ab5`.
It is a statically linked Linux x86-64 ELF, so `ldd` is expected to report
`not a dynamic executable`; package validation starts `--help` instead.

Samtools/HTSlib 1.23.1 is the newest release compatible with the exact
wfmash 0.14 Bioconda build. RaMAx uses only
`samtools faidx --fai-idx`, whose behavior is unchanged from 1.24.

## Prepare the Cactus helper

Download the precompiled helper into the repository and validate it before
packaging:

```bash
cd /mnt/d/code/RaMAx
mkdir -p bin
scp -P 50623 \
  zhangpinglu@bh.mmszxc.xin:/mnt/sda/tianqinzhong/software/cactus-bin-v2.9.9/bin/halAppendCactusSubtree \
  bin/halAppendCactusSubtree
chmod 0755 bin/halAppendCactusSubtree

export RAMAX_CACTUS_HELPER=/mnt/d/code/RaMAx/bin/halAppendCactusSubtree

test -x "${RAMAX_CACTUS_HELPER}"
file "${RAMAX_CACTUS_HELPER}"
ldd "${RAMAX_CACTUS_HELPER}"
readelf -d "${RAMAX_CACTUS_HELPER}"
objdump -T "${RAMAX_CACTUS_HELPER}" \
  | grep -o 'GLIBC_[0-9.]*' \
  | sort -Vu
sha256sum "${RAMAX_CACTUS_HELPER}"
```

The build fails if the helper has missing libraries, resolves a private
`/mnt/sda` library, or requires a GLIBC version newer than 2.17 for Conda.
A static helper must start and print its `USAGE:` help text. No Cactus library
closure is copied implicitly.

## Build the Conda package

```bash
export RAMAX_CONDA_SOURCE=path
export RAMAX_CACTUS_HELPER_SHA256="$(sha256sum /mnt/d/code/RaMAx/bin/halAppendCactusSubtree | awk '{print $1}')"

conda build /mnt/d/code/RaMAx/recipe \
  -c malab \
  -c conda-forge \
  -c bioconda \
  --output-folder /mnt/d/Result/RaMAx/package/conda-bld
```

Install and test the local package:

```bash
conda create -y \
  -p /mnt/d/Result/RaMAx/package/ramax-1.0.7-test \
  -c file:///mnt/d/Result/RaMAx/package/conda-bld \
  -c malab \
  -c conda-forge \
  -c bioconda \
  ramax=1.0.7

conda run -p /mnt/d/Result/RaMAx/package/ramax-1.0.7-test ramax --version
conda run -p /mnt/d/Result/RaMAx/package/ramax-1.0.7-test wfmash --version
conda run -p /mnt/d/Result/RaMAx/package/ramax-1.0.7-test samtools --version
conda run -p /mnt/d/Result/RaMAx/package/ramax-1.0.7-test minipoa -v
```

After validation, upload manually:

```bash
anaconda upload \
  --user malab \
  /mnt/d/Result/RaMAx/package/conda-bld/linux-64/ramax-1.0.7-*.conda
```

## Build the Docker image

The Dockerfile installs minipoa from the malab Conda channel and copies the
validated `bin/halAppendCactusSubtree` included in the release source tree.
Micromamba creates the complete environment directly at `/opt/ramax`; RaMAx,
Mash, minipoa, wfmash, Samtools, and the Cactus helper therefore all reside in
`/opt/ramax/bin`, while their required shared libraries remain in
`/opt/ramax/lib`. The image copies this complete prefix instead of copying
individual Conda executables without their libraries.

```bash
cd /mnt/d/code/RaMAx

export RAMAX_CACTUS_SHA256="$(sha256sum bin/halAppendCactusSubtree | awk '{print $1}')"

docker buildx build \
  --platform linux/amd64 \
  --build-arg BUILD_JOBS=16 \
  --build-arg RAMAX_VERSION=1.0.7 \
  --build-arg VCS_REF="$(git rev-parse HEAD)" \
  --build-arg CACTUS_HELPER_SHA256="${RAMAX_CACTUS_SHA256}" \
  -t ramax:1.0.7 \
  --load \
  .
```

Validate the runtime:

```bash
docker run --rm ramax:1.0.7 --version

docker run --rm --entrypoint bash ramax:1.0.7 -lc '
set -e
mash --version
wfmash --version
samtools --version | head -2
minipoa -v
command -v halAppendCactusSubtree
ldd "$(command -v halAppendCactusSubtree)" || true
halAppendCactusSubtree --help 2>&1 | grep -q 'USAGE:'
'
```

Push only after local validation:

```bash
docker buildx build \
  --platform linux/amd64 \
  --build-arg BUILD_JOBS=16 \
  --build-arg RAMAX_VERSION=1.0.7 \
  --build-arg VCS_REF="$(git rev-parse HEAD)" \
  --build-arg CACTUS_HELPER_SHA256="${RAMAX_CACTUS_SHA256}" \
  -t ghcr.io/pinglu-zhang/ramax:1.0.7 \
  -t ghcr.io/pinglu-zhang/ramax:latest \
  --push \
  .
```

## Final release tag

The current v1.0.7 tag predates the final packaging commit. Move it only
after both package formats pass validation. Do not force-push a branch.

```bash
git tag -d v1.0.7
git tag -a v1.0.7 -m "RaMAx 1.0.7" HEAD
git push origin dev-pan:dev-pan
git push --force origin refs/tags/v1.0.7
```
