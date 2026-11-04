# Packaging and releasing RaMAx 1.1.0

The release is built from committed `dev` code, merged into `master`, and
identified by the annotated `v1.1.0` tag. Experimental server `dev-pan` worktree
changes are excluded. Use a clean, independent checkout for release work.

## Supported distributions

- Conda: `malab/ramax`, version `1.1.0`, build number `0`, `linux-64`.
- Docker Hub: `pingluzhang/ramax:1.1.0`, `linux/amd64`.
- GitHub Release: source archives, the Conda package, `SHA256SUMS` and
  `BUILD_INFO.json`. There is no standalone binary archive.

The runtime tool versions remain Mash 2.3 (`hb105d93_9`), wfmash 0.14.0
(`h11f254b_0`), Samtools 1.23.1 (`ha83d96e_0`), HTSlib 1.23.1
(`h633afcb_0`) and minipoa 1.4.2 (`hd5d28ae_0`). HAL is linked into RaMAx;
no `halAppendCactusSubtree` binary is needed. Conda targets the glibc 2.17
sysroot. Both packages disable `-march=native`.

## Validation before tagging

The `Release validation and publication` workflow runs on pull requests to
`dev` and `master`. PR jobs never publish. They validate version consistency,
build the committed source, run release metadata and output-validator tests,
and exercise all four output formats on a deterministic three-genome fixture.
The old internal `qa/` directory and its CMake targets are removed. Release
validation helpers live under `ci/` and are not installed in the packages. GFA walks and HAL leaf DNA are compared with the complete original
contigs; MAF sequence rows and PAF coordinate/CIGAR consumption are checked.
An independent HAL reader checks the generated hierarchy.

Conda is built with flexible channel priority using `malab`, `bioconda` and
`conda-forge`. It is installed into a fresh environment before the integration
fixture is run. Docker is built without pushing, then tested through a bind
mount as the runner's non-root user. The tested image is retained with
`docker save`; publishing loads this exact artifact rather than rebuilding it.

The source version in `CMakeLists.txt`, `recipe/meta.yaml` and `Dockerfile`
must agree. Update the current installation documentation and changelog entry;
never rename versions in historical changelog entries.

For an authorized local packaging run in a disposable checkout:

```bash
RAMAX_CONDA_SOURCE=path conda build recipe \
  -c malab -c bioconda -c conda-forge --no-anaconda-upload \
  --output-folder /path/outside/the/source/conda-bld

docker build --build-arg RAMAX_VERSION=1.1.0 \
  --build-arg VCS_REF="$(git rev-parse HEAD)" \
  --build-arg BUILD_JOBS=2 -t ramax:1.1.0 .
```

Build and test before creating the release tag. Use the PR merge candidates to
check the final source tree before merging preparation into `dev`, then `dev`
into `master`. Preserve branch history; do not force-update branches or tags.

## Credentials and publication

Configure repository Secrets `DOCKERHUB_TOKEN` (write access for `pingluzhang`)
and `ANACONDA_API_TOKEN` (upload access to `malab`) before publication. Use access
tokens through the corresponding tools/Secrets, never in files or logs.
GitHub Release publication uses the workflow's `GITHUB_TOKEN` with contents
write permission. Package or branch-protection approvals, if configured on the
repository, must be satisfied through their normal GitHub controls.

After the final master commit is validated:

```bash
git tag -a v1.1.0 -m "RaMAx 1.1.0" <validated-master-commit>
git push origin v1.1.0
```

The tag workflow rebuilds and validates the frozen tagged source, retains both
package artifacts, creates a Draft Release and publishes versioned Conda and
Docker artifacts. It then performs fresh public-channel installation and a
Docker Hub pull, repeats the integration checks, attaches checksums and source/
image identity, advances Docker `1.2` and `latest`, and publishes the Release.
Existing versioned artifacts are reused only when their content identity
matches; different same-named artifacts are never overwritten.

## Resume after credentials or network failures

Keep the failed run and its artifacts. Select the `Resume a validated release`
workflow, give `v1.1.0` and the original tag-validation run ID. It verifies the
source identity, validation records and retained checksums, and publishes those
artifacts without rebuilding. The original three build/validation jobs must
have succeeded. A draft or partial publication is not a completed release.
If the run's retained artifacts have expired, stop and review the recovery
rather than silently replacing published versioned content.

## Completion

Check all public entry points:

```bash
conda create -n ramax120 -c malab -c conda-forge -c bioconda ramax=1.1.0
conda run -n ramax120 ramax --version
docker pull pingluzhang/ramax:1.1.0
docker run --rm pingluzhang/ramax:1.1.0 --version
```

Record the final master/dev commits, annotated tag, Release URL, Conda package
SHA-256, Docker digest and validation reports. All public installations must
pass before reporting the release complete. Packaging checks do not establish
new biological accuracy or performance claims.
