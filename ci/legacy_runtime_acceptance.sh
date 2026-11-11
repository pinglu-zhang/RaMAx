#!/usr/bin/env bash
# Run inside manylinux2014 with a previously installed, 2.17-targeted package.
set -euo pipefail
: "${RAMAX_PREFIX:?}" "${RAMAX_VERSION:?}" "${EVIDENCE:?}"
mkdir -p "$EVIDENCE"
exec > >(tee "$EVIDENCE/runtime.log") 2>&1
test "$(getconf GNU_LIBC_VERSION)" = 'glibc 2.17'
getconf GNU_LIBC_VERSION
uname -a
unset LD_PRELOAD LD_LIBRARY_PATH
export PATH="$RAMAX_PREFIX/bin:/opt/python/cp311-cp311/bin:/opt/rh/devtoolset-10/root/usr/bin:/usr/bin:/bin"
python=/opt/python/cp311-cp311/bin/python
test -x "$python"
files=()
for tool in ramax ramax-spawn-helper ramax-paf-fasta mash wfmash samtools minipoa; do
    test -x "$RAMAX_PREFIX/bin/$tool"
    files+=("$RAMAX_PREFIX/bin/$tool")
done
files+=("$RAMAX_PREFIX/libexec/ramax/ramax_release_hal_check")
"$python" /src/ci/check_glibc.py --maximum 2.17 --closure \
    --output "$EVIDENCE/elf-closure.json" "${files[@]}"
# Compile the production launcher modules with native old headers and libc.
"$python" /src/ci/external_msa_spawn_test.py --output "$EVIDENCE/spawn" \
    --minipoa "$RAMAX_PREFIX/bin/minipoa" \
    --cc /opt/rh/devtoolset-10/root/usr/bin/gcc \
    --cxx /opt/rh/devtoolset-10/root/usr/bin/g++
mash --version
wfmash --version
samtools --version
minipoa -v
"$python" /src/ci/release_smoke.py --ramax "$RAMAX_PREFIX/bin/ramax" \
    --hal-validator "$RAMAX_PREFIX/libexec/ramax/ramax_release_hal_check" \
    --version "$RAMAX_VERSION" --output "$EVIDENCE/full-program"
# This marker is reached only after four-format content and HAL validation.
printf 'PASS: glibc 2.17 runtime, dependency closure, launchers, full program\n' > "$EVIDENCE/PASS"
