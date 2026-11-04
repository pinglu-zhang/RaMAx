#!/usr/bin/env bash
set -euo pipefail
test "$(ramax --version)" = "RaMAx version ${PKG_VERSION}"
ramax --help >/dev/null
for executable in ramax ramax-paf-fasta mash wfmash samtools minipoa; do
  command -v "${executable}"
done
test "$(mash --version)" = '2.3'
test "$(wfmash --version 2>&1 | sed -n '1p')" = 'v0.14.0-0-g517e1bc'
test "$(samtools --version | sed -n '1p')" = 'samtools 1.23.1'
test "$(samtools --version | sed -n '2p')" = 'Using htslib 1.23.1'
test "$(minipoa -v)" = 'minipoa version 1.4.2'
for executable in ramax minipoa wfmash; do
  report="$(ldd "$(command -v "${executable}")")"
  printf '%s\n' "${report}"
  if printf '%s\n' "${report}" | grep -qE 'not found|/mnt/sda/|_build_env|_h_env'; then
    echo 'Unresolved or build-private runtime dependency' >&2
    exit 1
  fi
done
test -s "${PREFIX}/share/ramax/packaged-tools.tsv"
test -s "${PREFIX}/share/licenses/ramax/HAL-LICENSE.txt"
