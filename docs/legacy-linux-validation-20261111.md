# glibc 2.17 full-program acceptance — 2026-11-11

Status: PASS. Server evidence:
`/mnt/sda/zhangpinglu/exp/RaMAx/legacy-full-validation-20261111`.
The authoritative record is `acceptance.complete.json`, with per-run
`RESULT.json`, loader resolution, symbol reports and hashes.

The complete RaMAx 1.1.0 executable, cwd helper and HAL reader were built
with GCC 15.2 against the existing glibc 2.17 sysroot. Their ELF interpreter
points to ld-2.17.so. Runtime resolution confirms libc and pthread from this
sysroot. External tools were copied into this isolated validation directory
and their interpreter/RPATH fixed to the same runtime; installed tools were
not modified. This is an actual old libc test on the server's modern host
kernel, not a CentOS 7/kernel 3.10 or Slurm reproduction.

## Full-program checks

The existing deterministic release fixture contains three genomes, each with
a ~60 kb primary sequence and a 3,001 bp auxiliary sequence, including
substitutions and indels. Both runs use two threads.

| Route | Whole alignment/export wall time | Result |
|---|---:|---|
| Default Mash/wfmash routing | 0.920 s | MAF/HAL/PAF/GFA accepted |
| Native route, near-distance=0 | 0.321 s | MAF/HAL/PAF/GFA accepted |

These are tiny functional tests, not comparative performance benchmarks.
MAF rows were checked against input coordinates and DNA; PAF CIGAR consumption
and GFA complete input reconstruction were checked. The independent reader
called HAL validateAlignment and compared complete leaf DNA and contig counts.
The default run exported a 181,150-byte MAF and 161,169-byte HAL.

Neither full-program fixture invoked minipoa, so its production launcher was
also compiled and tested under the actual 2.17 runtime. Ten launcher
configurations passed, including temp-directory cleanup, concurrent calls,
memfd/tmpfs/file fallback, simulated tmpfs exhaustion, failure handling and
process-group timeout/cancellation. Eight executions of real minipoa on the
same synthetic indel input produced identical complete outputs.
Both the native-header and forced-legacy variants use real glibc 2.17 here.
Seven executable symbol audits passed. Runtime dependency resolution yielded
43 paths / 40 unique ELF files; their required GLIBC versions were <=2.17.

## Issues found and retained

- Nested sonLib/HAL Make builds did not inherit CMake's sysroot. This introduced
  __libc_single_threaded references and prevented old-libc linkage. CMake now
  forwards the sysroot and C/C++ flags to those builds.
- glibc 2.17 may report an exec failure through child exit 127 instead of
  returning ENOENT from posix_spawn. The regression now accepts both valid
  failure forms and still requires resource cleanup.
- Explicit command-line invocation of the old loader produced startup crashes
  for these binaries. Direct ELF execution with a fixed PT_INTERP passed.
  Failed entry-point experiments and earlier configuration/build logs remain.

No production binaries, original inputs or running genome tasks were replaced.
GitHub's full Conda/container release workflow has been configured but was not
run as part of this server acceptance.
