# Linux / glibc 2.17 support

The launchers support glibc 2.17 builds. This is separate from the compatibility
of a particular RaMAx binary, its external tools and shared libraries. A binary
built against a newer libc is not made portable by this source change.

## Child processes

All production posix_spawn calls use the same checked attribute implementation.
When compiled against glibc older than 2.24, it requests POSIX_SPAWN_USEVFORK.
When a process group is required, POSIX_SPAWN_SETPGROUP is combined with that flag.
Modern builds retain libc's default launch path. No application callback runs
inside vfork, and the parent never changes its working directory.

On glibc headers before 2.29, minipoa's cwd is set by the installed
`ramax-spawn-helper`. The helper executes chdir followed by execvp, preserving
the PID used to clean `minipoa_paths_<pid>.tmp`. Modern builds use the native
spawn chdir action. Both branches use the configured scratch directory even
when the program was launched elsewhere; signal/nonzero exits are cleaned too.
An unavailable helper is an explicit configuration error. Install it next to
`ramax`; embedding applications can set `RAMAX_SPAWN_HELPER` to its absolute
path. Copying only the ramax executable is not a complete legacy installation.

## Temporary input

Input uses memfd when supported by the running kernel and policy. Otherwise
the launcher tries an immediately unlinked file on /dev/shm, verifies it is
tmpfs, and passes the descriptor through exec. It writes through write(), not
mmap, so a full tmpfs can fall back without a SIGBUS. If both choices fail,
it creates a private mkostemp file in the scratch directory and removes it
after the invocation. File actions and attributes use scoped cleanup.

The summary retains `inputs(memfd/file)` and adds `tmpfs_inputs`. The three
counts represent input preparation attempts, excluding cache hits. A disk fallback is
reported once with its directory. A kernel or policy denial of memfd is cached.

For a cluster with shared work storage, select a **job-private local** directory
before starting RaMAx, for example:

```bash
export RAMAX_MINIPOA_TMPDIR="$SLURM_TMPDIR/ramax-minipoa-$SLURM_JOB_ID"
```

Only use this example when SLURM_TMPDIR is set to a valid node-local directory.
The override applies to minipoa's own temporary files as well as FASTA fallback;
it does not relocate final alignment outputs. Do not share this directory
between live RaMAx processes. It must remain available for the entire run.
Without the override, the existing work/minipoa_tmp location remains in use.
SIGKILL of the parent or loss of the node can still leave named temporary files;
normal completion and child failure cleanup do not imply crash-proof cleanup.

## Acceptance gates

1. BUILD_TESTING registers `ramax_external_spawn_test`: production modules with
   native and forced legacy compile paths, cwd/temp cleanup, memfd/tmpfs/file,
   redirection, attribute/action failures, concurrent calls and descendant
   cancellation. Forced legacy compilation is not an old runtime.
2. Conda tests inspect required GLIBC symbols for the main executable, helper,
   companion tool, minipoa, Mash, wfmash, samtools and the HAL acceptance reader.
   A requirement newer than 2.17 blocks acceptance.
3. The release workflow mounts a freshly installed package into a CentOS 7 /
   glibc 2.17 manylinux2014 container. It verifies the actual runtime, audits the
   resolved shared-library closure, builds/runs the launcher regressions there,
   and runs the complete synthetic MAF/HAL/PAF/GFA export and content validation.
   The HAL reader is built with the package's sysroot, not copied from Ubuntu.
   Image identity, exact package list and logs are retained. Publication depends
   on this job succeeding; workflow_dispatch performs validation only.

The internal reader can be installed with
`-DRAMAX_INSTALL_VALIDATION_TOOLS=ON` under libexec/ramax. It is enabled in
the Conda recipe. CI's manylinux image is documented by the
[manylinux project](https://github.com/pypa/manylinux#manylinux2014-centos-7-based-glibc-217).

Container testing uses the host kernel. It does not reproduce a partner's
kernel 3.10 backports, seccomp restrictions, Slurm limits, NUMA placement or
shared filesystem. Acceptance on their actual compute node is still required
before claiming that their entire environment is reproduced.

## Allocation profiling

No allocator, arena or malloc_trim setting is changed by this fix. If slowdown
persists after checking spawn mode and input counters, collect a bounded CPU
profile on the compute node and distinguish malloc/free, allocator locks,
malloc_trim and the actual alignment code. Record wall/user/system times,
peak RSS, minor faults and the loaded libc path. Only change the allocator after
a comparable workload establishes that allocation is a meaningful bottleneck.
