#!/usr/bin/env python3
"""Linux launcher regression; legacy headers are simulated, libc is native.

This does not replace real-window A/B testing under an actual old glibc loader.
Outputs use a new directory supplied by the caller and are retained on failure.
"""
import argparse
import pathlib
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--cxx", default="g++")
    args = parser.parse_args()
    root = pathlib.Path(__file__).resolve().parents[1]
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    source = (root / "src/align/external_msa_runner.cpp").read_text()
    # Inject only after all system headers: exercise our legacy compile branch
    # without changing declarations or pretending the host runtime is glibc 2.17.
    anchor = "extern char** environ;"
    assert source.count(anchor) == 1
    legacy = source.replace(anchor, anchor + "\n#undef __GLIBC_PREREQ\n"
        "#define __GLIBC_PREREQ(major, minor) ((major) < 2 || "
        "((major) == 2 && (minor) <= 17))\n")
    for mode, text in (("native", source), ("legacy", legacy)):
        unit = output / (mode + ".cpp")
        unit.write_text(text)
        binary = output / mode
        command = [args.cxx, "-std=c++20", "-O1", "-pthread",
                   "-I" + str(root / "include"),
                   "-I" + str(root / "third_party/spdlog/include"),
                   str(unit), str(root / "ci/external_msa_spawn_test.cpp"),
                   "-o", str(binary)]
        for symbol in ("posix_spawn", "posix_spawnp", "posix_spawnattr_init",
                       "posix_spawnattr_setflags", "posix_spawnattr_destroy"):
            command.append("-Wl,--wrap=" + symbol)
        with (output / (mode + ".build.log")).open("w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                           check=True, timeout=180)
        with (output / (mode + ".test.log")).open("w") as log:
            subprocess.run([str(binary), str(output / (mode + "-scratch")), mode],
                           stdout=log, stderr=subprocess.STDOUT,
                           check=True, timeout=60)
        print(mode, (output / (mode + ".test.log")).read_text().splitlines()[-1])


if __name__ == "__main__":
    main()
