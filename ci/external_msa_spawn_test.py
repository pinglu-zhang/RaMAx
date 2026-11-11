#!/usr/bin/env python3
"""Small Linux launcher tests: native and forced legacy code paths.

Forced legacy compilation is not a glibc 2.17 runtime test. CI also runs this
suite inside manylinux2014 (actual glibc 2.17). Every run retains its own evidence.
"""
import argparse
import pathlib
import os
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--cc", default="gcc")
    parser.add_argument("--minipoa", type=pathlib.Path)
    args = parser.parse_args()
    root = pathlib.Path(__file__).resolve().parents[1]
    args.output.mkdir(parents=True, exist_ok=True)
    output = pathlib.Path(tempfile.mkdtemp(prefix="run-", dir=args.output.resolve()))
    print("Evidence:", output, flush=True)
    real_outputs = []
    env = os.environ.copy()
    if args.minipoa:
        env["RAMAX_TEST_MINIPOA"] = str(args.minipoa.resolve(strict=True))
    subprocess.run([args.cc, "-O2", str(root / "src/align/spawn_helper.c"),
                    "-o", str(output / "ramax-spawn-helper")], check=True, timeout=60)
    for mode in ("native", "legacy"):
        for suite in ("msa", "tool"):
            binary = output / (mode + "-" + suite)
            command = [args.cxx, "-std=c++20", "-O1", "-pthread",
                       "-I" + str(root / "include"),
                       "-I" + str(root / "third_party/spdlog/include")]
            if mode == "legacy":
                command += ["-DRAMAX_TEST_LEGACY_SPAWN"]
            command += [str(root / ("src/align/external_msa_runner.cpp" if suite == "msa"
                                    else "src/align/external_tool.cpp")),
                        str(root / ("ci/external_msa_spawn_test.cpp" if suite == "msa"
                                    else "ci/external_tool_spawn_test.cpp")),
                        "-o", str(binary)]
            symbols = ["posix_spawn", "posix_spawnattr_init",
                       "posix_spawnattr_setflags", "posix_spawnattr_destroy",
                       "posix_spawn_file_actions_adddup2"]
            if suite == "msa":
                symbols += ["posix_spawnp", "syscall", "mkostemp", "write"]
            for symbol in symbols:
                command.append("-Wl,--wrap=" + symbol)
            with (output / (binary.name + ".build.log")).open("w") as log:
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                               check=True, timeout=180)
            for backend in (("memfd", "shm", "file", "fullshm") if suite == "msa" else ("default",)):
                name = binary.name + "-" + backend
                cwd = output / (name + "-launch")
                cwd.mkdir()
                with (output / (name + ".test.log")).open("w") as log:
                    subprocess.run([str(binary), str(output / (name + "-scratch")), mode, backend],
                                   cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT,
                                   check=True, timeout=60)
                print(name, (output / (name + ".test.log")).read_text().splitlines()[-1],
                      flush=True)
                if suite == "msa" and args.minipoa:
                    lines = (output / (name + ".test.log")).read_text().splitlines()
                    values = [line for line in lines if line.startswith("REAL_OUTPUT ")]
                    assert len(values) == 1
                    real_outputs.append(values[0])
    if args.minipoa:
        assert len(real_outputs) == 8 and len(set(real_outputs)) == 1, "minipoa outputs differ"
        (output / "real-minipoa-output.txt").write_text(real_outputs[0] + "\n")
    (output / "PASS").write_text("native + legacy launch regressions passed\n")


if __name__ == "__main__":
    main()
