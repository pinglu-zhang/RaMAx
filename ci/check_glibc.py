#!/usr/bin/env python3
"""Inspect required (not exported) GLIBC symbol versions in trusted ELF files."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess


def required_versions(text):
    needs = text.split("Version needs section", 1)
    if len(needs) == 1:
        return []
    return sorted(set(re.findall(r"Name: (GLIBC_\d+(?:\.\d+)+)\b", needs[1])),
                  key=lambda s: tuple(map(int, s[6:].split("."))))


def inspect(paths, maximum, closure):
    pending = [Path(p).resolve(strict=True) for p in paths]
    seen, records, errors = set(), [], []
    while pending:
        path = pending.pop()
        if path in seen:
            continue
        seen.add(path)
        version = subprocess.check_output(["readelf", "--version-info", str(path)],
                                         text=True, env={**os.environ, "LC_ALL": "C"})
        required = required_versions(version)
        bad = [v for v in required if tuple(map(int, v[6:].split("."))) > maximum]
        record = {"path": str(path), "required_glibc": required, "too_new": bad}
        if bad:
            errors.append(f"{path}: {', '.join(bad)}")
        if closure:
            result = subprocess.run(["ldd", str(path)], text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, env={**os.environ, "LC_ALL": "C"})
            record["ldd"] = result.stdout
            static = "statically linked" in result.stdout or "not a dynamic executable" in result.stdout
            if "not found" in result.stdout or (result.returncode and not static):
                errors.append(f"{path}: unresolved dependencies: {result.stdout.strip()}")
            for line in result.stdout.splitlines():
                match = re.match(r"\s*(?:\S+\s+=>\s+)?(/\S+)\s+\(", line)
                if match:
                    pending.append(Path(match[1]).resolve(strict=True))
        records.append(record)
    return {"status": "FAIL" if errors else "PASS", "maximum_glibc": ".".join(map(str, maximum)),
            "closure_checked": closure, "files": records, "errors": errors}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--maximum", default="2.17")
    p.add_argument("--closure", action="store_true",
                   help="Include resolved shared libraries; run in the target runtime")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("files", nargs="+")
    args = p.parse_args()
    report = inspect(args.files, tuple(map(int, args.maximum.split("."))), args.closure)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(report["status"], len(report["files"]), "ELF files;", args.output)
    if report["errors"]:
        raise SystemExit("\n".join(report["errors"]))


if __name__ == "__main__":
    main()
