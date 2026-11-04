#!/usr/bin/env python3
"""Run an isolated, resource-guarded RaMAx HAL/MAF A/B comparison."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import pathlib
import shutil
import signal
import socket
import subprocess
import sys
import time
from dataclasses import dataclass
from typing import Any

GIB = 1024**3
KIB = 1024
REPORT_SCHEMA = 1


class RunnerError(RuntimeError):
    pass


@dataclass(frozen=True)
class ProcessStat:
    pid: int
    ppid: int
    pgid: int
    start_ticks: int
    name: str
    rss_kib: int
    hwm_kib: int
    vm_kib: int


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stat_manifest(path: pathlib.Path, include_sha256: bool) -> dict[str, Any]:
    info = path.stat()
    result: dict[str, Any] = {
        "path": str(path),
        "size_bytes": info.st_size,
        "mtime_ns": info.st_mtime_ns,
    }
    if include_sha256:
        result["sha256"] = sha256_file(path)
    return result


def maf_manifest(path: pathlib.Path) -> dict[str, Any]:
    """Preserve every block byte and duplicate, ignoring only block order."""
    manifest = stat_manifest(path, False)
    raw_digest = hashlib.sha256()
    preamble_digest = hashlib.sha256()
    block_digest = None
    block_digests: list[bytes] = []
    with path.open("rb") as stream:
        for line in stream:
            raw_digest.update(line)
            if not line.strip():
                if block_digest is not None:
                    block_digests.append(block_digest.digest())
                    block_digest = None
                continue
            if line.startswith(b"a ") or line.rstrip(b"\r\n") == b"a":
                if block_digest is not None:
                    block_digests.append(block_digest.digest())
                block_digest = hashlib.sha256()
            if block_digest is None:
                preamble_digest.update(line)
            else:
                block_digest.update(line)
    if block_digest is not None:
        block_digests.append(block_digest.digest())
    block_digests.sort()
    multiset_digest = hashlib.sha256()
    for digest in block_digests:
        multiset_digest.update(digest)
    manifest.update(
        sha256=raw_digest.hexdigest(),
        preamble_sha256=preamble_digest.hexdigest(),
        block_multiset_sha256=multiset_digest.hexdigest(),
        block_count=len(block_digests),
    )
    return manifest


def parse_seqfile(path: pathlib.Path, hash_inputs: bool) -> dict[str, Any]:
    mappings: list[dict[str, Any]] = []
    tree: str | None = None
    with path.open("rt", encoding="utf-8") as stream:
        for raw_line in stream:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if tree is None and line.startswith("("):
                tree = line
                continue
            fields = line.split(None, 1)
            if len(fields) != 2:
                raise RunnerError(f"invalid seqfile mapping: {line!r}")
            species, location = fields
            entry: dict[str, Any] = {"species": species, "location": location}
            if "://" in location:
                entry["kind"] = "url"
                entry["note"] = "not fetched or hashed by benchmark runner"
            else:
                source = pathlib.Path(location).expanduser()
                if not source.is_absolute():
                    raise RunnerError(
                        "benchmark seqfile local paths must be absolute: "
                        f"{location}"
                    )
                source = source.resolve()
                if not source.is_file():
                    raise RunnerError(f"seqfile input is not a regular file: {source}")
                entry["kind"] = "local"
                entry["stat"] = stat_manifest(source, hash_inputs)
            mappings.append(entry)
    if tree is None:
        raise RunnerError("HAL benchmark seqfile has no Newick tree record")
    if not mappings:
        raise RunnerError("seqfile has no genome mappings")
    species = [entry["species"] for entry in mappings]
    if len(set(species)) != len(species):
        raise RunnerError("seqfile contains duplicate genome names")
    return {
        "seqfile": stat_manifest(path, True),
        "tree_record": tree,
        "mappings": mappings,
        "large_input_hashes_enabled": hash_inputs,
    }


def read_mem_available_kib() -> int:
    with open("/proc/meminfo", "rt", encoding="ascii") as stream:
        for line in stream:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1])
    raise RunnerError("/proc/meminfo does not expose MemAvailable")


def read_process_stat(pid: int) -> ProcessStat | None:
    proc = pathlib.Path("/proc") / str(pid)
    try:
        raw = (proc / "stat").read_text(encoding="ascii")
        right_paren = raw.rfind(")")
        if right_paren < 0:
            return None
        name = raw[raw.find("(") + 1 : right_paren]
        fields = raw[right_paren + 2 :].split()
        ppid = int(fields[1])
        pgid = int(fields[2])
        start_ticks = int(fields[19])
        status: dict[str, int] = {}
        with (proc / "status").open("rt", encoding="ascii") as stream:
            for line in stream:
                key, _, value = line.partition(":")
                if key in {"VmRSS", "VmHWM", "VmSize"}:
                    parts = value.split()
                    status[key] = int(parts[0]) if parts else 0
        return ProcessStat(
            pid=pid,
            ppid=ppid,
            pgid=pgid,
            start_ticks=start_ticks,
            name=name,
            rss_kib=status.get("VmRSS", 0),
            hwm_kib=status.get("VmHWM", 0),
            vm_kib=status.get("VmSize", 0),
        )
    except (FileNotFoundError, ProcessLookupError, PermissionError, ValueError):
        return None


def process_group_stats(pgid: int) -> list[ProcessStat]:
    result: list[ProcessStat] = []
    for entry in pathlib.Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        stat = read_process_stat(int(entry.name))
        if stat is not None and stat.pgid == pgid:
            result.append(stat)
    result.sort(key=lambda item: item.pid)
    return result

def process_tree_stats(root_pid: int) -> list[ProcessStat]:
    all_stats: list[ProcessStat] = []
    for entry in pathlib.Path("/proc").iterdir():
        if entry.name.isdigit():
            stat = read_process_stat(int(entry.name))
            if stat is not None:
                all_stats.append(stat)
    children: dict[int, list[ProcessStat]] = {}
    by_pid = {item.pid: item for item in all_stats}
    for item in all_stats:
        children.setdefault(item.ppid, []).append(item)
    result: list[ProcessStat] = []
    pending = [root_pid]
    seen: set[int] = set()
    while pending:
        pid = pending.pop()
        if pid in seen:
            continue
        seen.add(pid)
        item = by_pid.get(pid)
        if item is not None:
            result.append(item)
        pending.extend(child.pid for child in children.get(pid, []))
    result.sort(key=lambda item: item.pid)
    return result


def require_owned_group(
    process: subprocess.Popen[Any], pgid: int, leader_start_ticks: int
) -> None:
    if pgid != process.pid or pgid == os.getpgrp():
        raise RunnerError("refusing to signal a process group not created by this runner")
    leader = read_process_stat(process.pid)
    if leader is not None:
        if leader.pgid != pgid or leader.start_ticks != leader_start_ticks:
            raise RunnerError("refusing to signal a reused or foreign process group leader")
        return
    if not process_group_stats(pgid):
        return
    if process.poll() is None:
        raise RunnerError("cannot prove ownership of the live process-group leader")


def process_group_exists(pgid: int) -> bool:
    return bool(process_group_stats(pgid))


def stop_owned_group(
    process: subprocess.Popen[Any],
    pgid: int,
    leader_start_ticks: int,
    grace_seconds: float,
) -> dict[str, Any]:
    require_owned_group(process, pgid, leader_start_ticks)
    actions: list[str] = []
    if process_group_exists(pgid):
        try:
            os.killpg(pgid, signal.SIGTERM)
            actions.append("SIGTERM")
        except ProcessLookupError:
            pass
    deadline = time.monotonic() + grace_seconds
    while time.monotonic() < deadline:
        process.poll()
        if not process_group_exists(pgid):
            break
        time.sleep(0.2)
    if process_group_exists(pgid):
        require_owned_group(process, pgid, leader_start_ticks)
        try:
            os.killpg(pgid, signal.SIGKILL)
            actions.append("SIGKILL")
        except ProcessLookupError:
            pass
    try:
        process.wait(timeout=max(1.0, grace_seconds))
    except subprocess.TimeoutExpired as error:
        raise RunnerError("process-group leader did not exit after SIGKILL") from error
    return {"signals_sent": actions, "returncode": process.returncode}


def resolve_executable(value: str, label: str) -> pathlib.Path:
    candidate = shutil.which(value)
    path = pathlib.Path(candidate if candidate else value).expanduser().resolve()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise RunnerError(f"{label} is not executable: {value}")
    return path


def run_capture(command: list[str], label: str) -> str:
    completed = subprocess.run(
        command,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise RunnerError(f"{label} failed with status {completed.returncode}: {detail}")
    return completed.stdout.replace("\r\n", "\n").rstrip() + "\n"


def run_validate(tool: pathlib.Path, hal_path: pathlib.Path, log_path: pathlib.Path) -> None:
    command = [str(tool), str(hal_path)]
    completed = subprocess.run(
        command,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    log_path.write_bytes(completed.stdout)
    if completed.returncode != 0:
        raise RunnerError(
            f"halValidate failed with status {completed.returncode}; see {log_path}"
        )


def hal2fasta_digest(
    tool: pathlib.Path,
    hal_path: pathlib.Path,
    genome: str,
    stderr_path: pathlib.Path,
    grace_seconds: float,
) -> dict[str, Any]:
    command = [str(tool), str(hal_path), genome]
    with stderr_path.open("wb") as errors:
        process = subprocess.Popen(
            command,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=errors,
            start_new_session=True,
        )
        leader = read_process_stat(process.pid)
        if leader is None or leader.pgid != process.pid:
            process.kill()
            process.wait()
            raise RunnerError("could not establish an owned hal2fasta process group")
        assert process.stdout is not None
        records: list[dict[str, Any]] = []
        current_name: str | None = None
        current_digest = hashlib.sha256()
        current_length = 0

        def finish_record() -> None:
            nonlocal current_name, current_digest, current_length
            if current_name is None:
                return
            records.append(
                {
                    "name": current_name,
                    "length": current_length,
                    "dna_sha256": current_digest.hexdigest(),
                }
            )

        try:
            for raw_line in process.stdout:
                line = raw_line.rstrip(b"\r\n")
                if line.startswith(b">"):
                    finish_record()
                    current_name = line[1:].split(None, 1)[0].decode(
                        "utf-8", errors="strict"
                    )
                    current_digest = hashlib.sha256()
                    current_length = 0
                elif line:
                    if current_name is None:
                        raise RunnerError("hal2fasta emitted DNA before a FASTA header")
                    current_digest.update(line)
                    current_length += len(line)
            finish_record()
            returncode = process.wait()
        except BaseException:
            stop_owned_group(process, process.pid, leader.start_ticks, grace_seconds)
            raise
    if returncode != 0:
        raise RunnerError(
            f"hal2fasta failed for {genome} with status {returncode}; see {stderr_path}"
        )
    aggregate = hashlib.sha256()
    for record in records:
        aggregate.update(record["name"].encode("utf-8"))
        aggregate.update(b"\0")
        aggregate.update(str(record["length"]).encode("ascii"))
        aggregate.update(b"\0")
        aggregate.update(record["dna_sha256"].encode("ascii"))
        aggregate.update(b"\n")
    return {
        "genome": genome,
        "records": records,
        "total_bases": sum(record["length"] for record in records),
        "digest_sha256": aggregate.hexdigest(),
    }


def collect_hal_semantics(
    hal_path: pathlib.Path,
    leaf_names: list[str],
    hal_validate: pathlib.Path,
    hal_stats: pathlib.Path,
    hal2fasta: pathlib.Path,
    logs: pathlib.Path,
    grace_seconds: float,
) -> dict[str, Any]:
    run_validate(hal_validate, hal_path, logs / "halValidate.log")
    basic = run_capture([str(hal_stats), str(hal_path)], "halStats basic")
    tree = run_capture([str(hal_stats), "--tree", str(hal_path)], "halStats tree")
    genomes_text = run_capture(
        [str(hal_stats), "--genomes", str(hal_path)], "halStats genomes"
    )
    genomes = genomes_text.split()
    if len(set(genomes)) != len(genomes):
        raise RunnerError("halStats --genomes returned duplicate genome names")
    missing_leaves = sorted(set(leaf_names) - set(genomes))
    if missing_leaves:
        raise RunnerError(f"HAL is missing seqfile leaves: {missing_leaves}")

    metadata: dict[str, Any] = {}
    for genome in sorted(genomes):
        metadata[genome] = {
            "sequence_stats": run_capture(
                [str(hal_stats), "--sequenceStats", genome, str(hal_path)],
                f"halStats sequenceStats {genome}",
            ),
            "num_segments": run_capture(
                [str(hal_stats), "--numSegments", genome, str(hal_path)],
                f"halStats numSegments {genome}",
            ),
        }

    leaf_dna: dict[str, Any] = {}
    for leaf in sorted(leaf_names):
        leaf_dna[leaf] = hal2fasta_digest(
            hal2fasta,
            hal_path,
            leaf,
            logs / f"hal2fasta.{leaf}.stderr.log",
            grace_seconds,
        )
    return {
        "validated": True,
        "basic": basic,
        "tree": tree,
        "genomes": genomes,
        "genome_metadata": metadata,
        "leaf_dna": leaf_dna,
    }


def output_semantic_view(run: dict[str, Any]) -> dict[str, Any]:
    semantics = run["outputs"]["hal_semantics"]
    return {
        "tree": semantics["tree"],
        "genomes": semantics["genomes"],
        "genome_metadata": semantics["genome_metadata"],
        "leaf_dna": semantics["leaf_dna"],
    }


def write_json_atomic(path: pathlib.Path, value: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wt", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
    os.replace(temporary, path)


def split_returncode(returncode: int | None) -> tuple[int | None, int | None]:
    if returncode is None:
        return None, None
    if returncode < 0:
        return None, -returncode
    return returncode, None


def run_one(
    label: str,
    binary: pathlib.Path,
    args: argparse.Namespace,
    root_dir: pathlib.Path,
    input_manifest: dict[str, Any],
    tools: dict[str, pathlib.Path],
) -> dict[str, Any]:
    current_manifest = parse_seqfile(args.seqfile, args.hash_input_files)
    if current_manifest != input_manifest:
        raise RunnerError(
            "seqfile or a local input changed after benchmark preflight"
        )
    run_dir = root_dir / label
    work_dir = run_dir / "work"
    output_dir = run_dir / "output"
    logs_dir = run_dir / "logs"
    work_dir.mkdir(parents=True)
    output_dir.mkdir()
    logs_dir.mkdir()
    hal_path = output_dir / "alignment.hal"
    maf_path = output_dir / "alignment.maf"
    command = [
        str(binary),
        "--input",
        str(args.seqfile),
        "--output",
        str(maf_path),
        "--output",
        str(hal_path),
        "--workdir",
        str(work_dir),
        "--threads",
        str(args.threads),
        "--root",
        args.root,
        *args.ramax_args,
    ]
    result: dict[str, Any] = {
        "label": label,
        "binary": stat_manifest(binary, True),
        "command": command,
        "work_directory": str(work_dir),
        "started_at": utc_now(),
        "status": "running",
    }
    aggregate_path = logs_dir / "memory_samples.tsv"
    process_path = logs_dir / "process_samples.tsv"
    console_path = logs_dir / "console.log"
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(args.threads)
    start = time.monotonic()
    safety_trigger: dict[str, Any] | None = None
    stop_record: dict[str, Any] | None = None
    peak_tree_rss = 0
    peak_tree_hwm = 0
    min_available = read_mem_available_kib()
    min_free = shutil.disk_usage(root_dir).free

    with console_path.open("wb") as console, aggregate_path.open(
        "wt", encoding="utf-8", newline=""
    ) as aggregate_stream, process_path.open(
        "wt", encoding="utf-8", newline=""
    ) as process_stream:
        aggregate_writer = csv.writer(aggregate_stream, delimiter="\t")
        process_writer = csv.writer(process_stream, delimiter="\t")
        aggregate_writer.writerow(
            [
                "timestamp_utc",
                "elapsed_seconds",
                "pgid",
                "process_count",
                "tree_rss_kib",
                "tree_hwm_kib",
                "mem_available_kib",
                "output_free_bytes",
            ]
        )
        process_writer.writerow(
            [
                "timestamp_utc",
                "elapsed_seconds",
                "pid",
                "ppid",
                "pgid",
                "name",
                "rss_kib",
                "hwm_kib",
                "vm_kib",
            ]
        )
        process = subprocess.Popen(
            command,
            cwd=run_dir,
            env=environment,
            stdin=subprocess.DEVNULL,
            stdout=console,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        leader = read_process_stat(process.pid)
        if leader is None or leader.pgid != process.pid:
            process.kill()
            process.wait()
            raise RunnerError("could not establish an owned RaMAx process group")
        pgid = process.pid
        result["pid"] = process.pid
        result["pgid"] = pgid
        result["leader_start_ticks"] = leader.start_ticks
        try:
            while process.poll() is None:
                stamp = utc_now()
                elapsed = time.monotonic() - start
                members = process_tree_stats(process.pid)
                tree_rss = sum(item.rss_kib for item in members)
                tree_hwm = sum(item.hwm_kib for item in members)
                available = read_mem_available_kib()
                free_bytes = shutil.disk_usage(root_dir).free
                peak_tree_rss = max(peak_tree_rss, tree_rss)
                peak_tree_hwm = max(peak_tree_hwm, tree_hwm)
                min_available = min(min_available, available)
                min_free = min(min_free, free_bytes)
                aggregate_writer.writerow(
                    [
                        stamp,
                        f"{elapsed:.6f}",
                        pgid,
                        len(members),
                        tree_rss,
                        tree_hwm,
                        available,
                        free_bytes,
                    ]
                )
                for item in members:
                    process_writer.writerow(
                        [
                            stamp,
                            f"{elapsed:.6f}",
                            item.pid,
                            item.ppid,
                            item.pgid,
                            item.name,
                            item.rss_kib,
                            item.hwm_kib,
                            item.vm_kib,
                        ]
                    )
                aggregate_stream.flush()
                process_stream.flush()
                if tree_rss * KIB > args.max_rss_gib * GIB:
                    safety_trigger = {
                        "kind": "max_rss_gib",
                        "limit_gib": args.max_rss_gib,
                        "observed_tree_rss_kib": tree_rss,
                        "timestamp_utc": stamp,
                    }
                elif available * KIB < args.min_available_gib * GIB:
                    safety_trigger = {
                        "kind": "min_available_gib",
                        "limit_gib": args.min_available_gib,
                        "observed_mem_available_kib": available,
                        "timestamp_utc": stamp,
                    }
                elif free_bytes < args.min_free_gib * GIB:
                    safety_trigger = {
                        "kind": "min_free_gib",
                        "limit_gib": args.min_free_gib,
                        "observed_output_free_bytes": free_bytes,
                        "timestamp_utc": stamp,
                    }
                if safety_trigger is not None:
                    stop_record = stop_owned_group(
                        process, pgid, leader.start_ticks, args.termination_grace_seconds
                    )
                    break
                time.sleep(args.sample_seconds)
            if process.poll() is None:
                process.wait()
        except BaseException:
            if process.poll() is None or process_group_exists(pgid):
                stop_record = stop_owned_group(
                    process, pgid, leader.start_ticks, args.termination_grace_seconds
                )
            raise

    elapsed = time.monotonic() - start
    exit_code, terminating_signal = split_returncode(process.returncode)
    result.update(
        {
            "finished_at": utc_now(),
            "elapsed_seconds": elapsed,
            "returncode": process.returncode,
            "exit_code": exit_code,
            "terminating_signal": terminating_signal,
            "safety_trigger": safety_trigger,
            "stop_record": stop_record,
            "resource_summary": {
                "peak_process_tree_rss_kib": peak_tree_rss,
                "peak_sum_process_vmhwm_kib": peak_tree_hwm,
                "minimum_mem_available_kib": min_available,
                "minimum_output_free_bytes": min_free,
                "aggregate_samples_tsv": str(aggregate_path),
                "per_process_samples_tsv": str(process_path),
            },
        }
    )
    if safety_trigger is not None:
        result["status"] = "safety_stopped"
        return result
    if process.returncode != 0:
        result["status"] = "failed"
        return result
    if not hal_path.is_file() or hal_path.stat().st_size == 0:
        result["status"] = "failed"
        result["failure"] = "RaMAx exited zero without a nonempty HAL output"
        return result
    if not maf_path.is_file() or maf_path.stat().st_size == 0:
        result["status"] = "failed"
        result["failure"] = "RaMAx exited zero without a nonempty MAF output"
        return result

    result["outputs"] = {
        "hal": stat_manifest(hal_path, True),
        "maf": maf_manifest(maf_path),
    }
    try:
        semantics = collect_hal_semantics(
            hal_path,
            args.expected_leaf,
            tools["hal_validate"],
            tools["hal_stats"],
            tools["hal2fasta"],
            logs_dir,
            args.termination_grace_seconds,
        )
        result["outputs"]["hal_semantics"] = semantics
        result["status"] = "passed"
    except Exception as error:
        result["status"] = "validation_failed"
        result["failure"] = str(error)
    return result


def compare_runs(
    baseline: dict[str, Any],
    candidate: dict[str, Any],
    h5diff: pathlib.Path,
    comparison_log: pathlib.Path,
) -> dict[str, Any]:
    if baseline.get("status") != "passed" or candidate.get("status") != "passed":
        return {
            "performed": False,
            "reason": "both runs must pass execution and semantic validation",
        }
    baseline_hal = baseline["outputs"]["hal"]
    candidate_hal = candidate["outputs"]["hal"]
    baseline_maf = baseline["outputs"]["maf"]
    candidate_maf = candidate["outputs"]["maf"]
    with comparison_log.open("wb") as log:
        payload_comparison = subprocess.run(
            [str(h5diff), "-q", baseline_hal["path"], candidate_hal["path"]],
            stdin=subprocess.DEVNULL,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if payload_comparison.returncode not in (0, 1):
        raise RunnerError(
            f"h5diff failed with status {payload_comparison.returncode}; "
            f"see {comparison_log}"
        )
    payload_equal = payload_comparison.returncode == 0
    semantic_equal = payload_equal and (
        output_semantic_view(baseline) == output_semantic_view(candidate)
    )
    hal_bytes_equal = baseline_hal["sha256"] == candidate_hal["sha256"]
    maf_equal = baseline_maf["sha256"] == candidate_maf["sha256"]
    maf_blocks_equal = all(
        baseline_maf[key] == candidate_maf[key]
        for key in ("preamble_sha256", "block_multiset_sha256", "block_count")
    )
    if hal_bytes_equal:
        hal_difference = "identical_bytes"
    elif semantic_equal:
        hal_difference = "different_binary_layout_same_semantics"
    else:
        hal_difference = "semantic_mismatch"
    return {
        "performed": True,
        "recorded_input_manifest_rechecked": True,
        "hal_binary_sha256_equal": hal_bytes_equal,
        "hal_semantically_equal": semantic_equal,
        "hdf5_payload_equal": payload_equal,
        "h5diff_log": str(comparison_log),
        "hal_difference_classification": hal_difference,
        "maf_sha256_equal": maf_equal,
        "maf_block_multiset_equal": maf_blocks_equal,
        "passed": semantic_equal and maf_blocks_equal,
    }


def write_summary_tsv(path: pathlib.Path, report: dict[str, Any]) -> None:
    fields = [
        "label",
        "status",
        "binary_sha256",
        "elapsed_seconds",
        "exit_code",
        "terminating_signal",
        "safety_trigger",
        "peak_process_tree_rss_kib",
        "peak_sum_process_vmhwm_kib",
        "minimum_mem_available_kib",
        "minimum_output_free_bytes",
        "hal_validated",
        "hal_sha256",
        "maf_sha256",
        "hal_semantically_equal",
        "hdf5_payload_equal",
        "hal_difference_classification",
        "maf_sha256_equal",
        "maf_block_multiset_equal",
        "comparison_passed",
    ]
    with path.open("wt", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        comparison = report.get("comparison", {})
        for run in report.get("runs", []):
            resources = run.get("resource_summary", {})
            outputs = run.get("outputs", {})
            semantics = outputs.get("hal_semantics", {})
            writer.writerow(
                {
                    "label": run.get("label"),
                    "status": run.get("status"),
                    "binary_sha256": run.get("binary", {}).get("sha256"),
                    "elapsed_seconds": run.get("elapsed_seconds"),
                    "exit_code": run.get("exit_code"),
                    "terminating_signal": run.get("terminating_signal"),
                    "safety_trigger": json.dumps(
                        run.get("safety_trigger"), sort_keys=True
                    ),
                    "peak_process_tree_rss_kib": resources.get(
                        "peak_process_tree_rss_kib"
                    ),
                    "peak_sum_process_vmhwm_kib": resources.get(
                        "peak_sum_process_vmhwm_kib"
                    ),
                    "minimum_mem_available_kib": resources.get(
                        "minimum_mem_available_kib"
                    ),
                    "minimum_output_free_bytes": resources.get(
                        "minimum_output_free_bytes"
                    ),
                    "hal_validated": semantics.get("validated", False),
                    "hal_sha256": outputs.get("hal", {}).get("sha256"),
                    "maf_sha256": outputs.get("maf", {}).get("sha256"),
                    "hal_semantically_equal": comparison.get(
                        "hal_semantically_equal"
                    ),
                    "hdf5_payload_equal": comparison.get("hdf5_payload_equal"),
                    "hal_difference_classification": comparison.get(
                        "hal_difference_classification"
                    ),
                    "maf_sha256_equal": comparison.get("maf_sha256_equal"),
                    "maf_block_multiset_equal": comparison.get(
                        "maf_block_multiset_equal"
                    ),
                    "comparison_passed": comparison.get("passed"),
                }
            )


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run baseline then candidate in fresh directories, enforce memory/disk "
            "guards, and compare all HDF5 payloads plus exact MAF block multisets."
        )
    )
    parser.add_argument("--baseline", help="baseline RaMAx executable")
    parser.add_argument("--candidate", required=True, help="candidate RaMAx executable")
    parser.add_argument("--seqfile", required=True, type=pathlib.Path)
    parser.add_argument("--root", required=True, help="HAL root genome name")
    parser.add_argument(
        "--expected-leaf", action="append",
        help="expected leaf below --root; repeat for subtree runs using a larger seqfile",
    )
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    parser.add_argument(
        "--candidate-only",
        action="store_true",
        help="run only candidate; no semantic A/B claim is made",
    )
    parser.add_argument("--sample-seconds", type=float, default=1.0)
    parser.add_argument("--termination-grace-seconds", type=float, default=30.0)
    parser.add_argument("--min-available-gib", type=float, required=True)
    parser.add_argument("--max-rss-gib", type=float, required=True)
    parser.add_argument("--min-free-gib", type=float, required=True)
    parser.add_argument(
        "--hash-input-files",
        action="store_true",
        help="SHA-256 every local FASTA; default records size/mtime only",
    )
    parser.add_argument("--hal-validate", default="halValidate")
    parser.add_argument("--hal-stats", default="halStats")
    parser.add_argument("--hal2fasta", default="hal2fasta")
    parser.add_argument("--h5diff", default="h5diff")
    parser.add_argument(
        "ramax_args",
        nargs=argparse.REMAINDER,
        help="additional RaMAx arguments after --",
    )
    args = parser.parse_args(argv)
    if args.ramax_args and args.ramax_args[0] == "--":
        args.ramax_args = args.ramax_args[1:]
    protected = {
        "-i", "--input", "-o", "--output", "-w", "--workdir",
        "-t", "--threads", "--root", "--restart",
    }
    for token in args.ramax_args:
        option = token.split("=", 1)[0]
        if option in protected or (
            not option.startswith("--") and option[:2] in protected
        ):
            parser.error(
                f"{option} cannot be overridden after --; use the runner option"
            )
    if args.threads <= 0:
        parser.error("--threads must be positive")
    for name in (
        "sample_seconds",
        "termination_grace_seconds",
        "min_available_gib",
        "max_rss_gib",
        "min_free_gib",
    ):
        value = getattr(args, name)
        if not math.isfinite(value) or value <= 0:
            parser.error(f"--{name.replace('_', '-')} must be finite and positive")
    if args.candidate_only:
        if args.baseline:
            parser.error("--baseline and --candidate-only are mutually exclusive")
    elif not args.baseline:
        parser.error("--baseline is required unless --candidate-only is used")
    return args


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    if os.name != "posix" or not pathlib.Path("/proc/self/stat").is_file():
        raise RunnerError("this resource sampler requires Linux /proc")
    args.seqfile = args.seqfile.expanduser().resolve()
    args.output_dir = args.output_dir.expanduser().resolve()
    if not args.seqfile.is_file():
        raise RunnerError(f"seqfile is not a regular file: {args.seqfile}")
    if args.output_dir.exists():
        raise RunnerError(f"refusing to reuse output directory: {args.output_dir}")
    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    initial_free = shutil.disk_usage(args.output_dir.parent).free
    if initial_free < args.min_free_gib * GIB:
        raise RunnerError(
            f"output filesystem has {initial_free / GIB:.3f} GiB free, below guard"
        )
    initial_available = read_mem_available_kib()
    if initial_available * KIB < args.min_available_gib * GIB:
        raise RunnerError(
            f"MemAvailable is {initial_available * KIB / GIB:.3f} GiB, below guard"
        )
    candidate = resolve_executable(args.candidate, "candidate")
    baseline = (
        resolve_executable(args.baseline, "baseline") if args.baseline else None
    )
    tools = {
        "hal_validate": resolve_executable(args.hal_validate, "halValidate"),
        "hal_stats": resolve_executable(args.hal_stats, "halStats"),
        "hal2fasta": resolve_executable(args.hal2fasta, "hal2fasta"),
        "h5diff": resolve_executable(args.h5diff, "h5diff"),
    }
    input_manifest = parse_seqfile(args.seqfile, args.hash_input_files)
    mapped_species = {entry["species"] for entry in input_manifest["mappings"]}
    if args.expected_leaf is None:
        args.expected_leaf = sorted(mapped_species)
    if len(args.expected_leaf) != len(set(args.expected_leaf)):
        raise RunnerError("--expected-leaf contains duplicate names")
    if not set(args.expected_leaf).issubset(mapped_species):
        raise RunnerError("--expected-leaf contains names absent from the seqfile")
    args.output_dir.mkdir()
    report: dict[str, Any] = {
        "schema_version": REPORT_SCHEMA,
        "mode": "candidate-only" if args.candidate_only else "ab",
        "started_at": utc_now(),
        "host": socket.gethostname(),
        "configuration": {
            "root": args.root,
            "expected_leaves": args.expected_leaf,
            "threads": args.threads,
            "sample_seconds": args.sample_seconds,
            "termination_grace_seconds": args.termination_grace_seconds,
            "min_available_gib": args.min_available_gib,
            "max_rss_gib": args.max_rss_gib,
            "min_free_gib": args.min_free_gib,
            "extra_ramax_args": args.ramax_args,
            "tools": {
                name: stat_manifest(path, True) for name, path in tools.items()
            },
        },
        "input_manifest": input_manifest,
        "runs": [],
    }

    success = False
    try:
        if baseline is not None:
            baseline_run = run_one(
                "baseline", baseline, args, args.output_dir, input_manifest, tools
            )
            report["runs"].append(baseline_run)
            if baseline_run["status"] != "passed":
                report["comparison"] = {
                    "performed": False,
                    "reason": "baseline did not pass; candidate was not started",
                }
            else:
                candidate_run = run_one(
                    "candidate", candidate, args, args.output_dir, input_manifest, tools
                )
                report["runs"].append(candidate_run)
                report["comparison"] = compare_runs(
                    baseline_run, candidate_run, tools["h5diff"],
                    args.output_dir / "h5diff.log",
                )
                success = bool(report["comparison"].get("passed"))
        else:
            candidate_run = run_one(
                "candidate", candidate, args, args.output_dir, input_manifest, tools
            )
            report["runs"].append(candidate_run)
            report["comparison"] = {
                "performed": False,
                "reason": "candidate-only mode intentionally has no semantic control",
            }
            success = candidate_run["status"] == "passed"
    except BaseException as error:
        report["runner_error"] = f"{type(error).__name__}: {error}"
        if isinstance(error, KeyboardInterrupt):
            report["interrupted"] = True
    finally:
        report["finished_at"] = utc_now()
        report["passed"] = success
        write_json_atomic(args.output_dir / "report.json", report)
        write_summary_tsv(args.output_dir / "summary.tsv", report)
    return 0 if success else 1


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1:]))
    except RunnerError as error:
        print(f"run_hal_memory_ab.py: {error}", file=sys.stderr)
        sys.exit(2)
