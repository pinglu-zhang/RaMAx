from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import math
import os
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

from . import SCHEMA_VERSION, __version__


class EvalError(RuntimeError):
    """A user-facing, evidence-preserving evaluation error."""


def now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")


utc_now = now_iso


def canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def stable_hash(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path, include_hash: bool = False) -> dict[str, Any]:
    resolved = path.expanduser().resolve(strict=True)
    stat = resolved.stat()
    record: dict[str, Any] = {"path": str(resolved), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if include_hash:
        record["sha256"] = sha256_file(resolved)
    return record


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def _json_safe(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    return value


def atomic_json(path: Path, value: Any) -> None:
    atomic_text(path, json.dumps(_json_safe(value), ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def format_scalar(value: Any) -> Any:
    if isinstance(value, float):
        return "NA" if not math.isfinite(value) else repr(value)
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return ""
    return value


def atomic_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fieldnames)
        for row in rows:
            if isinstance(row, Mapping):
                writer.writerow([format_scalar(row.get(key, "")) for key in fieldnames])
            else:
                writer.writerow([format_scalar(value) for value in row])
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def wilson(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        return float("nan"), float("nan")
    p = successes / total
    denominator = 1.0 + z * z / total
    centre = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt((p * (1 - p) + z * z / (4 * total)) / total) / denominator
    return max(0.0, centre - half), min(1.0, centre + half)


def safe_div(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator else float("nan")


def resolve_tool(name: str, explicit: Optional[Any] = None, tool_dir: Optional[Any] = None) -> Path:
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(Path(explicit).expanduser())
    if tool_dir is not None:
        candidates.append(Path(tool_dir).expanduser() / name)
    found = shutil.which(name)
    if found:
        candidates.append(Path(found))
    for candidate in candidates:
        try:
            resolved = candidate.resolve(strict=True)
        except FileNotFoundError:
            continue
        if resolved.is_file() and os.access(resolved, os.X_OK):
            return resolved
    raise EvalError(f"missing executable {name}; use --tool-dir or an explicit tool option")


def tool_record(path: Any, version_args: Sequence[str] = ("--version",)) -> dict[str, Any]:
    resolved = path if isinstance(path, Path) else resolve_tool(str(path))
    version = "unknown"
    try:
        proc = subprocess.run([str(resolved), *version_args], capture_output=True, text=True, timeout=30, check=False)
        text = (proc.stdout or proc.stderr).strip()
        if text:
            version = text.splitlines()[0]
    except (OSError, subprocess.SubprocessError):
        pass
    return {"path": str(resolved), "sha256": sha256_file(resolved), "version": version}


@dataclass
class CommandResult:
    args: Sequence[str]
    returncode: int
    stdout: str
    stderr: str
    stdout_path: Optional[str] = None
    stderr_path: Optional[str] = None

    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)


def run_command(
    argv: Sequence[str], *, input_text: Optional[str] = None, log_path: Optional[Path] = None,
    log_dir: Optional[Path] = None, label: str = "command", stdout_path: Optional[Path] = None,
    cwd: Optional[Path] = None, env: Optional[Mapping[str, str]] = None, check: bool = True,
) -> CommandResult:
    proc = subprocess.run(
        list(argv), input=input_text, capture_output=True, text=True,
        cwd=str(cwd) if cwd else None, env=dict(env) if env else None, check=False,
    )
    stderr_path: Optional[Path] = None
    if log_dir is not None:
        log_dir.mkdir(parents=True, exist_ok=True)
        log_path = log_path or log_dir / f"{label}.log"
        stderr_path = log_dir / f"{label}.stderr.log"
        atomic_text(stderr_path, proc.stderr)
    if log_path is not None:
        atomic_text(log_path, f"command={shlex.join(argv)}\nreturncode={proc.returncode}\n--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}")
    if stdout_path is not None:
        atomic_text(stdout_path, proc.stdout)
    if check and proc.returncode != 0:
        suffix = f"; see {log_path}" if log_path else ""
        raise EvalError(f"command failed ({proc.returncode}): {shlex.join(argv)}{suffix}")
    return CommandResult(tuple(argv), proc.returncode, proc.stdout, proc.stderr,
                         str(stdout_path) if stdout_path else None,
                         str(stderr_path or log_path) if (stderr_path or log_path) else None)


def script_record() -> dict[str, Any]:
    package_dir = Path(__file__).resolve().parent
    records = [{"name": path.name, "sha256": sha256_file(path), "size": path.stat().st_size}
               for path in sorted(package_dir.glob("*.py"))]
    return {"version": __version__, "files": records, "digest": stable_hash(records)}


@dataclass
class RunContext:
    outdir: Path
    command: str
    settings: dict[str, Any]
    inputs: list[Any]
    check_only: Any = False
    hash_inputs: bool = False

    def __post_init__(self) -> None:
        self.outdir = Path(self.outdir).expanduser().resolve()
        self.tools: dict[str, Any] = {}
        if isinstance(self.check_only, dict):
            self.tools = self.check_only
            self.check_only = False
        self.inputs = [file_record(Path(item), self.hash_inputs) if isinstance(item, (str, Path)) else item
                       for item in self.inputs]
        self.logs, self.work, self.metrics = self.outdir / "logs", self.outdir / "work", self.outdir / "metrics"
        self.signature = stable_hash({
            "schema_version": SCHEMA_VERSION, "command": self.command,
            "settings": _json_safe(self.settings), "inputs": self.inputs,
            "tools": self.tools, "tool_source": script_record(),
        })
        self.manifest_path = self.outdir / "run_manifest.json"
        self.complete_path = self.outdir / "analysis.complete.json"
        self.lock_path = self.outdir / ".analysis.lock"
        self.started_at = now_iso()

    def begin(self) -> bool:
        if self.check_only:
            return False
        self.outdir.mkdir(parents=True, exist_ok=True)
        if self.lock_path.exists():
            raise EvalError(f"analysis lock exists: {self.lock_path}")
        if self.manifest_path.is_file():
            old = json.loads(self.manifest_path.read_text(encoding="utf-8"))
            if old.get("signature") != self.signature:
                raise EvalError(f"outdir contains a different run signature: {self.outdir}")
            if old.get("status") == "complete" and self.complete_path.is_file():
                return True
        atomic_text(self.lock_path, f"pid={os.getpid()}\nstarted_at={self.started_at}\n")
        atomic_json(self.manifest_path, self._manifest("running"))
        return False

    def _manifest(self, status: str, **extra: Any) -> dict[str, Any]:
        return {
            "schema_version": SCHEMA_VERSION, "package_version": __version__, "status": status,
            "signature": self.signature, "command": self.command, "started_at": self.started_at,
            "settings": _json_safe(self.settings), "inputs": self.inputs, "tools": self.tools,
            "tool_source": script_record(), **extra,
        }

    def complete(self, outputs: Sequence[Any], inputs_after: Optional[list[dict[str, Any]]] = None) -> None:
        after = inputs_after if inputs_after is not None else self.inputs
        if after != self.inputs:
            raise EvalError("one or more inputs changed during analysis")
        resolved = []
        for raw in outputs:
            path = Path(raw)
            resolved.append(path if path.is_absolute() else self.outdir / path)
        records = [file_record(path) for path in resolved]
        completed_at = now_iso()
        atomic_json(self.complete_path, {
            "schema_version": SCHEMA_VERSION, "package_version": __version__, "status": "complete",
            "signature": self.signature, "completed_at": completed_at,
            "inputs_unchanged": True, "outputs": records,
        })
        atomic_json(self.manifest_path, self._manifest(
            "complete", completed_at=completed_at, inputs_unchanged=True, outputs=records
        ))
        self.lock_path.unlink(missing_ok=True)

    def fail(self, exc: BaseException, inputs_after: Optional[list[dict[str, Any]]] = None) -> None:
        if self.check_only:
            return
        unchanged = inputs_after == self.inputs if inputs_after is not None else None
        atomic_json(self.manifest_path, self._manifest(
            "failed", failed_at=now_iso(), error=str(exc), inputs_unchanged=unchanged
        ))
        self.lock_path.unlink(missing_ok=True)


def main_guard(function: Any) -> int:
    try:
        return int(function())
    except KeyboardInterrupt:
        print("ERROR: interrupted by user", file=sys.stderr)
        return 130
    except EvalError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
