#!/usr/bin/env python3
"""Evaluate RaMAx phase-1 candidate windows against Alignathon truth.

Two evaluation modes are supported:

* direct MAF mode (recommended) streams the current candidate boundaries
  against an Alignathon truth MAF and does not require endpoint identity with
  an older RaMAx run;
* legacy replay mode consumes an existing breakpoint_classifications.tsv.gz
  file for backwards-compatible comparisons.

Neither mode modifies a HAL, MAF, or graph.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import math
import multiprocessing
import os
import shutil
import subprocess
import sys
import time
from array import array
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Iterable, Iterator, Sequence, TextIO


POSITIVE_LABEL = "over_fragmentation"
PROTECTED_LABELS = {
    "correct_truth_break",
    "avoided_cactus_overmerge",
    "truth_supported_break",
}
CLASSIFIED_LABELS = {POSITIVE_LABEL, *PROTECTED_LABELS}
DIRECT_TRUTH_VERSION = "phase1-direct-maf-v2"


def log(message: str, quiet: bool = False) -> None:
    if not quiet:
        print(f"[window-eval] {message}", file=sys.stderr, flush=True)


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV has no header: {path}")
        return list(reader)


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(path)


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    partial.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    partial.replace(path)


def as_int(row: dict[str, str], key: str, default: int = 0) -> int:
    value = row.get(key, "")
    if value in (None, ""):
        return default
    return int(value)


def as_bool(row: dict[str, str], key: str) -> bool:
    return row.get(key, "").strip().lower() in {"1", "true", "yes"}


def safe_ratio(numerator: int | float, denominator: int | float) -> float | None:
    return numerator / denominator if denominator else None


def fmt_ratio(value: float | None) -> str:
    return "NA" if value is None else f"{value:.6f}"


def percentile(values: Sequence[int], fraction: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return float(ordered[lower])
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def percentile_sorted(ordered: Sequence[int], fraction: float) -> float:
    if not ordered:
        return 0.0
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return float(ordered[lower])
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def n50(lengths: Sequence[int]) -> int:
    positive = sorted((value for value in lengths if value > 0), reverse=True)
    total = sum(positive)
    if total == 0:
        return 0
    threshold = total / 2
    cumulative = 0
    for value in positive:
        cumulative += value
        if cumulative >= threshold:
            return value
    return 0


@dataclass(frozen=True)
class MatchResult:
    truth: dict[str, str] | None
    method: str
    distance: int | None
    reason: str


class TruthMatcher:
    def __init__(self, rows: list[dict[str, str]], tolerance: int) -> None:
        self.rows = rows
        self.tolerance = tolerance
        # Sorted left-reference endpoint index per species pair/chromosome.
        # Tolerance matching first uses bisect on left_ref_start and then checks
        # the remaining seven endpoints. This turns the old O(E x T) scan into
        # O(E log T), where E is boundary evidence and T is truth breakpoints.
        self.buckets: dict[
            tuple[str, str, str], list[tuple[int, int]]
        ] = defaultdict(list)
        self.by_exact_reference: dict[
            tuple[str, str, str, int, int, int, int], list[dict[str, str]]
        ] = defaultdict(list)
        for index, row in enumerate(rows):
            pair_key = (
                row["reference_species"],
                row["target_species"],
                row["ref_chrom"],
            )
            exact_key = (
                *pair_key,
                as_int(row, "left_ref_start"),
                as_int(row, "left_ref_end"),
                as_int(row, "right_ref_start"),
                as_int(row, "right_ref_end"),
            )
            self.by_exact_reference[exact_key].append(row)
            self.buckets[pair_key].append(
                (as_int(row, "left_ref_start"), index)
            )
        for values in self.buckets.values():
            values.sort()

    @staticmethod
    def _target_compatible(evidence: dict[str, str], truth: dict[str, str], tolerance: int) -> bool:
        if evidence.get("left_target_chromosome") != truth.get("left_target_chrom"):
            return False
        if evidence.get("right_target_chromosome") != truth.get("right_target_chrom"):
            return False
        if evidence.get("left_strand") != truth.get("left_strand"):
            return False
        if evidence.get("right_strand") != truth.get("right_strand"):
            return False
        target_fields = (
            ("left_target_start", "left_target_start"),
            ("left_target_end", "left_target_end"),
            ("right_target_start", "right_target_start"),
            ("right_target_end", "right_target_end"),
        )
        return all(
            abs(as_int(evidence, evidence_key) - as_int(truth, truth_key)) <= tolerance
            for evidence_key, truth_key in target_fields
        )

    @staticmethod
    def _distance(boundary: dict[str, str], evidence: dict[str, str], truth: dict[str, str]) -> int:
        pairs = (
            (as_int(boundary, "left_ref_start"), as_int(truth, "left_ref_start")),
            (as_int(boundary, "left_ref_end"), as_int(truth, "left_ref_end")),
            (as_int(boundary, "right_ref_start"), as_int(truth, "right_ref_start")),
            (as_int(boundary, "right_ref_end"), as_int(truth, "right_ref_end")),
            (as_int(evidence, "left_target_start"), as_int(truth, "left_target_start")),
            (as_int(evidence, "left_target_end"), as_int(truth, "left_target_end")),
            (as_int(evidence, "right_target_start"), as_int(truth, "right_target_start")),
            (as_int(evidence, "right_target_end"), as_int(truth, "right_target_end")),
        )
        return sum(abs(left - right) for left, right in pairs)

    def match(self, boundary: dict[str, str], evidence: dict[str, str]) -> MatchResult:
        if boundary.get("coordinate_status") != "resolved":
            return MatchResult(None, "none", None, "reference_coordinate_unresolved")
        if evidence.get("coordinate_status") != "resolved":
            return MatchResult(None, "none", None, "target_coordinate_unresolved")
        if not as_bool(evidence, "unique_target"):
            return MatchResult(None, "none", None, "target_not_unique")

        pair_key = (
            boundary["reference_species"],
            evidence["target_species"],
            boundary["reference_chromosome"],
        )
        exact_key = (
            *pair_key,
            as_int(boundary, "left_ref_start"),
            as_int(boundary, "left_ref_end"),
            as_int(boundary, "right_ref_start"),
            as_int(boundary, "right_ref_end"),
        )
        exact_candidates = [
            row
            for row in self.by_exact_reference.get(exact_key, [])
            if self._target_compatible(evidence, row, 0)
        ]
        if len(exact_candidates) == 1:
            return MatchResult(exact_candidates[0], "exact", 0, "")
        if len(exact_candidates) > 1:
            labels = {row["classification"] for row in exact_candidates}
            if len(labels) == 1:
                chosen = min(exact_candidates, key=lambda row: int(row["breakpoint_id"]))
                return MatchResult(chosen, "exact_equivalent", 0, "")
            return MatchResult(None, "none", None, "ambiguous_exact_truth")

        if self.tolerance <= 0:
            return MatchResult(None, "none", None, "no_exact_match")

        bucket = self.buckets.get(pair_key, [])
        if not bucket:
            return MatchResult(None, "none", None, "no_tolerance_match")

        left_ref_start = as_int(boundary, "left_ref_start")
        left_ref_end = as_int(boundary, "left_ref_end")
        right_ref_start = as_int(boundary, "right_ref_start")
        right_ref_end = as_int(boundary, "right_ref_end")
        candidates: list[tuple[int, dict[str, str]]] = []
        lo = bisect.bisect_left(
            bucket, (left_ref_start - self.tolerance, -1)
        )
        hi = bisect.bisect_right(
            bucket, (left_ref_start + self.tolerance, len(self.rows))
        )
        for _, index in bucket[lo:hi]:
            truth = self.rows[index]
            reference_differences = (
                abs(left_ref_start - as_int(truth, "left_ref_start")),
                abs(left_ref_end - as_int(truth, "left_ref_end")),
                abs(right_ref_start - as_int(truth, "right_ref_start")),
                abs(right_ref_end - as_int(truth, "right_ref_end")),
            )
            if max(reference_differences) > self.tolerance:
                continue
            if not self._target_compatible(evidence, truth, self.tolerance):
                continue
            candidates.append((self._distance(boundary, evidence, truth), truth))

        if not candidates:
            return MatchResult(None, "none", None, "no_tolerance_match")
        candidates.sort(key=lambda item: (item[0], int(item[1]["breakpoint_id"])))
        best_distance = candidates[0][0]
        best = [truth for distance, truth in candidates if distance == best_distance]
        labels = {row["classification"] for row in best}
        if len(best) > 1 and len(labels) > 1:
            return MatchResult(None, "none", best_distance, "ambiguous_tolerance_truth")
        return MatchResult(best[0], "tolerance", best_distance, "")


@dataclass(frozen=True, slots=True)
class DirectQueryWindow:
    start: int
    end: int
    boundary_index: int
    side: int
    boundary_coordinate: int


@dataclass(frozen=True, slots=True)
class DirectMafRow:
    name: str
    species: str
    start: int
    size: int
    strand: str
    src_size: int
    sequence: bytes


@dataclass(frozen=True, slots=True)
class DirectCandidate:
    side: int
    ref_coord: int
    target_chrom: str
    target_coord: int
    strand: str
    distance: int


@dataclass(frozen=True, slots=True)
class DirectTruthResult:
    classification: str
    boundary_label: str
    reason: str
    left_candidates: tuple[DirectCandidate, ...]
    right_candidates: tuple[DirectCandidate, ...]
    left_truncated: int
    right_truncated: int
    bridge: tuple[DirectCandidate, DirectCandidate, int, int] | None
    compatible_pair_count: int
    ambiguous: bool


class DirectQueryIndex:
    """Compact interval index for current boundary flanks.

    Each boundary contributes two windows regardless of the number of target
    genomes. This is three times smaller than indexing every pairwise evidence
    row in the four-genome Alignathon run.
    """

    def __init__(
        self,
        boundaries: Sequence[dict[str, str]],
        flank: int,
    ) -> None:
        grouped: dict[str, list[DirectQueryWindow]] = defaultdict(list)
        for boundary_index, row in enumerate(boundaries):
            if row.get("coordinate_status") != "resolved":
                continue
            left_start = as_int(row, "left_ref_start")
            left_end = as_int(row, "left_ref_end")
            right_start = as_int(row, "right_ref_start")
            right_end = as_int(row, "right_ref_end")
            chromosome = row["reference_chromosome"]
            query_start = max(left_start, left_end - flank)
            if query_start < left_end:
                grouped[chromosome].append(
                    DirectQueryWindow(
                        query_start,
                        left_end,
                        boundary_index,
                        0,
                        left_end - 1,
                    )
                )
            query_end = min(right_end, right_start + flank)
            if right_start < query_end:
                grouped[chromosome].append(
                    DirectQueryWindow(
                        right_start,
                        query_end,
                        boundary_index,
                        1,
                        right_start,
                    )
                )
        self.windows: dict[str, list[DirectQueryWindow]] = {}
        self.starts: dict[str, list[int]] = {}
        self.flank = flank
        for chromosome, values in grouped.items():
            values.sort(
                key=lambda item: (
                    item.start,
                    item.end,
                    item.boundary_index,
                    item.side,
                )
            )
            self.windows[chromosome] = values
            self.starts[chromosome] = [item.start for item in values]

    @property
    def query_count(self) -> int:
        return sum(len(values) for values in self.windows.values())

    def overlaps(
        self, chromosome: str, start: int, end: int
    ) -> Sequence[DirectQueryWindow]:
        values = self.windows.get(chromosome, [])
        if not values:
            return ()
        starts = self.starts[chromosome]
        index = bisect.bisect_left(starts, start - self.flank)
        result: list[DirectQueryWindow] = []
        while index < len(values) and values[index].start < end:
            window = values[index]
            if window.end > start:
                result.append(window)
            index += 1
        return result


def canonical_sequence_name(sequence: str) -> str:
    fields = sequence.split(".")
    if len(fields) >= 3 and fields[0].startswith("sim") and fields[1] == fields[0]:
        return ".".join(fields[1:])
    return sequence


def sequence_species(sequence: str) -> str:
    return canonical_sequence_name(sequence).split(".", 1)[0]


def parse_direct_maf_s_line(raw: bytes, validate: bool) -> DirectMafRow:
    fields = raw.split()
    if len(fields) < 7:
        raise ValueError(f"Malformed MAF s line: {raw[:100]!r}")
    name = canonical_sequence_name(fields[1].decode("utf-8"))
    start = int(fields[2])
    size = int(fields[3])
    strand = fields[4].decode("ascii")
    src_size = int(fields[5])
    sequence = fields[6]
    if strand not in {"+", "-"}:
        raise ValueError(f"Invalid MAF strand for {name}: {strand}")
    if start < 0 or size < 0 or src_size <= 0 or start + size > src_size:
        raise ValueError(
            f"Invalid MAF coordinates: {name} {start} {size} {src_size}"
        )
    if validate and len(sequence) - sequence.count(b"-") != size:
        raise ValueError(f"MAF size mismatch for {name}: expected {size}")
    return DirectMafRow(
        name,
        sequence_species(name),
        start,
        size,
        strand,
        src_size,
        sequence,
    )


def direct_maf_s_line_name(raw: bytes) -> str:
    index = 1
    while index < len(raw) and raw[index] in b" \t":
        index += 1
    end = index
    while end < len(raw) and raw[end] not in b" \t\r\n":
        end += 1
    if end == index:
        raise ValueError(f"Malformed MAF s line: {raw[:100]!r}")
    return canonical_sequence_name(raw[index:end].decode("utf-8"))


def direct_row_interval(row: DirectMafRow) -> tuple[int, int]:
    if row.strand == "+":
        return row.start, row.start + row.size
    return row.src_size - row.start - row.size, row.src_size - row.start


def direct_coordinate_vector(row: DirectMafRow) -> list[int | None]:
    positions: list[int | None] = []
    consumed = 0
    if row.strand == "+":
        for base in row.sequence:
            if base == 45:
                positions.append(None)
            else:
                positions.append(row.start + consumed)
                consumed += 1
    else:
        first = row.src_size - row.start - 1
        for base in row.sequence:
            if base == 45:
                positions.append(None)
            else:
                positions.append(first - consumed)
                consumed += 1
    if consumed != row.size:
        raise ValueError(f"MAF coordinate mismatch for {row.name}")
    return positions


def _seek_to_complete_maf_block(handle: BinaryIO, start: int) -> None:
    if start <= 0:
        handle.seek(0)
        return
    handle.seek(start - 1)
    previous = handle.read(1)
    handle.seek(start)
    if previous != b"\n":
        handle.readline()
    position = handle.tell()
    line = handle.readline()
    if not line:
        return
    if not line.strip():
        return
    if line.startswith(b"a ") or line.startswith(b"a\t"):
        handle.seek(position)
        return
    while line and line.strip():
        line = handle.readline()


def iter_direct_maf_blocks_range(
    path: Path,
    start: int,
    end: int,
    max_blocks: int,
    validate: bool,
) -> Iterator[list[DirectMafRow]]:
    blocks = 0
    if path.suffix == ".gz":
        if start != 0 or end != -1:
            raise ValueError("gzip truth MAF supports only a full sequential scan")
        context = gzip.open(path, "rb")
    else:
        context = path.open("rb")
    with context as handle:
        if path.suffix != ".gz":
            _seek_to_complete_maf_block(handle, start)
        rows: list[DirectMafRow] = []
        while True:
            if not rows and end >= 0 and handle.tell() >= end:
                return
            raw = handle.readline()
            if not raw:
                if rows:
                    yield rows
                return
            if raw.startswith(b"s ") or raw.startswith(b"s\t"):
                name = direct_maf_s_line_name(raw)
                if sequence_species(name) in _DIRECT_INPUT_SPECIES:
                    rows.append(parse_direct_maf_s_line(raw, validate))
            elif not raw.strip() and rows:
                yield rows
                blocks += 1
                if max_blocks and blocks >= max_blocks:
                    return
                rows = []


_DIRECT_QUERY_INDEX: DirectQueryIndex | None = None
_DIRECT_INPUT_SPECIES: tuple[str, ...] = ()
_DIRECT_VALIDATE_MAF = False


def _scan_direct_truth_range(task: dict) -> dict:
    if _DIRECT_QUERY_INDEX is None:
        raise RuntimeError("Direct truth worker has no query index")
    truth_path = Path(task["truth_maf"])
    output_path = Path(task["output"])
    partial = output_path.with_name(output_path.name + ".partial")
    input_species = set(_DIRECT_INPUT_SPECIES)
    blocks = 0
    relevant_blocks = 0
    candidate_count = 0
    started = time.monotonic()
    with partial.open("w", encoding="utf-8", newline="") as destination:
        for rows in iter_direct_maf_blocks_range(
            truth_path,
            int(task["start"]),
            int(task["end"]),
            int(task["max_blocks"]),
            _DIRECT_VALIDATE_MAF,
        ):
            blocks += 1
            by_species: dict[str, list[DirectMafRow]] = defaultdict(list)
            for row in rows:
                if row.species in input_species:
                    by_species[row.species].append(row)
            block_relevant = False
            vectors: dict[int, list[int | None]] = {}
            for ref_species, ref_rows in by_species.items():
                for ref_row in ref_rows:
                    ref_start, ref_end = direct_row_interval(ref_row)
                    windows = _DIRECT_QUERY_INDEX.overlaps(
                        ref_row.name, ref_start, ref_end
                    )
                    if not windows:
                        continue
                    block_relevant = True
                    ref_key = id(ref_row)
                    ref_positions = vectors.get(ref_key)
                    if ref_positions is None:
                        ref_positions = direct_coordinate_vector(ref_row)
                        vectors[ref_key] = ref_positions
                    for target_species, target_rows in by_species.items():
                        if target_species == ref_species:
                            continue
                        for target_row in target_rows:
                            target_key = id(target_row)
                            target_positions = vectors.get(target_key)
                            if target_positions is None:
                                target_positions = direct_coordinate_vector(target_row)
                                vectors[target_key] = target_positions
                            if len(ref_positions) != len(target_positions):
                                raise ValueError(
                                    "MAF rows have unequal alignment lengths in one block"
                                )
                            paired = [
                                (ref_coord, target_coord)
                                for ref_coord, target_coord in zip(
                                    ref_positions, target_positions
                                )
                                if ref_coord is not None and target_coord is not None
                            ]
                            if not paired:
                                continue
                            if paired[0][0] > paired[-1][0]:
                                paired.reverse()
                            ref_coordinates = [value[0] for value in paired]
                            target_coordinates = [value[1] for value in paired]
                            relative_strand = (
                                "+" if ref_row.strand == target_row.strand else "-"
                            )
                            for window in windows:
                                if window.side == 0:
                                    position = (
                                        bisect.bisect_left(
                                            ref_coordinates, window.end
                                        )
                                        - 1
                                    )
                                    if (
                                        position < 0
                                        or ref_coordinates[position] < window.start
                                    ):
                                        continue
                                else:
                                    position = bisect.bisect_left(
                                        ref_coordinates, window.start
                                    )
                                    if (
                                        position >= len(ref_coordinates)
                                        or ref_coordinates[position] >= window.end
                                    ):
                                        continue
                                ref_coord = ref_coordinates[position]
                                target_coord = target_coordinates[position]
                                distance = abs(
                                    window.boundary_coordinate - ref_coord
                                )
                                destination.write(
                                    f"{window.boundary_index}\t{window.side}\t"
                                    f"{target_species}\t{target_row.name}\t"
                                    f"{relative_strand}\t{ref_coord}\t"
                                    f"{target_coord}\t{distance}\n"
                                )
                                candidate_count += 1
            if block_relevant:
                relevant_blocks += 1
            if blocks % 250000 == 0:
                print(
                    f"[window-eval worker {os.getpid()}] "
                    f"{blocks:,} MAF blocks; {candidate_count:,} candidates",
                    file=sys.stderr,
                    flush=True,
                )
    partial.replace(output_path)
    return {
        "output": str(output_path),
        "range_start": int(task["start"]),
        "range_end": int(task["end"]),
        "maf_blocks": blocks,
        "relevant_maf_blocks": relevant_blocks,
        "raw_candidate_count": candidate_count,
        "elapsed_seconds": time.monotonic() - started,
    }


def file_fingerprint(path: Path) -> dict[str, str | int]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def signature_for(payload: dict) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def cache_marker_valid(path: Path, marker: Path, signature: str) -> bool:
    if not path.is_file() or path.stat().st_size == 0 or not marker.is_file():
        return False
    try:
        record = json.loads(marker.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return False
    return (
        record.get("signature") == signature
        and record.get("output_size") == path.stat().st_size
    )


def write_cache_marker(
    path: Path, marker: Path, signature: str, payload: dict
) -> None:
    record = dict(payload)
    record.update(
        {
            "signature": signature,
            "output": str(path.resolve()),
            "output_size": path.stat().st_size,
            "created_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        }
    )
    write_json(marker, record)


DIRECT_CANDIDATE_FIELDS = (
    "boundary_index",
    "side",
    "target_species",
    "target_chromosome",
    "strand",
    "ref_coord",
    "target_coord",
    "distance",
    "truncated_candidates",
)


def reduce_direct_candidates(
    sorted_path: Path,
    output_path: Path,
    max_candidates: int,
) -> tuple[int, int]:
    partial = output_path.with_name(output_path.name + ".partial")
    retained = 0
    truncated = 0

    def flush_group(
        writer: csv.DictWriter,
        rows: list[list[str]],
        unique_count: int,
    ) -> tuple[int, int]:
        dropped = max(0, unique_count - len(rows))
        for fields in rows:
            writer.writerow(
                dict(
                    zip(
                        DIRECT_CANDIDATE_FIELDS,
                        [*fields, str(dropped)],
                    )
                )
            )
        return len(rows), dropped

    with sorted_path.open("r", encoding="utf-8", newline="") as source, gzip.open(
        partial, "wt", encoding="utf-8", newline="", compresslevel=3
    ) as destination:
        writer = csv.DictWriter(
            destination,
            fieldnames=DIRECT_CANDIDATE_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        current_group: tuple[str, str, str] | None = None
        kept: list[list[str]] = []
        seen: set[tuple[str, str, str, str]] = set()
        unique_count = 0
        for raw in source:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 8:
                raise ValueError(f"Malformed direct truth candidate: {raw[:120]}")
            group = (fields[0], fields[1], fields[2])
            if current_group is not None and group != current_group:
                added, dropped = flush_group(writer, kept, unique_count)
                retained += added
                truncated += dropped
                kept = []
                seen = set()
                unique_count = 0
            current_group = group
            candidate_key = (fields[3], fields[4], fields[5], fields[6])
            if candidate_key in seen:
                continue
            seen.add(candidate_key)
            unique_count += 1
            if len(kept) < max_candidates:
                kept.append(fields)
        if current_group is not None:
            added, dropped = flush_group(writer, kept, unique_count)
            retained += added
            truncated += dropped
    partial.replace(output_path)
    return retained, truncated


def ensure_direct_truth_candidates(
    truth_maf: Path,
    boundaries: Sequence[dict[str, str]],
    boundary_sources: Sequence[Path],
    input_species: Sequence[str],
    cache_dir: Path,
    flank: int,
    max_candidates: int,
    threads: int,
    sort_memory: str,
    validate_maf: bool,
    max_maf_blocks: int,
    keep_shards: bool,
    quiet: bool,
) -> tuple[Path, dict]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    signature_payload = {
        "version": DIRECT_TRUTH_VERSION,
        "truth_maf": file_fingerprint(truth_maf),
        "boundary_sources": [file_fingerprint(path) for path in boundary_sources],
        "input_species": list(input_species),
        "flank": flank,
        "max_candidates": max_candidates,
        "validate_maf": validate_maf,
        "max_maf_blocks": max_maf_blocks,
    }
    signature = signature_for(signature_payload)
    reduced_path = cache_dir / f"truth_candidates.{signature[:16]}.tsv.gz"
    reduced_marker = reduced_path.with_name(reduced_path.name + ".complete.json")
    if cache_marker_valid(reduced_path, reduced_marker, signature):
        record = json.loads(reduced_marker.read_text(encoding="utf-8"))
        log(f"reusing direct truth cache: {reduced_path}", quiet)
        return reduced_path, record

    query_index = DirectQueryIndex(boundaries, flank)
    query_count = query_index.query_count
    query_sequence_count = len(query_index.windows)
    log(
        f"direct truth index: {query_count:,} flank queries "
        f"across {query_sequence_count:,} reference sequences",
        quiet,
    )
    if query_count == 0:
        raise ValueError("No resolved boundary flanks are available for truth lookup")

    global _DIRECT_QUERY_INDEX, _DIRECT_INPUT_SPECIES, _DIRECT_VALIDATE_MAF
    _DIRECT_QUERY_INDEX = query_index
    _DIRECT_INPUT_SPECIES = tuple(input_species)
    _DIRECT_VALIDATE_MAF = validate_maf

    if truth_maf.suffix == ".gz":
        if threads > 1:
            log("gzip truth MAF cannot be byte-sharded; using one worker", quiet)
        ranges = [(0, -1)]
        threads = 1
    elif max_maf_blocks:
        if threads > 1:
            log("--max-maf-blocks forces a single truth worker", quiet)
        ranges = [(0, truth_maf.stat().st_size)]
        threads = 1
    else:
        file_size = truth_maf.stat().st_size
        worker_count = max(1, min(threads, file_size or 1))
        ranges = [
            (
                file_size * index // worker_count,
                file_size * (index + 1) // worker_count,
            )
            for index in range(worker_count)
        ]

    tasks: list[dict] = []
    raw_paths: list[Path] = []
    worker_stats: list[dict] = []
    for index, (start, end) in enumerate(ranges):
        raw_path = cache_dir / f"truth_candidates.{signature[:16]}.part-{index:03d}.tsv"
        marker = raw_path.with_name(raw_path.name + ".complete.json")
        shard_signature = signature_for(
            {
                "parent": signature,
                "range_start": start,
                "range_end": end,
                "part": index,
            }
        )
        raw_paths.append(raw_path)
        if cache_marker_valid(raw_path, marker, shard_signature):
            worker_stats.append(json.loads(marker.read_text(encoding="utf-8")))
            continue
        tasks.append(
            {
                "truth_maf": str(truth_maf),
                "output": str(raw_path),
                "marker": str(marker),
                "signature": shard_signature,
                "start": start,
                "end": end,
                "max_blocks": max_maf_blocks,
                "part": index,
            }
        )

    def record_worker_result(task: dict, result: dict) -> None:
        marker = Path(task["marker"])
        output = Path(result["output"])
        write_cache_marker(output, marker, str(task["signature"]), result)
        worker_stats.append(result)
        log(
            f"truth shard {task['part'] + 1}/{len(ranges)}: "
            f"{result['maf_blocks']:,} blocks, "
            f"{result['raw_candidate_count']:,} candidates, "
            f"{result['elapsed_seconds']:.1f}s",
            quiet,
        )

    if tasks:
        use_fork = "fork" in multiprocessing.get_all_start_methods()
        if threads > 1 and use_fork:
            context = multiprocessing.get_context("fork")
            with ProcessPoolExecutor(
                max_workers=min(threads, len(tasks)), mp_context=context
            ) as executor:
                future_to_task = {
                    executor.submit(_scan_direct_truth_range, task): task
                    for task in tasks
                }
                for future in as_completed(future_to_task):
                    task = future_to_task[future]
                    record_worker_result(task, future.result())
        else:
            if threads > 1 and not use_fork:
                log(
                    "fork multiprocessing is unavailable; using one truth worker",
                    quiet,
                )
            for task in tasks:
                record_worker_result(task, _scan_direct_truth_range(task))

    _DIRECT_QUERY_INDEX = None
    del query_index

    sort_executable = shutil.which("sort")
    if sort_executable is None:
        raise FileNotFoundError("GNU sort is required for direct truth evaluation")
    sorted_path = cache_dir / f"truth_candidates.{signature[:16]}.sorted.tsv"
    sorted_partial = sorted_path.with_name(sorted_path.name + ".partial")
    sort_tmp = cache_dir / f"sort-tmp-{signature[:12]}"
    sort_tmp.mkdir(parents=True, exist_ok=True)
    sort_command = [
        sort_executable,
        "--field-separator=\t",
        f"--parallel={max(1, threads)}",
        "--buffer-size",
        sort_memory,
        "--temporary-directory",
        str(sort_tmp),
        "-k1,1n",
        "-k2,2n",
        "-k3,3",
        "-k8,8n",
        "-k4,4",
        "-k5,5",
        "-k6,6n",
        "-k7,7n",
        *[str(path) for path in raw_paths],
        "--output",
        str(sorted_partial),
    ]
    environment = dict(os.environ)
    environment["LC_ALL"] = "C"
    log(
        f"sorting direct truth candidates with {max(1, threads)} workers",
        quiet,
    )
    completed = subprocess.run(sort_command, env=environment, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            f"Direct truth candidate sort failed with status {completed.returncode}"
        )
    sorted_partial.replace(sorted_path)

    retained, truncated = reduce_direct_candidates(
        sorted_path, reduced_path, max_candidates
    )
    stats = {
        "mode": "direct_maf",
        "version": DIRECT_TRUTH_VERSION,
        "signature_payload": signature_payload,
        "query_count": query_count,
        "query_sequence_count": query_sequence_count,
        "worker_count": len(ranges),
        "worker_stats": sorted(
            worker_stats, key=lambda row: int(row.get("range_start", 0))
        ),
        "maf_block_count": sum(
            int(row.get("maf_blocks", 0)) for row in worker_stats
        ),
        "relevant_maf_block_count": sum(
            int(row.get("relevant_maf_blocks", 0)) for row in worker_stats
        ),
        "raw_candidate_count": sum(
            int(row.get("raw_candidate_count", 0)) for row in worker_stats
        ),
        "retained_candidate_count": retained,
        "truncated_candidate_count": truncated,
    }
    write_cache_marker(reduced_path, reduced_marker, signature, stats)
    if not keep_shards:
        for path in raw_paths:
            path.unlink(missing_ok=True)
            path.with_name(path.name + ".complete.json").unlink(missing_ok=True)
        sorted_path.unlink(missing_ok=True)
        try:
            sort_tmp.rmdir()
        except OSError:
            pass
    log(
        f"direct truth cache complete: {retained:,} retained candidates; "
        f"{truncated:,} truncated",
        quiet,
    )
    return reduced_path, stats


def iter_direct_candidate_groups(
    path: Path,
) -> Iterator[
    tuple[
        int,
        dict[tuple[int, str], list[DirectCandidate]],
        dict[tuple[int, str], int],
    ]
]:
    current_index: int | None = None
    candidates: dict[tuple[int, str], list[DirectCandidate]] = defaultdict(list)
    truncated: dict[tuple[int, str], int] = defaultdict(int)
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            boundary_index = int(row["boundary_index"])
            if current_index is None:
                current_index = boundary_index
            if boundary_index != current_index:
                yield current_index, dict(candidates), dict(truncated)
                current_index = boundary_index
                candidates = defaultdict(list)
                truncated = defaultdict(int)
            key = (int(row["side"]), row["target_species"])
            candidates[key].append(
                DirectCandidate(
                    int(row["side"]),
                    int(row["ref_coord"]),
                    row["target_chromosome"],
                    int(row["target_coord"]),
                    row["strand"],
                    int(row["distance"]),
                )
            )
            truncated[key] = max(
                truncated[key], int(row["truncated_candidates"])
            )
    if current_index is not None:
        yield current_index, dict(candidates), dict(truncated)


def compatible_direct_truth_pair(
    left: DirectCandidate,
    right: DirectCandidate,
    max_gap: int,
) -> tuple[int, int] | None:
    if left.target_chrom != right.target_chrom or left.strand != right.strand:
        return None
    ref_gap = right.ref_coord - left.ref_coord - 1
    if left.strand == "+":
        target_gap = right.target_coord - left.target_coord - 1
    else:
        target_gap = left.target_coord - right.target_coord - 1
    if 0 <= ref_gap <= max_gap and 0 <= target_gap <= max_gap:
        return ref_gap, target_gap
    return None


def infer_direct_truth_break_reason(
    left: Sequence[DirectCandidate],
    right: Sequence[DirectCandidate],
    max_gap: int,
) -> str:
    if not left or not right:
        return "insufficient_truth_support"
    if not any(a.target_chrom == b.target_chrom for a in left for b in right):
        return "target_switch"
    if not any(
        a.target_chrom == b.target_chrom and a.strand == b.strand
        for a in left
        for b in right
    ):
        return "strand_switch"
    saw_order_break = False
    saw_large_gap = False
    for left_candidate in left:
        for right_candidate in right:
            if (
                left_candidate.target_chrom != right_candidate.target_chrom
                or left_candidate.strand != right_candidate.strand
            ):
                continue
            ref_gap = right_candidate.ref_coord - left_candidate.ref_coord - 1
            target_gap = (
                right_candidate.target_coord - left_candidate.target_coord - 1
                if left_candidate.strand == "+"
                else left_candidate.target_coord - right_candidate.target_coord - 1
            )
            if ref_gap < 0 or target_gap < 0:
                saw_order_break = True
            elif ref_gap > max_gap or target_gap > max_gap:
                saw_large_gap = True
    if saw_order_break:
        return "order_break"
    if saw_large_gap:
        return "large_gap_break"
    return "ambiguous_incompatibility"


def classify_direct_truth(
    target_species: str,
    candidates: dict[tuple[int, str], list[DirectCandidate]],
    truncated: dict[tuple[int, str], int],
    max_gap: int,
) -> DirectTruthResult:
    left = tuple(candidates.get((0, target_species), ()))
    right = tuple(candidates.get((1, target_species), ()))
    left_truncated = truncated.get((0, target_species), 0)
    right_truncated = truncated.get((1, target_species), 0)
    if not left or not right:
        return DirectTruthResult(
            "insufficient_truth_support",
            "UNRESOLVED_BOUNDARY",
            "insufficient_truth_support",
            left,
            right,
            left_truncated,
            right_truncated,
            None,
            0,
            False,
        )

    compatible: list[
        tuple[
            tuple[int, int, str, str, int, int],
            DirectCandidate,
            DirectCandidate,
            int,
            int,
        ]
    ] = []
    for left_candidate in left:
        for right_candidate in right:
            gaps = compatible_direct_truth_pair(
                left_candidate, right_candidate, max_gap
            )
            if gaps is None:
                continue
            ref_gap, target_gap = gaps
            score = (
                left_candidate.distance + right_candidate.distance,
                abs(ref_gap - target_gap),
                left_candidate.target_chrom,
                left_candidate.strand,
                left_candidate.target_coord,
                right_candidate.target_coord,
            )
            compatible.append(
                (
                    score,
                    left_candidate,
                    right_candidate,
                    ref_gap,
                    target_gap,
                )
            )
    if compatible:
        compatible.sort(key=lambda item: item[0])
        best = compatible[0]
        best_primary_score = best[0][:2]
        equivalent_loci = {
            (
                item[1].target_chrom,
                item[1].strand,
                item[1].target_coord,
                item[2].target_coord,
            )
            for item in compatible
            if item[0][:2] == best_primary_score
        }
        ambiguous = len(equivalent_loci) > 1
        if ambiguous:
            return DirectTruthResult(
                "copy_ambiguous",
                "UNRESOLVED_BOUNDARY",
                "multiple_equivalent_truth_bridges",
                left,
                right,
                left_truncated,
                right_truncated,
                None,
                len(compatible),
                True,
            )
        return DirectTruthResult(
            POSITIVE_LABEL,
            "MERGE_BOUNDARY",
            "continuous",
            left,
            right,
            left_truncated,
            right_truncated,
            (best[1], best[2], best[3], best[4]),
            len(compatible),
            False,
        )

    if left_truncated or right_truncated:
        return DirectTruthResult(
            "insufficient_truth_support",
            "UNRESOLVED_BOUNDARY",
            "truth_candidate_limit_reached",
            left,
            right,
            left_truncated,
            right_truncated,
            None,
            0,
            False,
        )
    return DirectTruthResult(
        "truth_supported_break",
        "KEEP_BREAK",
        infer_direct_truth_break_reason(left, right, max_gap),
        left,
        right,
        left_truncated,
        right_truncated,
        None,
        0,
        False,
    )


def discover_round_dirs(report_dir: Path, selection: str) -> list[Path]:
    round_dirs = sorted(
        path for path in report_dir.glob("round_*") if path.is_dir()
    )
    if not round_dirs:
        raise FileNotFoundError(f"No round_* directories found under {report_dir}")
    if selection == "all":
        return round_dirs
    if selection == "latest":
        return [round_dirs[-1]]
    normalized = selection if selection.startswith("round_") else f"round_{int(selection):03d}"
    selected = report_dir / normalized
    if selected not in round_dirs:
        raise FileNotFoundError(f"Requested round directory not found: {selected}")
    return [selected]


def load_reports(round_dirs: Sequence[Path]) -> tuple[
    list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]
]:
    windows: list[dict[str, str]] = []
    boundaries: list[dict[str, str]] = []
    boundary_genomes: list[dict[str, str]] = []
    blocks: list[dict[str, str]] = []
    for round_dir in round_dirs:
        required = {
            "windows.tsv": windows,
            "window_boundaries.tsv": boundaries,
            "boundary_genomes.tsv": boundary_genomes,
            "window_blocks.tsv": blocks,
        }
        for name, destination in required.items():
            path = round_dir / name
            if not path.is_file():
                raise FileNotFoundError(f"Missing report file: {path}")
            rows = read_tsv(path)
            for row in rows:
                row["source_round_dir"] = round_dir.name
            destination.extend(rows)
    return windows, boundaries, boundary_genomes, blocks


def load_direct_reports(
    round_dirs: Sequence[Path],
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[Path]]:
    windows: list[dict[str, str]] = []
    boundaries: list[dict[str, str]] = []
    boundary_sources: list[Path] = []
    for round_dir in round_dirs:
        windows_path = round_dir / "windows.tsv"
        boundaries_path = round_dir / "window_boundaries.tsv"
        evidence_path = round_dir / "boundary_genomes.tsv"
        for path in (windows_path, boundaries_path, evidence_path):
            if not path.is_file():
                raise FileNotFoundError(f"Missing report file: {path}")
        for row in read_tsv(windows_path):
            row["source_round_dir"] = round_dir.name
            windows.append(row)
        for row in read_tsv(boundaries_path):
            row["source_round_dir"] = round_dir.name
            boundaries.append(row)
        boundary_sources.append(boundaries_path)
    return windows, boundaries, boundary_sources


def iter_boundary_evidence_groups(
    round_dirs: Sequence[Path],
) -> Iterator[tuple[tuple[str, str, str], list[dict[str, str]]]]:
    current_key: tuple[str, str, str] | None = None
    current_rows: list[dict[str, str]] = []
    for round_dir in round_dirs:
        path = round_dir / "boundary_genomes.tsv"
        with open_text(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"TSV has no header: {path}")
            for row in reader:
                key = (round_dir.name, row["window_id"], row["boundary_id"])
                if current_key is None:
                    current_key = key
                if key != current_key:
                    yield current_key, current_rows
                    current_key = key
                    current_rows = []
                current_rows.append(row)
    if current_key is not None:
        yield current_key, current_rows


def direct_candidate_matches_prediction(
    candidate: DirectCandidate,
    evidence: dict[str, str],
    side: str,
    tolerance: int,
) -> bool:
    if evidence.get("coordinate_status") != "resolved":
        return False
    prefix = "left" if side == "left" else "right"
    if not as_bool(evidence, f"{prefix}_present"):
        return False
    chromosome = evidence.get(f"{prefix}_target_chromosome", "")
    strand = evidence.get(f"{prefix}_strand", "")
    start = as_int(evidence, f"{prefix}_target_start")
    end = as_int(evidence, f"{prefix}_target_end")
    return (
        candidate.target_chrom == chromosome
        and candidate.strand == strand
        and start - tolerance <= candidate.target_coord < end + tolerance
    )


DIRECT_BOUNDARY_FIELDS = [
    "window_id",
    "boundary_id",
    "source_round_dir",
    "detector_priority",
    "detector_recommends_bridge",
    "signal_type",
    "reference_species",
    "target_species",
    "reference_chromosome",
    "left_ref_start",
    "left_ref_end",
    "right_ref_start",
    "right_ref_end",
    "reference_gap",
    "min_adjacent_block_length",
    "truth_source",
    "truth_classification",
    "boundary_truth_label",
    "truth_break_reason",
    "truth_left_candidates",
    "truth_right_candidates",
    "truth_left_truncated",
    "truth_right_truncated",
    "truth_bridge",
    "truth_bridge_target_chromosome",
    "truth_bridge_strand",
    "truth_left_ref_coord",
    "truth_left_target_coord",
    "truth_right_ref_coord",
    "truth_right_target_coord",
    "truth_ref_gap",
    "truth_target_gap",
    "truth_compatible_pair_count",
    "truth_mapping_ambiguous",
    "truth_left_matches_prediction",
    "truth_right_matches_prediction",
    "truth_bridge_matches_prediction",
    "truth_bridge_relation",
    "replay_match_method",
    "replay_match_distance",
    "replay_unmatched_reason",
    "replay_truth_breakpoint_id",
    "replay_truth_classification",
    "cactus_bridges",
]


def direct_scope_metrics(
    name: str,
    tiers: set[str],
    tier_counts: dict[str, Counter],
    window_labels: Sequence[dict[str, str | int]],
) -> dict:
    counts: Counter = Counter()
    for tier in tiers:
        counts.update(tier_counts.get(tier, Counter()))
    classified = counts[POSITIVE_LABEL] + counts["truth_supported_break"]
    protected = counts["truth_supported_break"]
    total = sum(counts.values())
    selected_windows = [
        row for row in window_labels if str(row["priority_tier"]) in tiers
    ]
    mixed = sum(
        row["truth_window_label"] == "MIXED_WINDOW"
        for row in selected_windows
    )
    return {
        "scope": name,
        "tiers": sorted(tiers),
        "candidate_pairwise_boundary_count": total,
        "directly_classified_pairwise_boundary_count": classified,
        "over_fragmentation_detected": counts[POSITIVE_LABEL],
        "truth_supported_break_detected": protected,
        "insufficient_truth_support_count": counts[
            "insufficient_truth_support"
        ],
        "copy_ambiguous_count": counts["copy_ambiguous"],
        "classified_precision": safe_ratio(counts[POSITIVE_LABEL], classified),
        "over_fragmentation_recall": None,
        "recall_status": "not_evaluable_without_all_graph_boundaries",
        "protected_boundary_contamination": safe_ratio(protected, classified),
        "direct_truth_coverage": safe_ratio(classified, total),
        "selected_window_count": len(selected_windows),
        "mixed_window_count": mixed,
        "mixed_window_rate": safe_ratio(mixed, len(selected_windows)),
    }


def direct_acceptance(metrics: dict[str, dict]) -> dict:
    tier_a = metrics["tier_a"]
    tier_ab = metrics["tier_a_b"]

    def at_least(value: float | None, threshold: float) -> bool | None:
        return None if value is None else value >= threshold

    def at_most(value: float | None, threshold: float) -> bool | None:
        return None if value is None else value <= threshold

    checks: dict[str, bool | None] = {
        "tier_a_precision_ge_0_95": at_least(
            tier_a["classified_precision"], 0.95
        ),
        "tier_a_recall_ge_0_40": None,
        "tier_a_protected_contamination_le_0_05": at_most(
            tier_a["protected_boundary_contamination"], 0.05
        ),
        "tier_a_b_precision_ge_0_75": at_least(
            tier_ab["classified_precision"], 0.75
        ),
        "tier_a_b_recall_ge_0_75": None,
        "tier_a_b_mixed_window_rate_le_0_05": at_most(
            tier_ab["mixed_window_rate"], 0.05
        ),
    }
    evaluated = [value for value in checks.values() if value is not None]
    return {
        "checks": checks,
        "all_evaluable_checks_passed": bool(evaluated) and all(evaluated),
        "fully_evaluable": False,
        "not_evaluable": [name for name, value in checks.items() if value is None],
    }


def feature_rows_from_arrays(
    values: dict[str, tuple[array, array]],
) -> list[dict[str, str | int]]:
    output: list[dict[str, str | int]] = []
    for label, (lengths_array, gaps_array) in sorted(values.items()):
        lengths = sorted(lengths_array)
        gaps = sorted(gaps_array)
        output.append(
            {
                "truth_classification": label,
                "count": len(lengths),
                "min_block_p10": f"{percentile_sorted(lengths, 0.10):.3f}",
                "min_block_p50": f"{percentile_sorted(lengths, 0.50):.3f}",
                "min_block_p90": f"{percentile_sorted(lengths, 0.90):.3f}",
                "reference_gap_p10": f"{percentile_sorted(gaps, 0.10):.3f}",
                "reference_gap_p50": f"{percentile_sorted(gaps, 0.50):.3f}",
                "reference_gap_p90": f"{percentile_sorted(gaps, 0.90):.3f}",
            }
        )
    return output


def map_truth_label(classification: str) -> str:
    return {
        "over_fragmentation": "MERGE_BOUNDARY",
        "correct_truth_break": "KEEP_BREAK",
        "avoided_cactus_overmerge": "KEEP_CONSERVATIVE_BREAK",
        "insufficient_truth_support": "UNRESOLVED_BOUNDARY",
    }.get(classification, "UNRESOLVED_BOUNDARY")


def evaluate_boundaries(
    boundaries: list[dict[str, str]],
    boundary_genomes: list[dict[str, str]],
    matcher: TruthMatcher,
    quiet: bool = False,
) -> list[dict[str, str | int]]:
    boundary_by_key = {
        (row["window_id"], row["boundary_id"]): row for row in boundaries
    }
    output: list[dict[str, str | int]] = []
    seen_evidence: set[tuple[str, str, str]] = set()
    for row_index, evidence in enumerate(boundary_genomes, start=1):
        evidence_key = (
            evidence["window_id"],
            evidence["boundary_id"],
            evidence["target_species"],
        )
        if evidence_key in seen_evidence:
            continue
        seen_evidence.add(evidence_key)
        boundary = boundary_by_key.get(evidence_key[:2])
        if boundary is None:
            raise ValueError(f"boundary_genomes row has no parent boundary: {evidence_key}")
        result = matcher.match(boundary, evidence)
        truth = result.truth or {}
        classification = truth.get("classification", "unmatched")
        output.append(
            {
                "window_id": evidence["window_id"],
                "boundary_id": evidence["boundary_id"],
                "source_round_dir": boundary.get("source_round_dir", ""),
                "detector_priority": boundary["detector_priority"],
                "detector_recommends_bridge": boundary["detector_recommends_bridge"],
                "signal_type": boundary["signal_type"],
                "reference_species": boundary["reference_species"],
                "target_species": evidence["target_species"],
                "reference_chromosome": boundary["reference_chromosome"],
                "left_ref_start": boundary["left_ref_start"],
                "left_ref_end": boundary["left_ref_end"],
                "right_ref_start": boundary["right_ref_start"],
                "right_ref_end": boundary["right_ref_end"],
                "reference_gap": boundary["reference_gap"],
                "min_adjacent_block_length": boundary["min_adjacent_block_length"],
                "match_method": result.method,
                "match_distance": "" if result.distance is None else result.distance,
                "unmatched_reason": result.reason,
                "truth_breakpoint_id": truth.get("breakpoint_id", ""),
                "truth_classification": classification,
                "boundary_truth_label": map_truth_label(classification),
                "truth_break_reason": truth.get("truth_break_reason", ""),
                "cactus_bridges": truth.get("cactus_bridges", ""),
            }
        )
        if row_index % 100000 == 0:
            log(
                f"matched {row_index:,}/{len(boundary_genomes):,} "
                "boundary-genome rows",
                quiet,
            )
    output.sort(
        key=lambda row: (
            str(row["source_round_dir"]),
            str(row["window_id"]),
            str(row["boundary_id"]),
            str(row["target_species"]),
        )
    )
    log(f"matched {len(output):,} boundary-genome rows", quiet)
    return output


def classify_windows(
    windows: list[dict[str, str]],
    boundary_labels: list[dict[str, str | int]],
) -> list[dict[str, str | int]]:
    labels_by_window: dict[str, list[dict[str, str | int]]] = defaultdict(list)
    for row in boundary_labels:
        labels_by_window[str(row["window_id"])].append(row)

    output: list[dict[str, str | int]] = []
    for window in windows:
        evidence = labels_by_window.get(window["window_id"], [])
        classifications = [str(row["truth_classification"]) for row in evidence]
        positive_ids = {
            str(row["truth_breakpoint_id"])
            for row in evidence
            if row["truth_classification"] == POSITIVE_LABEL
            and row["truth_breakpoint_id"]
        }
        protected_ids = {
            str(row["truth_breakpoint_id"])
            for row in evidence
            if row["truth_classification"] in PROTECTED_LABELS
            and row["truth_breakpoint_id"]
        }
        if positive_ids and protected_ids:
            label = "MIXED_WINDOW"
        elif positive_ids:
            label = "TRUE_POSITIVE_WINDOW"
        elif protected_ids:
            label = "FALSE_POSITIVE_WINDOW"
        else:
            label = "UNRESOLVED_WINDOW"
        output.append(
            {
                "window_id": window["window_id"],
                "source_round_dir": window.get("source_round_dir", ""),
                "priority_tier": window["priority_tier"],
                "report_species": window["report_species"],
                "report_chromosome": window["report_chromosome"],
                "original_start": window["original_start"],
                "original_end": window["original_end"],
                "coordinate_status": window["coordinate_status"],
                "signal_list": window["signal_list"],
                "max_possible_k": window["max_possible_k"],
                "truth_window_label": label,
                "merge_boundary_count": len(positive_ids),
                "protected_boundary_count": len(protected_ids),
                "unresolved_evidence_count": sum(
                    value not in CLASSIFIED_LABELS for value in classifications
                ),
                "matched_breakpoint_ids": ",".join(
                    sorted(positive_ids | protected_ids, key=int)
                ),
            }
        )
    output.sort(
        key=lambda row: (
            str(row["source_round_dir"]),
            str(row["report_species"]),
            str(row["report_chromosome"]),
            int(str(row["original_start"] or 0)),
            str(row["window_id"]),
        )
    )
    return output


def scope_metrics(
    name: str,
    tiers: set[str],
    boundary_labels: list[dict[str, str | int]],
    window_labels: list[dict[str, str | int]],
    truth_counts: Counter,
) -> dict:
    selected = [row for row in boundary_labels if row["detector_priority"] in tiers]
    unique_truth: dict[str, str] = {}
    for row in selected:
        breakpoint_id = str(row["truth_breakpoint_id"])
        classification = str(row["truth_classification"])
        if breakpoint_id:
            unique_truth[breakpoint_id] = classification
    counts = Counter(unique_truth.values())
    classified = sum(counts[label] for label in CLASSIFIED_LABELS)
    protected = sum(counts[label] for label in PROTECTED_LABELS)
    selected_windows = [row for row in window_labels if row["priority_tier"] in tiers]
    mixed = sum(row["truth_window_label"] == "MIXED_WINDOW" for row in selected_windows)
    metrics = {
        "scope": name,
        "tiers": sorted(tiers),
        "candidate_pairwise_boundary_count": len(selected),
        "matched_unique_breakpoint_count": len(unique_truth),
        "over_fragmentation_detected": counts[POSITIVE_LABEL],
        "correct_truth_break_detected": counts["correct_truth_break"],
        "avoided_cactus_overmerge_detected": counts["avoided_cactus_overmerge"],
        "classified_precision": safe_ratio(counts[POSITIVE_LABEL], classified),
        "over_fragmentation_recall": safe_ratio(
            counts[POSITIVE_LABEL], truth_counts[POSITIVE_LABEL]
        ),
        "protected_boundary_contamination": safe_ratio(protected, classified),
        "selected_window_count": len(selected_windows),
        "mixed_window_count": mixed,
        "mixed_window_rate": safe_ratio(mixed, len(selected_windows)),
    }
    return metrics


def component_lengths(
    boundaries: Sequence[dict[str, str]], bridge_ids: set[str]
) -> tuple[list[int], list[int], int]:
    intervals: set[tuple[int, int]] = set()
    edges: list[tuple[tuple[int, int], tuple[int, int], str]] = []
    for boundary in boundaries:
        left = (as_int(boundary, "left_ref_start"), as_int(boundary, "left_ref_end"))
        right = (as_int(boundary, "right_ref_start"), as_int(boundary, "right_ref_end"))
        if left[1] <= left[0] or right[1] <= right[0]:
            continue
        intervals.add(left)
        intervals.add(right)
        edges.append((left, right, boundary["boundary_id"]))
    ordered = sorted(intervals)
    index = {interval: position for position, interval in enumerate(ordered)}
    parent = list(range(len(ordered)))

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    breaks_removed = 0
    for left, right, boundary_id in edges:
        if boundary_id in bridge_ids:
            union(index[left], index[right])
            breaks_removed += 1

    grouped: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for position, interval in enumerate(ordered):
        grouped[find(position)].append(interval)
    spans = [max(end for _, end in group) - min(start for start, _ in group) for group in grouped.values()]
    support = [sum(end - start for start, end in group) for group in grouped.values()]
    return spans, support, breaks_removed


def continuity_rows(
    windows: list[dict[str, str]],
    boundaries: list[dict[str, str]],
    boundary_labels: list[dict[str, str | int]],
    boundary_truth_bits: bytearray | None = None,
) -> tuple[list[dict], dict[str, dict[str, list[int]]]]:
    boundaries_by_window: dict[
        str, list[tuple[int, dict[str, str]]]
    ] = defaultdict(list)
    for boundary_index, row in enumerate(boundaries):
        boundaries_by_window[row["window_id"]].append((boundary_index, row))
    labels_by_boundary: dict[tuple[str, str], list[dict[str, str | int]]] = defaultdict(list)
    for row in boundary_labels:
        labels_by_boundary[(str(row["window_id"]), str(row["boundary_id"]))].append(row)

    output: list[dict] = []
    distributions: dict[str, dict[str, list[int]]] = {}
    for window in windows:
        window_id = window["window_id"]
        indexed_rows = sorted(
            boundaries_by_window.get(window_id, []),
            key=lambda item: (
                as_int(item[1], "left_ref_start"),
                as_int(item[1], "right_ref_start"),
                item[1]["boundary_id"],
            ),
        )
        rows = [item[1] for item in indexed_rows]
        if not rows:
            continue
        projected_bridges = {
            row["boundary_id"] for row in rows if as_bool(row, "detector_recommends_bridge")
        }
        oracle_bridges: set[str] = set()
        if boundary_truth_bits is not None:
            for boundary_index, row in indexed_rows:
                bits = boundary_truth_bits[boundary_index]
                if bits & 1 and not bits & 2:
                    oracle_bridges.add(row["boundary_id"])
        else:
            for row in rows:
                evidence = labels_by_boundary.get(
                    (window_id, row["boundary_id"]), []
                )
                labels = {
                    str(item["truth_classification"]) for item in evidence
                }
                if POSITIVE_LABEL in labels and not (labels & PROTECTED_LABELS):
                    oracle_bridges.add(row["boundary_id"])

        old_spans, old_support, _ = component_lengths(rows, set())
        projected_spans, projected_support, projected_removed = component_lengths(
            rows, projected_bridges
        )
        oracle_spans, oracle_support, oracle_removed = component_lengths(rows, oracle_bridges)
        distributions[window_id] = {
            "old": old_spans,
            "projected": projected_spans,
            "oracle": oracle_spans,
        }
        row: dict[str, str | int | float] = {
            "window_id": window_id,
            "source_round_dir": window.get("source_round_dir", ""),
            "priority_tier": window["priority_tier"],
            "report_species": window["report_species"],
            "report_chromosome": window["report_chromosome"],
            "original_start": window["original_start"],
            "original_end": window["original_end"],
            "old_fragment_count": len(old_spans),
            "projected_fragment_count": len(projected_spans),
            "oracle_fragment_count": len(oracle_spans),
            "breaks_removed_projected": projected_removed,
            "breaks_removed_oracle": oracle_removed,
            "old_max_span": max(old_spans, default=0),
            "projected_max_span": max(projected_spans, default=0),
            "oracle_max_span": max(oracle_spans, default=0),
            "old_local_n50": n50(old_spans),
            "projected_local_n50": n50(projected_spans),
            "oracle_local_n50": n50(oracle_spans),
            "n50_fold_change": fmt_ratio(
                safe_ratio(n50(projected_spans), n50(old_spans))
            ),
            "old_aligned_support_bases": sum(old_support),
            "projected_aligned_support_bases": sum(projected_support),
            "oracle_aligned_support_bases": sum(oracle_support),
        }
        for threshold in (100, 500, 1000):
            row[f"old_span_bases_ge_{threshold}"] = sum(
                length for length in old_spans if length >= threshold
            )
            row[f"projected_span_bases_ge_{threshold}"] = sum(
                length for length in projected_spans if length >= threshold
            )
            row[f"oracle_span_bases_ge_{threshold}"] = sum(
                length for length in oracle_spans if length >= threshold
            )
        output.append(row)
    return output, distributions


def non_overlapping_windows(
    windows: Sequence[dict[str, str | int]], allowed_tiers: set[str]
) -> list[dict[str, str | int]]:
    eligible = [
        row
        for row in windows
        if row["priority_tier"] in allowed_tiers
        and row["truth_window_label"] == "TRUE_POSITIVE_WINDOW"
        and row["coordinate_status"] == "resolved"
    ]
    eligible.sort(
        key=lambda row: (
            str(row["source_round_dir"]),
            str(row["report_species"]),
            str(row["report_chromosome"]),
            int(row["original_end"]),
            int(row["original_start"]),
            str(row["window_id"]),
        )
    )
    selected: list[dict[str, str | int]] = []
    last_end: dict[tuple[str, str, str], int] = {}
    for row in eligible:
        key = (
            str(row["source_round_dir"]),
            str(row["report_species"]),
            str(row["report_chromosome"]),
        )
        start = int(row["original_start"])
        end = int(row["original_end"])
        if start < last_end.get(key, -1):
            continue
        selected.append(row)
        last_end[key] = end
    return selected


def global_projection_rows(
    window_labels: list[dict[str, str | int]],
    distributions: dict[str, dict[str, list[int]]],
) -> list[dict]:
    scopes = (
        ("tier_a_conservative", {"A"}, "projected"),
        ("tier_a_b_extended", {"A", "B"}, "projected"),
        ("truth_oracle_upper_bound", {"A", "B", "C"}, "oracle"),
    )
    output: list[dict] = []
    for name, tiers, state in scopes:
        selected = non_overlapping_windows(window_labels, tiers)
        old_lengths: list[int] = []
        new_lengths: list[int] = []
        for window in selected:
            distribution = distributions.get(str(window["window_id"]))
            if not distribution:
                continue
            old_lengths.extend(distribution["old"])
            new_lengths.extend(distribution[state])
        output.append(
            {
                "scope": name,
                "projection_state": state,
                "selected_non_overlapping_window_count": len(selected),
                "old_fragment_count_in_selected_windows": len(old_lengths),
                "projected_fragment_count_in_selected_windows": len(new_lengths),
                "breaks_removed_in_selected_windows": len(old_lengths) - len(new_lengths),
                "old_window_scope_n50": n50(old_lengths),
                "projected_window_scope_n50": n50(new_lengths),
                "old_span_bases_ge_500": sum(value for value in old_lengths if value >= 500),
                "projected_span_bases_ge_500": sum(value for value in new_lengths if value >= 500),
                "old_span_bases_ge_1000": sum(value for value in old_lengths if value >= 1000),
                "projected_span_bases_ge_1000": sum(value for value in new_lengths if value >= 1000),
                "scope_note": "Metrics cover selected windows only; they are not whole-genome HAL N50.",
            }
        )
    return output


def feature_rows(
    boundary_labels: list[dict[str, str | int]],
) -> list[dict[str, str | int | float]]:
    grouped: dict[str, list[dict[str, str | int]]] = defaultdict(list)
    for row in boundary_labels:
        grouped[str(row["truth_classification"])].append(row)
    output: list[dict[str, str | int | float]] = []
    for label, rows in sorted(grouped.items()):
        lengths = [int(row["min_adjacent_block_length"]) for row in rows]
        gaps = [int(row["reference_gap"]) for row in rows]
        output.append(
            {
                "truth_classification": label,
                "count": len(rows),
                "min_block_p10": f"{percentile(lengths, 0.10):.3f}",
                "min_block_p50": f"{percentile(lengths, 0.50):.3f}",
                "min_block_p90": f"{percentile(lengths, 0.90):.3f}",
                "reference_gap_p10": f"{percentile(gaps, 0.10):.3f}",
                "reference_gap_p50": f"{percentile(gaps, 0.50):.3f}",
                "reference_gap_p90": f"{percentile(gaps, 0.90):.3f}",
            }
        )
    return output


def acceptance(metrics: dict[str, dict]) -> dict:
    tier_a = metrics["tier_a"]
    tier_ab = metrics["tier_a_b"]

    def at_least(value: float | None, threshold: float) -> bool:
        return value is not None and value >= threshold

    def at_most(value: float | None, threshold: float) -> bool:
        return value is not None and value <= threshold

    checks = {
        "tier_a_precision_ge_0_95": at_least(tier_a["classified_precision"], 0.95),
        "tier_a_recall_ge_0_40": at_least(tier_a["over_fragmentation_recall"], 0.40),
        "tier_a_protected_contamination_le_0_05": at_most(
            tier_a["protected_boundary_contamination"], 0.05
        ),
        "tier_a_b_precision_ge_0_75": at_least(tier_ab["classified_precision"], 0.75),
        "tier_a_b_recall_ge_0_75": at_least(tier_ab["over_fragmentation_recall"], 0.75),
        "tier_a_b_mixed_window_rate_le_0_05": at_most(
            tier_ab["mixed_window_rate"], 0.05
        ),
    }
    return {"checks": checks, "all_passed": all(checks.values())}


def markdown_report(
    report_dir: Path,
    breakpoint_path: Path,
    round_dirs: Sequence[Path],
    truth_counts: Counter,
    metrics: dict[str, dict],
    acceptance_result: dict,
    window_labels: list[dict[str, str | int]],
    global_rows: list[dict],
) -> str:
    window_counts = Counter(str(row["truth_window_label"]) for row in window_labels)
    lines = [
        "# RaMAx 第一阶段窗口识别评估",
        "",
        "## 输入",
        "",
        f"- 窗口报告：`{report_dir}`",
        f"- 断点 truth 分类：`{breakpoint_path}`",
        f"- 评估轮次：{', '.join(path.name for path in round_dirs)}",
        f"- truth 断点总数：{sum(truth_counts.values()):,}",
        f"- truth-confirmed over-fragmentation：{truth_counts[POSITIVE_LABEL]:,}",
        "",
        "## 检测结果",
        "",
        f"- TRUE_POSITIVE_WINDOW：{window_counts['TRUE_POSITIVE_WINDOW']:,}",
        f"- FALSE_POSITIVE_WINDOW：{window_counts['FALSE_POSITIVE_WINDOW']:,}",
        f"- MIXED_WINDOW：{window_counts['MIXED_WINDOW']:,}",
        f"- UNRESOLVED_WINDOW：{window_counts['UNRESOLVED_WINDOW']:,}",
        "",
        "## Tier 指标",
        "",
        "| 范围 | classified precision | over-fragmentation recall | 保护边界污染率 | mixed-window rate |",
        "|---|---:|---:|---:|---:|",
    ]
    for key in ("tier_a", "tier_a_b"):
        item = metrics[key]
        lines.append(
            "| {scope} | {precision} | {recall} | {contamination} | {mixed} |".format(
                scope=item["scope"],
                precision=fmt_ratio(item["classified_precision"]),
                recall=fmt_ratio(item["over_fragmentation_recall"]),
                contamination=fmt_ratio(item["protected_boundary_contamination"]),
                mixed=fmt_ratio(item["mixed_window_rate"]),
            )
        )
    lines.extend(["", "## 验收", ""])
    for name, passed in acceptance_result["checks"].items():
        lines.append(f"- {'PASS' if passed else 'FAIL'}：{name}")
    lines.extend(
        [
            "",
            "## 连续性投影",
            "",
            "以下 N50 只覆盖互不重叠的候选窗口，不是全基因组 HAL N50。",
            "",
            "| 范围 | 窗口数 | 移除断点 | OLD N50 | projected N50 | >=500 bp 增量 | >=1000 bp 增量 |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in global_rows:
        lines.append(
            f"| {row['scope']} | {row['selected_non_overlapping_window_count']} | "
            f"{row['breaks_removed_in_selected_windows']} | "
            f"{row['old_window_scope_n50']} | {row['projected_window_scope_n50']} | "
            f"{row['projected_span_bases_ge_500'] - row['old_span_bases_ge_500']} | "
            f"{row['projected_span_bases_ge_1000'] - row['old_span_bases_ge_1000']} |"
        )
    lines.extend(
        [
            "",
            "## 解释边界",
            "",
            "- 本报告复用已有 breakpoint classification；Cactus 只通过该分类提供辅助证据，不被当作 truth。",
            "- DETECTOR_PROJECTION 只连接检测器建议的边界；TRUTH_ORACLE 只连接 truth-confirmed over-fragmentation。",
            "- 第一阶段不删除 Block、不运行 minipoa，也不修改 HAL/MAF。",
            "- 精确错误 Block 与 missing-homology 碱基级枚举属于下一层局部 MAF/truth 对照。",
            "",
        ]
    )
    return "\n".join(lines)


def markdown_report_direct(
    report_dir: Path,
    truth_maf: Path,
    breakpoint_path: Path | None,
    round_dirs: Sequence[Path],
    classification_counts: Counter,
    metrics: dict[str, dict],
    acceptance_result: dict,
    window_labels: Sequence[dict[str, str | int]],
    global_rows: Sequence[dict],
    candidate_stats: dict,
) -> str:
    window_counts = Counter(
        str(row["truth_window_label"]) for row in window_labels
    )
    lines = [
        "# RaMAx 第一阶段窗口直接 Ground Truth 评估",
        "",
        "## 输入",
        "",
        f"- 窗口报告：`{report_dir}`",
        f"- Alignathon truth MAF：`{truth_maf}`",
        f"- 旧断点分类辅助输入：`{breakpoint_path or '未使用'}`",
        f"- 评估轮次：{', '.join(path.name for path in round_dirs)}",
        f"- truth MAF 扫描 block：{candidate_stats.get('maf_block_count', 0):,}",
        f"- 保留 truth 锚点候选：{candidate_stats.get('retained_candidate_count', 0):,}",
        "",
        "## 当前候选边界的直接 Truth 分类",
        "",
        f"- MERGE_BOUNDARY：{classification_counts[POSITIVE_LABEL]:,}",
        f"- KEEP_BREAK：{classification_counts['truth_supported_break']:,}",
        f"- insufficient truth support：{classification_counts['insufficient_truth_support']:,}",
        f"- copy ambiguous：{classification_counts['copy_ambiguous']:,}",
        "",
        "## 窗口分类",
        "",
        f"- TRUE_POSITIVE_WINDOW：{window_counts['TRUE_POSITIVE_WINDOW']:,}",
        f"- FALSE_POSITIVE_WINDOW：{window_counts['FALSE_POSITIVE_WINDOW']:,}",
        f"- MIXED_WINDOW：{window_counts['MIXED_WINDOW']:,}",
        f"- UNRESOLVED_WINDOW：{window_counts['UNRESOLVED_WINDOW']:,}",
        "",
        "## Tier 指标",
        "",
        "| 范围 | classified precision | direct truth coverage | 保护边界污染率 | mixed-window rate | recall |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for key in ("tier_a", "tier_a_b"):
        item = metrics[key]
        lines.append(
            "| {scope} | {precision} | {coverage} | {contamination} | "
            "{mixed} | NA |".format(
                scope=item["scope"],
                precision=fmt_ratio(item["classified_precision"]),
                coverage=fmt_ratio(item["direct_truth_coverage"]),
                contamination=fmt_ratio(
                    item["protected_boundary_contamination"]
                ),
                mixed=fmt_ratio(item["mixed_window_rate"]),
            )
        )
    lines.extend(["", "## 验收", ""])
    for name, passed in acceptance_result["checks"].items():
        state = "NOT_EVALUATED" if passed is None else ("PASS" if passed else "FAIL")
        lines.append(f"- {state}：{name}")
    lines.extend(
        [
            "",
            "## 连续性投影",
            "",
            "以下 N50 只覆盖互不重叠的候选窗口，不是全基因组 HAL N50。",
            "",
            "| 范围 | 窗口数 | 移除断点 | OLD N50 | projected N50 | >=500 bp 增量 | >=1000 bp 增量 |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in global_rows:
        lines.append(
            f"| {row['scope']} | {row['selected_non_overlapping_window_count']} | "
            f"{row['breaks_removed_in_selected_windows']} | "
            f"{row['old_window_scope_n50']} | {row['projected_window_scope_n50']} | "
            f"{row['projected_span_bases_ge_500'] - row['old_span_bases_ge_500']} | "
            f"{row['projected_span_bases_ge_1000'] - row['old_span_bases_ge_1000']} |"
        )
    lines.extend(
        [
            "",
            "## 解释边界",
            "",
            "- 本报告直接查询当前窗口在 truth MAF 中的左右锚点，不要求与旧 RaMAx 断点坐标相同。",
            "- classified precision 只描述已经进入候选集合的边界；由于报告不包含所有未入选图边界，正式 recall 不能由本次输出计算。",
            "- truth 支持连接表示应送入局部重比对，不表示可以不经 minipoa 直接修改图。",
            "- 第一阶段不删除 Block、不运行 minipoa，也不修改 HAL/MAF。",
            "",
        ]
    )
    return "\n".join(lines)


def evaluate_direct_maf(
    args: argparse.Namespace,
    report_dir: Path,
    output_dir: Path,
    manifest: dict,
    round_dirs: Sequence[Path],
) -> int:
    truth_maf = args.truth_maf.resolve()
    if not truth_maf.is_file():
        raise FileNotFoundError(truth_maf)
    breakpoint_path = (
        args.breakpoint_classifications.resolve()
        if args.breakpoint_classifications is not None
        else None
    )
    if breakpoint_path is not None and not breakpoint_path.is_file():
        raise FileNotFoundError(breakpoint_path)

    log(
        f"loading direct-mode reports from {', '.join(path.name for path in round_dirs)}",
        args.quiet,
    )
    windows, boundaries, boundary_sources = load_direct_reports(round_dirs)
    input_species = tuple(manifest.get("input_genomes", ()))
    if len(input_species) < 2:
        raise ValueError("run_manifest.json does not contain at least two input genomes")
    log(
        f"loaded {len(windows):,} windows and {len(boundaries):,} boundaries; "
        "boundary-genome evidence will be streamed",
        args.quiet,
    )

    cache_dir = (
        args.truth_cache_dir.resolve()
        if args.truth_cache_dir is not None
        else output_dir / "truth_maf_cache"
    )
    candidate_path, candidate_stats = ensure_direct_truth_candidates(
        truth_maf=truth_maf,
        boundaries=boundaries,
        boundary_sources=boundary_sources,
        input_species=input_species,
        cache_dir=cache_dir,
        flank=args.truth_flank_bp,
        max_candidates=args.truth_max_candidates,
        threads=args.threads,
        sort_memory=args.sort_memory,
        validate_maf=args.validate_maf,
        max_maf_blocks=args.max_maf_blocks,
        keep_shards=args.keep_truth_shards,
        quiet=args.quiet,
    )

    legacy_matcher: TruthMatcher | None = None
    replay_truth_rows: list[dict[str, str]] = []
    if breakpoint_path is not None:
        log("loading optional legacy breakpoint replay labels", args.quiet)
        replay_truth_rows = read_tsv(breakpoint_path)
        legacy_matcher = TruthMatcher(
            replay_truth_rows, args.coordinate_tolerance
        )

    window_index: dict[tuple[str, str], int] = {}
    for index, window in enumerate(windows):
        key = (window.get("source_round_dir", ""), window["window_id"])
        if key in window_index:
            raise ValueError(f"Duplicate window key in report: {key}")
        window_index[key] = index

    boundary_bits = bytearray(len(boundaries))
    window_merge_pairs = array("Q", [0]) * len(windows)
    window_keep_pairs = array("Q", [0]) * len(windows)
    window_unresolved_pairs = array("Q", [0]) * len(windows)
    window_ambiguous_pairs = array("Q", [0]) * len(windows)
    species_index = {species: index for index, species in enumerate(input_species)}
    star_slot_count = len(windows) * len(input_species)
    window_target_bits = bytearray(star_slot_count)
    window_target_merge_counts = array("I", [0]) * star_slot_count
    tier_counts: dict[str, Counter] = defaultdict(Counter)
    classification_counts: Counter = Counter()
    boundary_label_counts: Counter = Counter()
    bridge_relation_counts: Counter = Counter()
    replay_method_counts: Counter = Counter()
    feature_values: dict[str, tuple[array, array]] = {}

    candidate_groups = iter(iter_direct_candidate_groups(candidate_path))
    current_candidates = next(candidate_groups, None)
    evidence_groups = iter(iter_boundary_evidence_groups(round_dirs))
    current_evidence = next(evidence_groups, None)
    evidence_count = 0
    output_dir.mkdir(parents=True, exist_ok=True)
    boundary_path = output_dir / "boundary_truth_labels.tsv"
    boundary_partial = boundary_path.with_name(boundary_path.name + ".partial")
    with boundary_partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=DIRECT_BOUNDARY_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for boundary_index, boundary in enumerate(boundaries):
            while (
                current_candidates is not None
                and current_candidates[0] < boundary_index
            ):
                current_candidates = next(candidate_groups, None)
            if (
                current_candidates is not None
                and current_candidates[0] == boundary_index
            ):
                candidates = current_candidates[1]
                truncated = current_candidates[2]
                current_candidates = next(candidate_groups, None)
            else:
                candidates = {}
                truncated = {}

            expected_key = (
                boundary.get("source_round_dir", ""),
                boundary["window_id"],
                boundary["boundary_id"],
            )
            if current_evidence is None or current_evidence[0] != expected_key:
                observed = None if current_evidence is None else current_evidence[0]
                raise ValueError(
                    "window_boundaries.tsv and boundary_genomes.tsv are not in "
                    f"the same grouped order: expected {expected_key}, observed {observed}"
                )
            evidence_rows = current_evidence[1]
            current_evidence = next(evidence_groups, None)
            window_key = (
                boundary.get("source_round_dir", ""),
                boundary["window_id"],
            )
            window_position = window_index.get(window_key)
            if window_position is None:
                raise ValueError(f"Boundary has no parent window: {expected_key}")

            for evidence in evidence_rows:
                evidence_count += 1
                target_species = evidence["target_species"]
                result = classify_direct_truth(
                    target_species,
                    candidates,
                    truncated,
                    args.truth_max_gap_bp,
                )
                classification_counts[result.classification] += 1
                boundary_label_counts[result.boundary_label] += 1
                tier_counts[boundary["detector_priority"]][
                    result.classification
                ] += 1
                if result.boundary_label == "MERGE_BOUNDARY":
                    boundary_bits[boundary_index] |= 1
                    window_merge_pairs[window_position] += 1
                    target_position = species_index.get(target_species)
                    if target_position is not None:
                        star_slot = (
                            window_position * len(input_species) + target_position
                        )
                        window_target_bits[star_slot] |= 1
                        window_target_merge_counts[star_slot] += 1
                elif result.boundary_label == "KEEP_BREAK":
                    boundary_bits[boundary_index] |= 2
                    window_keep_pairs[window_position] += 1
                    target_position = species_index.get(target_species)
                    if target_position is not None:
                        window_target_bits[
                            window_position * len(input_species) + target_position
                        ] |= 2
                else:
                    boundary_bits[boundary_index] |= 4
                    window_unresolved_pairs[window_position] += 1
                    target_position = species_index.get(target_species)
                    if target_position is not None:
                        window_target_bits[
                            window_position * len(input_species) + target_position
                        ] |= 4
                    if result.ambiguous:
                        boundary_bits[boundary_index] |= 8
                        window_ambiguous_pairs[window_position] += 1

                values = feature_values.get(result.classification)
                if values is None:
                    values = (array("q"), array("q"))
                    feature_values[result.classification] = values
                values[0].append(as_int(boundary, "min_adjacent_block_length"))
                values[1].append(as_int(boundary, "reference_gap"))

                replay = (
                    legacy_matcher.match(boundary, evidence)
                    if legacy_matcher is not None
                    else MatchResult(None, "not_requested", None, "")
                )
                replay_truth = replay.truth or {}
                replay_method_counts[replay.method] += 1
                bridge = result.bridge
                bridge_left = bridge[0] if bridge is not None else None
                bridge_right = bridge[1] if bridge is not None else None
                left_prediction_match = any(
                    direct_candidate_matches_prediction(
                        candidate,
                        evidence,
                        "left",
                        args.prediction_coordinate_tolerance,
                    )
                    for candidate in result.left_candidates
                )
                right_prediction_match = any(
                    direct_candidate_matches_prediction(
                        candidate,
                        evidence,
                        "right",
                        args.prediction_coordinate_tolerance,
                    )
                    for candidate in result.right_candidates
                )
                bridge_prediction_match = bool(
                    bridge_left is not None
                    and bridge_right is not None
                    and direct_candidate_matches_prediction(
                        bridge_left,
                        evidence,
                        "left",
                        args.prediction_coordinate_tolerance,
                    )
                    and direct_candidate_matches_prediction(
                        bridge_right,
                        evidence,
                        "right",
                        args.prediction_coordinate_tolerance,
                    )
                )
                if bridge is None:
                    bridge_relation = ""
                elif bridge_prediction_match:
                    bridge_relation = "EXISTING_FLANKS_CONTINUOUS"
                elif left_prediction_match or right_prediction_match:
                    bridge_relation = "PARTIAL_ALTERNATIVE_MAPPING"
                else:
                    bridge_relation = "ALTERNATIVE_TARGET_MAPPING"
                if bridge_relation:
                    bridge_relation_counts[bridge_relation] += 1
                writer.writerow(
                    {
                        "window_id": boundary["window_id"],
                        "boundary_id": boundary["boundary_id"],
                        "source_round_dir": boundary.get(
                            "source_round_dir", ""
                        ),
                        "detector_priority": boundary["detector_priority"],
                        "detector_recommends_bridge": boundary[
                            "detector_recommends_bridge"
                        ],
                        "signal_type": boundary["signal_type"],
                        "reference_species": boundary["reference_species"],
                        "target_species": target_species,
                        "reference_chromosome": boundary[
                            "reference_chromosome"
                        ],
                        "left_ref_start": boundary["left_ref_start"],
                        "left_ref_end": boundary["left_ref_end"],
                        "right_ref_start": boundary["right_ref_start"],
                        "right_ref_end": boundary["right_ref_end"],
                        "reference_gap": boundary["reference_gap"],
                        "min_adjacent_block_length": boundary[
                            "min_adjacent_block_length"
                        ],
                        "truth_source": "direct_maf",
                        "truth_classification": result.classification,
                        "boundary_truth_label": result.boundary_label,
                        "truth_break_reason": result.reason,
                        "truth_left_candidates": len(result.left_candidates),
                        "truth_right_candidates": len(result.right_candidates),
                        "truth_left_truncated": result.left_truncated,
                        "truth_right_truncated": result.right_truncated,
                        "truth_bridge": int(bridge is not None),
                        "truth_bridge_target_chromosome": (
                            bridge_left.target_chrom if bridge_left else ""
                        ),
                        "truth_bridge_strand": (
                            bridge_left.strand if bridge_left else ""
                        ),
                        "truth_left_ref_coord": (
                            bridge_left.ref_coord if bridge_left else ""
                        ),
                        "truth_left_target_coord": (
                            bridge_left.target_coord if bridge_left else ""
                        ),
                        "truth_right_ref_coord": (
                            bridge_right.ref_coord if bridge_right else ""
                        ),
                        "truth_right_target_coord": (
                            bridge_right.target_coord if bridge_right else ""
                        ),
                        "truth_ref_gap": bridge[2] if bridge else "",
                        "truth_target_gap": bridge[3] if bridge else "",
                        "truth_compatible_pair_count": (
                            result.compatible_pair_count
                        ),
                        "truth_mapping_ambiguous": int(result.ambiguous),
                        "truth_left_matches_prediction": int(
                            left_prediction_match
                        ),
                        "truth_right_matches_prediction": int(
                            right_prediction_match
                        ),
                        "truth_bridge_matches_prediction": int(
                            bridge_prediction_match
                        ),
                        "truth_bridge_relation": bridge_relation,
                        "replay_match_method": replay.method,
                        "replay_match_distance": (
                            "" if replay.distance is None else replay.distance
                        ),
                        "replay_unmatched_reason": replay.reason,
                        "replay_truth_breakpoint_id": replay_truth.get(
                            "breakpoint_id", ""
                        ),
                        "replay_truth_classification": replay_truth.get(
                            "classification", ""
                        ),
                        "cactus_bridges": replay_truth.get(
                            "cactus_bridges", ""
                        ),
                    }
                )
                if evidence_count % 100000 == 0:
                    log(
                        f"directly classified {evidence_count:,} "
                        "boundary-genome rows",
                        args.quiet,
                    )
    if current_evidence is not None or next(evidence_groups, None) is not None:
        raise ValueError("boundary_genomes.tsv contains rows without boundaries")
    if current_candidates is not None or next(candidate_groups, None) is not None:
        raise ValueError("truth candidate cache contains an invalid boundary index")
    boundary_partial.replace(boundary_path)

    merge_boundaries = array("Q", [0]) * len(windows)
    keep_boundaries = array("Q", [0]) * len(windows)
    unresolved_boundaries = array("Q", [0]) * len(windows)
    for boundary_index, boundary in enumerate(boundaries):
        window_position = window_index[
            (boundary.get("source_round_dir", ""), boundary["window_id"])
        ]
        bits = boundary_bits[boundary_index]
        if bits & 1:
            merge_boundaries[window_position] += 1
        if bits & 2:
            keep_boundaries[window_position] += 1
        if not bits or bits & 4:
            unresolved_boundaries[window_position] += 1

    window_labels: list[dict[str, str | int]] = []
    k_comparison_counts: Counter = Counter()
    for index, window in enumerate(windows):
        has_merge = merge_boundaries[index] > 0
        has_keep = keep_boundaries[index] > 0
        if has_merge and has_keep:
            label = "MIXED_WINDOW"
        elif has_merge:
            label = "TRUE_POSITIVE_WINDOW"
        elif has_keep:
            label = "FALSE_POSITIVE_WINDOW"
        else:
            label = "UNRESOLVED_WINDOW"
        internal_boundary_count = as_int(window, "internal_boundary_count")
        reference_species = window["report_species"]
        compatible_species = [reference_species]
        excluded: list[str] = []
        for species_position, species in enumerate(input_species):
            if species == reference_species:
                continue
            star_slot = index * len(input_species) + species_position
            bits = window_target_bits[star_slot]
            merge_count = window_target_merge_counts[star_slot]
            if (
                internal_boundary_count > 0
                and bits == 1
                and merge_count == internal_boundary_count
            ):
                compatible_species.append(species)
            elif bits & 2:
                excluded.append(f"{species}:truth_break")
            elif bits & 4 or merge_count < internal_boundary_count:
                excluded.append(f"{species}:truth_unresolved")
            else:
                excluded.append(f"{species}:no_truth_evidence")
        truth_reference_star_k = (
            len(compatible_species) if internal_boundary_count > 0 else 0
        )
        detector_k = as_int(window, "max_possible_k")
        if truth_reference_star_k == 0:
            k_comparison = "K_UNRESOLVED"
        elif detector_k == truth_reference_star_k:
            k_comparison = "K_EXACT"
        elif detector_k > truth_reference_star_k:
            k_comparison = "K_OVERESTIMATED"
        else:
            k_comparison = "K_UNDERESTIMATED"
        k_comparison_counts[k_comparison] += 1
        window_labels.append(
            {
                "window_id": window["window_id"],
                "source_round_dir": window.get("source_round_dir", ""),
                "priority_tier": window["priority_tier"],
                "report_species": window["report_species"],
                "report_chromosome": window["report_chromosome"],
                "original_start": window["original_start"],
                "original_end": window["original_end"],
                "coordinate_status": window["coordinate_status"],
                "signal_list": window["signal_list"],
                "max_possible_k": window["max_possible_k"],
                "truth_window_label": label,
                "merge_boundary_count": merge_boundaries[index],
                "protected_boundary_count": keep_boundaries[index],
                "unresolved_boundary_count": unresolved_boundaries[index],
                "merge_pair_count": window_merge_pairs[index],
                "protected_pair_count": window_keep_pairs[index],
                "unresolved_evidence_count": window_unresolved_pairs[index],
                "copy_ambiguous_pair_count": window_ambiguous_pairs[index],
                "truth_reference_star_k": truth_reference_star_k,
                "truth_reference_star_genomes": ",".join(compatible_species),
                "truth_reference_star_excluded": ";".join(excluded),
                "k_comparison": k_comparison,
                "matched_breakpoint_ids": "",
                "truth_evaluation_mode": "direct_maf",
            }
        )

    metrics = {
        "tier_a": direct_scope_metrics(
            "Tier A", {"A"}, tier_counts, window_labels
        ),
        "tier_a_b": direct_scope_metrics(
            "Tier A+B", {"A", "B"}, tier_counts, window_labels
        ),
    }
    acceptance_result = direct_acceptance(metrics)
    log("calculating direct-truth continuity projections", args.quiet)
    continuity, distributions = continuity_rows(
        windows, boundaries, [], boundary_bits
    )
    global_rows = global_projection_rows(window_labels, distributions)
    features = feature_rows_from_arrays(feature_values)

    window_fields = [
        "window_id",
        "source_round_dir",
        "priority_tier",
        "report_species",
        "report_chromosome",
        "original_start",
        "original_end",
        "coordinate_status",
        "signal_list",
        "max_possible_k",
        "truth_window_label",
        "merge_boundary_count",
        "protected_boundary_count",
        "unresolved_boundary_count",
        "merge_pair_count",
        "protected_pair_count",
        "unresolved_evidence_count",
        "copy_ambiguous_pair_count",
        "truth_reference_star_k",
        "truth_reference_star_genomes",
        "truth_reference_star_excluded",
        "k_comparison",
        "matched_breakpoint_ids",
        "truth_evaluation_mode",
    ]
    continuity_fields = list(continuity[0].keys()) if continuity else ["window_id"]
    global_fields = list(global_rows[0].keys())
    feature_fields = (
        list(features[0].keys())
        if features
        else ["truth_classification", "count"]
    )
    write_tsv(output_dir / "window_truth_labels.tsv", window_fields, window_labels)
    write_tsv(output_dir / "continuity_gain.tsv", continuity_fields, continuity)
    write_tsv(
        output_dir / "global_continuity_projection.tsv",
        global_fields,
        global_rows,
    )
    write_tsv(
        output_dir / "feature_distributions.tsv", feature_fields, features
    )

    input_manifest = {
        "window_report_dir": str(report_dir),
        "window_report_manifest": manifest,
        "selected_rounds": [path.name for path in round_dirs],
        "truth_maf": str(truth_maf),
        "truth_maf_fingerprint": file_fingerprint(truth_maf),
        "truth_candidate_cache": str(candidate_path),
        "truth_candidate_cache_stats": candidate_stats,
        "breakpoint_classifications": (
            str(breakpoint_path) if breakpoint_path is not None else None
        ),
        "truth_flank_bp": args.truth_flank_bp,
        "truth_max_gap_bp": args.truth_max_gap_bp,
        "truth_max_candidates": args.truth_max_candidates,
        "prediction_coordinate_tolerance": (
            args.prediction_coordinate_tolerance
        ),
        "threads": args.threads,
        "window_count": len(windows),
        "boundary_count": len(boundaries),
        "boundary_genome_count": evidence_count,
        "window_block_row_count": None,
    }
    write_json(output_dir / "input_manifest.json", input_manifest)
    summary = {
        "status": "complete",
        "evaluation_version": DIRECT_TRUTH_VERSION,
        "input": input_manifest,
        "truth_classification_counts": dict(
            sorted(classification_counts.items())
        ),
        "boundary_truth_label_counts": dict(
            sorted(boundary_label_counts.items())
        ),
        "truth_bridge_relation_counts": dict(
            sorted(bridge_relation_counts.items())
        ),
        "k_comparison_counts": dict(sorted(k_comparison_counts.items())),
        "replay_match_method_counts": dict(
            sorted(replay_method_counts.items())
        ),
        "window_truth_label_counts": dict(
            sorted(
                Counter(
                    str(row["truth_window_label"]) for row in window_labels
                ).items()
            )
        ),
        "metrics": metrics,
        "acceptance": acceptance_result,
        "global_continuity_projection": global_rows,
        "limitations": [
            "Recall is not evaluable because the detector report does not contain all non-candidate graph boundaries.",
            "Direct MAF labels test boundary continuity; full per-base false-Block and missing-edge enumeration remains a later layer.",
            "Continuity N50 is calculated within selected windows, not over the whole HAL.",
            "A MERGE_BOUNDARY is a local-realignment candidate, not authorization to modify the graph without minipoa validation.",
        ],
    }
    write_json(output_dir / "summary.json", summary)
    report = markdown_report_direct(
        report_dir,
        truth_maf,
        breakpoint_path,
        round_dirs,
        classification_counts,
        metrics,
        acceptance_result,
        window_labels,
        global_rows,
        candidate_stats,
    )
    report_path = output_dir / "window_detection_evaluation.md"
    report_partial = report_path.with_name(report_path.name + ".partial")
    report_partial.write_text(report, encoding="utf-8")
    report_partial.replace(report_path)
    log(f"direct truth evaluation complete: {output_dir}", args.quiet)
    return 0


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate RaMAx phase-1 windows directly against an Alignathon "
            "truth MAF, or replay legacy breakpoint labels."
        )
    )
    parser.add_argument("--window-report-dir", type=Path, required=True)
    parser.add_argument(
        "--truth-maf",
        type=Path,
        help=(
            "Alignathon truth MAF for direct boundary evaluation "
            "(recommended)."
        ),
    )
    parser.add_argument(
        "--breakpoint-classifications",
        type=Path,
        help=(
            "Optional legacy breakpoint labels. In direct MAF mode they are "
            "reported as auxiliary replay evidence; without --truth-maf they "
            "select legacy replay mode."
        ),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--round",
        default="latest",
        help="latest, all, a numeric round id, or round_### (default: latest)",
    )
    parser.add_argument(
        "--coordinate-tolerance",
        type=int,
        default=5,
        help="Maximum endpoint difference for fallback matching (default: 5 bp)",
    )
    parser.add_argument(
        "--truth-flank-bp",
        type=int,
        default=100,
        help="Reference bases searched on each side of a boundary (default: 100).",
    )
    parser.add_argument(
        "--truth-max-gap-bp",
        type=int,
        default=1000,
        help="Maximum truth reference/target gap considered continuous (default: 1000).",
    )
    parser.add_argument(
        "--truth-max-candidates",
        type=int,
        default=8,
        help="Maximum truth mappings retained per boundary side/species (default: 8).",
    )
    parser.add_argument(
        "--prediction-coordinate-tolerance",
        type=int,
        default=5,
        help="Tolerance for checking whether truth supports current target blocks (default: 5).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=max(1, min(8, os.cpu_count() or 1)),
        help="Parallel MAF scan/sort workers (default: min(8, available CPUs)).",
    )
    parser.add_argument(
        "--truth-cache-dir",
        type=Path,
        help="Reusable direct-truth cache directory (default: OUTPUT/truth_maf_cache).",
    )
    parser.add_argument(
        "--sort-memory",
        default="25%",
        help="GNU sort memory budget used for truth candidates (default: 25%%).",
    )
    parser.add_argument(
        "--validate-maf",
        action="store_true",
        help="Validate every parsed MAF row; safer but slower.",
    )
    parser.add_argument(
        "--max-maf-blocks",
        type=int,
        default=0,
        help="Development-only limit on truth MAF blocks; 0 means all.",
    )
    parser.add_argument(
        "--keep-truth-shards",
        action="store_true",
        help="Keep raw/sorted truth-candidate shards after cache completion.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite evaluator-owned files in an existing output directory.",
    )
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)
    if args.truth_maf is None and args.breakpoint_classifications is None:
        parser.error(
            "one of --truth-maf or --breakpoint-classifications is required"
        )
    if args.coordinate_tolerance < 0:
        parser.error("--coordinate-tolerance must be >= 0")
    if args.truth_flank_bp <= 0:
        parser.error("--truth-flank-bp must be > 0")
    if args.truth_max_gap_bp < 0:
        parser.error("--truth-max-gap-bp must be >= 0")
    if args.truth_max_candidates <= 0:
        parser.error("--truth-max-candidates must be > 0")
    if args.prediction_coordinate_tolerance < 0:
        parser.error("--prediction-coordinate-tolerance must be >= 0")
    if args.threads <= 0:
        parser.error("--threads must be > 0")
    if args.max_maf_blocks < 0:
        parser.error("--max-maf-blocks must be >= 0")
    return args


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    report_dir = args.window_report_dir.resolve()
    output_dir = args.output_dir.resolve()

    manifest_path = report_dir / "run_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing run manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "complete":
        raise ValueError(f"Window report is not complete: {manifest_path}")
    owned_outputs = {
        "input_manifest.json",
        "boundary_truth_labels.tsv",
        "window_truth_labels.tsv",
        "continuity_gain.tsv",
        "global_continuity_projection.tsv",
        "feature_distributions.tsv",
        "summary.json",
        "window_detection_evaluation.md",
    }
    existing_owned = [name for name in owned_outputs if (output_dir / name).exists()]
    if existing_owned and not args.force:
        raise FileExistsError(
            f"Output directory already contains evaluator files: {output_dir}; "
            "choose a new directory or pass --force"
        )
    output_dir.mkdir(parents=True, exist_ok=True)

    round_dirs = discover_round_dirs(report_dir, args.round)
    if args.truth_maf is not None:
        return evaluate_direct_maf(
            args,
            report_dir,
            output_dir,
            manifest,
            round_dirs,
        )

    if args.breakpoint_classifications is None:
        raise ValueError("Legacy replay mode requires --breakpoint-classifications")
    breakpoint_path = args.breakpoint_classifications.resolve()
    if not breakpoint_path.is_file():
        raise FileNotFoundError(breakpoint_path)
    log(f"loading reports from {', '.join(path.name for path in round_dirs)}", args.quiet)
    windows, boundaries, boundary_genomes, blocks = load_reports(round_dirs)
    log(
        f"loaded {len(windows)} windows, {len(boundaries)} boundaries, "
        f"{len(boundary_genomes)} boundary-genome rows",
        args.quiet,
    )

    log("loading Alignathon breakpoint classifications", args.quiet)
    truth_rows = read_tsv(breakpoint_path)
    truth_counts = Counter(row["classification"] for row in truth_rows)
    matcher = TruthMatcher(truth_rows, args.coordinate_tolerance)

    log("matching candidate boundaries to truth labels", args.quiet)
    boundary_labels = evaluate_boundaries(
        boundaries, boundary_genomes, matcher, args.quiet
    )
    window_labels = classify_windows(windows, boundary_labels)

    metrics = {
        "tier_a": scope_metrics(
            "Tier A", {"A"}, boundary_labels, window_labels, truth_counts
        ),
        "tier_a_b": scope_metrics(
            "Tier A+B", {"A", "B"}, boundary_labels, window_labels, truth_counts
        ),
    }
    acceptance_result = acceptance(metrics)

    log("calculating local and non-overlapping continuity projections", args.quiet)
    continuity, distributions = continuity_rows(windows, boundaries, boundary_labels)
    global_rows = global_projection_rows(window_labels, distributions)
    features = feature_rows(boundary_labels)

    boundary_fields = [
        "window_id", "boundary_id", "source_round_dir", "detector_priority",
        "detector_recommends_bridge", "signal_type", "reference_species",
        "target_species", "reference_chromosome", "left_ref_start",
        "left_ref_end", "right_ref_start", "right_ref_end", "reference_gap",
        "min_adjacent_block_length", "match_method", "match_distance",
        "unmatched_reason", "truth_breakpoint_id", "truth_classification",
        "boundary_truth_label", "truth_break_reason", "cactus_bridges",
    ]
    window_fields = [
        "window_id", "source_round_dir", "priority_tier", "report_species",
        "report_chromosome", "original_start", "original_end", "coordinate_status",
        "signal_list", "max_possible_k", "truth_window_label",
        "merge_boundary_count", "protected_boundary_count",
        "unresolved_evidence_count", "matched_breakpoint_ids",
    ]
    continuity_fields = list(continuity[0].keys()) if continuity else ["window_id"]
    global_fields = list(global_rows[0].keys())
    feature_fields = list(features[0].keys()) if features else ["truth_classification", "count"]

    write_tsv(output_dir / "boundary_truth_labels.tsv", boundary_fields, boundary_labels)
    write_tsv(output_dir / "window_truth_labels.tsv", window_fields, window_labels)
    write_tsv(output_dir / "continuity_gain.tsv", continuity_fields, continuity)
    write_tsv(output_dir / "global_continuity_projection.tsv", global_fields, global_rows)
    write_tsv(output_dir / "feature_distributions.tsv", feature_fields, features)

    input_manifest = {
        "window_report_dir": str(report_dir),
        "window_report_manifest": manifest,
        "selected_rounds": [path.name for path in round_dirs],
        "breakpoint_classifications": str(breakpoint_path),
        "coordinate_tolerance": args.coordinate_tolerance,
        "window_count": len(windows),
        "boundary_count": len(boundaries),
        "boundary_genome_count": len(boundary_genomes),
        "window_block_row_count": len(blocks),
        "truth_breakpoint_count": len(truth_rows),
    }
    write_json(output_dir / "input_manifest.json", input_manifest)

    matched_methods = Counter(str(row["match_method"]) for row in boundary_labels)
    summary = {
        "status": "complete",
        "evaluation_version": "phase1-breakpoint-replay-v1",
        "input": input_manifest,
        "truth_classification_counts": dict(sorted(truth_counts.items())),
        "match_method_counts": dict(sorted(matched_methods.items())),
        "window_truth_label_counts": dict(
            sorted(Counter(str(row["truth_window_label"]) for row in window_labels).items())
        ),
        "metrics": metrics,
        "acceptance": acceptance_result,
        "global_continuity_projection": global_rows,
        "limitations": [
            "Continuity N50 is calculated within selected windows, not over the whole HAL.",
            "This replay does not yet enumerate per-base false Block edges or all missing truth edges.",
            "Cactus is auxiliary evidence embedded in the existing breakpoint classification, not truth.",
        ],
    }
    write_json(output_dir / "summary.json", summary)
    report = markdown_report(
        report_dir,
        breakpoint_path,
        round_dirs,
        truth_counts,
        metrics,
        acceptance_result,
        window_labels,
        global_rows,
    )
    (output_dir / "window_detection_evaluation.md").write_text(report, encoding="utf-8")
    log(f"evaluation complete: {output_dir}", args.quiet)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (
        FileNotFoundError,
        FileExistsError,
        ValueError,
        KeyError,
        RuntimeError,
    ) as error:
        print(f"[window-eval] ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
