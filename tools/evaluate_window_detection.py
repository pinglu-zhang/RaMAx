#!/usr/bin/env python3
"""Evaluate RaMAx phase-1 candidate windows against Alignathon labels.

This first-stage evaluator deliberately consumes the existing
breakpoint_classifications.tsv.gz file. That file already combines old RaMAx,
Cactus and Alignathon truth evidence. It does not modify a HAL, MAF or graph.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence, TextIO


POSITIVE_LABEL = "over_fragmentation"
PROTECTED_LABELS = {"correct_truth_break", "avoided_cactus_overmerge"}
CLASSIFIED_LABELS = {POSITIVE_LABEL, *PROTECTED_LABELS}


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
        self.by_pair: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
        self.by_exact_reference: dict[
            tuple[str, str, str, int, int, int, int], list[dict[str, str]]
        ] = defaultdict(list)
        for row in rows:
            pair_key = (
                row["reference_species"],
                row["target_species"],
                row["ref_chrom"],
            )
            self.by_pair[pair_key].append(row)
            exact_key = (
                *pair_key,
                as_int(row, "left_ref_start"),
                as_int(row, "left_ref_end"),
                as_int(row, "right_ref_start"),
                as_int(row, "right_ref_end"),
            )
            self.by_exact_reference[exact_key].append(row)

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

        candidates: list[tuple[int, dict[str, str]]] = []
        for truth in self.by_pair.get(pair_key, []):
            reference_differences = (
                abs(as_int(boundary, "left_ref_start") - as_int(truth, "left_ref_start")),
                abs(as_int(boundary, "left_ref_end") - as_int(truth, "left_ref_end")),
                abs(as_int(boundary, "right_ref_start") - as_int(truth, "right_ref_start")),
                abs(as_int(boundary, "right_ref_end") - as_int(truth, "right_ref_end")),
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
) -> list[dict[str, str | int]]:
    boundary_by_key = {
        (row["window_id"], row["boundary_id"]): row for row in boundaries
    }
    output: list[dict[str, str | int]] = []
    seen_evidence: set[tuple[str, str, str]] = set()
    for evidence in boundary_genomes:
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
    output.sort(
        key=lambda row: (
            str(row["source_round_dir"]),
            str(row["window_id"]),
            str(row["boundary_id"]),
            str(row["target_species"]),
        )
    )
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
) -> tuple[list[dict], dict[str, dict[str, list[int]]]]:
    boundaries_by_window: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in boundaries:
        boundaries_by_window[row["window_id"]].append(row)
    labels_by_boundary: dict[tuple[str, str], list[dict[str, str | int]]] = defaultdict(list)
    for row in boundary_labels:
        labels_by_boundary[(str(row["window_id"]), str(row["boundary_id"]))].append(row)

    output: list[dict] = []
    distributions: dict[str, dict[str, list[int]]] = {}
    for window in windows:
        window_id = window["window_id"]
        rows = sorted(
            boundaries_by_window.get(window_id, []),
            key=lambda row: (
                as_int(row, "left_ref_start"),
                as_int(row, "right_ref_start"),
                row["boundary_id"],
            ),
        )
        if not rows:
            continue
        projected_bridges = {
            row["boundary_id"] for row in rows if as_bool(row, "detector_recommends_bridge")
        }
        oracle_bridges: set[str] = set()
        for row in rows:
            evidence = labels_by_boundary.get((window_id, row["boundary_id"]), [])
            labels = {str(item["truth_classification"]) for item in evidence}
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


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate RaMAx phase-1 windows against Alignathon breakpoint labels."
    )
    parser.add_argument("--window-report-dir", type=Path, required=True)
    parser.add_argument("--breakpoint-classifications", type=Path, required=True)
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
        "--force",
        action="store_true",
        help="Overwrite evaluator-owned files in an existing output directory.",
    )
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)
    if args.coordinate_tolerance < 0:
        parser.error("--coordinate-tolerance must be >= 0")
    return args


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    report_dir = args.window_report_dir.resolve()
    breakpoint_path = args.breakpoint_classifications.resolve()
    output_dir = args.output_dir.resolve()

    manifest_path = report_dir / "run_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing run manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "complete":
        raise ValueError(f"Window report is not complete: {manifest_path}")
    if not breakpoint_path.is_file():
        raise FileNotFoundError(breakpoint_path)

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
    boundary_labels = evaluate_boundaries(boundaries, boundary_genomes, matcher)
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
    except (FileNotFoundError, FileExistsError, ValueError, KeyError) as error:
        print(f"[window-eval] ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
