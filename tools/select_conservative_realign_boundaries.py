#!/usr/bin/env python3
"""Build an evaluation-only allowlist for the safest realignment experiment.

The output is deliberately only a list of stable boundary IDs.  RaMAx still
rechecks every graph-side structural invariant and every minipoa QC threshold
before an ID can be applied.  Direct truth MAF is used here only to approve a
controlled Alignathon experiment; this file must not be used as a production
detector.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class BoundaryState:
    accepted: bool = True
    rows: int = 0
    targets: set[str] = field(default_factory=set)
    reasons: set[str] = field(default_factory=set)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description=(
            "Select boundaries whose direct-MAF truth rows all support the "
            "existing two flanks as one continuous relationship."
        )
    )
    result.add_argument("--boundary-truth-labels", required=True, type=Path)
    result.add_argument("--output", required=True, type=Path)
    result.add_argument(
        "--expected-target-count",
        type=int,
        default=0,
        help="Require this many distinct non-reference target genomes (0: no fixed count).",
    )
    result.add_argument(
        "--summary",
        type=Path,
        help="Optional JSON summary (default: <output>.summary.json).",
    )
    return result


def false_like(value: str) -> bool:
    return value.strip().lower() in {"", "0", "false", "no"}


def evaluate_row(row: dict[str, str]) -> list[str]:
    reasons: list[str] = []
    if row["truth_source"] != "direct_maf":
        reasons.append("NOT_DIRECT_MAF")
    if row["boundary_truth_label"] != "MERGE_BOUNDARY":
        reasons.append("NOT_MERGE")
    if row["truth_bridge_relation"] != "EXISTING_FLANKS_CONTINUOUS":
        reasons.append("NOT_EXISTING_FLANK_CONTINUITY")
    if row["detector_recommends_bridge"] != "1":
        reasons.append("DETECTOR_DID_NOT_RECOMMEND")
    if row["truth_bridge_matches_prediction"] != "1":
        reasons.append("TRUTH_DOES_NOT_MATCH_PREDICTED_FLANKS")
    if not false_like(row["truth_mapping_ambiguous"]):
        reasons.append("TRUTH_MAPPING_AMBIGUOUS")
    if not false_like(row["truth_left_truncated"]):
        reasons.append("LEFT_CANDIDATES_TRUNCATED")
    if not false_like(row["truth_right_truncated"]):
        reasons.append("RIGHT_CANDIDATES_TRUNCATED")
    return reasons


def atomic_write(path: Path, content_writer) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as output:
            content_writer(output)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def main() -> int:
    args = parser().parse_args()
    states: dict[str, BoundaryState] = {}
    with args.boundary_truth_labels.open("r", encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        required = {
            "boundary_id",
            "target_species",
            "truth_source",
            "boundary_truth_label",
            "truth_bridge_relation",
            "detector_recommends_bridge",
            "truth_bridge_matches_prediction",
            "truth_mapping_ambiguous",
            "truth_left_truncated",
            "truth_right_truncated",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            boundary_id = row["boundary_id"]
            state = states.setdefault(boundary_id, BoundaryState())
            state.rows += 1
            state.targets.add(row["target_species"])
            failures = evaluate_row(row)
            if failures:
                state.accepted = False
                state.reasons.update(failures)

    rejection_counts: dict[str, int] = {}
    approved: list[str] = []
    for boundary_id, state in states.items():
        if (
            args.expected_target_count
            and len(state.targets) != args.expected_target_count
        ):
            state.accepted = False
            state.reasons.add("TARGET_COUNT_MISMATCH")
        if state.accepted:
            approved.append(boundary_id)
        else:
            for reason in state.reasons:
                rejection_counts[reason] = rejection_counts.get(reason, 0) + 1
    approved.sort()

    def write_allowlist(output) -> None:
        output.write("boundary_id\n")
        for boundary_id in approved:
            output.write(boundary_id + "\n")

    atomic_write(args.output, write_allowlist)
    summary_path = args.summary or Path(str(args.output) + ".summary.json")
    payload = {
        "source": str(args.boundary_truth_labels.resolve()),
        "output": str(args.output.resolve()),
        "boundary_count": len(states),
        "approved_boundary_count": len(approved),
        "expected_target_count": args.expected_target_count,
        "rejection_counts": dict(sorted(rejection_counts.items())),
        "policy": {
            "truth_source": "direct_maf",
            "boundary_truth_label": "MERGE_BOUNDARY",
            "truth_bridge_relation": "EXISTING_FLANKS_CONTINUOUS",
            "detector_recommends_bridge": True,
            "truth_bridge_matches_prediction": True,
            "truth_mapping_ambiguous": False,
            "candidate_lists_truncated": False,
        },
    }

    def write_summary(output) -> None:
        json.dump(payload, output, indent=2, sort_keys=True)
        output.write("\n")

    atomic_write(summary_path, write_summary)
    print(
        f"approved {len(approved):,}/{len(states):,} boundaries; "
        f"wrote {args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
