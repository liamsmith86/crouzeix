#!/usr/bin/env python3
"""Audit one-delay covariance of the edge-deleted scalar volume.

For a completely delayed grade ``k``, remove the first left wandering
layer and call the retained partial isometry ``S_1``.  The live
all-grade statement is the associated coefficient identity

    [c^(2k)] mu_S^circ = [c^(2k-2)] mu_(S_1)^circ.

This checker constructs both sides in the same exact rational word
algebra and verifies equality after trace cyclicity through finite
grades.  It is an exact finite-grade audit, not the missing
arbitrary-grade proof.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

from repeated_crabb_cyclic_radial_volume import (
    DelayedQuotient,
    cyclic_reduce,
)
from repeated_crabb_delayed_slack_anticommutator import (
    E,
    F,
    IDENTITY,
    S,
    STAR,
)


@dataclass(frozen=True)
class VolumeDelayCovarianceRecord:
    """One exact associated-volume covariance audit."""

    grade: int
    full_face_degree: int
    tail_face_degree: int
    full_face_word_count: int
    embedded_tail_face_word_count: int
    quotient_difference_word_count: int
    cyclic_difference_word_count: int
    all_checks_passed: bool


def audit_grade(grade: int) -> VolumeDelayCovarianceRecord:
    """Audit the active coefficient after removing one clean layer."""

    if grade < 2:
        raise ValueError("one-delay covariance starts at grade two")

    quotient = DelayedQuotient(grade)
    _, _, full_face = quotient.mass_components()

    retained_identity = quotient.add(
        IDENTITY,
        quotient.scale(-1, F),
    )
    tail = quotient.multiply(S, retained_identity)
    tail_final = quotient.multiply(
        quotient.multiply(S, F),
        STAR,
    )
    tail_degree = 2 * (grade - 1)
    _, _, embedded_tail_face = quotient.mass_components_for(
        tail,
        E,
        tail_final,
        retained_identity,
        tail_degree,
    )

    difference = quotient.add(
        full_face,
        quotient.scale(-1, embedded_tail_face),
    )
    cyclic_difference = cyclic_reduce(
        difference,
        grade - 1,
    )
    verified = not cyclic_difference
    if not verified:
        raise RuntimeError(
            "the exact volume delay-covariance audit failed: "
            f"grade={grade}, cyclic={len(cyclic_difference)}"
        )

    return VolumeDelayCovarianceRecord(
        grade=grade,
        full_face_degree=2 * grade,
        tail_face_degree=tail_degree,
        full_face_word_count=len(full_face),
        embedded_tail_face_word_count=len(embedded_tail_face),
        quotient_difference_word_count=len(difference),
        cyclic_difference_word_count=len(cyclic_difference),
        all_checks_passed=verified,
    )


def exact_records(maximum_grade: int) -> list[VolumeDelayCovarianceRecord]:
    """Return exact audits from grade two through one maximum."""

    return [audit_grade(grade) for grade in range(2, maximum_grade + 1)]


def write_records(
    records: list[VolumeDelayCovarianceRecord],
    output: Path,
) -> str:
    """Write deterministic JSON Lines atomically and return its hash."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)
    return hashlib.sha256(output.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-grade", type=int, default=5)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_volume_delay_covariance_s70225.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact covariance records."""

    args = parse_args()
    records = exact_records(args.maximum_grade)
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
