#!/usr/bin/env python3
"""Audit the trace-only one-delay closed-return recursion."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
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
    Polynomial,
    S,
    STAR,
)


@dataclass(frozen=True)
class ClosedReturnRecursionRecord:
    """One exact one-delay operator/trace comparison."""

    grade: int
    full_face_word_count: int
    tail_face_word_count: int
    difference_word_count: int
    cyclic_difference_word_count: int
    deep_remainder_present: bool
    deep_remainder_depth: int
    radial_remainder_word_count: int
    all_earlier_faces_vanish: bool
    all_checks_passed: bool


def deep_divergence(grade: int) -> Polynomial:
    """Return the observed depth-``grade+1`` Stein divergence."""

    middle = (
        "s"
        + "a" * (grade + 2)
        + "s" * (grade + 2)
        + "a"
    )
    shifted = "s" + middle + "a"
    return {
        middle: Fraction(1),
        shifted: Fraction(-1),
    }


def is_radial_word(word: str) -> bool:
    """Return whether a word is a one-sided radial power."""

    if not word:
        return True
    for split in range(1, len(word)):
        if (
            word == "a" * split + "s" * split
            or word == "s" * split + "a" * split
        ):
            return True
    return False


def audit_grade(grade: int) -> ClosedReturnRecursionRecord:
    """Audit the full-versus-deflated first closed-return face."""

    if grade < 2:
        raise ValueError("one-delay recursion starts at grade two")
    quotient = DelayedQuotient(grade)
    full = quotient.closed_return_defect_for(
        S,
        E,
        F,
        IDENTITY,
        2 * grade,
    )

    retained = quotient.add(
        IDENTITY,
        quotient.scale(-1, F),
    )
    tail_partial = quotient.multiply(S, retained)
    tail_final = quotient.multiply(
        quotient.multiply(S, F),
        STAR,
    )
    tail = quotient.closed_return_defect_for(
        tail_partial,
        E,
        tail_final,
        retained,
        2 * (grade - 1),
    )
    difference = quotient.add(
        full[2 * grade],
        quotient.scale(-1, tail[2 * (grade - 1)]),
    )
    cyclic_difference = cyclic_reduce(
        difference,
        grade - 1,
    )
    deep = deep_divergence(grade)
    deep_present = all(
        difference.get(word) == coefficient
        for word, coefficient in deep.items()
    )
    radial_remainder = quotient.add(
        difference,
        quotient.scale(-1, deep),
    )
    radial_support = all(
        is_radial_word(word)
        for word in radial_remainder
    )
    earlier_vanish = bool(
        not any(full[: 2 * grade])
        and not any(tail[: 2 * (grade - 1)])
    )

    expected_grade_two = {
        "": Fraction(-2),
        "as": Fraction(5),
        "sa": Fraction(-1),
        "aass": Fraction(-3),
        "ssaa": Fraction(1),
        "saaaassssa": Fraction(1),
        "ssaaaassssaa": Fraction(-1),
    }
    grade_two_exact = (
        grade != 2
        or difference == expected_grade_two
    )
    verified = bool(
        earlier_vanish
        and not cyclic_difference
        and deep_present
        and radial_support
        and grade_two_exact
    )
    if not verified:
        raise RuntimeError(
            "the closed-return recursion audit failed: "
            f"grade={grade}, cyclic={len(cyclic_difference)}, "
            f"deep={deep_present}, radial={radial_support}"
        )
    return ClosedReturnRecursionRecord(
        grade=grade,
        full_face_word_count=len(full[2 * grade]),
        tail_face_word_count=len(tail[2 * (grade - 1)]),
        difference_word_count=len(difference),
        cyclic_difference_word_count=len(cyclic_difference),
        deep_remainder_present=deep_present,
        deep_remainder_depth=grade + 1,
        radial_remainder_word_count=len(radial_remainder),
        all_earlier_faces_vanish=earlier_vanish,
        all_checks_passed=verified,
    )


def exact_records(
    maximum_grade: int,
) -> list[ClosedReturnRecursionRecord]:
    """Return exact audits from grade two through one maximum."""

    return [
        audit_grade(grade)
        for grade in range(2, maximum_grade + 1)
    ]


def write_records(
    records: list[ClosedReturnRecursionRecord],
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
            "repeated_crabb_closed_return_recursion_s70225.jsonl",
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact recursion audits."""

    args = parse_args()
    records = exact_records(args.maximum_grade)
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
