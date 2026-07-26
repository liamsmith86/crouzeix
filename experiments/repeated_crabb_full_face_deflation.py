#!/usr/bin/env python3
"""Audit the algebraic deflation of the proposed full-face support form."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from repeated_crabb_delayed_slack_anticommutator import (
    IDENTITY,
    Polynomial,
    add,
    multiply,
    reduce_delay,
    scale,
)


@dataclass(frozen=True)
class FullFaceDeflationRecord:
    """One exact tail-transport audit."""

    grade: int
    radial_generator_checks: int
    right_orbit_checks: int
    remote_word_checks: int
    all_checks_passed: bool


def monomial(word: str) -> Polynomial:
    """Return one word with coefficient one."""

    return {word: Fraction(1)}


def product(*polynomials: Polynomial) -> Polynomial:
    """Multiply an ordered list of word polynomials."""

    result = IDENTITY
    for polynomial in polynomials:
        result = multiply(result, polynomial)
    return result


def power(polynomial: Polynomial, exponent: int) -> Polynomial:
    """Raise one word polynomial to a nonnegative power."""

    result = IDENTITY
    for _ in range(exponent):
        result = multiply(result, polynomial)
    return result


def reduced(polynomial: Polynomial, grade: int) -> Polynomial:
    """Reduce in the complete-delay quotient for one grade."""

    return reduce_delay(polynomial, grade - 1)


def q(grade: int, index: int) -> Polynomial:
    """Return ``Q_index=(S*)^index S^index``."""

    return reduced(monomial("a" * index + "s" * index), grade)


def r(grade: int, index: int) -> Polynomial:
    """Return ``R_index=S^index (S*)^index``."""

    return reduced(monomial("s" * index + "a" * index), grade)


def deep(grade: int) -> Polynomial:
    """Return ``G_grade``."""

    return reduced(
        monomial(
            "s"
            + "a" * (grade + 2)
            + "s" * (grade + 2)
            + "a"
        ),
        grade,
    )


def remote_words(grade: int) -> tuple[Polynomial, ...]:
    """Return the four proposed universal remote words."""

    return (
        reduced(
            monomial(
                "a"
                + "s" * (grade + 2)
                + "a" * (grade + 1)
            ),
            grade,
        ),
        reduced(
            monomial(
                "s" * (grade + 1)
                + "a" * (grade + 2)
                + "s"
            ),
            grade,
        ),
        r(grade, grade + 2),
        deep(grade),
    )


def audit_grade(grade: int) -> FullFaceDeflationRecord:
    """Verify every tail-to-full identity used by the reduction."""

    if grade < 2:
        raise ValueError("deflation starts at grade two")

    identity = IDENTITY
    tail_identity = r(grade, 1)
    tail_s = reduced(monomial("ssa"), grade)
    tail_a = reduced(monomial("saa"), grade)

    radial_checks: list[bool] = [tail_identity == r(grade, 1)]
    for index in range(1, grade):
        tail_q = reduced(
            product(
                power(tail_a, index),
                power(tail_s, index),
            ),
            grade,
        )
        expected = reduced(
            add(
                q(grade, index),
                r(grade, 1),
                scale(-1, identity),
            ),
            grade,
        )
        radial_checks.append(tail_q == expected)

    right_checks: list[bool] = []
    for index in range(1, grade + 1):
        tail_r = reduced(
            product(
                power(tail_s, index),
                power(tail_a, index),
            ),
            grade,
        )
        right_checks.append(tail_r == r(grade, index + 1))

    tail_grade = grade - 1
    tail_remote = (
        reduced(
            product(
                tail_a,
                power(tail_s, tail_grade + 2),
                power(tail_a, tail_grade + 1),
            ),
            grade,
        ),
        reduced(
            product(
                power(tail_s, tail_grade + 1),
                power(tail_a, tail_grade + 2),
                tail_s,
            ),
            grade,
        ),
        reduced(
            product(
                power(tail_s, tail_grade + 2),
                power(tail_a, tail_grade + 2),
            ),
            grade,
        ),
        reduced(
            product(
                tail_s,
                power(tail_a, tail_grade + 2),
                power(tail_s, tail_grade + 2),
                tail_a,
            ),
            grade,
        ),
    )
    full_remote = remote_words(grade)
    shifted_deep = reduced(
        product(monomial("s"), deep(grade), monomial("a")),
        grade,
    )
    remote_expected = (
        full_remote[0],
        full_remote[1],
        full_remote[2],
        shifted_deep,
    )
    remote_checks = [
        actual == expected
        for actual, expected in zip(
            tail_remote,
            remote_expected,
            strict=True,
        )
    ]

    verified = all(radial_checks + right_checks + remote_checks)
    if not verified:
        raise RuntimeError(f"full-face deflation failed at grade {grade}")
    return FullFaceDeflationRecord(
        grade=grade,
        radial_generator_checks=len(radial_checks),
        right_orbit_checks=len(right_checks),
        remote_word_checks=len(remote_checks),
        all_checks_passed=verified,
    )


def write_records(
    records: list[FullFaceDeflationRecord],
    output: Path,
) -> str:
    """Write deterministic records atomically and return their hash."""

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
    parser.add_argument("--maximum-grade", type=int, default=12)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_full_face_deflation_s70225.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the exact tail-transport audit."""

    args = parse_args()
    records = [
        audit_grade(grade)
        for grade in range(2, args.maximum_grade + 1)
    ]
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
