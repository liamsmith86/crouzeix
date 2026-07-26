#!/usr/bin/env python3
"""Audit the all-grade metric-inverse cancellation of the interior fan."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from repeated_crabb_cyclic_radial_volume import DelayedQuotient
from repeated_crabb_delayed_slack_anticommutator import (
    F,
    IDENTITY,
    Polynomial,
)


@dataclass(frozen=True)
class MetricFanTelescopeRecord:
    """One exact delayed metric-inverse support audit."""

    grade: int
    active_inverse_word_count: int
    interior_fan_word_count: int
    radial_remainder_word_count: int
    lower_radial_coefficients: int
    radial_product_checks: int
    split_convolutions: int
    all_checks_passed: bool


def monomial(word: str) -> Polynomial:
    """Return one word with coefficient one."""

    return {word: Fraction(1)}


def q_word(index: int) -> str:
    """Return the word for ``Q_index=(S*)^index S^index``."""

    return "a" * index + "s" * index


def r_word(index: int) -> str:
    """Return the word for ``R_index=S^index (S*)^index``."""

    return "s" * index + "a" * index


def is_radial_word(word: str) -> bool:
    """Return whether a word is an initial or final radial power."""

    if not word:
        return True
    if len(word) % 2:
        return False
    index = len(word) // 2
    return word in {q_word(index), r_word(index)}


def interior_pair(grade: int, split: int) -> Polynomial:
    """Return the two fan words assigned to one inverse convolution."""

    left_index = split + 1
    right_index = grade - split + 1
    return {
        (
            "a" * left_index
            + "s" * (grade + 2)
            + "a" * right_index
        ): Fraction(1),
        (
            "s" * left_index
            + "a" * (grade + 2)
            + "s" * right_index
        ): Fraction(1),
    }


def interior_fan(
    quotient: DelayedQuotient,
    grade: int,
) -> Polynomial:
    """Return all ``2 grade - 2`` interior fan words."""

    result: Polynomial = {}
    for split in range(1, grade):
        result = quotient.add(result, interior_pair(grade, split))
    return result


def retained_compression(
    quotient: DelayedQuotient,
    polynomial: Polynomial,
) -> Polynomial:
    """Compress one coefficient away from the final defect."""

    retained = quotient.add(
        IDENTITY,
        quotient.scale(-1, F),
    )
    return quotient.multiply(
        quotient.multiply(retained, polynomial),
        retained,
    )


def radial(polynomial: Polynomial) -> bool:
    """Return whether every word in a polynomial is radial."""

    return all(is_radial_word(word) for word in polynomial)


def audit_grade(grade: int) -> MetricFanTelescopeRecord:
    """Audit the metric-inverse fan theorem in one grade."""

    if grade < 2:
        raise ValueError("the interior-fan theorem starts at grade two")

    quotient = DelayedQuotient(grade)
    metric = quotient.metric_series()
    metric[2 * grade] = {}
    metric_inverse = quotient.inverse_series(metric, IDENTITY)

    radial_product_checks: list[bool] = []
    for left_index in range(1, grade + 1):
        for right_index in range(1, grade + 1):
            if left_index + right_index > grade + 1:
                continue
            left_top = quotient.add(
                monomial(q_word(left_index)),
                quotient.scale(-1, monomial(r_word(left_index))),
            )
            right_top = quotient.add(
                quotient.scale(-1, monomial(q_word(right_index))),
                monomial(r_word(right_index)),
            )
            minimum = min(left_index, right_index)
            expected = quotient.add(
                monomial(q_word(minimum)),
                monomial(r_word(minimum)),
                quotient.scale(-2, IDENTITY),
            )
            radial_product_checks.append(
                quotient.multiply(left_top, right_top) == expected
            )

    lower_checks: list[bool] = []
    for lower_grade in range(1, grade):
        coefficient = retained_compression(
            quotient,
            metric_inverse[2 * lower_grade],
        )
        top_index = lower_grade + 1
        lower_checks.append(
            radial(coefficient)
            and coefficient.get(q_word(top_index)) == -1
            and coefficient.get(r_word(top_index)) == 1
        )

    split_checks: list[bool] = []
    for split in range(1, grade):
        convolution = quotient.multiply(
            metric[2 * split],
            metric_inverse[2 * (grade - split)],
        )
        convolution = retained_compression(quotient, convolution)
        remainder = quotient.add(
            convolution,
            quotient.scale(-1, interior_pair(grade, split)),
        )
        split_checks.append(radial(remainder))

    active = retained_compression(
        quotient,
        metric_inverse[2 * grade],
    )
    fan = interior_fan(quotient, grade)
    radial_remainder = quotient.add(active, fan)
    expected_fan_count = 2 * grade - 2
    verified = bool(
        all(lower_checks)
        and all(radial_product_checks)
        and all(split_checks)
        and len(fan) == expected_fan_count
        and radial(radial_remainder)
    )
    if not verified:
        raise RuntimeError(
            "the metric fan telescope failed: "
            f"grade={grade}, lower={all(lower_checks)}, "
            f"products={all(radial_product_checks)}, "
            f"splits={all(split_checks)}, "
            f"fan={len(fan)}/{expected_fan_count}, "
            f"radial={radial(radial_remainder)}"
        )
    return MetricFanTelescopeRecord(
        grade=grade,
        active_inverse_word_count=len(active),
        interior_fan_word_count=len(fan),
        radial_remainder_word_count=len(radial_remainder),
        lower_radial_coefficients=len(lower_checks),
        radial_product_checks=len(radial_product_checks),
        split_convolutions=len(split_checks),
        all_checks_passed=verified,
    )


def write_records(
    records: list[MetricFanTelescopeRecord],
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
    parser.add_argument("--maximum-grade", type=int, default=16)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_metric_fan_telescope_s70225.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist every exact fan audit."""

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
