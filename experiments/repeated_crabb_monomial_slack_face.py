#!/usr/bin/env python3
"""Audit the all-length monomial boundary-slack face exactly.

On a length-``k`` repeated Crabb shift, compare L219's boundary metric
with L117's exact elliptic-axis metric.  Their first difference is at
``q^k`` and, in balanced coordinates, is ``3I-2E``.  Since the exact
axis has zero Stein Schur residual, the boundary metric's first
residual is consequently ``2E_1``, equal to L228's anticommutator on
the monomial channel.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path


ScalarSeries = list[Fraction]


@dataclass(frozen=True)
class MonomialSlackFaceRecord:
    """One exact all-length scalar-channel audit."""

    transfer_grade: int
    state_dimension: int
    maximum_earlier_metric_difference: str
    lower_metric_face: str
    interior_metric_face: str
    upper_metric_face: str
    balanced_face_matches_3i_minus_2e: bool
    slack_face_matches_e_plus_2e1: bool
    schur_face_matches_anticommutator: bool
    all_checks_passed: bool


def zero_series(maximum_degree: int) -> ScalarSeries:
    """Return a zero scalar series."""

    return [Fraction(0) for _ in range(maximum_degree + 1)]


def add(left: ScalarSeries, right: ScalarSeries) -> ScalarSeries:
    """Add scalar series."""

    return [
        left_value + right_value
        for left_value, right_value in zip(left, right, strict=True)
    ]


def multiply(left: ScalarSeries, right: ScalarSeries) -> ScalarSeries:
    """Multiply truncated scalar series."""

    maximum_degree = len(left) - 1
    result = zero_series(maximum_degree)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            degree = left_degree + right_degree
            if degree <= maximum_degree:
                result[degree] += left_value * right_value
    return result


def inverse(series: ScalarSeries) -> ScalarSeries:
    """Invert a scalar series with nonzero constant term."""

    maximum_degree = len(series) - 1
    result = zero_series(maximum_degree)
    result[0] = 1 / series[0]
    for degree in range(1, maximum_degree + 1):
        result[degree] = -sum(
            series[positive] * result[degree - positive]
            for positive in range(1, degree + 1)
        ) / series[0]
    return result


def shifted_sech(
    exponent: int,
    shift: int,
    maximum_degree: int,
) -> ScalarSeries:
    """Return ``c^(-shift) sech(exponent*(-log c))``."""

    absolute_exponent = abs(exponent)
    result = zero_series(maximum_degree)
    if absolute_exponent == 0:
        if shift:
            raise ValueError("a shifted zero exponent is singular")
        result[0] = Fraction(1)
        return result

    multiple = 0
    while True:
        degree = absolute_exponent * (2 * multiple + 1) - shift
        if degree > maximum_degree:
            break
        if degree >= 0:
            result[degree] += 2 * (-1) ** multiple
        multiple += 1
    return result


def axis_physical_weight(
    length: int,
    level: int,
    maximum_degree: int,
) -> ScalarSeries:
    """Return L117's exact physical weight as a ``c`` series."""

    numerator = zero_series(maximum_degree)
    maximum_alias = 2 + maximum_degree // (2 * length)
    for alias in range(-maximum_alias, maximum_alias + 1):
        exponent = level + 2 * length * alias
        numerator = add(
            numerator,
            shifted_sech(exponent, level, maximum_degree),
        )

    denominator = zero_series(maximum_degree)
    for alias in range(-maximum_alias, maximum_alias + 1):
        denominator = add(
            denominator,
            shifted_sech(
                2 * length * alias,
                0,
                maximum_degree,
            ),
        )
    return multiply(numerator, inverse(denominator))


def geometric_weight(
    exponent: int,
    scale: int,
    maximum_degree: int,
) -> ScalarSeries:
    """Return ``scale/(1+c^exponent)``."""

    result = zero_series(maximum_degree)
    multiple = 0
    while exponent * multiple <= maximum_degree:
        result[exponent * multiple] = scale * (-1) ** multiple
        multiple += 1
    return result


def boundary_physical_weight(
    length: int,
    level: int,
    maximum_degree: int,
) -> ScalarSeries:
    """Return L219's physical metric on a monomial shift."""

    if level == 0:
        result = zero_series(maximum_degree)
        result[0] = Fraction(1)
        result[2 * length] = Fraction(1)
        return result
    if level == length:
        return geometric_weight(2 * length, 4, maximum_degree)

    result = geometric_weight(2 * level, 2, maximum_degree)
    result[2 * (length - level)] += 2
    return result


def exact_record(grade: int) -> MonomialSlackFaceRecord:
    """Return one exact monomial face audit."""

    maximum_degree = 2 * grade
    differences: list[ScalarSeries] = []
    for level in range(grade + 1):
        differences.append([
            boundary - axis
            for boundary, axis in zip(
                boundary_physical_weight(
                    grade,
                    level,
                    maximum_degree,
                ),
                axis_physical_weight(
                    grade,
                    level,
                    maximum_degree,
                ),
                strict=True,
            )
        ])

    maximum_earlier = max(
        (
            abs(value)
            for difference in differences
            for value in difference[:maximum_degree]
        ),
        default=Fraction(0),
    )
    physical_face = [
        difference[maximum_degree]
        for difference in differences
    ]
    expected_physical = [
        Fraction(1),
        *(Fraction(6) for _ in range(grade - 1)),
        Fraction(12),
    ]
    balanced_face = [
        value / (1 if level == 0 else 4 if level == grade else 2)
        for level, value in enumerate(physical_face)
    ]
    expected_balanced = [
        Fraction(1),
        *(Fraction(3) for _ in range(grade)),
    ]

    slack_face = [
        balanced_face[0],
        *(
            balanced_face[level] - balanced_face[level - 1]
            for level in range(1, grade + 1)
        ),
    ]
    expected_slack = [
        Fraction(1),
        Fraction(2),
        *(Fraction(0) for _ in range(grade - 1)),
    ]
    schur_face = slack_face[1:]
    expected_schur = [
        Fraction(2),
        *(Fraction(0) for _ in range(grade - 1)),
    ]

    physical_matches = physical_face == expected_physical
    balanced_matches = balanced_face == expected_balanced
    slack_matches = slack_face == expected_slack
    schur_matches = schur_face == expected_schur
    verified = (
        maximum_earlier == 0
        and physical_matches
        and balanced_matches
        and slack_matches
        and schur_matches
    )
    if not verified:
        raise RuntimeError(
            "the monomial slack-face audit failed: "
            f"grade={grade}, physical={physical_face}, "
            f"balanced={balanced_face}, slack={slack_face}"
        )
    interior_face = (
        physical_face[1]
        if grade > 1
        else Fraction(0)
    )
    return MonomialSlackFaceRecord(
        transfer_grade=grade,
        state_dimension=grade + 1,
        maximum_earlier_metric_difference=str(maximum_earlier),
        lower_metric_face=str(physical_face[0]),
        interior_metric_face=str(interior_face),
        upper_metric_face=str(physical_face[-1]),
        balanced_face_matches_3i_minus_2e=balanced_matches,
        slack_face_matches_e_plus_2e1=slack_matches,
        schur_face_matches_anticommutator=schur_matches,
        all_checks_passed=verified,
    )


def exact_records(maximum_grade: int) -> list[MonomialSlackFaceRecord]:
    """Return exact audits through one transfer grade."""

    return [
        exact_record(grade)
        for grade in range(1, maximum_grade + 1)
    ]


def write_records(
    records: list[MonomialSlackFaceRecord],
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
    parser.add_argument("--maximum-grade", type=int, default=12)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_monomial_slack_face_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact monomial audits."""

    args = parse_args()
    records = exact_records(args.maximum_grade)
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
