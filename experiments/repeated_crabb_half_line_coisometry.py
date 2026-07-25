#!/usr/bin/env python3
"""Audit the coisometric normalization of the elliptic half-line.

L240 proves ``D - A* D A = d d*``.  The normalized defect
``e = D**(-1/2) d`` has norm one, so the similar operator
``D**(1/2) A D**(-1/2)`` is a partial isometry.  The analytic proof
shows that it is onto and hence a coisometry.  This script independently
checks the equivalent dual Stein identity in exact coisometric word
algebra and the scalar theta normalization coefficient by coefficient.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from repeated_crabb_half_line_stein_factor import (
    IDENTITY,
    Series,
    add,
    inverse_theta_three,
    metric_series,
    multiply,
    operator_series,
    scale,
    series_add,
    series_adjoint,
    series_multiply,
    zero_series,
)


@dataclass(frozen=True)
class HalfLineCoisometryRecord:
    """One exact coefficient audit."""

    degree: int
    inverse_metric_word_count: int
    dual_stein_residual_word_count: int
    defect_normalization_residual: str
    all_checks_passed: bool


def inverse_operator_series(series: Series) -> Series:
    """Return the formal inverse of a series with identity constant term."""

    maximum_degree = len(series) - 1
    result = zero_series(maximum_degree)
    result[0] = IDENTITY
    for degree in range(1, maximum_degree + 1):
        convolution = {}
        for positive_degree in range(1, degree + 1):
            convolution = add(
                convolution,
                multiply(
                    series[positive_degree],
                    result[degree - positive_degree],
                ),
            )
        result[degree] = scale(-1, convolution)
    return result


def scalar_convolution(
    left: list[Fraction],
    right: list[Fraction],
) -> list[Fraction]:
    """Multiply equal-length scalar series."""

    maximum_degree = len(left) - 1
    result = [Fraction(0) for _ in range(maximum_degree + 1)]
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            degree = left_degree + right_degree
            if degree <= maximum_degree:
                result[degree] += left_value * right_value
    return result


def defect_normalization_series(maximum_degree: int) -> list[Fraction]:
    """Return the series of ``d* D^{-1} d``."""

    inverse_theta = inverse_theta_three(maximum_degree)
    bracket = [Fraction(0) for _ in range(maximum_degree + 1)]
    bracket[0] = Fraction(1)
    for frequency in range(1, maximum_degree // 2 + 1):
        base_degree = 2 * frequency
        denominator_power = 0
        while base_degree + 4 * frequency * denominator_power <= maximum_degree:
            degree = base_degree + 4 * frequency * denominator_power
            bracket[degree] += 4 * (-1) ** denominator_power
            denominator_power += 1
    inverse_theta_square = scalar_convolution(
        inverse_theta,
        inverse_theta,
    )
    return scalar_convolution(inverse_theta_square, bracket)


def exact_records(maximum_degree: int) -> list[HalfLineCoisometryRecord]:
    """Return exact dual-Stein and normalization audits."""

    operator = operator_series(maximum_degree)
    inverse_metric = inverse_operator_series(
        metric_series(maximum_degree)
    )
    dual = series_add(
        series_multiply(
            series_multiply(operator, inverse_metric),
            series_adjoint(operator),
        ),
        [scale(-1, coefficient) for coefficient in inverse_metric],
    )
    normalization = defect_normalization_series(maximum_degree)

    records: list[HalfLineCoisometryRecord] = []
    for degree in range(maximum_degree + 1):
        normalization_residual = normalization[degree] - (
            Fraction(1) if degree == 0 else Fraction(0)
        )
        verified = (
            not dual[degree]
            and normalization_residual == 0
        )
        if not verified:
            raise RuntimeError(
                "the half-line coisometry audit failed: "
                f"degree={degree}, dual={dual[degree]}, "
                f"normalization={normalization_residual}"
            )
        records.append(
            HalfLineCoisometryRecord(
                degree=degree,
                inverse_metric_word_count=len(inverse_metric[degree]),
                dual_stein_residual_word_count=len(dual[degree]),
                defect_normalization_residual=str(
                    normalization_residual
                ),
                all_checks_passed=verified,
            )
        )
    return records


def write_records(
    records: list[HalfLineCoisometryRecord],
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
    parser.add_argument("--maximum-degree", type=int, default=8)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_half_line_coisometry_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact coisometry audits."""

    args = parse_args()
    records = exact_records(args.maximum_degree)
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
