#!/usr/bin/env python3
"""Disprove coefficientwise separation of the two ellipse orientations.

The direct ellipse map has scalar coefficients

    [w^(2*n+1)] phi_c(w) = sum_l alpha_(n,l) c^(n+2*l).

A tempting bookkeeping device replaces one such scalar monomial by
``a^(n+l) b^l`` and replaces the reverse edge ``c*J`` by ``a*J``;
formal adjoint then swaps ``a`` and ``b``.  On the physical diagonal
``a=b=c`` this recovers the original series.

This checker performs that lift in the exact grade-one delayed word
algebra.  The separate bidegrees do not vanish cyclically.  Only their
sum at fixed total degree has the physical cancellation.  Thus the
two orientations may not be proved or discarded independently.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from crabb_palindromic_elliptic_hessian import direct_map_coefficients
from repeated_crabb_cyclic_radial_volume import (
    DelayedQuotient,
    cyclic_reduce,
    radial_target,
)
from repeated_crabb_delayed_slack_anticommutator import (
    E,
    F,
    IDENTITY,
    S,
)


Bidegree = tuple[int, int]
BiSeries = dict[Bidegree, dict[str, Fraction]]


@dataclass(frozen=True)
class OrientationSplitObstructionRecord:
    """Exact grade-one obstruction data."""

    grade: int
    nonzero_degree_one_sectors: int
    nonzero_degree_two_sectors: int
    collapsed_degree_one_word_count: int
    collapsed_degree_two_target_difference_word_count: int
    operator_diagonal_recovery_passed: bool
    metric_diagonal_recovery_passed: bool
    all_checks_passed: bool


class BivariateGradeOne:
    """Sparse bivariate series in the exact grade-one word quotient."""

    maximum_degree = 2

    def __init__(self) -> None:
        self.quotient = DelayedQuotient(1)

    def add(
        self,
        left: BiSeries,
        right: BiSeries,
        coefficient: int = 1,
    ) -> BiSeries:
        """Add two bivariate series."""

        result = dict(left)
        for degree, polynomial in right.items():
            result[degree] = self.quotient.add(
                result.get(degree, {}),
                self.quotient.scale(coefficient, polynomial),
            )
            if not result[degree]:
                del result[degree]
        return result

    def multiply(self, left: BiSeries, right: BiSeries) -> BiSeries:
        """Multiply and truncate two bivariate series."""

        result: BiSeries = {}
        for (left_a, left_b), left_polynomial in left.items():
            for (right_a, right_b), right_polynomial in right.items():
                degree = (left_a + right_a, left_b + right_b)
                if sum(degree) > self.maximum_degree:
                    continue
                product = self.quotient.multiply(
                    left_polynomial,
                    right_polynomial,
                )
                result[degree] = self.quotient.add(
                    result.get(degree, {}),
                    product,
                )
                if not result[degree]:
                    del result[degree]
        return result

    def adjoint(self, series: BiSeries) -> BiSeries:
        """Take the formal adjoint and swap orientation degrees."""

        return {
            (right_degree, left_degree):
                self.quotient.polynomial_adjoint(polynomial)
            for (left_degree, right_degree), polynomial in series.items()
        }

    def power(self, series: BiSeries, exponent: int) -> BiSeries:
        """Raise a bivariate series to one nonnegative power."""

        result: BiSeries = {(0, 0): IDENTITY}
        for _ in range(exponent):
            result = self.multiply(result, series)
        return result

    def inverse(self, series: BiSeries) -> BiSeries:
        """Invert a series whose constant coefficient is the identity."""

        if series.get((0, 0)) != IDENTITY:
            raise ValueError("the constant coefficient must be the identity")
        result: BiSeries = {(0, 0): IDENTITY}
        for total_degree in range(1, self.maximum_degree + 1):
            for left_degree in range(total_degree + 1):
                degree = (left_degree, total_degree - left_degree)
                convolution: dict[str, Fraction] = {}
                for positive_degree, polynomial in series.items():
                    if positive_degree == (0, 0):
                        continue
                    remainder = (
                        degree[0] - positive_degree[0],
                        degree[1] - positive_degree[1],
                    )
                    if min(remainder) < 0 or remainder not in result:
                        continue
                    convolution = self.quotient.add(
                        convolution,
                        self.quotient.multiply(
                            polynomial,
                            result[remainder],
                        ),
                    )
                if convolution:
                    result[degree] = self.quotient.scale(
                        -1,
                        convolution,
                    )
        return result

    def operator(self) -> BiSeries:
        """Return the natural two-orientation lift of the direct map."""

        reverse = self.quotient.multiply(
            self.quotient.add(IDENTITY, F),
            self.quotient.multiply(
                self.quotient.polynomial_adjoint(S),
                self.quotient.add(IDENTITY, E),
            ),
        )
        pencil: BiSeries = {
            (0, 0): S,
            (1, 0): reverse,
        }
        result: BiSeries = {}
        coefficients = direct_map_coefficients(2, 3)
        for index, scalar_series in enumerate(coefficients):
            power = self.power(pencil, 2 * index + 1)
            for scalar_degree in range(3):
                scalar = scalar_series.coefficient(scalar_degree)
                if not scalar:
                    continue
                remainder = scalar_degree - index
                if remainder < 0 or remainder % 2:
                    raise RuntimeError("unexpected direct-map parity")
                paired_degree = remainder // 2
                for (left_degree, right_degree), polynomial in power.items():
                    degree = (
                        index + paired_degree + left_degree,
                        paired_degree + right_degree,
                    )
                    if sum(degree) > self.maximum_degree:
                        continue
                    result[degree] = self.quotient.add(
                        result.get(degree, {}),
                        self.quotient.scale(scalar, polynomial),
                    )
        return result

    def metric(self) -> BiSeries:
        """Return the physical metric with the active edge deleted."""

        metric = {
            (degree // 2, degree // 2): polynomial
            for degree, polynomial in enumerate(
                self.quotient.metric_series(),
            )
            if polynomial
        }
        metric.pop((1, 1), None)
        return metric

    def mass(self) -> tuple[BiSeries, BiSeries, BiSeries]:
        """Return the bivariate normalized mass, operator, and metric."""

        operator = self.operator()
        metric = self.metric()
        row_gram = self.multiply(
            self.multiply(self.adjoint(operator), {(0, 0): F}),
            operator,
        )
        denominator = self.add(metric, row_gram, -1)
        pulled_metric = self.multiply(
            self.multiply(self.adjoint(operator), metric),
            operator,
        )
        slack = self.add(metric, pulled_metric, -1)
        return (
            self.multiply(self.inverse(denominator), slack),
            operator,
            metric,
        )

    def collapse(self, series: BiSeries) -> list[dict[str, Fraction]]:
        """Set ``a=b=c`` and collect by total degree."""

        result = [{} for _ in range(self.maximum_degree + 1)]
        for degree, polynomial in series.items():
            total_degree = sum(degree)
            result[total_degree] = self.quotient.add(
                result[total_degree],
                polynomial,
            )
        return result


def exact_record() -> OrientationSplitObstructionRecord:
    """Construct and verify the exact grade-one obstruction."""

    algebra = BivariateGradeOne()
    mass, operator, metric = algebra.mass()
    operator_recovery = (
        algebra.collapse(operator)
        == algebra.quotient.operator_series()
    )
    expected_metric = algebra.quotient.metric_series()
    expected_metric[2] = {}
    metric_recovery = algebra.collapse(metric) == expected_metric

    cyclic_sectors = {
        degree: cyclic_reduce(polynomial, 0)
        for degree, polynomial in mass.items()
        if degree != (0, 0)
    }
    degree_one = {
        degree: polynomial
        for degree, polynomial in cyclic_sectors.items()
        if sum(degree) == 1 and polynomial
    }
    degree_two = {
        degree: polynomial
        for degree, polynomial in cyclic_sectors.items()
        if sum(degree) == 2 and polynomial
    }
    collapsed_mass = algebra.collapse(mass)
    collapsed_degree_one = cyclic_reduce(collapsed_mass[1], 0)
    active_difference = cyclic_reduce(
        algebra.quotient.add(
            collapsed_mass[2],
            algebra.quotient.scale(-1, radial_target(1)),
        ),
        0,
    )

    verified = (
        operator_recovery
        and metric_recovery
        and len(degree_one) == 2
        and len(degree_two) == 3
        and not collapsed_degree_one
        and not active_difference
    )
    if not verified:
        raise RuntimeError("the orientation-split obstruction failed")
    return OrientationSplitObstructionRecord(
        grade=1,
        nonzero_degree_one_sectors=len(degree_one),
        nonzero_degree_two_sectors=len(degree_two),
        collapsed_degree_one_word_count=len(collapsed_degree_one),
        collapsed_degree_two_target_difference_word_count=(
            len(active_difference)
        ),
        operator_diagonal_recovery_passed=operator_recovery,
        metric_diagonal_recovery_passed=metric_recovery,
        all_checks_passed=verified,
    )


def write_record(
    record: OrientationSplitObstructionRecord,
    output: Path,
) -> str:
    """Write one deterministic JSON record and return its hash."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(
        json.dumps(asdict(record), sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(output)
    return hashlib.sha256(output.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_orientation_split_obstruction_s70225.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact obstruction."""

    args = parse_args()
    record = exact_record()
    digest = write_record(record, args.output)
    print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
