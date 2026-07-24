#!/usr/bin/env python3
"""Regenerate the first canonical-metric jet after L176 recentering.

The exact disk path is

    H(s) = I/2 + s Z(z) + s**2 B(z wedge J conjugate(z)).

For L122's canonical rank-one Stein metric, let ``lambda_-`` and
``lambda_+`` be the two simple generalized endpoint levels.  The sign
of

    delta(s) = lambda_+(s) - 4 lambda_-(s)

is the sign of the condition-square excess over four.  L176--L177
force its terms through degree four to vanish.  This checker extracts
the next exact coefficients, including the larger exceptional zero
set of the sixth-order term.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Sequence

import sympy as sp

from crabb_circular_normal_quadratic_exact import (
    deterministic_disk_direction,
    disk_model_from_hermitian_series,
)
from crabb_disk_toeplitz_quartic import (
    extend,
    generalized_endpoint_polynomial_series,
    lyapunov_series,
)
from crabb_full_disk_correction_isometry import plucker_correction


@dataclass(frozen=True)
class FullDiskBaseJetRecord:
    """One exact recentered canonical-metric endpoint audit."""

    dimension: int
    length: int
    direction_kind: str
    direction: tuple[str, ...]
    series_order: int
    vanishing_through_degree: int
    leading_degree: int
    leading_delta: str
    identity_verified: bool


def toeplitz_direction(
    direction: Sequence[sp.Expr],
) -> sp.Matrix:
    """Return the Hermitian Toeplitz matrix encoded by ``direction``."""

    length = len(direction)
    matrix = sp.zeros(length)
    for offset in range(1, length):
        for row in range(length - offset):
            matrix[row, row + offset] = direction[offset]
            matrix[row + offset, row] = sp.conjugate(
                direction[offset]
            )
    return matrix


def endpoint_delta_series(
    direction: tuple[sp.Expr, ...],
    order: int,
) -> list[sp.Expr]:
    """Return ``lambda_+ - 4 lambda_-`` through ``order`` exactly."""

    length = len(direction)
    dimension = length + 1
    hermitian = [
        extend(sp.eye(length) / 2),
        extend(toeplitz_direction(direction)),
        extend(plucker_correction(direction)),
        *[sp.zeros(dimension) for _ in range(order - 2)],
    ]
    operator, disk_metric = disk_model_from_hermitian_series(
        hermitian
    )
    return endpoint_delta_from_disk_series(
        hermitian,
        disk_metric,
        operator,
    )


def endpoint_delta_from_disk_series(
    hermitian: Sequence[sp.Matrix],
    disk_metric: Sequence[sp.Matrix],
    operator: Sequence[sp.Matrix],
) -> list[sp.Expr]:
    """Return the endpoint excess for arbitrary equal-order disk series."""

    if not (
        len(hermitian) == len(disk_metric) == len(operator)
        and hermitian
    ):
        raise ValueError("disk series must have equal positive lengths")
    order = len(operator) - 1
    dimension = operator[0].rows
    defect = [coefficient[:, 0] for coefficient in hermitian]
    stein_metric = lyapunov_series(operator, defect, order)
    lower = generalized_endpoint_polynomial_series(
        stein_metric,
        disk_metric,
        endpoint=0,
        initial_value=sp.Rational(1, 2),
        order=order,
    )
    upper = generalized_endpoint_polynomial_series(
        stein_metric,
        disk_metric,
        endpoint=dimension - 1,
        initial_value=sp.Integer(2),
        order=order,
    )
    return [
        sp.factor(upper[degree] - 4 * lower[degree])
        for degree in range(order + 1)
    ]


def make_record(
    direction: tuple[sp.Expr, ...],
    direction_kind: str,
    order: int,
    expected_leading_degree: int,
    expected_leading_delta: sp.Expr | None = None,
) -> FullDiskBaseJetRecord:
    """Build and verify one exact base-jet record."""

    delta = endpoint_delta_series(direction, order)
    nonzero = [
        degree
        for degree, coefficient in enumerate(delta)
        if coefficient != 0
    ]
    if nonzero != [expected_leading_degree]:
        raise RuntimeError(
            f"unexpected delta support for {direction_kind}: {delta}"
        )
    leading = delta[expected_leading_degree]
    identity_verified = (
        expected_leading_delta is None
        or sp.expand(leading - expected_leading_delta) == 0
    )
    if not identity_verified:
        raise RuntimeError(
            f"leading identity failed for {direction_kind}: {leading}"
        )
    return FullDiskBaseJetRecord(
        dimension=len(direction) + 1,
        length=len(direction),
        direction_kind=direction_kind,
        direction=tuple(str(value) for value in direction[1:]),
        series_order=order,
        vanishing_through_degree=expected_leading_degree - 1,
        leading_degree=expected_leading_degree,
        leading_delta=str(leading),
        identity_verified=identity_verified,
    )


def symbolic_length_four_record() -> FullDiskBaseJetRecord:
    """Verify the factorized real length-four sixth-order face."""

    left, middle, right = sp.symbols(
        "a_left a_middle a_right",
        real=True,
    )
    direction = (sp.Integer(0), left, middle, right)
    expected = (
        -sp.Rational(128, 9)
        * middle**2
        * (left - right) ** 2
        * (2 * left**2 + middle**2)
    )
    return make_record(
        direction,
        "symbolic_real",
        order=6,
        expected_leading_degree=6,
        expected_leading_delta=expected,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-length", type=int, default=3)
    parser.add_argument("--maximum-length", type=int, default=6)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run the exact recentered base-jet audit."""

    args = parse_args()
    if (
        args.minimum_length < 3
        or args.maximum_length < args.minimum_length
    ):
        raise ValueError("require 3 <= minimum length <= maximum length")
    records = []
    for length in range(
        args.minimum_length,
        args.maximum_length + 1,
    ):
        leading_degree = 8 if length == 3 else 6
        records.append(
            make_record(
                deterministic_disk_direction(length),
                "deterministic_gaussian_rational",
                order=leading_degree,
                expected_leading_degree=leading_degree,
            )
        )
    if args.minimum_length <= 4 <= args.maximum_length:
        terminal_direction = (
            sp.Integer(0),
            sp.Rational(1, 5),
            sp.Integer(0),
            sp.Rational(1, 11),
        )
        records.append(
            make_record(
                terminal_direction,
                "terminal_only_rational",
                order=8,
                expected_leading_degree=8,
            )
        )
        records.append(symbolic_length_four_record())

    args.output.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        json.dumps(asdict(record), sort_keys=True)
        for record in records
    ]
    args.output.write_text(
        "".join(f"{line}\n" for line in lines),
        encoding="utf-8",
    )
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
