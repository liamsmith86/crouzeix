#!/usr/bin/env python3
"""Audit reciprocal-reversal duality for rank-one Crabb metrics.

There are two independent exact checks.

1.  For an arbitrary rational defect of the nilpotent persymmetric
    shift, reversal of the inverse Stein Gramian again satisfies a
    rank-one Stein equation.
2.  On the grade-two reflected Crabb path, the locally stationary
    defect through jet four makes

        J P(epsilon)^(-1) J = alpha(epsilon) P(epsilon)

    through the same four jets.  The first residual occurs at jet five,
    exactly where the omitted fifth optimizer coefficient can enter.

The second calculation checks a consequence of the uniqueness argument
in the accompanying proof; finite jets are not its proof.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Sequence

import sympy as sp

from crabb_circular_normal_series import (
    optimized_defect_jets,
    physical_reflected_path,
)
from general_crabb_weighted_series import inverse_riemann_series
from rank_one_stein_series import stein_gramian_series


DIMENSION = 5
GRADE = 2
JET_COUNT = 4


@dataclass(frozen=True)
class ReciprocalReversalRecord:
    """One exact reciprocal-reversal audit record."""

    rational_dual_rank: int
    rational_dual_stein_verified: bool
    dimension: int
    reflected_grade: int
    optimized_jet_count: int
    scalar_coefficients: tuple[str, ...]
    identity_through_optimized_jet: bool
    first_omitted_degree: int
    first_omitted_residual_nonzero_entries: int


def reversal_matrix(dimension: int) -> sp.Matrix:
    """Return coordinate reversal."""

    reversal = sp.zeros(dimension)
    for index in range(dimension):
        reversal[index, dimension - 1 - index] = 1
    return reversal


def inverse_matrix_series(
    coefficients: Sequence[sp.Matrix],
) -> list[sp.Matrix]:
    """Invert an ordinary matrix power series."""

    dimension = coefficients[0].rows
    inverse = [coefficients[0].inv()]
    for degree in range(1, len(coefficients)):
        convolution = sum(
            (
                coefficients[source_degree]
                * inverse[degree - source_degree]
                for source_degree in range(1, degree + 1)
            ),
            sp.zeros(dimension),
        )
        inverse.append(sp.simplify(-inverse[0] * convolution))
    return inverse


def scalar_series_quotient(
    numerator: Sequence[sp.Expr],
    denominator: Sequence[sp.Expr],
) -> list[sp.Expr]:
    """Divide two scalar power series with nonzero denominator constant."""

    quotient: list[sp.Expr] = []
    for degree in range(len(numerator)):
        known = sum(
            (
                quotient[source_degree]
                * denominator[degree - source_degree]
                for source_degree in range(degree)
            ),
            sp.Integer(0),
        )
        quotient.append(
            sp.simplify(
                (numerator[degree] - known) / denominator[0]
            )
        )
    return quotient


def rational_dual_audit() -> tuple[int, bool]:
    """Verify rank-one inverse duality away from the Crabb defect."""

    dimension = 4
    shift = sp.zeros(dimension)
    for index in range(dimension - 1):
        shift[index, index + 1] = 1
    defect = sp.Matrix(
        [
            1,
            sp.Rational(1, 3),
            sp.Rational(-2, 5),
            sp.Rational(1, 7),
        ]
    )
    gramian = sum(
        (
            (shift.T**power)
            * defect
            * defect.T
            * (shift**power)
            for power in range(dimension)
        ),
        sp.zeros(dimension),
    )
    inverse_defect = sp.simplify(
        gramian.inv()
        - shift * gramian.inv() * shift.T
    )
    reversal = reversal_matrix(dimension)
    dual_gramian = reversal * gramian.inv() * reversal
    dual_forcing = reversal * inverse_defect * reversal
    dual_residual = sp.simplify(
        dual_gramian
        - shift.T * dual_gramian * shift
        - dual_forcing
    )
    rank = inverse_defect.rank()
    verified = dual_residual == sp.zeros(dimension)
    if rank != 1 or not verified:
        raise AssertionError("the rational inverse-Stein duality failed")
    return rank, verified


def stationary_jet_audit() -> tuple[tuple[str, ...], int]:
    """Verify self-duality through every available optimizer jet."""

    epsilon = sp.symbols("epsilon", real=True)
    work_order = 2 * JET_COUNT
    path = physical_reflected_path(
        dimension=DIMENSION,
        equality_grade=GRADE,
        strong_parameter=sp.Integer(0),
        strong_direction=sp.zeros(DIMENSION),
        strong_degree=1,
        order=work_order,
    )
    _, operator = inverse_riemann_series(path, work_order)
    defect = optimized_defect_jets(
        operator,
        epsilon,
        jet_count=JET_COUNT,
    )

    comparison_order = JET_COUNT + 1
    gramian = stein_gramian_series(
        operator,
        defect,
        epsilon,
        comparison_order,
    )
    inverse = inverse_matrix_series(gramian)
    reversal = reversal_matrix(DIMENSION)
    reversed_inverse = [
        reversal * coefficient * reversal
        for coefficient in inverse
    ]
    scalar = scalar_series_quotient(
        [coefficient[0, 0] for coefficient in reversed_inverse],
        [coefficient[0, 0] for coefficient in gramian],
    )

    residuals = []
    for degree in range(comparison_order + 1):
        residual = sp.simplify(
            reversed_inverse[degree]
            - sum(
                (
                    scalar[source_degree]
                    * gramian[degree - source_degree]
                    for source_degree in range(degree + 1)
                ),
                sp.zeros(DIMENSION),
            )
        )
        residuals.append(residual)

    if any(
        residual != sp.zeros(DIMENSION)
        for residual in residuals[: JET_COUNT + 1]
    ):
        raise AssertionError(
            "reciprocal reversal failed within the optimized jet"
        )
    omitted_residual = residuals[JET_COUNT + 1]
    omitted_nonzero = sum(
        omitted_residual[row, column] != 0
        for row in range(DIMENSION)
        for column in range(DIMENSION)
    )
    if omitted_nonzero == 0:
        raise AssertionError(
            "the truncation boundary was not exposed by the checker"
        )
    return (
        tuple(str(sp.factor(value)) for value in scalar),
        omitted_nonzero,
    )


def make_record() -> ReciprocalReversalRecord:
    """Run both exact audits."""

    dual_rank, dual_verified = rational_dual_audit()
    scalar, omitted_nonzero = stationary_jet_audit()
    return ReciprocalReversalRecord(
        rational_dual_rank=dual_rank,
        rational_dual_stein_verified=dual_verified,
        dimension=DIMENSION,
        reflected_grade=GRADE,
        optimized_jet_count=JET_COUNT,
        scalar_coefficients=scalar,
        identity_through_optimized_jet=True,
        first_omitted_degree=JET_COUNT + 1,
        first_omitted_residual_nonzero_entries=omitted_nonzero,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run the deterministic exact audit."""

    args = parse_args()
    record = make_record()
    line = json.dumps(asdict(record), sort_keys=True)
    print(line, flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(line + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
