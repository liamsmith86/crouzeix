#!/usr/bin/env python3
"""Regenerate the first circular-normal completed square at a Crabb block.

Let ``L = dimension - 1``.  Give a phase-palindromic disk coefficient and
the ellipse parameter the same bookkeeping weight ``epsilon`` and give the
bottom circular-normal entry ``E[L, 0]`` weight ``epsilon**2``.  This script
derives, without floating-point arithmetic,

    D_strong [epsilon**4] Gamma = -8 (5 L - 1) / L.

Together with the already-proved L65 bottom square and the L134 flat face,
this leaves the strictly positive residual

    32 (L - 1) (2 L**2 - L + 3)
    --------------------------------.
       L (L**2 + 36 L - 13)

The finite regenerations audit the Riemann/Stein recurrence.  The all-size
step is the closed endpoint calculation recorded in the accompanying proof,
not extrapolation from the tested dimensions.
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
    real_circular_normal_direction,
)
from general_crabb_weighted_series import inverse_riemann_series
from rank_one_stein_series import (
    diagonal_gramian_condition_series,
    simple_diagonal_eigenvalue_coefficients,
    stein_gramian_series,
)


SERIES_ORDER = 4


@dataclass(frozen=True)
class CircularNormalFaceRecord:
    """One exact completed-square audit record."""

    dimension: int
    length: int
    flat_epsilon_four_coefficient: str
    bottom_cross_coefficient: str
    predicted_bottom_cross_coefficient: str
    lower_endpoint_cross_coefficient: str
    predicted_lower_endpoint_cross_coefficient: str
    upper_endpoint_cross_coefficient: str
    predicted_upper_endpoint_cross_coefficient: str
    bottom_hessian_coefficient: str
    completed_residual: str
    predicted_completed_residual: str
    real_normal_mode_cross_coefficients: dict[str, str] | None


def physical_weighted_path(
    dimension: int,
    strong_parameter: sp.Symbol,
    strong_direction: sp.Matrix | None = None,
) -> list[sp.Matrix]:
    """Construct the equality/ellipse/bottom-normal path through order four."""

    length = dimension - 1
    if strong_direction is None:
        strong_direction = sp.zeros(dimension)
        strong_direction[length, 0] = 1
    return physical_reflected_path(
        dimension=dimension,
        equality_grade=1,
        strong_parameter=strong_parameter,
        strong_direction=strong_direction,
        strong_degree=2,
        order=SERIES_ORDER,
    )


def optimized_flat_defect(
    operator: Sequence[sp.Matrix],
    epsilon: sp.Symbol,
) -> sp.Matrix:
    """Solve the first two defect stationarity equations on the flat path."""

    return optimized_defect_jets(
        operator,
        epsilon,
        jet_count=2,
    )


def endpoint_cross_coefficients(
    operator: Sequence[sp.Matrix],
    defect: sp.Matrix,
    epsilon: sp.Symbol,
    strong: sp.Symbol,
) -> tuple[sp.Expr, sp.Expr]:
    """Return the two fourth-order endpoint derivatives in ``strong``."""

    gramian = stein_gramian_series(
        operator,
        defect,
        epsilon,
        SERIES_ORDER,
    )
    lower = simple_diagonal_eigenvalue_coefficients(
        gramian,
        endpoint=0,
        base_eigenvalue=1,
    )
    upper = simple_diagonal_eigenvalue_coefficients(
        gramian,
        endpoint=gramian[0].rows - 1,
        base_eigenvalue=4,
    )
    for degree in range(4):
        if (
            sp.diff(lower[degree], strong) != 0
            or sp.diff(upper[degree], strong) != 0
        ):
            raise AssertionError(
                "the strong insertion appeared below fourth order"
            )
    return (
        sp.factor(sp.diff(lower[4], strong).subs(strong, 0)),
        sp.factor(sp.diff(upper[4], strong).subs(strong, 0)),
    )


def normal_mode_cross_coefficients(
    dimension: int,
    flat_path: Sequence[sp.Matrix],
    defect: sp.Matrix,
    epsilon: sp.Symbol,
    strong: sp.Symbol,
) -> dict[str, str]:
    """Audit every real coercive circular-normal character in one size."""

    crosses: dict[str, str] = {}
    for mode in range(3, dimension + 1):
        path = [coefficient.copy() for coefficient in flat_path]
        path[2] += (
            strong
            * real_circular_normal_direction(dimension, mode)
        )
        _, operator = inverse_riemann_series(path, SERIES_ORDER)
        condition = diagonal_gramian_condition_series(
            operator,
            defect,
            epsilon,
            SERIES_ORDER,
        )
        cross = sp.factor(
            sp.diff(
                sp.expand(condition).coeff(epsilon, 4),
                strong,
            ).subs(strong, 0)
        )
        predicted = (
            sp.Rational(-4 * (5 * (dimension - 1) - 1), dimension - 1)
            if mode == dimension
            else sp.Integer(0)
        )
        if cross != predicted:
            raise AssertionError(
                f"dimension {dimension}: mode {mode} cross changed"
            )
        crosses[str(mode)] = str(cross)
    return crosses


def make_record(
    dimension: int,
    audit_normal_modes: bool,
) -> CircularNormalFaceRecord:
    """Regenerate and verify one exact size."""

    length = dimension - 1
    epsilon, strong = sp.symbols(
        "epsilon strong",
        real=True,
    )
    path = physical_weighted_path(dimension, strong)
    _, operator = inverse_riemann_series(path, SERIES_ORDER)
    flat_operator = [
        coefficient.subs(strong, 0)
        for coefficient in operator
    ]
    defect = optimized_flat_defect(flat_operator, epsilon)
    condition = diagonal_gramian_condition_series(
        operator,
        defect,
        epsilon,
        SERIES_ORDER,
    )
    quartic = sp.expand(condition).coeff(epsilon, 4)
    flat = sp.factor(quartic.subs(strong, 0))
    cross = sp.factor(sp.diff(quartic, strong).subs(strong, 0))
    lower_cross, upper_cross = endpoint_cross_coefficients(
        operator,
        defect,
        epsilon,
        strong,
    )
    predicted_cross = sp.Rational(-8 * (5 * length - 1), length)
    if length == 2:
        predicted_lower_cross = 8 * (sp.sqrt(2) + 2)
        predicted_upper_cross = 4 * (7 + 8 * sp.sqrt(2))
    elif length == 3:
        predicted_lower_cross = 8 * (1 + sp.sqrt(2))
        predicted_upper_cross = sp.Rational(16, 3) * (
            -1 + 6 * sp.sqrt(2)
        )
    else:
        predicted_lower_cross = 8 * (sp.sqrt(2) + 2)
        predicted_upper_cross = (
            24 + sp.Rational(8, length) + 32 * sp.sqrt(2)
        )
    bottom_hessian = sp.Rational(
        length**2 + 36 * length - 13,
        6 * length,
    )
    residual = sp.factor(
        64 - cross**2 / (4 * bottom_hessian)
    )
    predicted_residual = sp.factor(
        sp.Rational(
            32
            * (length - 1)
            * (2 * length**2 - length + 3),
            length * (length**2 + 36 * length - 13),
        )
    )

    expected_flat = -80 if length == 2 else -64
    if flat != expected_flat:
        raise AssertionError("the flat weighted coefficient changed")
    if cross != predicted_cross:
        raise AssertionError("the bottom cross coefficient changed")
    if sp.simplify(lower_cross - predicted_lower_cross) != 0:
        raise AssertionError("the lower endpoint cross coefficient changed")
    if sp.simplify(upper_cross - predicted_upper_cross) != 0:
        raise AssertionError("the upper endpoint cross coefficient changed")
    if sp.simplify(upper_cross - 4 * lower_cross - cross) != 0:
        raise AssertionError("the endpoint cross coefficients do not combine")
    if residual != predicted_residual or residual <= 0:
        raise AssertionError("the completed residual is not strictly positive")

    mode_crosses = None
    if audit_normal_modes:
        flat_path = physical_weighted_path(
            dimension,
            sp.Integer(0),
        )
        mode_crosses = normal_mode_cross_coefficients(
            dimension,
            flat_path,
            defect,
            epsilon,
            strong,
        )

    return CircularNormalFaceRecord(
        dimension=dimension,
        length=length,
        flat_epsilon_four_coefficient=str(flat),
        bottom_cross_coefficient=str(cross),
        predicted_bottom_cross_coefficient=str(predicted_cross),
        lower_endpoint_cross_coefficient=str(lower_cross),
        predicted_lower_endpoint_cross_coefficient=str(
            predicted_lower_cross
        ),
        upper_endpoint_cross_coefficient=str(upper_cross),
        predicted_upper_endpoint_cross_coefficient=str(
            predicted_upper_cross
        ),
        bottom_hessian_coefficient=str(bottom_hessian),
        completed_residual=str(residual),
        predicted_completed_residual=str(predicted_residual),
        real_normal_mode_cross_coefficients=mode_crosses,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=3)
    parser.add_argument("--maximum-size", type=int, default=6)
    parser.add_argument(
        "--selection-maximum-size",
        type=int,
        default=5,
        help=(
            "audit every real circular-normal mode through this size; "
            "use 0 to disable"
        ),
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the exact regeneration."""

    args = parse_args()
    if args.minimum_size < 3 or args.maximum_size < args.minimum_size:
        raise ValueError("size range must satisfy 3 <= minimum <= maximum")
    if args.selection_maximum_size < 0:
        raise ValueError("selection maximum must be nonnegative")
    lines = []
    for dimension in range(
        args.minimum_size,
        args.maximum_size + 1,
    ):
        record = make_record(
            dimension,
            audit_normal_modes=(
                3 <= dimension <= args.selection_maximum_size
            ),
        )
        line = json.dumps(asdict(record), sort_keys=True)
        lines.append(line)
        print(line, flush=True)
    if args.output is not None:
        args.output.write_text(
            "".join(f"{line}\n" for line in lines),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
