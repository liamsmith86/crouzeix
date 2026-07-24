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


def inverse_square_root_series(
    constant: sp.Matrix,
    tangent: sp.Matrix,
    order: int,
) -> list[sp.Matrix]:
    """Solve ``G(epsilon) K(epsilon) G(epsilon) = I`` with symmetric ``G``."""

    dimension = constant.rows
    square_roots = [
        sp.sqrt(constant[index, index])
        for index in range(dimension)
    ]
    coefficients = [
        sp.diag(*[1 / value for value in square_roots])
    ]
    gramian = [
        constant,
        tangent,
        *[sp.zeros(dimension) for _ in range(order - 1)],
    ]

    for degree in range(1, order + 1):
        known = sp.zeros(dimension)
        for left_degree in range(degree + 1):
            for middle_degree in range(degree - left_degree + 1):
                right_degree = degree - left_degree - middle_degree
                is_left_unknown = (
                    left_degree == degree
                    and middle_degree == 0
                    and right_degree == 0
                )
                is_right_unknown = (
                    left_degree == 0
                    and middle_degree == 0
                    and right_degree == degree
                )
                if is_left_unknown or is_right_unknown:
                    continue
                if (
                    left_degree < len(coefficients)
                    and right_degree < len(coefficients)
                ):
                    known += (
                        coefficients[left_degree]
                        * gramian[middle_degree]
                        * coefficients[right_degree]
                    )

        coefficient = sp.zeros(dimension)
        for row in range(dimension):
            for column in range(dimension):
                coefficient[row, column] = sp.simplify(
                    -known[row, column]
                    / (square_roots[row] + square_roots[column])
                )
        coefficients.append(coefficient)

    for degree in range(order + 1):
        residual = triple_series_coefficient(
            coefficients,
            gramian,
            coefficients,
            degree,
        )
        target = sp.eye(dimension) if degree == 0 else sp.zeros(dimension)
        if residual != target:
            raise AssertionError(
                "the inverse-square-root recurrence did not regenerate"
            )
    return coefficients


def triple_series_coefficient(
    left: Sequence[sp.Matrix],
    middle: Sequence[sp.Matrix],
    right: Sequence[sp.Matrix],
    degree: int,
) -> sp.Matrix:
    """Return one coefficient of a product of three matrix series."""

    dimension = left[0].rows
    result = sp.zeros(dimension)
    for left_degree in range(degree + 1):
        for middle_degree in range(degree - left_degree + 1):
            right_degree = degree - left_degree - middle_degree
            if (
                left_degree < len(left)
                and middle_degree < len(middle)
                and right_degree < len(right)
            ):
                result += (
                    left[left_degree]
                    * middle[middle_degree]
                    * right[right_degree]
                )
    return sp.simplify(result)


def physical_weighted_path(
    dimension: int,
    strong_parameter: sp.Symbol,
    strong_direction: sp.Matrix | None = None,
) -> list[sp.Matrix]:
    """Construct the equality/ellipse/bottom-normal path through order four."""

    length = dimension - 1
    shift = sp.zeros(dimension)
    for index in range(length):
        shift[index, index + 1] = 1

    base_toeplitz = sp.zeros(dimension)
    tangent_toeplitz = sp.zeros(dimension)
    for index in range(length):
        base_toeplitz[index, index] = sp.Rational(1, 2)
    for offset in {1, length - 1}:
        for row in range(length - offset):
            tangent_toeplitz[row, row + offset] = 1
            tangent_toeplitz[row + offset, row] = 1

    base_coordinate_gramian = (
        base_toeplitz
        + shift.T * base_toeplitz * shift
    )
    tangent_coordinate_gramian = (
        tangent_toeplitz
        + shift.T * tangent_toeplitz * shift
    )
    inverse_square_root = inverse_square_root_series(
        base_coordinate_gramian,
        tangent_coordinate_gramian,
        SERIES_ORDER,
    )
    toeplitz_series = [
        base_toeplitz,
        tangent_toeplitz,
        *[sp.zeros(dimension) for _ in range(SERIES_ORDER - 1)],
    ]
    shifted_inverse = [
        shift * coefficient
        for coefficient in inverse_square_root
    ]
    disk_operator = [
        2
        * triple_series_coefficient(
            inverse_square_root,
            toeplitz_series,
            shifted_inverse,
            degree,
        )
        for degree in range(SERIES_ORDER + 1)
    ]

    if strong_direction is None:
        strong_direction = sp.zeros(dimension)
        strong_direction[length, 0] = 1
    if strong_direction.shape != (dimension, dimension):
        raise ValueError("the strong direction has the wrong shape")

    path = []
    for degree in range(SERIES_ORDER + 1):
        coefficient = disk_operator[degree]
        if degree >= 1:
            coefficient += disk_operator[degree - 1].T
        if degree == 2:
            coefficient += strong_parameter * strong_direction
        path.append(sp.simplify(coefficient))
    return path


def optimized_flat_defect(
    operator: Sequence[sp.Matrix],
    epsilon: sp.Symbol,
) -> sp.Matrix:
    """Solve the first two defect stationarity equations on the flat path."""

    dimension = operator[0].rows
    first_variables = sp.symbols(
        f"first_defect_1:{dimension}",
        real=True,
    )
    second_variables = sp.symbols(
        f"second_defect_1:{dimension}",
        real=True,
    )
    base = sp.zeros(dimension, 1)
    base[0] = 1
    first = sp.Matrix([0, *first_variables])
    second = sp.Matrix([0, *second_variables])
    defect = base + epsilon * first + epsilon**2 * second

    condition = diagonal_gramian_condition_series(
        operator,
        defect,
        epsilon,
        SERIES_ORDER,
    )
    quadratic = sp.expand(condition).coeff(epsilon, 2)
    first_solutions = sp.solve(
        [
            sp.diff(quadratic, variable)
            for variable in first_variables
        ],
        first_variables,
        dict=True,
    )
    if len(first_solutions) != 1:
        raise RuntimeError("the first defect jet was not uniquely stationary")

    quartic = sp.expand(
        condition.subs(first_solutions[0])
    ).coeff(epsilon, 4)
    second_solutions = sp.solve(
        [
            sp.diff(quartic, variable)
            for variable in second_variables
        ],
        second_variables,
        dict=True,
    )
    if len(second_solutions) != 1:
        raise RuntimeError("the second defect jet was not uniquely stationary")
    return sp.simplify(
        defect.subs(
            {
                **first_solutions[0],
                **second_solutions[0],
            }
        )
    )


def real_circular_normal_direction(
    dimension: int,
    mode: int,
) -> sp.Matrix:
    """Return the real Riesz representative of one disk-normal support mode."""

    if not 2 <= mode <= dimension:
        raise ValueError("a circular normal mode must lie between 2 and p")
    weights = [
        1 / sp.sqrt(2) if index in (0, dimension - 1) else sp.Integer(1)
        for index in range(dimension)
    ]
    direction = sp.zeros(dimension)
    for row in range(dimension):
        for column in range(dimension):
            grade = column - row - 1
            if grade in (mode, -mode):
                direction[row, column] = weights[row] * weights[column]
    return direction


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
