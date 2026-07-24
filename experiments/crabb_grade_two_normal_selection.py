#!/usr/bin/env python3
"""Regenerate the grade-two circular-normal leading cancellation.

Give the equality amplitude and ellipse parameter weight one.  The compact
grade-two coordinate then has weight three, as does its only eligible
coercive circular normal.  This checker proves that their apparent
weight-six cross is zero after exact rank-one defect optimization.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from crabb_circular_normal_series import (
    optimized_defect_jets,
    physical_reflected_path,
    real_circular_normal_direction,
)
from general_crabb_weighted_series import inverse_riemann_series
from rank_one_stein_series import (
    simple_diagonal_eigenvalue_coefficients,
    stein_gramian_series,
)


SERIES_ORDER = 6
EQUALITY_GRADE = 2


@dataclass(frozen=True)
class GradeTwoSelectionRecord:
    """One exact grade-two normal-selection audit."""

    dimension: int
    length: int
    independent_amplitudes: bool
    flat_epsilon_six_coefficient: str
    predicted_flat_epsilon_six_coefficient: str
    normal_cross_coefficient: str
    lower_endpoint_cross_coefficient: str
    upper_endpoint_cross_coefficient: str


def endpoint_condition_coefficients(
    operator: list[sp.Matrix],
    defect: sp.Matrix,
    epsilon: sp.Symbol,
) -> tuple[list[sp.Expr], list[sp.Expr], list[sp.Expr]]:
    """Return the two endpoint series and their quotient."""

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
    ratio = [sp.Integer(4), *[sp.Integer(0) for _ in range(SERIES_ORDER)]]
    for degree in range(1, SERIES_ORDER + 1):
        ratio[degree] = sp.expand(
            upper[degree]
            - sum(
                lower[source_degree] * ratio[degree - source_degree]
                for source_degree in range(1, degree + 1)
            )
        )
    return lower, upper, ratio


def make_record(
    dimension: int,
    independent_amplitudes: bool,
) -> GradeTwoSelectionRecord:
    """Regenerate and verify one exact dimension."""

    epsilon, strong = sp.symbols(
        "epsilon strong",
        real=True,
    )
    if independent_amplitudes:
        amplitude, ellipse = sp.symbols(
            "amplitude ellipse",
            real=True,
        )
    else:
        amplitude = ellipse = sp.Integer(1)

    strong_direction = real_circular_normal_direction(
        dimension,
        dimension - 1,
    )
    path = physical_reflected_path(
        dimension=dimension,
        equality_grade=EQUALITY_GRADE,
        strong_parameter=strong,
        strong_direction=strong_direction,
        strong_degree=3,
        order=SERIES_ORDER,
        amplitude=amplitude,
        ellipse=ellipse,
    )
    _, operator = inverse_riemann_series(path, SERIES_ORDER)
    flat_operator = [
        coefficient.subs(strong, 0)
        for coefficient in operator
    ]
    defect = optimized_defect_jets(
        flat_operator,
        epsilon,
        jet_count=3,
    )
    lower, upper, ratio = endpoint_condition_coefficients(
        operator,
        defect,
        epsilon,
    )

    for degree in range(SERIES_ORDER):
        if sp.diff(ratio[degree], strong) != 0:
            raise AssertionError(
                "the normal cross appeared below weight six"
            )

    flat = sp.factor(ratio[6].subs(strong, 0))
    predicted_flat = -64 * amplitude**2 * ellipse**4
    lower_cross = sp.factor(
        sp.diff(lower[6], strong).subs(strong, 0)
    )
    upper_cross = sp.factor(
        sp.diff(upper[6], strong).subs(strong, 0)
    )
    normal_cross = sp.factor(
        sp.diff(ratio[6], strong).subs(strong, 0)
    )
    length = dimension - 1

    if sp.simplify(flat - predicted_flat) != 0:
        raise AssertionError("the grade-two flat face changed")
    if normal_cross != 0:
        raise AssertionError("the grade-two leading normal cross is nonzero")
    if sp.simplify(upper_cross - 4 * lower_cross) != 0:
        raise AssertionError("the endpoint normal shifts do not cancel")

    if length >= 7 and (lower_cross != 0 or upper_cross != 0):
        raise AssertionError("the separated endpoint derivatives are nonzero")
    if length == 6:
        collision = -2 * ellipse * (
            amplitude**2
            + 2 * amplitude * ellipse
            - 2 * ellipse**2
        )
        if (
            sp.simplify(lower_cross - collision) != 0
            or sp.simplify(upper_cross - 4 * collision) != 0
        ):
            raise AssertionError("the length-six collision formula changed")

    return GradeTwoSelectionRecord(
        dimension=dimension,
        length=length,
        independent_amplitudes=independent_amplitudes,
        flat_epsilon_six_coefficient=str(flat),
        predicted_flat_epsilon_six_coefficient=str(
            sp.factor(predicted_flat)
        ),
        normal_cross_coefficient=str(normal_cross),
        lower_endpoint_cross_coefficient=str(lower_cross),
        upper_endpoint_cross_coefficient=str(upper_cross),
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=5)
    parser.add_argument("--maximum-size", type=int, default=9)
    parser.add_argument(
        "--bivariate-maximum-size",
        type=int,
        default=7,
        help=(
            "keep amplitude and ellipse symbolic through this size; "
            "use 0 to disable"
        ),
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the exact grade-two regeneration."""

    args = parse_args()
    if args.minimum_size < 5 or args.maximum_size < args.minimum_size:
        raise ValueError("size range must satisfy 5 <= minimum <= maximum")
    if args.bivariate_maximum_size < 0:
        raise ValueError("bivariate maximum must be nonnegative")

    output = None
    if args.output is not None:
        output = args.output.open("w", encoding="utf-8")
    try:
        for dimension in range(
            args.minimum_size,
            args.maximum_size + 1,
        ):
            record = make_record(
                dimension,
                independent_amplitudes=(
                    dimension <= args.bivariate_maximum_size
                ),
            )
            line = json.dumps(asdict(record), sort_keys=True)
            print(line, flush=True)
            if output is not None:
                output.write(line + "\n")
                output.flush()
    finally:
        if output is not None:
            output.close()


if __name__ == "__main__":
    main()
