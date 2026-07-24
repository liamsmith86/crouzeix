#!/usr/bin/env python3
"""Regenerate the grade-three circular-normal leading cancellation.

The eligible circular normal is mode ``p-2`` (not the grade-two mode
``p-1``).  The exact calculation retains the final rank-one-defect jet,
lifts both Stein endpoints through weight eight, and verifies that the
two endpoint shifts occur in the ratio ``1:4``.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from crabb_circular_normal_series import (
    endpoint_condition_coefficients,
    optimized_defect_jets,
    physical_reflected_path,
    real_circular_normal_direction,
)
from general_crabb_weighted_series import inverse_riemann_series


SERIES_ORDER = 8
EQUALITY_GRADE = 3


@dataclass(frozen=True)
class GradeThreeSelectionRecord:
    """One exact grade-three normal-selection audit."""

    dimension: int
    length: int
    independent_amplitudes: bool
    flat_epsilon_eight_coefficient: str
    predicted_flat_epsilon_eight_coefficient: str
    normal_cross_coefficient: str
    lower_endpoint_cross_coefficient: str
    upper_endpoint_cross_coefficient: str


def make_record(
    dimension: int,
    independent_amplitudes: bool,
) -> GradeThreeSelectionRecord:
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
        dimension - 2,
    )
    path = physical_reflected_path(
        dimension=dimension,
        equality_grade=EQUALITY_GRADE,
        strong_parameter=strong,
        strong_direction=strong_direction,
        strong_degree=4,
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
        jet_count=4,
    )
    lower, upper, ratio = endpoint_condition_coefficients(
        operator,
        defect,
        epsilon,
        SERIES_ORDER,
    )

    for degree in range(SERIES_ORDER):
        if sp.diff(ratio[degree], strong) != 0:
            raise AssertionError(
                "the normal cross appeared below weight eight"
            )

    flat = sp.factor(ratio[8].subs(strong, 0))
    predicted_flat = -64 * amplitude**2 * ellipse**6
    lower_cross = sp.factor(
        sp.diff(lower[8], strong).subs(strong, 0)
    )
    upper_cross = sp.factor(
        sp.diff(upper[8], strong).subs(strong, 0)
    )
    normal_cross = sp.factor(
        sp.diff(ratio[8], strong).subs(strong, 0)
    )

    if sp.simplify(flat - predicted_flat) != 0:
        raise AssertionError("the grade-three flat face changed")
    if normal_cross != 0:
        raise AssertionError("the grade-three leading cross is nonzero")
    if sp.simplify(upper_cross - 4 * lower_cross) != 0:
        raise AssertionError("the endpoint shifts do not cancel")

    return GradeThreeSelectionRecord(
        dimension=dimension,
        length=dimension - 1,
        independent_amplitudes=independent_amplitudes,
        flat_epsilon_eight_coefficient=str(flat),
        predicted_flat_epsilon_eight_coefficient=str(
            sp.factor(predicted_flat)
        ),
        normal_cross_coefficient=str(normal_cross),
        lower_endpoint_cross_coefficient=str(lower_cross),
        upper_endpoint_cross_coefficient=str(upper_cross),
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=7)
    parser.add_argument("--maximum-size", type=int, default=8)
    parser.add_argument(
        "--bivariate-maximum-size",
        type=int,
        default=7,
        help=(
            "keep amplitude and ellipse symbolic through this size; "
            "use 0 to disable"
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run the exact grade-three regeneration."""

    args = parse_args()
    if args.minimum_size < 7 or args.maximum_size < args.minimum_size:
        raise ValueError("size range must satisfy 7 <= minimum <= maximum")
    if args.bivariate_maximum_size < 0:
        raise ValueError("bivariate maximum must be nonnegative")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as output:
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
            output.write(line + "\n")


if __name__ == "__main__":
    main()
