#!/usr/bin/env python3
"""Regenerate the universal leading-Hardy-residual identities.

Let an analytic full-disk path start at the Crabb matrix, and suppose
its projected Hardy residual first appears in order ``m`` as a matrix
``F_m``.  This checker tests the order-independent identities

    -2 [s^(2m)] delta = 4 ||F_m||_F^2

and

    G_(L,k)^(m)
      = 4 (4k - 1) / L^2
        sum_(a+b=L-k) (b-a) skew(F_m)_(a,b).

It covers genuinely different strata: arbitrary first- and
second-order defects, the cubic recentered face, a weighted pair
fourth-order face, and an exactly canceled fourth-order face whose
first residual is fifth order.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from crabb_circular_normal_quadratic_exact import (
    characteristic_series,
    crabb_normal_directions,
    directional_gradient_series,
    disk_model_from_hermitian_series,
)
from crabb_disk_toeplitz_quartic import extend
from crabb_full_disk_base_jet_exact import (
    endpoint_delta_from_disk_series,
    toeplitz_direction,
)
from crabb_full_disk_correction_isometry import plucker_correction
from crabb_full_disk_cubic_normal_exact import normal_curvature
from crabb_full_disk_hardy_geometry import (
    skew_part,
    weighted_projection,
)
from crabb_full_disk_kernel_eighth import pair_direction
from crabb_full_disk_pair_blowup import (
    recentered_pair_blowup_series,
)
from crabb_full_disk_sixth_hardy_factor import (
    hardy_residual_from_disk_series,
    matrix_is_zero,
)


@dataclass(frozen=True)
class LeadingResidualRecord:
    """One exact first-residual audit."""

    dimension: int
    length: int
    path_kind: str
    leading_residual_order: int
    lower_residual_orders_zero: bool
    lower_endpoint_orders_zero: bool
    lower_response_orders_zero: bool
    residual_skew_symmetric: bool
    endpoint_deficit: str
    residual_frobenius_square: str
    endpoint_factor_verified: bool
    active_normal_modes: tuple[int, ...]
    response_projection_verified: bool
    conjugate_projection_verified: bool
    flux_gain: str
    actual_normal_gain: str
    flux_margin: str
    actual_schur_margin: str
    all_identities_verified: bool


def deterministic_hermitian(length: int, seed: int) -> sp.Matrix:
    """Return one dense exact Hermitian coefficient."""

    matrix = sp.zeros(length)
    for row in range(length):
        matrix[row, row] = sp.Rational(
            (seed + 2) * (row + 1),
            (101 + 2 * seed) * length,
        )
        for column in range(row + 1, length):
            value = (
                sp.Rational(
                    (seed + row + 1) * (column + 2),
                    (107 + 3 * seed) * length,
                )
                + sp.I
                * sp.Rational(
                    seed + row + column + 2,
                    (109 + 5 * seed) * length,
                )
            )
            matrix[row, column] = value
            matrix[column, row] = sp.conjugate(value)
    return matrix


def build_disk_series(
    coefficients: tuple[sp.Matrix, ...],
    order: int,
) -> tuple[
    tuple[sp.Matrix, ...],
    tuple[sp.Matrix, ...],
    tuple[sp.Matrix, ...],
]:
    """Extend an ``L``-square Hermitian series and build its disk model."""

    if len(coefficients) > order + 1:
        raise ValueError("too many Hermitian coefficients for the order")
    length = coefficients[0].rows
    dimension = length + 1
    hermitian = (
        *(extend(coefficient) for coefficient in coefficients),
        *[
            sp.zeros(dimension)
            for _ in range(order + 1 - len(coefficients))
        ],
    )
    operator, metric = disk_model_from_hermitian_series(
        list(hermitian)
    )
    return hermitian, tuple(metric), tuple(operator)


def audit_path(
    path_kind: str,
    hermitian: tuple[sp.Matrix, ...],
    metric: tuple[sp.Matrix, ...],
    operator: tuple[sp.Matrix, ...],
    leading_order: int,
) -> LeadingResidualRecord:
    """Audit the endpoint and every true-normal leading coefficient."""

    if len(operator) <= 2 * leading_order:
        raise ValueError("the path must be expanded through order 2m")
    length = operator[0].rows - 1
    residuals, support_verified = hardy_residual_from_disk_series(
        hermitian[: leading_order + 1],
        operator[: leading_order + 1],
    )
    lower_residual_zero = all(
        matrix_is_zero(residual)
        for residual in residuals[: leading_order - 1]
    )
    leading_residual = residuals[leading_order - 1]
    if matrix_is_zero(leading_residual):
        raise RuntimeError(f"{path_kind} has no residual in expected order")
    skew_residual = skew_part(leading_residual)

    endpoint = endpoint_delta_from_disk_series(
        hermitian,
        metric,
        operator,
    )
    lower_endpoint_zero = all(
        coefficient == 0
        for coefficient in endpoint[: 2 * leading_order]
    )
    endpoint_deficit = sp.factor(
        -2 * endpoint[2 * leading_order]
    )
    residual_norm = sp.factor(
        sum(
            (
                entry * sp.conjugate(entry)
                for entry in leading_residual
            ),
            sp.Integer(0),
        )
    )
    endpoint_factor = (
        sp.simplify(endpoint_deficit - 4 * residual_norm) == 0
    )

    characteristic = characteristic_series(
        list(operator[: leading_order + 1])
    )
    lower_response_zero = True
    response_verified = True
    conjugate_verified = True
    active_modes = []
    flux_gain = sp.Integer(0)
    actual_gain = sp.Integer(0)
    for mode in range(3, length + 2):
        responses = []
        for normal in crabb_normal_directions([metric[0]], mode):
            response = directional_gradient_series(
                list(operator[: leading_order + 1]),
                list(metric[: leading_order + 1]),
                characteristic,
                normal,
            )
            lower_response_zero &= all(
                sp.simplify(coefficient) == 0
                for coefficient in response[:leading_order]
            )
            responses.append(sp.factor(response[leading_order]))
        complex_response = sp.factor(
            responses[0] + sp.I * responses[1]
        )
        lower, upper = weighted_projection(
            skew_residual,
            length,
            mode,
        )
        predicted = sp.factor(
            sp.Rational(4 * (4 * mode - 1), length**2) * lower
        )
        response_verified &= (
            sp.simplify(complex_response - predicted) == 0
        )
        conjugate_verified &= (
            sp.simplify(upper - sp.conjugate(lower)) == 0
        )
        if complex_response != 0:
            active_modes.append(mode)
            actual_gain += sum(
                response**2 for response in responses
            ) / (4 * normal_curvature(length, mode))
        if 3 <= mode <= length - 3:
            flux_gain += (
                16
                * lower
                * sp.conjugate(lower)
                / sp.binomial(length - mode, 3)
            )

    flux_gain = sp.factor(sp.simplify(flux_gain))
    actual_gain = sp.factor(sp.simplify(actual_gain))
    flux_margin = sp.factor(endpoint_deficit - flux_gain)
    actual_margin = sp.factor(endpoint_deficit - actual_gain)
    verified = bool(
        support_verified
        and lower_residual_zero
        and lower_endpoint_zero
        and lower_response_zero
        and endpoint_factor
        and response_verified
        and conjugate_verified
        and flux_margin >= 0
        and actual_margin > 0
    )
    if not verified:
        raise RuntimeError(f"leading-residual audit failed for {path_kind}")
    return LeadingResidualRecord(
        dimension=length + 1,
        length=length,
        path_kind=path_kind,
        leading_residual_order=leading_order,
        lower_residual_orders_zero=lower_residual_zero,
        lower_endpoint_orders_zero=lower_endpoint_zero,
        lower_response_orders_zero=lower_response_zero,
        residual_skew_symmetric=matrix_is_zero(
            leading_residual + leading_residual.T
        ),
        endpoint_deficit=str(endpoint_deficit),
        residual_frobenius_square=str(residual_norm),
        endpoint_factor_verified=endpoint_factor,
        active_normal_modes=tuple(active_modes),
        response_projection_verified=response_verified,
        conjugate_projection_verified=conjugate_verified,
        flux_gain=str(flux_gain),
        actual_normal_gain=str(actual_gain),
        flux_margin=str(flux_margin),
        actual_schur_margin=str(actual_margin),
        all_identities_verified=verified,
    )


def audit_suite() -> list[LeadingResidualRecord]:
    """Build representative paths with first residual orders one to five."""

    length = 6
    base = sp.eye(length) / 2
    order_one = build_disk_series(
        (base, deterministic_hermitian(length, 1)),
        order=2,
    )

    toeplitz = (
        sp.Integer(0),
        *[
            sp.Rational(index + 1, 50 + 3 * index)
            + sp.I * sp.Rational(index, 70 + 2 * index)
            for index in range(1, length)
        ],
    )
    order_two = build_disk_series(
        (
            base,
            toeplitz_direction(toeplitz),
            deterministic_hermitian(length, 2),
        ),
        order=4,
    )
    order_three = build_disk_series(
        (
            base,
            toeplitz_direction(toeplitz),
            plucker_correction(toeplitz),
        ),
        order=6,
    )

    pair = pair_direction(
        length,
        0,
        sp.Rational(2, 19) + sp.I * sp.Rational(1, 23),
        sp.Rational(1, 17) - sp.I * sp.Rational(2, 29),
    )
    transverse = (
        sp.Integer(0),
        *[
            sp.Rational(2 * index + 1, 41 * length)
            + sp.I * sp.Rational(index + 2, 47 * length)
            for index in range(1, length)
        ],
    )
    order_four = recentered_pair_blowup_series(
        pair,
        transverse,
        order=8,
    )

    canceled_pair = pair_direction(
        length,
        0,
        sp.Rational(1, 10),
        sp.Rational(1, 20),
    )
    canceling_transverse = (
        sp.Integer(0),
        sp.Integer(0),
        sp.Rational(1, 120),
        sp.Integer(0),
        sp.Rational(1, 240),
        sp.Integer(0),
    )
    order_five = recentered_pair_blowup_series(
        canceled_pair,
        canceling_transverse,
        order=10,
    )

    paths = (
        ("generic_first_order", order_one, 1),
        ("toeplitz_then_generic_second_order", order_two, 2),
        ("plucker_recentered_cubic", order_three, 3),
        ("weighted_pair_fourth_order", order_four, 4),
        ("canceled_pair_fifth_order", order_five, 5),
    )
    records = []
    for path_kind, series, leading_order in paths:
        record = audit_path(
            path_kind,
            *series,
            leading_order=leading_order,
        )
        records.append(record)
        print(
            f"verified leading residual order {leading_order}: "
            f"{path_kind}",
            flush=True,
        )
    return records


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run and serialize the complete leading-residual audit."""

    args = parse_args()
    records = audit_suite()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    lines = [json.dumps(asdict(record), sort_keys=True) for record in records]
    args.output.write_text(
        "".join(f"{line}\n" for line in lines),
        encoding="utf-8",
    )
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
