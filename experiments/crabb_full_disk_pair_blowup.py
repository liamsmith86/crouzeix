#!/usr/bin/env python3
"""Regenerate the weighted blow-up of a sixth-kernel pair.

For a reversal-pair direction ``z0`` and an arbitrary transverse
direction ``v``, use the weighted path

    z(s) = s z0 + s**2 v,
    H(s) = I/2 + Z(z(s)) + B(z(s) wedge C z(s)).

For the checked directions, the first Hardy residual is order four.
This checker verifies that it is the pair fourth residual plus the
derivative of the cubic residual, that the eighth endpoint deficit is
four times its Frobenius norm, and that every quartic true-normal
response is its weighted anti-diagonal projection.  The same
Cauchy/null-lift argument as L184 proves nonnegativity of each checked
weighted eighth face.  Some special transverse directions cancel this
residual and must be continued to higher order; they are audited by
``crabb_full_disk_leading_residual.py``.
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
from crabb_full_disk_kernel_eighth import (
    closed_cubic_residual,
    pair_direction,
    predicted_pair_fourth_residual,
)
from crabb_full_disk_hardy_geometry import weighted_projection
from crabb_full_disk_sixth_hardy_factor import (
    hardy_residual_from_disk_series,
    matrix_is_zero,
)


@dataclass(frozen=True)
class PairBlowupRecord:
    """One exact weighted reversal-pair blow-up audit."""

    dimension: int
    length: int
    left_index: int
    right_index: int
    residual_lower_orders_zero: bool
    residual_fourth_skew_symmetric: bool
    weighted_residual_formula_verified: bool
    endpoint_eighth_deficit: str
    hardy_fourth_norm_square: str
    hardy_factor_verified: bool
    active_quartic_modes: tuple[int, ...]
    response_projection_verified: bool
    flux_only_gain: str
    actual_normal_gain: str
    actual_schur_margin: str
    strict_checked_direction: bool
    all_identities_verified: bool


def deterministic_transverse_direction(
    length: int,
    left: int,
) -> tuple[sp.Expr, ...]:
    """Return one dense exact transverse Toeplitz direction."""

    return (
        sp.Integer(0),
        *[
            sp.Rational(
                2 * offset + left + 1,
                (31 + left) * length,
            )
            + sp.I
            * sp.Rational(
                offset + 2 * left + 2,
                (37 + 2 * left) * length,
            )
            for offset in range(length - 1)
        ],
    )


def recentered_pair_blowup_series(
    pair: tuple[sp.Expr, ...],
    transverse: tuple[sp.Expr, ...],
    order: int,
) -> tuple[
    tuple[sp.Matrix, ...],
    tuple[sp.Matrix, ...],
    tuple[sp.Matrix, ...],
]:
    """Construct the exact weighted recentered disk-chart series."""

    if len(pair) != len(transverse):
        raise ValueError("pair and transverse directions must match")
    if order < 4:
        raise ValueError("the weighted blow-up requires order >= 4")
    length = len(pair)
    dimension = length + 1
    pair_correction = plucker_correction(pair)
    transverse_correction = plucker_correction(transverse)
    combined = tuple(
        pair[index] + transverse[index]
        for index in range(length)
    )
    mixed_correction = (
        plucker_correction(combined)
        - pair_correction
        - transverse_correction
    ).applyfunc(sp.expand)
    hermitian = (
        extend(sp.eye(length) / 2),
        extend(toeplitz_direction(pair)),
        extend(toeplitz_direction(transverse) + pair_correction),
        extend(mixed_correction),
        extend(transverse_correction),
        *[sp.zeros(dimension) for _ in range(order - 4)],
    )
    operator, metric = disk_model_from_hermitian_series(
        list(hermitian)
    )
    return hermitian, tuple(metric), tuple(operator)


def differentiated_cubic_residual(
    pair: tuple[sp.Expr, ...],
    transverse: tuple[sp.Expr, ...],
) -> sp.Matrix:
    """Differentiate the homogeneous cubic residual at the pair."""

    epsilon = sp.symbols("weighted_pair_epsilon", real=True)
    varied = tuple(
        pair[index] + epsilon * transverse[index]
        for index in range(len(pair))
    )
    return closed_cubic_residual(varied).applyfunc(
        lambda entry: sp.expand(entry).coeff(epsilon, 1)
    )


def audit_pair(length: int, left: int) -> PairBlowupRecord:
    """Audit one pair plane and one dense transverse direction."""

    right = length - 2 - left
    pair = pair_direction(
        length,
        left,
        sp.Rational(left + 2, 7)
        + sp.I * sp.Rational(left + 1, 11),
        sp.Rational(left + 3, 13)
        - sp.I * sp.Rational(left + 2, 17),
    )
    transverse = deterministic_transverse_direction(length, left)
    hermitian, metric, operator = recentered_pair_blowup_series(
        pair,
        transverse,
        order=8,
    )
    residuals, _ = hardy_residual_from_disk_series(
        hermitian[:5],
        operator[:5],
    )
    lower_zero = all(
        matrix_is_zero(residuals[degree])
        for degree in range(3)
    )
    fourth = residuals[3]
    fourth_skew = matrix_is_zero(fourth + fourth.T)
    predicted_fourth = (
        predicted_pair_fourth_residual(pair, left)
        + differentiated_cubic_residual(pair, transverse)
    )
    residual_formula = matrix_is_zero(fourth - predicted_fourth)

    endpoint = endpoint_delta_from_disk_series(
        hermitian,
        metric,
        operator,
    )
    base_deficit = sp.factor(-2 * endpoint[8])
    norm_square = sp.factor(
        sum(
            (
                entry * sp.conjugate(entry)
                for entry in fourth
            ),
            sp.Integer(0),
        )
    )
    hardy_factor = sp.simplify(base_deficit - 4 * norm_square) == 0

    characteristic = characteristic_series(list(operator[:5]))
    active_modes = []
    projection_verified = True
    actual_gain = sp.Integer(0)
    flux_gain = sp.Integer(0)
    for mode in range(3, length + 2):
        responses = []
        for normal in crabb_normal_directions([metric[0]], mode):
            response = directional_gradient_series(
                list(operator[:5]),
                list(metric[:5]),
                characteristic,
                normal,
            )
            responses.append(sp.factor(response[4]))
        complex_response = responses[0] + sp.I * responses[1]
        lower, upper = weighted_projection(fourth, length, mode)
        predicted_response = (
            sp.Rational(4 * (4 * mode - 1), length**2) * lower
        )
        projection_verified &= (
            sp.simplify(complex_response - predicted_response) == 0
            and sp.simplify(upper - sp.conjugate(lower)) == 0
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

    actual_gain = sp.factor(sp.simplify(actual_gain))
    flux_gain = sp.factor(sp.simplify(flux_gain))
    margin = sp.factor(sp.simplify(base_deficit - actual_gain))
    strict_face = bool(margin > 0)
    verified = bool(
        lower_zero
        and fourth_skew
        and residual_formula
        and hardy_factor
        and projection_verified
        and sp.simplify(base_deficit - flux_gain) >= 0
        and strict_face
    )
    if not verified:
        raise RuntimeError(
            f"weighted pair blow-up failed in length {length}, "
            f"pair {left},{right}"
        )
    return PairBlowupRecord(
        dimension=length + 1,
        length=length,
        left_index=left,
        right_index=right,
        residual_lower_orders_zero=lower_zero,
        residual_fourth_skew_symmetric=fourth_skew,
        weighted_residual_formula_verified=residual_formula,
        endpoint_eighth_deficit=str(base_deficit),
        hardy_fourth_norm_square=str(norm_square),
        hardy_factor_verified=hardy_factor,
        active_quartic_modes=tuple(active_modes),
        response_projection_verified=projection_verified,
        flux_only_gain=str(flux_gain),
        actual_normal_gain=str(actual_gain),
        actual_schur_margin=str(margin),
        strict_checked_direction=strict_face,
        all_identities_verified=verified,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-length", type=int, default=4)
    parser.add_argument("--maximum-length", type=int, default=8)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run and serialize every requested reversal-pair audit."""

    args = parse_args()
    if (
        args.minimum_length < 4
        or args.maximum_length < args.minimum_length
    ):
        raise ValueError("require 4 <= minimum length <= maximum length")
    records = []
    for length in range(args.minimum_length, args.maximum_length + 1):
        for left in range((length - 1) // 2):
            records.append(audit_pair(length, left))
            print(
                f"verified weighted pair length {length}, left {left}",
                flush=True,
            )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    lines = [json.dumps(asdict(record), sort_keys=True) for record in records]
    args.output.write_text(
        "".join(f"{line}\n" for line in lines),
        encoding="utf-8",
    )
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
