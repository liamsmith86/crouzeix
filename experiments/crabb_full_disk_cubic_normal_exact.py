#!/usr/bin/env python3
"""Audit the cubic true-normal response after L176 recentering.

For the exact disk path

    H(s) = I/2 + s Z(z) + s**2 B(z wedge J conjugate(z)),

L177 cancels the complete quadratic ambient response.  This checker
extends the independent characteristic-Blaschke/Riemann series engine
through cubic order, records the first surviving true circular-normal
mode, and compares its exact completed-square gain with the canonical
disk-base sixth-order deficit.
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
    deterministic_disk_direction,
    directional_gradient_series,
    disk_model_series,
)
from crabb_full_disk_correction_isometry import plucker_correction
from crabb_full_disk_base_jet_exact import endpoint_delta_series


@dataclass(frozen=True)
class CubicNormalRecord:
    """One exact dimension's recentered normal-response audit."""

    dimension: int
    length: int
    direction_index: int
    disk_direction: tuple[str, ...]
    normal_polarization_count: int
    active_cubic_modes: tuple[int, ...]
    cubic_responses: tuple[str, ...]
    quadratic_response_zero: bool
    cubic_response_zero: bool
    base_deficit_sixth_coefficient: str
    cubic_schur_gain_coefficient: str
    gain_over_base_deficit: str
    strict_sixth_order_face: bool
    length_six_formula_verified: bool


def normal_curvature(length: int, mode: int) -> sp.Expr:
    """Return L173's exact positive curvature on a support Riesz row."""

    remainder = length - mode
    flux = (
        (4 * mode - 1) ** 2
        * remainder
        * (remainder - 1)
        * (remainder - 2)
        / sp.Integer(24)
    )
    lifted_null = (
        2
        * mode
        * (mode - 1)
        * (mode - 2)
        * (sp.Rational(remainder) + sp.Rational(1, 4)) ** 2
        / 3
    )
    return sp.factor((flux + lifted_null) / length**4)


def length_six_cubic_response(
    direction: tuple[sp.Expr, ...],
) -> sp.Expr:
    """Return the first ``p=7`` cubic response in closed form."""

    coefficients = sp.Matrix(direction[1:])
    reversal = sp.Matrix(
        [
            sp.conjugate(coefficients[4 - index])
            for index in range(5)
        ]
    )
    plucker = (
        coefficients * reversal.T - reversal * coefficients.T
    )
    cubic = (
        6 * coefficients[3] * plucker[0, 3]
        - 5 * coefficients[4] * plucker[0, 2]
        + 2 * coefficients[3] * plucker[1, 2]
    )
    return sp.expand(-sp.Rational(22, 45) * cubic)


def audit_direction(
    direction: tuple[sp.Expr, ...],
    direction_index: int,
) -> CubicNormalRecord:
    """Audit every fixed Crabb true-normal polarization on one ray."""

    length = len(direction)
    correction = plucker_correction(direction)
    operator, metric = disk_model_series(
        direction,
        order=3,
        correction=correction,
    )
    characteristic = characteristic_series(operator)
    normal_count = 0
    quadratic_zero = True
    cubic_zero = True
    active_modes = []
    cubic_responses = []
    schur_gain = sp.Integer(0)
    for mode in range(3, length + 2):
        normal_directions = crabb_normal_directions([metric[0]], mode)
        mode_responses = []
        for normal in normal_directions:
            response = directional_gradient_series(
                operator,
                metric,
                characteristic,
                normal,
            )
            quadratic_zero &= sp.simplify(response[2]) == 0
            cubic_zero &= sp.simplify(response[3]) == 0
            mode_responses.append(sp.factor(response[3]))
            normal_count += 1
        if any(response != 0 for response in mode_responses):
            active_modes.append(mode)
            cubic_responses.extend(
                str(response) for response in mode_responses
            )
            curvature = normal_curvature(length, mode)
            schur_gain += sum(
                response**2 for response in mode_responses
            ) / (4 * curvature)

    if not quadratic_zero:
        raise RuntimeError(
            f"recentered normal response failed in length {length}: "
            f"quadratic={quadratic_zero}"
        )
    length_six_formula_verified = True
    if length == 6:
        actual_complex = (
            sp.sympify(cubic_responses[0])
            + sp.I * sp.sympify(cubic_responses[1])
        )
        length_six_formula_verified = (
            active_modes == [3]
            and sp.simplify(
                actual_complex - length_six_cubic_response(direction)
            )
            == 0
        )
        if not length_six_formula_verified:
            raise RuntimeError("the closed length-six response failed")
    base_delta = endpoint_delta_series(direction, order=6)[6]
    base_deficit = sp.factor(-2 * base_delta)
    schur_gain = sp.factor(schur_gain)
    if base_deficit == 0:
        ratio = sp.Integer(0) if schur_gain == 0 else sp.oo
    else:
        ratio = sp.factor(schur_gain / base_deficit)
    strict_face = bool(schur_gain < base_deficit) if base_deficit else bool(
        schur_gain == 0
    )
    if not strict_face:
        raise RuntimeError(
            f"the sixth-order Schur face failed in length {length}: "
            f"base={base_deficit}, gain={schur_gain}"
        )
    return CubicNormalRecord(
        dimension=length + 1,
        length=length,
        direction_index=direction_index,
        disk_direction=tuple(str(value) for value in direction[1:]),
        normal_polarization_count=normal_count,
        active_cubic_modes=tuple(active_modes),
        cubic_responses=tuple(cubic_responses),
        quadratic_response_zero=quadratic_zero,
        cubic_response_zero=cubic_zero,
        base_deficit_sixth_coefficient=str(base_deficit),
        cubic_schur_gain_coefficient=str(schur_gain),
        gain_over_base_deficit=str(ratio),
        strict_sixth_order_face=strict_face,
        length_six_formula_verified=length_six_formula_verified,
    )


def varied_disk_direction(
    length: int,
    direction_index: int,
) -> tuple[sp.Expr, ...]:
    """Return a deterministic family of exact complex directions."""

    if direction_index == 0:
        return deterministic_disk_direction(length)
    values = [sp.Integer(0)]
    for offset in range(1, length):
        real_sign = -1 if direction_index * offset % 2 else 1
        imaginary_sign = (
            -1 if (direction_index + offset) % 2 else 1
        )
        values.append(
            real_sign
            * sp.Rational(
                (direction_index + 2) * offset + 1,
                (17 + 2 * direction_index) * length,
            )
            + sp.I
            * imaginary_sign
            * sp.Rational(
                (2 * direction_index + 1) * offset
                + direction_index
                + 1,
                (23 + 3 * direction_index) * length,
            )
        )
    return tuple(values)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-length", type=int, default=3)
    parser.add_argument("--maximum-length", type=int, default=6)
    parser.add_argument("--direction-count", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run the exact recentered cubic-normal audit."""

    args = parse_args()
    if (
        args.minimum_length < 3
        or args.maximum_length < args.minimum_length
    ):
        raise ValueError("require 3 <= minimum length <= maximum length")
    if args.direction_count < 1:
        raise ValueError("direction count must be positive")
    records = [
        audit_direction(
            varied_disk_direction(length, direction_index),
            direction_index,
        )
        for length in range(
            args.minimum_length,
            args.maximum_length + 1,
        )
        for direction_index in range(args.direction_count)
    ]
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
