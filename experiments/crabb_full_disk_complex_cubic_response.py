#!/usr/bin/env python3
"""Verify the full symbolic complex cubic response in the first active size.

The ordinary finite-direction audit in
``crabb_full_disk_cubic_normal_exact.py`` checks several exact complex
rays.  This slower checker instead propagates five generic complex
Toeplitz coefficients through the complete characteristic/Riemann
recurrence.  It proves that mode three is the only cubic true-normal
mode for ``L=6`` and verifies its closed Pluecker formula coefficientwise.
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
    disk_model_series,
)
from crabb_full_disk_correction_isometry import plucker_correction
from crabb_full_disk_cubic_normal_exact import highest_mode_cubic_response


@dataclass(frozen=True)
class ComplexCubicResponseRecord:
    """Exact metadata for the generic complex length-six response."""

    dimension: int
    length: int
    symbolic_complex_variable_count: int
    normal_polarization_count: int
    active_cubic_modes: tuple[int, ...]
    quadratic_response_zero: bool
    inactive_cubic_responses_zero: bool
    mode_three_tensor_formula_verified: bool
    mode_three_closed_formula_verified: bool
    full_complex_response_verified: bool


def complex_expression_zero(expression: sp.Expr) -> bool:
    """Test an identity after expanding generic variables into real parts."""

    return sp.expand(sp.expand_complex(expression)) == 0


def closed_mode_three_response(
    direction: tuple[sp.Expr, ...],
) -> sp.Expr:
    """Return ``-22 C_6/45`` on the full complex length-six slice."""

    coefficients = direction[1:]
    reversal = tuple(sp.conjugate(coefficients[4 - index]) for index in range(5))

    def pluecker(left: int, right: int) -> sp.Expr:
        return (
            coefficients[left] * reversal[right] - coefficients[right] * reversal[left]
        )

    cubic = (
        6 * coefficients[3] * pluecker(0, 3)
        - 5 * coefficients[4] * pluecker(0, 2)
        + 2 * coefficients[3] * pluecker(1, 2)
    )
    return sp.expand(-sp.Rational(22, 45) * cubic)


def build_record() -> ComplexCubicResponseRecord:
    """Run the generic complex characteristic/Riemann audit."""

    variables = sp.symbols("z0:5")
    direction = (sp.Integer(0), *variables)
    correction = plucker_correction(direction)
    operator, metric = disk_model_series(
        direction,
        order=3,
        correction=correction,
    )
    characteristic = characteristic_series(operator)
    normal_count = 0
    active_modes = []
    quadratic_zero = True
    inactive_cubic_zero = True
    mode_three_responses = []
    for mode in range(3, len(direction) + 2):
        mode_responses = []
        for normal in crabb_normal_directions([metric[0]], mode):
            response = directional_gradient_series(
                operator,
                metric,
                characteristic,
                normal,
            )
            quadratic_zero &= complex_expression_zero(response[2])
            mode_responses.append(response[3])
            normal_count += 1
        mode_active = any(
            not complex_expression_zero(response) for response in mode_responses
        )
        if mode_active:
            active_modes.append(mode)
        if mode == 3:
            mode_three_responses = mode_responses
        else:
            inactive_cubic_zero &= not mode_active

    if len(mode_three_responses) != 2:
        raise RuntimeError("mode three did not have two polarizations")
    actual = mode_three_responses[0] + sp.I * mode_three_responses[1]
    tensor_candidate = highest_mode_cubic_response(direction)
    closed_candidate = closed_mode_three_response(direction)
    tensor_formula = complex_expression_zero(actual - tensor_candidate)
    closed_formula = complex_expression_zero(tensor_candidate - closed_candidate)
    verified = bool(
        quadratic_zero
        and inactive_cubic_zero
        and active_modes == [3]
        and tensor_formula
        and closed_formula
    )
    if not verified:
        raise RuntimeError("the full complex cubic-response audit failed")
    return ComplexCubicResponseRecord(
        dimension=7,
        length=6,
        symbolic_complex_variable_count=len(variables),
        normal_polarization_count=normal_count,
        active_cubic_modes=tuple(active_modes),
        quadratic_response_zero=quadratic_zero,
        inactive_cubic_responses_zero=inactive_cubic_zero,
        mode_three_tensor_formula_verified=tensor_formula,
        mode_three_closed_formula_verified=closed_formula,
        full_complex_response_verified=verified,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Regenerate and serialize the exact response metadata."""

    args = parse_args()
    record = build_record()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    line = json.dumps(asdict(record), sort_keys=True)
    args.output.write_text(f"{line}\n", encoding="utf-8")
    print(line, flush=True)


if __name__ == "__main__":
    main()
