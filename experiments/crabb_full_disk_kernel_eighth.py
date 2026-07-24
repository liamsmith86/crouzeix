#!/usr/bin/env python3
"""Regenerate the sixth-kernel formula and its eighth-order pair lift.

The cubic Hardy residual depends only on the low Wronskian
coefficients

    S_t = sum_(i+j=t) (j-i) z_i conjugate(z_(n-1-j)).

This checker verifies its closed all-index formula with unrestricted
complex symbols.  It then restricts to every reversal-pair kernel
plane, verifies the single-entry fourth Hardy residual and its exact
eighth-order endpoint deficit, and independently checks that every
true circular-normal response vanishes through quartic order.
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
from crabb_full_disk_base_jet_exact import endpoint_delta_series
from crabb_full_disk_correction_isometry import plucker_correction
from crabb_full_disk_sixth_hardy_factor import (
    hardy_residual_matrices,
    matrix_is_zero,
)


@dataclass(frozen=True)
class SymbolicKernelRecord:
    """One unrestricted-symbol cubic and pair-fourth audit."""

    record_kind: str
    dimension: int
    length: int
    variable_count: int
    reversal_pair_count: int
    cubic_closed_formula_verified: bool
    pair_lower_orders_zero: bool
    pair_fourth_formula_verified: bool
    all_identities_verified: bool


@dataclass(frozen=True)
class ExactPairRecord:
    """One exact reversal-pair eighth-order audit."""

    record_kind: str
    dimension: int
    length: int
    left_index: int
    right_index: int
    separation: int
    base_eighth_deficit: str
    predicted_eighth_deficit: str
    hardy_fourth_norm_square: str
    endpoint_formula_verified: bool
    hardy_factor_verified: bool
    normal_response_count: int
    responses_through_quartic_zero: bool
    all_identities_verified: bool


def wronskian_fluxes(
    direction: tuple[sp.Expr, ...],
) -> dict[int, sp.Expr]:
    """Return the low coefficients of ``p q' - p' q``."""

    length = len(direction)
    coefficient_count = length - 1
    coefficients = direction[1:]
    return {
        total: sp.expand(
            sum(
                (
                    (right - left)
                    * coefficients[left]
                    * sp.conjugate(
                        coefficients[
                            coefficient_count - 1 - right
                        ]
                    )
                    for left in range(coefficient_count)
                    for right in range(coefficient_count)
                    if left + right == total
                ),
                sp.Integer(0),
            )
        )
        for total in range(1, coefficient_count - 1)
    }


def closed_cubic_residual(
    direction: tuple[sp.Expr, ...],
) -> sp.Matrix:
    """Return the all-index cubic Hardy residual from its fluxes."""

    length = len(direction)
    coefficient_count = length - 1
    coefficients = direction[1:]
    fluxes = wronskian_fluxes(direction)
    residual = sp.zeros(coefficient_count)
    for row in range(coefficient_count):
        for column in range(row + 1, coefficient_count):
            value = sp.Integer(0)
            positive_flux = coefficient_count - 2 - column
            if positive_flux in fluxes:
                value += (
                    4
                    * coefficients[coefficient_count - 1 - row]
                    * fluxes[positive_flux]
                    / (coefficient_count - column)
                )
            negative_flux = coefficient_count - 2 - row
            if negative_flux in fluxes:
                value -= (
                    4
                    * coefficients[coefficient_count - 1 - column]
                    * fluxes[negative_flux]
                    / (coefficient_count - row)
                )
            reflected_positive = row - 1
            if reflected_positive in fluxes:
                value += (
                    4
                    * sp.conjugate(coefficients[column])
                    * sp.conjugate(fluxes[reflected_positive])
                    / (row + 1)
                )
            reflected_negative = column - 1
            if reflected_negative in fluxes:
                value -= (
                    4
                    * sp.conjugate(coefficients[row])
                    * sp.conjugate(fluxes[reflected_negative])
                    / (column + 1)
                )
            residual[row, column] = sp.expand(value)
            residual[column, row] = -residual[row, column]
    return residual


def pair_direction(
    length: int,
    left: int,
    left_value: sp.Expr,
    right_value: sp.Expr,
) -> tuple[sp.Expr, ...]:
    """Return a direction supported on one reversal pair."""

    coefficient_count = length - 1
    right = coefficient_count - 1 - left
    if not 0 <= left < right < coefficient_count:
        raise ValueError("indices do not form a strict reversal pair")
    coefficients = [sp.Integer(0)] * coefficient_count
    coefficients[left] = left_value
    coefficients[right] = right_value
    return (sp.Integer(0), *coefficients)


def predicted_pair_fourth_residual(
    direction: tuple[sp.Expr, ...],
    left: int,
) -> sp.Matrix:
    """Return the single central anti-diagonal fourth residual."""

    length = len(direction)
    coefficient_count = length - 1
    right = coefficient_count - 1 - left
    separation = right - left
    left_value = direction[left + 1]
    right_value = direction[right + 1]
    imbalance = (
        left_value * sp.conjugate(left_value)
        - right_value * sp.conjugate(right_value)
    )
    value = sp.expand(
        sp.Rational(8 * separation, length) * imbalance**2
    )
    residual = sp.zeros(coefficient_count)
    residual[left, right] = value
    residual[right, left] = -value
    return residual


def symbolic_audit(length: int) -> SymbolicKernelRecord:
    """Verify the cubic formula and every pair-fourth formula."""

    variables = sp.symbols(f"z0:{length - 1}")
    direction = (sp.Integer(0), *variables)
    cubic = hardy_residual_matrices(direction)[0][2]
    cubic_formula = matrix_is_zero(
        cubic - closed_cubic_residual(direction)
    )

    left_value, right_value = sp.symbols("pair_left pair_right")
    pair_count = (length - 1) // 2
    lower_zero = True
    fourth_formula = True
    for left in range(pair_count):
        pair = pair_direction(
            length,
            left,
            left_value,
            right_value,
        )
        residuals, _ = hardy_residual_matrices(pair, order=4)
        lower_zero &= all(
            matrix_is_zero(residuals[degree])
            for degree in range(3)
        )
        fourth_formula &= matrix_is_zero(
            residuals[3]
            - predicted_pair_fourth_residual(pair, left)
        )
    verified = bool(cubic_formula and lower_zero and fourth_formula)
    if not verified:
        raise RuntimeError(
            f"symbolic kernel/eighth audit failed in length {length}"
        )
    return SymbolicKernelRecord(
        record_kind="symbolic_kernel",
        dimension=length + 1,
        length=length,
        variable_count=2 * (length - 1),
        reversal_pair_count=pair_count,
        cubic_closed_formula_verified=cubic_formula,
        pair_lower_orders_zero=lower_zero,
        pair_fourth_formula_verified=fourth_formula,
        all_identities_verified=verified,
    )


def exact_pair_audit(length: int, left: int) -> ExactPairRecord:
    """Audit one exact pair endpoint and every true-normal response."""

    right = length - 2 - left
    left_value = (
        sp.Rational(left + 2, 7)
        + sp.I * sp.Rational(left + 1, 11)
    )
    right_value = (
        sp.Rational(left + 3, 13)
        - sp.I * sp.Rational(left + 2, 17)
    )
    direction = pair_direction(
        length,
        left,
        left_value,
        right_value,
    )
    residuals, _ = hardy_residual_matrices(direction, order=4)
    fourth = residuals[3]
    norm_square = sp.factor(
        sum(
            (
                entry * sp.conjugate(entry)
                for entry in fourth
            ),
            sp.Integer(0),
        )
    )
    base_deficit = sp.factor(
        -2 * endpoint_delta_series(direction, order=8)[8]
    )
    separation = right - left
    imbalance = (
        left_value * sp.conjugate(left_value)
        - right_value * sp.conjugate(right_value)
    )
    predicted = sp.factor(
        sp.Rational(512 * separation**2, length**2)
        * imbalance**4
    )
    endpoint_formula = sp.simplify(base_deficit - predicted) == 0
    hardy_factor = sp.simplify(base_deficit - 4 * norm_square) == 0

    correction = plucker_correction(direction)
    operator, metric = disk_model_series(
        direction,
        order=4,
        correction=correction,
    )
    characteristic = characteristic_series(operator)
    response_count = 0
    responses_zero = True
    for mode in range(3, length + 2):
        for normal in crabb_normal_directions([metric[0]], mode):
            response = directional_gradient_series(
                operator,
                metric,
                characteristic,
                normal,
            )
            responses_zero &= all(
                sp.simplify(response[degree]) == 0
                for degree in range(2, 5)
            )
            response_count += 1
    verified = bool(endpoint_formula and hardy_factor and responses_zero)
    if not verified:
        raise RuntimeError(
            f"exact pair eighth audit failed in length {length}, "
            f"pair {left},{right}"
        )
    return ExactPairRecord(
        record_kind="exact_pair",
        dimension=length + 1,
        length=length,
        left_index=left,
        right_index=right,
        separation=separation,
        base_eighth_deficit=str(base_deficit),
        predicted_eighth_deficit=str(predicted),
        hardy_fourth_norm_square=str(norm_square),
        endpoint_formula_verified=endpoint_formula,
        hardy_factor_verified=hardy_factor,
        normal_response_count=response_count,
        responses_through_quartic_zero=responses_zero,
        all_identities_verified=verified,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--symbolic-maximum-length", type=int, default=8)
    parser.add_argument("--exact-maximum-length", type=int, default=10)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run and serialize the symbolic and exact audits."""

    args = parse_args()
    if (
        args.symbolic_maximum_length < 4
        or args.exact_maximum_length < 4
    ):
        raise ValueError("maximum lengths are below four")

    records: list[SymbolicKernelRecord | ExactPairRecord] = []
    for length in range(4, args.symbolic_maximum_length + 1):
        records.append(symbolic_audit(length))
        print(f"verified symbolic kernel length {length}", flush=True)
    for length in range(4, args.exact_maximum_length + 1):
        for left in range((length - 1) // 2):
            records.append(exact_pair_audit(length, left))
            print(
                f"verified exact pair length {length}, left {left}",
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
