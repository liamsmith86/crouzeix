#!/usr/bin/env python3
"""Prove exact rigidity of the doubled terminal edge at the first delay.

Replace the balanced ellipse pencil

    S + c (I + F) S* (I + E)

by ``S + c (I + lambda*F) S* (I + E)`` while keeping the boundary
metric fixed.  On the length-three monomial channel, the transfer has
``B_1=0``.  The degree-two final-defect Schur residual is computed
symbolically and vanishes only at the physical value ``lambda=1``.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import sympy as sp

from crabb_palindromic_elliptic_hessian import direct_map_coefficients


Matrix = sp.Matrix


@dataclass(frozen=True)
class TerminalMultiplierRigidityRecord:
    """The exact first-delay terminal-multiplier calculation."""

    state_dimension: int
    delayed_transfer_grade: int
    first_diagonal_entry: str
    second_diagonal_entry: str
    physical_multiplier: str
    physical_residual_is_zero: bool
    multiplier_is_unique: bool
    all_checks_passed: bool


def boundary_metric_second_coefficient(
    partial: Matrix,
    initial: Matrix,
    final: Matrix,
) -> Matrix:
    """Return ``[c^2] P_bl`` from L219."""

    return partial * final * partial.T - partial.T * initial * partial


def exact_record() -> TerminalMultiplierRigidityRecord:
    """Return the exact symbolic length-three obstruction."""

    multiplier = sp.symbols("lambda", real=True)
    identity = sp.eye(3)
    partial = sp.Matrix(
        [
            [0, 1, 0],
            [0, 0, 1],
            [0, 0, 0],
        ]
    )
    initial = identity - partial.T * partial
    final = identity - partial * partial.T

    coefficients = direct_map_coefficients(2, 3)
    scalar_linear = coefficients[0]
    scalar_cubic = coefficients[1]
    scalar_quintic = coefficients[2]
    scalar_checks = (
        scalar_linear.coefficient(0) == 1
        and scalar_linear.coefficient(1) == 0
        and scalar_linear.coefficient(2) == 2
        and scalar_cubic.coefficient(0) == 0
        and scalar_cubic.coefficient(1) == -1
        and scalar_cubic.coefficient(2) == 0
        and scalar_quintic.coefficient(0) == 0
        and scalar_quintic.coefficient(1) == 0
        and scalar_quintic.coefficient(2) == 1
    )

    reverse_edge = (identity + multiplier * final) * partial.T * (identity + initial)
    operator_zero = partial
    operator_one = reverse_edge - partial**3
    operator_two = (
        2 * partial
        - reverse_edge * partial**2
        - partial * reverse_edge * partial
        - partial**2 * reverse_edge
        + partial**5
    )
    metric_two = boundary_metric_second_coefficient(
        partial,
        initial,
        final,
    )

    dual_zero = identity - operator_zero * operator_zero.T
    dual_one = -(operator_one * operator_zero.T + operator_zero * operator_one.T)
    dual_two = (
        -metric_two
        - operator_two * operator_zero.T
        - operator_zero * operator_two.T
        - operator_one * operator_one.T
        + operator_zero * metric_two * operator_zero.T
    )

    final_index = next(index for index in range(3) if final[index, index] == 1)
    complement_indices = [index for index in range(3) if index != final_index]

    def pivot(matrix: Matrix) -> Matrix:
        return matrix.extract([final_index], [final_index])

    def cross(matrix: Matrix) -> Matrix:
        return matrix.extract(complement_indices, [final_index])

    def corner(matrix: Matrix) -> Matrix:
        return matrix.extract(
            complement_indices,
            complement_indices,
        )

    pivot_zero = pivot(dual_zero)
    pivot_one = pivot(dual_one)
    pivot_inverse_zero = pivot_zero.inv()
    pivot_inverse_one = -pivot_inverse_zero * pivot_one * pivot_inverse_zero
    cross_zero = cross(dual_zero)
    cross_one = cross(dual_one)
    cross_two = cross(dual_two)
    square_two = (
        cross_zero * pivot_inverse_zero * cross_two.T
        + cross_zero * pivot_inverse_one * cross_one.T
        + cross_one * pivot_inverse_zero * cross_one.T
        + cross_one * pivot_inverse_one * cross_zero.T
        + cross_two * pivot_inverse_zero * cross_zero.T
    )
    residual_two = sp.simplify(corner(dual_two) - square_two)
    expected = sp.diag(
        1 - multiplier**2,
        2 * multiplier - 2,
    )

    physical_zero = residual_two.subs(multiplier, 1) == sp.zeros(2)
    unique = sp.solve(
        list(residual_two),
        [multiplier],
        dict=True,
    ) == [{multiplier: 1}]
    verified = scalar_checks and residual_two == expected and physical_zero and unique
    if not verified:
        raise RuntimeError("the exact terminal-multiplier rigidity calculation failed")
    return TerminalMultiplierRigidityRecord(
        state_dimension=3,
        delayed_transfer_grade=2,
        first_diagonal_entry=str(1 - multiplier**2),
        second_diagonal_entry=str(2 * multiplier - 2),
        physical_multiplier="1",
        physical_residual_is_zero=physical_zero,
        multiplier_is_unique=unique,
        all_checks_passed=verified,
    )


def write_records(
    records: list[TerminalMultiplierRigidityRecord],
    output: Path,
) -> str:
    """Write deterministic JSON Lines atomically and return its hash."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)
    return hashlib.sha256(output.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/repeated_crabb_terminal_multiplier_rigidity_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact rigidity record."""

    args = parse_args()
    records = [exact_record()]
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
