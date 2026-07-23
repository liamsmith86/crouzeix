#!/usr/bin/env python3
"""Audit the Newton-edge LDL factor of the pure-defect Hessian.

At the elliptic Crabb axis, the quadratic part in the free defect
tangent is independent of the equality direction.  This checker
reconstructs that universal Hessian exactly and performs an LDL*
factorization over rational truncated series.

The observed leading lower factor is the finite-path version of

    (1 + c z^2) / (1 - c z^2)
      = 1 + 2 sum_{r>=1} c^r z^(2r).

The terminal row has the expected generalized-endpoint weight 4/3.
The finite audit identifies the spectral-factor target; it is not an
all-size proof of the factorization.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path

from crabb_palindromic_elliptic_hessian import (
    Matrix,
    Series,
    defect_quadratic_data,
    one,
    operator_expansion_from_coefficients,
    zero,
)


DEFAULT_LENGTH = 8
DEFAULT_ORDER = 6


@dataclass(frozen=True)
class DefectHessianFactorRecord:
    """Serializable exact LDL-edge audit."""

    length: int
    dimension: int
    audit_order: int
    interior_diagonal_constant: str
    terminal_diagonal_constant: str
    lower_terms_vanish: bool
    interior_edge_coefficient: str
    terminal_edge_coefficient: str
    parity_cross_terms_vanish: bool
    checked_interior_edges: int
    checked_terminal_edges: int


def ldl_factor(matrix: Matrix) -> tuple[Matrix, list[Series]]:
    """Return a unit-lower LDL* factorization over the series ring."""

    dimension = len(matrix)
    order = matrix[0][0].order
    lower = [
        [zero(order) for _ in range(dimension)]
        for _ in range(dimension)
    ]
    diagonal = [zero(order) for _ in range(dimension)]
    for row in range(dimension):
        lower[row][row] = one(order)
        for column in range(row):
            residual = matrix[row][column] - sum(
                lower[row][index]
                * diagonal[index]
                * lower[column][index]
                for index in range(column)
            )
            lower[row][column] = residual / diagonal[column]
        diagonal[row] = matrix[row][row] - sum(
            lower[row][index] ** 2 * diagonal[index]
            for index in range(row)
        )
    return lower, diagonal


def make_record(
    length: int = DEFAULT_LENGTH,
    order: int = DEFAULT_ORDER,
) -> DefectHessianFactorRecord:
    """Construct and verify one exact universal-Hessian factor."""

    if length < 4 or order <= length // 2:
        raise ValueError(
            "require length >= 4 and order > floor(length / 2)"
        )
    operator, coefficients = operator_expansion_from_coefficients(
        length + 1,
        [0] * (length - 1),
        order,
    )
    _, _, hessian, _ = defect_quadratic_data(
        operator,
        coefficients,
        order,
    )
    lower, diagonal = ldl_factor(hessian)

    interior_edges = 0
    terminal_edges = 0
    lower_terms_vanish = True
    parity_cross_terms_vanish = True
    for row in range(length):
        for column in range(row):
            difference = row - column
            entry = lower[row][column]
            if difference % 2:
                parity_cross_terms_vanish &= (
                    entry.valuation() == order
                )
                continue
            grade = difference // 2
            lower_terms_vanish &= all(
                entry.coefficient(degree) == 0
                for degree in range(grade)
            )
            expected = (
                Fraction(4, 3)
                if row == length - 1
                else Fraction(2)
            )
            if entry.coefficient(grade) != expected:
                raise AssertionError(
                    "the LDL factor violated its Newton edge"
                )
            if row == length - 1:
                terminal_edges += 1
            else:
                interior_edges += 1

    if not lower_terms_vanish or not parity_cross_terms_vanish:
        raise AssertionError(
            "the LDL factor acquired a lower or wrong-parity term"
        )
    if any(
        entry.coefficient(0) != 4
        for entry in diagonal[:-1]
    ) or diagonal[-1].coefficient(0) != Fraction(8, 3):
        raise AssertionError("the disk Hessian diagonal changed")

    return DefectHessianFactorRecord(
        length=length,
        dimension=length + 1,
        audit_order=order,
        interior_diagonal_constant="4",
        terminal_diagonal_constant="8/3",
        lower_terms_vanish=lower_terms_vanish,
        interior_edge_coefficient="2",
        terminal_edge_coefficient="4/3",
        parity_cross_terms_vanish=parity_cross_terms_vanish,
        checked_interior_edges=interior_edges,
        checked_terminal_edges=terminal_edges,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length", type=int, default=DEFAULT_LENGTH)
    parser.add_argument("--order", type=int, default=DEFAULT_ORDER)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the exact factor audit and optionally persist it."""

    args = parse_args()
    record = make_record(args.length, args.order)
    line = json.dumps(asdict(record), sort_keys=True)
    print(line, flush=True)
    if args.output is not None:
        args.output.write_text(f"{line}\n", encoding="utf-8")


if __name__ == "__main__":
    main()
