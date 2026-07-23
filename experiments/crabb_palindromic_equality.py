#!/usr/bin/env python3
"""Regenerate the phase-palindromic Crabb disk equality theorem.

The script works in the diagonal gauge

    z_k = conjugate(z_(L-k)).

It avoids symbolic matrix square roots and inverses.  The coefficient
operator is constructed in its explicit companion form, and every
identity is checked after multiplication by the coordinate Gramian.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from crabb_disk_toeplitz_quartic import extend, shift


@dataclass(frozen=True)
class PalindromicEqualityRecord:
    dimension: int
    direction_kind: str
    companion_identity: bool
    stein_identity: bool
    endpoint_gram_identity: bool
    characteristic_identity: bool
    observability_identity: bool
    blaschke_identity: bool


def phase_palindromic_coefficients(
    length: int,
    symbolic: bool,
) -> list[sp.Expr]:
    """Return coefficients satisfying ``z_k=conjugate(z_(L-k))``."""

    count = length - 1
    coefficients: list[sp.Expr | None] = [None] * count
    for index in range(count):
        if coefficients[index] is not None:
            continue
        reverse_index = count - 1 - index
        if index == reverse_index:
            coefficient = (
                sp.symbols(f"ar{index + 1}", real=True)
                if symbolic
                else sp.Rational(index + 1, index + 4)
            )
            coefficients[index] = coefficient
            continue

        if symbolic:
            real, imaginary = sp.symbols(
                f"ar{index + 1} ai{index + 1}",
                real=True,
            )
            coefficient = real + sp.I * imaginary
        else:
            coefficient = (
                sp.Rational(index + 1, index + 5)
                + sp.I * sp.Rational(index + 2, index + 7)
            )
        coefficients[index] = coefficient
        coefficients[reverse_index] = sp.conjugate(coefficient)

    return [
        coefficient
        for coefficient in coefficients
        if coefficient is not None
    ]


def toeplitz_chart(
    length: int,
    coefficients: list[sp.Expr],
) -> sp.Matrix:
    """Return the extended Toeplitz matrix with diagonal one half."""

    matrix = sp.eye(length) / 2
    for offset, coefficient in enumerate(coefficients, start=1):
        for row in range(length - offset):
            matrix[row, row + offset] = coefficient
            matrix[row + offset, row] = sp.conjugate(coefficient)
    return extend(matrix)


def companion_operator(
    dimension: int,
    coefficients: list[sp.Expr],
) -> sp.Matrix:
    """Return the explicit coefficient operator from proof equation (8)."""

    operator = sp.zeros(dimension)
    operator[0, 1] = 2
    last = dimension - 1
    for column in range(2, dimension):
        coefficient = coefficients[column - 2]
        operator[column - 1, column] = 1
        operator[0, column] = 2 * coefficient
        operator[last, column] = -2 * coefficient
    return operator


def reversal(dimension: int) -> sp.Matrix:
    """Return coordinate reversal."""

    matrix = sp.zeros(dimension)
    for index in range(dimension):
        matrix[index, dimension - 1 - index] = 1
    return matrix


def evaluate_polynomial(
    expression: sp.Expr,
    variable: sp.Symbol,
    matrix: sp.Matrix,
) -> sp.Matrix:
    """Evaluate a scalar polynomial at a square matrix by Horner's rule."""

    result = sp.zeros(matrix.rows)
    for coefficient in sp.Poly(expression, variable).all_coeffs():
        result = result * matrix + coefficient * sp.eye(matrix.rows)
    return result


def matrix_is_zero(matrix: sp.Matrix) -> bool:
    """Test exact zero after expansion."""

    return all(
        sp.expand(entry) == 0
        for entry in matrix
    )


def audit_dimension(
    dimension: int,
    symbolic: bool,
) -> PalindromicEqualityRecord:
    """Audit every algebraic identity in one dimension."""

    length = dimension - 1
    coefficients = phase_palindromic_coefficients(length, symbolic)
    chart = toeplitz_chart(length, coefficients)
    nilpotent_shift = shift(dimension)
    coordinate_gramian = (
        chart + nilpotent_shift.T * chart * nilpotent_shift
    )
    operator = companion_operator(dimension, coefficients)
    companion_identity = matrix_is_zero(
        coordinate_gramian * operator
        - 2 * chart * nilpotent_shift
    )

    first = sp.eye(dimension)[:, 0]
    last = sp.eye(dimension)[:, -1]
    defect = chart * first
    reflected_defect = reversal(dimension) * sp.conjugate(defect)
    metric = (
        coordinate_gramian
        - defect * defect.conjugate().T
        + 2 * reflected_defect * reflected_defect.conjugate().T
    )
    stein_identity = matrix_is_zero(
        metric
        - operator.conjugate().T * metric * operator
        - defect * defect.conjugate().T
    )

    endpoint_gram_identity = (
        matrix_is_zero(coordinate_gramian * first - defect)
        and matrix_is_zero(
            coordinate_gramian * last - reflected_defect
        )
        and sp.expand(defect.conjugate().T * first)[0]
        == sp.Rational(1, 2)
        and sp.expand(reflected_defect.conjugate().T * last)[0]
        == sp.Rational(1, 2)
        and sp.expand(defect.conjugate().T * last)[0] == 0
    )

    variable = sp.symbols("xi")
    polynomial = (
        variable**length
        + 2
        * sum(
            coefficient * variable**offset
            for offset, coefficient in enumerate(
                coefficients,
                start=1,
            )
        )
    )
    characteristic_identity = (
        sp.expand(
            operator.charpoly(variable).as_expr()
            - variable * polynomial
        )
        == 0
    )

    observability = sp.Matrix.vstack(
        *[
            defect.conjugate().T * operator**power
            for power in range(dimension)
        ]
    )
    unextended_chart = chart[:length, :length]
    observability_identity = (
        sp.expand(
            observability.det()
            - 2 ** (length - 1) * unextended_chart.det()
        )
        == 0
    )

    reversed_polynomial = (
        1
        + 2
        * sum(
            coefficient * variable**offset
            for offset, coefficient in enumerate(
                coefficients,
                start=1,
            )
        )
    )
    polynomial_at_operator = evaluate_polynomial(
        polynomial,
        variable,
        operator,
    )
    reversed_at_operator = evaluate_polynomial(
        reversed_polynomial,
        variable,
        operator,
    )
    blaschke_identity = matrix_is_zero(
        polynomial_at_operator
        - (
            4
            * first
            * reflected_defect.conjugate().T
            * reversed_at_operator
        )
    )

    checks = (
        companion_identity,
        stein_identity,
        endpoint_gram_identity,
        characteristic_identity,
        observability_identity,
        blaschke_identity,
    )
    if not all(checks):
        raise RuntimeError(
            f"palindromic equality audit failed in size {dimension}: "
            f"companion={companion_identity}, stein={stein_identity}, "
            f"gram={endpoint_gram_identity}, "
            f"characteristic={characteristic_identity}, "
            f"observability={observability_identity}, "
            f"blaschke={blaschke_identity}"
        )

    return PalindromicEqualityRecord(
        dimension=dimension,
        direction_kind="symbolic" if symbolic else "exact_rational",
        companion_identity=companion_identity,
        stein_identity=stein_identity,
        endpoint_gram_identity=endpoint_gram_identity,
        characteristic_identity=characteristic_identity,
        observability_identity=observability_identity,
        blaschke_identity=blaschke_identity,
    )


def nongauged_coefficients(length: int) -> list[sp.Expr]:
    """Return small exact coefficients with phase ``omega=-1``."""

    count = length - 1
    coefficients: list[sp.Expr | None] = [None] * count
    denominator = 80 * max(count, 1)
    for index in range(count):
        if coefficients[index] is not None:
            continue
        reverse_index = count - 1 - index
        if index == reverse_index:
            coefficients[index] = (
                sp.I * sp.Rational(index + 1, denominator)
            )
            continue
        coefficient = (
            sp.Rational(index + 1, denominator)
            + sp.I * sp.Rational(index + 2, denominator + count)
        )
        coefficients[index] = coefficient
        coefficients[reverse_index] = -sp.conjugate(coefficient)
    return [
        coefficient
        for coefficient in coefficients
        if coefficient is not None
    ]


def audit_nongauged_phase(maximum_size: int = 5) -> dict[str, object]:
    """Check the invariant formulas without first removing the phase."""

    checked_dimensions: list[int] = []
    variable = sp.symbols("xi")
    for dimension in range(3, maximum_size + 1):
        length = dimension - 1
        coefficients = nongauged_coefficients(length)
        chart = toeplitz_chart(length, coefficients)
        unextended_chart = chart[:length, :length]
        if not all(
            sp.simplify(
                unextended_chart[:size, :size].det()
            ).is_positive
            for size in range(1, length + 1)
        ):
            raise RuntimeError("the nongauged audit chart is not positive")

        nilpotent_shift = shift(dimension)
        coordinate_gramian = (
            chart + nilpotent_shift.T * chart * nilpotent_shift
        )
        operator = (
            2
            * coordinate_gramian.inv()
            * chart
            * nilpotent_shift
        )
        first = sp.eye(dimension)[:, 0]
        defect = chart * first
        reflected_defect = reversal(dimension) * sp.conjugate(defect)
        metric = (
            coordinate_gramian
            - defect * defect.conjugate().T
            + 2
            * reflected_defect
            * reflected_defect.conjugate().T
        )
        if not matrix_is_zero(
            metric
            - operator.conjugate().T * metric * operator
            - defect * defect.conjugate().T
        ):
            raise RuntimeError("nongauged Stein identity failed")

        characteristic = operator.charpoly(variable).as_expr()
        polynomial = sp.cancel(characteristic / variable)
        coefficients_descending = sp.Poly(
            polynomial,
            variable,
        ).all_coeffs()
        reversed_polynomial = sum(
            sp.conjugate(coefficient) * variable**index
            for index, coefficient in enumerate(
                coefficients_descending
            )
        )
        polynomial_at_operator = evaluate_polynomial(
            polynomial,
            variable,
            operator,
        )
        reversed_at_operator = evaluate_polynomial(
            reversed_polynomial,
            variable,
            operator,
        )
        if not matrix_is_zero(
            polynomial_at_operator
            - (
                4
                * first
                * reflected_defect.conjugate().T
                * reversed_at_operator
            )
        ):
            raise RuntimeError("nongauged Blaschke identity failed")
        checked_dimensions.append(dimension)

    return {
        "phase": "-1",
        "nongauged_checked_dimensions": checked_dimensions,
        "stein_identity": True,
        "blaschke_identity": True,
    }


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=3)
    parser.add_argument("--maximum-size", type=int, default=10)
    parser.add_argument(
        "--symbolic-maximum-size",
        type=int,
        default=6,
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run exact audits and optionally persist their JSONL records."""

    args = parse_args()
    if args.minimum_size < 3 or args.maximum_size < args.minimum_size:
        raise ValueError("size range must satisfy 3 <= minimum <= maximum")
    records = [
        audit_dimension(
            dimension,
            symbolic=dimension <= args.symbolic_maximum_size,
        )
        for dimension in range(args.minimum_size, args.maximum_size + 1)
    ]
    lines = [
        *[
            json.dumps(asdict(record), sort_keys=True)
            for record in records
        ],
        json.dumps(audit_nongauged_phase(), sort_keys=True),
    ]
    for line in lines:
        print(line, flush=True)
    if args.output is not None:
        args.output.write_text(
            "".join(f"{line}\n" for line in lines),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
