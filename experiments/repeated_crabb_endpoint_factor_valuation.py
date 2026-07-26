#!/usr/bin/env python3
"""Audit analytic endpoint-factor divisibility in Smith coordinates."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import sympy as sp


@dataclass(frozen=True)
class EndpointFactorValuationRecord:
    """One exact Smith-coordinate factorization audit."""

    record_type: str
    output_dimension: int
    factor_domain_dimension: int
    generic_rank: int
    smith_exponents: str
    surviving_generic_kernel_dimension: int
    kernel_block_zero: bool
    minimum_valuation_surplus: int
    analytic_factor_expected: bool
    analytic_factor_obtained: bool
    exact_factorization_residual_count: int
    pointwise_flag_condition_including_origin: bool
    small_parameter_factor_norm: str
    all_checks_passed: bool


def smith_diagonal(
    parameter: sp.Symbol,
    exponents: tuple[int, ...],
    rows: int,
    columns: int,
) -> sp.Matrix:
    """Return a rectangular Smith diagonal."""

    diagonal = sp.zeros(rows, columns)
    for index, exponent in enumerate(exponents):
        diagonal[index, index] = parameter**exponent
    return diagonal


def polynomial_valuation(
    expression: sp.Expr,
    parameter: sp.Symbol,
) -> int | None:
    """Return the valuation of a polynomial, or None for zero."""

    expanded = sp.expand(expression)
    if expanded == 0:
        return None
    polynomial = sp.Poly(expanded, parameter)
    return min(monomial[0] for monomial, _ in polynomial.terms())


def valuation_conditions(
    matrix: sp.Matrix,
    exponents: tuple[int, ...],
    parameter: sp.Symbol,
) -> tuple[bool, int]:
    """Return kernel-block zero and the minimum required surplus."""

    rank = len(exponents)
    dimension = matrix.rows
    kernel_zero = all(
        sp.expand(matrix[row, column]) == 0
        for row in range(rank, dimension)
        for column in range(rank, dimension)
    )
    surpluses: list[int] = []
    for row in range(rank):
        for column in range(row, rank):
            valuation = polynomial_valuation(
                matrix[row, column],
                parameter,
            )
            if valuation is not None:
                surpluses.append(
                    valuation - min(exponents[row], exponents[column])
                )
        for column in range(rank, dimension):
            valuation = polynomial_valuation(
                matrix[row, column],
                parameter,
            )
            if valuation is not None:
                surpluses.append(valuation - exponents[row])
    minimum = min(surpluses, default=0)
    return kernel_zero, minimum


def reconstruct_factor(
    matrix: sp.Matrix,
    exponents: tuple[int, ...],
    domain_dimension: int,
    parameter: sp.Symbol,
) -> sp.Matrix:
    """Return the explicit Smith-coordinate factor, possibly Laurent."""

    rank = len(exponents)
    factor = sp.zeros(matrix.rows, domain_dimension)
    for row in range(rank):
        diagonal = parameter ** exponents[row]
        factor[row, row] = sp.cancel(
            matrix[row, row] / (2 * diagonal)
        )
        for column in range(row + 1, rank):
            if exponents[row] <= exponents[column]:
                factor[column, row] = sp.cancel(
                    matrix[column, row] / diagonal
                )
            else:
                other_diagonal = parameter ** exponents[column]
                factor[row, column] = sp.cancel(
                    matrix[row, column] / other_diagonal
                )
        for column in range(rank, matrix.rows):
            factor[column, row] = sp.cancel(
                matrix[column, row] / diagonal
            )
    return factor


def is_polynomial_matrix(
    matrix: sp.Matrix,
    parameter: sp.Symbol,
) -> bool:
    """Return whether every entry is polynomial in the parameter."""

    for entry in matrix:
        _, denominator = sp.fraction(sp.cancel(entry))
        if sp.degree(denominator, parameter) > 0:
            return False
    return True


def maximum_entry_norm(
    matrix: sp.Matrix,
    parameter: sp.Symbol,
    value: sp.Rational,
) -> str:
    """Return the maximum absolute entry after one exact substitution."""

    entries = [
        abs(sp.N(entry.subs(parameter, value), 16))
        for entry in matrix
    ]
    return str(max(entries, default=sp.Integer(0)))


def pointwise_flag_condition(
    matrix: sp.Matrix,
    exponents: tuple[int, ...],
    parameter: sp.Symbol,
) -> bool:
    """Check the generic kernel block and the enlarged origin flag."""

    rank = len(exponents)
    generic_kernel_zero = all(
        sp.expand(matrix[row, column]) == 0
        for row in range(rank, matrix.rows)
        for column in range(rank, matrix.rows)
    )
    origin_zero = all(
        sp.expand(entry.subs(parameter, 0)) == 0
        for entry in matrix
    )
    return generic_kernel_zero and origin_zero


def audit_matrix(
    record_type: str,
    matrix: sp.Matrix,
    exponents: tuple[int, ...],
    domain_dimension: int,
    expected: bool,
    parameter: sp.Symbol,
) -> EndpointFactorValuationRecord:
    """Audit one good or obstructed analytic factor problem."""

    diagonal = smith_diagonal(
        parameter,
        exponents,
        matrix.rows,
        domain_dimension,
    )
    kernel_zero, minimum_surplus = valuation_conditions(
        matrix,
        exponents,
        parameter,
    )
    analytic_expected = kernel_zero and minimum_surplus >= 0
    factor = reconstruct_factor(
        matrix,
        exponents,
        domain_dimension,
        parameter,
    )
    residual = (
        matrix
        - factor * diagonal.T
        - diagonal * factor.T
    ).applyfunc(sp.expand)
    residual_count = sum(entry != 0 for entry in residual)
    factor_is_analytic = is_polynomial_matrix(factor, parameter)
    pointwise = pointwise_flag_condition(
        matrix,
        exponents,
        parameter,
    )
    verified = bool(
        analytic_expected == expected
        and factor_is_analytic == expected
        and residual_count == 0
        and pointwise
    )
    if not verified:
        raise RuntimeError(
            "the endpoint-factor valuation audit failed: "
            f"type={record_type}, "
            f"expected={expected}, "
            f"criterion={analytic_expected}, "
            f"factor={factor_is_analytic}, "
            f"residual={residual_count}, "
            f"pointwise={pointwise}"
        )
    return EndpointFactorValuationRecord(
        record_type=record_type,
        output_dimension=matrix.rows,
        factor_domain_dimension=domain_dimension,
        generic_rank=len(exponents),
        smith_exponents=",".join(map(str, exponents)),
        surviving_generic_kernel_dimension=(
            matrix.rows - len(exponents)
        ),
        kernel_block_zero=kernel_zero,
        minimum_valuation_surplus=minimum_surplus,
        analytic_factor_expected=analytic_expected,
        analytic_factor_obtained=factor_is_analytic,
        exact_factorization_residual_count=residual_count,
        pointwise_flag_condition_including_origin=pointwise,
        small_parameter_factor_norm=maximum_entry_norm(
            factor,
            parameter,
            sp.Rational(1, 1000),
        ),
        all_checks_passed=verified,
    )


def standard_records() -> list[EndpointFactorValuationRecord]:
    """Return exact positive and pointwise-only obstruction audits."""

    parameter = sp.Symbol("s", real=True)

    diagonal = smith_diagonal(parameter, (1, 3), 3, 4)
    good_factor = sp.Matrix(
        [
            [1 + parameter, 2, 0, 0],
            [parameter**2, 1 - parameter, 0, 0],
            [3, parameter, 0, 0],
        ]
    )
    good_matrix = (
        good_factor * diagonal.T
        + diagonal * good_factor.T
    )

    full_diagonal = smith_diagonal(parameter, (1, 2, 4), 3, 3)
    full_factor = sp.Matrix(
        [
            [1, parameter, 2],
            [2 - parameter, 3, parameter],
            [parameter**2, 1, 4],
        ]
    )
    full_matrix = (
        full_factor * full_diagonal.T
        + full_diagonal * full_factor.T
    )

    bad_scalar = sp.Matrix([[parameter]])
    bad_cross = sp.zeros(3)
    bad_cross[1, 2] = parameter
    bad_cross[2, 1] = parameter

    return [
        audit_matrix(
            "analytic_active_kernel",
            good_matrix,
            (1, 3),
            4,
            True,
            parameter,
        ),
        audit_matrix(
            "analytic_full_rank",
            full_matrix,
            (1, 2, 4),
            3,
            True,
            parameter,
        ),
        audit_matrix(
            "pointwise_only_scalar",
            bad_scalar,
            (2,),
            1,
            False,
            parameter,
        ),
        audit_matrix(
            "pointwise_only_cross",
            bad_cross,
            (1, 3),
            3,
            False,
            parameter,
        ),
    ]


def write_records(
    records: list[EndpointFactorValuationRecord],
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
            "experiments/"
            "repeated_crabb_endpoint_factor_valuation_s70226.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the deterministic audit."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
