#!/usr/bin/env python3
"""Audit the Schur-parameter expansion of the transfer Markov gap."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import sympy as sp


Matrix = sp.Matrix


@dataclass(frozen=True)
class SchurMarkovLaplacianRecord:
    """One exact all-grade Markov-Laplacian audit."""

    length: int
    multiplicity: int
    active_schur_grades: tuple[int, ...]
    tangent_degrees: tuple[int, ...]
    hermitian_basis_size: int
    left_parseval_residual_count: int
    right_parseval_residual_count: int
    superoperator_residual_count: int
    dirichlet_residual_count: int
    skew_terminal_cancellation_residual_count: int
    positive_test_energy: str
    all_checks_passed: bool


def adjoint(matrix: Matrix) -> Matrix:
    """Return the conjugate transpose."""

    return matrix.conjugate().T


def zero_matrix(size: int) -> Matrix:
    """Return a square exact zero matrix."""

    return sp.zeros(size, size)


def matrix_sum(matrices: list[Matrix], size: int) -> Matrix:
    """Sum a finite list of equal-size matrices exactly."""

    total = zero_matrix(size)
    for matrix in matrices:
        total += matrix
    return total


def residual_count(matrix: Matrix) -> int:
    """Count entries that remain nonzero after exact simplification."""

    return sum(
        int(sp.simplify(entry) != 0)
        for entry in matrix
    )


def exact_matrix(
    size: int,
    grade: int,
    offset: int,
) -> Matrix:
    """Build a deterministic noncommuting Gaussian-rational matrix."""

    return Matrix(
        size,
        size,
        lambda row, column: (
            sp.Rational(
                (grade + 2) * (row + 1) - (column + offset),
                grade + offset + 2,
            )
            + sp.I
            * sp.Rational(
                (row + offset + 1) * (column + 1) - grade,
                grade + offset + 3,
            )
        ),
    )


def hermitian_basis(size: int) -> list[Matrix]:
    """Return the standard real Hermitian basis."""

    basis: list[Matrix] = []
    for row in range(size):
        diagonal = zero_matrix(size)
        diagonal[row, row] = 1
        basis.append(diagonal)
        for column in range(row + 1, size):
            real = zero_matrix(size)
            real[row, column] = 1
            real[column, row] = 1
            basis.append(real)

            imaginary = zero_matrix(size)
            imaginary[row, column] = sp.I
            imaginary[column, row] = -sp.I
            basis.append(imaginary)
    return basis


def commutator(left: Matrix, right: Matrix) -> Matrix:
    """Return ``left * right - right * left``."""

    return left * right - right * left


def frobenius_square(matrix: Matrix) -> sp.Expr:
    """Return the exact squared Frobenius norm."""

    return sp.simplify(sp.trace(adjoint(matrix) * matrix))


def expected_laplacian(
    test: Matrix,
    directions: list[Matrix],
) -> Matrix:
    """Return the double-commutator Laplacian from L293."""

    terms: list[Matrix] = []
    for direction in directions:
        direction_adjoint = adjoint(direction)
        terms.append(
            commutator(
                direction_adjoint,
                commutator(direction, test),
            )
            + commutator(
                direction,
                commutator(direction_adjoint, test),
            )
        )
    return matrix_sum(terms, test.rows)


def audit_case(
    length: int,
    multiplicity: int,
    active_grades: tuple[int, ...],
    offset: int,
) -> SchurMarkovLaplacianRecord:
    """Audit one exact formal Schur tangent."""

    directions = [
        exact_matrix(multiplicity, grade, offset)
        for grade in active_grades
    ]
    tangent_coefficients: list[Matrix] = []
    tangent_degrees: list[int] = []
    for grade, direction in zip(active_grades, directions, strict=True):
        tangent_degrees.extend((grade, 2 * length - grade))
        tangent_coefficients.extend((direction, -adjoint(direction)))

    radial = matrix_sum(
        [
            direction * adjoint(direction)
            + adjoint(direction) * direction
            for direction in directions
        ],
        multiplicity,
    )
    skew_seed = exact_matrix(multiplicity, length, offset + 4)
    skew_terminal = (skew_seed - adjoint(skew_seed)) / 2
    terminal_two = -radial / 2 + skew_terminal

    left_parseval = (
        terminal_two
        + adjoint(terminal_two)
        + matrix_sum(
            [
                coefficient * adjoint(coefficient)
                for coefficient in tangent_coefficients
            ],
            multiplicity,
        )
    )
    right_parseval = (
        terminal_two
        + adjoint(terminal_two)
        + matrix_sum(
            [
                adjoint(coefficient) * coefficient
                for coefficient in tangent_coefficients
            ],
            multiplicity,
        )
    )

    def channel_two(test: Matrix, terminal: Matrix) -> Matrix:
        return (
            terminal * test
            + test * adjoint(terminal)
            + matrix_sum(
                [
                    coefficient * test * adjoint(coefficient)
                    for coefficient in tangent_coefficients
                ],
                multiplicity,
            )
        )

    basis = hermitian_basis(multiplicity)
    superoperator_residuals = 0
    skew_residuals = 0
    dirichlet_residuals = 0
    for test in basis:
        channel = channel_two(test, terminal_two)
        channel_adjoint = channel_two(
            test,
            adjoint(terminal_two),
        )
        laplacian = -channel - channel_adjoint
        expected = expected_laplacian(test, directions)
        superoperator_residuals += residual_count(laplacian - expected)

        no_skew_channel = channel_two(test, -radial / 2)
        no_skew_adjoint = channel_two(test, -radial / 2)
        skew_residuals += residual_count(
            laplacian + no_skew_channel + no_skew_adjoint
        )

        energy = sp.simplify(sp.trace(test * laplacian))
        expected_energy = sp.simplify(
            2
            * sum(
                (
                    frobenius_square(
                        commutator(test, direction)
                    )
                    for direction in directions
                ),
                sp.S.Zero,
            )
        )
        dirichlet_residuals += int(
            sp.simplify(energy - expected_energy) != 0
        )

    positive_test = exact_matrix(multiplicity, length + 1, offset + 7)
    positive_test = positive_test + adjoint(positive_test)
    positive_energy = sp.simplify(
        2
        * sum(
            (
                frobenius_square(
                    commutator(positive_test, direction)
                )
                for direction in directions
            ),
            sp.S.Zero,
        )
    )

    degrees_are_valid = bool(
        len(tangent_degrees) == len(set(tangent_degrees))
        and length not in tangent_degrees
    )
    checks_passed = bool(
        residual_count(left_parseval) == 0
        and residual_count(right_parseval) == 0
        and superoperator_residuals == 0
        and dirichlet_residuals == 0
        and skew_residuals == 0
        and degrees_are_valid
        and positive_energy > 0
    )
    if not checks_passed:
        raise RuntimeError(
            "Schur Markov-Laplacian audit failed: "
            f"L={length}, m={multiplicity}, "
            f"left={residual_count(left_parseval)}, "
            f"right={residual_count(right_parseval)}, "
            f"super={superoperator_residuals}, "
            f"energy={dirichlet_residuals}, "
            f"skew={skew_residuals}, "
            f"degrees={degrees_are_valid}"
        )

    return SchurMarkovLaplacianRecord(
        length=length,
        multiplicity=multiplicity,
        active_schur_grades=active_grades,
        tangent_degrees=tuple(sorted(tangent_degrees)),
        hermitian_basis_size=len(basis),
        left_parseval_residual_count=residual_count(left_parseval),
        right_parseval_residual_count=residual_count(right_parseval),
        superoperator_residual_count=superoperator_residuals,
        dirichlet_residual_count=dirichlet_residuals,
        skew_terminal_cancellation_residual_count=skew_residuals,
        positive_test_energy=str(positive_energy),
        all_checks_passed=checks_passed,
    )


def standard_records() -> list[SchurMarkovLaplacianRecord]:
    """Return deterministic exact audits."""

    return [
        audit_case(3, 2, (1, 2), 1),
        audit_case(5, 2, (1, 3, 4), 3),
        audit_case(4, 3, (1, 2, 3), 5),
    ]


def write_records(
    records: list[SchurMarkovLaplacianRecord],
    output: Path,
) -> str:
    """Write deterministic JSON Lines and return its SHA-256."""

    output.parent.mkdir(parents=True, exist_ok=True)
    payload = "".join(
        json.dumps(asdict(record), sort_keys=True) + "\n"
        for record in records
    )
    output.write_text(payload, encoding="utf-8")
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the exact audit."""

    args = parse_args()
    records = standard_records()
    digest = None
    if args.output is not None:
        digest = write_records(records, args.output)
    print(
        json.dumps(
            {
                "all_checks_passed": all(
                    record.all_checks_passed for record in records
                ),
                "dataset_sha256": digest,
                "records": [asdict(record) for record in records],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
