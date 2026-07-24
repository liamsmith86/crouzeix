#!/usr/bin/env python3
"""Audit inverse-block-Toeplitz Hardy equality at repeated Crabb blocks.

Replace each scalar entry in the ``L``-level disk Gramian by an
``m x m`` copy-space block.  If ``B = H^{-1}``, block row elimination
gives the exact reflection identity

    P_I (A - S) = -(B_+ - B_-) (B_+ + B_-)^{-1} P_2.

Thus the finite Hardy residual vanishes locally exactly when ``B`` is
Hermitian block Toeplitz.  The normalized rank-``m`` Stein defect then
gives generalized endpoint eigenspaces at one and four.

This script checks the reflection identity off equality, the exact
diagonal-difference rank at the Crabb point, and explicit equality
certificates with genuinely noncommuting block-Toeplitz coefficients.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.linalg import eigh, solve_discrete_lyapunov


@dataclass(frozen=True)
class BlockDifferentialRecord:
    """One exact-rank block residual differential audit."""

    record_kind: str
    dimension: int
    length: int
    multiplicity: int
    hermitian_variable_count: int
    residual_real_rank: int
    expected_rank: int
    kernel_dimension: int
    expected_block_toeplitz_dimension: int
    all_identities_verified: bool


@dataclass(frozen=True)
class BlockReflectionRecord:
    """One nonlinear inverse-Gram reflection-factor audit."""

    record_kind: str
    dimension: int
    length: int
    multiplicity: int
    diagonal_difference_norm: str
    reflection_identity_error: str
    minimum_hermitian_eigenvalue: str
    all_checks_passed: bool


@dataclass(frozen=True)
class BlockEqualityRecord:
    """One explicit inverse-block-Toeplitz equality audit."""

    record_kind: str
    dimension: int
    length: int
    multiplicity: int
    residual_norm: str
    interior_shift_error: str
    inverse_toeplitz_error: str
    distance_from_crabb_gram: str
    minimum_hermitian_eigenvalue: str
    coefficient_commutator_norm: str
    krylov_tail_minimum_singular_value: str
    minimum_generalized_eigenvalue: str
    maximum_generalized_eigenvalue: str
    condition_square: str
    lower_endpoint_multiplicity: int
    upper_endpoint_multiplicity: int
    upper_sandwich_maximum_eigenvalue: str
    endpoint_relation_error: str
    all_checks_passed: bool


AuditRecord = (
    BlockDifferentialRecord | BlockReflectionRecord | BlockEqualityRecord
)


def format_float(value: float) -> str:
    """Format a deterministic binary64 diagnostic."""

    return f"{value:.12e}"


def block_shift(length: int, multiplicity: int) -> np.ndarray:
    """Return the level shift tensored with the copy identity."""

    shift = np.zeros(
        ((length + 1) * multiplicity,) * 2,
        dtype=complex,
    )
    identity = np.eye(multiplicity)
    for level in range(length):
        shift[
            level * multiplicity : (level + 1) * multiplicity,
            (level + 1) * multiplicity : (level + 2) * multiplicity,
        ] = identity
    return shift


def extend_hermitian(
    hermitian: np.ndarray,
    length: int,
    multiplicity: int,
) -> np.ndarray:
    """Append one zero level to a block Hermitian Gramian."""

    dimension = (length + 1) * multiplicity
    extended = np.zeros((dimension, dimension), dtype=complex)
    extended[: length * multiplicity, : length * multiplicity] = hermitian
    return extended


def block_disk_model(
    hermitian: np.ndarray,
    length: int,
    multiplicity: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Construct the coefficient-gauge block disk model."""

    extended = extend_hermitian(hermitian, length, multiplicity)
    shift = block_shift(length, multiplicity)
    coordinate = extended + shift.conj().T @ extended @ shift
    operator = np.linalg.solve(coordinate, 2 * extended @ shift)
    return extended, coordinate, operator, shift


def block_hardy_residual(
    hermitian: np.ndarray,
    length: int,
    multiplicity: int,
) -> np.ndarray:
    """Return the finite terminal-orbit Hardy residual."""

    _, _, operator, shift = block_disk_model(
        hermitian,
        length,
        multiplicity,
    )
    dimension = len(operator)
    endpoint = np.zeros((dimension, multiplicity), dtype=complex)
    endpoint[length * multiplicity :, :] = np.eye(multiplicity)
    residual = np.zeros(
        ((length - 1) * multiplicity,) * 2,
        dtype=complex,
    )
    power = np.eye(dimension, dtype=complex)
    for degree in range(1, length):
        previous = power
        power = power @ operator
        coefficient = power @ endpoint - shift @ previous @ endpoint
        for column in range(1, length):
            block = coefficient[
                column * multiplicity : (column + 1) * multiplicity,
                :,
            ]
            residual[
                (degree - 1)
                * multiplicity : degree
                * multiplicity,
                (column - 1)
                * multiplicity : column
                * multiplicity,
            ] = block.T
    return residual


def linearized_block_residual(
    coefficient: np.ndarray,
    length: int,
    multiplicity: int,
) -> np.ndarray:
    """Return the block diagonal-difference differential."""

    residual = np.zeros(
        ((length - 1) * multiplicity,) * 2,
        dtype=complex,
    )
    for row in range(length - 1):
        reflected = length - 2 - row
        for column in range(length - 1):
            upper = coefficient[
                (column + 1)
                * multiplicity : (column + 2)
                * multiplicity,
                (reflected + 1)
                * multiplicity : (reflected + 2)
                * multiplicity,
            ]
            lower = coefficient[
                column * multiplicity : (column + 1) * multiplicity,
                reflected
                * multiplicity : (reflected + 1)
                * multiplicity,
            ]
            residual[
                row * multiplicity : (row + 1) * multiplicity,
                column * multiplicity : (column + 1) * multiplicity,
            ] = (upper - lower).T
    return residual


def pack_hermitian(matrix: np.ndarray) -> np.ndarray:
    """Pack a complex Hermitian matrix into real coordinates."""

    size = len(matrix)
    values = [float(np.real(matrix[index, index])) for index in range(size)]
    for row in range(size):
        for column in range(row + 1, size):
            value = matrix[row, column]
            values.extend((float(np.real(value)), float(np.imag(value))))
    return np.asarray(values)


def pack_reversed_residual(
    residual: np.ndarray,
    length: int,
    multiplicity: int,
) -> np.ndarray:
    """Pack the Hermitian row reversal of a linearized residual."""

    reversal = np.fliplr(np.eye(length - 1))
    reversed_residual = (
        np.kron(reversal, np.eye(multiplicity)) @ residual
    )
    hermitian = (
        reversed_residual + reversed_residual.conj().T
    ) / 2
    return pack_hermitian(hermitian)


def hermitian_basis(size: int) -> Sequence[np.ndarray]:
    """Return the standard real basis of Hermitian matrices."""

    basis = []
    for index in range(size):
        matrix = np.zeros((size, size), dtype=complex)
        matrix[index, index] = 1
        basis.append(matrix)
    for row in range(size):
        for column in range(row + 1, size):
            real = np.zeros((size, size), dtype=complex)
            real[row, column] = 1
            real[column, row] = 1
            basis.append(real)
            imaginary = np.zeros((size, size), dtype=complex)
            imaginary[row, column] = 1j
            imaginary[column, row] = -1j
            basis.append(imaginary)
    return tuple(basis)


def inverse_gram_diagonal_data(
    hermitian: np.ndarray,
    multiplicity: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``B``, its adjacent principal difference, and their sum."""

    inverse = np.linalg.inv(hermitian)
    leading = inverse[:-multiplicity, :-multiplicity]
    trailing = inverse[multiplicity:, multiplicity:]
    return inverse, trailing - leading, trailing + leading


def tail_selector(length: int, multiplicity: int) -> np.ndarray:
    """Return the matrix extracting levels two through ``length``."""

    selector = np.zeros(
        ((length - 1) * multiplicity, (length + 1) * multiplicity),
        dtype=complex,
    )
    selector[:, 2 * multiplicity :] = np.eye(
        (length - 1) * multiplicity
    )
    return selector


def terminal_krylov_tail(
    operator: np.ndarray,
    length: int,
    multiplicity: int,
) -> np.ndarray:
    """Return the tail restriction of the terminal Krylov matrix."""

    endpoint = np.zeros(
        ((length + 1) * multiplicity, multiplicity),
        dtype=complex,
    )
    endpoint[length * multiplicity :, :] = np.eye(multiplicity)
    columns = []
    power = np.eye(len(operator), dtype=complex)
    for _ in range(length - 1):
        columns.append(power @ endpoint)
        power = power @ operator
    krylov = np.hstack(columns)
    return tail_selector(length, multiplicity) @ krylov


def differential_record(
    length: int,
    multiplicity: int,
) -> BlockDifferentialRecord:
    """Compute the integer rank of the block differential."""

    size = length * multiplicity
    columns = []
    for basis_element in hermitian_basis(size):
        residual = linearized_block_residual(
            basis_element,
            length,
            multiplicity,
        )
        columns.append(
            pack_reversed_residual(
                residual,
                length,
                multiplicity,
            )
        )
    differential = np.stack(columns, axis=1)
    rank = int(np.linalg.matrix_rank(differential, tol=1e-10))
    expected_rank = ((length - 1) * multiplicity) ** 2
    variable_count = size**2
    kernel_dimension = variable_count - rank
    expected_kernel = (2 * length - 1) * multiplicity**2
    verified = bool(
        rank == expected_rank
        and kernel_dimension == expected_kernel
    )
    if not verified:
        raise RuntimeError(
            "the block Hardy differential rank audit failed"
        )
    return BlockDifferentialRecord(
        record_kind="block_hardy_differential",
        dimension=(length + 1) * multiplicity,
        length=length,
        multiplicity=multiplicity,
        hermitian_variable_count=variable_count,
        residual_real_rank=rank,
        expected_rank=expected_rank,
        kernel_dimension=kernel_dimension,
        expected_block_toeplitz_dimension=expected_kernel,
        all_identities_verified=verified,
    )


def generic_hermitian(
    length: int,
    multiplicity: int,
) -> np.ndarray:
    """Return a deterministic positive non-Toeplitz Hermitian input."""

    size = length * multiplicity
    perturbation = np.zeros((size, size), dtype=complex)
    for row in range(size):
        perturbation[row, row] = (1 + 2 * row) / (20 * size)
        for column in range(row + 1, size):
            value = (
                (1 + row + 3 * column) / (80 * size)
                + 1j * (2 + 2 * row - column) / (100 * size)
            )
            perturbation[row, column] = value
            perturbation[column, row] = np.conjugate(value)
    return np.eye(size, dtype=complex) / 2 + 0.01 * perturbation


def reflection_record(
    length: int,
    multiplicity: int,
) -> BlockReflectionRecord:
    """Audit the inverse-Gram reflection identity off equality."""

    hermitian = generic_hermitian(length, multiplicity)
    _, difference, total = inverse_gram_diagonal_data(
        hermitian,
        multiplicity,
    )
    _, _, operator, shift = block_disk_model(
        hermitian,
        length,
        multiplicity,
    )
    interior = (operator - shift)[
        multiplicity : length * multiplicity,
        :,
    ]
    reflected = -difference @ np.linalg.solve(
        total,
        tail_selector(length, multiplicity),
    )
    identity_error = float(np.linalg.norm(interior - reflected))
    difference_norm = float(np.linalg.norm(difference))
    minimum_hermitian = float(np.linalg.eigvalsh(hermitian)[0])
    verified = bool(
        identity_error < 2e-12
        and difference_norm > 1e-6
        and minimum_hermitian > 0
    )
    if not verified:
        raise RuntimeError(
            "the inverse-Gram reflection identity audit failed"
        )
    return BlockReflectionRecord(
        record_kind="block_inverse_gram_reflection",
        dimension=(length + 1) * multiplicity,
        length=length,
        multiplicity=multiplicity,
        diagonal_difference_norm=format_float(difference_norm),
        reflection_identity_error=format_float(identity_error),
        minimum_hermitian_eigenvalue=format_float(minimum_hermitian),
        all_checks_passed=verified,
    )


def inverse_block_toeplitz(
    length: int,
    multiplicity: int,
) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
    """Return a positive inverse Gram with noncommuting Toeplitz legs."""

    inverse = np.eye(length * multiplicity, dtype=complex) * 2
    coefficients = []
    for offset in range(1, length):
        block = np.zeros((multiplicity, multiplicity), dtype=complex)
        for row in range(multiplicity):
            for column in range(multiplicity):
                block[row, column] = (
                    0.012
                    * (1 + offset + 2 * row + 3 * column)
                    / (length * multiplicity)
                    + 1j
                    * 0.008
                    * (2 + 2 * offset + row - column)
                    / (length * multiplicity)
                )
        coefficients.append(block)
        for level in range(length - offset):
            inverse[
                level
                * multiplicity : (level + 1)
                * multiplicity,
                (level + offset)
                * multiplicity : (level + offset + 1)
                * multiplicity,
            ] = block
            inverse[
                (level + offset)
                * multiplicity : (level + offset + 1)
                * multiplicity,
                level
                * multiplicity : (level + 1)
                * multiplicity,
            ] = block.conj().T
    return inverse, tuple(coefficients)


def equality_record(
    length: int,
    multiplicity: int,
) -> BlockEqualityRecord:
    """Construct and audit a noncommuting block equality point."""

    inverse, coefficients = inverse_block_toeplitz(
        length,
        multiplicity,
    )
    hermitian = np.linalg.inv(inverse)
    _, inverse_difference, _ = inverse_gram_diagonal_data(
        hermitian,
        multiplicity,
    )
    residual_norm = float(
        np.linalg.norm(
            block_hardy_residual(
                hermitian,
                length,
                multiplicity,
            )
        )
    )
    inverse_toeplitz_error = float(np.linalg.norm(inverse_difference))
    distance_from_crabb = float(
        np.linalg.norm(
            hermitian
            - np.eye(length * multiplicity, dtype=complex) / 2
        )
    )
    minimum_hermitian = float(np.linalg.eigvalsh(hermitian)[0])
    if len(coefficients) >= 2:
        commutator = (
            coefficients[0] @ coefficients[1]
            - coefficients[1] @ coefficients[0]
        )
    else:
        commutator = (
            coefficients[0] @ coefficients[0].conj().T
            - coefficients[0].conj().T @ coefficients[0]
        )
    commutator_norm = float(np.linalg.norm(commutator))

    extended, coordinate, operator, shift = block_disk_model(
        hermitian,
        length,
        multiplicity,
    )
    interior_shift_error = float(
        np.linalg.norm(
            (operator - shift)[
                multiplicity : length * multiplicity,
                :,
            ]
        )
    )
    krylov_minimum = float(
        np.linalg.svd(
            terminal_krylov_tail(operator, length, multiplicity),
            compute_uv=False,
        )[-1]
    )
    endpoint_zero = np.zeros(
        ((length + 1) * multiplicity, multiplicity),
        dtype=complex,
    )
    endpoint_zero[:multiplicity, :] = np.eye(multiplicity)
    endpoint_top = np.zeros_like(endpoint_zero)
    endpoint_top[length * multiplicity :, :] = np.eye(multiplicity)
    corner = endpoint_zero.conj().T @ extended @ endpoint_zero
    defect = (
        extended
        @ endpoint_zero
        @ np.linalg.inv(np.linalg.cholesky(corner)).conj().T
    )
    metric = solve_discrete_lyapunov(
        operator.conj().T,
        defect @ defect.conj().T,
    )
    generalized = eigh(metric, coordinate, eigvals_only=True)
    condition_square = float(generalized[-1] / generalized[0])
    lower_multiplicity = int(
        np.count_nonzero(np.abs(generalized - 1) < 2e-8)
    )
    upper_multiplicity = int(
        np.count_nonzero(np.abs(generalized - 4) < 2e-8)
    )
    sandwich = metric - 4 * coordinate
    sandwich_maximum = float(np.linalg.eigvalsh(sandwich)[-1])
    endpoint_error = float(
        np.linalg.norm(metric @ endpoint_top - 4 * coordinate @ endpoint_top)
    )
    verified = bool(
        residual_norm < 1e-10
        and interior_shift_error < 1e-10
        and inverse_toeplitz_error < 1e-10
        and distance_from_crabb > 1e-8
        and minimum_hermitian > 0
        and (multiplicity == 1 or commutator_norm > 1e-8)
        and krylov_minimum > 0.5
        and abs(generalized[0] - 1) < 2e-8
        and abs(generalized[-1] - 4) < 2e-8
        and abs(condition_square - 4) < 3e-8
        and lower_multiplicity == multiplicity
        and upper_multiplicity == multiplicity
        and sandwich_maximum < 2e-8
        and endpoint_error < 1e-9
    )
    if not verified:
        raise RuntimeError(
            "the inverse-block-Toeplitz equality audit failed for "
            f"length {length}, multiplicity {multiplicity}"
        )
    return BlockEqualityRecord(
        record_kind="inverse_block_toeplitz_equality",
        dimension=(length + 1) * multiplicity,
        length=length,
        multiplicity=multiplicity,
        residual_norm=format_float(residual_norm),
        interior_shift_error=format_float(interior_shift_error),
        inverse_toeplitz_error=format_float(inverse_toeplitz_error),
        distance_from_crabb_gram=format_float(distance_from_crabb),
        minimum_hermitian_eigenvalue=format_float(minimum_hermitian),
        coefficient_commutator_norm=format_float(commutator_norm),
        krylov_tail_minimum_singular_value=format_float(krylov_minimum),
        minimum_generalized_eigenvalue=format_float(generalized[0]),
        maximum_generalized_eigenvalue=format_float(generalized[-1]),
        condition_square=format_float(condition_square),
        lower_endpoint_multiplicity=lower_multiplicity,
        upper_endpoint_multiplicity=upper_multiplicity,
        upper_sandwich_maximum_eigenvalue=format_float(
            sandwich_maximum
        ),
        endpoint_relation_error=format_float(endpoint_error),
        all_checks_passed=verified,
    )


def write_records(path: Path, records: Sequence[AuditRecord]) -> None:
    """Write deterministic JSON Lines output atomically."""

    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(asdict(record), sort_keys=True))
            handle.write("\n")
    temporary.replace(path)


def main() -> None:
    """Run repeated block-Hardy equality audits."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--minimum-length", type=int, default=3)
    parser.add_argument("--maximum-length", type=int, default=4)
    parser.add_argument(
        "--multiplicities",
        type=int,
        nargs="+",
        default=(2, 3),
        help="copy multiplicities audited at every requested length",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/crabb_block_hardy_equality_s70224.jsonl"
        ),
    )
    args = parser.parse_args()

    records: list[AuditRecord] = []
    for multiplicity in args.multiplicities:
        for length in range(args.minimum_length, args.maximum_length + 1):
            differential = differential_record(length, multiplicity)
            records.append(differential)
            print(
                "verified block differential "
                f"length {length}, multiplicity {multiplicity}",
                flush=True,
            )
            reflection = reflection_record(length, multiplicity)
            records.append(reflection)
            print(
                "verified inverse-Gram reflection "
                f"length {length}, multiplicity {multiplicity}",
                flush=True,
            )
            equality = equality_record(length, multiplicity)
            records.append(equality)
            print(
                "verified block equality "
                f"length {length}, multiplicity {multiplicity}",
                flush=True,
            )
    write_records(args.output, records)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
