#!/usr/bin/env python3
"""Audit promotion of the first residual kernel in inverse-Gram coordinates.

The equality manifold is ``B = H^{-1}`` block Toeplitz.  After
subtracting block-diagonal means, the adjacent-principal difference
is an isomorphism on the transverse complement.  A common copy-space
kernel of the first residual therefore forces the leading transverse
inverse-Gram coefficient to live entirely on the orthogonal copy
complement.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Sequence

import numpy as np

from crabb_block_hardy_equality import (
    block_disk_model,
    format_float,
    linearized_block_residual,
)


@dataclass(frozen=True)
class KernelPromotionRecord:
    """One inverse-Gram residual-kernel promotion audit."""

    record_kind: str
    length: int
    multiplicity: int
    promoted_kernel_dimension: int
    diagonal_difference_reconstruction_error: str
    residual_factorization_error: str
    residual_gram_kernel_error: str
    transverse_kernel_support_error: str
    first_active_gram_eigenvalue: str
    all_checks_passed: bool


def block_toeplitz_projection(
    matrix: np.ndarray,
    length: int,
    multiplicity: int,
) -> np.ndarray:
    """Project a Hermitian block matrix onto constant block diagonals."""

    projection = np.zeros_like(matrix)
    for offset in range(length):
        blocks = [
            matrix[
                level
                * multiplicity : (level + 1)
                * multiplicity,
                (level + offset)
                * multiplicity : (level + offset + 1)
                * multiplicity,
            ]
            for level in range(length - offset)
        ]
        mean = sum(blocks, np.zeros_like(blocks[0])) / len(blocks)
        for level in range(length - offset):
            projection[
                level
                * multiplicity : (level + 1)
                * multiplicity,
                (level + offset)
                * multiplicity : (level + offset + 1)
                * multiplicity,
            ] = mean
            if offset:
                projection[
                    (level + offset)
                    * multiplicity : (level + offset + 1)
                    * multiplicity,
                    level
                    * multiplicity : (level + 1)
                    * multiplicity,
                ] = mean.conj().T
    return projection


def adjacent_principal_difference(
    matrix: np.ndarray,
    multiplicity: int,
) -> np.ndarray:
    """Return the trailing-minus-leading principal difference."""

    return (
        matrix[multiplicity:, multiplicity:]
        - matrix[:-multiplicity, :-multiplicity]
    )


def integrate_zero_mean_difference(
    difference: np.ndarray,
    length: int,
    multiplicity: int,
) -> np.ndarray:
    """Invert the adjacent difference on the zero-Toeplitz complement."""

    matrix = np.zeros(
        (length * multiplicity,) * 2,
        dtype=complex,
    )
    for offset in range(length):
        count = length - offset
        deltas = [
            difference[
                level
                * multiplicity : (level + 1)
                * multiplicity,
                (level + offset)
                * multiplicity : (level + offset + 1)
                * multiplicity,
            ]
            for level in range(count - 1)
        ]
        first = -sum(
            (
                (count - 1 - index) * delta
                for index, delta in enumerate(deltas)
            ),
            np.zeros((multiplicity, multiplicity), dtype=complex),
        ) / count
        blocks = [first]
        for delta in deltas:
            blocks.append(blocks[-1] + delta)
        for level, block in enumerate(blocks):
            matrix[
                level
                * multiplicity : (level + 1)
                * multiplicity,
                (level + offset)
                * multiplicity : (level + offset + 1)
                * multiplicity,
            ] = block
            if offset:
                matrix[
                    (level + offset)
                    * multiplicity : (level + offset + 1)
                    * multiplicity,
                    level
                    * multiplicity : (level + 1)
                    * multiplicity,
                ] = block.conj().T
    return (matrix + matrix.conj().T) / 2


def deterministic_copy_unitary(multiplicity: int) -> np.ndarray:
    """Return a deterministic dense copy-space unitary."""

    seed = np.zeros((multiplicity, multiplicity), dtype=complex)
    for row in range(multiplicity):
        for column in range(multiplicity):
            seed[row, column] = (
                1
                + 2 * row
                + 3 * column
                + 1j * (2 - row + 2 * column)
            )
            if row == column:
                seed[row, column] += 2 * multiplicity
    unitary, triangular = np.linalg.qr(seed)
    phases = np.diag(triangular).copy()
    phases /= np.abs(phases)
    return unitary @ np.diag(np.conjugate(phases))


def deterministic_transverse(
    length: int,
    multiplicity: int,
    kernel_dimension: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return a zero-Toeplitz transverse coefficient supported off a kernel."""

    unitary = deterministic_copy_unitary(multiplicity)
    kernel = unitary[:, :kernel_dimension]
    complement = unitary[:, kernel_dimension:]
    active_dimension = multiplicity - kernel_dimension
    active = np.zeros(
        (length * active_dimension,) * 2,
        dtype=complex,
    )
    for row in range(len(active)):
        active[row, row] = (1 + row) / (2 * len(active))
        for column in range(row + 1, len(active)):
            value = (
                (1 + row + 2 * column) / (3 * len(active))
                + 1j
                * (2 - row + column)
                / (4 * len(active))
            )
            active[row, column] = value
            active[column, row] = np.conjugate(value)
    active -= block_toeplitz_projection(
        active,
        length,
        active_dimension,
    )
    embedding = np.kron(np.eye(length), complement)
    transverse = embedding @ active @ embedding.conj().T
    return transverse, kernel


def terminal_krylov_tail(
    length: int,
    multiplicity: int,
) -> np.ndarray:
    """Return the base terminal Krylov tail."""

    base_hermitian = np.eye(length * multiplicity, dtype=complex) / 2
    _, _, operator, _ = block_disk_model(
        base_hermitian,
        length,
        multiplicity,
    )
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
    return krylov[2 * multiplicity :, :]


def arrange_actual_residual(
    actual: np.ndarray,
    length: int,
    multiplicity: int,
) -> np.ndarray:
    """Arrange output/degree blocks in the checker's transposed convention."""

    arranged = np.zeros_like(actual)
    for interior in range(length - 1):
        for degree in range(length - 1):
            block = actual[
                interior
                * multiplicity : (interior + 1)
                * multiplicity,
                degree
                * multiplicity : (degree + 1)
                * multiplicity,
            ]
            arranged[
                degree
                * multiplicity : (degree + 1)
                * multiplicity,
                interior
                * multiplicity : (interior + 1)
                * multiplicity,
            ] = block.T
    return arranged


def residual_copy_gram(
    actual: np.ndarray,
    length: int,
    multiplicity: int,
) -> np.ndarray:
    """Return the terminal copy Gram of all actual residual blocks."""

    gram = np.zeros((multiplicity, multiplicity), dtype=complex)
    for interior in range(length - 1):
        for degree in range(length - 1):
            block = actual[
                interior
                * multiplicity : (interior + 1)
                * multiplicity,
                degree
                * multiplicity : (degree + 1)
                * multiplicity,
            ]
            gram += block.conj().T @ block
    return (gram + gram.conj().T) / 2


def promotion_record(
    length: int,
    multiplicity: int,
    kernel_dimension: int,
) -> KernelPromotionRecord:
    """Audit leading residual-kernel promotion."""

    transverse, kernel = deterministic_transverse(
        length,
        multiplicity,
        kernel_dimension,
    )
    difference = adjacent_principal_difference(
        transverse,
        multiplicity,
    )
    recovered = integrate_zero_mean_difference(
        difference,
        length,
        multiplicity,
    )
    reconstruction_error = float(np.linalg.norm(recovered - transverse))

    krylov = terminal_krylov_tail(length, multiplicity)
    actual = -difference @ krylov / 4
    gram = residual_copy_gram(actual, length, multiplicity)
    gram_kernel_error = float(np.linalg.norm(gram @ kernel))
    level_kernel = np.kron(np.eye(length), kernel)
    support_error = float(
        np.linalg.norm(transverse @ level_kernel)
        + np.linalg.norm(level_kernel.conj().T @ transverse)
    )

    gram_eigenvalues = np.linalg.eigvalsh(gram)
    first_active = float(gram_eigenvalues[kernel_dimension])
    hermitian_direction = -transverse / 4
    arranged = linearized_block_residual(
        hermitian_direction,
        length,
        multiplicity,
    )
    factorization_error = float(
        np.linalg.norm(arranged - arrange_actual_residual(
            actual,
            length,
            multiplicity,
        ))
    )
    verified = bool(
        reconstruction_error < 2e-12
        and factorization_error < 2e-12
        and gram_kernel_error < 2e-12
        and support_error < 2e-12
        and np.count_nonzero(np.abs(gram_eigenvalues) < 2e-10)
        == kernel_dimension
        and first_active > 1e-8
    )
    if not verified:
        raise RuntimeError("the inverse-Gram kernel promotion audit failed")
    return KernelPromotionRecord(
        record_kind="inverse_gram_kernel_promotion",
        length=length,
        multiplicity=multiplicity,
        promoted_kernel_dimension=kernel_dimension,
        diagonal_difference_reconstruction_error=format_float(
            reconstruction_error
        ),
        residual_factorization_error=format_float(factorization_error),
        residual_gram_kernel_error=format_float(gram_kernel_error),
        transverse_kernel_support_error=format_float(support_error),
        first_active_gram_eigenvalue=format_float(first_active),
        all_checks_passed=verified,
    )


def write_records(
    path: Path,
    records: Sequence[KernelPromotionRecord],
) -> None:
    """Write deterministic JSON Lines output atomically."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(asdict(record), sort_keys=True))
            handle.write("\n")
    temporary.replace(path)


def main() -> None:
    """Run inverse-Gram kernel-promotion audits."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--minimum-length", type=int, default=2)
    parser.add_argument("--maximum-length", type=int, default=6)
    parser.add_argument(
        "--multiplicities",
        type=int,
        nargs="+",
        default=(2, 3, 4),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_inverse_gram_kernel_s70224.jsonl"
        ),
    )
    args = parser.parse_args()

    records = []
    for multiplicity in args.multiplicities:
        for kernel_dimension in range(1, multiplicity):
            for length in range(
                args.minimum_length,
                args.maximum_length + 1,
            ):
                record = promotion_record(
                    length,
                    multiplicity,
                    kernel_dimension,
                )
                records.append(record)
                print(
                    "verified inverse-Gram kernel promotion "
                    f"length {length}, multiplicity {multiplicity}, "
                    f"kernel {kernel_dimension}",
                    flush=True,
                )
    write_records(args.output, records)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
