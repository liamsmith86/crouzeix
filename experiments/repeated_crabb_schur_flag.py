#!/usr/bin/env python3
"""Audit the finite Schur-orthogonal flag of repeated upper endpoints."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Sequence

import numpy as np
import sympy as sp

from crabb_block_hardy_equality import format_float
from repeated_crabb_first_residual_endpoint import canonical_endpoint
from repeated_crabb_inverse_gram_kernel import (
    block_toeplitz_projection,
    deterministic_copy_unitary,
    deterministic_transverse,
)


@dataclass(frozen=True)
class SymbolicFlagRecord:
    """Exact two-stage analytic Gram/Schur audit."""

    record_kind: str
    first_reduced_valuation: int
    first_reduced_leading_rank: int
    second_reduced_valuation: int
    second_reduced_leading_coefficient: str
    all_checks_passed: bool


@dataclass(frozen=True)
class CanonicalFlagRecord:
    """One actual repeated-Hardy endpoint flag audit."""

    record_kind: str
    length: int
    multiplicity: int
    first_kernel_dimension: int
    step: str
    full_endpoint_minimum_eigenvalue: str
    active_scaled_minimum_eigenvalue: str
    reduced_scaled_minimum_eigenvalue: str
    reduced_scaled_norm: str
    all_checks_passed: bool


AuditRecord = SymbolicFlagRecord | CanonicalFlagRecord


def symbolic_flag_record() -> SymbolicFlagRecord:
    """Verify two exact Schur promotions in a polynomial Gram family."""

    parameter = sp.symbols("parameter", real=True)
    residual = sp.Matrix(
        [
            [parameter, 2 * parameter**2, 3 * parameter**3],
            [parameter**2, parameter**2, 2 * parameter**3],
            [0, parameter**3, parameter**3],
            [0, 0, parameter**4],
        ]
    )
    gram = residual.T * residual
    active = gram[:1, :1]
    cross = gram[:1, 1:]
    kernel = gram[1:, 1:]
    first_reduced = sp.simplify(
        kernel - cross.T * active.inv() * cross
    )
    first_leading = first_reduced.applyfunc(
        lambda entry: sp.limit(entry / parameter**4, parameter, 0)
    )

    second_active = first_reduced[:1, :1]
    second_cross = first_reduced[:1, 1:]
    second_kernel = first_reduced[1:, 1:]
    second_reduced = sp.factor(
        (
            second_kernel
            - second_cross.T * second_active.inv() * second_cross
        )[0]
    )
    second_leading = sp.limit(
        second_reduced / parameter**6,
        parameter,
        0,
    )
    verified = bool(
        first_leading == sp.diag(1, 0)
        and first_leading.rank() == 1
        and second_leading == 1
    )
    if not verified:
        raise RuntimeError("the symbolic Schur flag audit failed")
    return SymbolicFlagRecord(
        record_kind="symbolic_analytic_schur_flag",
        first_reduced_valuation=4,
        first_reduced_leading_rank=1,
        second_reduced_valuation=6,
        second_reduced_leading_coefficient=str(second_leading),
        all_checks_passed=verified,
    )


def deterministic_full_transverse(
    length: int,
    multiplicity: int,
) -> np.ndarray:
    """Return one dense zero-Toeplitz Hermitian coefficient."""

    size = length * multiplicity
    matrix = np.zeros((size, size), dtype=complex)
    for row in range(size):
        matrix[row, row] = (1 + row) / (5 * size)
        for column in range(row + 1, size):
            value = (
                (1 + row + 2 * column) / (4 * size)
                + 1j * (2 - row + column) / (6 * size)
            )
            matrix[row, column] = value
            matrix[column, row] = np.conjugate(value)
    return matrix - block_toeplitz_projection(
        matrix,
        length,
        multiplicity,
    )


def canonical_flag_record(
    length: int,
    multiplicity: int,
    kernel_dimension: int,
    step: float,
) -> CanonicalFlagRecord:
    """Audit the first Schur promotion on an actual Hardy endpoint."""

    first, _ = deterministic_transverse(
        length,
        multiplicity,
        kernel_dimension,
    )
    second = deterministic_full_transverse(length, multiplicity)
    inverse_gram = (
        2 * np.eye(length * multiplicity, dtype=complex)
        + step * first
        + step**2 * second
    )
    hermitian = np.linalg.inv(inverse_gram)
    positive_endpoint = -canonical_endpoint(
        hermitian,
        length,
        multiplicity,
    )

    unitary = deterministic_copy_unitary(multiplicity)
    active_dimension = multiplicity - kernel_dimension
    ordered = np.hstack(
        (
            unitary[:, kernel_dimension:],
            unitary[:, :kernel_dimension],
        )
    )
    rotated = ordered.conj().T @ positive_endpoint @ ordered
    active = rotated[:active_dimension, :active_dimension]
    cross = rotated[:active_dimension, active_dimension:]
    kernel = rotated[active_dimension:, active_dimension:]
    reduced = kernel - cross.conj().T @ np.linalg.solve(active, cross)

    full_minimum = float(np.linalg.eigvalsh(positive_endpoint)[0])
    active_minimum = float(np.linalg.eigvalsh(active)[0] / step**2)
    reduced_scaled = (reduced + reduced.conj().T) / (2 * step**4)
    reduced_minimum = float(np.linalg.eigvalsh(reduced_scaled)[0])
    reduced_norm = float(np.linalg.norm(reduced_scaled))
    verified = bool(
        full_minimum > -2e-10
        and active_minimum > 1e-5
        and reduced_minimum > -2e-4
        and reduced_norm > 1e-4
    )
    if not verified:
        raise RuntimeError("the canonical Schur flag audit failed")
    return CanonicalFlagRecord(
        record_kind="canonical_repeated_schur_flag",
        length=length,
        multiplicity=multiplicity,
        first_kernel_dimension=kernel_dimension,
        step=format_float(step),
        full_endpoint_minimum_eigenvalue=format_float(full_minimum),
        active_scaled_minimum_eigenvalue=format_float(active_minimum),
        reduced_scaled_minimum_eigenvalue=format_float(reduced_minimum),
        reduced_scaled_norm=format_float(reduced_norm),
        all_checks_passed=verified,
    )


def write_records(path: Path, records: Sequence[AuditRecord]) -> None:
    """Write deterministic JSON Lines output atomically."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(asdict(record), sort_keys=True))
            handle.write("\n")
    temporary.replace(path)


def main() -> None:
    """Run exact and canonical repeated Schur-flag audits."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--minimum-length", type=int, default=2)
    parser.add_argument("--maximum-length", type=int, default=5)
    parser.add_argument(
        "--multiplicities",
        type=int,
        nargs="+",
        default=(3, 4),
    )
    parser.add_argument("--step", type=float, default=3e-3)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/repeated_crabb_schur_flag_s70224.jsonl"
        ),
    )
    args = parser.parse_args()

    records: list[AuditRecord] = [symbolic_flag_record()]
    print("verified exact symbolic Schur flag", flush=True)
    for multiplicity in args.multiplicities:
        for kernel_dimension in range(1, multiplicity):
            for length in range(
                args.minimum_length,
                args.maximum_length + 1,
            ):
                record = canonical_flag_record(
                    length,
                    multiplicity,
                    kernel_dimension,
                    args.step,
                )
                records.append(record)
                print(
                    "verified canonical Schur flag "
                    f"length {length}, multiplicity {multiplicity}, "
                    f"kernel {kernel_dimension}",
                    flush=True,
                )
    write_records(args.output, records)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
