#!/usr/bin/env python3
"""Audit the genuine matrix-inner Faber coordinates at repeated Crabb.

L193's canonical Stein metric turns every inverse-block-Toeplitz equality
anchor into a finite pure partial isometry ``C``.  If ``V`` and ``W`` span
``ker(C)`` and ``ker(C*)``, respectively, then

    B(z) = W* (I - z C*)^{-1} V

is the matrix-inner transfer that replaces the false raw matrix quotient
disproved by L200.  It has ``B(0)=0`` and Taylor coefficients
``B_n = W* (C*)^n V``.  Scalar Dickson polynomials then give the exact
coefficientwise reflection

    F_c B(zeta + c/zeta) = B(zeta) + B(c/zeta).

This checker validates the partial-isometry defects, circle innerness,
Parseval Gram, Faber reflection, and the pure monomial transfer at the
Crabb apex.  It does not identify the reflected Gram with the full
elliptic similarity endpoint.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.linalg import null_space, solve_discrete_lyapunov, sqrtm

from crabb_block_hardy_equality import (
    block_disk_model,
    format_float,
    inverse_block_toeplitz,
)


@dataclass(frozen=True)
class InnerFaberRecord:
    """One genuine colligation-transfer audit."""

    length: int
    multiplicity: int
    toeplitz_strength: str
    spectral_radius: str
    right_defect_projection_error: str
    left_defect_projection_error: str
    defect_basis_error: str
    zero_coefficient_error: str
    circle_inner_error: str
    parseval_error: str
    faber_reflection_error: str
    reflected_gram_minimum_eigenvalue: str
    apex_monomial_error: str
    all_checks_passed: bool


@dataclass(frozen=True)
class TransferData:
    """Canonical partial isometry and its two defect bases."""

    contraction: np.ndarray
    right_defect_basis: np.ndarray
    left_defect_basis: np.ndarray
    right_defect_error: float
    left_defect_error: float
    basis_error: float


def strengthened_inverse_toeplitz(
    length: int,
    multiplicity: int,
    requested_strength: float,
) -> tuple[np.ndarray, float]:
    """Scale L193's deterministic noncommuting legs while preserving positivity."""

    _, coefficients = inverse_block_toeplitz(length, multiplicity)
    strength = requested_strength
    identity = np.eye(length * multiplicity, dtype=complex)
    while True:
        inverse = 2 * identity.copy()
        for offset, coefficient in enumerate(coefficients, start=1):
            block = strength * coefficient
            for level in range(length - offset):
                row = slice(
                    level * multiplicity,
                    (level + 1) * multiplicity,
                )
                column = slice(
                    (level + offset) * multiplicity,
                    (level + offset + 1) * multiplicity,
                )
                inverse[row, column] = block
                inverse[column, row] = block.conj().T
        if np.linalg.eigvalsh(inverse)[0] > 0.45:
            return inverse, strength
        strength *= 0.8


def canonical_transfer_data(
    hermitian: np.ndarray,
    length: int,
    multiplicity: int,
) -> TransferData:
    """Construct L193's canonical partial isometry and defect frames."""

    extended, _, operator, _ = block_disk_model(
        hermitian,
        length,
        multiplicity,
    )
    main_block = hermitian[:multiplicity, :multiplicity]
    main_root = np.asarray(sqrtm(main_block), dtype=complex)
    forcing = (
        extended[:, :multiplicity] @ np.linalg.inv(main_root)
    )
    gramian = solve_discrete_lyapunov(
        operator.conj().T,
        forcing @ forcing.conj().T,
    )
    gramian = (gramian + gramian.conj().T) / 2
    root = np.asarray(sqrtm(gramian), dtype=complex)
    inverse_root = np.linalg.inv(root)
    contraction = root @ operator @ inverse_root

    right_basis = inverse_root @ forcing
    right_basis = right_basis @ np.linalg.inv(
        np.asarray(
            sqrtm(right_basis.conj().T @ right_basis),
            dtype=complex,
        )
    )
    left_basis = null_space(contraction.conj().T)
    if left_basis.shape[1] != multiplicity:
        raise RuntimeError("the left defect has the wrong multiplicity")

    dimension = len(contraction)
    identity = np.eye(dimension, dtype=complex)
    right_projection = right_basis @ right_basis.conj().T
    left_projection = left_basis @ left_basis.conj().T
    right_defect = identity - contraction.conj().T @ contraction
    left_defect = identity - contraction @ contraction.conj().T
    right_error = float(np.linalg.norm(right_defect - right_projection))
    left_error = float(np.linalg.norm(left_defect - left_projection))
    basis_error = max(
        float(
            np.linalg.norm(
                right_basis.conj().T @ right_basis
                - np.eye(multiplicity)
            )
        ),
        float(
            np.linalg.norm(
                left_basis.conj().T @ left_basis
                - np.eye(multiplicity)
            )
        ),
    )
    return TransferData(
        contraction=contraction,
        right_defect_basis=right_basis,
        left_defect_basis=left_basis,
        right_defect_error=right_error,
        left_defect_error=left_error,
        basis_error=basis_error,
    )


def transfer_value(data: TransferData, value: complex) -> np.ndarray:
    """Evaluate ``B(z)=W*(I-zC*)^{-1}V``."""

    contraction = data.contraction
    identity = np.eye(len(contraction), dtype=complex)
    return (
        data.left_defect_basis.conj().T
        @ np.linalg.solve(
            identity - value * contraction.conj().T,
            data.right_defect_basis,
        )
    )


def transfer_coefficients(
    data: TransferData,
    count: int,
) -> tuple[np.ndarray, ...]:
    """Return Taylor coefficients ``W*(C*)^nV``."""

    coefficients = []
    orbit = data.right_defect_basis.copy()
    for _ in range(count):
        coefficients.append(data.left_defect_basis.conj().T @ orbit)
        orbit = data.contraction.conj().T @ orbit
    return tuple(coefficients)


def dickson_faber_sum(
    coefficients: Sequence[np.ndarray],
    boundary: complex,
    ellipse: complex,
) -> np.ndarray:
    """Evaluate the truncated matrix Dickson/Faber series."""

    if len(coefficients) < 2:
        raise ValueError("at least two transfer coefficients are required")
    previous = 2 + 0j
    current = boundary
    total = coefficients[1] * current
    for degree in range(2, len(coefficients)):
        following = boundary * current - ellipse * previous
        total += coefficients[degree] * following
        previous, current = current, following
    return total


def transfer_statistics(
    data: TransferData,
    length: int,
    multiplicity: int,
    ellipse: float,
) -> tuple[float, float, float, float, float]:
    """Return zero, inner, Parseval, and Faber errors."""

    zero_error = float(np.linalg.norm(transfer_value(data, 0)))
    identity = np.eye(multiplicity, dtype=complex)
    inner_error = 0.0
    for index in range(97):
        value = np.exp(2j * np.pi * (index + 0.23) / 97)
        transfer = transfer_value(data, value)
        inner_error = max(
            inner_error,
            float(np.linalg.norm(transfer.conj().T @ transfer - identity)),
        )

    coefficient_count = 240
    coefficients = transfer_coefficients(data, coefficient_count)
    parseval = sum(
        coefficient.conj().T @ coefficient
        for coefficient in coefficients
    )
    parseval_error = float(np.linalg.norm(parseval - identity))

    faber_error = 0.0
    for index in range(29):
        zeta = np.exp(2j * np.pi * (index + 0.31) / 29)
        boundary = zeta + ellipse / zeta
        transformed = dickson_faber_sum(
            coefficients,
            boundary,
            ellipse,
        )
        expected = transfer_value(data, zeta) + transfer_value(
            data,
            ellipse / zeta,
        )
        faber_error = max(
            faber_error,
            float(np.linalg.norm(transformed - expected)),
        )

    reflected_gram = sum(
        ellipse ** (2 * degree)
        * coefficient.conj().T
        @ coefficient
        for degree, coefficient in enumerate(coefficients)
        if degree > 0
    )
    reflected_minimum = float(np.linalg.eigvalsh(reflected_gram)[0])
    return (
        zero_error,
        inner_error,
        parseval_error,
        faber_error,
        reflected_minimum,
    )


def apex_monomial_error(
    length: int,
    multiplicity: int,
) -> float:
    """Check that the Crabb transfer is one unitary times ``z**L``."""

    hermitian = np.eye(length * multiplicity, dtype=complex) / 2
    data = canonical_transfer_data(hermitian, length, multiplicity)
    coefficients = transfer_coefficients(data, 2 * length + 3)
    dominant = coefficients[length]
    identity = np.eye(multiplicity, dtype=complex)
    error = float(
        np.linalg.norm(dominant.conj().T @ dominant - identity)
    )
    for degree, coefficient in enumerate(coefficients):
        if degree != length:
            error = max(error, float(np.linalg.norm(coefficient)))
    return error


def make_record(
    length: int,
    multiplicity: int,
    strength: float,
) -> InnerFaberRecord:
    """Build one complete genuine-transfer record."""

    inverse, actual_strength = strengthened_inverse_toeplitz(
        length,
        multiplicity,
        strength,
    )
    hermitian = np.linalg.inv(inverse)
    data = canonical_transfer_data(hermitian, length, multiplicity)
    (
        zero_error,
        inner_error,
        parseval_error,
        faber_error,
        reflected_minimum,
    ) = transfer_statistics(
        data,
        length,
        multiplicity,
        0.23,
    )
    monomial_error = apex_monomial_error(length, multiplicity)
    spectral_radius = float(
        np.max(np.abs(np.linalg.eigvals(data.contraction)))
    )
    verified = bool(
        spectral_radius < 0.95
        and data.right_defect_error < 2e-10
        and data.left_defect_error < 2e-10
        and data.basis_error < 2e-11
        and zero_error < 2e-11
        and inner_error < 2e-9
        and parseval_error < 2e-9
        and faber_error < 2e-8
        and reflected_minimum > 1e-12
        and monomial_error < 2e-10
    )
    if not verified:
        raise RuntimeError("the genuine inner-Faber transfer audit failed")
    return InnerFaberRecord(
        length=length,
        multiplicity=multiplicity,
        toeplitz_strength=format_float(actual_strength),
        spectral_radius=format_float(spectral_radius),
        right_defect_projection_error=format_float(
            data.right_defect_error
        ),
        left_defect_projection_error=format_float(
            data.left_defect_error
        ),
        defect_basis_error=format_float(data.basis_error),
        zero_coefficient_error=format_float(zero_error),
        circle_inner_error=format_float(inner_error),
        parseval_error=format_float(parseval_error),
        faber_reflection_error=format_float(faber_error),
        reflected_gram_minimum_eigenvalue=format_float(
            reflected_minimum
        ),
        apex_monomial_error=format_float(monomial_error),
        all_checks_passed=verified,
    )


def write_records(
    path: Path,
    records: Sequence[InnerFaberRecord],
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
    """Run all genuine inner-Faber transfer audits."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--minimum-length", type=int, default=2)
    parser.add_argument("--maximum-length", type=int, default=5)
    parser.add_argument(
        "--multiplicities",
        type=int,
        nargs="+",
        default=(2, 3),
    )
    parser.add_argument("--strength", type=float, default=18.0)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_inner_faber_transfer_s70224.jsonl"
        ),
    )
    args = parser.parse_args()

    records = [
        make_record(length, multiplicity, args.strength)
        for multiplicity in args.multiplicities
        for length in range(args.minimum_length, args.maximum_length + 1)
    ]
    write_records(args.output, records)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
