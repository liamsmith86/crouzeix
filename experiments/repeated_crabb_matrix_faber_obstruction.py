#!/usr/bin/env python3
"""Audit the matrix-Faber obstruction at repeated Crabb equality anchors.

The obstruction already occurs for repeated ``C_3`` (Crabb length two).
For a normalized positive block Toeplitz inverse Gram with main block
``2I`` and upper coefficient ``Z``, its equality operator has an exact
block-companion form.  Put ``C=Z/2``.  Its terminal matrix polynomial is

    G(z) = z**2 I - C* z.

The scalar proof uses the finite Blaschke quotient ``g/g#``.  A tempting
matrix lift replaces it by either ordered quotient of ``G`` and its
reciprocal adjoint ``G#``.  Both quotients fail to be unitary on the circle
for nonnormal copy coefficients, even though the coefficientwise
Dickson/Faber reflection identity remains true.

This checker verifies the exact companion and Faber identities, then
adversarially verifies the failure of both ordered inner quotients.  It also
checks that ``X + c X*`` still has the expected repeated elliptic support.
The calculations diagnose a false proof route; they do not test or
disprove the desired similarity bound.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.linalg import sqrtm

from crabb_block_hardy_equality import (
    block_disk_model,
    format_float,
)


@dataclass(frozen=True)
class MatrixFaberRecord:
    """One block-companion/Faber obstruction audit."""

    length: int
    multiplicity: int
    companion_error: str
    characteristic_error: str
    faber_reflection_error: str
    coefficient_nonnormality: str
    right_inner_error: str
    left_inner_error: str
    ellipse_support_error: str
    ellipse_support_multiplicity: int
    all_checks_passed: bool


def deterministic_coefficients(
    length: int,
    multiplicity: int,
    seed: int,
) -> tuple[np.ndarray, ...]:
    """Return bounded, genuinely nonnormal block Toeplitz coefficients."""

    rng = np.random.default_rng(seed)
    coefficients = []
    for offset in range(1, length):
        block = (
            rng.standard_normal((multiplicity, multiplicity))
            + 1j * rng.standard_normal((multiplicity, multiplicity))
        ) / np.sqrt(2 * length * multiplicity)
        if offset == 1:
            block += 0.7 * np.diag(
                np.ones(multiplicity - 1, dtype=complex),
                1,
            )
        coefficients.append(block)

    scale = 0.45
    while True:
        inverse = inverse_block_toeplitz(
            length,
            multiplicity,
            tuple(scale * block for block in coefficients),
        )
        if np.linalg.eigvalsh(inverse)[0] > 0.5:
            return tuple(scale * block for block in coefficients)
        scale *= 0.8


def inverse_block_toeplitz(
    length: int,
    multiplicity: int,
    coefficients: Sequence[np.ndarray],
) -> np.ndarray:
    """Build a Hermitian block Toeplitz matrix with main block ``2I``."""

    inverse = 2 * np.eye(length * multiplicity, dtype=complex)
    for offset, block in enumerate(coefficients, start=1):
        for level in range(length - offset):
            row = slice(level * multiplicity, (level + 1) * multiplicity)
            column = slice(
                (level + offset) * multiplicity,
                (level + offset + 1) * multiplicity,
            )
            inverse[row, column] = block
            inverse[column, row] = block.conj().T
    return inverse


def expected_companion(
    length: int,
    multiplicity: int,
    coefficients: Sequence[np.ndarray],
) -> np.ndarray:
    """Return the exact normalized repeated-``C_3`` equality companion."""

    if length != 2:
        raise ValueError("the exact companion audit is intentionally length two")

    dimension = (length + 1) * multiplicity
    companion = np.zeros((dimension, dimension), dtype=complex)
    identity = np.eye(multiplicity, dtype=complex)
    companion[:multiplicity, multiplicity : 2 * multiplicity] = 2 * identity
    for level in range(1, length):
        companion[
            level * multiplicity : (level + 1) * multiplicity,
            (level + 1) * multiplicity : (level + 2) * multiplicity,
        ] = identity
    for offset, block in enumerate(coefficients, start=1):
        companion[
            :multiplicity,
            (offset + 1) * multiplicity : (offset + 2) * multiplicity,
        ] = -block / 2
        column = length - offset + 1
        companion[
            length * multiplicity :,
            column * multiplicity : (column + 1) * multiplicity,
        ] = block.conj().T / 2
    return companion


def matrix_polynomial(
    value: complex,
    length: int,
    multiplicity: int,
    coefficients: Sequence[np.ndarray],
) -> np.ndarray:
    """Evaluate the terminal matrix polynomial ``G``."""

    result = value**length * np.eye(multiplicity, dtype=complex)
    for offset, block in enumerate(coefficients, start=1):
        result -= block.conj().T * value ** (length - offset) / 2
    return result


def reciprocal_polynomial(
    value: complex,
    multiplicity: int,
    coefficients: Sequence[np.ndarray],
) -> np.ndarray:
    """Evaluate ``G#(z)=z^L G(1/conj(z))*`` coefficientwise."""

    result = np.eye(multiplicity, dtype=complex)
    for offset, block in enumerate(coefficients, start=1):
        result -= block * value**offset / 2
    return result


def dickson_values(
    maximum_degree: int,
    value: complex,
    ellipse: complex,
) -> list[complex]:
    """Return monic Dickson polynomials through ``maximum_degree``."""

    values = [2 + 0j]
    if maximum_degree == 0:
        return values
    values.append(value)
    for _ in range(2, maximum_degree + 1):
        values.append(value * values[-1] - ellipse * values[-2])
    return values


def characteristic_error(
    operator: np.ndarray,
    length: int,
    multiplicity: int,
    coefficients: Sequence[np.ndarray],
) -> float:
    """Check ``det(zI-A)=z^m det G(z)`` at fixed safe points."""

    identity = np.eye(len(operator), dtype=complex)
    worst = 0.0
    for value in (0.63 + 0.41j, -0.72 + 0.29j, 1.17 - 0.23j):
        left = np.linalg.det(value * identity - operator)
        right = value**multiplicity * np.linalg.det(
            matrix_polynomial(
                value,
                length,
                multiplicity,
                coefficients,
            )
        )
        worst = max(worst, abs(left - right))
    return float(worst)


def faber_reflection_error(
    length: int,
    multiplicity: int,
    coefficients: Sequence[np.ndarray],
    ellipse: complex,
) -> float:
    """Check the valid coefficientwise Dickson/Faber identity."""

    identity = np.eye(multiplicity, dtype=complex)
    worst = 0.0
    for index in range(19):
        zeta = np.exp(2j * np.pi * (index + 0.31) / 19)
        boundary = zeta + ellipse / zeta
        dickson = dickson_values(length, boundary, ellipse)
        transformed = dickson[length] * identity
        positive = matrix_polynomial(
            zeta,
            length,
            multiplicity,
            coefficients,
        )
        negative = ellipse**length * zeta ** (-length) * identity
        for offset, block in enumerate(coefficients, start=1):
            degree = length - offset
            transformed -= block.conj().T * dickson[degree] / 2
            negative -= (
                block.conj().T
                * ellipse**degree
                * zeta ** (-degree)
                / 2
            )
        worst = max(worst, float(np.linalg.norm(transformed - positive - negative)))
    return worst


def inner_quotient_errors(
    length: int,
    multiplicity: int,
    coefficients: Sequence[np.ndarray],
) -> tuple[float, float]:
    """Return the maximum unitary defects of both scalar-like quotients."""

    identity = np.eye(multiplicity, dtype=complex)
    right_error = 0.0
    left_error = 0.0
    for index in range(23):
        value = np.exp(2j * np.pi * (index + 0.17) / 23)
        polynomial = matrix_polynomial(
            value,
            length,
            multiplicity,
            coefficients,
        )
        reciprocal = reciprocal_polynomial(
            value,
            multiplicity,
            coefficients,
        )
        right = polynomial @ np.linalg.inv(reciprocal)
        left = np.linalg.solve(reciprocal, polynomial)
        right_error = max(
            right_error,
            float(np.linalg.norm(right.conj().T @ right - identity)),
        )
        left_error = max(
            left_error,
            float(np.linalg.norm(left.conj().T @ left - identity)),
        )
    return right_error, left_error


def ellipse_support_statistics(
    operator: np.ndarray,
    multiplicity: int,
    ellipse: float,
) -> tuple[float, int]:
    """Check the exact elliptic support formula and top multiplicity."""

    elliptic = operator + ellipse * operator.conj().T
    worst = 0.0
    minimum_multiplicity = len(operator)
    for index in range(73):
        angle = 2 * np.pi * (index + 0.23) / 73
        phase = np.exp(-1j * angle)
        support = (
            phase * elliptic
            + np.conjugate(phase) * elliptic.conj().T
        ) / 2
        eigenvalues = np.linalg.eigvalsh(support)
        expected = abs(phase + ellipse * np.conjugate(phase))
        worst = max(worst, abs(float(eigenvalues[-1]) - expected))
        count = int(np.sum(abs(eigenvalues - expected) < 2e-9))
        minimum_multiplicity = min(minimum_multiplicity, count)
    return worst, minimum_multiplicity


def make_record(
    length: int,
    multiplicity: int,
    seed: int,
) -> MatrixFaberRecord:
    """Build one complete obstruction record."""

    if length != 2:
        raise ValueError("L200's obstruction record is intentionally length two")
    coefficients = deterministic_coefficients(length, multiplicity, seed)
    inverse = inverse_block_toeplitz(length, multiplicity, coefficients)
    hermitian = np.linalg.inv(inverse)
    _, coordinate, operator, _ = block_disk_model(
        hermitian,
        length,
        multiplicity,
    )
    expected = expected_companion(length, multiplicity, coefficients)
    companion_error = float(np.linalg.norm(operator - expected))
    determinant_error = characteristic_error(
        operator,
        length,
        multiplicity,
        coefficients,
    )
    reflection_error = faber_reflection_error(
        length,
        multiplicity,
        coefficients,
        0.19 + 0.07j,
    )
    first = coefficients[0]
    nonnormality = float(
        np.linalg.norm(
            first @ first.conj().T - first.conj().T @ first
        )
    )
    right_error, left_error = inner_quotient_errors(
        length,
        multiplicity,
        coefficients,
    )

    root = np.asarray(sqrtm(coordinate), dtype=complex)
    physical = root @ operator @ np.linalg.inv(root)
    support_error, support_multiplicity = ellipse_support_statistics(
        physical,
        multiplicity,
        0.27,
    )
    verified = bool(
        companion_error < 2e-11
        and determinant_error < 2e-10
        and reflection_error < 2e-10
        and nonnormality > 1e-3
        and right_error > 1e-3
        and left_error > 1e-3
        and support_error < 2e-10
        and support_multiplicity == multiplicity
    )
    if not verified:
        raise RuntimeError("the matrix-Faber obstruction audit failed")

    return MatrixFaberRecord(
        length=length,
        multiplicity=multiplicity,
        companion_error=format_float(companion_error),
        characteristic_error=format_float(determinant_error),
        faber_reflection_error=format_float(reflection_error),
        coefficient_nonnormality=format_float(nonnormality),
        right_inner_error=format_float(right_error),
        left_inner_error=format_float(left_error),
        ellipse_support_error=format_float(support_error),
        ellipse_support_multiplicity=support_multiplicity,
        all_checks_passed=verified,
    )


def write_records(
    path: Path,
    records: Sequence[MatrixFaberRecord],
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
    """Run the block companion/Faber obstruction audits."""

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--multiplicities",
        type=int,
        nargs="+",
        default=(2, 3, 4),
    )
    parser.add_argument("--seed", type=int, default=70224)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_matrix_faber_obstruction_s70224.jsonl"
        ),
    )
    args = parser.parse_args()

    records = [
        make_record(
            2,
            multiplicity,
            args.seed + 17 * multiplicity,
        )
        for multiplicity in args.multiplicities
    ]
    write_records(args.output, records)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
