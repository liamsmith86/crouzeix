#!/usr/bin/env python3
"""Audit elliptic first-jet stationarity at block-Hardy equality.

At every inverse-block-Toeplitz equality anchor, the normalized Hardy
metric has only the eigenvalues 1, 2, and 4.  After balancing by that
metric, the disk operator is a partial isometry with copy-dimensional
right and left defects.  This structure makes the Stein-defect derivative
of the normalized elliptic tangent ``T* - T^3`` vanish on the old defect
kernel.

The calculation is fully matrix ordered and the test anchors have
noncommuting block-Toeplitz coefficients.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve_discrete_lyapunov, sqrtm

from crabb_block_hardy_equality import (
    block_disk_model,
    format_float,
)
from repeated_crabb_inner_faber_transfer import (
    strengthened_inverse_toeplitz,
)


@dataclass(frozen=True)
class EllipticFirstJetRecord:
    """One block-Hardy first-jet audit."""

    length: int
    multiplicity: int
    toeplitz_strength: str
    endpoint_energy_identity_error: str
    metric_three_level_error: str
    right_partial_isometry_error: str
    left_partial_isometry_error: str
    elliptic_kernel_compression_error: str
    elliptic_left_kernel_compression_error: str
    defect_row_factorization_error: str
    coefficient_commutator_norm: str
    all_checks_passed: bool


def endpoint_embedding(
    length: int,
    multiplicity: int,
    level: int,
) -> np.ndarray:
    """Return the copy-space embedding at one level."""

    embedding = np.zeros(
        ((length + 1) * multiplicity, multiplicity),
        dtype=complex,
    )
    embedding[
        level * multiplicity : (level + 1) * multiplicity,
        :,
    ] = np.eye(multiplicity)
    return embedding


def normalized_columns(columns: np.ndarray) -> np.ndarray:
    """Normalize a full-rank column block to an isometry."""

    gram = columns.conj().T @ columns
    root = np.asarray(sqrtm(gram), dtype=complex)
    return columns @ np.linalg.inv(root)


def physical_hardy_data(
    length: int,
    multiplicity: int,
    requested_strength: float,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    float,
    float,
    float,
]:
    """Return physical Hardy data and two structural diagnostics."""

    inverse, actual_strength = strengthened_inverse_toeplitz(
        length,
        multiplicity,
        requested_strength,
    )
    hermitian = np.linalg.inv(inverse)
    extended, coordinate, operator, _ = block_disk_model(
        hermitian,
        length,
        multiplicity,
    )
    coordinate_root = np.asarray(sqrtm(coordinate), dtype=complex)
    coordinate_root_inverse = np.linalg.inv(coordinate_root)
    physical = coordinate_root @ operator @ coordinate_root_inverse

    initial = endpoint_embedding(length, multiplicity, 0)
    terminal = endpoint_embedding(length, multiplicity, length)
    corner = initial.conj().T @ extended @ initial
    forcing = (
        extended
        @ initial
        @ np.linalg.inv(np.asarray(sqrtm(corner), dtype=complex))
    )
    terminal_corner = terminal.conj().T @ coordinate @ terminal
    terminal_forcing = (
        coordinate
        @ terminal
        @ np.linalg.inv(
            np.asarray(sqrtm(terminal_corner), dtype=complex)
        )
    )
    endpoint_energy = (
        coordinate
        - operator.conj().T @ coordinate @ operator
        - forcing @ forcing.conj().T
        + operator.conj().T
        @ forcing
        @ forcing.conj().T
        @ operator
        / 2
        + terminal_forcing @ terminal_forcing.conj().T
    )
    endpoint_energy_error = float(np.linalg.norm(endpoint_energy))
    metric = solve_discrete_lyapunov(
        operator.conj().T,
        forcing @ forcing.conj().T,
    )
    physical_metric = (
        coordinate_root_inverse.conj().T
        @ metric
        @ coordinate_root_inverse
    )
    right_defect = normalized_columns(coordinate_root @ initial)
    left_defect = normalized_columns(coordinate_root @ terminal)

    first_coefficient = inverse[
        :multiplicity,
        multiplicity : 2 * multiplicity,
    ]
    if length >= 3:
        second_coefficient = inverse[
            :multiplicity,
            2 * multiplicity : 3 * multiplicity,
        ]
        commutator = (
            first_coefficient @ second_coefficient
            - second_coefficient @ first_coefficient
        )
    else:
        commutator = (
            first_coefficient @ first_coefficient.conj().T
            - first_coefficient.conj().T @ first_coefficient
        )
    return (
        physical,
        physical_metric,
        right_defect,
        left_defect,
        endpoint_energy_error,
        float(np.linalg.norm(commutator)),
        actual_strength,
    )


def defect_row_factorization(
    derivative: np.ndarray,
    defect: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Factor a Hermitian matrix supported on one defect row and column."""

    projection = defect @ defect.conj().T
    corner = defect.conj().T @ derivative @ defect
    column = (
        (np.eye(len(derivative)) - projection) @ derivative @ defect
        + defect @ corner / 2
    )
    reconstructed = defect @ column.conj().T + column @ defect.conj().T
    return column, float(np.linalg.norm(derivative - reconstructed))


def make_record(
    length: int,
    multiplicity: int,
    strength: float,
) -> EllipticFirstJetRecord:
    """Build one exact-structure first-jet audit."""

    (
        operator,
        metric,
        right,
        left,
        endpoint_energy_error,
        commutator,
        actual_strength,
    ) = physical_hardy_data(length, multiplicity, strength)
    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    identity = np.eye(len(operator), dtype=complex)
    expected_metric = (
        2 * identity - right_projection + 2 * left_projection
    )
    metric_error = float(np.linalg.norm(metric - expected_metric))

    metric_root = np.asarray(sqrtm(metric), dtype=complex)
    balanced = metric_root @ operator @ np.linalg.inv(metric_root)
    right_error = float(
        np.linalg.norm(
            identity - balanced.conj().T @ balanced - right_projection
        )
    )
    left_error = float(
        np.linalg.norm(
            identity - balanced @ balanced.conj().T - left_projection
        )
    )

    kernel_projection = identity - right_projection
    compression_error = 0.0
    left_compression_error = 0.0
    factorization_error = 0.0
    balanced_left_kernel = identity - left_projection
    for phase in (1 + 0j, 1j, np.exp(0.37j)):
        elliptic_tangent = (
            phase * operator.conj().T
            - np.conjugate(phase) * np.linalg.matrix_power(operator, 3)
        )
        defect_derivative = -(
            elliptic_tangent.conj().T @ metric @ operator
            + operator.conj().T @ metric @ elliptic_tangent
        )
        compression_error = max(
            compression_error,
            float(
                np.linalg.norm(
                    kernel_projection
                    @ defect_derivative
                    @ kernel_projection
                )
            ),
        )
        balanced_tangent = (
            phase * metric @ balanced.conj().T @ np.linalg.inv(metric)
            - np.conjugate(phase)
            * np.linalg.matrix_power(balanced, 3)
        )
        left_compression_error = max(
            left_compression_error,
            float(
                np.linalg.norm(
                    balanced_left_kernel
                    @ (
                        balanced_tangent @ balanced.conj().T
                        + balanced @ balanced_tangent.conj().T
                    )
                    @ balanced_left_kernel
                )
            ),
        )
        _, current_factorization_error = defect_row_factorization(
            defect_derivative,
            right,
        )
        factorization_error = max(
            factorization_error,
            current_factorization_error,
        )

    verified = bool(
        endpoint_energy_error < 2e-10
        and metric_error < 2e-10
        and right_error < 2e-10
        and left_error < 2e-10
        and compression_error < 2e-10
        and left_compression_error < 2e-10
        and factorization_error < 2e-10
        and (multiplicity == 1 or commutator > 1e-8)
    )
    if not verified:
        raise RuntimeError("the repeated elliptic first-jet audit failed")
    return EllipticFirstJetRecord(
        length=length,
        multiplicity=multiplicity,
        toeplitz_strength=format_float(actual_strength),
        endpoint_energy_identity_error=format_float(
            endpoint_energy_error
        ),
        metric_three_level_error=format_float(metric_error),
        right_partial_isometry_error=format_float(right_error),
        left_partial_isometry_error=format_float(left_error),
        elliptic_kernel_compression_error=format_float(compression_error),
        elliptic_left_kernel_compression_error=format_float(
            left_compression_error
        ),
        defect_row_factorization_error=format_float(factorization_error),
        coefficient_commutator_norm=format_float(commutator),
        all_checks_passed=verified,
    )


def write_records(
    records: list[EllipticFirstJetRecord],
    output: Path,
) -> None:
    """Write deterministic JSON Lines records."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)


def parse_args() -> argparse.Namespace:
    """Parse command-line options."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/repeated_crabb_elliptic_first_jet_s70224.jsonl"
        ),
    )
    parser.add_argument("--strength", type=float, default=18.0)
    return parser.parse_args()


def main() -> None:
    """Run the standard noncommutative audit suite."""

    args = parse_args()
    records = [
        make_record(length, multiplicity, args.strength)
        for multiplicity in (2, 3)
        for length in range(2, 6)
    ]
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
