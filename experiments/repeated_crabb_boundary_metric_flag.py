#!/usr/bin/env python3
"""Audit ordered partial flags of the L219 boundary metric.

The balanced upper gap has an exact transfer-factor expansion.  After
the state complement and every earlier active left-copy range have
been Schur eliminated, its first face on

    intersection(kernel(B_j*), j < k)

is ``B_k B_k*`` at order ``q**k``.  Physical coordinates multiply
this face by four.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.linalg import null_space

from crabb_block_hardy_equality import format_float
from repeated_crabb_delayed_jet import inverse_series, series_multiply
from repeated_crabb_matrix_schur_chart import schur_value
from repeated_crabb_schur_kernel_flag import schur_features
from repeated_crabb_transfer_channel_covariance import (
    transfer_coefficient,
)
from repeated_crabb_transfer_deflation import heterogeneous_shift


Matrix = np.ndarray


@dataclass(frozen=True)
class BoundaryMetricFlagRecord:
    """One ordered partial-flag face audit."""

    construction_kind: str
    multiplicity: int
    state_dimension: int
    grade: int
    flag_dimension: int
    active_rank: int
    earlier_flagged_transfer_error: str
    leading_valuation: int
    leading_face_error: str
    leading_positive_eigenvalue: str
    colligation_error: str
    all_checks_passed: bool


def lossless_schur_realization(
    parameters: list[Matrix],
    terminal: Matrix,
) -> tuple[Matrix, Matrix, Matrix, float]:
    """Realize a finite matrix-inner Schur recursion as a partial isometry."""

    length = len(parameters)
    multiplicity = len(terminal)
    model_dimension = length * multiplicity
    sources = []
    targets = []
    for sample in range(length + 1):
        value = 0.51 * np.exp(
            2j * np.pi * (sample + 0.173) / (length + 1)
        )
        features = np.vstack(
            schur_features(parameters, terminal, value)
        )
        function = schur_value(parameters, terminal, value)
        sources.append(
            np.vstack((np.eye(multiplicity), value * features))
        )
        targets.append(np.vstack((features, function)))

    source = np.hstack(sources)
    target = np.hstack(targets)
    colligation = np.linalg.solve(source.T, target.T).T
    state = colligation[:model_dimension, multiplicity:]
    input_column = colligation[:model_dimension, :multiplicity]
    output_row = colligation[model_dimension:, multiplicity:]
    feedthrough = colligation[model_dimension:, :multiplicity]

    dimension = model_dimension + multiplicity
    adjoint_partial = np.block(
        [
            [state, input_column],
            [
                np.zeros((multiplicity, model_dimension), dtype=complex),
                np.zeros((multiplicity, multiplicity), dtype=complex),
            ],
        ]
    )
    partial = adjoint_partial.conj().T
    right = np.vstack(
        (
            np.zeros((model_dimension, multiplicity), dtype=complex),
            np.eye(multiplicity, dtype=complex),
        )
    )
    left = np.vstack((output_row.conj().T, feedthrough.conj().T))

    identity = np.eye(dimension, dtype=complex)
    errors = (
        np.linalg.norm(
            colligation.conj().T @ colligation
            - np.eye(len(colligation), dtype=complex)
        ),
        np.linalg.norm(
            identity
            - partial.conj().T @ partial
            - right @ right.conj().T
        ),
        np.linalg.norm(
            identity
            - partial @ partial.conj().T
            - left @ left.conj().T
        ),
        np.linalg.norm(right.conj().T @ left),
    )
    return partial, right, left, float(max(errors))


def balanced_upper_gap_coefficients(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    degree: int,
) -> list[Matrix]:
    """Return ``4 P^-1 - P_bl(q)`` through the requested q-degree."""

    dimension = len(partial)
    identity = np.eye(dimension, dtype=complex)
    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    coefficients = [
        np.zeros((dimension, dimension), dtype=complex)
        for _ in range(degree + 1)
    ]
    coefficients[0] = 2 * right_projection + identity - left_projection

    partial_powers = [
        np.eye(dimension, dtype=complex)
    ]
    adjoint_powers = [
        np.eye(dimension, dtype=complex)
    ]
    for _ in range(degree):
        partial_powers.append(partial_powers[-1] @ partial)
        adjoint_powers.append(adjoint_powers[-1] @ partial.conj().T)

    right_orbits = [
        adjoint_powers[index]
        @ right_projection
        @ partial_powers[index]
        for index in range(degree + 1)
    ]
    left_orbits = [
        partial_powers[index]
        @ left_projection
        @ adjoint_powers[index]
        for index in range(degree + 1)
    ]
    for order in range(1, degree + 1):
        coefficients[order] -= left_orbits[order]
        for divisor in range(1, order + 1):
            if order % divisor == 0:
                coefficients[order] += (
                    (-1) ** (order // divisor - 1)
                    * right_orbits[divisor]
                )
    return coefficients


def full_endpoint_gap_series(
    gap_coefficients: list[Matrix],
    left: Matrix,
) -> tuple[list[Matrix], float]:
    """Short the balanced upper gap to the complete left-copy endpoint."""

    degree = len(gap_coefficients) - 1
    complement = null_space(left.conj().T)
    interior = [
        complement.conj().T @ coefficient @ complement
        for coefficient in gap_coefficients
    ]
    cross = [
        left.conj().T @ coefficient @ complement
        for coefficient in gap_coefficients
    ]
    endpoint = [
        left.conj().T @ coefficient @ left
        for coefficient in gap_coefficients
    ]
    cross_square = series_multiply(
        series_multiply(
            cross,
            inverse_series(interior),
            degree,
        ),
        [coefficient.conj().T for coefficient in cross],
        degree,
    )
    result = []
    for order in range(degree + 1):
        coefficient = endpoint[order] - cross_square[order]
        result.append((coefficient + coefficient.conj().T) / 2)

    constant_error = float(np.linalg.norm(result[0]))
    result[0] = np.zeros_like(result[0])
    return result, constant_error


def first_nonzero_coefficient(
    coefficients: list[Matrix],
    tolerance: float,
) -> tuple[int, Matrix]:
    """Return the valuation and leading coefficient of a matrix series."""

    for order, coefficient in enumerate(coefficients):
        if np.linalg.norm(coefficient) > tolerance:
            return order, coefficient
    raise RuntimeError("the truncated endpoint series vanished")


def schur_reduce_leading_kernel(
    coefficients: list[Matrix],
    tolerance: float,
) -> tuple[list[Matrix], Matrix, Matrix, int, Matrix]:
    """Split the leading range and return the formal kernel Schur series."""

    degree = len(coefficients) - 1
    valuation, leading = first_nonzero_coefficient(
        coefficients,
        tolerance,
    )
    eigenvalues, eigenvectors = np.linalg.eigh(leading)
    active = eigenvectors[:, eigenvalues > tolerance]
    kernel = eigenvectors[:, eigenvalues <= tolerance]
    if active.shape[1] == 0:
        raise RuntimeError("the leading endpoint face had zero rank")
    if kernel.shape[1] == 0:
        return [], active, kernel, valuation, leading

    active_series = [
        active.conj().T @ coefficient @ active
        for coefficient in coefficients
    ]
    cross_series = [
        active.conj().T @ coefficient @ kernel
        for coefficient in coefficients
    ]
    kernel_series = [
        kernel.conj().T @ coefficient @ kernel
        for coefficient in coefficients
    ]
    normalized_active = active_series[valuation:]
    normalized_inverse = inverse_series(normalized_active)
    product = series_multiply(
        series_multiply(
            [coefficient.conj().T for coefficient in cross_series],
            normalized_inverse,
            degree,
        ),
        cross_series,
        degree,
    )

    reliable_degree = degree - valuation
    reduced = []
    for order in range(reliable_degree + 1):
        coefficient = (
            kernel_series[order] - product[order + valuation]
        )
        reduced.append((coefficient + coefficient.conj().T) / 2)
    return reduced, active, kernel, valuation, leading


def rank_chain_case(
    multiplicity: int,
    seed: int,
    parameter_scale: float = 1.0,
) -> tuple[Matrix, Matrix, Matrix, float, int]:
    """Return a noncommuting Schur realization with one new left rank per grade."""

    generator = np.random.default_rng(seed)

    def random_unitary() -> Matrix:
        matrix = (
            generator.standard_normal((multiplicity, multiplicity))
            + 1j
            * generator.standard_normal((multiplicity, multiplicity))
        )
        unitary, triangular = np.linalg.qr(matrix)
        phases = np.diag(triangular)
        phases = np.where(np.abs(phases) > 0, phases / np.abs(phases), 1)
        return unitary @ np.diag(np.conjugate(phases))

    left_basis = random_unitary()
    right_basis = random_unitary()
    parameters = [
        np.zeros((multiplicity, multiplicity), dtype=complex)
    ]
    for grade in range(1, multiplicity):
        left_vector = (
            left_basis[:, grade - 1]
            + 0.31 * left_basis[:, grade % multiplicity]
            + 0.17 * left_basis[:, -1]
        )
        right_vector = (
            right_basis[:, (2 * grade - 1) % multiplicity]
            + 0.29 * right_basis[:, grade - 1]
        )
        left_vector /= np.linalg.norm(left_vector)
        right_vector /= np.linalg.norm(right_vector)
        radius = parameter_scale * (0.11 + 0.018 * grade)
        parameters.append(
            radius
            * np.outer(left_vector, np.conjugate(right_vector))
        )

    terminal = random_unitary()
    partial, right, left, error = lossless_schur_realization(
        parameters,
        terminal,
    )
    return partial, right, left, error, multiplicity


def audit_case(
    construction_kind: str,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    colligation_error: float,
    maximum_grade: int,
) -> list[BoundaryMetricFlagRecord]:
    """Audit every active member of one ordered transfer flag."""

    multiplicity = right.shape[1]
    degree = maximum_grade * (maximum_grade + 1) // 2 + 2
    gap = balanced_upper_gap_coefficients(
        partial,
        right,
        left,
        degree,
    )
    endpoint, constant_error = full_endpoint_gap_series(gap, left)
    endpoint = [4 * coefficient for coefficient in endpoint]
    transfers = [
        transfer_coefficient(partial, right, left, grade)
        for grade in range(maximum_grade + 1)
    ]

    coordinates = np.eye(multiplicity, dtype=complex)
    current = endpoint
    records = []
    tolerance = 2e-7
    for grade in range(1, maximum_grade + 1):
        earlier_error = max(
            (
                float(
                    np.linalg.norm(
                        coordinates.conj().T @ transfers[index]
                    )
                )
                for index in range(1, grade)
            ),
            default=0.0,
        )
        active_transfer = coordinates.conj().T @ transfers[grade]
        if np.linalg.norm(active_transfer) < tolerance:
            continue

        (
            reduced,
            active,
            kernel,
            valuation,
            leading,
        ) = schur_reduce_leading_kernel(current, tolerance)
        expected = 4 * active_transfer @ active_transfer.conj().T
        face_error = float(np.linalg.norm(leading - expected))
        positive_eigenvalues = np.linalg.eigvalsh(leading)
        minimum_positive = float(
            positive_eigenvalues[positive_eigenvalues > tolerance][0]
        )
        active_rank = int(np.linalg.matrix_rank(expected, tolerance))
        verified = bool(
            constant_error < 3e-8
            and colligation_error < 3e-8
            and earlier_error < 3e-8
            and valuation == grade
            and active.shape[1] == active_rank
            and face_error < 3e-7
            and minimum_positive > tolerance
        )
        if not verified:
            raise RuntimeError(
                "the boundary-metric flag audit failed: "
                f"kind={construction_kind}, grade={grade}, "
                f"valuation={valuation}, expected={grade}, "
                f"face={face_error:.3e}, earlier={earlier_error:.3e}"
            )
        records.append(
            BoundaryMetricFlagRecord(
                construction_kind=construction_kind,
                multiplicity=multiplicity,
                state_dimension=len(partial),
                grade=grade,
                flag_dimension=coordinates.shape[1],
                active_rank=active_rank,
                earlier_flagged_transfer_error=format_float(
                    earlier_error
                ),
                leading_valuation=valuation,
                leading_face_error=format_float(face_error),
                leading_positive_eigenvalue=format_float(
                    minimum_positive
                ),
                colligation_error=format_float(colligation_error),
                all_checks_passed=verified,
            )
        )
        if kernel.shape[1] == 0:
            coordinates = coordinates[:, :0]
            break
        coordinates = coordinates @ kernel
        current = reduced

    if coordinates.shape[1] != 0:
        raise RuntimeError(
            f"the transfer flag did not terminate for {construction_kind}"
        )
    return records


def standard_records() -> list[BoundaryMetricFlagRecord]:
    """Return noncommuting chains and rank-jumping shift flags."""

    records = []
    for multiplicity in range(2, 6):
        case = rank_chain_case(multiplicity, 108_200 + multiplicity)
        records.extend(
            audit_case(
                "noncommuting_schur_rank_chain",
                *case,
            )
        )

    for index, lengths in enumerate(
        (
            (1, 1, 3, 5),
            (2, 4, 4, 6),
            (3, 3, 3, 5, 7),
        )
    ):
        partial, right, left = heterogeneous_shift(lengths)
        records.extend(
            audit_case(
                f"rank_jumping_shift_{index}",
                partial,
                right,
                left,
                0.0,
                max(lengths),
            )
        )
    return records


def write_records(
    records: Sequence[BoundaryMetricFlagRecord],
    output: Path,
) -> None:
    """Write deterministic JSON Lines atomically."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_boundary_metric_flag_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the complete audit."""

    args = parse_args()
    records = standard_records()
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
