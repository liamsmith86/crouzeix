#!/usr/bin/env python3
"""Audit ordered right-copy flags of the L219 lower metric gap.

The balanced lower gap ``P_bl - P^-1`` has a transfer-factor expansion
dual to L224's upper expansion.  After the state complement and every
earlier active right-copy range are Schur eliminated, its first face on

    intersection(kernel(B_j), j < k)

is ``B_k* B_k`` at order ``q**k``.  Unlike the upper face, physical
coordinates introduce no scalar factor at the right defect.
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
from repeated_crabb_boundary_metric_flag import (
    rank_chain_case,
    schur_reduce_leading_kernel,
)
from repeated_crabb_delayed_jet import inverse_series, series_multiply
from repeated_crabb_transfer_channel_covariance import (
    transfer_coefficient,
)
from repeated_crabb_transfer_deflation import heterogeneous_shift


Matrix = np.ndarray


@dataclass(frozen=True)
class LowerMetricFlagRecord:
    """One ordered lower partial-flag face audit."""

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


def balanced_lower_gap_coefficients(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    degree: int,
) -> list[Matrix]:
    """Return ``P_bl(q) - P^-1`` through the requested q-degree."""

    dimension = len(partial)
    identity = np.eye(dimension, dtype=complex)
    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    metric = 2 * identity - right_projection + 2 * left_projection
    coefficients = [
        np.zeros((dimension, dimension), dtype=complex)
        for _ in range(degree + 1)
    ]
    coefficients[0] = identity - np.linalg.inv(metric)

    partial_powers = [np.eye(dimension, dtype=complex)]
    adjoint_powers = [np.eye(dimension, dtype=complex)]
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
        coefficients[order] += left_orbits[order]
        for divisor in range(1, order + 1):
            if order % divisor == 0:
                coefficients[order] += (
                    (-1) ** (order // divisor)
                    * right_orbits[divisor]
                )
    return coefficients


def full_lower_endpoint_series(
    gap_coefficients: list[Matrix],
    right: Matrix,
) -> tuple[list[Matrix], float]:
    """Short the balanced lower gap to the complete right-copy endpoint."""

    degree = len(gap_coefficients) - 1
    complement = null_space(right.conj().T)
    interior = [
        complement.conj().T @ coefficient @ complement
        for coefficient in gap_coefficients
    ]
    cross = [
        right.conj().T @ coefficient @ complement
        for coefficient in gap_coefficients
    ]
    endpoint = [
        right.conj().T @ coefficient @ right
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


def audit_case(
    construction_kind: str,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    colligation_error: float,
    maximum_grade: int,
) -> list[LowerMetricFlagRecord]:
    """Audit every active member of one ordered right-transfer flag."""

    multiplicity = right.shape[1]
    degree = maximum_grade * (maximum_grade + 1) // 2 + 2
    gap = balanced_lower_gap_coefficients(
        partial,
        right,
        left,
        degree,
    )
    endpoint, constant_error = full_lower_endpoint_series(gap, right)
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
                        transfers[index] @ coordinates
                    )
                )
                for index in range(1, grade)
            ),
            default=0.0,
        )
        active_transfer = transfers[grade] @ coordinates
        if np.linalg.norm(active_transfer) < tolerance:
            continue

        (
            reduced,
            active,
            kernel,
            valuation,
            leading,
        ) = schur_reduce_leading_kernel(current, tolerance)
        expected = active_transfer.conj().T @ active_transfer
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
                "the lower boundary-metric flag audit failed: "
                f"kind={construction_kind}, grade={grade}, "
                f"valuation={valuation}, expected={grade}, "
                f"face={face_error:.3e}, earlier={earlier_error:.3e}"
            )
        records.append(
            LowerMetricFlagRecord(
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
            f"the right transfer flag did not terminate for {construction_kind}"
        )
    return records


def standard_records() -> list[LowerMetricFlagRecord]:
    """Return noncommuting chains and rank-jumping shift flags."""

    records = []
    for multiplicity in range(2, 6):
        case = rank_chain_case(multiplicity, 109_200 + multiplicity)
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
    records: Sequence[LowerMetricFlagRecord],
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
            "experiments/repeated_crabb_lower_metric_flag_s70224.jsonl"
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
