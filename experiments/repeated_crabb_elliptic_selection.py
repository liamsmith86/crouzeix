#!/usr/bin/env python3
"""Audit the explicit analytic selection for the grade-one elliptic face.

For a balanced pure partial isometry ``S`` with orthogonal right and
left defect frames ``V,W``, put

    B_1 = W* S* V,    Q = I - V V*.

The L207 correction is the polynomial column

    C_hat = -(7/2) Q S W B_1.

In physical coordinates it is ``C = P^(1/2) C_hat``.  The checker
verifies that this column solves L204's complete endpoint equation,
both on deterministic inverse-block-Toeplitz equality anchors and on
unstructured partial isometries.  It also checks the finite
noncommutative expansion used in the proof.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import sqrtm

from repeated_crabb_elliptic_cokernel import (
    ReducedFaceData,
    build_reduced_face_data,
    endpoint_motion,
    stein_inverse,
)
from repeated_crabb_elliptic_second_face import (
    canonical_first_forcing_motion,
    second_ellipse_coefficients,
    second_forcing,
)
from crabb_block_hardy_equality import format_float


@dataclass(frozen=True)
class EllipticSelectionRecord:
    """One explicit-selection audit."""

    anchor_kind: str
    state_dimension: int
    defect_dimension: int
    length: int
    sample_index: int
    anchor_scale: str
    spectral_radius: str
    normal_corner_transfer_error: str
    perpendicular_column_error: str
    forcing_expansion_error: str
    range_equation_error: str
    full_second_stein_error: str
    lower_endpoint_error: str
    upper_endpoint_error: str
    correction_norm: str
    all_checks_passed: bool


def matrix_power(matrix: np.ndarray, exponent: int) -> np.ndarray:
    """Return a nonnegative integer matrix power."""

    return np.linalg.matrix_power(matrix, exponent)


def explicit_balanced_column(
    balanced: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return L207's balanced correction and first transfer coefficient."""

    projection = np.eye(len(balanced)) - right @ right.conj().T
    first_transfer = left.conj().T @ balanced.conj().T @ right
    column = (
        -3.5
        * projection
        @ balanced
        @ left
        @ first_transfer
    )
    return column, first_transfer


def expanded_balanced_forcing(balanced: np.ndarray) -> np.ndarray:
    """Return the reduced word expansion of L207's total forcing."""

    operator = balanced
    adjoint = balanced.conj().T
    identity = np.eye(len(operator))
    return (
        -4 * identity
        + 4 * operator @ adjoint
        - 2 * matrix_power(operator, 4)
        + 4 * matrix_power(adjoint, 2) @ matrix_power(operator, 2)
        - 2 * matrix_power(adjoint, 4)
        - 0.5
        * matrix_power(operator, 2)
        @ matrix_power(adjoint, 3)
        @ operator
        - operator
        @ matrix_power(adjoint, 3)
        @ matrix_power(operator, 2)
        + 2 * adjoint @ matrix_power(operator, 5)
        - 0.5
        * adjoint
        @ matrix_power(operator, 3)
        @ matrix_power(adjoint, 2)
        - matrix_power(adjoint, 2)
        @ matrix_power(operator, 3)
        @ adjoint
        - matrix_power(adjoint, 3) @ matrix_power(operator, 3)
        + 2 * matrix_power(adjoint, 5) @ operator
    )


def face_data_from_partial_isometry(
    balanced: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
) -> ReducedFaceData:
    """Construct L204 data from an unstructured partial isometry."""

    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    metric = (
        2 * np.eye(len(balanced))
        - right_projection
        + 2 * left_projection
    )
    metric_root = np.asarray(sqrtm(metric), dtype=complex)
    metric_root_inverse = np.linalg.inv(metric_root)
    operator = metric_root_inverse @ balanced @ metric_root
    first, second = second_ellipse_coefficients(operator)
    first_column = canonical_first_forcing_motion(
        operator,
        metric,
        right,
        first,
    )
    fixed_forcing = second_forcing(
        operator,
        metric,
        first,
        second,
        first_column,
    )
    lower_corner = right.conj().T @ fixed_forcing @ right
    lower_corner = (lower_corner + lower_corner.conj().T) / 2
    parallel_column = -right @ lower_corner / 2
    base_forcing = (
        fixed_forcing
        + right @ parallel_column.conj().T
        + parallel_column @ right.conj().T
    )
    base_metric = stein_inverse(operator, base_forcing)
    balanced_first = metric_root @ first @ metric_root_inverse
    normal_corner = left.conj().T @ balanced_first @ right
    target = -normal_corner @ normal_corner.conj().T
    base_endpoint = left.conj().T @ base_metric @ left
    target_gap = target - base_endpoint
    target_gap = (target_gap + target_gap.conj().T) / 2
    return ReducedFaceData(
        operator=operator,
        metric=metric,
        right_defect=right,
        left_defect=left,
        fixed_forcing=fixed_forcing,
        parallel_column=parallel_column,
        base_metric_direction=base_metric,
        normal_corner=normal_corner,
        target_endpoint=target,
        target_gap=target_gap,
        commutator_norm=0.0,
        actual_strength=0.0,
    )


def audit_data(
    data: ReducedFaceData,
    *,
    anchor_kind: str,
    length: int,
    sample_index: int,
    anchor_scale: float,
) -> EllipticSelectionRecord:
    """Audit the explicit column on one reduced face."""

    operator = data.operator
    metric = data.metric
    right = data.right_defect
    left = data.left_defect
    metric_root = np.asarray(sqrtm(metric), dtype=complex)
    metric_root_inverse = np.linalg.inv(metric_root)
    balanced = metric_root @ operator @ metric_root_inverse
    balanced_column, first_transfer = explicit_balanced_column(
        balanced,
        right,
        left,
    )
    column = metric_root @ balanced_column

    normal_error = float(
        np.linalg.norm(data.normal_corner - 4 * first_transfer)
    )
    perpendicular_error = float(
        np.linalg.norm(right.conj().T @ column)
    )
    range_error = float(
        np.linalg.norm(
            endpoint_motion(
                operator,
                right,
                left,
                column,
            )
            - data.target_gap
        )
    )

    second_column = data.parallel_column + column
    total_forcing = (
        data.fixed_forcing
        + right @ second_column.conj().T
        + second_column @ right.conj().T
    )
    metric_direction = stein_inverse(operator, total_forcing)
    reconstructed = (
        metric_direction
        - operator.conj().T @ metric_direction @ operator
    )
    full_stein_error = float(
        np.linalg.norm(reconstructed - total_forcing)
    )
    lower_error = float(
        np.linalg.norm(
            right.conj().T @ metric_direction @ right
        )
    )
    upper_error = float(
        np.linalg.norm(
            left.conj().T @ metric_direction @ left
            - data.target_endpoint
        )
    )

    balanced_forcing = (
        metric_root_inverse
        @ total_forcing
        @ metric_root_inverse
    )
    expansion_error = float(
        np.linalg.norm(
            balanced_forcing
            - expanded_balanced_forcing(balanced)
        )
    )
    spectral_radius = float(
        np.max(np.abs(np.linalg.eigvals(balanced)))
    )
    correction_norm = float(np.linalg.norm(column))
    verified = bool(
        spectral_radius < 1
        and normal_error < 2e-9
        and perpendicular_error < 2e-9
        and expansion_error < 2e-8
        and range_error < 2e-8
        and full_stein_error < 2e-8
        and lower_error < 2e-8
        and upper_error < 2e-8
    )
    if not verified:
        raise RuntimeError("the explicit elliptic selection audit failed")

    return EllipticSelectionRecord(
        anchor_kind=anchor_kind,
        state_dimension=len(operator),
        defect_dimension=right.shape[1],
        length=length,
        sample_index=sample_index,
        anchor_scale=format_float(anchor_scale),
        spectral_radius=format_float(spectral_radius),
        normal_corner_transfer_error=format_float(normal_error),
        perpendicular_column_error=format_float(perpendicular_error),
        forcing_expansion_error=format_float(expansion_error),
        range_equation_error=format_float(range_error),
        full_second_stein_error=format_float(full_stein_error),
        lower_endpoint_error=format_float(lower_error),
        upper_endpoint_error=format_float(upper_error),
        correction_norm=format_float(correction_norm),
        all_checks_passed=verified,
    )


def haar_unitary(size: int, rng: np.random.Generator) -> np.ndarray:
    """Return a deterministic-seed Haar unitary."""

    matrix = (
        rng.standard_normal((size, size))
        + 1j * rng.standard_normal((size, size))
    )
    unitary, triangular = np.linalg.qr(matrix)
    diagonal = np.diag(triangular)
    phases = np.conjugate(diagonal) / np.abs(diagonal)
    return unitary @ np.diag(phases)


def random_partial_isometry(
    state_dimension: int,
    defect_dimension: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a partial isometry with orthogonal fixed defect frames."""

    right = np.eye(state_dimension, dtype=complex)[:, :defect_dimension]
    left = np.eye(state_dimension, dtype=complex)[:, -defect_dimension:]
    domain = np.eye(state_dimension, dtype=complex)[:, defect_dimension:]
    range_frame = np.eye(state_dimension, dtype=complex)[
        :, : state_dimension - defect_dimension
    ]
    bridge = haar_unitary(
        state_dimension - defect_dimension,
        rng,
    )
    balanced = range_frame @ bridge @ domain.conj().T
    return balanced, right, left


def standard_records() -> list[EllipticSelectionRecord]:
    """Return deterministic equality and unstructured audit records."""

    records = [
        audit_data(
            build_reduced_face_data(length, multiplicity, strength),
            anchor_kind="block_toeplitz_equality",
            length=length,
            sample_index=0,
            anchor_scale=strength,
        )
        for multiplicity in (2, 3)
        for length in range(2, 6)
        for strength in (0.1, 8.0)
    ]

    rng = np.random.default_rng(70224)
    specifications = ((5, 1), (6, 2), (8, 2), (9, 3), (12, 3))
    for state_dimension, defect_dimension in specifications:
        accepted = 0
        while accepted < 4:
            balanced, right, left = random_partial_isometry(
                state_dimension,
                defect_dimension,
                rng,
            )
            spectral_radius = float(
                np.max(np.abs(np.linalg.eigvals(balanced)))
            )
            if spectral_radius >= 0.98:
                continue
            data = face_data_from_partial_isometry(
                balanced,
                right,
                left,
            )
            records.append(
                audit_data(
                    data,
                    anchor_kind="unstructured_partial_isometry",
                    length=-1,
                    sample_index=accepted,
                    anchor_scale=spectral_radius,
                )
            )
            accepted += 1
    return records


def write_records(
    records: list[EllipticSelectionRecord],
    output: Path,
) -> None:
    """Write deterministic JSON Lines records atomically."""

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
            "repeated_crabb_elliptic_selection_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the explicit analytic-selection audit."""

    args = parse_args()
    records = standard_records()
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
