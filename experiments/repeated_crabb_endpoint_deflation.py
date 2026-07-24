#!/usr/bin/env python3
"""Audit endpoint-response covariance under wandering-delay deflation.

If the first ``r`` transfer coefficients vanish, removing the state
columns ``W, SW, ..., S^(r-1)W`` promotes ``S^r W`` to the left defect
of a smaller partial isometry.  A Stein forcing supported on the
retained state space has exactly the same upper endpoint response in
the original and deflated systems.

The checker also verifies that L212's grade-``r+1`` column is the
literal pullback of L207's grade-one column after this deflation.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import null_space, sqrtm

from crabb_block_hardy_equality import format_float
from repeated_crabb_elliptic_cokernel import endpoint_motion, stein_inverse
from repeated_crabb_grade_two_face import gauge_case
from repeated_crabb_transfer_channel_covariance import (
    inflated_case,
    transfer_channel,
    transfer_coefficient,
)
from repeated_crabb_transfer_deflation import heterogeneous_shift


@dataclass(frozen=True)
class EndpointDeflationRecord:
    """One all-delay endpoint covariance audit."""

    construction_kind: str
    delay_depth: int
    state_dimension_before: int
    state_dimension_after: int
    defect_dimension: int
    spectral_radius: str
    maximum_earlier_transfer_norm: str
    retained_invariance_error: str
    deflated_right_defect_error: str
    deflated_left_defect_error: str
    balanced_stein_compression_error: str
    balanced_endpoint_covariance_error: str
    physical_endpoint_covariance_error: str
    l212_column_pullback_error: str
    l212_endpoint_covariance_error: str
    l212_endpoint_identity_error: str
    all_checks_passed: bool


@dataclass(frozen=True)
class DeflatedSystem:
    """The retained realization after removing a complete delay line."""

    basis: np.ndarray
    partial: np.ndarray
    right: np.ndarray
    left: np.ndarray


def deflate_complete_delay(
    partial: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
    delay_depth: int,
) -> DeflatedSystem:
    """Remove ``W, SW, ..., S^(delay_depth-1)W``."""

    removed = np.hstack(
        [
            np.linalg.matrix_power(partial, degree) @ left
            for degree in range(delay_depth)
        ]
    )
    basis = null_space(removed.conj().T)
    promoted_left = (
        np.linalg.matrix_power(partial, delay_depth) @ left
    )
    return DeflatedSystem(
        basis=basis,
        partial=basis.conj().T @ partial @ basis,
        right=basis.conj().T @ right,
        left=basis.conj().T @ promoted_left,
    )


def physical_data(
    partial: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the canonical physical metric root and operator."""

    identity = np.eye(len(partial), dtype=complex)
    metric = (
        2 * identity
        - right @ right.conj().T
        + 2 * left @ left.conj().T
    )
    root = np.asarray(sqrtm(metric), dtype=complex)
    operator = np.linalg.inv(root) @ partial @ root
    return root, operator


def audit_case(
    construction_kind: str,
    partial: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
    delay_depth: int,
    seed: int,
) -> EndpointDeflationRecord:
    """Audit arbitrary and L212 endpoint motions after deflation."""

    grade = delay_depth + 1
    dimension = len(partial)
    multiplicity = right.shape[1]
    identity = np.eye(dimension, dtype=complex)
    deflated = deflate_complete_delay(
        partial,
        right,
        left,
        delay_depth,
    )
    retained_dimension = len(deflated.partial)
    retained_identity = np.eye(retained_dimension, dtype=complex)

    earlier = [
        transfer_coefficient(partial, right, left, degree)
        for degree in range(1, grade)
    ]
    earlier_norm = max(float(np.linalg.norm(value)) for value in earlier)
    invariance_error = float(
        np.linalg.norm(
            (identity - deflated.basis @ deflated.basis.conj().T)
            @ partial
            @ deflated.basis
        )
    )
    right_defect_error = float(
        np.linalg.norm(
            retained_identity
            - deflated.partial.conj().T @ deflated.partial
            - deflated.right @ deflated.right.conj().T
        )
    )
    left_defect_error = float(
        np.linalg.norm(
            retained_identity
            - deflated.partial @ deflated.partial.conj().T
            - deflated.left @ deflated.left.conj().T
        )
    )

    generator = np.random.default_rng(seed)
    arbitrary = (
        generator.standard_normal((retained_dimension, multiplicity))
        + 1j
        * generator.standard_normal((retained_dimension, multiplicity))
    )
    arbitrary = (
        retained_identity
        - deflated.right @ deflated.right.conj().T
    ) @ arbitrary
    lifted_arbitrary = deflated.basis @ arbitrary
    forcing = (
        right @ lifted_arbitrary.conj().T
        + lifted_arbitrary @ right.conj().T
    )
    deflated_forcing = (
        deflated.right @ arbitrary.conj().T
        + arbitrary @ deflated.right.conj().T
    )
    response = stein_inverse(partial, forcing)
    deflated_response = stein_inverse(
        deflated.partial,
        deflated_forcing,
    )
    compression_error = float(
        np.linalg.norm(
            deflated.basis.conj().T
            @ response
            @ deflated.basis
            - deflated_response
        )
    )
    balanced_endpoint_error = float(
        np.linalg.norm(
            left.conj().T @ response @ left
            - deflated.left.conj().T
            @ deflated_response
            @ deflated.left
        )
    )

    root, operator = physical_data(partial, right, left)
    deflated_root, deflated_operator = physical_data(
        deflated.partial,
        deflated.right,
        deflated.left,
    )
    physical_endpoint_error = float(
        np.linalg.norm(
            endpoint_motion(
                operator,
                right,
                left,
                root @ lifted_arbitrary,
            )
            - endpoint_motion(
                deflated_operator,
                deflated.right,
                deflated.left,
                deflated_root @ arbitrary,
            )
        )
    )

    coefficient = transfer_coefficient(
        partial,
        right,
        left,
        grade,
    )
    deflated_first = transfer_coefficient(
        deflated.partial,
        deflated.right,
        deflated.left,
        1,
    )
    projection = identity - right @ right.conj().T
    correction = (
        -3.5
        * projection
        @ np.linalg.matrix_power(partial, grade)
        @ left
        @ coefficient
    )
    deflated_projection = (
        retained_identity
        - deflated.right @ deflated.right.conj().T
    )
    deflated_correction = (
        -3.5
        * deflated_projection
        @ deflated.partial
        @ deflated.left
        @ deflated_first
    )
    column_pullback_error = float(
        np.linalg.norm(
            deflated.basis.conj().T @ correction
            - deflated_correction
        )
    )
    correction_endpoint = endpoint_motion(
        operator,
        right,
        left,
        root @ correction,
    )
    deflated_correction_endpoint = endpoint_motion(
        deflated_operator,
        deflated.right,
        deflated.left,
        deflated_root @ deflated_correction,
    )
    correction_endpoint_error = float(
        np.linalg.norm(
            correction_endpoint - deflated_correction_endpoint
        )
    )
    expected = 28 * (
        transfer_channel(
            partial,
            right,
            left,
            coefficient.conj().T @ coefficient,
        )
        - coefficient @ coefficient.conj().T
    )
    correction_identity_error = float(
        np.linalg.norm(correction_endpoint - expected)
    )

    spectral_radius = float(
        np.max(np.abs(np.linalg.eigvals(partial)))
    )
    tolerance = 3e-9
    verified = bool(
        spectral_radius < 1
        and earlier_norm < tolerance
        and invariance_error < tolerance
        and right_defect_error < tolerance
        and left_defect_error < tolerance
        and compression_error < tolerance
        and balanced_endpoint_error < tolerance
        and physical_endpoint_error < tolerance
        and column_pullback_error < tolerance
        and correction_endpoint_error < tolerance
        and correction_identity_error < tolerance
    )
    if not verified:
        raise RuntimeError(
            "endpoint deflation covariance failed: "
            f"earlier={earlier_norm:.3e}, "
            f"invariance={invariance_error:.3e}, "
            f"compression={compression_error:.3e}, "
            f"balanced={balanced_endpoint_error:.3e}, "
            f"physical={physical_endpoint_error:.3e}, "
            f"column={column_pullback_error:.3e}, "
            f"endpoint={correction_endpoint_error:.3e}, "
            f"identity={correction_identity_error:.3e}"
        )

    return EndpointDeflationRecord(
        construction_kind=construction_kind,
        delay_depth=delay_depth,
        state_dimension_before=dimension,
        state_dimension_after=retained_dimension,
        defect_dimension=multiplicity,
        spectral_radius=format_float(spectral_radius),
        maximum_earlier_transfer_norm=format_float(earlier_norm),
        retained_invariance_error=format_float(invariance_error),
        deflated_right_defect_error=format_float(right_defect_error),
        deflated_left_defect_error=format_float(left_defect_error),
        balanced_stein_compression_error=format_float(compression_error),
        balanced_endpoint_covariance_error=format_float(
            balanced_endpoint_error
        ),
        physical_endpoint_covariance_error=format_float(
            physical_endpoint_error
        ),
        l212_column_pullback_error=format_float(column_pullback_error),
        l212_endpoint_covariance_error=format_float(
            correction_endpoint_error
        ),
        l212_endpoint_identity_error=format_float(
            correction_identity_error
        ),
        all_checks_passed=verified,
    )


def standard_records() -> list[EndpointDeflationRecord]:
    """Return deterministic unstructured, heterogeneous, and apex cases."""

    records: list[EndpointDeflationRecord] = []
    generator = np.random.default_rng(71624)
    for delay_depth in range(1, 6):
        for multiplicity in (1, 2, 3):
            records.append(
                audit_case(
                    "inflated_unstructured_partial_isometry",
                    *inflated_case(
                        8 + delay_depth + multiplicity,
                        multiplicity,
                        delay_depth + 1,
                        multiplicity,
                        716_000 + 10 * delay_depth + multiplicity,
                    )[:3],
                    delay_depth,
                    816_000 + 10 * delay_depth + multiplicity,
                )
            )

        lengths = (
            delay_depth + 1,
            delay_depth + 2,
            delay_depth + 4,
        )
        records.append(
            audit_case(
                "gauged_heterogeneous_shift",
                *gauge_case(
                    *heterogeneous_shift(lengths),
                    generator,
                ),
                delay_depth,
                916_000 + delay_depth,
            )
        )
        records.append(
            audit_case(
                "repeated_apex",
                *gauge_case(
                    *heterogeneous_shift(
                        (delay_depth + 1,) * 2
                    ),
                    generator,
                ),
                delay_depth,
                926_000 + delay_depth,
            )
        )
    return records


def write_records(
    records: list[EndpointDeflationRecord],
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
            "repeated_crabb_endpoint_deflation_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the endpoint-deflation audit."""

    args = parse_args()
    records = standard_records()
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
