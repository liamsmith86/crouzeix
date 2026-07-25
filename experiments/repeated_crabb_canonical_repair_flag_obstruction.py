#!/usr/bin/env python3
"""Falsify the canonical boundary repair on partial copy-space flags.

L227's canonical correction makes the elliptic Stein slack positive.
Its grade-one upper metric gap is the favorable matrix

    12 c^2 B_1 B_1*.

When ``B_1`` is rank deficient, however, the next coefficient on its
left kernel can be an indefinite cubic matrix.  This checker computes
that coefficient from the formal canonical-repair series.  The scaled
Schur rank chains approach the repeated monomial colligation as the
parameter scale tends to zero, so the obstruction is local rather
than a far-away artifact.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import sqrtm

from crabb_block_hardy_equality import format_float
from repeated_crabb_boundary_metric_flag import (
    full_endpoint_gap_series,
    rank_chain_case,
)
from repeated_crabb_boundary_slack_deflation import (
    slack_schur_coefficients,
)
from repeated_crabb_delayed_jet import ellipse_operator_coefficients
from repeated_crabb_elliptic_cokernel import stein_inverse
from repeated_crabb_one_image_generator import (
    boundary_metric_coefficient,
)
from repeated_crabb_transfer_channel_covariance import (
    transfer_coefficient,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class CanonicalRepairFlagObstructionRecord:
    """One cubic partial-flag obstruction."""

    multiplicity: int
    parameter_scale: str
    state_dimension: int
    colligation_error: str
    spectral_radius: str
    first_transfer_rank: int
    first_transfer_norm: str
    second_face_error: str
    cubic_trace_error: str
    cubic_kernel_dimension: int
    cubic_kernel_minimum_eigenvalue: str
    cubic_kernel_maximum_eigenvalue: str
    maximum_repair_stein_residual: str
    upper_bound_fails_for_small_positive_c: bool
    all_checks_passed: bool


def equality_metric_root(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
) -> tuple[Matrix, Matrix, Matrix]:
    """Return the repeated equality metric and its two square roots."""

    identity = np.eye(len(partial), dtype=complex)
    metric = (
        2 * identity
        - right @ right.conj().T
        + 2 * left @ left.conj().T
    )
    root = np.asarray(sqrtm(metric), dtype=complex)
    return metric, root, np.linalg.inv(root)


def canonical_repair_metric_series(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    degree: int,
) -> tuple[list[Matrix], float]:
    """Return ``4 P^-1 - (P_bl + X)`` through one ``c`` degree."""

    metric, root, inverse_root = equality_metric_root(
        partial,
        right,
        left,
    )
    physical_operator = inverse_root @ partial @ root
    physical_series = ellipse_operator_coefficients(
        physical_operator,
        degree,
    )
    balanced_operator = [
        root @ coefficient @ inverse_root
        for coefficient in physical_series
    ]
    residual = slack_schur_coefficients(
        partial,
        right,
        left,
        degree,
    )

    repair: list[Matrix] = []
    maximum_stein_residual = 0.0
    for order in range(degree + 1):
        forcing = -residual[order].copy()
        for left_degree in range(order + 1):
            for metric_degree in range(order):
                right_degree = order - left_degree - metric_degree
                if 0 <= right_degree <= degree:
                    forcing += (
                        balanced_operator[left_degree].conj().T
                        @ repair[metric_degree]
                        @ balanced_operator[right_degree]
                    )
        coefficient = stein_inverse(partial, forcing)
        repair.append(coefficient)
        stein_residual = (
            coefficient
            - partial.conj().T @ coefficient @ partial
            - forcing
        )
        maximum_stein_residual = max(
            maximum_stein_residual,
            float(np.linalg.norm(stein_residual)),
        )

    identity = np.eye(len(partial), dtype=complex)
    boundary = [identity] + [
        boundary_metric_coefficient(
            partial,
            right,
            left,
            order,
        )
        for order in range(1, degree + 1)
    ]
    upper_gap = []
    for order in range(degree + 1):
        base = 4 * np.linalg.inv(metric) if order == 0 else 0
        upper_gap.append(base - boundary[order] - repair[order])
    return upper_gap, maximum_stein_residual


def leading_left_kernel(matrix: Matrix) -> Matrix:
    """Return a stable orthonormal basis for a PSD matrix's kernel."""

    eigenvalues, eigenvectors = np.linalg.eigh(
        (matrix + matrix.conj().T) / 2
    )
    tolerance = max(1e-11, 1e-8 * float(eigenvalues[-1]))
    return eigenvectors[:, eigenvalues < tolerance]


def audit_case(
    multiplicity: int,
    parameter_scale: float,
    seed: int,
) -> CanonicalRepairFlagObstructionRecord:
    """Compute one canonical-repair cubic obstruction."""

    partial, right, left, colligation_error, _ = rank_chain_case(
        multiplicity,
        seed,
        parameter_scale,
    )
    upper_gap, stein_residual = canonical_repair_metric_series(
        partial,
        right,
        left,
        3,
    )
    endpoint, constant_error = full_endpoint_gap_series(
        upper_gap,
        left,
    )
    physical_endpoint = [4 * coefficient for coefficient in endpoint]

    first_transfer = transfer_coefficient(
        partial,
        right,
        left,
        1,
    )
    first_gram = first_transfer @ first_transfer.conj().T
    expected_second_face = 12 * first_gram
    second_face_error = float(
        np.linalg.norm(physical_endpoint[2] - expected_second_face)
    )

    kernel = leading_left_kernel(expected_second_face)
    cubic = (
        kernel.conj().T @ physical_endpoint[3] @ kernel
    )
    cubic = (cubic + cubic.conj().T) / 2
    cubic_eigenvalues = np.linalg.eigvalsh(cubic)
    cubic_minimum = float(cubic_eigenvalues[0])
    cubic_maximum = float(cubic_eigenvalues[-1])
    cubic_trace_error = abs(float(np.trace(physical_endpoint[3]).real))

    tolerance = 3e-8
    failure_detected = bool(cubic_minimum < -1e-9)
    verified = bool(
        constant_error < tolerance
        and colligation_error < tolerance
        and np.max(np.abs(np.linalg.eigvals(partial))) < 1
        and np.linalg.matrix_rank(first_transfer, 1e-8) == 1
        and kernel.shape[1] == multiplicity - 1
        and second_face_error < tolerance
        and cubic_trace_error < tolerance
        and stein_residual < tolerance
        and failure_detected
    )
    if not verified:
        raise RuntimeError(
            "the canonical-repair flag obstruction audit failed: "
            f"multiplicity={multiplicity}, scale={parameter_scale}, "
            f"second={second_face_error:.3e}, "
            f"cubic_min={cubic_minimum:.3e}, "
            f"Stein={stein_residual:.3e}"
        )

    return CanonicalRepairFlagObstructionRecord(
        multiplicity=multiplicity,
        parameter_scale=format_float(parameter_scale),
        state_dimension=len(partial),
        colligation_error=format_float(colligation_error),
        spectral_radius=format_float(
            float(np.max(np.abs(np.linalg.eigvals(partial))))
        ),
        first_transfer_rank=int(
            np.linalg.matrix_rank(first_transfer, 1e-8)
        ),
        first_transfer_norm=format_float(
            float(np.linalg.norm(first_transfer))
        ),
        second_face_error=format_float(second_face_error),
        cubic_trace_error=format_float(cubic_trace_error),
        cubic_kernel_dimension=kernel.shape[1],
        cubic_kernel_minimum_eigenvalue=format_float(cubic_minimum),
        cubic_kernel_maximum_eigenvalue=format_float(cubic_maximum),
        maximum_repair_stein_residual=format_float(stein_residual),
        upper_bound_fails_for_small_positive_c=failure_detected,
        all_checks_passed=verified,
    )


def standard_records() -> list[CanonicalRepairFlagObstructionRecord]:
    """Return scaled noncommuting rank chains approaching the apex."""

    records = []
    for multiplicity in (4, 5):
        for parameter_scale in (1.0, 0.5, 0.2, 0.1):
            records.append(
                audit_case(
                    multiplicity,
                    parameter_scale,
                    109_200 + multiplicity,
                )
            )
    return records


def write_records(
    records: list[CanonicalRepairFlagObstructionRecord],
    output: Path,
) -> str:
    """Write deterministic JSON Lines atomically and return its hash."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)
    return hashlib.sha256(output.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_canonical_repair_flag_obstruction_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the complete audit."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
