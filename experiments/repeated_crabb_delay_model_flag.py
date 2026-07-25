#!/usr/bin/env python3
"""Audit exact delay covariance of Schur features and state orbits.

If the first ``r`` transfer coefficients vanish, the matrix-inner
transfer factors as ``B(z)=z^r B_r(z)``.  The model kernel then splits
into ``r`` monomial layers and a shifted copy of the deflated kernel.
On the state side, the same statement is the wandering-line
decomposition from L209, including an exact compression formula for
L219's boundary-layer metric.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import null_space

from crabb_block_hardy_equality import format_float
from repeated_crabb_matrix_schur_chart import schur_value
from repeated_crabb_schur_kernel_flag import schur_features
from repeated_crabb_transfer_channel_covariance import inflated_case
from repeated_crabb_transfer_deflation import transfer_coefficient


Matrix = np.ndarray


@dataclass(frozen=True)
class DelayModelFlagRecord:
    """One audit of the kernel, feature, orbit, and metric shifts."""

    delay_length: int
    state_dimension_before: int
    state_dimension_after: int
    defect_dimension: int
    maximum_earlier_transfer_norm: str
    forward_orbit_split_error: str
    retained_boundary_metric_error: str
    transfer_factorization_error: str
    model_kernel_split_error: str
    schur_feature_shift_error: str
    all_checks_passed: bool


def hermitian_part(matrix: Matrix) -> Matrix:
    """Return the Hermitian part of ``matrix``."""

    return (matrix + matrix.conj().T) / 2


def boundary_metric(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    ellipse_parameter: float,
    maximum_terms: int = 2_000,
    tail_tolerance: float = 2e-16,
) -> Matrix:
    """Return L219's balanced boundary-layer metric."""

    identity = np.eye(len(operator), dtype=complex)
    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    nome = ellipse_parameter**2
    result = identity.copy()
    forward = identity.copy()
    backward = identity.copy()
    for index in range(1, maximum_terms + 1):
        forward = forward @ operator
        backward = backward @ operator.conj().T
        result -= (
            nome**index
            / (1 + nome**index)
            * backward
            @ right_projection
            @ forward
        )
        result += (
            nome**index
            * forward
            @ left_projection
            @ backward
        )
        if nome**index < tail_tolerance:
            break
    else:
        raise RuntimeError("the partial-isometry orbit did not decay")
    return hermitian_part(result)


def retained_boundary_metric(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    ellipse_parameter: float,
    delay_length: int,
    maximum_terms: int = 2_000,
    tail_tolerance: float = 2e-16,
) -> Matrix:
    """Return the predicted retained compression after ``r`` delays."""

    identity = np.eye(len(operator), dtype=complex)
    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    nome = ellipse_parameter**2
    result = identity + nome**delay_length * left_projection
    forward = identity.copy()
    backward = identity.copy()
    for index in range(1, maximum_terms + 1):
        forward = forward @ operator
        backward = backward @ operator.conj().T
        result -= (
            nome**index
            / (1 + nome**index)
            * backward
            @ right_projection
            @ forward
        )
        result += (
            nome ** (delay_length + index)
            * forward
            @ left_projection
            @ backward
        )
        if nome**index < tail_tolerance:
            break
    else:
        raise RuntimeError("the deflated orbit did not decay")
    return hermitian_part(result)


def synthetic_delayed_schur_data(
    delay_length: int,
    defect_dimension: int,
    seed: int,
) -> tuple[list[Matrix], Matrix]:
    """Return noncommuting strict Schur data with a forced delay."""

    generator = np.random.default_rng(seed)
    length = delay_length + 3
    parameters: list[Matrix] = []
    for index in range(length):
        if index <= delay_length:
            parameters.append(
                np.zeros(
                    (defect_dimension, defect_dimension),
                    dtype=complex,
                )
            )
            continue
        matrix = (
            generator.standard_normal(
                (defect_dimension, defect_dimension)
            )
            + 1j
            * generator.standard_normal(
                (defect_dimension, defect_dimension)
            )
        )
        parameters.append(0.18 * matrix / np.linalg.norm(matrix, ord=2))

    matrix = (
        generator.standard_normal(
            (defect_dimension, defect_dimension)
        )
        + 1j
        * generator.standard_normal(
            (defect_dimension, defect_dimension)
        )
    )
    terminal, triangular = np.linalg.qr(matrix)
    diagonal = np.diag(triangular)
    phases = np.ones(defect_dimension, dtype=complex)
    nonzero = np.abs(diagonal) > 0
    phases[nonzero] = (
        np.conjugate(diagonal[nonzero])
        / np.abs(diagonal[nonzero])
    )
    return parameters, terminal @ np.diag(phases)


def audit_case(
    delay_length: int,
    state_dimension: int,
    defect_dimension: int,
    seed: int,
) -> DelayModelFlagRecord:
    """Audit one fully delayed unstructured colligation."""

    grade = delay_length + 1
    operator, right, left, _ = inflated_case(
        state_dimension,
        defect_dimension,
        grade,
        defect_dimension,
        seed,
    )
    coefficients = [
        transfer_coefficient(operator, right, left, index)
        for index in range(3 * len(operator) + 1)
    ]
    earlier_norm = max(
        float(np.linalg.norm(coefficients[index]))
        for index in range(1, grade)
    )

    delay_columns = [
        np.linalg.matrix_power(operator, index) @ left
        for index in range(delay_length)
    ]
    removed = np.hstack(delay_columns)
    retained_basis = null_space(removed.conj().T)
    compressed = retained_basis.conj().T @ operator @ retained_basis
    compressed_right = retained_basis.conj().T @ right
    promoted_left = (
        retained_basis.conj().T
        @ np.linalg.matrix_power(operator, delay_length)
        @ left
    )

    maximum_orbit_error = 0.0
    for index in range(2 * grade + 5):
        full_orbit = (
            np.linalg.matrix_power(operator.conj().T, index)
            @ right
        )
        retained_orbit = (
            retained_basis
            @ np.linalg.matrix_power(
                compressed.conj().T,
                index,
            )
            @ compressed_right
        )
        removed_orbit = sum(
            delay_columns[delay_index]
            @ coefficients[index + delay_index]
            for delay_index in range(delay_length)
        )
        maximum_orbit_error = max(
            maximum_orbit_error,
            float(
                np.linalg.norm(
                    full_orbit - retained_orbit - removed_orbit
                )
            ),
        )

    ellipse_parameter = 0.14
    full_boundary = boundary_metric(
        operator,
        right,
        left,
        ellipse_parameter,
    )
    compressed_boundary = (
        retained_basis.conj().T
        @ full_boundary
        @ retained_basis
    )
    expected_boundary = retained_boundary_metric(
        compressed,
        compressed_right,
        promoted_left,
        ellipse_parameter,
        delay_length,
    )
    boundary_error = float(
        np.linalg.norm(compressed_boundary - expected_boundary)
    )

    transfer_factorization_error = 0.0
    kernel_split_error = 0.0
    feature_shift_error = 0.0
    parameters, terminal = synthetic_delayed_schur_data(
        delay_length,
        defect_dimension,
        seed + 500_000,
    )
    tail_parameters = parameters[delay_length:]
    identity = np.eye(defect_dimension, dtype=complex)
    for first_index in range(5):
        first = 0.67 * np.exp(
            2j * np.pi * (first_index + 0.19) / 5
        )
        full_first = schur_value(
            parameters,
            terminal,
            first,
        )
        tail_first = schur_value(
            tail_parameters,
            terminal,
            first,
        )
        transfer_factorization_error = max(
            transfer_factorization_error,
            float(
                np.linalg.norm(
                    full_first - first**delay_length * tail_first
                )
            ),
        )

        full_features = schur_features(
            parameters,
            terminal,
            first,
        )
        tail_features = schur_features(
            tail_parameters,
            terminal,
            first,
        )
        for index in range(delay_length):
            feature_shift_error = max(
                feature_shift_error,
                float(
                    np.linalg.norm(
                        full_features[index] - first**index * identity
                    )
                ),
            )
        for index, feature in enumerate(tail_features):
            feature_shift_error = max(
                feature_shift_error,
                float(
                    np.linalg.norm(
                        full_features[index + delay_length]
                        - first**delay_length * feature
                    )
                ),
            )

        for second_index in range(4):
            second = 0.63 * np.exp(
                2j * np.pi * (second_index + 0.27) / 4
            )
            full_second = schur_value(
                parameters,
                terminal,
                second,
            )
            tail_second = schur_value(
                tail_parameters,
                terminal,
                second,
            )
            direct_kernel = (
                identity - full_first.conj().T @ full_second
            ) / (1 - np.conjugate(first) * second)
            monomial_kernel = sum(
                (np.conjugate(first) * second) ** index * identity
                for index in range(delay_length)
            )
            tail_kernel = (
                identity - tail_first.conj().T @ tail_second
            ) / (1 - np.conjugate(first) * second)
            expected_kernel = (
                monomial_kernel
                + (np.conjugate(first) * second) ** delay_length
                * tail_kernel
            )
            kernel_split_error = max(
                kernel_split_error,
                float(np.linalg.norm(direct_kernel - expected_kernel)),
            )

    tolerance = 3e-8
    verified = bool(
        earlier_norm < tolerance
        and maximum_orbit_error < tolerance
        and boundary_error < tolerance
        and transfer_factorization_error < tolerance
        and kernel_split_error < tolerance
        and feature_shift_error < tolerance
    )
    if not verified:
        raise RuntimeError(
            "delay model flag audit failed: "
            f"earlier={earlier_norm:.3e}, "
            f"orbit={maximum_orbit_error:.3e}, "
            f"boundary={boundary_error:.3e}, "
            f"transfer={transfer_factorization_error:.3e}, "
            f"kernel={kernel_split_error:.3e}, "
            f"feature={feature_shift_error:.3e}"
        )

    return DelayModelFlagRecord(
        delay_length=delay_length,
        state_dimension_before=len(operator),
        state_dimension_after=len(compressed),
        defect_dimension=defect_dimension,
        maximum_earlier_transfer_norm=format_float(earlier_norm),
        forward_orbit_split_error=format_float(maximum_orbit_error),
        retained_boundary_metric_error=format_float(boundary_error),
        transfer_factorization_error=format_float(
            transfer_factorization_error
        ),
        model_kernel_split_error=format_float(kernel_split_error),
        schur_feature_shift_error=format_float(feature_shift_error),
        all_checks_passed=verified,
    )


def standard_records() -> list[DelayModelFlagRecord]:
    """Return deterministic delay lengths one through four."""

    records: list[DelayModelFlagRecord] = []
    for delay_length in range(1, 5):
        for defect_dimension in (1, 2):
            records.append(
                audit_case(
                    delay_length,
                    6 + defect_dimension,
                    defect_dimension,
                    105_000 + 10 * delay_length + defect_dimension,
                )
            )
    return records


def write_records(
    records: list[DelayModelFlagRecord],
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
            "experiments/repeated_crabb_delay_model_flag_s70224.jsonl"
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
