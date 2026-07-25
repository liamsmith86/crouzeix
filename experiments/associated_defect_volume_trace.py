#!/usr/bin/env python3
"""Audit the logarithmic-volume formula for a first defect-Schur trace.

For an analytic path through a partial isometry, the final-defect
Schur residual can contain very large corner and Schur-square
coefficients which cancel at the first active grade.  This checker
verifies that the cancellation is equivalently encoded by one scalar
log-determinant ratio, and by its final-row-whitened form.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from math import factorial
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import null_space

from crabb_block_hardy_equality import format_float
from repeated_crabb_boundary_slack_deflation import adjoint_series
from repeated_crabb_delayed_jet import inverse_series, series_multiply
from repeated_crabb_transfer_channel_covariance import (
    random_partial_isometry,
)


Matrix = np.ndarray
Series = list[Matrix]


@dataclass(frozen=True)
class DefectVolumeTraceRecord:
    """One prescribed first-active-grade volume audit."""

    first_active_degree: int
    state_dimension: int
    defect_dimension: int
    maximum_earlier_residual_norm: str
    schur_log_volume_error: str
    defect_determinant_ratio_error: str
    whitened_channel_ratio_error: str
    all_checks_passed: bool


def zero_series(dimension: int, maximum_degree: int) -> Series:
    """Return a zero square-matrix series."""

    return [
        np.zeros((dimension, dimension), dtype=complex)
        for _ in range(maximum_degree + 1)
    ]


def identity_series(dimension: int, maximum_degree: int) -> Series:
    """Return the constant identity series."""

    result = zero_series(dimension, maximum_degree)
    result[0] = np.eye(dimension, dtype=complex)
    return result


def exponential_series(
    generator: Matrix,
    maximum_degree: int,
) -> Series:
    """Return ``exp(c*generator)`` through one degree."""

    identity = np.eye(len(generator), dtype=complex)
    return [
        np.linalg.matrix_power(generator, degree) / factorial(degree)
        if degree
        else identity
        for degree in range(maximum_degree + 1)
    ]


def series_add(left: Series, right: Series) -> Series:
    """Add equal-length matrix series."""

    return [
        left_coefficient + right_coefficient
        for left_coefficient, right_coefficient in zip(
            left,
            right,
            strict=True,
        )
    ]


def series_left_right(
    left: Matrix,
    series: Series,
    right: Matrix,
) -> Series:
    """Apply constant left and right factors to a matrix series."""

    return [left @ coefficient @ right for coefficient in series]


def trace_log_coefficients(series: Series) -> list[complex]:
    """Return coefficients of ``log det(series)`` above its constant."""

    maximum_degree = len(series) - 1
    inverse_constant = np.linalg.inv(series[0])
    normalized = [
        inverse_constant @ coefficient
        for coefficient in series
    ]
    perturbation = normalized
    perturbation[0] = perturbation[0] - np.eye(len(series[0]))

    result = [0j for _ in range(maximum_degree + 1)]
    power = identity_series(len(series[0]), maximum_degree)
    for exponent in range(1, maximum_degree + 1):
        power = series_multiply(
            power,
            perturbation,
            maximum_degree,
        )
        sign = 1 if exponent % 2 else -1
        for degree in range(1, maximum_degree + 1):
            result[degree] += (
                sign * np.trace(power[degree]) / exponent
            )
    return result


def schur_residual(
    defect: Series,
    frame: Matrix,
) -> tuple[Series, Series]:
    """Return the fixed-frame Schur residual and pivot series."""

    maximum_degree = len(defect) - 1
    identity = np.eye(len(defect[0]), dtype=complex)
    projection = frame @ frame.conj().T
    complement = identity - projection
    pivot = [
        frame.conj().T @ coefficient @ frame
        for coefficient in defect
    ]
    cross = [
        complement @ coefficient @ frame
        for coefficient in defect
    ]
    corner = [
        complement @ coefficient @ complement
        for coefficient in defect
    ]
    square = series_multiply(
        series_multiply(
            cross,
            inverse_series(pivot),
            maximum_degree,
        ),
        adjoint_series(cross),
        maximum_degree,
    )
    return (
        [
            corner_coefficient - square_coefficient
            for corner_coefficient, square_coefficient in zip(
                corner,
                square,
                strict=True,
            )
        ],
        pivot,
    )


def audit_degree(
    first_active_degree: int,
    seed: int,
) -> DefectVolumeTraceRecord:
    """Audit one analytic path with a prescribed first defect grade."""

    dimension = 8
    multiplicity = 2
    generator = np.random.default_rng(seed)
    partial, _, final_frame = random_partial_isometry(
        dimension,
        multiplicity,
        generator,
    )

    left_random = (
        generator.standard_normal((dimension, dimension))
        + 1j * generator.standard_normal((dimension, dimension))
    )
    right_random = (
        generator.standard_normal((dimension, dimension))
        + 1j * generator.standard_normal((dimension, dimension))
    )
    left_generator = 0.06 * (
        left_random - left_random.conj().T
    )
    right_generator = 0.06 * (
        right_random - right_random.conj().T
    )
    left_unitary = exponential_series(
        left_generator,
        first_active_degree,
    )
    right_unitary = exponential_series(
        right_generator,
        first_active_degree,
    )
    constant_partial = [
        partial,
        *(
            np.zeros_like(partial)
            for _ in range(first_active_degree)
        ),
    ]
    operator = series_multiply(
        series_multiply(
            left_unitary,
            constant_partial,
            first_active_degree,
        ),
        right_unitary,
        first_active_degree,
    )
    operator[first_active_degree] += (
        generator.standard_normal((dimension, dimension))
        + 1j * generator.standard_normal((dimension, dimension))
    ) / (9 * dimension)

    final_product = series_multiply(
        operator,
        adjoint_series(operator),
        first_active_degree,
    )
    defect = identity_series(dimension, first_active_degree)
    defect = [
        identity - product
        for identity, product in zip(
            defect,
            final_product,
            strict=True,
        )
    ]
    residual, pivot = schur_residual(defect, final_frame)

    final_projection = final_frame @ final_frame.conj().T
    complement = np.eye(dimension, dtype=complex) - final_projection
    augmented_defect = defect.copy()
    augmented_defect[0] = augmented_defect[0] + complement
    augmented_log = trace_log_coefficients(augmented_defect)
    pivot_log = trace_log_coefficients(pivot)
    schur_log = [
        left - right
        for left, right in zip(
            augmented_log,
            pivot_log,
            strict=True,
        )
    ]

    half_weight = final_projection + 0.5 * complement
    final_gram = series_multiply(
        adjoint_series(operator),
        series_left_right(
            final_projection,
            operator,
            np.eye(dimension, dtype=complex),
        ),
        first_active_degree,
    )
    weighted_gram = series_multiply(
        adjoint_series(operator),
        series_left_right(
            half_weight,
            operator,
            np.eye(dimension, dtype=complex),
        ),
        first_active_degree,
    )
    final_denominator = [
        identity - gram
        for identity, gram in zip(
            identity_series(dimension, first_active_degree),
            final_gram,
            strict=True,
        )
    ]
    weighted_denominator = [
        identity - gram
        for identity, gram in zip(
            identity_series(dimension, first_active_degree),
            weighted_gram,
            strict=True,
        )
    ]
    determinant_ratio = [
        left - right
        for left, right in zip(
            trace_log_coefficients(weighted_denominator),
            trace_log_coefficients(final_denominator),
            strict=True,
        )
    ]

    complement_frame = null_space(final_frame.conj().T)
    final_row = [
        final_frame.conj().T @ coefficient
        for coefficient in operator
    ]
    main_row = [
        complement_frame.conj().T @ coefficient
        for coefficient in operator
    ]
    row_gram = series_multiply(
        adjoint_series(final_row),
        final_row,
        first_active_degree,
    )
    row_denominator = [
        identity - gram
        for identity, gram in zip(
            identity_series(dimension, first_active_degree),
            row_gram,
            strict=True,
        )
    ]
    whitened_gram = series_multiply(
        series_multiply(
            main_row,
            inverse_series(row_denominator),
            first_active_degree,
        ),
        adjoint_series(main_row),
        first_active_degree,
    )
    whitened_denominator = [
        identity - 0.5 * gram
        for identity, gram in zip(
            identity_series(
                dimension - multiplicity,
                first_active_degree,
            ),
            whitened_gram,
            strict=True,
        )
    ]
    whitened_ratio = trace_log_coefficients(
        whitened_denominator
    )

    active_trace = np.trace(residual[first_active_degree])
    schur_error = abs(
        active_trace - schur_log[first_active_degree]
    )
    determinant_error = max(
        abs(left - right)
        for left, right in zip(
            schur_log[1:],
            determinant_ratio[1:],
            strict=True,
        )
    )
    whitened_error = max(
        abs(left - right)
        for left, right in zip(
            determinant_ratio[1:],
            whitened_ratio[1:],
            strict=True,
        )
    )
    earlier = max(
        (
            float(np.linalg.norm(coefficient))
            for coefficient in residual[:first_active_degree]
        ),
        default=0.0,
    )

    tolerance = 2e-9
    verified = max(
        earlier,
        schur_error,
        determinant_error,
        whitened_error,
    ) < tolerance
    if not verified:
        raise RuntimeError(
            "defect volume-trace audit failed: "
            f"degree={first_active_degree}, "
            f"earlier={earlier:.3e}, "
            f"schur={schur_error:.3e}, "
            f"determinant={determinant_error:.3e}, "
            f"whitened={whitened_error:.3e}"
        )
    return DefectVolumeTraceRecord(
        first_active_degree=first_active_degree,
        state_dimension=dimension,
        defect_dimension=multiplicity,
        maximum_earlier_residual_norm=format_float(earlier),
        schur_log_volume_error=format_float(float(schur_error)),
        defect_determinant_ratio_error=format_float(
            float(determinant_error)
        ),
        whitened_channel_ratio_error=format_float(
            float(whitened_error)
        ),
        all_checks_passed=verified,
    )


def standard_records(
    maximum_degree: int,
) -> list[DefectVolumeTraceRecord]:
    """Return deterministic audits through one active degree."""

    return [
        audit_degree(degree, 248_000 + degree)
        for degree in range(1, maximum_degree + 1)
    ]


def write_records(
    records: list[DefectVolumeTraceRecord],
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
    parser.add_argument("--maximum-degree", type=int, default=6)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "associated_defect_volume_trace_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist all audits."""

    args = parse_args()
    records = standard_records(args.maximum_degree)
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
