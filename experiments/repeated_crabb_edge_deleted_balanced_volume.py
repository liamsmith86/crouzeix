#!/usr/bin/env python3
"""Audit the balanced-coordinate form of the edge-deleted volume.

On a completely delayed grade ``k``, delete the full degree-``2*k``
coefficient of L219's boundary metric.  Its final-defect row is then
constant and orthogonal to the retained metric through the target
degree.  Consequently L248's normalized initial-defect mass can be
computed without either metric square root:

    -m + tr((R - A* F A)^(-1) (R - A* R A)).

This checker compares that balanced expression with the literal
metric-normalized expression and with the final-defect Schur trace.
The observed value ``4 ||B_k||_F^2`` remains an audit of A194, not a
proof of that open coefficient.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from associated_defect_volume_trace import (
    schur_residual,
    trace_log_coefficients,
)
from crabb_block_hardy_equality import format_float
from repeated_crabb_boundary_slack_deflation import (
    adjoint_series,
    physical_metric_root,
)
from repeated_crabb_delayed_jet import (
    ellipse_operator_coefficients,
    inverse_series,
    series_multiply,
)
from repeated_crabb_one_image_generator import (
    boundary_metric_coefficient,
)
from repeated_crabb_transfer_channel_covariance import (
    inflated_case,
    random_partial_isometry,
    transfer_coefficient,
)


Matrix = np.ndarray
Series = list[Matrix]


@dataclass(frozen=True)
class EdgeDeletedBalancedVolumeRecord:
    """One deterministic balanced-volume audit."""

    grade: int
    state_dimension: int
    defect_dimension: int
    final_metric_block_error: str
    final_metric_cross_error: str
    balanced_normalized_mass_error: str
    balanced_determinant_error: str
    output_determinant_error: str
    balanced_schur_trace_error: str
    maximum_earlier_mass: str
    active_volume_error: str
    all_checks_passed: bool


def zero_series(dimension: int, maximum_degree: int) -> Series:
    """Return a zero square-matrix series."""

    return [
        np.zeros((dimension, dimension), dtype=complex)
        for _ in range(maximum_degree + 1)
    ]


def identity_series(dimension: int, maximum_degree: int) -> Series:
    """Return a constant identity series."""

    result = zero_series(dimension, maximum_degree)
    result[0] = np.eye(dimension, dtype=complex)
    return result


def subtract_series(left: Series, right: Series) -> Series:
    """Subtract equal-length matrix series."""

    return [
        left_coefficient - right_coefficient
        for left_coefficient, right_coefficient in zip(
            left,
            right,
            strict=True,
        )
    ]


def square_root_series(series: Series) -> Series:
    """Return the square root of a series with identity constant term."""

    maximum_degree = len(series) - 1
    dimension = len(series[0])
    result = identity_series(dimension, maximum_degree)
    for degree in range(1, maximum_degree + 1):
        convolution = np.zeros_like(series[0])
        for left_degree in range(1, degree):
            convolution += result[left_degree] @ result[degree - left_degree]
        result[degree] = 0.5 * (series[degree] - convolution)
    return result


def weighted_product(
    operator: Series,
    weight: Series,
    maximum_degree: int,
) -> Series:
    """Return ``operator* weight operator`` as a series."""

    return series_multiply(
        series_multiply(
            adjoint_series(operator),
            weight,
            maximum_degree,
        ),
        operator,
        maximum_degree,
    )


def traces(series: Series) -> list[complex]:
    """Return the coefficientwise traces."""

    return [np.trace(coefficient) for coefficient in series]


def delayed_case(
    grade: int,
) -> tuple[Matrix, Matrix, Matrix]:
    """Return one deterministic completely delayed colligation."""

    multiplicity = 2 if grade <= 4 else 1
    if grade == 1:
        return random_partial_isometry(
            8,
            multiplicity,
            np.random.default_rng(250_001),
        )
    partial, right, left, _ = inflated_case(
        grade + 3,
        multiplicity,
        grade,
        multiplicity,
        250_000 + grade,
    )
    return partial, right, left


def audit_grade(
    grade: int,
) -> EdgeDeletedBalancedVolumeRecord:
    """Audit one edge-deleted delayed grade."""

    partial, right, left = delayed_case(grade)
    dimension = len(partial)
    multiplicity = right.shape[1]
    maximum_degree = 2 * grade
    identity = np.eye(dimension, dtype=complex)
    final = left @ left.conj().T
    complement = identity - final

    physical_root, physical_root_inverse = physical_metric_root(
        partial,
        right,
        left,
    )
    physical_operator = physical_root_inverse @ partial @ physical_root
    physical_series = ellipse_operator_coefficients(
        physical_operator,
        maximum_degree,
    )
    operator = [
        physical_root @ coefficient @ physical_root_inverse
        for coefficient in physical_series
    ]

    metric = [
        identity,
        *(
            boundary_metric_coefficient(
                partial,
                right,
                left,
                degree,
            )
            for degree in range(1, maximum_degree + 1)
        ),
    ]
    metric[maximum_degree] = np.zeros_like(identity)

    final_block_error = max(
        (
            float(np.linalg.norm(left.conj().T @ coefficient @ left))
            for coefficient in metric[1:]
        ),
        default=0.0,
    )
    final_cross_error = max(
        (
            float(np.linalg.norm(complement @ coefficient @ final))
            for coefficient in metric[1:]
        ),
        default=0.0,
    )

    final_weight = zero_series(dimension, maximum_degree)
    final_weight[0] = final
    row_denominator = subtract_series(
        metric,
        weighted_product(
            operator,
            final_weight,
            maximum_degree,
        ),
    )
    stein_slack = subtract_series(
        metric,
        weighted_product(
            operator,
            metric,
            maximum_degree,
        ),
    )
    balanced_mass = series_multiply(
        inverse_series(row_denominator),
        stein_slack,
        maximum_degree,
    )
    balanced_trace = traces(balanced_mass)
    balanced_trace[0] -= multiplicity

    metric_root = square_root_series(metric)
    metric_root_inverse = inverse_series(metric_root)
    normalized_operator = series_multiply(
        series_multiply(
            metric_root,
            operator,
            maximum_degree,
        ),
        metric_root_inverse,
        maximum_degree,
    )
    normalized_initial_defect = subtract_series(
        identity_series(dimension, maximum_degree),
        series_multiply(
            adjoint_series(normalized_operator),
            normalized_operator,
            maximum_degree,
        ),
    )
    normalized_row_gram = weighted_product(
        normalized_operator,
        final_weight,
        maximum_degree,
    )
    normalized_row_denominator = subtract_series(
        identity_series(dimension, maximum_degree),
        normalized_row_gram,
    )
    normalized_mass = series_multiply(
        inverse_series(normalized_row_denominator),
        normalized_initial_defect,
        maximum_degree,
    )
    normalized_trace = traces(normalized_mass)
    normalized_trace[0] -= multiplicity

    balanced_half_weight = [
        0.5 * complement @ coefficient @ complement for coefficient in metric
    ]
    balanced_half_weight[0] += final
    balanced_numerator = subtract_series(
        metric,
        weighted_product(
            operator,
            balanced_half_weight,
            maximum_degree,
        ),
    )
    balanced_log_volume = [
        numerator - denominator
        for numerator, denominator in zip(
            trace_log_coefficients(balanced_numerator),
            trace_log_coefficients(row_denominator),
            strict=True,
        )
    ]
    output_gram = series_multiply(
        series_multiply(
            operator,
            inverse_series(metric),
            maximum_degree,
        ),
        adjoint_series(operator),
        maximum_degree,
    )
    output_numerator = subtract_series(
        identity_series(dimension, maximum_degree),
        series_multiply(
            balanced_half_weight,
            output_gram,
            maximum_degree,
        ),
    )
    output_denominator = subtract_series(
        identity_series(dimension, maximum_degree),
        series_multiply(
            final_weight,
            output_gram,
            maximum_degree,
        ),
    )
    output_log_volume = [
        numerator - denominator
        for numerator, denominator in zip(
            trace_log_coefficients(output_numerator),
            trace_log_coefficients(output_denominator),
            strict=True,
        )
    ]

    normalized_half_weight = zero_series(
        dimension,
        maximum_degree,
    )
    normalized_half_weight[0] = final + 0.5 * complement
    normalized_numerator = subtract_series(
        identity_series(dimension, maximum_degree),
        weighted_product(
            normalized_operator,
            normalized_half_weight,
            maximum_degree,
        ),
    )
    normalized_log_volume = [
        numerator - denominator
        for numerator, denominator in zip(
            trace_log_coefficients(normalized_numerator),
            trace_log_coefficients(normalized_row_denominator),
            strict=True,
        )
    ]

    final_product = series_multiply(
        normalized_operator,
        adjoint_series(normalized_operator),
        maximum_degree,
    )
    final_defect = subtract_series(
        identity_series(dimension, maximum_degree),
        final_product,
    )
    residual, _ = schur_residual(final_defect, left)
    residual_trace = traces(residual)

    balanced_normalized_error = max(
        abs(left_value - right_value)
        for left_value, right_value in zip(
            balanced_trace,
            normalized_trace,
            strict=True,
        )
    )
    balanced_determinant_error = max(
        abs(left_value - right_value)
        for left_value, right_value in zip(
            balanced_log_volume,
            normalized_log_volume,
            strict=True,
        )
    )
    output_determinant_error = max(
        abs(left_value - right_value)
        for left_value, right_value in zip(
            balanced_log_volume,
            output_log_volume,
            strict=True,
        )
    )
    balanced_schur_error = max(
        abs(left_value - right_value)
        for left_value, right_value in zip(
            balanced_trace,
            residual_trace,
            strict=True,
        )
    )
    earlier_mass = max(
        (abs(value) for value in balanced_trace[:maximum_degree]),
        default=0.0,
    )
    active = transfer_coefficient(
        partial,
        right,
        left,
        grade,
    )
    target = 4 * float(np.linalg.norm(active) ** 2)
    active_error = abs(float(balanced_trace[maximum_degree].real) - target)

    tolerance = 8e-7
    verified = (
        max(
            final_block_error,
            final_cross_error,
            balanced_normalized_error,
            balanced_determinant_error,
            output_determinant_error,
            balanced_schur_error,
            earlier_mass,
            active_error,
        )
        < tolerance
    )
    if not verified:
        raise RuntimeError(
            "the edge-deleted balanced-volume audit failed: "
            f"grade={grade}, block={final_block_error:.3e}, "
            f"cross={final_cross_error:.3e}, "
            f"normalized={balanced_normalized_error:.3e}, "
            f"determinant={balanced_determinant_error:.3e}, "
            f"output={output_determinant_error:.3e}, "
            f"schur={balanced_schur_error:.3e}, "
            f"earlier={earlier_mass:.3e}, "
            f"active={active_error:.3e}"
        )
    return EdgeDeletedBalancedVolumeRecord(
        grade=grade,
        state_dimension=dimension,
        defect_dimension=multiplicity,
        final_metric_block_error=format_float(final_block_error),
        final_metric_cross_error=format_float(final_cross_error),
        balanced_normalized_mass_error=format_float(float(balanced_normalized_error)),
        balanced_determinant_error=format_float(float(balanced_determinant_error)),
        output_determinant_error=format_float(float(output_determinant_error)),
        balanced_schur_trace_error=format_float(float(balanced_schur_error)),
        maximum_earlier_mass=format_float(float(earlier_mass)),
        active_volume_error=format_float(float(active_error)),
        all_checks_passed=bool(verified),
    )


def standard_records(
    maximum_grade: int,
) -> list[EdgeDeletedBalancedVolumeRecord]:
    """Return the deterministic delayed audits."""

    return [audit_grade(grade) for grade in range(1, maximum_grade + 1)]


def write_records(
    records: list[EdgeDeletedBalancedVolumeRecord],
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
    parser.add_argument("--maximum-grade", type=int, default=5)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/repeated_crabb_edge_deleted_balanced_volume_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the balanced-volume audits."""

    args = parse_args()
    records = standard_records(args.maximum_grade)
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
