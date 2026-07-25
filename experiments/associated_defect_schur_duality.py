#!/usr/bin/env python3
"""Audit associated-grade duality of initial/final defect Schur faces."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from math import factorial
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_boundary_slack_deflation import adjoint_series
from repeated_crabb_delayed_jet import inverse_series, series_multiply
from repeated_crabb_transfer_channel_covariance import (
    random_partial_isometry,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class DefectSchurDualityRecord:
    """One deterministic associated-grade duality audit."""

    first_active_degree: int
    state_dimension: int
    defect_dimension: int
    maximum_earlier_initial_residual_norm: str
    maximum_earlier_final_residual_norm: str
    normalized_face_conjugacy_error: str
    metric_face_conjugacy_error: str
    metric_trace_error: str
    all_checks_passed: bool


def zero_series(
    dimension: int,
    maximum_degree: int,
) -> list[Matrix]:
    """Return a zero square-matrix series."""

    return [
        np.zeros((dimension, dimension), dtype=complex)
        for _ in range(maximum_degree + 1)
    ]


def series_add(
    left: list[Matrix],
    right: list[Matrix],
) -> list[Matrix]:
    """Add equal-length matrix series."""

    return [
        left_coefficient + right_coefficient
        for left_coefficient, right_coefficient in zip(
            left,
            right,
            strict=True,
        )
    ]


def exponential_series(
    generator: Matrix,
    maximum_degree: int,
) -> list[Matrix]:
    """Return ``exp(c*generator)`` through one degree."""

    identity = np.eye(len(generator), dtype=complex)
    return [
        np.linalg.matrix_power(generator, degree) / factorial(degree)
        if degree
        else identity
        for degree in range(maximum_degree + 1)
    ]


def schur_residual(
    series: list[Matrix],
    frame: Matrix,
) -> list[Matrix]:
    """Return the fixed-frame Schur residual series."""

    maximum_degree = len(series) - 1
    identity = np.eye(len(series[0]), dtype=complex)
    projection = frame @ frame.conj().T
    complement = identity - projection
    pivot = [
        frame.conj().T @ coefficient @ frame
        for coefficient in series
    ]
    cross = [
        complement @ coefficient @ frame
        for coefficient in series
    ]
    corner = [
        complement @ coefficient @ complement
        for coefficient in series
    ]
    square = series_multiply(
        series_multiply(cross, inverse_series(pivot), maximum_degree),
        adjoint_series(cross),
        maximum_degree,
    )
    return [
        corner_coefficient - square_coefficient
        for corner_coefficient, square_coefficient in zip(
            corner,
            square,
            strict=True,
        )
    ]


def defect_series(
    operator: list[Matrix],
) -> tuple[list[Matrix], list[Matrix]]:
    """Return initial and final defect series."""

    maximum_degree = len(operator) - 1
    dimension = len(operator[0])
    identity_series = zero_series(dimension, maximum_degree)
    identity_series[0] = np.eye(dimension, dtype=complex)
    initial_product = series_multiply(
        adjoint_series(operator),
        operator,
        maximum_degree,
    )
    final_product = series_multiply(
        operator,
        adjoint_series(operator),
        maximum_degree,
    )
    return (
        [
            identity - product
            for identity, product in zip(
                identity_series,
                initial_product,
                strict=True,
            )
        ],
        [
            identity - product
            for identity, product in zip(
                identity_series,
                final_product,
                strict=True,
            )
        ],
    )


def audit_degree(
    first_active_degree: int,
    seed: int,
) -> DefectSchurDualityRecord:
    """Audit one prescribed first residual degree."""

    dimension = 8
    multiplicity = 2
    generator = np.random.default_rng(seed)
    partial, right, left = random_partial_isometry(
        dimension,
        multiplicity,
        generator,
    )

    random_left = (
        generator.standard_normal((dimension, dimension))
        + 1j * generator.standard_normal((dimension, dimension))
    )
    random_right = (
        generator.standard_normal((dimension, dimension))
        + 1j * generator.standard_normal((dimension, dimension))
    )
    left_generator = 0.08 * (random_left - random_left.conj().T)
    right_generator = 0.08 * (random_right - random_right.conj().T)
    left_unitary = exponential_series(
        left_generator,
        first_active_degree,
    )
    right_unitary = exponential_series(
        right_generator,
        first_active_degree,
    )
    base = series_multiply(
        series_multiply(
            left_unitary,
            [partial, *(
                np.zeros_like(partial)
                for _ in range(first_active_degree)
            )],
            first_active_degree,
        ),
        right_unitary,
        first_active_degree,
    )
    perturbation = (
        generator.standard_normal((dimension, dimension))
        + 1j * generator.standard_normal((dimension, dimension))
    ) / (7 * dimension)
    base[first_active_degree] += perturbation

    initial, final = defect_series(base)
    initial_residual = schur_residual(initial, right)
    final_residual = schur_residual(final, left)
    normalized_error = float(
        np.linalg.norm(
            final_residual[first_active_degree]
            - partial
            @ initial_residual[first_active_degree]
            @ partial.conj().T
        )
    )

    random_metric = (
        generator.standard_normal((dimension, dimension))
        + 1j * generator.standard_normal((dimension, dimension))
    )
    metric_generator = 0.05 * (
        random_metric + random_metric.conj().T
    )
    metric_root = exponential_series(
        metric_generator,
        first_active_degree,
    )
    inverse_root = exponential_series(
        -metric_generator,
        first_active_degree,
    )
    physical_operator = series_multiply(
        series_multiply(inverse_root, base, first_active_degree),
        metric_root,
        first_active_degree,
    )
    metric = series_multiply(
        metric_root,
        metric_root,
        first_active_degree,
    )
    inverse_metric = series_multiply(
        inverse_root,
        inverse_root,
        first_active_degree,
    )
    physical_initial = series_add(
        metric,
        [
            -coefficient
            for coefficient in series_multiply(
                series_multiply(
                    adjoint_series(physical_operator),
                    metric,
                    first_active_degree,
                ),
                physical_operator,
                first_active_degree,
            )
        ],
    )
    physical_final = series_add(
        inverse_metric,
        [
            -coefficient
            for coefficient in series_multiply(
                series_multiply(
                    physical_operator,
                    inverse_metric,
                    first_active_degree,
                ),
                adjoint_series(physical_operator),
                first_active_degree,
            )
        ],
    )
    metric_initial_residual = schur_residual(
        physical_initial,
        right,
    )
    metric_final_residual = schur_residual(
        physical_final,
        left,
    )
    metric_error = float(
        np.linalg.norm(
            metric_final_residual[first_active_degree]
            - partial
            @ metric_initial_residual[first_active_degree]
            @ partial.conj().T
        )
    )
    trace_error = abs(
        np.trace(metric_final_residual[first_active_degree])
        - np.trace(metric_initial_residual[first_active_degree])
    )
    earlier_initial = max(
        (
            float(np.linalg.norm(coefficient))
            for coefficient in metric_initial_residual[
                :first_active_degree
            ]
        ),
        default=0.0,
    )
    earlier_final = max(
        (
            float(np.linalg.norm(coefficient))
            for coefficient in metric_final_residual[
                :first_active_degree
            ]
        ),
        default=0.0,
    )

    tolerance = 3e-10
    verified = max(
        earlier_initial,
        earlier_final,
        normalized_error,
        metric_error,
        trace_error,
    ) < tolerance
    if not verified:
        raise RuntimeError(
            "associated defect-Schur duality audit failed: "
            f"degree={first_active_degree}, "
            f"normalized={normalized_error:.3e}, "
            f"metric={metric_error:.3e}"
        )
    return DefectSchurDualityRecord(
        first_active_degree=first_active_degree,
        state_dimension=dimension,
        defect_dimension=multiplicity,
        maximum_earlier_initial_residual_norm=format_float(
            earlier_initial
        ),
        maximum_earlier_final_residual_norm=format_float(
            earlier_final
        ),
        normalized_face_conjugacy_error=format_float(
            normalized_error
        ),
        metric_face_conjugacy_error=format_float(metric_error),
        metric_trace_error=format_float(float(trace_error)),
        all_checks_passed=verified,
    )


def standard_records(
    maximum_degree: int,
) -> list[DefectSchurDualityRecord]:
    """Return deterministic audits through one first-active degree."""

    return [
        audit_degree(degree, 247_000 + degree)
        for degree in range(1, maximum_degree + 1)
    ]


def write_records(
    records: list[DefectSchurDualityRecord],
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
            "associated_defect_schur_duality_s70224.jsonl"
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
