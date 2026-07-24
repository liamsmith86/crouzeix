#!/usr/bin/env python3
"""Audit the all-length repeated-Crabb analytic metric chart.

The metric is made exactly lower-tight by a Schur complement.  Its
range block is then selected by prescribing the Stein Schur
complement.  At a repeated Crabb block the derivative in that range
block is a triangular weighted diagonal recurrence, hence invertible
in every length and copy multiplicity.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.optimize import least_squares

from crabb_block_hardy_equality import (
    format_float,
    hermitian_basis,
    pack_hermitian,
)


@dataclass(frozen=True)
class ChartLinearizationRecord:
    """One all-block Jacobian and inverse-recurrence audit."""

    record_kind: str
    dimension: int
    length: int
    multiplicity: int
    hermitian_variable_count: int
    jacobian_rank: int
    inverse_recurrence_error: str
    all_checks_passed: bool


@dataclass(frozen=True)
class NonlinearChartRecord:
    """One nonlinear prescribed-Stein-slack chart audit."""

    record_kind: str
    dimension: int
    length: int
    multiplicity: int
    solver_success: bool
    lower_schur_error: str
    prescribed_slack_error: str
    minimum_lower_metric_eigenvalue: str
    minimum_stein_eigenvalue: str
    upper_endpoint_maximum_eigenvalue: str
    upper_feasibility_equivalent: bool
    all_checks_passed: bool


AuditRecord = ChartLinearizationRecord | NonlinearChartRecord


def crabb_weights(length: int) -> np.ndarray:
    """Return the canonical Crabb link weights."""

    if length < 2:
        raise ValueError("Crabb length must be at least two")
    weights = np.ones(length)
    weights[0] = np.sqrt(2)
    weights[-1] = np.sqrt(2)
    return weights


def repeated_crabb(length: int, multiplicity: int) -> np.ndarray:
    """Return the repeated canonical Crabb weighted shift."""

    dimension = (length + 1) * multiplicity
    operator = np.zeros((dimension, dimension), dtype=complex)
    identity = np.eye(multiplicity)
    for level, weight in enumerate(crabb_weights(length)):
        operator[
            level * multiplicity : (level + 1) * multiplicity,
            (level + 1) * multiplicity : (level + 2) * multiplicity,
        ] = weight * identity
    return operator


def base_metric(length: int, multiplicity: int) -> np.ndarray:
    """Return ``diag(1,2,...,2,4)`` tensored with copy identity."""

    level_values = [1, *([2] * (length - 1)), 4]
    return np.kron(np.diag(level_values), np.eye(multiplicity))


def unpack_hermitian(values: np.ndarray, size: int) -> np.ndarray:
    """Unpack real coordinates into a complex Hermitian matrix."""

    matrix = np.zeros((size, size), dtype=complex)
    cursor = 0
    for index in range(size):
        matrix[index, index] = values[cursor]
        cursor += 1
    for row in range(size):
        for column in range(row + 1, size):
            value = values[cursor] + 1j * values[cursor + 1]
            matrix[row, column] = value
            matrix[column, row] = np.conjugate(value)
            cursor += 2
    return matrix


def lower_tight_metric(
    range_metric: np.ndarray,
    cross_block: np.ndarray,
    multiplicity: int,
) -> np.ndarray:
    """Construct a metric whose lower Schur complement is exactly zero."""

    range_size = len(range_metric)
    shifted = range_metric - np.eye(range_size)
    corner = (
        np.eye(multiplicity)
        + cross_block @ np.linalg.solve(shifted, cross_block.conj().T)
    )
    metric = np.zeros(
        (multiplicity + range_size,) * 2,
        dtype=complex,
    )
    metric[:multiplicity, :multiplicity] = corner
    metric[:multiplicity, multiplicity:] = cross_block
    metric[multiplicity:, :multiplicity] = cross_block.conj().T
    metric[multiplicity:, multiplicity:] = range_metric
    return metric


def lower_schur_complement(
    metric: np.ndarray,
    multiplicity: int,
) -> np.ndarray:
    """Return the level-zero Schur complement of ``metric-I``."""

    shifted = metric - np.eye(len(metric))
    corner = shifted[:multiplicity, :multiplicity]
    cross = shifted[:multiplicity, multiplicity:]
    range_block = shifted[multiplicity:, multiplicity:]
    return corner - cross @ np.linalg.solve(
        range_block,
        cross.conj().T,
    )


def stein_schur_complement(
    operator: np.ndarray,
    metric: np.ndarray,
    multiplicity: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the Stein defect and its level-zero Schur complement."""

    defect = metric - operator.conj().T @ metric @ operator
    corner = defect[:multiplicity, :multiplicity]
    cross = defect[:multiplicity, multiplicity:]
    kernel = defect[multiplicity:, multiplicity:]
    schur = kernel - cross.conj().T @ np.linalg.solve(corner, cross)
    return defect, schur


def upper_endpoint(
    metric: np.ndarray,
    multiplicity: int,
) -> np.ndarray:
    """Return the final-level Schur endpoint for ``metric<=4I``."""

    split = len(metric) - multiplicity
    leading = metric[:split, :split]
    cross = metric[:split, split:]
    final = metric[split:, split:]
    return (
        final
        - 4 * np.eye(multiplicity)
        + cross.conj().T
        @ np.linalg.solve(
            4 * np.eye(split) - leading,
            cross,
        )
    )


def stein_range_linearization(
    direction: np.ndarray,
    length: int,
    multiplicity: int,
) -> np.ndarray:
    """Apply the triangular weighted diagonal recurrence."""

    weights = crabb_weights(length)
    result = direction.copy()
    for row in range(1, length):
        for column in range(1, length):
            result[
                row * multiplicity : (row + 1) * multiplicity,
                column * multiplicity : (column + 1) * multiplicity,
            ] -= (
                weights[row]
                * weights[column]
                * direction[
                    (row - 1)
                    * multiplicity : row
                    * multiplicity,
                    (column - 1)
                    * multiplicity : column
                    * multiplicity,
                ]
            )
    return result


def invert_stein_range_linearization(
    target: np.ndarray,
    length: int,
    multiplicity: int,
) -> np.ndarray:
    """Invert the weighted diagonal recurrence by forward propagation."""

    weights = crabb_weights(length)
    direction = target.copy()
    for row in range(1, length):
        for column in range(1, length):
            direction[
                row * multiplicity : (row + 1) * multiplicity,
                column * multiplicity : (column + 1) * multiplicity,
            ] += (
                weights[row]
                * weights[column]
                * direction[
                    (row - 1)
                    * multiplicity : row
                    * multiplicity,
                    (column - 1)
                    * multiplicity : column
                    * multiplicity,
                ]
            )
    return direction


def linearization_record(
    length: int,
    multiplicity: int,
) -> ChartLinearizationRecord:
    """Audit full rank and the explicit inverse on every Hermitian basis."""

    size = length * multiplicity
    basis = hermitian_basis(size)
    columns = [
        pack_hermitian(
            stein_range_linearization(
                element,
                length,
                multiplicity,
            )
        )
        for element in basis
    ]
    jacobian = np.stack(columns, axis=1)
    rank = int(np.linalg.matrix_rank(jacobian, tol=1e-10))
    worst_inverse_error = 0.0
    for element in basis:
        recovered = invert_stein_range_linearization(
            stein_range_linearization(
                element,
                length,
                multiplicity,
            ),
            length,
            multiplicity,
        )
        worst_inverse_error = max(
            worst_inverse_error,
            float(np.linalg.norm(recovered - element)),
        )
    expected = size**2
    verified = bool(rank == expected and worst_inverse_error < 1e-12)
    if not verified:
        raise RuntimeError("the repeated metric-chart Jacobian audit failed")
    return ChartLinearizationRecord(
        record_kind="repeated_metric_chart_linearization",
        dimension=(length + 1) * multiplicity,
        length=length,
        multiplicity=multiplicity,
        hermitian_variable_count=expected,
        jacobian_rank=rank,
        inverse_recurrence_error=format_float(worst_inverse_error),
        all_checks_passed=verified,
    )


def deterministic_complex(
    rows: int,
    columns: int,
    scale: float,
) -> np.ndarray:
    """Return a deterministic dense complex perturbation."""

    matrix = np.zeros((rows, columns), dtype=complex)
    for row in range(rows):
        for column in range(columns):
            matrix[row, column] = scale * (
                (1 + 2 * row + 3 * column) / (rows + columns)
                + 1j
                * (2 - row + 2 * column)
                / (2 * (rows + columns))
            )
    return matrix


def nonlinear_record(
    length: int,
    multiplicity: int,
) -> NonlinearChartRecord:
    """Solve one nearby prescribed-slack chart point."""

    dimension = (length + 1) * multiplicity
    range_size = length * multiplicity
    operator = repeated_crabb(length, multiplicity) + deterministic_complex(
        dimension,
        dimension,
        2e-5,
    )
    cross = deterministic_complex(
        multiplicity,
        range_size,
        3e-5,
    )
    slack_factor = deterministic_complex(
        range_size,
        multiplicity,
        3e-4,
    )
    prescribed_slack = slack_factor @ slack_factor.conj().T
    base_range = base_metric(length, multiplicity)[
        multiplicity:,
        multiplicity:,
    ]

    def equations(values: np.ndarray) -> np.ndarray:
        range_metric = unpack_hermitian(values, range_size)
        metric = lower_tight_metric(
            range_metric,
            cross,
            multiplicity,
        )
        _, schur = stein_schur_complement(
            operator,
            metric,
            multiplicity,
        )
        return pack_hermitian(schur - prescribed_slack)

    solution = least_squares(
        equations,
        pack_hermitian(base_range),
        xtol=1e-13,
        ftol=1e-13,
        gtol=1e-13,
        max_nfev=1_000,
    )
    range_metric = unpack_hermitian(solution.x, range_size)
    metric = lower_tight_metric(
        range_metric,
        cross,
        multiplicity,
    )
    defect, recovered_slack = stein_schur_complement(
        operator,
        metric,
        multiplicity,
    )
    lower_error = float(
        np.linalg.norm(lower_schur_complement(metric, multiplicity))
    )
    slack_error = float(np.linalg.norm(recovered_slack - prescribed_slack))
    lower_minimum = float(np.linalg.eigvalsh(metric - np.eye(dimension))[0])
    stein_minimum = float(np.linalg.eigvalsh(defect)[0])
    endpoint = upper_endpoint(metric, multiplicity)
    endpoint_maximum = float(np.linalg.eigvalsh(endpoint)[-1])
    metric_upper = bool(
        np.linalg.eigvalsh(4 * np.eye(dimension) - metric)[0] >= -1e-9
    )
    endpoint_upper = bool(endpoint_maximum <= 1e-9)
    upper_equivalent = metric_upper == endpoint_upper
    verified = bool(
        solution.success
        and lower_error < 1e-10
        and slack_error < 1e-10
        and lower_minimum > -1e-9
        and stein_minimum > -1e-9
        and upper_equivalent
    )
    if not verified:
        raise RuntimeError("the nonlinear repeated metric-chart audit failed")
    return NonlinearChartRecord(
        record_kind="repeated_metric_chart_nonlinear",
        dimension=dimension,
        length=length,
        multiplicity=multiplicity,
        solver_success=bool(solution.success),
        lower_schur_error=format_float(lower_error),
        prescribed_slack_error=format_float(slack_error),
        minimum_lower_metric_eigenvalue=format_float(lower_minimum),
        minimum_stein_eigenvalue=format_float(stein_minimum),
        upper_endpoint_maximum_eigenvalue=format_float(endpoint_maximum),
        upper_feasibility_equivalent=upper_equivalent,
        all_checks_passed=verified,
    )


def write_records(path: Path, records: Sequence[AuditRecord]) -> None:
    """Write deterministic JSON Lines output atomically."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(asdict(record), sort_keys=True))
            handle.write("\n")
    temporary.replace(path)


def main() -> None:
    """Run all-length repeated metric-chart audits."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--minimum-length", type=int, default=2)
    parser.add_argument("--maximum-length", type=int, default=5)
    parser.add_argument(
        "--multiplicities",
        type=int,
        nargs="+",
        default=(1, 2),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/repeated_crabb_exact_metric_chart_s70224.jsonl"
        ),
    )
    args = parser.parse_args()

    records: list[AuditRecord] = []
    for multiplicity in args.multiplicities:
        for length in range(args.minimum_length, args.maximum_length + 1):
            linear = linearization_record(length, multiplicity)
            records.append(linear)
            print(
                "verified metric-chart linearization "
                f"length {length}, multiplicity {multiplicity}",
                flush=True,
            )
            nonlinear = nonlinear_record(length, multiplicity)
            records.append(nonlinear)
            print(
                "verified nonlinear metric chart "
                f"length {length}, multiplicity {multiplicity}",
                flush=True,
            )
    write_records(args.output, records)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
