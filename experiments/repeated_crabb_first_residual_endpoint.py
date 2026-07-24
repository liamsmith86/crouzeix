#!/usr/bin/env python3
"""Audit the repeated-Crabb first-residual endpoint Gram square.

The normalized rank-``m`` Hardy metric belongs exactly to the
zero-Stein-slack branch of the all-length metric chart.  If the first
Hardy residual blocks are ``F[r,c]``, the chart's upper endpoint has
leading coefficient

    -4 * sum(F[r,c]^* F[r,c]).

The script checks both the exact zero-slack bridge and this
operator-valued quadratic coefficient.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.linalg import eigh, solve_discrete_lyapunov

from crabb_block_hardy_equality import (
    block_disk_model,
    format_float,
    generic_hermitian,
    linearized_block_residual,
)
from repeated_crabb_exact_metric_chart import (
    lower_schur_complement,
    stein_schur_complement,
    upper_endpoint,
)


@dataclass(frozen=True)
class ZeroSlackBridgeRecord:
    """One exact canonical-Hardy/metric-chart compatibility audit."""

    record_kind: str
    dimension: int
    length: int
    multiplicity: int
    lower_schur_error: str
    stein_schur_error: str
    upper_endpoint_maximum_eigenvalue: str
    all_checks_passed: bool


@dataclass(frozen=True)
class FirstResidualEndpointRecord:
    """One operator-valued first-residual endpoint audit."""

    record_kind: str
    dimension: int
    length: int
    multiplicity: int
    residual_gram_minimum_eigenvalue: str
    residual_gram_norm: str
    symmetric_second_coefficient_error: str
    relative_second_coefficient_error: str
    endpoint_second_coefficient_maximum_eigenvalue: str
    all_checks_passed: bool


AuditRecord = ZeroSlackBridgeRecord | FirstResidualEndpointRecord


def matrix_power(matrix: np.ndarray, exponent: float) -> np.ndarray:
    """Return a Hermitian matrix power."""

    values, vectors = eigh(matrix)
    return (vectors * (values**exponent)) @ vectors.conj().T


def canonical_physical_data(
    hermitian: np.ndarray,
    length: int,
    multiplicity: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the physical disk operator and normalized Hardy metric."""

    extended, coordinate, operator, _ = block_disk_model(
        hermitian,
        length,
        multiplicity,
    )
    dimension = (length + 1) * multiplicity
    endpoint = np.zeros((dimension, multiplicity), dtype=complex)
    endpoint[:multiplicity, :] = np.eye(multiplicity)
    corner = endpoint.conj().T @ extended @ endpoint
    defect = extended @ endpoint @ matrix_power(corner, -0.5)
    metric = solve_discrete_lyapunov(
        operator.conj().T,
        defect @ defect.conj().T,
    )
    coordinate_half = matrix_power(coordinate, 0.5)
    coordinate_inverse_half = matrix_power(coordinate, -0.5)
    physical_operator = (
        coordinate_half @ operator @ coordinate_inverse_half
    )
    physical_metric = (
        coordinate_inverse_half @ metric @ coordinate_inverse_half
    )
    return physical_operator, physical_metric


def actual_residual_gram(
    direction: np.ndarray,
    length: int,
    multiplicity: int,
) -> np.ndarray:
    """Return ``sum F_rc^* F_rc`` in terminal copy coordinates."""

    arranged = linearized_block_residual(
        direction,
        length,
        multiplicity,
    )
    gram = np.zeros((multiplicity, multiplicity), dtype=complex)
    for degree in range(length - 1):
        for interior in range(length - 1):
            arranged_block = arranged[
                degree
                * multiplicity : (degree + 1)
                * multiplicity,
                interior
                * multiplicity : (interior + 1)
                * multiplicity,
            ]
            coefficient = arranged_block.T
            gram += coefficient.conj().T @ coefficient
    return (gram + gram.conj().T) / 2


def deterministic_hermitian_direction(
    length: int,
    multiplicity: int,
) -> np.ndarray:
    """Return one normalized dense Hermitian Gram direction."""

    size = length * multiplicity
    matrix = np.zeros((size, size), dtype=complex)
    for row in range(size):
        matrix[row, row] = (1 + 2 * row) / (3 * size)
        for column in range(row + 1, size):
            value = (
                (2 + row + 3 * column) / (4 * size)
                + 1j * (1 - 2 * row + column) / (5 * size)
            )
            matrix[row, column] = value
            matrix[column, row] = np.conjugate(value)
    return matrix / np.linalg.norm(matrix)


def canonical_endpoint(
    hermitian: np.ndarray,
    length: int,
    multiplicity: int,
) -> np.ndarray:
    """Return the zero-slack chart endpoint of the Hardy metric."""

    _, metric = canonical_physical_data(
        hermitian,
        length,
        multiplicity,
    )
    endpoint = upper_endpoint(metric, multiplicity)
    return (endpoint + endpoint.conj().T) / 2


def bridge_record(
    length: int,
    multiplicity: int,
) -> ZeroSlackBridgeRecord:
    """Audit that the normalized Hardy metric is exactly zero-slack."""

    hermitian = generic_hermitian(length, multiplicity)
    operator, metric = canonical_physical_data(
        hermitian,
        length,
        multiplicity,
    )
    lower_error = float(
        np.linalg.norm(lower_schur_complement(metric, multiplicity))
    )
    _, stein_schur = stein_schur_complement(
        operator,
        metric,
        multiplicity,
    )
    stein_error = float(np.linalg.norm(stein_schur))
    endpoint_maximum = float(
        np.linalg.eigvalsh(upper_endpoint(metric, multiplicity))[-1]
    )
    verified = bool(
        lower_error < 2e-10
        and stein_error < 2e-10
        and endpoint_maximum < 2e-9
    )
    if not verified:
        raise RuntimeError("the zero-slack Hardy bridge audit failed")
    return ZeroSlackBridgeRecord(
        record_kind="canonical_hardy_zero_slack_bridge",
        dimension=(length + 1) * multiplicity,
        length=length,
        multiplicity=multiplicity,
        lower_schur_error=format_float(lower_error),
        stein_schur_error=format_float(stein_error),
        upper_endpoint_maximum_eigenvalue=format_float(endpoint_maximum),
        all_checks_passed=verified,
    )


def endpoint_record(
    length: int,
    multiplicity: int,
    step: float,
) -> FirstResidualEndpointRecord:
    """Audit the operator-valued quadratic endpoint coefficient."""

    size = length * multiplicity
    base = np.eye(size, dtype=complex) / 2
    direction = deterministic_hermitian_direction(length, multiplicity)
    gram = actual_residual_gram(direction, length, multiplicity)
    base_endpoint = canonical_endpoint(base, length, multiplicity)
    positive = canonical_endpoint(
        base + step * direction,
        length,
        multiplicity,
    )
    negative = canonical_endpoint(
        base - step * direction,
        length,
        multiplicity,
    )
    second_coefficient = (
        positive + negative - 2 * base_endpoint
    ) / (2 * step**2)
    expected = -4 * gram
    error = float(np.linalg.norm(second_coefficient - expected))
    expected_norm = float(np.linalg.norm(expected))
    relative_error = error / expected_norm
    gram_minimum = float(np.linalg.eigvalsh(gram)[0])
    endpoint_maximum = float(
        np.linalg.eigvalsh(second_coefficient)[-1]
    )
    verified = bool(
        gram_minimum > -1e-10
        and expected_norm > 1e-5
        and relative_error < 2e-5
        and endpoint_maximum < 2e-5
    )
    if not verified:
        raise RuntimeError(
            "the operator-valued first-residual endpoint audit failed"
        )
    return FirstResidualEndpointRecord(
        record_kind="operator_first_residual_endpoint",
        dimension=(length + 1) * multiplicity,
        length=length,
        multiplicity=multiplicity,
        residual_gram_minimum_eigenvalue=format_float(gram_minimum),
        residual_gram_norm=format_float(float(np.linalg.norm(gram))),
        symmetric_second_coefficient_error=format_float(error),
        relative_second_coefficient_error=format_float(relative_error),
        endpoint_second_coefficient_maximum_eigenvalue=format_float(
            endpoint_maximum
        ),
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
    """Run repeated first-residual endpoint audits."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--minimum-length", type=int, default=2)
    parser.add_argument("--maximum-length", type=int, default=5)
    parser.add_argument(
        "--multiplicities",
        type=int,
        nargs="+",
        default=(2, 3),
    )
    parser.add_argument("--step", type=float, default=3e-4)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_first_residual_endpoint_s70224.jsonl"
        ),
    )
    args = parser.parse_args()

    records: list[AuditRecord] = []
    for multiplicity in args.multiplicities:
        for length in range(args.minimum_length, args.maximum_length + 1):
            bridge = bridge_record(length, multiplicity)
            records.append(bridge)
            print(
                "verified zero-slack bridge "
                f"length {length}, multiplicity {multiplicity}",
                flush=True,
            )
            endpoint = endpoint_record(
                length,
                multiplicity,
                args.step,
            )
            records.append(endpoint)
            print(
                "verified first-residual endpoint "
                f"length {length}, multiplicity {multiplicity}",
                flush=True,
            )
    write_records(args.output, records)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
