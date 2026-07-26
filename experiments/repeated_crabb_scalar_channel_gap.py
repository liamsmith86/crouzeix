#!/usr/bin/env python3
"""Audit the quantitative scalar-channel gap at repeated disk anchors."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Sequence

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_inner_faber_transfer import (
    canonical_transfer_data,
    strengthened_inverse_toeplitz,
    transfer_coefficients,
)
from repeated_crabb_scalar_channel import (
    maximize_product_state_energy,
)


@dataclass(frozen=True)
class AbstractGapRecord:
    """A batch audit of the abstract contraction endpoint gap."""

    record_kind: str
    dimension: int
    defect_multiplicity: int
    case_count: int
    minimum_bound_margin: str
    maximum_contraction_norm: str
    all_checks_passed: bool


@dataclass(frozen=True)
class TransferGapRecord:
    """One actual matrix-inner transfer/channel-gap audit."""

    record_kind: str
    length: int
    multiplicity: int
    toeplitz_strength: str
    transfer_coefficient_count: int
    parseval_error: str
    apex_channel_error: str
    noncommuting_channel_score: str
    noncommuting_channel_gap: str
    sampled_polynomial_count: int
    minimum_direct_gap_margin: str
    minimum_transfer_score_margin: str
    all_checks_passed: bool


AuditRecord = AbstractGapRecord | TransferGapRecord


def positive_root(matrix: np.ndarray) -> np.ndarray:
    """Return the positive Hermitian square root."""

    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    return (
        eigenvectors
        @ np.diag(np.sqrt(eigenvalues))
        @ eigenvectors.conj().T
    )


def endpoint_metric(
    right: np.ndarray,
    left: np.ndarray,
) -> np.ndarray:
    """Return ``2I - VV* + 2WW*``."""

    identity = np.eye(right.shape[0], dtype=complex)
    return (
        2 * identity
        - right @ right.conj().T
        + 2 * left @ left.conj().T
    )


def contraction_from_seed(
    dimension: int,
    generator: np.random.Generator,
) -> np.ndarray:
    """Return a deterministic random strict contraction."""

    raw = (
        generator.standard_normal((dimension, dimension))
        + 1j * generator.standard_normal((dimension, dimension))
    )
    left, _, right = np.linalg.svd(raw)
    singular_values = generator.uniform(0.05, 0.999, dimension)
    return (left * singular_values) @ right


def abstract_gap_record(
    dimension: int,
    multiplicity: int,
    case_count: int,
    seed: int,
) -> AbstractGapRecord:
    """Test L320's abstract endpoint inequality."""

    generator = np.random.default_rng(seed)
    minimum_margin = float("inf")
    maximum_contraction_norm = 0.0
    for _ in range(case_count):
        frame_seed = (
            generator.standard_normal((dimension, 2 * multiplicity))
            + 1j
            * generator.standard_normal((dimension, 2 * multiplicity))
        )
        frame, _ = np.linalg.qr(frame_seed)
        right = frame[:, :multiplicity]
        left = frame[:, multiplicity : 2 * multiplicity]
        metric = endpoint_metric(right, left)
        root = positive_root(metric)
        inverse_root = np.linalg.inv(root)
        contraction = contraction_from_seed(dimension, generator)
        transformed = inverse_root @ contraction @ root
        corner = right.conj().T @ contraction @ left
        score = float(np.linalg.norm(corner, 2))
        bound = 4 - (1 - score) ** 2
        actual = float(np.linalg.norm(transformed, 2) ** 2)
        minimum_margin = min(minimum_margin, bound - actual)
        maximum_contraction_norm = max(
            maximum_contraction_norm,
            float(np.linalg.norm(contraction, 2)),
        )

    verified = bool(
        minimum_margin > -2e-11
        and maximum_contraction_norm <= 1 + 1e-12
    )
    if not verified:
        raise RuntimeError("the abstract scalar-channel gap audit failed")
    return AbstractGapRecord(
        record_kind="abstract_scalar_channel_gap",
        dimension=dimension,
        defect_multiplicity=multiplicity,
        case_count=case_count,
        minimum_bound_margin=format_float(minimum_margin),
        maximum_contraction_norm=format_float(maximum_contraction_norm),
        all_checks_passed=verified,
    )


def scalar_schur_polynomials(
    count: int,
    maximum_degree: int,
    seed: int,
) -> list[np.ndarray]:
    """Return deterministic polynomials with coefficient l1 norm at most one."""

    generator = np.random.default_rng(seed)
    polynomials = [np.array([1.0]), np.array([0.0, 1.0])]
    for _ in range(count - len(polynomials)):
        degree = int(generator.integers(1, maximum_degree + 1))
        coefficients = (
            generator.standard_normal(degree + 1)
            + 1j * generator.standard_normal(degree + 1)
        )
        coefficients /= max(1.0, float(np.sum(np.abs(coefficients))))
        polynomials.append(coefficients)
    return polynomials


def polynomial_value(
    coefficients: np.ndarray,
    matrix: np.ndarray,
) -> np.ndarray:
    """Evaluate a scalar polynomial by Horner's rule."""

    identity = np.eye(len(matrix), dtype=complex)
    value = np.zeros_like(matrix)
    for coefficient in coefficients[::-1]:
        value = value @ matrix + coefficient * identity
    return value


def transfer_gap_record(
    length: int,
    multiplicity: int,
    strength: float,
    coefficient_count: int,
    polynomial_count: int,
) -> TransferGapRecord:
    """Audit L320 on one noncommuting L193 equality anchor."""

    inverse, actual_strength = strengthened_inverse_toeplitz(
        length,
        multiplicity,
        strength,
    )
    data = canonical_transfer_data(
        np.linalg.inv(inverse),
        length,
        multiplicity,
    )
    coefficients = transfer_coefficients(data, coefficient_count)[1:]
    parseval = sum(
        coefficient.conj().T @ coefficient
        for coefficient in coefficients
    )
    parseval_error = float(
        np.linalg.norm(parseval - np.eye(multiplicity))
    )
    score, _ = maximize_product_state_energy(coefficients)

    apex = canonical_transfer_data(
        np.eye(length * multiplicity, dtype=complex) / 2,
        length,
        multiplicity,
    )
    apex_coefficients = transfer_coefficients(
        apex,
        coefficient_count,
    )[1:]
    apex_score, apex_leakage = maximize_product_state_energy(
        apex_coefficients
    )
    apex_error = max(abs(1 - apex_score), apex_leakage)

    metric = endpoint_metric(
        data.right_defect_basis,
        data.left_defect_basis,
    )
    root = positive_root(metric)
    inverse_root = np.linalg.inv(root)
    polynomials = scalar_schur_polynomials(
        polynomial_count,
        min(18, coefficient_count - 2),
        32000 + 10 * length + multiplicity,
    )
    minimum_direct_margin = float("inf")
    minimum_score_margin = float("inf")
    for polynomial in polynomials:
        value = polynomial_value(polynomial, data.contraction)
        transformed = inverse_root @ value @ root
        corner = (
            data.right_defect_basis.conj().T
            @ value
            @ data.left_defect_basis
        )
        corner_score = float(np.linalg.norm(corner, 2) ** 2)
        direct_bound = 4 - (1 - np.sqrt(corner_score)) ** 2
        actual = float(np.linalg.norm(transformed, 2) ** 2)
        minimum_direct_margin = min(
            minimum_direct_margin,
            direct_bound - actual,
        )
        minimum_score_margin = min(
            minimum_score_margin,
            score - corner_score,
        )

    verified = bool(
        parseval_error < 5e-9
        and apex_error < 5e-9
        and 1 - score > 1e-8
        and minimum_direct_margin > -2e-9
        and minimum_score_margin > -2e-8
    )
    if not verified:
        raise RuntimeError("the transfer scalar-channel gap audit failed")
    return TransferGapRecord(
        record_kind="noncommuting_transfer_channel_gap",
        length=length,
        multiplicity=multiplicity,
        toeplitz_strength=format_float(actual_strength),
        transfer_coefficient_count=len(coefficients),
        parseval_error=format_float(parseval_error),
        apex_channel_error=format_float(apex_error),
        noncommuting_channel_score=format_float(score),
        noncommuting_channel_gap=format_float(1 - score),
        sampled_polynomial_count=len(polynomials),
        minimum_direct_gap_margin=format_float(minimum_direct_margin),
        minimum_transfer_score_margin=format_float(minimum_score_margin),
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
    """Run the abstract and actual-transfer gap audits."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/repeated_crabb_scalar_channel_gap_s70224.jsonl"
        ),
    )
    parser.add_argument("--case-count", type=int, default=80)
    parser.add_argument("--strength", type=float, default=18.0)
    parser.add_argument("--coefficient-count", type=int, default=100)
    parser.add_argument("--polynomial-count", type=int, default=36)
    args = parser.parse_args()

    records: list[AuditRecord] = []
    for multiplicity in (1, 2, 3):
        dimension = 2 * multiplicity + 3
        records.append(
            abstract_gap_record(
                dimension,
                multiplicity,
                args.case_count,
                70224 + multiplicity,
            )
        )
    for multiplicity in (2, 3):
        for length in range(2, 6):
            records.append(
                transfer_gap_record(
                    length,
                    multiplicity,
                    args.strength,
                    args.coefficient_count,
                    args.polynomial_count,
                )
            )

    write_records(args.output, records)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
