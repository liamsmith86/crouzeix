#!/usr/bin/env python3
"""Audit the sharp linear scalar-channel gap on repeated disk anchors."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

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
from repeated_crabb_scalar_channel_gap import (
    contraction_from_seed,
    endpoint_metric,
    polynomial_value,
    positive_root,
    scalar_schur_polynomials,
)


@dataclass(frozen=True)
class LinearGapRecord:
    """One batch of sharp, random, or transfer gap audits."""

    record_kind: str
    dimension: int
    defect_multiplicity: int
    case_count: int
    minimum_exact_bound_margin: str
    minimum_linear_bound_margin: str
    maximum_auxiliary_error: str
    all_checks_passed: bool


def endpoint_angle_bound(corner_norm: float) -> float:
    """Return the sharp abstract bound (4)."""

    return (
        4
        + corner_norm**2
        + corner_norm * np.sqrt(corner_norm**2 + 8)
    ) / 2


def channel_score_bound(score: float) -> float:
    """Return (1) from the squared product-state score."""

    return (4 + score + np.sqrt(score**2 + 8 * score)) / 2


def sharp_rotation_record() -> LinearGapRecord:
    """Verify equality on the three-level unitary rotation family."""

    metric = np.diag([1.0, 2.0, 4.0]).astype(complex)
    root = positive_root(metric)
    inverse_root = np.linalg.inv(root)
    exact_margin = float("inf")
    linear_margin = float("inf")
    unitary_error = 0.0
    samples = np.linspace(0.0, 1.0, 41)
    for corner in samples:
        complement = np.sqrt(1 - corner**2)
        contraction = np.array(
            [
                [0, -complement, corner],
                [0, corner, complement],
                [1, 0, 0],
            ],
            dtype=complex,
        )
        actual = float(
            np.linalg.norm(
                inverse_root @ contraction @ root,
                ord=2,
            )
            ** 2
        )
        exact_margin = min(
            exact_margin,
            endpoint_angle_bound(corner) - actual,
        )
        linear_margin = min(
            linear_margin,
            4 - 4 / 3 * (1 - corner**2) - actual,
        )
        unitary_error = max(
            unitary_error,
            float(
                np.linalg.norm(
                    contraction.conj().T @ contraction
                    - np.eye(3)
                )
            ),
        )

    verified = bool(
        abs(exact_margin) < 3e-12
        and linear_margin > -3e-12
        and unitary_error < 3e-12
    )
    if not verified:
        raise RuntimeError("the sharp endpoint-rotation audit failed")
    return LinearGapRecord(
        record_kind="sharp_unitary_endpoint_rotation",
        dimension=3,
        defect_multiplicity=1,
        case_count=len(samples),
        minimum_exact_bound_margin=format_float(exact_margin),
        minimum_linear_bound_margin=format_float(linear_margin),
        maximum_auxiliary_error=format_float(unitary_error),
        all_checks_passed=verified,
    )


def random_contraction_record(
    multiplicity: int,
    case_count: int,
    seed: int,
) -> LinearGapRecord:
    """Test (5) on deterministic complex contractions."""

    dimension = 2 * multiplicity + 4
    generator = np.random.default_rng(seed)
    exact_margin = float("inf")
    linear_margin = float("inf")
    contraction_error = 0.0
    for _ in range(case_count):
        seed_frame = (
            generator.standard_normal((dimension, 2 * multiplicity))
            + 1j
            * generator.standard_normal((dimension, 2 * multiplicity))
        )
        frame, _ = np.linalg.qr(seed_frame)
        right = frame[:, :multiplicity]
        left = frame[:, multiplicity:]
        metric = endpoint_metric(right, left)
        root = positive_root(metric)
        inverse_root = np.linalg.inv(root)
        contraction = contraction_from_seed(dimension, generator)
        corner = float(
            np.linalg.norm(
                right.conj().T @ contraction @ left,
                ord=2,
            )
        )
        actual = float(
            np.linalg.norm(
                inverse_root @ contraction @ root,
                ord=2,
            )
            ** 2
        )
        exact_margin = min(
            exact_margin,
            endpoint_angle_bound(corner) - actual,
        )
        linear_margin = min(
            linear_margin,
            4 - 4 / 3 * (1 - corner**2) - actual,
        )
        contraction_error = max(
            contraction_error,
            max(0.0, float(np.linalg.norm(contraction, ord=2)) - 1),
        )

    verified = bool(
        exact_margin > -3e-11
        and linear_margin > -3e-11
        and contraction_error < 3e-12
    )
    if not verified:
        raise RuntimeError("the random endpoint-angle audit failed")
    return LinearGapRecord(
        record_kind="random_complex_contractions",
        dimension=dimension,
        defect_multiplicity=multiplicity,
        case_count=case_count,
        minimum_exact_bound_margin=format_float(exact_margin),
        minimum_linear_bound_margin=format_float(linear_margin),
        maximum_auxiliary_error=format_float(contraction_error),
        all_checks_passed=verified,
    )


def transfer_record(
    length: int,
    multiplicity: int,
    strength: float,
    polynomial_count: int,
) -> LinearGapRecord:
    """Test (1) on one noncommuting L193 equality anchor."""

    inverse, _ = strengthened_inverse_toeplitz(
        length,
        multiplicity,
        strength,
    )
    data = canonical_transfer_data(
        np.linalg.inv(inverse),
        length,
        multiplicity,
    )
    coefficients = transfer_coefficients(data, 120)[1:]
    score, _ = maximize_product_state_energy(
        coefficients,
        restarts=72,
        iterations=100,
    )
    metric = endpoint_metric(
        data.right_defect_basis,
        data.left_defect_basis,
    )
    root = positive_root(metric)
    inverse_root = np.linalg.inv(root)
    polynomials = scalar_schur_polynomials(
        polynomial_count,
        min(18, len(data.contraction) - 1),
        322_000 + 10 * length + multiplicity,
    )

    exact_margin = float("inf")
    linear_margin = float("inf")
    corner_score_error = 0.0
    for polynomial in polynomials:
        value = polynomial_value(polynomial, data.contraction)
        corner = (
            data.right_defect_basis.conj().T
            @ value
            @ data.left_defect_basis
        )
        corner_norm = float(np.linalg.norm(corner, ord=2))
        actual = float(
            np.linalg.norm(
                inverse_root @ value @ root,
                ord=2,
            )
            ** 2
        )
        exact_margin = min(
            exact_margin,
            endpoint_angle_bound(corner_norm) - actual,
        )
        linear_margin = min(
            linear_margin,
            channel_score_bound(score) - actual,
        )
        corner_score_error = max(
            corner_score_error,
            max(0.0, corner_norm**2 - score),
        )

    parseval = sum(
        coefficient.conj().T @ coefficient
        for coefficient in coefficients
    )
    parseval_error = float(
        np.linalg.norm(parseval - np.eye(multiplicity))
    )
    auxiliary_error = max(corner_score_error, parseval_error)
    verified = bool(
        0 <= score <= 1 + 3e-10
        and exact_margin > -3e-9
        and linear_margin > -3e-9
        and auxiliary_error < 3e-8
    )
    if not verified:
        raise RuntimeError("the matrix-inner transfer gap audit failed")
    return LinearGapRecord(
        record_kind=f"matrix_inner_transfer_length_{length}",
        dimension=len(data.contraction),
        defect_multiplicity=multiplicity,
        case_count=len(polynomials),
        minimum_exact_bound_margin=format_float(exact_margin),
        minimum_linear_bound_margin=format_float(linear_margin),
        maximum_auxiliary_error=format_float(auxiliary_error),
        all_checks_passed=verified,
    )


def linear_simplification_record() -> LinearGapRecord:
    """Audit (13) on a dense scalar-channel grid."""

    exact_margin = float("inf")
    for score in np.linspace(0.0, 1.0, 10_001):
        exact_margin = min(
            exact_margin,
            4
            - 4 / 3 * (1 - score)
            - channel_score_bound(float(score)),
        )
    verified = bool(exact_margin > -3e-13)
    if not verified:
        raise RuntimeError("the linear simplification audit failed")
    return LinearGapRecord(
        record_kind="global_linear_simplification",
        dimension=0,
        defect_multiplicity=0,
        case_count=10_001,
        minimum_exact_bound_margin=format_float(exact_margin),
        minimum_linear_bound_margin=format_float(exact_margin),
        maximum_auxiliary_error=format_float(0.0),
        all_checks_passed=verified,
    )


def standard_records(
    case_count: int,
    strength: float,
    polynomial_count: int,
) -> list[LinearGapRecord]:
    """Return all standard deterministic audits."""

    records = [
        sharp_rotation_record(),
        linear_simplification_record(),
    ]
    records.extend(
        random_contraction_record(
            multiplicity,
            case_count,
            322_270 + multiplicity,
        )
        for multiplicity in (1, 2, 3)
    )
    records.extend(
        transfer_record(
            length,
            multiplicity,
            strength,
            polynomial_count,
        )
        for multiplicity in (2, 3)
        for length in range(2, 6)
    )
    return records


def write_records(
    records: list[LinearGapRecord],
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
            "repeated_crabb_scalar_channel_linear_gap_s70224.jsonl"
        ),
    )
    parser.add_argument("--case-count", type=int, default=120)
    parser.add_argument("--strength", type=float, default=18.0)
    parser.add_argument("--polynomial-count", type=int, default=48)
    return parser.parse_args()


def main() -> None:
    """Run the sharp scalar-channel linear-gap audit."""

    args = parse_args()
    records = standard_records(
        args.case_count,
        args.strength,
        args.polynomial_count,
    )
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(f"sha256={digest}")


if __name__ == "__main__":
    main()
