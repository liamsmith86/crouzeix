#!/usr/bin/env python3
"""Audit the exact delayed model-kernel resolvent formula.

After the right endpoint correction is absorbed into the bulk resolvent,
the remaining finite-delay reflection is a rank-m Krein correction.  Its
small denominator is the characteristic/model kernel of the delayed
transfer ``C(u) = u**r B(u)``.  This checker verifies every state-level
factor in that formula on deterministic noncommuting colligations.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_multidelay_terminal_resolvent import (
    delay_decomposition,
    reflected_self_energy,
)
from repeated_crabb_scattering_schur_collapse import (
    adjoint_transfer,
    transfer,
)
from repeated_crabb_transfer_channel_covariance import inflated_case


Matrix = np.ndarray


@dataclass(frozen=True)
class ModelKernelResolventRecord:
    """One deterministic state-lift audit."""

    delay: int
    state_dimension: int
    defect_dimension: int
    tail_dimension: int
    maximum_right_bulk_error: str
    maximum_right_compression_error: str
    maximum_endpoint_factor_error: str
    maximum_model_kernel_error: str
    maximum_krein_resolvent_error: str
    all_checks_passed: bool


def audit_case(
    delay: int,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
) -> ModelKernelResolventRecord:
    """Audit the model-kernel formula for one delayed colligation."""

    dimension, multiplicity = right.shape
    _, _, tail, tail_right, tail_left = delay_decomposition(
        partial,
        right,
        left,
        delay,
    )
    tail_dimension = len(tail)
    identity = np.eye(tail_dimension, dtype=complex)
    copy_identity = np.eye(multiplicity, dtype=complex)
    right_projection = tail_right @ tail_right.conj().T
    left_projection = tail_left @ tail_left.conj().T

    right_bulk_error = 0.0
    right_compression_error = 0.0
    endpoint_factor_error = 0.0
    model_kernel_error = 0.0
    krein_resolvent_error = 0.0

    for c in (0.03, 0.17, 0.41):
        for zeta in (
            1.7 * np.exp(0.23j),
            1.9 * np.exp(1.17j),
            2.2 * np.exp(-0.71j),
        ):
            x = 1 / zeta
            rho = c / zeta
            t = c / zeta**2
            z = zeta + rho
            gamma = reflected_self_energy(delay, z, c)
            delta = gamma - rho

            bulk = (
                (zeta * identity - tail)
                @ (identity - rho * tail.conj().T)
            )
            right_bulk_denominator = (
                bulk
                - c * tail.conj().T @ right_projection
            )
            right_bulk = np.linalg.inv(right_bulk_denominator)

            reflected = transfer(
                tail,
                tail_right,
                tail_left,
                rho,
            )
            positive = adjoint_transfer(
                tail,
                tail_right,
                tail_left,
                x,
            )
            predicted_compression = (
                x
                / (1 - t)
                * (copy_identity + reflected @ positive)
            )
            actual_compression = (
                tail_left.conj().T
                @ right_bulk
                @ tail_left
            )
            right_compression_error = max(
                right_compression_error,
                float(
                    np.linalg.norm(
                        actual_compression - predicted_compression
                    )
                ),
            )

            direct_right_bulk = np.linalg.inv(
                z * identity
                - tail
                - c
                * tail.conj().T
                @ (identity + right_projection)
                - rho * left_projection
            )
            right_bulk_error = max(
                right_bulk_error,
                float(
                    np.linalg.norm(
                        right_bulk - direct_right_bulk
                    )
                ),
            )

            endpoint = (
                copy_identity - delta * actual_compression
            )
            delayed_transfer = rho**delay * reflected
            delayed_positive = x**delay * positive
            numerator = (
                copy_identity
                - delayed_transfer @ delayed_positive
            )
            predicted_endpoint = numerator / (1 + t**delay)
            endpoint_factor_error = max(
                endpoint_factor_error,
                float(
                    np.linalg.norm(endpoint - predicted_endpoint)
                ),
            )

            model_kernel = numerator / (1 - t)
            split_kernel = sum(
                t**index * copy_identity
                for index in range(delay)
            ) + t**delay * (
                copy_identity - reflected @ positive
            ) / (1 - t)
            model_kernel_error = max(
                model_kernel_error,
                float(np.linalg.norm(model_kernel - split_kernel)),
            )

            direct_denominator = (
                right_bulk_denominator - delta * left_projection
            )
            direct_resolvent = np.linalg.inv(direct_denominator)
            predicted_resolvent = (
                right_bulk
                + zeta
                * t**delay
                * right_bulk
                @ tail_left
                @ np.linalg.inv(model_kernel)
                @ tail_left.conj().T
                @ right_bulk
            )
            krein_resolvent_error = max(
                krein_resolvent_error,
                float(
                    np.linalg.norm(
                        direct_resolvent - predicted_resolvent
                    )
                ),
            )

    tolerance = 8e-11
    verified = max(
        right_bulk_error,
        right_compression_error,
        endpoint_factor_error,
        model_kernel_error,
        krein_resolvent_error,
    ) < tolerance
    if not verified:
        raise RuntimeError(
            "the model-kernel resolvent audit failed: "
            f"delay={delay}, n={dimension}"
        )
    return ModelKernelResolventRecord(
        delay=delay,
        state_dimension=dimension,
        defect_dimension=multiplicity,
        tail_dimension=tail_dimension,
        maximum_right_bulk_error=format_float(right_bulk_error),
        maximum_right_compression_error=format_float(
            right_compression_error
        ),
        maximum_endpoint_factor_error=format_float(
            endpoint_factor_error
        ),
        maximum_model_kernel_error=format_float(model_kernel_error),
        maximum_krein_resolvent_error=format_float(
            krein_resolvent_error
        ),
        all_checks_passed=verified,
    )


def standard_records() -> list[ModelKernelResolventRecord]:
    """Return deterministic noncommuting audits."""

    records: list[ModelKernelResolventRecord] = []
    for delay, multiplicity in (
        (1, 2),
        (2, 2),
        (3, 3),
        (4, 2),
        (5, 2),
        (6, 1),
    ):
        base_dimension = max(3 * multiplicity + 1, 7)
        partial, right, left, _ = inflated_case(
            base_dimension,
            multiplicity,
            delay + 1,
            multiplicity,
            74000 + 10 * delay + multiplicity,
        )
        records.append(
            audit_case(delay, partial, right, left)
        )
    return records


def write_records(
    records: list[ModelKernelResolventRecord],
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
            "repeated_crabb_model_kernel_resolvent_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the model-kernel state-lift audits."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
