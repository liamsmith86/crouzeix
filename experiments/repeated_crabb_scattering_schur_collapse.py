#!/usr/bin/env python3
"""Audit the exact Schur collapse of the two-defect scattering matrix.

L238 reduces a completely delayed retained resolvent to a ``2m`` endpoint
scattering matrix.  The lower-right bulk compression is a characteristic
kernel of the deflated transfer.  After eliminating the right endpoint,
the remaining left reflection has the exact scalar weight

    lambda_r(t) = t^r / (1 + t^r),  t = c / zeta^2.

This checker verifies the characteristic-kernel identity, the resulting
Schur complement, and the exact continued-fraction/Joukowski collapse on
deterministic noncommuting tails.
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
from repeated_crabb_transfer_channel_covariance import inflated_case


Matrix = np.ndarray


@dataclass(frozen=True)
class ScatteringSchurCollapseRecord:
    """One audit of the endpoint Schur collapse."""

    delay: int
    state_dimension: int
    defect_dimension: int
    tail_dimension: int
    maximum_kernel_error: str
    maximum_scattering_error: str
    maximum_scalar_weight_error: str
    maximum_schur_collapse_error: str
    all_checks_passed: bool


def transfer(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    parameter: complex,
) -> Matrix:
    """Return ``W*(I-parameter*T*)^-1 V``."""

    identity = np.eye(len(partial), dtype=complex)
    return (
        left.conj().T
        @ np.linalg.inv(identity - parameter * partial.conj().T)
        @ right
    )


def adjoint_transfer(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    parameter: complex,
) -> Matrix:
    """Return the coefficientwise-adjoint transfer ``B^sharp(parameter)``."""

    identity = np.eye(len(partial), dtype=complex)
    return (
        right.conj().T
        @ np.linalg.inv(identity - parameter * partial)
        @ left
    )


def audit_case(
    delay: int,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
) -> ScatteringSchurCollapseRecord:
    """Audit one completely delayed colligation."""

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

    kernel_error = 0.0
    scattering_error = 0.0
    scalar_weight_error = 0.0
    schur_collapse_error = 0.0
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

            left_resolvent = np.linalg.inv(zeta * identity - tail)
            right_resolvent = np.linalg.inv(
                identity - rho * tail.conj().T
            )
            bulk_inverse = right_resolvent @ left_resolvent

            positive = adjoint_transfer(
                tail,
                tail_right,
                tail_left,
                x,
            )
            reflected = transfer(
                tail,
                tail_right,
                tail_left,
                rho,
            )
            lower_right = (
                tail_left.conj().T
                @ bulk_inverse
                @ tail_left
            )
            predicted_lower_right = (
                x
                / (1 - t)
                * (
                    copy_identity
                    - t * reflected @ positive
                )
            )
            kernel_error = max(
                kernel_error,
                float(
                    np.linalg.norm(
                        lower_right - predicted_lower_right
                    )
                ),
            )

            columns = np.hstack(
                (
                    c * tail.conj().T @ tail_right,
                    delta * tail_left,
                )
            )
            rows = np.vstack(
                (
                    tail_right.conj().T,
                    tail_left.conj().T,
                )
            )
            scattering = rows @ bulk_inverse @ columns
            predicted_scattering = np.block(
                [
                    [
                        t * copy_identity,
                        delta * x * positive,
                    ],
                    [
                        (1 + t) * reflected,
                        delta * predicted_lower_right,
                    ],
                ]
            )
            scattering_error = max(
                scattering_error,
                float(
                    np.linalg.norm(
                        scattering - predicted_scattering
                    )
                ),
            )

            scalar_weight = delta * x / (1 - t)
            predicted_scalar_weight = t**delay / (1 + t**delay)
            scalar_weight_error = max(
                scalar_weight_error,
                abs(scalar_weight - predicted_scalar_weight),
            )

            upper_left = scattering[:multiplicity, :multiplicity]
            upper_right = scattering[:multiplicity, multiplicity:]
            lower_left = scattering[multiplicity:, :multiplicity]
            lower_scattering = scattering[multiplicity:, multiplicity:]
            endpoint_schur = (
                copy_identity
                - lower_scattering
                - lower_left
                @ np.linalg.inv(copy_identity - upper_left)
                @ upper_right
            )
            predicted_schur = (
                copy_identity
                - predicted_scalar_weight
                * (copy_identity + reflected @ positive)
            )
            schur_collapse_error = max(
                schur_collapse_error,
                float(
                    np.linalg.norm(
                        endpoint_schur - predicted_schur
                    )
                ),
            )

    tolerance = 8e-11
    verified = max(
        kernel_error,
        scattering_error,
        scalar_weight_error,
        schur_collapse_error,
    ) < tolerance
    if not verified:
        raise RuntimeError(
            "the scattering Schur-collapse audit failed: "
            f"delay={delay}, n={dimension}"
        )
    return ScatteringSchurCollapseRecord(
        delay=delay,
        state_dimension=dimension,
        defect_dimension=multiplicity,
        tail_dimension=tail_dimension,
        maximum_kernel_error=format_float(kernel_error),
        maximum_scattering_error=format_float(scattering_error),
        maximum_scalar_weight_error=format_float(
            scalar_weight_error
        ),
        maximum_schur_collapse_error=format_float(
            schur_collapse_error
        ),
        all_checks_passed=verified,
    )


def standard_records() -> list[ScatteringSchurCollapseRecord]:
    """Return deterministic noncommuting audits."""

    records: list[ScatteringSchurCollapseRecord] = []
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
            73900 + 10 * delay + multiplicity,
        )
        records.append(
            audit_case(
                delay,
                partial,
                right,
                left,
            )
        )
    return records


def write_records(
    records: list[ScatteringSchurCollapseRecord],
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
            "repeated_crabb_scattering_schur_collapse_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact-collapse audits."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
