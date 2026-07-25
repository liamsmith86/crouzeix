#!/usr/bin/env python3
"""Audit the two-defect scattering reduction after complete delay.

L237 eliminates the clean wandering chain into a scalar self-energy.
On the Joukowski contour, the remaining bulk resolvent factors
exactly.  The two endpoint corrections then enter through a ``2m``
Woodbury matrix whose off-diagonal blocks are evaluations of the
deflated transfer.  This checker verifies that reduction on
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
from repeated_crabb_one_delay_block import direct_pencil
from repeated_crabb_transfer_channel_covariance import inflated_case


Matrix = np.ndarray


@dataclass(frozen=True)
class TwoDefectScatteringRecord:
    """One audit of the factored bulk and Woodbury scattering matrix."""

    delay: int
    state_dimension: int
    defect_dimension: int
    tail_dimension: int
    maximum_bulk_factor_error: str
    maximum_right_diagonal_error: str
    maximum_positive_transfer_error: str
    maximum_reflected_transfer_error: str
    maximum_woodbury_resolvent_error: str
    all_checks_passed: bool


def transfer_value(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    parameter: complex,
) -> Matrix:
    """Return ``W*(I-parameter*S*)^-1 V``."""

    identity = np.eye(len(partial), dtype=complex)
    return (
        left.conj().T
        @ np.linalg.inv(identity - parameter * partial.conj().T)
        @ right
    )


def audit_case(
    delay: int,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
) -> TwoDefectScatteringRecord:
    """Audit one arbitrary complete-delay scattering reduction."""

    dimension, multiplicity = right.shape
    _, retained, tail, tail_right, tail_left = delay_decomposition(
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

    bulk_factor_error = 0.0
    right_diagonal_error = 0.0
    positive_transfer_error = 0.0
    reflected_transfer_error = 0.0
    woodbury_error = 0.0
    for parameter in (0.03, 0.17, 0.41):
        full_pencil = direct_pencil(
            partial,
            right,
            left,
            parameter,
        )
        for zeta in (
            1.7 * np.exp(0.23j),
            1.9 * np.exp(1.17j),
            2.2 * np.exp(-0.71j),
        ):
            reflected_parameter = parameter / zeta
            spectral_parameter = zeta + reflected_parameter
            gamma = reflected_self_energy(
                delay,
                spectral_parameter,
                parameter,
            )
            delta = gamma - reflected_parameter

            left_factor = zeta * identity - tail
            right_factor = (
                identity
                - reflected_parameter * tail.conj().T
            )
            bulk = left_factor @ right_factor
            direct_denominator = (
                spectral_parameter * identity
                - (
                    tail
                    + parameter
                    * tail.conj().T
                    @ (identity + right_projection)
                )
                - gamma * left_projection
            )
            predicted_denominator = (
                bulk
                - parameter
                * tail.conj().T
                @ right_projection
                - delta * left_projection
            )
            bulk_factor_error = max(
                bulk_factor_error,
                float(
                    np.linalg.norm(
                        direct_denominator - predicted_denominator
                    )
                ),
            )

            bulk_inverse = np.linalg.inv(bulk)
            columns = np.hstack(
                (
                    parameter * tail.conj().T @ tail_right,
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

            right_diagonal_error = max(
                right_diagonal_error,
                float(
                    np.linalg.norm(
                        scattering[:multiplicity, :multiplicity]
                        - (parameter / zeta**2) * copy_identity
                    )
                ),
            )

            positive_transfer = (
                tail_right.conj().T
                @ np.linalg.inv(zeta * identity - tail)
                @ tail_left
            )
            expected_positive = delta * positive_transfer
            positive_transfer_error = max(
                positive_transfer_error,
                float(
                    np.linalg.norm(
                        scattering[
                            :multiplicity,
                            multiplicity:,
                        ]
                        - expected_positive
                    )
                ),
            )

            reflected_transfer = transfer_value(
                tail,
                tail_right,
                tail_left,
                reflected_parameter,
            )
            expected_reflected = (
                spectral_parameter / zeta
            ) * reflected_transfer
            reflected_transfer_error = max(
                reflected_transfer_error,
                float(
                    np.linalg.norm(
                        scattering[
                            multiplicity:,
                            :multiplicity,
                        ]
                        - expected_reflected
                    )
                ),
            )

            woodbury_resolvent = (
                bulk_inverse
                + bulk_inverse
                @ columns
                @ np.linalg.inv(
                    np.eye(2 * multiplicity, dtype=complex)
                    - scattering
                )
                @ rows
                @ bulk_inverse
            )
            direct_resolvent = np.linalg.inv(direct_denominator)
            retained_resolvent = (
                retained.conj().T
                @ np.linalg.inv(
                    spectral_parameter
                    * np.eye(dimension, dtype=complex)
                    - full_pencil
                )
                @ retained
            )
            woodbury_error = max(
                woodbury_error,
                float(
                    np.linalg.norm(
                        direct_resolvent - woodbury_resolvent
                    )
                ),
                float(
                    np.linalg.norm(
                        retained_resolvent - woodbury_resolvent
                    )
                ),
            )

    tolerance = 8e-11
    verified = max(
        bulk_factor_error,
        right_diagonal_error,
        positive_transfer_error,
        reflected_transfer_error,
        woodbury_error,
    ) < tolerance
    if not verified:
        raise RuntimeError(
            "the two-defect scattering audit failed: "
            f"delay={delay}, n={dimension}"
        )
    return TwoDefectScatteringRecord(
        delay=delay,
        state_dimension=dimension,
        defect_dimension=multiplicity,
        tail_dimension=tail_dimension,
        maximum_bulk_factor_error=format_float(bulk_factor_error),
        maximum_right_diagonal_error=format_float(
            right_diagonal_error
        ),
        maximum_positive_transfer_error=format_float(
            positive_transfer_error
        ),
        maximum_reflected_transfer_error=format_float(
            reflected_transfer_error
        ),
        maximum_woodbury_resolvent_error=format_float(
            woodbury_error
        ),
        all_checks_passed=verified,
    )


def standard_records() -> list[TwoDefectScatteringRecord]:
    """Return deterministic noncommuting scattering audits."""

    records: list[TwoDefectScatteringRecord] = []
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
            73800 + 10 * delay + multiplicity,
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
    records: list[TwoDefectScatteringRecord],
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
            "repeated_crabb_two_defect_scattering_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the scattering audits."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
