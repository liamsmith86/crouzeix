#!/usr/bin/env python3
"""Audit the causal Toeplitz-window form of L252's radial energy."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_transfer_channel_covariance import (
    inflated_case,
    random_partial_isometry,
    transfer_coefficient,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class ToeplitzWindowEnergyRecord:
    """One deterministic finite-window identity audit."""

    window_grade: int
    complete_delay: bool
    state_dimension: int
    defect_dimension: int
    radial_toeplitz_error: str
    coefficient_sum_error: str
    model_compression_error: str
    delayed_single_cell_error: str
    all_checks_passed: bool


def transfer_coefficients(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    maximum_degree: int,
) -> list[Matrix]:
    """Return transfer coefficients from zero through one degree."""

    zero = np.zeros((right.shape[1], right.shape[1]), dtype=complex)
    return [
        zero,
        *(
            transfer_coefficient(
                partial,
                right,
                left,
                degree,
            )
            for degree in range(1, maximum_degree + 1)
        ),
    ]


def toeplitz_window(
    coefficients: list[Matrix],
    grade: int,
) -> Matrix:
    """Return the causal Toeplitz compression on rows zero to grade."""

    multiplicity = len(coefficients[0])
    window = np.zeros(
        ((grade + 1) * multiplicity, (grade + 1) * multiplicity),
        dtype=complex,
    )
    for row in range(grade + 1):
        for column in range(row + 1):
            degree = row - column
            if degree >= len(coefficients):
                continue
            window[
                row * multiplicity : (row + 1) * multiplicity,
                column * multiplicity : (column + 1) * multiplicity,
            ] = coefficients[degree]
    return window


def model_kernel_window(
    coefficients: list[Matrix],
    grade: int,
) -> Matrix:
    """Return the model projection compressed to one Hardy window."""

    toeplitz = toeplitz_window(coefficients, grade)
    identity = np.eye(len(toeplitz), dtype=complex)
    return identity - toeplitz @ toeplitz.conj().T


def audit_case(
    grade: int,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    complete_delay: bool,
) -> ToeplitzWindowEnergyRecord:
    """Audit one state/Hardy finite-window identity."""

    coefficients = transfer_coefficients(
        partial,
        right,
        left,
        grade,
    )
    toeplitz = toeplitz_window(coefficients, grade)
    toeplitz_energy = float(np.linalg.norm(toeplitz) ** 2)
    coefficient_energy = sum(
        (grade + 1 - degree)
        * float(np.linalg.norm(coefficients[degree]) ** 2)
        for degree in range(1, grade + 1)
    )

    dimension = len(partial)
    multiplicity = right.shape[1]
    identity = np.eye(dimension, dtype=complex)
    q_one = partial.conj().T @ partial
    q_top = (
        np.linalg.matrix_power(partial.conj().T, grade + 2)
        @ np.linalg.matrix_power(partial, grade + 2)
    )
    radial_energy = float(
        np.trace(
            q_top
            - (grade + 2) * q_one
            + (grade + 1) * identity
        ).real
    )

    kernel = model_kernel_window(coefficients, grade)
    model_energy = (grade + 1) * multiplicity - float(
        np.trace(kernel).real
    )

    if complete_delay:
        delayed_energy = float(
            np.linalg.norm(coefficients[grade]) ** 2
        )
    else:
        delayed_energy = toeplitz_energy

    radial_error = abs(radial_energy - toeplitz_energy)
    coefficient_error = abs(coefficient_energy - toeplitz_energy)
    model_error = abs(model_energy - toeplitz_energy)
    delayed_error = abs(delayed_energy - toeplitz_energy)
    tolerance = 2e-9
    verified = bool(
        radial_error < tolerance
        and coefficient_error < tolerance
        and model_error < tolerance
        and delayed_error < tolerance
    )
    if not verified:
        raise RuntimeError(
            "Toeplitz-window energy audit failed: "
            f"grade={grade}, delayed={complete_delay}, "
            f"radial={radial_error:.3e}, "
            f"coefficients={coefficient_error:.3e}, "
            f"model={model_error:.3e}, "
            f"single={delayed_error:.3e}"
        )
    return ToeplitzWindowEnergyRecord(
        window_grade=grade,
        complete_delay=complete_delay,
        state_dimension=dimension,
        defect_dimension=multiplicity,
        radial_toeplitz_error=format_float(radial_error),
        coefficient_sum_error=format_float(coefficient_error),
        model_compression_error=format_float(model_error),
        delayed_single_cell_error=format_float(delayed_error),
        all_checks_passed=verified,
    )


def standard_records() -> list[ToeplitzWindowEnergyRecord]:
    """Return general and completely delayed deterministic audits."""

    records: list[ToeplitzWindowEnergyRecord] = []
    for grade in range(1, 7):
        multiplicity = 2 if grade <= 4 else 1
        general = random_partial_isometry(
            grade + 7,
            multiplicity,
            np.random.default_rng(253_000 + grade),
        )
        records.append(
            audit_case(grade, *general, complete_delay=False)
        )

        if grade == 1:
            delayed = random_partial_isometry(
                8,
                multiplicity,
                np.random.default_rng(253_100 + grade),
            )
        else:
            delayed = inflated_case(
                grade + 3,
                multiplicity,
                grade,
                multiplicity,
                253_100 + grade,
            )[:3]
        records.append(
            audit_case(grade, *delayed, complete_delay=True)
        )
    return records


def write_records(
    records: list[ToeplitzWindowEnergyRecord],
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
            "repeated_crabb_toeplitz_window_energy_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the deterministic audits."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
