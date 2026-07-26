#!/usr/bin/env python3
"""Audit the fixed-half finite-flag margin assembly.

L318 uses the exact endpoint faces

    lower retained coefficient = 1 - theta,
    upper-gap direct coefficient = 12 + 4 theta,

at ``theta=1/2``.  L317 puts every other finite-jet block behind one
extra power of ``c``.  This checker revalidates the endpoint direction,
tests the global middle-space Weyl bound on noncommuting Hermitian
remainders, and attacks the terminal analytic-tail estimate.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_all_grade_preimage import (
    stable_random_partial_isometry,
)
from repeated_crabb_oriented_retightening_transport import (
    oriented_retightening_direction,
)
from repeated_crabb_transfer_flag import transfer_coefficient


Matrix = np.ndarray
THETA = 0.5
LOWER_DIRECT = 1 - THETA
UPPER_DIRECT = 12 + 4 * THETA


@dataclass(frozen=True)
class FixedHalfFlagRecord:
    """One endpoint-face, middle-bound, and tail audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    grade: int
    lower_direction_error: str
    physical_upper_direction_error: str
    retained_lower_coefficient: str
    retained_upper_gap_coefficient: str
    remainder_norm: str
    ellipse_parameter: str
    minimum_middle_eigenvalue: str
    weyl_middle_lower_bound: str
    terminal_minimum_singular_value: str
    analytic_tail_constant: str
    minimum_full_endpoint_eigenvalue: str
    terminal_lower_bound: str
    all_checks_passed: bool


def hermitian_with_norm(
    dimension: int,
    norm: float,
    generator: np.random.Generator,
) -> Matrix:
    """Return a deterministic noncommuting Hermitian matrix."""

    raw = (
        generator.standard_normal((dimension, dimension))
        + 1j * generator.standard_normal((dimension, dimension))
    )
    hermitian = (raw + raw.conj().T) / 2
    return norm * hermitian / np.linalg.norm(hermitian, ord=2)


def random_unitary(
    dimension: int,
    generator: np.random.Generator,
) -> Matrix:
    """Return a deterministic Haar-style unitary."""

    raw = (
        generator.standard_normal((dimension, dimension))
        + 1j * generator.standard_normal((dimension, dimension))
    )
    unitary, triangular = np.linalg.qr(raw)
    diagonal = np.diag(triangular)
    phases = np.where(
        np.abs(diagonal) > 0,
        diagonal / np.abs(diagonal),
        1,
    )
    return unitary @ np.diag(np.conjugate(phases))


def audit_case(
    construction_kind: str,
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    grade: int,
    seed: int,
) -> FixedHalfFlagRecord:
    """Audit one physical direction and one adversarial block margin."""

    direction = oriented_retightening_direction(
        operator,
        right,
        left,
        grade,
    )
    transfer = transfer_coefficient(
        operator,
        right,
        left,
        grade,
    )
    right_gram = transfer.conj().T @ transfer
    left_gram = transfer @ transfer.conj().T
    lower_error = float(
        np.linalg.norm(
            right.conj().T @ direction.metric @ right + right_gram
        )
    )
    physical_upper_error = float(
        np.linalg.norm(
            4 * left.conj().T @ direction.metric @ left
            + 4 * left_gram
        )
    )

    generator = np.random.default_rng(seed)
    multiplicity = right.shape[1]
    block_dimension = grade * multiplicity
    remainder_norm = 2.0 + 0.4 * grade
    remainder = hermitian_with_norm(
        block_dimension,
        remainder_norm,
        generator,
    )
    parameter = min(
        0.19,
        0.91 * (UPPER_DIRECT / 2) / remainder_norm,
    )
    middle = (
        UPPER_DIRECT * np.eye(block_dimension, dtype=complex)
        + parameter * remainder
    )
    minimum_middle = float(np.min(np.linalg.eigvalsh(middle)))
    weyl_lower = UPPER_DIRECT - parameter * remainder_norm

    terminal = random_unitary(multiplicity, generator)
    terminal += (
        0.07
        * (
            generator.standard_normal((multiplicity, multiplicity))
            + 1j
            * generator.standard_normal(
                (multiplicity, multiplicity)
            )
        )
        / np.sqrt(multiplicity)
    )
    terminal_minimum = float(
        np.min(np.linalg.svd(terminal, compute_uv=False))
    )
    channels = [
        0.08
        * (
            generator.standard_normal((multiplicity, multiplicity))
            + 1j
            * generator.standard_normal(
                (multiplicity, multiplicity)
            )
        )
        / np.sqrt(multiplicity)
        for _ in range(grade - 1)
    ]
    channels.append(terminal)
    weighted = np.hstack(
        [
            parameter**index * channel
            for index, channel in enumerate(channels, start=1)
        ]
    )
    finite_endpoint = weighted @ middle @ weighted.conj().T

    tail_constant = 3.0 + 0.3 * grade
    tail_parameter_limit = (
        7 * terminal_minimum**2 / (2 * tail_constant)
    )
    if parameter > tail_parameter_limit:
        parameter = 0.91 * tail_parameter_limit
        middle = (
            UPPER_DIRECT * np.eye(block_dimension, dtype=complex)
            + parameter * remainder
        )
        minimum_middle = float(np.min(np.linalg.eigvalsh(middle)))
        weyl_lower = UPPER_DIRECT - parameter * remainder_norm
        weighted = np.hstack(
            [
                parameter**index * channel
                for index, channel in enumerate(channels, start=1)
            ]
        )
        finite_endpoint = weighted @ middle @ weighted.conj().T

    tail = -tail_constant * parameter ** (2 * grade + 1) * np.eye(
        multiplicity,
        dtype=complex,
    )
    full_endpoint = finite_endpoint + tail
    minimum_full = float(np.min(np.linalg.eigvalsh(full_endpoint)))
    terminal_lower = (
        parameter ** (2 * grade)
        * (
            7 * terminal_minimum**2
            - tail_constant * parameter
        )
    )

    verified = bool(
        lower_error < 2e-9
        and physical_upper_error < 2e-9
        and LOWER_DIRECT == 0.5
        and UPPER_DIRECT == 14
        and minimum_middle >= weyl_lower - 1e-10
        and weyl_lower >= 7
        and minimum_full >= terminal_lower - 1e-10
        and terminal_lower > 0
    )
    return FixedHalfFlagRecord(
        construction_kind=construction_kind,
        state_dimension=len(operator),
        defect_dimension=multiplicity,
        grade=grade,
        lower_direction_error=format_float(lower_error),
        physical_upper_direction_error=format_float(
            physical_upper_error
        ),
        retained_lower_coefficient=format_float(LOWER_DIRECT),
        retained_upper_gap_coefficient=format_float(UPPER_DIRECT),
        remainder_norm=format_float(remainder_norm),
        ellipse_parameter=format_float(parameter),
        minimum_middle_eigenvalue=format_float(minimum_middle),
        weyl_middle_lower_bound=format_float(weyl_lower),
        terminal_minimum_singular_value=format_float(terminal_minimum),
        analytic_tail_constant=format_float(tail_constant),
        minimum_full_endpoint_eigenvalue=format_float(minimum_full),
        terminal_lower_bound=format_float(terminal_lower),
        all_checks_passed=verified,
    )


def standard_records() -> list[FixedHalfFlagRecord]:
    """Return deterministic unstructured audits through grade six."""

    records: list[FixedHalfFlagRecord] = []
    generator = np.random.default_rng(318_000)
    cases = [
        (
            "unstructured",
            *stable_random_partial_isometry(
                dimension,
                multiplicity,
                generator,
            ),
        )
        for dimension, multiplicity in ((7, 2), (9, 3), (12, 3))
    ]
    for case_index, case in enumerate(cases):
        construction_kind, operator, right, left = case
        for grade in range(1, 7):
            records.append(
                audit_case(
                    construction_kind,
                    operator,
                    right,
                    left,
                    grade,
                    318_100 + 10 * case_index + grade,
                )
            )
    return records


def write_records(
    records: list[FixedHalfFlagRecord],
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
            "experiments/repeated_crabb_fixed_half_flag_s70226.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist every fixed-half flag audit."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(
        json.dumps(
            {"records": len(records), "sha256": digest},
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
