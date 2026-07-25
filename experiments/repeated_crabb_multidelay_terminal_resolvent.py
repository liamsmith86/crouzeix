#!/usr/bin/env python3
"""Audit the exact terminal resolvent after an arbitrary complete delay.

After ``r`` complete left-wandering stages are removed, the balanced
ellipse pencil has a scalar ``r``-step Jacobi corner.  Eliminating that
corner inserts one scalar continued fraction times the promoted left
defect.  This checker verifies the block pencil, the retained resolvent,
the boundary-metric compression, and the exact first reflected
coefficient on deterministic noncommuting colligations.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from math import comb
from pathlib import Path

import numpy as np
from scipy.linalg import null_space

from crabb_block_hardy_equality import format_float
from repeated_crabb_one_delay_block import direct_pencil
from repeated_crabb_one_image_generator import (
    boundary_metric_coefficient,
)
from repeated_crabb_transfer_channel_covariance import (
    inflated_case,
    transfer_coefficient,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class MultiDelayTerminalRecord:
    """One arbitrary-delay terminal-resolvent audit."""

    delay: int
    state_dimension: int
    defect_dimension: int
    tail_dimension: int
    maximum_earlier_transfer_norm: str
    block_operator_error: str
    maximum_pencil_block_error: str
    maximum_closed_gamma_error: str
    maximum_resolvent_error: str
    maximum_boundary_compression_error: str
    catalan_prefix_exact: bool
    first_reflection_coefficient: str
    all_checks_passed: bool


def delay_decomposition(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    delay: int,
) -> tuple[Matrix, Matrix, Matrix, Matrix, Matrix]:
    """Remove ``delay`` complete left-wandering layers."""

    removed = np.hstack(
        [
            np.linalg.matrix_power(partial, degree) @ left
            for degree in range(delay)
        ]
    )
    retained = null_space(removed.conj().T)
    tail = retained.conj().T @ partial @ retained
    tail_right = retained.conj().T @ right
    tail_left = (
        retained.conj().T
        @ np.linalg.matrix_power(partial, delay)
        @ left
    )
    return removed, retained, tail, tail_right, tail_left


def delay_chain_geometry(
    delay: int,
    multiplicity: int,
) -> tuple[Matrix, Matrix, Matrix]:
    """Return the forward chain, left balance, and last inclusion."""

    chain_dimension = delay * multiplicity
    forward = np.zeros(
        (chain_dimension, chain_dimension),
        dtype=complex,
    )
    identity = np.eye(multiplicity, dtype=complex)
    for level in range(delay - 1):
        row = slice(
            (level + 1) * multiplicity,
            (level + 2) * multiplicity,
        )
        column = slice(level * multiplicity, (level + 1) * multiplicity)
        forward[row, column] = identity

    left_balance = np.eye(chain_dimension, dtype=complex)
    left_balance[:multiplicity, :multiplicity] *= 2

    last = np.zeros((chain_dimension, multiplicity), dtype=complex)
    last[-multiplicity:, :] = identity
    return forward, left_balance, last


def chain_pencil(
    delay: int,
    multiplicity: int,
    parameter: float,
) -> tuple[Matrix, Matrix, Matrix]:
    """Return the clean delay pencil and its last-layer inclusion."""

    forward, left_balance, last = delay_chain_geometry(
        delay,
        multiplicity,
    )
    pencil = forward + parameter * left_balance @ forward.conj().T
    return pencil, left_balance, last


def continuant(
    delay: int,
    spectral_parameter: complex,
    parameter: float,
) -> complex:
    """Return the scalar determinant of the delay Jacobi corner."""

    if delay == 0:
        return 1.0 + 0.0j
    previous_previous = 1.0 + 0.0j
    previous = spectral_parameter
    if delay == 1:
        return previous
    current = spectral_parameter**2 - 2 * parameter
    for _ in range(3, delay + 1):
        previous_previous, previous = previous, current
        current = (
            spectral_parameter * previous
            - parameter * previous_previous
        )
    return current


def reflected_self_energy(
    delay: int,
    spectral_parameter: complex,
    parameter: float,
) -> complex:
    """Return the closed terminal self-energy ``gamma_r``."""

    if delay == 1:
        return 2 * parameter / spectral_parameter
    return (
        parameter
        * continuant(delay - 1, spectral_parameter, parameter)
        / continuant(delay, spectral_parameter, parameter)
    )


def inverse_scalar_series(
    series: list[Fraction],
) -> list[Fraction]:
    """Invert a scalar series with unit constant coefficient."""

    result = [Fraction(0) for _ in series]
    result[0] = 1 / series[0]
    for degree in range(1, len(series)):
        result[degree] = -result[0] * sum(
            series[index] * result[degree - index]
            for index in range(1, degree + 1)
        )
    return result


def normalized_self_energy_series(delay: int) -> list[Fraction]:
    """Return ``gamma_r / (c/z)`` as a series in ``t=c/z^2``."""

    order = delay + 1
    result = [Fraction(0) for _ in range(order)]
    result[0] = 2
    for _ in range(2, delay + 1):
        denominator = [Fraction(1)] + [
            -result[index - 1]
            for index in range(1, order)
        ]
        result = inverse_scalar_series(denominator)
    return result


def catalan_prefix(delay: int) -> tuple[bool, Fraction]:
    """Check the exact first departure from the Catalan fixed point."""

    self_energy = normalized_self_energy_series(delay)
    catalan = [
        Fraction(comb(2 * degree, degree), degree + 1)
        for degree in range(delay + 1)
    ]
    first_degree = delay - 1
    prefix_exact = all(
        self_energy[degree] == catalan[degree]
        for degree in range(first_degree)
    )
    return (
        prefix_exact,
        self_energy[first_degree] - catalan[first_degree],
    )


def predicted_boundary_compression(
    tail: Matrix,
    tail_right: Matrix,
    tail_left: Matrix,
    delay: int,
    half_order: int,
) -> Matrix:
    """Return L221's retained coefficient at ``q**half_order``."""

    right_projection = tail_right @ tail_right.conj().T
    result = np.zeros_like(tail)
    for divisor in range(1, half_order + 1):
        if half_order % divisor:
            continue
        result += (
            (-1) ** (half_order // divisor)
            * np.linalg.matrix_power(tail.conj().T, divisor)
            @ right_projection
            @ np.linalg.matrix_power(tail, divisor)
        )
    if half_order >= delay:
        orbit_degree = half_order - delay
        left_projection = tail_left @ tail_left.conj().T
        result += (
            np.linalg.matrix_power(tail, orbit_degree)
            @ left_projection
            @ np.linalg.matrix_power(
                tail.conj().T,
                orbit_degree,
            )
        )
    return result


def audit_case(
    delay: int,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
) -> MultiDelayTerminalRecord:
    """Audit one complete delay of arbitrary length."""

    dimension, multiplicity = right.shape
    removed, retained, tail, tail_right, tail_left = (
        delay_decomposition(
            partial,
            right,
            left,
            delay,
        )
    )
    basis = np.hstack((removed, retained))
    chain_dimension = delay * multiplicity
    chain_operator, _, last = delay_chain_geometry(
        delay,
        multiplicity,
    )
    identity_copy = np.eye(multiplicity, dtype=complex)
    predicted_operator = np.block(
        [
            [
                chain_operator,
                np.zeros(
                    (chain_dimension, len(tail)),
                    dtype=complex,
                ),
            ],
            [tail_left @ last.conj().T, tail],
        ]
    )
    block_operator_error = float(
        np.linalg.norm(
            basis.conj().T @ partial @ basis
            - predicted_operator
        )
    )

    earlier_transfer_norm = max(
        (
            float(
                np.linalg.norm(
                    transfer_coefficient(
                        partial,
                        right,
                        left,
                        degree,
                    )
                )
            )
            for degree in range(1, delay + 1)
        ),
        default=0.0,
    )

    tail_identity = np.eye(len(tail), dtype=complex)
    tail_right_projection = tail_right @ tail_right.conj().T
    tail_left_projection = tail_left @ tail_left.conj().T
    pencil_block_error = 0.0
    gamma_error = 0.0
    resolvent_error = 0.0
    for parameter in (0.03, 0.17, 0.41):
        clean_chain, left_balance, last = chain_pencil(
            delay,
            multiplicity,
            parameter,
        )
        core_pencil = (
            tail
            + parameter
            * tail.conj().T
            @ (tail_identity + tail_right_projection)
        )
        top_right = (
            parameter
            * left_balance
            @ last
            @ tail_left.conj().T
        )
        predicted_pencil = np.block(
            [
                [clean_chain, top_right],
                [tail_left @ last.conj().T, core_pencil],
            ]
        )
        full_pencil = direct_pencil(
            partial,
            right,
            left,
            parameter,
        )
        pencil_block_error = max(
            pencil_block_error,
            float(
                np.linalg.norm(
                    basis.conj().T @ full_pencil @ basis
                    - predicted_pencil
                )
            ),
        )

        for spectral_parameter in (2.1 + 0.2j, -2.4 + 0.3j):
            chain_resolvent = np.linalg.inv(
                spectral_parameter
                * np.eye(chain_dimension, dtype=complex)
                - clean_chain
            )
            numeric_gamma_matrix = (
                parameter
                * last.conj().T
                @ chain_resolvent
                @ left_balance
                @ last
            )
            closed_gamma = reflected_self_energy(
                delay,
                spectral_parameter,
                parameter,
            )
            gamma_error = max(
                gamma_error,
                float(
                    np.linalg.norm(
                        numeric_gamma_matrix
                        - closed_gamma * identity_copy
                    )
                ),
            )

            full_resolvent = np.linalg.inv(
                spectral_parameter
                * np.eye(dimension, dtype=complex)
                - full_pencil
            )
            retained_resolvent = (
                retained.conj().T @ full_resolvent @ retained
            )
            predicted_resolvent = np.linalg.inv(
                spectral_parameter * tail_identity
                - core_pencil
                - closed_gamma * tail_left_projection
            )
            resolvent_error = max(
                resolvent_error,
                float(
                    np.linalg.norm(
                        retained_resolvent - predicted_resolvent
                    )
                ),
            )

    boundary_error = 0.0
    for half_order in range(1, delay + 5):
        full_coefficient = boundary_metric_coefficient(
            partial,
            right,
            left,
            2 * half_order,
        )
        predicted = predicted_boundary_compression(
            tail,
            tail_right,
            tail_left,
            delay,
            half_order,
        )
        boundary_error = max(
            boundary_error,
            float(
                np.linalg.norm(
                    retained.conj().T
                    @ full_coefficient
                    @ retained
                    - predicted
                )
            ),
        )

    prefix_exact, first_reflection = catalan_prefix(delay)
    tolerance = 5e-11
    verified = (
        max(
            earlier_transfer_norm,
            block_operator_error,
            pencil_block_error,
            gamma_error,
            resolvent_error,
            boundary_error,
        )
        < tolerance
        and prefix_exact
        and first_reflection == 1
    )
    if not verified:
        raise RuntimeError(
            "the arbitrary-delay terminal audit failed: "
            f"delay={delay}, n={dimension}"
        )
    return MultiDelayTerminalRecord(
        delay=delay,
        state_dimension=dimension,
        defect_dimension=multiplicity,
        tail_dimension=len(tail),
        maximum_earlier_transfer_norm=format_float(
            earlier_transfer_norm
        ),
        block_operator_error=format_float(block_operator_error),
        maximum_pencil_block_error=format_float(
            pencil_block_error
        ),
        maximum_closed_gamma_error=format_float(gamma_error),
        maximum_resolvent_error=format_float(resolvent_error),
        maximum_boundary_compression_error=format_float(
            boundary_error
        ),
        catalan_prefix_exact=prefix_exact,
        first_reflection_coefficient=str(first_reflection),
        all_checks_passed=verified,
    )


def standard_records() -> list[MultiDelayTerminalRecord]:
    """Return deterministic noncommuting delay audits."""

    records: list[MultiDelayTerminalRecord] = []
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
            73600 + 10 * delay + multiplicity,
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
    records: list[MultiDelayTerminalRecord],
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
            "repeated_crabb_multidelay_terminal_resolvent_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the arbitrary-delay audits."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
