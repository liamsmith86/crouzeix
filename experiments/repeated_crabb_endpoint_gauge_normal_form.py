#!/usr/bin/env python3
"""Audit the endpoint-gauge normal form for two Schur states."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import null_space

from crabb_block_hardy_equality import format_float


Matrix = np.ndarray


@dataclass(frozen=True)
class EndpointGaugeRecord:
    """One audit of gauge-plus-hereditary endpoint transport."""

    record_type: str
    endpoint_dimension: int
    interior_dimension: int
    transfer_row_rank: int
    surviving_flag_dimension: int
    old_factorization_error: str
    new_factorization_error: str
    interior_congruence_error: str
    quotient_fixing_error: str
    gauge_normal_form_error: str
    endpoint_difference_error: str
    flag_compression_error: str
    converse_factorization_error: str
    all_checks_passed: bool


def random_matrix(
    generator: np.random.Generator,
    rows: int,
    columns: int,
) -> Matrix:
    """Return one deterministic complex matrix."""

    return (
        generator.standard_normal((rows, columns))
        + 1j * generator.standard_normal((rows, columns))
    )


def positive_matrix(
    generator: np.random.Generator,
    dimension: int,
    floor: float,
) -> Matrix:
    """Return one deterministic positive-definite Hermitian matrix."""

    raw = random_matrix(generator, dimension, dimension)
    return raw.conj().T @ raw + floor * np.eye(dimension)


def hermitian_power(matrix: Matrix, exponent: float) -> Matrix:
    """Return a Hermitian spectral power of a positive matrix."""

    values, vectors = np.linalg.eigh(matrix)
    if values.min() <= 0:
        raise ValueError("the matrix must be positive definite")
    return (vectors * (values**exponent)) @ vectors.conj().T


def assemble_state(
    endpoint: Matrix,
    cross: Matrix,
    interior: Matrix,
) -> Matrix:
    """Assemble a state with the requested Schur endpoint."""

    corner = (
        endpoint
        + cross @ np.linalg.solve(interior, cross.conj().T)
    )
    return np.block(
        [
            [corner, cross],
            [cross.conj().T, interior],
        ]
    )


def schur_endpoint(
    state: Matrix,
    endpoint_dimension: int,
) -> Matrix:
    """Return the leading Schur endpoint."""

    size = endpoint_dimension
    corner = state[:size, :size]
    cross = state[:size, size:]
    interior = state[size:, size:]
    return corner - cross @ np.linalg.solve(
        interior,
        cross.conj().T,
    )


def elimination_factor(
    state: Matrix,
    endpoint_dimension: int,
) -> Matrix:
    """Return the unit lower factor in the Schur decomposition."""

    size = endpoint_dimension
    cross = state[:size, size:]
    interior = state[size:, size:]
    factor = np.eye(len(state), dtype=complex)
    factor[size:, :size] = np.linalg.solve(
        interior,
        cross.conj().T,
    )
    return factor


def interior_congruence(old: Matrix, new: Matrix) -> Matrix:
    """Return the canonical positive C with C* old C = new."""

    old_half = hermitian_power(old, 0.5)
    old_inverse_half = hermitian_power(old, -0.5)
    middle = old_half @ new @ old_half
    return (
        old_inverse_half
        @ hermitian_power(middle, 0.5)
        @ old_inverse_half
    )


def audit_case(
    endpoint_dimension: int,
    interior_dimension: int,
    transfer_rank: int,
    seed: int,
) -> EndpointGaugeRecord:
    """Audit the exact block formulas on one complex case."""

    generator = np.random.default_rng(seed)
    endpoint = positive_matrix(generator, endpoint_dimension, 1.5)
    interior = positive_matrix(generator, interior_dimension, 1.0)
    new_interior = positive_matrix(
        generator,
        interior_dimension,
        1.2,
    )
    cross = 0.12 * random_matrix(
        generator,
        endpoint_dimension,
        interior_dimension,
    )
    new_cross = 0.12 * random_matrix(
        generator,
        endpoint_dimension,
        interior_dimension,
    )

    transfer_left = random_matrix(
        generator,
        endpoint_dimension,
        transfer_rank,
    )
    transfer_right = random_matrix(
        generator,
        transfer_rank,
        2 * endpoint_dimension,
    )
    transfer_row = transfer_left @ transfer_right
    multiplier = 0.015 * random_matrix(
        generator,
        endpoint_dimension,
        2 * endpoint_dimension,
    )
    endpoint_change = (
        multiplier @ transfer_row.conj().T
        + transfer_row @ multiplier.conj().T
    )

    old_state = assemble_state(endpoint, cross, interior)
    new_endpoint = endpoint + endpoint_change
    new_state = assemble_state(
        new_endpoint,
        new_cross,
        new_interior,
    )

    old_factor = elimination_factor(
        old_state,
        endpoint_dimension,
    )
    new_factor = elimination_factor(
        new_state,
        endpoint_dimension,
    )
    old_diagonal = np.block(
        [
            [
                endpoint,
                np.zeros(
                    (endpoint_dimension, interior_dimension),
                    dtype=complex,
                ),
            ],
            [
                np.zeros(
                    (interior_dimension, endpoint_dimension),
                    dtype=complex,
                ),
                interior,
            ],
        ]
    )
    new_diagonal = np.block(
        [
            [
                new_endpoint,
                np.zeros(
                    (endpoint_dimension, interior_dimension),
                    dtype=complex,
                ),
            ],
            [
                np.zeros(
                    (interior_dimension, endpoint_dimension),
                    dtype=complex,
                ),
                new_interior,
            ],
        ]
    )
    old_factor_error = float(
        np.linalg.norm(
            old_state
            - old_factor.conj().T @ old_diagonal @ old_factor
        )
    )
    new_factor_error = float(
        np.linalg.norm(
            new_state
            - new_factor.conj().T @ new_diagonal @ new_factor
        )
    )

    congruence = interior_congruence(interior, new_interior)
    block_congruence = np.block(
        [
            [
                np.eye(endpoint_dimension),
                np.zeros(
                    (endpoint_dimension, interior_dimension),
                    dtype=complex,
                ),
            ],
            [
                np.zeros(
                    (interior_dimension, endpoint_dimension),
                    dtype=complex,
                ),
                congruence,
            ],
        ]
    )
    gauge = (
        np.linalg.solve(old_factor, block_congruence)
        @ new_factor
    )
    injected_endpoint = np.zeros_like(old_state)
    injected_endpoint[
        :endpoint_dimension,
        :endpoint_dimension,
    ] = endpoint_change

    interior_error = float(
        np.linalg.norm(
            congruence.conj().T
            @ interior
            @ congruence
            - new_interior
        )
    )
    quotient_error = float(
        np.linalg.norm(
            gauge[:endpoint_dimension, :endpoint_dimension]
            - np.eye(endpoint_dimension)
        )
        + np.linalg.norm(
            gauge[:endpoint_dimension, endpoint_dimension:]
        )
    )
    normal_form_error = float(
        np.linalg.norm(
            new_state
            - gauge.conj().T @ old_state @ gauge
            - injected_endpoint
        )
    )
    endpoint_error = float(
        np.linalg.norm(
            schur_endpoint(new_state, endpoint_dimension)
            - schur_endpoint(old_state, endpoint_dimension)
            - endpoint_change
        )
    )

    flag = null_space(transfer_row.conj().T)
    flag_error = float(
        np.linalg.norm(
            flag.conj().T @ endpoint_change @ flag
        )
    )

    transfer_pseudoinverse = np.linalg.pinv(transfer_row)
    range_projection = transfer_row @ transfer_pseudoinverse
    reconstructed_multiplier = (
        (
            0.5
            * range_projection
            @ endpoint_change
            @ range_projection
            + (
                np.eye(endpoint_dimension) - range_projection
            )
            @ endpoint_change
            @ range_projection
        )
        @ transfer_pseudoinverse.conj().T
    )
    converse_error = float(
        np.linalg.norm(
            endpoint_change
            - reconstructed_multiplier @ transfer_row.conj().T
            - transfer_row @ reconstructed_multiplier.conj().T
        )
    )

    tolerance = 2e-9
    verified = bool(
        old_factor_error < tolerance
        and new_factor_error < tolerance
        and interior_error < tolerance
        and quotient_error < tolerance
        and normal_form_error < tolerance
        and endpoint_error < tolerance
        and flag_error < tolerance
        and converse_error < tolerance
    )
    if not verified:
        raise RuntimeError(
            "the endpoint-gauge normal-form audit failed: "
            f"factor=({old_factor_error:.3e},"
            f"{new_factor_error:.3e}), "
            f"interior={interior_error:.3e}, "
            f"quotient={quotient_error:.3e}, "
            f"normal={normal_form_error:.3e}, "
            f"endpoint={endpoint_error:.3e}, "
            f"flag={flag_error:.3e}, "
            f"converse={converse_error:.3e}"
        )

    return EndpointGaugeRecord(
        record_type=f"rank_{transfer_rank}",
        endpoint_dimension=endpoint_dimension,
        interior_dimension=interior_dimension,
        transfer_row_rank=int(
            np.linalg.matrix_rank(transfer_row)
        ),
        surviving_flag_dimension=flag.shape[1],
        old_factorization_error=format_float(old_factor_error),
        new_factorization_error=format_float(new_factor_error),
        interior_congruence_error=format_float(interior_error),
        quotient_fixing_error=format_float(quotient_error),
        gauge_normal_form_error=format_float(normal_form_error),
        endpoint_difference_error=format_float(endpoint_error),
        flag_compression_error=format_float(flag_error),
        converse_factorization_error=format_float(converse_error),
        all_checks_passed=verified,
    )


def standard_records() -> list[EndpointGaugeRecord]:
    """Return deterministic rank-changing normal-form audits."""

    return [
        audit_case(3, 4, 1, 702_301),
        audit_case(4, 5, 2, 702_302),
        audit_case(5, 6, 3, 702_303),
        audit_case(5, 7, 4, 702_304),
    ]


def write_records(
    records: list[EndpointGaugeRecord],
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
            "repeated_crabb_endpoint_gauge_normal_form_s70226.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the deterministic audit."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
