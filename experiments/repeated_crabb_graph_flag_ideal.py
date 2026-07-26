#!/usr/bin/env python3
"""Audit hereditary transfer ideals under exact two-graph transport."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import null_space

from crabb_block_hardy_equality import format_float
from repeated_crabb_boundary_metric_flag import rank_chain_case
from repeated_crabb_canonical_quintic_preimage import joint_left_kernel
from repeated_crabb_transfer_channel_covariance import (
    random_partial_isometry,
    transfer_coefficient,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class GraphFlagIdealRecord:
    """One matrix audit of hereditary-ideal Schur transport."""

    record_type: str
    state_dimension: int
    defect_dimension: int
    delay_depth: int
    surviving_flag_dimension: int
    old_graph_equation_error: str
    new_graph_equation_error: str
    old_graph_endpoint_normalization_error: str
    new_graph_endpoint_normalization_error: str
    two_graph_endpoint_error: str
    hereditary_factorization_error: str
    surviving_flag_compression_error: str
    endpoint_difference_norm: str
    factor_bound_surplus: str
    all_checks_passed: bool


def random_hermitian(
    generator: np.random.Generator,
    dimension: int,
) -> Matrix:
    """Return one deterministic complex Hermitian matrix."""

    raw = (
        generator.standard_normal((dimension, dimension))
        + 1j * generator.standard_normal((dimension, dimension))
    )
    return (raw + raw.conj().T) / 2


def endpoint_graph(
    matrix: Matrix,
    endpoint: Matrix,
) -> tuple[Matrix, Matrix]:
    """Return the Schur graph and endpoint relative to one frame."""

    complement = null_space(endpoint.conj().T)
    interior = complement.conj().T @ matrix @ complement
    graph = endpoint - (
        complement
        @ np.linalg.solve(
            interior,
            complement.conj().T @ matrix @ endpoint,
        )
    )
    schur = endpoint.conj().T @ matrix @ graph
    return graph, schur


def audit_case(
    construction_kind: str,
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    delay_depth: int,
    seed: int,
) -> GraphFlagIdealRecord:
    """Audit the exact hereditary factorization on one colligation."""

    generator = np.random.default_rng(seed)
    state_dimension = len(operator)
    defect_dimension = right.shape[1]
    identity = np.eye(state_dimension, dtype=complex)
    initial_projection = right @ right.conj().T
    final_projection = left @ left.conj().T

    old = 3.0 * identity + (
        0.035 * random_hermitian(generator, state_dimension)
    )
    multipliers = [
        0.025 * (
            generator.standard_normal(
                (state_dimension, state_dimension)
            )
            + 1j
            * generator.standard_normal(
                (state_dimension, state_dimension)
            )
        )
        for _ in range(delay_depth)
    ]
    delay_channels = [
        initial_projection
        @ np.linalg.matrix_power(operator, degree)
        @ final_projection
        for degree in range(1, delay_depth + 1)
    ]
    perturbation = np.zeros_like(operator)
    for multiplier, channel in zip(
        multipliers,
        delay_channels,
        strict=True,
    ):
        term = multiplier @ channel
        perturbation += term + term.conj().T
    new = old + perturbation

    old_graph, old_schur = endpoint_graph(old, left)
    new_graph, new_schur = endpoint_graph(new, left)
    endpoint_difference = new_schur - old_schur
    mixed = new_graph.conj().T @ perturbation @ old_graph

    transfers = [
        transfer_coefficient(operator, right, left, degree)
        for degree in range(1, delay_depth + 1)
    ]
    predicted = np.zeros_like(endpoint_difference)
    factor_norm_sum = 0.0
    for multiplier, transfer in zip(
        multipliers,
        transfers,
        strict=True,
    ):
        new_factor = new_graph.conj().T @ multiplier @ right
        old_factor = old_graph.conj().T @ multiplier @ right
        predicted += (
            new_factor @ transfer.conj().T
            + transfer @ old_factor.conj().T
        )
        factor_norm_sum += (
            np.linalg.norm(new_factor, ord=2)
            + np.linalg.norm(old_factor, ord=2)
        ) * np.linalg.norm(transfer, ord=2)

    flag = joint_left_kernel(*transfers)
    flag_error = float(
        np.linalg.norm(
            flag.conj().T @ endpoint_difference @ flag
        )
    )
    complement = identity - final_projection
    old_graph_error = float(
        np.linalg.norm(complement @ old @ old_graph)
    )
    new_graph_error = float(
        np.linalg.norm(complement @ new @ new_graph)
    )
    old_normalization_error = float(
        np.linalg.norm(
            left.conj().T @ old_graph
            - np.eye(defect_dimension)
        )
    )
    new_normalization_error = float(
        np.linalg.norm(
            left.conj().T @ new_graph
            - np.eye(defect_dimension)
        )
    )
    two_graph_error = float(
        np.linalg.norm(endpoint_difference - mixed)
    )
    factor_error = float(np.linalg.norm(mixed - predicted))
    endpoint_norm = float(np.linalg.norm(endpoint_difference, ord=2))
    bound_surplus = factor_norm_sum - endpoint_norm

    tolerance = 3e-10
    verified = bool(
        old_graph_error < tolerance
        and new_graph_error < tolerance
        and old_normalization_error < tolerance
        and new_normalization_error < tolerance
        and two_graph_error < tolerance
        and factor_error < tolerance
        and flag_error < tolerance
        and bound_surplus > -tolerance
    )
    if not verified:
        raise RuntimeError(
            "the hereditary graph-flag audit failed: "
            f"kind={construction_kind}, "
            f"old_graph={old_graph_error:.3e}, "
            f"new_graph={new_graph_error:.3e}, "
            f"two_graph={two_graph_error:.3e}, "
            f"factor={factor_error:.3e}, "
            f"flag={flag_error:.3e}, "
            f"bound={bound_surplus:.3e}"
        )
    return GraphFlagIdealRecord(
        record_type=construction_kind,
        state_dimension=state_dimension,
        defect_dimension=defect_dimension,
        delay_depth=delay_depth,
        surviving_flag_dimension=flag.shape[1],
        old_graph_equation_error=format_float(old_graph_error),
        new_graph_equation_error=format_float(new_graph_error),
        old_graph_endpoint_normalization_error=format_float(
            old_normalization_error
        ),
        new_graph_endpoint_normalization_error=format_float(
            new_normalization_error
        ),
        two_graph_endpoint_error=format_float(two_graph_error),
        hereditary_factorization_error=format_float(factor_error),
        surviving_flag_compression_error=format_float(flag_error),
        endpoint_difference_norm=format_float(endpoint_norm),
        factor_bound_surplus=format_float(bound_surplus),
        all_checks_passed=verified,
    )


def standard_records() -> list[GraphFlagIdealRecord]:
    """Return unstructured and rank-changing hereditary-ideal audits."""

    records = []
    operator, right, left = random_partial_isometry(
        12,
        2,
        np.random.default_rng(702_291),
    )
    records.append(
        audit_case(
            "unstructured_depth_two",
            operator,
            right,
            left,
            2,
            702_292,
        )
    )
    for multiplicity, delay_depth, seed in (
        (3, 1, 702_293),
        (4, 2, 702_294),
        (5, 3, 702_295),
    ):
        operator, right, left, _, _ = rank_chain_case(
            multiplicity,
            seed,
            0.8,
        )
        records.append(
            audit_case(
                f"rank_changing_{multiplicity}",
                operator,
                right,
                left,
                delay_depth,
                seed + 100,
            )
        )
    return records


def write_records(
    records: list[GraphFlagIdealRecord],
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
            "repeated_crabb_graph_flag_ideal_s70226.jsonl"
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
