#!/usr/bin/env python3
"""Audit the ideal-valued normalization of every physical rooted endpoint.

L310 converts a rooted transfer-channel value to a literal hereditary
endpoint factor, but its standalone polarized response need not vanish
on a complete delay.  For a physical root ``p G_j q``, the multiplier
``K = W* q p V`` has the state realization

    X = W K V* = F q p E.

The relative commutator ``[G_j, X] + [G_j, X]*`` has endpoint equal to
the transfer-channel value minus the literal factor.  Expanding both
defect projections in ``X`` turns it into eight L305 rooted terms whose
roots cancel pairwise.  This checker validates the resulting explicit
ideal-valued response, its delay vanishing, and a uniform length bound.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_all_grade_preimage import generalized_grade_column
from repeated_crabb_boundary_metric_flag import rank_chain_case
from repeated_crabb_elliptic_cokernel import stein_inverse
from repeated_crabb_endpoint_null_gauge import (
    delayed_random_partial_isometry,
)
from repeated_crabb_rooted_bridge_flux import (
    response_endpoint,
    rooted_bridge_split,
    transfer_channel,
    word_value,
)
from repeated_crabb_transfer_channel_covariance import (
    random_partial_isometry,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class RootedEndpointNormalizationRecord:
    """One audit of the physical relative-commutator normalization."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    delay: int
    prefix: str
    suffix: str
    transfer_norm: str
    multiplier_norm: str
    raw_polarized_column_norm: str
    normalized_response_column_norm: str
    normalized_column_bound_slack: str
    defect_projection_expansion_error: str
    relative_forcing_expansion_error: str
    rooted_quotient_cancellation_error: str
    relative_endpoint_identity_error: str
    final_endpoint_split_error: str
    ideal_lift_reconstruction_error: str
    complete_delay_normalization_error: str
    all_checks_passed: bool


def transfer_coefficient(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    delay: int,
) -> Matrix:
    """Return ``B_j = W* (S*)^j V``."""

    return (
        left.conj().T
        @ np.linalg.matrix_power(operator.conj().T, delay)
        @ right
    )


def projected_cross_word_expansion(
    prefix: str,
    suffix: str,
) -> tuple[tuple[float, str], ...]:
    """Expand ``F q p E`` using ``F=I-SS*`` and ``E=I-S*S``."""

    word = suffix + prefix
    return (
        (1.0, word),
        (-1.0, word + "as"),
        (-1.0, "sa" + word),
        (1.0, "sa" + word + "as"),
    )


def ideal_relative_response(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    delay: int,
    prefix: str,
    suffix: str,
) -> tuple[Matrix, Matrix, Matrix, Matrix, float]:
    """Return the expanded relative forcing, column, lift, root error.

    For every monomial ``x`` in ``F q p E``, L305 is applied to
    ``G_j x - x G_j``.  The two terms have exactly the same rooted
    quotient, so their difference is an ideal-valued response.
    """

    dimension = len(operator)
    forcing = np.zeros((dimension, dimension), dtype=complex)
    column = np.zeros_like(right)
    lifted_column = np.zeros((dimension, dimension), dtype=complex)
    root_error = 0.0

    for coefficient, word in projected_cross_word_expansion(
        prefix,
        suffix,
    ):
        (
            right_term,
            right_root,
            right_column,
            _,
            right_lift,
        ) = rooted_bridge_split(
            operator,
            right,
            left,
            delay,
            "",
            word,
        )
        (
            left_term,
            left_root,
            left_column,
            _,
            left_lift,
        ) = rooted_bridge_split(
            operator,
            right,
            left,
            delay,
            word,
            "",
        )
        commutator = right_term - left_term
        forcing += coefficient * (
            commutator + commutator.conj().T
        )
        column += coefficient * (right_column - left_column)
        lifted_column += coefficient * (right_lift - left_lift)
        root_error = max(
            root_error,
            float(np.linalg.norm(right_root - left_root)),
        )

    return forcing, column, lifted_column, root_error


def audit_case(
    construction_kind: str,
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    delay: int,
    prefix: str,
    suffix: str,
    *,
    require_complete_delay: bool = False,
) -> RootedEndpointNormalizationRecord:
    """Audit one physical rooted endpoint normalization."""

    dimension = len(operator)
    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T

    prefix_value = word_value(prefix, operator)
    suffix_value = word_value(suffix, operator)
    transfer = transfer_coefficient(operator, right, left, delay)
    multiplier = (
        left.conj().T @ suffix_value @ prefix_value @ right
    )
    bridge = (
        right @ transfer.conj().T @ left.conj().T
    )

    projected_cross = (
        left_projection
        @ suffix_value
        @ prefix_value
        @ right_projection
    )
    expanded_cross = np.zeros_like(operator)
    for coefficient, word in projected_cross_word_expansion(
        prefix,
        suffix,
    ):
        expanded_cross += coefficient * word_value(word, operator)

    direct_relative_forcing = bridge @ projected_cross
    direct_relative_forcing -= projected_cross @ bridge
    direct_relative_forcing += direct_relative_forcing.conj().T
    (
        expanded_relative_forcing,
        relative_column,
        relative_lift,
        root_error,
    ) = ideal_relative_response(
        operator,
        right,
        left,
        delay,
        prefix,
        suffix,
    )

    channel_root = (
        transfer.conj().T @ multiplier
        + multiplier.conj().T @ transfer
    )
    channel_value = transfer_channel(
        channel_root,
        operator,
        right,
        left,
    )
    endpoint_factor = (
        multiplier @ transfer.conj().T
        + transfer @ multiplier.conj().T
    )
    relative_endpoint = (
        left.conj().T
        @ stein_inverse(operator, direct_relative_forcing)
        @ left
    )
    selected_relative_endpoint = response_endpoint(
        relative_column,
        operator,
        right,
        left,
    )

    (
        rooted_term,
        _,
        rooted_column,
        _,
        rooted_lift,
    ) = rooted_bridge_split(
        operator,
        right,
        left,
        delay,
        prefix,
        suffix,
    )
    actual_endpoint = (
        left.conj().T
        @ stein_inverse(
            operator,
            rooted_term + rooted_term.conj().T,
        )
        @ left
    )
    normalized_column = rooted_column + relative_column
    normalized_lift = rooted_lift + relative_lift
    predicted_endpoint = endpoint_factor + response_endpoint(
        normalized_column,
        operator,
        right,
        left,
    )

    raw_polarized_column, _ = generalized_grade_column(
        operator,
        right,
        left,
        delay,
        multiplier,
    )

    transfer_norm = float(np.linalg.norm(transfer))
    normalized_column_norm = float(np.linalg.norm(normalized_column))
    word_length = len(prefix) + len(suffix)
    bound = (
        5 * len(prefix) + 4 * len(suffix) + 17
    ) * transfer_norm

    projection_error = float(
        np.linalg.norm(projected_cross - expanded_cross)
    )
    forcing_error = float(
        np.linalg.norm(
            direct_relative_forcing
            - expanded_relative_forcing
        )
    )
    relative_endpoint_error = float(
        max(
            np.linalg.norm(
                relative_endpoint
                - (channel_value - endpoint_factor)
            ),
            np.linalg.norm(
                selected_relative_endpoint
                - (channel_value - endpoint_factor)
            ),
        )
    )
    final_endpoint_error = float(
        np.linalg.norm(actual_endpoint - predicted_endpoint)
    )
    ideal_lift_error = float(
        np.linalg.norm(
            normalized_column @ right.conj().T
            - normalized_lift
        )
    )
    perpendicular_error = float(
        np.linalg.norm(right.conj().T @ normalized_column)
    )
    complete_delay_error = (
        max(
            transfer_norm,
            normalized_column_norm,
            float(np.linalg.norm(endpoint_factor)),
        )
        if require_complete_delay
        else 0.0
    )
    bound_slack = bound - normalized_column_norm

    tolerance = 4e-8
    verified = bool(
        projection_error < tolerance
        and forcing_error < tolerance
        and root_error < tolerance
        and relative_endpoint_error < tolerance
        and final_endpoint_error < tolerance
        and ideal_lift_error < tolerance
        and perpendicular_error < tolerance
        and bound_slack > -tolerance
        and complete_delay_error < tolerance
    )
    if not verified:
        raise RuntimeError(
            "the rooted endpoint normalization audit failed: "
            f"{construction_kind=}, {delay=}, {prefix=}, {suffix=}, "
            f"{projection_error=}, {forcing_error=}, {root_error=}, "
            f"{relative_endpoint_error=}, {final_endpoint_error=}, "
            f"{ideal_lift_error=}, {perpendicular_error=}, "
            f"{bound_slack=}, {complete_delay_error=}, "
            f"{word_length=}"
        )

    return RootedEndpointNormalizationRecord(
        construction_kind=construction_kind,
        state_dimension=dimension,
        defect_dimension=right.shape[1],
        delay=delay,
        prefix=prefix or "1",
        suffix=suffix or "1",
        transfer_norm=format_float(transfer_norm),
        multiplier_norm=format_float(
            float(np.linalg.norm(multiplier))
        ),
        raw_polarized_column_norm=format_float(
            float(np.linalg.norm(raw_polarized_column))
        ),
        normalized_response_column_norm=format_float(
            normalized_column_norm
        ),
        normalized_column_bound_slack=format_float(bound_slack),
        defect_projection_expansion_error=format_float(
            projection_error
        ),
        relative_forcing_expansion_error=format_float(
            forcing_error
        ),
        rooted_quotient_cancellation_error=format_float(root_error),
        relative_endpoint_identity_error=format_float(
            relative_endpoint_error
        ),
        final_endpoint_split_error=format_float(final_endpoint_error),
        ideal_lift_reconstruction_error=format_float(ideal_lift_error),
        complete_delay_normalization_error=format_float(
            complete_delay_error
        ),
        all_checks_passed=verified,
    )


def standard_records() -> list[RootedEndpointNormalizationRecord]:
    """Return deterministic general, rank-chain, and delayed audits."""

    word_cases = (
        (1, "", ""),
        (1, "a", "a"),
        (2, "sa", "aa"),
        (3, "asa", "aas"),
    )
    records = []

    for multiplicity in (2, 3, 4):
        operator, right, left = random_partial_isometry(
            3 * multiplicity + 6,
            multiplicity,
            np.random.default_rng(142_000 + multiplicity),
        )
        for delay, prefix, suffix in word_cases:
            records.append(
                audit_case(
                    "unstructured",
                    operator,
                    right,
                    left,
                    delay,
                    prefix,
                    suffix,
                )
            )

    for multiplicity in (2, 3, 4):
        operator, right, left, _, _ = rank_chain_case(
            multiplicity,
            142_100 + multiplicity,
            0.12,
        )
        for delay, prefix, suffix in word_cases[:3]:
            records.append(
                audit_case(
                    "rank_chain",
                    operator,
                    right,
                    left,
                    delay,
                    prefix,
                    suffix,
                )
            )

    generator = np.random.default_rng(142_200)
    for multiplicity in (2, 3, 4):
        operator, right, left = delayed_random_partial_isometry(
            4 * multiplicity,
            multiplicity,
            generator,
        )
        records.append(
            audit_case(
                "complete_first_delay",
                operator,
                right,
                left,
                1,
                "aa",
                "",
                require_complete_delay=True,
            )
        )

    return records


def write_records(
    records: list[RootedEndpointNormalizationRecord],
    output: Path,
) -> str:
    """Write deterministic JSON Lines atomically and return its hash."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(
                json.dumps(asdict(record), sort_keys=True) + "\n"
            )
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
            "repeated_crabb_rooted_endpoint_normalization_s70226.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the rooted endpoint normalization audit."""

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
