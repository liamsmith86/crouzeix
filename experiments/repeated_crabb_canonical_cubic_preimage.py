#!/usr/bin/env python3
"""Audit the polynomial preimage of the canonical cubic endpoint.

The cubic upper-gap endpoint of A172--L229 is cancelled by the
perpendicular polynomial column

    C_3 = 3 (S*)^2 V (B_1* B_1).

The companion note proves the identity by an exact four-word Stein
coboundary.  This checker regenerates that word identity and audits
the endpoint formula on unstructured and rank-changing colligations.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_boundary_metric_flag import (
    full_endpoint_gap_series,
    rank_chain_case,
)
from repeated_crabb_canonical_repair_flag_obstruction import (
    canonical_repair_metric_series,
)
from repeated_crabb_delayed_slack_anticommutator import (
    E,
    F,
    IDENTITY,
    Polynomial,
    S,
    STAR,
    add,
    adjoint,
    multiply,
    scale,
    schur_residual_series,
)
from repeated_crabb_elliptic_cokernel import (
    endpoint_motion,
    stein_inverse,
)
from repeated_crabb_transfer_channel_covariance import (
    random_partial_isometry,
    transfer_coefficient,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class CanonicalCubicPreimageRecord:
    """One exact-structure and numerical endpoint audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    parameter_scale: str
    colligation_error: str
    spectral_radius: str
    first_transfer_norm: str
    cubic_endpoint_norm: str
    polynomial_column_norm: str
    cubic_preimage_error: str
    perpendicular_column_error: str
    lower_endpoint_error: str
    metric_stein_residual: str
    exact_coboundary_residual_word_count: int
    exact_upper_compression_word_count: int
    all_checks_passed: bool


def canonical_cubic_forcing() -> Polynomial:
    """Return the balanced forcing defining the cubic repair metric."""

    residual = schur_residual_series(3)
    compact = residual[2]
    cubic = residual[3]
    reflected = multiply(
        add(IDENTITY, F),
        multiply(STAR, add(IDENTITY, E)),
    )
    first_operator = add(
        reflected,
        scale(-1, {"sss": 1}),
    )
    return add(
        scale(-1, cubic),
        scale(
            -1,
            multiply(
                multiply(adjoint(first_operator), compact),
                S,
            ),
        ),
        scale(
            -1,
            multiply(
                multiply(STAR, compact),
                first_operator,
            ),
        ),
    )


def polynomial_preimage_forcing() -> Polynomial:
    """Return ``V C_3* + C_3 V*`` in state-word form."""

    lifted_right_gram = multiply(
        multiply(
            multiply(
                multiply(E, S),
                F,
            ),
            STAR,
        ),
        E,
    )
    return scale(
        3,
        add(
            multiply(lifted_right_gram, {"ss": 1}),
            multiply({"aa": 1}, lifted_right_gram),
        ),
    )


def coboundary_witness() -> Polynomial:
    """Return a four-word witness whose upper compression vanishes."""

    return {
        "aaas": 1,
        "asss": 1,
        "aaasssaa": -1,
        "ssaaasss": -1,
    }


def exact_residuals() -> tuple[int, int]:
    """Return word counts in the two exact proof residuals."""

    witness = coboundary_witness()
    coboundary = add(
        witness,
        scale(
            -1,
            multiply(
                multiply(STAR, witness),
                S,
            ),
        ),
    )
    residual = add(
        canonical_cubic_forcing(),
        polynomial_preimage_forcing(),
        scale(-1, coboundary),
    )
    upper_compression = multiply(
        multiply(F, witness),
        F,
    )
    return len(residual), len(upper_compression)


def polynomial_column(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
) -> Matrix:
    """Return ``3 (S*)^2 V (B_1* B_1)``."""

    first_transfer = transfer_coefficient(
        partial,
        right,
        left,
        1,
    )
    right_gram = first_transfer.conj().T @ first_transfer
    return (
        3
        * np.linalg.matrix_power(partial.conj().T, 2)
        @ right
        @ right_gram
    )


def audit_case(
    construction_kind: str,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    parameter_scale: float,
) -> CanonicalCubicPreimageRecord:
    """Audit one cubic endpoint and its polynomial preimage."""

    upper_gap, _ = canonical_repair_metric_series(
        partial,
        right,
        left,
        3,
    )
    endpoint, constant_error = full_endpoint_gap_series(
        upper_gap,
        left,
    )
    cubic_endpoint = 4 * endpoint[3]

    column = polynomial_column(partial, right, left)
    forcing = (
        right @ column.conj().T
        + column @ right.conj().T
    )
    metric_response = stein_inverse(partial, forcing)
    response = 4 * endpoint_motion(
        partial,
        right,
        left,
        column,
    )
    preimage_error = float(
        np.linalg.norm(cubic_endpoint - response)
    )
    perpendicular_error = float(
        np.linalg.norm(right.conj().T @ column)
    )
    lower_endpoint_error = float(
        np.linalg.norm(
            right.conj().T @ metric_response @ right
        )
    )
    metric_stein_residual = float(
        np.linalg.norm(
            metric_response
            - partial.conj().T @ metric_response @ partial
            - forcing
        )
    )
    identity = np.eye(len(partial), dtype=complex)
    colligation_error = max(
        float(
            np.linalg.norm(
                partial.conj().T @ partial
                - (identity - right @ right.conj().T)
            )
        ),
        float(
            np.linalg.norm(
                partial @ partial.conj().T
                - (identity - left @ left.conj().T)
            )
        ),
        float(np.linalg.norm(right.conj().T @ left)),
    )
    spectral_radius = float(
        np.max(np.abs(np.linalg.eigvals(partial)))
    )
    exact_residual, exact_upper = exact_residuals()
    first_transfer = transfer_coefficient(
        partial,
        right,
        left,
        1,
    )

    tolerance = 3e-8
    verified = bool(
        constant_error < tolerance
        and colligation_error < tolerance
        and spectral_radius < 1
        and preimage_error < tolerance
        and perpendicular_error < tolerance
        and lower_endpoint_error < tolerance
        and metric_stein_residual < tolerance
        and exact_residual == 0
        and exact_upper == 0
    )
    if not verified:
        raise RuntimeError(
            "the canonical cubic polynomial preimage audit failed: "
            f"kind={construction_kind}, "
            f"preimage={preimage_error:.3e}, "
            f"perpendicular={perpendicular_error:.3e}, "
            f"lower={lower_endpoint_error:.3e}, "
            f"Stein={metric_stein_residual:.3e}, "
            f"word_residual={exact_residual}, "
            f"upper_words={exact_upper}"
        )

    return CanonicalCubicPreimageRecord(
        construction_kind=construction_kind,
        state_dimension=len(partial),
        defect_dimension=right.shape[1],
        parameter_scale=format_float(parameter_scale),
        colligation_error=format_float(colligation_error),
        spectral_radius=format_float(spectral_radius),
        first_transfer_norm=format_float(
            float(np.linalg.norm(first_transfer))
        ),
        cubic_endpoint_norm=format_float(
            float(np.linalg.norm(cubic_endpoint))
        ),
        polynomial_column_norm=format_float(
            float(np.linalg.norm(column))
        ),
        cubic_preimage_error=format_float(preimage_error),
        perpendicular_column_error=format_float(
            perpendicular_error
        ),
        lower_endpoint_error=format_float(
            lower_endpoint_error
        ),
        metric_stein_residual=format_float(
            metric_stein_residual
        ),
        exact_coboundary_residual_word_count=exact_residual,
        exact_upper_compression_word_count=exact_upper,
        all_checks_passed=verified,
    )


def standard_records() -> list[CanonicalCubicPreimageRecord]:
    """Return unstructured and rank-changing deterministic records."""

    records: list[CanonicalCubicPreimageRecord] = []
    for defect_dimension in range(1, 5):
        for repetition in range(3):
            state_dimension = (
                3 * defect_dimension + 2 + repetition
            )
            partial, right, left = random_partial_isometry(
                state_dimension,
                defect_dimension,
                np.random.default_rng(
                    110_000 + 100 * defect_dimension + repetition
                ),
            )
            records.append(
                audit_case(
                    "unstructured",
                    partial,
                    right,
                    left,
                    1,
                )
            )

    for defect_dimension in range(2, 7):
        for parameter_scale in (1.0, 0.3):
            partial, right, left, _, _ = rank_chain_case(
                defect_dimension,
                109_200 + defect_dimension,
                parameter_scale,
            )
            records.append(
                audit_case(
                    "rank_chain",
                    partial,
                    right,
                    left,
                    parameter_scale,
                )
            )
    return records


def write_records(
    records: list[CanonicalCubicPreimageRecord],
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
            "repeated_crabb_canonical_cubic_preimage_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the complete audit."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
