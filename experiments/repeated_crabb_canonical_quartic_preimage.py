#!/usr/bin/env python3
"""Audit an explicit polynomial preparation of the quartic endpoint.

The ten-term perpendicular column in ``quartic_column_lift`` changes
L231's complete quartic upper gap to the positive matrix

    12 B_2 B_2*
    + 32 B_1 B_1*
    + 56 B_1 (B_1* B_1) B_1*.

The proof is a finite rational word identity with a 29-term Hermitian
Stein witness.  This checker regenerates both exact identities and
audits the physical endpoint on unstructured, rank-changing, and
completely delayed colligations.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_boundary_metric_flag import rank_chain_case
from repeated_crabb_canonical_quartic_trace import (
    exact_quartic_components,
    prepared_endpoint_series,
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
)
from repeated_crabb_elliptic_cokernel import (
    endpoint_motion,
    stein_inverse,
)
from repeated_crabb_endpoint_null_gauge import (
    delayed_random_partial_isometry,
)
from repeated_crabb_transfer_channel_covariance import (
    random_partial_isometry,
    transfer_coefficient,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class CanonicalQuarticPreimageRecord:
    """One exact-structure and numerical quartic preparation audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    parameter_scale: str
    colligation_error: str
    spectral_radius: str
    first_transfer_norm: str
    second_transfer_norm: str
    quartic_endpoint_norm: str
    positive_target_norm: str
    positive_target_minimum_eigenvalue: str
    polynomial_column_norm: str
    quartic_preimage_error: str
    perpendicular_column_error: str
    lower_endpoint_motion_error: str
    exact_coboundary_residual_word_count: int
    exact_endpoint_residual_word_count: int
    exact_witness_hermitian_residual_word_count: int
    all_checks_passed: bool


def quartic_column_bracket() -> Polynomial:
    """Return the ten-term polynomial inside the left projection."""

    return add(
        scale(10, {"ssaa": 1}),
        scale(-1, {"aassssaa": 1}),
        {"ssaaaaaa": 1},
        scale(-1, {"ssssssaa": 1}),
        scale(-7, {"aaaaasssaa": 1}),
        scale(3, {"aaaasssaaa": 1}),
        {"aaasssaaaa": 1},
        {"aasssaaaaa": 1},
        scale(-1, {"saaassssaa": 1}),
        {"ssaaasssaa": 1},
    )


def quartic_column_lift() -> Polynomial:
    """Return ``C_4 V*`` as an exact state polynomial."""

    initial_projection = add(IDENTITY, scale(-1, E))
    return scale(
        Fraction(1, 2),
        multiply(
            multiply(
                initial_projection,
                quartic_column_bracket(),
            ),
            E,
        ),
    )


def quartic_coboundary_witness() -> Polynomial:
    """Return the Hermitian 29-term quartic Stein witness."""

    terms: tuple[tuple[Fraction | int, Polynomial], ...] = (
        (-53, {"": 1}),
        (31, {"as": 1}),
        (21, {"sa": 1}),
        (2, {"aaaa": 1, "ssss": 1}),
        (19, {"aass": 1}),
        (33, {"ssaa": 1}),
        (3, {"aaasss": 1}),
        (-14, {"aasssa": 1, "saaass": 1}),
        (-23, {"asssaa": 1, "ssaaas": 1}),
        (-1, {"sssaaa": 1}),
        (-1, {"aaassssa": 1, "saaaasss": 1}),
        (
            Fraction(-5, 2),
            {"aassssaa": 1, "ssaaaass": 1},
        ),
        (15, {"asssaaas": 1}),
        (12, {"saaasssa": 1}),
        (
            Fraction(1, 2),
            {"aaaaaasssa": 1, "saaassssss": 1},
        ),
        (
            Fraction(1, 2),
            {"aaaaasssaa": 1, "ssaaasssss": 1},
        ),
        (
            Fraction(-1, 2),
            {"aaaasssaaa": 1, "sssaaassss": 1},
        ),
        (
            Fraction(1, 2),
            {"aaasssaaaa": 1, "ssssaaasss": 1},
        ),
        (
            Fraction(1, 2),
            {"aasssaaaaa": 1, "sssssaaass": 1},
        ),
        (10, {"aasssaaass": 1}),
        (
            Fraction(3, 2),
            {"aassssaaas": 1, "asssaaaass": 1},
        ),
        (
            Fraction(-1, 2),
            {"aasssssssa": 1, "saaaaaaass": 1},
        ),
        (
            Fraction(1, 2),
            {"asssaaaaaa": 1, "ssssssaaas": 1},
        ),
        (
            Fraction(-1, 2),
            {"asssssssaa": 1, "ssaaaaaaas": 1},
        ),
        (
            Fraction(3, 2),
            {"saaassssaa": 1, "ssaaaasssa": 1},
        ),
        (16, {"ssaaasssaa": 1}),
        (
            Fraction(-13, 2),
            {"aasssaaasssa": 1, "saaasssaaass": 1},
        ),
        (
            Fraction(-3, 2),
            {"asssaaaasssa": 1, "saaassssaaas": 1},
        ),
        (
            Fraction(-15, 2),
            {"asssaaasssaa": 1, "ssaaasssaaas": 1},
        ),
    )
    result: Polynomial = {}
    for coefficient, polynomial in terms:
        result = add(result, scale(coefficient, polynomial))
    return result


def positive_target_lift() -> Polynomial:
    """Return one quarter of the positive physical endpoint lift."""

    first_left_gram = multiply(
        multiply(
            multiply(
                multiply(F, STAR),
                E,
            ),
            S,
        ),
        F,
    )
    second_left_gram = multiply(
        multiply(
            multiply(
                multiply(F, {"aa": 1}),
                E,
            ),
            {"ss": 1},
        ),
        F,
    )
    first_right_gram = multiply(
        multiply(
            multiply(
                multiply(E, S),
                F,
            ),
            STAR,
        ),
        E,
    )
    nonlinear_left_gram = multiply(
        multiply(
            multiply(
                multiply(F, STAR),
                first_right_gram,
            ),
            S,
        ),
        F,
    )
    return add(
        scale(8, first_left_gram),
        scale(3, second_left_gram),
        scale(14, nonlinear_left_gram),
    )


@lru_cache(maxsize=1)
def exact_residuals() -> tuple[int, int, int]:
    """Return word counts in the exact quartic certificate residuals."""

    forcing, _, upper_schur_square = exact_quartic_components()
    column_lift = quartic_column_lift()
    witness = quartic_coboundary_witness()
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
    forcing_residual = add(
        forcing,
        column_lift,
        adjoint(column_lift),
        scale(-1, coboundary),
    )
    endpoint_residual = add(
        multiply(multiply(F, witness), F),
        multiply(F, upper_schur_square),
        positive_target_lift(),
    )
    hermitian_residual = add(
        witness,
        scale(-1, adjoint(witness)),
    )
    return (
        len(forcing_residual),
        len(endpoint_residual),
        len(hermitian_residual),
    )


def evaluate_polynomial(
    polynomial: Polynomial,
    partial: Matrix,
) -> Matrix:
    """Evaluate one reduced state-word polynomial."""

    result = np.zeros_like(partial)
    adjoint_partial = partial.conj().T
    for word, coefficient in polynomial.items():
        value = np.eye(len(partial), dtype=complex)
        for letter in word:
            value = value @ (
                partial if letter == "s" else adjoint_partial
            )
        result += float(coefficient) * value
    return result


def quartic_column(
    partial: Matrix,
    right: Matrix,
) -> Matrix:
    """Return the explicit perpendicular quartic column."""

    return evaluate_polynomial(
        quartic_column_lift(),
        partial,
    ) @ right


def positive_quartic_target(
    first: Matrix,
    second: Matrix,
) -> Matrix:
    """Return the three-Gram positive quartic target."""

    first_right_gram = first.conj().T @ first
    return (
        12 * second @ second.conj().T
        + 32 * first @ first.conj().T
        + 56 * first @ first_right_gram @ first.conj().T
    )


def audit_case(
    construction_kind: str,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    parameter_scale: float,
) -> CanonicalQuarticPreimageRecord:
    """Audit one exact quartic polynomial preparation."""

    upper, _ = prepared_endpoint_series(partial, right, left)
    first = transfer_coefficient(partial, right, left, 1)
    second = transfer_coefficient(partial, right, left, 2)
    target = positive_quartic_target(first, second)
    target = (target + target.conj().T) / 2

    column = quartic_column(partial, right)
    response = 4 * endpoint_motion(
        partial,
        right,
        left,
        column,
    )
    preimage_error = float(
        np.linalg.norm(upper[4] - response - target)
    )
    perpendicular_error = float(
        np.linalg.norm(right.conj().T @ column)
    )
    forcing = (
        right @ column.conj().T
        + column @ right.conj().T
    )
    metric_response = stein_inverse(partial, forcing)
    lower_endpoint_error = float(
        np.linalg.norm(
            right.conj().T @ metric_response @ right
        )
    )

    identity = np.eye(len(partial), dtype=complex)
    colligation_error = max(
        float(
            np.linalg.norm(
                identity
                - partial.conj().T @ partial
                - right @ right.conj().T
            )
        ),
        float(
            np.linalg.norm(
                identity
                - partial @ partial.conj().T
                - left @ left.conj().T
            )
        ),
        float(np.linalg.norm(right.conj().T @ left)),
    )
    spectral_radius = float(
        np.max(np.abs(np.linalg.eigvals(partial)))
    )
    target_minimum = float(np.linalg.eigvalsh(target)[0])
    exact_forcing, exact_endpoint, exact_hermitian = exact_residuals()

    tolerance = 3e-8
    verified = bool(
        colligation_error < tolerance
        and spectral_radius < 1
        and target_minimum > -tolerance
        and preimage_error < tolerance
        and perpendicular_error < tolerance
        and lower_endpoint_error < tolerance
        and exact_forcing == 0
        and exact_endpoint == 0
        and exact_hermitian == 0
    )
    if not verified:
        raise RuntimeError(
            "the canonical quartic preimage audit failed: "
            f"kind={construction_kind}, "
            f"preimage={preimage_error:.3e}, "
            f"perpendicular={perpendicular_error:.3e}, "
            f"lower={lower_endpoint_error:.3e}, "
            f"exact=({exact_forcing},{exact_endpoint},"
            f"{exact_hermitian})"
        )

    return CanonicalQuarticPreimageRecord(
        construction_kind=construction_kind,
        state_dimension=len(partial),
        defect_dimension=right.shape[1],
        parameter_scale=format_float(parameter_scale),
        colligation_error=format_float(colligation_error),
        spectral_radius=format_float(spectral_radius),
        first_transfer_norm=format_float(float(np.linalg.norm(first))),
        second_transfer_norm=format_float(float(np.linalg.norm(second))),
        quartic_endpoint_norm=format_float(
            float(np.linalg.norm(upper[4]))
        ),
        positive_target_norm=format_float(float(np.linalg.norm(target))),
        positive_target_minimum_eigenvalue=format_float(target_minimum),
        polynomial_column_norm=format_float(float(np.linalg.norm(column))),
        quartic_preimage_error=format_float(preimage_error),
        perpendicular_column_error=format_float(perpendicular_error),
        lower_endpoint_motion_error=format_float(lower_endpoint_error),
        exact_coboundary_residual_word_count=exact_forcing,
        exact_endpoint_residual_word_count=exact_endpoint,
        exact_witness_hermitian_residual_word_count=exact_hermitian,
        all_checks_passed=verified,
    )


def standard_records() -> list[CanonicalQuarticPreimageRecord]:
    """Return unstructured, rank-changing, and delayed records."""

    records: list[CanonicalQuarticPreimageRecord] = []
    for defect_dimension in range(1, 5):
        for repetition in range(2):
            state_dimension = (
                3 * defect_dimension + 3 + repetition
            )
            partial, right, left = random_partial_isometry(
                state_dimension,
                defect_dimension,
                np.random.default_rng(
                    114_000 + 100 * defect_dimension + repetition
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

    for defect_dimension in range(3, 7):
        for parameter_scale in (1.0, 0.5, 0.2, 0.05):
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

    generator = np.random.default_rng(114_800)
    for defect_dimension in (2, 3, 4):
        partial, right, left = delayed_random_partial_isometry(
            4 * defect_dimension,
            defect_dimension,
            generator,
        )
        records.append(
            audit_case(
                "complete_delay",
                partial,
                right,
                left,
                1,
            )
        )
    return records


def write_records(
    records: list[CanonicalQuarticPreimageRecord],
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
            "repeated_crabb_canonical_quartic_preimage_s70224.jsonl"
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
