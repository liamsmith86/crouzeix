#!/usr/bin/env python3
"""Audit pointwise and scaled selection for the canonical cubic flag.

The cubic upper-gap obstruction exposed in A172 has an exact
ten-word state forcing.  Cyclic word reduction proves that it pairs
to zero with every reducing-commutant cokernel.  Numerically, the
minimum-norm endpoint preimage on scaled irreducible Schur chains is
``O(scale**2)`` while the target is ``O(scale**3)``.

The range statement is proved in the companion note.  The scaling
test is evidence for, not a proof of, bounded analytic divisibility.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import sqrtm

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
    reduce_word,
    scale,
    schur_residual_series,
)
from repeated_crabb_elliptic_cokernel import (
    endpoint_motion,
    projected_column_basis,
)
from repeated_crabb_elliptic_second_face import (
    hermitian_coordinates,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class CanonicalCubicSelectionRecord:
    """One scaled endpoint-preimage audit."""

    multiplicity: int
    parameter_scale: str
    state_dimension: int
    endpoint_map_rank: int
    endpoint_cokernel_dimension: int
    minimum_nonzero_singular_value: str
    cubic_target_norm: str
    preimage_column_norm: str
    preimage_over_scale_squared: str
    range_equation_error: str
    perpendicular_column_error: str
    exact_cubic_word_count: int
    remaining_cyclic_trace_classes: int
    all_checks_passed: bool


def cubic_forcing_words() -> Polynomial:
    """Return the exact state forcing for the cubic canonical repair."""

    residual = schur_residual_series(3)
    compact = residual[2]
    cubic = residual[3]
    j_operator = multiply(
        add(IDENTITY, F),
        multiply(STAR, add(IDENTITY, E)),
    )
    first_operator = add(
        j_operator,
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


def cyclic_normal_form(word: str) -> str:
    """Return the shortest reduced cyclic representative of one word."""

    if not word:
        return word
    candidates = []
    for index in range(len(word)):
        rotated = word[index:] + word[:index]
        reduced = reduce_word(rotated)
        if len(reduced) == 1 and next(iter(reduced.values())) == 1:
            candidates.append(next(iter(reduced)))
    if not candidates:
        raise RuntimeError(
            f"the cubic trace word did not have a monomial reduction: {word}"
        )
    return min(candidates, key=lambda item: (len(item), item))


def cyclic_trace_classes(
    polynomial: Polynomial,
) -> Polynomial:
    """Collect a word polynomial modulo cyclic trace equivalence."""

    classes: Polynomial = {}
    for word, coefficient in polynomial.items():
        representative = cyclic_normal_form(word)
        classes[representative] = (
            classes.get(representative, 0) + coefficient
        )
    return {
        word: coefficient
        for word, coefficient in classes.items()
        if coefficient
    }


def physical_equality_operator(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
) -> Matrix:
    """Return the physical disk operator associated with ``partial``."""

    identity = np.eye(len(partial), dtype=complex)
    metric = (
        2 * identity
        - right @ right.conj().T
        + 2 * left @ left.conj().T
    )
    root = np.asarray(sqrtm(metric), dtype=complex)
    return np.linalg.inv(root) @ partial @ root


def minimum_endpoint_preimage(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    target: Matrix,
) -> tuple[Matrix, int, float, float]:
    """Return a minimum-coordinate preimage of one Hermitian target."""

    column_basis = projected_column_basis(right)
    endpoint_map = np.stack(
        [
            hermitian_coordinates(
                endpoint_motion(
                    operator,
                    right,
                    left,
                    column,
                )
            )
            for column in column_basis
        ],
        axis=1,
    )
    target_coordinates = hermitian_coordinates(target)
    solution, _, rank, singular_values = np.linalg.lstsq(
        endpoint_map,
        target_coordinates,
        rcond=1e-10,
    )
    column = sum(
        coefficient * basis_column
        for coefficient, basis_column in zip(
            solution,
            column_basis,
            strict=True,
        )
    )
    error = float(
        np.linalg.norm(endpoint_map @ solution - target_coordinates)
    )
    minimum_singular_value = float(singular_values[rank - 1])
    return column, int(rank), minimum_singular_value, error


def audit_case(
    multiplicity: int,
    parameter_scale: float,
    seed: int,
    exact_word_count: int,
    cyclic_class_count: int,
) -> CanonicalCubicSelectionRecord:
    """Audit one cubic endpoint selection."""

    partial, right, left, colligation_error, _ = rank_chain_case(
        multiplicity,
        seed,
        parameter_scale,
    )
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
    target = 4 * endpoint[3]

    operator = physical_equality_operator(partial, right, left)
    column, rank, minimum_singular_value, range_error = (
        minimum_endpoint_preimage(
            operator,
            right,
            left,
            target,
        )
    )
    perpendicular_error = float(
        np.linalg.norm(right.conj().T @ column)
    )
    column_norm = float(np.linalg.norm(column))
    target_norm = float(np.linalg.norm(target))

    tolerance = 3e-8
    verified = bool(
        colligation_error < tolerance
        and constant_error < tolerance
        and rank == multiplicity**2 - 1
        and range_error < tolerance
        and perpendicular_error < tolerance
        and exact_word_count == 10
        and cyclic_class_count == 0
        and column_norm / parameter_scale**2 < 0.03
    )
    if not verified:
        raise RuntimeError(
            "the canonical cubic selection audit failed: "
            f"multiplicity={multiplicity}, scale={parameter_scale}, "
            f"rank={rank}, range={range_error:.3e}, "
            f"column/scale^2={column_norm / parameter_scale**2:.3e}"
        )

    return CanonicalCubicSelectionRecord(
        multiplicity=multiplicity,
        parameter_scale=format_float(parameter_scale),
        state_dimension=len(partial),
        endpoint_map_rank=rank,
        endpoint_cokernel_dimension=multiplicity**2 - rank,
        minimum_nonzero_singular_value=format_float(
            minimum_singular_value
        ),
        cubic_target_norm=format_float(target_norm),
        preimage_column_norm=format_float(column_norm),
        preimage_over_scale_squared=format_float(
            column_norm / parameter_scale**2
        ),
        range_equation_error=format_float(range_error),
        perpendicular_column_error=format_float(perpendicular_error),
        exact_cubic_word_count=exact_word_count,
        remaining_cyclic_trace_classes=cyclic_class_count,
        all_checks_passed=verified,
    )


def standard_records() -> list[CanonicalCubicSelectionRecord]:
    """Return exact-word and scaled numerical selection records."""

    forcing = cubic_forcing_words()
    cyclic_classes = cyclic_trace_classes(forcing)
    records = []
    for multiplicity in (4, 5):
        for parameter_scale in (0.5, 0.2, 0.1, 0.05):
            records.append(
                audit_case(
                    multiplicity,
                    parameter_scale,
                    109_200 + multiplicity,
                    len(forcing),
                    len(cyclic_classes),
                )
            )
    return records


def write_records(
    records: list[CanonicalCubicSelectionRecord],
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
            "repeated_crabb_canonical_cubic_selection_s70224.jsonl"
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
