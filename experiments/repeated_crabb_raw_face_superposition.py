#!/usr/bin/env python3
"""Falsify naive superposition of all raw delayed-face columns.

The one-image gauge and L212 close each *first active* delayed face.
It is tempting to insert every such even-order column into a single
frame series and expect the full endpoint to be
``-16 sum c^(2k) B_k B_k*``.  This checker shows that the tempting
identity already fails at order three on generic noncommuting
colligations.  Later faces must be Schur-orthogonalized first.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_delayed_jet import (
    construct_delayed_metric_jet,
    endpoint_coefficient,
)
from repeated_crabb_one_image_generator import canonical_adjustment
from repeated_crabb_transfer_channel_covariance import (
    random_partial_isometry,
    transfer_coefficient,
)


@dataclass(frozen=True)
class RawFaceSuperpositionRecord:
    """One robust counterexample to the raw superposition identity."""

    state_dimension: int
    defect_dimension: int
    leading_face_error: str
    first_mismatch_order: int
    first_mismatch_norm: str
    maximum_lower_endpoint_error: str
    all_checks_passed: bool


def naive_all_face_adjustment(maximum_degree: int):
    """Return the candidate gauge plus every unflagged face column."""

    one_image = canonical_adjustment(maximum_degree)

    def adjustment(
        order: int,
        partial: np.ndarray,
        right: np.ndarray,
        left: np.ndarray,
    ) -> np.ndarray:
        result = one_image(order, partial, right, left)
        if order % 2:
            return result

        grade = order // 2
        identity = np.eye(len(partial), dtype=complex)
        complement = identity - right @ right.conj().T
        coefficients = [
            transfer_coefficient(
                partial,
                right,
                left,
                index,
            )
            for index in range(grade + 1)
        ]
        column = (
            np.linalg.matrix_power(partial, grade)
            @ left
            @ coefficients[grade]
        )
        for earlier in range(1, grade):
            column += (
                np.linalg.matrix_power(
                    partial.conj().T,
                    grade - earlier,
                )
                @ right
                @ coefficients[grade].conj().T
                @ coefficients[earlier]
            )
        result += -3.5 * complement @ column
        if grade % 2 == 0:
            result += (
                4
                * complement
                @ np.linalg.matrix_power(partial, grade)
                @ left
                @ coefficients[grade]
            )
        return result

    return adjustment


def make_record(
    state_dimension: int,
    defect_dimension: int,
    seed: int,
    maximum_degree: int,
) -> RawFaceSuperpositionRecord:
    """Return one generic noncommuting falsification record."""

    partial, right, left = random_partial_isometry(
        state_dimension,
        defect_dimension,
        np.random.default_rng(seed),
    )
    jet = construct_delayed_metric_jet(
        partial,
        right,
        left,
        maximum_degree,
        naive_all_face_adjustment(maximum_degree),
    )
    upper = [
        endpoint_coefficient(jet, order, upper=True)
        for order in range(1, maximum_degree + 1)
    ]
    lower = [
        endpoint_coefficient(jet, order, upper=False)
        for order in range(1, maximum_degree + 1)
    ]

    first_coefficient = transfer_coefficient(
        partial,
        right,
        left,
        1,
    )
    leading_target = (
        -16 * first_coefficient @ first_coefficient.conj().T
    )
    leading_error = float(np.linalg.norm(upper[1] - leading_target))

    first_mismatch_order = 0
    first_mismatch_norm = 0.0
    for order, endpoint in enumerate(upper, start=1):
        if order % 2:
            target = np.zeros_like(endpoint)
        else:
            grade = order // 2
            coefficient = transfer_coefficient(
                partial,
                right,
                left,
                grade,
            )
            target = -16 * coefficient @ coefficient.conj().T
        error = float(np.linalg.norm(endpoint - target))
        if order > 2 and error > 1e-3:
            first_mismatch_order = order
            first_mismatch_norm = error
            break

    maximum_lower_error = max(
        float(np.linalg.norm(endpoint))
        for endpoint in lower
    )
    verified = bool(
        leading_error < 2e-8
        and first_mismatch_order == 3
        and first_mismatch_norm > 1e-3
        and maximum_lower_error < 2e-8
    )
    if not verified:
        raise RuntimeError(
            "raw face superposition falsification was not reproduced: "
            f"leading={leading_error:.3e}, "
            f"order={first_mismatch_order}, "
            f"mismatch={first_mismatch_norm:.3e}, "
            f"lower={maximum_lower_error:.3e}"
        )
    return RawFaceSuperpositionRecord(
        state_dimension=state_dimension,
        defect_dimension=defect_dimension,
        leading_face_error=format_float(leading_error),
        first_mismatch_order=first_mismatch_order,
        first_mismatch_norm=format_float(first_mismatch_norm),
        maximum_lower_endpoint_error=format_float(
            maximum_lower_error
        ),
        all_checks_passed=verified,
    )


def standard_records(
    maximum_degree: int,
) -> list[RawFaceSuperpositionRecord]:
    """Return deterministic generic counterexamples."""

    return [
        make_record(7, 2, 803_101, maximum_degree),
        make_record(8, 3, 803_102, maximum_degree),
    ]


def write_records(
    records: list[RawFaceSuperpositionRecord],
    output: Path,
) -> None:
    """Write deterministic JSON Lines atomically."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-degree", type=int, default=6)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_raw_face_superposition_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the complete falsification audit."""

    args = parse_args()
    if not 4 <= args.maximum_degree <= 10:
        raise ValueError("maximum-degree must lie between four and ten")
    records = standard_records(args.maximum_degree)
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
