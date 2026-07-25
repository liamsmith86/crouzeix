#!/usr/bin/env python3
"""Audit the canonical repair of the grade-one boundary Stein slack.

The grade-one Schur residual of L219's boundary metric collapses to

    K_2 = 2F - Q_2 F - F Q_2,
    Q_2 = (S*)^2 S^2.

Its Stein endpoint is exactly ``2 B_1 B_1*``.  Consequently the
canonical repair ``X_2 = -G_S(K_2)`` preserves the leading lower face
and changes the physical upper face by ``-8 B_1 B_1*``.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_boundary_slack_deflation import (
    slack_schur_coefficients,
)
from repeated_crabb_elliptic_cokernel import stein_inverse
from repeated_crabb_one_image_generator import (
    boundary_metric_coefficient,
)
from repeated_crabb_transfer_channel_covariance import (
    random_partial_isometry,
    transfer_coefficient,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class BoundarySlackRepairRecord:
    """One grade-one compact-residual and repair audit."""

    state_dimension: int
    defect_dimension: int
    compact_residual_error: str
    endpoint_response_error: str
    repair_stein_error: str
    repaired_lower_face_error: str
    repaired_upper_face_error: str
    repaired_upper_maximum_eigenvalue: str
    all_checks_passed: bool


def compact_residual(
    partial: Matrix,
    left: Matrix,
) -> Matrix:
    """Return the closed form for the grade-one slack residual."""

    left_projection = left @ left.conj().T
    second_initial_projection = (
        np.linalg.matrix_power(partial.conj().T, 2)
        @ np.linalg.matrix_power(partial, 2)
    )
    return (
        2 * left_projection
        - second_initial_projection @ left_projection
        - left_projection @ second_initial_projection
    )


def audit_case(
    state_dimension: int,
    defect_dimension: int,
    seed: int,
) -> BoundarySlackRepairRecord:
    """Audit one unstructured pure partial isometry."""

    partial, right, left = random_partial_isometry(
        state_dimension,
        defect_dimension,
        np.random.default_rng(seed),
    )
    transfer = transfer_coefficient(partial, right, left, 1)
    left_gram = transfer @ transfer.conj().T
    right_gram = transfer.conj().T @ transfer

    direct_residual = slack_schur_coefficients(
        partial,
        right,
        left,
        2,
    )[2]
    residual = compact_residual(partial, left)
    residual_error = float(np.linalg.norm(direct_residual - residual))

    response = stein_inverse(partial, residual)
    endpoint = left.conj().T @ response @ left
    endpoint_error = float(np.linalg.norm(endpoint - 2 * left_gram))

    repair = -response
    repair_forcing = repair - partial.conj().T @ repair @ partial
    stein_error = float(np.linalg.norm(repair_forcing + residual))

    boundary = boundary_metric_coefficient(
        partial,
        right,
        left,
        2,
    )
    repaired_metric = boundary + repair
    lower_face = right.conj().T @ repaired_metric @ right
    upper_face = 4 * left.conj().T @ repaired_metric @ left
    lower_error = float(np.linalg.norm(lower_face - right_gram))
    upper_error = float(np.linalg.norm(upper_face + 12 * left_gram))
    upper_maximum = float(
        np.linalg.eigvalsh(
            (upper_face + upper_face.conj().T) / 2
        )[-1]
    )

    tolerance = 3e-8
    verified = bool(
        residual_error < tolerance
        and endpoint_error < tolerance
        and stein_error < tolerance
        and lower_error < tolerance
        and upper_error < tolerance
        and upper_maximum < tolerance
    )
    if not verified:
        raise RuntimeError(
            "the boundary-slack repair audit failed: "
            f"dimension={state_dimension}, defect={defect_dimension}, "
            f"residual={residual_error:.3e}, "
            f"endpoint={endpoint_error:.3e}, "
            f"Stein={stein_error:.3e}, lower={lower_error:.3e}, "
            f"upper={upper_error:.3e}"
        )
    return BoundarySlackRepairRecord(
        state_dimension=state_dimension,
        defect_dimension=defect_dimension,
        compact_residual_error=format_float(residual_error),
        endpoint_response_error=format_float(endpoint_error),
        repair_stein_error=format_float(stein_error),
        repaired_lower_face_error=format_float(lower_error),
        repaired_upper_face_error=format_float(upper_error),
        repaired_upper_maximum_eigenvalue=format_float(upper_maximum),
        all_checks_passed=verified,
    )


def standard_records() -> list[BoundarySlackRepairRecord]:
    """Return deterministic scalar and noncommuting audits."""

    records = []
    for defect_dimension in (1, 2, 3):
        for repetition in range(4):
            state_dimension = 3 * defect_dimension + 2 + repetition
            records.append(
                audit_case(
                    state_dimension,
                    defect_dimension,
                    109_000 + 100 * defect_dimension + repetition,
                )
            )
    return records


def write_records(
    records: list[BoundarySlackRepairRecord],
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
            "repeated_crabb_boundary_slack_repair_s70224.jsonl"
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
