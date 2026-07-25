#!/usr/bin/env python3
"""Audit the exact boundary-layer metric sandwich.

The candidate one-image construction contains an explicit metric

    P_bl = I - sum_j a_j (S*)^j E S^j
             + sum_j b_j S^j F (S*)^j,

where ``a_j=q^j/(1+q^j)``, ``b_j=q^j``, and ``q=c^2``.  Two orbit
resolutions turn both metric gaps into positive sums.  This script
checks those decompositions on deterministic partial isometries and,
separately, records that the same metric need not have positive Stein
slack for the elliptic pullback.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import sqrtm

from crabb_block_hardy_equality import format_float
from crabb_full_equality_elliptic_merger import ellipse_pullback
from repeated_crabb_transfer_channel_covariance import (
    inflated_case,
    random_partial_isometry,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class BoundaryMetricRecord:
    """One numerical audit of the two positive gap formulas."""

    construction_kind: str
    ellipse_parameter: str
    state_dimension: int
    defect_dimension: int
    spectral_radius: str
    right_orbit_resolution_error: str
    left_orbit_resolution_error: str
    lower_gap_decomposition_error: str
    upper_gap_decomposition_error: str
    lower_physical_gap_minimum: str
    upper_physical_gap_minimum: str
    elliptic_stein_slack_minimum: str
    metric_sandwich_verified: bool


def hermitian_part(matrix: Matrix) -> Matrix:
    """Return the Hermitian part of ``matrix``."""

    return (matrix + matrix.conj().T) / 2


def orbit_data(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    ellipse_parameter: float,
    maximum_terms: int = 2_000,
    tail_tolerance: float = 2e-15,
) -> tuple[Matrix, Matrix, Matrix, Matrix, Matrix]:
    """Return the boundary metric and its four positive orbit sums."""

    dimension = len(partial)
    identity = np.eye(dimension, dtype=complex)
    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    nome = ellipse_parameter**2

    boundary_metric = identity.copy()
    right_resolution = right_projection.copy()
    left_resolution = left_projection.copy()
    lower_positive = left_projection / 4
    upper_positive = 2 * right_projection

    forward = np.eye(dimension, dtype=complex)
    backward = np.eye(dimension, dtype=complex)
    for index in range(1, maximum_terms + 1):
        forward = forward @ partial
        backward = backward @ partial.conj().T
        right_orbit = backward @ right_projection @ forward
        left_orbit = forward @ left_projection @ backward
        right_resolution += right_orbit
        left_resolution += left_orbit

        right_weight = nome**index / (1 + nome**index)
        left_weight = nome**index
        boundary_metric -= right_weight * right_orbit
        boundary_metric += left_weight * left_orbit
        lower_positive += (0.5 - right_weight) * right_orbit
        lower_positive += left_weight * left_orbit
        upper_positive += right_weight * right_orbit
        upper_positive += (1 - left_weight) * left_orbit

        if (
            np.linalg.norm(forward, ord=2) < tail_tolerance
            and np.linalg.norm(backward, ord=2) < tail_tolerance
        ):
            break
    else:
        raise RuntimeError("the pure-partial-isometry orbit did not decay")

    return (
        hermitian_part(boundary_metric),
        hermitian_part(right_resolution),
        hermitian_part(left_resolution),
        hermitian_part(lower_positive),
        hermitian_part(upper_positive),
    )


def audit_case(
    construction_kind: str,
    ellipse_parameter: float,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
) -> BoundaryMetricRecord:
    """Audit one partial isometry."""

    dimension = len(partial)
    identity = np.eye(dimension, dtype=complex)
    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    base_metric = 2 * identity - right_projection + 2 * left_projection
    inverse_base_metric = np.linalg.inv(base_metric)
    metric_root = np.asarray(sqrtm(base_metric), dtype=complex)

    (
        boundary_metric,
        right_resolution,
        left_resolution,
        lower_positive,
        upper_positive,
    ) = orbit_data(
        partial,
        right,
        left,
        ellipse_parameter,
    )
    lower_gap = boundary_metric - inverse_base_metric
    upper_gap = 4 * inverse_base_metric - boundary_metric
    physical_metric = metric_root @ boundary_metric @ metric_root

    physical_operator = (
        np.linalg.inv(metric_root) @ partial @ metric_root
    )
    elliptic_operator = ellipse_pullback(
        physical_operator,
        ellipse_parameter,
        0.0,
    )
    balanced_elliptic_operator = (
        metric_root @ elliptic_operator @ np.linalg.inv(metric_root)
    )
    stein_slack = hermitian_part(
        boundary_metric
        - balanced_elliptic_operator.conj().T
        @ boundary_metric
        @ balanced_elliptic_operator
    )

    right_resolution_error = float(
        np.linalg.norm(right_resolution - identity)
    )
    left_resolution_error = float(
        np.linalg.norm(left_resolution - identity)
    )
    lower_decomposition_error = float(
        np.linalg.norm(lower_gap - lower_positive)
    )
    upper_decomposition_error = float(
        np.linalg.norm(upper_gap - upper_positive)
    )
    lower_minimum = float(
        np.linalg.eigvalsh(
            hermitian_part(physical_metric - identity)
        )[0]
    )
    upper_minimum = float(
        np.linalg.eigvalsh(
            hermitian_part(4 * identity - physical_metric)
        )[0]
    )
    slack_minimum = float(np.linalg.eigvalsh(stein_slack)[0])
    tolerance = 2e-10
    verified = bool(
        right_resolution_error < tolerance
        and left_resolution_error < tolerance
        and lower_decomposition_error < tolerance
        and upper_decomposition_error < tolerance
        and lower_minimum > -tolerance
        and upper_minimum > -tolerance
    )
    if not verified:
        raise RuntimeError(
            "the boundary metric audit failed: "
            f"right={right_resolution_error:.3e}, "
            f"left={left_resolution_error:.3e}, "
            f"lower={lower_decomposition_error:.3e}, "
            f"upper={upper_decomposition_error:.3e}, "
            f"lower_min={lower_minimum:.3e}, "
            f"upper_min={upper_minimum:.3e}"
        )

    return BoundaryMetricRecord(
        construction_kind=construction_kind,
        ellipse_parameter=format_float(ellipse_parameter),
        state_dimension=dimension,
        defect_dimension=right.shape[1],
        spectral_radius=format_float(
            float(np.max(np.abs(np.linalg.eigvals(partial))))
        ),
        right_orbit_resolution_error=format_float(
            right_resolution_error
        ),
        left_orbit_resolution_error=format_float(
            left_resolution_error
        ),
        lower_gap_decomposition_error=format_float(
            lower_decomposition_error
        ),
        upper_gap_decomposition_error=format_float(
            upper_decomposition_error
        ),
        lower_physical_gap_minimum=format_float(lower_minimum),
        upper_physical_gap_minimum=format_float(upper_minimum),
        elliptic_stein_slack_minimum=format_float(slack_minimum),
        metric_sandwich_verified=verified,
    )


def standard_records() -> list[BoundaryMetricRecord]:
    """Return deterministic general and fully delayed audits."""

    records: list[BoundaryMetricRecord] = []
    for repetition in range(3):
        partial, right, left = random_partial_isometry(
            7 + repetition,
            2,
            np.random.default_rng(104_100 + repetition),
        )
        for ellipse_parameter in (0.04, 0.16):
            records.append(
                audit_case(
                    "unstructured_partial_isometry",
                    ellipse_parameter,
                    partial,
                    right,
                    left,
                )
            )

    for grade in range(2, 6):
        partial, right, left, _ = inflated_case(
            7,
            2,
            grade,
            2,
            104_200 + grade,
        )
        for ellipse_parameter in (0.06, 0.18):
            records.append(
                audit_case(
                    f"fully_delayed_grade_{grade}",
                    ellipse_parameter,
                    partial,
                    right,
                    left,
                )
            )
    return records


def write_records(
    records: list[BoundaryMetricRecord],
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
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_boundary_metric_sandwich_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the complete audit."""

    args = parse_args()
    records = standard_records()
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
