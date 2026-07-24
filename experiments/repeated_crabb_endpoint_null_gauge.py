#!/usr/bin/env python3
"""Audit endpoint-null gauges in the repeated elliptic selection.

For a balanced pure partial isometry ``S`` with orthogonal defect
frames ``V,W``, the balanced second-frame column ``2W`` produces the
Stein response

    2(VW* + WV*).

Its physical response is ``4(VW* + WV*)`` and both endpoint
compressions vanish.  Adding this column to L207 therefore leaves its
proved endpoint face unchanged.

There is a stronger delayed-face normalization.  If
``B_1=W*S*V=0``, then the balanced metric coefficient

    X_axis = S W W* S* - S* V V* S

is compatible with the complete second elliptic Stein equation.  Its
two endpoint compressions vanish, and on every repeated monomial
shift it is the exact second coefficient of the elliptic-axis metric.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import sqrtm

from crabb_block_hardy_equality import format_float
from repeated_crabb_all_grade_preimage import (
    stable_random_partial_isometry,
)
from repeated_crabb_elliptic_cokernel import (
    endpoint_motion,
    stein_inverse,
)
from repeated_crabb_elliptic_selection import (
    explicit_balanced_column,
    face_data_from_partial_isometry,
    haar_unitary,
)
from repeated_crabb_transfer_deflation import heterogeneous_shift


@dataclass(frozen=True)
class EndpointNullGaugeRecord:
    """One endpoint-null and optional axis-normalization audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    chain_length: int
    spectral_radius: str
    balanced_stein_error: str
    physical_stein_error: str
    lower_endpoint_error: str
    upper_endpoint_error: str
    l207_endpoint_invariance_error: str
    first_transfer_norm: str
    delayed_gauge_applied: bool
    ordered_expansion_error: str
    delayed_row_support_error: str
    delayed_column_formula_error: str
    delayed_reconstruction_error: str
    delayed_lower_endpoint_error: str
    delayed_upper_endpoint_error: str
    axis_second_metric_error: str
    all_checks_passed: bool


def audit_case(
    construction_kind: str,
    balanced: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
    chain_length: int = 0,
) -> EndpointNullGaugeRecord:
    """Audit the null column and, for shifts, the exact axis gauge."""

    data = face_data_from_partial_isometry(balanced, right, left)
    metric = data.metric
    metric_root = np.asarray(sqrtm(metric), dtype=complex)
    metric_root_inverse = np.linalg.inv(metric_root)
    operator = data.operator

    balanced_null_column = 2 * left
    physical_null_column = metric_root @ balanced_null_column
    balanced_null_forcing = (
        right @ balanced_null_column.conj().T
        + balanced_null_column @ right.conj().T
    )
    balanced_null_metric = 2 * (
        right @ left.conj().T + left @ right.conj().T
    )
    balanced_stein_error = float(
        np.linalg.norm(
            balanced_null_metric
            - balanced.conj().T
            @ balanced_null_metric
            @ balanced
            - balanced_null_forcing
        )
    )

    physical_null_forcing = (
        right @ physical_null_column.conj().T
        + physical_null_column @ right.conj().T
    )
    physical_null_metric = 4 * (
        right @ left.conj().T + left @ right.conj().T
    )
    physical_stein_error = float(
        np.linalg.norm(
            physical_null_metric
            - operator.conj().T
            @ physical_null_metric
            @ operator
            - physical_null_forcing
        )
    )
    lower_error = float(
        np.linalg.norm(
            right.conj().T @ physical_null_metric @ right
        )
    )
    upper_error = float(
        np.linalg.norm(
            left.conj().T @ physical_null_metric @ left
        )
    )

    original_balanced_column, _ = explicit_balanced_column(
        balanced,
        right,
        left,
    )
    original_physical_column = metric_root @ original_balanced_column
    gauged_physical_column = (
        metric_root
        @ (original_balanced_column + balanced_null_column)
    )
    endpoint_invariance_error = float(
        np.linalg.norm(
            endpoint_motion(
                operator,
                right,
                left,
                gauged_physical_column,
            )
            - endpoint_motion(
                operator,
                right,
                left,
                original_physical_column,
            )
        )
    )

    first_transfer = (
        left.conj().T @ balanced.conj().T @ right
    )
    delayed_gauge_applied = bool(
        construction_kind
        in {
            "delayed_unstructured_partial_isometry",
            "repeated_monomial_shift",
        }
    )
    ordered_expansion_error = 0.0
    delayed_row_error = 0.0
    delayed_column_formula_error = 0.0
    delayed_reconstruction_error = 0.0
    delayed_lower_error = 0.0
    delayed_upper_error = 0.0
    axis_error = 0.0
    if delayed_gauge_applied:
        right_projection = right @ right.conj().T
        left_projection = left @ left.conj().T
        projection = np.eye(len(balanced)) - right_projection
        expected_axis_metric = (
            balanced
            @ left_projection
            @ balanced.conj().T
            - balanced.conj().T
            @ right_projection
            @ balanced
        )
        balanced_fixed_forcing = (
            metric_root_inverse
            @ data.fixed_forcing
            @ metric_root_inverse
        )
        reduced_fixed_forcing = (
            balanced_fixed_forcing
            - right_projection
            @ balanced_fixed_forcing
            @ right_projection
        )
        predicted_reduced_forcing = (
            -left_projection
            + balanced
            @ left_projection
            @ balanced.conj().T
            - balanced.conj().T
            @ right_projection
            @ balanced
            - 2
            * right_projection
            @ np.linalg.matrix_power(balanced, 4)
            + np.linalg.matrix_power(balanced.conj().T, 2)
            @ right_projection
            @ np.linalg.matrix_power(balanced, 2)
            - 2
            * np.linalg.matrix_power(balanced.conj().T, 4)
            @ right_projection
        )
        ordered_expansion_error = float(
            np.linalg.norm(
                reduced_fixed_forcing
                - predicted_reduced_forcing
            )
        )
        required_row_forcing = (
            expected_axis_metric
            - balanced.conj().T
            @ expected_axis_metric
            @ balanced
            - balanced_fixed_forcing
        )
        delayed_row_error = float(
            np.linalg.norm(
                projection
                @ required_row_forcing
                @ projection
            )
        )
        delayed_balanced_column = (
            projection @ required_row_forcing @ right
            + 0.5
            * right
            @ (
                right.conj().T
                @ required_row_forcing
                @ right
            )
        )
        explicit_delayed_column = (
            -0.5
            * right
            @ (
                right.conj().T
                @ balanced_fixed_forcing
                @ right
            )
            + 2
            * np.linalg.matrix_power(
                balanced.conj().T,
                4,
            )
            @ right
        )
        delayed_column_formula_error = float(
            np.linalg.norm(
                delayed_balanced_column - explicit_delayed_column
            )
        )
        delayed_physical_column = (
            metric_root @ delayed_balanced_column
        )
        delayed_total_forcing = (
            data.fixed_forcing
            + right @ delayed_physical_column.conj().T
            + delayed_physical_column @ right.conj().T
        )
        delayed_physical_metric = stein_inverse(
            operator,
            delayed_total_forcing,
        )
        expected_physical_metric = (
            metric_root @ expected_axis_metric @ metric_root
        )
        delayed_reconstruction_error = float(
            np.linalg.norm(
                delayed_physical_metric - expected_physical_metric
            )
        )
        delayed_lower_error = float(
            np.linalg.norm(
                right.conj().T
                @ delayed_physical_metric
                @ right
            )
        )
        delayed_upper_error = float(
            np.linalg.norm(
                left.conj().T
                @ delayed_physical_metric
                @ left
            )
        )
        if construction_kind == "repeated_monomial_shift":
            axis_error = float(
                np.linalg.norm(
                    metric_root_inverse
                    @ delayed_physical_metric
                    @ metric_root_inverse
                    - expected_axis_metric
                )
            )

    spectral_radius = float(
        np.max(np.abs(np.linalg.eigvals(balanced)))
    )
    tolerance = 2e-9
    verified = bool(
        spectral_radius < 1
        and balanced_stein_error < tolerance
        and physical_stein_error < tolerance
        and lower_error < tolerance
        and upper_error < tolerance
        and endpoint_invariance_error < tolerance
        and (
            not delayed_gauge_applied
            or (
                np.linalg.norm(first_transfer) < tolerance
                and ordered_expansion_error < tolerance
                and delayed_row_error < tolerance
                and delayed_column_formula_error < tolerance
                and delayed_reconstruction_error < tolerance
                and delayed_lower_error < tolerance
                and delayed_upper_error < tolerance
            )
        )
        and axis_error < tolerance
    )
    if not verified:
        raise RuntimeError(
            "the endpoint-null gauge audit failed: "
            f"balanced={balanced_stein_error:.3e}, "
            f"physical={physical_stein_error:.3e}, "
            f"lower={lower_error:.3e}, upper={upper_error:.3e}, "
            f"invariance={endpoint_invariance_error:.3e}, "
            f"ordered={ordered_expansion_error:.3e}, "
            f"delayed-row={delayed_row_error:.3e}, "
            f"delayed-column={delayed_column_formula_error:.3e}, "
            f"delayed-reconstruction={delayed_reconstruction_error:.3e}, "
            f"axis={axis_error:.3e}"
        )

    return EndpointNullGaugeRecord(
        construction_kind=construction_kind,
        state_dimension=len(balanced),
        defect_dimension=right.shape[1],
        chain_length=chain_length,
        spectral_radius=format_float(spectral_radius),
        balanced_stein_error=format_float(balanced_stein_error),
        physical_stein_error=format_float(physical_stein_error),
        lower_endpoint_error=format_float(lower_error),
        upper_endpoint_error=format_float(upper_error),
        l207_endpoint_invariance_error=format_float(
            endpoint_invariance_error
        ),
        first_transfer_norm=format_float(
            float(np.linalg.norm(first_transfer))
        ),
        delayed_gauge_applied=delayed_gauge_applied,
        ordered_expansion_error=format_float(
            ordered_expansion_error
        ),
        delayed_row_support_error=format_float(delayed_row_error),
        delayed_column_formula_error=format_float(
            delayed_column_formula_error
        ),
        delayed_reconstruction_error=format_float(
            delayed_reconstruction_error
        ),
        delayed_lower_endpoint_error=format_float(
            delayed_lower_error
        ),
        delayed_upper_endpoint_error=format_float(
            delayed_upper_error
        ),
        axis_second_metric_error=format_float(axis_error),
        all_checks_passed=verified,
    )


def delayed_random_partial_isometry(
    dimension: int,
    multiplicity: int,
    generator: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return an unstructured stable partial isometry with ``B_1=0``."""

    middle_dimension = dimension - 2 * multiplicity
    if middle_dimension < multiplicity:
        raise ValueError("the delayed construction needs n >= 3m")

    identity = np.eye(dimension, dtype=complex)
    right = identity[:, :multiplicity]
    middle = identity[
        :,
        multiplicity : multiplicity + middle_dimension,
    ]
    left = identity[:, multiplicity + middle_dimension :]
    domain = np.hstack((middle, left))
    range_frame = np.hstack((right, middle))

    while True:
        input_gauge = haar_unitary(middle_dimension, generator)
        output_gauge = haar_unitary(middle_dimension, generator)
        residual_gauge = haar_unitary(
            middle_dimension - multiplicity,
            generator,
        )
        upper_left = (
            input_gauge[:, :multiplicity].conj().T
        )
        lower_right = output_gauge[:, :multiplicity]
        lower_left = (
            output_gauge[:, multiplicity:]
            @ residual_gauge
            @ input_gauge[:, multiplicity:].conj().T
        )
        bridge = np.block(
            [
                [
                    upper_left,
                    np.zeros(
                        (multiplicity, multiplicity),
                        dtype=complex,
                    ),
                ],
                [lower_left, lower_right],
            ]
        )
        balanced = range_frame @ bridge @ domain.conj().T
        spectral_radius = float(
            np.max(np.abs(np.linalg.eigvals(balanced)))
        )
        if spectral_radius < 0.97:
            return balanced, right, left


def standard_records() -> list[EndpointNullGaugeRecord]:
    """Return deterministic unstructured and exact-axis audits."""

    records: list[EndpointNullGaugeRecord] = []
    generator = np.random.default_rng(71324)
    for dimension, multiplicity in ((7, 2), (9, 3), (12, 3), (14, 4)):
        for _ in range(2):
            records.append(
                audit_case(
                    "unstructured_partial_isometry",
                    *stable_random_partial_isometry(
                        dimension,
                        multiplicity,
                        generator,
                    ),
                )
            )

    for dimension, multiplicity in ((8, 2), (10, 2), (12, 3), (15, 3)):
        for _ in range(2):
            records.append(
                audit_case(
                    "delayed_unstructured_partial_isometry",
                    *delayed_random_partial_isometry(
                        dimension,
                        multiplicity,
                        generator,
                    ),
                )
            )

    for length in range(2, 9):
        for multiplicity in (1, 2, 3):
            records.append(
                audit_case(
                    "repeated_monomial_shift",
                    *heterogeneous_shift((length,) * multiplicity),
                    chain_length=length,
                )
            )
    return records


def write_records(
    records: list[EndpointNullGaugeRecord],
    output: Path,
) -> None:
    """Write deterministic JSON Lines records atomically."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)


def parse_args() -> argparse.Namespace:
    """Parse command-line options."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_endpoint_null_gauge_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the endpoint-null gauge audit."""

    args = parse_args()
    records = standard_records()
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
