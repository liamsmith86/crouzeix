#!/usr/bin/env python3
"""Audit the prepared fourth-order face on the delayed ``B_1=0`` stratum.

L213 fixes the axis-compatible second metric gauge.  In that gauge,
the fourth-order upper Schur coefficient has one additional explicit
polynomial preparation:

    C_hat[prep] = Q (
        4 S^2 W B_2 + 2 S* V B_3* B_2
    ).

This checker builds the complete Riemann/Stein jet through order four,
including both endpoint Schur complements.  It verifies that the
preparation changes the raw face to

    12 B_2 B_2* - 28 channel(B_2* B_2).

Adding L212's grade-two channel coboundary then leaves the coercive
endpoint ``-16 B_2 B_2*``.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import null_space, sqrtm

from crabb_block_hardy_equality import format_float
from repeated_crabb_elliptic_cokernel import (
    endpoint_motion,
    stein_inverse,
)
from repeated_crabb_endpoint_null_gauge import (
    delayed_random_partial_isometry,
)
from repeated_crabb_elliptic_selection import haar_unitary
from repeated_crabb_transfer_channel_covariance import (
    inflated_case,
    transfer_channel,
    transfer_coefficient,
)
from repeated_crabb_transfer_deflation import heterogeneous_shift


@dataclass(frozen=True)
class GradeTwoFaceRecord:
    """One complete fourth-order delayed-face audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    spectral_radius: str
    first_transfer_norm: str
    second_transfer_norm: str
    canonical_first_column_error: str
    delayed_second_metric_error: str
    maximum_lower_schur_error: str
    maximum_pre_fourth_upper_schur_error: str
    raw_endpoint_identity_error: str
    prepared_endpoint_identity_error: str
    final_endpoint_identity_error: str
    apex_total_correction_norm: str
    all_checks_passed: bool


def series_multiply(
    left: list[np.ndarray],
    right: list[np.ndarray],
    degree: int,
) -> list[np.ndarray]:
    """Multiply two square matrix series through ``degree``."""

    dimension = len(left[0])
    result = [
        np.zeros((dimension, dimension), dtype=complex)
        for _ in range(degree + 1)
    ]
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            if left_degree + right_degree <= degree:
                result[left_degree + right_degree] += (
                    left_coefficient @ right_coefficient
                )
    return result


def series_power(
    series: list[np.ndarray],
    power: int,
    degree: int,
) -> list[np.ndarray]:
    """Raise a square matrix series to a nonnegative integer power."""

    dimension = len(series[0])
    result = [np.eye(dimension, dtype=complex)] + [
        np.zeros((dimension, dimension), dtype=complex)
        for _ in range(degree)
    ]
    for _ in range(power):
        result = series_multiply(result, series, degree)
    return result


def ellipse_operator_jets(operator: np.ndarray) -> list[np.ndarray]:
    """Return the exact direct ellipse pullback through order four.

    The centered direct map is

        phi_c(w)
          = w - c w^3 + c^2(2w+w^5)
            - c^3(3w^3+w^7)
            + c^4(w+5w^5+w^9) + O(c^5).
    """

    degree = 4
    pencil = [operator, operator.conj().T] + [
        np.zeros_like(operator) for _ in range(degree - 1)
    ]
    scalar_coefficients = (
        {1: 1},
        {3: -1},
        {1: 2, 5: 1},
        {3: -3, 7: -1},
        {1: 1, 5: 5, 9: 1},
    )
    result = [np.zeros_like(operator) for _ in range(degree + 1)]
    for ellipse_degree, polynomial in enumerate(scalar_coefficients):
        for power, scalar in polynomial.items():
            power_series = series_power(
                pencil,
                power,
                degree - ellipse_degree,
            )
            for pencil_degree, coefficient in enumerate(power_series):
                total_degree = ellipse_degree + pencil_degree
                if total_degree <= degree:
                    result[total_degree] += scalar * coefficient
    return result


def inverse_series(series: list[np.ndarray]) -> list[np.ndarray]:
    """Invert a square matrix series with invertible constant term."""

    result = [np.linalg.inv(series[0])]
    for degree in range(1, len(series)):
        convolution = np.zeros_like(series[0])
        for positive_degree in range(1, degree + 1):
            convolution += (
                series[positive_degree]
                @ result[degree - positive_degree]
            )
        result.append(-result[0] @ convolution)
    return result


def schur_cross_coefficient(
    metric_jets: list[np.ndarray | None],
    endpoint: np.ndarray,
    endpoint_value: float,
    degree: int,
) -> np.ndarray:
    """Return ``[c^degree] B(c)D(c)^(-1)B(c)*`` at one endpoint."""

    complement = null_space(endpoint.conj().T)
    zero = np.zeros_like(metric_jets[0])
    cross: list[np.ndarray] = []
    denominator: list[np.ndarray] = []
    for index in range(degree + 1):
        metric = (
            metric_jets[index]
            if index < len(metric_jets)
            and metric_jets[index] is not None
            else zero
        )
        cross.append(endpoint.conj().T @ metric @ complement)
        denominator.append(
            complement.conj().T @ metric @ complement
        )
    denominator[0] -= endpoint_value * np.eye(complement.shape[1])
    denominator_inverse = inverse_series(denominator)
    left_product = rectangular_series_multiply(
        cross,
        denominator_inverse,
        degree,
    )
    adjoint_cross = [
        coefficient.conj().T for coefficient in cross
    ]
    return rectangular_series_multiply(
        left_product,
        adjoint_cross,
        degree,
    )[degree]


def rectangular_series_multiply(
    left: list[np.ndarray],
    right: list[np.ndarray],
    degree: int,
) -> list[np.ndarray]:
    """Multiply compatible rectangular matrix series."""

    result = [
        np.zeros(
            (left[0].shape[0], right[0].shape[1]),
            dtype=complex,
        )
        for _ in range(degree + 1)
    ]
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            if left_degree + right_degree <= degree:
                result[left_degree + right_degree] += (
                    left_coefficient @ right_coefficient
                )
    return result


def axis_perpendicular_column(
    partial: np.ndarray,
    right: np.ndarray,
    degree: int,
) -> np.ndarray:
    """Return the delayed zero-reflection frame coefficient.

    Through order four these are the perpendicular coefficients of

        (1+2c^2)^(-1)
        (V + 2 sum_{j>=1} (-c)^j (S*)^(2j)V).

    The parallel component is fixed separately by the lower Schur
    condition.
    """

    result = np.zeros_like(right)
    for orbit_degree in range(1, degree + 1):
        remainder = degree - orbit_degree
        if remainder % 2:
            continue
        scalar = (
            2
            * (-1) ** orbit_degree
            * (-2) ** (remainder // 2)
        )
        result += (
            scalar
            * np.linalg.matrix_power(
                partial.conj().T,
                2 * orbit_degree,
            )
            @ right
        )
    return result


def fourth_order_jet(
    operator: np.ndarray,
    metric: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
    partial: np.ndarray,
    metric_root: np.ndarray,
) -> tuple[
    list[np.ndarray],
    list[np.ndarray],
    list[np.ndarray],
    np.ndarray,
]:
    """Construct the lower-tight metric and raw upper fourth face."""

    degree = 4
    operator_jets = ellipse_operator_jets(operator)
    metric_jets: list[np.ndarray | None] = [metric] + [None] * degree
    defect_jets: list[np.ndarray | None] = [right] + [None] * degree

    for index in range(1, degree + 1):
        fixed = np.zeros_like(operator)
        for left_degree in range(index + 1):
            for metric_degree in range(index + 1 - left_degree):
                right_degree = index - left_degree - metric_degree
                if metric_degree >= index:
                    continue
                fixed += (
                    operator_jets[left_degree].conj().T
                    @ metric_jets[metric_degree]
                    @ operator_jets[right_degree]
                )
        for left_degree in range(1, index):
            fixed += (
                defect_jets[left_degree]
                @ defect_jets[index - left_degree].conj().T
            )

        lower_target = schur_cross_coefficient(
            metric_jets,
            right,
            1,
            index,
        )
        parallel = (
            lower_target - right.conj().T @ fixed @ right
        ) / 2
        defect_jets[index] = (
            right @ parallel
            + metric_root
            @ axis_perpendicular_column(partial, right, index)
        )
        forcing = (
            fixed
            + right @ defect_jets[index].conj().T
            + defect_jets[index] @ right.conj().T
        )
        metric_jets[index] = stein_inverse(operator, forcing)

    raw_upper = (
        left.conj().T @ metric_jets[degree] @ left
        - schur_cross_coefficient(
            metric_jets,
            left,
            4,
            degree,
        )
    )
    return (
        operator_jets,
        list(metric_jets),
        list(defect_jets),
        raw_upper,
    )


def gauge_case(
    partial: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
    generator: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply independent state and endpoint unitary gauges."""

    state_dimension = len(partial)
    multiplicity = right.shape[1]
    state_unitary = haar_unitary(state_dimension, generator)
    right_unitary = haar_unitary(multiplicity, generator)
    left_unitary = haar_unitary(multiplicity, generator)
    return (
        state_unitary @ partial @ state_unitary.conj().T,
        state_unitary @ right @ right_unitary,
        state_unitary @ left @ left_unitary,
    )


def audit_case(
    construction_kind: str,
    partial: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
) -> GradeTwoFaceRecord:
    """Audit one delayed partial-isometry colligation."""

    dimension = len(partial)
    identity = np.eye(dimension, dtype=complex)
    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    complement = identity - right_projection
    metric = 2 * identity - right_projection + 2 * left_projection
    metric_root = np.asarray(sqrtm(metric), dtype=complex)
    metric_root_inverse = np.linalg.inv(metric_root)
    operator = metric_root_inverse @ partial @ metric_root

    (
        operator_jets,
        metric_jets,
        defect_jets,
        raw_upper,
    ) = fourth_order_jet(
        operator,
        metric,
        right,
        left,
        partial,
        metric_root,
    )

    first = transfer_coefficient(partial, right, left, 1)
    second = transfer_coefficient(partial, right, left, 2)
    third = transfer_coefficient(partial, right, left, 3)
    left_gram = second @ second.conj().T
    channel_gram = transfer_channel(
        partial,
        right,
        left,
        second.conj().T @ second,
    )

    expected_first_column = (
        metric_root
        @ (
            -2
            * np.linalg.matrix_power(partial.conj().T, 2)
            @ right
        )
    )
    first_column_error = float(
        np.linalg.norm(defect_jets[1] - expected_first_column)
    )

    balanced_second_metric = (
        metric_root_inverse
        @ metric_jets[2]
        @ metric_root_inverse
    )
    expected_second_metric = (
        partial
        @ left_projection
        @ partial.conj().T
        - partial.conj().T
        @ right_projection
        @ partial
    )
    second_metric_error = float(
        np.linalg.norm(
            balanced_second_metric - expected_second_metric
        )
    )

    lower_errors = []
    upper_lower_errors = []
    for degree in range(1, 5):
        lower_schur = (
            right.conj().T @ metric_jets[degree] @ right
            - schur_cross_coefficient(
                metric_jets,
                right,
                1,
                degree,
            )
        )
        lower_errors.append(float(np.linalg.norm(lower_schur)))
        if degree < 4:
            upper_schur = (
                left.conj().T @ metric_jets[degree] @ left
                - schur_cross_coefficient(
                    metric_jets,
                    left,
                    4,
                    degree,
                )
            )
            upper_lower_errors.append(
                float(np.linalg.norm(upper_schur))
            )

    tail_column = (
        complement
        @ partial.conj().T
        @ right
        @ third.conj().T
        @ second
    )
    tail_endpoint = endpoint_motion(
        operator,
        right,
        left,
        metric_root @ tail_column,
    )
    raw_expected = (
        4 * channel_gram
        - 20 * left_gram
        - 2 * tail_endpoint
    )
    raw_error = float(np.linalg.norm(raw_upper - raw_expected))

    preparation_column = complement @ (
        4
        * np.linalg.matrix_power(partial, 2)
        @ left
        @ second
        + 2
        * partial.conj().T
        @ right
        @ third.conj().T
        @ second
    )
    preparation_endpoint = endpoint_motion(
        operator,
        right,
        left,
        metric_root @ preparation_column,
    )
    prepared_upper = raw_upper + preparation_endpoint
    prepared_expected = 12 * left_gram - 28 * channel_gram
    prepared_error = float(
        np.linalg.norm(prepared_upper - prepared_expected)
    )

    l212_column = (
        -3.5
        * complement
        @ np.linalg.matrix_power(partial, 2)
        @ left
        @ second
    )
    l212_endpoint = endpoint_motion(
        operator,
        right,
        left,
        metric_root @ l212_column,
    )
    final_upper = prepared_upper + l212_endpoint
    final_expected = -16 * left_gram
    final_error = float(np.linalg.norm(final_upper - final_expected))

    apex_correction_norm = 0.0
    if construction_kind == "repeated_grade_two_apex":
        apex_correction_norm = float(
            np.linalg.norm(preparation_column + l212_column)
        )

    spectral_radius = float(
        np.max(np.abs(np.linalg.eigvals(partial)))
    )
    tolerance = 3e-9
    verified = bool(
        spectral_radius < 1
        and np.linalg.norm(first) < tolerance
        and first_column_error < tolerance
        and second_metric_error < tolerance
        and max(lower_errors) < tolerance
        and max(upper_lower_errors) < tolerance
        and raw_error < tolerance
        and prepared_error < tolerance
        and final_error < tolerance
        and apex_correction_norm < tolerance
    )
    if not verified:
        raise RuntimeError(
            "the grade-two face audit failed: "
            f"B1={np.linalg.norm(first):.3e}, "
            f"C1={first_column_error:.3e}, X2={second_metric_error:.3e}, "
            f"lower={max(lower_errors):.3e}, "
            f"upper<4={max(upper_lower_errors):.3e}, "
            f"raw={raw_error:.3e}, prepared={prepared_error:.3e}, "
            f"final={final_error:.3e}, apex={apex_correction_norm:.3e}"
        )

    return GradeTwoFaceRecord(
        construction_kind=construction_kind,
        state_dimension=dimension,
        defect_dimension=right.shape[1],
        spectral_radius=format_float(spectral_radius),
        first_transfer_norm=format_float(float(np.linalg.norm(first))),
        second_transfer_norm=format_float(float(np.linalg.norm(second))),
        canonical_first_column_error=format_float(first_column_error),
        delayed_second_metric_error=format_float(second_metric_error),
        maximum_lower_schur_error=format_float(max(lower_errors)),
        maximum_pre_fourth_upper_schur_error=format_float(
            max(upper_lower_errors)
        ),
        raw_endpoint_identity_error=format_float(raw_error),
        prepared_endpoint_identity_error=format_float(prepared_error),
        final_endpoint_identity_error=format_float(final_error),
        apex_total_correction_norm=format_float(apex_correction_norm),
        all_checks_passed=verified,
    )


def standard_records() -> list[GradeTwoFaceRecord]:
    """Return deterministic unstructured, delayed, and apex audits."""

    records: list[GradeTwoFaceRecord] = []
    generator = np.random.default_rng(71424)

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

    for index, (dimension, multiplicity) in enumerate(
        ((8, 2), (10, 3), (12, 3), (14, 4))
    ):
        records.append(
            audit_case(
                "inflated_unstructured_partial_isometry",
                *inflated_case(
                    dimension,
                    multiplicity,
                    2,
                    multiplicity,
                    71_500 + index,
                )[:3],
            )
        )

    for index, lengths in enumerate(
        ((2, 3), (2, 4, 5), (2, 2, 5), (3, 4, 6))
    ):
        records.append(
            audit_case(
                "gauged_heterogeneous_shift",
                *gauge_case(
                    *heterogeneous_shift(lengths),
                    generator,
                ),
            )
        )

    for multiplicity in (1, 2, 3):
        records.append(
            audit_case(
                "repeated_grade_two_apex",
                *gauge_case(
                    *heterogeneous_shift((2,) * multiplicity),
                    generator,
                ),
            )
        )
    return records


def write_records(
    records: list[GradeTwoFaceRecord],
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
            "experiments/repeated_crabb_grade_two_face_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the complete fourth-order audit."""

    args = parse_args()
    records = standard_records()
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
