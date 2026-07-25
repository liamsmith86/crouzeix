#!/usr/bin/env python3
"""Audit the canonical zero/one-image gauge through transfer grade six.

The exact scalar axis has a normalized ``nd`` Fourier frame.  Before
two reflected legs can occur, its lift to a partial-isometry
colligation has two pieces:

* forward even orbits ``(S*)^(2j)V``;
* one terminal image
  ``S^a F (S*)^(a+2j)V = S^a W B_(a+2j)``.

On ``B_1=...=B_(k-1)=0`` this gauge should reproduce the explicit
boundary-layer metric through order ``2k-1``.  Its order-``2k`` upper
face should be the universal prepared base, apart from one even-grade
endpoint column.  This checker tests those identities at grades one
through six; it is evidence for the all-grade reduction, not its
ordered-algebra proof.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_delayed_jet import (
    axis_perpendicular_column,
    construct_delayed_metric_jet,
    endpoint_coefficient,
)
from repeated_crabb_elliptic_cokernel import endpoint_motion
from repeated_crabb_grade_two_face import gauge_case
from repeated_crabb_transfer_channel_covariance import (
    inflated_case,
    random_partial_isometry,
    transfer_channel,
    transfer_coefficient,
)
from repeated_crabb_transfer_deflation import heterogeneous_shift


Matrix = np.ndarray


@dataclass(frozen=True)
class OneImageRecord:
    """One leading-face audit in the canonical image gauge."""

    construction_kind: str
    grade: int
    state_dimension: int
    defect_dimension: int
    spectral_radius: str
    maximum_earlier_transfer_norm: str
    active_transfer_norm: str
    maximum_lower_schur_error: str
    maximum_preface_upper_schur_error: str
    maximum_preface_metric_error: str
    raw_face_error: str
    prepared_face_error: str
    final_face_error: str
    all_checks_passed: bool


def inverse_theta_three_coefficients(degree: int) -> tuple[int, ...]:
    """Return ``1/theta_3(c^2)`` through ``degree`` exactly."""

    theta = [0] * (degree + 1)
    theta[0] = 1
    index = 1
    while 2 * index * index <= degree:
        theta[2 * index * index] = 2
        index += 1

    inverse = [0] * (degree + 1)
    inverse[0] = 1
    for order in range(1, degree + 1):
        inverse[order] = -sum(
            theta[positive] * inverse[order - positive]
            for positive in range(1, order + 1)
        )
    return tuple(inverse)


def divided_shifted_coefficient(
    inverse_theta: tuple[int, ...],
    order: int,
    shift: int,
    denominator_step: int,
) -> int:
    """Return ``[c^order] c^shift/(theta_3(c^2)(1+c^step))``."""

    coefficient = 0
    denominator_power = 0
    while shift + denominator_power * denominator_step <= order:
        remainder = order - shift - denominator_power * denominator_step
        coefficient += (-1) ** denominator_power * inverse_theta[remainder]
        denominator_power += 1
    return coefficient


def canonical_perpendicular_column(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    order: int,
    maximum_degree: int,
) -> Matrix:
    """Return the zero/one-image balanced frame coefficient.

    The forward part is the exact normalized ``nd`` Fourier series

        2/theta_3(c^2) sum_j
          (-c)^j/(1+c^(4j)) (S*)^(2j)V.

    The image with delay ``a`` and forward frequency ``j`` is

        2 (-1)^(a+j) c^(2a+j)
          / (theta_3(c^2)(1+c^(4(a+j))))
          Q S^a F (S*)^(a+2j)V.

    Only finitely many terms can contribute to ``maximum_degree``.
    """

    inverse_theta = inverse_theta_three_coefficients(maximum_degree)
    result = np.zeros_like(right)

    for frequency in range(1, maximum_degree + 1):
        scalar = (
            2
            * (-1) ** frequency
            * divided_shifted_coefficient(
                inverse_theta,
                order,
                frequency,
                4 * frequency,
            )
        )
        if scalar:
            result += (
                scalar
                * np.linalg.matrix_power(
                    partial.conj().T,
                    2 * frequency,
                )
                @ right
            )

    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    complement = np.eye(len(partial), dtype=complex) - right_projection
    for delay in range(1, maximum_degree + 1):
        for frequency in range(1, maximum_degree + 1):
            shift = 2 * delay + frequency
            if shift > order:
                break
            scalar = (
                2
                * (-1) ** (delay + frequency)
                * divided_shifted_coefficient(
                    inverse_theta,
                    order,
                    shift,
                    4 * (delay + frequency),
                )
            )
            if scalar:
                result += (
                    scalar
                    * complement
                    @ np.linalg.matrix_power(partial, delay)
                    @ left_projection
                    @ np.linalg.matrix_power(
                        partial.conj().T,
                        delay + 2 * frequency,
                    )
                    @ right
                )
    return result


def canonical_adjustment(maximum_degree: int):
    """Return the shared-engine adjustment for the image gauge."""

    def adjustment(
        order: int,
        partial: Matrix,
        right: Matrix,
        left: Matrix,
    ) -> Matrix:
        return canonical_perpendicular_column(
            partial,
            right,
            left,
            order,
            maximum_degree,
        ) - axis_perpendicular_column(partial, right, order)

    return adjustment


def boundary_metric_coefficient(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    order: int,
) -> Matrix:
    """Return the explicit preface boundary-layer metric coefficient."""

    if order % 2:
        return np.zeros_like(partial)

    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    half_order = order // 2
    result = (
        np.linalg.matrix_power(partial, half_order)
        @ left_projection
        @ np.linalg.matrix_power(partial.conj().T, half_order)
    )
    for divisor in range(1, half_order + 1):
        if half_order % divisor:
            continue
        result += (
            (-1) ** (half_order // divisor)
            * np.linalg.matrix_power(partial.conj().T, divisor)
            @ right_projection
            @ np.linalg.matrix_power(partial, divisor)
        )
    return result


def audit_case(
    construction_kind: str,
    grade: int,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    maximum_degree: int,
) -> OneImageRecord:
    """Audit one fully delayed colligation at its first active face."""

    face_order = 2 * grade
    jet = construct_delayed_metric_jet(
        partial,
        right,
        left,
        maximum_degree,
        canonical_adjustment(maximum_degree),
    )
    coefficients = [
        transfer_coefficient(partial, right, left, index)
        for index in range(maximum_degree + grade + 2)
    ]
    active = coefficients[grade]
    left_gram = active @ active.conj().T
    channel = transfer_channel(
        partial,
        right,
        left,
        active.conj().T @ active,
    )

    metric_root_inverse = np.linalg.inv(jet.metric_root)
    metric_errors = []
    for order in range(1, face_order):
        balanced_metric = (
            metric_root_inverse
            @ jet.metric_coefficients[order]
            @ metric_root_inverse
        )
        metric_errors.append(
            float(
                np.linalg.norm(
                    balanced_metric
                    - boundary_metric_coefficient(
                        partial,
                        right,
                        left,
                        order,
                    )
                )
            )
        )

    lower_errors = [
        float(
            np.linalg.norm(
                endpoint_coefficient(jet, order, upper=False)
            )
        )
        for order in range(1, face_order + 1)
    ]
    earlier_upper_errors = [
        float(
            np.linalg.norm(
                endpoint_coefficient(jet, order, upper=True)
            )
        )
        for order in range(1, face_order)
    ]

    complement = np.eye(len(partial), dtype=complex) - right @ right.conj().T
    channel_column = (
        complement
        @ np.linalg.matrix_power(partial, grade)
        @ left
        @ active
    )
    channel_endpoint = endpoint_motion(
        jet.operator,
        right,
        left,
        jet.metric_root @ channel_column,
    )

    base = 12 * left_gram - 28 * channel
    raw = endpoint_coefficient(jet, face_order, upper=True)
    raw_expected = (
        base - 4 * channel_endpoint
        if grade % 2 == 0
        else base
    )
    raw_error = float(np.linalg.norm(raw - raw_expected))

    prepared = raw
    if grade % 2 == 0:
        prepared = prepared + 4 * channel_endpoint
    prepared_error = float(np.linalg.norm(prepared - base))

    l212_column = -3.5 * channel_column
    final = prepared + endpoint_motion(
        jet.operator,
        right,
        left,
        jet.metric_root @ l212_column,
    )
    final_error = float(np.linalg.norm(final + 16 * left_gram))

    earlier_transfer_norm = max(
        (
            float(np.linalg.norm(coefficients[index]))
            for index in range(1, grade)
        ),
        default=0.0,
    )
    spectral_radius = float(
        np.max(np.abs(np.linalg.eigvals(partial)))
    )
    tolerance = 2e-8
    verified = bool(
        spectral_radius < 1
        and earlier_transfer_norm < tolerance
        and max(lower_errors) < tolerance
        and max(earlier_upper_errors, default=0.0) < tolerance
        and max(metric_errors, default=0.0) < tolerance
        and raw_error < tolerance
        and prepared_error < tolerance
        and final_error < tolerance
    )
    if not verified:
        raise RuntimeError(
            "the one-image audit failed: "
            f"grade={grade}, earlier={earlier_transfer_norm:.3e}, "
            f"lower={max(lower_errors):.3e}, "
            f"upper={max(earlier_upper_errors, default=0.0):.3e}, "
            f"metric={max(metric_errors, default=0.0):.3e}, "
            f"raw={raw_error:.3e}, prepared={prepared_error:.3e}, "
            f"final={final_error:.3e}"
        )

    return OneImageRecord(
        construction_kind=construction_kind,
        grade=grade,
        state_dimension=len(partial),
        defect_dimension=right.shape[1],
        spectral_radius=format_float(spectral_radius),
        maximum_earlier_transfer_norm=format_float(
            earlier_transfer_norm
        ),
        active_transfer_norm=format_float(
            float(np.linalg.norm(active))
        ),
        maximum_lower_schur_error=format_float(max(lower_errors)),
        maximum_preface_upper_schur_error=format_float(
            max(earlier_upper_errors, default=0.0)
        ),
        maximum_preface_metric_error=format_float(
            max(metric_errors, default=0.0)
        ),
        raw_face_error=format_float(raw_error),
        prepared_face_error=format_float(prepared_error),
        final_face_error=format_float(final_error),
        all_checks_passed=verified,
    )


def standard_records(maximum_grade: int) -> list[OneImageRecord]:
    """Return unstructured, rank-changing, and apex audits."""

    maximum_degree = 2 * maximum_grade
    records: list[OneImageRecord] = []
    generator = np.random.default_rng(103_024)

    for grade in range(1, maximum_grade + 1):
        if grade == 1:
            for repetition in range(2):
                partial, right, left = random_partial_isometry(
                    7 + repetition,
                    2,
                    np.random.default_rng(103_100 + repetition),
                )
                records.append(
                    audit_case(
                        "unstructured_partial_isometry",
                        grade,
                        partial,
                        right,
                        left,
                        maximum_degree,
                    )
                )
        else:
            for multiplicity in (1, 2):
                partial, right, left, _ = inflated_case(
                    6 + multiplicity,
                    multiplicity,
                    grade,
                    multiplicity,
                    103_200 + 10 * grade + multiplicity,
                )
                records.append(
                    audit_case(
                        "inflated_unstructured_partial_isometry",
                        grade,
                        partial,
                        right,
                        left,
                        maximum_degree,
                    )
                )

        lengths = (grade, grade + 1, grade + 3)
        records.append(
            audit_case(
                "rank_deficient_heterogeneous_shift",
                grade,
                *gauge_case(
                    *heterogeneous_shift(lengths),
                    generator,
                ),
                maximum_degree,
            )
        )
        records.append(
            audit_case(
                "rank_zero_heterogeneous_shift",
                grade,
                *gauge_case(
                    *heterogeneous_shift(
                        (grade + 1, grade + 2)
                    ),
                    generator,
                ),
                maximum_degree,
            )
        )
        apex_multiplicity = 2
        records.append(
            audit_case(
                "repeated_grade_apex",
                grade,
                *gauge_case(
                    *heterogeneous_shift(
                        (grade,) * apex_multiplicity
                    ),
                    generator,
                ),
                maximum_degree,
            )
        )
    return records


def write_records(
    records: list[OneImageRecord],
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
    parser.add_argument("--maximum-grade", type=int, default=6)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_one_image_generator_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the complete one-image audit."""

    args = parse_args()
    if not 1 <= args.maximum_grade <= 6:
        raise ValueError("maximum-grade must lie between one and six")
    records = standard_records(args.maximum_grade)
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
