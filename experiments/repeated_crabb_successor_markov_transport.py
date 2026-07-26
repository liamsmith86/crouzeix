#!/usr/bin/env python3
"""Audit the gap-free Markov transport of one Stein direction."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import block_diag
import sympy as sp

from crabb_block_hardy_equality import format_float
from repeated_crabb_all_grade_preimage import (
    gauged_shift,
    stable_random_partial_isometry,
)
from repeated_crabb_cubic_markov_flux import (
    dual_commutator_coefficient,
    minimum_response_column,
    simplified_cubic_defect,
)
from repeated_crabb_elliptic_cokernel import (
    dual_stein_inverse,
    stein_inverse,
)
from repeated_crabb_markov_response import (
    transfer_channel_adjoint,
)
from repeated_crabb_moving_retightening_defect import (
    ellipse_tangent,
    exact_stein_inverse,
)
from repeated_crabb_transfer_flag import transfer_channel


Matrix = np.ndarray


@dataclass(frozen=True)
class SuccessorTransportRecord:
    """One exact or numerical successor-transport audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    direction_metric_norm: str
    direction_frame_norm: str
    successor_endpoint_norm: str
    decomposition_error: str
    endpoint_trace_absolute: str
    dual_pairing_error: str
    dirichlet_identity_error: str
    flux_bound_slack: str
    response_synthesis_error: str
    response_column_norm: str
    response_column_bound_slack: str
    all_checks_passed: bool


def hermitian_part(matrix: Matrix) -> Matrix:
    """Return the Hermitian part of a square matrix."""

    return (matrix + matrix.conj().T) / 2


def successor_defect(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    metric_direction: Matrix,
    frame_direction: Matrix,
) -> Matrix:
    """Return the first moving successor of one fixed-base direction."""

    operator_tangent = ellipse_tangent(operator, right, left)
    frame_tangent = (
        -2
        * operator.conj().T
        @ operator.conj().T
        @ right
    )
    successor = -(
        operator_tangent.conj().T
        @ metric_direction
        @ operator
        + operator.conj().T
        @ metric_direction
        @ operator_tangent
        + frame_tangent @ frame_direction.conj().T
        + frame_direction @ frame_tangent.conj().T
    )
    return hermitian_part(successor)


def successor_decomposition(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    metric_direction: Matrix,
    frame_direction: Matrix,
) -> Matrix:
    """Return L301's commutator-plus-remainder decomposition."""

    adjoint = operator.conj().T
    universal = simplified_cubic_defect(
        operator,
        right,
        left,
        metric_direction,
    )
    remainder = -2 * (
        adjoint
        @ adjoint
        @ frame_direction
        @ right.conj().T
        + right
        @ frame_direction.conj().T
        @ operator
        @ operator
    )
    return hermitian_part(universal + remainder)


def successor_endpoint(
    operator: Matrix,
    left: Matrix,
    successor: Matrix,
) -> Matrix:
    """Return the balanced upper endpoint of a successor forcing."""

    endpoint = (
        left.conj().T
        @ stein_inverse(operator, successor)
        @ left
    )
    return hermitian_part(endpoint)


def random_direction(
    operator: Matrix,
    right: Matrix,
    generator: np.random.Generator,
) -> tuple[Matrix, Matrix]:
    """Return one deterministic Stein direction with both frame parts."""

    frame = (
        generator.standard_normal(right.shape)
        + 1j * generator.standard_normal(right.shape)
    )
    frame /= max(1.0, float(np.linalg.norm(frame)))
    metric = stein_inverse(
        operator,
        right @ frame.conj().T + frame @ right.conj().T,
    )
    return metric, frame


def numerical_record(
    construction_kind: str,
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    seed: int,
    *,
    require_collapsed_response: bool = False,
) -> SuccessorTransportRecord:
    """Audit one numerical colligation and direction."""

    generator = np.random.default_rng(seed)
    metric, frame = random_direction(operator, right, generator)
    successor = successor_defect(
        operator,
        right,
        left,
        metric,
        frame,
    )
    decomposed = successor_decomposition(
        operator,
        right,
        left,
        metric,
        frame,
    )
    decomposition_error = float(
        np.linalg.norm(successor - decomposed)
    )
    endpoint = successor_endpoint(operator, left, successor)
    trace_absolute = abs(float(np.trace(endpoint).real))

    multiplicity = right.shape[1]
    raw_test = (
        generator.standard_normal((multiplicity, multiplicity))
        + 1j
        * generator.standard_normal((multiplicity, multiplicity))
    )
    endpoint_test = hermitian_part(raw_test)
    observability = dual_stein_inverse(
        operator,
        left @ endpoint_test @ left.conj().T,
    )
    projection = (
        np.eye(len(operator), dtype=complex)
        - right @ right.conj().T
    )
    defect_column = projection @ observability @ right
    commutator_coefficient = dual_commutator_coefficient(
        operator,
        right,
        left,
        defect_column,
    )
    predicted_pairing = (
        float(
            np.trace(
                metric @ commutator_coefficient
            ).real
        )
        - 4
        * float(
            np.trace(
                defect_column.conj().T
                @ operator.conj().T
                @ operator.conj().T
                @ frame
            ).real
        )
    )
    direct_pairing = float(
        np.trace(observability @ successor).real
    )
    endpoint_pairing = float(
        np.trace(endpoint_test @ endpoint).real
    )
    pairing_error = max(
        abs(direct_pairing - endpoint_pairing),
        abs(direct_pairing - predicted_pairing),
    )

    channel_adjoint = transfer_channel_adjoint(
        operator,
        right,
        left,
        endpoint_test,
    )
    markov_image = transfer_channel(
        operator,
        right,
        left,
        channel_adjoint,
    )
    dirichlet = float(
        np.trace(
            endpoint_test @ (endpoint_test - markov_image)
        ).real
    )
    dirichlet_error = abs(
        dirichlet - 2 * float(np.linalg.norm(defect_column) ** 2)
    )
    bound = (
        (
            20 * float(np.linalg.norm(metric))
            + 4 * float(np.linalg.norm(frame))
        )
        * float(np.linalg.norm(defect_column))
    )
    flux_slack = bound - abs(endpoint_pairing)

    if (
        require_collapsed_response
        and float(np.linalg.norm(endpoint)) < 2e-8
    ):
        response_column = np.zeros_like(frame)
        synthesis_error = float(np.linalg.norm(endpoint))
    else:
        response_column, synthesis_error = minimum_response_column(
            operator,
            right,
            left,
            endpoint,
        )
    response_norm = float(np.linalg.norm(response_column))
    response_bound = (
        10 * float(np.linalg.norm(metric))
        + 2 * float(np.linalg.norm(frame))
    )
    response_slack = response_bound - response_norm
    tolerance = 4e-8
    verified = bool(
        decomposition_error < tolerance
        and trace_absolute < tolerance
        and pairing_error < tolerance
        and dirichlet_error < tolerance
        and flux_slack > -tolerance
        and synthesis_error < tolerance
        and response_slack > -tolerance
        and (
            not require_collapsed_response
            or float(np.linalg.norm(endpoint)) < tolerance
        )
    )
    if not verified:
        raise RuntimeError(
            "the successor Markov-transport audit failed: "
            f"{construction_kind=}, {decomposition_error=}, "
            f"{trace_absolute=}, {pairing_error=}, "
            f"{synthesis_error=}"
        )
    return SuccessorTransportRecord(
        construction_kind=construction_kind,
        state_dimension=len(operator),
        defect_dimension=multiplicity,
        direction_metric_norm=format_float(float(np.linalg.norm(metric))),
        direction_frame_norm=format_float(float(np.linalg.norm(frame))),
        successor_endpoint_norm=format_float(
            float(np.linalg.norm(endpoint))
        ),
        decomposition_error=format_float(decomposition_error),
        endpoint_trace_absolute=format_float(trace_absolute),
        dual_pairing_error=format_float(pairing_error),
        dirichlet_identity_error=format_float(dirichlet_error),
        flux_bound_slack=format_float(flux_slack),
        response_synthesis_error=format_float(synthesis_error),
        response_column_norm=format_float(response_norm),
        response_column_bound_slack=format_float(response_slack),
        all_checks_passed=verified,
    )


def exact_rational_record() -> SuccessorTransportRecord:
    """Audit the successor formulas with exact rational arithmetic."""

    dimension = 6
    multiplicity = 2
    rotation_one = sp.eye(4)
    rotation_one[0, 0] = rotation_one[1, 1] = sp.Rational(3, 5)
    rotation_one[0, 1] = -sp.Rational(4, 5)
    rotation_one[1, 0] = sp.Rational(4, 5)
    rotation_two = sp.eye(4)
    rotation_two[1, 1] = rotation_two[2, 2] = sp.Rational(5, 13)
    rotation_two[1, 2] = -sp.Rational(12, 13)
    rotation_two[2, 1] = sp.Rational(12, 13)
    bridge = rotation_two * rotation_one

    identity = sp.eye(dimension)
    right = identity[:, :multiplicity]
    left = identity[:, 4:]
    operator = (
        identity[:, :4]
        * bridge
        * identity[:, 2:].T
    )
    right_projection = right * right.T
    left_projection = left * left.T
    projection = identity - right_projection
    frame = projection * sp.Matrix(
        [
            [1, 0],
            [0, 1],
            [1, 2],
            [-1, 1],
            [2, -1],
            [1, 3],
        ]
    )
    metric = exact_stein_inverse(
        operator,
        right * frame.T + frame * right.T,
    )
    reflection = (
        (identity + left_projection)
        * operator.T
        * (identity + right_projection)
    )
    operator_tangent = reflection - operator**3
    frame_tangent = -2 * operator.T**2 * right
    successor = -(
        operator_tangent.T * metric * operator
        + operator.T * metric * operator_tangent
        + frame_tangent * frame.T
        + frame * frame_tangent.T
    )

    leading = (
        (identity + right_projection)
        * operator
        * (identity + left_projection)
        * metric
        * operator
    )
    universal = (
        2 * (operator.T**2 * metric + metric * operator**2)
        - (
            operator.T**3 * metric * operator
            + operator.T * metric * operator**3
        )
        - leading
        - leading.T
    )
    remainder = -2 * (
        operator.T**2 * frame * right.T
        + right * frame.T * operator**2
    )

    endpoint_test = sp.Matrix([[2, 1], [1, -3]])
    observability = exact_stein_inverse(
        operator.T,
        left * endpoint_test * left.T,
    )
    defect_column = projection * observability * right
    commutator = operator * defect_column * right.T
    right_defect_commutator = (
        right * defect_column.T - defect_column * right.T
    )
    adjoint_commutator = (
        -right * defect_column.T * operator.T
    )
    reflected_commutator = (
        (identity + left_projection)
        * (
            operator.T * right_defect_commutator
            + adjoint_commutator * (identity + right_projection)
        )
    )
    coefficient = (
        2 * (operator * commutator + commutator * operator)
        - commutator * operator.T**3
        - (
            operator**2 * commutator
            + operator * commutator * operator
            + commutator * operator**2
        )
        * operator.T
        - commutator
        * (identity + right_projection)
        * operator
        * (identity + left_projection)
        - reflected_commutator * operator.T
    )
    endpoint = (
        left.T
        * exact_stein_inverse(operator, successor)
        * left
    )
    predicted_pairing = (
        sp.trace(metric * coefficient)
        - 4
        * sp.trace(
            defect_column.T * operator.T**2 * frame
        )
    )
    exact_checks = (
        successor == universal + remainder,
        sp.trace(endpoint) == 0,
        sp.trace(observability * successor)
        == sp.trace(endpoint_test * endpoint),
        sp.trace(observability * successor)
        == predicted_pairing,
    )
    if not all(exact_checks):
        raise RuntimeError(
            "the exact successor Markov-transport guard failed"
        )

    return SuccessorTransportRecord(
        construction_kind="exact_rational_noncommuting_direction",
        state_dimension=dimension,
        defect_dimension=multiplicity,
        direction_metric_norm="not_applicable",
        direction_frame_norm="not_applicable",
        successor_endpoint_norm=format_float(
            float(
                np.linalg.norm(
                    np.asarray(endpoint, dtype=float)
                )
            )
        ),
        decomposition_error="0",
        endpoint_trace_absolute="0",
        dual_pairing_error="0",
        dirichlet_identity_error="not_applicable",
        flux_bound_slack="exact_nonnegative_by_theorem",
        response_synthesis_error="not_applicable",
        response_column_norm="not_applicable",
        response_column_bound_slack="not_applicable",
        all_checks_passed=True,
    )


def direct_sum_case(seed: int) -> tuple[Matrix, Matrix, Matrix]:
    """Return a reducible colligation with nonscalar commutant."""

    first = stable_random_partial_isometry(
        7,
        2,
        np.random.default_rng(seed),
    )
    second = stable_random_partial_isometry(
        8,
        2,
        np.random.default_rng(seed + 1),
    )
    return (
        block_diag(first[0], second[0]),
        block_diag(first[1], second[1]),
        block_diag(first[2], second[2]),
    )


def standard_records() -> list[SuccessorTransportRecord]:
    """Return exact, irreducible, reducible, and apex audits."""

    records = [exact_rational_record()]
    for multiplicity in range(1, 5):
        operator, right, left = stable_random_partial_isometry(
            3 * multiplicity + 5,
            multiplicity,
            np.random.default_rng(790_000 + multiplicity),
        )
        records.append(
            numerical_record(
                "unstructured_partial_isometry",
                operator,
                right,
                left,
                791_000 + multiplicity,
            )
        )
    records.append(
        numerical_record(
            "reducible_two_block",
            *direct_sum_case(792_000),
            793_000,
        )
    )
    for length in (1, 2, 3, 4):
        operator, right, left = gauged_shift(
            (length,) * 2,
            794_000 + length,
        )
        records.append(
            numerical_record(
                "repeated_monomial_apex",
                operator,
                right,
                left,
                795_000 + length,
                require_collapsed_response=True,
            )
        )
    return records


def write_records(
    records: list[SuccessorTransportRecord],
    output: Path,
) -> str:
    """Write deterministic JSONL and return its SHA-256."""

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(
                json.dumps(
                    asdict(record),
                    sort_keys=True,
                    separators=(",", ":"),
                )
                + "\n"
            )
    return hashlib.sha256(output.read_bytes()).hexdigest()


def main() -> None:
    """Run the deterministic audit."""

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_successor_markov_transport_s70226.jsonl"
        ),
    )
    arguments = parser.parse_args()
    records = standard_records()
    digest = write_records(records, arguments.output)
    print(
        json.dumps(
            {
                "all_checks_passed": all(
                    record.all_checks_passed for record in records
                ),
                "output": str(arguments.output),
                "records": len(records),
                "sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
