#!/usr/bin/env python3
"""Audit the gap-free Markov bound for the L299 cubic defect."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import block_diag, null_space
import sympy as sp

from crabb_block_hardy_equality import format_float
from repeated_crabb_all_grade_preimage import (
    stable_random_partial_isometry,
)
from repeated_crabb_boundary_metric_flag import rank_chain_case
from repeated_crabb_elliptic_cokernel import (
    dual_stein_inverse,
    stein_inverse,
)
from repeated_crabb_markov_response import (
    transfer_channel_adjoint,
)
from repeated_crabb_moving_retightening_defect import (
    cubic_moving_defect,
    exact_stein_inverse,
)
from repeated_crabb_oriented_retightening_transport import (
    oriented_retightening_direction,
)
from repeated_crabb_transfer_flag import transfer_channel


Matrix = np.ndarray


@dataclass(frozen=True)
class CubicFluxRecord:
    """One exact or numerical cubic-flux audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    cubic_endpoint_norm: str
    simplified_defect_error: str
    scalar_coboundary_error: str
    commutator_formula_error: str
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


def simplified_cubic_defect(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    metric_direction: Matrix,
) -> Matrix:
    """Return L300's frame-free formula for the cubic defect."""

    identity = np.eye(len(operator), dtype=complex)
    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    adjoint = operator.conj().T
    leading = (
        (identity + right_projection)
        @ operator
        @ (identity + left_projection)
        @ metric_direction
        @ operator
    )
    return (
        2
        * (
            adjoint @ adjoint @ metric_direction
            + metric_direction @ operator @ operator
        )
        - (
            adjoint @ adjoint @ adjoint @ metric_direction @ operator
            + adjoint
            @ metric_direction
            @ operator
            @ operator
            @ operator
        )
        - leading
        - leading.conj().T
    )


def scalar_coboundary(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
) -> Matrix:
    """Return the operator multiplying a scalar dual test."""

    identity = np.eye(len(operator), dtype=complex)
    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    adjoint = operator.conj().T
    return (
        2 * adjoint @ adjoint
        + 2 * operator @ operator
        - operator @ adjoint @ adjoint @ adjoint
        - operator @ operator @ operator @ adjoint
        - operator
        @ (identity + right_projection)
        @ operator
        @ (identity + left_projection)
        - (identity + left_projection)
        @ adjoint
        @ (identity + right_projection)
        @ adjoint
    )


def dual_commutator_coefficient(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    defect_column: Matrix,
) -> Matrix:
    """Return the coefficient ``K_Y`` in ``tr(R K_Y)``."""

    identity = np.eye(len(operator), dtype=complex)
    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    adjoint = operator.conj().T
    commutator = operator @ defect_column @ right.conj().T
    right_defect_commutator = (
        right @ defect_column.conj().T
        - defect_column @ right.conj().T
    )
    adjoint_commutator = (
        -right @ defect_column.conj().T @ adjoint
    )
    reflected_commutator = (
        (identity + left_projection)
        @ (
            adjoint @ right_defect_commutator
            + adjoint_commutator @ (identity + right_projection)
        )
    )
    return (
        2
        * (
            operator @ commutator
            + commutator @ operator
        )
        - commutator @ adjoint @ adjoint @ adjoint
        - (
            operator @ operator @ commutator
            + operator @ commutator @ operator
            + commutator @ operator @ operator
        )
        @ adjoint
        - commutator
        @ (identity + right_projection)
        @ operator
        @ (identity + left_projection)
        - reflected_commutator @ adjoint
    )


def balanced_endpoint_response(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    column: Matrix,
) -> Matrix:
    """Apply the balanced free-column endpoint response."""

    forcing = (
        right @ column.conj().T
        + column @ right.conj().T
    )
    return (
        left.conj().T
        @ stein_inverse(operator, forcing)
        @ left
    )


def hermitian_basis(size: int) -> list[Matrix]:
    """Return an orthonormal real basis of Hermitian matrices."""

    basis = []
    for row in range(size):
        diagonal = np.zeros((size, size), dtype=complex)
        diagonal[row, row] = 1
        basis.append(diagonal)
    for row in range(size):
        for column in range(row + 1, size):
            symmetric = np.zeros((size, size), dtype=complex)
            symmetric[row, column] = 1 / np.sqrt(2)
            symmetric[column, row] = 1 / np.sqrt(2)
            basis.append(symmetric)
            skew = np.zeros((size, size), dtype=complex)
            skew[row, column] = 1j / np.sqrt(2)
            skew[column, row] = -1j / np.sqrt(2)
            basis.append(skew)
    return basis


def hermitian_coordinates(
    matrix: Matrix,
    basis: list[Matrix],
) -> Matrix:
    """Return real Hilbert--Schmidt coordinates."""

    return np.asarray(
        [
            float(np.trace(element @ matrix).real)
            for element in basis
        ]
    )


def minimum_response_column(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    target: Matrix,
) -> tuple[Matrix, float]:
    """Return a numerical minimum-norm perpendicular response column."""

    multiplicity = right.shape[1]
    complement = null_space(right.conj().T, rcond=1e-11)
    endpoint_basis = hermitian_basis(multiplicity)
    columns = []
    state_columns = []
    for state_index in range(complement.shape[1]):
        for copy_index in range(multiplicity):
            for phase in (1, 1j):
                column = np.zeros(
                    (len(operator), multiplicity),
                    dtype=complex,
                )
                column[:, copy_index] = (
                    phase * complement[:, state_index]
                )
                state_columns.append(column)
                response = balanced_endpoint_response(
                    operator,
                    right,
                    left,
                    column,
                )
                columns.append(
                    hermitian_coordinates(response, endpoint_basis)
                )
    response_matrix = np.column_stack(columns)
    target_coordinates = hermitian_coordinates(
        target,
        endpoint_basis,
    )
    coefficients = np.linalg.lstsq(
        response_matrix,
        target_coordinates,
        rcond=1e-11,
    )[0]
    selected = sum(
        coefficient * column
        for coefficient, column in zip(
            coefficients,
            state_columns,
            strict=True,
        )
    )
    residual = float(
        np.linalg.norm(
            response_matrix @ coefficients - target_coordinates
        )
    )
    return selected, residual


def numerical_record(
    construction_kind: str,
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    seed: int,
) -> CubicFluxRecord:
    """Audit one numerical colligation."""

    dimension = len(operator)
    multiplicity = right.shape[1]
    identity = np.eye(dimension, dtype=complex)
    right_projection = right @ right.conj().T
    generator = np.random.default_rng(seed)

    state_defect, endpoint, _ = cubic_moving_defect(
        operator,
        right,
        left,
    )
    direction = oriented_retightening_direction(
        operator,
        right,
        left,
        1,
    )
    simplified = simplified_cubic_defect(
        operator,
        right,
        left,
        direction.metric,
    )
    simplified_error = float(
        np.linalg.norm(state_defect - simplified)
    )
    scalar_error = float(
        np.linalg.norm(
            scalar_coboundary(operator, right, left)
        )
    )

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
    defect_column = (
        (identity - right_projection)
        @ observability
        @ right
    )
    coefficient = dual_commutator_coefficient(
        operator,
        right,
        left,
        defect_column,
    )

    adjoint = operator.conj().T
    right_defect = right_projection
    left_defect = left @ left.conj().T
    raw_coefficient = (
        2 * observability @ adjoint @ adjoint
        + 2 * operator @ operator @ observability
        - operator @ observability @ adjoint @ adjoint @ adjoint
        - operator
        @ operator
        @ operator
        @ observability
        @ adjoint
        - operator
        @ observability
        @ (identity + right_defect)
        @ operator
        @ (identity + left_defect)
        - (identity + left_defect)
        @ adjoint
        @ (identity + right_defect)
        @ observability
        @ adjoint
    )
    commutator_error = float(
        np.linalg.norm(raw_coefficient - coefficient)
    )
    direct_pairing = float(
        np.trace(observability @ state_defect).real
    )
    endpoint_pairing = float(
        np.trace(endpoint_test @ endpoint).real
    )
    commutator_pairing = float(
        np.trace(direction.metric @ coefficient).real
    )
    pairing_error = max(
        abs(direct_pairing - endpoint_pairing),
        abs(direct_pairing - commutator_pairing),
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
    flux_bound = (
        20
        * float(np.linalg.norm(direction.metric))
        * float(np.linalg.norm(defect_column))
    )
    flux_slack = flux_bound - abs(endpoint_pairing)

    response_column, synthesis_error = minimum_response_column(
        operator,
        right,
        left,
        endpoint,
    )
    response_column_norm = float(np.linalg.norm(response_column))
    response_column_slack = (
        10 * float(np.linalg.norm(direction.metric))
        - response_column_norm
    )
    tolerance = 2e-8
    verified = bool(
        simplified_error < tolerance
        and scalar_error < tolerance
        and commutator_error < tolerance
        and pairing_error < tolerance
        and dirichlet_error < tolerance
        and flux_slack > -tolerance
        and synthesis_error < tolerance
        and response_column_slack > -tolerance
    )
    if not verified:
        raise RuntimeError(
            "the cubic Markov-flux audit failed: "
            f"{construction_kind=}, {simplified_error=}, "
            f"{commutator_error=}, {pairing_error=}, "
            f"{synthesis_error=}"
        )
    return CubicFluxRecord(
        construction_kind=construction_kind,
        state_dimension=dimension,
        defect_dimension=multiplicity,
        cubic_endpoint_norm=format_float(float(np.linalg.norm(endpoint))),
        simplified_defect_error=format_float(simplified_error),
        scalar_coboundary_error=format_float(scalar_error),
        commutator_formula_error=format_float(commutator_error),
        dual_pairing_error=format_float(pairing_error),
        dirichlet_identity_error=format_float(dirichlet_error),
        flux_bound_slack=format_float(flux_slack),
        response_synthesis_error=format_float(synthesis_error),
        response_column_norm=format_float(response_column_norm),
        response_column_bound_slack=format_float(response_column_slack),
        all_checks_passed=verified,
    )


def exact_rational_record() -> CubicFluxRecord:
    """Audit every algebraic identity on one rational colligation."""

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
    domain = identity[:, 2:]
    range_frame = identity[:, :4]
    operator = range_frame * bridge * domain.T
    right_projection = right * right.T
    left_projection = left * left.T
    projection = identity - right_projection

    first_transfer = left.T * operator.T * right
    right_gram = first_transfer.T * first_transfer
    perpendicular = (
        -sp.Rational(1, 2)
        * projection
        * operator
        * left
        * first_transfer
    )
    orbit = exact_stein_inverse(
        operator,
        right * right_gram * right.T,
    )
    perpendicular_metric = exact_stein_inverse(
        operator,
        right * perpendicular.T + perpendicular * right.T,
    )
    metric_direction = -orbit + perpendicular_metric
    frame_direction = -right * right_gram / 2 + perpendicular

    reflection = (
        (identity + left_projection)
        * operator.T
        * (identity + right_projection)
    )
    operator_tangent = reflection - operator**3
    frame_tangent = -2 * operator.T**2 * right
    state_defect = -(
        operator_tangent.T * metric_direction * operator
        + operator.T * metric_direction * operator_tangent
        + frame_tangent * frame_direction.T
        + frame_direction * frame_tangent.T
    )

    leading = (
        (identity + right_projection)
        * operator
        * (identity + left_projection)
        * metric_direction
        * operator
    )
    simplified = (
        2 * (operator.T**2 * metric_direction + metric_direction * operator**2)
        - (
            operator.T**3 * metric_direction * operator
            + operator.T * metric_direction * operator**3
        )
        - leading
        - leading.T
    )
    scalar = (
        2 * operator.T**2
        + 2 * operator**2
        - operator * operator.T**3
        - operator**3 * operator.T
        - operator
        * (identity + right_projection)
        * operator
        * (identity + left_projection)
        - (identity + left_projection)
        * operator.T
        * (identity + right_projection)
        * operator.T
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
    raw_coefficient = (
        2 * observability * operator.T**2
        + 2 * operator**2 * observability
        - operator * observability * operator.T**3
        - operator**3 * observability * operator.T
        - operator
        * observability
        * (identity + right_projection)
        * operator
        * (identity + left_projection)
        - (identity + left_projection)
        * operator.T
        * (identity + right_projection)
        * observability
        * operator.T
    )
    endpoint = (
        left.T
        * exact_stein_inverse(operator, state_defect)
        * left
    )
    exact_checks = (
        state_defect == simplified,
        scalar == sp.zeros(dimension),
        raw_coefficient == coefficient,
        sp.trace(observability * state_defect)
        == sp.trace(endpoint_test * endpoint),
        sp.trace(observability * state_defect)
        == sp.trace(metric_direction * coefficient),
        sp.trace(endpoint) == 0,
    )
    if not all(exact_checks):
        raise RuntimeError("the exact cubic Markov-flux guard failed")

    endpoint_norm = float(
        np.linalg.norm(np.asarray(endpoint, dtype=float))
    )
    return CubicFluxRecord(
        construction_kind="exact_rational_rank_one_flag",
        state_dimension=dimension,
        defect_dimension=multiplicity,
        cubic_endpoint_norm=format_float(endpoint_norm),
        simplified_defect_error="0",
        scalar_coboundary_error="0",
        commutator_formula_error="0",
        dual_pairing_error="0",
        dirichlet_identity_error="not_applicable",
        flux_bound_slack="exact_nonnegative_by_theorem",
        response_synthesis_error="not_applicable",
        response_column_norm="not_applicable",
        response_column_bound_slack="not_applicable",
        all_checks_passed=True,
    )


def direct_sum_case(seed: int) -> tuple[Matrix, Matrix, Matrix]:
    """Return a reducible colligation with a nonscalar cokernel."""

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
    operator = block_diag(first[0], second[0])
    right = block_diag(first[1], second[1])
    left = block_diag(first[2], second[2])
    return operator, right, left


def standard_records() -> list[CubicFluxRecord]:
    """Return the exact guard and deterministic numerical suite."""

    records = [exact_rational_record()]
    for multiplicity in range(1, 5):
        for index in range(2):
            operator, right, left = stable_random_partial_isometry(
                3 * multiplicity + 4 + index,
                multiplicity,
                np.random.default_rng(
                    730_000 + 100 * multiplicity + index
                ),
            )
            records.append(
                numerical_record(
                    "unstructured_partial_isometry",
                    operator,
                    right,
                    left,
                    740_000 + 100 * multiplicity + index,
                )
            )
    for multiplicity in range(2, 5):
        operator, right, left, _, _ = rank_chain_case(
            multiplicity,
            750_000 + multiplicity,
            0.37,
        )
        records.append(
            numerical_record(
                "noncommuting_rank_chain",
                operator,
                right,
                left,
                760_000 + multiplicity,
            )
        )
    records.append(
        numerical_record(
            "reducible_two_block",
            *direct_sum_case(770_000),
            780_000,
        )
    )
    return records


def write_records(
    records: list[CubicFluxRecord],
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
            "repeated_crabb_cubic_markov_flux_s70226.jsonl"
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
