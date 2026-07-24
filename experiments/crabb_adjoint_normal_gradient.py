#!/usr/bin/env python3
"""Regenerate the adjoint-gradient reduction of circular-normal faces.

For the locally stationary defect, the derivative of the logarithmic
condition number can be evaluated without differentiating either
endpoint:

    Z - T Z T^* = v_+v_+^*/lambda_+ - v_-v_-^*/lambda_-,
    D log(kappa)[dot T] = 2 Re tr(Z T^* P dot T).

This checker derives every series in that identity exactly.  It verifies
the nonzero grade-one exception, the termwise-zero grade-two face, and
the first genuinely covariant grade-three cancellation, where five
nonzero pairings telescope to zero.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Sequence

import sympy as sp

from crabb_circular_normal_series import (
    endpoint_condition_coefficients,
    physical_reflected_path,
    real_circular_normal_direction,
    reciprocal_reversal_defect_jets,
)
from general_crabb_weighted_series import inverse_riemann_series
from rank_one_stein_series import (
    simple_diagonal_eigenpair_coefficients,
    stein_gramian_series,
    stein_series_from_forcing,
)


CASES = (
    (4, 1, True, False),
    (5, 2, True, False),
    (7, 3, False, True),
    (9, 4, False, False),
)


@dataclass(frozen=True)
class AdjointNormalRecord:
    """One exact adjoint normal-gradient record."""

    dimension: int
    length: int
    reflected_grade: int
    target_weight: int
    normal_mode: int
    independent_amplitudes: bool
    direct_endpoint_cross_recomputed: bool
    direct_condition_derivative: str | None
    adjoint_log_derivative: str
    direct_equals_four_times_adjoint: bool | None
    defect_stationarity_through_face_degree: bool
    gradient_persymmetric_through_face_degree: bool
    endpoint_commutator_verified: bool
    gradient_degree_pairings: tuple[str, ...]
    commutator_boundary_pairings: tuple[str, ...]
    companion_slice_pairings: tuple[str, ...]
    commutator_boundary_sum: str
    companion_slice_sum: str
    companion_decomposition_verified: bool
    every_target_pairing_zero: bool


def scalar_series_product(
    left: Sequence[sp.Expr],
    right: Sequence[sp.Expr],
) -> list[sp.Expr]:
    """Multiply two equally truncated scalar series."""

    order = min(len(left), len(right)) - 1
    return [
        sp.simplify(
            sum(
                (
                    left[source_degree]
                    * right[degree - source_degree]
                    for source_degree in range(degree + 1)
                ),
                sp.Integer(0),
            )
        )
        for degree in range(order + 1)
    ]


def scalar_series_inverse(
    coefficients: Sequence[sp.Expr],
) -> list[sp.Expr]:
    """Invert a scalar series."""

    inverse = [sp.simplify(1 / coefficients[0])]
    for degree in range(1, len(coefficients)):
        inverse.append(
            sp.simplify(
                -sum(
                    (
                        coefficients[source_degree]
                        * inverse[degree - source_degree]
                        for source_degree in range(1, degree + 1)
                    ),
                    sp.Integer(0),
                )
                / coefficients[0]
            )
        )
    return inverse


def normalized_projector_over_eigenvalue(
    eigenvalues: Sequence[sp.Expr],
    eigenvectors: Sequence[sp.Matrix],
) -> list[sp.Matrix]:
    """Return the series ``xx^T / ((x^Tx) lambda)``."""

    order = len(eigenvalues) - 1
    dimension = eigenvectors[0].rows
    outer = []
    norm = []
    for degree in range(order + 1):
        outer.append(
            sp.simplify(
                sum(
                    (
                        eigenvectors[source_degree]
                        * eigenvectors[degree - source_degree].T
                        for source_degree in range(degree + 1)
                    ),
                    sp.zeros(dimension),
                )
            )
        )
        norm.append(
            sp.simplify(
                sum(
                    (
                        (
                            eigenvectors[source_degree].T
                            * eigenvectors[degree - source_degree]
                        )[0]
                        for source_degree in range(degree + 1)
                    ),
                    sp.Integer(0),
                )
            )
        )

    reciprocal = scalar_series_inverse(
        scalar_series_product(norm, eigenvalues)
    )
    return [
        sp.simplify(
            sum(
                (
                    reciprocal[source_degree]
                    * outer[degree - source_degree]
                    for source_degree in range(degree + 1)
                ),
                sp.zeros(dimension),
            )
        )
        for degree in range(order + 1)
    ]


def ordinary_triple_series(
    left: Sequence[sp.Matrix],
    middle: Sequence[sp.Matrix],
    right: Sequence[sp.Matrix],
) -> list[sp.Matrix]:
    """Multiply three equally truncated ordinary matrix series."""

    order = min(len(left), len(middle), len(right)) - 1
    dimension = left[0].rows
    result = []
    for degree in range(order + 1):
        coefficient = sp.zeros(dimension)
        for left_degree in range(degree + 1):
            for middle_degree in range(degree - left_degree + 1):
                right_degree = degree - left_degree - middle_degree
                coefficient += (
                    left[left_degree]
                    * middle[middle_degree]
                    * right[right_degree]
                )
        result.append(sp.simplify(coefficient))
    return result


def coefficient_vector(
    vector: sp.Matrix,
    parameter: sp.Symbol,
    order: int,
) -> list[sp.Matrix]:
    """Extract an ordinary vector power series."""

    return [
        vector.applyfunc(
            lambda entry: sp.expand(entry).coeff(parameter, degree)
        )
        for degree in range(order + 1)
    ]


def reversal_matrix(dimension: int) -> sp.Matrix:
    """Return coordinate reversal."""

    reversal = sp.zeros(dimension)
    for index in range(dimension):
        reversal[index, dimension - 1 - index] = 1
    return reversal


def crabb_companion_decomposition(
    matrix: sp.Matrix,
    crabb: sp.Matrix,
) -> tuple[sp.Matrix, sp.Matrix]:
    """Write ``matrix=[crabb,X]+H`` with first row ``X=0`` and bottom-row ``H``."""

    dimension = matrix.rows
    length = dimension - 1
    weights = [crabb[index, index + 1] for index in range(length)]
    conjugator = sp.zeros(dimension)
    companion = sp.zeros(dimension)
    for row in range(length):
        for column in range(dimension):
            previous = (
                weights[column - 1] * conjugator[row, column - 1]
                if column > 0
                else sp.Integer(0)
            )
            conjugator[row + 1, column] = sp.expand(
                (matrix[row, column] + previous) / weights[row]
            )
    for column in range(dimension):
        previous = (
            weights[column - 1]
            * conjugator[length, column - 1]
            if column > 0
            else sp.Integer(0)
        )
        companion[length, column] = sp.expand(
            matrix[length, column] + previous
        )
    if (
        sp.simplify(
            crabb * conjugator
            - conjugator * crabb
            + companion
            - matrix
        )
        != sp.zeros(dimension)
    ):
        raise AssertionError("the Crabb companion decomposition failed")
    return conjugator, companion


def companion_decomposition_series(
    operator: Sequence[sp.Matrix],
    variation: Sequence[sp.Matrix],
) -> tuple[list[sp.Matrix], list[sp.Matrix]]:
    """Decompose a variation series as ``[T,X]+H`` recursively."""

    conjugators: list[sp.Matrix] = []
    companions: list[sp.Matrix] = []
    for degree, coefficient in enumerate(variation):
        residual = coefficient.copy()
        for source_degree in range(1, degree + 1):
            residual -= (
                operator[source_degree]
                * conjugators[degree - source_degree]
                - conjugators[degree - source_degree]
                * operator[source_degree]
            )
        conjugator, companion = crabb_companion_decomposition(
            sp.simplify(residual),
            operator[0],
        )
        conjugators.append(conjugator)
        companions.append(companion)

        regenerated = companion.copy()
        for source_degree in range(degree + 1):
            regenerated += (
                operator[source_degree]
                * conjugators[degree - source_degree]
                - conjugators[degree - source_degree]
                * operator[source_degree]
            )
        if sp.simplify(regenerated - coefficient) != sp.zeros(
            coefficient.rows
        ):
            raise AssertionError(
                "the companion series did not regenerate the variation"
            )
    return conjugators, companions


def make_record(
    dimension: int,
    grade: int,
    recompute_endpoint_cross: bool = True,
    independent_amplitudes: bool = False,
) -> AdjointNormalRecord:
    """Derive one exact adjoint-gradient face."""

    length = dimension - 1
    face_degree = grade + 1
    target_weight = 2 * face_degree
    normal_mode = length + 2 - grade
    epsilon, strong = sp.symbols("epsilon strong", real=True)
    if independent_amplitudes:
        amplitude, ellipse = sp.symbols(
            "amplitude ellipse",
            real=True,
        )
    else:
        amplitude = ellipse = sp.Integer(1)

    normal = real_circular_normal_direction(dimension, normal_mode)
    if recompute_endpoint_cross:
        path = physical_reflected_path(
            dimension=dimension,
            equality_grade=grade,
            strong_parameter=strong,
            strong_direction=normal,
            strong_degree=face_degree,
            order=target_weight,
            amplitude=amplitude,
            ellipse=ellipse,
        )
        _, strong_operator = inverse_riemann_series(
            path,
            target_weight,
        )
        operator = [
            coefficient.subs(strong, 0)
            for coefficient in strong_operator
        ]
        operator_derivative = [
            sp.simplify(sp.diff(coefficient, strong).subs(strong, 0))
            for coefficient in strong_operator
        ]
    else:
        # The strong insertion starts at ``face_degree``.  Through twice
        # that weight every coefficient is at most quadratic in it, so the
        # central difference at +/-1 is the exact first derivative.
        operators = []
        for strong_value in (0, 1, -1):
            path = physical_reflected_path(
                dimension=dimension,
                equality_grade=grade,
                strong_parameter=sp.Integer(strong_value),
                strong_direction=normal,
                strong_degree=face_degree,
                order=target_weight,
                amplitude=amplitude,
                ellipse=ellipse,
            )
            _, value = inverse_riemann_series(path, target_weight)
            operators.append(value)
        operator = operators[0]
        operator_derivative = [
            sp.simplify((plus - minus) / 2)
            for plus, minus in zip(
                operators[1],
                operators[2],
                strict=True,
            )
        ]
        strong_operator = None

    defect = reciprocal_reversal_defect_jets(
        operator,
        epsilon,
        jet_count=face_degree,
    )
    gramian = stein_gramian_series(
        operator,
        defect,
        epsilon,
        face_degree,
    )
    lower_values, lower_vectors = (
        simple_diagonal_eigenpair_coefficients(
            gramian,
            endpoint=0,
            base_eigenvalue=1,
        )
    )
    upper_values, upper_vectors = (
        simple_diagonal_eigenpair_coefficients(
            gramian,
            endpoint=dimension - 1,
            base_eigenvalue=4,
        )
    )
    lower_projector = normalized_projector_over_eigenvalue(
        lower_values,
        lower_vectors,
    )
    upper_projector = normalized_projector_over_eigenvalue(
        upper_values,
        upper_vectors,
    )
    endpoint_gradient = [
        sp.simplify(upper - lower)
        for upper, lower in zip(
            upper_projector,
            lower_projector,
            strict=True,
        )
    ]

    adjoint_operator = [coefficient.T for coefficient in operator]
    adjoint_metric = stein_series_from_forcing(
        adjoint_operator,
        endpoint_gradient,
        expand_only=True,
    )
    defect_coefficients = coefficient_vector(
        defect,
        epsilon,
        target_weight,
    )
    stationarity = []
    for degree in range(face_degree + 1):
        stationarity.append(
            sp.simplify(
                sum(
                    (
                        adjoint_metric[source_degree]
                        * defect_coefficients[degree - source_degree]
                        for source_degree in range(degree + 1)
                    ),
                    sp.zeros(dimension, 1),
                )
            )
        )
    stationary = all(
        value == sp.zeros(dimension, 1)
        for value in stationarity
    )
    if not stationary:
        raise AssertionError("the adjoint metric did not annihilate the defect")

    operator_gradient = ordinary_triple_series(
        gramian,
        operator,
        adjoint_metric,
    )
    reversal = reversal_matrix(dimension)
    gradient_persymmetric = all(
        sp.simplify(
            reversal * operator_gradient[degree].T * reversal
            - operator_gradient[degree]
        )
        == sp.zeros(dimension)
        for degree in range(face_degree + 1)
    )
    if not gradient_persymmetric:
        raise AssertionError("reciprocal reversal did not fix the gradient")

    identity_series = [
        sp.eye(dimension),
        *[sp.zeros(dimension) for _ in range(face_degree)],
    ]
    operator_adjoint = [
        coefficient.T
        for coefficient in operator[: face_degree + 1]
    ]
    adjoint_times_gradient = ordinary_triple_series(
        operator_adjoint,
        operator_gradient,
        identity_series,
    )
    gradient_times_adjoint = ordinary_triple_series(
        operator_gradient,
        operator_adjoint,
        identity_series,
    )
    metric_times_endpoint = ordinary_triple_series(
        gramian,
        endpoint_gradient,
        identity_series,
    )
    commutator_verified = all(
        sp.simplify(
            adjoint_times_gradient[degree]
            - gradient_times_adjoint[degree]
            - metric_times_endpoint[degree]
        )
        == sp.zeros(dimension)
        for degree in range(face_degree + 1)
    )
    if not commutator_verified:
        raise AssertionError("the endpoint commutator law failed")

    pairings = [
        sp.factor(
            sp.trace(
                operator_gradient[degree].T
                * operator_derivative[target_weight - degree]
            )
        )
        for degree in range(face_degree + 1)
    ]
    adjoint_derivative = sp.factor(2 * sum(pairings))

    unity = [
        sp.Integer(1),
        *[sp.Integer(0) for _ in range(face_degree)],
    ]
    lower_unit_projector = normalized_projector_over_eigenvalue(
        unity,
        lower_vectors,
    )
    upper_unit_projector = normalized_projector_over_eigenvalue(
        unity,
        upper_vectors,
    )
    endpoint_projector_difference = [
        sp.simplify(upper - lower)
        for upper, lower in zip(
            upper_unit_projector,
            lower_unit_projector,
            strict=True,
        )
    ]
    conjugators, companions = companion_decomposition_series(
        operator,
        operator_derivative,
    )
    boundary_pairings = [
        sp.factor(
            sp.trace(
                endpoint_projector_difference[degree]
                * conjugators[target_weight - degree]
            )
        )
        for degree in range(face_degree + 1)
    ]
    companion_pairings = [
        sp.factor(
            sp.trace(
                operator_gradient[degree].T
                * companions[target_weight - degree]
            )
        )
        for degree in range(face_degree + 1)
    ]
    boundary_sum = sp.factor(sum(boundary_pairings))
    companion_sum = sp.factor(sum(companion_pairings))
    decomposition_verified = sp.simplify(
        sum(pairings)
        - boundary_sum
        - companion_sum
    ) == 0
    if not decomposition_verified:
        raise AssertionError(
            "the commutator/companion pairing balance failed"
        )

    direct_derivative: sp.Expr | None = None
    matches: bool | None = None
    if recompute_endpoint_cross:
        if strong_operator is None:
            raise AssertionError("the direct strong series is unavailable")
        _, _, ratio = endpoint_condition_coefficients(
            strong_operator,
            defect,
            epsilon,
            target_weight,
        )
        direct_derivatives = [
            sp.factor(
                sp.diff(coefficient, strong).subs(strong, 0)
            )
            for coefficient in ratio
        ]
        if any(
            value != 0
            for value in direct_derivatives[:target_weight]
        ):
            raise AssertionError(
                "the normal derivative appeared below its face"
            )
        direct_derivative = direct_derivatives[target_weight]
        matches = sp.simplify(
            direct_derivative - 4 * adjoint_derivative
        ) == 0
        if not matches:
            raise AssertionError(
                "the adjoint derivative did not match endpoints"
            )

    every_pairing_zero = all(value == 0 for value in pairings)
    if grade == 1 and adjoint_derivative == 0:
        raise AssertionError("the grade-one exception disappeared")
    if grade >= 2 and adjoint_derivative != 0:
        raise AssertionError("the higher-grade adjoint selection changed")
    if grade == 2 and not every_pairing_zero:
        raise AssertionError("the grade-two termwise selection changed")
    if grade == 3 and every_pairing_zero:
        raise AssertionError("the grade-three covariant telescope disappeared")

    return AdjointNormalRecord(
        dimension=dimension,
        length=length,
        reflected_grade=grade,
        target_weight=target_weight,
        normal_mode=normal_mode,
        independent_amplitudes=independent_amplitudes,
        direct_endpoint_cross_recomputed=recompute_endpoint_cross,
        direct_condition_derivative=(
            str(direct_derivative)
            if direct_derivative is not None
            else None
        ),
        adjoint_log_derivative=str(adjoint_derivative),
        direct_equals_four_times_adjoint=matches,
        defect_stationarity_through_face_degree=stationary,
        gradient_persymmetric_through_face_degree=(
            gradient_persymmetric
        ),
        endpoint_commutator_verified=commutator_verified,
        gradient_degree_pairings=tuple(
            str(value) for value in pairings
        ),
        commutator_boundary_pairings=tuple(
            str(value) for value in boundary_pairings
        ),
        companion_slice_pairings=tuple(
            str(value) for value in companion_pairings
        ),
        commutator_boundary_sum=str(boundary_sum),
        companion_slice_sum=str(companion_sum),
        companion_decomposition_verified=decomposition_verified,
        every_target_pairing_zero=every_pairing_zero,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run the deterministic exact audit."""

    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as output:
        for (
            dimension,
            grade,
            recompute_endpoint_cross,
            independent_amplitudes,
        ) in CASES:
            record = make_record(
                dimension,
                grade,
                recompute_endpoint_cross,
                independent_amplitudes,
            )
            line = json.dumps(asdict(record), sort_keys=True)
            print(line, flush=True)
            output.write(line + "\n")


if __name__ == "__main__":
    main()
