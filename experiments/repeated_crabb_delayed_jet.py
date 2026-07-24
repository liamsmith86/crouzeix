"""Shared finite-jet machinery for delayed repeated-Crabb faces.

The module constructs the exact ellipse pullback, a coefficientwise
lower-Schur-tight Stein metric, and both generalized metric endpoints.
Individual grade audits supply only their additional perpendicular
defect-frame gauge.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from scipy.linalg import null_space, sqrtm

from crabb_palindromic_elliptic_hessian import direct_map_coefficients
from repeated_crabb_elliptic_cokernel import stein_inverse


Matrix = np.ndarray
GaugeAdjustment = Callable[
    [int, Matrix, Matrix, Matrix],
    Matrix,
]


@dataclass(frozen=True)
class DelayedMetricJet:
    """One lower-tight delayed metric jet."""

    partial: Matrix
    right: Matrix
    left: Matrix
    metric: Matrix
    metric_root: Matrix
    operator: Matrix
    operator_coefficients: tuple[Matrix, ...]
    metric_coefficients: tuple[Matrix, ...]
    frame_coefficients: tuple[Matrix, ...]


def series_multiply(
    left: list[Matrix],
    right: list[Matrix],
    degree: int,
) -> list[Matrix]:
    """Multiply compatible matrix series through ``degree``."""

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


def series_power(
    series: list[Matrix],
    power: int,
    degree: int,
) -> list[Matrix]:
    """Raise a square matrix series to an integer power."""

    dimension = len(series[0])
    result = [np.eye(dimension, dtype=complex)] + [
        np.zeros_like(series[0]) for _ in range(degree)
    ]
    for _ in range(power):
        result = series_multiply(result, series, degree)
    return result


def ellipse_operator_coefficients(
    operator: Matrix,
    degree: int,
) -> list[Matrix]:
    """Return ``phi_c(T+cT*)`` through ``degree`` exactly."""

    pencil = [operator, operator.conj().T] + [
        np.zeros_like(operator) for _ in range(degree - 1)
    ]
    result = [np.zeros_like(operator) for _ in range(degree + 1)]
    for index, coefficient in enumerate(
        direct_map_coefficients(degree, degree + 1)
    ):
        power = series_power(pencil, 2 * index + 1, degree)
        for scalar_degree in range(degree + 1):
            scalar = complex(coefficient.coefficient(scalar_degree))
            if scalar == 0:
                continue
            for matrix_degree, matrix in enumerate(power):
                total_degree = scalar_degree + matrix_degree
                if total_degree <= degree:
                    result[total_degree] += scalar * matrix
    return result


def inverse_series(series: list[Matrix]) -> list[Matrix]:
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
    metric_coefficients: list[Matrix | None],
    endpoint: Matrix,
    endpoint_value: float,
    degree: int,
) -> Matrix:
    """Return the Schur cross-square coefficient at one endpoint."""

    complement = null_space(endpoint.conj().T)
    zero = np.zeros_like(metric_coefficients[0])
    cross: list[Matrix] = []
    denominator: list[Matrix] = []
    for index in range(degree + 1):
        metric = (
            metric_coefficients[index]
            if index < len(metric_coefficients)
            and metric_coefficients[index] is not None
            else zero
        )
        cross.append(endpoint.conj().T @ metric @ complement)
        denominator.append(
            complement.conj().T @ metric @ complement
        )
    denominator[0] -= endpoint_value * np.eye(complement.shape[1])
    denominator_inverse = inverse_series(denominator)
    return series_multiply(
        series_multiply(cross, denominator_inverse, degree),
        [coefficient.conj().T for coefficient in cross],
        degree,
    )[degree]


def axis_perpendicular_column(
    partial: Matrix,
    right: Matrix,
    degree: int,
) -> Matrix:
    """Return the zero-reflection perpendicular frame coefficient."""

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


def construct_delayed_metric_jet(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    degree: int,
    adjustment: GaugeAdjustment | None = None,
) -> DelayedMetricJet:
    """Construct a lower-Schur-tight metric through ``degree``."""

    dimension = len(partial)
    identity = np.eye(dimension, dtype=complex)
    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    metric = 2 * identity - right_projection + 2 * left_projection
    metric_root = np.asarray(sqrtm(metric), dtype=complex)
    operator = np.linalg.inv(metric_root) @ partial @ metric_root
    operator_coefficients = ellipse_operator_coefficients(
        operator,
        degree,
    )
    metric_coefficients: list[Matrix | None] = (
        [metric] + [None] * degree
    )
    frame_coefficients: list[Matrix | None] = (
        [right] + [None] * degree
    )

    for index in range(1, degree + 1):
        fixed = np.zeros_like(operator)
        for left_degree in range(index + 1):
            for metric_degree in range(index + 1 - left_degree):
                right_degree = index - left_degree - metric_degree
                if metric_degree >= index:
                    continue
                fixed += (
                    operator_coefficients[left_degree].conj().T
                    @ metric_coefficients[metric_degree]
                    @ operator_coefficients[right_degree]
                )
        for left_degree in range(1, index):
            fixed += (
                frame_coefficients[left_degree]
                @ frame_coefficients[
                    index - left_degree
                ].conj().T
            )

        lower_target = schur_cross_coefficient(
            metric_coefficients,
            right,
            1,
            index,
        )
        parallel = (
            lower_target - right.conj().T @ fixed @ right
        ) / 2
        perpendicular = axis_perpendicular_column(
            partial,
            right,
            index,
        )
        if adjustment is not None:
            perpendicular += adjustment(
                index,
                partial,
                right,
                left,
            )
        frame_coefficients[index] = (
            right @ parallel + metric_root @ perpendicular
        )
        forcing = (
            fixed
            + right @ frame_coefficients[index].conj().T
            + frame_coefficients[index] @ right.conj().T
        )
        metric_coefficients[index] = stein_inverse(
            operator,
            forcing,
        )

    return DelayedMetricJet(
        partial=partial,
        right=right,
        left=left,
        metric=metric,
        metric_root=metric_root,
        operator=operator,
        operator_coefficients=tuple(operator_coefficients),
        metric_coefficients=tuple(metric_coefficients),
        frame_coefficients=tuple(frame_coefficients),
    )


def endpoint_coefficient(
    jet: DelayedMetricJet,
    degree: int,
    *,
    upper: bool,
) -> Matrix:
    """Return one generalized endpoint Schur coefficient."""

    endpoint = jet.left if upper else jet.right
    endpoint_value = 4 if upper else 1
    return (
        endpoint.conj().T
        @ jet.metric_coefficients[degree]
        @ endpoint
        - schur_cross_coefficient(
            list(jet.metric_coefficients),
            endpoint,
            endpoint_value,
            degree,
        )
    )
