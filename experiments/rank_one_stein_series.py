#!/usr/bin/env python3
"""Exact rank-one-defect Stein and eigenvalue power-series helpers."""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp


def matrix_power_series(
    coefficients: Sequence[sp.Matrix],
    epsilon: sp.Symbol,
) -> sp.Matrix:
    """Assemble a matrix power series as a symbolic matrix."""

    value = sp.zeros(coefficients[0].rows)
    for degree, coefficient in enumerate(coefficients):
        value += epsilon**degree * coefficient
    return value


def stein_gramian_series(
    operator: Sequence[sp.Matrix],
    defect: sp.Matrix,
    epsilon: sp.Symbol,
    order: int,
) -> list[sp.Matrix]:
    """Solve ``P-T^*P*T=defect*defect^*`` through ``order``."""

    dimension = operator[0].rows
    base = operator[0]
    real_path = all(
        entry.is_real is True
        for coefficient in operator
        for entry in coefficient
    ) and all(entry.is_real is True for entry in defect)

    def adjoint(matrix: sp.Matrix) -> sp.Matrix:
        return matrix.T if real_path else matrix.conjugate().T

    defect_square = sp.expand(defect * adjoint(defect))
    forcing = [
        defect_square.applyfunc(
            lambda entry: sp.expand(entry).coeff(epsilon, degree)
        )
        for degree in range(order + 1)
    ]
    gramian: list[sp.Matrix] = []
    for degree in range(order + 1):
        right_hand_side = forcing[degree]
        for left_degree in range(degree + 1):
            for right_degree in range(degree + 1 - left_degree):
                metric_degree = degree - left_degree - right_degree
                if metric_degree < len(gramian):
                    right_hand_side += (
                        adjoint(operator[left_degree])
                        * gramian[metric_degree]
                        * operator[right_degree]
                    )
        coefficient = sp.zeros(dimension)
        for power in range(dimension):
            coefficient += (
                (adjoint(base) ** power)
                * right_hand_side
                * (base**power)
            )
        gramian.append(sp.simplify(coefficient))
    return gramian


def simple_eigenvalue_series(
    matrix: sp.Matrix,
    base_eigenvalue: int,
    epsilon: sp.Symbol,
    order: int,
) -> sp.Expr:
    """Lift a simple eigenvalue of the constant matrix through ``order``."""

    dimension = matrix.rows
    spectral_parameter = sp.symbols("spectral_parameter")
    series = sp.Integer(base_eigenvalue)
    characteristic = (
        spectral_parameter * sp.eye(dimension) - matrix
    ).det()
    for degree in range(1, order + 1):
        unknown = sp.symbols(f"eigenvalue_{base_eigenvalue}_{degree}")
        equation = sp.expand(
            sp.series(
                characteristic.subs(
                    spectral_parameter,
                    series + unknown * epsilon**degree,
                ),
                epsilon,
                0,
                degree + 1,
            ).removeO()
        ).coeff(epsilon, degree)
        solution = sp.solve(equation, unknown)
        if len(solution) != 1:
            raise AssertionError("a simple Gramian eigenvalue did not lift uniquely")
        series = sp.expand(series + solution[0] * epsilon**degree)
    return series


def simple_diagonal_eigenvalue_coefficients(
    matrix_coefficients: Sequence[sp.Matrix],
    endpoint: int,
    base_eigenvalue: int,
) -> list[sp.Expr]:
    """Lift a simple eigenvalue when the constant matrix is diagonal.

    This is the Rayleigh--Schrödinger recurrence in intermediate
    normalization.  Unlike :func:`simple_eigenvalue_series`, it never forms
    a symbolic determinant, so it remains practical when the matrix entries
    contain several auxiliary symbols.
    """

    order = len(matrix_coefficients) - 1
    dimension = matrix_coefficients[0].rows
    constant = matrix_coefficients[0]
    if constant != sp.diag(*constant.diagonal()):
        raise ValueError("the constant matrix must be diagonal")
    if any(
        constant[index, index] == base_eigenvalue
        for index in range(dimension)
        if index != endpoint
    ):
        raise ValueError("the requested base eigenvalue must be simple")

    eigenvectors = [sp.eye(dimension)[:, endpoint]]
    eigenvalues = [sp.Integer(base_eigenvalue)]
    for degree in range(1, order + 1):
        image = sum(
            (
                matrix_coefficients[source_degree]
                * eigenvectors[degree - source_degree]
                for source_degree in range(1, degree + 1)
            ),
            sp.zeros(dimension, 1),
        )
        eigenvalue = sp.simplify(image[endpoint])
        eigenvalues.append(eigenvalue)

        right_hand_side = sp.zeros(dimension, 1)
        for source_degree in range(1, degree + 1):
            right_hand_side -= (
                matrix_coefficients[source_degree]
                - eigenvalues[source_degree] * sp.eye(dimension)
            ) * eigenvectors[degree - source_degree]

        eigenvector = sp.zeros(dimension, 1)
        for index in range(dimension):
            if index == endpoint:
                continue
            eigenvector[index] = sp.simplify(
                right_hand_side[index]
                / (constant[index, index] - base_eigenvalue)
            )
        eigenvectors.append(eigenvector)
    return eigenvalues


def diagonal_gramian_condition_series(
    operator: Sequence[sp.Matrix],
    defect: sp.Matrix,
    epsilon: sp.Symbol,
    order: int,
) -> sp.Expr:
    """Return the endpoint ratio using the diagonal-eigenpair recurrence."""

    gramian = stein_gramian_series(operator, defect, epsilon, order)
    lower = simple_diagonal_eigenvalue_coefficients(
        gramian,
        endpoint=0,
        base_eigenvalue=1,
    )
    upper = simple_diagonal_eigenvalue_coefficients(
        gramian,
        endpoint=gramian[0].rows - 1,
        base_eigenvalue=4,
    )
    ratio = [sp.Integer(4), *[sp.Integer(0) for _ in range(order)]]
    for degree in range(1, order + 1):
        ratio[degree] = sp.expand(
            upper[degree]
            - sum(
                lower[source_degree] * ratio[degree - source_degree]
                for source_degree in range(1, degree + 1)
            )
        )
    return sum(
        ratio[degree] * epsilon**degree
        for degree in range(order + 1)
    )


def gramian_condition_series(
    operator: Sequence[sp.Matrix],
    defect: sp.Matrix,
    epsilon: sp.Symbol,
    order: int,
) -> sp.Expr:
    """Return the extreme-eigenvalue ratio of a Stein Gramian series."""

    gramian = stein_gramian_series(operator, defect, epsilon, order)
    full_gramian = matrix_power_series(gramian, epsilon)
    lower = simple_eigenvalue_series(full_gramian, 1, epsilon, order)
    upper = simple_eigenvalue_series(full_gramian, 4, epsilon, order)
    return sp.expand(
        sp.series(upper / lower, epsilon, 0, order + 1).removeO()
    )
