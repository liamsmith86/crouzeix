#!/usr/bin/env python3
"""Sparse exact series engine for the canonical p=3 Crabb block.

This specialized engine makes high-order certificate regeneration practical.
The older general symbolic helpers remain as an independent low-order audit.
"""

from __future__ import annotations

from itertools import permutations
from typing import TypeAlias

import sympy as sp


Laurent: TypeAlias = dict[int, sp.Expr]
Series: TypeAlias = list[Laurent]


def clean(value: sp.Expr) -> sp.Expr:
    return sp.expand(value)


def lp_add(left: Laurent, right: Laurent) -> Laurent:
    result = dict(left)
    for degree, coefficient in right.items():
        value = clean(result.get(degree, 0) + coefficient)
        if value == 0:
            result.pop(degree, None)
        else:
            result[degree] = value
    return result


def lp_scale(value: Laurent, scalar: sp.Expr) -> Laurent:
    if scalar == 0:
        return {}
    return {
        degree: coefficient * scalar
        for degree, coefficient in value.items()
        if coefficient * scalar != 0
    }


def lp_shift(value: Laurent, shift: int) -> Laurent:
    return {degree + shift: coefficient for degree, coefficient in value.items()}


def lp_mul(left: Laurent, right: Laurent) -> Laurent:
    result: Laurent = {}
    for left_degree, left_coefficient in left.items():
        for right_degree, right_coefficient in right.items():
            degree = left_degree + right_degree
            result[degree] = result.get(degree, 0) + left_coefficient * right_coefficient
    return {
        degree: clean(coefficient)
        for degree, coefficient in result.items()
        if clean(coefficient) != 0
    }


def lp_star(value: Laurent) -> Laurent:
    return {
        -degree: sp.conjugate(coefficient)
        for degree, coefficient in value.items()
    }


def series_mul(left: Series, right: Series, order: int) -> Series:
    result = [{} for _ in range(order + 1)]
    for degree in range(order + 1):
        for left_degree in range(degree + 1):
            right_degree = degree - left_degree
            if left_degree < len(left) and right_degree < len(right):
                result[degree] = lp_add(
                    result[degree],
                    lp_mul(left[left_degree], right[right_degree]),
                )
    return result


def series_exp(value: Series, order: int) -> Series:
    """Exponentiate a series with zero constant coefficient."""

    result: Series = [{0: sp.Integer(1)}, *[{} for _ in range(order)]]
    for degree in range(1, order + 1):
        coefficient: Laurent = {}
        for source_degree in range(1, degree + 1):
            if source_degree < len(value):
                product = lp_mul(value[source_degree], result[degree - source_degree])
                coefficient = lp_add(
                    coefficient,
                    lp_scale(product, sp.Rational(source_degree, degree)),
                )
        result[degree] = coefficient
    return result


def scalar_matrix_series_product(
    left: list[sp.Matrix],
    right: list[sp.Matrix],
    order: int,
) -> list[sp.Matrix]:
    dimension = left[0].rows
    result = [sp.zeros(dimension) for _ in range(order + 1)]
    for degree in range(order + 1):
        for left_degree in range(degree + 1):
            right_degree = degree - left_degree
            if left_degree < len(left) and right_degree < len(right):
                result[degree] += left[left_degree] * right[right_degree]
    return [matrix.applyfunc(clean) for matrix in result]


def polynomial_at_matrix_series(
    polynomial: Laurent,
    matrix: list[sp.Matrix],
    order: int,
) -> list[sp.Matrix]:
    if any(power < 0 for power in polynomial):
        raise ValueError("analytic polynomial expected")
    dimension = matrix[0].rows
    result = [sp.zeros(dimension) for _ in range(order + 1)]
    identity_series = [sp.eye(dimension), *[sp.zeros(dimension) for _ in range(order)]]
    powers = {0: identity_series}
    maximum = max(polynomial, default=0)
    for power in range(1, maximum + 1):
        powers[power] = scalar_matrix_series_product(powers[power - 1], matrix, order)
    for power, coefficient in polynomial.items():
        for degree in range(order + 1):
            result[degree] += coefficient * powers[power][degree]
    return [entry.applyfunc(clean) for entry in result]


def determinant_series(matrix: list[list[Series]], order: int) -> Series:
    result: Series = [{} for _ in range(order + 1)]
    for permutation in permutations(range(3)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(3)
            for j in range(i + 1, 3)
        )
        term: Series = [{0: sp.Integer(1)}, *[{} for _ in range(order)]]
        for row, column in enumerate(permutation):
            term = series_mul(term, matrix[row][column], order)
        sign = -1 if inversions % 2 else 1
        for degree in range(order + 1):
            result[degree] = lp_add(result[degree], lp_scale(term[degree], sign))
    return result


def inverse_riemann_path_sparse(
    path: list[sp.Matrix],
    order: int,
    *,
    progress: bool = False,
) -> tuple[list[Laurent], list[sp.Matrix]]:
    """Return exact inverse-map and pullback series for a real p=3 path."""

    root_two = sp.sqrt(2)
    canonical_base = sp.Matrix(
        [[0, root_two, 0], [0, 0, root_two], [0, 0, 0]]
    )
    if not path or any(coefficient.shape != (3, 3) for coefficient in path):
        raise ValueError("this exact engine requires a nonempty 3x3 path")
    if path[0] != canonical_base:
        raise ValueError("the base must be the canonical p=3 Crabb block")
    if any(
        entry.is_real is not True
        for coefficient in path
        for entry in coefficient
    ):
        raise ValueError("this exact engine requires a real symbolic path")

    zero: Series = [{} for _ in range(order + 1)]
    support: list[list[Series]] = [[list(zero) for _ in range(3)] for _ in range(3)]
    for row in range(3):
        for column in range(3):
            entry = [{} for _ in range(order + 1)]
            for degree in range(min(len(path), order + 1)):
                forward = path[degree][row, column]
                reverse = path[degree][column, row]
                if forward != 0:
                    entry[degree][-1] = entry[degree].get(-1, 0) + forward / 2
                if reverse != 0:
                    entry[degree][1] = entry[degree].get(1, 0) + reverse / 2
            support[row][column] = entry

    support_value: Series = [{0: sp.Integer(1)}, *[{} for _ in range(order)]]
    for degree in range(1, order + 1):
        characteristic_matrix: list[list[Series]] = []
        for row in range(3):
            characteristic_row = []
            for column in range(3):
                entry = [lp_scale(value, -1) for value in support[row][column]]
                if row == column:
                    entry = [
                        lp_add(entry[index], support_value[index])
                        for index in range(order + 1)
                    ]
                characteristic_row.append(entry)
            characteristic_matrix.append(characteristic_row)
        residual = determinant_series(characteristic_matrix, degree)[degree]
        support_value[degree] = lp_scale(residual, -sp.Rational(1, 2))
        if progress:
            print("support", degree, len(support_value[degree]), flush=True)

    boundary: Series = []
    for coefficient in support_value:
        boundary.append(
            {
                mode + 1: clean((1 - mode) * value)
                for mode, value in coefficient.items()
                if clean((1 - mode) * value) != 0
            }
        )

    angle_shift: Series = [{}, *[{} for _ in range(order)]]
    inverse_map: list[Laurent] = [{1: sp.Integer(1)}]

    def reparameterized_coefficient(degree: int) -> Laurent:
        result: Laurent = {}
        modes = {mode for coefficient in boundary for mode in coefficient}
        for mode in modes:
            exponent = [lp_scale(value, sp.I * mode) for value in angle_shift]
            phase = series_exp(exponent, degree)
            for boundary_degree in range(degree + 1):
                value = boundary[boundary_degree].get(mode, 0)
                if value != 0:
                    term = lp_shift(
                        lp_scale(phase[degree - boundary_degree], value),
                        mode,
                    )
                    result = lp_add(result, term)
        return result

    for degree in range(1, order + 1):
        unmatched = reparameterized_coefficient(degree)
        correction: Laurent = {}
        for mode, coefficient in unmatched.items():
            if mode <= 0:
                analytic_half = {mode - 1: sp.I * coefficient}
                correction = lp_add(correction, analytic_half)
                correction = lp_add(correction, lp_star(analytic_half))
        angle_shift[degree] = correction
        analytic = reparameterized_coefficient(degree)
        bad = {mode: value for mode, value in analytic.items() if mode <= 0 and value != 0}
        if bad:
            raise AssertionError((degree, bad))
        inverse_map.append(analytic)
        if progress:
            print("map", degree, len(analytic), flush=True)

    operator = [path[0]]
    for degree in range(1, order + 1):
        trial = [*operator, sp.zeros(3)]
        known = sp.zeros(3)
        for map_degree in range(1, degree + 1):
            evaluated = polynomial_at_matrix_series(
                inverse_map[map_degree],
                trial,
                degree - map_degree,
            )
            known += evaluated[degree - map_degree]
        target = path[degree] if degree < len(path) else sp.zeros(3)
        operator.append((target - known).applyfunc(clean))
        if progress:
            print("operator", degree, flush=True)

    return inverse_map, operator


def stein_gramian_coefficients(
    operator: list[sp.Matrix],
    defect_coefficients: list[sp.Matrix],
    order: int,
) -> list[sp.Matrix]:
    dimension = operator[0].rows
    base = operator[0]
    forcing = [sp.zeros(dimension) for _ in range(order + 1)]
    for degree in range(order + 1):
        for left_degree in range(degree + 1):
            right_degree = degree - left_degree
            if left_degree < len(defect_coefficients) and right_degree < len(defect_coefficients):
                forcing[degree] += (
                    defect_coefficients[left_degree]
                    * defect_coefficients[right_degree].T
                )
        forcing[degree] = forcing[degree].applyfunc(clean)

    gramian: list[sp.Matrix] = []
    for degree in range(order + 1):
        right_hand_side = forcing[degree]
        for left_degree in range(degree + 1):
            for right_degree in range(degree + 1 - left_degree):
                metric_degree = degree - left_degree - right_degree
                if metric_degree < len(gramian):
                    right_hand_side += (
                        operator[left_degree].T
                        * gramian[metric_degree]
                        * operator[right_degree]
                    )
        coefficient = sp.zeros(dimension)
        for power in range(dimension):
            coefficient += (
                (base.T**power)
                * right_hand_side
                * (base**power)
            )
        gramian.append(coefficient.applyfunc(clean))
    return gramian


def scalar_series_mul(left: list[sp.Expr], right: list[sp.Expr], order: int) -> list[sp.Expr]:
    result = [sp.Integer(0) for _ in range(order + 1)]
    for degree in range(order + 1):
        result[degree] = clean(
            sum(
                left[left_degree] * right[degree - left_degree]
                for left_degree in range(degree + 1)
            )
        )
    return result


def determinant_coefficients(matrix: list[list[list[sp.Expr]]], order: int) -> list[sp.Expr]:
    result = [sp.Integer(0) for _ in range(order + 1)]
    for permutation in permutations(range(3)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(3)
            for j in range(i + 1, 3)
        )
        term = [sp.Integer(1), *[sp.Integer(0) for _ in range(order)]]
        for row, column in enumerate(permutation):
            term = scalar_series_mul(term, matrix[row][column], order)
        sign = -1 if inversions % 2 else 1
        for degree in range(order + 1):
            result[degree] += sign * term[degree]
    return [clean(value) for value in result]


def eigenvalue_coefficients(
    matrix_coefficients: list[sp.Matrix],
    base_eigenvalue: int,
    derivative: int,
    order: int,
) -> list[sp.Expr]:
    eigenvalue = [sp.Integer(base_eigenvalue), *[sp.Integer(0) for _ in range(order)]]
    for degree in range(1, order + 1):
        characteristic_matrix: list[list[list[sp.Expr]]] = []
        for row in range(3):
            characteristic_row = []
            for column in range(3):
                entry = [
                    -matrix_coefficients[index][row, column]
                    for index in range(order + 1)
                ]
                if row == column:
                    entry = [
                        clean(entry[index] + eigenvalue[index])
                        for index in range(order + 1)
                    ]
                characteristic_row.append(entry)
            characteristic_matrix.append(characteristic_row)
        residual = determinant_coefficients(characteristic_matrix, degree)[degree]
        eigenvalue[degree] = clean(-residual / derivative)
    return eigenvalue


def gramian_condition_coefficients(
    operator: list[sp.Matrix],
    defect_coefficients: list[sp.Matrix],
    order: int,
) -> list[sp.Expr]:
    """Return the endpoint-eigenvalue ratio of the rank-one Stein Gramian."""

    if len(operator) < order + 1 or len(defect_coefficients) < order + 1:
        raise ValueError("operator and defect series must include the requested order")
    gramian = stein_gramian_coefficients(operator, defect_coefficients, order)
    lower = eigenvalue_coefficients(gramian, 1, 3, order)
    upper = eigenvalue_coefficients(gramian, 4, 6, order)
    ratio = [sp.Integer(4), *[sp.Integer(0) for _ in range(order)]]
    for degree in range(1, order + 1):
        ratio[degree] = clean(
            upper[degree]
            - sum(lower[index] * ratio[degree - index] for index in range(1, degree + 1))
        )
    return ratio


def optimized_rank_one_condition(
    path: list[sp.Matrix],
    order: int,
) -> tuple[tuple[sp.Expr, ...], dict[sp.Symbol, sp.Expr]]:
    """Optimize the analytic rank-one defect through ``order``."""

    _, operator = inverse_riemann_path_sparse(path, order)
    correction_count = order // 2
    a = sp.symbols(f"a1:{correction_count + 1}", real=True)
    b = sp.symbols(f"b1:{correction_count + 1}", real=True)
    defect = [sp.Matrix([1, 0, 0])]
    for degree in range(1, order + 1):
        if degree <= correction_count:
            defect.append(sp.Matrix([0, a[degree - 1], b[degree - 1]]))
        else:
            defect.append(sp.zeros(3, 1))
    condition = gramian_condition_coefficients(operator, defect, order)

    substitutions: dict[sp.Symbol, sp.Expr] = {}
    for degree in range(2, 2 * correction_count + 1, 2):
        variables = (a[degree // 2 - 1], b[degree // 2 - 1])
        coefficient = sp.factor(condition[degree].subs(substitutions))
        solutions = sp.solve(
            [sp.diff(coefficient, variable) for variable in variables],
            variables,
            dict=True,
        )
        if len(solutions) != 1:
            raise AssertionError("the formal defect correction was not unique")
        substitutions.update(solutions[0])

    optimized = tuple(
        sp.factor(coefficient.subs(substitutions))
        for coefficient in condition
    )
    return optimized, substitutions
