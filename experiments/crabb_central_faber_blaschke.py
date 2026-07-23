#!/usr/bin/env python3
"""Regenerate the central corrected Faber--Blaschke face.

For ``L=2k`` the corrected A99 polynomial factors through the degree-``k``
Dickson map:

    P_(2k) + 2 a (1+c^k) P_k
      = Q_(a,r)(P_k),       r=c^k,
    Q_(a,r)(y) = y^2 + 2 a (1+r)y - 2r.

The corresponding root Blaschke product therefore descends to the
three-dimensional outer block from L126.  This checker performs the
remaining size-three calculation over an exact polynomial ring.  It
uses the elliptic Riemann map through total degree four and proves

    ||b_(a,r)(T_(a,r))||^2
      = 4 - 64 a^2 r^2 - 16 r^4 + O_total(5).
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


TOTAL_ORDER = 4


@dataclass(frozen=True)
class CentralFaberBlaschkeRecord:
    """Serializable exact central-face audit."""

    total_order: int
    dickson_composition_identity: bool
    lower_eigenvalue_coefficients: tuple[str, ...]
    fourth_eigenvalue_coefficient: str
    mixed_coefficient: str
    predicted_mixed_coefficient: str
    axis_coefficient: str
    predicted_axis_coefficient: str
    root_sum_series: str
    root_product_series: str


def truncate(expression: sp.Expr, scale: sp.Symbol) -> sp.Expr:
    """Truncate a scalar expression after ``TOTAL_ORDER`` in ``scale``."""

    return sp.expand(
        sp.series(
            expression,
            scale,
            0,
            TOTAL_ORDER + 1,
        ).removeO()
    )


def truncate_matrix(matrix: sp.Matrix, scale: sp.Symbol) -> sp.Matrix:
    """Apply :func:`truncate` entrywise."""

    return matrix.applyfunc(lambda entry: truncate(entry, scale))


def dickson_polynomial(
    degree: int,
    parameter: sp.Expr,
    variable: sp.Symbol,
) -> sp.Expr:
    """Return the Dickson polynomial ``P_degree``."""

    if degree == 0:
        return sp.Integer(2)
    previous_previous = sp.Integer(2)
    previous = variable
    for _ in range(2, degree + 1):
        current = sp.expand(
            variable * previous - parameter * previous_previous
        )
        previous_previous, previous = previous, current
    return previous


def verify_dickson_composition() -> bool:
    """Check ``P_(2k)=P_k^2-2c^k`` on several symbolic degrees."""

    variable, parameter = sp.symbols("z c")
    return all(
        sp.expand(
            dickson_polynomial(2 * degree, parameter, variable)
            - dickson_polynomial(degree, parameter, variable) ** 2
            + 2 * parameter**degree
        )
        == 0
        for degree in range(1, 9)
    )


def elliptic_direct_map(
    matrix: sp.Matrix,
    parameter: sp.Expr,
) -> sp.Matrix:
    """Return the part of ``phi_parameter(matrix)`` relevant to degree four.

    L125's recurrence gives, through total parameter degree four,

        phi_r(z)
          =(1+2r^2+r^4)z-(r+3r^3)z^3
           +(r^2+5r^4)z^5-r^3z^7+r^4z^9+O(r^5).
    """

    return (
        (1 + 2 * parameter**2 + parameter**4) * matrix
        - (parameter + 3 * parameter**3) * matrix**3
        + (parameter**2 + 5 * parameter**4) * matrix**5
        - parameter**3 * matrix**7
        + parameter**4 * matrix**9
    )


def root_symmetric_series(
    amplitude: sp.Expr,
    parameter: sp.Expr,
    scale: sp.Symbol,
    variable: sp.Symbol,
) -> tuple[sp.Expr, sp.Expr]:
    """Return the sum and product of the mapped roots of ``Q_(a,r)``."""

    linear = 2 * amplitude * (1 + parameter)
    constant = -2 * parameter
    power_sums: dict[int, sp.Expr] = {
        0: sp.Integer(2),
        1: -linear,
    }

    def power_sum(degree: int) -> sp.Expr:
        while degree not in power_sums:
            index = max(power_sums) + 1
            power_sums[index] = sp.expand(
                -linear * power_sums[index - 1]
                - constant * power_sums[index - 2]
            )
        return power_sums[degree]

    mapped = (
        (1 + 2 * parameter**2 + parameter**4) * variable
        - (parameter + 3 * parameter**3) * variable**3
        + (parameter**2 + 5 * parameter**4) * variable**5
        - parameter**3 * variable**7
        + parameter**4 * variable**9
    )
    root_sum = truncate(
        (1 + 2 * parameter**2 + parameter**4) * power_sum(1)
        - (parameter + 3 * parameter**3) * power_sum(3)
        + (parameter**2 + 5 * parameter**4) * power_sum(5)
        - parameter**3 * power_sum(7)
        + parameter**4 * power_sum(9),
        scale,
    )
    squared_sum = sp.Integer(0)
    for term in sp.Add.make_args(sp.expand(mapped**2)):
        degree = int(term.as_powers_dict().get(variable, 0))
        squared_sum += term / variable**degree * power_sum(degree)
    root_product = truncate(
        (root_sum**2 - squared_sum) / 2,
        scale,
    )
    return root_sum, root_product


def inverse_matrix_series(
    matrix: sp.Matrix,
    scale: sp.Symbol,
) -> sp.Matrix:
    """Invert a matrix with identity constant term through total degree four."""

    coefficients = [
        matrix.applyfunc(
            lambda entry, degree=degree: sp.expand(entry).coeff(
                scale,
                degree,
            )
        )
        for degree in range(TOTAL_ORDER + 1)
    ]
    if coefficients[0] != sp.eye(matrix.rows):
        raise AssertionError("the denominator did not have identity constant term")
    inverse_coefficients = [sp.eye(matrix.rows)]
    for degree in range(1, TOTAL_ORDER + 1):
        inverse_coefficients.append(
            -sum(
                (
                    coefficients[index]
                    * inverse_coefficients[degree - index]
                    for index in range(1, degree + 1)
                ),
                sp.zeros(matrix.rows),
            )
        )
    return truncate_matrix(
        sum(
            (
                scale**degree * coefficient
                for degree, coefficient in enumerate(inverse_coefficients)
            ),
            sp.zeros(matrix.rows),
        ),
        scale,
    )


def top_generalized_eigenvalue_series(
    gramian: sp.Matrix,
    coordinate_metric: sp.Matrix,
    scale: sp.Symbol,
) -> tuple[sp.Expr, ...]:
    """Lift the simple top generalized eigenvalue from its Crabb endpoint."""

    eigenvalue = sp.Integer(4)
    eigenvector = sp.Matrix([0, 0, 1])
    coefficients: list[sp.Expr] = []
    for degree in range(1, TOTAL_ORDER + 1):
        eigenvalue_coefficient = sp.symbols(f"lambda_{degree}")
        first, second = sp.symbols(f"x_{degree}_0 x_{degree}_1")
        trial_eigenvalue = (
            eigenvalue + eigenvalue_coefficient * scale**degree
        )
        trial_eigenvector = (
            eigenvector + sp.Matrix([first, second, 0]) * scale**degree
        )
        residual = (
            (gramian - trial_eigenvalue * coordinate_metric)
            * trial_eigenvector
        )
        equations = [
            truncate(entry, scale).coeff(scale, degree)
            for entry in residual
        ]
        solutions = sp.solve(
            equations,
            [first, second, eigenvalue_coefficient],
            dict=True,
            simplify=False,
        )
        if len(solutions) != 1:
            raise AssertionError("the simple endpoint eigenpair did not lift")
        solution = solutions[0]
        coefficient = sp.factor(solution[eigenvalue_coefficient])
        coefficients.append(coefficient)
        eigenvalue = sp.expand(
            eigenvalue + coefficient * scale**degree
        )
        eigenvector = sp.simplify(
            eigenvector
            + sp.Matrix(
                [solution[first], solution[second], 0]
            )
            * scale**degree
        )
    return tuple(coefficients)


def make_record() -> CentralFaberBlaschkeRecord:
    """Run the exact associated central calculation."""

    scale, amplitude_symbol, parameter_symbol, variable = sp.symbols(
        "t A R z",
        real=True,
    )
    amplitude = amplitude_symbol * scale
    parameter = parameter_symbol * scale

    coefficient_operator = sp.Matrix(
        [
            [0, 2, 0],
            [parameter, 0, 1],
            [0, 2 * parameter, 0],
        ]
    ) + amplitude * sp.Matrix(
        [
            [-2 * parameter, 0, 2],
            [0, 0, 0],
            [2 * parameter, 0, -2],
        ]
    )
    coordinate_metric = sp.Matrix(
        [
            [sp.Rational(1, 2), amplitude, 0],
            [amplitude, 1, amplitude],
            [0, amplitude, sp.Rational(1, 2)],
        ]
    )
    operator = truncate_matrix(
        elliptic_direct_map(coefficient_operator, parameter),
        scale,
    )
    root_sum, root_product = root_symmetric_series(
        amplitude,
        parameter,
        scale,
        variable,
    )

    identity = sp.eye(3)
    numerator = truncate_matrix(
        operator**2 - root_sum * operator + root_product * identity,
        scale,
    )
    denominator = truncate_matrix(
        identity
        - root_sum * operator
        + root_product * operator**2,
        scale,
    )
    blaschke_image = truncate_matrix(
        numerator * inverse_matrix_series(denominator, scale),
        scale,
    )
    gramian = truncate_matrix(
        blaschke_image.T * coordinate_metric * blaschke_image,
        scale,
    )
    eigenvalue_coefficients = top_generalized_eigenvalue_series(
        gramian,
        coordinate_metric,
        scale,
    )
    if any(coefficient != 0 for coefficient in eigenvalue_coefficients[:3]):
        raise AssertionError("a term appeared below the central Newton face")
    fourth_coefficient = sp.factor(eigenvalue_coefficients[3])
    mixed_coefficient = sp.expand(fourth_coefficient).coeff(
        amplitude_symbol,
        2,
    ).coeff(parameter_symbol, 2)
    axis_coefficient = sp.expand(fourth_coefficient).coeff(
        amplitude_symbol,
        0,
    ).coeff(parameter_symbol, 4)
    if mixed_coefficient != -64:
        raise AssertionError("the central corrected coefficient was not -64")
    if axis_coefficient != -16:
        raise AssertionError("the elliptic-axis coefficient was not -16")

    return CentralFaberBlaschkeRecord(
        total_order=TOTAL_ORDER,
        dickson_composition_identity=verify_dickson_composition(),
        lower_eigenvalue_coefficients=tuple(
            str(coefficient)
            for coefficient in eigenvalue_coefficients[:3]
        ),
        fourth_eigenvalue_coefficient=str(fourth_coefficient),
        mixed_coefficient=str(mixed_coefficient),
        predicted_mixed_coefficient="-64",
        axis_coefficient=str(axis_coefficient),
        predicted_axis_coefficient="-16",
        root_sum_series=str(sp.factor(root_sum)),
        root_product_series=str(sp.factor(root_product)),
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the checker and optionally persist its JSON record."""

    args = parse_args()
    record = make_record()
    line = json.dumps(asdict(record), sort_keys=True)
    print(line, flush=True)
    if args.output is not None:
        args.output.write_text(f"{line}\n", encoding="utf-8")


if __name__ == "__main__":
    main()
