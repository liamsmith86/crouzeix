#!/usr/bin/env python3
"""Derive and certify the optimized transverse Hessian at Gau--Wu models."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp

from gau_wu_disk_model import parameter_from_schwarz


@dataclass(frozen=True)
class GauWuHessianRecord:
    """One exact specialization of the symbolic sign certificate."""

    schwarz_parameter: str
    model_parameter: str
    tangent_determinant: str
    zero_velocity_minor: str
    zero_velocity_determinant: str
    real_normal_minors: tuple[str, ...]
    imaginary_normal_minors: tuple[str, ...]
    all_checks_passed: bool


def rational_residue(
    expression: sp.Expr,
    variable: sp.Symbol,
    pole: sp.Expr,
) -> sp.Expr:
    """Return a rational residue by a short local quotient recurrence."""

    expression = sp.cancel(expression)
    numerator, denominator = sp.fraction(expression)
    numerator_poly = sp.Poly(numerator, variable, domain="EX")
    denominator_poly = sp.Poly(denominator, variable, domain="EX")
    factor = sp.Poly(variable - pole, variable, domain="EX")
    order = 0
    while True:
        quotient, remainder = denominator_poly.div(factor)
        if not remainder.is_zero:
            break
        denominator_poly = quotient
        order += 1
    if order == 0:
        return sp.Integer(0)

    numerator_coefficients = [
        sp.cancel(
            sp.diff(numerator_poly.as_expr(), variable, degree).subs(
                variable,
                pole,
            )
            / sp.factorial(degree)
        )
        for degree in range(order)
    ]
    reduced_denominator = denominator_poly.as_expr()
    denominator_coefficients = [
        sp.cancel(
            sp.diff(reduced_denominator, variable, degree).subs(
                variable,
                pole,
            )
            / sp.factorial(degree)
        )
        for degree in range(order)
    ]
    quotient_coefficients: list[sp.Expr] = []
    for degree in range(order):
        value = numerator_coefficients[degree] - sum(
            denominator_coefficients[index]
            * quotient_coefficients[degree - index]
            for index in range(1, degree + 1)
        )
        quotient_coefficients.append(
            sp.cancel(value / denominator_coefficients[0])
        )
    return quotient_coefficients[-1]


def circle_integral(
    expression: sp.Expr,
    variable: sp.Symbol,
    poles: list[sp.Expr],
) -> sp.Expr:
    """Sum all prescribed interior residues."""

    return sp.cancel(
        sum(rational_residue(expression, variable, pole) for pole in poles)
    )


def matrix_function_from_values(
    matrix: sp.Matrix,
    parameter: sp.Expr,
    value_zero: sp.Expr,
    derivative_zero: sp.Expr,
    value_parameter: sp.Expr,
) -> sp.Matrix:
    """Evaluate an analytic function using ``z^2(z-a)``."""

    quadratic = sp.cancel(
        (
            value_parameter
            - value_zero
            - parameter * derivative_zero
        )
        / parameter**2
    )
    identity = sp.eye(matrix.rows)
    return (
        value_zero * identity
        + derivative_zero * matrix
        + quadratic * matrix**2
    ).applyfunc(sp.cancel)


def function_data(
    function: sp.Expr,
    variable: sp.Symbol,
    parameter: sp.Expr,
) -> dict[str, sp.Expr]:
    """Return the Hermite data needed for a Frechet derivative."""

    return {
        "value_zero": sp.cancel(function.subs(variable, 0)),
        "derivative_zero": sp.cancel(
            sp.diff(function, variable).subs(variable, 0)
        ),
        "second_zero": sp.cancel(
            sp.diff(function, variable, 2).subs(variable, 0)
        ),
        "third_zero": sp.cancel(
            sp.diff(function, variable, 3).subs(variable, 0)
        ),
        "value_parameter": sp.cancel(function.subs(variable, parameter)),
        "derivative_parameter": sp.cancel(
            sp.diff(function, variable).subs(variable, parameter)
        ),
    }


def frechet_from_data(
    matrix: sp.Matrix,
    direction: sp.Matrix,
    parameter: sp.Expr,
    data: dict[str, sp.Expr],
) -> sp.Matrix:
    """Evaluate the Frechet derivative by a degree-five Hermite polynomial."""

    coefficients = [
        data["value_zero"],
        data["derivative_zero"],
        data["second_zero"] / 2,
        data["third_zero"] / 6,
    ]
    remainder_value = sp.cancel(
        data["value_parameter"]
        - sum(
            coefficients[degree] * parameter**degree
            for degree in range(4)
        )
    )
    remainder_derivative = sp.cancel(
        data["derivative_parameter"]
        - sum(
            degree
            * coefficients[degree]
            * parameter ** (degree - 1)
            for degree in range(1, 4)
        )
    )
    coefficients.extend(
        (
            sp.cancel(
                5 * remainder_value / parameter**4
                - remainder_derivative / parameter**3
            ),
            sp.cancel(
                remainder_derivative / parameter**4
                - 4 * remainder_value / parameter**5
            ),
        )
    )

    size = matrix.rows
    block = sp.zeros(2 * size)
    block[:size, :size] = matrix
    block[:size, size:] = direction
    block[size:, size:] = matrix
    power = sp.eye(2 * size)
    value = sp.zeros(2 * size)
    for coefficient in coefficients:
        value += coefficient * power
        power *= block
    return value[:size, size:].applyfunc(sp.cancel)


def jet_multiply(
    left: list[sp.Matrix],
    right: list[sp.Matrix],
) -> list[sp.Matrix]:
    """Multiply two matrix jets through second order."""

    size = left[0].rows
    return [
        sum(
            (
                left[index] * right[degree - index]
                for index in range(degree + 1)
            ),
            sp.zeros(size),
        )
        for degree in range(3)
    ]


def jet_inverse(value: list[sp.Matrix]) -> list[sp.Matrix]:
    """Invert a matrix jet through second order."""

    inverse_zero = value[0].inv()
    inverse_one = -inverse_zero * value[1] * inverse_zero
    inverse_two = -inverse_zero * (
        value[1] * inverse_one + value[2] * inverse_zero
    )
    return [inverse_zero, inverse_one, inverse_two]


def expected_optimized_hessian(q: sp.Symbol) -> sp.Matrix:
    """Return the closed five-normal Hessian printed in the proof note."""

    d = q**6 + 2 * q**4 - 4 * q**2 + 3
    matrix = sp.zeros(5)
    matrix[0, 0] = (
        q**10
        + 3 * q**8
        - 11 * q**6
        + 17 * q**4
        - 29 * q**2
        + 21
    ) / (16 * (q**2 - 1))
    matrix[0, 2] = -sp.sqrt(2) * q * d / 16
    matrix[0, 4] = -q**2 * d / (16 * (q**2 - 1))
    matrix[1, 1] = (
        q**10
        - 3 * q**8
        + 5 * q**6
        - q**4
        - 13 * q**2
        - 21
    ) / (16 * (q**2 + 1))
    matrix[1, 3] = (
        -sp.sqrt(2)
        * q
        * (q**4 + 1)
        * (q**4 - q**2 + 3)
        / (16 * (q**2 + 1))
    )
    matrix[2, 2] = (
        q**2
        * (q**12 - 3 * q**8 + 3 * q**4 - 64 * q**2 - 1)
        / (8 * (q**2 - 1) ** 2 * (q**2 + 1) ** 2)
    )
    matrix[2, 4] = (
        sp.sqrt(2)
        * q**3
        * (q**8 + 2 * q**6 - 2 * q**2 + 31)
        / (16 * (q**2 - 1) * (q**2 + 1) ** 2)
    )
    matrix[3, 3] = (
        q**2
        * (q**4 + 1) ** 2
        / (8 * (q**2 - 1) * (q**2 + 1))
    )
    matrix[4, 4] = (
        q**2
        * (q**8 + 3 * q**6 + 3 * q**4 - 15 * q**2 + 16)
        / (16 * (q**2 - 1) * (q**2 + 1) ** 2)
    )
    for row, column in ((0, 2), (0, 4), (1, 3), (2, 4)):
        matrix[column, row] = matrix[row, column]
    return matrix.applyfunc(sp.factor)


def derive_hessian() -> tuple[sp.Symbol, sp.Expr, sp.Matrix, sp.Matrix]:
    """Derive the zero-motion and optimized normal Hessians exactly."""

    q = sp.symbols("q", real=True, nonzero=True)
    phase, variable = sp.symbols("z w", nonzero=True)
    parameter = parameter_from_schwarz(q)
    edge = sp.sqrt(2) * (1 - q**2) / (1 + q**2)
    matrix = sp.Matrix(
        [[0, edge, -2 * parameter], [0, parameter, edge], [0, 0, 0]]
    )
    identity = sp.eye(3)

    normal_basis: list[sp.Matrix] = []
    for row, column, imaginary in (
        (2, 0, False),
        (2, 0, True),
        (2, 1, False),
        (2, 1, True),
        (2, 2, False),
    ):
        direction = sp.zeros(3)
        direction[row, column] = sp.I if imaginary else 1
        normal_basis.append(direction)

    support_matrix = (matrix / phase + phase * matrix.T) / 2
    middle = parameter * (phase + 1 / phase) / 2
    top_projection = (
        (support_matrix + identity)
        * (support_matrix - middle * identity)
        / (2 * (1 - middle))
    ).applyfunc(sp.cancel)
    bottom_projection = (
        (support_matrix - identity)
        * (support_matrix - middle * identity)
        / (2 * (1 + middle))
    )
    middle_projection = (
        (support_matrix**2 - identity) / (middle**2 - 1)
    )
    reduced_resolvent = (
        bottom_projection / 2 + middle_projection / (1 - middle)
    ).applyfunc(sp.cancel)
    support_directions = [
        (
            direction / phase
            + phase * direction.conjugate().T
        )
        / 2
        for direction in normal_basis
    ]

    first_support: list[sp.Expr] = []
    first_map: list[sp.Expr] = []
    first_data: list[dict[str, sp.Expr]] = []
    angle_shift: list[sp.Expr] = []
    for support_direction in support_directions:
        support = sp.cancel(sp.trace(top_projection * support_direction))
        integrand = (
            (phase + variable)
            * support
            / ((phase - variable) * phase)
        )
        schwarz = sp.cancel(
            sp.residue(integrand, phase, 0)
            + sp.limit((phase - q) * integrand, phase, q)
            + sp.limit(
                (phase - variable) * integrand,
                phase,
                variable,
            )
        )
        conformal_map = sp.cancel(variable * schwarz)
        shift = sp.cancel(
            -sp.I
            * (
                schwarz.subs(variable, phase)
                - support
                + phase * sp.diff(support, phase)
            )
        )
        first_support.append(support)
        first_map.append(conformal_map)
        first_data.append(
            function_data(conformal_map, variable, parameter)
        )
        angle_shift.append(shift)

    second_map_data: dict[tuple[int, int], tuple[sp.Expr, sp.Expr]] = {}
    for row in range(5):
        for column in range(row, 5):
            second_support = sp.cancel(
                sp.trace(
                    top_projection
                    * (
                        support_directions[row]
                        * reduced_resolvent
                        * support_directions[column]
                        + support_directions[column]
                        * reduced_resolvent
                        * support_directions[row]
                    )
                )
                / 2
            )
            second_normal = sp.cancel(
                second_support
                - angle_shift[row] * angle_shift[column] / 2
            )
            _, second_denominator = sp.fraction(second_normal)
            if sp.simplify(second_denominator.subs(phase, -q)) == 0:
                raise RuntimeError("unlisted interior pole at -q")
            derivative_zero = circle_integral(
                second_normal / phase,
                phase,
                [sp.Integer(0), q],
            )
            weighted = (
                (phase + parameter)
                * second_normal
                / ((phase - parameter) * phase)
            )
            value_parameter = sp.cancel(
                parameter
                * circle_integral(
                    weighted,
                    phase,
                    [sp.Integer(0), q, parameter],
                )
            )
            second_map_data[row, column] = (
                derivative_zero,
                value_parameter,
            )

    def evaluate_function(function: sp.Expr) -> sp.Matrix:
        return matrix_function_from_values(
            matrix,
            parameter,
            sp.cancel(function.subs(variable, 0)),
            sp.cancel(sp.diff(function, variable).subs(variable, 0)),
            sp.cancel(function.subs(variable, parameter)),
        )

    map_at_matrix = [evaluate_function(function) for function in first_map]
    tangent = [
        (normal_basis[index] - map_at_matrix[index]).applyfunc(sp.cancel)
        for index in range(5)
    ]
    frechet = [
        [
            frechet_from_data(
                matrix,
                normal_basis[column],
                parameter,
                first_data[row],
            )
            for column in range(5)
        ]
        for row in range(5)
    ]
    second_operator: dict[tuple[int, int], sp.Matrix] = {}
    for row in range(5):
        for column in range(row, 5):
            product = sp.cancel(
                (
                    sp.diff(first_map[row], variable)
                    * first_map[column]
                    + sp.diff(first_map[column], variable)
                    * first_map[row]
                )
                / 2
            )
            derivative_zero, value_parameter = second_map_data[row, column]
            second_conformal = matrix_function_from_values(
                matrix,
                parameter,
                sp.Integer(0),
                derivative_zero,
                value_parameter,
            )
            second_operator[row, column] = (
                -(frechet[row][column] + frechet[column][row]) / 2
                + evaluate_function(product)
                - second_conformal
            ).applyfunc(sp.cancel)

    normal_variables = sp.symbols("e0:5", real=True)
    zero_variables = sp.symbols("u0:4", real=True)
    first_operator = sp.zeros(3)
    second_operator_total = sp.zeros(3)
    for row in range(5):
        first_operator += normal_variables[row] * tangent[row]
        second_operator_total += (
            normal_variables[row] ** 2 * second_operator[row, row]
        )
        for column in range(row + 1, 5):
            second_operator_total += (
                2
                * normal_variables[row]
                * normal_variables[column]
                * second_operator[row, column]
            )

    def blaschke_jet(
        zero: sp.Expr,
        velocity: sp.Expr,
    ) -> list[sp.Matrix]:
        numerator = [
            matrix - zero * identity,
            first_operator - velocity * identity,
            second_operator_total,
        ]
        denominator = [
            identity - sp.conjugate(zero) * matrix,
            -sp.conjugate(zero) * first_operator
            - sp.conjugate(velocity) * matrix,
            -sp.conjugate(zero) * second_operator_total
            - sp.conjugate(velocity) * first_operator,
        ]
        return jet_multiply(numerator, jet_inverse(denominator))

    first_zero_velocity = zero_variables[0] + sp.I * zero_variables[1]
    second_zero_velocity = zero_variables[2] + sp.I * zero_variables[3]
    image_jet = jet_multiply(
        blaschke_jet(sp.Integer(0), first_zero_velocity),
        blaschke_jet(parameter, second_zero_velocity),
    )
    first_image, second_image = image_jet[1], image_jet[2]
    first_norm = sp.expand_complex(sp.re(first_image[0, 2]))
    if sp.simplify(first_norm) != 0:
        raise RuntimeError("joint first variation did not vanish")
    second_norm = (
        4 * sp.re(second_image[0, 2]).expand(complex=True)
        + sum(
            sp.conjugate(first_image[row, 2]) * first_image[row, 2]
            for row in range(3)
        )
        + sum(
            sp.conjugate(first_image[0, column])
            * first_image[0, column]
            for column in range(2)
        )
    ) / 4
    second_norm = sp.expand_complex(second_norm)

    variables = list(normal_variables) + list(zero_variables)
    joint = sp.zeros(9)
    for row in range(9):
        for column in range(row, 9):
            entry = sp.cancel(
                sp.diff(
                    second_norm,
                    variables[row],
                    variables[column],
                )
                / 2
            )
            joint[row, column] = entry
            joint[column, row] = entry
    zero_block = joint[5:, 5:]
    optimized = (
        joint[:5, :5]
        - joint[:5, 5:] * zero_block.inv() * joint[5:, :5]
    ).applyfunc(sp.factor)
    return q, parameter, zero_block, optimized


def real_vectorize(matrix: sp.Matrix) -> sp.Matrix:
    """Vectorize a complex matrix into interlaced real coordinates."""

    return sp.Matrix(
        [
            coordinate
            for entry in matrix
            for coordinate in (
                sp.re(entry).expand(complex=True),
                sp.im(entry).expand(complex=True),
            )
        ]
    )


def tangent_determinant(q: sp.Symbol) -> sp.Expr:
    """Return a determinant proving the equality tangent has dimension 13."""

    parameter = parameter_from_schwarz(q)
    edge = sp.sqrt(2) * (1 - q**2) / (1 + q**2)
    matrix = sp.Matrix(
        [[0, edge, -2 * parameter], [0, parameter, edge], [0, 0, 0]]
    )
    tangents: list[sp.Matrix] = []
    for index in (0, 1):
        generator = sp.zeros(3)
        generator[index, index] = sp.I
        tangents.append(generator * matrix - matrix * generator)
    for row in range(3):
        for column in range(row + 1, 3):
            real_generator = sp.zeros(3)
            real_generator[row, column] = 1
            real_generator[column, row] = -1
            tangents.append(
                real_generator * matrix - matrix * real_generator
            )
            imaginary_generator = sp.zeros(3)
            imaginary_generator[row, column] = sp.I
            imaginary_generator[column, row] = sp.I
            tangents.append(
                imaginary_generator * matrix
                - matrix * imaginary_generator
            )
    identity = sp.eye(3)
    tangents.extend((identity, sp.I * identity, matrix, sp.I * matrix))
    tangents.append(sp.diff(matrix, q))

    full_columns = [real_vectorize(tangent) for tangent in tangents]
    for direction in (
        (2, 0, False),
        (2, 0, True),
        (2, 1, False),
        (2, 1, True),
        (2, 2, False),
    ):
        row, column, imaginary = direction
        complement = sp.zeros(3)
        complement[row, column] = sp.I if imaginary else 1
        full_columns.append(real_vectorize(complement))
    return sp.factor(sp.Matrix.hstack(*full_columns).det())


def bernstein_coefficients(
    polynomial: sp.Expr,
    variable: sp.Symbol,
) -> tuple[sp.Expr, ...]:
    """Convert a power polynomial to Bernstein coefficients on ``[0,1]``."""

    degree = sp.degree(polynomial, variable)
    power = sp.Poly(polynomial, variable).all_coeffs()[::-1]
    return tuple(
        sp.simplify(
            sum(
                power[index]
                * sp.binomial(order, index)
                / sp.binomial(degree, index)
                for index in range(order + 1)
            )
        )
        for order in range(degree + 1)
    )


def audit_records() -> list[GauWuHessianRecord]:
    """Regenerate the Hessian and return exact rational audit records."""

    q, parameter, zero_block, optimized = derive_hessian()
    expected = expected_optimized_hessian(q)
    if any(sp.simplify(entry) != 0 for entry in optimized - expected):
        raise RuntimeError("optimized Hessian identity failed")

    repeated_zero_block = sp.Matrix(
        [
            [
                -(q**4 + 4 * q**2 + 1) / (q**2 + 1) ** 2,
                -1,
            ],
            [
                -1,
                -(q**2 + 1) ** 2 / (q**2 - 1) ** 2,
            ],
        ]
    )
    expected_zero = sp.Matrix(
        [
            [repeated_zero_block[0, 0], 0, -1, 0],
            [0, repeated_zero_block[0, 0], 0, -1],
            [-1, 0, repeated_zero_block[1, 1], 0],
            [0, -1, 0, repeated_zero_block[1, 1]],
        ]
    )
    if any(sp.simplify(entry) != 0 for entry in zero_block - expected_zero):
        raise RuntimeError("zero-velocity Hessian identity failed")

    determinant = tangent_determinant(q)
    expected_determinant = (
        -2048
        * q**2
        * (q - 1)
        * (q + 1)
        * (q**4 + 1) ** 2
        / (q**2 + 1) ** 8
    )
    if sp.simplify(determinant - expected_determinant) != 0:
        raise RuntimeError("equality-tangent determinant failed")

    square = sp.symbols("t", real=True)
    sign_polynomials = (
        square**5
        + 3 * square**4
        - 11 * square**3
        + 17 * square**2
        - 29 * square
        + 21,
        square**9
        - 5 * square**8
        + 12 * square**7
        + 12 * square**6
        + 22 * square**5
        - 150 * square**4
        + 284 * square**3
        - 484 * square**2
        + 337 * square
        + 3,
        square**3 - 3 * square**2 + 7 * square + 3,
        -square**5
        + 3 * square**4
        - 5 * square**3
        + square**2
        + 13 * square
        + 21,
    )
    for polynomial in sign_polynomials:
        coefficients = bernstein_coefficients(polynomial, square)
        if not all(coefficient > 0 for coefficient in coefficients):
            raise RuntimeError("Bernstein positivity certificate failed")

    real_block = optimized.extract([0, 2, 4], [0, 2, 4])
    imaginary_block = optimized.extract([1, 3], [1, 3])
    records: list[GauWuHessianRecord] = []
    for value in (
        Fraction(1, 5),
        Fraction(1, 3),
        Fraction(1, 2),
        Fraction(2, 3),
        Fraction(4, 5),
    ):
        exact_q = sp.Rational(value.numerator, value.denominator)
        real_minors = tuple(
            sp.factor((-real_block[:size, :size]).det()).subs(q, exact_q)
            for size in range(1, 4)
        )
        imaginary_minors = tuple(
            sp.factor((-imaginary_block[:size, :size]).det()).subs(
                q,
                exact_q,
            )
            for size in range(1, 3)
        )
        zero_minor = sp.factor(-repeated_zero_block[0, 0]).subs(
            q,
            exact_q,
        )
        zero_determinant = sp.factor(repeated_zero_block.det()).subs(
            q,
            exact_q,
        )
        tangent_value = determinant.subs(q, exact_q)
        checks = (
            tangent_value != 0
            and zero_minor > 0
            and zero_determinant > 0
            and all(minor > 0 for minor in real_minors)
            and all(minor > 0 for minor in imaginary_minors)
        )
        if not checks:
            raise RuntimeError(f"specialized sign check failed at q={value}")
        records.append(
            GauWuHessianRecord(
                schwarz_parameter=str(value),
                model_parameter=str(parameter.subs(q, exact_q)),
                tangent_determinant=str(tangent_value),
                zero_velocity_minor=str(zero_minor),
                zero_velocity_determinant=str(zero_determinant),
                real_normal_minors=tuple(map(str, real_minors)),
                imaginary_normal_minors=tuple(map(str, imaginary_minors)),
                all_checks_passed=True,
            )
        )
    return records


def write_records(
    records: list[GauWuHessianRecord],
    output: Path,
) -> str:
    """Write JSON Lines atomically and return its SHA-256."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)
    return hashlib.sha256(output.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/gau_wu_transverse_hessian_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the exact derivation and report its data hash."""

    arguments = parse_args()
    records = audit_records()
    digest = write_records(records, arguments.output)
    print(
        json.dumps(
            {
                "all_checks_passed": all(
                    record.all_checks_passed for record in records
                ),
                "record_count": len(records),
                "sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
