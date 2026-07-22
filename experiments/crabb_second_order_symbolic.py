"""Exact symbolic helpers for the low-order Crabb second variation.

The public entry point constructs the physical-gauge coefficients ``G,H``
from exact support eigenvectors.  Keeping this derivation shared prevents the
``p=3`` and ``p=4`` proof certificates from silently using different Fourier
or functional-calculus conventions.
"""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp


def symbolic_perturbation(
    dimension: int,
) -> tuple[sp.Matrix, tuple[sp.Symbol, ...], tuple[sp.Symbol, ...]]:
    """Return a general complex matrix with independent real coordinates."""

    real_parts = sp.symbols(f"x0:{dimension**2}", real=True)
    imaginary_parts = sp.symbols(f"y0:{dimension**2}", real=True)
    perturbation = sp.Matrix(
        dimension,
        dimension,
        lambda row, column: (
            real_parts[dimension * row + column]
            + sp.I * imaginary_parts[dimension * row + column]
        ),
    )
    return perturbation, real_parts, imaginary_parts


def conformal_coefficients(
    base: sp.Matrix,
    perturbation: sp.Matrix,
    top_vector: sp.Matrix,
    lower_spectral_data: Sequence[tuple[sp.Expr, sp.Matrix]],
    z: sp.Symbol,
) -> tuple[sp.Matrix, sp.Matrix]:
    """Derive the exact physical-gauge operator coefficients ``G,H``.

    The vectors must be normalized eigenvectors of
    ``Re(z**-1 * base)``.  ``lower_spectral_data`` stores each non-top
    eigenvalue together with its vector; the top eigenvalue is normalized to
    one for every canonical Crabb block.
    """

    dimension = base.rows

    def star(expression: sp.Expr) -> sp.Expr:
        return sp.expand(sp.conjugate(expression).subs(sp.conjugate(z), 1 / z))

    def vector_star(vector: sp.Matrix) -> sp.Matrix:
        return sp.Matrix([[star(vector[index]) for index in range(vector.rows)]])

    def coefficient(expression: sp.Expr, degree: int) -> sp.Expr:
        return sp.expand(expression).coeff(z, degree)

    support_base = (base / z + z * base.conjugate().T) / 2
    spectral_data = [(sp.Integer(1), top_vector), *lower_spectral_data]
    for eigenvalue, eigenvector in spectral_data:
        norm_residual = sp.simplify((vector_star(eigenvector) * eigenvector)[0] - 1)
        eigen_residual = support_base * eigenvector - eigenvalue * eigenvector
        if norm_residual != 0 or any(
            sp.simplify(entry) != 0 for entry in eigen_residual
        ):
            raise AssertionError("invalid exact support eigenpair")

    support_perturbation = (perturbation / z + z * perturbation.conjugate().T) / 2
    first_support = sp.expand(
        (vector_star(top_vector) * support_perturbation * top_vector)[0]
    )
    second_support = sp.Integer(0)
    for eigenvalue, eigenvector in lower_spectral_data:
        coupling = sp.expand(
            (vector_star(eigenvector) * support_perturbation * top_vector)[0]
        )
        second_support += star(coupling) * coupling / (1 - eigenvalue)
    second_support = sp.expand(second_support)

    boundary_schwarz = coefficient(first_support, 0) + sum(
        2 * coefficient(first_support, degree) * z**degree
        for degree in range(1, dimension + 1)
    )
    imaginary_boundary = sp.expand(
        (boundary_schwarz - star(boundary_schwarz)) / (2 * sp.I)
    )
    support_derivative = sum(
        sp.I * degree * coefficient(first_support, degree) * z**degree
        for degree in range(-dimension, dimension + 1)
    )
    normal_angle_shift = sp.expand(imaginary_boundary - support_derivative)
    second_normal_data = sp.expand(second_support - normal_angle_shift**2 / 2)

    # F needs degrees through 2p-1 for DF(A)[E]; K(A) only needs degrees
    # below p.  Frequencies outside the finite Laurent support vanish exactly.
    first_coefficients = [sp.Integer(0)] * (2 * dimension)
    first_coefficients[1] = coefficient(first_support, 0)
    for degree in range(2, 2 * dimension):
        first_coefficients[degree] = 2 * coefficient(
            first_support,
            degree - 1,
        )
    second_coefficients = [sp.Integer(0)] * dimension
    second_coefficients[1] = coefficient(second_normal_data, 0)
    for degree in range(2, dimension):
        second_coefficients[degree] = 2 * coefficient(
            second_normal_data,
            degree - 1,
        )

    def polynomial_at_base(coefficients: Sequence[sp.Expr]) -> sp.Matrix:
        value = sp.zeros(dimension)
        power = sp.eye(dimension)
        for scalar in coefficients:
            value += scalar * power
            power *= base
        return value

    tangent = sp.simplify(
        perturbation - polynomial_at_base(first_coefficients[:dimension])
    )
    powers = [sp.eye(dimension)]
    for _ in range(1, dimension):
        powers.append(powers[-1] * base)
    frechet_term = sp.zeros(dimension)
    for degree in range(1, 2 * dimension):
        for left_degree in range(degree):
            right_degree = degree - 1 - left_degree
            if left_degree < dimension and right_degree < dimension:
                frechet_term += first_coefficients[degree] * (
                    powers[left_degree] * perturbation * powers[right_degree]
                )

    derivative_product = [sp.Integer(0)] * dimension
    for left_degree in range(1, 2 * dimension):
        for right_degree in range(1, 2 * dimension):
            product_degree = left_degree - 1 + right_degree
            if product_degree < dimension:
                derivative_product[product_degree] += (
                    left_degree
                    * first_coefficients[left_degree]
                    * first_coefficients[right_degree]
                )
    second_order = sp.simplify(
        -frechet_term
        + polynomial_at_base(derivative_product)
        - polynomial_at_base(second_coefficients)
    )
    return tangent, second_order
