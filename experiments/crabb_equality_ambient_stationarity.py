#!/usr/bin/env python3
"""Audit ambient first-order stationarity on the Crabb equality ridge.

The exact part checks the residue, transfer-function, and boundary-kernel
identities used in L162 on deterministic Gaussian-rational anchors.  The
numerical part independently constructs the first Riemann correction and
checks every real and imaginary matrix unit.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, solve_discrete_lyapunov
import sympy as sp


FFT_SIZE = 4096
SCHWARZ_DEGREE = 512


@dataclass(frozen=True)
class ExactRecord:
    """One exact algebraic regeneration."""

    dimension: int
    length: int
    stein_identity: bool
    transfer_identity: bool
    boundary_kernel_identity: bool
    residue_identity: bool


@dataclass(frozen=True)
class NumericalRecord:
    """One full ambient first-variation audit."""

    dimension: int
    length: int
    condition: float
    spectral_radius: float
    raw_gradient_norm: float
    pulled_gradient_maximum: float
    pulled_gradient_norm: float
    pulled_lower_endpoint_maximum: float
    pulled_upper_endpoint_maximum: float


def phase_palindromic_coefficients(
    length: int,
    *,
    exact: bool,
) -> list[complex | sp.Expr]:
    """Return a small deterministic phase-one positive anchor."""

    zero = sp.Integer(0) if exact else 0j
    coefficients: list[complex | sp.Expr] = [
        zero for _ in range(length)
    ]
    for offset in range(1, length):
        reflected = length - offset
        if offset < reflected:
            if exact:
                value = (
                    sp.Rational(offset + 1, 100 * length)
                    + sp.I
                    * sp.Rational(offset + 2, 140 * length)
                )
            else:
                value = (
                    (offset + 1) / (100 * length)
                    + 1j * (offset + 2) / (140 * length)
                )
            coefficients[offset] = value
            coefficients[reflected] = (
                sp.conjugate(value) if exact else np.conj(value)
            )
        elif offset == reflected:
            value = (
                sp.Rational(offset + 1, 100 * length)
                if exact
                else (offset + 1) / (100 * length)
            )
            coefficients[offset] = value
    return coefficients


def exact_polynomial_value(
    operator: sp.Matrix,
    polynomial: sp.Expr,
    variable: sp.Symbol,
) -> sp.Matrix:
    """Evaluate a scalar polynomial by Horner's rule."""

    value = sp.zeros(operator.rows)
    identity = sp.eye(operator.rows)
    for coefficient in sp.Poly(polynomial, variable).all_coeffs():
        value = value * operator + coefficient * identity
    return value


def exact_model(
    dimension: int,
) -> tuple[
    list[sp.Expr],
    sp.Matrix,
    sp.Matrix,
    sp.Matrix,
    sp.Matrix,
]:
    """Build the coefficient-coordinate L123 model exactly."""

    length = dimension - 1
    coefficients = phase_palindromic_coefficients(
        length,
        exact=True,
    )
    gram_block = sp.zeros(dimension)
    for index in range(length):
        gram_block[index, index] = sp.Rational(1, 2)
    for offset in range(1, length):
        for row in range(length - offset):
            gram_block[row, row + offset] = coefficients[offset]
            gram_block[row + offset, row] = sp.conjugate(
                coefficients[offset]
            )
    shift = sp.zeros(dimension)
    for index in range(length):
        shift[index, index + 1] = 1
    coordinate_metric = (
        gram_block
        + shift.conjugate().T * gram_block * shift
    )
    operator = sp.simplify(
        2
        * coordinate_metric.inv()
        * gram_block
        * shift
    )
    defect = gram_block[:, 0]
    reflected_defect = sp.Matrix(
        [
            sp.conjugate(defect[length - index])
            for index in range(dimension)
        ]
    )
    stein_metric = (
        coordinate_metric
        - defect * defect.conjugate().T
        + 2
        * reflected_defect
        * reflected_defect.conjugate().T
    )
    return (
        coefficients,
        coordinate_metric,
        operator,
        defect,
        stein_metric,
    )


def exact_record(dimension: int) -> ExactRecord:
    """Check all finite algebraic shadows of the L162 proof."""

    length = dimension - 1
    (
        coefficients,
        coordinate_metric,
        operator,
        defect,
        stein_metric,
    ) = exact_model(dimension)
    variable, transfer, boundary, reciprocal = sp.symbols(
        "variable transfer boundary reciprocal"
    )
    g_sharp = 1 + 2 * sum(
        coefficients[offset] * variable**offset
        for offset in range(1, length)
    )
    g = variable**length + 2 * sum(
        coefficients[offset] * variable**offset
        for offset in range(1, length)
    )

    stein_residual = sp.simplify(
        stein_metric
        - operator.conjugate().T * stein_metric * operator
        - defect * defect.conjugate().T
    )
    if stein_residual != sp.zeros(dimension):
        raise AssertionError("the exact Stein identity changed")

    blaschke_derivative = sp.cancel(
        sp.diff(g / g_sharp, variable)
    )
    numerator, denominator = sp.fraction(blaschke_derivative)
    derivative_at_operator = (
        exact_polynomial_value(operator, numerator, variable)
        * exact_polynomial_value(operator, denominator, variable).inv()
    )
    endpoint = sp.eye(dimension)[:, length]
    transfer_left = sp.factor(
        (
            defect.conjugate().T
            * derivative_at_operator
            * operator
            * (sp.eye(dimension) - transfer * operator).inv()
            * endpoint
        )[0]
    )
    transfer_denominator = 1 + 2 * sum(
        sp.conjugate(coefficients[offset]) * transfer**offset
        for offset in range(1, length)
    )
    transfer_numerator = length + 2 * sum(
        (length - offset)
        * sp.conjugate(coefficients[offset])
        * transfer**offset
        for offset in range(1, length)
    )
    if sp.cancel(
        transfer_left
        - transfer_numerator / transfer_denominator
    ) != 0:
        raise AssertionError("the endpoint transfer identity changed")

    boundary_g_sharp = g_sharp.subs(variable, boundary)
    boundary_g = g.subs(variable, boundary)
    boundary_numerator = length + 2 * sum(
        (length - offset)
        * coefficients[offset]
        * boundary**offset
        for offset in range(1, length)
    )
    conjugate_numerator = length + 2 * sum(
        (length - offset)
        * sp.conjugate(coefficients[offset])
        * boundary ** (-offset)
        for offset in range(1, length)
    )
    conjugate_g_sharp = boundary ** (-length) * boundary_g
    support_vector = sp.Matrix(
        [boundary**index for index in range(dimension)]
    )
    support_row = sp.Matrix(
        [[boundary ** (-index) for index in range(dimension)]]
    )
    support_denominator = sp.expand(
        (support_row * coordinate_metric * support_vector)[0]
    )
    kernel_residual = sp.cancel(
        boundary_numerator / boundary_g_sharp
        + conjugate_numerator / conjugate_g_sharp
        - length
        - support_denominator
        / (boundary_g_sharp * conjugate_g_sharp)
    )
    if kernel_residual != 0:
        raise AssertionError("the boundary kernel identity changed")

    polynomial_row = sp.Matrix(
        [
            [
                boundary ** (length - index)
                for index in range(dimension)
            ]
        ]
    )
    first_endpoint = sp.eye(dimension)[:, 0]
    last_endpoint = sp.eye(dimension)[:, length]
    c_vector = support_vector + boundary_g_sharp * first_endpoint
    d_row = polynomial_row + boundary_g_sharp * last_endpoint.T
    residue_numerator = sp.simplify(
        c_vector * d_row
        - support_vector * polynomial_row
    )
    expected_numerator = (
        first_endpoint * polynomial_row
        + support_vector * last_endpoint.T
        + boundary_g_sharp
        * first_endpoint
        * last_endpoint.T
    )
    if sp.simplify(
        residue_numerator / boundary_g_sharp
        - expected_numerator
    ) != sp.zeros(dimension):
        raise AssertionError("the residue numerator did not telescope")
    for row in range(dimension):
        for column in range(dimension):
            integrand = sp.cancel(
                residue_numerator[row, column]
                / (
                    boundary**2
                    * boundary_g
                    * boundary_g_sharp
                )
            )
            residue_at_infinity = sp.residue(
                -integrand.subs(boundary, 1 / reciprocal)
                / reciprocal**2,
                reciprocal,
                0,
            )
            if sp.simplify(residue_at_infinity) != 0:
                raise AssertionError(
                    "the boundary residue identity changed"
                )

    return ExactRecord(
        dimension=dimension,
        length=length,
        stein_identity=True,
        transfer_identity=True,
        boundary_kernel_identity=True,
        residue_identity=True,
    )


def numerical_model(
    dimension: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build a floating-point equality anchor and its Stein metric."""

    length = dimension - 1
    coefficients = phase_palindromic_coefficients(
        length,
        exact=False,
    )
    gram_block = np.zeros((dimension, dimension), dtype=complex)
    gram_block[:length, :length] = np.eye(length) / 2
    for offset in range(1, length):
        for row in range(length - offset):
            gram_block[row, row + offset] = coefficients[offset]
            gram_block[row + offset, row] = np.conj(
                coefficients[offset]
            )
    shift = np.zeros_like(gram_block)
    for index in range(length):
        shift[index, index + 1] = 1
    coordinate_metric = (
        gram_block
        + shift.conj().T @ gram_block @ shift
    )
    operator = 2 * np.linalg.solve(
        coordinate_metric,
        gram_block @ shift,
    )
    defect = gram_block[:, 0]
    stein_metric = solve_discrete_lyapunov(
        operator.conj().T,
        np.outer(defect, defect.conj()),
    )
    return operator, coordinate_metric, stein_metric


def matrix_polynomial(
    operator: np.ndarray,
    coefficients: np.ndarray,
) -> np.ndarray:
    """Evaluate an increasing-coefficient polynomial."""

    value = np.zeros_like(operator)
    identity = np.eye(len(operator), dtype=complex)
    for coefficient in coefficients[::-1]:
        value = value @ operator + coefficient * identity
    return value


def first_riemann_correction(
    operator: np.ndarray,
    coordinate_metric: np.ndarray,
    perturbation: np.ndarray,
) -> np.ndarray:
    """Return h_E(A) from the boundary support derivative."""

    dimension = len(operator)
    indices = np.arange(dimension)
    support_derivative = np.empty(FFT_SIZE)
    for sample in range(FFT_SIZE):
        boundary = np.exp(2j * np.pi * sample / FFT_SIZE)
        vector = boundary**indices
        denominator = np.real(
            np.vdot(vector, coordinate_metric @ vector)
        )
        support_derivative[sample] = np.real(
            np.conj(boundary)
            * np.vdot(
                vector,
                coordinate_metric @ perturbation @ vector,
            )
            / denominator
        )
    fourier = np.fft.fft(support_derivative) / FFT_SIZE
    coefficients = np.zeros(SCHWARZ_DEGREE + 2, dtype=complex)
    coefficients[1] = fourier[0]
    coefficients[2:] = 2 * fourier[1 : SCHWARZ_DEGREE + 1]
    return matrix_polynomial(operator, coefficients)


def endpoint_derivatives(
    operator: np.ndarray,
    stein_metric: np.ndarray,
    perturbation: np.ndarray,
    eigenvectors: np.ndarray,
) -> tuple[float, float]:
    """Differentiate the two generalized metric endpoints."""

    forcing = (
        perturbation.conj().T @ stein_metric @ operator
        + operator.conj().T @ stein_metric @ perturbation
    )
    metric_derivative = solve_discrete_lyapunov(
        operator.conj().T,
        forcing,
    )
    lower_vector = eigenvectors[:, 0]
    upper_vector = eigenvectors[:, -1]
    lower = float(
        np.real(np.vdot(lower_vector, metric_derivative @ lower_vector))
    )
    upper = float(
        np.real(np.vdot(upper_vector, metric_derivative @ upper_vector))
    )
    return lower, upper


def numerical_record(dimension: int) -> NumericalRecord:
    """Check every real ambient coordinate after Riemann pullback."""

    operator, coordinate_metric, stein_metric = numerical_model(
        dimension
    )
    eigenvalues, eigenvectors = eigh(
        stein_metric,
        coordinate_metric,
    )
    lower_value = float(eigenvalues[0])
    upper_value = float(eigenvalues[-1])
    raw_derivatives: list[float] = []
    pulled_derivatives: list[float] = []
    pulled_lower: list[float] = []
    pulled_upper: list[float] = []
    for scalar in (1, 1j):
        for row in range(dimension):
            for column in range(dimension):
                perturbation = np.zeros_like(operator)
                perturbation[row, column] = scalar
                raw_lower, raw_upper = endpoint_derivatives(
                    operator,
                    stein_metric,
                    perturbation,
                    eigenvectors,
                )
                raw_derivatives.append(
                    raw_upper / lower_value
                    - upper_value
                    * raw_lower
                    / lower_value**2
                )
                pulled = (
                    perturbation
                    - first_riemann_correction(
                        operator,
                        coordinate_metric,
                        perturbation,
                    )
                )
                lower, upper = endpoint_derivatives(
                    operator,
                    stein_metric,
                    pulled,
                    eigenvectors,
                )
                pulled_lower.append(lower)
                pulled_upper.append(upper)
                pulled_derivatives.append(
                    upper / lower_value
                    - upper_value * lower / lower_value**2
                )
    raw_array = np.asarray(raw_derivatives)
    pulled_array = np.asarray(pulled_derivatives)
    return NumericalRecord(
        dimension=dimension,
        length=dimension - 1,
        condition=upper_value / lower_value,
        spectral_radius=float(
            np.max(np.abs(np.linalg.eigvals(operator)))
        ),
        raw_gradient_norm=float(np.linalg.norm(raw_array)),
        pulled_gradient_maximum=float(
            np.max(np.abs(pulled_array))
        ),
        pulled_gradient_norm=float(np.linalg.norm(pulled_array)),
        pulled_lower_endpoint_maximum=float(
            np.max(np.abs(pulled_lower))
        ),
        pulled_upper_endpoint_maximum=float(
            np.max(np.abs(pulled_upper))
        ),
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=3)
    parser.add_argument("--maximum-exact-size", type=int, default=6)
    parser.add_argument("--maximum-numerical-size", type=int, default=8)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run exact and independent floating-point audits."""

    args = parse_args()
    if args.minimum_size < 3:
        raise ValueError("minimum size must be at least three")
    maximum_size = max(
        args.maximum_exact_size,
        args.maximum_numerical_size,
    )
    if maximum_size < args.minimum_size:
        raise ValueError("at least one maximum must reach the minimum")

    output = (
        args.output.open("w", encoding="utf-8")
        if args.output is not None
        else None
    )
    try:
        for dimension in range(args.minimum_size, maximum_size + 1):
            if dimension <= args.maximum_exact_size:
                exact = exact_record(dimension)
                line = json.dumps(
                    {"kind": "exact", **asdict(exact)},
                    sort_keys=True,
                )
                print(line, flush=True)
                if output is not None:
                    output.write(line + "\n")
                    output.flush()
            if dimension <= args.maximum_numerical_size:
                numerical = numerical_record(dimension)
                line = json.dumps(
                    {"kind": "numerical", **asdict(numerical)},
                    sort_keys=True,
                )
                print(line, flush=True)
                if output is not None:
                    output.write(line + "\n")
                    output.flush()
    finally:
        if output is not None:
            output.close()


if __name__ == "__main__":
    main()
