"""Shared exact Pluecker-tensor Gram helpers for the full-disk sixth face."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


TensorIndex = tuple[int, int, int]


@dataclass(frozen=True)
class RankOneFactor:
    """One nonnegative rational rank-one Gram factor."""

    indices: tuple[int, ...]
    scale: sp.Rational
    coefficients: tuple[sp.Rational, ...]


def gram_from_factors(
    factors: tuple[RankOneFactor, ...],
    size: int,
) -> sp.Matrix:
    """Assemble a rational Gram matrix from rank-one factors."""

    gram = sp.zeros(size)
    for factor in factors:
        vector = sp.zeros(size, 1)
        for index, coefficient in zip(
            factor.indices,
            factor.coefficients,
            strict=True,
        ):
            vector[index] = coefficient
        gram += factor.scale * vector * vector.T
    return gram


def pluecker_tensor_features(
    mapping: tuple[TensorIndex, ...],
    variables: tuple[sp.Symbol, ...],
    conjugate_variables: tuple[sp.Symbol, ...],
) -> tuple[sp.Matrix, sp.Matrix]:
    """Return selected ``-z_a W_ij`` coordinates and their conjugates."""

    coefficient_count = len(variables)
    if len(conjugate_variables) != coefficient_count:
        raise ValueError("variable tuples must have equal length")
    features = []
    conjugate_features = []
    for outer, left, right in mapping:
        pluecker = (
            variables[left] * conjugate_variables[coefficient_count - 1 - right]
            - variables[right] * conjugate_variables[coefficient_count - 1 - left]
        )
        conjugate_pluecker = (
            conjugate_variables[left] * variables[coefficient_count - 1 - right]
            - conjugate_variables[right] * variables[coefficient_count - 1 - left]
        )
        features.append(-variables[outer] * pluecker)
        conjugate_features.append(-conjugate_variables[outer] * conjugate_pluecker)
    return sp.Matrix(features), sp.Matrix(conjugate_features)


def positive_factors_verified(
    factors: tuple[RankOneFactor, ...],
) -> bool:
    """Check the structural conditions making every factor nonnegative."""

    return all(
        factor.scale > 0 and len(factor.indices) == len(factor.coefficients)
        for factor in factors
    )
