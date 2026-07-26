"""Shared exact formulas for the three-dimensional Gau--Wu disk model."""

from __future__ import annotations

import sympy as sp


def gau_wu_matrix(parameter: sp.Expr) -> sp.Matrix:
    """Return the exact three-dimensional model with zeros 0 and parameter."""

    edge = sp.sqrt(2 * (1 - parameter**2))
    return sp.Matrix(
        [
            [0, edge, -2 * parameter],
            [0, parameter, edge],
            [0, 0, 0],
        ]
    )


def blaschke_value(matrix: sp.Matrix, parameter: sp.Expr) -> sp.Matrix:
    """Evaluate z(z-a)/(1-az) at the model matrix."""

    identity = sp.eye(matrix.rows)
    return sp.simplify(
        matrix
        * (matrix - parameter * identity)
        * (identity - parameter * matrix).inv()
    )


def parameter_from_schwarz(schwarz_parameter: sp.Expr) -> sp.Expr:
    """Rationalize a by a=2q/(1+q^2)."""

    return 2 * schwarz_parameter / (1 + schwarz_parameter**2)
