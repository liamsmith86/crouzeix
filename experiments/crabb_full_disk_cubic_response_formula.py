#!/usr/bin/env python3
"""Closed all-size cubic true-normal response after disk recentering.

The Toeplitz direction convention used throughout the Crabb disk
experiments is ``direction=(0,z_0,...,z_{L-2})``.  L176's quadratic
correction is built from

    W = z wedge J conjugate(z).

This module evaluates the closed response formula in terms of the
weighted anti-diagonal fluxes ``S_t``.  It contains no characteristic
polynomial or fitted data and is therefore suitable as the candidate
side of the independent exact audits.
"""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp


def weighted_plucker_antidiagonal(
    coefficients: Sequence[sp.Expr],
    total: int,
) -> sp.Expr:
    """Return ``S_total`` for ``z wedge J conjugate(z)``.

    The ordered form

    ``S_t = sum_{i+j=t} (j-i) z_i conjugate(z_{n-1-j})``

    is algebraically identical to the intrinsic ``i<j`` definition and
    automatically handles the valid endpoint range.
    """

    coefficient_count = len(coefficients)
    if total < 1:
        return sp.Integer(0)
    return sp.expand(
        sum(
            (
                (total - 2 * left)
                * coefficients[left]
                * sp.conjugate(
                    coefficients[coefficient_count - 1 - (total - left)]
                )
                for left in range(coefficient_count)
                if 0 <= total - left < coefficient_count
            ),
            sp.Integer(0),
        )
    )

def cubic_normal_response(
    direction: Sequence[sp.Expr],
    mode: int,
) -> sp.Expr:
    """Return the complex cubic response in one circular-normal mode.

    Modes ``3,...,L-3`` are the only possibly active modes.  The first
    sum has circle grade ``mode`` and the conjugated second sum has the
    same grade.  Coefficients outside the intrinsic L176 anti-diagonal
    range are deliberately absent rather than represented through
    exterior-algebra syzygies.
    """

    length = len(direction)
    if length < 1:
        raise ValueError("the direction must contain its zero constant term")
    if mode < 3 or mode > length - 3:
        return sp.Integer(0)

    coefficients = tuple(direction[1:])
    coefficient_count = len(coefficients)
    response = sp.Integer(0)

    for outer in range(mode, coefficient_count):
        total = length + mode - 3 - outer
        flux = weighted_plucker_antidiagonal(coefficients, total)
        response += (
            sp.Rational(outer - total - 1, total + 2)
            * coefficients[outer]
            * flux
        )

    conjugate_total = length - mode - 3
    for outer in range(conjugate_total):
        total = conjugate_total - outer
        flux = weighted_plucker_antidiagonal(coefficients, total)
        response += (
            sp.Rational(outer - total - 1, total + 2)
            * sp.conjugate(coefficients[outer] * flux)
        )

    return sp.expand(
        sp.Rational(16 * (4 * mode - 1), length**2) * response
    )
