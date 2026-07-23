#!/usr/bin/env python3
"""Audit the finite algebra used in the local p=3 disk-center reduction."""

from __future__ import annotations

from itertools import product

import sympy as sp


def defect_hessian() -> sp.Expr:
    """Return the second rank-one-defect coefficient at the Crabb block."""

    root_two = sp.sqrt(2)
    operator = sp.Matrix(
        [[0, root_two, 0], [0, 0, root_two], [0, 0, 0]]
    )
    parameter = sp.symbols("parameter", real=True)
    x_real, x_imaginary, y_real, y_imaginary = sp.symbols(
        "x_real x_imaginary y_real y_imaginary",
        real=True,
    )
    defect = sp.Matrix(
        [
            1,
            parameter * (x_real + sp.I * x_imaginary),
            parameter * (y_real + sp.I * y_imaginary),
        ]
    )
    gramian = sp.zeros(3)
    iterate = defect
    for _ in range(3):
        gramian += iterate * iterate.conjugate().T
        iterate = operator.conjugate().T * iterate

    base_metric = sp.diag(1, 2, 4)
    linear = gramian.diff(parameter).subs(parameter, 0)
    quadratic = gramian.diff(parameter, 2).subs(parameter, 0) / 2
    lower_quadratic = quadratic[0, 0] + sum(
        linear[index, 0].conjugate() * linear[index, 0]
        / (1 - base_metric[index, index])
        for index in (1, 2)
    )
    upper_quadratic = quadratic[2, 2] + sum(
        linear[index, 2].conjugate() * linear[index, 2]
        / (4 - base_metric[index, index])
        for index in (0, 1)
    )
    return sp.factor(upper_quadratic - 4 * lower_quadratic)


def invariant_patterns(total_degree: int) -> set[tuple[int, int, int, int]]:
    """Enumerate rotation-invariant z^a zbar^b U^c Ubar^d patterns."""

    patterns = set()
    for exponents in product(range(total_degree + 1), repeat=4):
        z_degree, zbar_degree, soft_degree, softbar_degree = exponents
        if sum(exponents) != total_degree:
            continue
        if soft_degree + softbar_degree < 2:
            continue
        weight = (
            z_degree
            - zbar_degree
            + 2 * (soft_degree - softbar_degree)
        )
        if weight == 0:
            patterns.add(exponents)
    return patterns


def main() -> None:
    x_real, x_imaginary, y_real, y_imaginary = sp.symbols(
        "x_real x_imaginary y_real y_imaginary",
        real=True,
    )
    expected_defect_hessian = sp.Rational(8, 3) * (
        3 * x_real**2
        + 3 * x_imaginary**2
        + y_real**2
        + y_imaginary**2
    )
    if sp.expand(defect_hessian() - expected_defect_hessian) != 0:
        raise AssertionError("the rank-one-defect Hessian lost positivity")

    expected_patterns = {
        2: {(0, 0, 1, 1)},
        3: set(),
        4: {(0, 0, 2, 2), (1, 1, 1, 1)},
        5: {(0, 2, 2, 1), (2, 0, 1, 2)},
    }
    actual_patterns = {
        degree: invariant_patterns(degree) for degree in range(2, 6)
    }
    if actual_patterns != expected_patterns:
        raise AssertionError(f"unexpected low-order invariant patterns: {actual_patterns}")

    root_two = sp.sqrt(2)
    soft_real, strong_real, transverse_real = sp.symbols(
        "soft_real strong_real transverse_real",
        real=True,
    )
    soft_center = 3 * root_two / 64
    strong_center = -sp.Rational(9, 64)
    leading = (
        -sp.Rational(171, 4096)
        - sp.Rational(31, 8) * soft_real**2
        - sp.Rational(123, 256) * root_two * soft_real
        - sp.Rational(117, 128) * strong_real
        - 6 * root_two * soft_real * strong_real
        - sp.Rational(21, 4) * strong_real**2
        - 8 * transverse_real**2
    )
    strong_optimizer = sp.solve(
        sp.diff(leading, strong_real),
        strong_real,
    )[0]
    reduced = sp.factor(
        leading.subs(
            {
                soft_real: soft_center + soft_real,
                strong_real: strong_optimizer.subs(
                    soft_real,
                    soft_center + soft_real,
                ),
                transverse_real: 0,
            },
            simultaneous=True,
        )
    )
    expected_reduced = -sp.Rational(25, 56) * soft_real**2
    if sp.expand(reduced - expected_reduced) != 0:
        raise AssertionError("strong-variable elimination changed the soft coefficient")
    if strong_optimizer.subs(soft_real, soft_center) != strong_center:
        raise AssertionError("the leading strong optimizer missed the disk center")

    print("PASS p=3 disk-center Morse-Bott algebra")
    print(f"defect Hessian = {expected_defect_hessian}")
    print("no rotation-invariant degree-three soft term exists")
    print("degree four contains only |U|^4 and |z|^2|U|^2")
    print("strong elimination leaves -(25/56)|U|^2 at weighted order six")


if __name__ == "__main__":
    main()
