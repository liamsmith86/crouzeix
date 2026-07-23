#!/usr/bin/env python3
"""Audit stationarity of the L70 rank-one certificate on its disk center."""

from __future__ import annotations

import sympy as sp


def main() -> None:
    root_two = sp.sqrt(2)
    disk_eigenvalue, off_diagonal = sp.symbols(
        "disk_eigenvalue off_diagonal",
        real=True,
    )
    operator = sp.Matrix(
        [
            [0, off_diagonal, -2 * disk_eigenvalue],
            [0, disk_eigenvalue, off_diagonal],
            [0, 0, 0],
        ]
    )
    metric = sp.diag(1, 2, 4)
    defect = sp.Matrix([1, 0, 0])

    def vanishes_on_disk(expression: sp.Expr) -> bool:
        """Reduce a rational identity modulo a^2 = 2(1-l^2)."""
        numerator = sp.together(expression).as_numer_denom()[0]
        relation = off_diagonal**2 - 2 * (1 - disk_eigenvalue**2)
        remainder = sp.rem(
            sp.Poly(numerator, off_diagonal),
            sp.Poly(relation, off_diagonal),
        ).as_expr()
        return sp.factor(remainder) == 0

    stein_defect = metric - operator.T * metric * operator - defect * defect.T
    if not all(vanishes_on_disk(entry) for entry in stein_defect):
        raise AssertionError("the canonical disk metric lost its rank-one defect")

    real_entries = sp.symbols("real_entry_0:9", real=True)
    real_direction = sp.Matrix(3, 3, real_entries)
    defect_one, defect_two = sp.symbols("defect_one defect_two", real=True)
    defect_tangent = sp.Matrix([0, defect_one, defect_two])
    p00, p01, p02, p11, p12, p22 = sp.symbols(
        "p00 p01 p02 p11 p12 p22",
        real=True,
    )
    metric_tangent = sp.Matrix(
        [[p00, p01, p02], [p01, p11, p12], [p02, p12, p22]]
    )
    forcing = (
        real_direction.T * metric * operator
        + operator.T * metric * real_direction
        + defect * defect_tangent.T
        + defect_tangent * defect.T
    )
    stein_tangent = sp.expand(
        metric_tangent - operator.T * metric_tangent * operator - forcing
    )
    unknowns = (p00, p01, p02, p11, p12, p22)
    solutions = sp.solve(list(stein_tangent), unknowns, dict=True)
    if len(solutions) != 1:
        raise AssertionError("the real Stein tangent did not solve uniquely")
    condition_tangent = (p22 - 4 * p00).subs(solutions[0])

    def endpoint_functional(matrix: sp.Matrix) -> sp.Expr:
        return (
            off_diagonal * (matrix[0, 1] + matrix[1, 2])
            + 2 * disk_eigenvalue * matrix[1, 1]
            - disk_eigenvalue * matrix[0, 2]
        )

    expected_condition_tangent = 4 * endpoint_functional(real_direction)
    if not vanishes_on_disk(condition_tangent - expected_condition_tangent):
        raise AssertionError("unexpected real endpoint-eigenvalue variation")

    imaginary_entries = sp.symbols("imaginary_entry_0:9", real=True)
    imaginary_direction = sp.I * sp.Matrix(3, 3, imaginary_entries)
    imaginary_defect_one, imaginary_defect_two = sp.symbols(
        "imaginary_defect_one imaginary_defect_two",
        real=True,
    )
    imaginary_defect_tangent = sp.I * sp.Matrix(
        [0, imaginary_defect_one, imaginary_defect_two]
    )
    u00, u01, u02, u11, u12, u22 = sp.symbols(
        "u00 u01 u02 u11 u12 u22",
        real=True,
    )
    v01, v02, v12 = sp.symbols("v01 v02 v12", real=True)
    hermitian_tangent = sp.Matrix(
        [
            [u00, u01 + sp.I * v01, u02 + sp.I * v02],
            [u01 - sp.I * v01, u11, u12 + sp.I * v12],
            [u02 - sp.I * v02, u12 - sp.I * v12, u22],
        ]
    )
    imaginary_forcing = (
        imaginary_direction.conjugate().T * metric * operator
        + operator.T * metric * imaginary_direction
        + defect * imaginary_defect_tangent.conjugate().T
        + imaginary_defect_tangent * defect.T
    )
    imaginary_stein = sp.expand(
        hermitian_tangent
        - operator.T * hermitian_tangent * operator
        - imaginary_forcing
    )
    imaginary_unknowns = (
        u00,
        u01,
        u02,
        u11,
        u12,
        u22,
        v01,
        v02,
        v12,
    )
    imaginary_equations: list[sp.Expr] = []
    for entry in imaginary_stein:
        imaginary_equations.extend((sp.re(entry), sp.im(entry)))
    imaginary_solutions = sp.solve(
        imaginary_equations,
        imaginary_unknowns,
        dict=True,
    )
    if len(imaginary_solutions) != 1:
        raise AssertionError("the imaginary Stein tangent did not solve uniquely")
    imaginary_condition_tangent = (u22 - 4 * u00).subs(
        imaginary_solutions[0]
    )
    if not vanishes_on_disk(imaginary_condition_tangent):
        raise AssertionError("an imaginary operator tangent changed the endpoint ratio")

    schur_parameter = sp.symbols("schur_parameter", real=True, nonzero=True)
    rational_eigenvalue = 2 * schur_parameter / (1 + schur_parameter**2)
    rational_off_diagonal = (
        root_two
        * (1 - schur_parameter**2)
        / (1 + schur_parameter**2)
    )
    boundary_variable = sp.symbols("boundary_variable", nonzero=True)
    numerator = sp.Matrix(
        [
            1 / boundary_variable - rational_eigenvalue,
            rational_off_diagonal,
            boundary_variable - rational_eigenvalue,
        ]
    )
    conjugate_numerator = sp.Matrix(
        [
            boundary_variable - rational_eigenvalue,
            rational_off_diagonal,
            1 / boundary_variable - rational_eigenvalue,
        ]
    )
    normalization = 4 * (
        1
        - rational_eigenvalue
        * (boundary_variable + 1 / boundary_variable)
        / 2
    )
    rational_operator = operator.subs(
        {
            disk_eigenvalue: rational_eigenvalue,
            off_diagonal: rational_off_diagonal,
        }
    )
    support_operator = (
        rational_operator / boundary_variable
        + boundary_variable * rational_operator.T
    ) / 2
    if any(
        sp.simplify(entry) != 0
        for entry in (support_operator - sp.eye(3)) * numerator
    ):
        raise AssertionError("the disk support vector lost its top eigenvalue")
    if sp.factor((conjugate_numerator.T * numerator)[0] - normalization) != 0:
        raise AssertionError("the disk support vector lost its normalization")

    def positive_support_value(direction: sp.Matrix) -> sp.Expr:
        support = sp.cancel(
            (
                conjugate_numerator.T
                * (
                    direction / boundary_variable
                    + boundary_variable * direction.conjugate().T
                )
                * numerator
            )[0]
            / (2 * normalization)
        )
        integrand = sp.cancel(
            support / (boundary_variable - rational_eigenvalue)
        )
        residue_zero = sp.factor(
            sp.diff(
                sp.cancel(boundary_variable**2 * integrand),
                boundary_variable,
            ).subs(boundary_variable, 0)
        )
        residue_inner = sp.factor(
            sp.cancel(
                (boundary_variable - schur_parameter) * integrand
            ).subs(boundary_variable, schur_parameter)
        )
        residue_eigenvalue = sp.factor(
            support.subs(boundary_variable, rational_eigenvalue)
        )
        return sp.factor(residue_zero + residue_inner + residue_eigenvalue)

    rational_real_direction = real_direction
    rational_endpoint = endpoint_functional(real_direction).subs(
        {
            disk_eigenvalue: rational_eigenvalue,
            off_diagonal: rational_off_diagonal,
        }
    )
    real_positive_part = positive_support_value(rational_real_direction)
    if sp.factor(real_positive_part - rational_endpoint / 4) != 0:
        raise AssertionError(
            "the real support residues did not reproduce the endpoint functional"
        )
    # Expanding all nine imaginary symbols at once creates a needlessly large
    # intermediate expression.  Real linearity makes the nine matrix units an
    # exact and much faster audit of the same statement.
    imaginary_basis = []
    for row in range(3):
        for column in range(3):
            basis_direction = sp.zeros(3)
            basis_direction[row, column] = sp.I
            imaginary_basis.append(basis_direction)
    if any(
        sp.simplify(sp.re(positive_support_value(direction))) != 0
        for direction in imaginary_basis
    ):
        raise AssertionError("the imaginary support residues had a real endpoint component")

    power_differences = [
        endpoint_functional(operator) - 4,
        endpoint_functional(operator**2) - 2 * disk_eigenvalue,
        *(operator**3 - disk_eigenvalue * operator**2),
    ]
    if not all(vanishes_on_disk(entry) for entry in power_differences):
        raise AssertionError("the inverse-map functional-calculus moments failed")

    print("PASS p=3 disk-center tangent cancellation")
    print("P - T* P T = e1 e1*")
    print("delta(condition) = 4*endpoint_functional(Re G)")
    print("support residues reproduce the endpoint functional")
    print("the Riemann-map correction cancels every complex first variation")


if __name__ == "__main__":
    main()
