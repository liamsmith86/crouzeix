#!/usr/bin/env python3
"""Exact diagnostics for the reconstructed scalar power lemma.

The accompanying proof establishes all powers and dimensions analytically.
This script checks its algebra on a nonnormal rational example, and checks
that dropping either essential hypothesis admits an explicit counterexample.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


def require_zero(value: sp.Expr | sp.MatrixBase, label: str) -> None:
    """Fail rather than round a purported exact identity."""

    if isinstance(value, sp.MatrixBase):
        passed = all(sp.cancel(entry) == 0 for entry in value)
    else:
        passed = sp.cancel(value) == 0
    if not passed:
        raise AssertionError(f"{label}: nonzero exact residual {value}")


def audit(maximum_power: int) -> dict[str, object]:
    """Check the lemma algebra with no floating-point tolerance."""

    identity = sp.eye(2)
    nilpotent = sp.Matrix([[0, 1], [0, 0]])
    c, d = sp.Rational(1, 3), sp.Rational(4, 5)
    contraction = c * identity + d * nilpotent
    operator = c * identity + 2 * d * nilpotent
    kappa = sp.Rational(5, 3)
    unnormalized_vector = sp.Matrix([1, 5])
    normalization_square = 26
    curvature = kappa * (kappa - 1)

    def pairing(left: sp.MatrixBase, right: sp.MatrixBase) -> sp.Expr:
        # All entries are real rational in this exact diagnostic.
        return (right.T * left)[0] / normalization_square

    def squared_norm(vector: sp.MatrixBase) -> sp.Expr:
        return pairing(vector, vector)

    contraction_slack = identity - contraction.T * contraction
    assert contraction_slack[0, 0] > 0
    assert contraction_slack.det() > 0
    require_zero(
        operator.T * operator * unnormalized_vector - kappa**2 * unnormalized_vector,
        "top singular vector",
    )
    other_singular_square = sp.trace(operator.T * operator) - kappa**2
    assert 0 <= other_singular_square <= kappa**2

    defect_vector = (contraction.T * operator - kappa * identity) * unnormalized_vector
    defect_square = squared_norm(defect_vector)

    moments = {}
    for power in range(1, maximum_power + 2):
        correction = 2 * contraction.T**power - operator.T**power
        require_zero(correction - c**power * identity, f"E_{power}")
        require_zero(
            correction * operator - operator * correction,
            f"commutation {power}",
        )
        moments[power] = pairing(
            correction * operator**power * unnormalized_vector,
            unnormalized_vector,
        )

    records = []
    accumulated_bound = 0
    for power in range(1, maximum_power + 1):
        response = (
            2 * contraction.T ** (power + 1) * operator
            - 2 * kappa * contraction.T**power
        ) * unnormalized_vector
        adjoint_orbit = operator.T**power * unnormalized_vector
        recurrence = kappa * moments[power] - moments[power + 1]
        algebraic_expression = curvature * squared_norm(adjoint_orbit) - pairing(
            response, adjoint_orbit
        )
        completed_square = curvature * squared_norm(
            adjoint_orbit - response / (2 * curvature)
        ) - squared_norm(response) / (4 * curvature)
        require_zero(recurrence - algebraic_expression, "recurrence")
        require_zero(recurrence - completed_square, "square completion")
        require_zero(
            response - 2 * contraction.T**power * defect_vector,
            "common-contraction factorization",
        )
        assert squared_norm(response) <= 4 * defect_square
        assert recurrence >= -defect_square / curvature

        accumulated_bound -= kappa ** (-power) * defect_square / curvature
        finite_lower_bound = kappa ** (-power) * moments[power + 1] + accumulated_bound
        assert moments[1] >= finite_lower_bound
        records.append(
            {
                "power": power,
                "recurrence": str(recurrence),
                "square_completion_residual": "0",
                "response_square": str(squared_norm(response)),
                "finite_iteration_margin": str(moments[1] - finite_lower_bound),
            }
        )

    defect_upper_bound = 2 * kappa**2 - kappa * moments[1] - kappa**3
    assert defect_square <= defect_upper_bound
    assert moments[1] >= -defect_square / (kappa * (kappa - 1) ** 2)
    assert defect_square * (1 - 1 / (kappa - 1) ** 2) <= kappa**2 * (2 - kappa)

    # These explicitly violate the removed hypothesis, not the lemma.
    large_nilpotent = 3 * nilpotent
    bounded_noncommuting_correction = -large_nilpotent.T
    commutator = (
        bounded_noncommuting_correction * large_nilpotent
        - large_nilpotent * bounded_noncommuting_correction
    )
    assert commutator != sp.zeros(2)
    require_zero(large_nilpotent**2, "bounded nilpotent correction sequence")
    require_zero(
        large_nilpotent.T * large_nilpotent - sp.diag(0, 9),
        "noncommuting counterexample norm square",
    )

    # Symbolic all-power scalar countermodel.  The exact recurrence and
    # initial magnitude show exponential unboundedness analytically.
    scalar_operator = sp.Integer(3)
    symbolic_power = sp.Symbol("j", integer=True, positive=True)
    scalar_correction = -(scalar_operator**symbolic_power)
    assert scalar_operator > 2
    require_zero(
        2 * sp.Integer(0) ** symbolic_power
        - scalar_operator**symbolic_power
        - scalar_correction,
        "unbounded scalar correction definition",
    )
    require_zero(
        scalar_correction * scalar_operator - scalar_operator * scalar_correction,
        "unbounded scalar correction commutation",
    )
    require_zero(
        scalar_correction.subs(symbolic_power, symbolic_power + 1)
        - 3 * scalar_correction,
        "unbounded scalar correction growth recurrence",
    )
    require_zero(
        scalar_correction.subs(symbolic_power, 1) + 3,
        "unbounded scalar correction initial value",
    )

    return {
        "arithmetic": "exact SymPy rationals; no numerical tolerances",
        "maximum_power": maximum_power,
        "operator_norm": str(kappa),
        "contraction_slack_first_minor": str(contraction_slack[0, 0]),
        "contraction_slack_determinant": str(contraction_slack.det()),
        "all_power_correction_formula": "E_j = (1/3)^j I (binomial; N^2=0)",
        "defect_square": str(defect_square),
        "defect_expansion_margin": str(defect_upper_bound - defect_square),
        "power_records": records,
        "removed_hypothesis_counterexamples": {
            "no_commutation": {
                "T": [[0, 3], [0, 0]],
                "Q": "0",
                "corrections": "E_1=-T*, E_j=0 for j>=2",
                "commutator": str(commutator),
            },
            "no_uniform_bound": {
                "T": "3I",
                "Q": "0",
                "corrections": "E_j=-3^j I, commuting but unbounded",
                "symbolically_verified": "E_1=-3, E_(j+1)=3 E_j, [E_j,T]=0",
                "unboundedness_argument": "|E_1|=3 and growth ratio 3>1 give |E_j|=3^j -> infinity",
            },
        },
        "all_checks_passed": True,
        "scope": "algebra diagnostic; general proof is in the audit note",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-power", type=int, default=8)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    if arguments.maximum_power < 1:
        parser.error("--maximum-power must be positive")
    result = audit(arguments.maximum_power)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.output is not None:
        arguments.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
