#!/usr/bin/env python3
"""Numerically probe the real ``p=7`` sixth-face SOS constant.

The exact base polynomial ``P_6`` and cubic row ``C_6`` are built with
SymPy.  A degree-three Gram SDP maximizes ``alpha`` subject to

    P_6 - alpha C_6**2 = m_3.T Q m_3,   Q >= 0.

This is a discovery calculation only.  Floating SDP output, especially
an ``optimal_inaccurate`` status, is not an exact sum-of-squares proof.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import itertools
import json
from pathlib import Path

import cvxpy as cp
import numpy as np
import sympy as sp

from crabb_full_disk_base_jet_exact import endpoint_delta_series


@dataclass(frozen=True)
class RealSosProbeRecord:
    """One numerical degree-three Gram optimization."""

    dimension: int
    length: int
    solver: str
    solver_status: str
    gram_size: int
    optimized_constant: float
    conjectured_constant: float
    required_schur_constant: float
    minimum_gram_eigenvalue: float


def degree_exponents(
    variable_count: int,
    degree: int,
) -> list[tuple[int, ...]]:
    """Return exponent tuples of one homogeneous degree."""

    return [
        exponent
        for exponent in itertools.product(
            range(degree + 1),
            repeat=variable_count,
        )
        if sum(exponent) == degree
    ]


def exact_polynomials() -> tuple[
    tuple[sp.Symbol, ...],
    sp.Poly,
    sp.Poly,
]:
    """Return the exact real ``P_6`` and ``C_6`` polynomials."""

    variables = sp.symbols("a1:6", real=True)
    delta = endpoint_delta_series(
        (sp.Integer(0), *variables),
        order=6,
    )[6]
    base = sp.Poly(
        sp.expand(-sp.Rational(225, 32) * delta),
        *variables,
    )
    left, second, middle, fourth, right = variables
    cubic = sp.Poly(
        6 * left * second * fourth
        - 5 * left * middle * right
        + 2 * second * middle * fourth
        - 2 * middle * fourth**2
        + 5 * middle * right**2
        - 6 * fourth**2 * right,
        *variables,
    )
    return variables, base, cubic


def solve_gram_sdp(solver: str) -> RealSosProbeRecord:
    """Build and solve the homogeneous Gram SDP."""

    variables, base, cubic = exact_polynomials()
    monomials = degree_exponents(len(variables), 3)
    degree_six = degree_exponents(len(variables), 6)
    base_coefficients = base.as_dict()
    cubic_square_coefficients = sp.Poly(
        cubic.as_expr() ** 2,
        *variables,
    ).as_dict()

    gram = cp.Variable(
        (len(monomials), len(monomials)),
        symmetric=True,
    )
    alpha = cp.Variable()
    constraints = [gram >> 0]
    for target in degree_six:
        entries = [
            gram[left, right]
            for left, left_exponent in enumerate(monomials)
            for right, right_exponent in enumerate(monomials)
            if tuple(
                left_exponent[index] + right_exponent[index]
                for index in range(len(variables))
            )
            == target
        ]
        coefficient = (
            float(base_coefficients.get(target, 0))
            - alpha
            * float(cubic_square_coefficients.get(target, 0))
        )
        constraints.append(cp.sum(cp.hstack(entries)) == coefficient)

    problem = cp.Problem(cp.Maximize(alpha), constraints)
    problem.solve(solver=solver, verbose=False)
    if gram.value is None or alpha.value is None:
        raise RuntimeError(f"{solver} did not return a Gram matrix")
    minimum_eigenvalue = float(
        np.linalg.eigvalsh(np.asarray(gram.value))[
            0
        ]
    )
    return RealSosProbeRecord(
        dimension=7,
        length=6,
        solver=solver,
        solver_status=str(problem.status),
        gram_size=len(monomials),
        optimized_constant=float(alpha.value),
        conjectured_constant=9.0,
        required_schur_constant=float(sp.Rational(1089, 290)),
        minimum_gram_eigenvalue=minimum_eigenvalue,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", default="CLARABEL")
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run the numerical real-SOS probe."""

    args = parse_args()
    if args.solver not in cp.installed_solvers():
        raise ValueError(f"solver {args.solver!r} is not installed")
    record = solve_gram_sdp(args.solver)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    line = json.dumps(asdict(record), sort_keys=True)
    args.output.write_text(f"{line}\n", encoding="utf-8")
    print(line, flush=True)


if __name__ == "__main__":
    main()
