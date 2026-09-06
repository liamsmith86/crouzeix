#!/usr/bin/env python3
"""Non-load-bearing exact and numerical audits of the all-powers proof.

Reference: Lorist--Schwenninger, arXiv:2608.03841v2, Lemma 1/Theorem 3.
The contour tests use positively oriented convex enclosing polygons. Neither
quadrature nor a finite collection of powers certifies the general theorem.
The output distinguishes exact identities from floating-point diagnostics.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
import sympy as sp

from crouzeix import nr_support, poly_A, poly_z


TOLERANCE = 2e-9
# Flat numerical ranges can put a pole close to a long polygon edge. Both
# orders resolve the segment fixture; 48 nodes failed the unchanged mass gate.
QUADRATURE_ORDERS = (192, 384)
DEFAULT_OUTPUT = Path(__file__).with_suffix(".json")


def require(condition: bool, message: str) -> None:
    """Keep audit gates active even when Python assertions are disabled."""

    if not condition:
        raise RuntimeError(message)


def exact_zero(expression: sp.Expr | sp.MatrixBase, label: str) -> None:
    """Reject a symbolic residual unless simplification proves it zero."""

    entries = (
        list(expression) if isinstance(expression, sp.MatrixBase) else [expression]
    )
    require(all(sp.simplify(entry) == 0 for entry in entries), label)


def symbolic_identities() -> dict[str, object]:
    """Check the recurrence before imposing any dilation estimates."""

    kappa, sigma = sp.symbols("kappa sigma", real=True)
    rotation = sp.Matrix([[3, -4], [4, 3]]) / 5
    operator = rotation * sp.diag(kappa, sigma)
    vector = sp.Matrix([1, 0])
    exact_zero(operator.T * operator * vector - kappa**2 * vector, "singular vector")
    # Each independently parameterized defect commutes with T. S_n* is then
    # defined by the defect identity; no norm or positivity claim is needed.
    parameters = sp.symbols("a1:5 b1:5", real=True)
    defects = [
        parameters[n] * sp.eye(2) + parameters[n + 4] * operator for n in range(4)
    ]
    for power in range(1, 4):
        current, following = defects[power - 1 : power + 1]
        m_current = (vector.T * current * operator**power * vector)[0]
        m_next = (vector.T * following * operator ** (power + 1) * vector)[0]
        s_current = current + operator.T**power
        s_next = following + operator.T ** (power + 1)
        y = (s_next * operator - kappa * s_current) * vector
        v = operator.T**power * vector
        exact_zero(
            kappa * m_current
            - m_next
            - kappa * (kappa - 1) * (v.T * v)[0]
            + (y.T * v)[0],
            f"symbolic recurrence n={power}",
        )

    count = 7
    moments = sp.symbols(f"m1:{count + 2}")
    telescoped = sum(
        kappa ** (1 - n) * (kappa * moments[n - 1] - moments[n])
        for n in range(1, count + 1)
    )
    exact_zero(
        telescoped - kappa * moments[0] + kappa ** (1 - count) * moments[count],
        "finite telescoping including terminal term",
    )
    a, v2, z2, cross = sp.symbols("a v2 z2 cross", real=True, nonzero=True)
    exact_zero(
        a * v2 - 2 * cross - (a * (v2 - 2 * cross / a + z2 / a**2) - z2 / a),
        "square completion",
    )
    d2, m1 = sp.symbols("d2 m1", real=True)
    # Substituting the lower bound for kappa*m1 produces exactly the final
    # contradiction expression, including the crucial (kappa-1)^2 factor.
    upper_residual = d2 - (2 * kappa**2 - kappa * m1 - kappa**3)
    exact_zero(
        upper_residual.subs(m1, -d2 / (kappa * (kappa - 1) ** 2))
        - (d2 * (1 - 1 / (kappa - 1) ** 2) - kappa**2 * (2 - kappa)),
        "final scalar elimination",
    )
    return {
        "passed": True,
        "symbolic_recurrence_powers": [1, 2, 3],
        "finite_telescoping_length": count,
        "terminal_term_retained": True,
        "square_completion": "exact",
        "final_scalar_elimination": "exact",
        "scope": "Algebraic identities, not an automated proof of the limit or theorem.",
    }


def exact_hypothesis_checks() -> dict[str, object]:
    """Give all-power formulas, not finite-power guesses, for countermodels."""

    shift = sp.Matrix([[0, 1], [0, 0]])
    identity = sp.eye(2)
    zero = sp.zeros(2)
    sharp = 2 * shift
    exact_zero(sharp.T * sharp - sp.diag(0, 4), "sharp norm")
    a, b = sp.symbols("a b", real=True)
    exact_zero((a**2 + b**2) ** 2 - 4 * a**2 * b**2 - (a**2 - b**2) ** 2, "radius")
    # Unit vectors (1,z)/sqrt(2), |z|=1, attain every boundary point.
    z = sp.symbols("z")
    boundary_value = (sp.Matrix([[1, 1 / z]]) * sharp * sp.Matrix([1, z]))[0] / 2
    exact_zero(boundary_value - z, "sharp numerical-range boundary")
    exact_zero(sharp**2, "sharp nilpotence")
    exact_zero(2 * shift.T - sharp.T, "sharp exact dilation")

    models = [
        ("commutation", 3 * shift, zero, identity, "E_1=-3 J*, E_n=0 for n>=2"),
        ("boundedness", 3 * identity, zero, identity, "E_n=-3^n I for every n>=1"),
        ("isometry", 8 * shift, shift, 2 * identity, "E_n=0 for every n>=1"),
        (
            "contraction",
            3 * shift,
            sp.Rational(3, 2) * shift,
            identity,
            "E_n=0 for every n>=1",
        ),
    ]
    records = []
    for dropped, operator, contraction, embedding, formula in models:
        defect = 2 * embedding.T * contraction.T * embedding - operator.T
        isometry = embedding.T * embedding == identity
        contractive = (
            identity - contraction.T * contraction
        ).is_positive_semidefinite is True
        commutes = defect * operator == operator * defect
        bounded = dropped != "boundedness"
        conditions = {
            "commutation": commutes,
            "boundedness": bounded,
            "isometry": isometry,
            "contraction": contractive,
        }
        require(
            all(value == (name != dropped) for name, value in conditions.items()),
            dropped,
        )
        if dropped == "boundedness":
            n = sp.symbols("n", integer=True, positive=True)
            exact_zero(operator**n - 3**n * identity, "unbounded exact powers")
            exact_zero(contraction, "zero dilation powers")
        else:
            # T^2=Q^2=0 proves the displayed tail for every n>=2.
            exact_zero(operator**2, f"{dropped}: nilpotent T")
            exact_zero(contraction**2, f"{dropped}: nilpotent Q")
            if dropped != "commutation":
                exact_zero(defect, f"{dropped}: zero first defect")
        norm_squared = max((operator.T * operator).eigenvals())
        require(norm_squared > 4, f"{dropped}: failure of conclusion")
        records.append(
            {
                "dropped_hypothesis": dropped,
                "hypotheses_satisfied": conditions,
                "operator": [
                    [str(entry) for entry in row] for row in operator.tolist()
                ],
                "Q": [[str(entry) for entry in row] for row in contraction.tolist()],
                "V": [[str(entry) for entry in row] for row in embedding.tolist()],
                "operator_norm_squared": str(norm_squared),
                "all_powers_defect_formula": formula,
                "passed": True,
            }
        )
    return {
        "sharpness": {
            "operator": "2 J, J=[[0,1],[0,0]]",
            "norm": 2,
            "numerical_range": "closed unit disk",
            "polynomial": "z",
            "Q": "J",
            "V": "I",
            "defects": "zero for every n>=1",
            "passed": True,
        },
        "one_hypothesis_countermodels": records,
    }


def enclosing_polygon(matrix: np.ndarray, sides: int = 12) -> np.ndarray:
    """Intersect adjacent positively separated numerical-range support lines."""

    angles = np.arange(sides) * (2 * np.pi / sides)
    supports, _ = nr_support(matrix, angles)
    margin = 0.08 * max(1.0, float(np.linalg.norm(matrix, 2)))
    supports += margin
    normals = np.column_stack((np.cos(angles), np.sin(angles)))
    vertices = np.empty(sides, dtype=complex)
    for index in range(sides):
        indices = [index, (index + 1) % sides]
        point = np.linalg.solve(normals[indices], supports[indices])
        vertices[index] = point[0] + 1j * point[1]
    coordinates = np.column_stack((vertices.real, vertices.imag))
    require(
        float(np.max(normals @ coordinates.T - supports[:, None])) < TOLERANCE,
        "polygon support feasibility",
    )
    area_twice = np.imag(np.sum(vertices.conj() * np.roll(vertices, -1)))
    require(area_twice > 0, "polygon positive orientation")
    return vertices


def contour_layer(
    matrix: np.ndarray, vertices: np.ndarray, order: int
) -> dict[str, np.ndarray | float]:
    """Construct Cauchy weights C_j and positive Gram weights (C_j+C_j*)/2."""

    nodes, weights = leggauss(order)
    starts = vertices[:, None]
    edges = (np.roll(vertices, -1) - vertices)[:, None]
    points = (starts + edges * (nodes + 1) / 2).ravel()
    differentials = (edges * weights / 2).ravel()
    identity = np.eye(len(matrix), dtype=complex)
    resolvents = np.linalg.inv(points[:, None, None] * identity - matrix)
    cauchy = differentials[:, None, None] * resolvents / (2j * np.pi)
    gram = (cauchy + cauchy.conj().transpose(0, 2, 1)) / 2
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    require(float(np.min(eigenvalues)) > 0, "strict contour positivity")
    square_roots = (
        eigenvectors * np.sqrt(eigenvalues)[:, None, :]
    ) @ eigenvectors.conj().transpose(0, 2, 1)
    mass_error = float(np.linalg.norm(np.sum(gram, axis=0) - identity, 2))
    cauchy_mass_error = float(np.linalg.norm(np.sum(cauchy, axis=0) - identity, 2))
    require(
        max(mass_error, cauchy_mass_error) < TOLERANCE,
        f"contour mass at order {order}: positive={mass_error}, Cauchy={cauchy_mass_error}",
    )
    return {
        "points": points,
        "cauchy": cauchy,
        "gram": gram,
        "square_roots": square_roots,
        "mass_error": mass_error,
        "cauchy_mass_error": cauchy_mass_error,
        "minimum_gram_eigenvalue": float(np.min(eigenvalues)),
    }


def polynomial_families(dimension: int) -> dict[str, np.ndarray]:
    """Include zero/constant cases and a polynomial above the matrix dimension."""

    high_degree = np.zeros(dimension + 4, dtype=complex)
    high_degree[[0, 1, -1]] = [0.02 + 0.01j, 0.9, 0.04 - 0.03j]
    return {
        "zero": np.array([0.0]),
        "constant": np.array([1.0]),
        "linear": np.array([0.0, 1.0]),
        "high_degree": high_degree,
    }


def contour_case(
    label: str,
    matrix: np.ndarray,
    vertices: np.ndarray,
    family: str,
    coefficients: np.ndarray,
    layer: dict[str, np.ndarray | float],
    order: int,
    maximum_power: int,
) -> dict[str, object]:
    """Check the application, recurrence, and finite telescoping independently."""

    points = layer["points"]
    cauchy, gram, roots = (layer[name] for name in ("cauchy", "gram", "square_roots"))
    dimension = len(matrix)
    center = np.trace(matrix) / dimension
    radius = float(np.max(abs(vertices - center)))
    # The entire polygon lies in this disk. Triangle inequality, not a grid
    # maximum, ensures |f|<=1 everywhere after coefficient normalization.
    coefficient_bound = float(np.sum(abs(coefficients)))
    normalized = coefficients / coefficient_bound if coefficient_bound else coefficients
    values = poly_z(normalized, (points - center) / radius)
    operator = poly_A(normalized, (matrix - center * np.eye(dimension)) / radius)
    require(float(np.max(abs(values))) <= 1 + TOLERANCE, "polynomial contraction")
    _, singular_values, vh = np.linalg.svd(operator)
    kappa = float(singular_values[0])
    vector = vh[0].conj()
    vx = roots @ vector
    vtx = roots @ (operator @ vector)
    u = values.conj()[:, None] * vtx - kappa * vx
    d2 = float(np.vdot(u, u).real)

    powers = [np.eye(dimension, dtype=complex)]
    moments = [0.0]
    errors: dict[str, list[float]] = {
        name: []
        for name in (
            "cauchy_power",
            "companion_identity",
            "commutator",
            "recurrence",
            "completed_square",
        )
    }
    defects = []
    for power in range(1, maximum_power + 2):
        powers.append(powers[-1] @ operator)
        moment = np.einsum("j,jab->ab", values.conj() ** power, gram)
        defect = 2 * moment - powers[-1].conj().T
        companion = np.einsum("j,jab->ab", values.conj() ** power, cauchy)
        cauchy_power = np.einsum("j,jab->ab", values**power, cauchy)
        errors["cauchy_power"].append(
            float(np.linalg.norm(cauchy_power - powers[-1], 2))
        )
        errors["companion_identity"].append(
            float(np.linalg.norm(defect - companion, 2))
        )
        errors["commutator"].append(
            float(np.linalg.norm(defect @ operator - operator @ defect, 2))
        )
        defects.append(defect)
        moments.append(float(np.vdot(vector, defect @ powers[-1] @ vector).real))

    recurrence_slacks = []
    for power in range(1, maximum_power + 1):
        v = powers[power].conj().T @ vector
        z = np.einsum("jab,jb->a", roots, values.conj()[:, None] ** power * u)
        lhs = kappa * moments[power] - moments[power + 1]
        a = kappa * (kappa - 1)
        rhs = a * float(np.vdot(v, v).real) - 2 * float(np.vdot(z, v).real)
        errors["recurrence"].append(abs(lhs - rhs))
        if kappa > 1 + 1e-8:
            completed = (
                a * float(np.vdot(v - z / a, v - z / a).real)
                - float(np.vdot(z, z).real) / a
            )
            errors["completed_square"].append(abs(lhs - completed))
            recurrence_slacks.append(lhs + d2 / a)

    telescope_error = None
    telescope_tail = None
    if kappa > 1 + 1e-8:
        finite_sum = sum(
            kappa ** (1 - power) * (kappa * moments[power] - moments[power + 1])
            for power in range(1, maximum_power + 1)
        )
        telescope_tail = kappa ** (1 - maximum_power) * moments[maximum_power + 1]
        telescope_error = abs(kappa * moments[1] - finite_sum - telescope_tail)
        require(telescope_error < TOLERANCE, "finite numerical telescope")
        require(min(recurrence_slacks) > -TOLERANCE, "recurrence lower bound")
    norm_bound_slack = 2 * kappa**2 - kappa * moments[1] - kappa**3 - d2
    maxima = {name: max(values, default=0.0) for name, values in errors.items()}
    require(max(maxima.values()) < TOLERANCE, f"{label}/{family}: {maxima}")
    require(norm_bound_slack > -TOLERANCE, "one-step error-vector norm bound")
    companion_bound = float(np.sum(np.linalg.norm(cauchy, ord=2, axis=(1, 2))))
    max_defect = max(float(np.linalg.norm(defect, 2)) for defect in defects)
    require(max_defect <= companion_bound + TOLERANCE, "discrete companion bound")
    return {
        "matrix": label,
        "dimension": dimension,
        "polynomial_family": family,
        "polynomial_degree": len(coefficients) - 1,
        "quadrature_order_per_edge": order,
        "powers_tested": maximum_power + 1,
        "operator_norm": kappa,
        "normalization": "l1 coefficient bound on an enclosing disk",
        "mass_error": layer["mass_error"],
        "cauchy_mass_error": layer["cauchy_mass_error"],
        "minimum_gram_eigenvalue": layer["minimum_gram_eigenvalue"],
        "maximum_residuals": maxima,
        "one_step_norm_bound_slack": norm_bound_slack,
        "minimum_recurrence_lower_bound_slack": min(recurrence_slacks, default=None),
        "finite_telescope_residual": telescope_error,
        "finite_telescope_terminal_term": telescope_tail,
        "discrete_companion_uniform_bound": companion_bound,
        "maximum_tested_defect_norm": max_defect,
        "passed": True,
    }


def matrix_families(seed: int) -> dict[str, np.ndarray]:
    """Fixed structural edge cases plus seeded genuinely complex dense cases."""

    generator = np.random.default_rng(seed)
    jordan = np.diag(np.full(3, 1.2), 1).astype(complex) + (0.2 + 0.1j) * np.eye(4)
    output = {
        "scalar1": np.array([[0.3 + 0.4j]]),
        "zero2": np.zeros((2, 2), dtype=complex),
        "sharp_jordan2": np.array([[0, 2], [0, 0]], dtype=complex),
        "hermitian_segment3": np.diag([-1.0, 0.0, 2.0]).astype(complex),
        "normal_polygon3": np.diag([0.0, 1.0, 1.0j]),
        "large_nilpotent3": np.diag([1000.0, 1000.0], 1).astype(complex),
        "defective_jordan4": jordan,
        "singular_triangular3": np.array(
            [[0, 1.1, 0.3j], [0, 0.2, 0.8], [0, 0, -0.3]], dtype=complex
        ),
    }
    for dimension in (3, 5):
        output[f"dense{dimension}"] = (
            generator.standard_normal((dimension, dimension))
            + 1j * generator.standard_normal((dimension, dimension))
        ) / np.sqrt(dimension)
    return output


def run_audit(seed: int, maximum_power: int) -> dict[str, object]:
    """Execute every exact gate and the two-resolution contour suite."""

    exact = {"identities": symbolic_identities(), **exact_hypothesis_checks()}
    records = []
    for label, matrix in matrix_families(seed).items():
        vertices = enclosing_polygon(matrix)
        for order in QUADRATURE_ORDERS:
            try:
                layer = contour_layer(matrix, vertices, order)
            except RuntimeError as error:
                raise RuntimeError(f"{label}: {error}") from error
            for family, coefficients in polynomial_families(len(matrix)).items():
                records.append(
                    contour_case(
                        label,
                        matrix,
                        vertices,
                        family,
                        coefficients,
                        layer,
                        order,
                        maximum_power,
                    )
                )
    require(
        any(record["operator_norm"] > 1.2 for record in records),
        "nontrivial norm regime coverage",
    )
    return {
        "schema_version": 1,
        "reference": "https://arxiv.org/html/2608.03841v2",
        "scope": "Non-load-bearing regression tests; finite floating-point quadrature is not a proof.",
        "seed": seed,
        "tolerance": TOLERANCE,
        "exact": exact,
        "numerical_records": records,
        "passed": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--seed", type=int, default=20260905)
    parser.add_argument("--maximum-power", type=int, default=16)
    args = parser.parse_args()
    if args.maximum_power < 1:
        parser.error("--maximum-power must be positive")
    result = run_audit(args.seed, args.maximum_power)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print(
        f"Passed exact gates and {len(result['numerical_records'])} numerical cases; {args.output}"
    )


if __name__ == "__main__":
    main()
