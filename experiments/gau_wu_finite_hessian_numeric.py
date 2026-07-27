#!/usr/bin/env python3
"""Falsify the optimized Hessian sign at random finite Gau--Wu models."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from gau_wu_finite_model import (
    blaschke_at_matrix,
    disk_model_residuals,
    expected_tangent_rank,
    extremal_zeros,
    gau_wu_model,
    normal_basis,
    pack_zeros,
    random_interior_zeros,
    real_matrix,
    unpack_zeros,
)
from theodorsen import GeneralPullback, theodorsen_map


@dataclass(frozen=True)
class FiniteGauWuHessianRecord:
    """One random finite-model quotient-Hessian falsification record."""

    dimension: int
    sample: int
    seed: int
    interior_zeros: tuple[str, ...]
    tangent_rank: int
    expected_tangent_rank: int
    normal_dimension: int
    tested_dimension: int
    full_normal_hessian: bool
    support_error: float
    functional_error: float
    minimum_support_gap: float
    base_optimized_value: float
    maximum_map_residual: float
    finite_difference_step: float
    eigenvalues: tuple[float, ...]
    maximum_eigenvalue: float
    sign_classification: str
    refined_witness_coefficients: tuple[float, ...]
    all_checks_passed: bool


def normalized_operator(
    matrix: np.ndarray,
    resolution: int,
) -> tuple[np.ndarray, float]:
    """Compute the Riemann-pulled operator by Theodorsen and Cauchy."""

    boundary, derivative, residual = theodorsen_map(
        matrix,
        N=resolution,
        inflate=0.0,
        iters=400,
        tol=3e-14,
    )
    pullback = GeneralPullback(boundary, derivative)
    resolvents = pullback.resolvent_stack(matrix)
    operator = pullback.calc(pullback.w, matrix, resolvents)
    return operator, float(residual)


def optimize_zeros(
    operator: np.ndarray,
    initial_parameters: np.ndarray,
) -> tuple[float, np.ndarray]:
    """Optimize all zeros of the finite Blaschke extremal locally."""

    def objective(parameters: np.ndarray) -> float:
        zeros = unpack_zeros(parameters)
        return -float(np.linalg.norm(blaschke_at_matrix(operator, zeros), 2))

    result = minimize(
        objective,
        initial_parameters,
        method="BFGS",
        options={"maxiter": 500, "gtol": 3e-9},
    )
    if not result.success or np.linalg.norm(result.jac) > 5e-5:
        result = minimize(
            objective,
            result.x,
            method="Nelder-Mead",
            options={
                "maxiter": 4000,
                "xatol": 3e-10,
                "fatol": 3e-13,
            },
        )
    return -float(result.fun), np.asarray(result.x)


def optimal_value(
    matrix: np.ndarray,
    initial_parameters: np.ndarray,
    resolution: int,
) -> tuple[float, np.ndarray, float]:
    """Normalize the numerical range and optimize the Blaschke zeros."""

    operator, residual = normalized_operator(matrix, resolution)
    value, parameters = optimize_zeros(operator, initial_parameters)
    return value, parameters, residual


def directional_coefficient(
    matrix: np.ndarray,
    direction: np.ndarray,
    step: float,
    base_value: float,
    initial_parameters: np.ndarray,
    resolution: int,
) -> tuple[float, float]:
    """Return the symmetric second coefficient in one physical direction."""

    plus, _, plus_residual = optimal_value(
        matrix + step * direction,
        initial_parameters,
        resolution,
    )
    minus, _, minus_residual = optimal_value(
        matrix - step * direction,
        initial_parameters,
        resolution,
    )
    coefficient = (plus + minus - 2 * base_value) / (2 * step**2)
    return coefficient, max(plus_residual, minus_residual)


def probe_basis(
    normal: np.ndarray,
    maximum_dimension: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, bool]:
    """Return either the full normal basis or a random orthonormal subspace."""

    normal_dimension = normal.shape[1]
    if normal_dimension <= maximum_dimension:
        return normal, True
    random_coordinates = rng.normal(
        size=(normal_dimension, maximum_dimension)
    )
    coordinate_basis, _ = np.linalg.qr(random_coordinates)
    return normal @ coordinate_basis[:, :maximum_dimension], False


def finite_difference_hessian(
    matrix: np.ndarray,
    basis: np.ndarray,
    dimension: int,
    step: float,
    base_value: float,
    initial_parameters: np.ndarray,
    resolution: int,
) -> tuple[np.ndarray, float]:
    """Reconstruct the Hessian on the supplied real physical subspace."""

    count = basis.shape[1]
    hessian = np.zeros((count, count))
    maximum_residual = 0.0
    directions = [
        real_matrix(basis[:, index], dimension) for index in range(count)
    ]
    for row, direction in enumerate(directions):
        coefficient, residual = directional_coefficient(
            matrix,
            direction,
            step,
            base_value,
            initial_parameters,
            resolution,
        )
        hessian[row, row] = coefficient
        maximum_residual = max(maximum_residual, residual)

    for row in range(count):
        for column in range(row + 1, count):
            coefficient, residual = directional_coefficient(
                matrix,
                directions[row] + directions[column],
                step,
                base_value,
                initial_parameters,
                resolution,
            )
            entry = (
                coefficient - hessian[row, row] - hessian[column, column]
            ) / 2
            hessian[row, column] = entry
            hessian[column, row] = entry
            maximum_residual = max(maximum_residual, residual)
    return hessian, maximum_residual


def refined_witness(
    matrix: np.ndarray,
    basis: np.ndarray,
    eigenvector: np.ndarray,
    step: float,
    base_value: float,
    initial_parameters: np.ndarray,
    resolution: int,
) -> tuple[float, ...]:
    """Recheck a possible positive mode at smaller steps and higher resolution."""

    direction_vector = basis @ eigenvector
    direction = real_matrix(direction_vector, matrix.shape[0])
    coefficients: list[float] = []
    for divisor in (1, 2):
        coefficient, _ = directional_coefficient(
            matrix,
            direction,
            step / divisor,
            base_value,
            initial_parameters,
            2 * resolution,
        )
        coefficients.append(coefficient)
    return tuple(coefficients)


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    maximum_probe_dimension: int,
    step: float,
    resolution: int,
) -> FiniteGauWuHessianRecord:
    """Generate and audit one random finite Gau--Wu equality model."""

    rng = np.random.default_rng(seed)
    interior = random_interior_zeros(dimension, rng)
    matrix = gau_wu_model(interior)
    zeros = extremal_zeros(interior)
    support_error, functional_error, minimum_gap = disk_model_residuals(
        matrix,
        zeros,
    )

    normal, tangent_rank, _ = normal_basis(interior)
    expected_rank = expected_tangent_rank(dimension)
    basis, is_full = probe_basis(normal, maximum_probe_dimension, rng)

    initial_parameters = pack_zeros(zeros)
    base_value, optimized_parameters, base_residual = optimal_value(
        matrix,
        initial_parameters,
        resolution,
    )
    hessian, maximum_residual = finite_difference_hessian(
        matrix,
        basis,
        dimension,
        step,
        base_value,
        optimized_parameters,
        resolution,
    )
    maximum_residual = max(maximum_residual, base_residual)
    eigenvalues, eigenvectors = np.linalg.eigh(hessian)
    maximum_eigenvalue = float(eigenvalues[-1])

    if maximum_eigenvalue > 2e-3:
        refinement = refined_witness(
            matrix,
            basis,
            eigenvectors[:, -1],
            step,
            base_value,
            optimized_parameters,
            resolution,
        )
        if min(refinement) > 5e-4:
            classification = "positive witness"
        elif max(refinement) < -1e-5:
            classification = "rejected reconstructed positive mode"
        else:
            classification = "numerically unresolved"
    elif maximum_eigenvalue < -5e-4:
        classification = "strictly negative on tested subspace"
        refinement = ()
    else:
        classification = "numerically unresolved"
        refinement = refined_witness(
            matrix,
            basis,
            eigenvectors[:, -1],
            step,
            base_value,
            optimized_parameters,
            resolution,
        )

    checks = (
        tangent_rank == expected_rank
        and support_error < 2e-12
        and functional_error < 2e-12
        and minimum_gap > 1e-3
        and abs(base_value - 2) < 2e-8
        and maximum_residual < 2e-8
        and (
            classification != "positive witness"
            or min(refinement) > 5e-4
        )
    )
    if not checks:
        raise RuntimeError(
            "finite Gau--Wu Hessian audit failed diagnostics: "
            f"n={dimension}, rank={tangent_rank}/{expected_rank}, "
            f"support={support_error}, functional={functional_error}, "
            f"gap={minimum_gap}, base={base_value}, "
            f"map={maximum_residual}, class={classification}, "
            f"refinement={refinement}"
        )

    return FiniteGauWuHessianRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        interior_zeros=tuple(
            f"{zero.real:+.12f}{zero.imag:+.12f}j" for zero in interior
        ),
        tangent_rank=tangent_rank,
        expected_tangent_rank=expected_rank,
        normal_dimension=normal.shape[1],
        tested_dimension=basis.shape[1],
        full_normal_hessian=is_full,
        support_error=support_error,
        functional_error=functional_error,
        minimum_support_gap=minimum_gap,
        base_optimized_value=base_value,
        maximum_map_residual=maximum_residual,
        finite_difference_step=step,
        eigenvalues=tuple(map(float, eigenvalues)),
        maximum_eigenvalue=maximum_eigenvalue,
        sign_classification=classification,
        refined_witness_coefficients=refinement,
        all_checks_passed=True,
    )


def write_records(
    records: list[FiniteGauWuHessianRecord],
    output: Path,
) -> str:
    """Write deterministic JSON Lines atomically and return its SHA-256."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)
    return hashlib.sha256(output.read_bytes()).hexdigest()


def parse_dimensions(value: str) -> tuple[int, ...]:
    """Parse a comma-separated dimension list."""

    dimensions = tuple(int(item) for item in value.split(","))
    if not dimensions or any(dimension < 3 for dimension in dimensions):
        raise argparse.ArgumentTypeError("dimensions must be integers at least 3")
    return dimensions


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dimensions", type=parse_dimensions, default=(4,))
    parser.add_argument("--samples", type=int, default=1)
    parser.add_argument("--seed", type=int, default=70224)
    parser.add_argument("--maximum-probe-dimension", type=int, default=12)
    parser.add_argument("--step", type=float, default=0.003)
    parser.add_argument("--resolution", type=int, default=192)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/gau_wu_finite_hessian_numeric_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the requested deterministic finite-model falsification grid."""

    arguments = parse_args()
    records: list[FiniteGauWuHessianRecord] = []
    for dimension in arguments.dimensions:
        for sample in range(arguments.samples):
            seed = arguments.seed + 1009 * dimension + sample
            record = audit_model(
                dimension,
                sample,
                seed,
                arguments.maximum_probe_dimension,
                arguments.step,
                arguments.resolution,
            )
            records.append(record)
            print(
                json.dumps(
                    {
                        "dimension": dimension,
                        "sample": sample,
                        "maximum_eigenvalue": record.maximum_eigenvalue,
                        "classification": record.sign_classification,
                    },
                    sort_keys=True,
                ),
                flush=True,
            )
    digest = write_records(records, arguments.output)
    print(
        json.dumps(
            {
                "all_checks_passed": all(
                    record.all_checks_passed for record in records
                ),
                "record_count": len(records),
                "sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
