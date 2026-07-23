#!/usr/bin/env python3
"""Audit the exact central Crabb metric lift from L130.

The computation deliberately avoids optimizing a defect vector in the full
``(2k+1)``-dimensional problem.  It optimizes only the descended size-three
problem, factors the outer critical polynomial of the degree-``k`` Blaschke
product, and uses L130's explicit multiplier to construct the full defect.
It then compares that vector with the independently transferred dual kernel
and constructs the resulting full rank-one Stein metric.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import null_space, solve_discrete_lyapunov

from blaschke_stein_composition import (
    dual_transfer,
    model_functions,
    real_blaschke_critical_factor,
)
from crabb_central_metric_descent import central_coordinate_gramian
from crabb_palindromic_elliptic_face import (
    disk_matrix,
    ellipse_pullback,
    palindromic_direction,
    rank_one_envelope,
)
from crabb_touching_gradient import chebyshev_blaschke_zeros


DEFAULT_DEGREES = (2, 3, 4, 5)
DEFAULT_ELLIPSE_PARAMETERS = (0.15, 0.3, 0.5)
DEFAULT_AMPLITUDE = 0.08


@dataclass(frozen=True)
class CentralDualLiftRecord:
    degree: int
    dimension: int
    ellipse_parameter: float
    descended_parameter: float
    amplitude: float
    outer_image_residual: float
    critical_factorization_residual: float
    model_multiplier_residual: float
    scalar_compression_residual: float
    dual_smallest_eigenvalue: float
    dual_second_eigenvalue: float
    explicit_dual_kernel_residual: float
    explicit_kernel_angle_residual: float
    weighted_reconstruction_residual: float
    kernel_alignment_residual: float
    forcing_cross_residual: float
    random_dual_kernel_alignment_residual: float
    random_dual_forcing_cross_residual: float
    metric_cross_residual: float
    outer_metric_residual: float
    constructed_condition: float
    outer_condition: float
    condition_gap: float


def central_operator(
    degree: int,
    ellipse_parameter: float,
    amplitude: float,
) -> np.ndarray:
    """Return the physical central operator of size ``2*degree+1``."""

    dimension = 2 * degree + 1
    coefficients = palindromic_direction(2 * degree, degree)
    disk_operator = disk_matrix(dimension, amplitude, coefficients)
    return np.asarray(
        ellipse_pullback(disk_operator, ellipse_parameter),
        dtype=complex,
    )


def optimized_rank_one_metric(
    operator: np.ndarray,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Return the optimized condition, defect, and rank-one Stein metric."""

    condition, tail = rank_one_envelope(np.real_if_close(operator).real)
    defect = np.concatenate(([1.0], tail)).astype(complex)
    metric = solve_discrete_lyapunov(
        operator.conj().T,
        np.outer(defect, defect.conj()),
    )
    return condition, defect, metric


def outer_isometry(
    degree: int,
    amplitude: float,
) -> np.ndarray:
    """Map the exact outer coefficient coordinates into physical space."""

    coordinate_gramian = central_coordinate_gramian(degree, amplitude)
    eigenvalues, eigenvectors = np.linalg.eigh(coordinate_gramian)
    square_root = (
        eigenvectors * np.sqrt(eigenvalues)
    ) @ eigenvectors.T
    outer_indices = [0, degree, 2 * degree]
    embedding = np.eye(2 * degree + 1)[:, outer_indices]
    outer_gramian = coordinate_gramian[
        np.ix_(outer_indices, outer_indices)
    ]
    outer_values, outer_vectors = np.linalg.eigh(outer_gramian)
    outer_inverse_square_root = (
        outer_vectors * (1 / np.sqrt(outer_values))
    ) @ outer_vectors.T
    isometry = square_root @ embedding @ outer_inverse_square_root
    if np.linalg.norm(isometry.conj().T @ isometry - np.eye(3), 2) > 2e-13:
        raise AssertionError("the outer coordinate map is not isometric")
    return np.asarray(isometry, dtype=complex)


def normalized_residual(value: np.ndarray, scale: float = 1.0) -> float:
    """Return a spectral-norm residual with a nonvanishing scale."""

    return float(np.linalg.norm(value, 2) / max(1.0, abs(scale)))


def polynomial_matrix_value(
    coefficients: np.ndarray,
    matrix: np.ndarray,
) -> np.ndarray:
    """Evaluate an ascending-coefficient polynomial by Horner's rule."""

    identity = np.eye(matrix.shape[0], dtype=complex)
    value = np.zeros_like(matrix, dtype=complex)
    for coefficient in coefficients[::-1]:
        value = value @ matrix + coefficient * identity
    return value


def critical_multiplier(
    operator: np.ndarray,
    zeros: tuple[complex, ...],
    functions: list[np.ndarray],
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Construct ``F=Q/D`` from the outside critical factor of ``B``."""

    factorization = real_blaschke_critical_factor(zeros)
    outer_factor = factorization.outer_factor
    denominator = factorization.denominator

    numerator_value = polynomial_matrix_value(outer_factor, operator)
    denominator_value = polynomial_matrix_value(denominator, operator)
    multiplier = numerator_value @ np.linalg.inv(denominator_value)
    model_matrix = np.column_stack(
        [function.reshape(-1) for function in functions]
    )
    coefficients = np.linalg.lstsq(
        model_matrix,
        multiplier.reshape(-1),
        rcond=None,
    )[0]
    reconstructed = sum(
        (
            coefficient * function
            for coefficient, function in zip(coefficients, functions)
        ),
        start=np.zeros_like(operator, dtype=complex),
    )
    model_residual = normalized_residual(
        multiplier - reconstructed,
        np.linalg.norm(multiplier, 2),
    )
    return (
        multiplier,
        coefficients,
        factorization.factorization_residual,
        model_residual,
    )


def fiber_cross_residuals(
    functions: list[np.ndarray],
    outer_map: np.ndarray,
    outer_kernel: np.ndarray,
    lifted_kernel: np.ndarray,
    inner_map: np.ndarray,
) -> tuple[float, float, np.ndarray, float]:
    """Return kernel alignment, cross forcing, forcing, and outer mass."""

    dimension = lifted_kernel.size
    forcing = np.zeros((dimension, dimension), dtype=complex)
    coefficient_sum = 0.0
    maximum_alignment = 0.0
    for function in functions:
        fiber = function.conj().T @ lifted_kernel
        outer_fiber = outer_map.conj().T @ fiber
        coefficient = np.vdot(outer_kernel, outer_fiber)
        inner_fiber = fiber - coefficient * outer_map @ outer_kernel
        maximum_alignment = max(
            maximum_alignment,
            float(
                np.linalg.norm(
                    outer_fiber - coefficient * outer_kernel
                )
            ),
        )
        forcing += np.outer(fiber, fiber.conj())
        coefficient_sum += abs(coefficient) ** 2
    forcing_cross = outer_map.conj().T @ forcing @ inner_map
    return (
        maximum_alignment,
        normalized_residual(forcing_cross, np.linalg.norm(forcing, 2)),
        forcing,
        coefficient_sum,
    )


def random_dual_audit(
    functions: list[np.ndarray],
    outer_map: np.ndarray,
    inner_map: np.ndarray,
    seed: int,
    trials: int = 3,
) -> tuple[float, float]:
    """Test cancellation for unrelated deterministic rank-two duals."""

    rng = np.random.default_rng(seed)
    maximum_alignment = 0.0
    maximum_cross = 0.0
    for _ in range(trials):
        factor = (
            rng.standard_normal((3, 2))
            + 1j * rng.standard_normal((3, 2))
        )
        outer_dual = factor @ factor.conj().T
        outer_values, outer_vectors = np.linalg.eigh(outer_dual)
        outer_kernel = outer_vectors[:, 0]
        lifted_dual = dual_transfer(
            functions,
            outer_map @ outer_dual @ outer_map.conj().T,
        )
        lifted_values, lifted_vectors = np.linalg.eigh(
            (lifted_dual + lifted_dual.conj().T) / 2
        )
        if abs(lifted_values[0]) > 2e-10 or lifted_values[1] <= 1e-8:
            raise AssertionError(
                "an unrelated rank-two dual did not lift to corank one"
            )
        alignment, cross, _, _ = fiber_cross_residuals(
            functions,
            outer_map,
            outer_kernel,
            lifted_vectors[:, 0],
            inner_map,
        )
        maximum_alignment = max(maximum_alignment, alignment)
        maximum_cross = max(maximum_cross, cross)
    return maximum_alignment, maximum_cross


def make_record(
    degree: int,
    ellipse_parameter: float,
    amplitude: float,
) -> CentralDualLiftRecord:
    """Construct and audit one deterministic central dual lift."""

    dimension = 2 * degree + 1
    descended_parameter = ellipse_parameter**degree
    operator = central_operator(degree, ellipse_parameter, amplitude)
    outer_operator = central_operator(1, descended_parameter, amplitude)
    outer_condition, _, outer_metric = optimized_rank_one_metric(
        outer_operator
    )

    zeros = tuple(
        complex(zero)
        for zero in chebyshev_blaschke_zeros(
            degree + 1,
            ellipse_parameter,
        )
    )
    functions, descended_operator = model_functions(operator, zeros)
    outer_map = outer_isometry(degree, amplitude)
    outer_image = outer_map.conj().T @ descended_operator @ outer_map
    phase = np.vdot(outer_operator, outer_image) / np.vdot(
        outer_operator, outer_operator
    )
    outer_image_residual = normalized_residual(
        outer_image - phase * outer_operator,
        np.linalg.norm(outer_operator, 2),
    )
    (
        multiplier,
        multiplier_coefficients,
        critical_factorization_residual,
        model_multiplier_residual,
    ) = critical_multiplier(operator, zeros, functions)
    multiplier_norm_square = float(
        np.vdot(multiplier_coefficients, multiplier_coefficients).real
    )
    explicit_lift = np.linalg.solve(
        multiplier.conj().T,
        outer_map,
    )
    scalar_compression_residual = max(
        normalized_residual(
            outer_map.conj().T
            @ function.conj().T
            @ explicit_lift
            - (
                np.conjugate(coefficient) / multiplier_norm_square
            )
            * np.eye(3),
        )
        for function, coefficient in zip(
            functions,
            multiplier_coefficients,
        )
    )

    metric_values, metric_vectors = np.linalg.eigh(outer_metric)
    lower_vector = metric_vectors[:, 0]
    upper_vector = metric_vectors[:, -1]
    dual_difference = (
        np.outer(upper_vector, upper_vector.conj())
        - outer_condition
        * np.outer(lower_vector, lower_vector.conj())
    )
    outer_dual = solve_discrete_lyapunov(
        outer_operator,
        dual_difference,
    )
    outer_dual = (outer_dual + outer_dual.conj().T) / 2
    outer_dual_values, outer_dual_vectors = np.linalg.eigh(outer_dual)
    outer_kernel = outer_dual_vectors[:, 0]
    if outer_dual_values[1] <= 1e-8:
        raise AssertionError("the outer dual did not have numerical rank two")

    embedded_dual = outer_map @ outer_dual @ outer_map.conj().T
    lifted_dual = dual_transfer(functions, embedded_dual)
    lifted_dual = (lifted_dual + lifted_dual.conj().T) / 2
    dual_values, dual_vectors = np.linalg.eigh(lifted_dual)
    explicit_defect = explicit_lift @ outer_kernel
    explicit_defect /= np.linalg.norm(explicit_defect)
    explicit_dual_kernel_residual = normalized_residual(
        lifted_dual @ explicit_defect,
        np.linalg.norm(lifted_dual, 2),
    )
    explicit_kernel_angle_residual = float(
        np.sqrt(
            max(
                0.0,
                1.0
                - abs(np.vdot(dual_vectors[:, 0], explicit_defect)) ** 2,
            )
        )
    )
    defect = explicit_defect
    weighted_reconstruction = sum(
        (
            np.conjugate(coefficient)
            * function.conj().T
            @ defect
            for coefficient, function in zip(
                multiplier_coefficients,
                functions,
            )
        ),
        start=np.zeros(dimension, dtype=complex),
    )
    explicit_scale = np.vdot(
        explicit_lift @ outer_kernel,
        defect,
    ) / np.vdot(
        explicit_lift @ outer_kernel,
        explicit_lift @ outer_kernel,
    )
    expected_reconstruction = explicit_scale * outer_map @ outer_kernel
    weighted_reconstruction_residual = normalized_residual(
        weighted_reconstruction - expected_reconstruction,
    )

    inner_map = null_space(outer_map.conj().T)
    outer_defect_metric = solve_discrete_lyapunov(
        outer_operator.conj().T,
        np.outer(outer_kernel, outer_kernel.conj()),
    )
    maximum_alignment, forcing_cross_residual, forcing, coefficient_sum = (
        fiber_cross_residuals(
            functions,
            outer_map,
            outer_kernel,
            defect,
            inner_map,
        )
    )
    random_alignment, random_cross = random_dual_audit(
        functions,
        outer_map,
        inner_map,
        seed=70_223 + degree,
    )
    outer_projection = outer_map @ outer_map.conj().T
    metric = solve_discrete_lyapunov(
        operator.conj().T,
        np.outer(defect, defect.conj()),
    )
    metric_cross = outer_map.conj().T @ metric @ inner_map
    outer_metric = outer_map.conj().T @ metric @ outer_map
    expected_outer_metric = coefficient_sum * outer_defect_metric
    constructed_values = np.linalg.eigvalsh(metric)
    constructed_condition = float(
        constructed_values[-1] / constructed_values[0]
    )

    record = CentralDualLiftRecord(
        degree=degree,
        dimension=dimension,
        ellipse_parameter=ellipse_parameter,
        descended_parameter=descended_parameter,
        amplitude=amplitude,
        outer_image_residual=outer_image_residual,
        critical_factorization_residual=(
            critical_factorization_residual
        ),
        model_multiplier_residual=model_multiplier_residual,
        scalar_compression_residual=scalar_compression_residual,
        dual_smallest_eigenvalue=float(abs(dual_values[0])),
        dual_second_eigenvalue=float(dual_values[1]),
        explicit_dual_kernel_residual=explicit_dual_kernel_residual,
        explicit_kernel_angle_residual=explicit_kernel_angle_residual,
        weighted_reconstruction_residual=(
            weighted_reconstruction_residual
        ),
        kernel_alignment_residual=maximum_alignment,
        forcing_cross_residual=forcing_cross_residual,
        random_dual_kernel_alignment_residual=random_alignment,
        random_dual_forcing_cross_residual=random_cross,
        metric_cross_residual=normalized_residual(
            metric_cross, np.linalg.norm(metric, 2)
        ),
        outer_metric_residual=normalized_residual(
            outer_metric - expected_outer_metric,
            np.linalg.norm(expected_outer_metric, 2),
        ),
        constructed_condition=constructed_condition,
        outer_condition=outer_condition,
        condition_gap=constructed_condition - outer_condition,
    )

    if record.outer_image_residual > 2e-11:
        raise AssertionError("the Blaschke image missed the outer block")
    if record.critical_factorization_residual > 2e-10:
        raise AssertionError("the outer critical factorization failed")
    if record.model_multiplier_residual > 2e-10:
        raise AssertionError("the critical multiplier left the model space")
    if record.scalar_compression_residual > 2e-10:
        raise AssertionError("the constant fiber trace identity failed")
    if record.dual_smallest_eigenvalue > 2e-10:
        raise AssertionError("the transferred dual lost its kernel")
    if record.dual_second_eigenvalue <= 1e-8:
        raise AssertionError("the transferred dual gained extra kernel")
    if record.explicit_dual_kernel_residual > 2e-10:
        raise AssertionError("the explicit defect missed the dual kernel")
    if record.explicit_kernel_angle_residual > 2e-7:
        raise AssertionError("the explicit and numerical kernels separated")
    if record.weighted_reconstruction_residual > 2e-10:
        raise AssertionError("the multiplier reconstruction identity failed")
    if record.kernel_alignment_residual > 2e-10:
        raise AssertionError("the model-space fibers lost kernel alignment")
    if record.forcing_cross_residual > 2e-10:
        raise AssertionError("the cross-fiber cancellation failed")
    if record.random_dual_kernel_alignment_residual > 2e-10:
        raise AssertionError("an unrelated dual lost kernel alignment")
    if record.random_dual_forcing_cross_residual > 2e-10:
        raise AssertionError("an unrelated dual lost cross-fiber cancellation")
    if record.metric_cross_residual > 2e-10:
        raise AssertionError("the lifted metric did not reduce the outer block")
    if record.outer_metric_residual > 2e-9:
        raise AssertionError("the outer metric was not the descended metric")
    if abs(record.condition_gap) > 2e-8:
        raise AssertionError("the constructed and descended conditions separated")

    # Keep this explicit: the complement is defined as the orthogonal
    # complement of the outer isometry, not by coordinate index assumptions.
    if np.linalg.norm(outer_projection @ inner_map, 2) > 2e-13:
        raise AssertionError("the computed inner basis is not orthogonal")
    return record


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--degrees",
        type=int,
        nargs="+",
        default=DEFAULT_DEGREES,
    )
    parser.add_argument(
        "--ellipse-parameters",
        type=float,
        nargs="+",
        default=DEFAULT_ELLIPSE_PARAMETERS,
    )
    parser.add_argument("--amplitude", type=float, default=DEFAULT_AMPLITUDE)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the deterministic dual-lift audit."""

    args = parse_args()
    records = [
        make_record(degree, ellipse_parameter, args.amplitude)
        for degree in args.degrees
        for ellipse_parameter in args.ellipse_parameters
    ]
    lines = [
        json.dumps(asdict(record), sort_keys=True)
        for record in records
    ]
    for line in lines:
        print(line, flush=True)
    if args.output is not None:
        args.output.write_text(
            "".join(f"{line}\n" for line in lines),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
