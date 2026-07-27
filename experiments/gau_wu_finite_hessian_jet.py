#!/usr/bin/env python3
"""Compute arbitrary finite Gau--Wu Hessians from the analytic Riemann jet."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from gau_wu_finite_model import (
    disk_model_residuals,
    expected_tangent_rank,
    extremal_zeros,
    gau_wu_model,
    normal_basis,
    random_interior_zeros,
    real_matrix,
)
from theodorsen import hilbert_periodic


@dataclass(frozen=True)
class FiniteGauWuJetRecord:
    """One support/Riemann/Blaschke jet Hessian record."""

    dimension: int
    sample: int
    seed: int
    interior_zeros: tuple[str, ...]
    angle_count: int
    tangent_rank: int
    normal_dimension: int
    support_error: float
    functional_error: float
    minimum_support_gap: float
    maximum_joint_first_variation: float
    maximum_zero_hessian_eigenvalue: float
    minimum_zero_hessian_eigenvalue: float
    zero_hessian_condition_number: float
    optimized_eigenvalues: tuple[float, ...]
    maximum_optimized_eigenvalue: float
    sign_classification: str
    all_checks_passed: bool


def matrix_polynomial(
    coefficients: np.ndarray,
    matrix: np.ndarray,
) -> np.ndarray:
    """Evaluate a scalar polynomial at a matrix by Horner's rule."""

    identity = np.eye(matrix.shape[0])
    value = coefficients[-1] * identity
    for coefficient in coefficients[-2::-1]:
        value = value @ matrix + coefficient * identity
    return value


def frechet_polynomial(
    coefficients: np.ndarray,
    matrix: np.ndarray,
    direction: np.ndarray,
) -> np.ndarray:
    """Evaluate a polynomial Frechet derivative with one block calculation."""

    dimension = matrix.shape[0]
    block = np.zeros((2 * dimension, 2 * dimension), dtype=complex)
    block[:dimension, :dimension] = matrix
    block[:dimension, dimension:] = direction
    block[dimension:, dimension:] = matrix
    return matrix_polynomial(coefficients, block)[:dimension, dimension:]


def schwarz_map_coefficients(boundary_values: np.ndarray) -> np.ndarray:
    """Return coefficients of w times the Schwarz transform of real data."""

    count = len(boundary_values)
    fourier = np.fft.fft(boundary_values) / count
    maximum_degree = count // 2
    coefficients = np.zeros(maximum_degree + 1, dtype=complex)
    coefficients[1] = fourier[0]
    coefficients[2:] = 2 * fourier[1:maximum_degree]
    return coefficients


def derivative_coefficients(coefficients: np.ndarray) -> np.ndarray:
    """Differentiate scalar polynomial coefficients."""

    return np.asarray(
        [
            degree * coefficients[degree]
            for degree in range(1, len(coefficients))
        ]
    )


def support_jet(
    matrix: np.ndarray,
    directions: list[np.ndarray],
    angle_count: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Return first and polarized second support coefficients."""

    direction_count = len(directions)
    first = np.empty((direction_count, angle_count))
    second = np.empty((direction_count, direction_count, angle_count))
    minimum_gap = np.inf
    for angle_index, angle in enumerate(
        np.linspace(0, 2 * np.pi, angle_count, endpoint=False)
    ):
        phase = np.exp(-1j * angle)
        support = (phase * matrix + np.conj(phase) * matrix.conj().T) / 2
        eigenvalues, eigenvectors = np.linalg.eigh(support)
        top_vector = eigenvectors[:, -1]
        gaps = 1 - eigenvalues[:-1]
        minimum_gap = min(minimum_gap, gaps[-1])

        projected = np.empty(
            (matrix.shape[0] - 1, direction_count),
            dtype=complex,
        )
        for index, direction in enumerate(directions):
            support_direction = (
                phase * direction
                + np.conj(phase) * direction.conj().T
            ) / 2
            image = support_direction @ top_vector
            first[index, angle_index] = np.vdot(
                top_vector,
                image,
            ).real
            projected[:, index] = (
                eigenvectors[:, :-1].conj().T @ image / np.sqrt(gaps)
            )
        second[:, :, angle_index] = (
            projected.conj().T @ projected
        ).real
    return first, second, float(minimum_gap)


def normalized_operator_jet(
    matrix: np.ndarray,
    directions: list[np.ndarray],
    angle_count: int,
) -> tuple[list[np.ndarray], np.ndarray, float]:
    """Return first and polarized second inverse-Riemann operator jets."""

    first_support, second_support, minimum_gap = support_jet(
        matrix,
        directions,
        angle_count,
    )
    direction_count = len(directions)
    first_coefficients = [
        schwarz_map_coefficients(values) for values in first_support
    ]
    first_values = [
        matrix_polynomial(coefficients, matrix)
        for coefficients in first_coefficients
    ]
    first_derivatives = [
        matrix_polynomial(derivative_coefficients(coefficients), matrix)
        for coefficients in first_coefficients
    ]
    first_operator = [
        direction - value
        for direction, value in zip(directions, first_values, strict=True)
    ]

    angle_shift = np.empty_like(first_support)
    frequencies = np.fft.fftfreq(angle_count, d=1 / angle_count)
    for index, values in enumerate(first_support):
        derivative = np.fft.ifft(
            1j * frequencies * np.fft.fft(values)
        ).real
        angle_shift[index] = hilbert_periodic(values) - derivative

    frechet = [
        [
            frechet_polynomial(coefficients, matrix, direction)
            for direction in directions
        ]
        for coefficients in first_coefficients
    ]
    second_operator = np.empty(
        (
            direction_count,
            direction_count,
            matrix.shape[0],
            matrix.shape[1],
        ),
        dtype=complex,
    )
    for row in range(direction_count):
        for column in range(row, direction_count):
            normal_data = (
                second_support[row, column]
                - angle_shift[row] * angle_shift[column] / 2
            )
            second_conformal = matrix_polynomial(
                schwarz_map_coefficients(normal_data),
                matrix,
            )
            product = (
                first_derivatives[row] @ first_values[column]
                + first_derivatives[column] @ first_values[row]
            ) / 2
            value = (
                -(frechet[row][column] + frechet[column][row]) / 2
                + product
                - second_conformal
            )
            second_operator[row, column] = value
            second_operator[column, row] = value
    return first_operator, second_operator, minimum_gap


def jet_multiply(
    left: list[np.ndarray],
    right: list[np.ndarray],
) -> list[np.ndarray]:
    """Multiply two matrix jets through second order."""

    return [
        sum(
            (
                left[index] @ right[degree - index]
                for index in range(degree + 1)
            ),
            np.zeros_like(left[0]),
        )
        for degree in range(3)
    ]


def jet_inverse(value: list[np.ndarray]) -> list[np.ndarray]:
    """Invert a matrix jet through second order."""

    inverse_zero = np.linalg.inv(value[0])
    inverse_one = -inverse_zero @ value[1] @ inverse_zero
    inverse_two = -inverse_zero @ (
        value[1] @ inverse_one + value[2] @ inverse_zero
    )
    return [inverse_zero, inverse_one, inverse_two]


def blaschke_image_jet(
    matrix: np.ndarray,
    first_operator: np.ndarray,
    second_operator: np.ndarray,
    zeros: np.ndarray,
    velocities: np.ndarray,
) -> list[np.ndarray]:
    """Return the Blaschke image jet for moving operator and zeros."""

    identity = np.eye(matrix.shape[0])
    image = [
        identity.astype(complex),
        np.zeros_like(matrix),
        np.zeros_like(matrix),
    ]
    for zero, velocity in zip(zeros, velocities, strict=True):
        numerator = [
            matrix - zero * identity,
            first_operator - velocity * identity,
            second_operator,
        ]
        denominator = [
            identity - np.conj(zero) * matrix,
            -np.conj(zero) * first_operator
            - np.conj(velocity) * matrix,
            -np.conj(zero) * second_operator
            - np.conj(velocity) * first_operator,
        ]
        image = jet_multiply(
            image,
            jet_multiply(numerator, jet_inverse(denominator)),
        )
    return image


def norm_coefficient(image_jet: list[np.ndarray]) -> tuple[float, float]:
    """Return first and second top-singular-value coefficients."""

    first_image = image_jet[1]
    second_image = image_jet[2]
    first = float(np.real(first_image[0, -1]))
    second = (
        4 * np.real(second_image[0, -1])
        + np.sum(abs(first_image[:, -1]) ** 2)
        + np.sum(abs(first_image[0, :-1]) ** 2)
    ) / 4
    return first, float(np.real(second))


def joint_hessian(
    matrix: np.ndarray,
    physical_first: list[np.ndarray],
    physical_second: np.ndarray,
    zeros: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Build the joint physical/zero-velocity quadratic coefficient."""

    physical_dimension = len(physical_first)
    zero_count = len(zeros)
    variable_count = physical_dimension + 2 * zero_count

    def coefficient(variable: np.ndarray) -> tuple[float, float]:
        physical = variable[:physical_dimension]
        zero_real = variable[
            physical_dimension : physical_dimension + zero_count
        ]
        zero_imaginary = variable[physical_dimension + zero_count :]
        velocities = zero_real + 1j * zero_imaginary
        first_operator = sum(
            (
                physical[index] * physical_first[index]
                for index in range(physical_dimension)
            ),
            np.zeros_like(matrix),
        )
        second_operator = np.einsum(
            "i,j,ijab->ab",
            physical,
            physical,
            physical_second,
        )
        return norm_coefficient(
            blaschke_image_jet(
                matrix,
                first_operator,
                second_operator,
                zeros,
                velocities,
            )
        )

    basis = np.eye(variable_count)
    diagonal = np.empty(variable_count)
    maximum_first = 0.0
    for index in range(variable_count):
        first, diagonal[index] = coefficient(basis[index])
        maximum_first = max(maximum_first, abs(first))

    hessian = np.diag(diagonal)
    for row in range(variable_count):
        for column in range(row + 1, variable_count):
            first, value = coefficient(basis[row] + basis[column])
            maximum_first = max(maximum_first, abs(first))
            entry = (value - diagonal[row] - diagonal[column]) / 2
            hessian[row, column] = entry
            hessian[column, row] = entry
    return hessian, maximum_first


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> FiniteGauWuJetRecord:
    """Compute one complete optimized normal Hessian from analytic jets."""

    rng = np.random.default_rng(seed)
    interior = random_interior_zeros(dimension, rng)
    matrix = gau_wu_model(interior)
    zeros = extremal_zeros(interior)
    support_error, functional_error, sampled_gap = disk_model_residuals(
        matrix,
        zeros,
    )
    normal, tangent_rank, _ = normal_basis(interior)
    directions = [
        real_matrix(normal[:, index], dimension)
        for index in range(normal.shape[1])
    ]
    physical_first, physical_second, jet_gap = normalized_operator_jet(
        matrix,
        directions,
        angle_count,
    )
    joint, maximum_first = joint_hessian(
        matrix,
        physical_first,
        physical_second,
        zeros,
    )
    physical_dimension = len(directions)
    zero_block = joint[physical_dimension:, physical_dimension:]
    zero_eigenvalues = np.linalg.eigvalsh(zero_block)
    zero_condition = float(
        abs(zero_eigenvalues[0] / zero_eigenvalues[-1])
    )
    optimized = (
        joint[:physical_dimension, :physical_dimension]
        - joint[:physical_dimension, physical_dimension:]
        @ np.linalg.solve(
            zero_block,
            joint[physical_dimension:, :physical_dimension],
        )
    )
    optimized = (optimized + optimized.T) / 2
    optimized_eigenvalues = np.linalg.eigvalsh(optimized)
    maximum_eigenvalue = float(optimized_eigenvalues[-1])
    if maximum_eigenvalue > 2e-7:
        classification = "positive witness"
    elif maximum_eigenvalue < -2e-7:
        classification = "strictly negative"
    else:
        classification = "numerically unresolved"

    checks = (
        tangent_rank == expected_tangent_rank(dimension)
        and support_error < 2e-12
        and functional_error < 2e-12
        and min(sampled_gap, jet_gap) > 1e-3
        and maximum_first < 2e-8
        and zero_eigenvalues[-1] < -1e-10
        and zero_condition < 1e9
    )
    if not checks:
        raise RuntimeError(
            "finite Gau--Wu jet audit failed diagnostics: "
            f"n={dimension}, rank={tangent_rank}, "
            f"support={support_error}, functional={functional_error}, "
            f"gaps={sampled_gap}/{jet_gap}, first={maximum_first}, "
            f"zero_eigenvalues={zero_eigenvalues}, "
            f"zero_condition={zero_condition}"
        )

    return FiniteGauWuJetRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        interior_zeros=tuple(
            f"{zero.real:+.12f}{zero.imag:+.12f}j" for zero in interior
        ),
        angle_count=angle_count,
        tangent_rank=tangent_rank,
        normal_dimension=normal.shape[1],
        support_error=support_error,
        functional_error=functional_error,
        minimum_support_gap=min(sampled_gap, jet_gap),
        maximum_joint_first_variation=maximum_first,
        maximum_zero_hessian_eigenvalue=float(zero_eigenvalues[-1]),
        minimum_zero_hessian_eigenvalue=float(zero_eigenvalues[0]),
        zero_hessian_condition_number=zero_condition,
        optimized_eigenvalues=tuple(map(float, optimized_eigenvalues)),
        maximum_optimized_eigenvalue=maximum_eigenvalue,
        sign_classification=classification,
        all_checks_passed=True,
    )


def write_records(
    records: list[FiniteGauWuJetRecord],
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
    parser.add_argument("--dimensions", type=parse_dimensions, default=(3, 4))
    parser.add_argument("--samples", type=int, default=1)
    parser.add_argument("--seed", type=int, default=70224)
    parser.add_argument("--angle-count", type=int, default=512)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/gau_wu_finite_hessian_jet_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the finite-model analytic-jet Hessian grid."""

    arguments = parse_args()
    records: list[FiniteGauWuJetRecord] = []
    for dimension in arguments.dimensions:
        for sample in range(arguments.samples):
            seed = arguments.seed + 1009 * dimension + sample
            record = audit_model(
                dimension,
                sample,
                seed,
                arguments.angle_count,
            )
            records.append(record)
            print(
                json.dumps(
                    {
                        "dimension": dimension,
                        "sample": sample,
                        "maximum_eigenvalue": (
                            record.maximum_optimized_eigenvalue
                        ),
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
