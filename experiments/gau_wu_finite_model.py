"""Numerical utilities for arbitrary finite Gau--Wu disk models."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np


def gau_wu_model(interior_zeros: Sequence[complex]) -> np.ndarray:
    """Return the Gau--Wu model whose ordered phi-zeros are 0, b_1, ..., 0."""

    zeros = np.asarray((0, *interior_zeros, 0), dtype=complex)
    dimension = len(zeros)
    model = np.zeros((dimension, dimension), dtype=complex)
    np.fill_diagonal(model, zeros)
    for row in range(dimension):
        for column in range(row + 1, dimension):
            middle = zeros[row + 1 : column]
            product = np.prod(np.conj(middle)) if len(middle) else 1
            model[row, column] = (
                (-1) ** (column - row - 1)
                * product
                * np.sqrt(
                    (1 - abs(zeros[row]) ** 2)
                    * (1 - abs(zeros[column]) ** 2)
                )
            )

    scaling = np.ones(dimension)
    scaling[0] = np.sqrt(2)
    scaling[-1] = 1 / np.sqrt(2)
    return scaling[:, None] * model / scaling[None, :]


def extremal_zeros(interior_zeros: Sequence[complex]) -> np.ndarray:
    """Return the zeros of f when the model inner function is phi=z f."""

    return np.asarray((*interior_zeros, 0), dtype=complex)


def blaschke_at_matrix(
    matrix: np.ndarray,
    zeros: Sequence[complex],
) -> np.ndarray:
    """Evaluate the finite Blaschke product with the given zeros."""

    identity = np.eye(matrix.shape[0])
    value = np.eye(matrix.shape[0], dtype=complex)
    for zero in zeros:
        value = (
            value
            @ (matrix - zero * identity)
            @ np.linalg.inv(identity - np.conj(zero) * matrix)
        )
    return value


def pack_zeros(zeros: np.ndarray) -> np.ndarray:
    """Map disk coordinates to unconstrained real optimizer coordinates."""

    points = zeros / np.sqrt(1 - abs(zeros) ** 2)
    return np.concatenate((points.real, points.imag))


def unpack_zeros(parameters: np.ndarray) -> np.ndarray:
    """Map unconstrained real optimizer coordinates into the unit disk."""

    degree = len(parameters) // 2
    points = parameters[:degree] + 1j * parameters[degree:]
    return np.asarray(
        [point / np.sqrt(1 + abs(point) ** 2) for point in points]
    )


def real_vectorize(matrix: np.ndarray) -> np.ndarray:
    """Vectorize a complex matrix into interlaced real coordinates."""

    return np.asarray(
        [
            coordinate
            for entry in matrix.flat
            for coordinate in (entry.real, entry.imag)
        ]
    )


def real_matrix(vector: np.ndarray, dimension: int) -> np.ndarray:
    """Invert :func:`real_vectorize` at the requested matrix dimension."""

    entries = [
        vector[2 * index] + 1j * vector[2 * index + 1]
        for index in range(dimension**2)
    ]
    return np.asarray(entries).reshape(dimension, dimension)


def equality_tangent(
    interior_zeros: Sequence[complex],
    derivative_step: float = 1e-6,
) -> np.ndarray:
    """Return a spanning matrix for affine, unitary, and model-zero tangents."""

    zeros = np.asarray(interior_zeros, dtype=complex)
    matrix = gau_wu_model(zeros)
    dimension = matrix.shape[0]
    tangents: list[np.ndarray] = []

    for index in range(dimension - 1):
        generator = np.zeros_like(matrix)
        generator[index, index] = 1j
        tangents.append(generator @ matrix - matrix @ generator)
    for row in range(dimension):
        for column in range(row + 1, dimension):
            real_generator = np.zeros_like(matrix)
            real_generator[row, column] = 1
            real_generator[column, row] = -1
            tangents.append(
                real_generator @ matrix - matrix @ real_generator
            )

            imaginary_generator = np.zeros_like(matrix)
            imaginary_generator[row, column] = 1j
            imaginary_generator[column, row] = 1j
            tangents.append(
                imaginary_generator @ matrix - matrix @ imaginary_generator
            )

    identity = np.eye(dimension)
    tangents.extend((identity, 1j * identity, matrix, 1j * matrix))

    for index in range(len(zeros)):
        for increment in (derivative_step, 1j * derivative_step):
            plus = zeros.copy()
            minus = zeros.copy()
            plus[index] += increment
            minus[index] -= increment
            tangents.append(
                (gau_wu_model(plus) - gau_wu_model(minus))
                / (2 * derivative_step)
            )

    return np.stack(
        [real_vectorize(tangent) for tangent in tangents],
        axis=1,
    )


def expected_tangent_rank(dimension: int) -> int:
    """Return the generic real dimension of the equality tangent."""

    return dimension**2 + 2 * dimension - 2


def normal_basis(
    interior_zeros: Sequence[complex],
) -> tuple[np.ndarray, int, np.ndarray]:
    """Return an orthonormal normal basis, tangent rank, and singular values."""

    tangent = equality_tangent(interior_zeros)
    left, singular_values, _ = np.linalg.svd(tangent, full_matrices=True)
    threshold = 1e-8 * singular_values[0]
    rank = int(np.sum(singular_values > threshold))
    return left[:, rank:], rank, singular_values


def disk_model_residuals(
    matrix: np.ndarray,
    zeros: Sequence[complex],
    angle_count: int = 720,
) -> tuple[float, float, float]:
    """Return support, functional-value, and support-gap diagnostics."""

    support_error = 0.0
    minimum_gap = np.inf
    for angle in np.linspace(0, 2 * np.pi, angle_count, endpoint=False):
        support = (
            np.exp(-1j * angle) * matrix
            + np.exp(1j * angle) * matrix.conj().T
        ) / 2
        eigenvalues = np.linalg.eigvalsh(support)
        support_error = max(support_error, abs(eigenvalues[-1] - 1))
        minimum_gap = min(minimum_gap, eigenvalues[-1] - eigenvalues[-2])

    expected = np.zeros_like(matrix)
    expected[0, -1] = 2
    functional_error = float(
        np.linalg.norm(blaschke_at_matrix(matrix, zeros) - expected, 2)
    )
    return float(support_error), functional_error, float(minimum_gap)


def random_interior_zeros(
    dimension: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Generate separated, nondegenerate interior model zeros."""

    zeros: list[complex] = []
    maximum_attempts = 10_000 * dimension
    for _ in range(maximum_attempts):
        if len(zeros) == dimension - 2:
            return np.asarray(zeros)
        radius = rng.uniform(0.12, 0.62)
        angle = rng.uniform(-np.pi, np.pi)
        candidate = radius * np.exp(1j * angle)
        if all(abs(candidate - zero) > 0.12 for zero in zeros):
            zeros.append(candidate)
    raise RuntimeError(
        f"could not generate {dimension - 2} separated interior zeros"
    )
