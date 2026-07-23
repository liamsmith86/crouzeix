#!/usr/bin/env python3
"""Audit the seven-real-dimensional affine-unitary slice at C_3."""

from __future__ import annotations

import numpy as np

from crabb_second_order_equality import (
    affine_unitary_generators,
    real_vector,
)
from crouzeix import crabb_matrix
from general_crabb_second_order_modes import full_quadratic_form


def slice_generators() -> np.ndarray:
    """Return the four residual and three negative slice generators."""

    mode_one_real = np.diag([-1 / 3, 2 / 3, -1 / 3]).astype(complex)
    mode_one_real[0, 2] = 1
    mode_one_imaginary = 1j * np.diag([-1 / 3, 2 / 3, -1 / 3])
    mode_one_imaginary[0, 2] = -1j

    mode_two_real = np.zeros((3, 3), dtype=complex)
    mode_two_real[1, 0] = 1
    mode_two_real[2, 1] = 1
    mode_two_imaginary = 1j * mode_two_real

    superdiagonal_difference = np.zeros((3, 3), dtype=complex)
    superdiagonal_difference[0, 1] = 1
    superdiagonal_difference[1, 2] = -1
    bottom_real = np.zeros((3, 3), dtype=complex)
    bottom_real[2, 0] = 1
    bottom_imaginary = 1j * bottom_real

    matrices = (
        mode_one_real,
        mode_one_imaginary,
        mode_two_real,
        mode_two_imaginary,
        superdiagonal_difference,
        bottom_real,
        bottom_imaginary,
    )
    return np.stack([real_vector(matrix) for matrix in matrices], axis=1)


def main() -> None:
    base = crabb_matrix(2)
    orbit = affine_unitary_generators(base)
    slice_basis = slice_generators()
    orbit_rank = np.linalg.matrix_rank(orbit, tol=1e-10)
    slice_rank = np.linalg.matrix_rank(slice_basis, tol=1e-10)
    combined_rank = np.linalg.matrix_rank(
        np.column_stack((orbit, slice_basis)),
        tol=1e-10,
    )
    orthogonality_error = np.linalg.norm(orbit.T @ slice_basis, 2)
    if (orbit_rank, slice_rank, combined_rank) != (11, 7, 18):
        raise AssertionError("the orbit and slice dimensions did not close")
    if orthogonality_error > 1e-12:
        raise AssertionError("the proposed slice was not orbit-orthogonal")

    quadratic = full_quadratic_form(3, 64)
    restricted = slice_basis.T @ quadratic @ slice_basis
    eigenvalues = np.linalg.eigvalsh((restricted + restricted.T) / 2)
    if np.count_nonzero(eigenvalues < -1e-8) != 3:
        raise AssertionError("the slice did not have three negative directions")
    if np.max(eigenvalues) > 1e-10:
        raise AssertionError("the slice quadratic form had a positive direction")

    print("PASS p=3 affine-unitary local-slice audit")
    print(f"orbit/slice/combined ranks = {orbit_rank}/{slice_rank}/{combined_rank}")
    print(f"orthogonality error = {orthogonality_error:.3e}")
    print(f"restricted eigenvalues = {eigenvalues.tolist()}")


if __name__ == "__main__":
    main()
