#!/usr/bin/env python3
"""Falsify a tempting kernel-Gram lower bound for the Gau--Wu shape form."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from gau_wu_conformal_shape_phase import complex_blocks, conformal_shape_map
from gau_wu_finite_hessian_jet import (
    joint_hessian,
    normalized_operator_jet,
    support_jet,
)
from gau_wu_finite_model import (
    extremal_zeros,
    gau_wu_model,
    normal_basis,
    real_matrix,
)


@dataclass(frozen=True)
class KernelGramFalsificationRecord:
    """One resolution-independent counterexample record."""

    dimension: int
    interior_zeros: tuple[str, ...]
    angle_counts: tuple[int, ...]
    first_pivot_errors: tuple[float, ...]
    minimum_candidate_gap_eigenvalues: tuple[float, ...]
    maximum_candidate_gap_eigenvalues: tuple[float, ...]
    minimum_shape_eigenvalues: tuple[float, ...]
    maximum_resolution_drift: float
    all_checks_passed: bool


def model_kernel_gram(interior_zeros: np.ndarray) -> np.ndarray:
    """Return the monomial coefficient Gram of the model kernel."""

    degree = len(interior_zeros)
    numerator = np.asarray([1 + 0j])
    denominator = np.asarray([1 + 0j])
    for zero in interior_zeros:
        numerator = np.convolve(numerator, (-zero, 1))
        denominator = np.convolve(denominator, (1, -np.conj(zero)))

    coefficients = np.zeros(3 * degree + 3, dtype=complex)
    for index in range(len(coefficients)):
        numerator_value = (
            numerator[index] if index < len(numerator) else 0
        )
        previous = sum(
            denominator[offset] * coefficients[index - offset]
            for offset in range(1, min(index + 1, len(denominator)))
        )
        coefficients[index] = numerator_value - previous

    gram = np.empty((degree, degree), dtype=complex)
    for row in range(degree):
        for column in range(degree):
            convolution = sum(
                coefficients[row - offset]
                * np.conj(coefficients[column - offset])
                for offset in range(min(row, column) + 1)
            )
            gram[row, column] = (row == column) - convolution
    return gram


def negative_shape_matrix(
    interior_zeros: np.ndarray,
    angle_count: int,
) -> np.ndarray:
    """Return the candidate negative Hermitian conformal-shape matrix."""

    matrix = gau_wu_model(interior_zeros)
    zeros = extremal_zeros(interior_zeros)
    normal, _, _ = normal_basis(interior_zeros)
    directions = [
        real_matrix(normal[:, index], len(zeros) + 1)
        for index in range(normal.shape[1])
    ]
    first_operator, second_operator, _ = normalized_operator_jet(
        matrix,
        directions,
        angle_count,
    )
    joint, _ = joint_hessian(
        matrix,
        first_operator,
        second_operator,
        zeros,
    )
    physical_dimension = len(directions)
    zero_block = joint[physical_dimension:, physical_dimension:]
    scalar_schur = (
        joint[:physical_dimension, :physical_dimension]
        - joint[:physical_dimension, physical_dimension:]
        @ np.linalg.solve(
            zero_block,
            joint[physical_dimension:, :physical_dimension],
        )
    )
    similarity_form = 4 * scalar_schur

    first_support, _, _ = support_jet(matrix, directions, angle_count)
    shape_map = conformal_shape_map(first_support, matrix.shape[0])
    inverse_shape_form = (
        shape_map @ np.linalg.solve(similarity_form, shape_map.T)
    )
    shape_form = np.linalg.inv(inverse_shape_form)
    shape_form = (shape_form + shape_form.T) / 2
    hermitian, _ = complex_blocks(shape_form, matrix.shape[0] - 1)
    return -(hermitian + hermitian.conj().T) / 2


def audit_witness(angle_counts: tuple[int, ...]) -> KernelGramFalsificationRecord:
    """Evaluate the fixed rational witness at several Fourier resolutions."""

    interior_zeros = np.asarray((3 / 8 + 1j / 4, -1 / 2 - 1j / 4))
    degree = len(interior_zeros)
    derivative = np.diag(np.arange(1, degree + 1))
    kernel_gram = model_kernel_gram(interior_zeros)
    proposed_lower_bound = (
        3 * derivative @ np.linalg.inv(kernel_gram) @ derivative
    )
    expected_first_pivot = 16 * abs(np.prod(interior_zeros)) ** 2

    pivots: list[float] = []
    minimum_gaps: list[float] = []
    maximum_gaps: list[float] = []
    minimum_shapes: list[float] = []
    shapes: list[np.ndarray] = []
    for angle_count in angle_counts:
        shape = negative_shape_matrix(interior_zeros, angle_count)
        shapes.append(shape)
        pivots.append(float(abs(shape[0, 0] - expected_first_pivot)))
        remainder = (
            shape[1:, 1:]
            - np.outer(shape[1:, 0], shape[0, 1:]) / shape[0, 0]
        )
        candidate_gap = remainder - proposed_lower_bound
        candidate_gap = (candidate_gap + candidate_gap.conj().T) / 2
        gap_eigenvalues = np.linalg.eigvalsh(candidate_gap)
        minimum_gaps.append(float(gap_eigenvalues[0]))
        maximum_gaps.append(float(gap_eigenvalues[-1]))
        minimum_shapes.append(float(np.linalg.eigvalsh(shape)[0]))

    reference = shapes[-1]
    maximum_drift = max(
        float(np.linalg.norm(shape - reference, 2)) for shape in shapes
    )
    checks = (
        max(pivots) < 2e-10
        and max(minimum_gaps) < -0.1
        and min(maximum_gaps) > 1
        and min(minimum_shapes) > 0.05
        and maximum_drift < 2e-9
    )
    if not checks:
        raise RuntimeError(
            "kernel-Gram falsification witness was not reproduced: "
            f"pivots={pivots}, minimum_gaps={minimum_gaps}, "
            f"maximum_gaps={maximum_gaps}, "
            f"minimum_shapes={minimum_shapes}, drift={maximum_drift}"
        )

    return KernelGramFalsificationRecord(
        dimension=4,
        interior_zeros=tuple(
            f"{zero.real:+.12f}{zero.imag:+.12f}j"
            for zero in interior_zeros
        ),
        angle_counts=angle_counts,
        first_pivot_errors=tuple(pivots),
        minimum_candidate_gap_eigenvalues=tuple(minimum_gaps),
        maximum_candidate_gap_eigenvalues=tuple(maximum_gaps),
        minimum_shape_eigenvalues=tuple(minimum_shapes),
        maximum_resolution_drift=maximum_drift,
        all_checks_passed=True,
    )


def write_record(
    record: KernelGramFalsificationRecord,
    output: Path,
) -> str:
    """Write the deterministic record atomically and return its SHA-256."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(
        json.dumps(asdict(record), sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(output)
    return hashlib.sha256(output.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/gau_wu_kernel_gram_falsification_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the fixed rational counterexample."""

    arguments = parse_args()
    record = audit_witness((256, 512, 1024))
    digest = write_record(record, arguments.output)
    print(
        json.dumps(
            {
                "all_checks_passed": record.all_checks_passed,
                "minimum_gap": min(
                    record.minimum_candidate_gap_eigenvalues
                ),
                "maximum_resolution_drift": (
                    record.maximum_resolution_drift
                ),
                "sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
