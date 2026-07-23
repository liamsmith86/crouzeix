#!/usr/bin/env python3
"""Probe the mixed Faber--Blaschke dual square on coprime grades.

For a noncentral one-grade equality direction, form

    G = P_L + 2 a (P_k + P_(L-k))
        + lambda a c^k P_(L-k),

map the roots of ``G`` from the ellipse to the disk, and use them as
the zeros of a finite Blaschke product ``B``.  This interpolates between
the characteristic Blaschke product on the disk equality family and the
Chebyshev--Blaschke product on the elliptic Crabb axis.

The observed associated loss is

    (4 - ||B(T(a,c))||^2) / (16 a^2 c^(2k))
        -> 4 + |lambda - 2|^2 / 4.

Thus the explicit mixed correction ``lambda=2`` improves the natural
Faber-root loss from coefficient 80 to the conjecturally sharp 64.
The checker is numerical evidence for A99, not a proof of the limit.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

from mpmath import mp
import numpy as np

from crabb_elliptic_axis import elliptic_modulus_from_nome
from crabb_palindromic_elliptic_face import (
    disk_matrix,
    ellipse_pullback,
    palindromic_direction,
)


DEFAULT_CASES = ((5, 2), (7, 3), (9, 4))
DEFAULT_CORRECTIONS = (-2.0, 0.0, 1.0, 2.0, 3.0, 4.0)


@dataclass(frozen=True)
class FaberBlaschkeDualRecord:
    """One finite-amplitude dual-square regression."""

    length: int
    dimension: int
    grade: int
    ellipse_parameter: float
    amplitude: float
    maximum_zero_modulus: float
    correction_values: tuple[float, ...]
    normalized_losses: tuple[float, ...]
    limiting_predictions: tuple[float, ...]
    maximum_limit_residual: float
    natural_loss_coefficient: float
    corrected_loss_coefficient: float
    predicted_natural_coefficient: float
    predicted_corrected_coefficient: float


def add_polynomials(
    left: np.ndarray,
    right: np.ndarray,
) -> np.ndarray:
    """Add ascending-order coefficient arrays."""

    result = np.zeros(max(len(left), len(right)), dtype=complex)
    result[: len(left)] += left
    result[: len(right)] += right
    return result


def dickson_polynomials(
    length: int,
    ellipse_parameter: float,
) -> tuple[np.ndarray, ...]:
    """Return ``P_0,...,P_L`` in ascending coefficient order."""

    polynomials = [
        np.array([2.0], dtype=complex),
        np.array([0.0, 1.0], dtype=complex),
    ]
    for _ in range(2, length + 1):
        polynomials.append(
            add_polynomials(
                np.concatenate(([0.0], polynomials[-1])),
                -ellipse_parameter * polynomials[-2],
            )
        )
    return tuple(polynomials)


def ellipse_to_disk(value: complex, ellipse_parameter: float) -> complex:
    """Evaluate the normalized elliptic Riemann map at one point."""

    modulus = elliptic_modulus_from_nome(ellipse_parameter**2)
    elliptic_parameter = modulus**2
    quarter_period = mp.ellipk(elliptic_parameter)
    argument = (
        2
        * quarter_period
        / mp.pi
        * mp.asin(complex(value) / (2 * np.sqrt(ellipse_parameter)))
    )
    return complex(
        mp.sqrt(modulus)
        * mp.ellipfun("sn")(argument, elliptic_parameter)
    )


def evaluate_matrix_polynomial(
    coefficients: np.ndarray,
    matrix: np.ndarray,
) -> np.ndarray:
    """Evaluate an ascending-order polynomial by Horner's rule."""

    result = np.zeros_like(matrix, dtype=complex)
    identity = np.eye(matrix.shape[0])
    for coefficient in coefficients[::-1]:
        result = result @ matrix + coefficient * identity
    return result


def faber_polynomial(
    length: int,
    grade: int,
    amplitude: float,
    ellipse_parameter: float,
    correction: float,
) -> np.ndarray:
    """Return the corrected Faber characteristic polynomial ``G``."""

    polynomials = dickson_polynomials(length, ellipse_parameter)
    result = polynomials[length].copy()
    result = add_polynomials(
        result,
        2
        * amplitude
        * add_polynomials(
            polynomials[grade],
            polynomials[length - grade],
        ),
    )
    result = add_polynomials(
        result,
        (
            correction
            * amplitude
            * ellipse_parameter**grade
            * polynomials[length - grade]
        ),
    )
    return result


def blaschke_image(
    operator: np.ndarray,
    length: int,
    grade: int,
    amplitude: float,
    ellipse_parameter: float,
    correction: float,
) -> tuple[np.ndarray, float]:
    """Evaluate the root-defined Blaschke product at ``operator``."""

    polynomial = faber_polynomial(
        length,
        grade,
        amplitude,
        ellipse_parameter,
        correction,
    )
    ellipse_roots = np.roots(polynomial[::-1])
    disk_zeros = np.array(
        [
            ellipse_to_disk(root, ellipse_parameter)
            for root in ellipse_roots
        ]
    )
    numerator = np.poly(disk_zeros)[::-1]
    denominator = np.array([1.0], dtype=complex)
    for zero in disk_zeros:
        denominator = np.convolve(
            denominator,
            np.array([1.0, -np.conj(zero)]),
        )
    image = np.linalg.solve(
        evaluate_matrix_polynomial(denominator, operator),
        evaluate_matrix_polynomial(numerator, operator),
    )
    return image, float(np.max(np.abs(disk_zeros)))


def squared_operator_norm(matrix: np.ndarray) -> float:
    """Return the square of the largest singular value."""

    largest = np.linalg.svd(matrix, compute_uv=False)[0]
    return float(largest**2)


def make_record(
    length: int,
    grade: int,
    ellipse_parameter: float,
    amplitude: float,
    corrections: tuple[float, ...],
) -> FaberBlaschkeDualRecord:
    """Construct one correction scan and validate its qualitative square."""

    if not 1 <= grade < length / 2:
        raise ValueError("the probe requires a noncentral grade")
    coefficients = palindromic_direction(length, grade)
    axis = ellipse_pullback(
        disk_matrix(length + 1, 0.0, coefficients),
        ellipse_parameter,
    )
    axis_image, maximum_zero_modulus = blaschke_image(
        axis,
        length,
        grade,
        0.0,
        ellipse_parameter,
        0.0,
    )
    axis_value = squared_operator_norm(axis_image)

    losses = []
    maximum_zero_moduli = [maximum_zero_modulus]
    for correction in corrections:
        values = []
        for sign in (-1.0, 1.0):
            signed_amplitude = sign * amplitude
            operator = ellipse_pullback(
                disk_matrix(
                    length + 1,
                    signed_amplitude,
                    coefficients,
                ),
                ellipse_parameter,
            )
            image, zero_modulus = blaschke_image(
                operator,
                length,
                grade,
                signed_amplitude,
                ellipse_parameter,
                correction,
            )
            maximum_zero_moduli.append(zero_modulus)
            values.append(squared_operator_norm(image))
        loss = (
            2 * axis_value - values[0] - values[1]
        ) / (2 * amplitude**2)
        losses.append(loss)

    scale = 16 * ellipse_parameter ** (2 * grade)
    normalized = tuple(loss / scale for loss in losses)
    predictions = tuple(
        4 + abs(correction - 2) ** 2 / 4
        for correction in corrections
    )
    maximum_residual = max(
        abs(actual - predicted)
        for actual, predicted in zip(normalized, predictions)
    )
    natural_index = corrections.index(0.0)
    corrected_index = corrections.index(2.0)
    natural = losses[natural_index]
    corrected = losses[corrected_index]
    predicted_natural = 80 * ellipse_parameter ** (2 * grade)
    predicted_corrected = 64 * ellipse_parameter ** (2 * grade)

    if max(maximum_zero_moduli) >= 1:
        raise AssertionError("a Faber root left the disk after mapping")
    if not corrected < natural:
        raise AssertionError("the mixed correction did not improve the dual")
    if maximum_residual > 0.25:
        raise AssertionError(
            "the finite scan left the predicted dual square"
        )

    return FaberBlaschkeDualRecord(
        length=length,
        dimension=length + 1,
        grade=grade,
        ellipse_parameter=ellipse_parameter,
        amplitude=amplitude,
        maximum_zero_modulus=max(maximum_zero_moduli),
        correction_values=corrections,
        normalized_losses=normalized,
        limiting_predictions=predictions,
        maximum_limit_residual=maximum_residual,
        natural_loss_coefficient=natural,
        corrected_loss_coefficient=corrected,
        predicted_natural_coefficient=predicted_natural,
        predicted_corrected_coefficient=predicted_corrected,
    )


def parse_cases(values: list[str]) -> tuple[tuple[int, int], ...]:
    """Parse ``L:k`` case specifications."""

    result = []
    for value in values:
        length_text, grade_text = value.split(":", maxsplit=1)
        result.append((int(length_text), int(grade_text)))
    return tuple(result)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cases",
        nargs="+",
        default=[f"{length}:{grade}" for length, grade in DEFAULT_CASES],
    )
    parser.add_argument("--ellipse-parameter", type=float, default=0.18)
    parser.add_argument("--amplitude", type=float, default=0.01)
    parser.add_argument(
        "--corrections",
        nargs="+",
        type=float,
        default=DEFAULT_CORRECTIONS,
    )
    parser.add_argument("--precision", type=int, default=60)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the deterministic Faber--Blaschke correction scan."""

    args = parse_args()
    cases = parse_cases(args.cases)
    corrections = tuple(args.corrections)
    if 0.0 not in corrections or 2.0 not in corrections:
        raise ValueError("the correction scan must contain 0 and 2")
    if not 0 < args.ellipse_parameter < 1:
        raise ValueError("the ellipse parameter must lie in (0,1)")
    if not 0 < args.amplitude < 0.05:
        raise ValueError("the amplitude must lie in (0,.05)")
    if args.precision < 40:
        raise ValueError("use at least 40 digits")
    mp.dps = args.precision

    records = [
        make_record(
            length,
            grade,
            args.ellipse_parameter,
            args.amplitude,
            corrections,
        )
        for length, grade in cases
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
