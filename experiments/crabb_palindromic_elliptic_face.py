#!/usr/bin/env python3
"""Probe the elliptic Newton face along palindromic disk equality branches.

For a phase-one palindromic Toeplitz direction ``z`` and small real ``a``,
construct the exact disk matrix ``X(a z)`` from L122.  Its elliptic image

    X(a z) + c X(a z)^*

has numerical range with semiaxes ``1+c`` and ``1-c``.  This script applies
the exact scalar ellipse-to-disk map by eigenvalue functional calculus and
compares the locally optimized rank-one Stein envelope with the full L21
similarity SDP.

The conjectured joint leading face is

    Gamma = Gamma_axis - 64 a^2 sum_k |z_k|^2 c^(2k) + higher terms.

The output is exploratory evidence and a regression target, not a proof of
that all-size asymptotic formula.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

from mpmath import mp
import numpy as np
from scipy.linalg import solve_discrete_lyapunov
from scipy.optimize import minimize

from crabb_disk_toeplitz_quartic import shift
from crabb_elliptic_axis import elliptic_modulus_from_nome
from slice_cb_sdp import solve_similarity_sdp


DEFAULT_ELLIPSE_PARAMETERS = (0.12, 0.20)
DEFAULT_AMPLITUDES = (0.03, 0.06)


@dataclass(frozen=True)
class EllipticFaceRecord:
    dimension: int
    first_offset: int
    ellipse_parameter: float
    smaller_amplitude: float
    larger_amplitude: float
    axis_bound: float
    larger_amplitude_bound: float
    estimated_quadratic_coefficient: float
    predicted_leading_coefficient: float
    prediction_ratio: float
    quadratic_scaling_ratio: float
    expected_scaling_ratio: float
    maximum_rank_one_sdp_gap: float


def palindromic_direction(
    length: int,
    first_offset: int,
) -> np.ndarray:
    """Return a real phase-one direction supported on one reversed pair."""

    if not 1 <= first_offset < length:
        raise ValueError("the offset must satisfy 1 <= offset < length")
    reverse_offset = length - first_offset
    coefficients = np.zeros(length - 1)
    coefficients[first_offset - 1] = 1
    coefficients[reverse_offset - 1] = 1
    return coefficients


def disk_matrix(
    dimension: int,
    amplitude: float,
    coefficients: np.ndarray,
) -> np.ndarray:
    """Construct the exact normalized disk chart matrix."""

    length = dimension - 1
    toeplitz = np.eye(length) / 2
    for offset, coefficient in enumerate(coefficients, start=1):
        for row in range(length - offset):
            toeplitz[row, row + offset] = amplitude * coefficient
            toeplitz[row + offset, row] = amplitude * coefficient

    chart = np.asarray(extend_matrix(toeplitz))
    nilpotent_shift = np.asarray(shift(dimension), dtype=float)
    coordinate_gramian = (
        chart + nilpotent_shift.T @ chart @ nilpotent_shift
    )
    eigenvalues, eigenvectors = np.linalg.eigh(coordinate_gramian)
    if eigenvalues[0] <= 0:
        raise ValueError("the requested Toeplitz chart left the positive cone")
    inverse_square_root = (
        eigenvectors * (1 / np.sqrt(eigenvalues))
    ) @ eigenvectors.T
    return (
        2
        * inverse_square_root
        @ chart
        @ nilpotent_shift
        @ inverse_square_root
    )


def extend_matrix(matrix: np.ndarray) -> np.ndarray:
    """Append one zero row and column without converting through SymPy."""

    result = np.zeros((matrix.shape[0] + 1, matrix.shape[1] + 1))
    result[:-1, :-1] = matrix
    return result


def ellipse_pullback(
    disk_operator: np.ndarray,
    ellipse_parameter: float,
) -> np.ndarray:
    """Apply the normalized ellipse-to-disk map by spectral calculus."""

    parameter = ellipse_parameter
    ellipse_operator = (
        disk_operator + parameter * disk_operator.T
    )
    modulus = elliptic_modulus_from_nome(parameter**2)
    elliptic_parameter = modulus**2
    quarter_period = mp.ellipk(elliptic_parameter)
    sn = mp.ellipfun("sn")

    eigenvalues, eigenvectors = np.linalg.eig(ellipse_operator)
    mapped = []
    for eigenvalue in eigenvalues:
        argument = (
            (2 * quarter_period / mp.pi)
            * mp.asin(complex(eigenvalue) / (2 * np.sqrt(parameter)))
        )
        mapped.append(
            complex(mp.sqrt(modulus) * sn(argument, elliptic_parameter))
        )
    operator = (
        eigenvectors
        @ np.diag(mapped)
        @ np.linalg.inv(eigenvectors)
    )
    return np.real_if_close(operator, tol=1_000).real


def rank_one_envelope(
    operator: np.ndarray,
    initial: np.ndarray | None = None,
) -> tuple[float, np.ndarray]:
    """Minimize the real rank-one Stein branch with first defect entry one."""

    dimension = operator.shape[0]
    if initial is None:
        initial = np.zeros(dimension - 1)

    def objective(tail: np.ndarray) -> float:
        defect = np.concatenate(([1.0], tail))
        gramian = solve_discrete_lyapunov(
            operator.T,
            np.outer(defect, defect),
        )
        eigenvalues = np.linalg.eigvalsh(gramian)
        return float(eigenvalues[-1] / eigenvalues[0])

    result = minimize(
        objective,
        initial,
        method="BFGS",
        options={"gtol": 1e-10, "maxiter": 1_000},
    )
    if not np.isfinite(result.fun):
        raise RuntimeError("rank-one envelope optimization failed")
    return float(result.fun), np.asarray(result.x)


def make_record(
    dimension: int,
    first_offset: int,
    ellipse_parameter: float,
    amplitudes: tuple[float, float],
) -> EllipticFaceRecord:
    """Measure one equality branch and compare the two optimizations."""

    length = dimension - 1
    coefficients = palindromic_direction(length, first_offset)
    sample_amplitudes = (0.0, *amplitudes)
    rank_one_bounds: list[float] = []
    sdp_bounds: list[float] = []
    initial: np.ndarray | None = None
    for amplitude in sample_amplitudes:
        disk_operator = disk_matrix(
            dimension,
            amplitude,
            coefficients,
        )
        operator = ellipse_pullback(disk_operator, ellipse_parameter)
        rank_one_bound, initial = rank_one_envelope(operator, initial)
        rank_one_bounds.append(rank_one_bound)
        sdp_bounds.append(solve_similarity_sdp(operator).bound)

    smaller_amplitude, larger_amplitude = amplitudes
    estimated_coefficient = (
        rank_one_bounds[0] - rank_one_bounds[2]
    ) / larger_amplitude**2
    predicted_coefficient = 64 * sum(
        abs(coefficient) ** 2 * ellipse_parameter ** (2 * offset)
        for offset, coefficient in enumerate(coefficients, start=1)
    )
    quadratic_scaling_ratio = (
        (rank_one_bounds[0] - rank_one_bounds[1])
        / (rank_one_bounds[0] - rank_one_bounds[2])
    )
    expected_scaling_ratio = (
        smaller_amplitude / larger_amplitude
    ) ** 2
    maximum_gap = max(
        abs(rank_one - sdp)
        for rank_one, sdp in zip(rank_one_bounds, sdp_bounds)
    )

    if estimated_coefficient <= 0:
        raise RuntimeError("the sampled elliptic equality face was not strict")
    if abs(quadratic_scaling_ratio - expected_scaling_ratio) > 0.02:
        raise RuntimeError("the sampled equality face was not quadratic")
    if maximum_gap > 2e-6:
        raise RuntimeError("rank-one and full SDP branches separated")

    return EllipticFaceRecord(
        dimension=dimension,
        first_offset=first_offset,
        ellipse_parameter=ellipse_parameter,
        smaller_amplitude=smaller_amplitude,
        larger_amplitude=larger_amplitude,
        axis_bound=rank_one_bounds[0],
        larger_amplitude_bound=rank_one_bounds[2],
        estimated_quadratic_coefficient=estimated_coefficient,
        predicted_leading_coefficient=predicted_coefficient,
        prediction_ratio=estimated_coefficient / predicted_coefficient,
        quadratic_scaling_ratio=quadratic_scaling_ratio,
        expected_scaling_ratio=expected_scaling_ratio,
        maximum_rank_one_sdp_gap=maximum_gap,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=3)
    parser.add_argument("--maximum-size", type=int, default=7)
    parser.add_argument(
        "--ellipse-parameters",
        type=float,
        nargs="+",
        default=DEFAULT_ELLIPSE_PARAMETERS,
    )
    parser.add_argument(
        "--amplitudes",
        type=float,
        nargs=2,
        default=DEFAULT_AMPLITUDES,
    )
    parser.add_argument("--precision", type=int, default=50)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the scaling probes and optionally persist JSONL output."""

    args = parse_args()
    if args.minimum_size < 3 or args.maximum_size < args.minimum_size:
        raise ValueError("size range must satisfy 3 <= minimum <= maximum")
    if not 0 < args.amplitudes[0] < args.amplitudes[1]:
        raise ValueError("amplitudes must be positive and increasing")
    mp.dps = args.precision

    records = [
        make_record(
            dimension,
            first_offset,
            ellipse_parameter,
            tuple(args.amplitudes),
        )
        for dimension in range(args.minimum_size, args.maximum_size + 1)
        for first_offset in range(1, (dimension - 1) // 2 + 1)
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
