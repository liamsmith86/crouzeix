#!/usr/bin/env python3
"""Audit the elliptic lift of the full Hardy equality manifold.

L187 replaces the old phase-palindromic equality cone by the analytic
manifold ``Psi(H)=0`` in the general Hermitian disk chart.  This
checker tests the two new interfaces needed by the final local merger:

1. the nonconstant characteristic coefficients are local coordinates
   on that manifold and obey the unrestricted Faber reflection law;
2. a first Hardy residual is orthogonal to the leading reflected
   elliptic leg in the optimized rank-one Stein envelope.

The first interface is checked on exact formal equality jets.  The
second is an exact weighted-series calculation, including the Riemann
map and optimized defect jets.  A floating nonlinear check separately
tests genuinely non-Toeplitz equality anchors and several ellipse
phases.  The finite checks are adversarial regeneration, not the
all-size proof recorded in the accompanying proof note.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Sequence

from mpmath import mp
import numpy as np
from scipy.linalg import solve_discrete_lyapunov
from scipy.optimize import least_squares, minimize
import sympy as sp

from crabb_circular_normal_quadratic_exact import characteristic_series
from crabb_circular_normal_series import (
    inverse_square_root_from_series,
    optimized_defect_jets,
    ordinary_triple_series_coefficient,
)
from crabb_elliptic_axis import elliptic_modulus_from_nome
from crabb_disk_one_reflection_jets import (
    cleared_resolvent_gradient_at,
)
from crabb_full_disk_equality_manifold import (
    nonlinear_equations,
    numerical_disk_model,
    pack_hermitian,
    toeplitz_coordinates,
    unpack_hermitian,
)
from crabb_full_disk_hardy_geometry import (
    formal_equality_jet,
    linearized_hardy_residual,
)
from general_crabb_weighted_series import inverse_riemann_series
from rank_one_stein_series import diagonal_gramian_condition_series


@dataclass(frozen=True)
class CharacteristicChartRecord:
    """One exact formal characteristic-coordinate audit."""

    record_kind: str
    dimension: int
    length: int
    jet_order: int
    characteristic_linearization_verified: bool
    characteristic_constant_zero_through_order: int
    faber_reflection_identity_verified: bool
    all_identities_verified: bool


@dataclass(frozen=True)
class LeadingOrthogonalityRecord:
    """One exact residual/reflection leading-face audit."""

    record_kind: str
    dimension: int
    length: int
    reflected_grade: int
    residual_real_parameter: str
    hardy_residual_norm_square: str
    optimized_face: str
    predicted_face: str
    mixed_coefficient_zero: bool
    all_identities_verified: bool


@dataclass(frozen=True)
class EndpointGradientRecord:
    """One exact full-Gram endpoint first-jet audit."""

    record_kind: str
    dimension: int
    length: int
    residual_value_zero: bool
    gradient_endpoint_support_only: bool
    endpoint_values_match: bool
    all_identities_verified: bool


@dataclass(frozen=True)
class NonlinearEllipticRecord:
    """One floating full-equality elliptic phase audit."""

    record_kind: str
    dimension: int
    length: int
    equality_residual_norm: str
    transverse_nontoeplitz_norm: str
    characteristic_grade_one: str
    predicted_quadratic_descent: str
    fitted_quadratic_descents: tuple[str, ...]
    maximum_phase_spread: str
    maximum_prediction_error: str
    all_checks_passed: bool


def format_float(value: float) -> str:
    """Format a deterministic binary64 diagnostic."""

    return f"{value:.12e}"


def characteristic_direction(length: int) -> tuple[sp.Expr, ...]:
    """Return one generic exact Hermitian-Toeplitz tangent."""

    return (
        sp.Integer(0),
        *[
            sp.Rational(offset + 1, 43 * length)
            + sp.I
            * (-1) ** offset
            * sp.Rational(offset + 2, 47 * length)
            for offset in range(1, length)
        ],
    )


def dickson_polynomial(
    degree: int,
    variable: sp.Symbol,
    ellipse: sp.Symbol,
) -> sp.Expr:
    """Return the monic Dickson/Faber polynomial of one degree."""

    if degree == 0:
        return sp.Integer(2)
    previous, current = sp.Integer(2), variable
    for _ in range(2, degree + 1):
        previous, current = (
            current,
            sp.expand(variable * current - ellipse * previous),
        )
    return current


def faber_identity(length: int) -> bool:
    """Verify the unrestricted characteristic reflection identity."""

    variable, boundary, ellipse = sp.symbols(
        "variable boundary ellipse",
        nonzero=True,
    )
    coefficients = sp.symbols(f"faber_u_1:{length}")
    prepared = dickson_polynomial(length, variable, ellipse)
    prepared += 2 * sum(
        (
            coefficients[degree - 1]
            * dickson_polynomial(degree, variable, ellipse)
            for degree in range(1, length)
        ),
        sp.Integer(0),
    )
    positive = boundary**length + 2 * sum(
        (
            coefficients[degree - 1] * boundary**degree
            for degree in range(1, length)
        ),
        sp.Integer(0),
    )
    reflected = ellipse**length * boundary ** (-length)
    reflected += 2 * sum(
        (
            coefficients[degree - 1]
            * ellipse**degree
            * boundary ** (-degree)
            for degree in range(1, length)
        ),
        sp.Integer(0),
    )
    residual = sp.together(
        prepared.subs(variable, boundary + ellipse / boundary)
        - positive
        - reflected
    )
    return sp.expand(residual * boundary**length) == 0


def characteristic_record(
    length: int,
    order: int,
) -> CharacteristicChartRecord:
    """Audit the characteristic chart on one exact equality jet."""

    direction = characteristic_direction(length)
    _, _, operator = formal_equality_jet(direction, order)
    characteristic = characteristic_series(list(operator))
    linearization = all(
        sp.simplify(
            characteristic[offset][1]
            - 2 * sp.conjugate(direction[offset])
        )
        == 0
        for offset in range(1, length)
    )
    factor_constant_zero = all(
        coefficient == 0 for coefficient in characteristic[length]
    )
    determinant_constant_zero = all(
        coefficient == 0 for coefficient in characteristic[length + 1]
    )
    faber_verified = faber_identity(length)
    verified = bool(
        linearization
        and factor_constant_zero
        and determinant_constant_zero
        and faber_verified
    )
    if not verified:
        raise RuntimeError(
            f"characteristic-chart audit failed in length {length}"
        )
    return CharacteristicChartRecord(
        record_kind="exact_characteristic_chart",
        dimension=length + 1,
        length=length,
        jet_order=order,
        characteristic_linearization_verified=linearization,
        characteristic_constant_zero_through_order=order,
        faber_reflection_identity_verified=faber_verified,
        all_identities_verified=verified,
    )


def deterministic_residual_direction(length: int) -> sp.Matrix:
    """Return a zero-Toeplitz real Hermitian transverse direction."""

    matrix = sp.zeros(length + 1)
    matrix[0, 1] = 1
    matrix[1, 0] = 1
    matrix[1, 2] = -1
    matrix[2, 1] = -1
    return matrix


def physical_ellipse_path(
    length: int,
    grade: int,
    residual_parameter: sp.Symbol,
    order: int,
    residual_direction: sp.Matrix | None = None,
) -> tuple[list[sp.Matrix], sp.Matrix]:
    """Build the exact equality/ellipse/residual physical path."""

    dimension = length + 1
    weighted_degree = grade + 1
    direction = [sp.Integer(0)] * length
    direction[length - grade] = sp.Rational(1, 10)
    hermitian, _, _ = formal_equality_jet(direction, order)
    hermitian = list(hermitian)
    if residual_direction is None:
        residual_direction = deterministic_residual_direction(length)
    if residual_direction.shape != (dimension, dimension):
        raise ValueError("the residual direction has the wrong shape")
    hermitian[weighted_degree] += (
        residual_parameter * residual_direction
    )

    shift = sp.zeros(dimension)
    for index in range(length):
        shift[index, index + 1] = 1
    coordinate = [
        sp.simplify(coefficient + shift.T * coefficient * shift)
        for coefficient in hermitian
    ]
    inverse_square_root = inverse_square_root_from_series(coordinate)
    shifted_hermitian = [
        coefficient * shift for coefficient in hermitian
    ]
    disk_operator = [
        2
        * ordinary_triple_series_coefficient(
            inverse_square_root,
            shifted_hermitian,
            inverse_square_root,
            degree,
        )
        for degree in range(order + 1)
    ]
    ellipse_path = [
        sp.simplify(
            disk_operator[degree]
            + (
                disk_operator[degree - 1].T
                if degree
                else sp.zeros(dimension)
            )
        )
        for degree in range(order + 1)
    ]
    return ellipse_path, residual_direction[:-1, :-1]


def leading_orthogonality_record(
    length: int,
    grade: int,
) -> LeadingOrthogonalityRecord:
    """Regenerate one exact residual/reflection completed face."""

    weighted_degree = grade + 1
    order = 2 * weighted_degree
    epsilon = sp.symbols("elliptic_epsilon", real=True)
    residual_parameter = sp.symbols(
        f"hardy_residual_{length}_{grade}",
        real=True,
    )
    path, residual_direction = physical_ellipse_path(
        length,
        grade,
        residual_parameter,
        order,
    )
    _, operator = inverse_riemann_series(path, order)
    defect = optimized_defect_jets(
        operator,
        epsilon,
        weighted_degree,
    )
    condition = diagonal_gramian_condition_series(
        operator,
        defect,
        epsilon,
        order,
    )
    face = sp.factor(sp.expand(condition).coeff(epsilon, order))
    hardy_residual = linearized_hardy_residual(residual_direction)
    residual_norm_square = sp.simplify(
        sum(entry * sp.conjugate(entry) for entry in hardy_residual)
    )
    predicted = sp.factor(
        -sp.Rational(16, 25)
        - 4 * residual_norm_square * residual_parameter**2
    )
    mixed_zero = sp.diff(face, residual_parameter).subs(
        residual_parameter,
        0,
    ) == 0
    verified = bool(sp.expand(face - predicted) == 0 and mixed_zero)
    if not verified:
        raise RuntimeError(
            "the residual/reflection leading face failed for "
            f"length {length}, grade {grade}"
        )
    return LeadingOrthogonalityRecord(
        record_kind="exact_leading_orthogonality",
        dimension=length + 1,
        length=length,
        reflected_grade=grade,
        residual_real_parameter=str(residual_parameter),
        hardy_residual_norm_square=str(residual_norm_square),
        optimized_face=str(face),
        predicted_face=str(predicted),
        mixed_coefficient_zero=mixed_zero,
        all_identities_verified=verified,
    )


def endpoint_gradient_record(length: int) -> EndpointGradientRecord:
    """Verify the unrestricted Crabb endpoint-residual gradient."""

    dimension = length + 1
    hermitian = sp.zeros(dimension)
    for index in range(length):
        hermitian[index, index] = sp.Rational(1, 2)
    variable = sp.symbols(
        f"endpoint_variable_{length}",
        nonzero=True,
    )
    residual, transfer, gradient = cleared_resolvent_gradient_at(
        hermitian,
        variable,
    )
    expected = sp.zeros(length)
    expected[0, 0] = transfer
    expected[-1, -1] = transfer
    gradient_block = gradient[:length, :length]
    support_only = all(
        sp.simplify(entry) == 0
        for entry in gradient_block - expected
    )
    endpoint_match = sp.simplify(
        transfer - variable ** (-(length + 1))
    ) == 0
    verified = bool(
        sp.simplify(residual) == 0
        and support_only
        and endpoint_match
    )
    if not verified:
        raise RuntimeError(
            f"the endpoint-gradient audit failed in length {length}"
        )
    return EndpointGradientRecord(
        record_kind="exact_endpoint_full_gram_gradient",
        dimension=dimension,
        length=length,
        residual_value_zero=True,
        gradient_endpoint_support_only=support_only,
        endpoint_values_match=endpoint_match,
        all_identities_verified=verified,
    )


def solve_nonlinear_equality(length: int) -> np.ndarray:
    """Return one genuinely non-Toeplitz numerical equality Gram block."""

    initial = np.eye(length, dtype=complex) / 2
    for offset in range(1, length):
        value = (
            0.018 * (offset + 1) / length
            + 1j * 0.011 * (length - offset) / length
        )
        for row in range(length - offset):
            initial[row, row + offset] = value
            initial[row + offset, row] = np.conjugate(value)
    target = toeplitz_coordinates(initial)
    solution = least_squares(
        nonlinear_equations,
        pack_hermitian(initial),
        args=(length, target),
        xtol=1e-13,
        ftol=1e-13,
        gtol=1e-13,
        max_nfev=300,
    )
    if not solution.success:
        raise RuntimeError("the nonlinear equality solve failed")
    return unpack_hermitian(solution.x, length)


def physical_disk_operator(hermitian: np.ndarray) -> np.ndarray:
    """Transform one coefficient-gauge disk matrix to Euclidean gauge."""

    _, coordinate, operator, _ = numerical_disk_model(hermitian)
    eigenvalues, eigenvectors = np.linalg.eigh(coordinate)
    square_root = (
        eigenvectors * np.sqrt(eigenvalues)
    ) @ eigenvectors.conj().T
    inverse_square_root = (
        eigenvectors * (1 / np.sqrt(eigenvalues))
    ) @ eigenvectors.conj().T
    return square_root @ operator @ inverse_square_root


def ellipse_pullback(
    disk_operator: np.ndarray,
    ellipse_parameter: float,
    phase: float,
) -> np.ndarray:
    """Apply the exact elliptic scalar map in one rotated phase."""

    rotated = np.exp(-0.5j * phase) * disk_operator
    ellipse_operator = (
        rotated + ellipse_parameter * rotated.conj().T
    )
    modulus = elliptic_modulus_from_nome(ellipse_parameter**2)
    elliptic_parameter = modulus**2
    quarter_period = mp.ellipk(elliptic_parameter)
    sn = mp.ellipfun("sn")
    eigenvalues, eigenvectors = np.linalg.eig(ellipse_operator)
    mapped = []
    for eigenvalue in eigenvalues:
        argument = (
            (2 * quarter_period / mp.pi)
            * mp.asin(
                complex(eigenvalue)
                / (2 * np.sqrt(ellipse_parameter))
            )
        )
        mapped.append(
            complex(
                mp.sqrt(modulus)
                * sn(argument, elliptic_parameter)
            )
        )
    return (
        eigenvectors
        @ np.diag(mapped)
        @ np.linalg.inv(eigenvectors)
    )


def rank_one_envelope(
    operator: np.ndarray,
    initial: np.ndarray | None = None,
) -> tuple[float, np.ndarray]:
    """Optimize the complex rank-one Stein branch."""

    dimension = len(operator)
    if initial is None:
        initial = np.zeros(2 * (dimension - 1))

    def objective(values: np.ndarray) -> float:
        tail = (
            values[: dimension - 1]
            + 1j * values[dimension - 1 :]
        )
        defect = np.concatenate(([1.0], tail))
        gramian = solve_discrete_lyapunov(
            operator.conj().T,
            np.outer(defect, defect.conj()),
        )
        eigenvalues = np.linalg.eigvalsh(gramian)
        return float(eigenvalues[-1] / eigenvalues[0])

    result = minimize(
        objective,
        initial,
        method="BFGS",
        options={"gtol": 1e-10, "maxiter": 2_000},
    )
    if not np.isfinite(result.fun):
        raise RuntimeError("the complex rank-one envelope failed")
    return float(result.fun), np.asarray(result.x)


def toeplitz_projection(hermitian: np.ndarray) -> np.ndarray:
    """Return the Hermitian-Toeplitz projection of one Gram block."""

    length = len(hermitian)
    projection = np.zeros_like(hermitian)
    for offset in range(length):
        mean = np.mean(np.diag(hermitian, offset))
        for row in range(length - offset):
            projection[row, row + offset] = mean
            projection[row + offset, row] = np.conjugate(mean)
    return projection


def nonlinear_elliptic_record(length: int) -> NonlinearEllipticRecord:
    """Test the elliptic descent on one nonlinear equality anchor."""

    hermitian = solve_nonlinear_equality(length)
    residual_norm = float(
        np.linalg.norm(
            nonlinear_equations(
                pack_hermitian(hermitian),
                length,
                toeplitz_coordinates(hermitian),
            )[: (length - 1) ** 2]
        )
    )
    nontoeplitz = float(
        np.linalg.norm(hermitian - toeplitz_projection(hermitian))
    )
    disk_operator = physical_disk_operator(hermitian)
    characteristic_factor = np.poly(disk_operator)[:-1]
    grade_one = characteristic_factor[-2] / 2
    predicted = 64 * abs(grade_one) ** 2

    parameters = np.asarray((0.004, 0.008, 0.016))
    phases = (0.0, np.pi / 2, np.pi)
    fitted = []
    for phase in phases:
        values = []
        initial = None
        for parameter in parameters:
            operator = ellipse_pullback(
                disk_operator,
                float(parameter),
                phase,
            )
            value, initial = rank_one_envelope(operator, initial)
            values.append(value - 4)
        ratios = -np.asarray(values) / parameters**2
        fitted.append(
            float(np.polyval(np.polyfit(parameters, ratios, 1), 0))
        )
    phase_spread = max(fitted) - min(fitted)
    prediction_error = max(abs(value - predicted) for value in fitted)
    tolerance = max(2e-5, 0.01 * predicted)
    verified = bool(
        residual_norm < 1e-10
        and nontoeplitz > 1e-7
        and phase_spread < tolerance
        and prediction_error < tolerance
    )
    if not verified:
        raise RuntimeError(
            f"the nonlinear elliptic audit failed in length {length}"
        )
    return NonlinearEllipticRecord(
        record_kind="nonlinear_equality_elliptic_phase",
        dimension=length + 1,
        length=length,
        equality_residual_norm=format_float(residual_norm),
        transverse_nontoeplitz_norm=format_float(nontoeplitz),
        characteristic_grade_one=(
            f"{grade_one.real:.12e}{grade_one.imag:+.12e}j"
        ),
        predicted_quadratic_descent=format_float(predicted),
        fitted_quadratic_descents=tuple(
            format_float(value) for value in fitted
        ),
        maximum_phase_spread=format_float(phase_spread),
        maximum_prediction_error=format_float(prediction_error),
        all_checks_passed=verified,
    )


def write_records(
    path: Path,
    records: Sequence[
        CharacteristicChartRecord
        | EndpointGradientRecord
        | LeadingOrthogonalityRecord
        | NonlinearEllipticRecord
    ],
) -> None:
    """Write deterministic JSON Lines output atomically."""

    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(asdict(record), sort_keys=True))
            handle.write("\n")
    temporary.replace(path)


def main() -> None:
    """Run exact and floating full-equality elliptic audits."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--minimum-length", type=int, default=3)
    parser.add_argument("--maximum-length", type=int, default=6)
    parser.add_argument("--formal-order", type=int, default=4)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "crabb_full_equality_elliptic_merger_s70224.jsonl"
        ),
    )
    args = parser.parse_args()

    records: list[
        CharacteristicChartRecord
        | EndpointGradientRecord
        | LeadingOrthogonalityRecord
        | NonlinearEllipticRecord
    ] = []
    for length in range(args.minimum_length, args.maximum_length + 1):
        record = characteristic_record(length, args.formal_order)
        records.append(record)
        print(f"verified characteristic chart length {length}", flush=True)
        endpoint_record = endpoint_gradient_record(length)
        records.append(endpoint_record)
        print(f"verified endpoint gradient length {length}", flush=True)
    for length in range(args.minimum_length, args.maximum_length + 1):
        record = leading_orthogonality_record(length, grade=1)
        records.append(record)
        print(f"verified leading grade one length {length}", flush=True)
    if args.minimum_length <= 4 <= args.maximum_length:
        record = leading_orthogonality_record(4, grade=2)
        records.append(record)
        print("verified leading grade two length 4", flush=True)
    for length in range(args.minimum_length, args.maximum_length + 1):
        record = nonlinear_elliptic_record(length)
        records.append(record)
        print(f"verified nonlinear elliptic length {length}", flush=True)
    write_records(args.output, records)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
