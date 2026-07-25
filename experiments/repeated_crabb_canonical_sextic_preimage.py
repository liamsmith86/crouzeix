#!/usr/bin/env python3
"""Audit the polynomial preparation of the canonical sextic flag.

After the cubic, quartic, and quintic defect-factor columns have been
inserted, this checker verifies an explicit bounded sixth column whose
direct face on

    ker(B_1*) intersect ker(B_2*)

is exactly ``12 B_3 B_3*``.  It then includes the fifth-order
off-diagonal Schur cost and proves a uniform positive remaining
coefficient.  The large discovered certificate is stored separately
as exact rational JSON and is fully rechecked here.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_boundary_metric_flag import (
    lossless_schur_realization,
    rank_chain_case,
)
from repeated_crabb_canonical_quintic_preimage import (
    cubic_coboundary_witness,
    cubic_column_lift,
    evaluate_polynomial,
    exact_canonical_factor_lifts,
    prepared_endpoint_series,
    quintic_coboundary_witness,
    quintic_column_lift,
)
from repeated_crabb_canonical_quartic_preimage import (
    quartic_coboundary_witness,
    quartic_column_lift,
)
from repeated_crabb_delayed_slack_anticommutator import (
    E,
    F,
    IDENTITY,
    Polynomial,
    S,
    STAR,
    add,
    adjoint,
    boundary_metric_series,
    multiply,
    operator_series,
    scale,
    schur_residual_series,
)
from repeated_crabb_elliptic_cokernel import stein_inverse
from repeated_crabb_transfer_channel_covariance import (
    inflated_case,
    random_partial_isometry,
    transfer_coefficient,
)


Matrix = np.ndarray
DEGREE = 6
CERTIFICATE_PATH = Path(__file__).with_name(
    "repeated_crabb_canonical_sextic_certificate.json"
)
UNIFORM_MIDDLE_MARGIN = Fraction(215, 22)


@dataclass(frozen=True)
class CanonicalSexticPreimageRecord:
    """One exact-structure and numerical sextic preparation audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    parameter_scale: str
    colligation_error: str
    spectral_radius: str
    first_transfer_norm: str
    second_transfer_norm: str
    third_transfer_norm: str
    sextic_endpoint_norm: str
    endpoint_reconstruction_error: str
    left_flag_dimension: int
    direct_flag_face_error: str
    effective_sextic_minimum_eigenvalue: str
    analytic_lower_bound_surplus_minimum_eigenvalue: str
    analytic_middle_weight_minimum_eigenvalue: str
    polynomial_column_norm: str
    perpendicular_column_error: str
    lower_endpoint_motion_error: str
    canonical_factor_error: str
    exact_coboundary_residual_word_count: int
    exact_flag_factor_residual_word_count: int
    exact_witness_hermitian_residual_word_count: int
    exact_column_perpendicular_residual_word_count: int
    exact_quintic_channel_residual_word_count: int
    exact_right_gram_residual_word_count: int
    all_checks_passed: bool


def decode_fraction(encoded: list[int]) -> Fraction:
    """Decode one ``[numerator, denominator]`` pair."""

    return Fraction(encoded[0], encoded[1])


@lru_cache(maxsize=1)
def load_certificate() -> dict[str, object]:
    """Load the tracked exact rational sextic certificate."""

    certificate = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
    if certificate.get("format_version") != 1:
        raise RuntimeError("unsupported sextic certificate format")
    expected_counts = {
        "column_bracket": 53,
        "witness": 139,
        "flag_multipliers": 28,
    }
    for key, expected in expected_counts.items():
        value = certificate.get(key)
        if not isinstance(value, list) or len(value) != expected:
            raise RuntimeError(
                f"invalid sextic certificate section {key!r}"
            )
    return certificate


def sextic_column_bracket() -> Polynomial:
    """Return the 53-term polynomial inside the left projection."""

    result: Polynomial = {}
    entries = load_certificate()["column_bracket"]
    if not isinstance(entries, list):
        raise RuntimeError("invalid sextic column certificate")
    for numerator, denominator, word in entries:
        result = add(
            result,
            scale(
                Fraction(numerator, denominator),
                {word: 1},
            ),
        )
    return result


def sextic_column_lift() -> Polynomial:
    """Return the exact perpendicular state lift ``C_6 V*``."""

    initial_complement = add(IDENTITY, scale(-1, E))
    return multiply(
        multiply(
            initial_complement,
            sextic_column_bracket(),
        ),
        E,
    )


def sextic_column(
    partial: Matrix,
    right: Matrix,
) -> Matrix:
    """Evaluate the explicit sixth defect-factor column."""

    return evaluate_polynomial(sextic_column_lift(), partial) @ right


def sextic_coboundary_witness() -> Polynomial:
    """Return the exact Hermitian sixth Stein witness."""

    result: Polynomial = {}
    entries = load_certificate()["witness"]
    if not isinstance(entries, list):
        raise RuntimeError("invalid sextic witness certificate")
    for numerator, denominator, encoded_polynomial in entries:
        polynomial = {
            word: decode_fraction(encoded)
            for word, encoded in encoded_polynomial.items()
        }
        result = add(
            result,
            scale(
                Fraction(numerator, denominator),
                polynomial,
            ),
        )
    return result


def sextic_flag_multipliers() -> tuple[Polynomial, Polynomial]:
    """Return the two exact right-ideal endpoint multipliers."""

    multipliers: list[Polynomial] = [{}, {}]
    entries = load_certificate()["flag_multipliers"]
    if not isinstance(entries, list):
        raise RuntimeError("invalid sextic flag certificate")
    for numerator, denominator, word, channel in entries:
        multipliers[channel - 1] = add(
            multipliers[channel - 1],
            scale(
                Fraction(numerator, denominator),
                {word: 1},
            ),
        )
    return multipliers[0], multipliers[1]


def prepared_exact_metrics() -> list[Polynomial]:
    """Return completed metric coefficients through degree five."""

    boundary = boundary_metric_series(DEGREE)
    residual = schur_residual_series(DEGREE)
    metric = [{} for _ in range(DEGREE + 1)]
    metric[0] = IDENTITY
    metric[2] = add(
        boundary[2],
        scale(-1, residual[2]),
    )
    metric[3] = cubic_coboundary_witness()
    metric[4] = quartic_coboundary_witness()
    metric[5] = quintic_coboundary_witness()
    return metric


@lru_cache(maxsize=1)
def exact_sixth_components() -> tuple[Polynomial, Polynomial]:
    """Return the sixth Stein forcing and upper Schur square."""

    factors = exact_canonical_factor_lifts(DEGREE)
    factors[3] = add(factors[3], cubic_column_lift())
    factors[4] = add(factors[4], quartic_column_lift())
    factors[5] = add(factors[5], quintic_column_lift())
    prepared_slack: list[Polynomial] = []
    for order in range(DEGREE + 1):
        coefficient: Polynomial = {}
        for left_degree in range(order + 1):
            coefficient = add(
                coefficient,
                multiply(
                    factors[left_degree],
                    adjoint(factors[order - left_degree]),
                ),
            )
        prepared_slack.append(coefficient)

    operator = operator_series(DEGREE)
    metric = prepared_exact_metrics()
    forcing = prepared_slack[DEGREE]
    for left_degree in range(DEGREE + 1):
        for metric_degree in range(DEGREE):
            right_degree = DEGREE - left_degree - metric_degree
            if 0 <= right_degree <= DEGREE:
                forcing = add(
                    forcing,
                    multiply(
                        multiply(
                            adjoint(operator[left_degree]),
                            metric[metric_degree],
                        ),
                        operator[right_degree],
                    ),
                )

    inverse = add(
        IDENTITY,
        scale(-1, F),
        scale(Fraction(-2, 3), E),
    )
    schur_square: Polynomial = {}
    for left_degree, right_degree in ((2, 4), (4, 2), (3, 3)):
        schur_square = add(
            schur_square,
            multiply(
                multiply(
                    multiply(
                        metric[left_degree],
                        inverse,
                    ),
                    metric[right_degree],
                ),
                F,
            ),
        )
    schur_square = add(
        schur_square,
        multiply(
            multiply(
                multiply(
                    multiply(
                        multiply(metric[2], inverse),
                        metric[2],
                    ),
                    inverse,
                ),
                metric[2],
            ),
            F,
        ),
    )
    return forcing, schur_square


def channel_lift(order: int) -> Polynomial:
    """Return ``E S^order F``."""

    return multiply(
        multiply(E, {"s" * order: 1}),
        F,
    )


def grade_three_gram_lift() -> Polynomial:
    """Return the state lift of ``B_3 B_3*``."""

    third = channel_lift(3)
    return multiply(adjoint(third), third)


def sextic_endpoint_lift() -> Polynomial:
    """Return one minus-quarter of the physical sixth endpoint."""

    _, schur_square = exact_sixth_components()
    witness = sextic_coboundary_witness()
    return add(
        multiply(multiply(F, witness), F),
        multiply(F, schur_square),
    )


def quintic_flag_multipliers() -> tuple[Polynomial, Polynomial]:
    """Return L233's two fifth-endpoint multipliers."""

    first = multiply(
        F,
        add(
            scale(Fraction(7, 2), {"aaaasssaa": 1}),
            scale(Fraction(5, 2), {"aasssaaaa": 1}),
        ),
    )
    second = multiply(F, {"aassssaa": 1})
    return first, second


def quintic_third_channel_multiplier() -> Polynomial:
    """Return the multiplier controlling L233's flag cross term."""

    return scale(
        Fraction(7, 2),
        add(
            IDENTITY,
            scale(-1, {"ssaa": 1}),
            {"asssaa": 1},
        ),
    )


@lru_cache(maxsize=1)
def exact_residuals() -> tuple[int, int, int, int, int, int]:
    """Return word counts in every load-bearing exact identity."""

    forcing, _ = exact_sixth_components()
    column = sextic_column_lift()
    witness = sextic_coboundary_witness()
    coboundary = add(
        witness,
        scale(
            -1,
            multiply(multiply(STAR, witness), S),
        ),
    )
    forcing_residual = add(
        forcing,
        column,
        adjoint(column),
        scale(-1, coboundary),
    )

    first_flag, second_flag = sextic_flag_multipliers()
    first_channel = channel_lift(1)
    second_channel = channel_lift(2)
    flag_remainder = add(
        multiply(first_flag, first_channel),
        adjoint(multiply(first_flag, first_channel)),
        multiply(second_flag, second_channel),
        adjoint(multiply(second_flag, second_channel)),
    )
    flag_residual = add(
        sextic_endpoint_lift(),
        scale(3, grade_three_gram_lift()),
        scale(-1, flag_remainder),
    )

    first_quintic, second_quintic = quintic_flag_multipliers()
    third_multiplier = quintic_third_channel_multiplier()
    first_earlier_multiplier = add(
        scale(Fraction(-5, 2), {"ssssaa": 1}),
        scale(Fraction(5, 2), {"asssssaa": 1}),
    )
    second_earlier_multiplier = add(
        scale(-1, {"ssaaa": 1}),
        {"asssaaa": 1},
    )
    first_channel_residual = add(
        multiply(multiply(E, adjoint(first_quintic)), F),
        scale(
            -1,
            add(
                multiply(third_multiplier, channel_lift(3)),
                multiply(first_earlier_multiplier, first_channel),
            ),
        ),
    )
    second_channel_residual = add(
        multiply(multiply(E, adjoint(second_quintic)), F),
        scale(
            -1,
            multiply(second_earlier_multiplier, first_channel),
        ),
    )
    first_right_gram = multiply(
        multiply(
            multiply(
                multiply(E, S),
                F,
            ),
            STAR,
        ),
        E,
    )
    right_gram_residual = add(
        multiply(
            multiply(E, third_multiplier),
            E,
        ),
        scale(Fraction(-7, 2), first_right_gram),
    )
    return (
        len(forcing_residual),
        len(flag_residual),
        len(add(witness, scale(-1, adjoint(witness)))),
        len(multiply(E, column)),
        len(first_channel_residual) + len(second_channel_residual),
        len(right_gram_residual),
    )


def effective_sextic_face(
    upper: list[Matrix],
) -> tuple[Matrix, Matrix]:
    """Return the quartic-kernel sixth Schur face and its frame."""

    quartic = (upper[4] + upper[4].conj().T) / 2
    eigenvalues, eigenvectors = np.linalg.eigh(quartic)
    tolerance = max(1e-9, 1e-8 * eigenvalues[-1])
    active = eigenvectors[:, eigenvalues > tolerance]
    kernel = eigenvectors[:, eigenvalues <= tolerance]
    if kernel.shape[1] == 0:
        return np.zeros((0, 0), dtype=complex), kernel
    direct = kernel.conj().T @ upper[6] @ kernel
    if active.shape[1]:
        leading = active.conj().T @ quartic @ active
        cross = active.conj().T @ upper[5] @ kernel
        direct -= cross.conj().T @ np.linalg.solve(leading, cross)
    return (direct + direct.conj().T) / 2, kernel


def analytic_middle_weight(first: Matrix) -> Matrix:
    """Return the uniform positive middle weight after the Schur cost."""

    right_gram = first.conj().T @ first
    identity = np.eye(len(right_gram), dtype=complex)
    correction = (
        float(Fraction(49, 2))
        * right_gram
        @ np.linalg.solve(
            4 * identity + 7 * right_gram,
            right_gram,
        )
    )
    return 12 * identity - correction


def deterministic_unitary(
    generator: np.random.Generator,
    dimension: int,
) -> Matrix:
    """Return a deterministic Haar-style unitary."""

    matrix = (
        generator.standard_normal((dimension, dimension))
        + 1j * generator.standard_normal((dimension, dimension))
    )
    unitary, triangular = np.linalg.qr(matrix)
    diagonal = np.diag(triangular)
    phases = np.where(
        np.abs(diagonal) > 0,
        diagonal / np.abs(diagonal),
        1,
    )
    return unitary @ np.diag(np.conjugate(phases))


def low_rank_schur_parameter(
    generator: np.random.Generator,
    dimension: int,
    rank: int,
    norm: float,
) -> Matrix:
    """Return one deterministic Schur parameter of prescribed rank."""

    left = deterministic_unitary(generator, dimension)[:, :rank]
    right = deterministic_unitary(generator, dimension)[:rank, :]
    return norm * left @ right


def high_amplitude_flag_case(
    repetition: int,
) -> tuple[Matrix, Matrix, Matrix, float]:
    """Return a valid rank-two flag near the rejected-candidate edge."""

    dimension = 5
    generator = np.random.default_rng(180_500 + repetition)
    ranks = (2, 1, 2)
    scale_value = 0.8
    parameters = [np.zeros((dimension, dimension), dtype=complex)]
    parameters.extend(
        low_rank_schur_parameter(
            generator,
            dimension,
            rank,
            scale_value * (1 + 0.1 * index),
        )
        for index, rank in enumerate(ranks)
    )
    terminal = deterministic_unitary(generator, dimension)
    return lossless_schur_realization(parameters, terminal)


def audit_case(
    construction_kind: str,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    parameter_scale: float,
) -> CanonicalSexticPreimageRecord:
    """Audit one exact sextic polynomial preparation."""

    column = sextic_column(partial, right)
    upper, _, factor_error = prepared_endpoint_series(
        partial,
        right,
        left,
        DEGREE,
        additional_factor_columns={DEGREE: column},
    )
    first = transfer_coefficient(partial, right, left, 1)
    second = transfer_coefficient(partial, right, left, 2)
    third = transfer_coefficient(partial, right, left, 3)

    endpoint_lift = evaluate_polynomial(
        sextic_endpoint_lift(),
        partial,
    )
    predicted_endpoint = (
        -4 * left.conj().T @ endpoint_lift @ left
    )
    endpoint_error = float(
        np.linalg.norm(upper[6] - predicted_endpoint)
    )

    effective, flag = effective_sextic_face(upper)
    direct = flag.conj().T @ upper[6] @ flag
    expected_direct = (
        12
        * flag.conj().T
        @ third
        @ third.conj().T
        @ flag
    )
    direct_error = float(np.linalg.norm(direct - expected_direct))
    middle_weight = analytic_middle_weight(first)
    analytic_bound = (
        flag.conj().T
        @ third
        @ middle_weight
        @ third.conj().T
        @ flag
    )
    surplus = (effective - analytic_bound)
    surplus = (surplus + surplus.conj().T) / 2
    effective_minimum = (
        float(np.linalg.eigvalsh(effective)[0])
        if len(effective)
        else 0.0
    )
    surplus_minimum = (
        float(np.linalg.eigvalsh(surplus)[0])
        if len(surplus)
        else 0.0
    )
    middle_minimum = float(
        np.linalg.eigvalsh(middle_weight)[0]
    )

    perpendicular_error = float(
        np.linalg.norm(right.conj().T @ column)
    )
    column_forcing = (
        right @ column.conj().T
        + column @ right.conj().T
    )
    column_response = stein_inverse(partial, column_forcing)
    lower_motion_error = float(
        np.linalg.norm(
            right.conj().T @ column_response @ right
        )
    )

    identity = np.eye(len(partial), dtype=complex)
    colligation_error = max(
        float(
            np.linalg.norm(
                partial.conj().T @ partial
                - (identity - right @ right.conj().T)
            )
        ),
        float(
            np.linalg.norm(
                partial @ partial.conj().T
                - (identity - left @ left.conj().T)
            )
        ),
        float(np.linalg.norm(right.conj().T @ left)),
    )
    spectral_radius = float(
        np.max(np.abs(np.linalg.eigvals(partial)))
    )
    exact = exact_residuals()

    tolerance = 4e-8
    verified = bool(
        colligation_error < tolerance
        and spectral_radius < 1
        and endpoint_error < tolerance
        and direct_error < tolerance
        and effective_minimum > -tolerance
        and surplus_minimum > -tolerance
        and middle_minimum
        > float(UNIFORM_MIDDLE_MARGIN) - tolerance
        and np.linalg.norm(column) <= 181 + tolerance
        and perpendicular_error < tolerance
        and lower_motion_error < tolerance
        and factor_error < tolerance
        and all(value == 0 for value in exact)
    )
    if not verified:
        raise RuntimeError(
            "the canonical sextic preimage audit failed: "
            f"kind={construction_kind}, "
            f"endpoint={endpoint_error:.3e}, "
            f"direct={direct_error:.3e}, "
            f"effective={effective_minimum:.3e}, "
            f"surplus={surplus_minimum:.3e}, "
            f"middle={middle_minimum:.3e}, "
            f"factor={factor_error:.3e}, exact={exact}"
        )

    return CanonicalSexticPreimageRecord(
        construction_kind=construction_kind,
        state_dimension=len(partial),
        defect_dimension=right.shape[1],
        parameter_scale=format_float(parameter_scale),
        colligation_error=format_float(colligation_error),
        spectral_radius=format_float(spectral_radius),
        first_transfer_norm=format_float(float(np.linalg.norm(first))),
        second_transfer_norm=format_float(float(np.linalg.norm(second))),
        third_transfer_norm=format_float(float(np.linalg.norm(third))),
        sextic_endpoint_norm=format_float(
            float(np.linalg.norm(upper[6]))
        ),
        endpoint_reconstruction_error=format_float(endpoint_error),
        left_flag_dimension=flag.shape[1],
        direct_flag_face_error=format_float(direct_error),
        effective_sextic_minimum_eigenvalue=format_float(
            effective_minimum
        ),
        analytic_lower_bound_surplus_minimum_eigenvalue=format_float(
            surplus_minimum
        ),
        analytic_middle_weight_minimum_eigenvalue=format_float(
            middle_minimum
        ),
        polynomial_column_norm=format_float(float(np.linalg.norm(column))),
        perpendicular_column_error=format_float(perpendicular_error),
        lower_endpoint_motion_error=format_float(lower_motion_error),
        canonical_factor_error=format_float(factor_error),
        exact_coboundary_residual_word_count=exact[0],
        exact_flag_factor_residual_word_count=exact[1],
        exact_witness_hermitian_residual_word_count=exact[2],
        exact_column_perpendicular_residual_word_count=exact[3],
        exact_quintic_channel_residual_word_count=exact[4],
        exact_right_gram_residual_word_count=exact[5],
        all_checks_passed=verified,
    )


def standard_records() -> list[CanonicalSexticPreimageRecord]:
    """Return unstructured, rank-changing, delayed, and edge records."""

    records: list[CanonicalSexticPreimageRecord] = []
    for defect_dimension in range(1, 4):
        for repetition in range(2):
            state_dimension = 5 * defect_dimension + 3 + repetition
            partial, right, left = random_partial_isometry(
                state_dimension,
                defect_dimension,
                np.random.default_rng(
                    118_000 + 100 * defect_dimension + repetition
                ),
            )
            records.append(
                audit_case(
                    "unstructured",
                    partial,
                    right,
                    left,
                    1,
                )
            )

    for defect_dimension in range(3, 7):
        for parameter_scale in (1.0, 0.5, 0.2, 0.05):
            partial, right, left, _, _ = rank_chain_case(
                defect_dimension,
                109_200 + defect_dimension,
                parameter_scale,
            )
            records.append(
                audit_case(
                    "rank_chain",
                    partial,
                    right,
                    left,
                    parameter_scale,
                )
            )

    for defect_dimension in (2, 3, 4):
        partial, right, left, _ = inflated_case(
            4 * defect_dimension + 3,
            defect_dimension,
            3,
            defect_dimension,
            118_800 + defect_dimension,
        )
        records.append(
            audit_case(
                "complete_double_delay",
                partial,
                right,
                left,
                1,
            )
        )

    for repetition in range(3):
        partial, right, left, _ = high_amplitude_flag_case(repetition)
        records.append(
            audit_case(
                "high_amplitude_rank_two_flag",
                partial,
                right,
                left,
                0.8,
            )
        )
    return records


def write_records(
    records: list[CanonicalSexticPreimageRecord],
    output: Path,
) -> str:
    """Write deterministic JSON Lines atomically and return its hash."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)
    return hashlib.sha256(output.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_canonical_sextic_preimage_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the complete audit."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(
        json.dumps(
            {
                "certificate_sha256": hashlib.sha256(
                    CERTIFICATE_PATH.read_bytes()
                ).hexdigest(),
                "dataset_sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
