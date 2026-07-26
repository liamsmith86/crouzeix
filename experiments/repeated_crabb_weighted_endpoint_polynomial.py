#!/usr/bin/env python3
"""Audit weighted transfer factorization of causal endpoint polynomials.

Put ``J=(I+F)S*(I+E)``, the physical reverse edge.  If a word is
multiplied by ``c**d`` and contains at most ``d`` copies of ``J``, a
wandering-chain recursion factors its endpoint without losing a power
of ``c``:

    W* P(c,S,J) V = [c B_1, ..., c**R B_R] C(c).

This checker tests random causal word polynomials and the exact finite
jets of the direct elliptic map ``phi_c(S+cJ)``.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from crabb_palindromic_elliptic_hessian import direct_map_coefficients
from repeated_crabb_boundary_metric_flag import rank_chain_case
from repeated_crabb_delayed_jet import ellipse_operator_coefficients
from repeated_crabb_endpoint_word_gram import transfer_coefficient
from repeated_crabb_endpoint_null_gauge import (
    delayed_random_partial_isometry,
)
from repeated_crabb_transfer_channel_covariance import (
    random_partial_isometry,
)


Matrix = np.ndarray
Term = tuple[complex, int, str]
MAXIMUM_DEGREE = 6


@dataclass(frozen=True)
class WeightedEndpointPolynomialRecord:
    """One weighted endpoint-factor audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    maximum_degree: int
    term_count: int
    ellipse_parameter: str
    maximum_reverse_excess: int
    factor_norm: str
    theoretical_factor_bound: str
    factor_bound_ratio: str
    reconstruction_error: str
    minimum_gram_slack_eigenvalue: str
    individual_words_checked: int
    maximum_individual_reconstruction_error: str
    maximum_individual_bound_ratio: str
    direct_operator_jet_error: str
    all_checks_passed: bool


def causal_random_terms(
    maximum_degree: int,
    generator: np.random.Generator,
) -> list[Term]:
    """Return a deterministic collection of causal word terms."""

    terms: list[Term] = []
    for degree in range(1, maximum_degree + 1):
        for _ in range(5):
            length = int(generator.integers(1, 2 * maximum_degree + 2))
            maximum_reverses = min(degree, length)
            reverse_count = int(
                generator.integers(0, maximum_reverses + 1)
            )
            letters = ["j"] * reverse_count + ["s"] * (
                length - reverse_count
            )
            generator.shuffle(letters)
            coefficient = (
                generator.standard_normal()
                + 1j * generator.standard_normal()
            ) / (degree + length)
            terms.append((complex(coefficient), degree, "".join(letters)))
    return terms


def direct_map_terms(maximum_degree: int) -> list[Term]:
    """Expand the finite direct-map jet as causal word terms."""

    terms: list[Term] = []
    coefficients = direct_map_coefficients(
        maximum_degree,
        maximum_degree + 1,
    )
    for power_index, scalar_series in enumerate(coefficients):
        word_length = 2 * power_index + 1
        for scalar_degree in range(maximum_degree + 1):
            scalar = scalar_series.coefficient(scalar_degree)
            if scalar == 0:
                continue
            for reverse_count in range(word_length + 1):
                total_degree = scalar_degree + reverse_count
                if not 1 <= total_degree <= maximum_degree:
                    continue
                for reverse_positions in itertools.combinations(
                    range(word_length),
                    reverse_count,
                ):
                    positions = set(reverse_positions)
                    word = "".join(
                        "j" if index in positions else "s"
                        for index in range(word_length)
                    )
                    terms.append(
                        (
                            complex(Fraction(scalar)),
                            total_degree,
                            word,
                        )
                    )
    return terms


def physical_reverse(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
) -> Matrix:
    """Return ``J=(I+F)S*(I+E)``."""

    identity = np.eye(len(operator), dtype=complex)
    initial_defect = right @ right.conj().T
    final_defect = left @ left.conj().T
    return (
        (identity + final_defect)
        @ operator.conj().T
        @ (identity + initial_defect)
    )


def physical_word_value(
    word: str,
    operator: Matrix,
    reverse: Matrix,
) -> Matrix:
    """Evaluate a word in the physical letters ``S,J``."""

    result = np.eye(len(operator), dtype=complex)
    for letter in word:
        result = result @ (operator if letter == "s" else reverse)
    return result


def physical_word_factor_coefficients(
    word: str,
    operator: Matrix,
    right: Matrix,
    left: Matrix,
) -> tuple[Matrix, ...]:
    """Construct ``W*w(S,J)V=sum B_j C_j`` by chain propagation."""

    reverse_count = word.count("j")
    multiplicity = right.shape[1]
    coefficients_adjoint = [
        np.zeros((multiplicity, multiplicity), dtype=complex)
        for _ in range(reverse_count)
    ]
    if reverse_count == 0:
        return tuple(coefficients_adjoint)

    reverse = physical_reverse(operator, right, left)
    adjoint_letters = [
        operator.conj().T if letter == "s" else reverse.conj().T
        for letter in word
    ]

    height = 0
    principal_coefficient = 1.0
    principal_alive = True
    for index, letter in enumerate(word):
        if not principal_alive:
            break

        error_index = 0
        error_coefficient = 0.0
        if letter == "s":
            if height == 0:
                principal_alive = False
                continue
            error_index = height - 1
            error_coefficient = -principal_coefficient
            height -= 1
        elif height == 0:
            error_index = 1
            error_coefficient = 2 * principal_coefficient
            principal_coefficient *= 2
            height = 1
        else:
            error_index = height + 1
            error_coefficient = principal_coefficient
            height += 1

        if error_index == 0:
            continue
        remaining = np.eye(len(operator), dtype=complex)
        for adjoint_letter in adjoint_letters[index + 1 :]:
            remaining = adjoint_letter @ remaining
        coefficients_adjoint[error_index - 1] += (
            error_coefficient
            * right.conj().T
            @ remaining
            @ right
        )

    if principal_alive and height > 0:
        coefficients_adjoint[height - 1] += (
            principal_coefficient
            * np.eye(multiplicity, dtype=complex)
        )
    return tuple(
        coefficient.conj().T
        for coefficient in coefficients_adjoint
    )


def weighted_factor(
    terms: list[Term],
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    parameter: float,
    maximum_degree: int,
) -> tuple[Matrix, float, int]:
    """Return the weighted factor, its bound, and reverse-edge excess."""

    multiplicity = right.shape[1]
    coefficients = [
        np.zeros((multiplicity, multiplicity), dtype=complex)
        for _ in range(maximum_degree)
    ]
    bound = 0.0
    maximum_excess = -maximum_degree

    for scalar, degree, word in terms:
        reverse_count = word.count("j")
        maximum_excess = max(maximum_excess, reverse_count - degree)
        normal = physical_word_factor_coefficients(
            word,
            operator,
            right,
            left,
        )
        for index, coefficient in enumerate(normal, start=1):
            if index <= maximum_degree:
                coefficients[index - 1] += (
                    scalar
                    * parameter ** (degree - index)
                    * coefficient
                )
        bound += (
            abs(scalar)
            * (len(word) + 1)
            * 4**reverse_count
            * parameter ** (degree - reverse_count)
        )

    return np.vstack(coefficients), bound, maximum_excess


def exhaustive_word_audit(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    maximum_length: int = 8,
) -> tuple[int, float, float]:
    """Audit the constructive factor on every short physical word."""

    words_checked = 0
    maximum_error = 0.0
    maximum_ratio = 0.0
    reverse = physical_reverse(operator, right, left)
    for length in range(1, maximum_length + 1):
        for letters in itertools.product("sj", repeat=length):
            word = "".join(letters)
            reverse_count = word.count("j")
            endpoint = (
                left.conj().T
                @ physical_word_value(word, operator, reverse)
                @ right
            )
            if reverse_count == 0:
                reconstruction = np.zeros_like(endpoint)
                factor_norm = 0.0
            else:
                blocks = [
                    transfer_coefficient(
                        operator,
                        right,
                        left,
                        degree,
                    )
                    for degree in range(1, reverse_count + 1)
                ]
                factor = np.vstack(
                    physical_word_factor_coefficients(
                        word,
                        operator,
                        right,
                        left,
                    )
                )
                reconstruction = np.hstack(blocks) @ factor
                factor_norm = float(np.linalg.norm(factor, ord=2))
            maximum_error = max(
                maximum_error,
                float(np.linalg.norm(reconstruction - endpoint)),
            )
            maximum_ratio = max(
                maximum_ratio,
                factor_norm / ((length + 1) * 4**reverse_count),
            )
            words_checked += 1
    return words_checked, maximum_error, maximum_ratio


def direct_operator_jet_error(
    terms: list[Term],
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    maximum_degree: int,
) -> float:
    """Compare the physical word jet with the balanced implementation."""

    dimension = len(operator)
    identity = np.eye(dimension, dtype=complex)
    initial_defect = right @ right.conj().T
    final_defect = left @ left.conj().T
    square_root_two = np.sqrt(2.0)
    metric_root = (
        square_root_two * identity
        + (1 - square_root_two) * initial_defect
        + (2 - square_root_two) * final_defect
    )
    inverse_root = (
        identity / square_root_two
        + (1 - 1 / square_root_two) * initial_defect
        + (0.5 - 1 / square_root_two) * final_defect
    )
    physical_operator = inverse_root @ operator @ metric_root
    implemented = ellipse_operator_coefficients(
        physical_operator,
        maximum_degree,
    )
    implemented = [
        metric_root @ coefficient @ inverse_root
        for coefficient in implemented
    ]

    reverse = physical_reverse(operator, right, left)
    expanded = [
        np.zeros_like(operator) for _ in range(maximum_degree + 1)
    ]
    expanded[0] = operator.copy()
    for scalar, degree, word in terms:
        expanded[degree] += (
            scalar * physical_word_value(word, operator, reverse)
        )
    return max(
        float(np.linalg.norm(actual - expected))
        for actual, expected in zip(expanded, implemented, strict=True)
    )


def audit_terms(
    construction_kind: str,
    terms: list[Term],
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    parameter: float,
    maximum_degree: int,
) -> WeightedEndpointPolynomialRecord:
    """Audit one weighted polynomial."""

    factor, bound, maximum_excess = weighted_factor(
        terms,
        operator,
        right,
        left,
        parameter,
        maximum_degree,
    )
    multiplicity = right.shape[1]
    reverse = physical_reverse(operator, right, left)
    endpoint = np.zeros((multiplicity, multiplicity), dtype=complex)
    for scalar, degree, word in terms:
        endpoint += (
            scalar
            * parameter**degree
            * left.conj().T
            @ physical_word_value(word, operator, reverse)
            @ right
        )

    weighted_blocks = [
        parameter**degree
        * transfer_coefficient(operator, right, left, degree)
        for degree in range(1, maximum_degree + 1)
    ]
    weighted_row = np.hstack(weighted_blocks)
    reconstruction_error = float(
        np.linalg.norm(weighted_row @ factor - endpoint)
    )
    factor_norm = float(np.linalg.norm(factor, ord=2))
    gram_slack = (
        bound**2 * weighted_row @ weighted_row.conj().T
        - endpoint @ endpoint.conj().T
    )
    gram_slack = (gram_slack + gram_slack.conj().T) / 2
    minimum_slack = float(np.linalg.eigvalsh(gram_slack)[0])
    (
        individual_words_checked,
        individual_error,
        individual_bound_ratio,
    ) = exhaustive_word_audit(operator, right, left)
    operator_jet_error = (
        direct_operator_jet_error(
            terms,
            operator,
            right,
            left,
            maximum_degree,
        )
        if construction_kind.endswith("direct_elliptic_map_jet")
        else 0.0
    )
    tolerance = 3e-9
    verified = bool(
        maximum_excess <= 0
        and factor_norm <= bound + tolerance
        and reconstruction_error < tolerance
        and minimum_slack > -tolerance
        and individual_error < tolerance
        and individual_bound_ratio <= 1 + tolerance
        and operator_jet_error < tolerance
    )
    if not verified:
        raise RuntimeError(
            "weighted endpoint-polynomial audit failed: "
            f"{construction_kind=}, {maximum_excess=}, "
            f"{factor_norm=}, {bound=}, "
            f"{reconstruction_error=}, {minimum_slack=}, "
            f"{individual_error=}, {individual_bound_ratio=}, "
            f"{operator_jet_error=}"
        )
    return WeightedEndpointPolynomialRecord(
        construction_kind=construction_kind,
        state_dimension=len(operator),
        defect_dimension=multiplicity,
        maximum_degree=maximum_degree,
        term_count=len(terms),
        ellipse_parameter=format_float(parameter),
        maximum_reverse_excess=maximum_excess,
        factor_norm=format_float(factor_norm),
        theoretical_factor_bound=format_float(bound),
        factor_bound_ratio=format_float(factor_norm / bound),
        reconstruction_error=format_float(reconstruction_error),
        minimum_gram_slack_eigenvalue=format_float(minimum_slack),
        individual_words_checked=individual_words_checked,
        maximum_individual_reconstruction_error=format_float(
            individual_error
        ),
        maximum_individual_bound_ratio=format_float(
            individual_bound_ratio
        ),
        direct_operator_jet_error=format_float(operator_jet_error),
        all_checks_passed=verified,
    )


def standard_records() -> list[WeightedEndpointPolynomialRecord]:
    """Return deterministic weighted-polynomial audits."""

    records: list[WeightedEndpointPolynomialRecord] = []
    direct_terms = direct_map_terms(MAXIMUM_DEGREE)
    cases: list[tuple[str, Matrix, Matrix, Matrix]] = []
    for defect_dimension in (1, 2, 3):
        generator = np.random.default_rng(151_300 + defect_dimension)
        operator, right, left = random_partial_isometry(
            4 * defect_dimension + 9,
            defect_dimension,
            generator,
        )
        cases.append(("unstructured", operator, right, left))

    for defect_dimension in (2, 3):
        operator, right, left, _, _ = rank_chain_case(
            defect_dimension,
            151_400 + defect_dimension,
            1e-6,
        )
        cases.append(("rank_chain", operator, right, left))

    delayed_generator = np.random.default_rng(151_500)
    for defect_dimension in (2, 3):
        operator, right, left = delayed_random_partial_isometry(
            4 * defect_dimension,
            defect_dimension,
            delayed_generator,
        )
        cases.append(("complete_first_delay", operator, right, left))

    for case_index, (case_kind, operator, right, left) in enumerate(cases):
        generator = np.random.default_rng(151_600 + case_index)
        for parameter in (0.13, 0.41, 0.73):
            records.append(
                audit_terms(
                    f"{case_kind}_random_causal_polynomial",
                    causal_random_terms(MAXIMUM_DEGREE, generator),
                    operator,
                    right,
                    left,
                    parameter,
                    MAXIMUM_DEGREE,
                )
            )
            records.append(
                audit_terms(
                    f"{case_kind}_direct_elliptic_map_jet",
                    direct_terms,
                    operator,
                    right,
                    left,
                    parameter,
                    MAXIMUM_DEGREE,
                )
            )
    return records


def write_records(
    records: list[WeightedEndpointPolynomialRecord],
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
            "repeated_crabb_weighted_endpoint_polynomial_s70226.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the weighted endpoint-polynomial audits."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(
        json.dumps(
            {"records": len(records), "sha256": digest},
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
