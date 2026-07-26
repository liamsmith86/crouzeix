#!/usr/bin/env python3
"""Audit the gap-free relative flux bound for L303's quartic response.

At the fixed partial-retightening scale theta=1/2, the complete
lower-neutral quartic residual is an exact sum of 28 commutators.
Every commutator contains either

    G = E S F = V B_1* W*

or its adjoint.  Pairing with an observability Gramian therefore
supplies both factors needed for a rank-stable estimate:
``||B_1||`` and the L281 observability defect ``||Z_Y||``.

The only nonlocal Green term is handled by the off-commutant error in
L302's dual Stein telescope.  This checker regenerates the exact
relative-commutator certificate and audits the complete numerical
pairing, telescope error, flux bounds, and response synthesis.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_boundary_metric_flag import rank_chain_case
from repeated_crabb_canonical_quintic_preimage import (
    canonical_metric_slack_and_operator,
)
from repeated_crabb_cubic_markov_flux import minimum_response_column
from repeated_crabb_delayed_slack_anticommutator import (
    E,
    F,
    Polynomial,
    S,
    STAR,
    add,
    multiply,
    scale,
)
from repeated_crabb_elliptic_cokernel import dual_stein_inverse
from repeated_crabb_endpoint_null_gauge import (
    delayed_random_partial_isometry,
)
from repeated_crabb_retightening_full_quartic import (
    exact_full_quartic_polynomials,
    full_quartic_deltas,
)
from repeated_crabb_retightening_quartic_invariant import (
    numerical_quartic_components,
    reducing_sum_case,
)
from repeated_crabb_transfer_channel_covariance import (
    random_partial_isometry,
)


Matrix = np.ndarray
THETA = Fraction(1, 2)
RELATIVE_CERTIFICATE_NORM = Fraction(177, 16)


@dataclass(frozen=True)
class RelativeCommutatorTerm:
    """One term ``coefficient [letter, prefix bridge suffix]``."""

    coefficient: Fraction
    letter: str
    bridge: str
    prefix: str
    suffix: str


@dataclass(frozen=True)
class QuarticFluxRecord:
    """One numerical audit of the complete quartic flux bound."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    parameter_scale: str
    first_transfer_norm: str
    stein_inverse_norm: str
    observability_defect_norm: str
    actual_residual_pairing: str
    relative_commutator_pairing: str
    telescope_error_pairing: str
    complete_pairing_error: str
    relative_flux_bound_slack: str
    telescope_flux_bound_slack: str
    response_synthesis_error: str
    response_column_norm: str
    response_column_bound_slack: str
    exact_certificate_term_count: int
    exact_certificate_residual_words: int
    all_checks_passed: bool


RELATIVE_COMMUTATOR_TERMS = (
    RelativeCommutatorTerm(Fraction(1, 8), "a", "h", "", "ss"),
    RelativeCommutatorTerm(Fraction(-1, 2), "s", "g", "", "aa"),
    RelativeCommutatorTerm(Fraction(-3, 8), "a", "h", "ss", ""),
    RelativeCommutatorTerm(Fraction(1, 4), "s", "h", "", "ssss"),
    RelativeCommutatorTerm(Fraction(9, 16), "s", "h", "", "ssaa"),
    RelativeCommutatorTerm(Fraction(-1, 2), "a", "g", "a", "aaa"),
    RelativeCommutatorTerm(Fraction(-1, 4), "s", "h", "ss", "ss"),
    RelativeCommutatorTerm(Fraction(-1, 2), "a", "g", "aa", "aa"),
    RelativeCommutatorTerm(Fraction(1, 4), "s", "g", "saaa", ""),
    RelativeCommutatorTerm(Fraction(5, 16), "s", "h", "aass", ""),
    RelativeCommutatorTerm(Fraction(1, 4), "a", "g", "aaaa", ""),
    RelativeCommutatorTerm(Fraction(1, 4), "s", "g", "", "aassss"),
    RelativeCommutatorTerm(
        Fraction(21, 16),
        "a",
        "g",
        "",
        "aasssa",
    ),
    RelativeCommutatorTerm(Fraction(1, 4), "s", "g", "a", "aaass"),
    RelativeCommutatorTerm(Fraction(1, 4), "s", "h", "ass", "sss"),
    RelativeCommutatorTerm(Fraction(1, 4), "a", "g", "aaa", "aas"),
    RelativeCommutatorTerm(Fraction(-1, 4), "a", "g", "saaa", "aa"),
    RelativeCommutatorTerm(Fraction(1, 4), "s", "h", "asss", "ss"),
    RelativeCommutatorTerm(Fraction(1, 4), "a", "h", "aaaass", ""),
    RelativeCommutatorTerm(
        Fraction(1, 4),
        "a",
        "h",
        "",
        "sssaaass",
    ),
    RelativeCommutatorTerm(
        Fraction(1, 4),
        "a",
        "h",
        "",
        "ssaaasss",
    ),
    RelativeCommutatorTerm(
        Fraction(-11, 16),
        "s",
        "g",
        "",
        "aasssaaa",
    ),
    RelativeCommutatorTerm(
        Fraction(1, 4),
        "s",
        "g",
        "",
        "aaasssaa",
    ),
    RelativeCommutatorTerm(
        Fraction(1, 4),
        "a",
        "g",
        "saaa",
        "aaas",
    ),
    RelativeCommutatorTerm(
        Fraction(1, 4),
        "a",
        "h",
        "saaasss",
        "s",
    ),
    RelativeCommutatorTerm(
        Fraction(3, 8),
        "a",
        "h",
        "sssaaass",
        "",
    ),
    RelativeCommutatorTerm(
        Fraction(1, 4),
        "a",
        "h",
        "ssaaasss",
        "",
    ),
    RelativeCommutatorTerm(
        Fraction(-25, 16),
        "s",
        "g",
        "aaasssaa",
        "",
    ),
)


def monomial(word: str) -> Polynomial:
    """Return one exact word monomial."""

    return {word: Fraction(1)}


def commutator(left: Polynomial, right: Polynomial) -> Polynomial:
    """Return the exact commutator ``left right - right left``."""

    return add(
        multiply(left, right),
        scale(-1, multiply(right, left)),
    )


def bridge_polynomials() -> dict[str, Polynomial]:
    """Return ``G=ESF`` and ``G*=FS*E``."""

    bridge = multiply(multiply(E, S), F)
    return {
        "g": bridge,
        "h": multiply(multiply(F, STAR), E),
    }


def relative_commutator_certificate() -> Polynomial:
    """Reconstruct the 28-term exact relative commutator."""

    letters = {"s": S, "a": STAR}
    bridges = bridge_polynomials()
    terms = []
    for term in RELATIVE_COMMUTATOR_TERMS:
        interior = multiply(
            multiply(
                monomial(term.prefix),
                bridges[term.bridge],
            ),
            monomial(term.suffix),
        )
        terms.append(
            scale(
                term.coefficient,
                commutator(letters[term.letter], interior),
            )
        )
    return add(*terms)


def exact_certificate() -> tuple[int, Fraction]:
    """Return the exact residual count and coefficient l1 norm."""

    target = exact_full_quartic_polynomials(
        THETA
    ).actual_neutralized_upper_residual
    residual = add(
        relative_commutator_certificate(),
        scale(-1, target),
    )
    coefficient_norm = sum(
        abs(term.coefficient)
        for term in RELATIVE_COMMUTATOR_TERMS
    )
    return len(residual), coefficient_norm


def evaluate_polynomial(polynomial: Polynomial, operator: Matrix) -> Matrix:
    """Evaluate one exact word polynomial numerically."""

    value = np.zeros_like(operator, dtype=complex)
    adjoint = operator.conj().T
    for word, coefficient in polynomial.items():
        product = np.eye(len(operator), dtype=complex)
        for letter in word:
            product = product @ (operator if letter == "s" else adjoint)
        value += float(coefficient) * product
    return value


def stein_inverse_norm(operator: Matrix) -> float:
    """Return the Frobenius operator norm of the stable Stein inverse."""

    dimension = len(operator)
    stein_matrix = np.eye(dimension * dimension, dtype=complex) - np.kron(
        operator.T,
        operator.conj().T,
    )
    smallest = float(
        np.linalg.svd(stein_matrix, compute_uv=False)[-1]
    )
    return 1 / smallest


def random_hermitian(size: int, seed: int) -> Matrix:
    """Return one deterministic Hermitian endpoint test."""

    generator = np.random.default_rng(seed)
    raw = (
        generator.standard_normal((size, size))
        + 1j * generator.standard_normal((size, size))
    )
    return (raw + raw.conj().T) / 2


def audit_case(
    construction_kind: str,
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    parameter_scale: float,
    seed: int,
) -> QuarticFluxRecord:
    """Audit the exact and numerical relative-flux identities."""

    theta = float(THETA)
    (
        _,
        _,
        neutral_upper_delta,
        _,
        first_transfer,
        _,
    ) = full_quartic_deltas(
        operator,
        right,
        left,
        theta,
    )
    left_gram = first_transfer @ first_transfer.conj().T
    retained_class = (
        -5 * left_gram / 4
        - 3 * left_gram @ left_gram / 4
    )
    response_target = neutral_upper_delta[4] / 4 - retained_class

    endpoint_test = random_hermitian(right.shape[1], seed)
    observability = dual_stein_inverse(
        operator,
        left @ endpoint_test @ left.conj().T,
    )
    identity = np.eye(len(operator), dtype=complex)
    right_projection = right @ right.conj().T
    defect_column = (
        (identity - right_projection)
        @ observability
        @ right
    )
    actual_pairing = float(
        np.trace(endpoint_test @ response_target).real
    )

    exact_target = exact_full_quartic_polynomials(
        THETA
    ).actual_neutralized_upper_residual
    relative_pairing = float(
        np.trace(
            observability
            @ evaluate_polynomial(exact_target, operator)
        ).real
    )

    components = numerical_quartic_components(
        operator,
        right,
        left,
    )
    _, _, operator_series, _ = canonical_metric_slack_and_operator(
        operator,
        right,
        left,
        4,
    )
    first_operator = operator_series[1]
    cubic_metric = components.metric_three
    cubic_forcing = -(
        cubic_metric
        - operator.conj().T @ cubic_metric @ operator
    )
    dual_solution = 2 * (
        operator @ operator
        + operator.conj().T @ operator.conj().T
    )
    original_telescope_term = -(
        first_operator.conj().T @ cubic_metric @ operator
        + operator.conj().T @ cubic_metric @ first_operator
    )
    effective_telescope_term = cubic_forcing @ dual_solution
    telescope_pairing = float(
        np.trace(
            observability
            @ (
                original_telescope_term
                - effective_telescope_term
            )
        ).real
    )
    complete_error = abs(
        actual_pairing
        - relative_pairing
        - theta * telescope_pairing
    )

    first_norm = float(np.linalg.norm(first_transfer))
    defect_norm = float(np.linalg.norm(defect_column))
    cubic_norm = float(np.linalg.norm(cubic_forcing))
    gamma = stein_inverse_norm(operator)
    relative_bound = (
        float(RELATIVE_CERTIFICATE_NORM)
        * first_norm
        * defect_norm
    )
    telescope_bound = (
        theta
        * 33
        * gamma
        * cubic_norm
        * defect_norm
    )
    relative_slack = relative_bound - abs(relative_pairing)
    telescope_slack = telescope_bound - abs(
        theta * telescope_pairing
    )

    if first_norm < 1e-10:
        response_column_norm = 0.0
        synthesis_error = float(np.linalg.norm(response_target))
    else:
        response_column, synthesis_error = minimum_response_column(
            operator,
            right,
            left,
            response_target,
        )
        response_column_norm = float(np.linalg.norm(response_column))
    proven_coefficient = (
        float(RELATIVE_CERTIFICATE_NORM)
        + 330 * gamma * gamma
        + 33 * gamma
    ) / 2
    response_bound_slack = (
        proven_coefficient * first_norm - response_column_norm
    )

    exact_residual, exact_norm = exact_certificate()
    tolerance = 2e-7
    verified = bool(
        exact_residual == 0
        and exact_norm == RELATIVE_CERTIFICATE_NORM
        and complete_error < tolerance
        and relative_slack > -tolerance
        and telescope_slack > -tolerance
        and synthesis_error < tolerance
        and response_bound_slack > -tolerance
    )
    if not verified:
        raise RuntimeError(
            "the quartic relative-flux audit failed: "
            f"{construction_kind=}, {complete_error=}, "
            f"{relative_slack=}, {telescope_slack=}, "
            f"{synthesis_error=}, {response_bound_slack=}"
        )
    return QuarticFluxRecord(
        construction_kind=construction_kind,
        state_dimension=len(operator),
        defect_dimension=right.shape[1],
        parameter_scale=format_float(parameter_scale),
        first_transfer_norm=format_float(first_norm),
        stein_inverse_norm=format_float(gamma),
        observability_defect_norm=format_float(defect_norm),
        actual_residual_pairing=format_float(actual_pairing),
        relative_commutator_pairing=format_float(relative_pairing),
        telescope_error_pairing=format_float(
            theta * telescope_pairing
        ),
        complete_pairing_error=format_float(complete_error),
        relative_flux_bound_slack=format_float(relative_slack),
        telescope_flux_bound_slack=format_float(telescope_slack),
        response_synthesis_error=format_float(synthesis_error),
        response_column_norm=format_float(response_column_norm),
        response_column_bound_slack=format_float(
            response_bound_slack
        ),
        exact_certificate_term_count=len(
            RELATIVE_COMMUTATOR_TERMS
        ),
        exact_certificate_residual_words=exact_residual,
        all_checks_passed=verified,
    )


def standard_records() -> list[QuarticFluxRecord]:
    """Return deterministic general, rank-changing, and delayed cases."""

    records = []
    for multiplicity in range(1, 5):
        operator, right, left = random_partial_isometry(
            3 * multiplicity + 3,
            multiplicity,
            np.random.default_rng(130_000 + multiplicity),
        )
        records.append(
            audit_case(
                "unstructured",
                operator,
                right,
                left,
                1,
                130_100 + multiplicity,
            )
        )

    for multiplicity in range(2, 5):
        for parameter_scale in (1.0, 0.1, 0.01):
            operator, right, left, _, _ = rank_chain_case(
                multiplicity,
                130_200 + multiplicity,
                parameter_scale,
            )
            records.append(
                audit_case(
                    "rank_chain",
                    operator,
                    right,
                    left,
                    parameter_scale,
                    130_300 + 10 * multiplicity,
                )
            )

    generator = np.random.default_rng(130_500)
    for multiplicity in (2, 3, 4):
        operator, right, left = delayed_random_partial_isometry(
            4 * multiplicity,
            multiplicity,
            generator,
        )
        records.append(
            audit_case(
                "complete_first_delay",
                operator,
                right,
                left,
                1,
                130_600 + multiplicity,
            )
        )

    for index in range(3):
        operator, right, left, _ = reducing_sum_case(
            130_700 + 2 * index
        )
        records.append(
            audit_case(
                "reducing_direct_sum",
                operator,
                right,
                left,
                1,
                130_800 + index,
            )
        )
    return records


def write_records(
    records: list[QuarticFluxRecord],
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
            "repeated_crabb_retightening_quartic_flux_s70226.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the relative-flux audit."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
