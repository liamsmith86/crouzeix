#!/usr/bin/env python3
"""Audit the two-ended selection-independent quartic invariant.

Insert L298's grade-one retightening into the moving canonical pair,
cancel its cubic residual with the convenient metric-only Green
representative, and retain the resulting quartic endpoint.  Its
pairing against every L206 reducing-commutant weight is

    (15/4) tr(Y B_1 B_1*) + (1/2) tr(Y (B_1 B_1*)^2).

The symbolic audit proves the scalar identity by exact finite word
reduction.  It also proves that the same forcing has lower corner

    6 B_1*B_1 + (9/4)(B_1*B_1)^2.

After the compulsory parallel lower neutralization, its upper
commutant class is the prior-flag cost

    (9/4) B_1B_1* + (7/4)(B_1B_1*)^2.

The numerical audit independently reconstructs the moving canonical
pair and checks scalar and weighted reducing sums.
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
from scipy.linalg import block_diag

from crabb_block_hardy_equality import format_float
from repeated_crabb_boundary_metric_flag import rank_chain_case
from repeated_crabb_canonical_cubic_selection import cyclic_trace_classes
from repeated_crabb_canonical_quintic_preimage import (
    canonical_factor_series,
    canonical_metric_slack_and_operator,
    exact_canonical_factor_lifts,
)
from repeated_crabb_delayed_slack_anticommutator import (
    E,
    F,
    Polynomial,
    S,
    STAR,
    add,
    adjoint,
    multiply,
    operator_series,
    scale,
)
from repeated_crabb_elliptic_cokernel import (
    dual_stein_inverse,
    stein_inverse,
)
from repeated_crabb_endpoint_null_gauge import (
    delayed_random_partial_isometry,
)
from repeated_crabb_oriented_retightening_transport import (
    oriented_retightening_direction,
)
from repeated_crabb_transfer_channel_covariance import (
    random_partial_isometry,
    transfer_coefficient,
)


Matrix = np.ndarray
DEGREE = 4


@dataclass(frozen=True)
class QuarticInvariantRecord:
    """One scalar or reducing-weight quartic audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    parameter_scale: str
    first_transfer_norm: str
    cubic_forcing_norm: str
    quartic_endpoint_norm: str
    weighted_quartic_pairing: str
    predicted_weighted_pairing: str
    weighted_pairing_error: str
    lower_corner_error: str
    lower_neutral_upper_pairing: str
    predicted_lower_neutral_upper_pairing: str
    lower_neutral_pairing_error: str
    observability_commutator_error: str
    canonical_factor_error: str
    exact_retightening_stein_residual_words: int
    exact_cubic_hermitian_residual_words: int
    exact_dual_telescope_residual_words: int
    exact_quartic_cyclic_residual_classes: int
    exact_lower_corner_residual_words: int
    all_checks_passed: bool


def finite_stein_sum(
    forcing: Polynomial,
    maximum_steps: int = 8,
) -> Polynomial:
    """Return a terminating exact Stein sum."""

    result: Polynomial = {}
    term = forcing
    for _ in range(maximum_steps):
        result = add(result, term)
        term = multiply(multiply(STAR, term), S)
        if not term:
            return result
    raise RuntimeError("the exact retightening Stein sum did not terminate")


def left_gram_lift() -> Polynomial:
    """Return the state lift ``W(B_1B_1*)W*``."""

    return multiply(
        multiply(
            multiply(
                multiply(F, STAR),
                E,
            ),
            S,
        ),
        F,
    )


def retightening_frame_lift() -> Polynomial:
    """Return the simplified L298 lift ``F_1 V*``."""

    return scale(
        Fraction(-1, 2),
        multiply(
            multiply(
                multiply(S, F),
                STAR,
            ),
            E,
        ),
    )


@lru_cache(maxsize=1)
def exact_certificate() -> tuple[int, int, int, int, int]:
    """Return exact residual counts for the quartic trace certificate."""

    operator = operator_series(DEGREE)
    canonical_factor = exact_canonical_factor_lifts(DEGREE)
    frame = retightening_frame_lift()
    frame_forcing = add(frame, adjoint(frame))
    metric = finite_stein_sum(frame_forcing)
    metric_stein_residual = add(
        metric,
        scale(-1, multiply(multiply(STAR, metric), S)),
        scale(-1, frame_forcing),
    )

    cubic = scale(
        -1,
        add(
            multiply(multiply(adjoint(operator[1]), metric), S),
            multiply(multiply(STAR, metric), operator[1]),
            multiply(canonical_factor[1], adjoint(frame)),
            multiply(frame, adjoint(canonical_factor[1])),
        ),
    )
    cubic_hermitian_residual = add(
        cubic,
        scale(-1, adjoint(cubic)),
    )

    # If X_0=-G_S(N_3), its quartic cross pairs through the dual
    # Stein solution 2(S^2+(S*)^2).  This finite telescope is what
    # removes the apparently nonlocal Green term from the trace.
    dual_solution = scale(
        2,
        add(
            multiply(S, S),
            multiply(STAR, STAR),
        ),
    )
    dual_forcing = add(
        multiply(S, adjoint(operator[1])),
        multiply(operator[1], STAR),
    )
    dual_residual = add(
        dual_solution,
        scale(
            -1,
            multiply(
                multiply(S, dual_solution),
                STAR,
            ),
        ),
        scale(-1, dual_forcing),
    )

    finite_quartic = scale(
        -1,
        add(
            multiply(multiply(adjoint(operator[2]), metric), S),
            multiply(multiply(STAR, metric), operator[2]),
            multiply(
                multiply(adjoint(operator[1]), metric),
                operator[1],
            ),
            multiply(canonical_factor[2], adjoint(frame)),
            multiply(frame, adjoint(canonical_factor[2])),
            multiply(frame, adjoint(frame)),
        ),
    )
    effective_trace_lift = add(
        finite_quartic,
        multiply(cubic, dual_solution),
    )
    first_left_gram = left_gram_lift()
    predicted_left_lift = add(
        scale(Fraction(15, 4), first_left_gram),
        scale(
            Fraction(1, 2),
            multiply(first_left_gram, first_left_gram),
        ),
    )
    cyclic_residual = cyclic_trace_classes(
        add(
            effective_trace_lift,
            scale(-1, predicted_left_lift),
        )
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
    lower_corner = scale(
        -1,
        multiply(
            multiply(
                E,
                add(
                    multiply(
                        multiply(adjoint(operator[2]), metric),
                        S,
                    ),
                    multiply(
                        multiply(STAR, metric),
                        operator[2],
                    ),
                    multiply(
                        multiply(adjoint(operator[1]), metric),
                        operator[1],
                    ),
                    multiply(
                        canonical_factor[2],
                        adjoint(frame),
                    ),
                    multiply(
                        frame,
                        adjoint(canonical_factor[2]),
                    ),
                    multiply(frame, adjoint(frame)),
                ),
            ),
            E,
        ),
    )
    predicted_lower_corner = add(
        scale(6, first_right_gram),
        scale(
            Fraction(9, 4),
            multiply(first_right_gram, first_right_gram),
        ),
    )
    lower_corner_residual = add(
        lower_corner,
        scale(-1, predicted_lower_corner),
    )
    return (
        len(metric_stein_residual),
        len(cubic_hermitian_residual),
        len(dual_residual),
        len(cyclic_residual),
        len(lower_corner_residual),
    )


def quartic_endpoint(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
) -> tuple[Matrix, Matrix, Matrix, float, float]:
    """Return the retained quartic endpoint and its first transfer."""

    _, slack, operator, _ = canonical_metric_slack_and_operator(
        partial,
        right,
        left,
        DEGREE,
    )
    canonical_factor, factor_error = canonical_factor_series(
        slack,
        right,
    )
    direction = oriented_retightening_direction(
        partial,
        right,
        left,
        1,
    )
    metric = direction.metric
    frame = direction.frame

    cubic = -(
        operator[1].conj().T @ metric @ partial
        + partial.conj().T @ metric @ operator[1]
        + canonical_factor[1] @ frame.conj().T
        + frame @ canonical_factor[1].conj().T
    )
    cubic_metric = -stein_inverse(partial, cubic)
    quartic = -(
        operator[2].conj().T @ metric @ partial
        + partial.conj().T @ metric @ operator[2]
        + operator[1].conj().T @ metric @ operator[1]
        + canonical_factor[2] @ frame.conj().T
        + frame @ canonical_factor[2].conj().T
        + frame @ frame.conj().T
        + operator[1].conj().T @ cubic_metric @ partial
        + partial.conj().T @ cubic_metric @ operator[1]
    )
    endpoint = (
        left.conj().T
        @ stein_inverse(partial, quartic)
        @ left
    )
    endpoint = (endpoint + endpoint.conj().T) / 2
    lower_corner = right.conj().T @ quartic @ right
    lower_corner = (lower_corner + lower_corner.conj().T) / 2
    first = transfer_coefficient(partial, right, left, 1)
    return (
        endpoint,
        lower_corner,
        first,
        float(np.linalg.norm(cubic)),
        factor_error,
    )


def audit_case(
    construction_kind: str,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    parameter_scale: float,
    weight: Matrix | None = None,
) -> QuarticInvariantRecord:
    """Audit one scalar or reducing-commutant pairing."""

    endpoint, lower_corner, first, cubic_norm, factor_error = (
        quartic_endpoint(
        partial,
        right,
        left,
        )
    )
    if weight is None:
        weight = np.eye(right.shape[1], dtype=complex)
    observability = dual_stein_inverse(
        partial,
        left @ weight @ left.conj().T,
    )
    commutator_error = float(
        np.linalg.norm(
            observability @ partial - partial @ observability
        )
    )

    left_gram = first @ first.conj().T
    positive_representative = (
        15 * left_gram / 4
        + left_gram @ left_gram / 2
    )
    actual = float(np.trace(weight @ endpoint).real)
    predicted = float(
        np.trace(weight @ positive_representative).real
    )
    pairing_error = abs(actual - predicted)
    right_gram = first.conj().T @ first
    predicted_lower_corner = (
        6 * right_gram
        + 9 * right_gram @ right_gram / 4
    )
    lower_corner_error = float(
        np.linalg.norm(lower_corner - predicted_lower_corner)
    )
    lower_channel = (
        left.conj().T
        @ stein_inverse(
            partial,
            right @ lower_corner @ right.conj().T,
        )
        @ left
    )
    lower_neutral_endpoint = -endpoint + lower_channel
    lower_neutral_actual = float(
        np.trace(weight @ lower_neutral_endpoint).real
    )
    predicted_lower_neutral = (
        9 * left_gram / 4
        + 7 * left_gram @ left_gram / 4
    )
    lower_neutral_predicted = float(
        np.trace(weight @ predicted_lower_neutral).real
    )
    lower_neutral_error = abs(
        lower_neutral_actual - lower_neutral_predicted
    )

    exact_residuals = exact_certificate()
    tolerance = 8e-8
    verified = bool(
        factor_error < tolerance
        and pairing_error < tolerance
        and lower_corner_error < tolerance
        and lower_neutral_error < tolerance
        and commutator_error < tolerance
        and exact_residuals == (0, 0, 0, 0, 0)
    )
    if not verified:
        raise RuntimeError(
            "the quartic-invariant audit failed: "
            f"{construction_kind=}, {factor_error=}, "
            f"{pairing_error=}, {lower_corner_error=}, "
            f"{lower_neutral_error=}, {commutator_error=}, "
            f"{exact_residuals=}"
        )
    return QuarticInvariantRecord(
        construction_kind=construction_kind,
        state_dimension=len(partial),
        defect_dimension=right.shape[1],
        parameter_scale=format_float(parameter_scale),
        first_transfer_norm=format_float(float(np.linalg.norm(first))),
        cubic_forcing_norm=format_float(cubic_norm),
        quartic_endpoint_norm=format_float(
            float(np.linalg.norm(endpoint))
        ),
        weighted_quartic_pairing=format_float(actual),
        predicted_weighted_pairing=format_float(predicted),
        weighted_pairing_error=format_float(pairing_error),
        lower_corner_error=format_float(lower_corner_error),
        lower_neutral_upper_pairing=format_float(
            lower_neutral_actual
        ),
        predicted_lower_neutral_upper_pairing=format_float(
            lower_neutral_predicted
        ),
        lower_neutral_pairing_error=format_float(
            lower_neutral_error
        ),
        observability_commutator_error=format_float(commutator_error),
        canonical_factor_error=format_float(factor_error),
        exact_retightening_stein_residual_words=exact_residuals[0],
        exact_cubic_hermitian_residual_words=exact_residuals[1],
        exact_dual_telescope_residual_words=exact_residuals[2],
        exact_quartic_cyclic_residual_classes=exact_residuals[3],
        exact_lower_corner_residual_words=exact_residuals[4],
        all_checks_passed=verified,
    )


def reducing_sum_case(seed: int) -> tuple[Matrix, Matrix, Matrix, Matrix]:
    """Return a two-block colligation and a nonscalar cokernel weight."""

    first_partial, first_right, first_left = random_partial_isometry(
        7,
        1,
        np.random.default_rng(seed),
    )
    second_partial, second_right, second_left = random_partial_isometry(
        11,
        2,
        np.random.default_rng(seed + 1),
    )
    partial = block_diag(first_partial, second_partial)
    right = block_diag(first_right, second_right)
    left = block_diag(first_left, second_left)
    weight = block_diag(
        np.asarray([[2 / 3]], dtype=complex),
        np.eye(2, dtype=complex) * (7 / 4),
    )
    return partial, right, left, weight


def standard_records() -> list[QuarticInvariantRecord]:
    """Return deterministic unstructured, flag, and reducing records."""

    records: list[QuarticInvariantRecord] = []
    for multiplicity in range(1, 5):
        partial, right, left = random_partial_isometry(
            3 * multiplicity + 3,
            multiplicity,
            np.random.default_rng(120_000 + multiplicity),
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

    for multiplicity in range(2, 6):
        for parameter_scale in (1.0, 0.3, 0.1):
            partial, right, left, _, _ = rank_chain_case(
                multiplicity,
                120_200 + multiplicity,
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

    generator = np.random.default_rng(120_500)
    for multiplicity in (2, 3, 4):
        partial, right, left = delayed_random_partial_isometry(
            4 * multiplicity,
            multiplicity,
            generator,
        )
        records.append(
            audit_case(
                "complete_first_delay",
                partial,
                right,
                left,
                1,
            )
        )

    for index in range(3):
        partial, right, left, weight = reducing_sum_case(
            120_700 + 2 * index
        )
        records.append(
            audit_case(
                "weighted_reducing_sum",
                partial,
                right,
                left,
                1,
                weight,
            )
        )
    return records


def write_records(
    records: list[QuarticInvariantRecord],
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
            "repeated_crabb_retightening_quartic_invariant_s70226.jsonl"
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
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
