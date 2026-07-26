#!/usr/bin/env python3
"""Audit the canonical transfer-flux/Markov-energy inequality."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import sympy as sp


Matrix = sp.Matrix


@dataclass(frozen=True)
class CanonicalFluxRecord:
    """One exact bistochastic-channel flux audit."""

    multiplicity: int
    channel_kind: str
    active_row: int
    active_column: int
    test_kind: str
    kraus_count: int
    left_unital_residual_count: int
    right_unital_residual_count: int
    pairing_residual: str
    dirichlet_residual: str
    cauchy_schwarz_slack: str
    flagged_pairing_residual: str
    flagged_cauchy_schwarz_slack: str
    flagged_candidate_face_pairing: str
    flagged_candidate_positive_part: str
    flagged_candidate_bound_square_slack: str
    flagged_correction_residual_count: int
    candidate_face_pairing: str
    candidate_positive_part: str
    candidate_bound_square_slack: str
    correction_residual_count: int
    all_checks_passed: bool


def adjoint(matrix: Matrix) -> Matrix:
    """Return the conjugate transpose."""

    return matrix.conjugate().T


def matrix_sum(matrices: list[Matrix], size: int) -> Matrix:
    """Sum a finite exact matrix list."""

    total = sp.zeros(size, size)
    for matrix in matrices:
        total += matrix
    return total


def residual_count(matrix: Matrix) -> int:
    """Count entries that remain nonzero after exact simplification."""

    return sum(
        int(sp.simplify(entry) != 0)
        for entry in matrix
    )


def frobenius_square(matrix: Matrix) -> sp.Expr:
    """Return the exact squared Frobenius norm."""

    return sp.simplify(sp.trace(adjoint(matrix) * matrix))


def matrix_unit_channel(size: int) -> list[Matrix]:
    """Return Kraus matrices for the exact depolarizing channel."""

    scale = sp.sqrt(size)
    kraus: list[Matrix] = []
    for row in range(size):
        for column in range(size):
            unit = sp.zeros(size, size)
            unit[row, column] = 1 / scale
            kraus.append(unit)
    return kraus


def orthogonal_mixture_channel(size: int) -> list[Matrix]:
    """Return an exact noncommuting mixture of three orthogonal maps."""

    identity = sp.eye(size)
    cycle = sp.zeros(size, size)
    for column in range(size):
        cycle[(column + 1) % size, column] = 1
    sign = sp.diag(*([-1] + [1] * (size - 1)))
    scale = sp.sqrt(3)
    return [identity / scale, cycle / scale, sign / scale]


def phi(kraus: list[Matrix], matrix: Matrix) -> Matrix:
    """Apply the bistochastic channel."""

    return matrix_sum(
        [
            coefficient * matrix * adjoint(coefficient)
            for coefficient in kraus
        ],
        matrix.rows,
    )


def phi_adjoint(kraus: list[Matrix], matrix: Matrix) -> Matrix:
    """Apply the Hilbert--Schmidt adjoint channel."""

    return matrix_sum(
        [
            adjoint(coefficient) * matrix * coefficient
            for coefficient in kraus
        ],
        matrix.rows,
    )


def positive_test(size: int) -> Matrix:
    """Return a deterministic positive non-diagonal exact test."""

    generator = Matrix(
        size,
        size,
        lambda row, column: (
            sp.Rational((row + 1) * (column + 2), size + 3)
            + sp.I
            * sp.Rational(row - column, size + 4)
        ),
    )
    coordinate = sp.zeros(size, size)
    coordinate[0, 0] = 2
    return generator * adjoint(generator) + coordinate


def audit_case(
    size: int,
    active_row: int,
    active_column: int,
    *,
    coordinate_test: bool,
    channel_kind: str = "matrix_unit",
) -> CanonicalFluxRecord:
    """Audit one active Kraus coefficient exactly."""

    if channel_kind == "matrix_unit":
        kraus = matrix_unit_channel(size)
        active_index = active_row * size + active_column
    elif channel_kind == "orthogonal_mixture":
        kraus = orthogonal_mixture_channel(size)
        active_index = active_column
    else:
        raise ValueError(f"unknown channel kind: {channel_kind}")
    active = kraus[active_index]
    identity = sp.eye(size)

    left_unital = matrix_sum(
        [
            coefficient * adjoint(coefficient)
            for coefficient in kraus
        ],
        size,
    ) - identity
    right_unital = matrix_sum(
        [
            adjoint(coefficient) * coefficient
            for coefficient in kraus
        ],
        size,
    ) - identity

    if coordinate_test:
        test = sp.zeros(size, size)
        test[active_row, active_row] = 1
        test_kind = "active_coordinate"
    else:
        test = positive_test(size)
        test_kind = "dense_positive"
    adjoint_state = phi_adjoint(kraus, test)
    defect = test * active - active * adjoint_state
    flux = phi(kraus, adjoint(active) * active) - (
        active * adjoint(active)
    )
    pairing = sp.simplify(sp.trace(test * flux))
    expected_pairing = sp.simplify(
        -sp.re(sp.trace(adjoint(active) * defect))
    )

    markov = test - phi(kraus, phi_adjoint(kraus, test))
    dirichlet = sp.simplify(sp.trace(test * markov))
    defect_energy = sp.simplify(
        sum(
            (
                frobenius_square(
                    test * coefficient
                    - coefficient * adjoint_state
                )
                for coefficient in kraus
            ),
            sp.S.Zero,
        )
    )
    active_energy = frobenius_square(active)
    cauchy_slack = sp.simplify(
        active_energy * dirichlet - pairing**2
    )

    flag = sp.zeros(size, size)
    flag[active_row, active_row] = 1
    flagged_test = flag * test * flag
    flagged_adjoint_state = phi_adjoint(kraus, flagged_test)
    flagged_active = flag * active
    flagged_defect = flag * (
        flagged_test * active
        - active * flagged_adjoint_state
    )
    flagged_flux = flag * (
        phi(
            kraus,
            adjoint(active) * flag * active,
        )
        - active * adjoint(active)
    ) * flag
    flagged_pairing = sp.simplify(
        sp.trace(flagged_test * flagged_flux)
    )
    flagged_expected_pairing = sp.simplify(
        -sp.re(
            sp.trace(
                adjoint(flagged_active) * flagged_defect
            )
        )
    )
    flagged_markov = flagged_test - phi(
        kraus,
        phi_adjoint(kraus, flagged_test),
    )
    flagged_dirichlet = sp.simplify(
        sp.trace(flagged_test * flagged_markov)
    )
    flagged_active_energy = frobenius_square(flagged_active)
    flagged_pairing_residual = sp.simplify(
        flagged_pairing - flagged_expected_pairing
    )
    flagged_cauchy_slack = sp.simplify(
        flagged_active_energy * flagged_dirichlet
        - flagged_pairing**2
    )
    flagged_active_gram = flag * active * adjoint(active) * flag
    flagged_candidate = (
        12 * flagged_active_gram
        - 28
        * flag
        * phi(
            kraus,
            adjoint(active) * flag * active,
        )
        * flag
    )
    flagged_candidate_pairing = sp.simplify(
        sp.trace(flagged_test * flagged_candidate)
    )
    flagged_candidate_positive = sp.Max(flagged_candidate_pairing, 0)
    flagged_candidate_bound_square_slack = sp.simplify(
        28**2 * flagged_active_energy * flagged_dirichlet
        - flagged_candidate_positive**2
    )
    flagged_correction_residual = (
        flagged_candidate
        + 28 * flagged_flux
        + 16 * flagged_active_gram
    )

    active_gram = active * adjoint(active)
    candidate = 12 * active_gram - 28 * phi(
        kraus,
        adjoint(active) * active,
    )
    candidate_pairing = sp.simplify(sp.trace(test * candidate))
    candidate_positive = sp.Max(candidate_pairing, 0)
    candidate_bound_square_slack = sp.simplify(
        28**2 * active_energy * dirichlet
        - candidate_positive**2
    )
    correction_residual = candidate + 28 * flux + 16 * active_gram

    pairing_residual = sp.simplify(pairing - expected_pairing)
    dirichlet_residual = sp.simplify(dirichlet - defect_energy)
    checks_passed = bool(
        residual_count(left_unital) == 0
        and residual_count(right_unital) == 0
        and pairing_residual == 0
        and dirichlet_residual == 0
        and cauchy_slack >= 0
        and flagged_pairing_residual == 0
        and flagged_cauchy_slack >= 0
        and flagged_candidate_bound_square_slack >= 0
        and residual_count(flagged_correction_residual) == 0
        and candidate_bound_square_slack >= 0
        and residual_count(correction_residual) == 0
    )
    if not checks_passed:
        raise RuntimeError(
            "canonical flux audit failed: "
            f"m={size}, row={active_row}, column={active_column}, "
            f"pairing={pairing_residual}, "
            f"Dirichlet={dirichlet_residual}, "
            f"CS={cauchy_slack}, "
            f"flag_pairing={flagged_pairing_residual}, "
            f"flag_CS={flagged_cauchy_slack}, "
            f"flag_face={flagged_candidate_bound_square_slack}, "
            "flag_correction="
            f"{residual_count(flagged_correction_residual)}, "
            f"face={candidate_bound_square_slack}, "
            f"correction={residual_count(correction_residual)}"
        )

    return CanonicalFluxRecord(
        multiplicity=size,
        channel_kind=channel_kind,
        active_row=active_row,
        active_column=active_column,
        test_kind=test_kind,
        kraus_count=len(kraus),
        left_unital_residual_count=residual_count(left_unital),
        right_unital_residual_count=residual_count(right_unital),
        pairing_residual=str(pairing_residual),
        dirichlet_residual=str(dirichlet_residual),
        cauchy_schwarz_slack=str(cauchy_slack),
        flagged_pairing_residual=str(flagged_pairing_residual),
        flagged_cauchy_schwarz_slack=str(flagged_cauchy_slack),
        flagged_candidate_face_pairing=str(flagged_candidate_pairing),
        flagged_candidate_positive_part=str(
            flagged_candidate_positive
        ),
        flagged_candidate_bound_square_slack=str(
            flagged_candidate_bound_square_slack
        ),
        flagged_correction_residual_count=residual_count(
            flagged_correction_residual
        ),
        candidate_face_pairing=str(candidate_pairing),
        candidate_positive_part=str(candidate_positive),
        candidate_bound_square_slack=str(candidate_bound_square_slack),
        correction_residual_count=residual_count(correction_residual),
        all_checks_passed=checks_passed,
    )


def standard_records() -> list[CanonicalFluxRecord]:
    """Return exact noncommuting channel audits."""

    return [
        audit_case(2, 0, 1, coordinate_test=False),
        audit_case(3, 0, 2, coordinate_test=True),
        audit_case(4, 0, 3, coordinate_test=True),
        audit_case(4, 2, 1, coordinate_test=True),
        audit_case(
            2,
            1,
            1,
            coordinate_test=False,
            channel_kind="orthogonal_mixture",
        ),
        audit_case(
            3,
            1,
            2,
            coordinate_test=False,
            channel_kind="orthogonal_mixture",
        ),
    ]


def write_records(
    records: list[CanonicalFluxRecord],
    output: Path,
) -> str:
    """Write deterministic JSON Lines and return its SHA-256."""

    output.parent.mkdir(parents=True, exist_ok=True)
    payload = "".join(
        json.dumps(asdict(record), sort_keys=True) + "\n"
        for record in records
    )
    output.write_text(payload, encoding="utf-8")
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the exact audit."""

    args = parse_args()
    records = standard_records()
    digest = None
    if args.output is not None:
        digest = write_records(records, args.output)
    print(
        json.dumps(
            {
                "all_checks_passed": all(
                    record.all_checks_passed for record in records
                ),
                "dataset_sha256": digest,
                "records": [asdict(record) for record in records],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
