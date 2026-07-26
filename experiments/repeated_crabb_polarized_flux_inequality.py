#!/usr/bin/env python3
"""Audit the polarized transfer-flux/Markov-energy inequality."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import sympy as sp

from repeated_crabb_canonical_flux_inequality import (
    adjoint,
    frobenius_square,
    matrix_sum,
    matrix_unit_channel,
    orthogonal_mixture_channel,
    phi,
    phi_adjoint,
    positive_test,
)


Matrix = sp.Matrix


@dataclass(frozen=True)
class PolarizedFluxRecord:
    """One exact polarized-response audit."""

    multiplicity: int
    channel_kind: str
    kraus_count: int
    single_pairing_residual: str
    single_cauchy_schwarz_slack: str
    summed_pairing_residual: str
    summed_cauchy_schwarz_slack: str
    flagged_pairing_residual: str
    flagged_cauchy_schwarz_slack: str
    dirichlet_residual: str
    all_checks_passed: bool


def symmetrize(matrix: Matrix) -> Matrix:
    """Return the Hermitian part."""

    return (matrix + adjoint(matrix)) / 2


def polarization(size: int, seed: int) -> Matrix:
    """Return one deterministic exact non-Hermitian polarization."""

    return Matrix(
        size,
        size,
        lambda row, column: (
            sp.Rational(
                (row + 1) * (column + seed + 2),
                size + seed + 3,
            )
            + sp.I
            * sp.Rational(
                (seed + 1) * (row - 2 * column + 1),
                size + seed + 5,
            )
        ),
    )


def flux(
    kraus: list[Matrix],
    active: Matrix,
    coefficient: Matrix,
) -> Matrix:
    """Return one polarized channel flux."""

    return (
        symmetrize(active * adjoint(coefficient))
        - phi(
            kraus,
            symmetrize(adjoint(coefficient) * active),
        )
    )


def audit_case(
    size: int,
    channel_kind: str,
) -> PolarizedFluxRecord:
    """Audit one exact bistochastic channel."""

    if channel_kind == "matrix_unit":
        kraus = matrix_unit_channel(size)
    elif channel_kind == "orthogonal_mixture":
        kraus = orthogonal_mixture_channel(size)
    else:
        raise ValueError(f"unknown channel kind: {channel_kind}")

    test = positive_test(size)
    adjoint_state = phi_adjoint(kraus, test)
    defects = [
        test * active - active * adjoint_state
        for active in kraus
    ]
    dirichlet = sp.simplify(
        sp.trace(
            test
            * (
                test
                - phi(kraus, adjoint_state)
            )
        )
    )
    defect_energy = sp.simplify(
        sum(
            (frobenius_square(defect) for defect in defects),
            sp.S.Zero,
        )
    )

    active = kraus[min(1, len(kraus) - 1)]
    active_index = min(1, len(kraus) - 1)
    coefficient = polarization(size, 1)
    single_flux = flux(kraus, active, coefficient)
    single_pairing = sp.simplify(sp.trace(test * single_flux))
    single_expected = sp.simplify(
        sp.re(
            sp.trace(
                adjoint(coefficient) * defects[active_index]
            )
        )
    )
    single_slack = sp.simplify(
        frobenius_square(coefficient) * dirichlet
        - single_pairing**2
    )

    coefficients = [
        polarization(size, index + 2)
        for index in range(len(kraus))
    ]
    summed_flux = matrix_sum(
        [
            flux(kraus, active_kraus, coefficient_k)
            for active_kraus, coefficient_k in zip(
                kraus,
                coefficients,
                strict=True,
            )
        ],
        size,
    )
    summed_pairing = sp.simplify(sp.trace(test * summed_flux))
    summed_expected = sp.simplify(
        sum(
            (
                sp.re(
                    sp.trace(
                        adjoint(coefficient_k) * defect_k
                    )
                )
                for coefficient_k, defect_k in zip(
                    coefficients,
                    defects,
                    strict=True,
                )
            ),
            sp.S.Zero,
        )
    )
    coefficient_energy = sp.simplify(
        sum(
            (
                frobenius_square(coefficient_k)
                for coefficient_k in coefficients
            ),
            sp.S.Zero,
        )
    )
    summed_slack = sp.simplify(
        coefficient_energy * dirichlet - summed_pairing**2
    )

    flag = sp.zeros(size, size)
    flag[size - 1, size - 1] = 1
    flagged_test = flag * test * flag
    flagged_adjoint_state = phi_adjoint(kraus, flagged_test)
    flagged_defects = [
        flagged_test * active_kraus
        - active_kraus * flagged_adjoint_state
        for active_kraus in kraus
    ]
    flagged_dirichlet = sp.simplify(
        sp.trace(
            flagged_test
            * (
                flagged_test
                - phi(kraus, flagged_adjoint_state)
            )
        )
    )
    flagged_coefficient = flag * coefficient
    flagged_flux = flux(kraus, active, flagged_coefficient)
    flagged_pairing = sp.simplify(
        sp.trace(flagged_test * flagged_flux)
    )
    flagged_expected = sp.simplify(
        sp.re(
            sp.trace(
                adjoint(flagged_coefficient)
                * flagged_defects[active_index]
            )
        )
    )
    flagged_slack = sp.simplify(
        frobenius_square(flagged_coefficient)
        * flagged_dirichlet
        - flagged_pairing**2
    )

    residuals = (
        sp.simplify(single_pairing - single_expected),
        sp.simplify(summed_pairing - summed_expected),
        sp.simplify(flagged_pairing - flagged_expected),
        sp.simplify(dirichlet - defect_energy),
    )
    checks_passed = bool(
        all(residual == 0 for residual in residuals)
        and single_slack >= 0
        and summed_slack >= 0
        and flagged_slack >= 0
    )
    if not checks_passed:
        raise RuntimeError(
            "polarized flux audit failed: "
            f"m={size}, kind={channel_kind}, "
            f"residuals={residuals}, "
            f"slacks={(single_slack, summed_slack, flagged_slack)}"
        )

    return PolarizedFluxRecord(
        multiplicity=size,
        channel_kind=channel_kind,
        kraus_count=len(kraus),
        single_pairing_residual=str(residuals[0]),
        single_cauchy_schwarz_slack=str(single_slack),
        summed_pairing_residual=str(residuals[1]),
        summed_cauchy_schwarz_slack=str(summed_slack),
        flagged_pairing_residual=str(residuals[2]),
        flagged_cauchy_schwarz_slack=str(flagged_slack),
        dirichlet_residual=str(residuals[3]),
        all_checks_passed=checks_passed,
    )


def standard_records() -> list[PolarizedFluxRecord]:
    """Return exact audits for two noncommuting channel families."""

    return [
        audit_case(2, "matrix_unit"),
        audit_case(3, "matrix_unit"),
        audit_case(2, "orthogonal_mixture"),
        audit_case(3, "orthogonal_mixture"),
    ]


def write_records(
    records: list[PolarizedFluxRecord],
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
