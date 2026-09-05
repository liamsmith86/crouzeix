# Crouzeix campaign — current state

Updated 2026-09-05. **Resolution is not yet verified in this repository.**

## Objective and scope

`goal.txt` asks for a rigorous proof for every finite complex square matrix and
every scalar polynomial, or a certified strict counterexample. Equality at two
is not a counterexample. The stronger completely bounded statement is not a
completion requirement.

## Current priority: audit newly located global proofs

The July literature baseline is stale. A primary manuscript now claims a global
scalar solution: Emiel Lorist and Felix Schwenninger,
[A solution to Crouzeix's conjecture, v2, 2026-08-17](https://arxiv.org/html/2608.03841v2).
Its key ingredient is a bounded commuting-error lemma for all powers of a
compressed contraction. We are reconstructing and adversarially checking the
argument, not treating the title or search results as proof.

Parallel tasks:

1. Verify source versions, dates, scope, and other recent claims.
2. Independently reconstruct the global all-powers argument and boundary passage.
3. Audit the latest local work and its claimed relationship to the full objective.
4. Produce a self-contained finite-matrix/polynomial proof and reproducible
   algebra/regression checks if the new argument survives.

## Preserved local frontier (not a global reduction)

The latest local work is L355/A304. At a fixed finite nondegenerate Gau–Wu
equality model, the recorded metric/zero elimination is

`e_phi(C) - 4 J_C(tau) = (1/2) ||T_phi bar(tau) + Delta_phi(C)||²_(D_+)`.

The remaining local obligations are (a) complex-linear shape factorization
`Delta_phi(C)=B_phi A_phi(C)`, (b) type-(1,1) of the disk-eliminated common
physical coefficient, and (c) its nonpositive sign. Those are not yet proved.
Even proving all three would not alone settle arbitrary matrices: degenerate
strata, local-to-global control, and a global extremal reduction remain separate.
The old phrase “exactly one current mathematical gate” meant the selected local
subproblem and must not be read as the only debt in `goal.txt`.

Read the compact chain:

- `proof/gau_wu_canonical_residual_normal_form.md` (L354).
- `proof/gau_wu_residual_parallel_sum.md` (L355).
- `proof/gau_wu_residual_complex_covariance.md` (L352 scope).
- `proof/gau_wu_first_jet_loss.md` (L353 alternative; not the common term).

Do not drop the first-curvature reserve, split the endpoint signs, infer exact
phase covariance from floating-point pairing, or retry the disproved fixed
inverse-kernel lower bounds. See the retained pitfalls in `APPROACH_LEDGER.md`.

## Evidence and navigation

- `LEMMA_LEDGER.md`: latest local statements; full history linked there.
- `LITERATURE_LEDGER.md`: dated source record; refresh currently in progress.
- `COUNTEREXAMPLE_SEARCH.md`: historical search and certification gaps; no
  certified violation. The best reported original direct-search ratio is 1.9747.
- `proof/`: preserved derivations and audits. Historical claims are not silently
  promoted by this cleanup.
- `experiments/`: existing reproducible scripts/data; dependencies in
  `requirements.txt`; use `.venv/bin/python` and single-threaded BLAS for audits.
- `checkpoints/RESTART_PACKET.md`: short current restart instruction.
- `archive/2026-07-27/`: byte-identical snapshots of the four bloated state files,
  with SHA-256 manifest. No mathematical history or experiment has been deleted.

## Completion gate

Before marking the goal complete: inspect the full proof; verify all imported
facts, every-dimensional scope, degenerate cases, equality and arbitrary degree;
obtain independent reconstruction; run exact algebra and adversarial tests;
record a requirement-by-requirement audit and attribution. Numerical tests or a
published claim alone do not satisfy this gate.
