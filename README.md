# Crouzeix research notebook

An AI-assisted mathematics experiment exploring Crouzeix's conjecture through
proof attempts, numerical searches, symbolic calculations and independent
checks. The research used models including GPT-5.6 Sol and Claude Opus.

## The question

A matrix has an associated region in the complex plane called its **numerical
range**. Crouzeix's conjecture says that applying a polynomial to the matrix
can amplify vectors by at most twice the largest absolute value of that
polynomial on this region.

The experiment's goal was to find a general proof or an exact counterexample,
using AI agents to develop different approaches and check each other's work.
The [original task](goal.txt) records the research setup.

## What we tried

- **Revisit the existing proof machinery.** We examined where the earlier
  bound of about 2.414 lost information and looked for ways to recover two.
- **Search for counterexamples.** Small nonnormal matrices and polynomial
  functions were optimized numerically, with known equality cases used as
  reference points.
- **Solve structured families.** Elliptic numerical ranges and weighted
  matrices made it possible to derive explicit formulas and check whole
  parameter ranges.
- **Study perturbations of equality cases.** We investigated how the bound
  changes when matrices already attaining two are slightly altered.
- **Construct matrix certificates.** Similarity transformations, positive
  kernels and semidefinite optimization provided alternative ways to establish
  bounds. Some of these approaches addressed the stronger version in which
  polynomial coefficients can themselves be matrices.

## Findings in the notes

| Finding | What it describes | Notes |
| --- | --- | --- |
| Exact elliptic-family constant | A closed formula for the optimal constant for the elliptic Crabb family in every matrix size, with an explicit transformation attaining it. | [L117: elliptic Crabb axis](proof/crabb_elliptic_axis.md) |
| Weighted four-by-four family | A computer-assisted proof of the complete bound of two for a family with three independent positive weights and an elliptic numerical range. | [L20/L59: similarity and certificates](proof/slice_similarity_duality.md) |
| Local second-order behaviour | Formulas expressing the change near a single Crabb equality block as negative quadratic terms, in arbitrary size. | [L65: second variation](proof/general_crabb_second_variation.md) |
| A three-dimensional neighbourhood | A local argument covering small complex perturbations of the three-by-three Crabb block, including directions where the quadratic term vanishes. | [L73: local theorem](proof/p3_crabb_local_theorem.md) |
| Equality and residual structure | Connections between equality cases, Toeplitz matrices and explicit square identities that simplify the local calculations. | [L187: equality geometry](proof/crabb_full_disk_leading_residual_tube.md), [L355: residual square](proof/gau_wu_residual_parallel_sum.md) |

The [literature ledger](LITERATURE_LEDGER.md) compares these derivations with
earlier work, including Kenan Li's elliptic-family calculations and classical
results on disk numerical ranges. The [approach ledger](APPROACH_LEDGER.md)
also records failed shortcuts and the examples that identified their limits.

## The later general proofs

After the July research, the repository added a finite-matrix reconstruction
of the general proof by Emiel Lorist and Felix Schwenninger. Their argument
combines identities for all powers of an operator. Shanmu Jin's separate proof
uses a positive-kernel sampling construction.

- [Lorist–Schwenninger manuscript](https://arxiv.org/html/2608.03841v2)
- [Jin manuscript](https://www.preprints.org/manuscript/202607.1919)
- [Our reconstruction](proof/CROUZEIX_PROOF.md) and
  [independent audit](proof/audit_global_proof_20260905.md)

## Repository guide

- `proof/` — derivations, special-case results and proof audits.
- `experiments/` — Python programs, exact checks and recorded experiment data.
- `LEMMA_LEDGER.md` — statement IDs, scope and supporting evidence.
- `LITERATURE_LEDGER.md` — sources and comparisons with existing work.
- `APPROACH_LEDGER.md` — research decisions and retained unsuccessful approaches.
- `COUNTEREXAMPLE_SEARCH.md` — search methods, results and coverage.
- `archive/` — the full July lemma and approach ledgers.

## Run selected checks

The recorded environment uses Python 3.14. Dependencies are pinned in
`requirements.txt`. Run commands from the repository root:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt

# Reconstruct the exact three-dimensional second-order identity.
.venv/bin/python experiments/p3_second_order_identity.py

# Check the elliptic-family formulas at selected sizes and parameters.
.venv/bin/python experiments/crabb_elliptic_axis_theorem.py \
  --minimum-size 3 --maximum-size 8 \
  --ellipse-parameters 0.15 0.6 --precision 100 \
  --output tmp/elliptic-axis.jsonl

# Check the algebra used in the later all-powers proof.
.venv/bin/python experiments/global_power_lemma_exact_audit_20260905.py \
  --output tmp/all-powers-exact.json
```

Individual proof notes give the corresponding reproduction commands and
describe whether each check uses exact arithmetic or numerical approximation.
