# Literature refresh — 2026-09-05

Targeted primary-source update to the July literature baseline. This records
source statements and scope, not a substitute for a manuscript-level proof
audit. Search snippets were used only to locate sources. No counterexample to
the scalar conjecture was identified in this refresh.

## 1. Material change: the scalar problem now has public proofs

**Shanmu Jin, _The Numerical Range Is a 2-Spectral Set_.** The latest source is
version 4, submitted August 6 and posted August 7, 2026; the original posting
was July 27. The page explicitly labels the manuscript **not peer-reviewed**.
Theorem 1 states the sharp scalar polynomial inequality for every finite
complex matrix, and the holomorphic-neighborhood extension. Theorem 2 is a
mass-parameterized positive-real completion result using an auxiliary
eigenbasis; the numerical-range layer supplies mass two. Simple-spectrum
approximation removes the auxiliary diagonalizability restriction.
[Primary manuscript and version record](https://www.preprints.org/manuscript/202607.1919).

There is unusually strong independent public corroboration, distinct from
formal peer review: Alex Townsend and Anne Greenbaum's August 15 report says
that they and Michel Crouzeix checked Jin's proof and consider it correct.
Their report also identifies an independent subsequent proof.
[Expert-authored report, pp. 1–3](https://alextownsend.net/essays/SIAMNews_CrouzeixConjecture.pdf).

**Emiel Lorist and Felix L. Schwenninger, _A solution to Crouzeix's
conjecture_, arXiv:2608.03841.** First submitted August 4; latest version 2
is August 17, seven pages. Theorem 3 covers rational scalar functions of
bounded Hilbert-space operators. Its reusable Lemma 1 says that, for a
contraction `Q` and isometry `V`, uniformly bounded defects

`E_n = 2 V* Q*^n V - T*^n`

commuting with `T` for every positive integer `n` imply `||T|| <= 2`.
The double-layer representation realizes these defects as
`alpha(f^n)(A)`. The key new ingredient is telescoping over **all powers**,
not nonnegativity of the single extremal correction. Section 3(ii) explicitly
allows negative corrections, including in dimension three. Section 3(iv)
warns that matrix amplification destroys the required commutation; this is
**not** a general complete-2 theorem. Version 2 extends the lemma to general
Hilbert spaces using an ultrapower; a finite-matrix proof does not need that
step. [Version history](https://arxiv.org/abs/2608.03841),
[primary v2, Lemma 1, Theorem 3, and §3](https://arxiv.org/html/2608.03841v2).

**Project consequence (assessment):** the July statement that `1+sqrt(2)` is
the best known scalar universal bound is obsolete. A direct audit of the
short all-powers argument is a substantially better next action than
continuing a local Hessian route under an unchanged open-problem premise.
Any adaptation must credit the source and must not claim independent
discovery or automatically extend scalar statements to complete ones.

## 2. September abstract extension

**Catalin Badea, Ryan O'Loughlin, Jani Virtanen,
_On the abstract approach to spectral constants: a proof of the
Clouâtre–Ostermann–Ransford conjecture_, arXiv:2609.03637v1.**
The manuscript is dated September 1; arXiv submission is September 3, within
this refresh's cutoff. Theorem 2.1 assumes a bounded unital homomorphism
`theta` from a uniform algebra, a bounded unital antilinear map `alpha`, and
contractivity of
`Phi(h) = (theta(h) + theta(alpha(h))*)/2`.
It concludes `||theta|| <= 2`, for arbitrary Hilbert-space targets.
Contractivity of `alpha` is not required. Theorem 3.3 isolates a positive-real
layer of mass two whose remainder commutes with `theta(f)*` on every power.
The proof uses Herglotz–Naimark dilation and the Lorist–Schwenninger lemma.
This is another **scalar norm** conclusion, not a dimension-free cb bound.
[Submission record](https://arxiv.org/abs/2609.03637),
[primary theorem statements and proof](https://arxiv.org/html/2609.03637v1).

## 3. Complete inequality: dimension-three preprint, not all dimensions

**Per Åhag, Rafał Czyż, Jani Virtanen, _Square Functions and the Complete
Crouzeix Conjecture in Dimension Three_, arXiv:2608.27346v2.**
First submitted August 27; latest revision September 1. Use the v2 author
list: search metadata can still display Antti Perälä from v1.

Theorem 1.1 gives the complete-2 inequality for every matrix of order at most
three. It further reduces any possible complete counterexample to a
rectangular matrix level with `r,s >= 2` and `r+s <= n`; dimension four
therefore leaves only level `2 x 2`. Theorem 1.2 gives an abstract
uniform-algebra counterpart in target dimension at most three. Theorem 1.3
provides sharp row and column square-function estimates. Theorem 1.4 equates
scalar and complete attainment of two in these dimensions. Section 12
certifies the rational-domain containment for the `3 x 3` KMS example using
exact rational Sturm sequences. The introduction explicitly retains
`1+sqrt(2)` as the best general complete bound.

These are very recent preprint claims: the full multi-variable kernel and
rectangular-extremal arguments were **not independently audited here**.
[Version record](https://arxiv.org/abs/2608.27346),
[primary v2](https://arxiv.org/html/2608.27346v2).

## 4. Precise older baselines and corrections

### Configuration constants

Malman–Mashreghi–O'Loughlin–Ransford,
_Double-layer potentials, configuration constants and applications to
numerical ranges_, published in IMRN 2025(8), rnaf084:
Theorem 2 gives `a(Omega) < 1` for **every compact convex domain with
nonempty interior**, without a smooth-boundary restriction. Theorem 3 gives
the scalar shape-dependent bound `1+sqrt(1+a(W))`; Corollary 4 gives a
strict fixed-dimension improvement of `1+sqrt(2)`. Thin quadrilaterals have
`a(Omega)` arbitrarily close to one, so this particular estimate cannot
produce a better universal constant. Theorem 1 identifies the real and
complex configuration constants, with `a(Omega) <= c(Omega)`. Theorem 6
gives `c(Omega) <= 1 - (1/(2*pi))*integral(ds/R_Omega)`, using the
containing-tangent-disk radius. Thus the ledger's smooth-domain question
should be replaced by the actual nonsmooth-domain theorem, while keeping
the analytic and ordinary configuration constants distinct.
[Published record](https://academic.oup.com/imrn/article/2025/8/rnaf084/8109682),
[primary theorem text](https://arxiv.org/html/2407.19049v2).

### Weighted shifts

Crouzeix–Greenbaum, arXiv:2508.12768v1, first submitted August 18, 2025:
the class is precisely `M = P_d diag(alpha_1,...,alpha_d)` for a cyclic
permutation `P_d` and arbitrary complex scalar weights, including zeros.
The paper proves `psi_cb(M) = psi(M) <= 2` for its numerical range.
Lemma 1 additionally identifies both disk-calculus constants with
`max_{0 <= k <= d-1} ||M^k||` when `|product alpha_j| <= 1`; otherwise both
are infinite. This does not cover arbitrary block-weighted shifts or all
weighted tridiagonal matrices. The arXiv submission record has only v1;
an August 24, 2026 date displayed inside the generated HTML is not evidence
of a new revision.
[Primary paper](https://arxiv.org/html/2508.12768v1),
[submission history](https://arxiv.org/abs/2508.12768).

### “KLS matrices (2025–26 work)” is not a verified historical theorem

The matched source is Crouzeix–Greenbaum–Li,
_Numerical bounds on the Crouzeix ratio for a class of matrices_,
arXiv:2311.13890v2 (December 5, 2023). The body calls these **KMS** matrices
and studies the strict-upper-triangle all-ones family. Its near-optimal
conformal-map estimates are explicitly numerical. The proposed `3 x 3`
inner-domain construction still estimates derivative suprema from sampled
differences; those samples alone are not certified global derivative
bounds. Consequently this source must not be cited as a rigorously audited
2025–26 theorem for a general “KLS” class. The new dimension-three preprint
above separately addresses exact containment certification.
[Primary original manuscript, §§1–2](https://arxiv.org/html/2311.13890v2).

## 5. Retrieval and completion cautions

- Search snippets can omit `2 x 2` from statements of the new Schwarz–Jack
  proof, falsely suggesting a general proof. That article's publisher
  listing also gives a November 2026 issue date; do not use that issue date
  as evidence of a result available after this September cutoff. The
  actually opened page describes a new proof of the **two-dimensional**
  theorem, not a general result.
  [Publisher scope](https://doi.org/10.1016/j.jmaa.2026.130786).
- A public preprint, an expert endorsement, and a locally audited proof are
  different evidence levels. The present file verifies availability and
  scope; it does not by itself complete `goal.txt` or validate the campaign's
  hundreds of local lemmas.
- Historical artifacts can remain archived as a research record, but any
  active summary must distinguish the scalar theorem from the still
  unrestricted complete problem and flag the change in literature date.
