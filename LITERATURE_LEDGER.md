# LITERATURE_LEDGER.md

Papers, theorems, assumptions, possible gaps. Status verified via web 2026-07-23.

## Core chain
| Ref | Result | Method | Notes / gaps to exploit |
|---|---|---|---|
| Crouzeix 2004 ("Bounds for analytical functions of matrices") | Conjecture stated; proved for 2×2 (constant 2, sharp) | Direct estimates, conformal maps | 2×2 case: W(A) is an ellipse; sharp via Jordan block [[0,2],[0,0]] |
| Crouzeix 2007 | Universal constant 11.08 | Cauchy integral + conformal splitting | Superseded |
| Crouzeix–Palencia 2017 (SIMAX/FoCM) | W(A) is a (1+√2)-spectral set | f(A)+g(A)* via positive double-layer kernel; abstract lemma | THE target to improve. g = Cauchy transform of f̄ conjugate |
| Ransford–Schwenninger 2018 (SIMAX 1708.08633) | The abstract lemma (‖f(A)+g(A)*‖≤2‖f‖, g=C(f̄)) cannot yield better than 1+√2 in abstract setting | Explicit abstract construction | KEY OBSTRUCTION: need concrete structure of the pair (f(A), g(A)) beyond the abstract hypotheses |
| Malman–Mashreghi–O'Loughlin–Ransford 2024 | For each fixed N: C_N < 1+√2 | Compactness + strict inequality analysis of C–P equality conditions | Non-constructive; no uniform bound. Their equality analysis of C–P may reveal extremal structure |
| MMOR 2025 | Configuration constants via Neumann–Poincaré operator; domain-dependent bounds | Double-layer potential spectral analysis | Domain-dependent c(Ω) < 1+√2 for smooth domains? Check exact statement |
| Schwenninger–de Vries, arXiv:2409.15954 "DLP for spectral constants revisited" | §6: ρ(f₀,x₀) := ∫Re(K_Ω(f₀)*f₀)dμ₀ = ReC; **K² + ρ ≤ 2K** (Thm 6.1: K ≤ 1+√(1−ρ)); ρ ≥ −1 (convex); ρ ≥ 0 ⟹ Crouzeix; disk ⟹ (K−1)ρ = 0; \|ρ\| ≤ a(Ω) (Rem 6.2) | Extremal pair + DLP | READ & VERIFIED (pp. 15–17). Closest baseline to campaign's P1; our increment = β-subtraction (q < 2). Also: Appendix A smooth-approximation lemma (useful for Ω ↓ W(A) rigor); Prop A.2 spectral-constant equivalence |
| arXiv:2410.10678 | "The algebraic numerical range as a spectral set in Banach algebras" | — | CORRECTED ID (audit): earlier ledger entry misattributed this ID to a parameterized C–P extension |
| O'Loughlin–Rani, arXiv:2603.15536 (2026) "q-Numerical Ranges and Spectral Sets" | Bound K = t + √(t²+a(Ω)), t = 1+γ(1)/2 | q-numerical range | CAUTION (audit): sentence after Thm 3.3 claiming K<2-type conclusion may be typo/gap — verify before use |
| Hartz–McCarthy, arXiv:2606.02922 (2026) "From Clouâtre–Ostermann–Ransford to Okubo–Ando" | For a continuous unital operator-algebra homomorphism `theta` and every scalar functional `beta`, `||theta^(n)|| <= max(1,||(theta+beta I)^(n)||)` at every matrix level; recovers the sharp Okubo–Ando similarity theorem | Potapov–Möbius variational argument at matrix levels | AUDITED FOR L21: when `||theta||cb>1`, their inequality plus `beta=0` gives `inf_beta ||theta+beta I||cb=||theta||cb`, so scalar-shift search is an exact restatement, not a relaxation. Crouzeix–Palencia supplies an operator-valued correction `R(f)=theta(alpha(f))*`, not the required scalar range except in the disk case. A resolution/cross-solver SDP probe also numerically excludes every positive-state scalarization on one dense 3×3 test (`proof/general_similarity_probe.md`). Do not claim a general CP-to-HM bridge. |

## Special cases proved (constant 2)
- Normal matrices (trivial: R ≤ 1; equality-2 impossible). von Neumann: W(A) ⊇ spectrum.
- All 2×2 matrices (Crouzeix 2004).
- n×n with W(A) a disk (Okubo--Ando: numerical radius at most one implies
  similarity to a contraction with condition at most two, so the disk is a complete
  `2`-spectral set).  Campaign L123 is not a new proof of this upper bound: it identifies
  a concrete all-size phase-palindromic stratified family on which the upper constant is attained,
  with an explicit rank-one metric and matching characteristic-polynomial Blaschke product.
  Targeted searches for disk-matrix equality classifications and Blaschke extremals found no
  matching parametrized equality theorem; treat the classification as apparently new pending
  publication-level review.
- Nearly Jordan blocks (Choi–Greenbaum).
- Weighted shift matrices (Choi; complete 2-spectral via arXiv:2508.12768 Aug 2025).
- Certain 3×3: tridiagonal with elliptic W(A) centered at eigenvalue; 3×3 KLS matrices (2025-26 work).
- The Crabb-derived elliptic family `A=(C+bC*)/(2sqrt(b))`, where `C` has the fixed
  superdiagonal weights `(sqrt(2),1,...,1,sqrt(2))`: Kenan Li's 2021 thesis, Chapter 3,
  proposes an explicit contraction similarity with squared condition
  `k(b^(2n-2))/b^(n-1)<4`. **Audit correction (2026-07-22):** the thesis explicitly says
  its load-bearing identities (3.5)--(3.7) were proved only for `n=2,...,6` and tested
  numerically beyond that; the displayed all-dimensional conclusion therefore has an
  unclosed proof gap unless another source supplies those identities.  Do not cite it as
  an established every-dimensional theorem.
- **Campaign L59 (2026-07-22):** the arbitrary-weight 4×4 elliptic family
  `A=S_a+cS_a^T`, with three independent positive weights in `S_a`, has a condition-two
  contraction similarity and hence a complete 2-spectral numerical-range ellipse.  The
  closest-source audit found no prior arbitrary-weight theorem: Li's thesis fixes the Crabb
  weights and explicitly leaves the general elliptical problem open.  Treat this as apparently
  new pending a publication-level literature audit; it is not a result for arbitrary 4×4
  matrices.  Proof chain L20/L21/L27--L59; `proof/slice_core_projective_reduction.md`.
- **Campaign L73 (2026-07-22):** every complex `3×3` matrix in a full neighbourhood of the
  single Crabb block has its numerical range as a complete `2`-spectral set.  This is much wider
  locally than the weighted-shift/perturbed-Jordan classes, but it is only a neighbourhood of one
  equality block.  Targeted searches for "neighbourhood", Crabb perturbations, and complete
  spectral sets found no matching theorem; the closest variational paper proves first-order
  scalar stationarity rather than a neighbourhood inequality.  Treat L73 as apparently new
  pending publication-level review.  `proof/p3_crabb_local_theorem.md`.
- Compressions of the shift / model space operators with certain Blaschke conditions (Bickel–Gorkin school).
- Matrices attaining ‖A‖ = numerical radius conditions (Crabb-type equality cases).

## Numerical evidence
- Greenbaum–Overton (2017-2018, "Numerical investigation of Crouzeix's conjecture"): nonsmooth optimization (BFGS/GRANSO) over A and p; global max of R appears to be exactly 2, attained on Crabb-type / Jordan-block configurations. No R > 2 ever found.
- AIM workshop reports (aimath.org/pastworkshops/crouzeixrep.pdf): consensus directions.

## Known equality structure (R = 2 attained/approached)
- A = [[0,2],[0,0]], p = z. R = 2 exactly. W(A) = disk radius 1.
- Crabb matrix family (nilpotent Jordan-like with specific superdiagonal weights √2,1,...,1,√2? — verify): ‖p(A)‖ = 2 with p = z^k, W(A) = unit disk. Attains 2 for higher powers.
- General principle: equality cases known are nilpotent + monomial + disk. Conjecturally all extremals are these (Greenbaum–Overton observation).

## Equivalent / stronger formulations (Track F cautions)
- Completely bounded version: "W(A) is a complete 2-spectral set" — open, possibly strictly stronger; Paulsen theory: cb-version ⇔ similarity to operator with dilation... (K-spectral ⇒ complete K'-spectral with K' possibly larger). Do NOT conflate.
- Equivalent: enough to prove for f holomorphic on int W(A), continuous on closure (Mergelyan/convexity).
- Equivalent: enough for A with W(A) having smooth boundary (perturb A ↦ A + εB dense case? verify continuity argument), and p with ‖p‖_{W(A)} = 1.
- Reduction: can assume p has degree < n (Cayley–Hamilton changes sup — INVALID as direct reduction; but valid: sup over W(A) of the reduced rep may increase, so inequality for reduced rep does NOT imply original. Only the direction ‖·‖ side is fixed.) — see APPROACH_LEDGER pitfalls.

## Read & digested (Epochs 1–3)
- [x] Ransford–Schwenninger arXiv:1708.08633 (full) — abstract lemma, sharpness (non-unital, Ω ⊅ W(T)), Question 4.1.
- [x] Schwenninger–de Vries arXiv:2302.05389 (full) — extremal pairs/measures, Theorem 5, Prop 8/9, Remarks 12–13.
- [x] MMOR arXiv:2407.19049 (pp. 1–8) — configuration constants c_R = c_C, a(Ω) < 1, ellipse formula, thin-domain obstruction, Theorem 6 curvature bound.
- [x] BGG+20 arXiv:2006.04901 (pp. 4–12) — Blaschke extremal structure (Crouzeix 04 Thm 2.1), cancellation Thm 4.1, extremal measure Thm 4.5, compressed shifts.
- [x] Schwenninger–de Vries arXiv:2409.15954 (pp. 15–20) — see table (P1 baseline).
- [x] NIST DLMF §23.8.1 — Weierstrass Fourier expansion used to prove the positive series in
  EL4; combined with `wp(x)=1/sn²(x|m)−(1+m)/3` for half-periods `(K,iK′)`.
- [x] NIST DLMF §22.11.6 — reciprocal Jacobi `nd=1/dn` Fourier series.  After multiplication
  by `k'`, DCT-I sampling at L117's elliptic Lobatto nodes gives L137's exact axis-defect
  Newton edge `d_(2r)=4(-1)^r c^r+O(c^(r+2))`.
- [x] Mashreghi--Moucha--O'Loughlin--Ransford--Roth, arXiv:2506.23831 — Schwarz--Jack
  convexity: the inverse Riemann map of a bi-circular domain (hence an ellipse) is convex on the
  positive radius, so the direct map is concave. Used only for `r≤p` in L24; no quantitative
  two-node estimate is imported.
- [x] Kanas--Sugawa (2006), *On conformal representations of the interior of an ellipse* — the
  centered inverse ellipse map has positive odd Taylor coefficients; an alternative precise
  justification of `r≤p`. L26 uses this coefficient positivity to obtain the cubic minorant
  that closes the `(p−r)/c` upper-block bound; the subsequent Möbius barrier and theta
  certificate are campaign-derived.
- [x] Crouzeix--Greenbaum, arXiv:2508.12768 (2025/26) — for scalar cyclic weighted shifts the
  disk complete bound equals the largest consecutive weight product via an explicit diagonal
  similarity. L25 uses the elementary nilpotent specialization on the `c→0` face; the paper's
  general theorem does not cover the interior block-weighted slice.
- [x] Lewis--Overton, arXiv:1901.00050 / SIMAX 2020, *Partial Smoothness of the Numerical
  Radius at Matrices Whose Fields of Values are Disks* — proves that disk matrices form a local
  analytic manifold and the numerical radius is partly smooth at every nonzero scalar
  superdiagonal matrix, including a single Crabb block; gives explicit tangent/normal spaces
  under a simple support-eigenvalue hypothesis. L115 now uses their exact codimension-`2p`
  centered-disk tangent theorem (Theorems 6.4, 7.4, 8.7), adjoins the two translation
  coordinates, and intersects the resulting codimension-`2p−2` circular-range tangent with
  L65. It does **not** analyze repeated support multiplicity, the conformally pulled-back L21
  similarity square, or this kernel intersection, so no direct overlap with the new
  second-order reduction was found. L122 rewrites their polynomial support-vector
  parametrization as the concrete chart
  `X(H)=2(H+R*HR)^−1/2 H R (H+R*HR)^−1/2` and derives the rank-one similarity quartic on
  its Toeplitz flat submanifold. The disk factorization is an equivalent concrete use of
  their/Crouzeix's parametrization; the Stein quartic and its phase-palindromic null are the
  campaign-derived additions.
- [x] Greenbaum--Lewis--Overton, *Variational Analysis of the Crouzeix Ratio*, Math.
  Programming 2017 — proves Clarke regularity and nonnegative directional derivative of the
  **scalar Crouzeix ratio** at the Crabb/monomial candidate `(p(z)=z^(n−1), A=Crabb⊕0)`, while
  explicitly leaving second-order analysis on the active manifold open. This is the closest
  conceptual precedent for L61. It does not treat repeated Crabb support multiplicity, the
  domain-dependent Riemann pullback, the completely bounded similarity optimum `t_*`, or
  L61's density-matrix/Jensen closed form. Phrase novelty narrowly: L61 appears to be a new
  first-order theorem for the stronger L21 quantity, not the first variational stationarity
  result at a Crabb configuration. Campaign L62 now gives a finite second-order
  conformal/metric-SDP reduction for that stronger quantity. L63–L65 now prove its
  second-order sign at every single Crabb block, using an arbitrary-size path-Green/mode
  factorization. No matching completely bounded second-variation result was found in the
  targeted audit; retain the narrow novelty wording pending a broader search and do not present
  this local theorem as a proof of the scalar conjecture.
- [x] Rodin, *Behavior of the Riemann Mapping Function under Complex Analytic Deformations of
  the Domain*, Complex Variables 5 (1986), 189--195, DOI
  `10.1080/17476938608814139`; Wu, *Analytic Dependence of Riemann Mappings for Bounded
  Domains and Minimal Surfaces*, CPAM 46 (1993), 1303--1326, DOI
  `10.1002/cpa.3160461002` -- normalized Riemann maps vary analytically under analytic
  boundary deformations.  L73 needs only the elementary strictly-convex near-circle case,
  for which the boundary equation plus Fourier splitting gives a direct Banach-IFT proof.
- [x] Arlinskii--Golinskii--Tsekanovskii, arXiv:math/0611439, *Contractions with rank one
  defect operators and truncated CMV matrices* -- Theorem 6.4 models every completely
  nonunitary contraction with one-dimensional defect spaces by a truncated CMV matrix whose
  coefficients are the Schur parameters of its characteristic function.  This is a
  potentially useful sparse coordinate system for L118's rank-one-defect optimizer.
  Proposition 6.2, however, gives unitary equivalence only under multiplication of **all**
  Schur parameters by one common phase.  It does not supply an independent phase torus by
  grade and therefore does not prove A95's grade-block Hessian conjecture.  Any use here must
  explicitly derive the Crabb equality tangent-to-Schur-parameter Jacobian and the optimized
  endpoint quadratic; first-order sparsity of the CMV matrix alone is insufficient.
- [x] Bruns--Vetter, *Determinantal Rings* (Springer LNM 1327, 1988), chapters on maximal
  minors and their powers -- the ideal of maximal minors of a generic matrix has coincident
  symbolic and ordinary powers.  L155 uses only the elementary low-degree consequence for a
  generic `2 x m` matrix: a polynomial whose value and gradient vanish on the rank-one cone
  lies in `I^(2)=I^2`, so it has degree at least four.  This imported commutative-algebra fact
  removes the formerly separate Crabb-apex jet obligation; the application to the
  phase-palindromic cone is campaign-derived.
- [x] O'Loughlin--Rani, arXiv:2603.15536, *q-Numerical Ranges and Spectral Sets* -- a
  March-2026 extension of the Crouzeix--Palencia framework to scaled `q`-numerical ranges.
  It still treats constant two for the ordinary numerical range as a conjecture and does not
  supply a rank-one-defect Hessian, Crabb equality-stratum theorem, or complete-`2` result
  relevant to A94--A98.  Logged as a current-status cross-check, not as an ingredient.

## To read next
- [ ] GKL arXiv:1701.01365 (Glader–Kurula–Lindström 2018, 3×3 tridiagonal elliptic W(A)) —
  novelty calibration vs our sym3 ρ-proof (their mechanism is cb/dilation-based; ours is new)
- [x] Kenan Li, 2021 UW thesis, Chapters 2.3 and 3 — Chapter 2.3 identifies three extremal
  phases for a proportional 3×3 elliptic family: centered degree 1, shifted real Möbius degree
  1, and even degree 2. Chapter 3 derives the candidate sharp contraction similarity for the
  all-dimensional one-parameter family `A=(C+bC*)/(2sqrt(b))`, with condition
  `sqrt(k(b^(2n-2))/b^(n-1))<=2`. Its 4×4 member has the fixed weight ratio
  `(sqrt(2),1,sqrt(2))`; it does not cover arbitrary positive weights `(a1,a2,a3)`.
  Crucially, page 46 states that identities (3.5)--(3.7), used to prove the similarity for
  general parity, were proved only for `n=2,...,6` and tested numerically for larger `n`.
  No later source closing that gap was found in targeted exact-formula and citation searches.
  Campaign L117 now supplies a different all-size proof through a rank-one Szegő kernel,
  DCT-I diagonalization, and a periodized-sech/Poisson formula.  Treat this as an apparently
  new closure pending publication-level review, not as an imported theorem.  L20 separately
  remains the arbitrary-weight `4 x 4` result.
- [x] Ming-Xi Wang, 2011 ETH thesis, *Rational Points and Transcendental Points*,
  §§2.3 and 3.1 — constructs normalized Chebyshev--Blaschke products in every degree,
  proves their interval-preimage and nesting properties, and derives the Jacobi-`cd`
  multiplication formula by descent of elliptic isogenies.  Applied at the elliptic
  Lobatto points, this gives the exact alternation identity used in L116.  It does not
  construct a contraction similarity for the Crabb matrix.
- [x] Ng--Tsang, JCAM 277 (2015), *Chebyshev--Blaschke products: Solutions to certain
  approximation problems and differential equations*, DOI `10.1016/j.cam.2014.08.028`
  — proves the degree-`n` Chebyshev--Blaschke product has least deviation from zero on
  the corresponding elliptic interval.  Müller-Hermes--Szehr,
  arXiv:1405.4031, Lemma 2, records the exact value `sqrt(k(q^n))` and its
  model-operator role.  Combined with Wang's nodal formula and the DCT reversal,
  this proves L116's exact lower bound
  `t_*(phi(C_p+cC_p*))>=k(c^(2p−2))/c^(p−1)`.  These sources do not supply
  the matching all-size upper similarity.
- [x] Campaign L117 literature re-audit — searches for Li's exact modulus formula, elliptic
  Crabb similarities, and post-thesis citations found the thesis itself, Li's 2020 extremal-
  function paper, the unrelated 2023 KLS numerical-bounds paper, and the known `3 x 3`
  Glader--Kurula--Lindström theorem, but no all-size proof of the thesis's missing identities
  or of the equivalent periodized-sech metric.  O'Loughlin--Virtanen (LAA 697, 2024,
  DOI `10.1016/j.laa.2023.12.008`) cites the thesis but develops different truncated-Toeplitz
  closure classes.  This supports only the narrow phrase "apparently new closure of Li's
  fixed-weight all-size gap"; it is not a comprehensive priority claim.
- [ ] de Vries thesis (SV24 companion) — extremal-pair machinery details
- [ ] MMOR 2024 (C_N < 1+√2) — equality-case analysis of C–P (relevant to P2 stability)
- [ ] SV24 full §§1–5 (their Prop 2.4, 2.8, 3.5 used in §6; Berger–Stampfli generalization)
- [ ] R–S Question 4.1 citation trail — is unitality known to be insufficient? (still undetermined)
- [ ] BCG+23 arXiv:2312.04537 (Crouzeix, compressions of shifts, nilpotent classes)

## Internal load-bearing audit (2026-07-21)
- L17 and L21 were independently re-derived after the rapid Epoch-6 proof sequence.  L17's
  disconjugacy/positivity step is the exact ground-state transform with positive solution
  `1/sqrt(g')`; the comparison direction is `E_Q≤E_−1` for `Q≥−1`.  L21 has a genuine Slater
  point after scaling the Lyapunov series, and parity reduction must average the linear dual
  triple before re-optimizing its invariant ray.  No gap was found.  See
  `proof/load_bearing_audit.md`.
