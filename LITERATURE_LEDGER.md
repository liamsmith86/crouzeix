# LITERATURE_LEDGER.md

Papers, theorems, assumptions, possible gaps. **Latest primary-source refresh:
2026-09-05**, in [the versioned refresh report](proof/literature_refresh_20260905.md).
The historical notes below were last checked in July unless stated otherwise.

## September update — scalar proof verification

Lorist–Schwenninger, [arXiv:2608.03841v2](https://arxiv.org/html/2608.03841v2)
(August 17), proves the scalar constant two with an all-powers dilation lemma.
Our finite-matrix reconstruction is `proof/CROUZEIX_PROOF.md`, with an independent
proof audit in `proof/audit_global_proof_20260905.md`. Credit belongs to the
source authors. Jin's separate July/August proof and the September abstract
extension are documented in the refresh report; neither is needed by our
reconstruction.

The old universal scalar baseline `1+sqrt(2)` is superseded. The unrestricted
**complete** inequality remains a different problem. The new complete result
through dimension three is a preprint claim recorded in the refresh, not a
load-bearing ingredient or an independently audited result here.

Historical future-work labels below describe that earlier campaign. In
particular, old one-step abstract obstructions do not forbid the new argument,
which retains the identities for every power and their commuting errors.

## Core chain
| Ref | Result | Method | Notes / gaps to exploit |
|---|---|---|---|
| Crouzeix 2004 ("Bounds for analytical functions of matrices") | Conjecture stated; proved for 2×2 (constant 2, sharp) | Direct estimates, conformal maps | 2×2 case: W(A) is an ellipse; sharp via Jordan block [[0,2],[0,0]] |
| Crouzeix 2007, *Numerical range and functional calculus in Hilbert space*, JFA 244 | Universal constant 11.08; Theorem 2.1 gives an extremal `B∘phi` with `B` a finite Blaschke product of degree at most `n−1` | Cauchy integral + conformal splitting; finite Schur--Pick extremality | The numerical constant is superseded, but the finite-extremal reduction is load-bearing for L329 and is also restated in BGG+20 and Li 2020/2021 |
| Crouzeix–Palencia 2017 (SIMAX/FoCM) | W(A) is a (1+√2)-spectral set | f(A)+g(A)* via positive double-layer kernel; abstract lemma | Historical one-step baseline, superseded for scalar polynomials by the all-powers proof; g is the Cauchy transform of f̄ |
| Ransford–Schwenninger 2018 (SIMAX 1708.08633) | The abstract lemma (‖f(A)+g(A)*‖≤2‖f‖, g=C(f̄)) cannot yield better than 1+√2 in abstract setting | Explicit abstract construction | KEY OBSTRUCTION: need concrete structure of the pair (f(A), g(A)) beyond the abstract hypotheses |
| Malman–Mashreghi–O'Loughlin–Ransford 2024 | For each fixed N: C_N < 1+√2 | Compactness + strict inequality analysis of C–P equality conditions | Non-constructive; no uniform bound. Their equality analysis of C–P may reveal extremal structure |
| MMOR 2025 | Configuration constants via Neumann–Poincaré operator; domain-dependent bounds | Double-layer potential spectral analysis | Corrected 2026-09-05: for every compact convex domain with nonempty interior, analytic configuration constant a(Ω)<1 and scalar bound 1+sqrt(1+a(Ω)); no smoothness restriction. See refresh report |
| Schwenninger–de Vries, arXiv:2409.15954 "DLP for spectral constants revisited" | §6: ρ(f₀,x₀) := ∫Re(K_Ω(f₀)*f₀)dμ₀ = ReC; **K² + ρ ≤ 2K** (Thm 6.1: K ≤ 1+√(1−ρ)); ρ ≥ −1 (convex); ρ ≥ 0 ⟹ Crouzeix; disk ⟹ (K−1)ρ = 0; \|ρ\| ≤ a(Ω) (Rem 6.2) | Extremal pair + DLP | READ & VERIFIED (pp. 15–17). Closest baseline to campaign's P1; our increment = β-subtraction (q < 2). Also: Appendix A smooth-approximation lemma (useful for Ω ↓ W(A) rigor); Prop A.2 spectral-constant equivalence |
| arXiv:2410.10678 | "The algebraic numerical range as a spectral set in Banach algebras" | — | CORRECTED ID (audit): earlier ledger entry misattributed this ID to a parameterized C–P extension |
| O'Loughlin–Rani, arXiv:2603.15536 (2026) "q-Numerical Ranges and Spectral Sets" | Bound K = t + √(t²+a(Ω)), t = 1+γ(1)/2 | q-numerical range | CAUTION (audit): sentence after Thm 3.3 claiming K<2-type conclusion may be typo/gap — verify before use |
| Hartz–McCarthy, arXiv:2606.02922 (2026) "From Clouâtre–Ostermann–Ransford to Okubo–Ando" | For a continuous unital operator-algebra homomorphism `theta` and every scalar functional `beta`, `||theta^(n)|| <= max(1,||(theta+beta I)^(n)||)` at every matrix level; recovers the sharp Okubo–Ando similarity theorem | Potapov–Möbius variational argument at matrix levels | AUDITED FOR L21: when `||theta||cb>1`, their inequality plus `beta=0` gives `inf_beta ||theta+beta I||cb=||theta||cb`, so scalar-shift search is an exact restatement, not a relaxation. Crouzeix–Palencia supplies an operator-valued correction `R(f)=theta(alpha(f))*`, not the required scalar range except in the disk case. A resolution/cross-solver SDP probe also numerically excludes every positive-state scalarization on one dense 3×3 test (`proof/general_similarity_probe.md`). Do not claim a general CP-to-HM bridge. |
| Lewis–Overton 2020, "Partial Smoothness of the Numerical Radius at Matrices Whose Fields of Values Are Disks" | Under the simple polynomial boundary-eigenvector hypothesis, disk matrices form an analytic manifold of codimension `2n`; gives Crouzeix's global disk parametrization | Variational/partly-smooth numerical-radius geometry | PRIMARY SOURCE CHECKED.  This supports the ambient disk-manifold chart used in L115/L122 but does not classify which disk matrices attain scalar Crouzeix constant two or similarity square four.  Targeted equality searches still found no counterpart of L187's residual-zero submanifold. |

## Special cases proved (constant 2)
- Normal matrices (trivial: R ≤ 1; equality-2 impossible). von Neumann: W(A) ⊇ spectrum.
- All 2×2 matrices (Crouzeix 2004).
- n×n with W(A) a disk (Okubo--Ando: numerical radius at most one implies
  similarity to a contraction with condition at most two, so the disk is a complete
  `2`-spectral set).  Campaign L123 is not a new proof of this upper bound: it identifies
  a concrete all-size phase-palindromic stratified family on which the upper constant is attained.
  Campaign L187 substantially enlarges this to the analytic residual-zero manifold in L122's
  complete general-H disk chart, of physical dimension `2L−2`, with an explicit rank-one metric,
  matching characteristic-polynomial Blaschke product, and full normalized ambient stationarity.
  Targeted searches for disk-matrix equality classifications and Blaschke extremals found no
  matching parametrized equality theorem; treat the classification as apparently new pending
  publication-level review. `proof/crabb_full_disk_leading_residual_tube.md`.
- Nearly Jordan blocks (Choi–Greenbaum).
- Weighted shift matrices (Choi; complete 2-spectral via arXiv:2508.12768 Aug 2025).
- Certain 3×3: tridiagonal with elliptic W(A) centered at an eigenvalue. The former “KLS matrices (2025–26)” entry was unsupported: the 2023 source studies KMS matrices and its conformal estimates were numerical. A September 2026 preprint gives exact containment certification; see the refresh report for scope and audit status.
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
- **Campaign L342 (2026-07-26):** at every finite nondegenerate
  Gau--Wu equality model, the second-order Paulsen/L21 similarity SDP
  reduces to one universal boundary-row quadratic and equals four
  times the optimized scalar Crouzeix Hessian.  This is an exact
  local osculation statement, not the global assertion `t_*<=4` and
  not a proof of either scalar or completely bounded Crouzeix in a
  neighbourhood.  Okubo--Ando supplies the base disk bound and
  Paulsen supplies the similarity equivalence, but neither source
  currently identified in the ledger supplies this second-variation
  identity.  Treat its novelty as unclassified pending a targeted
  publication audit; do not advertise priority from the current
  search record alone.  `proof/gau_wu_similarity_osculation.md`.
- **Campaign L343 (2026-07-26):** the polynomial support frame and
  general-H disk chart are concrete forms of the classical
  Crouzeix/Lewis--Overton disk-matrix parametrization, while
  Okubo--Ando supplies the complete-2 bound on the chart.  The
  inverse-Toeplitz equality coordinates also overlap the classical
  Toeplitz-moment/Levinson substrate already recorded for L187.
  Therefore no novelty is claimed for those ingredients.  The
  campaign-specific contribution is their exact recentering at every
  finite Gau--Wu equality model and the identification of the
  common L342 Hessian's strict disk block with the kernel of the
  `2n−2` non-affine weighted-support quotient.  No matching
  second-variation quotient theorem has been identified, but its
  novelty remains unclassified pending a targeted publication audit.
  The quotient sign itself is still open.
  `proof/gau_wu_disk_chart_recenter.md`.
- **Campaign L344 (2026-07-26):** equality of two support-factor
  Grams, the polar-decomposition range isometry, and replacement of
  a minimal Moore--Penrose solve by a fixed spectral factor are
  classical operator/systems arguments.  No novelty is claimed for
  those ingredients.  Poon--Spitkovsky--Woerdeman, *Factorization of
  Singular Matrix Polynomials and Matrices with Circular Higher Rank
  Numerical Ranges* (SIMAX 43 (2022), DOI
  `10.1137/22M1475934`), is a particularly close modern source for
  the singular Hermitian trigonometric-polynomial factorization
  substrate.  Their use to transport the Gau--Wu
  second-support reserve into L343's exact inverse-Toeplitz
  coefficient chart is campaign-specific bookkeeping.  It remains a
  reduction only: the final endpoint-output domination on the
  `2n−2` shape quotient has not been proved.
  `proof/gau_wu_toeplitz_port_transport.md`.
- **Campaign A292 / Grunsky caution (2026-07-26):** the exact
  conformal-Fourier chart is elementary finite Toeplitz positivity,
  while the observed phase-covariant Hermitian form resembles the
  classical complex-symmetric contraction forms arising from
  Grunsky coefficients.  Classical Grunsky inequalities say that
  the weighted Grunsky matrix of a univalent map is contractive.
  The campaign has not identified its Gau--Wu Hermitian matrix with
  a Grunsky matrix, and the finite Blaschke extremal is not itself a
  univalent map in degree greater than one.  Treat this only as a
  search cue for a Hardy/Grunsky factor, not as a citation that
  proves phase covariance or the sign.
  `proof/gau_wu_conformal_shape_phase.md`.
- Compressions of the shift / model space operators with certain Blaschke conditions (Bickel–Gorkin school).
- Matrices attaining ‖A‖ = numerical radius conditions (Crabb-type equality cases).

## Numerical evidence
- Greenbaum–Overton (2017-2018, "Numerical investigation of Crouzeix's conjecture"): nonsmooth optimization (BFGS/GRANSO) over A and p; global max of R appears to be exactly 2, attained on Crabb-type / Jordan-block configurations. No R > 2 ever found.
- AIM workshop reports (aimath.org/pastworkshops/crouzeixrep.pdf): consensus directions.

## Structured-matrix novelty audit for L187/L193/L218 (2026-07-24)

Classification here is about the **mathematical core**, not whether the
campaign independently re-proved it.  This was a targeted source audit,
not an exhaustive priority search.

Revalidation after L248 checked the listed DOI metadata against
Crossref and the primary publisher records; the three novelty
classifications below are unchanged.

- **Classical scalar substrate.**  Gohberg--Semencul (1972), *The
  inversion of finite Toeplitz matrices and their continual analogues*,
  Mat. Issled. 7(2), 201--223, gives the classical inverse-Toeplitz
  formula.  Kailath--Vieira--Morf, SIAM Review 20 (1978), 106--119,
  DOI `10.1137/1020006`, reviews the equivalence with
  innovations/orthogonal-polynomial structure.  Kailath--Kung--Morf,
  JMAA 68 (1979), 395--407, DOI
  `10.1016/0022-247X(79)90124-0`, places Toeplitz inverses and related
  matrices in the classical low-displacement-rank framework.
- **Classical block substrate.**  Block Toeplitz inversion predates the
  campaign (for example Akaike 1973).  Labahn--Choi--Cabay, SIAM J.
  Comput. 19 (1990), 98--123, DOI `10.1137/0219006`, explicitly gives
  block Hankel/Toeplitz inverse formulae generalizing the scalar
  Gohberg--Semencul formula.  Thus noncommuting matrix blocks do not by
  themselves establish novelty.
- **Classical Schur/Levinson substrate.**  Delsarte--Genin--Kamp,
  SIAM J. Appl. Math. 36 (1979), 34--46, DOI `10.1137/0136004`,
  gives the matrix Schur--Szegő parameterization of positive
  block-Toeplitz systems.  Wiggins--Robinson, JGR 70 (1965),
  1885--1891, DOI `10.1029/JZ070i008p01885`, is an early
  multichannel/block Levinson recursion.  Ackner--Lev-Ari--Kailath,
  SIMAX 15 (1994), 140--150, DOI
  `10.1137/S0895479891218922`, explicitly treats the Schur algorithm
  for matrix-valued functions.

Lemma-by-lemma scope:

1. **L187 = class (ii), known object in new coordinates.**  The cone
   `H^{-1}` Hermitian Toeplitz and its inverse/displacement structure
   are classical.  What is not supplied by those sources is the
   campaign's exact equivalence
   `Hardy residual = 0 iff H^{-1} is Toeplitz`, nor its identification
   with the local Crouzeix/similarity equality manifold and endpoint
   certificate.  Retain only that problem-specific connection as
   apparently new pending publication-level review.
2. **L193 = class (ii), known block object in new coordinates.**  The
   block Toeplitz cone, its inverses, and noncommuting block formulas
   are classical.  The problem-specific Hardy-residual equivalence,
   complete-`2` upper endpoint, and Crouzeix equality-manifold
   interpretation are the campaign contribution.  Do not advertise
   “arbitrary noncommuting block Toeplitz coefficients” alone as new.
3. **L218 = class (i) for its core.**  Ordered matrix Schur parameters,
   termination for finite rational inner functions, and their
   Toeplitz/Levinson role belong to classical matrix
   Schur--Potapov/Levinson theory.  The forced-zero Crabb
   specialization, dimension match to L193, and displayed reflected
   monomial tangent are useful campaign coordinates, but are not a new
   matrix Schur algorithm.  Novelty, if claimed, must be restricted to
   their later use in the Crouzeix endpoint construction.

## Numerical-contraction equality audit for L326 (2026-07-26)

- **Primary source:** H.-L. Gau and P. Y. Wu, *Inner functions of
  numerical contractions*, Linear Algebra Appl. **430** (2009),
  2182--2191, DOI `10.1016/j.laa.2008.11.020`.
- **Theorem 8 is exactly the missing terminal structure.**  If
  `w(T)<=1`, `T` has no singular unitary part,
  `||f||_infinity<=1`, and `||f(T)x||=2`, then `f` is inner,
  `f(0)=0`, and `T` has the orthogonal reducing summand
  `X_phi S(phi) X_phi^−1`, `phi=zf`, with endpoint scaling
  `sqrt(2) direct-sum I direct-sum 1/sqrt(2)`.
- **Theorem 5 identifies that summand.**  It is cyclic and
  irreducible, has no unitary part, has numerical range the disk, and
  attains two for `f`.  Thus the reducing scalar model itself and the
  equality mechanism are classical; do not advertise L326 as a new
  equality classification.
- **Campaign-only bookkeeping.**  Near `C_p tensor I_m`, compactness
  and irreducibility imply every limiting reducing projection is
  `I_p tensor q`.  Hence the Gau--Wu summand has dimension `p r` and
  its complement lowers repeated-block copy multiplicity.  L328 then
  uses Gau--Wu Corollary 3's explicit matrix and an elementary power
  invariant to prove `r=1` locally.  These specialize the classical
  theorem to L325's partial canonical endpoint; they are local
  applications, not a new model theorem.
- **Scope correction.**  General L193 disk points need not be
  nilpotent.  Therefore Gau--Wu's inner `f` cannot generally be
  replaced by a monomial.  Monomial/Crabb rigidity is valid only on
  the explicit nilpotent subcase.  `proof/repeated_crabb_scalar_model_split.md`.
- **Post-L329 coverage correction (L330).**  Corollary 3 already
  gives explicit nonnilpotent sharp disk matrices outside the
  repeated-Crabb stratum.  For zeros `0,a`, the specialized
  `3 x 3` matrix `G_a` has support eigenvalues
  `1,-1,a cos(theta)` and satisfies
  `f_a(G_a)=2E_13`, but its spectrum `{0,0,a}` rules out affine
  Crabb equivalence when `a!=0`.  Thus any globalization must first
  cover arbitrary Gau--Wu models; L192/L329 cover only their local
  monomial collision sector.  This correction is a direct
  application of the classical source, not a novelty claim.
- **L331 stationarity scope.**  Gau--Wu prove the exact disk model
  and equality but do not state the ambient numerical-range/Riemann
  first-variation identity used here.  L331's Poisson-weighted
  support-projection residue cancellation appears campaign-specific
  pending a broader variational-literature audit.  It is only
  stationarity and must not be advertised as a local inequality.
- **L332--L333 second-order/local scope.**  Targeted searches checked
  Greenbaum--Lewis--Overton, *Variational Analysis of the Crouzeix
  Ratio* (Math. Programming 164, 2017), Lewis--Overton, *Partial
  Smoothness of the Numerical Radius at Matrices Whose Fields of
  Values are Disks* (SIAM J. Matrix Anal. Appl. 41, 2020), Li,
  *On the Uniqueness of Functions that Maximize the Crouzeix Ratio*
  (arXiv:2002.01027, revised 2026), and Overton, *Local Minimizers of
  the Crouzeix Ratio* (arXiv:2105.14176).  The first-order paper
  establishes nonsmooth stationarity at classical candidate pairs;
  Lewis--Overton identify the centered `3 x 3` disk-matrix manifold
  and its dimension; Li proves uniqueness in selected classes and
  gives nonunique `3 x 3` elliptic examples; Overton gives numerical
  local-minimizer evidence.  None states L332's optimized
  support/Riemann/Blaschke five-normal Hessian or L333's scalar
  constant-two operator neighbourhood at the nonnilpotent Gau--Wu
  family.  Gau--Wu equality structure and Crouzeix's finite-extremal
  reduction remain classical inputs.  The Hessian and local theorem
  therefore appear new pending publication-level audit; do not
  promote this to a global `3 x 3` or cb theorem.
- **L334 arbitrary-degree scope.**  L334 is a campaign numerical
  falsification protocol, not a literature theorem or a novelty
  claim.  Its general support/Riemann/zero-motion jet finds negative
  complete normal Hessians in dimensions `4..8`; no checked source
  supplies the proposed arbitrary-degree defect-one negative Gram.
  That identity remains unproved and must not inherit L332's exact
  status from numerical agreement.
- **L335 support-Gram scope.**  The ingredients—compressed-shift
  model calculus, Blaschke logarithmic derivatives, and Schwarz
  transforms—are classical.  Their assembly as the exact
  second-support contribution to this optimized Crouzeix Hessian is
  campaign-specific, but no novelty claim is made before a targeted
  publication audit.  L335 signs only this component; the Hardy and
  zero-motion remainder remains open.

## Fixed repeated-Crabb neighbourhood audit for L329 (2026-07-26)

- **Finite extremals are classical.**  Crouzeix 2007, Theorem 2.1,
  proves that an extremal for an `n x n` matrix can be written
  `B∘phi`, with `B` a finite Blaschke product of degree at most
  `n−1`.  Bickel--Gorkin--Greenbaum--Li--Overton--Ransford--
  Schwenninger--Wegert, *Crouzeix's Conjecture and Related Problems*
  (arXiv:2006.04901), and Kenan Li, arXiv:2002.01027, explicitly
  restate this fact.  L329 uses it only for compact finite norming
  data and makes no novelty claim for that reduction.
- **Closest local literature is weaker or differently scoped.**
  Lewis--Overton prove partial smoothness of the disk-matrix
  manifold; Greenbaum--Overton study nonsmooth stationarity/local
  minimizers; BGG+20 Theorem 3.8 gives a neighbourhood of one Crabb
  matrix in which extremal Blaschke products have maximal degree.
  None of these sources states the Crouzeix inequality on a full
  neighbourhood of `C_p tensor I_m`, and the BGG+20 neighbourhood
  conclusion is about extremal degree, not the constant-two bound.
- **Targeted search outcome.**  Searches for “repeated Crabb,”
  “Crabb block neighbourhood,” and local Crouzeix theorems found no
  matching repeated-multiplicity result.  The local scalar theorem
  L329 therefore appears new, but this is not a priority proof and
  must remain qualified as pending publication-level audit.
- **Scope guard.**  L329 is fixed `(p,m)`, local, and scalar.  It
  neither proves a radius uniform in size nor the completely bounded
  conjecture, and it does not classify arbitrary global sharp
  sequences.  L330 additionally shows that it does not cover all
  already known Gau--Wu equality points away from the monomial
  collision.

## Known equality structure (R = 2 attained/approached)
- A = [[0,2],[0,0]], p = z. R = 2 exactly. W(A) = disk radius 1.
- Crabb matrix family (nilpotent Jordan-like with specific superdiagonal weights √2,1,...,1,√2? — verify): ‖p(A)‖ = 2 with p = z^k, W(A) = unit disk. Attains 2 for higher powers.
- Gau--Wu 2009 shows more generally that scalar equality for a numerical
  contraction splits a disk-model block `X_phi S(zf)X_phi^−1` with
  inner `f(0)=0`; equality is not restricted to nilpotent monomials.
  The Crabb family is its monomial specialization.
- **L336 novelty scope.**  Its ingredients are classical:
  simple-singular-value perturbation, the model-space kernel identity,
  and the right/left finite-Blaschke Stein identities.  The targeted
  search found Crouzeix's numerical observation that sharp disk
  matrices appear to be local maxima and Lewis--Overton's
  partial-smoothness geometry of disk matrices, but no source proving
  the coupled second-order endpoint-flux identity or its sign at
  arbitrary Gau--Wu models.  Treat the identity's assembly as
  campaign-specific and the still-open positivity as having no
  novelty claim until a publication-level audit.
- **L337 novelty scope.**  The canonical conjugation on a scalar
  model space and the formula `h(S_phi)1=P_(K_phi)h` are classical.
  L337 is their direct application to the two endpoint defects, so
  record it as a classical identity in campaign coordinates, not an
  original theorem.  Its role is to expose the exact Hilbert metric
  for the still-open mixed L336 calculation.
- **L338 novelty scope.**  The tangent parametrization
  `h=(zf)k−J_fk`, Szegő kernel Gram, and real Hilbert-space square
  completion use classical model-space ingredients.  Their assembly
  eliminates the campaign's moving Blaschke zeros, but it is best
  classified as a new coordinate reduction rather than a new
  theorem of independent scope.  No novelty claim is made for the
  still-open physical Riesz-energy inequality.
- **L339 novelty scope.**  Ando's numerical-contraction
  factorization, the unitary completion of a defect-one partial
  isometry, and boundary model-kernel norm
  `||k_zeta^f||²=D_f(zeta)` are classical.  L339's identification of
  L335's exact Hessian reserve with that support-port energy is a
  campaign-specific bridge.  Treat it as a new reduction assembled
  from classical identities, with no claim that the still-open port
  contraction is known or original.
- **L340 novelty scope.**  L340 is only a numerical structural
  diagnostic for the still-open physical loss: PSD to roundoff with
  rank `6n−14` and nullity `(n−4)²` in twelve models.  It is not a
  theorem and carries no novelty claim.  If an exact response
  factorization is later proved, its relationship to classical
  unitary-colligation, Stein, and model-space realization theory must
  be audited before making any originality statement.
- **L341 novelty scope.**  The primal/dual transferred-frame Stein
  equations are direct applications of L127's classical
  model-kernel identity, and the moving-kernel square is elementary
  finite-dimensional perturbation theory.  The identification of
  the L340 loss kernel with
  `(s_E,Y_1^opt q,(Y_1^opt)*p)` is presently numerical, not a theorem,
  and carries no novelty claim.  Any eventual response metric must
  be compared with classical conservative-realization and
  observability-Gramian theory before originality is assessed.

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
- [x] Y. Arlinskii, arXiv:1109.4020, *Schur parameters, Toeplitz matrices, and Krein
  shorted operators* -- develops the classical conservative-system/Schur-parameter
  framework and expresses shorted Toeplitz defect operators through Schur-parameter
  products.  Campaign L267's unitary Redheffer initial/final defect congruences belong to
  this classical lossless-feedback circle and carry no novelty claim.  The campaign-specific
  open question is narrower: whether the L243/L251/L258 repeated-Crabb first reflected
  coefficient admits one copy-scalar lossless feedback realization through associated order.
  Neither the cited shorted-defect theory nor the abstract L267 identity establishes that
  physical factorization.
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
