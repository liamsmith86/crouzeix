# RESEARCH_STATE.md — Crouzeix Conjecture Campaign

**Last updated:** 2026-07-22 (Epoch 6 — L56/L58 compact transfer frontier)

## NEWEST (2026-07-22): bridge closed; only the positive compact tail remains
- **L56 closes the former `.005→.01` gap.** Four shared-nome boxes certify both
  determinant charts and all three final-minor charts for both signs. Together with L55 and
  L54, (RT) is now proved continuously on `0≤c≤.020736`.
- **L57 advances the compact interval rigorously, not by fixed-nome sampling.** The negative
  sign passes through `.63`, overlapping L42 because `.63>2^(−2/3)`. The positive sign passes
  through the exact rational `c+=.5566585294072849…`. Therefore the sole missing interval for
  the full elliptic `4×4` slice is the positive-sign tail
  `c+<c<2^(−2/3)≈.6299605`.
- **L58 identifies why raw subdivision stalls there.** At `a=b=1`, each positive determinant
  chart is `81c^9` times explicit squares and a fourth power. In the inward deficits
  `A=1−a,B=1−b`, the two linear coefficients coincide and are nonnegative; their only
  non-square factor is `2+3k+3Zk−8Zk²≥0`. Exact regeneration checks all 197,563 records in
  each chart. Higher total deficit degree has exact collected l1 norm `41,235,531,913`, but
  that uncorrelated bound is too coarse. The next certificate should keep those higher
  coefficients correlated in the nome and the three surviving projective variables, rather
  than add a tolerance or continue raw binary64 subdivision.

## AUDIT-GATE OUTCOME (2026-07-21, after L48 and external steering review)
- **L51 compresses the live interior theorem to two projective polynomial
  inequalities.**  Strict positivity of the first diagonal block leaves one leading `3×3`
  minor and the full determinant.  In `P=p²,t=1−o`, their exact orders at the zero-node/
  orientation intersection are 3 and 4, so two largest-coordinate charts remove that
  singularity.  Full cubic-envelope Bernstein tests pass at every fixed `c` tested from
  `.001` to `.629`; this is strong evidence, not yet a continuous-`c` certificate.  Directed
  Taylor boxes locate further boundary intersections rather than a negative determinant.
  L52 then closes both extreme-orientation faces by an elementary scalar factorization and
  proves that the apparent square at `P=1,o=0` has no moving interior zero.  Its factor is
  strictly negative for every nondegenerate nome.  `proof/slice_core_projective_reduction.md`.
- **L53 resolves the nested minor degeneracies exactly:** after the main order-3 chart, the
  intersections `(v,1−|b|)=(0,0)` and then `(h,1−a)=(0,0)` each have exact order one.
  A new deterministic generator rebuilds all 49,448/197,563 envelope records and all chart
  orders from the original `4×4` core.  The earlier continuous-nome prototype also had a real
  implementation bug: it allowed the zeroth Taylor term to vanish.  Correcting that and
  keeping the Taylor coordinate correlated led to L54: an 80-digit Arb plus directed
  Bernstein proof of the full transfer theorem on `c∈[.01,.020736]`, all final charts and both
  signs, in four rational boxes. L55 normalizes the low-nome algebra (determinant order `c^9`,
  minor order `c^2`) and now factors every exceptional face exposed by the remainder charts.
  Exact rational terms through order 15 plus parity-aware Arb tails close det0's
  main-orientation chart and, for the positive sign on `[0,.01]`, four widened charts at its
  `a=0` corner. Det1 has two more exact blow-ups ending in the positive quadratic (23); two
  59-chart directed runs certify its complete main chart for both signs on `[0,.005]`
  (84 minutes positive, 76 minutes negative). Det0's former `U`-axis line now ends in two
  coefficient-positive exact forms; the complete 26-chart hierarchy passes at the root for
  both signs on `[0,.005]` (the negative regeneration took 3033 seconds, 11.02 GB, no swap).
  The new exact minor corner and midpoint forms (28)--(33), including the positive secondary
  `a=2c` ridge and tertiary-one ratio-zero face, now support complete local-plus-global
  certificates for all three final minors and both signs on `[0,.005]`. Positive global leaf
  counts are `81/32/184`; negative counts are `71/32/188`. Independent full regeneration took
  10m19s and 3m52s respectively. Scale-free nome/main arm charts now remove the det0
  Cartesian cutoff. Exact extraction confirms that the apparent normalized `S=1/2` feature was
  only a centered-model seam: (25) is uniformly transverse there. The four negative widened
  charts pass at the root; the complete positive main-dominant chart and a cap-`1/2` negative
  equality tube certify; and the global complements close in `235/399` leaves at depths `8/9`.
  Consequently L55 proves the complete polynomial core, every chart and both signs, on
  `[0,.005]`. The later L56 bridge connects it to L54, while L57 supplies the current compact
  frontiers stated above.
- **L49 finds the sharp nome estimate required by the zero-node ridge:**
  `k≤4c/(1+4c²)` for `c≤1/2`.  Two Jacobi-product factors suffice, and the
  remaining degree-23 polynomial has 24 positive exact Bernstein coefficients.  On `p=0`,
  L50 reduces BE to one scalar inequality and closes the complete face.  The high half uses
  125 exact Bernstein coefficients; the low half uses a finite ridge-centered/blow-up chart
  cover around the exact equality mechanism `b=1,a=2c`.  No floating-point sign decisions
  enter either certificate.
- **L47 proves the new sharp block tradeoff `||B|| ||C||≤2`.** The two matrix
  traces and determinants compress to one orientation scalar; a universal 2×2 singular-value
  majorant reduces the claim to a four-variable rational polynomial. Elementary nome bounds,
  four singular-corner blow-up charts, and exact integer Bernstein coefficients close the whole
  parameter box. The checker regenerates every coefficient in about 13 seconds. Consequently
  L48 closes both full transfer axes `a=0` and `b=0`; the remaining rank-one/rank-one obstruction
  is genuinely two-parameter with `ab≠0`. `proof/slice_coupled_defects.md` §6.3.
- **L17 and L21 survive independent re-derivation.** L17's load-bearing step is the exact
  ground-state identity with positive solution `1/sqrt(g')`; L21 has a strict Lyapunov-series
  Slater point, and its parity restriction is valid only by averaging the linear dual triple and
  then rescaling the invariant ray. `proof/load_bearing_audit.md` records the full sign audit.
- **The similarity route survives its first general-matrix falsification gate.** A seeded sweep
  of 120 varied complex matrices (`n=3..8`) accepted 102 through independent map, Cauchy,
  double-layer, resolution, and SDP primal/dual gates. None exceeded four; the maximum was
  `t*=3.8399296`. The sharp case rose from `3.68956` to `3.97783` as the outer offset shrank from
  `.02` to `.00125`, always from below. Seventeen nearly normal cases had unresolved polygonal
  map discretization and one SDP had a bad duality gap, so they were rejected rather than counted.
  This makes L21 a credible general attack, not a proof. `proof/general_similarity_probe.md`.
- **L44's exact sharp boundary passes analytically and at 80 digits.** At `p=a=b=0` its energy is
  `k(1+c²)/(4c)`, independently of the modal angle, and L23 bounds it by `1/(1+c²)<1`. A scaled
  80-decimal scan down to `c=10^-20` finds margin `~3c²` and no excess. The old
  `p=p*+c⁴x` coordinate belongs to L29's different ridge and must not be imported into L44.
  This proves only the exact boundary (L46); a finite neighbourhood scan is still not L44.
- **All three L29 certificate scripts regenerated cleanly and independently.** The centered
  completion rebuilt its 309,479- and 565,425-term expansions and passed after 85 minutes; the
  small-edge and compact-range checkers also exited exactly with no unresolved boxes. The
  detached job ended with `PASS 2026-07-21T15:10:52-07:00`. The L29 audit gate is complete.
- **L45's one-variable convex branch is genuinely present.** In 2,000,000 random quadratic-form
  probes on exact modal slice matrices, 772,130 had `q2>0` with the minimizing `a` inside
  `(-1,1)`; none had a negative discriminant (smallest sampled margin `5.48e-5`). Thus the new
  identity does not collapse merely by concavity. The discriminant route remains well supported,
  but it must retain the conformal coupling.
- The exact frontier is still (RT), but L47/L48 remove two complete one-parameter sections of
  its square and provide a new coupled invariant likely useful in the interior.

## NEWEST (2026-07-21, Epoch 6) — L20 reduced to an explicit trace-cone inequality
- **L21 PROVED (dimension-independent):** for every strictly stable matrix `T`, the least
  similarity-square in `I≤P≤tI`, `T*PT≤P` is exactly
  `max(1, sup_{Z≥0} tr(Z−TZT*)_-/tr(Z−TZT*)_+)`. The proof is an explicit Slater/SDP-dual
  calculation followed by positive-part minimization; it is not a numerical inference.
- **L22 PROVED (slice-specific):** diagonal symmetrization plus the SVD of the 2×2 bidiagonal
  block puts every elliptic-slice `T=φ(A)` in a three-real-parameter form `(c,r,u)`, with
  `tan(v)=r tan(u)` and the two conformal nodes explicit in Jacobi `sn`. The metric problem is
  exactly four coupled 2×2 modal LMIs. Independent nodal reconstruction agrees to `8.9e-15`.
- Chiral symmetry permits both primal metrics and sharp dual certificates to be parity-block
  diagonal. Thus **L20 is now exactly** `tr(Z−TZT*)_- ≤ 4 tr(Z−TZT*)_+` for two coupled
  2×2 positive blocks and the explicit modal `T`; no optimization or phase classification remains
  in the statement.
- Numerical stress test only: all `15×9×9=1215` deterministic modal-grid cases pass; largest
  `t=3.999771308` at `(c,r,u)=(.001,.97,.03)`. Primal and dual values agree to `1.1e-9` on the
  default cases. The near-four singular corner shows that a proof must be sharp and uniform.
  Proof and reproducer: `proof/slice_similarity_duality.md`,
  `experiments/slice_similarity_duality.py`.
- Failed construction audit: diagonal metrics, short observability Gramians, and forcing one
  modal contraction inequality to equality all fail before the true optimum reaches four. The
  next attack should characterize extreme dual block pairs (low rank/boundary numerically) and
  prove their trace ratio directly from the conformal coupling of the two nodes.
- **Boundary/modal progress (L23–L26 PROVED):** Jacobi's product gives the sharp focus-node bound
  `k/c ≤ 4/(1+c²)²`. Using the reciprocal-quadratic dependence on `tan²u`, this proves
  `||C||≤2` for one entire off-diagonal modal block. On the singular `c→0` face, `T` becomes
  a nilpotent scalar weighted shift; every consecutive weight product is ≤2, so an explicit
  diagonal similarity proves `t*≤4`, sharply at the Crabb weights `(√2,1,√2)`.
  `proof/slice_boundary_theorems.md`; regression `experiments/slice_boundary_check.py`.
- **L26 closes the other modal block:** `||B||≤2` follows from an exact determinant factorization
  and a Möbius barrier for the normalized inverse ellipse map. Positive odd Taylor coefficients
  give a cubic minorant; three remaining scalar theta inequalities are proved by exact rational
  bounds near `c=0` and 13,500 outward-rounded algebraic interval boxes. Proof and certificate:
  `proof/slice_upper_block_theorem.md`, `experiments/slice_upper_block_certificate.py`.
- Both block norms ≤2 close every one-block dual phase, but do not handle the observed coupled
  phases where `t*>max(||B||²,||C||²)`. The trace route's sole slice obstruction is now genuinely
  coupled.
- **L27 eliminates the coupled contraction LMIs exactly.** If `Q_o,Q_e≥0` are their Stein
  defects and `K_ij=(1−τ_i²τ_j²)^{-1}`, then every modal contraction metric is
  `H_o=K∘(Q_o+cΣQ_eΣ)`, `H_e=K∘(Q_e+c^{-1}ΣQ_oΣ)`. Complementary slackness gives
  `rank Z_i+rank Q_i≤2`. Since both dual blocks cannot be full and one-block phases are already
  closed, any hypothetical KKT optimum above four lies on a rank-one/rank-one or rank-one/full
  coupled face. Hence `Q_o=aa^T`, `Q_e=bb^T` (one vector may vanish), leaving two directions and
  one relative scale. This is an exact reduction, not a numerical rank guess. Arbitrary choices
  of these defects can be badly conditioned; the live task is to exclude an *optimal KKT point*
  above four using `tan(v)=r tan(u)` and `r=H(p)`. `proof/slice_coupled_defects.md`.
- **L28 identifies the generic rank-one/full face.** A sign-separated KKT boundary is an
  orthogonal colligation, and its exact dual value is
  `sup_{a∈[-1,1]} ||B R_a(CB)||²`, where this matrix is the upper block of the odd Blaschke
  product `f_a(T)=T(T²−aI)(I−aT²)^{-1}`. The parity-swapped lower block is now proved ≤2 by
  nodal-value convexity plus a reciprocal-orientation estimate. For the upper block, orientation
  also disappears exactly, leaving one determinant `N(c,p,H(p),a)≥0`. The rigorous envelope
  `s₀p+a₃p³≤H(p)≤s₀p+(1−s₀)p³` contains the target and passes global searches; the unproved
  sharp ridge is `c→0, p→1/2, a∼3c`. `proof/slice_odd_block_reduction.md`.
- **L30 removes the Blaschke parameter from L29.** At reciprocal orientation the upper odd
  block has equal diagonal, so `||M||≤2` is exactly `det(M)+2|M₁₂+M₂₁|≤4`. In the inner-node
  value coordinate `t`, the two signed residuals are quadratics. The minus-sign quadratic is
  concave and is therefore closed by its `t=±1` endpoints (L26). The plus-sign quadratic can
  fail only on `A>0, |B|<2A`, where the complete remaining condition is `4AC−B²≥0` at the lower
  and upper cubic-envelope values of `r`. Thus both orientation and `t` are gone: L29 is now a
  two-variable `(c,p)` theta inequality with a sharply localized branch.
- **L31 factors the final discriminant exactly.** With
  `Q=2g(1+c²p)−d(kp+4c)` and `S=−2g(c²+p)+d(kp+4c)`, the remaining numerator is
  `16rg²p(1−c²)²(1−d²)−(Qr−S)²`, where
  `1−d²=(1−k²)(1−k²p⁴)/(1−k²p²)²`. This is an exact SymPy-audited identity, not a fit.
  It turns L29 into the sharp distortion inequality
  `|Qr−S|≤4g(1−c²)sqrt(rp(1−d²))`. The numerator is concave in `r`, so only the two
  cubic-envelope endpoints remain. The small-nome ridge is still the live analytic obstruction.
- **L32/L33 sharpen and localize the certificate.** A second exact factorization writes the
  discriminant as `4[rF₁F₂−g²D²]`, with
  `D=(r−p)−c²(1−pr)+λ(1−r)`; direct endpoint formulas for `p−r` avoid catastrophic
  cancellation. Separately, if `cℓ≥1`, convexity of the full nodal rectangle plus an exact
  opposite-sign vertex factorization proves every upper block is ≤2. L26's product bound gives
  `cℓ≥1` for `c≥12599/20000` via a five-factor exact polynomial audit. Thus L29's live square
  certificate is restricted to `0<c<12599/20000`, with the cancellation-free form preferred for
  intervals. L34 further factors `F₁` through a quadratic with exact vertex
  `p*=(g²−4c²)/(2g²(1−2c²g))=1/2+3c²/2+O(c⁶)`. This explains the sharp
  small-nome ridge and supplies the stable coordinate `p=p*+c⁴x`; the identity itself makes
  no unproved sign assumption about its residual theta expression `H`.
- **L29 is now PROVED: the generic rank-one/full KKT face is closed.** L35 first certifies the
  singular core. L37 deepens this to the entire centered band
  `|(p−p*)/c²|≤50` for `c≤1/20`, and also supplies the bridge tube
  `|p−p*|≤4c⁴` through `c=1/12`. Its exact checker regenerates 309,479- and
  565,425-term polynomials, using Bernstein tube bounds and a multiscale annulus estimate.
- **L38 closes every remaining small-nome branch point.** Exact branch numerators exclude
  `p<c/4` and `p≥3/4`; coefficient domination handles the boundary regions; exact Bernstein
  bounds handle the middle whenever `|2p−1|≥96c²`. The identity
  `(p*−1/2)/c²=(g³−2)/(g²(1−2c²g))` puts the complement inside L37. The worst
  regular correction ratio is `0.907895<1`.
- **L36/L39 finish the compact ranges.** L39 uses exact ridge-complement coordinates and
  outward-rounded Taylor intervals to close `1/20≤c≤1/12` outside L37's tube (15,267
  bisections, none unresolved). L36 supplies `1/12≤c≤12599/20000` (32,641
  bisections), and L33 supplies the rest. Thus L30–L32 prove the upper odd block ≤2 for every
  parameter. The lower odd block was already analytic, so L28's rank-one/full face is done.
  **Only L27's rank-one/rank-one coupled face remains for the L20 elliptic-slice theorem.**
- **L40 now reduces that last face to one matrix-valued inner theorem.** Two spectral
  row-Gram rotations give the transfer
  `R_{a,b}(T)=−S⁻¹(T−A)(I−AT)⁻¹S`; every rank-one/rank-one trace ratio is at most
  `||R_{a,b}(T)||²`. A cancellation-free block formula extends to the closed parameter square.
  In symmetric modal coordinates its node function is a `2×2` rational inner function with
  determinant `(ab−z²)/(1−abz²)`. Proving `||R_{a,b}(T)||≤2` uniformly would prove L20.
  The boundary reduces to the proved scalar even sector, but numerical maxima can occur in the
  genuinely matrix-valued interior, so a scalar-only argument is insufficient.
- **L41 removes the transfer resolvents.** The exact operator-ball defect identity makes
  `||R_{a,b}(T)||≤2` equivalent to a quartic polynomial `4×4` LMI. L24/L26 and two bilinear
  scalar estimates prove both `2×2` diagonal blocks positive semidefinite for every `a,b`.
  The sole content left is its explicit off-diagonal Schur inequality, with coupling
  `3[a(1−b²)B+b(1−a²)C*]`. This is the preferred form for the conformal-node attack.
- **L42 closes the full high-nome range `c≥2^(−2/3)`.** In symmetric modal coordinates
  the transfer is a contraction because its two node values come from a rational inner matrix.
  Conjugating to physical coordinates costs exactly `c^(−3/2)≤2`. The remaining transfer
  theorem is confined to `0<c<2^(−2/3)`.
- **L43 square-completes the remaining coupling.** The polynomial core is PSD exactly when
  an explicit four-block small-gain matrix has norm at most one. Its only resolvents are those
  of `BC,CB` at contracted parameters `u=3a/(4−a²)`, `v=3b/(4−b²)`. The stronger scalar
  target L44, `Σ||S_ij||²≤1`, implies the theorem and survives global searches even over the
  full rigorous cubic envelope for `r=H(p)`, approaching equality only on the known
  `c→0,a,b→0` boundary. The softer strip `s₀p≤r≤p` is false (value `1.3067`), so the cubic
  conformal input is essential. L44 remains numerical, not proved.
- **L46 certifies L44's exact small-nome boundary.** At `p=a=b=0`, direct modal substitution
  gives `Σ||S_ij||²=k(1+c²)/(4c)`, with complete cancellation of the free modal angle. L23 bounds
  this by `1/(1+c²)<1`; an 80-decimal scaled scan confirms margin `~3c²` through `c=10^-20`.
  This removes the sharp-edge falsification concern but does not prove a neighbourhood or L44.
- **L45 splits the polynomial core into two exact residual squares plus one central defect.**
  Its quadratic form is
  `3(1−b²)||x−aBy||²+3(1−a²)||y−bCx||²`
  `+(1−a²)(1−b²)(||x||²−||Cx||²+||y||²−||By||²)`.
  This proves the full parameter boundary and the center directly and isolates the only possible
  interior loss as `diag(I−C*C,I−B*B)`. The next analytic attack should exploit the separate
  quadratic dependence on `a` or `b`, reducing any convex vertex branch to a discriminant as in
  L30 rather than attempting another free block-norm bound.
- **L47/L48 add a sharp coupled block invariant and close both transfer axes.** The common
  modal matrix compresses both block traces to one scalar `w`. The majorant
  `tr−det²/tr` turns `||B||||C||≤2` into a rational scalar inequality, and the elementary
  relaxations `k≤min(1,4c/(1+c²)²)` and
  `H(p)≥p(1−c⁴)²/(1+2c²−c⁴)²` suffice. Four blow-up charts resolve the only equality corner
  with exact nonpositive integer Bernstein coefficients. Factoring `R_{a,0}` through an
  orthogonal block, then swapping parity, proves (RT) whenever `ab=0`.
- **L49/L50 close the complete zero-node boundary.**  The strengthened nome bound
  `k≤4c/(1+4c²)` is exact-certified on `c≤1/2`.  At `p=0`, every small-gain block is rank one;
  BE becomes (55), and `k≤1` plus a 125-coefficient exact Bernstein certificate proves it for
  `c≥1/2`.  For `c≤1/2`, exact ridge-centered and largest-coordinate charts certify the
  resulting polynomial even on the curve `b=1,a=2c` and at the degenerate origin.  The next
  obstruction therefore has `p>0` as well as `ab≠0`.

## Previous Epoch-6 milestone — EL4 PROVED
- **EL4 is now an analytic theorem**, not a grid conjecture. New L17: `SG ≥ 0` on the real
  node interval implies D2, by converting to hyperbolic coordinates and comparing the associated
  Dirichlet form against the constant Schwarzian `−2` Möbius model.
- For the squared-ellipse map `G(k sn²U)=sin²(sU)`, `SG≥0` reduces to the exact inequality
  `1/sn²(x|m)−s²/sin²(sx)≥(1+m−s²)/3`; the remainder is the positive Weierstrass series
  `8s² Σ n q^(2n)/(1−q^(2n))(1−cos(2nsx))`. Therefore `Theta≤1` and `rho≥0` throughout the
  nondegenerate even midpoint phase. Proof: `proof/el4_schwarzian_theorem.md`; 70-dps regression:
  `experiments/el4_schwarzian_check.py`.
- **Scope audit:** this does not yet prove the whole elliptic 4×4 slice. L16 proves midpoint
  stationarity and conditional `q1=q2=1/2`. **The globality debt is now closed by L18**:
  `sup ||F(u1)Q1+F(u2)Q2||=max(1,d||Q1−Q2||)`, and for norm>1 the midpoint automorphism is the
  unique nonconstant global even maximizer up to phase (`proof/even_pick_globality.md`). Thus the
  complete even sector of the elliptic slice has rho≥0. Odd and degree-one phases remain.
- **Immediate correction (same epoch): parity is false on the elliptic slice.** Exact nodal
  search at weights `(.8,2.4,1.1,.3)` finds shifted `b_β`, β=.65306547, K=1.569762, above odd
  numerical 1.511514 and even global 1.489136; diag=1.8e-8 and rho=+0.15504. This mirrors Kenan
  Li's 3×3 Region-II Möbius phase. Reproducer: `experiments/slice_phase_audit.py`.
- **New bypass route L20:** solve `I≤P≤tI`, `T*PT≤P` for `T=φ(A)`. If `t≤4` uniformly, then
  `P^(1/2)TP^(−1/2)` is a contraction with similarity condition≤2, so von Neumann proves the
  complete bound for the entire elliptic slice without classifying phases. SDP: 80/80 random
  cases pass, worst `t=3.92642` (`sqrt(t)=1.98152`). Analytic construction of P is now the
  highest-leverage slice target. `experiments/slice_cb_sdp.py`.

## Previous frontier (2026-07-21 late session) — level-4 theory + D2 landscape
- **L15 PROVED**: level-4 nodal closed form K² = (T+√(T²−4δ²F₁²F₂²))/2 (frame invariants
  p_ij, δ only; verified 1e-16); odd-phase stationarity law explicit (1e-10); Hellmann–Feynman
  frame identities. Extremal problem = 4-point Pick problem with explicit objective.
- **L16 PROVED**: even phase = deg-1 Möbius in collapsed variable; midpoint law (symmetry
  proof) + q = 1/2 (involution identity) — DOMAIN-GENERAL, verified off-slice (2e-7).
- **EL4 (central target)**: even-phase ρ = 1 − Θ(k,U₂) explicit elliptic formula; ≤ 1 verified
  (50×50 grid, ε⁴-edges at dps 80); U₂→0 edge = Landen theorem exactly. Θ-form validated
  vs exact machinery (6-7 digits) and on a bi-conic case (sign/order).
- **Phases on slice**: odd {0,±α} / even {±α} / Möbius (= odd family's α→1 endpoint);
  family-restricted searches can miss the global phase — always cross-check best_extremal.
- **D2 landscape mapped, soft routes ALL DEAD (certified)**: free-convex false (z+tz² exact,
  Herglotz-random ≤1.45); odd+G'-increasing FALSE (linearized step-deficit closed form +
  certified counterexample V−1 = +9.8e-8 at 40 dps, entire-poly step, Noshiro–Warschawski);
  symmetric-node case PROVED (one line); wedge family exactly solvable (γ ≤ 1 ⟺ pass);
  convex trace bound |d/ds log g'| ≤ 2 proved (Poisson). Ellipse/sq-ellipse pass with FREE
  nodes (focus-pinning not needed — mechanism is the map class regularity).
- **Surviving proof route for EL4**: deformation path — V ≡ 1 at disk, show dV/dk ≤ 0 with
  the explicit Möbius-Jacobian kernel functional applied to the elliptic velocity field
  ∂ψ_k/∂k (linear in velocity at each k). See D2_landscape.md end.

## Objective & status
Resolve Crouzeix's conjecture. OPEN in literature (verified 2026-07). No counterexample found by
campaign (all candidate violations were certified numerical artifacts). Campaign has reduced the
conjecture to a single sharply-supported positivity conjecture and proved it for n = 2:

**H-r: ρ := Re⟨(f₀·Φ(f₀))(A)x₀,x₀⟩ ≥ 0 at extremal pairs for Ω = int W(A).**
H-r ⟹ Crouzeix (via SV24 Thm 6.1 K ≤ 1+√(1−ρ), + shrinking; proof/strategy_S.md).
Posed-but-unattacked in literature (SV24 §6 remark); no refutation exists (searched).

## Proved by campaign (see proof/)
- **2×2 theorem**: ρ = 1 − π/(2K(m)) ≥ 0 in closed form at critical (confocal-ellipse) domains;
  equality iff disk. Geometric form φ'(0) ≥ φ(1). [Modulo α=0 symmetry step — numerically certain,
  write-up pending.] Also re-derives K ≤ 2 for 2×2 via K = √k(h+√(1+h²)) ≤ 2. (rho_2x2_theorem.md)
- P1: K² + ρ ≤ Kq, q = √(W²−|β|²) (phase-preserving; refines SV24's K²+ρ≤2K by β-subtraction).
- Ceiling: K² + 2ρ + G² ≤ 4 ⟹ at K=2: ρ ≤ −G²/2; with H-r ⟹ Crabb rigidity at K=2.
- Hadamard variation of Φ (single-layer form); Λ-density formula for frozen dρ; ∮Λds = 0 at disks
  (verified numerically to 6 digits).
- DLP contact degeneracy at Ω = W(A): ker P(σ_θ) = span((σ_θ−A)x_θ) (derivation done, write-up pending).
- Scalar reduction L10 DISPROVED (odd-symmetric mechanism) — operator structure essential.

## Key structural picture (Epoch 3–4)
- Criticality = confocality: Kippenhahn curve of ∂W(A) has foci at σ(A); the 2×2 mechanism
  (confocal ellipse ⟹ ρ ≥ 0 via elliptic integrals) is the n=2 case of a general
  "Pick-problem at eigenvalue images + confocal domain" structure.
- Extremal problem = n-point Pick-boundary problem at wⱼ = φ(σ(A)); extremal zeros ≈ critical
  points of ∏b_{wⱼ} (EXACT at n=2 — hyperbolic midpoint; close at n=3, deviation = non-Hermitian
  weighting). (zero_geometry.py)
- ρ-forms: ρ = K·Re⟨g₀(A)u₀,x₀⟩ = 1 − ReΣrᵢ⟨hᵢ(A)x₀,x₀⟩ ("capacity partition", terms ≈ harmonic
  split; Σ = 1 exactly on disk).
- Empirics at critical domains: ρ ≥ 0 always; ρ ↑ as Ω ↓ W(A) (offsets, Minkowski; textbook-monotone);
  adversarial min-ρ floors: +2e-4 (n=3), +1.7e-2 (n=4), +1.6e-2 (n=5, weak search); n=6 pending.
  ρ = 0 attained at disk-like configs. Hadamard: pair-response dominates frozen term (heavy route).

## MILESTONES (Epoch 5, 2026-07-21) — THE LANDEN THEOREM
- **sym3 ρ-positivity now CLASSICAL-COMPLETE closed form** (proof/landen_theorem.md):
  W(A) ellipse explicit (semi-axes from support function; foci ±e proved); z₁ = e/√2 EXACT
  (quarter-period u₁ = K/2, sn(K/2) = 1/√(1+k′)); r₁ = eπ/(2√2 K k);
  **ρ = 1 − π/((1+k′)K(k)) = 1 − π/(2K(k₁)) ≥ 0** (Landen descent k₁ = (1−k′)/(1+k′)),
  equality iff disk. Verified: independent mpmath 30 dps + Landen identity 1e-31.
- ζ = z²-collapse ≡ Landen transformation ⟹ **Landen-tower conjecture** for the symmetric
  tridiagonal family (induction over levels; needs elliptic range at each level).
- **Elliptic sym4 slice found analytically: b_j = c·a_j** (= c-deformed weighted shifts!);
  W(A) exact ellipse (verified 3e-15). [CORRECTED later: foci = ±e₁ = OUTER eigenvalue pair,
  τ₁ = √k — see proof/landen_theorem.md; an earlier foci² = e₁²+e₂² note was an algebra slip.]
  ρ > 0 on slice; no pure-modulus Landen formula (marked point breaks second descent) — the
  correct closed form is EL4's ρ = 1 − Θ(k,U₂) (see NEWEST above).
- Adversarial H-r floors (all positive): n=3 +2e-4, n=4 +1.7e-2/+7.2e-2, n=5 +1.6e-2/+1.8e-2
  (n=6 dense: NOT OBTAINED — uncapped acceptance loop + solver stall at diag ~1e-3, job killed 2026-07-21; see pitfall P7; redesign on structured families). Direct ratio
  searches DONE: n=6 best 1.148, n=7 best 1.465 — nothing near 2.

## MILESTONES (Epoch 4, 2026-07-20)
- **sym3 (= GKL 2018 class, real slice) H-r PROVED-by-reduction**: exact ζ = z² transform onto the
  canonical 2×2 confocal configuration; stationarity = pseudo-hyperbolic-midpoint law; π₀ = 1/2
  exact; ρ = (α²/2)(g₀(e)−g₀(0)) = 1 − S ≥ 0 via the 2×2 elliptic theorem. NEW MECHANISM,
  known class (novelty calibrated vs arXiv:1701.01365). proof/sym3_reduction.md is the master note.
- Rigidity of the inequality: S ≤ 1 FALSE for free odd convex maps (1.019), false on partial-focal
  synthetic loci (up to 1.24) — full criticality essential; no soft proof exists.
- L13 Clark-type transition proved; L14 three-factor formula; collapse theorem v = α².
- **sym4 = FIRST BEYOND-LITERATURE TERRITORY: ρ > 0 confirmed numerically (3 cases); extremal
  Blaschke is ODD ({0,±α} zeros) ⟹ f₀ = z·F(z²) collapse exists; squared boundary NOT elliptic
  (resid 1e-4–1e-3).** Next: close sym4 analytically = first new Crouzeix class in campaign.

## Current next actions (Epoch 6, refreshed 2026-07-22)
1. **Close L57's positive compact tail.** Keep L58's exact `a=b=1` face and positive inward
   first variation intact, and certify only the quadratic-and-higher `(1−a,1−b)` remainder
   with a shared Bernstein nome coordinate and the surviving `X,R,Y` variables. The raw global
   coefficient l1 norm is rigorous but too coarse. Do not resume tolerance-based subdivision;
   the target is only `c+<c<2^(−2/3)`. This proves L20 for the complete elliptic `4×4` slice.
2. **Shifted Möbius phase**: derive its exact stationarity/rho formula (Kenan-Li quartic analog)
   and prove rho≥0 or K≤2. Definite parity is false.
3. **Bi-conic Schwarzian test**: compute `SG` for the off-slice collapsed map on the critical
   real interval. `SG≥0` would extend L17 immediately; otherwise test the weaker Sturm-potential
   comparison that the proof actually needs.
4. **Odd phase positivity**: ρ = B₁g₁q₁+B₂g₂q₂ ≥ 0 given the L15 stationarity law
   (3-parameter; interlacing τ₂ < α < τ₁; term-1 dominance observed). Try the same
   deformation/kernel machinery.
5. Rigor debts: n=6 structured floor; 2×2 α=0; contact degeneracy; novelty audit.
   General-n work must include symmetry-breaking phases; the former parity-collapse induction
   remains valid only inside a chosen parity sector.
Keep committing+pushing after each task (user instruction).

## Files map (handoff-ready, 2026-07-22)
proof/ — read in this order for the current frontier:
  rho_positivity_program.md (MASTER program), el4_schwarzian_theorem.md (EL4 proof + L17),
  even_pick_globality.md (L18), slice_closed_form.md (level-4 theory: L15/L16/EL4/D2-crit), D2_landscape.md (soft-class
  falsifications + Schwarzian route), landen_theorem.md
  (sym3 closed form), sym3_reduction.md (collapse theorem), sym4_program.md (family setup),
  rho_2x2_theorem.md (n=2 base), graded_induction_skeleton.md (general-n plan);
  background: strategy_S.md, P2_target.md, refined_master_inequality.md,
  epoch2_extremal_structure.md, track_A_crouzeix_palencia.md.
experiments/: crouzeix.py (basics: poly_A, nr_support, ratio_inner/outer, crabb_matrix),
  extremal_pullback.py (ψ-domain exact machinery), theodorsen.py (+GeneralPullback: arbitrary
  convex domains, unitality certificate), minkowski_test.py (best_extremal — THE solver; always
  use it to cross-check phases), ellipse_sandbox.py (2×2 exact), slice_exact.py (elliptic sym4
  slice, 25 dps; family-restricted — see pitfall P5), slice_invariants.py (L15 invariant
  formulas + verifications), sym3_structure/sweep/analysis.py, sym4_probe/sweep.py (OTHER-branch
  taus/f0e fields BUGGY — pitfall P5), Hr_test.py / Hr_adversarial.py (certified min-ρ search),
  zero_geometry.py (phi_of_points); historical: L10/L12/scalar/adversarial_L7/search.py.
  Epoch-6: el4_schwarzian_check.py, even_pick_globality_check.py, slice_phase_audit.py,
  slice_cb_sdp.py (requires cvxpy; exploratory similarity SDP),
  slice_similarity_duality.py (exact modal reconstruction + primal/dual regression + grid),
  slice_boundary_check.py (L23–L26 stress regression),
  slice_upper_block_certificate.py (exact factorization + finite scalar certificate for L26),
  slice_coupled_defects.py (L27 reconstruction and KKT rank regression),
  slice_odd_block_check.py (L28 identities and L29 floating-point stress test),
  slice_block_product_certificate.py (L47 exact integer Bernstein certificate),
  slice_sharp_nome_certificate.py (L49/L50 exact zero-node certificates),
  slice_positive_face_audit.py (L58 exact positive-face/first-deficit audit).
Data: sym3_sweep_s51.jsonl (40 rec), sym4_sweep_s61.jsonl (20 rec, ρ column trustworthy).
Ledgers: LEMMA_LEDGER.md, APPROACH_LEDGER.md (pitfalls P1–P8 — READ BEFORE ANY SEARCH),
LITERATURE_LEDGER.md, COUNTEREXAMPLE_SEARCH.md. Audit: chatgpt/FABLE_RESEARCH_AUDIT.md
(reconciled 2026-07-20). Restart: checkpoints/RESTART_PACKET.md (paste-ready instruction).
