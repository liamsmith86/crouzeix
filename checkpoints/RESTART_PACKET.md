# CROUZEIX RESEARCH CHECKPOINT — NOT A FINAL RESULT
(Restart packet per goal.txt forced-stop protocol; written proactively 2026-07-20 while work continues.)

## Current general Crabb frontier (2026-07-22): L115

L114 closes the full repeated-`C3` local neighbourhood for every fixed copy multiplicity.
L115 then returns to arbitrary block size and intersects L65's equality kernel with the
Lewis--Overton circular-range tangent manifold.  Modulo affine-unitary motions, the tangent
part has dimension `2p−4` out of L66's total `2p−2`; only one complex soft normal survives.
It is mode two, represented by `C_p*`, and its exact central family satisfies
`W(C_p+cC_p*)={z+c conjugate(z):|z|<=1}`.  Hence the next general-size task is a
dimension-independent condition-four certificate for this Chebyshev--Lobatto elliptic family,
then a tubular absorption using L65's coercive normal complement.  The pure family is already
proved at `p=3` (L68--L73) and `p=4` (L20).  Do not compute all residual jets and do not begin
a generic `5 x 5` slice grind.  See `proof/crabb_disk_tangent_intersection.md`.

## Completed repeated-C3 frontier (2026-07-22): L114

L73 proves a complete-2 neighbourhood of one `C3`; L74--L93 reduce arbitrary repeated
copies to weighted terminal blocks, with all fixed directions closed.  L94--L100 give
quantitative terminal gaps and locate the unique sharp normal center
`w=3conj(d)^2/(8sqrt(2))`.  L101 proves the crossing-stable conformal tangent.

L102 is a **support-only** identity: the fourth Feshbach coefficient is scalar at the center
and has zero full transverse derivative.  An audit caught that the default full metric still
has derivative `(25d^3/8)[[0,1],[1,0]]`.  L103 now cancels it exactly using the free third-metric
block
`aU`, `U=5d^2[[0,-1],[1,0]]/(16sqrt(2))`; arbitrary retained conformal coefficients cancel.
Thus the total fourth transverse remainder is quadratic and absorbable by L100's cubic gap.
L104 proves that the same adjacent-order free metric map is onto every symmetry-allowed
traceless linear edge coefficient at **all** finite orders.  The live task is therefore not
another coefficient.  L105 replaces all forced lower/Stein coefficients by an exact analytic
metric chart with one free cross block and one upper endpoint.  L106 freezes the actual
perturbed Riemann map and proves a uniform Cauchy expansion in the transverse operator, avoiding
all differentiation of the nonsmooth domain map.  L107 retains the frozen normal
certificate's PSD Stein Schur slack as an exact chart parameter, and L108 proves that it may
be tightened safely to zero on the normal anchor.  L109 then computes the joint weighted jet:
with a retained slack fraction `theta`,
`E2(N)=−5(1−theta)a²I/8`, `E2(A)=5theta a²I/8`, and
`E3(N)=E3(A)=−16(1−theta)m3I`.  Thus zero slack exactly cancels the normal/transverse
quadratics and preserves L100's negative
`−delta(delta²+2r²)I/4` cubic.  L110 closes the analytic remainder: the whole-normal support
gap is `O(delta(r+delta))`, so `h=f_A∘f_N^−1` and its induced Stein slack have weighted
degree two.  In the difference between the actual zero-slack endpoint and the inherited normal
endpoint, every post-cubic Taylor monomial therefore has bound `C delta(r+delta)^3`; L100's
cubic dominates uniformly.  The complete two-copy terminal tube is proved.
L111 removes the first arbitrary-copy obstruction: at `Z_ir tensor I_k`, every tangent
preserving the scalar-support relations is a unitary commutator plus the same three
two-dimensional parameters on each multiplicity coordinate.  There is no hidden
large-multiplicity terminal modulus.  L112 closes its normal boundary:
at `s diag(Ip,−Iq)`, the only unitary-invariant equality tangent is the rectangular edge
`B−C*`, whose SVD is a direct sum of L110 two-copy pairs and normal leftovers.
L113 closes proper common-top spaces: the flag forms are coercive off the blockwise tangent
kernel, while a block-sign unitary makes every active diagonal remainder even in the cross
coupling.  Thus `−c s²||R||²` absorbs `C s³||R||²`, and flag rank collapse only moves to
another member of a finite stratification.  The complete flat core is locally proved for each
fixed multiplicity.  L114 performs the final merge: L61 is strict off the zero-Jensen face;
L86 is a complete sum of negative strong forms on it; L87's kernels are exactly the L113 flat
blocks; and L72 ambient stationarity plus L105--L110 analytic remainders gives a stratified
Morse--Bott tube.  Therefore every fixed `I_m tensor C3` has a full complete-`2`
neighbourhood.  This does **not** solve the general conjecture.  Next re-audit the global L21
equality set and decide which non-`C3` local model or global mechanism is genuinely next;
do not start a blind larger Schur grind.  See `RESEARCH_STATE.md` and
`proof/repeated_p3_local_neighbourhood.md`.

## Latest completion (2026-07-22): elliptic 4×4 slice proved

L59 closes the final positive-sign interval with a full deficit-centered
Arb/Bernstein certificate. Ten exactly adjacent rational boxes cover `c+→.63`; both determinant
and all three final-minor charts pass. The clean run forced regeneration of both 197,563-record
tables and both 207-term corner squares, passed 10/10, and exited zero. Provenance artifact:
`experiments/positive_tail_full_20260722.log`, SHA-256
`fb79b2dfc652062307d69d26a00d82ff20c4133a0043eb6077458ea7cf70928c`, certified code commit
`3dd51884…`. Because `.63³−1/4=47/10^6>0`, it overlaps L42. Hence L20 is proved for every
`A=S_a+cS_a^T` with arbitrary positive three-weight `S_a`, `0<c<1`: its numerical-range ellipse
is a complete 2-spectral set. This is apparently new versus the closest audited fixed-Crabb-weight
theorem, but it is **not** arbitrary 4×4 and not the general conjecture.

The campaign now returns to the general L21/H-r frontier. Do not start a 5×5 slice grind.
Hartz--McCarthy scalar shifts are an exact restatement of the cb target, not a CP bridge; a
reproducible order-three Toeplitz SDP also numerically rules out all positive-state
scalarizations on one dense 3×3 case. See `proof/general_similarity_probe.md`.

**General gate update:** the repeated-Crabb experiment is now complete. It tests block sizes
3/4, multiplicities 2/3, dense/cross/operator-weight perturbations, and accepts 118/120
perturbed records with no value above four (max `3.999844312562`). The two rejects have bad
primal/dual gaps. Two ultralocal cross directions remain below four to `delta=1e-5`, with
stable first-order drops `1.557/1.637`. This supports, but does not prove, local maximality.
Reproducer and data are in `experiments/general_similarity_equality_probe.py` and the three
`general_similarity_equality*_s9173401.jsonl` files.

**Local first-order theorem:** L61 computes the exact one-sided support derivative at a
repeated Crabb block and, because the block is nilpotent, reduces the conformal operator
tangent to only `p−1` support Fourier coefficients. Its active-kernel SDP dual collapses to one
density matrix and has exact value `−4` times a support-compression Jensen gap, hence is
nonpositive for every perturbation. A strict feasible-metric lift proves the same upper Dini
bound under a uniform conformal expansion. Clarabel/SCS and the closed form agree on 12
directions within `3.17e-8`; nonlinear slopes agree to relative error `≤1.20e-4`. This is only
first order, not a full neighbourhood theorem. Next classify zero-gap directions, calculate
second order where needed, and close nonsmooth map regularity. See
`proof/general_similarity_tangent.md`.

**Second-order reduction:** for a single Crabb block the L61 Jensen gap is identically zero.
L62 uses first-order complementarity and a second-order PSD Schur lemma to reduce the next
coefficient to three finite affine block LMIs in the metric coefficients. A two-resolution,
two-solver run on 12 full/structured `p=3,4` directions gives only negative coefficients,
`−3.94535170≤e≤−.00284309`; analytic and independently fitted gauges agree within `3.57e-7`.
The second support variation, normal-angle correction, and Schwarz integral now give `H` in
finite analytic form. This is still a general-`p` reduction, not a universal sign theorem. See
`proof/general_similarity_second_order.md`.

**Exact `3×3` and `4×4` second variations:** L63 gives
`e₃(E)=−2(Re(E01−E12))²−21|E20|²/4≤0`; L64 gives the rank-eight five-square identity
(28). Regenerating SymPy checkers start from all 18/32 real coordinates, verify the exact support
eigenpairs, rebuild both conformal coefficients, eliminate the metric variables, and prove the
residuals identically zero. Thus the upper second-order Dini change of `t*` is nonpositive at
both low-order Crabb blocks, strict off their displayed equality spaces. Next seek the
arbitrary-`p` factorization and classify those equality directions before going to third order.
See `proof/general_similarity_second_order.md` §§8–9 and `experiments/p3_second_order_identity.py`,
`experiments/p4_second_order_identity.py`.

**Arbitrary-size second variation:** L65 now subsumes those signs for every single Crabb block.
The support defect is an exactly factored path Laplacian. Circle grading collapses each nonzero
mode to `−t*ᵊ t`, with `ᵊ^{-1}=(I−H_r)/4+c_kqq*` for `k≥3` and PSD limits for
`k=1,2`; grade zero is an exact weighted shift controlled by the numerical-radius power
inequality. The form has rank `p(p−2)`. Checks cover every mode through `p=30`, while L63–L64
remain exact symbolic base cases. The live task is no longer the sign: classify the
`p(p+2)`-dimensional equality space and extend to repeated common-maximizer faces. See
`proof/general_crabb_second_variation.md` and
`experiments/general_crabb_second_order_modes.py`.

**Equality quotient:** L66 removes exact symmetries from the L65 kernel. The raw dimension is
`p(p+2)`, while unitary conjugation, translation, and real scaling span an exact flat orbit of
dimension `p²+2`; only `2p−2` real directions remain. Grade by grade this is one complex mode-1
direction, two complex mode-2 directions, and one complex direction in each higher paired mode.
At `p=3`, numerical fits suggest quartic decrease `−4|w|²(|w|²+|z|²)` and sixth-order
decrease on `w=0`, but those coefficients are not proved. The live task is their exact
derivation, not another raw-kernel scan. See `proof/crabb_second_order_equality.md`.

## Immediate audit gate (2026-07-21, post-L45)
**Load-bearing audit passed:** L17's exact ground-state identity with
`y=1/sqrt(g')` fixes the Dirichlet-form direction, and L21's parity restriction follows by
averaging its linear dual triple before rescaling the invariant ray.  See
`proof/load_bearing_audit.md`.

**General similarity gate passed numerically:** 102/120 certificate-gated general complex
matrices (`n=3..8`) had `t*(phi(A))<=3.8399296`; none exceeded four.  On the sharp near-Crabb
case, shrinking the outer offset `.02→.00125` raised `t*` from `3.68956→3.97783`, always below
four.  Seventeen polygonal near-normal maps and one inaccurate SDP were rejected rather than
counted.  This makes the complete-bounded route credible, not proved; see
`proof/general_similarity_probe.md`.

**L44 sharp boundary passed exactly:** its true boundary is `p=a=b=0` (the modal angle is flat),
not L29's stale `p=p*+c^4 x` ridge.  There the energy is `k(1+c²)/(4c)<=1/(1+c²)<1`; an
80-decimal scaled scan through `c=10^-20` finds margin `~3c²` and no excess.  This is L46 and
does not prove a neighbourhood.  All three L29 exact certificate scripts have now regenerated
cleanly; the 85-minute centered completion rebuilt both deep expansions and the final detached
job ended with `PASS 2026-07-21T15:10:52-07:00`.

**New exact milestone:** L47 proves `||B||||C||≤2` for the two elliptic-slice parity blocks.
One orientation scalar gives both traces and determinants; elementary nome relaxations reduce
the claim to exact integer Bernstein polynomials. Four blow-up charts close the singular corner
without floating-point margins. L48 then factors the transfer through an orthogonal block and
proves (RT) on both complete axes `a=0` and `b=0`. The remaining transfer obstruction has
`ab≠0`; see `proof/slice_coupled_defects.md` §6.3 and
`experiments/slice_block_product_certificate.py`.

**New zero-node milestone:** L49 proves the sharper nome estimate
`k≤4c/(1+4c²)` on `c≤1/2` by retaining two Jacobi-product factors; the residual has 24 positive
exact Bernstein coefficients.  L50 then proves BE/RT on the complete `p=0` face after an exact
rank-one reduction.  The high half uses 125 nonnegative coefficients.  The low half is closed
by a finite exact cover centered on `b=1,a=2c`, with largest-coordinate blow-ups at the origin.
The next live chart is an interior neighbourhood with `p>0`, still centered on that ridge.

**New full-core reduction:** L51 uses `E_o≻0` to reduce the remaining PSD theorem to one
leading `3×3` minor and the full determinant.  Set `P=p²,t=1−o`, retain the complete cubic
envelope through `r=p(s+γP)`, and use `(P,t)=(w,wv)` or `(wv,w)`.  Exact coefficient
collection gives common orders 3 and 4, respectively.  Fixed-`c` Bernstein tests pass for
the complete envelope at every tested `c=.001,...,.629`, but this is not yet a theorem over
continuous `c`.  An exact order-one secondary blow-up resolves `v=0,|b|=1` at fixed `c`.
L52 proves both extreme-orientation faces analytically: the apparent square at `P=1,o=0`
has scalar factor `G_x<0` throughout the nondegenerate range, so there is no moving interior
zero curve.  Use these strict factors to normalize the next directed Taylor certificate; see
`proof/slice_core_projective_reduction.md`.  Frobenius and L47's rational block-norm majorant
were falsified and should not be retried.  (Only the claim that BE is maximized at orientation
endpoints was falsified; the exact transfer theorem on those faces is now L52.)

L53 adds a tertiary exact order-one blow-up at `(h,1−a)=(0,0)` in the first secondary
chart.  `experiments/slice_projective_core.py` now regenerates every core record and chart
order, so do not use the old `/tmp` polynomial cache as a source of truth.  The first Taylor
prototype also incorrectly replaced the zeroth term by the interval between zero and its
value.  After retaining that term and making `c` a common Bernstein axis, determinant and
final minor charts pass continuous boxes `[.01,.0101]` and `[.1,.101]`.  Next implement the
directed, memory-bounded finite cover.  This implementation is now committed as
`slice_projective_interval_certificate.py`, and L54 rigorously closes the complete interval
`[.01,.020736]` with four boxes (all five final charts, both signs). Treat `c→0` separately after
normalizing `k/c`, `gamma_-/c²`, and `(gamma_+−gamma_-)/c⁴`; two further determinant powers
cancel there.
L55 performs the exact normalization: determinant order `c^9`, final-minor order `c^2`.  At
`main_w=0, ratio=1, a=1, b=0` its first two exceptional coefficients are
`1536(A+3U)` and `192[8A2+24U2+3V²+18(B−2σ)²]`; the sole surviving positive-sign ridge
`B=b/c=2` starts with `46080c^4`. Further exact products close the main-orientation remainder
and identify the det1 hierarchy (20)--(23) and det0 `a=0` forms (24)--(27) in
`proof/slice_core_projective_reduction.md`.
The det0 main-orientation chart and all four widened corner charts certify for both signs on
`[0,.005]`. Its former equality line ends in coefficient-positive forms (26)--(27), and all 26
local hierarchy charts close both signs; the negative run took 3033 seconds and every chart
passed at the root. Det1's endpoint hierarchy ends in an exact positive-definite quadratic.
Two 59-chart directed runs certify the full det1 main chart for both signs on `[0,.005]` (about
84 and 76 minutes). Exact minor forms (28)--(33) and full integrated runs now certify all three
final minors for both signs on `[0,.005]` (global leaves `81/32/184` positive and
`71/32/188` negative). Scale-free nome/main arm charts remove the det0 cutoff; exact extraction
shows that the apparent `S=1/2` transition was only a centered-model seam. The det0 global
complements close in `235/399` leaves at depths `8/9`. Thus L55 proves the complete core on
`[0,.005]`, every chart and both signs. L56 subsequently closes `.005→.01` in four directed
boxes, so L55+L56+L54 now prove (RT) continuously through `.020736`. L57 certifies the compact
negative sign through `.63` (overlapping L42) and the positive sign through the exact rational
`c+=.5566585294072849…`. At that stage only the positive tail
`c+<c<2^(−2/3)` remained; L59 now closes it.

L58 factors the sharp positive determinant face. At `a=b=1`, chart zero is
`81c^9(1−k)^2(1−kX)^2[R+(1−RX)C^2]^4`; chart one has the analogous
`(1−kRX)^2[1+R(1−X)Ctilde^2]^4`. The coefficients of `A=1−a` and `B=1−b`
coincide and are nonnegative, with residual factor `2+3k+3Zk−8Zk²`. The exact checker
`experiments/slice_positive_face_audit.py` regenerates all 197,563 records per chart and the
higher-deficit l1 norm `41,235,531,913`. That global norm is too coarse: retain correlation in
`c,X,R,Y` for the quadratic-and-higher remainder; do not use a tolerance or more unbounded raw
subdivision.

L45's convex-in-`a` branch is common (772,130 of 2,000,000 random quadratic-form probes) but its
discriminant remained positive in the sample; no concavity-only shortcut exists.

## Exact current frontier
Campaign has reduced Crouzeix's conjecture to **H-r: ρ = Re∫f₀·Φ(f₀)dμ ≥ 0 at extremal pairs on
the critical domain Ω = int W(A)** (⟹ conjecture via SV24 Thm 6.1 + shrinking; proof/strategy_S.md,
proof/rho_positivity_program.md). H-r: PROVED for all 2×2 (closed form ρ = 1 − π/(2K(m)), modulo
the α=0 symmetry step); PROVED-by-reduction for the sym3/GKL 3×3 elliptic class (proof/sym3_reduction.md
— ζ = z² collapse, midpoint law, π₀ = 1/2); NUMERICALLY CONFIRMED beyond all known classes (sym4;
adversarial floors positive for n=3,4,5). No counterexample to the conjecture found anywhere
(all "violations" ever seen were certified numerical artifacts — APPROACH_LEDGER pitfalls P1–P8).

## Epoch-6 addition (2026-07-21): EL4 proved
`proof/el4_schwarzian_theorem.md` proves the explicit elliptic inequality `Theta≤1`. The new
sharp lemma is `SG≥0 ⇒ D2`, proved by Schwarzian/Dirichlet-form comparison; the squared-ellipse
map has `SG≥0` because its gap is a positive Weierstrass Fourier series. The 70-dps regression is
`experiments/el4_schwarzian_check.py`. Do not redo the modulus-deformation derivative.

**Historical EL4-only status (superseded by L59):** this was not yet the 4×4 theorem.
Midpoint globality is proved by L18 in
`proof/even_pick_globality.md`, but parity/symmetry-breaking and the odd/Möbius phases remain.
Symmetry of the objective does not by itself prove that every maximizer has definite parity.

**Correction:** parity is actually false on the slice. `experiments/slice_phase_audit.py` finds
a robust shifted degree-one Möbius phase at weights `(.8,2.4,1.1,.3)` that beats both parity
sectors. Its rho is positive. Do not attempt to prove definite parity.

**Best new bypass:** `experiments/slice_cb_sdp.py` searches for `I≤P≤tI`, `T*PT≤P`, `T=φ(A)`.
All 80 random slice cases gave `t<4` (worst 3.92642). A uniform analytic metric proves the
complete-2 bound for the entire elliptic slice via von Neumann and avoids all phase bookkeeping.

## Epoch-6 addition (2026-07-21): exact similarity duality and modal reduction
`proof/slice_similarity_duality.md` proves two new exact lemmas. L21 gives, for strictly stable
`T`,

`t*(T) = max(1, sup_{Z≥0} tr(Z−TZT*)_-/tr(Z−TZT*)_+)`.

L22 puts every elliptic-slice `T=φ(A)` into a three-parameter modal form `(c,r,u)` with
`tan(v)=r tan(u)` and reduces the primal to coupled 2×2 metric LMIs. Chiral averaging reduces
the dual to two 2×2 PSD blocks. Therefore L20 is exactly the concrete trace inequality
`tr D_-≤4 tr D_+`. A 1215-point deterministic grid passes (max `t=3.999771308` at the singular
corner `(.001,.97,.03)`), but this is not a proof. Reproducer:
`experiments/slice_similarity_duality.py`. Simple diagonal/Gramian/equality metrics have already
failed; attack extreme dual block pairs rather than repeating blind primal ansatzes.

**Next proved pieces (L23–L25):** `proof/slice_boundary_theorems.md` proves the sharp nome bound
`k/c≤4/(1+c²)²`, uses it to prove one full modal-block estimate `||C||≤2`, and proves L20
on the singular `c→0` face by reducing to a nilpotent scalar weighted shift. The boundary constant
is sharp at the Crabb weights `(√2,1,√2)`.

**L26 now proves the opposite block:** `proof/slice_upper_block_theorem.md` proves `||B||≤2` by
factoring the reciprocal-orientation determinant, converting its upper root to a Möbius barrier,
and minorizing the inverse ellipse map by its positive-coefficient cubic. Three scalar theta
inequalities are certified by exact rational bounds plus 13,500 outward-rounded algebraic
interval boxes (`experiments/slice_upper_block_certificate.py`). Thus all one-block dual phases
are closed; genuinely coupled witnesses remain.

**L27 now parameterizes the coupled obstruction exactly:** with modal Stein defects `Q_o,Q_e≥0`,
every contraction metric is
`H_o=K∘(Q_o+cΣQ_eΣ)`, `H_e=K∘(Q_e+c^{-1}ΣQ_oΣ)`, where
`K_ij=(1−τ_i²τ_j²)^{-1}`. Complementary slackness gives
`rank Z_i+rank Q_i≤2`. Because both dual blocks cannot be full and the one-block faces are
already closed, a hypothetical KKT optimum above four must have
`Q_o=aa^T`, `Q_e=bb^T`, allowing one vector to vanish. This leaves three projective parameters.
Do not try to bound arbitrary defect choices: some are badly conditioned. Exclude a KKT minimizer
above four using the sandwich activity and the conformal relations. See
`proof/slice_coupled_defects.md` and `experiments/slice_coupled_defects.py`.

**L28 reduces the generic rank-one/full face to the odd Blaschke sector.** On its sign-separated
boundary an orthogonal-colligation identity gives the exact value
`sup_{a∈[-1,1]} ||B R_a(CB)||²`; this is the upper block of
`f_a(T)=T(T²−aI)(I−aT²)^{-1}`. The parity-swapped lower block is proved ≤2. For the upper block,
the reciprocal-quadratic orientation identity leaves one explicit determinant numerator
`N(c,p,H(p),a)`. Kanas--Sugawa gives the rigorous envelope
`s₀p+a₃p³≤H(p)≤s₀p+(1−s₀)p³`. L30–L39 now prove the upper block as well, including the sharp
ridge `c→0, p→1/2, a∼3c`; hence this entire face is closed. See
`proof/slice_odd_block_reduction.md` and its three certificate scripts.

**L30 further removes the Blaschke parameter.** The reciprocal upper block has equal diagonal,
so its norm-two condition is exactly `det(M)+2|M12+M21|≤4`. In the inner-node value coordinate
`t`, both signed residuals are quadratics. One is concave and is completely closed by L26 at
`t=±1`. The other only needs `4AC−B²≥0` on `A>0, |B|<2A`, at the two cubic-envelope endpoints.
This is a two-variable `(c,p)` theta inequality; the old `(c,p,a)` determinant is obsolete as
the primary certificate target. Equations (13m)–(13p) contain the exact coefficients.

**L31 factors that discriminant into a square.** Put
`Q=2g(1+c²p)−d(kp+4c)` and `S=−2g(c²+p)+d(kp+4c)`. Exact expansion gives
`c²(1+r)²(4AC−B²)=16rg²p(1−c²)²(1−d²)−(Qr−S)²`, with
`1−d²=(1−k²)(1−k²p⁴)/(1−k²p²)²`. The new primary target is therefore
`|Qr−S|≤4g(1−c²)sqrt(rp(1−d²))` at `r=L(p),U(p)` on the branch `A>0, |B|<2A`.
The expression is concave in `r`, so endpoint reduction is automatic. The exact audit is in
`experiments/slice_odd_block_check.py`.

**L32/L33 remove cancellation and the high-nome range.** With
`λ=d(kp+4c)/(2g)` and `D=(r−p)−c²(1−pr)+λ(1−r)`, exact factorization gives
`c²(1+r)²(4AC−B²)=4[rF₁F₂−g²D²]`; evaluate `p−r` from its endpoint formula rather
than subtracting. If `cℓ≥1`, convexity reduces the entire nodal rectangle to four vertices and
an exact opposite-sign determinant factor proves all of them norm≤2. The existing L26 product
five-factor product certificate gives `cℓ≥1` for `c≥12599/20000`. Only `0<c<12599/20000` remains for
intervals/asymptotics. L34 factors `F₁` through a quadratic whose exact vertex is
`p*=(g²−4c²)/(2g²(1−2c²g))`; near zero use `p=p*+c⁴x`, not axis-aligned `(c,p)` boxes.
**L29 is now proved.** L37 deepens the centered expansion to
`|(p−p*)/c²|≤50` for `c≤1/20` and supplies `|p−p*|≤4c⁴` through `c=1/12`.
L38 excludes the two branch boundary strips and exactly covers every regular small-nome point;
the worst correction ratio is `0.907895`. L39 interval-certifies the bridge complement
`1/20≤c≤1/12` with 15,267 bisections and no unresolved boxes. L36 covers
`1/12≤c≤12599/20000`, and L33 covers the range above. Therefore L30–L32 prove the upper odd
block ≤2 for every parameter. Together with the analytic lower block, the generic rank-one/full
KKT face is closed; at this stage only L27's rank-one/rank-one face remained for L20. L59 later
closes it.

**L40 reduces the final rank-one/rank-one face to a matrix-valued inner transfer.** If
`Z_o=xx^T` and `Z_e=yy^T`, spectral positive/negative factors and two row-Gram rotations give
`tr D_-/tr D_+≤||R_{a,b}(T)||²`, where
`R_{a,b}(T)=−S⁻¹(T−A)(I−AT)⁻¹S`. Its cancellation-free block formula is valid on the full
closed square `(a,b)∈[−1,1]²`. In symmetric modal coordinates the two nodes see
`K_{a,b}(z)=[[a−bz²,−stz],[-stz,b−az²]]/(1−abz²)`, a rational inner matrix with determinant
`(ab−z²)/(1−abz²)`. The boundary is in the proved even sector; the interior is genuinely
matrix-valued and can beat all scalar lines. At this stage the remaining L20 theorem was
`||R_{a,b}(T)||≤2`; L59 later proves it. Exact audit and sharp default-case regression:
`experiments/slice_rank_one_transfer.py`; proof: `proof/slice_coupled_defects.md` §4–7.
L41 further removes both transfer resolvents: the norm test is equivalent to a quartic
polynomial `4×4` LMI. Its two diagonal `2×2` blocks are already PSD by L24/L26 and elementary
bilinear estimates, leaving only one explicit off-diagonal Schur inequality.
L42 also closes `c≥2^(−2/3)`: the symmetric-modal transfer is a contraction and the
physical similarity has condition `c^(−3/2)≤2`. Only the lower-nome Schur inequality remains.
L43 then square-completes that Schur coupling: it is equivalent to norm at most one for an
explicit four-block small-gain matrix with contracted parameters `u=3a/(4−a²)`,
`v=3b/(4−b²)`. The stronger scalar block-energy target `Σ||S_ij||²≤1` (L44) implies the
result and is search-supported on the entire rigorous cubic conformal envelope, but is not yet
proved. The softer strip `s₀p≤r≤p` fails, so retain the cubic envelope.
L45 gives a complementary exact identity: the polynomial core quadratic form is
`3(1−b²)||x−aBy||²+3(1−a²)||y−bCx||²`
`+(1−a²)(1−b²)(||x||²−||Cx||²+||y||²−||By||²)`.
It closes all `a,b` edges and the center and isolates the interior loss. For fixed `b` it is a
quadratic in `a`; use its two square endpoint values and reduce only the convex vertex branch to
a discriminant, mirroring L30.

## Strongest proved lemmas
P1 (K²+ρ ≤ Kq), L13 (Clark-type transition ⟨q(A)x₀,u₀⟩ = K∫q·conj(f₀)dμ), L14 + collapse theorem
(sym3: v = α², ρ = (α²/2)(g₀(e)−g₀(0))), **Landen theorem** (sym3 ρ = 1 − π/(2K(k₁))),
**L15** (level-4 nodal/Pick closed form + odd stationarity law + frame identities),
**L16** (even-phase midpoint law + q = 1/2, domain-general), ceiling K²+2ρ+G² ≤ 4,
**L29** (generic rank-one/full elliptic-slice KKT face), **L93** (arbitrary-copy fixed flat
core via the metric flag and two-dimensional Clifford terminal blocks), **L98** (every bounded
weighted chart at a nonnormal two-copy terminal block), 2×2 closed form; D2 partial results
(symmetric-node case, wedge family, convex trace bound).
See LEMMA_LEDGER.md.

## Failed approaches — do not repeat
Scalar reduction L10 (disproved, odd-symmetric mechanism); global-operator conditions like
∮λ_min(P)ds (not ⊕-stable); soft/free-map versions of the confocal inequality S ≤ 1 (false:
1.019 free convex, 1.24 partial-focal); Löwner/K-monotonicity constraints (vacuous);
L12-strong annihilation (false); frozen-Hadamard route (pair-response dominates);
frame-free level-4 laws + PSLQ hunts (frame-coupled — L15 explains); soft D2 classes:
free univalent (Koebe 1+δ²), convex (z+tz² exact), odd+G'-increasing (certified 40-dps
counterexample) — see proof/D2_landscape.md; pure-modulus double-Landen formula for the slice
(marked point breaks second descent); arithmetic-in-U midpoint identity (false, 1e-2 off).

## Unresolved candidate lemmas
H-r general; global optimizer/parity classification at level 4; odd- and Möbius-phase level-4
positivity; bi-conic Schwarzian sign; 2×2 α=0 step; ellipse-squared analytic proof (elementary);
contact-degeneracy write-up; D2 sharp-class question (open classical problem, NOT needed).

## Best numerical assets (experiments/)
extremal_pullback.py + theodorsen.py + minkowski_test.best_extremal = exact extremal-pair machinery
(certified; best_extremal is the phase-safe solver — pitfall P5). ellipse_sandbox.py = exact 2×2.
slice_exact.py = exact elliptic sym4 slice (25 dps; family-restricted — cross-check phases).
slice_invariants.py = L15 invariant formulas (K closed form 1e-16, law 1e-10, frame identities).
sym3_sweep_s51.jsonl (40 rec); sym4_sweep_s61.jsonl (20 rec; OTHER-branch taus/f0e BUGGY, ρ/K/qs
fine). Hr_adversarial floors: n3 +2e-4, n4 +1.7e-2, n5 +1.6e-2 (n6 dense: NOT OBTAINED — job killed
as futile, see pitfall P7 + COUNTEREXAMPLE_SEARCH; use structured families instead).
Direct ratio searches: n=6 → 1.148, n=7 → 1.465.

## Epoch-5 late findings (2026-07-21, after Landen theorem)
- Landen theorem BANKED (proof/landen_theorem.md): sym3 ρ = 1 − π/(2K(k₁)) classical-complete.
- Elliptic sym4 slice b_j = c·a_j: W(A) exact ellipse CONFOCAL WITH OUTER EIGENVALUES ±e₁
  (τ₁ = √k focus-image law persists); inner pair ±e₂ interior nodes.
- Level-4 stationarity: c₁·∂_αB(τ₁)/e₁ + c₂·∂_αB(τ₂)/e₂ = 0 with frame coefficients
  c_j = (h_jx₀)⟨g_j,y₀⟩; NEAR-balance c₁ ≈ −c₂ (⟺ ⟨Ax₀,u₀⟩ ≈ 0 via L13 moment condition,
  0.2–3%) but NOT exact — frame-free law loses 2e-5–6e-4 in K. Exact closed form must carry the
  algebraic frame system (heavier than sym3 but finite).
- Tools: experiments/slice_exact.py (exact elliptic slice machinery, 25 dps).

## Epoch-5 LATE additions (2026-07-21, second work block)
- L15 (nodal/Pick closed form, odd law) is proved; L16 gives the even-phase midpoint stationary
  point + conditional q=1/2 identity. EL4 = ρ = 1−Θ(k,U₂) is now proved separately by the
  Schwarzian theorem. `proof/slice_closed_form.md` is the structural master note.
- Slice phase geography: odd/even/Möbius(=odd α→1 endpoint). Sweep OTHER-branch taus/f0e
  fields are BUGGY — recompute, don't trust.
- D2 soft-class hunt CLOSED (all falsified, incl. odd+G'↑ via certified 40-dps
  counterexample); symmetric-node case + wedge family + convex trace bound PROVED.
  The former deformation route was superseded by the Schwarzian proof of EL4; see
  proof/D2_landscape.md and proof/el4_schwarzian_theorem.md.

## Next five concrete actions (refreshed 2026-07-22, Epoch 6)
1. **Finish normal-face uniformity and merge L93 into L86.** L94--L98 close every bounded
   weighted chart at a nonnormal terminal block. Treat `a->0` onto L88's exact normal
   direct-sum manifold, then lift the transverse estimate through the metric flag and couple
   it to the negative Gram/strong-mode and losing-sector terms. Do not grind larger Schur
   matrices.
2. **Full CP-correction moments.** Do not retry scalar shifts or positive-state scalarizations.
   Test whether the L21 trace inequality follows from the block-Toeplitz positivity already
   supplied by the full operator-valued Crouzeix--Palencia correction.
3. **Shifted degree-one Möbius H-r phase**: derive its stationarity/rho formula and prove
   positivity (or K≤2) independently as the scalar-conjecture fallback.
4. **Bi-conic Schwarzian route**: compute `SG` for the off-slice collapsed map on its critical
   real interval. If nonnegative, L17 gives D2 immediately; otherwise compare its Sturm potential
   directly with `−1`.
5. **Odd-phase level-4 positivity**: ρ = B₁g₁q₁+B₂g₂q₂ ≥ 0 under the L15 stationarity law
   (slice_closed_form.md §3; interlacing τ₂ < α < τ₁, dominant positive outer term observed).
   Try the same kernel/deformation machinery; the Möbius-equality structure should persist.

## Paste-ready continuation instruction
"Continue the Crouzeix campaign in /home/liam/Downloads/crouzeix (git repo; commit+push after each
task). Read RESEARCH_STATE.md (NEWEST section first), then proof/slice_similarity_duality.md,
proof/slice_boundary_theorems.md,
proof/el4_schwarzian_theorem.md,
proof/even_pick_globality.md, proof/slice_closed_form.md, and proof/D2_landscape.md; master program in
proof/rho_positivity_program.md. Resume at restart-packet action 1: turn the completed repeated-
Crabb fixed-direction/metric-flag theory (L61--L93) into a uniform repeated-`C3`
neighbourhood theorem. L98 already closes every bounded weighted chart off the normal terminal
face; resume at `a->0` and the flag/strong-gap merger. Do not start a larger copy-space Schur
grind: L93 reduces every terminal irreducible flat block to size at most two. The complete elliptic 4×4 slice is
proved by L59; do not start a 5×5 slice grind or redo its SDP duality, modal norms, projective
charts, or certificate. EL4 and even-sector midpoint globality are proved; definite parity is
false. Respect
APPROACH_LEDGER.md pitfalls P1–P8: every
numerical claim needs the certificate battery; treat any apparent violation as artifact until it
survives strict re-evaluation and an independent implementation; cross-check extremal phases with
best_extremal (P5); verify analyticity/univalence of probe map families (P6). Do not re-open dead
ends listed in the ledgers (esp. soft D2 classes — all falsified with certificates). The
stop-condition remains: rigorous general proof or certified counterexample."
