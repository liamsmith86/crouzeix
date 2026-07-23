# RESEARCH_STATE.md — Crouzeix Conjecture Campaign

**Last updated:** 2026-07-22 (Epoch 6 — L100 weighted normal center)

## NEWEST (2026-07-22): L100 locates the sharp weighted normal center
- On the normal terminal face `a=z=D1=0`, the diagonal traceless third support has coefficients
  `[q³]h=d(3d²−8sqrt(2)conj(w))/128` and
  `[q¹]h=−conj(d)(3d²−8sqrt(2)conj(w))/128`.
- Thus its unique center for `d!=0` is
  `w=3conj(d)²/(8sqrt(2))`.  An initial `d²` guess failed exact regeneration; the corrected
  conjugated formula is now checked by full symbolic substitution.
- At the center the cubic support vanishes on the exact normal direct-sum manifold, but every
  transverse nonnormal edge still has
  `E3≤−a(a²+2|d|²)I/4`, independent of the weighted coordinates.
- This supplies the correct normal/tangential chart for the remaining compactness argument.
  It does not yet control analytic remainders when the transverse edge appears at a later
  asymptotic scale. `proof/repeated_p3_flat_two_copy_weighted.md`.

## NEWEST (2026-07-22): L99 removes cubic degeneration at the normal face
- L94 used an off-diagonal third-support coefficient and therefore weakened as the Schur edge
  `a->0`.  The missing normal-face coefficient is diagonal:
  `[q³](Q3)00=3d³/128`.
- Combining both coefficients gives
  `m3≥(2|d|²+|a|²)^(3/2)/(192sqrt(3))`, hence a cubic endpoint at most
  `−(2|d|²+|a|²)^(3/2)/(12sqrt(3))` on every nonzero trace-zero block.
- Thus the unweighted terminal cubic is uniform even at `a=0`.  There the gain is caused by the
  convex hull of the two oppositely deformed normal summands, consistent with L88.
  The remaining normal-face task is only the weighted recentering/compactness patch when a
  later nonnormal transverse scale is compared with the exact normal manifold.
  `proof/repeated_p3_flat_two_copy.md`.

## NEWEST (2026-07-22): L98 closes the bounded weighted terminal chart
- Start with a trace-zero nonnormal terminal block `D` at order `epsilon`, then let an arbitrary
  tracezero tangent `D1`, scalar trace `z`, and common flat mode `w` enter at order `epsilon²`.
  The complete order-two effective support remains scalar.
- Exact order-three propagation, including the inverse-map Fréchet term and arbitrary first
  conformal coefficient, collapses to
  `16(mean Q3−mean(lambda_max(Q3))I)`.
- If `z=0`, the mean vanishes and L94's uncancellable cross mode makes the endpoint strict.
  If `z!=0`, `mean Q3` and its mode-two coefficient have commutator
  `−15z²[D*,D]/4096`, so a nonnormal `D` cannot have a common top line.
- Differentiating the canonical metric supplies the free second-metric block
  `−3sqrt(2)D1/8`, and the same endpoint identity survives exactly.  Therefore every bounded
  second-order flat recentering of a nonnormal two-copy terminal block is strictly descending.
  Only the degeneration `a->0` onto L88's exact normal direct-sum stratum remains.
  `proof/repeated_p3_flat_two_copy_weighted.md`.

## NEWEST (2026-07-22): L97 rules out trace/common cancellation
- The trace-driven two-copy Pauli support is even under `q->−q`; the common-`w` support is odd.
  Pairing the two boundary points makes the joint mean norm dominate each component, so adding
  `w` can never weaken L95's trace-splitting Jensen gap.
- Combining that fact with L96 in two quantitative regimes gives a joint endpoint bound
  proportional to
  `−|w|a³(4|d|²+a²)/(2|d|²+a²)²`, uniformly in the scalar trace coordinate `z`.
- L94--L97 now control every leading terminal interaction (`z`, `w`, and the cubic block) with
  common zero only on the exact normal stratum.  The next task is no longer another leading
  coefficient: it is analytic remainder absorption/recentering onto L88's normal direct-sum
  manifold, then lifting the estimate through the L93 flag.  `proof/repeated_p3_flat_two_copy.md`.

## NEWEST (2026-07-22): L96 quantifies the common-mode terminal transition
- On the trace-zero terminal block, the non-scalar common-mode support has zero mean and upper
  cross entry `sqrt(2)a(q^−1w−q³conj(w))/16`.
- Its mean top eigenvalue dominates either Fourier coefficient, giving the exact canonical
  endpoint bound `lambda_max(E2,w)≤−sqrt(2)|aw|`.
- L94--L96 now give coercive margins for the terminal cubic, trace split, and common mode, all
  with common zero only at the exact normal face.  The remaining two-copy issue is simultaneous
  nonzero `z,w`, where even/odd support modes can interact, and then uniform analytic remainder
  absorption.  `proof/repeated_p3_flat_two_copy.md`.

## NEWEST (2026-07-22): L95 quantifies the second/cubic terminal transition
- Write a two-copy matrix as `Z=zI+[[d,a],[0,−d]]`.  The non-scalar second
  support is a traceless Hermitian `2 x 2` curve, hence a Euclidean
  Pauli-vector curve `B+C(q)`.
- Exact Laurent averaging factors its angular energy perpendicular to the mean direction as
  `9|z|⁴a²(4|d|²+a²)/(8192b0)`.  The elementary Euclidean triangle-defect identity turns this
  into an explicit negative second endpoint whenever both `z` and `a` are nonzero.
- L95 controls departure from the trace-zero face; L94 controls the trace-zero face cubically.
  Their common zero `a=0` is exactly the normal stratum already controlled by L88/L73.
  The next missing quantitative coordinate is the common flat mode `w`, followed by analytic
  remainder absorption.  `proof/repeated_p3_flat_two_copy.md`.

## NEWEST (2026-07-22): L94 quantifies the two-copy terminal descent
- On L93's nonnormal terminal block `Z=[[d,a],[0,−d]]`, the L92 third-support
  cross entry has an uncancellable Fourier coefficient
  `a(2|a|²+4|d|²)/128`.
- Mean top eigenvalue dominates that coefficient, so the exact cubic metric endpoint obeys
  `lambda_max(E3)≤−|a|(|a|²+2|d|²)/4`.
- Thus the cubic margin is uniformly coercive away from the exact normal face `a=0`.
  The next weighted step is to combine this with the second-order flag/Jensen margin as a
  block approaches the normal or reducible strata.  `proof/repeated_p3_flat_two_copy.md`.

## NEWEST (2026-07-22): L93 closes every fixed flat-copy direction in arbitrary multiplicity
- The free L80 metric tangent extends from a pure star to an arbitrary copy projection `P`.
  Its exact endpoint compression is the negative Gram sum
  `−5sqrt(2)PZ(I−P)Z*P−15sqrt(2)PZ*(I−P)ZP/4`; all internal complement
  blocks cancel.
- Iterating the derivative kernel produces a descending flag.  On an irreducible copy block it
  either ends at zero, giving strict second-order descent, or fills the block.  In the latter
  case the full effective support is scalar:
  `Z²=alpha I`, `ZZ*+Z*Z=beta I`.
- The scalar-support relations form a two-generator Clifford algebra.  The nilpotent case is
  a direct sum of equal square-zero pairs; the invertible case is a direct sum of one- or
  two-dimensional anticommuting-symmetry blocks.  Hence every irreducible residual has size at
  most two and is already closed by L88 or L92.
- This removes the proposed `4 x 4`, `5 x 5`, ... Schur grind and closes all **fixed** directions
  in the flat copy core.  The remaining issue is genuinely uniform: weighted sequences where
  Jensen/flag gaps, `w`, strong normals, losing gaps, and the cubic pair margin collapse
  together.  `proof/repeated_p3_flat_metric_flag.md`.

## NEWEST (2026-07-22): L92 closes every fixed two-copy flat direction
- A nonnormal `2 x 2` Schur matrix can be second-order flat at `w=0` only when `tr Z=0`.
  Then `Z²` and `ZZ*+Z*Z` are scalar, so the second support is scalar.
- The third effective support is traceless with an explicit cross Laurent polynomial whose
  middle coefficient contains `2|a|²+4|d|²`; it is nonzero whenever the nonnormal edge is.
  Refactored exact third-metric propagation gives strict endpoint
  `−16 mean(lambda_max(Q3))I`.  This extends L77 to arbitrary opposite diagonal loops.
- Normal cases are exact by L88; nonnormal `w!=0` cases are strict by L91; nonzero-trace
  `w=0` cases have strict Jensen gap.  Thus no fixed two-copy flat direction remains
  unclassified.  Only weighted uniformity and larger nonnormal copy matrices remain.
  `proof/repeated_p3_flat_two_copy.md`.

## NEWEST (2026-07-22): L91 eliminates irreducible equality when `w!=0`
- The unique non-scalar Fourier modes `+3,-3` in L90 are nonzero multiples of `Z,Z*`.
  Any vector that is an eigenvector of the support family for every angle must therefore reduce
  `Z`.
- A higher-dimensional unitarily irreducible nonnormal block with `w!=0` has strict
  matrix-Jensen descent.  Minimal reducing blocks that can carry equality are one-dimensional,
  hence are precisely nearby single `C3` blocks controlled by L73 and domain monotonicity.
- The common flat mode creates no new irreducible equality mechanism.  The remaining issue is
  the weighted transition `w->0`, centered on the `w=0` nonnormal copy matrix treated at three
  copies by L89.  `proof/repeated_p3_common_mode_rigidity.md`.

## NEWEST (2026-07-22): L90 inserts the common flat mode exactly
- With common `W(w)=w(E10+E21)`, the full effective support is
  `Q_Z(q)+r_w(q)I+ell_w(q)Z+conj(ell_w(q))Z*`, with explicit Laurent scalars.
  The first support compression remains scalar.
- The `r_w I` term cannot change any eigenspace.  Every interaction between the common flat
  mode, Schur diagonal, and nonnormal upper part is now the finite Hermitian pencil
  `ell_w Z+conj(ell_w)Z*` plus L88's quadratic four-term polynomial.
- This makes the next target an exact simultaneous-eigenvector classification for the Laurent
  coefficients of a `3 x 3` upper-triangular `Z`, extending L89 beyond zero diagonal.  No
  numerical Riemann-map fit is needed.  `proof/repeated_p3_flat_common_support.md`.

## NEWEST (2026-07-22): L89 closes the first nonnormal copy-matrix stratum
- For `N=[[0,a,b],[0,0,c],[0,0,0]]`, the angular coefficients in L88 force any common top
  vector (when `ac!=0`) to be the middle coordinate; its constant coefficient is compatible
  only when `b=0`.  The complement top defect then has determinant `16|ac|²>0`.
- Thus endpoint equality requires `abc=0`.  The length-two path `b=0,ac!=0` is exactly L80's
  rank-two star and becomes strict at optimized second order.  The square-zero face `ac=0`
  unitarily reduces to one L77 pure pair and is strict at cubic order.  If `abc!=0`, the
  matrix-Jensen term is already strict.
- This completely closes the zero-diagonal, three-copy nonnormal flat core.  The next Schur
  target is to allow nonzero diagonal/common `w`, then prove the same path/star recursion in
  arbitrary copy size.  `proof/repeated_p3_flat_three_copy.md`.

## NEWEST (2026-07-22): L88 packages the flat core into one copy matrix
- Every generator-zero/two edge and generator-zero loop on the L87 center is exactly
  `Z* tensor X0+Z tensor Y0` for one arbitrary copy-space matrix `Z`.  Its effective support is
  `5(ZZ*+Z*Z)/128−3(q²Z²+q^−2(Z*)²)/128`.
- If `Z` is normal, a copy unitary turns the full repeated matrix (including the common flat
  `w` mode) into a direct sum of nearby single `C3` blocks.  L73 controls each block, and
  complete-spectral-set monotonicity transfers the bound to the direct sum's convex-hull
  numerical range.  This is an exact arbitrary-multiplicity local theorem on the normal stratum.
- Schur form now reduces the unresolved center to a strictly upper-triangular, generator-zero
  copy matrix.  The next target is to prove that its common-top strata reduce to L77's pair
  cubic or L80's rank-two star improvement; otherwise (4)'s matrix-Jensen gap is already strict.
  `proof/repeated_p3_flat_copy_matrix.md`.

## NEWEST (2026-07-22): L87 isolates the true higher-order center
- The L86 endpoint is a sum of six explicit negative PSD terms.  A null vector must be a
  common top eigenvector of every `Q(q)`, lie in `ker(H1)` and every `ker(Aj*)`, and have
  common strong coordinate `v=0`.
- After a copy-unitary, the selected branch has no first-order losing coupling and no
  generator-one coupling.  Only its generator-zero/two flat star and the common flat
  single-block modes remain.  These are exactly the previously proved L77/L80 and L67--L73
  boundary mechanisms.
- The uniform blow-up therefore needs to center only on this intersection stratum.  Generic
  graphs, winner--loser edges, generator-one data, and the common strong mode already carry a
  strict first/second-order margin and should be treated as normal variables, not expanded
  blindly to higher order.  `proof/repeated_p3_equality_reduction.md`.

## NEWEST (2026-07-22): L86 closes every fixed repeated-`C3` direction through second order
- Promote all common kernel vectors into a maximal winner sector.  The complementary losing
  mean compression is then strictly negative, so L61/L82 provide first-order Slater slack for
  every losing internal block, even with pointwise contacts.
- A two-winner/one-loser exact audit polarizes L81 into
  `−25A0A0*/8−8A1A1*−50A2A2*/9`.  Adding L85 gives the total winner endpoint as its matrix
  Jensen term, the `H1` square, the common `v` square, and these three negative Gram matrices.
- Consequently a positive Jensen gap descends at first order, while every fixed direction on
  its zero face is nonpositive at second order.  The remaining repeated-`C3` problem is uniform
  higher-order control when all these second-order terms and the losing mean gap collapse
  together—not an unclassified fixed first-order direction.
  `proof/repeated_p3_complete_fixed_direction.md`.

## NEWEST (2026-07-22): L85 closes the full tied-winner sector at second order
- Global affine normalization plus L84 leaves only two new common diagonal modes:
  `w(E10+E21)` and `vE20`.  The first inverse-map coefficient vanishes at `C3` but its
  Frechet derivative does not; retaining it is essential.
- The exact complete endpoint is
  `16(mean Q−mean(lambda_max Q)I)−8H1²−21|v|²I/4<=0`.
  The common flat `w` mode cancels against the normal-angle correction, while the `v`
  coefficient agrees with the independently proved single-block L63 curvature.
- L74, L84, and the two common modes exhaust the fixed tied-winner first-order quotient.
  What remains is not another missing direction class: it is higher-order optimization on the
  endpoint equality set and a uniform merger with L82's first-order losing sectors.
  `proof/repeated_p3_complete_winner_sign.md`.

## NEWEST (2026-07-22): L84 adds every relative diagonal winner motion
- The diagonal zero-support condition has real rank seven.  Its 11-dimensional kernel is the
  rank-eight within-copy unitary orbit plus exactly one complex generator-zero loop and one
  real generator-one loop.  This is a complete quotient classification, not an ansatz.
- Adding arbitrary canonical loops to every vertex of L83's full graph leaves the exact endpoint
  identity unchanged: `E=16(mean Q−mean(lambda_max Q)I)−8H1²<=0`, with each real loop placed
  on `diag(H1)`.  The symbolic triangle includes all loop/edge mixed paths and retains the
  first conformal mode until it cancels.
- The important remaining distinction is a **common** diagonal motion shared by all copies.
  It cannot be discarded as gauge: L84 controls the relative zero-support remainder after a
  reference block is removed, not its coupling to that shared single-block motion.  L85 now
  closes that coupling; graph equality optimization and weighted winner/loser merging remain.
  `proof/repeated_p3_winner_diagonal_sign.md`.

## NEWEST (2026-07-22): L83 proves the full tied-winner graph sign
- Put arbitrary L74 quotient data on every edge of a tied winner sector and let `Q(q)` be its
  complete second effective support matrix.  Summing the pairwise L76 metric blocks and
  propagating the full second metric gives
  `E=16(mean Q−mean(lambda_max Q)I)−8H1²`, where `H1` is the Hermitian copy matrix of
  generator-1 edge coefficients.  Matrix Jensen order and `H1²>=0` prove `E<=0` in every
  multiplicity; the possibly nonzero first conformal mode cancels exactly.
- The symbolic three-copy triangle exposed a false simplification before it entered the proof:
  on cyclic graphs, `mean Q` is not the naive weighted sum of the three coefficient-matrix
  squares.  Oriented generator-0/2 mixed paths survive.  L83 uses the actual support mean, so
  the sign is unaffected; the shortcut is explicitly logged as invalid.
- Five random nonlinear full-triangle cases have stable strict quadratic coefficients from
  about `-1.48` to `-4.20`, but equality of the explicit endpoint is not yet the optimized-SDP
  classification.  Tied diagonal residuals and uniform winner/loser merging also remain.
  `proof/repeated_p3_winner_graph_sign.md`;
  `experiments/repeated_p3_winner_graph_sign.py`.

## NEWEST (2026-07-22): L82 extends radial gaps to every unique-winner face
- Normalize by the selected copy's first boundary motion.  The losing copy-space support
  compression is then a continuous matrix function `B_-(q)<=0`.  If it has no vector in every
  pointwise kernel, its mean is negative definite.  L61's exact tangent dual plus its Slater
  construction therefore provide strict first-order metric slack on the whole losing sector,
  even when the pointwise gap touches zero at isolated angles.
- Copy-sector parity makes the selected second metric additive: its diagonal contribution is the
  single-block coefficient `e_selected<=0` from L63/L65, and its winner/loser star contribution
  is exactly L81's universal
  `−25||a0||²/8−8||a1||²−50||a2||²/9`.  Internal losing-space blocks are absorbed by the first
  slack and do not enter the selected endpoint.
- Hence every fixed unique-winner direction with nonzero star coupling descends strictly at
  second order, for arbitrary nonconstant gap shape.  The scalar touching gap `1-cos(theta)`
  numerically converges to all three exact L81 coefficients.  The remaining repeated-face issue
  is multiple common winners and a uniform estimate as losing gaps collapse.
  `proof/repeated_p3_unique_winner.md`.

## NEWEST (2026-07-22): L81 proves the first diagonal--cross repeated theorem
- In the clean strict common-maximizer model, the selected copy grows as `(1+d epsilon)C3` and
  the other copies shrink as `(1-d epsilon)C3`, while the off-diagonal blocks are arbitrary L74
  star directions.  A diagonal metric tangent `diag(d,0,-d)` makes each losing copy strict at
  first order, with Stein slack `diag(6d,15d)`.
- The selected second-order problem has two free cross-metric scalars.  Exact elimination gives
  independent quadratics whose minimizers are `x=sqrt(2)/4`, `y=4sqrt(2)/3`, producing
  `e=−25||a0||²/8−8||a1||²−50||a2||²/9<0`.  This matches nonlinear limits on all three pure
  generators and extends by direct summation to arbitrary star multiplicity.
- The coefficient is independent of fixed `d>0`, but the valid asymptotic neighbourhood shrinks
  as `d->0`; that singular weighted transition must be joined to L77/L80 rather than treated as
  a uniform gap theorem.  General diagonal shapes, internal losing-copy blocks, and multiple
  first-order winners remain.  `proof/repeated_p3_radial_gap.md`;
  `experiments/repeated_p3_radial_gap.py`.

## NEWEST (2026-07-22): L80 closes every nonzero pure star ray
- The nonlinear rank-two probe did more than reject fourth order: it identified a missing free
  first-metric tangent.  On the active span of independent `a0,a2`, use
  `U_tau=[[0,(-3sqrt(2)/8+tau)a0^T],[(1/sqrt(2)+tau)conj(a2),0]]`.
- Exact full metric propagation gives
  `E(0)=diag(0,-5 adj(H)/72)` and
  `E'(0)=(5sqrt(2)/3)diag(-trace K,K)`, where
  `H=9 conj(a0)a0^T+16 conj(a2)a2^T` and
  `K=3 conj(a0)a0^T+4 conj(a2)a2^T`.  Rank two makes the orthogonal block at zero strictly
  negative, while the selected derivative is strictly negative.  Hence `E(tau)<0` for all
  sufficiently small positive `tau`.
- Combined star classification: `a1!=0` is strict at second order by L78/L79; `a1=0` with
  rank-two `(a0,a2)` is strict at second order by L80; rank one is copy-unitarily L77 and strict
  at cubic order.  Thus every nonzero pure star ray at arbitrary multiplicity descends.  This is
  still not a repeated-block neighbourhood because diagonal and internal `y^perp` blocks and
  weighted mixtures remain.  `proof/repeated_p3_star_second_sign.md` §6;
  `experiments/repeated_p3_star_second_sign.py`.

## NEWEST (2026-07-22): L79 classifies only certificate equality; fourth-order route rejected
- The L78 effective support has block form `diag(s(q),R(q))` and the exact trace reversal
  `s(q)=trace R(-q)`.  Combining this with the endpoint's two NSD summands proves that its top
  eigenvalue is zero **exactly** when the generator-1 coefficient vector `a1` vanishes; every
  `a1!=0` star direction is already strict at second order.
- When `a1=0`, `s(q)=trace R(q)>=lambda_max R(q)`, so the selected copy is a common top branch.
  A fixed orthogonal-copy branch ties it exactly when `span{a0,a2}` has rank at most one.  That
  rank-one case is copy-unitarily the L77 pure pair and descends cubically.
- Genuine rank-two `(a0,a2)` data require at least two orthogonal copies and have only the
  selected common top support branch.  However, this classifies equality of the **particular
  L78 metric**, not the optimized second-order SDP.  A nonlinear probe on orthogonal unit data
  gave `(t_*-4)/epsilon² -> approximately -.566`, rejecting the fourth-order inference; L80 now
  proves the missing strict second-order tangent exactly.  `proof/repeated_p3_star_second_sign.md`.

## NEWEST (2026-07-22): L78 closes arbitrary star coupling at second order
- For `m` repeated `C3` copies, collect the three L74 generator coefficients across the
  `m-1` orthogonal copies into vectors `a0,a1,a2`.  If `Q(q)` is the resulting second effective
  support matrix, exact averaging gives `Qbar=diag(mu,G)`, where
  `G=5 conj(a0)a0^T/128+conj(a1)a1^T/4+5 conj(a2)a2^T/72` and `mu=trace G`.
- A falsification check caught that `lambda_max Q(q)` is generally not pi-periodic once `m>2`;
  its first conformal Fourier mode can be nonzero.  Retaining that `A0²` correction in the full
  second metric shows that it cancels exactly.  The upper endpoint is
  `16(Qbar−mean(lambda_max Q) I)−8 diag(||a1||²,conj(a1)a1^T)`.
- Pointwise matrix order `Q(q)<=lambda_max(Q(q))I`, followed by averaging, makes the first term
  NSD; the second is a negative Gram block.  Thus simultaneous star coupling has nonpositive
  second-order change for every multiplicity.  Higher order on its equality set and all
  diagonal/internal-orthogonal-copy mixtures remain.  The exact two-symbolic-copy polarization
  audit runs in under two seconds.  `proof/repeated_p3_star_second_sign.md`;
  `experiments/repeated_p3_star_second_sign.py`.

## NEWEST (2026-07-22): L77 closes the pure-pair flat plane cubically
- On L76's equality plane `alpha1=0`, the second effective support matrix remains scalar.  The
  exact third effective matrix is off-diagonal with entry `c(q)`, an explicit odd Laurent
  polynomial, so the top third support branch is `|c(q)|`.  Its pi-periodicity removes the first
  Fourier mode, leaving only `m3 A0` after the third Schwarz map is evaluated at `C3`.
- The full third metric is constructed recursively with the cubic lower, Stein, and upper Schur
  penalties.  All second-map, Frechet, and lower-metric terms cancel at the upper endpoint:
  `e3=-16 mean_{|q|=1}|c(q)|`.  The factorization of `c` shows that this is strictly negative
  for every nonzero `(alpha0,alpha2)`.  On the pure axes it is exactly `-|alpha0|³/4` and
  `-16|alpha2|³/27`.
- Therefore every nonzero pure multiplicity-two cross-pair direction descends, quadratically
  off the L76 plane and cubically on it.  L78 now controls simultaneous orthogonal-copy
  directions at second order, but not their equality set, diagonal/cross mixtures, or the
  common-maximizer top-order inequality.
  `proof/repeated_p3_third_sign.md`; `experiments/repeated_p3_third_sign.py`.

## NEWEST (2026-07-22): L76 proves the pure cross-pair second-order sign
- On the three-complex-parameter L74 quotient for one selected/orthogonal copy pair, an explicit
  Hermitian first metric tangent annihilates all lower, upper, and Stein active compressions.
  The complete second metric is then built by the two-level Crabb Stein recurrence, including
  both endpoint penalties and the contraction Schur penalty.
- The upper endpoint collapses exactly to a scalar.  Substituting L75's conformal mean gives the
  feasible coefficient
  `e=-8|alpha1|²-(4sqrt(2)/(3pi))|3alpha0 conj(alpha1)+4alpha1 conj(alpha2)|<=0`.
  Thus every pure cross-pair direction is nonincreasing at second order, without a discretized
  boundary or numerical SDP.  The full symbolic checker reconstructs the second metric and
  verifies all three Schur complements exactly in under one second.
- Equality is exactly `alpha1=0`, leaving a two-complex-dimensional flat plane; L77 now proves
  strict cubic descent on it.  This is not a repeated-block neighbourhood theorem: mixtures with
  diagonal single-copy directions, common-maximizer order inequalities, and simultaneous
  multiplicity directions remain.  `proof/repeated_p3_stein_sign.md`;
  `experiments/repeated_p3_stein_sign.py`.

## NEWEST (2026-07-22): L75 removes the repeated-block boundary calculation
- For one multiplicity-two L74 cross pair, the first support compression vanishes.  The exact
  reduced support resolvent `I-H/4-3H²/4` makes the second effective copy-space matrix diagonal.
  Its largest eigenvalue is
  `kappa(q)=M(q)+(sqrt(2)/24)|Re(q d)|`, where `M` has only Fourier modes zero and two and
  `d=3 alpha0 conj(alpha1)+4 alpha1 conj(alpha2)`.  The absolute cosine is the complete
  nonsmooth repeated-eigenvalue effect.
- `kappa` is pi-periodic, so its first Fourier coefficient is zero.  Since `C3³=0`, every
  nonconstant Fourier mode then disappears under functional calculus.  The pulled operator is
  exactly `T_e=A0+eE-e² kappa_hat(0) A0+o(e²)`, with
  `kappa_hat(0)=5|alpha0|²/128+|alpha1|²/4+5|alpha2|²/72
  +sqrt(2)|d|/(12pi)`.
- Thus no boundary discretization remains.  L76 now solves the finite `6x6` Stein sign exactly:
  generator 1 decreases, while the complete `(alpha0,alpha2)` plane is second-order flat.
  `proof/repeated_p3_second_support.md`; `experiments/repeated_p3_second_support.py`.

## NEWEST (2026-07-22): L74 reduces the repeated-block exceptional cross face
- On L61's zero-Jensen face choose the common maximizing copy vector `y`.  For each
  `xi perpendicular to y`, the cross blocks `(X,Y)=(E_{xi y},E_{y xi})` obey one exact
  Laurent identity against the Crabb top support vector.  Its real constraint rank is 14, so
  the kernel has dimension 22.
- Infinitesimal cross-copy unitary mixing has rank 16 inside that kernel.  Exact
  Hilbert--Schmidt quotienting leaves only **three complex directions per orthogonal copy**,
  with sparse canonical generators.  Thus multiplicity `m` contributes `3(m-1)` complex
  cross parameters, not two arbitrary `3x3` blocks.
- The quotient is nonzero, so repeated blocks do not reduce to L73 at first order.  The active
  task is their second-order effective support/Stein sign, first when `y` has a strict compressed
  top gap and then at nonsmooth ties.  Proof: `proof/repeated_p3_common_maximizer.md`; exact audit:
  `experiments/repeated_p3_common_maximizer.py`.

## NEWEST (2026-07-22): L73 closes a full neighbourhood of the single `C3` block
- The local rank-one Stein condition is an honest real-analytic feasible certificate: its defect
  Hessian at `C3` is `8|x|²+(8/3)|y|²>0`, so the implicit-function theorem gives a
  unique analytic locally optimized defect.  Analyticity of the numerical-range Riemann map
  follows here from the simple uniform support eigenvalue and the standard near-circle boundary
  Fourier/implicit-function argument (consistent with Rodin 1986 and Wu 1993).
- In the L69 slice, L71--L72 make the exact disk curve an ambient critical manifold.  The
  `(s,V)` Hessian is strictly negative, so those three real normal variables can be maximized out
  analytically.  The last complex soft germ is equivariant under `(z,U)->(e^{it}z,e^{2it}U)`.
  Stationarity removes soft degree zero and one; symmetry leaves only `|z|²|U|²` and
  `|U|⁴` at degree four.  L70 and L68 give their exact coefficients `-25/56` and `-4`.
  Every higher allowed monomial is an absorbable small multiple of these negative terms.
- Therefore the certificate is at most four on a full seven-real-dimensional slice neighbourhood;
  L69 lifts this to every complex `3×3` matrix in a full neighbourhood of `C3`.  Hence its
  numerical range is a **complete `2`-spectral set** there.  This is not a repeated-block,
  larger-size, or general theorem.  Targeted literature searches found no prior full-neighbourhood
  result; call it apparently new pending publication-level review.
- Exact audit: `experiments/p3_disk_morse_bott.py`; proof:
  `proof/p3_crabb_local_theorem.md`.  A non-load-bearing 18-case mixed/superweighted map probe
  reached at most `3.999999999999` with diagnostics below `9.8e-13`.

## NEWEST (2026-07-22): L72 proves exact stationarity on the disk center
- The L71 curve has an exact canonical Schur form
  `T_l=[[0,a,-2l],[0,l,a],[0,0,0]]`, `a²=2(1-l²)`.  This is not inferred from
  the disk property alone: a new exact identity shows that the normalized product of its two
  nonzero squared singular values is four, which selects the symmetric Schur subfamily.
- `P=diag(1,2,4)` satisfies the rank-one Stein identity `P-T_l* P T_l=e1 e1*` on
  the whole curve.  Exact linearization makes the condition derivative `4ℓ(Re G)`, independent
  of the rank-one-defect adjustment.  An exact three-pole support residue calculation proves that
  the first Riemann-map correction has precisely the same `ℓ` value for every complex ambient
  perturbation.  Therefore every first variation cancels.
- This establishes the critical-manifold half of the proposed weighted Morse--Bott argument.
  It does **not** establish the local inequality: the active task is a uniform negative normal
  Hessian (and analytic optimized-defect selection) near the curve.  The exact audit runs in
  about one second: `experiments/p3_disk_center_tangent.py`; proof:
  `proof/p3_disk_center_tangent.md`.

## NEWEST (2026-07-22): L71 identifies and closes L70's hidden center
- The recentered L70 jet is not an accidental sequence of cancellations.  It is the analytic
  root through `R=1` of `81ε⁴R²+(1152ε²−4096)R+4096=0`, with
  `u=(3√2/64)R`, `v=−(9/64)R`, `s=0`.  Its expansion starts
  `R=1+9ε²/32+405ε⁴/4096+...`, exactly reproducing the independent center shift.
- Exact Kippenhahn reduction proves that every matrix on this curve has a circular numerical
  range centered at a double eigenvalue: the homogeneous polynomial is
  `(z−2cx)((z+cx)²−r²(x²+y²))`, and the isolated point lies inside the circle near `C3`.
  Berger--Okubo--Ando therefore gives the complete L21 bound `t*≤4` on the whole curve.
- This closes the center itself but not yet a full neighbourhood.  The remaining `p=3` task is
  a uniform normal estimate for sequences approaching the disk curve faster than L70's leading
  weighted scale, preferably via an analytic weighted Morse--Bott/splitting argument rather
  than still higher jets.
- The disk theorem is classical.  The new campaign contribution is the exact identification
  of the hidden weighted center with that classical locus.

## NEWEST (2026-07-22): L70 closes the weighted leading sign at `C3`
- L69 gives an exact seven-real-dimensional affine-unitary normal slice.  In its sharp chart
  `C3+εR1(1)+ε²R2(u)+ε³(sY+vV)`, L70 derives the complete rank-one feasible-certificate
  coefficient
  `H=−8s²−(21/4)|v+39/448+(4√2/7)u|²−(25/56)|u−3√2/64|²`.
  Thus the feared weighted coupling is never positive and has only the center
  `(u,v,s)=(3√2/64,−9/64,0)`.
- The missing exact inputs are the ordinary fifth mode coefficient `−123√2/256`, the bottom
  linear coefficient `−117/128`, and the mode/bottom coupling `−6√2`.  A new sparse exact
  Riemann/Stein engine regenerates the whole leading certificate in about ten seconds, agrees
  with the older implementation through order four, reproduces L67--L68, and cross-checks the
  fifth coefficient through the older Stein engine.
- The last center is genuinely subtle: its fixed-center order-ten descent is cancelled exactly
  by the common recentering factor `1+9ε²/32`; the corrected certificate is flat through order
  twelve, while the remaining real transverse coordinate contributes `−8s²`.  This suggests
  the disk-matrix curve now identified exactly by L71.  A uniform normal estimate is still
  needed; do not infer a punctured-neighbourhood theorem from a finite jet.
- Scope remains the stronger L21 complete-similarity route near one `3×3` Crabb block.  It is
  neither a proof of the general scalar conjecture nor a repeated-block theorem.

## NEWEST (2026-07-22): L20/L59 prove the complete elliptic 4×4 slice
- **The last positive compact tail is closed rigorously.** L59 collects the complete
  determinant in `A=1−a,B=1−b`, retains its shared correlations with `c,X,R,Y`, and uses
  order-0--9 Taylor/Arb order-10 remainder bounds followed by outward Bernstein conversion.
  Ten exactly adjacent ratio-`81/80` rational boxes certify both determinant and all three
  final-minor charts from the old `c+` frontier through `.63`. Bounded recentering rebuilds
  exact rational physical subboxes and certifies them independently; no tolerance is used.
- **The clean forced-regeneration run passed 10/10 and exited zero.** It first regenerated and
  audited both 197,563-record determinant tables, both 207-term corner squares, and the
  Bernstein/recenter machinery. Provenance: commit `3dd51884…`, Python 3.14.6,
  python-flint 0.9.0, NumPy 2.5.1; 130-line log
  `experiments/positive_tail_full_20260722.log`, SHA-256
  `fb79b2dfc652062307d69d26a00d82ff20c4133a0043eb6077458ea7cf70928c`.
- **This proves the theorem for the slice.** Since `.63³−1/4=47/10^6>0`, L59 overlaps L42.
  L27--L29 close the other KKT faces, so `t*(φ(A))≤4` and a condition-two contraction
  similarity hold for every `A=S_a+cS_a^T`, arbitrary positive weights
  `(a1,a2,a3)`, `0<c<1`. Thus its numerical-range ellipse is a complete 2-spectral set.
- **Scope and novelty:** this is not all 4×4 matrices and not the general Crouzeix conjecture.
  The closest-source audit found only Kenan Li's all-dimensional Crabb-derived family with
  fixed weights; no prior arbitrary-weight 4×4 theorem was found. Call L59 apparently new
  pending a publication-level novelty audit.
- **The repeated-Crabb general gate also survives.** A new harness tests block sizes 3/4,
  multiplicities 2/3, three transverse perturbation types, and a `1e-4→1e-2` ladder. Of 120
  perturbed records, 118 pass every map/SDP/support-gap gate and none exceeds four; max
  `3.999844312562`. The two rejected records have bad primal/dual gaps. The strongest cross
  directions remain below four down to `delta=1e-5`, with stable first-order drops
  `(4−t*)/delta≈1.557,1.637`. This is numerical evidence, not a local theorem.
- **General frontier remains active.** Do not begin a 5×5 slice grind. Hartz--McCarthy scalar
  shifts are an exact restatement, not a shortcut; every positive-state scalarization of the CP
  correction is numerically ruled out on one dense 3×3 example. The next analytic target is to
  retain the full correction moments and derive L21's trace inequality from their block-Toeplitz
  positivity, or extract a local inequality from the equality-locus first variation.

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
  `[0,.005]`. The later L56 bridge connects it to L54, L57 supplies the former compact
  frontiers, and L59 completes the cover.
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
- At the L47/L48 stage the exact frontier was (RT); those lemmas removed two complete
  one-parameter sections and supplied the coupled invariant later used by L59.

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
  **At this stage only L27's rank-one/rank-one coupled face remained; L59 later closes it.**
- **L40 now reduces that last face to one matrix-valued inner theorem.** Two spectral
  row-Gram rotations give the transfer
  `R_{a,b}(T)=−S⁻¹(T−A)(I−AT)⁻¹S`; every rank-one/rank-one trace ratio is at most
  `||R_{a,b}(T)||²`. A cancellation-free block formula extends to the closed parameter square.
  In symmetric modal coordinates its node function is a `2×2` rational inner function with
  determinant `(ab−z²)/(1−abz²)`. L59 later proves `||R_{a,b}(T)||≤2` uniformly and hence L20.
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
- **Historical scope audit (superseded independently by L59):** EL4 alone did not prove the
  whole elliptic 4×4 slice. L16 proves midpoint
  stationarity and conditional `q1=q2=1/2`. **The globality debt is now closed by L18**:
  `sup ||F(u1)Q1+F(u2)Q2||=max(1,d||Q1−Q2||)`, and for norm>1 the midpoint automorphism is the
  unique nonconstant global even maximizer up to phase (`proof/even_pick_globality.md`). Thus the
  complete even sector of the elliptic slice has rho≥0. Odd and degree-one phases remain only
  for a phase-by-phase H-r proof; L59 already proves the slice inequality by similarity.
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
  (resid 1e-4–1e-3).** This was the initial signal; L59 now proves the complete slice by the
  independent contraction-similarity route, without resolving each H-r phase separately.
- **GENERAL EQUALITY-LOCUS FIRST ORDER (L61):** at every repeated Crabb disk block, the
  conformal tangent is finite by nilpotence and the L21 tangent SDP has exact value `−4` times
  a support-compression Jensen gap. Hence its value is nonpositive in every direction; under a
  uniform conformal expansion, the upper Dini derivative of `t*` is nonpositive. This is a
  local first-order theorem, not a neighbourhood theorem; equality directions need second
  order and nonsmooth crossings need a regularity argument.
- **SINGLE-CRABB SECOND-ORDER REDUCTION (L62):** first-order complementarity and a general
  second-order PSD Schur lemma reduce the zero-Jensen-gap face to three finite affine block
  LMIs. The second support variation and Schwarz integral give the conformal coefficient in
  finite form. All 12 `p=3,4` test directions have negative quadratic coefficient, stable across
  solvers, map resolutions, analytic/fitted gauges, and fit step. L65 subsequently closes the
  universal single-block sign; repeated common-maximizer faces remain separate.
- **EXACT `3×3` SECOND VARIATION (L63):** eliminating the L62 metric variables gives
  `e₃(E)=−2(Re(E01−E12))²−21|E20|²/4≤0` for every complex perturbation. A regenerating
  18-real-variable symbolic audit proves the identity exactly. Thus the upper second-order Dini
  change of the stronger L21 quantity is nonpositive at the `3×3` Crabb block. The displayed
  equality space still needs third order or an exact-orbit classification.
- **EXACT `4×4` SECOND VARIATION (L64):** the same elimination gives the rank-eight,
  five-square identity (28), so `e₄(E)≤0` for every complex perturbation. A regenerating
  32-real-variable symbolic audit proves it exactly, and all six saved `p=4` SDP cases agree
  within `3.17e-8`. The low-order pattern is now proved at `p=3,4`, but no uniform-`p` theorem or
  complete local-neighbourhood theorem on the equality spaces has been established.
- **ARBITRARY-SIZE CRABB SECOND VARIATION (L65):** the L62 optimum satisfies `e_p(E)≤0` for
  every single Crabb block size. Circle modes reduce to explicit negative kernels whose inverse
  is `(I−H_r)/4+c_kqq*`; the two singular low modes are PSD limits, the bottom modes are
  negative squares, and grade zero is controlled exactly as a weighted shift. The quadratic
  rank is `p(p−2)`. This closes the single-block second-order sign, not the equality directions,
  repeated-block nonsmooth face, or the general conjecture.
- **EQUALITY QUOTIENT (L66):** the `p(p+2)`-dimensional L65 kernel contains an exact
  affine-unitary orbit of dimension `p²+2`; only `2p−2` real directions survive. Circle grading
  makes the residual problem one- or two-complex-dimensional per mode.
- **PURE `p=3` RESIDUAL MODE (L67):** exact formal Riemann-map coefficients through order six,
  followed by a rank-one Stein Gramian, give the feasible condition number
  `4−171ε⁶/4096+O(ε⁷)` on the pure residual mode-one ray. This proves strict local descent on
  that ray without assuming the Gramian is globally optimal.
- **FULL `p=3` QUOTIENT, RAYWISE (L68):** circle-invariant exact calculations give quartic
  `−4|w|⁴−31|z|²|w|²/8` whenever the mode-two coordinate is nonzero. Together with L67, every
  fixed nontrivial canonical straight ray strictly descends. This corrects the preliminary
  approximate coefficient `−4|w|²(|w|²+|z|²)`. The transition `w=O(εz)` still needs a weighted
  blow-up, and an exact affine-unitary local slice must uniformly couple the residual modes to
  the negative second-order directions; both are needed for a neighbourhood theorem.
- **EXACT `p=3` LOCAL SLICE (L69):** every nearby matrix is affine-unitarily reducible to
  `C₃+E(z,w,s,v)` with seven real coordinates. The slice is exactly orbit-orthogonal and L63
  becomes `−8s²−21|v|²/4`, leaving precisely the four residual `(z,w)` coordinates. Thus the
  full local question is now one finite weighted coupling problem; there is no additional hidden
  flat direction.

## Current next actions (Epoch 6, refreshed 2026-07-22)
1. **Finish normal-face uniformity and merge L93 back into L86.** L94--L98 now close every
   bounded weighted chart at a nonnormal terminal block.  Treat `a->0`, where the block
   approaches L88's exact normal direct-sum manifold; lift that transverse estimate through
   the metric flag and couple it to the negative Gram terms, common `v` curvature, and L82's
   losing mean gap.  Do not grind larger Schur matrices.  The target is a full repeated-`C3`
   neighbourhood theorem.
2. **Retain the full CP correction.** Hartz--McCarthy cannot scalarize it. Test whether L21's
   trace inequality can instead be derived from the block-Toeplitz positivity of the full
   operator-valued correction moments. Do not retry trace/positive-state scalarizations.
3. **Shifted Möbius H-r fallback**: derive its exact stationarity/rho formula (Kenan-Li quartic
   analog) and prove rho≥0 or K≤2. Definite parity is false.
4. **Bi-conic Schwarzian test**: compute `SG` for the off-slice collapsed map on the critical
   real interval. `SG≥0` would extend L17 immediately; otherwise test the weaker Sturm-potential
   comparison that the proof actually needs.
5. **Odd phase positivity**: ρ = B₁g₁q₁+B₂g₂q₂ ≥ 0 given the L15 stationarity law
   (3-parameter; interlacing τ₂ < α < τ₁; term-1 dominance observed). Try the same
   deformation/kernel machinery.
6. Rigor debts: n=6 structured floor; 2×2 α=0; contact degeneracy; L59 publication-level
   novelty audit.
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
  slice_positive_face_audit.py (L58 exact positive-face/first-deficit audit),
  slice_positive_corner_audit.py, slice_positive_face_certificate.py,
  slice_positive_tail_certificate.py (L59 full regenerating tail certificate),
  general_similarity_sdp.py (general L21 probe),
  general_similarity_scalarization_probe.py (CP-to-HM state-scalarization falsification),
  general_similarity_equality_probe.py (repeated-Crabb general gate),
  general_similarity_tangent_probe.py (L61 finite conformal/tangent-SDP reduction),
  general_similarity_second_order_probe.py (L62 single-Crabb second-order SDP),
  crabb_second_order_symbolic.py (shared exact conformal derivation),
  p3_second_order_identity.py (L63 exact 18-variable regeneration),
  p4_second_order_identity.py (L64 exact 32-variable regeneration),
  general_crabb_second_order_modes.py (L65 arbitrary-size mode regeneration),
  crabb_second_order_equality.py (L66 exact-orbit/kernel quotient audit),
  formal_riemann_series.py + rank_one_stein_series.py (exact higher-order helpers),
  p3_crabb_sixth_order.py + p3_crabb_quartic.py (L67/L68 exact certificates),
  p3_crabb_local_slice.py (L69 orbit-normal slice audit),
  p3_sparse_series.py + p3_crabb_weighted_slice.py + p3_crabb_center_jet.py
  (L70 sparse exact weighted certificate), p3_crabb_disk_center.py
  (L71 disk factorization and canonical singular product), p3_disk_center_tangent.py
  (L72 exact ambient-stationarity audit), p3_disk_morse_bott.py
  (L73 defect Hessian/invariant audit), p3_local_theorem_probe.py
  (non-load-bearing L73 numerical smoke test), repeated_p3_common_maximizer.py
  (L74 exact repeated-block cross quotient), repeated_p3_second_support.py
  (L75 exact effective support and conformal collapse), repeated_p3_stein_sign.py
  (L76 exact full second-metric certificate and sign), repeated_p3_third_sign.py
  (L77 exact third support split and strict cubic certificate),
  repeated_p3_star_second_sign.py (L78--L80 arbitrary-multiplicity star theorems),
  repeated_p3_radial_gap.py (L81 strict radial diagonal--cross coupling),
  repeated_p3_winner_graph_sign.py (L83--L85 complete tied-winner endpoint),
  repeated_p3_multiwinner_gap.py (L86 full winner/loser Gram penalty),
  repeated_p3_flat_copy_matrix.py (L88/L90 flat/common copy-matrix support),
  repeated_p3_flat_three_copy.py (L89 three-copy Schur equality split),
  repeated_p3_third_metric.py (shared L77/L92 third-metric propagation),
  repeated_p3_flat_two_copy_third.py (L92 trace-zero cubic theorem),
  repeated_p3_flat_two_copy_jensen.py (L95--L97 terminal Jensen coercivity),
  repeated_p3_flat_two_copy_weighted.py (L98 weighted terminal theorem),
  repeated_p3_flat_metric_flag.py (L93 arbitrary-copy metric-flag derivative).
Proof artifact: experiments/positive_tail_full_20260722.log (L59 clean 10-box run).
Data: sym3_sweep_s51.jsonl (40 rec), sym4_sweep_s61.jsonl (20 rec, ρ column trustworthy),
general_similarity_equality_s9173401.jsonl plus its `sensitivity` and `ultralocal` companions
(L21 equality gate); general_similarity_tangent_s9173401.jsonl and
general_similarity_tangent_all_s9173401.jsonl (L61 cross-solver tangent data);
general_similarity_second_order_all_s9173401.jsonl and its `halfstep` companion (L62 data);
general_crabb_second_order_modes_s70221.jsonl (L65 sizes `3..30` mode/rank audit);
crabb_second_order_equality_s70221.jsonl (L66 sizes `3..8` quotient audit).
Ledgers: LEMMA_LEDGER.md, APPROACH_LEDGER.md (pitfalls P1–P8 — READ BEFORE ANY SEARCH),
LITERATURE_LEDGER.md, COUNTEREXAMPLE_SEARCH.md. Audit: chatgpt/FABLE_RESEARCH_AUDIT.md
(reconciled 2026-07-20). Restart: checkpoints/RESTART_PACKET.md (paste-ready instruction).
