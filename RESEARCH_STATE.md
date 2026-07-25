# RESEARCH_STATE.md — Crouzeix Conjecture Campaign

**Last updated:** 2026-07-24 (Epoch 6 — repeated circular/elliptic merger)

## NEWEST (2026-07-24): L248 turns the open flux into one volume coefficient
- For any first active final-defect Schur face, its trace is the
  corresponding coefficient of
  `log det(H+I−F)−log det(FHF)`.  This identity is exact before
  coefficient extraction and absorbs the large corner/Schur-square
  cancellation into one scalar determinant.
- After metric normalization, with `H=I−CC*`, the determinant is
  `2^p det(I−C*(F+(I−F)/2)C)/det(I−C*FC)`.  Equivalently it is the
  determinant of the main output row after whitening the final row by
  `(I−C*FC)^−1`.
- Thus L247's open operator flux is now the scalar target
  `[c^(2k)] log V=4||B_k||_F²`.  In generic delays 1--6 the separate
  corner and Schur-square traces grow to thousands of times
  `||B_k||²` while the whitened difference remains exactly four, so
  those terms must not be estimated separately.  The value four is
  still open; insert L243--L245 directly into the whitened channel.
  `proof/associated_defect_volume_trace.md`;
  `experiments/associated_defect_volume_trace.py`.

## NEWEST (2026-07-24): L247 isolates the only open delayed trace flux
- In the left-dual formulation from L246, deleting only the active
  boundary-metric coefficient changes the grade-`k` Schur face by
  `(I−F)(−X_k+SX_kS*)(I−F)`.
- Its trace is `tr((F−E)X_k)=−2||B_k||_F²`, exactly, from L223/L225's
  two metric endpoint faces.
- Therefore the candidate total trace `+2||B_k||_F²` is equivalent to
  one metric-edge-deleted reflection flux `+4||B_k||_F²`.  L245's
  remote Green column begins with amplitude two, so the normalization
  is explained; the coisometric boundary-energy identity itself is
  still open and must cancel the later `B#(x)` terms explicitly.
  `proof/repeated_crabb_dual_metric_trace_split.md`;
  `experiments/repeated_crabb_dual_metric_trace_split.py`.

## NEWEST (2026-07-24): L246 makes the first defect face endpoint-dual
- For any analytic perturbation of a partial isometry, the
  right-initial and left-final defect Schur quotients obey an exact
  graph intertwining.  If their first nonzero degree is `d`, then
  `K_left[d]=S K_right[d] S*`; analytic metric congruence does not
  alter that first face.
- Thus L228's matrix face, trace, rank, and inertia may be calculated
  at L245's remote clean-chain left endpoint and transported back.
  Its proposed right face becomes exactly `EF_k+F_kE`, the cleaner
  Hardy cell `(0,k)` at L244's zero-th coisometric boundary row.
  This does not yet calculate the delayed coefficient, but it removes
  the need to force the calculation through the less natural right
  pivot.
  `proof/associated_defect_schur_duality.md`;
  `experiments/associated_defect_schur_duality.py`.

## ACTIVE CANDIDATE (2026-07-24): trace-only endpoint filtration
- If the post-Schur trace at degree `d` contains only closed transfer
  words of total endpoint weight at most `d`, complete delay and
  independent defect-frame covariance force the grade-`k` face to be
  `alpha_k ||B_k||_F²`; L241 gives `alpha_k=2` in every grade.
- This would prove the pointwise all-grade separator sign needed by
  L222 without proving L228's full matrix anticommutator.  The complete
  pipeline passes on generic tails through grade eight.
- **Not proved:** `B#(x)` in L243 has no visible `c` valuation, so later
  transfer coefficients occur in the pre-Schur zero/one sectors.  The
  missing theorem is that L244's coisometric Gram quotient cancels
  those terms and restores symmetric endpoint valuation.  Do not bank
  the trace law until that identity is explicit.
  `proof/repeated_crabb_delayed_trace_law.md`;
  `experiments/repeated_crabb_delayed_trace_law.py`.

## NEWEST (2026-07-24): L245 closes every exiting-chain Green column
- For a clean delay `r`, both exterior chain factors in the full block
  inverse have exact continuant entries:
  the last row is `Delta_j/Delta_r`, while the incoming column is
  `2c^r/Delta_r` at level zero and
  `c^(r−j)Delta_j/Delta_r` thereafter.
- On the Joukowski contour these become explicit powers of `x` or
  `rho` times `(1+t^j)/(1+t^r)`.  At distance `d`, the incoming
  entry is exactly
  `rho^d+zeta^d lambda_r(t)(1−t^d)`: the half-line column plus the
  same single reflected weight as L239/L243.  Its endpoint pairing is
  exactly `gamma_r`.
  Hence the “exiting-chain columns” left open by L242 are no longer
  unknown state quantities; together with L243 they give the entire
  full resolvent at every delay.
- L228 is now a finite Laurent-residue/metric/Stein calculation with
  explicit scalar chain features, the two inverse-kernel sectors from
  L243, and L244's coisometric Gram quotient.
  `proof/repeated_crabb_chain_green_columns.md`;
  `experiments/repeated_crabb_chain_green_columns.py`.

## NEWEST (2026-07-24): L244 upgrades the half-line baseline to a coisometry
- L240's defect column is exactly normalized in the inverse metric:
  `d_inf* D_R^−1 d_inf=I`.  This is Jacobi's Lambert-series identity
  `theta_3(q)^2=1+4 sum_(j>=1) q^j/(1+q^(2j))`.
- Therefore `Atilde_inf=D_R^(1/2)A_inf D_R^(−1/2)` is a partial
  isometry.  The half-line ellipse pencil is explicitly surjective
  for `|c|<1`, and the odd Riemann-map factor does not change its
  kernel/range; hence `Atilde_inf` is a coisometry and
  `A_inf D_R^−1 A_inf*=D_R^−1`.
- This identifies the structural cancellation needed after L243:
  a reflected state column's quadratic Gram belongs to the updated
  coisometric defect and is removed by the right-defect Schur square.
  The remaining L228 content is the linear cross-cell response in
  the finite zero/one inverse-kernel sectors.
  `proof/repeated_crabb_half_line_coisometry.md`;
  `experiments/repeated_crabb_half_line_coisometry.py`.

## NEWEST (2026-07-24): L243 identifies the exact delayed model-kernel lift
- After absorbing the right endpoint, the retained bulk inverse is
  `R_R=(G-cT*E)^−1` and its left compression is exactly
  `x(1−t)^−1(I+B(rho)B#(x))`.
- For the delayed transfer `C(u)=u^r B(u)`, the remaining endpoint
  denominator is
  `(I−C(rho)C#(x))/(1+t^r)`.  Consequently the full retained
  resolvent correction is
  `zeta t^r R_R W K_C(rho,x)^−1 W*R_R`: both exterior state columns
  survive, while the middle factor is precisely the inverse delayed
  model kernel.
- The pointwise identity
  `K_C=sum_(j<r)t^j I+t^r K_B` makes L239's physical continued
  fraction identical to L221's abstract model-space delay splitting.
  Expanding the inverse proves that only its zeroth and first
  nonconstant powers can reach the direct-map target window; powers
  two and higher start after `c^(2r+2)`.
- L228 still needs the finite two-sector contour/metric/Stein
  assembly, but there is no longer an infinite reflection tail.
  `proof/repeated_crabb_model_kernel_resolvent.md`;
  `experiments/repeated_crabb_model_kernel_resolvent.py`.

## NEWEST (2026-07-24): L242 fixes the retained reflected jet in all delays
- Relative to L240's unreflected half-line, the retained direct map
  after a finite delay `r>=2` has the exact form
  `(-1)^r(c^(2r)D_0+c^(2r+1)D_1+c^(2r+2)D_2)+...`,
  with fixed ordered tail polynomials of `2,6,23` words.
- A Laurent-cost bound proves that two terminal reflections start at
  `c^(4r)`, after the requested order for every `r>=2`.  The remaining
  one-reflection enumeration has at most five background letters.
  L125 and L237 make every delay-dependent scalar coefficient cancel
  symbolically for `r>=3`; `r=2` gives the same result separately.
- This is the complete **retained direct block**, not the complete
  slack.  The exiting-chain direct columns, L219 metric terms, and
  right-defect Stein Schur square remain.  They must be varied against
  L240's exact rank-`m` baseline before L228 can be promoted.  L241
  separately fixes the eventual scalar multiplier on monomial cells.
  `proof/repeated_crabb_reflected_direct_jet.md`;
  `experiments/repeated_crabb_reflected_direct_jet.py`.

## NEWEST (2026-07-24): L241 fixes the delayed multiplier in every grade
- On the length-`k` monomial channel, L219's physical boundary metric
  and L117's exact axis metric agree below `q^k`; their first
  difference is exactly `diag(1,6,...,6,12)`.  In balanced
  coordinates this is `3I−2E`.
- Because the exact axis metric has zero Stein Schur residual, the
  boundary metric's first slack difference is
  `(3I−2E)−S*(3I−2E)S=E+2E_1`; Schur-compressing away from `E` leaves
  `2E_1`.
- On the monomial channel `F_(k−1)=E_1`, so this is exactly
  `E_1F_(k−1)+F_(k−1)E_1`.  The scalar multiplier in L228 is
  therefore one for every grade, proved from the all-size axis rather
  than extrapolated from finite jets.
- The remaining L228 issue is purely the general ordered-cell lift:
  prove that a nonunitary/noncommuting first active transfer cell
  creates no additional terms beyond the Hermitian lift already
  identified in L236.
  `proof/repeated_crabb_monomial_slack_face.md`;
  `experiments/repeated_crabb_monomial_slack_face.py`.

## NEWEST (2026-07-24): L240 closes the half-line bulk slack exactly
- On the multiplicity-`m` backward-shift half-line, the right-balanced
  ellipse pencil has boundary metric
  `D_R=diag(1,(1+q)^−1,(1+q²)^−1,...)` and exact Stein factor
  `D_R−A_inf*D_R A_inf=d_inf d_inf*`.
- The defect column is the normalized Jacobi `nd` Fourier column:
  its level-zero block is `theta_3(c²)^−1 I`, its level `2j` block is
  `2(−c)^j/[theta_3(c²)(1+c^(4j))] I`, and its odd blocks vanish.
  Hence the right-defect Schur residual of the full half-line slack
  is identically zero.
- The proof is all-order coefficientwise finite-path stabilization of
  L117's exact arbitrary-size elliptic Crabb-axis metric; it makes no
  unproved boundary infinite-operator calculus claim.  An independent
  coisometric word audit regenerates the full theta/ODE identity
  exactly through degree eight.
- This removes the hidden-flux ambiguity in retained-block
  calculations.  The only remaining L228 content is the first active
  coefficient of the exact L239-reflected transform of this zero
  Schur residual, including its state-lift columns.  Do not silently
  truncate to the term linear in `lambda_r`; the exact small endpoint
  inverse must justify which reflection powers contribute.
  `proof/repeated_crabb_half_line_stein_factor.md`;
  `experiments/repeated_crabb_half_line_stein_factor.py`.

## NEWEST (2026-07-24): L239 makes the delay weight exact on the contour
- In L238's Joukowski coordinate `t=c/zeta²`, the terminal
  continuants collapse to `Delta_j=zeta^j(1+t^j)`.  Hence the
  normalized reflection is exactly
  `lambda_r=(gamma_r−c/zeta)/(zeta(1−t))=t^r/(1+t^r)`, not merely
  `t^r+...`.  This is the same rational weight family that occurs in
  L219's boundary metric.
- The formerly unsimplified lower-right bulk compression is the
  characteristic kernel
  `x(1−t)^−1(I−t B(rho)B#(x))`.  Schur-eliminating the right endpoint
  from L238's `2m` scattering matrix therefore leaves exactly
  `I−lambda_r(I+B(rho)B#(x))`, with multiplication order preserved.
- This removes all higher delay corrections from the live gate.  It
  does not yet prove L228: the remaining calculation must lift the
  contour resolvent back to the state space, combine it with L236's
  Hardy-frame metric, and include the Stein right-defect Schur square.
  `proof/repeated_crabb_scattering_schur_collapse.md`;
  `experiments/repeated_crabb_scattering_schur_collapse.py`.

## NEWEST (2026-07-24): L238 compresses the two endpoints before expansion
- On the Joukowski contour `z=zeta+c/zeta`, put `rho=c/zeta` and
  `delta_r=gamma_r−rho`.  After L237's complete-delay elimination, the
  retained denominator factors exactly as
  `(zeta I−T)(I−rho T*)−cT*E−delta_rF`.
- Woodbury reduces the full retained resolvent to one `2m x 2m`
  endpoint scattering matrix.  Its off-diagonal blocks are exactly
  `delta_r V*(zeta I−T)^−1W` and `(z/zeta)B(rho)`.  Thus the two
  orientations of the tail transfer appear without reordering, and
  arbitrary delay remains confined to the scalar `delta_r`.
- Six deterministic noncommuting audits through delay six reproduce
  the original full retained resolvent below `1.7e−15`.
- This is an exact bounded-size normal form, not L228 itself.  The live
  gate is now to evaluate the **full** L125 theta/ODE contour map and
  L219 boundary metric/right-defect Schur square on this scattering
  matrix.  That calculation must establish the scalar multiplier and
  the final anticommutator; no further isolated free-word grade should
  be computed.
  `proof/repeated_crabb_two_defect_scattering.md`;
  `experiments/repeated_crabb_two_defect_scattering.py`.

## NEWEST (2026-07-24): L237 resums every complete delay at once
- If `B_1=...=B_r=0`, the removed left-wandering chain is the clean
  scalar Jacobi corner `L_r=N_r+c(I+F_0)N_r*`.  Eliminating the whole
  chain gives the exact retained resolvent
  `(zI-Xi_(r,−)-gamma_r F_r)^−1`; no noncommutative word list grows
  with `r`.
- The scalar self-energy is
  `gamma_1=2c/z` and `gamma_r=c Delta_(r−1)/Delta_r`, where
  `Delta_2=z²−2c` and
  `Delta_j=zDelta_(j−1)−cDelta_(j−2)`.
- Compare this with the unreflected half-line Catalan self-energy
  `eta=(z−sqrt(z²−4c))/2`.  An exact continued-fraction induction
  proves
  `gamma_r−eta=c^r/z^(2r−1)+O(c^(r+1))`.
  The first reflected endpoint coefficient is therefore exactly one
  in every grade.  This is the first genuinely arbitrary-grade
  recurrence on the A178 line.
- L221's retained boundary metric has the matching factor `q^r` on
  the promoted left orbits.  Together with L236's two-frame/Hankel
  normal form, the live gate is now a single
  two-boundary coefficient calculation: with the **full** L125
  theta/ODE scalar map, pair L237's left reflection with the right
  defect and include the Schur square to derive L228's
  `E_1F_r+F_rE_1`.
  `proof/repeated_crabb_multidelay_terminal_resolvent.md`;
  `experiments/repeated_crabb_multidelay_terminal_resolvent.py`.

## NEWEST (2026-07-24): L236 puts the target in two Hardy frames
- The right and left defect-orbit analysis maps `O_R,O_L` are exact
  isometries into vector-valued Hardy coefficient space and intertwine
  `S,S*` with the backward shift.
- Their cross Gram is the transfer Hankel matrix:
  `(O_R O_L*)_(n,j)=B_(n+j)*`.  Thus a complete delay is literally a
  zero prefix of Hankel anti-diagonals; the first active matrix block
  is `B_k*`.
- L219's boundary metric is exactly
  `O_R*D_R(q)O_R+O_L*D_L(q)O_L`, with universal diagonal scalar
  weights `D_R=diag(1,(1+q)^−1,(1+q²)^−1,...)` and
  `D_L=diag(0,q,q²,...)`.
- L228's proposed `E_1F_(k−1)+F_(k−1)E_1` is precisely the Hermitian
  state-space lift of the active Hankel cell `(1,k−1)`.  The remaining
  issue is no longer matrix ordering or target identification: it is
  proving that the full grouped physical ellipse/Schur operation is
  causal for these anti-diagonals and has scalar multiplier one.
  `proof/repeated_crabb_hardy_two_frame.md`;
  `experiments/repeated_crabb_hardy_two_frame.py`.

## NEWEST (2026-07-24): L235 isolates one delay as one terminal insertion
- After removing the first left wandering layer, the balanced ellipse
  pencil has the exact arrowhead form
  `[[0,2cWtilde*],[Wtilde,S_1+cS_1*(I+Etilde)]]`.
  Its retained resolvent is
  `(zI-Xi_--2cz^−1Ftilde)^−1`.  Thus every terminal round trip is
  generated by the single rank-`m` insertion `2cz^−1Ftilde`.
- The boundary-layer metric has a parallel exact `2 x 2` block
  formula in `B_d` and the retained right orbits.  These two formulas
  implement L228's requested zero/one/two-crossing grouping before
  expanding any scalar Riemann coefficient.
- The independently balanced tail pencil differs from `Xi_-` by
  `cFtilde S_1*(I+Etilde)`.  Accordingly, the tempting whole-series identity
  `J*K_SJ=c²K_(S_1)` is false: eight complete-delay tests reproduce
  its leading face within `1.5e−13` but miss at the next even
  coefficient by `0.866`--`2.016`.
- L125's Newton bottom edge is not sufficient by itself either.
  Dropping the true scalar term `a_0(c)=1+2c²+...` leaves an exact
  spurious `4S*S` already in the grade-one residual; the omitted
  `2c²S` contributes precisely `−4S*S`.  Use the full theta/ODE
  scalar recurrence inside the crossing sectors.
- The correct remaining target is only the associated-graded
  first-active cancellation after imposing the tail delay ideal.
  Expand L125's scalar coefficients inside the grouped resolvent
  sectors, not in the original free-word basis.
  `proof/repeated_crabb_one_delay_terminal_block.md`;
  `experiments/repeated_crabb_one_delay_block.py`.

## NEWEST (2026-07-24): L234 closes the prepared sextic third flag
- A universal 53-term perpendicular sixth column, of coefficient
  l1-norm `181`, has an exact 139-term Stein witness and 28-term
  right-ideal endpoint certificate.  On
  `ker B_1* intersect ker B_2*`, its direct sixth face is exactly
  `12B_3B_3*`.
- The fifth off-diagonal block cannot be ignored.  L233's exact
  factorization reduces its surviving action to
  `-14B_1 R B_3*`, `R=B_1*B_1`; L232's quartic range has weight
  `8(4I+7R)` in the same channel.
- Completing the square leaves
  `12I-(49/2)R(4I+7R)^-1R`.  Since `0<=R<=I`, this is uniformly at
  least `215I/22`.  The effective sixth face is therefore positive
  through every partial rank change and passes to the next flag
  exactly when `B_3` also vanishes.
- A first complete-delay-only certificate was not accepted: after
  passing 1,584 moderate tests it failed on valid high-amplitude
  rank-two flags with minimum `-0.0352898`.  The final right-ideal
  certificate fixes precisely that stale-risk gap.
- This closes the explicit preparation through transfer grade three,
  not the arbitrary-grade induction.  The next attack is to combine
  the odd-cancellation/even-Gram/Schur-budget pattern with L228's
  delayed anticommutator recursion while controlling coefficient
  growth and convergence in L194's chart.
- **Recursion status (checkpoint): no arbitrary-grade formula is
  currently known.**  The preparation columns were solved separately:
  `C_3` has one contraction term (bound `3`), `C_4` has ten terms
  (coefficient l1-bound `27/2`), `C_5` has twenty terms (bound `48`),
  and `C_6` has 53 terms (bound `181`).  These results rigorously close
  the first three transfer flags, but the observed growth does not
  constitute an induction or a convergent all-series construction.
  Closing another isolated grade will count as evidence only if it
  exposes a uniform right-ideal recurrence and a summable/analytic
  coefficient bound.  A171/L228's one-delay anticommutator recursion
  remains open and is **still the intended structural route** to
  arbitrary grade; it is neither parked nor superseded by L230--L234.
  `proof/repeated_crabb_canonical_sextic_preimage.md`;
  `experiments/repeated_crabb_canonical_sextic_preimage.py`.

## NEWEST (2026-07-24): L233 clears the complete quintic second flag
- After inserting L230's cubic and L232's quartic columns into the
  exact canonical factor, the fifth Stein forcing has an explicit
  bounded polynomial preparation `C_5=Q Ccal_5(S,S*)V`.
  `Ccal_5` has twenty half-integral terms and coefficient l1-norm
  `48`, so no rank-dependent inverse or selection is present.
- A 55-term Hermitian witness proves the corrected forcing is an
  exact Stein coboundary.  The resulting upper endpoint lift then
  factors much more compactly as
  `H_5=R_1ESF+R_2ES²F+h.c.`.
- Since `ES^jWx=VB_j*x`, the whole fifth endpoint compresses to zero
  on `ker B_1* intersect ker B_2*`, including partial rank changes.
  Exact rational residuals are zero; 27 independent matrix
  reconstructions agree below `1.4e-12`.
- This removes the necessary fifth odd obstruction but does not yet
  close sextic positivity.  On partial flags the off-diagonal fifth
  block incurs a sextic Schur cost against L232's positive quartic
  range.  On complete double delays, the sextic coefficient must be
  recomputed because `C_4` and `C_5` are endpoint-null but nonzero
  gauges.
  `proof/repeated_crabb_canonical_quintic_preimage.md`;
  `experiments/repeated_crabb_canonical_quintic_preimage.py`.

## NEWEST (2026-07-24): L232 makes the full quartic endpoint positive
- L231's positive trace can be redistributed by one universal
  ten-term perpendicular polynomial column.  Its norm is bounded by
  `27/2`, independently of transfer ranks and state dimension.
- A 29-term Hermitian Stein witness proves, with zero exact rational
  word residuals, that the final physical quartic upper gap is
  `12B_2B_2*+32B_1B_1*+56B_1(B_1*B_1)B_1* >= 0`.
- The correction has zero lower endpoint response.  It uses no
  inverse, singular vector, flag projection, or pseudoinverse, so
  the partial rank-changing quartic obstruction is completely
  removed.
- On complete delays the polynomial column is bounded and
  endpoint-null, but generally nonzero.  This is a legitimate gauge
  at quartic order and a warning that quintic/sextic coefficients
  must be recomputed after both preparations.
  `proof/repeated_crabb_canonical_quartic_preimage.md`;
  `experiments/repeated_crabb_canonical_quartic_preimage.py`.

## NEWEST (2026-07-24): L231 proves the prepared quartic trace is positive
- Complete L230 by perturbing the exact canonical Stein-slack factor
  by `c³C_3` and applying the variable Stein inverse.  This preserves
  exact contraction and makes the full cubic endpoint zero.
- Exact word reduction gives the physical quartic upper-gap trace
  `12||B_2||_F²+32||B_1||_F²+56tr((B_1*B_1)²)`.
  The only residual is
  `8tr(−I+2S*S−(S*)²S²)=0`, from defect ranks
  `n,n−m,n−2m`.
- On `B_1=0`, the complete quartic matrix is exactly
  `12B_2B_2*`; the cubic correction vanishes there.  Thus every
  L222 reducing separator on the first partial flag sees positive
  trace unless it also belongs to the next delay.
- Positive trace is not positive semidefiniteness.  The prepared
  quartic compression is slightly indefinite in some multiplicity
  five/six rank chains, so a bounded polynomial fourth-column
  redistribution is still required.
  `proof/repeated_crabb_canonical_quartic_trace.md`;
  `experiments/repeated_crabb_canonical_quartic_trace.py`.

## NEWEST (2026-07-24): L230 cancels the cubic flag polynomially
- L229's pointwise cubic range statement now has an explicit
  all-size selection:
  `C_3 = 3(S*)² V(B_1*B_1)`.
- If `F_3=VC_3*+C_3V*`, the canonical ten-word forcing satisfies
  `L_3+F_3=Z−S*ZS` for a four-word `Z`, and `FZF=0`.  Stein
  telescoping therefore proves
  `[c³]U_can=4M_S(C_3)` with the physical sign and factor intact.
- The column is perpendicular to `V` and obeys
  `||C_3||<=3||B_1||²`.  It is a bounded basis-free real-analytic
  selection through every transfer-rank jump, not merely a
  pointwise Fredholm preimage.
- L194 inserts the corresponding metric row while preserving exact
  lower tightness and nonnegative Stein slack, so A172's full cubic
  endpoint is cancelled.
- This is not yet the metric sandwich.  The L194 completion changes
  quartic and later coefficients; recompute the prepared quartic
  flag before combining it with L228's grade-two even face.
  `proof/repeated_crabb_canonical_cubic_preimage.md`;
  `experiments/repeated_crabb_canonical_cubic_preimage.py`.

## NEWEST (2026-07-24): L229 removes the cubic pointwise range obstruction
- The bad cubic face of the canonical repair has an exact ten-word
  state forcing `L_3`; its physical endpoint is
  `-4 W*G_S(L_3)W`.
- The ten words cancel in four cyclic trace classes.  Consequently
  `tr(HL_3)=0` for every self-adjoint colligation commutant `H`.
  L206's cokernel theorem and Stein adjointness therefore put the
  complete cubic endpoint in `ran M_T` at every fixed equality
  colligation, including reducible ones.
- Scaled irreducible Schur chains give minimum-norm preimages
  `||C_3||=O(lambda^2)` while the endpoint map loses one factor
  `lambda` and the target is `O(lambda^3)`.  The ratios stabilize,
  strongly indicating removable analytic divisibility.
- Pointwise cancellation is now proved; bounded analytic selection
  is not.  Derive a polynomial cubic column or prove the
  `O(lambda^2)` divisibility on every analytic path before moving to
  the delayed even face.
  `proof/repeated_crabb_canonical_cubic_range.md`;
  `experiments/repeated_crabb_canonical_cubic_selection.py`.

## NEWEST (2026-07-24): the canonical repair fails on partial flags
- L227's exact Stein repair is not by itself a condition-number-four
  metric.  On rank-one `B_1` faces its upper gap begins with the
  favorable `12c^2 B_1B_1*`, but the compression of its cubic
  coefficient to `ker B_1*` can be indefinite.
- A multiplicity-four example has cubic kernel eigenvalues
  `(-0.0119641, 0.0000766, 0.0090452)`.  Since the active/kernel
  Schur cross-square starts at order four, the negative cubic
  eigenvalue forces the repaired physical metric above `4I` for all
  sufficiently small positive `c`.
- Scaled noncommuting Schur chains retain the obstruction down to
  parameter scale `0.1`, so it occurs arbitrarily near the repeated
  monomial apex.  This does not contradict L228: the cubic vanishes
  on the exact `B_1=0` stratum.
- The corrected next target is to cancel the odd mixed-flag endpoint
  by a bounded analytic L204 free-row/defect-frame correction before
  using L228's delayed even Gram.  Complete-delay covariance remains
  useful but cannot by itself finish the partial-flag sandwich.
  `proof/repeated_crabb_canonical_repair_flag_obstruction.md`;
  `experiments/repeated_crabb_canonical_repair_flag_obstruction.py`.

## NEWEST (2026-07-24): L228 exposes the delayed slack as one anticommutator
- L227's compact grade-one residual simplifies once more to
  `K_2=E_1F+FE_1`, where `E_1=S*ES`.
- Therefore L225's unresolved all-delay covariance is equivalent to
  the concrete formula
  `[c^(2k)]K_S=E_1F_(k-1)+F_(k-1)E_1` under
  `B_1=...=B_(k-1)=0`.  This is still a candidate arbitrary-grade
  identity, not a proved theorem.
- Exact noncommutative rational word reduction verifies the formula
  and all earlier vanishings through grade five.  The degree-ten face
  has 3277 words before imposing the delay ideal, strengthening the
  earlier floating audit without turning it into an induction.
- The candidate's consequences are now proved in every grade:
  its Stein endpoint is exactly `2B_kB_k*` and its trace is
  `2||B_k||_F^2`.  Hence proving the coefficient identity would give
  the canonical repair the faces `+B_k*B_k` and `-12B_kB_k*`.
- The live proof debt is a one-delay associated-graded recursion.
  Group terminal-crossing words before expanding their scalar ellipse
  coefficients; another fixed-grade jet would not close the theorem.
  `proof/repeated_crabb_delayed_slack_anticommutator.md`;
  `experiments/repeated_crabb_delayed_slack_anticommutator.py`.

## NEWEST (2026-07-24): L227 gives an exact canonical Stein repair
- L225's grade-one boundary-slack Schur residual collapses from a
  long ellipse-jet expression to
  `K_2=2F−(S*)²S²F−F(S*)²S²`.
- Its complete upper Stein response is the matrix identity
  `W*G_S(K_2)W=2B_1B_1*`, strengthening L225's scalar trace law.
- More generally, if `H` is the full boundary slack and
  `G=V*HV`, then
  `X=−G_T(H−HVG^−1V*H)` is analytic and makes the repaired slack
  exactly `HVG^−1V*H>=0`.  Thus a bounded exact contraction repair
  exists canonically; endpoint condition geometry remains to prove.
- At grade one this repair preserves the lower face
  `+B_1*B_1` and changes the upper face from `−4B_1B_1*` to
  `−12B_1B_1*`.  If L225's delayed covariance is proved, L216
  transports this full matrix response to every completely delayed
  grade.
  `proof/repeated_crabb_boundary_slack_repair.md`;
  `experiments/repeated_crabb_boundary_slack_repair.py`.

## NEWEST (2026-07-24): L226 gives the matching ordered lower flag
- On the surviving right-copy flag
  `R_(k-1)=intersection_(j<k) ker B_j`, the boundary metric's first
  lower physical Schur face is exactly
  `+c^(2k)(B_k K)*(B_k K)`.
- The proof is the right-oriented counterpart of L224.  The complete
  lower endpoint factors into direct right transfer Grams and
  interactions containing two transfer legs; all earlier active
  directions can be Schur-eliminated without a pseudoinverse.
- L224 and L226 now expose both endpoint budgets through arbitrary
  noncommuting rank changes.  L225's conditional repair accounting
  is therefore flag-compatible on both sides.
- This still does not make the boundary metric contractive.  The
  live proof gate remains the associated-graded delay covariance of
  its Stein-slack Schur residual, followed by bounded repair
  selection.
  `proof/repeated_crabb_lower_metric_flag.md`;
  `experiments/repeated_crabb_lower_metric_flag.py`.

## NEWEST (2026-07-24): L225 reduces the delayed trace law to one slack covariance
- L223's boundary metric has not only the upper face
  `-4c^(2k) B_k B_k*`, but the exact lower face
  `+c^(2k) B_k* B_k` on a fully delayed grade.
- For the boundary metric's Stein slack, take the Schur complement
  away from the right defect.  At grade one its second coefficient
  has exact trace `2 ||B_1||_F^2`; the proof is an order-preserving
  partial-isometry trace reduction.
- If the metric is repaired to rank `m` and the lower face is
  re-tightened, the dual orbit identity forces the repaired upper
  trace to be `-16 ||B_k||_F^2` provided the first slack Schur trace
  remains `2 ||B_k||_F^2`.
- The remaining identity is now precise: on a complete delay, the
  first slack Schur coefficient should be the embedded grade-one
  coefficient of the deflated colligation.  This stronger matrix
  covariance passes unstructured grades one through five, but is
  still **numerical**, not proved.  Prove its one-delay associated-
  graded recursion; then L222's pointwise trace sign closes.
  `proof/repeated_crabb_boundary_slack_deflation.md`;
  `experiments/repeated_crabb_boundary_slack_deflation.py`.

## NEWEST (2026-07-24): L224 exposes the exact partial-flag upper budget
- Let `K_(k-1)=intersection_(j<k) ker B_j*`.  After eliminating the
  state complement and every earlier active copy direction, L219's
  physical boundary metric has exact first face
  `-4c^(2k)(U*B_k)(U*B_k)*` on an isometry `U` onto `K_(k-1)`.
- The proof uses the exact full endpoint expansion
  `H=sum a_jB_jB_j*−sum a_ja_l B_jC_j*A(q)^−1C_lB_l*`.
  The direct grade-`k` Gram survives; all two-leg interactions and the
  singular earlier-block Schur square start at least one order later.
- This extends L223 through arbitrary noncommuting rank changes and
  removes the need to guess a partial-flag whitening for the endpoint
  budget.  It remains an endpoint theorem, not a contraction theorem:
  the boundary metric's Stein slack can still be indefinite.
- The live quantitative gate is to construct a structured
  Stein-positive repair whose cost on every L224 flag is strictly
  below the available factor four.
  `proof/repeated_crabb_boundary_metric_flag.md`;
  `experiments/repeated_crabb_boundary_metric_flag.py`.

## NEWEST (2026-07-24): L223 gives every delayed flag an exact upper budget
- For L219's explicit boundary-layer metric, if
  `B_1=...=B_(k-1)=0`, the first upper physical Schur face is exactly
  `-4c^(2k)B_kB_k*`.
- The proof is elementary and all-grade.  In the positive upper gap,
  only the `k`-th right-defect orbit has a leading compression to
  `W`; every left orbit annihilates `W`, and the cross Schur square
  starts at order `c^(4k)`.
- Thus the correct left orientation, reflected weight, rank-changing
  kernel, and terminal coercivity are already present before any
  Stein factorization.  At the repeated apex the face is
  `-4c^(2L)I`.
- This does not solve contraction: the same boundary metric can have
  indefinite elliptic Stein slack.  The remaining task is now
  quantitative—repair the slack while spending strictly less than
  L223's explicit upper budget (or construct a stronger endpoint),
  uniformly through the L222 rank flags.
  `proof/repeated_crabb_boundary_metric_face.md`;
  `experiments/repeated_crabb_boundary_metric_face.py`.

## NEWEST (2026-07-24): L222 reduces partial flags to delayed trace signs
- For any left-copy flag `U`, compress L204's homogeneous endpoint
  response to `C -> U* M_T(C) U`.  Finite-dimensional semidefinite
  separation gives an exact alternative: a face `E` can be made
  strictly negative iff `tr(YE)<0` for every nonzero positive
  annihilator `Y`.
- L206 identifies those annihilators without a rank assumption:
  the observability Gramian of `U Y U*` must commute with the balanced
  partial isometry.  Thus every genuine obstruction is a reducing
  colligation summand, not a new noncommutative partial-flag term.
- If `U*B_1=...=U*B_(k-1)=0`, each positive reducing separator lives
  on a completely delayed summand.  Therefore the all-grade repeated
  elliptic flag's **pointwise range obstruction** needs only the trace
  sign of the fully delayed effective face; exact identification of
  every partial face with `-16B_kB_k*` is stronger than necessary.
- The remaining physical target is now
  `tr E_(2k,eff)=-16||B_k||_F^2` on fully delayed summands.
  L203 proves it at grade one and L214--L215 prove stronger matrix
  identities at grades two and three.  The all-grade trace law and a
  bounded pathwise/Schur selection through commutant rank changes are
  still open.
  `proof/repeated_crabb_flagged_endpoint_alternative.md`;
  `experiments/repeated_crabb_flagged_endpoint_alternative.py`.

## NEWEST (2026-07-24): L221 makes delay an exact Schur-model tail
- If `B_1=...=B_r=0`, then `B=z^r Btilde`, and its model kernel is
  exactly the orthogonal sum of `r` monomial layers and the shifted
  deflated kernel `(zbar w)^r K_Btilde`.
- L220's ordered features shift literally:
  `Phi_j=z^jI` below the tail and
  `Phi_(r+j)=z^r Phitilde_j` inside it.
- The matching state identity
  `(S*)^nV=J(S_r*)^nV_r+sum_(j<r)S^jWB_(n+j)`
  exposes every early future-transfer image in one orthogonal block.
  L219's retained boundary metric also has an exact weighted-tail
  formula.
- Combined with L216, every transfer, feature, state-orbit, metric,
  and linear endpoint object now respects delay removal.  The only
  unproved interface is invariant: show that the prepared physical
  two-reflection quotient is `−16` times the Schur-feature Gram.
  `proof/repeated_crabb_delay_model_flag.md`;
  `experiments/repeated_crabb_delay_model_flag.py`.

## NEWEST (2026-07-24): L220 gives the exact Schur-orthogonal model flag
- L218's one-step kernel identity splits the model kernel into one
  rank-`m` feature layer plus the shifted kernel of the next Schur
  iterate.
- Iteration gives
  `K_B(z,w)=sum_j Phi_j(z)*Phi_j(w)` with the noncommutative order
  `Phi_j(w)=w^j A_j(w)...A_0(w)`.  These features span the complete
  `Lm`-dimensional model space.
- At the repeated apex `Phi_j=w^jI`, so ordinary Hardy/Fourier grades
  are the associated graded of an exact triangular whitening at every
  nearby noncommuting equality anchor.
- This makes the next physical statement precise: prove that the
  first two-reflection upper endpoint is `−16` times the Gram of the
  corresponding feature coefficients.  L216 can then deflate the
  first active layer to the proved grade-one response.
- A stronger raw shortcut is false.  Simultaneously inserting every
  gradewise L212 column gives the correct order-two face but generic
  order-three upper residuals of norms `7.37` and `1.64`, while lower
  tightening remains exact.  Later faces must be Schur-orthogonalized,
  not superposed.
  `proof/repeated_crabb_schur_kernel_flag.md`;
  `experiments/repeated_crabb_schur_kernel_flag.py`;
  `experiments/repeated_crabb_raw_face_superposition.py`.

## NEWEST (2026-07-24): L219 isolates the boundary metric from the Stein gate
- The explicit one-image boundary-layer metric has two exact positive
  orbit decompositions.  After returning to physical coordinates they
  prove `I<P_phys(c)<4I` for every finite pure partial isometry and
  every `0<c<1`, without any delayed-transfer hypothesis.
- Thus the candidate's condition-number geometry is already an
  all-grade theorem; it need not be recovered from separate endpoint
  jets.
- The tempting direct shortcut is false: the same metric's elliptic
  Stein slack is indefinite on deterministic general and delayed
  colligations.  Grade-eight testing independently continued to match
  the one-image face within `9e−10`, but remains numerical evidence.
- The remaining all-grade task is now exactly a Stein repair or
  factorization inside the two explicit positive orbit budgets.  This
  is weaker than preserving a rank-`m` Stein defect and supplies a
  second route alongside the ordered zero/one-image proof.
  `proof/repeated_crabb_boundary_metric_sandwich.md`;
  `experiments/repeated_crabb_boundary_metric_sandwich.py`.

## NEWEST (2026-07-24): L218 gives the ordered matrix Schur generator
- Every square rational inner transfer near `z^L U` has a unique
  ordered matrix Schur recursion.  Its one-step kernel identity and
  inverse preserve multiplication order, while determinant winding
  proves termination after exactly `L` steps.
- L201 forces the first parameter to vanish.  The remaining
  `Gamma_1,...,Gamma_(L−1)` and terminal unitary give
  `(2L−1)m²` real coordinates, exactly matching L193's equality
  dimension.
- At the monomial, the grade-`j` tangent is
  `z^j Delta−z^(2L−j)U Delta*U`.  The late adjoint coefficient is
  forced by innerness and supplies the structural source of L217's
  early future reflections.
- Separately, `proof/repeated_crabb_one_image_generator.md` records a
  coherent theta-weighted zero/one-image physical frame.  Its direct
  jets give the desired `−16B_kB_k*` endpoint through grades one to
  six, including noncommuting and rank-changing cases, but it remains
  a candidate rather than a lemma.
- The immediate proof gate is now precise: prove the candidate's
  ordered zero/one-reflection word identity, then show that its first
  two-reflection quotient is the L216-deflated grade-one face.
  L218 alone does not identify the physical Riemann/metric endpoint.
  `proof/repeated_crabb_matrix_schur_chart.md`;
  `experiments/repeated_crabb_matrix_schur_chart.py`.

## NEWEST (2026-07-24): L217 rules out the naïve delayed-frame induction
- The exact repeated length-five elliptic-axis defect frame differs
  from the universal zero-reflection frame already at order four:
  `−2c^4 S W B_5+2c^5 S^3 W B_5+O(c^6)`.
- Thus a future terminal coefficient can change the lower-tight gauge
  before its own transfer grade.  “Keep the zero-reflection frame
  through order `k−1`, then add a `B_k` term” is false.
- The result is exact: the periodized-sech metric, theta/ODE Riemann
  pullback, and rank-one Stein factor are regenerated symbolically
  through degree ten.
- This explains why a blind grade-four continuation was unstable.
  The next induction must group the complete exact-axis
  normalization before applying L216, or construct its matrix
  Schur/Levinson generating recursion.
  `proof/repeated_crabb_axis_gauge_obstruction.md`;
  `experiments/repeated_crabb_axis_gauge_obstruction.py`.

## NEWEST (2026-07-24): L216 proves all-delay endpoint-response covariance
- If `B_1=...=B_r=0`, remove the wandering states
  `W,SW,...,S^(r−1)W`; the retained partial isometry has left defect
  `S^rW`.
- Every Stein forcing supported on the retained space has exactly the
  same retained Stein solution and upper endpoint before and after
  deflation.  After restoring the two physical equality metrics, the
  L204 endpoint maps are still identical.
- L212's grade-`r+1` column is therefore the literal lift of L207's
  grade-one column for the smaller colligation, including its complete
  channel-coboundary endpoint.
- This closes the linear physical-response half of the proposed
  induction in every grade.  It does not transport the nonlinear raw
  Riemann/metric forcing or its endpoint Schur squares; that prepared
  base covariance is now the sole delay-induction gate.
  `proof/repeated_crabb_endpoint_deflation.md`;
  `experiments/repeated_crabb_endpoint_deflation.py`.

## NEWEST (2026-07-24): L215 closes the twice-delayed physical face
- On `B_1=B_2=0`, add `+2QSWB_3` to the perpendicular third-frame
  coefficient and its negative at frame order five.  These motions
  are endpoint-null at their own orders but essential at order six.
- The complete sixth-order upper Schur coefficient becomes
  `12B_3B_3*−28C(B_3*B_3)` plus two explicit future-row endpoint
  terms involving `B_4*B_3` and `B_5*B_3`.
- One polynomial sixth-frame preparation cancels both future rows.
  L212 then leaves the coercive endpoint `−16B_3B_3*`.
- The gauge matches the exact length-three elliptic-axis defect
  frame.  No transfer inverse, kernel projection, or pseudoinverse
  occurs, and singular/rank-zero `B_3` cases are included.
- Grades two and three now exhibit the same prepared base.  The next
  task is to derive the finite future-row/gauge pattern from one
  generating colligation equation, not to hand-expand grade four.
  `proof/repeated_crabb_grade_three_face.md`;
  `experiments/repeated_crabb_grade_three_face.py`.

## NEWEST (2026-07-24): L214 closes the first delayed physical face
- After L213's axis-compatible gauge, the complete fourth-order upper
  Schur coefficient on `B_1=0` reduces exactly to
  `4 C(B_2*B_2)−20 B_2B_2*−2 M(QS*V B_3*B_2)`.
- The polynomial fourth-frame preparation
  `Q{4S^2WB_2+2S*V B_3*B_2}` changes this to the physical base
  `12B_2B_2*−28C(B_2*B_2)`.
- L212 then adds its grade-two channel coboundary and leaves the
  coercive endpoint `−16B_2B_2*`.
- No transfer inverse, kernel projection, or pseudoinverse is used;
  the formulas remain analytic when `B_2` changes rank and vanish at
  the repeated length-two apex.
- The next gate is to identify and organize the higher axis/null
  gauges before the `B_3` sixth face, then seek the all-grade
  generating formula rather than hand-expand every grade.
  `proof/repeated_crabb_grade_two_face.md`;
  `experiments/repeated_crabb_grade_two_face.py`.

## NEWEST (2026-07-24): L213 fixes the delayed elliptic gauge
- L207's endpoint solution is not unique; later coefficients cannot
  be inferred from its minimum-norm representative.
- On the delayed face `B_1=0`, the exact polynomial normalization is
  `X_hat_ax=SWW*S*−S*VV*S` with full second frame column
  `C_hat_2=−V(V*F_hat_2V)/2+2(S*)^4V`.
- An ordered partial-isometry reduction proves the complete second
  Stein equation and both endpoints vanish.  At every repeated
  monomial apex this is exactly the second coefficient of the known
  all-size elliptic-axis metric.
- Separately, `2W` is a universal endpoint-null column, but it equals
  the delayed axis correction only at the length-four apex.  This
  explains why the first grade-two minimum-norm computation produced
  a false base mismatch.
- Next recompute the fourth-order `B_2` endpoint after this
  normalization, then compare it with L212's
  `12L_2−28C(R_2)` base.
  `proof/repeated_crabb_delayed_axis_gauge.md`;
  `experiments/repeated_crabb_endpoint_null_gauge.py`.

## NEWEST (2026-07-24): L212 removes every higher elliptic range obstruction
- Retaining L208's lower-grade contamination gives the explicit
  all-grade column
  `C_hat_k=−(7/2)Q{S^kWB_k+sum_(j<k)(S*)^(k−j)V B_k*B_j}`.
- It maps exactly to
  `28{C(B_k*B_k)−B_kB_k*}` for every grade, with no copy
  projection, pseudoinverse, or rank assumption.
- Summing with the exact Faber weights `|c|^(2k)` gives a bounded
  analytic preimage for the complete right-channel-minus-left-Gram
  target.  The correction vanishes at the repeated Crabb apex.
- Thus the flag/Schur issue is no longer a **range** issue.  The sole
  elliptic endpoint debt is to derive the prepared physical base
  `12L_c−28C(R_c)` (or a one-sided endpoint below it).  L212 would
  then turn it into the coercive `−16L_c`, and L201 supplies terminal
  positivity.
  `proof/repeated_crabb_all_grade_preimage.md`;
  `experiments/repeated_crabb_all_grade_preimage.py`.

## NEWEST (2026-07-24): L211 closes the flagged CP-channel covariance
- The bottom block of the deflated transfer channel is exactly the
  original channel compressed to the active flag:
  `J_U*C_(S_def)(K)J_U=U*C_S(K)U`.
- Together with L210's promoted left Gram, this makes L208's flagged
  correction exactly the bottom compression of a deflated grade-one
  channel coboundary.  The dual and state-space mechanisms now agree
  term for term.
- An important overclaim was excluded: the full unreduced L207 base
  uses `Btilde_1*Btilde_1`, whose complementary top-row contribution
  is generally nonzero.  It cannot simply be discarded.
- The remaining physical gate is therefore exact: prove that prior
  metric/least-squares Schur elimination replaces the full right Gram
  by its active bottom compression.  No transfer or multiplication
  order remains ambiguous.
  `proof/repeated_crabb_transfer_channel_covariance.md`.

## NEWEST (2026-07-24): L210 closes the transfer-level weighted pullback
- L209's smaller colligation shifts every transfer coefficient, not
  only the first:
  `Btilde_n=[U_perp*B_n;U*B_(n+k−1)]`.
- Thus the surviving original transfer row is exactly
  `z^(k−1)` times the bottom row of the deflated matrix-inner
  transfer, before and after Faber reflection at `c/z`.
- Its left-oriented reflected Gram is exactly `|c|^(2k−2)` times
  the bottom deflated Gram.  The leading term is therefore
  `|c|^(2k)U*B_kB_k*U`, with no ordering guess.
- The remaining higher elliptic interface is now specifically
  **physical endpoint covariance**: prove that the prepared
  model-complement/Stein endpoint respects the same monomial shift.
  That is the rank-`m` version of L149--L150, not another transfer
  calculation.
  `proof/repeated_crabb_transfer_deflation.md`.

## NEWEST (2026-07-24): L209 geometrically deflates every later transfer grade
- If `B_j*U=0` before grade `k`, the state columns
  `WU,SWU,...,S^(k−1)WU` form an orthonormal lossless delay line.
- Removing its first `k−1` stages leaves an invariant smaller state
  space.  The compressed operator is again a partial isometry, keeps
  right defect `V`, and has left frame
  `[WU_perp,S^(k−1)WU]`.
- The first coefficient of the smaller transfer is
  `[U_perp*B_1; U*B_k]`.  Thus the surviving grade-`k` row becomes a
  genuine grade-one coefficient, providing a geometric complement
  to L208's dual preimage formula.
- This still does not identify the prepared physical ellipse jet
  with the weighted pullback of L207.  That precise weighted-jet
  identity is now the common remaining interface for both L208 and
  L209.
  `proof/repeated_crabb_transfer_deflation.md`;
  `experiments/repeated_crabb_transfer_deflation.py`.

## NEWEST (2026-07-24): L208 identifies the higher flagged range mechanism
- For every genuine transfer coefficient, an arbitrary dual Gramian
  obeys the exact ordered identity
  `YB_k−B_kA=W*(S*)^kR+sum_(j<k)B_jR*(S*)^(k−j)V`.
- On a copy projection `P` with `PB_j=0` for all earlier grades, every
  contamination term vanishes.
- Therefore every compressed channel coboundary
  `eta P{Phi(B_k*P B_k)−B_kB_k*}P` has the explicit polynomial
  preimage `C_hat=−(eta/8)Q S^k W P B_k`.
- Grade one with `P=I`, `eta=28` recovers L207 exactly.  The higher
  range/cokernel mechanism is therefore no longer mysterious.
- This is conditional progress, not the higher elliptic theorem:
  derive the actual target and coefficient from the complete
  Faber/Riemann jet, then replace discontinuous exact kernels by
  analytic Schur flags along arcs.
  `proof/repeated_crabb_transfer_flag.md`;
  `experiments/repeated_crabb_transfer_flag.py`.

## NEWEST (2026-07-24): L207 removes the elliptic rank jump
- L204's endpoint equation has the explicit solution
  `C_hat=−(7/2)QSWB_1`, `C=P^(1/2)C_hat`.
- Substitution turns the complete balanced second forcing into a
  finite ordered polynomial in `S,S*`.  The endpoint Stein functional
  converts its words into weighted transfer correlations.
- The diagonal correlations telescope to `−4B_1B_1*`; the only
  off-diagonal survivors are fourth Fourier autocorrelations of
  L201's square matrix-inner transfer, so they vanish exactly.
- This correction is polynomial, bounded, and vanishes at the Crabb
  apex.  The rank-changing pseudoinverse is not needed.
- L194 now lifts the grade-one jet jointly analytically, with the
  correctly oriented upper loss `−16B_1B_1*`.
- Next iterate the same associated-graded mechanism on
  `ker B_1*` through `B_2,...,B_L`, then merge it with L197/L199.
  `proof/repeated_crabb_elliptic_selection.md`;
  `experiments/repeated_crabb_elliptic_selection.py`.

## NEWEST (2026-07-24): L206 closes elliptic cokernel compatibility
- Normalize L204's dual observability Gramian by
  `H_Y=P^(1/2) Z_Y P^(1/2)/4`.  Its adjoint-kernel equation is
  equivalent to the exact state commutator `H_Y S=S H_Y`.
- Hence `H_YW=WY`, `H_YV=VA`, and
  `YB_n=B_nA` for every coefficient of L201's genuine matrix-inner
  transfer.  Conversely, transfer intertwining plus the telescoping
  observability identity recovers the state commutant.
- Spectral projections of `H_Y` reduce the complete colligation and
  both elliptic jets.  L203's trace identity therefore applies on each
  reducing block.  Weighting and summing proves `tr(YD_T)=0` for
  every cokernel direction, not only the scalar one.
- Finite-dimensional Fredholm now proves the oriented matrix face
  `V*XV=0`, `W*XW=−ZZ*` at every fixed equality anchor.
- The remaining gate is no longer compatibility: it is a bounded
  real-analytic selection through the rank-jumping Crabb apex, followed
  by the same construction on grades `B_2,...,B_L`.
  `proof/repeated_crabb_elliptic_commutant.md`;
  `experiments/repeated_crabb_elliptic_commutant.py`.

## NEWEST (2026-07-24): L205 separates scalar from complete equality
- Every L193 anchor has complete similarity square four, but this does
  not mean a scalar Schur function attains norm two.
- Equality in the condition-two similarity chain forces
  `f(C)Wu=Vv` between the upper and lower defect spaces.
- With L201's coefficients `B_n=W*(C*)^nV`, scalar and matrix
  Parseval turn that condition into
  `B_n v=b_n u` for every `n`.  Therefore
  `B_H(z)v=g(z)u` for a scalar inner `g` with `g(0)=0`.
- Conversely, any such constant scalar channel supplies the sharp
  scalar function `f=g`.
- Boundary unitarity makes the transfer block diagonal
  `B_H=g direct-sum B'`; the characteristic-kernel model then splits
  `C`, the exact metric, and the physical disk operator into a scalar
  full-Hardy equality block and a complementary matrix-inner block.
- This opens a weaker route tailored to the scalar conjecture:
  channel-free noncommuting anchors have a scalar gap and need not
  satisfy L204's stronger complete matrix endpoint.  The next gate is
  a uniform associated-face estimate combining channel leakage,
  L190 reflection, and L195 transverse residuals.
  `proof/repeated_crabb_scalar_channel_rigidity.md`;
  `experiments/repeated_crabb_scalar_channel.py`.

## NEWEST (2026-07-24): L204 isolates the exact elliptic cokernel
- L203's full second Stein system reduces exactly to one real-linear
  copy endpoint map.  If `G_T` is the stable Stein inverse, eliminate
  the forced parallel defect motion with
  `K=V*F_2V`, `X_0=G_T(F_2−V K V*)`.  The remaining column
  `C perpendicular V` acts by
  `M_T(C)=W*G_T(VC*+CV*)W`.
- The oriented face is now the single equation
  `M_T(C)=D_T`, where `D_T=−ZZ*−W*X_0W`.
- Its adjoint is explicit.  If
  `Z_Y−TZ_YT*=WYW*`, then
  `M_T*(Y)=2(I−VV*)Z_YV`.  Solvability is therefore equivalent to
  `tr(YD_T)=0` for every such cokernel vector.
- L203 proves the universal scalar condition; L206 now proves every
  additional reducible cokernel condition.
- A 48-case scaling audit through lengths five and multiplicities
  three finds `D_T=O(s²)` and a minimum correction `C=O(s)` as the
  equality amplitude `s` tends to zero.  Every additional cokernel
  condition vanishes.  This is strong bounded-divisibility evidence,
  not an analytic proof.
- The remaining issue is bounded analytic selection through the rank
  jump, not pointwise compatibility.
  `proof/repeated_crabb_elliptic_cokernel.md`;
  `experiments/repeated_crabb_elliptic_cokernel.py`.

## NEWEST (2026-07-24): L203 proves the scalar second elliptic face
- The second ellipse pullback coefficient is
  `E_2=2T−(T²T*+TT*T+T*T²)+T⁵`.  Together with L202's first
  defect-frame motion it gives a complete finite second Stein equation.
- The endpoint orientation has been corrected.  If
  `Z=W*dot(S)V=4B_1`, then `Z*Z=16B_1*B_1` acts on the lower/right
  defect, whereas the physical upper endpoint acts on the left defect.
  Its candidate matrix loss is therefore `−ZZ*=−16B_1B_1*`.
- The unnormalized equality dual
  `Z_0=4(P^−1−VV*)` satisfies
  `Z_0−TZ_0T*=WW*−4VV*`.  Pairing the second Stein equation with it
  proves the exact scalar identity
  `tr(W*XW)−4tr(V*XV)=−||Z||_F²`.
- A direct partial-isometry trace calculation proves the last negative
  square; it is not inferred from the numerical solver.
- The stronger matrix boundary problem
  `V*XV=0`, `W*XW=−ZZ*` was numerically solvable in every tested
  noncommuting anchor, and L206 now proves pointwise solvability at
  every equality anchor.  The linear map still changes rank at the
  Crabb apex, so a uniformly bounded real-analytic selection must be
  constructed before L194 can lift the jet.
  `proof/repeated_crabb_elliptic_second_face.md`;
  `experiments/repeated_crabb_elliptic_second_face.py`.

## NEWEST (2026-07-24): L202 identifies the grade-one elliptic normal Gram
- Block-Toeplitz endpoint elimination gives a stronger exact form for
  L193's physical metric:
  `P=2I−VV*+2WW*`, with `V,W` the initial and terminal copy
  columns.  Its spectrum is exactly `1,2,4` everywhere on the
  noncommuting equality manifold.
- The balanced disk operator is a partial isometry with right defect
  `VV*` and left defect `WW*`.
- The centered ellipse pullback has tangent
  `E_gamma=gamma T*−conj(gamma)T³`.  In balanced coordinates, the
  exact three-eigenvalue formula makes the compression of the
  Stein-defect derivative to `V^perp` vanish.
- Hence the entire first defect jet is `VC*+CV*` and can be absorbed
  by moving the rank-`m` defect frame.  No first-order metric or
  upper-endpoint motion is needed at any equality anchor.
- The adjoint calculation gives the matching left active-defect
  identity, but an immediate tangent-space audit caught an omitted
  block: these two Gram identities do **not** force full
  partial-isometry tangency.
- The sole remaining normal block is the defect-to-defect corner
  `Z_gamma=W*dot(S)V`.  It is generally nonzero and exactly
  `4gamma B_1` after the unitary identification with L201's genuine
  transfer.  Thus its oriented square is
  `16|gamma|²B_1*B_1`, the correct grade-one reflected Gram.
- This supplies the physical grade-one bridge to L201 without the
  false matrix-polynomial quotient.  The next load-bearing calculation
  is converting this normal Gram into the negative prepared endpoint,
  then lifting the argument to every `B_n`.
  `proof/repeated_crabb_elliptic_first_jet.md`;
  `experiments/repeated_crabb_elliptic_first_jet.py`.

## NEWEST (2026-07-24): L201 constructs the correct matrix-inner reflected coordinates
- L193's canonical equality metric makes
  `C=M^(1/2) A M^(−1/2)` a pure partial isometry with rank-`m`
  right and left defect projections.
- In orthonormal defect frames, the genuine transfer is
  `B_H(z)=W*(I−zC*)^−1V`.  Its characteristic function is
  `zB_H(z)`, so `B_H` is matrix inner.  The special disk chain
  `AE1=2E0` also puts `V` in `ran C`, forcing `B_H(0)=0`.
- Writing `B_H=sum B_n z^n`, the convergent matrix Dickson/Faber
  identity is exact:
  `F_cB_H(zeta+c/zeta)=B_H(zeta)+B_H(c/zeta)`.
  The negative legs have the oriented Gram
  `sum |c|^(2n) B_n*B_n`.
- At the Crabb apex `B_H=Uz^L`; nearby `B_L` stays invertible.
  Hence the reflected Gram retains a strictly positive terminal
  `|c|^(2L)B_L*B_L` leg despite arbitrary noncommuting equality
  coefficients.
- Remaining debt: identify this inner-model Gram with the actual
  fully prepared rank-`m` similarity endpoint by lifting L149's
  one-reflection and L150's orbit-complement arguments.
  `proof/repeated_crabb_inner_faber_transfer.md`;
  `experiments/repeated_crabb_inner_faber_transfer.py`.

## NEWEST (2026-07-24): L200 blocks the naïve matrix-Faber lift
- At repeated `C3`, the normalized inverse-block-Toeplitz equality
  anchor has an exact block companion and terminal matrix polynomial
  `G(z)=z²I−(Z*/2)z`.
- Dickson/Faber reflection still works coefficientwise.  That fact
  alone is not enough: the scalar proof needs `g/g#` to be inner.
- A square-zero coefficient gives an exact counterexample.  At one
  circle point, the first-column norm squared of `G(G#)^−1` is
  `1+t²`; reversing the order gives `1−t²+t⁴`.  Neither quotient is
  unitary.
- This invalidates a proof route, not the elliptic bound.  L114
  already proves the full repeated-`C3` neighbourhood by a different
  metric-stratification argument.  For an all-length replacement,
  the correct next object is L193's genuine rank-`m`
  unitary-colligation transfer, equivalently the matrix
  Schur/Levinson recursion associated with the positive block
  Toeplitz inverse Gram.
  `proof/repeated_crabb_matrix_faber_obstruction.md`;
  `experiments/repeated_crabb_matrix_faber_obstruction.py`.

## AUDIT CORRECTION (2026-07-24): L199 closes the first raw normal face, not the whole flag
- L61 handles noncommuting true-normal copy coefficients before the
  residual face: unless their first support-compression Jensen gap
  vanishes, the similarity endpoint is already strict.
- On a zero-Jensen maximal winner, Fourier independence of L115's
  normal modes forces every surviving normal coefficient to be
  scalar and removes its winner/loser cross blocks.
- Test L195's block Hardy residual on a copy vector.  Its scalar
  compression retains Hardy reflection, its scalar Frobenius energy
  is bounded by the copy Gram, and L188/L173 applies without any
  noncommutative ordering ambiguity.
- Equality kills every normal amplitude and both the rows and columns
  of the first raw residual blocks on the kernel.  This independently
  recovers the compatibility needed by L196's first promotion.
- The previous wording incorrectly iterated this raw-reflection
  argument through all of L197.  Later flag residuals are
  least-squares orthogonalized against earlier Hardy ranges and need
  not automatically retain the reflected form used by the statewise
  proof.  The transported analytic normal graph may also mix those
  earlier ranges into a later quotient.
- The remaining repeated gates are therefore (i) this later
  circular-normal flag lift and (ii) the elliptic soft coordinate and
  its operator-valued marked/Faber face over noncommuting
  inverse-block-Toeplitz equality anchors.
  `proof/repeated_crabb_circular_jensen.md`;
  `experiments/repeated_crabb_circular_jensen.py`.

## CORRECTION (2026-07-24): L198's stronger Jordan lift is not a dependency
- A square-zero copy block rigorously disproves the naïve one-sided
  gain `16P*P/C`.
- The symmetric inequality `8(P*P+PP*)/C` is algebraically absorbed
  by the reflected Hardy Gram, but the identification of that
  symmetric form with the actual repeated endpoint was not
  independently derived.  Scalar L173/L188 data cannot determine
  every noncommutative product order.
- The Jordan note is therefore retained only as a conditional guard.
  L199 supplies the valid, weaker first-face route by scalarizing
  normals on the zero-Jensen winner.

## NEWEST (2026-07-24): L197 closes the finite repeated disk flag
- Along any analytic block-disk path, the negative upper endpoint
  `S=−E` is an analytic positive semidefinite copy matrix.
- Split the positive range of its first nonzero coefficient from its
  kernel.  The active block is invertible after removing its even
  power of the path parameter; the cross block gains an order.
  Therefore exact triangular Schur congruence is analytic.
- The reduced kernel endpoint is again analytic PSD and has strictly
  higher valuation.  Repeating lowers dimension, so the flag
  terminates after at most `m` active steps.  A terminal zero block is
  an exact analytic kernel bundle.
- The first step is L195/L196.  Later steps are exactly residual
  columns orthogonalized against all earlier active Hardy ranges.
  This closes the repeated **disk-only** residual induction.
- The remaining repeated gates are now the circular-normal merger and
  elliptic matrix-Jensen/support rigidity.
  `proof/repeated_crabb_schur_flag.md`;
  `experiments/repeated_crabb_schur_flag.py`.

## NEWEST (2026-07-24): L196 promotes the first repeated residual kernel
- Use `B=H^−1` and subtract its block-Toeplitz diagonal means.  The
  remainder `N` is an exact linear transverse coordinate to L193's
  equality manifold.
- Adjacent-principal difference is an isomorphism on this zero-mean
  complement.  Its inverse is explicit integration along each block
  diagonal, so it preserves arbitrary copy-space compressions.
- At the Crabb base, the leading actual Hardy residual is
  `−(1/4)Delta(N_q)(J tensor I_m)`.  Therefore a common copy kernel of
  L195's residual Gram forces both the columns and rows of `N_q` to
  vanish on the same level-by-copy subspace.
- After recentering at the inverse-block-Toeplitz equality anchor, the
  transverse valuation on the kernel and its cross blocks strictly
  rises.  This is the first exact flag promotion.
- L197 now supplies the complete finite Schur-orthogonal induction.
  `proof/repeated_crabb_inverse_gram_kernel.md`;
  `experiments/repeated_crabb_inverse_gram_kernel.py`.

## NEWEST (2026-07-24): L195 identifies the repeated first-residual metric flag
- L193's normalized rank-`m` Hardy metric lies exactly on L194's
  zero-Stein-slack branch.  Its physical Stein defect is rank `m`
  with invertible level-zero block, hence zero range Schur complement;
  its exact lower endpoint similarly gives zero lower Schur
  complement.
- If the first terminal Hardy residual consists of copy blocks
  `s^q F_(r,c)`, the final upper endpoint is
  `−4s^(2q) sum F_(r,c)*F_(r,c)+O(s^(2q+1))`.
  This is the operator-valued form of L188's scalar Frobenius square.
- The first face is negative semidefinite, and its kernel is exactly
  the common right kernel of every residual block.  Active copy
  directions close immediately; only this common kernel advances to
  the next metric flag.
- L196 now promotes the first common kernel in exact inverse-Gram
  coordinates.  Iterating this after Schur orthogonalization, then
  merging circular-normal and elliptic support effects, remains.
  `proof/repeated_crabb_first_residual_endpoint.md`;
  `experiments/repeated_crabb_first_residual_endpoint.py`.

## NEWEST (2026-07-24): L194 gives the all-length repeated metric chart
- L105's exact repeated-`C3` metric chart extends to every Crabb
  length.  Keep the entire level-zero/range metric row `B` free and
  set `P00=I+B(C−I)^−1B*`; this makes the lower constraint exactly
  tight.
- The derivative of the Stein range Schur complement in `C` is the
  weighted diagonal recurrence
  `(LX)_(ij)=X_(ij)−a_(i−1)a_(j−1)X_(i−1,j−1)`.
  It is triangular along block diagonals and has an explicit finite
  forward inverse in every length and copy multiplicity.
- The analytic IFT therefore gives a unique convergent
  `P=P(T,B,R)` for arbitrary nearby operator, free row, and prescribed
  Stein slack `R`.  Lower feasibility is automatic; Stein feasibility
  is equivalent to `R>=0`.
- The entire upper constraint is one final-level `m x m` Schur
  endpoint.  Thus divergent forced metric coefficients are not an
  obstruction.  L195 now identifies the canonical zero-slack choice
  and its first residual endpoint; stable-kernel promotion remains.
  `proof/repeated_crabb_exact_metric_chart.md`;
  `experiments/repeated_crabb_exact_metric_chart.py`.

## NEWEST (2026-07-24): L193 explicitly classifies the repeated block-Hardy equality manifold
- L122/L183's disk factorization is operator-valued: replace every
  scalar Gram entry by an `m x m` copy block and tensor the level shift
  with `I_m`.  For `D=E0*HE0`, the normalized rank-`m` Stein defect is
  `Q=HE0 D^(−1/2)`.  Its physical kernel columns are orthonormal, so
  Berger dilation and operator Bessel give the exact sandwich
  `K<=M<=4K` near the repeated Crabb point.
- If `B=H^−1`, exact noncommutative block-row elimination gives
  `P_I(A−S)=−(B_+−B_-)(B_++B_-)^−1P_2`.  The terminal Krylov tail is
  invertible near Crabb.  Therefore the Hardy residual vanishes
  exactly when `B_+=B_-`, i.e. when **the inverse Gram is Hermitian
  block Toeplitz**.  This explicitly parameterizes the manifold,
  rather than merely invoking IFT.
- Its differential is the same block diagonal-difference map.  Its
  real rank is `((L−1)m)^2`, and the manifold has dimension
  `(2L−1)m²`.
- On residual zero, the level-zero copy space is an exact generalized
  eigenspace at `1` and the last-level copy space is an exact
  generalized eigenspace at `4`.  The rank-`m` metric therefore has
  condition square four.
- Crucially, the inverse-block-Toeplitz coefficients are arbitrary
  copy matrices.  Explicit constructions at multiplicities two and
  three have **noncommuting** coefficient blocks, condition square
  exactly four, and `M<=4K`.  The repeated equality stratum is
  therefore larger than direct sums of scalar equality anchors.
- This also repairs a subtle issue in L187's original presentation:
  the row-reversed raw residual is Hermitian only to first order, so
  a direct Hermitian-codomain IFT was not justified.  The exact
  inverse-Toeplitz factorization proves L187's scalar manifold and
  L193's block manifold without that shortcut.
- Together with L194 this suggests a cleaner repeated strategy: use
  the block Hardy residual Gram as the transverse form, promote its
  copy-space kernel into a smaller-multiplicity equality block, and
  iterate a metric flag.  The operator-valued elliptic/normal merger
  and convergence of that flag remain open.
  `proof/crabb_block_hardy_equality.md`;
  `experiments/crabb_block_hardy_equality.py`.

## NEWEST (2026-07-24): L190--L192 close the single-Crabb local chart
- L187's full equality manifold has a simpler intrinsic coordinate
  system than the old phase-palindromic section.  If
  `det(xi I−A)=xi g_u(xi)`, then
  `g_u=xi^L+2 sum_(j=1)^(L−1) u_j xi^j`, and the `L−1` complex
  coefficients `u_j` are local analytic coordinates on the equality
  manifold.  The unrestricted Faber identity has negative Hardy legs
  `r=c^L` and `w_j=c^j u_j`.
- The reflected model-complement construction extends invariantly
  over every full-Hardy equality anchor.  Defect-one colligation
  innerness kills the complete one-reflection row; the compact face is
  `−16|r|²−64 sum|w_j|²`.  Convergence gives a uniform negative
  elliptic/marked tube over the whole equality manifold, not only the
  phase-palindromic Toeplitz family.
- The missing residual/reflection block is exactly orthogonal at the
  Crabb apex.  If `F` is the first Hardy residual and `(rho,omega)` is
  the first reflected vector at the same valuation, the optimized
  initial form is
  `−4||F||_F²−16|rho|²−64 sum|omega_j|²`, with no mixed row.
  The proof uses the full normalized Gram gradient of L156's cleared
  endpoint residual: before endpoint normalization it is supported
  only on the two endpoint diagonals, which are exactly canceled by
  defect-line motion.  Positivity removes the model-gap row.
- Adding true circular normals spends no curvature twice.  L188's
  residual responses use modes `3,...,L−3`; L160's sole nonzero
  compact row uses the disjoint bottom mode and leaves a strict
  square; L163 kills every grade `>=2` compact row; L120 handles the
  pure axis.  The complete associated-graded form is therefore strict.
- Analytic curve selection in L115's disk/elliptic/normal tubular
  chart proves a full neighbourhood theorem for every fixed
  **single** Crabb block: the explicit rank-one Stein certificate has
  condition square at most four, with equality only on the full
  circular-range equality manifold modulo symmetries.
- This is not the global conjecture.  L193 now replaces the naive
  direct-sum repeated stratum by a larger operator-valued Hardy
  equality manifold; its metric-flag tube is the live frontier.
  `proof/crabb_full_local_chart_merger.md`;
  `experiments/crabb_full_equality_elliptic_merger.py`.

## NEWEST (2026-07-24): L187--L189 prove the first-residual tube
- The exact disk-chart Hardy residual has a finite
  `(L−1) x (L−1)` coordinate matrix `Psi(H)`.  At the Crabb point,
  `D Psi(E)_(r,c)=E_(c+1,L−1−r)−E_(c,L−2−r)`.  After row reversal
  this is the Hermitian diagonal-difference map, with real rank
  `(L−1)^2` and kernel exactly the Hermitian Toeplitz space.
- The exact inverse-Gram reflection factorization now identifies this
  manifold as `H^−1` Hermitian Toeplitz.  It is curved in the original
  `H` coordinates and has physical dimension `2L−2`.  Residual zero
  makes the canonical Hardy upper endpoint exactly four times the
  lower endpoint.  The defect-one model's characteristic Blaschke
  product simultaneously has norm exactly two, so these are genuine
  scalar and similarity equality points.  L123's phase-palindromic
  family is an explicit lower-dimensional section.
- The complete pulled ambient derivative vanishes on this whole
  manifold.  Independent nonlinear solves in lengths `3,...,6`
  produce equality points with non-Toeplitz `H` (but Toeplitz
  `H^−1`), condition square and
  characteristic-Blaschke norm square both `4`, top-vector overlap
  `1`, and all real/imaginary matrix-unit derivatives at roundoff.
- If an arbitrary analytic full-disk path first leaves the equality
  manifold in order `m`, with Hardy residual coefficient `F`, then
  the canonical base deficit is universally `4||F||_F^2`.  Every
  leading circular-normal response is the weighted anti-diagonal
  projection of the transpose-skew part `(F−F^T)/2`.
- Anti-diagonal Cauchy--Schwarz exactly absorbs those responses using
  L173's flux curvature; its positive null lift makes the actual
  first-residual face strict.  This holds in **every order**, not only
  the sixth/eighth faces.  Exact audits cover qualitatively different
  first residual orders `1,...,5`, including the adversarial weighted
  pair whose fourth residual was canceled exactly.
- Analytic curve selection now proves the nonlinear tube over L122's
  complete general-`H` disk chart and the transported true
  circular-normal fibres.  L190--L192 subsequently merge the
  elliptic and marked/compact charts and close the complete
  single-Crabb local quotient.
  `proof/crabb_full_disk_leading_residual_tube.md`.

## NEWEST (2026-07-24): L185--L186 classify and lift the sixth kernel
- The cubic Hardy residual has a closed all-index formula in terms of
  L182's Wronskian fluxes
  `S_t=sum_(i+j=t)(j−i)z_i conjugate(z_(n−1−j))`.
  A two-form argument proves `R_L=0` exactly when every low `S_t`
  vanishes.
- If `p=sum z_i x^i` and
  `q=sum conjugate(z_(n−1−i))x^i`, those equations say
  `p q'−p' q=c x^(n−2)`.  After removing a monomial gcd, the rational
  map `p/q` has only the critical points zero and infinity.
  Riemann--Hurwitz therefore classifies the kernel as the union of
  the phase-palindromic equality cone and the coordinate planes
  supported on one reversal pair `{j,n−1−j}`.
- On a strict reversal pair `z_j=a,z_(n−1−j)=b`, with
  `d=n−1−2j`, the fourth Hardy residual has the single skew pair
  `8d(|a|²−|b|²)²/L`.  Hence
  `D_(8,L)=512d²(|a|²−|b|²)^4/L²`.
- Every true circular-normal response vanishes through quartic order
  on these planes: the fourth residual lies on the central
  anti-diagonal, while nonzero normal modes use the shifted
  anti-diagonals.  The eighth face is therefore strictly positive
  away from `|a|=|b|`, exactly the phase-palindromic intersection.
- This closes the homogeneous kernel fallback.  L187--L189
  subsequently replace its finite-order ladder by the full Hardy
  equality manifold, universal first-residual absorption, and the
  nonlinear general-`H` disk/circular-normal tube.
  `proof/crabb_full_disk_kernel_eighth.md`.

## NEWEST (2026-07-24): L183--L184 prove the all-size sixth face
- L183 gives a sharp kernel-observability theorem.  If `w(T)≤1`,
  `Tv=0`, and `P_v=sum (T*)^jvv*T^j`, Berger's unitary
  `2`-dilation and Bessel's inequality give
  `P_v≤4||v||²I`.
- L122's general-`H` disk chart has an explicit Hardy realization of
  that dilation.  Along L176's recentered path, its orbit-complement
  residual vanishes through degree two and its cubic coefficient is a
  skew matrix `R_L(z)`.
- The complete canonical base deficit is exactly
  `D_(6,L)=8||R_L||²`.  L182's flux coordinate `T_(L,k)` is exactly
  one quarter of each of the two weighted grade-`±k`
  anti-diagonals of `R_L`.
- Cauchy--Schwarz on the disjoint anti-diagonals proves the sharper
  inequality
  `D_(6,L)≥256 sum_k |T_(L,k)|²/binom(L−k,3)`.
  The right side is the completed cubic gain using only L173's flux
  curvature; the actual curvature also has its strictly positive
  null-lift term.  Therefore the complete sixth-order Schur face is
  nonnegative in every size.
- At `L=6` this proves the formerly conjectural sharp complex
  `P_6≥9|C_6|²`; no large SOS is needed.  The sixth kernel is exactly
  `R_L(z)=0`.  The terminal two-coefficient edge saturates the stronger
  flux inequality but is strict for the actual face because of the
  null lift.
- Unrestricted-symbol transfer checks pass through length eight;
  exact endpoint-factor checks pass through length ten; terminal
  equality checks pass through length twelve.
- L185--L186 subsequently classify `R_L=0` and prove its eighth-order
  fallback.  The remaining task is the nonlinear singular blow-up.
  `proof/crabb_full_disk_sixth_hardy_factor.md`.

## NEWEST (2026-07-24): L182 proves the all-size cubic response
- Put `S_t=sum_(i<t−i)(t−2i)W_(i,t−i)`.  Every cubic
  true-normal response after L176 recentering is one triangular
  interval-flux transform:
  `G_(L,k)=16(4k−1)/L²` times the circle-selected sums of
  `(a−t−1)z_aS_t/(t+2)` and their conjugates.
- The coefficient is derived from the degree-three
  characteristic/reversed-Horner/inverse-Riemann recurrence.  After
  inserting L176's mean-zero pulse, its intrinsic Pluecker coefficient
  has constant first difference and solves to
  `(a−t−1)(t−2i)/(t+2)`.
- This proves simultaneously that modes `3,...,L−3` are the complete
  active range and specializes exactly to L178's highest-mode formula
  and both L181 formulas.
- A separate 14-direction exact recovery at `L=9` has full flux ranks
  `11,10,9,8` and reconstructs every coefficient without supplying the
  formula.  The independent complete series checker passes two dense
  complex directions in each length `6,...,11`.
- This closed the response half of the arbitrary-size sixth face.
  L183--L184 subsequently factor the endpoint base deficit through a
  Hardy residual and prove positivity; L185--L186 then close its
  kernel lift.  Only the nonlinear tube remains open.
  `proof/crabb_full_disk_cubic_response.md`.

## NEWEST (2026-07-24): L181 proves the complex `p=8` sixth face
- The complete length-seven Schur residual is
  `−2δ_6−|G_3|²/(4b_3)−|G_4|²/(4b_4)`, with exact curvatures
  `b_3=773/9604` and `b_4=901/9604`.
- A dense generic-complex characteristic expansion became
  computationally pathological.  Exact sparse polarization replaces
  it: 364 one-/two-/three-coordinate points have cubic evaluation rank
  `364` and quadratic rank `78`.
- All 4,368 exact true-normal polarization evaluations prove that the
  quadratic response vanishes, modes three and four are the only cubic
  modes, and both proposed Pluecker formulas hold coefficientwise.
- On 37 selected coordinates of
  `z tensor (z wedge J conjugate(z))`, the residual is exactly nine
  positive rational squares in the real parts plus six in the
  imaginary parts.  Independent endpoint regeneration matches all 299
  polarized sextic coefficients; the Gram ranks are `9+6`.
- L180--L181 now close the first two active complex sizes.  Their
  literal shared tail factors point to an all-size interval-flux LDL
  recurrence, but this is not yet an induction theorem.  Derive that
  recurrence rather than starting another raw fixed-size SOS.
  `proof/crabb_full_disk_complex_L7_sixth_certificate.md`.

## NEWEST (2026-07-24): L180 proves the complex `p=7` sixth face
- On A126's first active complex slice, the exact required inequality
  is `P_6−1089|C_6|²/290>=0`.
- Twenty selected coordinates `q_r=−z_a(z wedge J conjugate(z))_ij`
  give an exact rational identity with six positive squares in
  `Re q` and four positive squares in `Im q`.  The corresponding Gram
  ranks are six and four.
- A separate generic-complex characteristic/Riemann audit verifies all
  ten true-normal polarizations: the quadratic response vanishes,
  mode three is the sole cubic mode, and it equals `−22C_6/45`.
- The checker independently reconstructs `P_6` from the endpoint
  recurrence, verifies conjugation invariance, and matches all 137
  polarized sextic coefficients.  No floating solver remains.
- A formal Hermitian lift of L179's real Gram is false: the separate
  global-phase sectors in the ten-square identity are essential.
- This closes the complete complex sixth face at `p=7`; it does not
  prove the stronger constant `9`, the arbitrary-size block theorem,
  the higher-order kernel lift, or the nonlinear tube.
  `proof/crabb_full_disk_complex_sixth_certificate.md`.

## NEWEST (2026-07-24): L179 proves the real `p=7` sixth face
- On A126's real first-active slice, the exact required inequality is
  `P_6−1089C_6²/290>=0`.
- The real phase-palindromic union has a 21-dimensional homogeneous
  cubic equality ideal.  On an explicit binomial basis `q`, an exact
  rational rank-seven matrix gives
  `P_6−1089C_6²/290=q^TQq`.
- The Gram matrix is represented as `Q=FS^(−1)F^T`.  All seven leading
  principal minors of the symmetric pivot core `S` are explicitly
  positive, so the certificate is PSD with no numerical-solver
  premise.  Exact expansion checks every polynomial coefficient.
- This closes all real directions at `p=7`; L180 subsequently closes
  the complex slice.  Neither result proves the stronger conjectural
  constant `9`, the all-size sixth block, or the nonlinear tube.
  `proof/crabb_full_disk_real_sixth_certificate.md`.

## NEWEST (2026-07-24): L178 closes the terminal sixth-order edge
- On the recentered full-disk path with only its last two Toeplitz
  coefficients nonzero, put `k=L−3`, `z_k=a`, and `z_(k+1)=b`.
  The exact canonical sixth-order base deficit is
  `256k²|a|⁴|b|²/(k+2)²`.
- The sole cubic true-normal response is in the highest active mode:
  `16k(4k−1)a²conj(b)/(L²(k+2))`.  L173's exact `r=3` curvature
  therefore makes its completed-square gain/base ratio
  `6(4k−1)²/[6(4k−1)²+169k(k−1)(k−2)]`.
- The ratio is strictly below one and decreases from `121/290` at
  `L=6` toward zero.  Thus the terminal edge that makes A126's stronger
  real `p=7` inequality sharp is not an equality edge of the actual
  Schur problem; the positive `169` null-lift term retains a strict
  margin.
- Exact endpoint and characteristic/Riemann recurrences check the
  complex formulas through `L=15`.  The full sixth-order tensor block,
  its kernel, and the nonlinear tube remain open.
  `proof/crabb_full_disk_terminal_sixth_face.md`.

## NEWEST (2026-07-24): A126 isolates the recentered sixth-order Schur face
- The initial `p=5,6` finite data did not extrapolate: after L176's
  exact disk recentering, the cubic true-normal response first survives
  at `p=7`, in support mode three.  Thus the next universal face is
  sixth order, not generically eighth order.
- For `p=7`, with `W=z wedge J conjugate(z)`, the first complex response
  is exactly checked as
  `−22(6z_3W_03−5z_4W_02+2z_3W_12)/45` on five unrelated rational
  complex rays.  Its L65 curvature is `145/2592`.
- The canonical disk-base endpoint excess vanishes through degree five
  and is generically negative at degree six.  A real `p=5` symbolic
  slice factors as
  `−128 a_2²(a_1−a_3)²(2a_1²+a_2²)/9`, exposing additional
  terminal-only strata whose first base term is eighth order.
- On the `p=7` slice, the needed Schur inequality is
  `P_6>=1089|C|²/290`.  L179 proves it over the real slice and L180
  over the full complex slice.  A Gram SDP and 2,000 complex random
  rays plus 20 BFGS searches independently locate the much stronger
  sharp-looking constant `9`; equality occurs when only the last two
  Toeplitz coefficients remain.  That stronger bound remains discovery
  evidence, not an exact SOS theorem.
- Exact canonical ratios remain favorable through `p=10`; five exact
  `p=7` rays have maximum about `0.277`.  A separate 30-record
  optimized nonlinear `p=7` probe has positive residual throughout and
  maximum finite-scale ratio `0.3375`.
- The load-bearing target is now the all-size block positivity of the
  sixth face on `z tensor (z wedge J conjugate(z))`, followed by an
  eighth-order analysis on its kernel.
  `proof/crabb_full_disk_sixth_face.md`.

## NEWEST (2026-07-24): L177 proves the complete full-disk response identity
- The unprojected degree-two characteristic/Riemann response in every
  active paired circle mode is the interval vector
  `x_i=16(chi_i−(t−2i)q_0/L)` paired with L65's reduced ambient
  coordinate.  All inactive nonzero modes vanish by circle character.
- L176's exact curvature calculation already gives
  `K_(m,k)s_i=x_i`.  The paired ambient response is therefore twice
  the polarized curvature of L176's explicit disk-tangent correction.
- The terminal diagonal row is the unique covector satisfying
  `A_L^Tq=4 Delta^TDelta d` and `q^Tbeta=0`; this is likewise twice the
  grade-zero weighted-shift curvature.
- Hence `g_2=2 C_pE_z` on the complete ambient matrix space and the
  L175 face completes exactly to
  `−<Y−E_z,C_p(Y−E_z)>≤0`.  This proves the homogeneous full-H face in
  every size.
- An independent exact characteristic/Riemann engine checks every real
  and imaginary matrix unit through `p=7`; a second run is
  byte-identical.  The nonlinear full-circular-range tube is still
  open. `proof/crabb_full_disk_response_identity.md`.

## NEWEST (2026-07-24): L176 proves the full-disk correction isometry
- L175's numerically selected correction now has an explicit all-size
  formula.  On each anti-diagonal of `h wedge J conj(h)`, it is the
  universal interval pulse
  `2 1_[i+1,t−i]−2(t−2i)/(t+2)`, placed on the corresponding
  Hermitian disk-chart offset.
- Differentiating L122's physical chart sends these pulses to nested
  path-flux vectors in L65.  Solving L65's singular/rank-one path
  inverse gives the exact Gram matrix `64I`; the terminal diagonal
  weighted-shift face gives `32I`.
- Accounting for reflected versus terminal Pluecker coordinates proves
  exactly in every size that the correction has curvature energy
  `32Q`.  Exact SymPy regeneration passes through length 30 and repeats
  byte-identically.
- This is the curvature-energy half of L175.  L177 subsequently proves
  the ambient identity `g_2=2 C_p D X[B_2(h)]` and completes the
  homogeneous face. `proof/crabb_full_disk_correction_isometry.md`.

## NEWEST (2026-07-24): L175 isolates the exact-looking full-disk recentering
- The full-H leading face is not strict.  For the quadratic ambient
  response `g_2(h)`, L65 curvature `C_p`, and the physical tangent
  `T_p` of L122's general-H disk chart, 16 structured/random records
  in `p=4,...,7` give
  `g_2 perpendicular ker(C_p)`,
  `g_2 in range(C_p T_p)`, and
  `(1/4)g_2^T C_p^dagger g_2=32Q` to `1.4e−9`.
- Thus completing every strong direction appears to cancel the entire
  Toeplitz quartic, with the maximizing class tangent to the exact disk
  manifold.  This sharpens L160's earlier finite-difference observation
  and is compatible with L173, whose smaller true-normal quotient has a
  strict margin.
- After recentering by the corresponding Hermitian correction, the
  optimized rank-one base deficit is approximately sixth order in
  `p=5,6` (eighth order in `p=4`), while the genuine circular-normal
  Schur gain is approximately eighth order.  All 36 multidirection
  residuals are positive; the largest ratio is `0.223` at scale `0.15`
  and it decays toward zero with scale.
- Status discipline: L176--L177 prove the tight homogeneous full-H
  identity.  The finite-scale exponents and nonlinear full-disk tube
  are **not proved**.  The result should now feed a tubular splitting
  over the complete circular-range manifold, not another strict
  full-strong quartic attempt.
  `proof/crabb_full_disk_weighted_face.md`.

## NEWEST (2026-07-24): L174 proves the nonlinear disk/circular-normal slice tube
- Over L122's exact Toeplitz disk chart, transport the `2p−4` true
  coercive circular-normal representatives by the disk metric square
  root and evaluate L118's analytic optimized rank-one envelope.
- L65's negative Hessian gives uniform fibre concavity and an analytic
  maximizing graph `y_*(z)`.  L123's exact equality metric together
  with L162's ambient stationarity pins `y_*(u)=0`, `H(u)=0`, and
  `DH(u)=0` on the entire phase-palindromic cone.
- The pure disk quartic is at most `−32Q`.  A111/L157 identify L173's
  quadratic Pluecker response as the optimized-envelope response, and
  L173's strict Schur margin leaves `H_4<=−c_LQ`.
- L155's determinantal-ideal argument removes subquartic terms.
  L124's best-phase splitting then supplies a two-regime blow-up:
  ordinary radial domination away from the cone and normal Taylor
  domination near it.  This yields the nonlinear estimate
  `Gamma<=−a_LQ−b_L||y−y_*(z)||²`.
- Scope is important: this is a rigorous Toeplitz-disk × coercive-normal
  slice theorem, not yet a tube over every non-Toeplitz point of the
  Lewis--Overton circular-range manifold.  The elliptic soft normal and
  compact reflected variables also remain to be merged.
  `proof/crabb_disk_circular_normal_tube.md`.

## NEWEST (2026-07-24): L173 proves the sharp Crabb-apex normal Schur face
- The tempting stronger claim that the quadratic disk-to-circular-normal
  gradient vanishes is false from `L=6` onward.  Its first exact value is
  `-253/45000`.  The corrected checker preserves this failed shortcut
  rather than hiding it.
- The complete quadratic response has the exact exterior-square form
  `G_(L,k)=8(4k-1)/L² sum_(i+j=L+k)(j-i)
  (h_i conj(h_(L-j))-h_j conj(h_(L-i)))`.
  Thus it factors through `h wedge J conj(h)`, and different support
  modes occupy orthogonal Plücker anti-diagonals.
- In Takagi coordinates `Q=4||p||²`.  Either real polarization of mode
  `k`, with `r=L-k`, has squared row norm
  `128(4k-1)² binom(r,3)/L⁴`.
- Solving L65's singular path kernel by cumulative flux gives the exact
  physical-mode curvature
  `b=L^(-4)[(4k-1)²r(r-1)(r-2)/24
  +(2/3)k(k-1)(k-2)(r+1/4)²]`.
- The first curvature term is exactly `||m||²/512`; the positive
  rank-one null-lift term makes the inequality strict.  Orthogonality
  therefore makes the full completed normal gain strictly smaller
  than the disk deficit `32Q=128||p||²` in every fixed size.
- Audit correction: coefficient-gauge support directions must be
  conjugated by the Crabb metric square root before evaluation in
  L65's physical Hessian.  Omitting this gives a close but incorrect
  curvature formula.  The corrected formula retains the theorem,
  now with the strict margin supplied by the null lift.
- This proves the sharp homogeneous apex face, not yet the nonlinear
  tubular patch.  The remaining step is analytic: transport the normal
  frame over the Lewis--Overton circular-range manifold, combine the
  strict apex margin with L162's exact positive-ridge stationarity,
  then merge L160 and L163/L172.
- Exact Gaussian-rational response regeneration through `L=7` and an
  independent all-mode L65/finite-difference regeneration through
  `L=12` both pass; repeated runs are byte-identical.
  `proof/crabb_circular_normal_plucker_schur.md`.

## NEWEST (2026-07-23): L172 closes L163 by an explicit endpoint residue
- The gap found in the L171 audit is now closed without the discarded
  slogan that every positive Hardy mode is automatically projected
  out.
- Cauchy functional calculus and L156 identify the remaining direct
  endpoint/cofactor row with the `z^(-1)` coefficient of
  `dot(B)_(w_k) dot(R)_s`.
- L140 gives exactly two reflected Blaschke powers: `L-k` and `L+k`.
  L168--L170 give the relative endpoint discrepancy at power
  `m=L+2-k`; multiplication by `F_0=2z^(-L-1)` makes its absolute
  power `1-k`.
- The only possible product powers are therefore
  `L-2k+1` and `L+1`.  For `2<=k<=L/2`, both are at least one, so
  neither is the Cauchy residue power `-1`.  The central fold has
  first power exactly one.
- Together with L171's zero-reflection and sparse
  singular-Hessian reductions, this proves L163 in every size.
  Grade one is untouched: its coupling-index argument and Schur
  tangent are different, and L160/L165 remains nonzero.
- The associated checker regenerates all 36 `(L,k)` faces through
  `L=14`, including the residue powers; a second run is byte-identical.
  A separate end-to-end bilinear-jet checker reconstructs the full
  frozen prepared Blaschke norm: grades two and three are exactly zero,
  while three grade-one controls reproduce `-4(5L-1)/L`.
  `proof/crabb_dual_endpoint_residue.md`.

## NEWEST (2026-07-23): L171 reduced the all-grade row to one endpoint kernel
- At this checkpoint L163 was **not yet proved** all-size.  The earlier
  “positive Hardy mode is annihilated” sentence was not established
  by L156 and has been removed rather than promoted into a lemma.
- The rigorous part is substantial.  L149 reflection count and L162
  kill the zero-reflection sector; L166 transfers the one-reflection
  row to the prepared Blaschke side.
- The ordinary dual singular Hessian is exactly sparse.  L143's
  reflected tangent couples the top endpoint only through `L-k`;
  `D(C^L)[E_d]` couples it only through `k-2`.  Since
  `L-k>=k>k-2`, both the direct quadratic and eigenvector-coupling
  terms vanish in every size.
- L168--L170 split the only remaining direct endpoint response into
  a logarithmic-inner tangent and an optimized-defect relative mode.
  At this stage the precise remaining gate was to evaluate
  L149/L156's full preparation/endpoint functional on that latter
  mode.  L172 now does so by explicit Cauchy-residue selection.
- Grade one remains the mandatory discriminator: `m=L+1` is the
  first terminal/feedthrough alias and L160/L165 is nonzero.
- Exact grades two through four and floating grades through seven
  satisfy the target.  The independent grade-four run has six
  nonzero terms whose total cancels, ruling out termwise support.
  L172 supplies the all-size endpoint proof.
  `proof/crabb_all_grade_normal_bridge.md`.

## NEWEST (2026-07-23): L170 removes the moving feedthrough without root tracking
- For L169's degree-`p` characteristic inner function, put
  `d=theta(0)` and take the first Schur iterate
  `B=(theta-d)/(z(1-conj(d)theta))`.
- This is a real-analytic finite Blaschke product of degree `p-1`.
  When `d=0` it is exactly L149's prepared factor `theta/z`.
- At the Crabb point, L167 gives a sharp discriminator.  The
  grade-one tangent `-1+z^(2p)` maps to `dot B=0`.  Every grade
  `k>=2` maps to
  `dot B=-kz^(k-2)+kz^(2p-k)`, hence
  `dot B/B=-k(z^(-m)-z^m)`.
- This independently recovers the negative of L168's endpoint inner
  tangent and explains why grade one has no matching degree-`p-1`
  Blaschke motion.
- Unlike naive division by `z`, the exact Schur formula remains
  analytic when A115's later characteristic constant becomes
  nonzero.  No eigenvalue or zero is selected.
- L171 uses this leading discriminator and sparse dual support to
  reduce L163 to one differentiated endpoint-functional identity.
  The all-order formula remains the guard against reintroducing a
  fixed-zero assumption.
  `proof/crabb_characteristic_schur_step.md`.

## NEWEST (2026-07-23): L169 packages every mixed characteristic jet into one inner transfer
- Any rank-one Stein pair
  `P-T*PT=qq*`, `P^(-1)-TP^(-1)T*=rr*` balances to a
  contraction with one-dimensional left and right defects.  Completing
  those columns gives a unitary colligation.
- Its scalar transfer has the exact determinant form
  `theta(z)=d+zq*(I-zT)^(-1)r
  =omega det(zI-T*)/det(I-zT)`.  Hence it is inner and every analytic
  variation satisfies `Re(dot(theta)/theta)=0` on the circle.
- L164 makes this canonical on the optimized persymmetric branch:
  `r=sqrt(alpha)Jq`.  No second defect series is needed.
- The scalar feedthrough `d` is essential: it carries the moving
  constant characteristic term detected by A115.  Thus this exact
  identity survives precisely where the fixed-zero shortcut failed.
- L168 is the first Crabb tangent of the colligation formula.  L171
  combines it with L162 and a sparse dual Hessian calculation; L172
  subsequently closes the remaining direct endpoint functional.
  `proof/crabb_colligation_transfer.md`.

## NEWEST (2026-07-23): L168 exposes the leading inner endpoint tangent
- Put `m=L+2-k`.  For the first inverse-Riemann pullback `E` of the
  eligible grade-`k` circular normal, exact weighted-path counting
  gives
  `D log(e0*(zI-C)^(-1)eL)[E]=k z^(-m)+beta_(L,k)z^m`.
- The positive coefficient is `0` for `k=2`, `(L-9)/(4L)` for
  `k=3`, and `k-3-2(k-2)^2/L` for `k>=4`.
- The symmetrized defect endpoint factor is `1` for `k=2` and
  `sqrt(2)` thereafter.  Substitution of L163's independently derived
  normal/defect coefficient gives the exact all-size identity
  `beta+k=delta gamma_(L,k)`.
- Therefore the normal response minus its optimized defect response
  is exactly `k(z^(-m)-z^m)`: an anti-self-reciprocal
  logarithmic-inner tangent, purely imaginary on the circle and with
  zero real mean.  This is the first concrete scalar bridge from
  L163 to L149 rather than a finite-grade cancellation.
- The result is deliberately only leading order.  L171 shows that it
  is the only endpoint response left after reflection-count separation
  and sparse dual support; L172 subsequently evaluates the complete
  differentiated L156 norming functional on it.
- A115 remains an important guard: later characteristic factors do
  move, but their weights do not enter L171's one-reflection face.
  `proof/crabb_leading_endpoint_transfer.md`.

## NEWEST (2026-07-23): L166 transfers the compact normal row to the sharp dual
- Let `U` be L118's optimized Stein upper envelope and `R` the
  prepared Blaschke norm square, frozen in the strong-normal variable.
  For every fixed nonzero weighted scale, `Delta=U-R>=0` and
  `Delta=0` on the exact elliptic axis.
- The amplitude/normal Hessian of `Delta` is therefore PSD.
  L142--L145 make its pure grade-`k` reflected diagonal
  `o(epsilon^(2k+2))`; a normal inserted at weight `k+1` makes its
  strong diagonal `O(epsilon^(2k+2))`.
- PSD Cauchy--Schwarz forces the mixed gap to be
  `o(epsilon^(2k+2))`.  Hence the leading compact
  one-reflection/normal coefficient of the optimized envelope equals
  the corresponding sharp prepared-Blaschke coefficient.
- This is an all-size theorem and removes all optimizer jets from the
  one-reflection part of L163.  L171 proves that the common dual row
  has no ordinary singular-Hessian contribution; L172 closes its
  direct endpoint/cofactor term.
- Audit warning: do not implement that bridge by freezing the zero
  characteristic root.  Grade two preserves the determinant through
  its face, but exact grade three has nonzero normal determinant
  derivatives already at weights five, seven, and eight while its
  condition cross still vanishes.  The characteristic polynomial must
  move inside the transfer identity.
- At the associated leading strong weight there is a clean all-size
  discriminator: for grade `k`,
  `D det(zI-C)[E_d]=-k z^(k-1)`.  It follows directly from the
  nilpotent adjugate and the single trace-closing subdiagonal.  Thus
  grades at least two retain one zero factor on the compact face,
  whereas grade one changes the constant term immediately.  Any inner
  bridge may use this leading fact, but not promote it to later jets.
  `proof/crabb_dual_normal_transfer.md`.

## NEWEST (2026-07-23): L165 reduces L163 to one half-order adjoint telescope
- For the optimized rank-one metric, set
  `G=v_+v_+*/lambda_+−v_-v_-*/lambda_-` and solve
  `Z−TZT*=G`.  Exact Stein adjointness gives
  `D log kappa[E]=2 Re tr(ZT*PE)`.
- Defect stationarity and scale invariance give `Zq=0`, so no
  differentiated optimizer jets occur in the normal derivative.
  L164 also makes the operator gradient `PTZ` persymmetric.
- The two Stein equations imply the sharper endpoint-flux identity
  `T^*(PTZ)−(PTZ)T^*=v_+v_+^*−v_-v_-^*`.  Thus every interior
  contribution is a commutator divergence with only two endpoint
  sources, providing the likely algebraic origin of the telescope.
- Canonically decomposing the normal series as `dot T=[T,X]+H`,
  with `X` first-row normalized and `H` supported on the bottom
  companion row, splits the target into an endpoint-basis term and a
  characteristic-polynomial term.  This decomposition is exact and
  coefficientwise, not a numerical fit.
- For reflected grade `k` and `d=k+1`, the weight-`2d` normal
  coefficient is the `d+1`-term convolution
  `2 Re sum_(j=0)^d <[epsilon^j]PTZ,
  [epsilon^(2d-j)]dot T>`.  The adjoint state is needed only through
  weight `d`, rather than two endpoint states through `2d`.
- The exact checker distinguishes all three current regimes:
  grade one is nonzero (`D kappa=-56/3` in size four), every grade-two
  summand vanishes separately, and grade three has five nonzero
  terms in `Q(sqrt(2))` which telescope exactly to zero.
- In grade three the endpoint-basis and characteristic sums are
  respectively `-(-1060+683sqrt(2))/12` and its negative.  Grade one
  has the nonzero sum `-7/3`.  Therefore neither support
  disjointness nor endpoint flux alone can prove the rule.
- With equality and ellipse amplitudes independent, grade three's
  `a^3c` and `ac^3` monomials each cancel separately.  The first
  central grade-four case is now also exact: its six nonzero adjoint
  pairings sum to zero, with boundary and characteristic totals
  `2(-607+439sqrt(2))` and its negative.
- L171 splits this convolution into marked sectors and removes all
  ordinary dual Hessian terms, but the all-grade direct endpoint
  balance remains open.  The exact records remain independent
  regressions, and grade one remains the mandatory discriminator.
  `proof/crabb_adjoint_normal_recurrence.md`.

## NEWEST (2026-07-23): L164 proves reciprocal reversal of the optimized metric
- For every real persymmetric normalized operator near a Crabb block,
  the locally unique L118 rank-one metric satisfies
  `J P^(-1) J = alpha P`, with
  `alpha=(det P)^(-2/p)`.
- The proof is abstract and all-size.  Rank-one inverse-Stein duality
  makes the reversed inverse another rank-one metric for the same
  operator and condition number; projective defect normalization is
  an analytic involution, so uniqueness of the optimized defect makes
  it a fixed point.
- The first numerical residuals were not failures: a defect series
  optimized through jet `j` satisfies the identity through jet `j`
  and generally fails first at the omitted jet `j+1`.  Solving one
  more exact jet removes the entire previous residual.
- This pairs the two endpoint eigenvalue series reciprocally and
  reduces L163 to one endpoint, or to a single adjoint-gradient
  coefficient.  It does **not** prove L163: self-duality controls the
  sum of relative endpoint derivatives, not the difference that must
  vanish.
- The fixed-point equation is also a half-order optimizer solver.
  Defect jet `j` follows from self-duality at degree `j`; its
  fixed-point linearization is nonsingular at the Crabb block.  The
  first four size-seven grade-three jets take seconds and agree with
  the direct endpoint optimizer where independently compared.
  `proof/crabb_reciprocal_reversal.md`.

## NEWEST (2026-07-23): A111 falsifies the `O(Q)` ambient-gradient shortcut
- The full target is exactly false in size four:
  `Q=25281x^4/15625000000` while the normalized
  characteristic-dual derivative in `E_(1,2)` is
  `-424x²/15625`.
- This survives projection onto the transported support-mode
  covectors `3,...,p`: their computed rank is exactly `2p-4`, and the
  projected norm divided by `sqrt(Q)` tends to a nonzero constant.
- Dimensions `p=4,...,7`, four transverse scales, and FFT resolutions
  `1024,...,8192` are stable; the equality controls are `~1e-14`.
- L157 still makes the optimized-upper/dual gradient difference
  `O(Q)`.  It therefore cannot cancel the leading `sqrt(Q)` term.
  The sharp covariant Schur cancellation in L163 is necessary rather
  than a technical detour.
  `proof/crabb_disk_ambient_gradient_division.md`.

## NEWEST (2026-07-23): L163 isolates the all-grade circular-normal proof gate
- For reflected grade `k>=2`, the only character-eligible true
  circular normal has mode `m=L+2-k`.
- At the Crabb point, the first normalized normal/defect-column
  coefficient is `2` for `k=2`,
  `sqrt(2)(13L-9)/(8L)` for `k=3`, and
  `sqrt(2)((2k-3)/2-(k-2)^2/L)` for `k>=4`.
  The complete residual Stein chain also has a closed all-size
  formula, regenerated exactly through `L=12`.
- Audit correction: this first column does not include the mixed
  strong coefficients at weights `d+1,...,2d` created by the inverse
  Riemann, Stein, and endpoint recurrences.  The nonzero individual
  grade-three endpoint shifts show those terms are real and must be
  retained.  The earlier all-grade PROVED label was premature.
- The obstruction is exact, not just diagnostic: for grade three the
  ordinary L65 Hessian cross is strictly positive (generically
  `6m(m-1)(m-2)(2-sqrt(2))/(m+1)^2`, with one explicit
  `L=10` terminal fold), although the complete checked physical cross
  vanishes.  Higher covariant jets must cancel this nonzero base term.
- L161 proves the complete statement for grade two.  The corrected
  grade-three central collision and next length cancel exactly, and
  floating scans through grade seven support the all-grade rule.
  These are finite evidence, not an all-size proof.
- A stronger floating audit leaves the final defect jet free: through
  grades two to five, the complete normal derivative is exactly the
  same closed `gamma_(L,k)` times the matching final-defect derivative.
  This identifies the likely theorem; defect stationarity would kill
  it immediately once the covariant recurrence is proved.
- L171 subsequently reduced the coefficient further by
  reflection-count separation and sparse dual support.  The direct
  differentiated endpoint functional was then closed by L172; only
  the separate uniform tubular lift remains.
  `proof/crabb_all_grade_normal_selection.md`.

## NEWEST (2026-07-23): L162 proves all-size ambient stationarity on the equality ridge
- At every positive phase-palindromic disk-equality anchor, fix L123's
  coefficient metric `K` and defect `q`, continue its rank-one Stein
  solution under an arbitrary ambient perturbation, and apply the first
  Schwarz/Riemann correction.  The resulting certificate condition has
  zero derivative in every complex matrix direction.
- The proof is exact.  L123's finite Blaschke lower bound touches the
  Stein upper certificate at four, identifying the derivative as
  `8 Re(q* DB(A)[E]e_L)`.  L156's endpoint resolvents turn this into a
  boundary integral with density `1/|g#|²`.  A finite endpoint transfer
  function and a reflected-polynomial kernel identity show that the
  Schwarz correction has exactly the same real derivative.
- Exact nonreal rational audits reach `p=6`; generic symbolic parameters
  reach length four; every real/imaginary matrix unit and additional
  near-boundary random anchors pass through `p=9`.  Non-palindromic disk
  controls are not stationary.
- This supplies the critical-ridge half of a tubular/Morse--Bott
  argument.  It does not prove a negative normal Hessian or the final
  neighborhood theorem.  L171--L172 now close L163's all-grade mixed
  selection; the prepared uniform full-disk tubular chart remains.
  `proof/crabb_equality_ambient_stationarity.md`.

## NEWEST (2026-07-23): L160--L161 isolate the first circular-normal faces
- Optimizing all of L118's strong variables over the Toeplitz disk chart
  is false: part of that space recenters onto the larger exact
  circular-range manifold.  The correct L115 quotient has `2p-4` real
  coercive circular normals after the one complex elliptic soft normal.
- In the offset-one reflected sector, Fourier selection leaves only the
  bottom normal.  Its exact cross is `-8(5L-1)/L`; completing against
  L65 leaves the strict residual
  `32(L-1)(2L²-L+3)/(L(L²+36L-13))`.
- For grade two, the only character-allowed normal cross vanishes in
  every size.  Separated endpoint derivatives vanish individually; the
  unique short collision has upper/lower ratio `4:1` and cancels in the
  condition number.
- These are leading weighted faces, not yet a uniform all-grade
  circular-normal lift.  Raw single-grade division is the wrong
  coordinate because L159-old found nonlinear grade aliasing.
  `proof/crabb_circular_normal_face.md`;
  `proof/crabb_grade_two_normal_selection.md`.

## NEWEST (2026-07-23): A106 falsifies the two-sided raw Rees comparison
- Ordinary invertibility of the characteristic coefficient map does
  not survive the singular anisotropic pullback.
- For `L=6` with the sole raw coefficient `z_1=t`, exact elimination
  gives
  `a_3=32t^7/((6t²−1)(16t⁴−14t²+1))`.  With `c=t^5`, the proposed raw
  norm has order `t^52`, while the prepared norm contains a grade-three
  term of order `1024t^44`; the ratio diverges.
- This nonlinear high-to-low grade alias is favorable for L158: it
  adds negative prepared descent, so the complete disk-flat sign still
  pulls back physically.  What fails is only the attempted two-sided
  quantitative comparison.
- For the strong-variable merger, either prove the one-sided lower
  estimate `N_prep>=c_LN_raw` or retain `N_prep` directly.  Do not use
  an ordinary inverse-function theorem across `c=0`.
  `proof/crabb_disk_flat_marked_merger.md`; exact coefficient model.

## NEWEST (2026-07-23): L158 closes the mixed disk/reflection model gap
- In L149's independent real reflected-Rees polydisk, the model gap
  `Delta=U-R` is genuinely nonnegative, not merely a formal series.
  L157 gives its zero-reflection value `O(Q²)`.
- A uniform one-dimensional Taylor test for a nonnegative function
  gives `||D_eta Delta(z,0)||²<=C Delta(z,0)`.  Hence the entire mixed
  linear reflected row is `O(Q)` without another endpoint recurrence.
- L145--L146 make the complete compact reflected Hessian of `Delta`
  zero at the apex.  Continuity and Taylor expansion therefore give
  `Delta<=CQ²+CQ||eta||+epsilon||eta||²`.
- On the dual side, L152 supplies `-aQ`, L154 controls its linear
  reflected row by `CQ||eta||`, and L144 supplies a uniformly negative
  reflected Hessian.  Shrinking absorbs the model gap and proves
  `U-4<=-a_1Q-b_1(|r|²+sum|w_k|²)`.
- Thus the disk-normal/reflected tube is complete in the independent
  marked chart without exact all-disk complementarity and without a
  mixed coefficient induction.  Its sign survives physical pullback.
  Remaining: absorb L118's strong variables, retaining the prepared
  norm or proving only the one-sided raw comparison.
  `proof/crabb_disk_flat_marked_merger.md`.

## NEWEST (2026-07-23): L157 pushes the model gap above the disk face
- Let `T=q_can*adj(xi I-A)e_L` and `D=g#`.  L156 plus L155 gives
  coefficientwise `T-D=O(Q)`.
- If `E D=1 mod chi`, Laurent moment extraction writes every failure
  of `q_can` to annihilate the model orbit
  `D(A)^(-1)span{e_L,...,A^(L-1)e_L}` as a finite linear functional
  of `T-D`.  The observability system is uniformly invertible, so the
  normalized model defect obeys `d_model-q_can=O(Q)`.
- For the model Stein Gramian, the lower Rayleigh residual on `e_0`
  is exactly `(d_model-q_can)/2`.  The equality spectrum
  `{1/2,1,...,1,2}` has a uniform lower gap, so Temple's bound squares
  the error: `lambda_min=1/2-O(Q²)`.  Complementarity fixes
  `lambda_max=||B(A)||_K²/2`.
- Therefore the nonnegative failure of exact model complementarity
  satisfies
  `0<=Delta(z,0)<=C_LQ(z)^2`.  It is strictly above the quartic
  disk-normal face and cannot spoil L152.  Full real coordinate
  first jets of the normalized defect are exact through length ten.
  `proof/crabb_disk_flat_marked_merger.md`;
  `experiments/crabb_disk_one_reflection_jets.py`.

## NEWEST (2026-07-23): L156 closes the all-size endpoint first jet
- At a phase-one equality point, an arbitrary Hermitian Toeplitz disk
  tangent `h` satisfies an inverse-free three-piece formula for
  `K dotA`: a direct top row, a reverse-conjugate bottom row, and one
  rank-one interior correction.
- Two explicit companion endpoint vectors collapse that tangent to
  the scalar identity
  `d(K dotA)c=xi(P_h g-2g# P_h^flat)`.  At the roots of `g`, rank-one
  adjugate factorization then gives
  `D_h det(xi I-A)=2xi P_h^flat`.  Polynomial continuation removes
  the generic simple-root assumption.
- Direct resolvent differentiation consequently gives
  `D_h(q*adj(xi I-A)e_L)=2P_h` and proves that
  `q*(xi I-A)^(-1)e_L-g#/(xi g)` has zero full real disk-coordinate
  first jet on every phase-palindromic equality branch.
- L149 transfers this identity to every marked one-reflection Hardy
  coefficient.  L155's cone division now proves L154:
  `|one-leg_k|<=C_L|c|^kQ` uniformly through the Crabb apex.  Thus the
  dual one-leg part of the disk-flat merger is closed; the remaining
  gates are the nonnegative primal/model gap, its mixed-face
  polarization, the triangular marked-coordinate estimate, and then
  L118's strong variables.
- The checker independently verifies all four inverse-free
  recurrences at generic complex anchors through length ten, in
  addition to its stronger full ambient Gram-gradient audit.
  `proof/crabb_disk_one_reflection.md`;
  `experiments/crabb_disk_one_reflection_jets.py`.

## NEWEST (2026-07-23): L155 proves uniform division by the disk quartic
- Let `E={z=omega J conjugate(z)}` and
  `Q=||z||^4-|z^T Jz|²`.  If a real-analytic germ `F` starts in
  ordinary degree four and both its value and full gradient vanish on
  `E`, then `|F|<=C_L Q` locally.
- The proof is short and uniform through the singular apex.  L124's
  best-phase split gives `z=u+v`, `Q=4||u||²||v||²`, with
  `||v||<=||u||`.  Fourth-order apex vanishing implies
  `D²F(u)=O(||u||²)` and `D³F(u+theta v)=O(||u||)`.
  Since `F(u)=DF(u)=0`, normal Taylor expansion is therefore bounded
  by `C||u||²||v||²`.
- Literal scalar divisibility by `Q` is not required; the quartic
  face can be any quadratic form in the determinantal minors.
- Applied grade-by-grade, L155 says a one-reflection coefficient with
  these jets is `O(Q)`.  Its physical factor `c^k` is then absorbed by
  L152's `-a_L Q`.  L156 now supplies the full-normal stationarity;
  the strengthened symbolic-power argument makes a separate
  fourth-order apex hypothesis unnecessary.  A shared exact
  truncated-series engine verifies the conclusion through length ten
  without root tracking.
  `proof/crabb_quartic_cone_division.md`;
  `experiments/crabb_disk_one_reflection_jets.py`.

## NEWEST (2026-07-23): A101 endpoint-determinant route (closed by L156)
- The root-free checker now evaluates L149's endpoint-resolvent
  residual
  `q_z*(xi I-A_z)^(-1)e_L-g_z#/(xi g_z)` by a cleared determinant.
  At one generic complex phase-palindromic anchor in every length
  three through ten, its value vanishes and its full Gram-block
  gradient is supported only on the two endpoint diagonal directions.
  These directions are absent from the normalized Toeplitz disk chart.
- The stronger experimental statement is exact: the checker forms
  the complete gradient matrix, not one sampled tangent, and evaluates
  it at `L+1` rational resolvent points.  Since the cleared numerator
  has degree at most `L`, this proves the polynomial first-jet identity
  at each recorded anchor.
- Clearing the common denominator turns the residual into the single
  bordered-determinant identity
  `e0*K adj(xi K-2HR) K eL - det(K-2xi R*H)`.
  L156 proves the required Toeplitz-coordinate gradient by a shorter
  three-boundary-piece recurrence.  The stronger ambient Hermitian
  cofactor-gradient formula remains a finite-anchor observation, but
  it is no longer a proof gate.
- L149's logarithmic-inner proof therefore remains valid modulo the
  square of the equality-normal ideal, and L155 converts the jet
  statement to the required `O(Q)` bound.
  `proof/crabb_disk_one_reflection.md`;
  `experiments/crabb_disk_one_reflection_jets.py`.

## NEWEST (2026-07-23): L155 no longer needs a separate apex jet
- First-order vanishing on the complete phase-palindromic equality
  cone automatically removes Taylor degrees zero through three.
  After a real-orthogonal change, the cone is the rank-one locus of a
  generic `2 x (L-1)` matrix.  Value-plus-gradient vanishing puts each
  homogeneous Taylor term in the second symbolic power of its
  maximal-minor ideal; this equals the ordinary square and begins in
  degree four.
- Therefore L155 now needs only `F|E=0` and `DF|E=0`; the formerly
  separate assumption `F=O(||z||^4)` follows.  For the one-reflection
  merger, L149 supplies the value and L156 now supplies the gradient.
  The exact apex records remain an independent audit rather than a
  separate proof obligation.
  `proof/crabb_quartic_cone_division.md`.

## NEWEST (2026-07-23): A101 finds the general-disk one-reflection obstruction
- L149's one-reflection stationarity is special to the exact
  phase-palindromic equality cone.  It does not extend identically
  over the whole Toeplitz disk chart.
- The first tangent can be evaluated exactly without roots:
  `T'=A^dagger_K-A³`, while monic preparation gives
  `n1=rem_g[w³g'-(g'-g'(0))/w+a_(L-1)w^(L-1)]`.
  Rational Fréchet differentiation then gives the simple top
  Blaschke-norm derivative.
- At `(z1,z2)=(1/20,1/30)` this derivative is the explicit nonzero
  rational number recorded in `proof/crabb_disk_one_reflection.md`,
  about `-8.48e-6`.  Thus the second overly strong shortcut in the
  first L153 draft is also false.
- The corrected target survives exact adversarial tests.  Around
  rational equality anchors in lengths three through five,
  `partial_c R(u+epsilon v,0)/Q(u+epsilon v)` remains bounded and
  converges as `epsilon->0`; the odd normal part starts cubically.
  The viable lemma is therefore `|partial_c R(z,0)|<=C_L Q(z)`, not
  exact zero.  Such a term is absorbed by L152 after multiplying by
  its physical reflected grade `c^k`.
- L156 subsequently proved the required full first jet, and L155
  showed that it automatically includes the fourth-order apex start.
  Thus the corrected `O(Q)` target is now L154.
  `proof/crabb_disk_one_reflection.md`;
  `experiments/crabb_disk_one_reflection.py`.

## NEWEST (2026-07-23): A100 falsifies exact all-disk model complementarity
- The first L153 merger draft incorrectly promoted L145's one-line
  complementary relation to equality of the full condition number and
  characteristic Blaschke norm at every Toeplitz disk point.  Rank one
  only makes the top singular line a generalized eigenline; equality
  also requires the left singular line to be the bottom eigenline.
- An exact rational `4 x 4` counterexample with Toeplitz coefficients
  `(1/20,1/30)` has rank-one `B(A)` and exact model-kernel annihilation,
  but the proposed lower endpoint has a nonzero rational generalized
  eigenvector residual.  Therefore the orbit-complement condition is
  strictly larger than `||B(A)||_K²`.
- The old binary64 grid hid the failure because it sampled amplitudes
  below `.07`, where the gap is extremely high order and at most about
  `1e-7`.  The corrected checker starts with the exact counterexample
  and then records visible positive gaps on general complex samples at
  amplitudes `.12,...,.22`.
- This does not falsify the disk-flat theorem.  It removes a shortcut.
  L154 now controls the dual one-reflection sector and L157 proves
  that the nonnegative model gap is `O(Q²)` on the disk face.  L158
  closes the mixed gap by positivity.  A106 subsequently disproves a
  stronger two-sided raw/prepared norm comparison, which is not
  needed for the disk-flat sign.
  `proof/crabb_disk_flat_marked_merger.md`;
  `experiments/crabb_disk_model_complement.py`.

## NEWEST (2026-07-23): L152 proves the uniform disk-normal tube
- L122's canonical rank-one disk Stein metric now satisfies the
  uniform local inequality
  `kappa_K(M)-4 <= -a_L Q(z)`, not merely the directional expansion
  `-32Q+O(||z||^5)`.
- The two missing structural inputs are exact.  L123 gives value four
  on every phase-palindromic equality branch.  At each such point, a
  fixed characteristic Blaschke lower bound and the canonical Stein
  upper bound touch at four; the disk spectral-set theorem makes the
  lower bound locally maximal, so the canonical metric has zero first
  derivative in every disk-normal direction.
- Use L124's best-phase split `z=u+v`, for which
  `Q=4||u||²||v||²`.  If `v` is comparable to `u`, L122's apex
  quartic dominates its fifth-order remainder.  If `v<<u`, Taylor
  expansion normal to the exact equality branch starts with
  `-128||u||²||v||²`; exact branch vanishing and stationarity force
  every remainder to retain enough `u` and `v` factors.  Compactness
  of the phase/equality sphere makes both estimates uniform through
  the singular apex.
- Thus the pure disk-normal anchor is closed in every fixed size.
  Together with L151's weighted raw elliptic Hessian and L117's axis,
  the remaining single-Crabb gate is only the marked mixed remainder
  and then L118's already-coercive strong variables.
  `proof/crabb_disk_normal_tube.md`.

## NEWEST (2026-07-23): L151 proves the raw disk-flat elliptic face
- The phase-palindromic companion pencil is invalid for a general
  Toeplitz disk coefficient.  The new checker instead expands the full
  coefficient gauge
  `S(a,c)=2K(a)^(-1)(H(a)R+cR*H(a))` before applying the ellipse map.
- The all-size raw Faber endpoint lemma identifies the correct reflected
  grade `k=L-j`.  For coefficient
  phase `zeta`, the disk characteristic derivative is
  `2 conjugate(zeta) xi^(L-j+1)`.  If `Y_j` is the full ellipse-pencil
  tangent and `P_m` are the Dickson polynomials, then
  ```
  e_0^*(DP_L[Y_j]+2 conjugate(zeta)P_(L-j))
      =4 conjugate(zeta)e_(L-j)^*,
  e_L^*(DP_L[Y_j]+2 conjugate(zeta)P_(L-j))
      =4 conjugate(zeta)c^(L-j)e_j^*.
  ```
  An explicit two-path recurrence proves this polynomially for every
  `L,j`.
- The sharp raw scalar factor is
  `G=P_L+2a(1+c^j)P_(L-j)`.  Its prepared Blaschke tangent telescopes
  along one translated shift; the `c^j` copy cancels every premature
  fold.  The linear top column and first singular coupling vanish
  through `c^(2k)`, while
  `[c^(2k)](B_2)_(0L)=-16`.  Hence the dual square loses exactly
  `64a²c^(2k)`.
- L145's orbit-complement defect extends to the raw gauge.  Its
  triangular orthogonality recurrence makes the feasible Stein
  condition square agree with the dual norm through `a²c^(2k)`.
  Distinct reflected grades are orthogonal on their first face; the
  PSD zero-diagonal argument from L146 closes every mixed primal
  coefficient.
- Therefore the complete optimized disk-flat amplitude Hessian obeys
  ```
  Q_L(z;c) = -64 sum_(j=1)^(L-1) |z_j|² c^(2(L-j))
             + terms strictly above the diagonal face.
  ```
  After scaling by `diag(c^(L-j))` it tends to `-64I`, so the raw
  Hessian is uniformly negative for each fixed size and small
  nonzero `c`.
- Three independent exact audits cover the full-gauge optimized
  Hessian, the real/imaginary endpoint identity through length 14,
  and the prepared dual/model-complement face on every offset through
  length five.  The older formal and complementarity grids regenerate
  byte-for-byte after the shared-engine refactor.
  `proof/crabb_disk_flat_elliptic_face.md`;
  `experiments/crabb_disk_flat_elliptic_face.py`;
  `experiments/crabb_raw_faber_endpoint.py`;
  `experiments/crabb_raw_faber_blaschke.py`.

## NEWEST (2026-07-23): A99 finds the corrected Faber--Blaschke dual square
- For a noncentral grade `k`, define
  `G=P_L+2a(P_k+P_(L-k))+lambda*a*c^k*P_(L-k)`, map its
  roots from the ellipse into the disk, and use them as the zeros of
  a degree-`L` Blaschke product.
- This inner function specializes **exactly** to L123's
  characteristic Blaschke product on `c=0` and to L116's
  Chebyshev--Blaschke product on `a=0`.  A fixed Chebyshev product
  fails the first test and has a spurious unweighted `a^2` loss.
- Coprime exact-anchor probes expose the associated scalar square
  `(4-||B(T)||^2)/(16a^2c^(2k))
  ->4+|lambda-2|^2/4`.
  Thus the natural root interpolation (`lambda=0`) gives coefficient
  80, while the single explicit correction `lambda=2` gives the sharp
  coefficient 64.
- Equivalently, the corrected Faber characteristic factor is
  `P_L+2aP_k+2a(1+c^k)P_(L-k)`.  Unrestricted numerical optimization
  of all mixed Faber coefficients independently returns
  `lambda_(L-k)=1.99...`; other coefficients affect only higher
  finite-`c` orders.
- L140 proves that `lambda=2` is the unique Hardy reflection of
  L131's negative frequency, so the center is no longer numerical.
- L141 also proves the entire mismatch curvature.  With
  `epsilon=(lambda-2)ac^k`, the associated Blaschke perturbation obeys
  the exact nilpotent identity
  `B_epsilon(C)=(1-|epsilon|^2)C^L+epsilon*C^(L-k)`.
  Its only top singular coupling is a two-column block, whose largest
  eigenvalue is `4-4|epsilon|^2+O(|epsilon|^4)`.  Hence the
  `4|lambda-2|^2a^2c^(2k)` part of the observed parabola is exact.
- At that stage the remaining scalar content was only the corrected base
  loss `4-||B_(lambda=2)(T)||^2=64a^2c^(2k)+o(...)`:
  derive it by inner--outer/endpoint singular-value perturbation and
  establish complementary slackness with L139's defect correction.
- L142 then proves that corrected base loss on **every central
  collision** `L=2k`.  The Dickson identity factors the corrected
  polynomial as `Q_(a,c^k)(P_k)`, so its root inner is a degree-two
  corrected inner composed with the degree-`k` Chebyshev--Blaschke
  map.  L126/L129 leave one active physical size-three block.  An
  exact polynomial-ring generalized singular calculation has no
  mixed term below `a^2c^(2k)` and gives coefficient `-64` there.
- This narrowed the open scalar step to proving that a
  **noncentral** grade has the same associated active block as its
  central `2k+1` model.  This must be a filtered localization, not an
  all-`c` norm comparison (A97 already falsified the latter).
- L143 closes that noncentral step by formal inner--outer preparation.
  The exact first mixed operator coefficient is a sparse translated
  shift, and the top singular-vector coupling begins strictly above
  grade `k`, so its Schur square cannot reach `c^(2k)`.
  The sole terminal fold is at `d=L-k`: its `+8c^d` numerator
  contribution cancels the coordinate metric's `-8c^d`, including
  when `L=3k` puts it on the target face.  The remaining top Gram
  coefficient is `-32`; the terminal coordinate weight `1/2` gives
  `-64`.
- Hence L140--L143 prove the complete one-grade dual loss parabola.
  L144 then polarizes it over every complex grade.  Different grades
  occupy distinct shift/Fourier modes; the only possible mixed folds
  are `L=2k+ell` and `L=k+2ell`, and both cancel by the same
  `+8` numerator / `-8` coordinate-metric mechanism.
- The complete dual principal face is therefore
  `-64 sum_k |u_k|^2 c^(2k)`.  The next gates are a uniform remainder
  and complementary slackness with L139's primal defect square.
- L145 supplies that complement canonically.  If `B=N/D` and `x` is
  its top right singular vector, choose the Stein defect orthogonal to
  `D(T)^(-1)span{x,Tx,...,T^(L-1)x}`.  The model-kernel identity makes
  the primal metric and dual inner exactly complementary on `x`;
  their one-grade Hessians agree through the face.
- L146 then closes every mixed **primal** coefficient without another
  coordinate calculation.  The optimized primal Hessian minus the
  dual Hessian is PSD.  After grade scaling, L145 makes its diagonal
  zero, so positivity forces the entire matrix to vanish.
- Therefore A98/A99's complete primal/dual Newton face is proved.
  L131/L120 also kill the entire linear equality gradient exactly.
- **Audit correction:** the first L147 draft promoted L143's
  amplitude-degree-two recurrence to an arbitrary-amplitude support
  theorem without writing the required recentered induction.  The
  uniform equality tube is therefore still conditional.  Its precise
  remaining gate is that Weierstrass preparation and the full Stein
  critical recurrence preserve the square of the reflected ideal
  generated by `r=c^L` and `w_k=c^ku_k` after recentering at L123's
  exact equality metric.  If this holds, ordinary convergent Newton
  domination closes the tube immediately.
- Binary64 continuation and higher-amplitude fits on one- and
  mixed-grade equality branches through `L=7` found no violation and
  put visible cubic/quartic terms strictly above the proposed face.
  The persisted 180-record probe has minimum half-face domination
  ratio `1.8577`.  This is adversarial evidence, not the missing proof.
- L148 removes the grade-one all-amplitude obstruction.
  At every nonlinear L123 equality anchor, not only at the Crabb
  apex, the first ellipse-pullback derivative of the sharp Blaschke
  norm and touching Stein envelope is exactly zero.  The proof splits
  inner-function stationarity from the all-size companion endpoint
  identity for `A^dagger_K-A^3`.
- L149 closes the all-grade extension correctly.  L131's exact Faber
  boundary identity supplies a convergent marked algebra in the
  reflected variables `r=c^L`, `w_k=c^ku_k`.  The companion resolvent
  identity `q*(zI-A)^(-1)e_L=1/(zB(z))` turns every one-reflection top
  norm derivative into the real mean of the logarithmic tangent of a
  fixed-degree inner function, hence zero for every grade, phase, and
  central fold.
- L150 supplies the primal lift without assuming that L118's optimizer
  preserves the marking.  L145's canonical orbit-complement defect is
  analytic at arbitrary amplitude.  Its condition gap above the
  Blaschke norm is nonnegative and zero on the equality divisor, so it
  belongs to the reflected ideal square automatically.  L145--L146's
  PSD zero-diagonal argument gives the same compact face.
- Therefore L147 is proved: the explicit model upper certificate is
  at most `4-8|c^L|²-32sum|c^ku_k|²` in a uniform equality-stratum
  neighbourhood.  The Faber/Blaschke remainder is no longer open.
- A new exact finite-amplitude dual checker substitutes rational
  equality coefficients before expanding in `c`, so it retains every
  amplitude order.  In seven cases through dimension nine, with
  minimum grades one through three and mixed higher grades, every
  coefficient below `c^(2k)` vanishes exactly and the `c^(2k)`
  coefficient is negative.  This supports the all-grade
  Hardy/model-space stationarity mechanism on the dual side but does
  not supply the primal reflected-Rees lift.
  `proof/crabb_faber_blaschke_dual.md`;
  `experiments/crabb_faber_blaschke_dual.py`;
  `experiments/crabb_central_faber_blaschke.py`;
  `experiments/crabb_faber_blaschke_formal.py`;
  `experiments/crabb_faber_blaschke_mixed.py`;
  `experiments/crabb_faber_blaschke_complementarity.py`;
  `proof/crabb_uniform_weighted_remainder.md`.

## NEWEST (2026-07-23): L139 proves the all-size defect-Hessian LDL edge
- L138's endpoint Schur formula depends only on the zeroth and last
  rows of `B=D^(-1)U diag(s)U^*D`.  DCT-I endpoint
  product-to-sum therefore makes the pure-defect Hessian exactly
  diagonal by scalar cosine mode.
- The corresponding modal weights have leading values
  `4c^(-j)` for interior modes and `(8/3)c^(-L)` for the terminal
  mode.  These are the all-size disk constants, not a finite fit.
- Multiplication by `dn/k'` transfers coefficient coordinate `i` to
  mode `j=i+2r` with associated coefficient `2c^r`: the two is one
  half of Jacobi's Fourier coefficient four.  Only even Fourier
  modes occur, so the two parity blocks separate exactly.
- Hence the pure-defect quadratic matrix has
  `H_(i,j)=8c^r+...` when `j=i+2r<L` and
  `H_(i,L)=(16/3)c^r+...`.  The exact LDL recursion gives interior
  lower coefficient `2c^r`, terminal coefficient `(4/3)c^r`, and
  diagonal constants `4,...,4,8/3`.
- This promotes A98's formerly finite Hessian clue to an all-size
  theorem and proves the outer factor
  `(I+cS^2)/(I-cS^2)`.  The remaining gate is now solely to insert
  the operator/Faber linear and constant terms in the same modal
  endpoint formula and prove the completed-square orientation.
  `proof/crabb_defect_hessian_edge.md`;
  `experiments/crabb_defect_hessian_factor.py`.

## NEWEST (2026-07-23): L138 restores the homogeneous defect coordinate exactly
- For an arbitrary rank-one forcing `q`, put `y=R^*q`, divide its
  spectral-node values by L117's axis defect `beta`, and call the
  resulting multiplier `h`.  The complete physical Stein metric is
  then exactly
  `K_0^(-1/2)M(q)K_0^(-1/2)=D^(-1)A_h W A_h^*D^(-1)`,
  where `A_h=U diag(h)U^*`.
- Thus defect optimization is multiplication by a sampled scalar
  function in the DCT-I node basis.  The missing scale direction is
  `h=constant`; it is no longer hidden by the exact-series gauge
  `x_0=0`.
- Linearizing `h=1+as` gives the exact congruence
  `P(a)=(I+aB)P_0(I+aB)^*`.  The equal squared DCT endpoint rows make
  the condition number stationary in every homogeneous defect
  direction.  Its full quadratic term is now two explicit endpoint
  Schur sums in the Toeplitz-plus-Hankel matrix
  `H_s=U diag(s)U^*`.
- This proves A98's structural normal form.  L139 subsequently derives
  the all-size Newton-edge LDL factor and terminal fold; the
  linear/constant coupling to L131's Faber row remains.
- A direct Stein-solve regression over sizes 4, 6, and 9 and three
  ellipse parameters confirms the normal form and the endpoint
  Hessian formula.
  `proof/crabb_homogeneous_defect_normal_form.md`;
  `experiments/crabb_homogeneous_defect_normal_form.py`.

## NEWEST (2026-07-23): L137 proves the axis-defect half of the square
- DLMF's reciprocal-`dn` Fourier series applied at L117's Lobatto
  nodes gives only even DCT-I modes.  Conjugation back through L135's
  eigenvector matrix yields, for `1<=r<L/2`,
  `d_(2r)=4(-1)^r c^r+O(c^(r+2))`; odd defect coordinates vanish.
- Hence the all-size associated spatial series is exactly
  `d_edge=(I-3cS^2)/(I+cS^2)`.  This promotes one of A98's two fitted
  factors to a theorem.
- Combining it with L139's subsequently proved Hessian LDL edge
  `(I+cS^2)/(I-cS^2)` gives the whitened transported tail
  `2-4sum_(r>=1)c^rS^(2r)`, explaining the magnitude-four collision
  with L131's reflected Faber row.
- L137 alone does not prove the LDL factor; L139 now supplies that
  factor.  The opposite normal orientation, Faber coupling, and
  uniform remainder remain open.
  `proof/crabb_axis_defect_fourier.md`;
  `experiments/crabb_axis_defect_fourier.py`.

## NEWEST (2026-07-23): A98 sharpens the missing Schur step to one square
- The exact Hessian frontier now includes the first new coprime
  noncentral grade-four pair `(L,k)=(9,4)`: all coefficients below
  `c^8` vanish and the leading coefficient is exactly `-64`.
  Its endpoint metric pair is `(48,128)`, again giving the
  gauge-invariant combination `128-4(48)=-64`.
- The optimized defect agrees below `c^4` with L131's natural
  transported defect `2U(u)d_c`.  At the reflected grade their
  difference is exactly the one-coordinate correction `-8c^4e_8`.
  The disk defect Hessian charges that correction by
  `4*8^2=256`; hence the transported certificate has face `+192`
  while optimization changes it to `-64`.
- The earlier `(L,k)=(7,3)` record has the identical completed-square
  constants, with correction `+8c^3e_6`.
- The raw correction is not itself universal: exact low-order scans
  retain the grade-three `+8` correction through lengths eight to ten
  but find no raw grade-four correction from lengths ten to twelve.
  Terminal transport can absorb the square center, so the theorem
  must be stated in the homogeneous defect quotient.
- This identifies a precise conjectured associated-graded normal
  form.  For L131's reflected row
  `r=4sum_j u_jc^je_(L-j)*`, restore homogeneous defect scale and
  seek an isometry `J_def` such that
  `Q_face=4||eta-2J_def r||^2-4||r||^2`.
  Minimization gives the required
  `-64sum_j|u_j|^2c^(2j)`, and Faber-row orthogonality kills mixed
  grades at the same time.
- L139 proves that the universal pure-defect Hessian has interior
  lower factor
  `I+2sum_(r>=1)c^r shift^(2r)=(I+c shift^2)/(I-c shift^2)`,
  with terminal coefficient `4/3` and diagonal endpoint weight
  `8/3`.  This is the all-size analytic outer factor of the Newton
  edge of L117's Szegő weight.
- The coefficient arithmetic is now explained.  L137's axis-defect
  edge is `(I-3cS^2)/(I+cS^2)`.  Multiplying the transported defect
  by the LDL outer factor gives
  `2(I-3cS^2)/(I-cS^2)=2-4sum_(r>=1)c^rS^(2r)`.
  Its first reflected tail has magnitude four, opposite L131's Faber
  row of magnitude four; the gap eight costs `4*8^2=256`, while the
  negative row energy is `-4*4^2=-64`.  The constants and pure-defect
  factor are no longer fitted; the open work is the Faber/reversal
  coupling.
- The completed square is still a target, not a lemma.  The remaining
  proof must put L131's operator and coefficient derivatives into
  L138's modal endpoint coordinates, identify their signs, and show
  every omitted terminal fold has strictly higher weight.
  `proof/crabb_principal_face_completed_square.md`;
  `experiments/crabb_principal_face_locality.py`.

## NEWEST (2026-07-23): L136 reduces every one-grade problem to a coprime pair
- Write `L=dq`, `k=ds`.  Degree-`d` Dickson descent preserves
  residues modulo `d` and reduces the residue-zero polynomial pencil
  exactly to the complete size-`q+1` grade-`s` equality/ellipse family
  at parameter `c^d`.
- The Toeplitz coordinate Gramian reduces the same space, every
  subcritical Dickson compression vanishes, and the polynomial is
  exactly affine in the amplitude.
- L132's outer-critical-factor defect reconstruction therefore applies
  unchanged:
  `t_*(T_(q,s)(a,c^d)) <= t_*(T_(L,k)(a,c))
  <= Gamma_(q+1)(T_(q,s)(a,c^d))`.
- With `d=gcd(L,k)`, the reduced pair `(q,s)` is coprime.  Hence the
  unresolved nondivisor diagonal face can now assume
  `gcd(L,k)=1`; all noncoprime cases lift from a smaller coprime one.
- This is an exact rank-one sandwich, not equality of the two full
  similarity optima in general.  It narrows but does not close the
  unequal-residue gate.
- `proof/crabb_gcd_dickson_descent.md`;
  `experiments/crabb_gcd_dickson_descent.py`.

## NEWEST (2026-07-23): A95 suggests exact grade diagonality
- The distinct-grade cancellation is much deeper than the first
  Newton face: the exact size-five `(1,2)` cross jet is zero through
  `c^17`, and all three exact size-seven cross jets are zero through
  `c^11`, including same-parity `(1,3)`.
- Finite-amplitude probes at `c=.12,.25,.4` also put the mixed Hessian
  at optimizer tolerance.  The residual falls by a factor four when
  the amplitude step is halved, so it is quartic finite-difference
  contamination rather than a quadratic cross term.
- New stronger conjecture: L118's entire amplitude Hessian on a fixed
  phase-palindromic equality/ellipse branch is block diagonal in
  distinct Faber grades for each fixed `c`.
- Mixed real/imaginary phase probes support the same cross-grade
  cancellation.  However the real and imaginary diagonal values
  differ at fixed nonzero `c`, so the conjecture is grade-block
  diagonality, not all-`c` phase isotropy.  L134 supplies the needed
  phase isotropy only on the first face.
- The likely proof basis is L117's DCT-I diagonalization of the exact
  Toeplitz-plus-Hankel elliptic Szegő kernel.  The remaining hard step
  is to transform the optimized defect quadratic itself; Faber
  endpoint Parseval does not do this.
- L135 now proves the exact coordinate bridge.  If
  `R=K_0^(−1/2)DU` is the DCT-I right eigenvector matrix of the axis,
  then the grade-`j` companion row satisfies
  `(v_j^*R)_n=-2sqrt(2/L)eps_n c^((j+1)/2)
  sin(j theta_n)sin(theta_n)`.  Hence the coefficient grades are
  genuine DST-I modes in this spectral leg, with no terminal alias.
- The remaining A95 identity is correspondingly sharp: conjugate
  L118's defect Schur kernel into these coordinates, remove the common
  sine factors, and prove the kernel is Toeplitz-plus-Hankel (DST-I
  diagonal).  Its grade eigenvalues must start at `-64c^(2j)`.
- A95 remains conjectural despite L135.  The weaker all-size
  first-face statement in A94 remains the theorem-level gate.

## NEWEST (2026-07-23): A94 isolates the distinct-grade proof gate
- Exact polarization of L118's optimized amplitude Hessian for two
  different real phase-palindromic grades gives zero through and
  including the first allowed mixed weight `c^(k+l)` for every pair
  through size seven.
- Selected deeper runs find the complete recorded cross jet zero
  through order 11 (size five) and order nine (sizes six/seven).
  The same-parity pair `(1,3)` also vanishes, so this is not merely a
  parity selection rule.
- This is not yet L135.  The precise target is
  `B_(k,l)(c)=O(c^(k+l+1))` in every size.
- L125 puts the first possible interaction at weight `k+l`; L131
  makes the two associated Faber/reflected rows orthogonal; L65 makes
  the Crabb metric Hessian circle-mode diagonal.  The remaining proof
  is that L118's defect-variable Schur complement commutes with this
  associated-graded projection.  Any terminal fold/alias must be
  shown to gain strictly higher weight.
- `proof/crabb_mixed_grade_face.md`;
  `experiments/crabb_mixed_grade_face.py`.

## NEWEST (2026-07-23): L134 proves complex one-grade phase isotropy
- For the pure-imaginary first-offset coefficient
  `(u_1,u_(L-1))=(i,-i)`, exact rational real/imaginary-pair Stein
  arithmetic gives defect jet
  `i(-2e_1+2e_(L-1)+8ce_3)` and the same endpoint pair `(48,128)`.
  The optimized amplitude Hessian is again `-64c²+O(c³)`.
- The first-face coefficient is a real quadratic form in
  `u=x+iy`.  Entrywise conjugation maps `y` to `-y` while preserving
  the condition number, so the mixed `xy` term vanishes.  L133 and
  the imaginary calculation therefore prove the arbitrary-phase law
  `-64|u|²c²`.
- L132 complexifies without extra condition cost: Dickson paths and
  the Hermitian coordinate Gramian reduce the one-pair coefficients
  `(u,conjugate(u))`, while the critical-factor defect reconstruction
  is complex-linear in the outer defect.  Thus every divisor grade
  has rank-one upper face `-64|u|²c^(2k)`.
- Exact complex offset-one records cover sizes four through twelve;
  stationarity is independently reconstructed through size seven.
  A separate `u=i` Dickson/Gramian checker covers `k<=5`,
  `3<=q<=6`, through dimension 31.
- Remaining coefficient gates: cross terms between **distinct**
  grades and unequal-residue localization when `k` does not divide
  `L`.  Uniform Newton-remainder control remains after those.
  `proof/crabb_offset_one_complex_face.md`;
  `experiments/crabb_offset_one_complex_face.py`;
  `experiments/crabb_divisible_complex_descent.py`.

## NEWEST (2026-07-23): L133 proves the universal real offset-one face
- For the real phase-one first-offset equality direction in every size
  `L+1>=3`, L118's optimized rank-one amplitude Hessian is now proved
  to begin `-64c²+O(c³)`.
- The disk defect tangent is `2 sum_(j in {1,L-1})e_j`.  Expanding the
  exact defect stationarity equation one elliptic order further gives
  the dimension-independent correction `-8ce_3`, with the obvious
  omission/collisions on the shortest chains.
- At second amplitude order, the coefficient-`c²` Stein forcing
  stabilizes for `L>=7` at
  `(48,-96,-32,16,0,...,0,16,-32,64)`.
  The short lengths `2,...,6` have different collided diagonals but
  the same endpoint result.  Stein inversion gives lower/upper
  coordinate-metric coefficients `(48,128)`, and the generalized
  condition coefficient is `128-4(48)=-64`.
- The explicit defect-gradient residual is `O(c²)`.  L118's invertible
  analytic defect Hessian therefore shows that the true optimizer
  differs only at order `c²`, which cannot alter this coefficient.
- Combining L133 with L132 transfers the real phase-one upper face
  `-64a²c^(2k)` to every divisor grade `k|L`.
- L134 subsequently closes complex one-grade phase.  Mixed-grade
  polarization, unequal-residue localization for `k` not dividing
  `L`, and a dimension-uniform analytic remainder are still open.
  `proof/crabb_offset_one_face.md`;
  `experiments/crabb_offset_one_face.py`.

## NEWEST (2026-07-23): L132 extends central descent to every divisor grade
- If `L=qk`, the residue-zero space
  `V=span{e0,e_k,...,e_(qk)}` exactly reduces the degree-`k`
  Dickson image of the phase-one equality/ellipse pencil.
  Its block is the full size-`q+1` first-offset family at parameter
  `r=c^k`; all polynomial cross blocks vanish and the result is exactly
  affine in amplitude.
- The Toeplitz coordinate Gramian reduces the same space and restricts
  to the exact size-`q+1` Gramian.  Every subcritical Dickson grade
  compresses to zero on `V`.  The other residue chains are unweighted
  nilpotent shifts at the apex and are strictly inactive locally.
- L130's outside-critical-factor trace works for outer dimension
  `q+1`: it lifts L118's outer rank-one defect to a full defect, makes
  the metric reduce `V`, preserves its outer condition exactly, and
  traps all inner levels between the active endpoints.
- The precise conclusion is the sandwich
  `t_*(T_(q,1)) <= t_*(T_(L,k)) <= Gamma_(q+1)(T_(q,1))`.
  Do **not** replace the right side by `t_*` for `q>2`; L118's
  rank-one branch has not been proved globally SDP-optimal off axis.
- The exact checker covers `k<=6,q<=6` (up to dimension 37).
  An independent full physical/critical-factor metric checker gives
  residuals below `8e-13` for degrees two/three, quotients three/four,
  and two ellipse parameters, with strict inner endpoint margins.
  `proof/crabb_divisible_dickson_descent.md`;
  `experiments/crabb_divisible_dickson_descent.py`;
  `experiments/crabb_divisible_metric_lift.py`.

## NEWEST (2026-07-23): L131 identifies the complete reflected Hardy vector
- For `S0=C+cJCJ`, the physical-adjoint equality tangent attached to
  coefficient grade `j` is
  `H_j=E_j+cJE_(L-j)J=2(e0-eL)(e_(j+1)*-c e_(j-1)*)`.
  Exact differentiated Dickson recurrence proves
  `e0*DP_L(S0)[H_j]=eL*DP_L(S0)[H_j]=0` in every size.
- The intermediate endpoint paths have an explicit folded formula:
  one enters from grade `j+1`, one from grade `j-1` with weight `-c`,
  and at Dickson grade `L` they become the identical
  `c^j e_(L-j)*` row and cancel.
- Therefore the Faber transform of
  `g_a(xi)=xi^L+2a sum_j u_jxi^j` has amplitude endpoint rows
  `4sum_j u_je_j*` and `4sum_j u_jc^je_(L-j)*`.
  The second row has exact energy
  `16sum_j|u_j|²c^(2j)`.  The scalar Joukowski/Fourier form gives the
  same Parseval identity and rules out mixed circle grades.
- This is the algebraic source of A84's candidate Hardy face,
  including complex phase-one directions.  It is not yet the metric
  theorem: show that all other Faber rows and the model-space lift are
  inactive/removable to strictly higher Newton weight.  L130 supplies
  each central one-coordinate normal and its strict inner gap.
  `proof/crabb_faber_reflection.md`;
  `experiments/crabb_faber_reflection.py`.

## NEWEST (2026-07-23): the exact Hessian checker now has a correct internal map guard
- Extending A85 to `p=8,k=3` exposed a checker defect: high scalar
  coefficients of the reverted elliptic map were being computed at the
  same truncation order as the requested matrix jet.  Their repeated
  valuation-one divisions polluted the terminal coefficients and
  produced a spurious huge rational instead of `-64`.
- The map generator now computes only the required scalar degrees at
  internal order `output+2*maximum_degree+4`, as required by L125, and
  truncates afterward.  It regenerates every saved `p=3,...,7` record.
  The formerly failing `p=8,k=3` jet is exactly zero through degree
  five and exactly `-64` at degree six; its endpoint metric
  coefficients are exactly `(-16,-128)`.
- This strengthens the finite audit but does not replace L131's
  all-size algebra or the remaining metric-normal-form proof.

## NEWEST (2026-07-23): L130 proves the exact central metric lift
- L127 now includes the companion positive dual transfer
  `R(Z)=sum_j f_j(T)Zf_j(T)*`, whose Stein difference is exactly
  `Z-B(T)ZB(T)*`.  The independent general-matrix audit verifies both
  primal and dual identities through dimensions three to seven and
  Blaschke degrees one to four.
- For real `B=N/D`, factor its critical Wronskian as
  `N'D−ND'=kappa Q Q#`, with `Q` carrying the outside critical points,
  and set `F=Q/D`.  A residue calculation proves
  `(1/k)sum_(B(t)=u)f(t)/F(t)=<f,F>/||F||²`, independently of the
  fiber value.
- The central pencil has the matching exact filtration:
  `Pi P_m(S(a,c))Pi=0` for all Dickson grades `1<=m<k`, including
  every amplitude power.  Polynomial division transfers the scalar
  fiber trace to the physical outer block.
- If `F=sum gamma_jf_j`, then `J=F(T)^−*O` maps the size-three defect
  to a full defect and gives `sum gamma_jw_j=0` algebraically.  The
  resulting full metric reduces the outer space and restricts to the
  exact size-three optimal metric.  Its apex inner levels are all two,
  strictly between endpoints one and four, so continuity and L129
  prove `t_*(T_k(a,c))=t_*(T_1(a,c^k))` locally for every fixed `k`.
- This closes the central `−64a²c^(2k)` event.  It does not yet localize
  a noncentral offset inside a larger block or polarize different
  coefficient grades.
  `proof/crabb_central_dual_lift.md`;
  `experiments/crabb_central_dual_lift.py`.

## NEWEST (2026-07-23): the central descent is exact in amplitude, and its metric target is now explicit
- L126 is stronger than its first version.  In size `2k+1`,
  `P_k(S+aE)=P_k(S)+aDP_k(S)[E]` **exactly**, not only modulo `a²`.
  The outer `span{e0,e_k,e_(2k)}` block is the complete size-three
  equality/ellipse family at parameter `r=c^k`.
- The central Toeplitz coordinate Gramian reduces the same space and
  restricts to the exact size-three Gramian.  Thus the proper-map image
  has an orthogonal physical size-three summand for all admissible
  amplitudes.
- A new 18-case grid (`k=2,3,4`, three `c` values, two amplitudes)
  finds
  `t_rank1(full)=t_rank1(size3)` within `1.1e-13`; unrestricted SDP
  values agree within solver tolerance.  More strongly, the optimized
  full metric has outer/inner cross residual below `2.3e-7`, its outer
  compression equals the size-three metric below `1.6e-7`, and every
  inner generalized level is strictly between the outer endpoints.
- L127 proves the exact finite-Blaschke Stein transfer
  `L(X)-T*L(X)T=X-B(T)*XB(T)` and Gramian composition
  `P_T(q)=L(P_(B(T))(q))`.  Hence upper feasibility really does lift;
  the only missing special fact is condition control.  Prove that this
  transfer reduces the outer space and traps all inner metric levels.
  That would give the exact identity
  `t_*(T_k(a,c))=t_*(T_1(a,c^k))`, not merely its Hessian.
  `proof/crabb_central_dickson_descent.md`;
  `proof/blaschke_stein_composition.md`;
  `experiments/crabb_central_metric_descent.py`.

## NEWEST (2026-07-23): L128 collapses every inactive central fiber to one involution
- Order the inner coordinates in the two layers
  `(e1,...,e_(k−1))` and `(e_(k+1),...,e_(2k−1))`.  With
  `H=J diag(c^(k−1),...,c)` and `r=c^k`, the full descended pencil is
  ```
  [[H,I],[rI,H]] + 2a[[-rI,H],[H,-I]],
  ```
  while its coordinate Gramian is `[[I,2aI],[2aI,I]]`.
- The decisive identity is `H²=rI`.  A diagonal gauge sends `H` to
  `sqrt(r)J`, so algebraically the whole inactive complement contains
  only two `2×2` fiber types.  The gauge is nonunitary in the physical
  Gramian, so this is not yet the missing condition estimate.
- The next calculation should formulate L127's lift as a `2×2`
  operator Schur complement with the sole relation `H²=rI`, and show
  that its generalized inner spectrum remains between the exact
  size-three endpoints.  No induction over individual coordinates is
  now justified.
  `proof/crabb_central_inner_fibers.md`.

## NEWEST (2026-07-23): L129 proves the inner fibers cannot become active
- At the apex L128's inner block is
  `N=[[0,I],[0,0]]`.  The fixed metric `diag(I/2,I)` has exact Stein
  slack `I/2` and condition square two.
- Since `||H||=c` and `r=c^k<=c²`, the full inner disk operator is an
  `O(|a|+c)` perturbation uniformly in `k`; its coordinate Gramian is
  also uniformly close to identity.  The fixed metric therefore stays
  feasible with condition below three.
- The outer size-three value stays above three near the apex.  Hence
  the similarity optimum of the descended direct sum is **exactly**
  its size-three outer value.
- This is not yet an upper bound for the original `T`: the remaining
  gate is precisely the condition-number cost in L127's positive
  Blaschke lift, at the `a²c^(2k)` endpoint order.
  `proof/crabb_central_inner_inactivity.md`.

## NEWEST (2026-07-23): L126 finds the all-size source of the `-64`
- In size `p=2k+1`, put `S=C+cJCJ` and take the central equality
  direction
  `E=2(e0-e_(2k))(e_(k+1)^*-c e_(k-1)^*)`.
  The degree-`k` Dickson polynomial and its derivative both reduce
  `span{e0,e_k,e_(2k)}` exactly.
- On that outer space they are
  ```
  P_k(S)       = [[0,2,0],[c^k,0,1],[0,2c^k,0]],
  DP_k(S)[E]   = [[-2c^k,0,2],[0,0,0],[2c^k,0,-2]].
  ```
  These are precisely the size-three ellipse and its central
  phase-palindromic equality tangent at descended parameter `r=c^k`.
  Every polynomial cross block to the inner fibers is zero.
- Dickson first-reversal path cancellation proves this in arbitrary
  `k`; the symbolic checker regenerates it through `k=12`.
- This explains why A85's first central Hessian term is the size-three
  event `-64(c^k)^2`.  Exact high-order algebra also finds
  `H_(5,2)(c)=H_(3,1)(c²) mod c^14`.
- **Remaining direction matters:** the inner-function identity gives a
  lower comparison automatically, but A84 needs an upper metric.
  Lift the active size-three rank-one metric through L117's strictly
  interior fiber levels and prove their Schur feedback is above
  `c^(2k)`.  Then localize a noncentral offset to its `2k+1` window
  and polarize grades.
  `proof/crabb_central_dickson_descent.md`;
  `experiments/crabb_central_dickson_descent.py`.

## NEWEST (2026-07-23): A85 makes the candidate `-64` face exact in finite sizes
- A new rational formal-series engine expands the complete locally
  optimized rank-one Stein calculation, not merely the scalar Riemann map:
  ```
  T=T0+aT1+a²T2,  M=M0+aM1+a²M2,  K=K0+aK1.
  ```
  It reconstructs the quadratic in the defect tangent, minimizes it by
  exact series Gaussian elimination, and perturbs both generalized metric
  endpoints.
- In every first-offset case through `p=7`, all coefficients below
  `c^(2k)` vanish exactly and the first coefficient is `-64`.  Thus the
  floating A84 law is now an exact finite recurrence, not a precision
  artifact.
- The optimizer itself exposes the cancellation: it starts at L123's
  equality defect and introduces a forced `-8c` entry two coordinates
  later.  For offset one the first endpoint metric coefficients are
  `(48,128)`; for tested offsets two and three they are `(-16,-128)`.
  In either case the condition ratio combines them as
  `M2_LL-4M2_00=-64`.
- **Do not promote this to L126.**  The remaining proof debt is an
  arbitrary-`k` induction for the two-step stationarity/path recurrence,
  followed by polarization of distinct coefficient grades and a uniform
  analytic remainder.  The exact finite checker identifies those tasks
  but does not replace them.
  `proof/crabb_palindromic_elliptic_hessian.md`;
  `experiments/crabb_palindromic_elliptic_hessian.py`;
  `experiments/exact_truncated_series.py`.

## NEWEST (2026-07-23): L125 proves the scalar grade-selection rule in A84
- If the normalized disk-to-ellipse map is
  `Psi_c(z)=sum_n b_n(c)z^(2n+1)`, then
  ```
  b_n(c)=c^n(Catalan_n+O(c²)).
  ```
  In particular its lowest bivariate edge is
  `z+cz³+2c²z⁵+5c³z⁷+14c⁴z⁹+...`.
- The proof differentiates the exact elliptic formula and uses
  ```
  (Psi')²(k-z²)(1-kz²)=alpha²(4c-Psi²).
  ```
  Theta identities make `k/c`, `alpha`, and the linear coefficient
  analytic in `c²`.  The Newton blow-up `Psi=zH(c,cz²)` is therefore a
  regular even-`c` recursion.
- On its lowest face the recursion becomes
  `(H+2yH')²(1-4y)=H²`; its normalized solution is the Catalan
  generating function `2/(1+sqrt(1-4y))`.
- This rigorously proves that an offset-`k` equality coefficient cannot
  enter through the scalar Riemann map below elliptic grade `c^k`.
  It does **not** yet prove A84's candidate `-64` optimized Stein
  coefficient.  Resume by linearizing the companion-coordinate
  operator and solving the defect-vector recurrence at this first
  permitted grade.
  `proof/ellipse_inverse_newton_edge.md`;
  `experiments/ellipse_inverse_newton_edge.py`.

## NEWEST (2026-07-23): L124 resolves the singular equality cone exactly
- Put `Cz=J conjugate(z)` on the Toeplitz coefficient space and
  `E={z:z=omega Cz}`.  Choosing the phase of `z^T Jz` gives the
  real-orthogonal splitting
  ```
  z=u+v,   T u=u,   T v=-v,   T=eta C.
  ```
- It obeys the exact identities
  ```
  ||u||²=(||z||²+|z^TJz|)/2,
  ||v||²=(||z||²-|z^TJz|)/2,
  dist(z,E)=||v||,
  Q(z)=4||u||²||v||².
  ```
  Thus, with `d=dist(z,E)`,
  `2||z||²d²<=Q(z)<=4||z||²d²`.
- This is the correct stratified replacement for a Morse--Bott chart:
  away from the apex L122 is quadratically coercive in the disk normal,
  with coefficient `||u||²`; at the apex that coefficient must vanish.
- The local merger should now use variables `(u,v,c)`:
  equality amplitude, non-palindromic disk normal, and elliptic normal.
  The live calculation is the elliptic deficit generated by nonzero `u`,
  which numerical probing indicates has a circle-grade hierarchy
  beginning with `-64 sum_k |u_k|² |c|^(2k)`.
  `proof/crabb_palindromic_normal_form.md`;
  `experiments/crabb_palindromic_normal_form.py`.

## NEWEST (2026-07-23): L123 identifies the entire quartic null exactly
- L122's phase-palindromic cone is not merely fourth-order flat.  Every
  positive point on it is an exact equality point in every size:
  ```
  sup_(||f||_D<=1) ||f(X(H))|| = 2,   t_*(X(H))=4.
  ```
- With `K=H+R*HR`, `A=2K^-1HR`, `q=He_0`, and
  `r=J conjugate(q)`, a diagonal phase gauge exposes the companion form
  ```
  Ae_0=0,  Ae_1=2e_0,
  Ae_j=e_(j-1)+2z_(j-1)(e_0-e_L).
  ```
- The exact metric
  `M=K-qq*+2rr*` satisfies `M-A*MA=qq*`.  Its generalized spectrum
  relative to `K` is `{1/2,1,...,1,2}`, so it gives the complete
  similarity upper square four.
- Write `det(xi I-A)=xi g(xi)`.  The observability determinant is
  `2^(L-1)det H`, so the roots of `g` are strictly in the disk.
  Consequently `B=g/g^sharp` is a finite Blaschke product.
  Cayley--Hamilton and the one-dimensional left/right kernels give
  `B(A)=4e_0r*`, whose `K`-norm is exactly two.  This is the matching
  scalar lower bound—not just an upper certificate.
- **Course correction completed:** there is no hidden higher-order pure
  disk descent on the quartic null.  Treat the phase-palindromic family
  as a stratified equality set: its fixed-phase branches meet singularly
  at the Crabb point.  Construct a uniform stratified normal form.
  The next gate is strict descent in the non-palindromic disk normal
  together with the elliptic normal, including their mixed remainder.
  `proof/crabb_palindromic_equality.md`;
  `experiments/crabb_palindromic_equality.py`.

## NEWEST (2026-07-23): L122 disproves the coercive disk-flat quartic
- Crouzeix/Lewis--Overton's polynomial support certificate yields an exact
  normalized disk chart.  For positive Hermitian `L x L` data `H`, extend by
  one zero coordinate and put
  ```
  K=H+R*HR,   X(H)=2K^(-1/2) H R K^(-1/2).
  ```
  The identity
  `K-(conj(w)HR+wR*H)=(I-wR*)H(I-conj(w)R)>=0`
  proves `W(X(H))=closed unit disk`.
- `H=I/2` gives `C_p`.  Hermitian Toeplitz curves
  `H=I/2+sZ(z)` form an exact `2p-4` dimensional disk submanifold tangent
  to the whole L115 disk-flat quotient.
- In coefficient coordinates, take `A=2K^-1HR`, `q=He_0`, and solve the
  rank-one Stein equation `M-A*MA=qq*`.  The endpoint generalized
  eigenvalues of `(M,K)` are
  ```
  lambda_-=1/2,
  lambda_+=2-16s^4 Q_L(z)+O(s^5),
  Q_L(z)=||z||^4-|z^T Jz|^2.
  ```
  Hence the feasible condition square is `4-32s^4Q_L(z)+O(s^5)`.
- `Q_L>=0` by Cauchy--Schwarz, but it is **not coercive**:
  `Q_L=0` iff `z=omega J conjugate(z)`.  In size four it is only
  `(|z_1|^2-|z_2|^2)^2`.
- An independent exact defect-vector expansion proves that optimizing
  L118's rank-one branch in size four does not fill this null:
  the correction is
  `16(|u_1|^2+|u_2|^2)+(32/3)|u_3|^2`, minimized at `u=0`.
- Exact noninfinitesimal palindromic samples through `p=8` retain
  generalized spectrum `{1/2,1,...,1,2}`, and SDPs return four, but this
  was only finite evidence at L122; L123 subsequently proves the all-size
  equality theorem and lower bound.
- **Course correction:** do not seek `-a||d||^4`.  L123 classifies the
  phase-palindromic null stratum as an exact `t_*=4` stratified equality
  family; analyze the elliptic and disk normals uniformly across it.
  `proof/crabb_disk_toeplitz_quartic.md`;
  `experiments/crabb_disk_toeplitz_quartic.py`.

## NEWEST (2026-07-23): L121 raises every disk-flat mixed linear term
- For L120's endpoint functional
  `F_c(Y)=(DP[Y])_(L0)-c^L(DP[Y])_(0L)`, root-of-unity filtering gives
  every coefficient explicitly by the circle grade `q=m-j-1`.
- The unrestricted first power is `c^floor(L/2)`.  If `L` is odd, its
  coefficient is exactly `Y_(L0)`; if `L` is even, it is exactly
  `sqrt(2)(Y_(L-1,0)+Y_(L,1))`.
- These are precisely L65's two bottom-mode obstructions.  They vanish on
  the full L65 equality space for `p>=4`; in the exceptional `p=3` case the
  circular-tangent constraint kills the same combination.
- Hence every L115 disk-flat direction satisfies
  ```
  F_c(Y)=O(c^(floor(L/2)+1)),
  D Gamma_p(C+cC*)[Y]=O(c^(L+floor(L/2)+1)).
  ```
- Writing `m=floor(L/2)+1`, a hypothetical coercive disk-flat quartic
  `-a||d||^4` absorbs the linear mixed term by Young because
  `4(L+m)/3>2L`; the remainder is then smaller than L117's axis margin.
- This proves the strict exponent for the **linear** mixed term, not the
  required disk-flat quartic or the higher mixed remainder.  Exact symbolic
  regeneration covers every matrix entry for `p=3,...,16`.
  `proof/crabb_flat_endpoint_selection.md`;
  `experiments/crabb_flat_endpoint_selection.py`.

## NEWEST (2026-07-23): L120 closes the pure elliptic strong tube in every size
- Let `L=p−1`, `r=c^L`, and
  `P_(L,c)(z)=2c^(L/2)T_L(z/(2sqrt(c)))`.  The outer endpoint compression is
  `Q=Pi P(A_c)Pi=[[0,2],[2r,0]]`.
- The polynomial fibre trace `E h=(1/L)sum_(P(z)=w)h(z)` obeys the exact
  conditional expectation `Pi h(A_c)Pi=(E h)(Q)`.
- In the exterior ellipse coordinate, the Crabb support vectors are
  `u(zeta)_j=d_j zeta^j/sqrt(L)`, independent of `c`.  Root-of-unity filtering
  proves the arbitrary-direction differential quadrature
  ```
  v(xi)* Pi DP(A_c)[Y] Pi v(xi)
    = (1/L) sum_(zeta_j^L=xi) P'(z_j) u(zeta_j)*Y u(zeta_j).
  ```
  Hence first numerical-boundary variations, and therefore first Riemann shape
  derivatives, commute with the fibre trace (up to a harmless disk automorphism).
- L119's full gradient consequently descends **exactly** to the classical
  `2 x 2` elliptic gradient:
  ```
  D Gamma_p(A_c)[Y]
    = (tau'(r)/2) Re((DP[Y])_(L0) − r(DP[Y])_(0L)),
  tau(r)=k(r²)/r.
  ```
- Dickson recurrence paths force the bracket to be
  `O(c^floor(L/2))`; Jacobi's product gives `tau'(c^L)=O(c^L)`.
  Thus the full gradient is `O(c^(L+floor(L/2)))=o(c^L)`.
- L118 now applies: optimizing every L65 coercive direction gains only
  `o(c^(2L))`, while L117 supplies the strict `−16c^(2L)` axis margin.
  Therefore `H_p(0,c)<0` for all sufficiently small nonzero `c`, in every size.
- An exact checker regenerates the reversal, conditional expectation through
  monomial degree `3L`, arbitrary-symbolic-direction quadrature, and sharp first
  path power for `p=3,...,10`.  It also matches the independent coefficients
  `−64c^5`, `−16c^7`, and `−32c^10`.
- **Next:** only L115's `2p−4` disk-flat variables and their weighted mixed merger
  with `c` remain in the single-Crabb local theorem.  Strong transverse variables
  are finished.
  `proof/crabb_descent_gradient.md`;
  `experiments/crabb_descent_gradient.py`.

## NEWEST (2026-07-23): L119 turns the transverse gradient into one scalar endpoint
- Fix an L117 axis point `T_c`, its rank-one defect `q_c`, and L116's
  Chebyshev--Blaschke product `B_c`.  The two differentiable functions
  `R_c(T)=||B_c(T)||²` and `U_c(T)=kappa(P(T,q_c))` satisfy
  `R_c<=t_*<=U_c` everywhere nearby and agree at `T_c`.
- Therefore `D R_c(T_c)=D U_c(T_c)` in **every** matrix direction.  By the
  envelope theorem this is also the derivative of L118's `Gamma_p+4`.
- Since
  `B_c(T_c)=sqrt(k(c^(2L))) diag(c^(j/2)) J diag(c^(−j/2))`
  has the simple top singular pair `e_0,e_L`, the full gradient is
  ```
  D Gamma_p(A_c)[Y]
    = 2 sqrt(k(c^(2L))) c^(−L/2)
      Re <e_0,D_A[B_c(phi_A(A))]_(A_c)[Y]e_L>.
  ```
  The prefactor tends to four.  Thus L118's all-size `o(c^L)` gate is exactly
  one endpoint functional-calculus shape derivative.
- Jacobi multiplication gives the additional exact descent
  `B_c o phi_c = phi_(c^L) o P_(L,c)`, where
  `P_(L,c)(z)=2c^(L/2)T_L(z/(2sqrt(c)))`.
  Moreover `P_(L,c)(A_c)=2c^(L/2)D J D^−1`; after permutation this is a
  nested direct sum of `2 x 2` elliptic reversal blocks whose outer numerical
  range is precisely `E_(c^L)`.
- A direct checker differentiates the rational Blaschke product and the Stein
  equation independently in 192 deterministic complex directions
  (`p=3,...,10`), with worst binary64 derivative discrepancy `9.2e-11`;
  it regenerates the polynomial descent to `1.8e-15`.
- **Completed by L120:** the endpoint derivative is
  `O(c^(L+floor(L/2)))` for arbitrary `p`, via the polynomial fibre trace,
  support quadrature, and Dickson path counting.
  `proof/crabb_touching_gradient.md`;
  `experiments/crabb_touching_gradient.py`.

## NEWEST (2026-07-22): L118 puts L65 and L117 on one analytic envelope
- For a stable pullback `T`, the rank-one Stein Gramian
  `P(T,q)=sum_n(T*)^nqq*T^n` is a feasible similarity metric.  At `C_p,e_0`,
  its condition number has positive defect-vector Hessian
  `8 sum_(j<L)|q_j|²+(8/3)|q_L|²`.
- The analytic IFT therefore selects a unique local minimizing defect vector `q_*(A)`;
  its condition number minus four is a real-analytic feasible envelope `Gamma_p(A)`.
- Equality in L62's Stein Schur complement proves that the matrix Hessian of `Gamma_p`
  is exactly L65's arbitrary-size nonpositive form.  L117's globally optimal rank-one
  metric tends to `P_0`, so it lies on this same branch and gives
  `Gamma_p(C_p+cC_p*)=k(c^(2p−2))/c^(p−1)−4`.
- Analytic maximization in L65's strong space leaves exactly L66's `2p−2` quotient
  coordinates; L115 splits these as `2p−4` disk-tangent variables and the one complex
  elliptic normal.  No independent high-dimensional metric chart remains.
- On the pure elliptic face, it is sufficient to prove that the strong gradient is
  `o(c^(p−1))`: completing the strong square then changes the envelope by
  `o(c^(2p−2))`, while L117 supplies `−16c^(2p−2)+o(c^(2p−2))`.
- A new sparse exact support-resolvent/Riemann/Stein engine proves that every dangerous
  weighted gradient coefficient vanishes for all real directions at `p=3,4,5`, and for
  selected grade-four/grade-six directions at `p=6`.  It exposes the first later terms
  `−16c^7` for `c^3E_30` at `p=4` and `−32c^10` for
  `c^4(E_30+E_41)/sqrt(2)` at `p=5`, matching independent numerical slopes.
- **Completed by L120:** the arbitrary-`p` cancellation and pure elliptic strong
  tube are proved.  Add L115's disk-flat mixed variables next.
  `proof/crabb_rank_one_envelope.md`;
  `experiments/general_crabb_weighted_series.py`.

## NEWEST (2026-07-22): L117 proves the exact elliptic Crabb axis in every size
- For `T_c=phi(C_p+cC_p*)`, `L=p−1`, L116's Chebyshev--Blaschke alternation gives
  `t_*(T_c)>=k(c^(2L))/c^L`.
- At the elliptic Lobatto nodes, the rank-one Szegő Gramian
  `G_ij=beta_i beta_j/(1−x_i x_j)` obeys `G−XGX=beta beta*`.  A Jacobi addition
  identity turns `G` into a Toeplitz-plus-Hankel kernel exactly diagonalized by DCT-I.
- Its normalized diagonal weights are
  `w_m=S_m/S_0`, `S_m=sum_(n in Z) sech((m+2Ln)(−log c))`.  Poisson summation
  identifies them as complementary-modulus `dn` values and proves
  `w_mw_(L−m)=k(c^(2L))`.
- Therefore `P=diag(w_m/c^m)` has a rank-one positive Stein defect.  The elementary
  termwise inequality `S_(m+1)>cS_m` proves
  `I<=P<=k(c^(2L))c^(−L)I`, matching L116's lower bound:
  `t_*(T_c)=k(c^(2L))/c^L<4`.
- This closes the all-size proof gap explicitly recorded on page 46 of Kenan Li's thesis.
  The novelty claim is narrow and pending publication-level review: the theorem concerns
  the fixed-weight one-parameter axis, while L20 covers arbitrary weights only in `4 x 4`.
- The SDP/eigensolver regression passes through `p=10`.  A separate 260-decimal explicit
  DCT/Jacobi checker audits every identity through `p=30` without an SDP.
- **Completed by L118--L120 on the pure elliptic face:** the weighted coercive
  absorption is proved.  L115's disk-flat mixed merger remains.
  `proof/crabb_elliptic_axis.md`; `experiments/crabb_elliptic_axis_theorem.py`.

## COMPLETED INPUT (2026-07-22): L116 supplies the sharp axis lower bound
- The degree-`L` Chebyshev--Blaschke product alternates on the elliptic Lobatto spectrum,
  and DCT reversal gives
  `B_L(T_c)=+-sqrt(k(c^(2L)))D J D^−1`.
- Hence `t_*(T_c)>=k(c^(2L))/c^L`; Jacobi's product makes this strictly below four.
  L117 supplies the matching all-size upper metric.

## NEWEST (2026-07-22): L115 leaves one complex soft normal at every Crabb size
- Lewis--Overton's local theorem says centered disk matrices near `C_p` form an analytic
  codimension-`2p` manifold.  Adding arbitrary centers leaves the circular-range manifold
  of codimension `2p−2`; its tangent equations are exactly the vanishing of support Fourier
  modes `2,...,p`.
- Intersecting these equations with L65's nonpositive second variation gives
  `dim((ker e_p intersect T_disk)/O_p)=2p−4`.  Since L66's full equality quotient has
  dimension `2p−2`, only two real normal coordinates remain soft.
- Grade by grade, mode one and every residual mode at least three are tangent to the
  circular-range manifold.  Mode two has exactly one complex normal quotient, represented
  by `C_p*`.
- The central normal family is exact:
  `W(C_p+c C_p*)={z+c conjugate(z): |z|<=1}`.  Thus the survivor is the elliptic deformation,
  while every circular anchor is already complete-`2` by Berger--Okubo--Ando.
- The larger-size campaign is therefore one Chebyshev--Lobatto elliptic family plus strong
  transverse directions, not `2p−2` independent higher-order jets.  The pure family is covered
  at `p=3` by L68--L73 and at `p=4` by L20.
- **Literature audit correction:** Kenan Li's thesis gives the candidate sharp diagonal
  similarity `t=k(c^(2p-2))/c^(p-1)<4`, but explicitly says the identities supporting its
  all-size construction were proved only for sizes `2,...,6` and numerically tested beyond.
  No later closure was found.  General `p` therefore remains a real proof target, followed by
  tubular absorption. `proof/crabb_disk_tangent_intersection.md`.

## NEWEST (2026-07-22): L114 proves a full repeated-`C3` neighbourhood
- Fix any finite multiplicity `m`.  L61 gives strict first-order descent whenever the
  support-compression Jensen gap is positive.
- On the zero-Jensen face, L86's complete endpoint is a sum of the negative matrix-Jensen
  defect, `−8H1²`, common-strong curvature, and three winner--loser Gram forms.  L87 identifies
  every residual kernel as a reducing flat block, now covered uniformly by L113.
- The only flat-core points with no endpoint margin are direct sums of L71's disk critical
  manifolds (up to unitary/affine coordinates).  L72 proves ambient stationarity there.
  Therefore on a complementary strong slice L86 is the first nonzero term and is coercive
  after all reducing flat kernels are split off.
- L105--L110 turn the finite jets into exact zero-slack metrics with strong-variable remainder
  `C(s||y||²+||y||³)`.  L82 supplies strict losing-sector feasibility; when a losing,
  winner, common-top, or flag gap collapses, its kernel is promoted to one of finitely many
  adjacent rank strata.
- Compactness of the normalized slice modulo copy unitaries gives a finite cover.  Hence for
  every fixed `m`, all matrices sufficiently close to `I_m tensor C3` have their numerical
  range as a complete `2`-spectral set.  The radius may depend on `m`.
- This is a genuine repeated-block local theorem, but not the general Crouzeix conjecture and
  not yet a theorem near larger `C_p` equality blocks.  The campaign now returns to the
  general equality-set/global attack.
  `proof/repeated_p3_local_neighbourhood.md`.

## NEWEST (2026-07-22): L113 closes the arbitrary-copy flat core
- For fixed copy multiplicity, iterate L93's common-top flag.  Its exact leakage form on each
  layer has kernel equal to the next layer; the stable layer reduces `Z` and has scalar
  support, hence splits into the one-/two-dimensional blocks closed by L110--L112.
- L111--L112 identify the entire tangent kernel at every stable stratum as unitary/blockwise
  motion.  On a complementary local slice, a finite hierarchical combination of L93's flag
  forms is therefore coercive: `−c s²||R||²`.
- The block-sign unitary fixes the block-diagonal base and sends each cross coupling
  `R` to `−R`.  Hence active diagonal endpoints are even in `R`; there is no dangerous
  linear transverse analytic term.  With a fixed normalized complementary support gap,
  the whole-domain support enlargement and Stein slack are `O(s²||R||²)`, and L105's
  post-leading remainder is `O(s³||R||²)`.
- If the support gap collapses, its kernel is promoted into the common-top space, moving to
  another member of the same finite flag stratification.  Compactness gives a finite cover for
  every fixed multiplicity (the radius need not be dimension-free).
- Thus the complete L87 flat copy core now has a genuine local complete-`2` theorem.
  The sole repeated-`C3` gate is to glue it through L86's generator-one, common-strong,
  winner--loser, and positive first-Jensen variables.
  `proof/repeated_p3_flat_flag_tube.md`.

## NEWEST (2026-07-22): L112 splits every normal multiplicity collision
- At a normal full-common-top scalar-support point, rotate and scale to
  `Z0=s diag(Ip,−Iq)`.  For a tangent `X=[[A,B],[C,D]]`, the two linearized scalar
  relations force exactly `A=aIp,D=−aIq`; the rectangular blocks `B,C` remain free.
- A skew-Hermitian copy-orbit generator with upper block `L=C*/(2s)` kills the lower
  rectangle and leaves the single invariant edge `R=B−C*`.
- The SVD of `R` splits the tangent into independent blocks
  `[[s+ta,t sigma_j],[0,−s−ta]]` plus unmatched normal coordinates.  These are exactly
  L110's trace-zero two-copy terminal tubes and L88's normal stratum.
- Exact real-linear ranks for every `1<=p,q<=4` give dimensions `4pq+2` full,
  `4pq` homogeneous, and `2pq` unitary orbit; a symbolic `2 x 3` calculation verifies
  the gauge identity entry by entry.
- Together L111--L112 close all full-common-top multiplicity strata, including their normal
  boundary.  The remaining flat-core gate is a proper common-top space with collapsing L93
  flag ranks.
  `proof/repeated_p3_normal_collision.md`;
  `experiments/repeated_p3_normal_collision.py`.

## NEWEST (2026-07-22): L111 removes hidden multiplicity tangents
- L93's full common-top blocks satisfy `Z²=alpha I` and
  `ZZ*+Z*Z=beta I`; every nonnormal irreducible is `2 x 2`.  The possible remaining
  obstruction was a multiplicity-`k` direct sum of equivalent irreducibles.
- Linearizing both relations at a square-zero block gives
  `X=[[A,B],[cI,−A]]` with `B+B*=bI`.  The arbitrary `A` and skew part of `B` are a
  unitary commutator; the real `b` and complex `c` are exactly the same three parameters on
  every two-dimensional summand.
- At an invertible nonnormal block, a Pauli decomposition makes the homogeneous relation
  tangent `3k²`-dimensional.  The unitary commutator orbit also has dimension
  `4k²−k²=3k²`, so the spaces coincide.  Allowing `dot alpha` and `dot beta` again adds only
  three repeated single-block parameters.
- Exact symbolic real-linear ranks verify both strata for multiplicities `1..4`.
  Thus there is no new large-multiplicity terminal family beyond unitary basis motion and
  L110's two-copy tube.  The next subgate is the collision with the normal stratum and proper
  common-top metric flags.
  `proof/repeated_p3_scalar_support_rigidity.md`;
  `experiments/repeated_p3_scalar_support_rigidity.py`.

## NEWEST (2026-07-22): L110 closes the uniform two-copy terminal tube
- Pinch a terminal matrix `A=N+E` to its block-diagonal normal retraction.  Since every flat
  edge has zero compression to the base top support cluster, uniform cluster perturbation gives
  `0<=h_A−h_N<=C delta(r+delta)`; the lower bound is the exact inclusion `W(N)⊂W(A)`.
- Compare the Riemann maps of the **whole** normal and actual numerical ranges.
  The inclusion-induced disk map `h=f_A∘f_N^−1` fixes zero, has
  `h'(0)>=1−C delta(r+delta)`, and Schwarz--Pick makes it
  `O(delta(r+delta))` from the identity on compact subdisks.
- Apply `h` to L108's zero-slack normal certificate.  The same metric still contracts
  `f_A(N)=h(f_N(N))` and acquires a block-diagonal PSD Stein Schur slack
  `H_delta=O(delta(r+delta))`.  This whole-normal slack vanishes with the transverse edge,
  unlike a per-block slack that may already be present on the normal stratum.
- In the analytic endpoint difference between the actual zero-slack chart and that inherited
  normal chart, every Taylor monomial contains a transverse variable or `H_delta`.  Giving
  `H_delta` weight two, the post-cubic remainder is `C delta(r+delta)^3`.
  L109 cancels the quadratic exactly and L100 leaves
  `−delta(delta²+2r²)I/4`, which dominates uniformly because
  `delta²+2r²>=(2/3)(r+delta)²`.
- L95--L100 cover the compact complement of the sharp weighted center.  Hence the entire
  two-copy terminal chart is now locally complete-`2`.  The next live gate is uniformizing
  L93's arbitrary-copy metric flag and merging the resulting tube with L86's losing-space
  gaps.
  `proof/repeated_p3_terminal_tube.md`;
  `experiments/repeated_p3_slack_transfer.py`.

## NEWEST (2026-07-22): L109 proves the joint normal/transverse cubic jet
- The previous warning was essential: neither the frozen normal endpoint nor the transverse
  increment may be discarded.  At L100's real weighted center, the actual inherited
  frozen-normal Stein slack begins with
  `H2=diag(4Delta m,8Delta m,4Delta m,8Delta m)`, `Delta m=5a²/128`, and an exact mixed
  cubic coefficient `gamma=15sqrt(2)a²d/256` in `H3`.
- Retaining a fraction `theta` of this slack gives
  `E2(N)=−5(1−theta)a²I/8` but `E2(A)=5theta a²I/8`.  Thus freezing the slack creates the
  observed positive actual quadratic endpoint, while zero slack makes the negative normal
  quadratic pay for the positive transverse quadratic exactly.
- Both frozen-normal and actual cubic endpoints are
  `E3=−16(1−theta)m3I`.  At zero slack this preserves L100's physical bound
  `−delta(delta²+2r²)I/4`; full slack is flat through cubic order.
- The checker also proves that full inherited slack reproduces the unperturbed normal metric
  through third order.  This closes the finite-jet accounting error.  The remaining terminal
  gate is a genuinely analytic one: combine L105's Taylor majorant with L106's explicit
  transverse factor to make the exact remainder `o(delta(r²+delta²))`.
  `proof/repeated_p3_exact_metric_chart.md`;
  `experiments/repeated_p3_slack_transfer.py`.

## NEWEST (2026-07-22): L108 selects the sharp zero-slack branch
- L107 makes inherited Stein slack legal, but freezing all of it is not the sharp transverse
  certificate.  Exact-chart probes at six terminal scales made the frozen-slack upper endpoint
  positive, while `H=0` matched the full SDP deficit.  This is diagnostic evidence, not the proof.
- On the block-diagonal normal anchor, the chart decouples copy by copy.  For one PSD slack
  block `H`, the scalar upper endpoint has base derivative
  `D_H E=h22+2h11>=tr H>=||H||`.
- Uniform analyticity therefore gives `E(H)−E(0)>=(1−Cs)||H||>=0` nearby.  Since L106's
  inherited endpoint `E(HN)` is nonpositive, the zero-slack endpoint is also nonpositive.
  Direct sums make this dimension-free on the normal stratum.
- L109 subsequently computes the complete weighted slack jet and proves the required
  normal/transverse quadratic cancellation.

## NEWEST (2026-07-22): L107 embeds the inherited frozen-normal slack
- Extend L105's Stein Schur equation from zero to a prescribed Hermitian right side `H`.
  The implicit-function Jacobian is unchanged, and Stein positivity is now exactly `H>=0`.
- Scale each L106 single-block metric to minimum eigenvalue one.  Near the Crabb metric the
  minimum is simple and the physical-level range block of `Pj−I` is positive, so singular
  positivity forces its lower Schur complement to vanish.
- Its Stein Schur complement `Hj` is PSD even when the intervening disk self-map is strict.
  Direct sums therefore give the exact identity `PN=P(TN,BN,HN)` in the extended chart.
- Keeping `HN` fixed while moving from `TN=f(N)` to `T=f(A)` preserves lower and Stein
  feasibility exactly, but need not preserve the upper bound.  L108 proves that one can instead
  tighten the normal anchor safely to `H=0`, the branch used by L98--L104.
  `proof/repeated_p3_exact_metric_chart.md`.

## NEWEST (2026-07-22): L106 bypasses nonsmooth domain differentiation
- For the actual perturbed matrix `A=N+E`, freeze its own Riemann map `f=phi_A` and compare
  `f(A)` with `f(N)`.  Every diagonal block `Bj` of `N` is a compression of `A`, so
  `W(Bj)⊂W(A)`.
- L73 supplies a metric contracting `Sj=phi_j(Bj)`.  Since
  `f(Bj)=(f∘phi_j^−1)(Sj)` and the composition is a disk self-map, the same metric contracts
  `f(Bj)` by von Neumann.  Their direct sum is an exact normal certificate for `f(N)` using
  the actual perturbed domain.
- A fixed circle lies inside every nearby `W(A)` and surrounds both spectra.  With `|f|≤1`,
  the resolvent Neumann series expands `f(N+E)−f(N)` into homogeneous terms with bounds
  `C1(C0||E||)^k`.  Every remainder has an explicit transverse factor, uniformly over the
  nonsmoothly varying maps.
- A non-automorphic intervening disk self-map can create Stein slack.  L107 embeds its
  nonnegative Stein Schur complement as an exact chart parameter; L108 safely tightens it and
  L110 subsequently proves the quantitative two-regime endpoint estimate.
  `proof/repeated_p3_frozen_domain.md`.

## NEWEST (2026-07-22): L105 replaces forced metric series by one analytic chart
- Fix the whole level-zero/range cross block `B` of the metric.  The exact formula
  `P00=I+B(C−I)^−1B*` makes `P−I` tight.
- The tight Stein Schur-complement equation has derivative
  `(X11,X12,X22)->(X11,X12,X22−2X11)` in the remaining range block `C`, an invertible
  real-linear map for every copy multiplicity.  The analytic implicit-function theorem
  therefore supplies a unique `C=C(T,B)`.
- Lower and Stein positivity are then automatic, and `P≤4I` is equivalent to one analytic
  upper endpoint `E(T,B)≤0`.  Thus the forced metric coefficients converge for every analytic
  choice of `T,B`; only the uniform free-block selection and scalar endpoint sign remain.
  `proof/repeated_p3_exact_metric_chart.md`;
  `experiments/repeated_p3_exact_metric_chart.py`.

## NEWEST (2026-07-22): L104 removes every finite-order traceless recurrence
- At any metric order `k>=2`, the homogeneous free block
  `Delta Pk=[[0,U,0],[U*,0,2U],[0,2U*,0]]` preserves the active lower and
  Stein kernel compressions and does not change the order-`k` upper endpoint.
- At order `k+1`, coefficient convolution pairs it only with the fixed first-order normal
  data.  Thus its endpoint transfer is independent of `k` and equals
  `5sqrt(2)d[[0,U01−conj(U10)],[conj(U01)−U10,0]]` in the real canonical chart.
- On a `J`-equivariant recursive branch, copy sign symmetry forces every coefficient linear
  in the transverse Schur edge to be Hermitian off-diagonal.  Since the transfer is onto that
  space for fixed `d!=0`, every such traceless term can be canceled at every finite order.
- This replaces an infinite coefficient grind by one precise analytic gate: prove convergence
  or direct transverse factorization while retaining the negative nonsmooth scalar Jensen term
  uniformly as `d->0`.  Products of a free block recur at later orders and the inverse grows
  like `1/|d|`; L104 alone is a formal fixed-`d` result, not a neighbourhood theorem.
  `proof/repeated_p3_flat_two_copy_weighted.md`;
  `experiments/repeated_p3_flat_two_copy_fourth_metric.py`.

## NEWEST (2026-07-22): L103 cancels the full transverse fourth endpoint
- An independent audit caught an important scope gap in L102: scalarity of the fourth
  Feshbach **support** coefficient does not by itself control the full similarity metric.
  With the default tight metric, exact propagation produces the linear transverse term
  `(25d³/8)[[0,1],[1,0]]`.
- The relevant freedom lies in the level-zero/level-one block `aU` of the **third** metric.
  The full derivative is
  `(5d/8){5d²+8sqrt(2)(U01−U10)}[[0,1],[1,0]]`.
  Choosing `U=5d²[[0,−1],[1,0]]/(16sqrt(2))` cancels it exactly.
- The identity retains arbitrary third conformal-response coefficients and the surviving
  normal fourth coefficients; all cancel.  The normal endpoint remains scalar, so the free
  block, which vanishes on the normal face, is compatible with L88.
- Therefore the **total** fourth endpoint relative to the exact normal stratum starts at
  `O(delta²(r+delta)²)` after simultaneous metric selection and is absorbable by L100's
  `delta(r²+delta²)` cubic gap.  L110 subsequently factors the full analytic remainder; the
  remaining debt is the arbitrary-copy L93/L86 flag lift.
  `proof/repeated_p3_flat_two_copy_weighted.md`;
  `experiments/repeated_p3_flat_two_copy_fourth_metric.py`.

## NEWEST (2026-07-22): L102 isolates the fourth support coefficient
- At L100's center, the full fourth Feshbach coefficient includes an energy-dependent
  `−M2 N R² N` correction in addition to the five ordinary resolvent words.
- Exact compression proves that its traceless copy part is zero.  Differentiating the complete
  coefficient in the nonnormal Schur edge gives zero **full** derivative.
- Hence the weighted normal chart is copy-scalar through fourth order at the support level,
  and its raw quartic traceless splitting is at least `O(delta²(r+delta)²)`.
  L102 alone makes no claim about the simultaneous similarity metric; L103 supplies that step.
- A separate nonlinear probe at map resolutions 2048/4096 found
  `(4−t*)/[delta(r²+delta²)]` between `0.637` and `0.671` on seven two-scale cases, with no
  sign reversal.  This is supporting evidence only; the exact result is the Feshbach identity.
  `experiments/repeated_p3_normal_center_probe.py`.
- `proof/repeated_p3_flat_two_copy_weighted.md`.

## NEWEST (2026-07-22): L101 closes the support-crossing regularity debt
- If convex domains near the disk have
  `h_epsilon=1+epsilon s+o(epsilon)` with merely continuous `s`, their radial functions have
  the same first variation.  Applying the Schwarz integral to
  `log(Psi_epsilon(z)/z)` then gives the normalized Riemann-map tangent used in L61.
- The argument is uniform for compact families of continuous profiles and needs neither a
  differentiable largest-eigenvalue branch nor strict convexity of the perturbed domain.
- At repeated Crabb blocks, finite-dimensional degenerate perturbation theory gives a uniform
  `O(epsilon²)` support remainder for exact linear paths (and paths with an `O(epsilon²)`
  matrix remainder); a general `o(epsilon)` matrix remainder gives the `o(epsilon)` support
  remainder the theorem needs.  The `lambda_max` profile is uniformly Lipschitz, so L61's
  operator tangent and Dini bound hold even at support crossings.
- The radial/logarithmic argument can be reapplied after lower-order analytic recentering, so
  conformal regularity is no longer the obstruction in the weighted normal-face patch.  The
  remaining debt is to organize and absorb the finite hierarchy of endpoint remainders.
  `proof/general_similarity_nonsmooth_tangent.md`.

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
  The closest-source audit found Kenan Li's all-dimensional candidate formula for the
  fixed-weight Crabb-derived family (whose every-size proof gap is now closed by L117);
  no prior arbitrary-weight 4×4 theorem was found. Call L59 apparently new pending a
  publication-level novelty audit.
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

## Current next actions (Epoch 6, refreshed 2026-07-24)
1. **A178: extract an arbitrary-grade recurrence from L230--L234.**
   The cubic-through-sextic columns are four separate finite
   certificates, not an induction.  Seek a uniform recurrence in the
   partial right ideal that simultaneously produces the next odd
   cancellation, the next even transfer Gram, and a uniformly positive
   Schur-budget margin.  It must also control column-norm growth well
   enough to converge in L194's analytic chart.  Do not compute another
   isolated grade unless that computation exposes this recurrence.
2. **Use, rather than silently replace, A171/L228.**  The proposed
   one-delay associated-graded anticommutator identity remains open and
   is still the intended structural route to arbitrary grade.  Determine
   whether L237's unit first-reflection coefficient, acting on L236's
   active Hankel cell, generates the right-ideal identities and
   even-face budget observed in L230--L234.  L238 has compressed the
   full deflated resolvent to a `2m x 2m` endpoint scattering matrix
   whose off-diagonal blocks are the two transfer orientations.
   Evaluate the full L125 theta/ODE contour coefficients and L219
   metric/right-defect Schur square on that matrix to prove the
   first-active pairing.  The whole-series shift is disproved.
   L220/L221 supply exact model delay covariance, but not this final
   physical scalar multiplier.
3. **Prove an all-series analytic bound.**  Any recurrence must replace
   the observed column bounds `3, 27/2, 48, 181` by an explicit
   grade-`k` estimate that gives a genuine local analytic metric, not
   merely formal finite jets.  Keep L194's exact positivity and lower
   endpoint throughout.
4. **Finish the repeated normal/elliptic merger after A178.** L197 closes the
   disk Schur flag and L199 closes only its first raw circular-normal
   face.  Prove the later Schur-orthogonal normal response, then merge
   it with the all-grade elliptic flag without spending a negative
   square twice.
5. **Retain the full CP correction as fallback.** Test whether L21's trace
   inequality follows from block-Toeplitz positivity of the complete
   operator-valued Crouzeix--Palencia correction.  Do not retry
   scalar shifts or positive-state scalarizations.
6. **H-r/analytic-flag fallback and guardrails.** Continue the shifted Möbius and odd-phase
   level-four positivity attacks only after the local Crabb merger is
   banked; definite parity is false.  Reconcile L208's exact kernel
   projections with L197's pathwise analytic Schur flag, and never
   insert a discontinuous pseudoinverse into the metric.

Standing negative constraints: do not superpose the raw L212
representatives (A179); do not retry scalar shifts or positive-state
scalarizations; do not put a discontinuous pseudoinverse into the
metric; and do not assume definite parity.
Keep committing+pushing after each task (user instruction).

## Files map (handoff-ready, 2026-07-22)
proof/ — read in this order for the current frontier:
  crabb_disk_normal_tube.md (L152 uniform disk-normal anchor),
  crabb_disk_flat_elliptic_face.md (L151 weighted raw face),
  crabb_reflected_hardy_lift.md (L149--L150 equality tube),
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
  crabb_disk_tangent_intersection.py (L115 arbitrary-size circular-tangent intersection),
  crabb_elliptic_axis.py (L116/L117 SDP and explicit-metric regression),
  crabb_elliptic_axis_theorem.py (L117 high-precision DCT/Jacobi identity audit),
  general_crabb_weighted_series.py (L118 arbitrary-size sparse support/Riemann engine),
  crabb_transverse_weighted_gradient.py (L118 exact low-size gradient cancellations),
  crabb_touching_gradient.py (L119 upper/lower derivative and polynomial-descent audit),
  crabb_descent_gradient.py (L120 exact fibre trace/quadrature/path audit),
  crabb_flat_endpoint_selection.py (L121 exact disk-flat bottom-mode selection),
  crabb_disk_toeplitz_quartic.py (L122 exact disk chart/noncoercive quartic),
  crabb_palindromic_equality.py (L123 exact disk equality family),
  crabb_palindromic_normal_form.py (L124 exact stratified normal form),
  crabb_palindromic_elliptic_face.py (candidate graded elliptic Newton face),
  crabb_palindromic_elliptic_hessian.py (guarded exact amplitude/defect Hessian engine),
  crabb_mixed_grade_face.py (A94 exact finite polarization),
  crabb_spectral_sine_modes.py (L135 DCT/DST bridge),
  crabb_gcd_dickson_descent.py (L136 common-divisor reduction),
  crabb_principal_face_locality.py (A98 focused coprime grade-four/completed-square audit),
  crabb_transport_correction_scan.py (A98 raw-coordinate falsification guard),
  crabb_defect_hessian_factor.py (A98 universal-Hessian LDL edge audit),
  crabb_axis_defect_fourier.py (L137 reciprocal-dn defect edge),
  crabb_homogeneous_defect_normal_form.py (L138 exact DCT defect congruence),
  crabb_faber_blaschke_dual.py (A99 corrected scalar dual square),
  crabb_central_faber_blaschke.py (L142 central corrected dual face),
  crabb_faber_blaschke_formal.py (L143 noncentral formal dual face),
  crabb_faber_blaschke_mixed.py (L144 mixed-grade polarization),
  crabb_faber_blaschke_complementarity.py (L145 primal/dual defect),
  crabb_equality_normal_stationarity.py (L148 all-anchor normal derivative),
  crabb_finite_amplitude_dual_filtration.py (all-amplitude dual support probe),
  crabb_uniform_remainder_probe.py (candidate L147 falsification probe),
  crabb_disk_flat_elliptic_face.py (L151 direct full-gauge Hessian audit),
  crabb_raw_faber_endpoint.py (L151 raw characteristic/endpoint audit),
  crabb_raw_faber_blaschke.py (L151 dual/model-complement face audit),
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
  repeated_p3_flat_two_copy_fourth.py + repeated_p3_flat_two_copy_fourth_metric.py
  (L102 support persistence + L103 simultaneous fourth-metric cancellation),
  repeated_p3_exact_metric_chart.py (L105 exact tight lower/Stein chart Jacobian),
  repeated_p3_slack_transfer.py (L109 weighted normal/transverse slack jet),
  repeated_p3_scalar_support_rigidity.py (L111 exact multiplicity-tangent ranks),
  repeated_p3_normal_collision.py (L112 rectangular normal-collision reduction),
  repeated_p3_normal_center_probe.py (supporting two-scale regression only),
  repeated_p3_flat_metric_flag.py (L93 arbitrary-copy metric-flag derivative).
Proof artifact: experiments/positive_tail_full_20260722.log (L59 clean 10-box run).
Data: sym3_sweep_s51.jsonl (40 rec), sym4_sweep_s61.jsonl (20 rec, ρ column trustworthy),
general_similarity_equality_s9173401.jsonl plus its `sensitivity` and `ultralocal` companions
(L21 equality gate); general_similarity_tangent_s9173401.jsonl and
general_similarity_tangent_all_s9173401.jsonl (L61 cross-solver tangent data);
general_similarity_second_order_all_s9173401.jsonl and its `halfstep` companion (L62 data);
general_crabb_second_order_modes_s70221.jsonl (L65 sizes `3..30` mode/rank audit);
crabb_second_order_equality_s70221.jsonl (L66 sizes `3..8` quotient audit);
crabb_disk_tangent_intersection_s70222.jsonl (L115 sizes `3..7`);
crabb_elliptic_axis_s70222.jsonl (L116/L117 sizes `3..10`, four ellipse parameters);
crabb_elliptic_axis_theorem_s70222.jsonl (L117 sizes `3..30`, six ellipse parameters).
Ledgers: LEMMA_LEDGER.md, APPROACH_LEDGER.md (pitfalls P1–P8 — READ BEFORE ANY SEARCH),
LITERATURE_LEDGER.md, COUNTEREXAMPLE_SEARCH.md. Audit: chatgpt/FABLE_RESEARCH_AUDIT.md
(reconciled 2026-07-20). Restart: checkpoints/RESTART_PACKET.md (paste-ready instruction).
