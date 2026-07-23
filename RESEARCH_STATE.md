# RESEARCH_STATE.md — Crouzeix Conjecture Campaign

**Last updated:** 2026-07-23 (Epoch 6 — Faber endpoint localization)

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

## Current next actions (Epoch 6, refreshed 2026-07-23)
1. **Localize L130 from central windows to every palindromic offset.**
   For an offset `k` in a larger `p=L+1` block, isolate the minimal
   `2k+1` endpoint/offset path window, prove that unused tails are strictly
   inactive in the transferred metric, and quantify the first feedback order.
   Then polarize distinct circle grades.  The target is A84's full leading face
   `-64a² sum_k |u_k|²c^(2k)` with a uniform weighted remainder.  Do not
   redo the central coefficient or infer tail inactivity from numerics.
2. **Build the L115/L117 transverse tube.**  Anchor on the exact elliptic-axis metric,
   and use L120's completed all-size absorption of L65's coercive normal complement.
   What remains is the weighted merger of the elliptic coordinate with L115's
   `2p−4` disk-flat coordinates, using the exact circular-range anchors rather than
   estimating them as generic flat Taylor directions.  Do not reopen the strong-gradient
   or differentiated-Stein calculations.
   L121 raises the linear disk-flat/elliptic coupling by one full power, and L122 proves
   that the hoped-for coercive quartic is false.  L123 now classifies the full
   phase-palindromic null as an exact `t_*=4` stratified equality family.  L124 gives
   the exact best-phase normal form `Q=4||u||²||v||²`.  Prove the observed graded
   elliptic margin along `u`, then combine it with the disk-normal `v` margin uniformly
   where all branches meet at the Crabb point.  Do not recompute the endpoint derivative
   or search for pure disk descent along the equality tangent.
   Do not rely on a fixed positive margin or compute the old `2p−2` residual jets.
3. **Retain the full CP correction.** Hartz--McCarthy cannot scalarize it. Test whether L21's
   trace inequality can instead be derived from the block-Toeplitz positivity of the full
   operator-valued correction moments. Do not retry trace/positive-state scalarizations.
4. **Shifted Möbius H-r fallback**: derive its exact stationarity/rho formula (Kenan-Li quartic
   analog) and prove rho≥0 or K≤2. Definite parity is false.
5. **Bi-conic Schwarzian test**: compute `SG` for the off-slice collapsed map on the critical
   real interval. `SG≥0` would extend L17 immediately; otherwise test the weaker Sturm-potential
   comparison that the proof actually needs.
6. **Odd phase positivity**: ρ = B₁g₁q₁+B₂g₂q₂ ≥ 0 given the L15 stationarity law
   (3-parameter; interlacing τ₂ < α < τ₁; term-1 dominance observed). Try the same
   deformation/kernel machinery.
7. Rigor debts: n=6 structured floor; 2×2 α=0; contact degeneracy; L59/L117 publication-level
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
