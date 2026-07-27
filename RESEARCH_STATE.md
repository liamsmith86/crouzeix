# RESEARCH_STATE.md — Crouzeix Conjecture Campaign

**Last updated:** 2026-07-26 (Epoch 6 — kernel-Gram route falsified after L350)

## CANONICAL LIVE GATE (2026-07-26, after L350/A298)

There is exactly one current mathematical gate:
**sign the `2n−2` boundary-shape Schur form at every finite
nondegenerate Gau--Wu equality model.**  L342 first reduces the full
second-order L21 SDP exactly to

`e_phi(C)=min_(x in C^(n−1)) F_C(x)`,

where `x` is one first-metric boundary row and its free curvature is
the model-independent diagonal form
`F_0(x)=8||x_mid||²+(8/3)|x_L|²`.  If `J_C(h)` is the scalar norm
Hessian for a moving inner tangent, L342 proves the exact square gap

`F_C(x)−4J_C(h)
 =||r_-||²_(diag(1,...,1,3))
  +||r_+||²_(diag(3,1,...,1))`.

The endpoint residual map in `(x,h)` is bijective, hence
`min_x F_C(x)=4 max_h J_C(h)`.

L343 now removes the quadratic-size bulk from this common form.
The L339 polynomial null state gives an exact support frame which
places every such model on L187's inverse-Toeplitz disk chart, without
a near-Crabb restriction.  General-H disk directions transverse to
the equality locus form a `(n−2)^2`-dimensional subspace
`D_phi`; Okubo--Ando and L183's linearized Hardy residual make the
common Hessian strictly negative on it.  The canonical shape map
(with `delta(t)=det(I−tS)` and
`omega_phi=|delta(conj(zeta))|²D_f`)

`Sigma_phi(C)=[omega_phi s_C] mod omega_phi span_R{1,cos(theta),sin(theta)}`

has rank `2n−2` and kernel exactly `D_phi`.  Define

`Hhat_phi(sigma)=max_(Sigma_phi(C)=sigma) e_phi(C)`.

The live theorem is precisely `Hhat_phi(sigma)<=0`.  A292 supplies
the exact, better-conditioned conformal coordinate

`A_phi(C)=2(hat(s_C)(2),...,hat(s_C)(n)) in C^(n−1)`.

It is an isomorphism on L343's shape quotient.  In fifteen complete
models the Schur form is phase covariant in this coordinate within
`3.85e−12`, hence is the realification of one negative Hermitian
`(n−1)`-square.  The phase covariance and sign are evidence, not
theorems.  L345 replaces the phase question by one exact dual
criterion.  If `J_phi` is the uneliminated joint scalar Hessian,
`S=[S_x;S_y]` is the conformal shape map, and

`W_phi=J_phi^−1S_y*+iJ_phi^−1S_x*`,

then, whenever `J_phi` is nondegenerate, phase covariance is
equivalent exactly to
`W_phi^T J_phi W_phi=0`.  Numerically, the first Blaschke image of
this dual range has zero lower triangle and the two L336 endpoint
forms are separately bilinearly isotropic there.  L346 proves that
the full lower response is generated exactly by its `n` diagonal
entries and southwest corner.  Prove only

`diag Y_1(W_phi)=0` and `(Y_1(W_phi))_(L0)=0`

from the polynomial support frame and L344's exact transport

`eta_C(zeta)=Hhat^(dagger/2)(I−zeta R*)^−1 Ccal*
 (H_C(zeta)y_zeta−s_C(zeta)y_zeta)`

and polarize L336 on that flag.  Then identify the positive Gram for
the negative Hermitian matrix.  This is a
fixed inverse-Toeplitz Hardy norm, so do not return to L339's moving
pointwise pseudoinverse, re-expand the already strict disk block, or
infer phase covariance merely from the nearly paired spectrum in
the weighted coordinate.  A universal-diagonal-plus-rank-one
inverse ansatz was falsified after resolving its dominant soft
component; do not reopen it.  A299 also falsifies every post-pivot
estimate `S_phi>=cD G_g^−1D` with fixed `c>=1`: weight three fails
at moderate rational roots, and weight one already fails at
separated boundary roots.  The full shape matrices remain positive,
so this rejects only the isolated model-kernel comparisons.  Keep
L342's two endpoint residuals coupled to L344's fixed port and the
disk-fibre elimination.
L347 now identifies the cyclic condition exactly:

`(Y_1)_(L0)=4f'(0) hat(sigma_W)(−n)/gamma`,

where `sigma_W=omega_phi s_W` and `gamma` is the nonzero endpoint
normalization of L343's polynomial frame (`gamma=1` in the standard
frame).  Pure zero motion has no cyclic corner.  Therefore derive
vanishing of this last negative Hardy coefficient from L344's fixed
port; do not try to repair it by moving a Blaschke zero.
L348 organizes every remaining diagonal condition.  Once the cyclic
coefficient vanishes, the `n−2` interior equations say that the
moving zeros track the simple spectral velocities
`mu_j=tr(P_jG)`, the sum of the endpoint equations says that the
last zero tracks `tr(P_0G)/2`, and only the endpoint difference
`(Y_1)_(00)−(Y_1)_(LL)` remains.  These tracking identities are
observed on `W_phi` but are not yet proved there.
L349 closes the endpoint sum exactly.  Postcomposition by a disk
automorphism gives pure zero tangents `k_1,k_i` for which

`tau=(Y_1)_(00)+(Y_1)_(LL)
 =(4/3)(k_1^T J_phi+i k_i^T J_phi)`.

The L345 source is zero on pure zero motion, so `tau(W_phi)=0`.
After the cyclic coefficient is closed, L348's repeated-root mean
tracking follows automatically.  Do not carry the endpoint sum as
an open Euler debt.
L350 closes the endpoint difference exactly.  The canonical
model-space conjugation exchanges the endpoint vectors, fixes the
Gau--Wu operator under `T -> C T* C`, and preserves both the joint
Hessian and the complete shape source.  It therefore fixes
`W_phi`, making `Y_1(W_phi)` complex symmetric under endpoint
exchange.  Thus its two endpoint diagonal entries agree; L349 says
their sum is zero, so both entries vanish.  The only live
lower-flag debts are now the cyclic coefficient and the `n−2`
simple-root equations.
L341's paired frames and L340's rank law remain equivalent response
diagnostics, not separate concurrent gates.  A proposed factor must
still reproduce their endpoint rank `4n−8`, support increment
`2n−6`, total rank `6n−14`, and nullity `(n−4)²`.  Keep both
endpoints coupled: the one-sided signs are false.  Do not replace the
exact quotient by another dimension-by-dimension SDP grind.  This is
current next action 1 below.  Every later
occurrence of “historical gate,”
“historical frontier,” or “remaining historical task” records what
was unresolved at that dated checkpoint; none is a concurrent
frontier.  The CP/H-r routes remain parked fallbacks.

## L350/A298 EXACT CANONICAL-CONJUGATION ENDPOINT CLOSURE (2026-07-26)

- The canonical model-space conjugation exchanges `p=f` and `q=1`.
  Since it also sends the Gau--Wu weight `X` to `X^−1`, the
  equality operator satisfies `C A* C=A`.
- The involution `rho(T)=C T* C` preserves the numerical range,
  functional calculus, and norm.  On the joint tangent it therefore
  gives `R^T J_phi R=J_phi` and `S_phi R=S_phi`.
- Hence `J_phi^−1S_phi*`, and in particular `W_phi`, is fixed by
  `R`.  Its first image satisfies
  `Y_1(W_phi)=C Y_1(W_phi)* C`, so the two endpoint diagonal entries
  are equal.
- L349 gives their sum zero.  Both endpoint entries therefore vanish
  exactly, without phase covariance or a sign assumption.
- The live lower-flag debt is now only one cyclic equation and
  `n−2` simple-root tracking equations.
- Twelve audits through `n=8` give image-symmetry and endpoint
  residuals below `3.16e−12` and `4.49e−14`; dataset SHA-256
  `8e317c98ea8216aae05b979cff0fda8b6766b25980e19b9c3e39e34ee96fed8b`.

## A299 NUMERICAL FIXED MODEL-KERNEL LOWER-BOUND FALSIFICATION (2026-07-26)

- After the first shape pivot, every proposed sufficient estimate
  `S_phi>=cD G_g^−1D` with fixed `c>=1` is numerically false.
- At the separated rational `n=4` roots
  `3/8+i/4,−1/2−i/4`, its gap eigenvalues are approximately
  `−0.11585777,1.17891488`, while the complete negative shape matrix
  is positive definite with minimum eigenvalue `0.06403464`.
- The first pivot matches `16|b_1b_2|²=1.015625`; independent
  resolutions 512, 1024, and 2048 agree to numerical precision.
- At the separated boundary roots
  `−41/64−47i/64,−34/64−52i/64`, the weight-one gap already has
  eigenvalues approximately `−275.39428827,246481.1030077`, while
  the complete shape matrix remains positive definite.
- This falsifies only a stronger Gram shortcut.  It leaves the live
  phase/lower-flag identities and the conjecture untouched.
- Dataset SHA-256
  `41352e20a5af43b1d4d1d3d6bfcb6ad193e055abfbb93c902928be08f2cd9d7a`.

## L349/A297 EXACT POSTCOMPOSITION EULER IDENTITY (2026-07-26)

- Postcomposing `f` with
  `beta_(tw)(z)=(z−tw)/(1−t conj(w)z)` moves each zero with velocity
  `w/f'(a_j)` and has first image `−wI`, because `Y_0²=0`.
- If `tau(v)=(Y_1(v))_(00)+(Y_1(v))_(LL)`, direct
  simple-singular-value polarization gives
  `B_J(v,k_w)=(3/4)Re(conj(w)tau(v))`.
- Since L345's source has no zero-motion component,
  `tau(W_phi)=0` exactly.  This proves the endpoint-sum Euler
  equation, independently of phase covariance or the Hessian sign.
- At the L349 stage the lower-flag debt was one cyclic equation,
  `n−2` simple-root tracking equations, and one endpoint-difference
  equation; L350 subsequently closed that endpoint difference.
- Twelve audits through `n=8` give Euler-row residual below
  `9.97e−12`; dataset SHA-256
  `e572460e3f284788bcf5a121e93e49aaa0fb7f0439fa7e2d02c232fde483d1d4`.

## L348/A296 EXACT ROOT-TRACKING FLAG NORMAL FORM (2026-07-26)

- A southwest commutator gauge writes
  `G=U+c e_L e_0*+[A,X]`, with `U` upper triangular and
  `X_(L0)=0`.  Functional calculus gives
  `Df(A)[[A,X]]=[f(A),X]`.
- For the simple spectral projections `P_j` and repeated-zero
  projection `P_0`, exact formulas are
  `(Y_1)_(jj)=f'(b_j)(tr(P_jG)−v_j)+kappa_j(Y_1)_(L0)` and
  `(Y_1)_(00)+(Y_1)_(LL)
   =f'(0)(tr(P_0G)−2v_0)+kappa_0(Y_1)_(L0)`.
- Hence on the cyclic-free face the diagonal debt is exactly
  `n−2` simple-root tracking equations, one mean repeated-root
  tracking equation, and one endpoint-splitting equation.  L348
  does not prove that the L345 dual lift satisfies them.
- Twelve full joint-basis audits through `n=8` give response
  residuals below `1.87e−15`; dataset SHA-256
  `d87d5cc883aefe7ad77ab0a7c4ec81c8e8643d9779dfa2c4160b4c2e3875a04c`.

## L347/A295 EXACT CYCLIC-CORNER/SHAPE IDENTITY (2026-07-26)

- If `sigma_E=omega_phi s_E`, direct coefficient extraction in
  L343's polynomial frame gives
  `hat(sigma_E)(−n)=(Ccal*E Ccal)_(L0)/4`.
- The endpoint columns are `Ccal e_0=gamma p`,
  `Ccal e_L=q`.  The inverse-Riemann correction and every
  zero-motion correction are upper triangular, while the endpoint
  divided difference is `f'(0)`.  Hence
  `(Y_1)_(L0)=4f'(0)hat(sigma_E)(−n)/gamma`.
- The cyclic flag debt is therefore purely physical and is exactly
  the final negative coefficient of the weighted shape polynomial.
  L347 locates the coefficient but does not prove its vanishing on
  `W_phi`.
- Twelve audits through `n=8` give highest-mode and corner residuals
  below `6.36e−15`; dataset SHA-256
  `a4bad7eb40ff2aa41d65c8b7cee515acb4f8ae21b1ea4b0caa3eb353f8646829`.

## L346/A294 EXACT LOWER-FLAG RECURRENCE (2026-07-26)

- Differentiating `T_eY_e=Y_eT_e` gives
  `[S,Y_1]=Y_0C−CY_0`.  Since `Y_0=e_0e_L*`, the right side has no
  strict lower entries.
- The strict lower commutator equation is a southwest recurrence.
  Every coefficient `lambda_i−lambda_j` is nonzero except at the
  two repeated zero endpoints `(L,0)`, so the corner is the only
  free strict-lower entry.  Together with the diagonal, `n+1`
  scalars generate the lower response.
- Therefore A293's observed `tril Y_1(W)=0` is equivalent exactly to
  `diag Y_1(W)=0` plus `Y_1(W)_(L0)=0`.  Prove these moving-root and
  cyclic-closing Euler equations only; do not expand the other
  lower entries.
- Interior-zero collisions require a confluent recurrence and are
  outside this nondegenerate chart.
- Twelve audits through `n=8` give rank `n+1`, commutator residual
  `2.28e−15`, and row-space reconstruction residual `4.92e−9`.
- Dataset SHA-256
  `537bb48cebeb8c79ff3ba84df566f97a6fe01b5efeb1cd56ce81f70ece7476f3`.

## L345/A293 EXACT DUAL CRITERION + NUMERICAL HARDY POLARIZATION (2026-07-26)

- Before eliminating zero velocities and the L343 disk fibre, the
  joint scalar Hessian `J_phi` has exactly `n²` real variables.
  With the extended shape map `S=[S_x;S_y]`, its dual compliance is
  `K=S J_phi^−1 S*`.
- Whenever `J_phi` is nondegenerate, put
  `V_x=J_phi^−1S_x*`, `V_y=J_phi^−1S_y*`, and
  `W=V_y+iV_x`.  Direct block expansion gives
  `W^T J_phi W=K_yy−K_xx+i(K_yx+K_xy)`.  Thus total bilinear
  isotropy is exactly equivalent to phase covariance of `K` and of
  the final Schur form.  This equivalence is L345.
- Numerically `tril Y_1(W)=0`, including the diagonal: the conformal
  positive-frequency source selects a purely upper-Schur-flag first
  response.  The right and left L336 endpoint forms are each
  bilinearly isotropic on this complex range.  These are not real
  one-sided sign claims and remain unproved.
- Fifteen models in dimensions `4..8` survive.  The worst forward
  residual `4.29e−6` occurs at condition number `1.61e11`, with
  backward solve residual `9.05e−17`; well-conditioned models are
  many orders tighter.
- Prove the Euler response's lower-triangular vanishing in L344
  coordinates, then derive endpoint isotropy from L336.
- The soft rank-one dominance of the inverse form produced a false
  universal-diagonal-plus-rank-one ansatz; resolved corrections are
  model-dependent and can be indefinite.
- Dataset SHA-256
  `d02f3734280a2a5594fba11cdcc1514ad6b760c6b3fb112c1b76bf5b733ce8db`.

## A292 EXACT CONFORMAL CHART + NUMERICAL HERMITIAN COLLAPSE (2026-07-26)

- The ordinary support modes
  `2(hat(s)(2),...,hat(s)(n))` give an exact complex coordinate on
  `P_n^R/(omega_phi P_1^R)`.  If the low modes of `sigma` are gauged
  away and these modes of `s=sigma/omega_phi` vanish, then
  `integral sigma²/omega_phi=integral sigma s=0`, so the quotient
  class is zero.
- In these coordinates the live real Schur form commutes numerically
  with multiplication by `i`; its complex-symmetric block vanishes.
  Thus the exact-looking target is one Hermitian `(n−1)`-square,
  not an arbitrary real `(2n−2)`-square.
- This explains but corrects the preliminary paired-spectrum clue:
  the weighted L343 coordinates do not carry the obvious complex
  structure, and their pairs are not exactly equal.
- Fifteen complete models, three per dimension `4..8`, give maximum
  phase residual `3.85e−12`; all Hermitian matrices are negative,
  with softest maximum eigenvalue `−3.51e−7`.  Neither identity is
  promoted beyond the exact coordinate chart.
- Dataset SHA-256
  `41e0bec1e7a37a5ac981a59821191a7ebac44969593187203ffcffa4bc1a0ad6`.

## L344/A291 EXACT FIXED-TOEPLITZ PORT TRANSPORT (2026-07-26)

- L343 supplies two factors of the same support slack:
  `F_zeta Ccal` and `Hhat^(1/2)(I−conj(zeta)R)`.
  Their Grams and ranks agree exactly, so polar decomposition gives
  an isometry between their ranges.
- Therefore L339's minimal response to `r⊥y_zeta` has exactly the
  same norm as
  `Hhat^(dagger/2)(I−zeta R*)^−1 Ccal* r`.
  The only pseudoinverse is fixed, with constant kernel `span{e_L}`; the
  other factor is one elementary nilpotent triangular resolvent.
- For `r=H_Cy−s_Cy`, integrating this norm is exactly `P_phi(C)`.
  Multiplication by the L343 support polynomial gives forcing
  `Ccal*H_C Ccal f_zeta−s_C Kf_zeta`, whose scalar part is the
  finite shape coordinate `omega_phi s_C`.
- The next step is now sharply finite: eliminate the strict
  general-H disk direction in this coefficient system and compare
  the resulting endpoint output with the retained fixed Hardy
  energy.
- Twelve fully polarized audits through dimension eight agree with
  the original port and L339 Gram within `1.44e−13`; dataset
  SHA-256
  `0e800d68f9e233738a51028b72ca9ce728f1990bb55286984e8125bbb187d874`.

## L343/A290 EXACT DISK-CHART RECENTERING AND SHAPE QUOTIENT (2026-07-26)

- With `delta(t)=det(I−tS)`, the L339 null state gives the explicit
  polynomial support frame
  `C=[Lr_(n−1),...,Lr_0]`,
  `r_j=sum_(k<=j)delta_k S^(j−k)q`.
- Congruencing the support pencil by `C` gives exactly
  `K=Hhat+R*Hhat R` and `K A_c=2Hhat R`.  Gau--Wu equality,
  terminal cyclicity, L183's Hardy defect, and L187's residual
  factorization force the leading `H^−1` to be Hermitian Toeplitz.
  Thus every finite nondegenerate model is a global recentering of
  the same full disk chart, not merely one near the Crabb collision.
- General-H disk directions transverse to the inverse-Toeplitz
  equality locus have dimension `(n−2)^2` and strict negative common
  L342 curvature.
- The non-affine weighted support motion
  `Sigma_phi=[omega_phi s] mod omega_phi span{1,cos,sin}` has rank
  `2n−2` and
  kernel exactly that disk block.  The sole live form is its Schur
  maximum.
- Five audits in dimensions `4..8` match the chart within
  `1.10e−14`, inverse Toeplitz within `2.23e−14`, and disk/shape
  kernel projectors within `2.27e−8`.  The open quotient remains
  negative numerically, with softest tested eigenvalue
  `−2.16e−9`; this sign is evidence only.
- Dataset SHA-256
  `01dabed0b2aad2b584210928396051bd74a585cb5832fec49b3a83036ee10103`.

## L342/A289 EXACT SIMILARITY/SCALAR OSCULATION (2026-07-26)

- Every finite Gau--Wu equality model has the same sharp metric
  `P_0=diag(1,2,...,2,4)` and positive embedded dual weight
  `W=diag(0,2,...,2,1)`, with the exact balance
  `W−AWA*=e_L e_L*−4e_0 e_0*`.
- The full second-order metric SDP reduces exactly to a minimization
  over one complex boundary row.  Its homogeneous Hessian is
  `8||x_mid||²+(8/3)|x_L|²`, independent of the model zeros and
  dimension.
- The difference from four times the joint scalar Hessian is exactly
  two positive endpoint-matching squares.  Strictness of the pure
  metric and pure zero blocks makes the square residual map
  bijective, proving exact optimized osculation:
  `e_phi(C)=4 max_h J_C(h)`.
- This closes the strategic question of whether the completely
  bounded route has a distinct quadratic obstruction.  It does not;
  the sign of the common form remains open for `n>=4`.
- Twelve audits through dimension eight reproduce the full gap factor
  within `1.27e−13`; dataset SHA-256
  `30f38471eea4fb3728fe3a10b154053fcc44cbbc9c13b5673f125a27fb6a8ae1`.

## L341/A288 EXACT PAIRED FRAMES AND NUMERICAL OBSERVATION LAW (2026-07-26)

- For any moving orthonormal basis of `K_(f_e)`, the positive primal
  and dual L127 frames obey exact Stein equations whose right sides
  are `qq*−Y_eqq*Y_e*` and `pp*−Y_e*pp*Y_e`.
  Their traces are exactly L336's two endpoint defects.
- At the Gau--Wu base the frames are `I−pp*` and `I−qq*`.
  Differentiating gives stable Stein equations forced by the
  normalized operator direction and the endpoint vectors
  `u=Y_1q`, `v=Y_1*p`.  Their moving rank-`n−1` kernels give exact
  second-order squares.
- After L338's zero completion define
  `O_phi(E)=(s_E,Y_1^opt q,(Y_1^opt)*p)`.  In twelve complete models
  through dimension nine, `(u,v)` has rank `4n−8`, the support
  coordinate adds `2n−6`, and
  `ker O_phi=ker G_phi` with total rank `6n−14` and nullity
  `(n−4)²`.  The full optimized first image vanishes on the common
  kernel.
- The rank/kernel identification is numerical only.  The live task
  is to derive the response metric `M_phi` from the exact paired
  frames and prove its L339 port domination.
- Dataset SHA-256
  `73b7d48403d8321ec89068ae6f492b5517c81153cf8b5383d37c99c32dda8e56`.

## L340/A287 NUMERICAL BOUNDARY-RANK LOSS LAW (2026-07-26)

- After L338's exact zero completion, set
  `H_phys=Q(C,0)−2||kappa_C||²` and
  `G_phi=2I+2P_phi−H_phys`.
- Twelve complete generic models, two per dimension `4..9`, give
  `G_phi>=0` to roundoff, rank exactly `6n−14`, and nullity exactly
  `(n−4)²`.  The identity
  `(n−1)²+1=(6n−14)+(n−4)²` accounts for the full physical normal
  quotient.
- This is a numerical structural diagnostic only.  It suggests that
  the quadratic-size quotient has an automatic bulk and a
  boundary-size response.  It does not prove either positivity or
  the rank formula.
- The live proof target is to derive
  `G_phi=L_phi*L_phi` from the paired L127 Stein identities and
  L339's unitary colligation, then prove
  `||L_phi C||²<=2||C||²+2||gamma_C||_(L²)²`.
- Dataset SHA-256
  `5611d1090336303e28f7a423e91114ad53174e3350c62f4350acfa9f127e5572`.

## L339/A286 UNWEIGHTED ANDO SUPPORT-PORT ENERGY (2026-07-26)

- The defect-one completion `U=S+qp*` is unitary.  With
  `D=0 direct-sum I/sqrt(2) direct-sum 1` and
  `E=1 direct-sum I/sqrt(2) direct-sum 0`,
  `A=2EUD` and the pencil `F_zeta=UD−zeta E` satisfies
  `F_zeta*F_zeta=I−Re(conj(zeta)A)`.
- The null state is
  `y_zeta=(1 direct-sum sqrt(2)I direct-sum 1)
  (I−conj(zeta)S)^−1q`, and
  `||y_zeta||²=2D_f(zeta)`.
- Therefore the minimal support response
  `gamma_C=(F_zeta*)^dagger H_C y_zeta` obeys
  `||gamma_C||²=2D_f t_C`, and L335 becomes the unweighted identity
  `P_phi(C)=integral ||gamma_C||²`.
- The live gate is now to identify L338's remaining physical loss
  and `kappa_C` as outputs of this same unitary colligation and prove
  their energy is bounded by the physical input plus the support
  port.  The contraction is not yet proved.
- Twelve polarized support-port audits pass within `1.60e−14`;
  dataset SHA-256
  `cae577fc8c856279b734f6e6a7481d0eefc42f531af5ccea80e4863a2371a2d3`.

## L338/A285 EXACT MIXED HARDY COMPLETION (2026-07-26)

- Modulo phase, every degree-preserving finite-Blaschke tangent is uniquely
  `h_k=(zf)k−J_f k`, `k in K_f`.  If the zeros `a_j` move with
  velocities `v_j`, then
  `k=sum_j conj(v_j)/(1−conj(a_j)z)`.
- Model calculus gives `h_k(S)q=−J_f k` and
  `h_k(S)*p=−zk`, so the pure zero block is the exact Szegő Gram
  `2||k||²`.
- The mixed physical/zero term is one real functional of `k`.  Its
  Riesz vector `kappa_C` gives the exact completion
  `r_2+l_2=2||k−kappa_C||²+Q(C,0)−2||kappa_C||²`.
  Thus the optimizing zero tangent is `k=kappa_C` and every root
  velocity has been eliminated without a numerical Schur matrix.
- The sole remaining gate is the physical inequality
  `Q(C,0)>=2||kappa_C||²`.  Insert L335's support Gram and express
  the Riesz response through the contraction/characteristic
  residuals.
- Twelve block-Fréchet audits pass within `1.21e−13`; dataset
  SHA-256
  `e600a466d6ad1fd78e692f9c370a8a73b6e6fe2101e582f98e232991de83d67c`.

## L337/A284 EXACT ZERO-MOTION PROJECTION SQUARE (2026-07-26)

- At fixed `S=S(zf)`, for every inner `g`,
  `1−||g(S)q||²=1−||g(S)*p||²=||P_(zf H²)g||²`.
  The equality of the endpoints is the canonical conjugation
  `J_(zf)u=zf zbar ubar=f ubar`, not a numerical symmetry.
- Therefore an inner tangent `g_e=f+e h+...` contributes exactly
  `r_2=l_2=||P_(zf H²)h||²`.  The zero acceleration drops out and
  the complete pure zero block is a Hardy projection square.
- The historical next step was to complete the mixed square.  L338
  now does so exactly; L337 alone did not control the mixed form.
- Twelve independent finite-model/Fourier checks pass within
  `1.43e−15`; dataset SHA-256
  `233e72245d54a1d71ad1442332d0d26d9dbc5bff207d43c2e42f35390cad8b4f`.

## L336/A283 TWO-SIDED ENDPOINT-DEFECT FLUX (2026-07-26)

- Exact singular-value algebra at `XY_0X^−1=2pq*` gives
  `J=−(r_2+l_2)/2−(|Y_1(qq)|²+|Y_1(pp)|²)/4`.
  This is dimension-free and uses the complete moving-zero/Riemann
  jet, not a numerical fit.
- L127 and its left-handed companion turn `r_2+l_2` into a symmetric
  scalar flux of the initial and final contraction defects through
  the finite model space of `f`.
- Cheap falsification is decisive: right and left fluxes separately
  are indefinite.  Every standard model in dimensions `4..8` has
  two negative modes at each endpoint; the holdouts confirm the
  failure and generically the same index.  Their sum is positive
  definite to numerical precision in 12 standard and 42 independent
  holdout models.  The softest holdout margin is `3.43e−11`, positive
  but deliberately treated as numerical.
- The linear residual pair
  `R=S*C+C*S` and
  `Z=D(zf)(S)[C]+(z dot f)(S)` has combined rank `n²` on every
  tested joint quotient; the individual ranks overlap by exactly
  two.  This explains both the full Hessian rank and the two missing
  directions in the raw one-defect shortcut, but the rank law is not
  yet proved.
- Standard dataset SHA-256:
  `023ef06315c2524d1dcc0e818698c6ddc1ab1fea4562b3272955b93b5f802cbb`.

## L335/A282 SECOND-SUPPORT NEGATIVE GRAM (2026-07-26)

- For `A=XS(zf)X^−1`, the endpoint model vectors satisfy the exact
  transition `<p,h(A)q>=2<f,h>_(H²)`.
- If `t_E=<H_Ev,R H_Ev>` is the polarized second support
  coefficient and `K_t=wH_t`, the inverse-Riemann term `−K_t(A)`
  contributes exactly
  `−2 Re<f,f'K_t>=−2 integral D_f t_E dm`.
- The weight
  `D_f(zeta)=zeta f'(zeta)/f(zeta)` is the strictly positive sum of
  the Poisson kernels of all zeros of `f`.  Thus this entire term is
  a dimension-free negative Gram, not numerical evidence.
- The remaining sign is precisely
  `max_u H_rest(E,u)<=P_phi(E)`.  In every L334 sample the optimized
  remainder has positive index `2(n−2)`; all other directions are
  already favorable.  The inertia pattern is numerical and must be
  derived before use.
- The Ando support factorization remains useful, but its raw
  rank-one partial-isometry defect was explicitly kernel-tested and
  misses two Hessian directions.  Do not use it as the complete
  residual.
- Twelve full-normal audits in dimensions `3..8` agree with the
  exact transition/Schwarz formula within `5.33e−15`; dataset
  SHA-256:
  `cb068337f93287075e54239265d493e73108305bd9f98082354b01de16663e93`.

## L334/A281 FINITE GAU--WU HESSIAN FALSIFICATION (2026-07-26)

- For ordered `phi`-zeros `(0,b_1,...,b_(n−2),0)`, the explicit
  Gau--Wu matrix generator verifies `W(A)=D` and `f(A)=2E_(1n)` to
  machine precision before any Hessian is accepted.
- The generic affine/unitary/all-zero tangent has real rank
  `n²+2n−2`; its full orthogonal normal quotient has dimension
  `(n−1)²+1`.
- A general analytic jet samples the simple support projection and
  reduced resolvent, forms the two inverse-Riemann coefficients by
  Schwarz transforms, moves all `n−1` Blaschke zeros, and takes the
  finite-dimensional zero-motion Schur complement.
- At a complex `3 x 3` model the general code agrees with L332's
  closed exact spectrum within `1.76e−13`.
- Ten complete normal Hessians, two in every dimension `4..8`, have
  no positive eigenvalue.  The softest value is `−1.6048e−8` in an
  `n=8` near-root-collision sample and is stable from 256 through
  2048 Fourier nodes; it remains conservatively labelled numerical,
  not proved.
- A soft `n=5` jet direction agrees with independent true-map
  Theodorsen normalization and nonlinear optimization of all zeros
  across three finite-difference steps.
- Complete-jet dataset SHA-256:
  `0caa8336d1ae210d4f2df3d6489e6acd77628834c75d2f3e90faaebd5cb79936`.
  Independent nonlinear dataset SHA-256:
  `acc979a2172cafc8832bd81d271a1b14bc0cb630e279d315bca4dbdd400d8435`.
- This closes the requested cheap falsification only.  It does not
  prove the arbitrary-degree Hessian or a local theorem beyond L333.
  The next ordered task is the defect-one/model-space negative Gram
  and quotient injectivity stated in the canonical gate above.

## L332--L333/A280 STRICT GAU--WU HESSIAN AND LOCAL THEOREM (2026-07-26)

- The exact affine-unitary Gau--Wu equality tangent at a fixed
  nonzero `G_a` has real dimension `13`, not the preliminary count
  `12`.  A closed `18 x 18` determinant gives the five bottom-row
  normal coordinates.
- The support spectrum `1,-1,a cos(theta)` yields exact first and
  second Riemann jets.  Rationalizing `a=2q/(1+q²)` reduces every
  Schwarz mean to residues at `0,q,a`.
- Because `G_a` is nonnilpotent, the Frechet term must use the
  degree-five Hermite calculus of the doubled minimal polynomial
  `z²(z-a)`.  Reusing the old nilpotent power truncation gives a
  false Hessian and was rejected before banking.
- The four real velocities of the two Blaschke zeros have a strict
  negative Hessian.  Its Schur complement is a sparse five-normal
  form with real `3 x 3` and imaginary `2 x 2` blocks.
- All five Sylvester minors of the negative form factor into positive
  rational denominators and four polynomials in `q²`.  Their complete
  Bernstein coefficient lists are strictly positive on `[0,1]`.
  Thus the optimized transverse Hessian is negative definite for
  every fixed `0<|a|<1`.
- L331 supplies stationarity along the equality manifold.  Gau--Wu
  irreducibility makes `{0,a}` the unique equality zero pair at the
  base, while Crouzeix's finite-extremal theorem separates all other
  degree-at-most-two branches.  The strict Hessian therefore proves
  the scalar constant-two inequality in a full operator neighbourhood
  of each fixed nonzero `G_a` and every affine-unitary image.
- This is local, scalar, and order three.  It degenerates at the
  Crabb collision `a->0` and support collision `|a|->1`; it is not a
  global `3 x 3` or completely bounded theorem.
- Five exact rational specializations regenerate after the free-symbol
  residue/Hermite derivation; dataset SHA-256:
  `8d75b8541f19aba118bc73840b2f0c428058b693f582b7c8ab42a017b2999752`.
  An independent Theodorsen/optimized-Blaschke finite difference
  matches the exact normal Hessian within `2.93e-5`; dataset SHA-256:
  `7ea05baeef3557d64761b1e9d731f21a2dddcc670d35ee747f46c0a36c2d85b7`.

## L331/A279 AMBIENT STATIONARITY AT `G_a` (2026-07-26)

- For an arbitrary complex direction `E`, let `s_E` be the first
  support variation and `H_s` its analytic Schwarz transform.  L101
  gives normalized direction `Gcal_E=E−G_aH_s(G_a)`.
- Differentiating `f_a(X)=X(X−aI)(I−aX)^−1` gives the exact scalar
  functional
  `L_a(E)=r_aE_12−aE_13+2aE_22+r_aE_23`.
- The minimal polynomial gives
  `L_a(G_a)=4` and `L_a(G_a^(k+1))=2a^k`.
- The top support projection is an explicit quadratic polynomial in
  `Re(e^(−i theta)G_a)`.  After rationalizing
  `a=2q/(1+q²)`, the Poisson-weighted support mean has only the
  interior poles `0,q,a`; their nine entrywise residues reproduce
  exactly the Frechet-gradient matrix.
- Therefore `Re L_a(Gcal_E)=0`, so the squared sharp scalar norm has
  zero first derivative in every ambient direction.
- This proved stationarity only at the L331 checkpoint.  The then-live
  optimized transverse Hessian was closed by L332--L333.
- Ninety-five exact direction audits and the free-symbol residue
  identity pass; dataset SHA-256:
  `23630cbbf378a5a6de520b15c67db5de7c1d822b3a7c890c81272d433a154293`.

## L330/A278 NONCRABB GAU--WU SHARP STRATUM (2026-07-26)

- For real `0<|a|<1`, let
  `G_a=[[0,sqrt(2(1−a²)),−2a],
        [0,a,sqrt(2(1−a²))],[0,0,0]]` and
  `f_a(z)=z(z−a)/(1−az)`.
- The support characteristic polynomial is exactly
  `(lambda−1)(lambda+1)(lambda−a cos(theta))`.
  Therefore `W(G_a)=D`.
- Exact rational functional calculus gives `f_a(G_a)=2E_13`, so
  `G_a` is a scalar equality matrix.
- Its spectrum is `{0,0,a}`.  For `a!=0` it cannot be an affine image
  of the single-eigenvalue Crabb block.  This is the smallest exact
  certificate that L192/L329 do not cover all known equality points.
- The construction is Gau--Wu Corollary 3 specialized to zeros
  `0,a`; its equality content is classical.  The campaign result is
  the corrected scope and ordered frontier.
- Six rational parameter audits and the free-symbol identities pass;
  dataset SHA-256:
  `185c75e1ed68ad4634a366472d0130e8d2ff43143d657d41479d1470b8c1d31a`.

## L329/A277 FIXED REPEATED-BLOCK SCALAR NEIGHBOURHOOD (2026-07-26)

- For every fixed `p=L+1>=3` and multiplicity `m`, there is a full
  operator neighbourhood of `C_p tensor I_m` on which the scalar
  Crouzeix inequality holds.
- Crouzeix's 2007 finite extremal theorem reduces scalar norming data
  to Blaschke products of degree at most `pm−1`, so the direction,
  function, support, transfer-rank, and finite-disk-flag data admit
  compact ramified strata.
- L61/L199 remove positive-Jensen and nonscalar normal branches.
  L325 supplies the same-direction linear joint disk/channel defect;
  L324 stacks every scalar normal/reflected response; L327 proves its
  square coordinates are quantitatively equivalent to that exact
  defect.  One vector completion therefore costs only the square of
  a retained linear reserve.
- L317/L318 retain the grouped elliptic direct Grams and exact
  analytic tail.  The final associated-graded audit keeps the
  principal elliptic block inside that exact completion and puts only
  marked responses through L324, so no metric or curvature is spent
  twice.
- L197 makes the disk descent finite.  At exact equality, Gau--Wu
  (L326) splits a scalar disk model and L328 proves it removes exactly
  one `p`-dimensional copy; induction terminates at L192.
- The theorem is scalar, local, and fixed-dimensional.  It does not
  classify all global sharp pairs and does not prove the completely
  bounded conjecture.
- Forty-eight exact master-square audits pass; dataset SHA-256:
  `3267bb8d19a30edb18abb343b0a5ddca1d83fd90f7ce8ef5017d580e440e814d`.

## L328/A276 ONE-COPY EQUALITY-MODEL RIGIDITY (2026-07-26)

- L326 initially permits a Gau--Wu model dimension `(L+1)r`.
  L328 proves that `r=1` in a sufficiently small repeated-Crabb
  neighbourhood.
- In a hypothetical `r>=2` sequence, every zero of the finite
  Blaschke model tends to zero because it is an eigenvalue of a
  reducing summand of a matrix tending to the nilpotent apex.
- Gau--Wu Corollary 3's explicit triangular matrix then tends to the
  one long Crabb chain `C_((L+1)r)`, whereas reducing-projection
  compactness makes the same summand tend to
  `C_(L+1) tensor I_r`.
- Their `L`th powers have squared norms `2` and `4`, respectively,
  so they cannot be unitarily equivalent.  Thus the equality model
  is one nearby scalar disk block and its complement has multiplicity
  `m−1`.
- This does not restore the false claim that the nearby model is
  nilpotent or monomial; only its dimension is rigid.
- Ten exact power-separation audits pass; dataset SHA-256:
  `8182c3f7bd1564b1ab4c43115d33d8b854b8a3f02e6d7e93c19d542c4e94133d`.

## L327/A275 SAME-DIRECTION JOINT DEFECT COMPARISON (2026-07-26)

- Retain L325's pre-minimization defect on the scalar norming
  direction:
  `Delta_k=<D_Xk,k>+1−sigma_X(k)`.  Paying a response by
  `min_k Delta_k` from another direction would be invalid.
- Observability Parseval identifies `1−sigma_X(k)` with squared
  distance of `O_C k` to the rank-one product-line set.  At L324's
  terminal equality retraction, matrix-inner Parseval identifies
  `delta` with distance to the same set.
- Undoing L197's finite least-squares flag gives
  `O_C k−O_eq(Wu)=sum_j G_j eta_j`, with bounded analytic graph
  columns, while `<D_Xk,k> asymp sum_j||eta_j||² asymp d`.
- The one-Lipschitz distance inequality and Young give
  `Delta_k asymp delta+d`.  Hence L324's
  `O((delta+d)²)` stacked completion is a higher-order cost paid by a
  retained fraction of L325's linear directionwise reserve.
- The statement is curvewise through finite ramified rank strata and
  uses neither a global maximizing frame nor a pseudoinverse.
- Twenty-eight exact distance/Young audits pass; dataset SHA-256:
  `8eea0b6f8a2a721155634d02a47f1335f5b075307d67b09c66bd7c44338c74ba`.

## L326/A274 CLASSICAL TERMINAL SCALAR MODEL SPLIT (2026-07-26)

- Gau--Wu, Linear Algebra Appl. 430 (2009), Theorems 5 and 8,
  already prove the missing equality structure: if a numerical
  contraction satisfies `||f(T)x||=2`, then `f` is inner with
  `f(0)=0` and `T` has an orthogonal reducing summand
  `X_phi S(zf) X_phi^−1`.
- That summand is cyclic, irreducible, has no unitary part, has
  numerical range the disk, and itself attains two.  This is a
  classical theorem, not an original campaign breakthrough.
- A sufficiently small repeated-Crabb neighbourhood has no reducing
  unitary summand.  Compactness of reducing projections and
  irreducibility of one Crabb block show that the Gau--Wu summand has
  dimension `(L+1)r`; its complement is a numerical contraction near
  multiplicity `m−r`.
- L328 sharpens this to `r=1` by comparing the explicit Gau--Wu
  zero-collision limit with the incompatible long-chain power norm.
- The first specialization attempt incorrectly assumed every L193
  disk point was nilpotent.  The proof was corrected before banking:
  general equality retains the finite Blaschke model.  Only an
  explicitly nilpotent point forces the monomial `z^L` and one
  literal Crabb block.
- Thus L325's terminal partial-endpoint zero face is closed, and the
  remaining work is the quantitative finite merger rather than a new
  model-kernel proof.
- Seven exact nilpotent-model audits pass; dataset SHA-256:
  `6fb00450f8dd2222d5bf9c0d3466d4832fda40d047e82e984da84617af4dd376`.

## L325/A273 JOINT ACTUAL-DISK SCALAR GAP (2026-07-26)

- A first attempted final assembly exposed another invalid shortcut:
  L322's channel reserve at an equality anchor cannot simply be added
  to an independently chosen L197 disk compression.
- For the actual canonical disk metric, use its top spectral cluster
  `X`, disk defect `D_X=X*(4P^−1−I)X`, and scalar Hardy output-line
  score `sigma_X(k)` on the same unit direction.
- The exact three-loss identity retains input upper loss, contraction
  loss, and output lower loss.  If the scalar row angle is `t`,
  contraction forces its complementary row below `sqrt(1−t²)`.
- The resulting `2x2` form has trace `q+mu/2` and determinant
  `q*mu*(1−t²)/2`.  The other half of the complement gap pays the
  disk defect.  Thus every scalar Schur function has squared norm at
  most `4−c*Delta`, where
  `Delta=min_k{<D_Xk,k>+1−sigma_X(k)}`.
- `D_X` is analytically congruent to the complete L197 endpoint, so
  this is one joint reserve through all disk rank changes.  On L193's
  equality manifold it reduces to the older channel leakage.
- At a terminal partial endpoint, scalar equality forces
  `Pk=4k`, `V*C^n k=b_n v`, and `f(C)k=Vv`.  L326 now applies
  Gau--Wu's classical equality theorem to split the scalar disk-model
  summand even away from the full L193 equality manifold.
- Two hundred twenty-five exact angle/allocation records pass; dataset
  SHA-256:
  `bd7e63c2ffaa8ba1bdbe5e6c61cfae1efbb722c107e36eef92b9dfb07cd0c967`.

## L324/A272 JOINT CHANNEL/DISK RESPONSE FACTOR (2026-07-26)

- The first post-L323 adversarial audit found a genuine scope issue:
  transfer leakage `delta` alone does not measure off-channel disk
  residuals inside a later L197 response.  L323 remains correct for
  transfer-only responses; its broader joint-response application has
  been narrowed in the ledgers.
- Along a ramified analytic zero-Jensen arc, transport the selected
  scalar state through every L197 congruence.  The accumulated active
  Hardy quotient columns `eta_j` satisfy
  `d asymp sum_j ||eta_j||²`, where `d` is the whole scalar disk loss.
- One complement-frame sign change flips both L323's transfer cross
  blocks and every off-channel disk quotient block.  After subtracting
  L192's channel-preserving scalar response, physical invariance removes
  total cross degrees zero and one.
- Hence `|R_full−R_split|<=C(delta+d)`, including mixed
  `sqrt(delta*d)` terms.  Completing all finite scalar responses at
  once costs `O((delta+d)²)`.  L327 now proves
  `delta+d asymp Delta_k` on the actual L325 norming direction, so a
  retained fraction of L325's linear reserve absorbs that cost after
  shrinking.
- This is the endpoint-specific scalar information missing from
  L319's arbitrary Schur germs.  It does not assert positivity of the
  complete matrix endpoint.
- Fifteen exact invariant-monomial and completion audits pass; dataset
  SHA-256:
  `c06279887606fcff4954bdc059b8b7fd8a118fbdff4ea6b4e25c7d795eb35472`.

## L323/A271 TWO-SIDED CHANNEL DEFECT AND RESPONSE FACTOR (2026-07-26)

- In channel frames `B=[[b,r],[c,D]]`, matrix-inner boundary
  unitarity gives exactly
  `||r||_H2²=||c||_H2²=delta`, where
  `delta=1−sum|u*B_nv|²`.
- Pointwise, `1−|b|²=||r||²=||c||²` is nonnegative.  Every Fourier
  coefficient of this defect and every bounded cross-Hardy pairing
  is therefore at most `delta`.
- Flipping both complement defect frames sends `(r,c)->(−r,−c)` but
  changes no physical scalar response.  The finite matrix-Schur
  chart retracts analytically to the block-diagonal inner locus with
  only quadratic diagonal error.
- Hence every transfer-only scalar response which vanishes on the
  split-inner subchart is `O(delta)`, and its normal completion costs
  `O(delta²)`.
- The original downstream sentence also applied that bound to joint
  disk/normal responses.  That scope was too broad because
  off-channel disk residuals are not measured by `delta`.  L324
  corrects it with the accumulated disk loss `d`.
- Use L324, not this section alone, in the finite curve-selected flag
  across L197, L321, and L318.
- Nine deterministic exact matrix-inner path audits pass; dataset
  SHA-256:
  `beca54b375f93da923b105d4a5a1070d1c51a0da6e9f7bba8cee526f235474e2`.

## L322/A270 SHARP LINEAR SCALAR-CHANNEL GAP (2026-07-26)

- For `P=2I−VV*+2WW*`, every contraction `F` obeys the sharp
  endpoint-angle bound
  `||P^(−1/2)FP^(1/2)||²
   <=(4+s²+s sqrt(s²+8))/2`,
  where `s=||V*FW||`.
- The proof reduces the endpoint geometry to a `2x2` matrix with
  trace `4+s²` and determinant `4`.  A unitary rotation through the
  middle eigenspace of `P` attains equality for every `s`.
- L320's Hardy corner estimate `s²<=sigma(B)` therefore gives
  `||f(T)||²<=(4+sigma+sqrt(sigma²+8sigma))/2
   <=4−(4/3)(1−sigma)`.
- Combined with L321, the certified reserve begins at
  `(4/3)Lambda(D)s^(2q)`, quadratic in channel-breaking amplitude.
  L320's older quartic reserve was valid but non-sharp.
- L323 proves transfer-only scalar transport is `O(1−sigma)`.
  L324 adds the accumulated disk loss and proves the actual joint
  completion costs `O((delta+d)²)`.  Retain the linear reserve for
  that finite flag assembly.
- Thirteen sharp/random/actual-transfer audit batches pass; dataset
  SHA-256:
  `52c16681fcddfcb21c08f4483093ea2f9e8a175e4dc8f22d6dfa2e90a46835b2`.

## L321/A269 LEADING SCALAR-CHANNEL VALUATION (2026-07-26)

- Along an analytic matrix-inner curve through `B_0(z)=z^L I_m`,
  suppose the first transfer jet has order `q`, with off-monomial
  coefficients `D_n`.
- The fixed-input product score is the top eigenvalue of the trace-one
  covariance `sum B_n vv* B_n*`.  Uniform leading-eigenvalue
  separation gives
  `1−sigma(B_s)=s^(2q)Lambda(D)+O(s^(2q+1))`, where
  `Lambda=min_v sum_(n!=L)||(I−vv*)D_nv||²`.
- `Lambda=0` exactly when the off-monomial jets have a common
  eigenvector.  The monomial coefficient jet only rotates the common
  output line and contributes no first leakage.
- L320's older estimate alone supplied only
  `(Lambda²/4)s^(4q)`, but L322 proves the sharp certified reserve is
  at least `(4/3)Lambda s^(2q)`.
- On `Lambda=0`, carry the common eigenline to the next nonzero jet;
  if it persists exactly, L205 splits it and L192 handles the scalar
  block.
- Four deterministic exact matrix-Schur path audits pass; dataset
  SHA-256:
  `72b24fe92a8c2ec97148015cc540fba0f3a7bd7df3455bd8c886c4576c2fc035`.

## L320/A268 QUANTITATIVE SCALAR-CHANNEL RESERVE (2026-07-26)

- This estimate remains valid, but L322 quantitatively supersedes it.
  Use L320 only for its Hardy corner bound
  `||V*f(C)W||²<=sigma(B)`.
- For `P=2I−VV*+2WW*` and any contraction `F`, the exact input,
  contraction, and output loss decomposition gives
  `||P^(−1/2)FP^(1/2)||²<=4−(1−||V*FW||)²`.
- For `F=f(C)`, scalar Hardy Parseval and L201's transfer coefficients
  bound the endpoint corner by `sqrt(sigma(B))`, where
  `sigma(B)=sup_(u,v)sum|u*B_nv|²`.
- Hence
  `||f(T)||²<=4−(1−sqrt(sigma(B)))²`.
  L205 identifies `sigma=1` exactly with a constant scalar inner
  channel.
- Every compact channel-free part of L193's noncommuting equality
  manifold therefore has a uniform strict scalar tube.  Only the
  scalar-channel strata remain sharp.  Exact channels split off as
  single-copy full-Hardy blocks by L205.
- The live proof uses L322's linear reserve in the finite channel
  induction; do not revert to this section's quartic estimate.
- Eleven deterministic audit records pass; dataset SHA-256:
  `973f669a8ae58492c65b73492fbb8b3c5651764ba6aefa0db94adb323c05c40b`.

## L319/A267 SCHUR-ONLY LATER-FLAG OBSTRUCTION (2026-07-26)

- The proposed inference “the old response is the cross block, so
  the exact Schur square leaves a fresh raw L199 face” is false
  without an identity for the actual endpoint's kernel block.
- The exact germs
  `[[s²,s³],[s³,s⁴+s⁶]]` and
  `[[s²,s³],[s³,s⁴−s⁶]]` share the positive first face
  `diag(1,0)`, its reducing kernel, the one-order cross gain, and the
  transported square `s⁴`, but their Schur quotients are `+s⁶` and
  `−s⁶`.
- This does not challenge L197's disk flag, L199's first raw normal
  absorption, L318's elliptic metric, or Crouzeix.  It prevents
  stacking the merger on an unproved assertion.
- The parked complete-similarity fallback is exact: derive L194's
  complete later normal critical graph/Gram after disk
  orthogonalization, including its kernel-block term.  L322's scalar
  channel route was the canonical live gate at this historical
  checkpoint.
- The exact checker passes and regenerates with SHA-256
  `d272005024eceb7d34da8212871a3af27c6778dc01629bba48c85f9897971952`.

## L318/A266 FIXED-HALF ELLIPTIC FLAG (2026-07-26)

- In physical upper-gap orientation, L227/L283 supplies direct
  coefficient `12` in every grade and L298/L303 adds `4theta`.
  At `theta=1/2`, the zero-slack middle block is therefore `14I`.
- L317 gives
  `U_<=2L=Bcal(c){14I+cR_L(c)}Bcal(c)*` after the bounded
  response is canceled.  Every affine leftover and the complete
  L307 square are retained inside the bounded Hermitian `R_L`.
- On a compact repeated-block neighbourhood, `||R_L||<=M_L`.
  Hence `14I+cR_L>=7I` for `c<=7/M_L`, uniformly through transfer-
  rank changes.
- L201 keeps terminal `B_L` uniformly invertible.  Its
  `7c^(2L)B_LB_L*` margin absorbs L289's uniform analytic
  `O(c^(2L+1))` tail.  L194 supplies `P>=I` and
  `P-T*PT>=0` exactly for the selected free row.
- Therefore `I<=P<=4I` and `T*PT<=P` on every fixed repeated
  elliptic chart for sufficiently small `c`.  This is not yet a
  full operator neighbourhood: transverse disk/circular-normal
  residuals remain to be merged through L197.
- Eighteen endpoint-direction, noncommuting Weyl, and terminal-tail
  audits pass.  Dataset SHA-256:
  `2938008cb4b25fc86123940945418f14ea726365899dbecd633affe75a4243fe`.

## L317/A265 SHARP TWO-CHANNEL SUPPORT (2026-07-26)

- L313's wandering-chain factorization now starts at an arbitrary
  paid root height: a word with `r` physical reverse edges next to
  `B_k` exposes only `B_a`, `a<=k+r`.
- Stable Stein terms are not factored termwise.  L243/L245/L251 and
  L255/L272 first close them into a Toeplitz leakage return; L266
  kills every nonconstant shortest return, L261 keeps the diagonal,
  and L277--L283 include the exceptional deep boundary divergence.
- Therefore every retained pair at elliptic degree `d` satisfies
  `j+k<=d`.  Equality forces `j=k=d/2` and is exactly the direct
  L283/L298 Gram.  Every other block has one extra power of `c`.
- Through a fixed terminal jet this gives
  `U=Bcal(c){D_dir+cR(c)}Bcal(c)*+M_S(C_rsp)`, with
  bounded polynomial/rank-stable `R` and L311-normalized response.
  The statement is deliberately grouped: L217's raw order-four
  `B_5` frame term can form an apparent `B_1/B_5` equality cross
  before physical closure.
- Six noncommuting audits reconstruct all 6,120 marked words through
  length seven and separately test continuant-port valuation,
  noncommuting Potapov leakage-shift annihilation, and Rees
  regrouping.  Dataset SHA-256:
  `f3f24239f45beb756b4da21fd5ce7d7ce6fc25e662cefbfd344c395fde52ed51`.

## CONSOLIDATION REVIEW (2026-07-26, after L316/A264)

- All 75 lemma status strings containing `OPEN` were adjudicated.
  Fifty-four stale debts now name the later lemma that superseded
  them, L297's direct effective-face comparison is obsolete after the
  L303--L307 partial-scale route, and 20 genuinely live rows point
  explicitly to ordered next action 1, 2, 3, or 4.  There are now no
  unclassified `OPEN` strings in the lemma status column.
- Every exactly disproved row has its proof-file certificate and an
  explicit downstream guard.  In particular, L297 corrects L295:
  the `+4I` apex gap is a raw-versus-lower-tight normalization change,
  not a branch obstruction.  L19 was downgraded to numerical
  counterevidence because it has no exact certificate and is unused.
- The canonical gate is the one stated above.  All other former
  “live gate/frontier” phrases in this file are marked as historical
  at their dated checkpoint.
- A same-file repeated-block scan of all 253 proof notes found no
  remaining regenerated duplicate derivation, conflict marker, or
  duplicate exact section heading.  Numeric lemma IDs are exactly the
  312 intended handles through L316, and approach IDs are exactly
  A1--A264.
- The tracked endpoint-word dataset still has SHA-256
  `5b5dac6046295c1176d2878365e1627961603282b368bb5e8c14fcd20a56becf`.
  Repository-wide Ruff now passes all 303 experiment modules after
  formatting old scripts and removing only unused names/imports;
  all 303 modules also pass `py_compile`.
- The L290--L316 proof notes and lemma ledger now carry the explicit
  fixed finite repeated-Crabb-neighbourhood scope disclaimer.  A259
  remains conjectural and unused.  The two-agent incident was not
  re-audited in this pass.

## SINGLE-AGENT RECONCILIATION (rechecked 2026-07-26 after L312/A260): concurrent work deconflicted
- Two agents unintentionally shared this repository while commits
  `87c5107`, `4c7a152`, `2f82d15`, and `dcab861` were made.  The Git
  history stayed linear, each commit recorded an explicit coherent
  file set, and no overwritten proof, code, or unresolved merge was
  found.
- L248 was independently regenerated after the overlap.  All six
  audits passed and the result was byte-identical to the tracked
  dataset with SHA-256
  `cf49efe0025b9336a84068e3c92802b5ec23a6bb40d34d4d15028771c7b5ecaa`.
- The only incoming uncommitted artifacts were the L249 note and
  checker.  The interrupted conditional-expression artifact was
  removed; the construction, Schur orientation, exact residual, and
  uniqueness at `lambda=1` were independently recomputed.  The clean
  checker passes and regenerates SHA-256
  `114ba20c89b0d3e1e69db45ed510836a6306caca0cf1b02f0d89aa8f911ec642`.
- The later stopped-agent handoff was reconciled against the linear
  history through `646ade8`.  Its proved L248 contribution was already
  present, its preliminary L249 warning had already received the
  independent cleanup above, and no competing branch, staged file, or
  unresolved working-tree artifact remains.  The sole new in-progress
  result was independently tightened and is banked below as L272/A219.
- At the next single-agent checkpoint the history and working tree
  were reconciled again through L281/A228.  The subsequent
  L273--L281 work is one coherent primary-agent line; no stopped-agent
  artifact, alternate branch, staged overlap, or competing next-action
  instruction remains.
- The 2026-07-26 L303 checkpoint repeated the concurrency audit after
  the user confirmed that the overlap had lasted roughly two hours.
  All four handoff commits (`87c5107`, `4c7a152`, `2f82d15`,
  `dcab861`) are ancestors of the single `master` branch; there is one
  worktree, no staged overlap, and no alternate branch.  L248 and
  L249 still have their independently regenerated hashes
  `cf49efe...b5ecaa` and `114ba20c...ec642`, respectively, and the
  interrupted `if False` artifact is absent.  The only uncommitted
  files at reconciliation were the explicitly scoped L303 proof,
  checker, dataset, and its backward-compatible L302 helper refactor.
  No mathematical result or next-action instruction from the stopped
  agent remains unaccounted for.
- The L312/A260 checkpoint repeated the audit from a clean worktree.
  All four overlap commits are ancestors of `master`; `master` is the
  only local branch and worktree and agrees with `origin/master`.
  L248 and L249 were regenerated again in `/tmp`, byte-identically,
  with the two hashes above; Ruff, `py_compile`, and the interrupted-
  artifact check pass.  The stopped-agent recommendation to continue
  L248 at the physical multiplier is historical: L277--L279 have
  since closed that volume coefficient.  the historical frontier at that checkpoint is
  L312's still-missing **grouped weighted** L306/L307 recurrence,
  followed by the fixed-`theta=1/2` finite-flag induction.  A259
  remains conjectural and is not a competing task.
- L313--L316 were then derived on the sole-agent line.  L313 closes
  weighted finite physical words, L314 closes the compulsory frame
  cross, and L315 preserves displayed bridge-root valuation through
  the complete affine recurrence.  L316 now proves that a grouped
  closed endpoint is even under simultaneous sign reversal of every
  transfer channel.  Thus a single exposed L315 root cannot survive
  alone outside the response; the historical debt at that checkpoint is the sharp
  **two-channel** weight/energy estimate and its margin assembly.
- **Historical pre-L279 route record (superseded):** A194 at the complete
  physical terminal balance, now in A196/L250's balanced output
  coordinates, A197/L251's paired analytic ports, A198/L252's cyclic
  radial quotient, A199/L253's Hardy-window energy target, and
  A200/L254's shifted-left-Hardy complement, now represented
  physically by A201/L255's doubled-Hardy pencil and anchored at
  relative grade one by A202/L256.  A203/L257 removes the remaining
  realization choice and states the historical step at that checkpoint as one canonical
  model-space coefficient recursion.  A204/L258 now eliminates the
  final output row exactly and rewrites its whitening as a closed-return
  renewal, so the historical calculation at that checkpoint contains no one-sided port term.
  A205/L259 identifies every apparently premature future transfer
  coefficient as an off-diagonal block of the first active causal
  Toeplitz-leakage row; the remaining question is whether the physical
  return closes that complete row in the ordinary Hardy metric at the
  active face.  A206/L260 blocks a false orientation-by-orientation
  proof, while A207/L261 proves that the full row norm then collapses
  to its diagonal energy by matrix-inner Parseval.  A208/L262 proves
  that every genuine normalized half-line Wold sandwich uses exactly
  this ordinary Hardy metric.  The remaining historical issue at that checkpoint is physical
  channel placement: prove that L258's first nonconstant closed return
  is that sandwich after the two ellipse orientations have been summed.
  A209/L263 additionally proves that the fixed physical Schur row and
  the moving Wold/defect graph cannot alter a first new face; only the
  still-unidentified physical channel numerator remains.  A210/L264
  supplies the final legal gauge choice: restore the universal
  right-half-line metric while continuing to delete the active left
  orbit, final corner, and cross row.  A211/L265 rewrites the target
  itself as four times one finite right-half-line Stein divergence.
  A212/L266 proves that any copy-scalar Laurent return—and any
  unilateral shift polynomial whose boundary depth is at most the
  delay—contributes only its bilateral-symbol constant between the
  delayed leakage projections.  A213 records the required synthesis.
  L272/A219 now closes its first stop condition exactly: shifted
  Hankel intertwining and the model-space complement normal-order
  every closed copy-dependent cycle into the two-sided leakage ideal,
  while L243's product valuation leaves only one leakage selection at
  the active delayed coefficient.  This is not a proof of A194:
  copy-scalar boundary depth/support and the lower-face structure
  remain open.  L277 later closes the odd part exactly.  A214/L267
  supplies the classical exact
  defect-transport identity for an abstract unitary Redheffer
  feedback.  L268/A215 has now rejected its hoped-for one-shot
  physical use: on the exact grade-two monomial, the full active
  closed-return quotient is `diag(2,2)` while the deflated tail face
  is `diag(0,4)`.  Their traces agree, but ranks two and one cannot be
  related by the first face of a single analytic entrance carrying
  only the tail defect.  A larger network would require an additional
  defect channel and therefore would not make A194 automatic.  The
  historical calculation at that checkpoint returned to A213's
  trace-locality conditions; no
  operator congruence can replace them.  L269/A216 then corrects the
  second condition: the full face does contain a nonzero
  boundary-depth-`k+1` term, already at grade two, but it appears as
  `G_k−SG_kS*` with zero delayed trace.  The all-grade target is
  therefore a shallow L266 return **plus** explicitly trace-null deep
  Stein divergences, not a wholly shallow operator face.
  L270/A217 removes the separate radial-coefficient calculation:
  conditional on L269's structural normal form, arbitrarily long
  exact monomial axes force `p_k(1)=0` and `p_k'(1)=−1`; L252 then
  proves complete one-delay trace covariance for every multiplicity.
  The historical A194 gate at that checkpoint was to derive that
  normal form, including its lower vanishings, in arbitrary grade from
  L243/L251/L258.  The finite grade-`2..5` pattern is still not an
  induction.  L271/A218 reduces that two-system theorem to one
  simpler full-face support formula:
  `D_k=U_k(Q)−R_1+R_k−4R_(k+1)+H_k`, where `u_k(1)=0` and `H_k`
  is one explicit four-word remote packet.  Exact deflation algebra
  then produces L269 automatically.  Proving this single-face
  support formula and its lower vanishings is the unique historical step at that checkpoint.
  L272 removes the separate unpaired-`B#`/ideal-placement obligation:
  the remaining calculation is the unilateral copy-scalar return
  between the two leakage projections.  L273/A220 now proves the
  retained-metric half of that return in arbitrary grade: it cancels
  every interior word of the universal fan with coefficient `−1`.
  L274/A221 now also fixes all four remote coefficients of the
  unweighted full-minus-tail recursion in arbitrary grade.  Each
  crossed coefficient is the physical cancellation
  `1+4(−1)^k−4(−1)^k=1`; `G_k` comes only from the full direct Gram
  and `−SG_kS*` only from the subtracted tail direct Gram.  The unique
  nonradial gate is therefore no longer endpoint normalization.
  L275/A222 now completes the coefficient problem for the entire
  `2k+2`-word fan.  After the final-row loops have been resummed, the
  inherited fan is the coefficient-shifted tail fan; L274 supplies
  the two new crossed endpoints and `G_k−SG_kS*`, and exact tail
  embedding propagates coefficient one to every remaining word.  The
  raw direct/immediate/later pieces do not obey this recursion
  separately.  L276/A223 now removes the unnecessarily strong exact-
  support debt.  It is enough that lower even faces be two-sided
  radial (with odd faces zero) and that the active
  associated difference be cyclically radial through index `k+2`.
  Long exact axes kill the lower radial operator and the first two
  active trace moments; the active monomial/full-tail comparison
  kills the top energy coefficient.  Thus possible words outside the
  fan need only be controlled modulo trace, not literally excluded.
  L277/A224 now closes every lower odd face: before the next delay is
  imposed, the first odd frontier has a universal twelve-word form and
  factors through the two-sided ideal generated by that next delay
  defect and its adjoint.  The even retained metric therefore cannot
  reintroduce an odd weighted face.  L278/A225 now closes every lower
  even face as well: with the physical metric intact, the first even
  frontier is exactly
  `(ES^hF)(S*)^h+S^hF(S*)^hE`, so the next delay annihilates it.
  All lower coefficients are zero operator-wise.  L279/A226 then
  linearizes L258's exact closed renewal under deletion of the active
  metric coefficient.  The response and L278's intact face are each
  cyclically the same two-unit radial flux, so the edge-deleted
  log-volume coefficient is exactly `4||B_k||²`, the delayed
  dual-Schur trace is `2||B_k||²`, and the effective separator trace
  is `−16||B_k||²` in every grade.  This closes A194 and L222's
  pointwise range obstruction.  The live repeated-block gate is now
  bounded selection through commutant-rank changes.  L280/A227 now
  rewrites the whole response range as the coboundaries of the
  bistochastic transfer-channel Markov operator
  `I−Phi Phi*`, with an explicit analytic state preimage and exact
  Dirichlet form.  L281/A228 then replaces that long column by the
  exact observability defect column `2R_H`.  Its squared norm is twice
  the Markov Dirichlet energy, so a closing spectral gap causes no
  separate state-synthesis blow-up.  L282/A229 dualizes its least
  Dirichlet cost exactly.  L283/A230 then closes the full matrix face
  on every completely delayed stratum: L258's first-face similarity
  and the opposite L247/L279 metric responses cancel operator-wise,
  leaving L278's intact frontier.  L246 conjugates it to L228's former
  all-grade anticommutator, and L227's canonical analytic repair has
  exact complete-delay faces `+B_k*B_k` and `−12B_kB_k*`.
  This does **not** close bounded repeated-elliptic selection.  A172
  proves that mixed odd coefficients of the unmodified repair can be
  indefinite when earlier transfer rows are nonzero but singular.
  The live elliptic frontier remains A178's arbitrary-grade
  right-ideal preparation, with L282's uniform flux inequality as an
  alternative route.  L284/A231 now removes one finite-gauge
  obstruction: each stored order-four--six column has a globally
  endpoint-null normalization that vanishes on its matching complete
  delay.  L285/A232 derives the exact triangular transport of each
  normalization and disproves automatic delay-ideal preservation at
  its first successor.  After the compulsory lower correction, every
  successor is exactly one L204/L280 endpoint homology equation.
  L286/A233 solves the first such equation exactly: a bounded
  double-delay-divisible correction transports the normalized
  quartic gauge through quintic order without changing L232--L233's
  endpoints.  L287/A234 then proves by an exact scalar trace
  obstruction that the next even metric coefficient cannot generally
  be made endpoint-null.  The live rule must alternate odd
  homological cancellation with control of a nonzero even physical
  face.  L288/A235 now solves that finite even successor under the
  correct weaker condition: the recomputed sextic column is
  triple-delay divisible, and its complete upper-endpoint motion,
  including the nonlinear quartic Schur cross, factors through
  `ESF`.  It therefore vanishes on the surviving partial flag and
  leaves L234's effective sextic budget unchanged.  The nonlinear
  lower Schur coefficient is not globally invariant.  L289/A236 now
  packages every nonlinear endpoint cross in the exact mixed-graph
  identity
  `U(Htilde)−U(H)=J_Htilde*(Htilde−H)J_H`.
  It also corrects the convergence target: at fixed repeated length
  `L`, L201's invertible terminal `B_L` and L194's analytic tail mean
  that a positive prepared jet through order `2L` is sufficient.
  No infinite summable selected-column series is required.
  L290/A237 then identifies the exact graph-compatible sufficient
  invariant: a Hermitian state motion ending in exposed channels
  `ES^jF` is carried by the full nonlinear Schur shorting into
  endpoint factors `X_jB_j*+B_jY_j*`, so it vanishes on the surviving
  transfer flag despite arbitrary graph motion.  This is a
  one-sided/hereditary condition, not L285's false raw two-sided
  ideal claim.  L291/A238 then proves that state heredity is stronger
  than necessary: every positive-interior state transport is a
  locally bounded endpoint-quotient-fixing graph congruence plus its
  literal endpoint residual, and flag-zero endpoint residuals are
  exactly `X[B_1 ... B_r]*+[B_1 ... B_r]X*`.  Extracting bounded
  physical endpoint factors directly—without the pointwise
  pseudoinverse used only to prove equivalence—is reduced by
  L292/A239 to finite Smith-coordinate valuation inequalities along
  L197's analytic failure arcs.  Deriving those physical
  arbitrary-grade divisibilities and proving uniform margins through
  grade `L` is one live certificate.  L293/A240 now sharpens the
  equivalent L282 route: in L218's independent matrix-Schur
  coordinates, the first Markov gap is exactly the simultaneous
  commutator Laplacian of the first Schur jets.  Its exponent and
  square root are therefore controlled.  L294/A241 now supplies the
  exact nonlinear numerator pairing for L212's canonical
  channel-minus-Gram target and proves L282 with constant `28`,
  including flags, and more generally controls every
  `−aB_kB_k*−bG_k`, `a>=0`, with constant `|b|`.  L295/A242 catches
  a cross-normalization error in the proposed continuation:
  L227/L285's **raw upper** repair face is `−12I` at the monomial
  apex, while L212's **lower-tight** target is `−16I`, so that
  mismatched comparison is false.  L297/A244 now identifies the
  entire discrepancy as lower retightening: the endpoint ledger is
  `(lower,upper)=(+I,−4I)` for the boundary metric,
  `(+I,−12I)` for the raw repair, and `(0,−16I)` for the
  lower-tight/effective face.  Thus the raw-to-effective increment is
  `(-I,-4I)`, so the correctly normalized effective comparison with
  L212 is not obstructed.  L298/A245 then supplies the
  noncommutative transport directly: the negative right-orbit Gram
  contributes `(-B_k*B_k,-4C(B_k*B_k))`, and one seventh of L212
  returns `+4(C(B_k*B_k)-B_kB_k*)` upstairs.  The channel terms
  cancel, giving the exact Stein-compatible pair
  `(-B_k*B_k,-4B_kB_k*)` in every grade.  Added to L283, this closes
  the completely delayed lower-tight face as `(0,-16B_kB_k*)` in
  arbitrary grade; weighted sums close the fixed-`S`
  base-colligation transport.  The live effective-route debt is now
  the nonlinear defect produced by the moving elliptic operator and
  prior prepared grades.  L299/A246 reduces that defect exactly to
  six quadratic terms and proves it begins cubically.  A rational
  rank-one flag example has nonzero cubic compression on
  `ker B_1*`, so the no-extra-correction insertion is false: the
  defect is not automatically a prior-flag factor.  L300/A247 now
  closes the resulting full cubic range and bounded-energy problem:
  after eliminating the moving frame, every dual pairing is
  controlled by the exact L281 observability defect, and a response
  column exists with
  `||C_3||<=10||R_1||<=20 Gamma_S||B_1||`.  The frame elimination
  uses L298's special `(S*)²F_1=0`, so it is not automatically
  recursive.  L301/A248 now closes the generic first successor:
  every direction `X−S*XS=VC*+CV*` produces a moving successor in
  the Markov range with preimage norm at most `10||X||+2||C||`.
  Thus the retained quartic commutant pairing is independent of the
  chosen L300 preimage.  L302/A249 computes both ends of that
  invariant.  Its upper class is
  `P_4=(15/4)B_1B_1*+(1/2)(B_1B_1*)²`, but its exact lower corner is
  `K_4=6B_1*B_1+(9/4)(B_1*B_1)²`.  Thus the tempting metric-only
  cancellation has harmful lower motion `−K_4` and is not a
  completed lower-tight correction.  After the compulsory parallel
  lower neutralization, the upper class is the prior-flag cost
  `(9/4)B_1B_1*+(7/4)(B_1B_1*)²`, which vanishes on `ker B_1*`.
  L303/A250 subsequently incorporates every moving Schur-graph cross
  at fixed partial scale.  It restores the quartic lower face and
  proves compatible lower/upper margins for every `0<theta<1`.
  L304/A251 closes the remaining gap-free analytic
  ordered-flag/Smith realization at `theta=1/2` by factoring the
  complete residual into relative commutators containing `ESF` and
  bounding the off-cokernel dual-telescope error.  L305/A252 then
  supplies the arbitrary-word mechanism without assuming abstract
  cyclic excision: every rooted bridge word is one explicit retained
  channel plus an ideal-valued bounded response.  The response,
  lower elimination, and subsequent transport preserve the
  cumulative bridge module.  L306/A253 now evaluates the complete
  rooted class without a growing word expansion.  Wold resolution and
  response-cokernel covariance collapse all six L299 terms to the
  single initial-copy energy matrix
  `q=−Psi(Delta*XS+S*XDelta+Delta*XDelta)
  −(C*H+H*C+C*C)`, modulo another gap-free module-preserving
  response.  L307/A254 proves this is the actual affine recurrence,
  not only its first homogeneous step: a correction canceling a
  preceding forcing leaves the same six-term successor.  Scaling by
  `theta` adds the exact positive state/copy terms
  `theta(1−theta)CC*` and `theta(1−theta)C*C` relative to the scaled
  successor.  For an affine forcing `R!=0`, however, the predecessor
  also leaves `−(1−theta)R`; the Gram term is not alone a recursive
  margin.  L308/A255 then falsifies the one-shot raw-parity shortcut:
  the unprepared weighted L298 copy energy has exact nonzero quintic
  scalar trace `33264/15625`, even though its cubic is a response.
  L309/A256 then performs that triangular preparation exactly on the
  same scalar-copy colligation and disproves the stronger target too:
  after cubic cancellation, quartic lower neutralization, and
  quintic lower neutralization, the physical upper quintic is
  `66384/15625`, while the scalar response space is zero.  Since
  `B_1!=0`, this is supported on an already active transfer range and
  does not obstruct the finite flag.  the next historical step at that checkpoint is therefore
  to formulate the prepared affine recurrence modulo bounded
  hereditary cumulative-transfer factors, not global odd-response
  parity.  L310/A257 now closes that endpoint factor termwise:
  polarized L280 converts every L305 rooted channel exactly to
  `KB_j*+B_jK*` plus a bounded response, with automatic L292
  valuations.  Its converting column can remain nonzero and
  endpoint-null when `B_j=0`, but L311/A258 now closes that
  normalization debt: for the physical root, the relative
  commutator with `X=FqpE` has exactly the channel-minus-factor
  endpoint, while its four-word defect expansion has pairwise
  identical L305 roots.  The resulting response column remains in
  the bridge ideal and is `O(||B_j||)`, so it vanishes on every
  complete delay.  the historical step at that checkpoint is now only the prepared finite-jet
  sign/margin assembly: absorb the bounded hereditary odd factors
  and prior-transfer even costs using the margins retained at
  `theta=1/2`.  Do not factor the raw `q`, return to a growing graph
  expansion, or start an unrelated isolated-grade march.
  The
  alternate raw route retains its lower
  budget, separates L283's favorable `−12` Gram, and expresses only
  L285/L289's mixed transported remainder as a controlled flux.
  L296/A243
  removes the unnecessary restriction to the special quadratic
  `G_k`: every polarized L280 response has an exact observability-
  defect pairing, and a family of polarization columns has L282
  energy controlled by the sum of their squared Frobenius norms.
  The live repair-branch target is now to discard L290's bounded
  prior-flag factors, derive polarization columns only for the
  surviving compressed remainder, and prove their total energy is
  `O(||U*B_k||²)`.  Do not mix the raw and lower-tight endpoint
  ledgers or treat L298's fixed-`S` weighted sum as an already
  completed moving metric series; L299 explicitly disproves that
  shortcut.
  The circular-normal merger is downstream.
  A195/L249 is only a guardrail against detaching the doubled terminal
  edge, not a competing approach.

## NEWEST (2026-07-26): L315 proves the full affine recurrence preserves elliptic bridge valuation
- Give a state coefficient at elliptic degree `d` the filtered
  decomposition `sum_(2k<=d) Z_(d,k)`, with
  `Z_(d,k)` in the two-sided bridge ideal `I_k` generated by
  `G_k=ES^kF` and its adjoint.  Each component must carry a displayed
  bounded `pG_kq/pG_k*q` presentation; abstract ideal membership is
  not accepted because the generated algebra can be full.  Apply the
  same definition to `C_dV*` for columns.
- L314 sharpens the L298 seed to
  `X_k,C_kV*,VC_k* in I_k` exactly at degree `2k`.
  Two-sided multiplication, adjoint, products, and stable Stein
  inversion preserve each root.  Since `Delta=A−S` and `H=D−V`
  start at degree one, every linear moving successor gains a degree.
- L285's lower parallel column has lift
  `C_parallel V*=−(1/2)E R E`, so lower elimination also preserves
  the filtered ideal.  L305/L311's physical normalized response has
  its lift in the same `I_k`.  These facts close the induction under
  every finite prepared affine step.
- Consequently transfer grade `k` can never be exposed before degree
  `2k`.  L311 rewrites a rooted term `c^d pG_kq` as a literal factor
  through `c^kB_k` plus an ideal-valued response of norm
  `O(c^(d−k)||c^kB_k||)`.  This proves the rank-stable L292
  valuations without separately telescoping the dual Stein closure.
- Root weighting is now closed.  The remaining elliptic gate is
  narrower but still substantive: factor the exposed multipliers
  through already active weighted transfer rows (using L313), then
  show their square-completion cost is paid by retained earlier
  margins at fixed `theta=1/2`.  This is sign/energy bookkeeping, not
  another root-causality problem.
- Twelve mixed, rank-collapse, and complete-delay records pass.
  Dataset SHA-256:
  `e98a59d70874201afb32272a3e2408e38bda2bc00019f8fdf7cef5ca1ec4b3fb`.

## NEWEST (2026-07-26): L314 proves the complete raw-frame cross is weighted bridge-causal
- The canonical repaired frame has the exact graph factor
  `D=Hcal V(V*Hcal V)^(-1/2)`, with
  `Hcal=P_bl−A*P_bl A`.  Thus every finite frame coefficient is a
  finite physical word polynomial; the stable Stein repair never has
  to be expanded inside the frame.
- L298's grade-`k` column has the exact lift
  `C_kV*=−(1/2)G_kG_k*
    −(1/2)Q{S^kG_k*+sum_(j<k)(S*)^(k−j)G_kG_j*}`,
  where `G_k=ES^kF`.  It is therefore rooted in the physical
  grade-`k` bridge ideal.
- Every coefficient of `(D−V)C*+C(D−V)*` is consequently a finite
  rooted physical word carrying degree at least `2k`.  L311 converts
  it to a bounded literal factor through `c^kB_k` plus an
  ideal-valued response of size
  `O(c^k||c^kB_k||)`, which vanishes whenever `B_k=0`.
- This corrects L313's deliberately conservative scope guard.
  L217 still forbids factoring the raw frame by itself, but the frame
  never appears by itself in L306.  Its compulsory pairing with
  `C_k` already restores the weight, so no cancellation against
  `Psi_S(Delta*XS+S*XDelta)` is required.
- The remaining weighted debt is now only the moving-operator copy
  closure
  `Psi_S(Delta*XS+S*XDelta+Delta*XDelta)`, followed by the
  fixed-`theta=1/2` finite-jet margin assembly.  Keep `−C*C` intact.
- Fifteen unstructured/rank-chain/complete-delay records pass.
  Dataset SHA-256:
  `c0883301a48aef9211fd62e1c7ee09d418c4a05bdac09eaf48c2232fe476c514`.

## NEWEST (2026-07-26): L313 proves weighted causality for every finite physical reverse-edge polynomial
- Keep the enhanced physical reverse edge
  `J=(I+F)S*(I+E)` intact.  On the forward defect chain,
  `S*S^hW=S^(h−1)W−VB_(h−1)*`,
  `J*W=2SW+2VB_1*`, and
  `J*S^hW=S^(h+1)W+VB_(h+1)*`.
- Following the unique chain-valued path through a word with `r`
  copies of `J` emits only `VB_j*` with `j<=r`.  Passing an emission
  through the remaining suffix gives the explicit polynomial
  factor
  `W*w(S,J)V=sum_(j<=r)B_jC_j(w)`, with stacked norm at most
  `(|w|+1)4^r`.
- Therefore every causal physical polynomial
  `P=sum alpha c^d w`, `d>=#J(w)`, factors through the weighted row
  `[cB_1,c²B_2,...]`, with a displayed rank-stable bound.  L125's
  scalar filtration shows that every finite jet of the actual direct
  map `phi_c(S+cJ)−S` is causal in this sense.
- This is the first genuinely weighted version of L312 and closes the
  finite conformal-map part of L306's multiplier problem.  It also
  explains why expanding `E,F` inside `J` and counting raw `S*`
  letters gives the wrong filtration.
- Scope guard: L217 proves that the raw Hermitian-gauge defect frame
  contains premature future-transfer images.  L313 must not be
  applied to that frame separately.  The remaining grouped weight
  debt is now the specific cancellation of the noncausal part of
  `C*H+H*C` against
  `Psi_S(Delta*XS+S*XDelta)`.  After that pairing, L313 and the
  favorable quadratic squares can be used.
- Forty-two general/rank-chain/complete-delay records pass; each also
  checks all 510 physical words through length eight.  Dataset
  SHA-256:
  `5e20fe4db61c1b80582436faa3c2c98b2d6699fd329d6b006f895196523e0252`.

## NEWEST (2026-07-26): L312 gives every finite endpoint word a constructive transfer-prefix normal form
- For a length-`d` word `w`, repeatedly reduce its first reversal
  `a^p s v` by `a^p s=a^(p−1)(I−VV*)`.  Each step emits the exact
  term `−B_(p−1)V*vV` and deletes two letters; a terminal `a^r`
  emits `B_r`.
- Hence `K_w=W*wV=sum_(j<=d)B_jC_j(w)`.  There are at most
  `ceil(d/2)` contraction emissions, so the stacked factor has norm
  at most `ceil(d/2)` and
  `K_wK_w*<=ceil(d/2)² sum_(j<=d)B_jB_j*`.
- This is exact, polynomial, rank-stable, and dimension-free.  It
  proves a weaker fallback to A259 without a pseudoinverse.  A259's
  constant-one bound remains conjectural.
- Each emitted transfer index is no larger than the number of
  `S*` letters.  This is useful for finite elliptic words, but a
  termwise application after an infinite Stein/Hardy expansion can
  lose the coisometric cancellations that enforce physical
  valuation.  the historical gate at that checkpoint is therefore a **grouped weighted**
  normal form for L306/L307 in L220's Schur-orthogonal coordinates.
- The enlarged tracked dataset regenerates with SHA-256
  `5b5dac6046295c1176d2878365e1627961603282b368bb5e8c14fcd20a56becf`.

## NEWEST (2026-07-26): A259 isolates a strong unweighted endpoint-word conjecture, but its proof is open
- For every length-`d` word `w` in `S,S*`, exhaustive tests support
  `K_wK_w* <= sum_(j=1)^d B_jB_j*`, where
  `K_w=W*w(S,S*)V`.
- Defect telescoping proves that this is exactly
  `||V*w*Wy||²+||S^(d+1)Wy||²<=||y||²`.  Hence a proof would give a
  constant-one, dimension-free Douglas factor for every L311 word
  multiplier through the first `d` transfer blocks.
- Ten deterministic general/rank-chain/delayed colligations pass all
  2,046 words through length ten; separate exhaustive tests through
  length sixteen and random tests through length twenty also pass.
  The enlarged L312 tracked dataset has SHA-256
  `5b5dac6046295c1176d2878365e1627961603282b368bb5e8c14fcd20a56becf`.
- This is **not L312**.  The tempting Julia-colligation circuit proof
  has an unproved residual-block step after mixed gate orientations.
  The nonlinear word `(S*)²S³(S*)²` already contains
  `−(I−B_1B_1*)B_1`, so naive Dyck-path straightening is inadequate.
- Do not stack the conjecture as a lemma.  Either prove its explicit
  finite-horizon residual cancellation, or use it only as a guide.
  More importantly, even a proof is unweighted: the live physical
  gate remains causal elliptic weight matching in the complete
  L306/L307 series.

## NEWEST (2026-07-26): L311 makes every physical endpoint split ideal-valued and delay-normalized
- For a physical root `T=pG_jq`, with
  `K=W*qpV`, its copy multiplier has the exact state realization
  `X=WKV*=FqpE`.
- The Hermitian relative commutator
  `[G_j,X]+[G_j,X]*` has upper endpoint
  `Phi(B_j*K+K*B_j)−(KB_j*+B_jK*)`.  It therefore replaces L310's
  standalone polarized converting column without changing the
  literal hereditary factor.
- Expanding `FqpE` as
  `r−rS*S−SS*r+SS*rS*S`, `r=qp`, gives four differences
  `G_jx−xG_j`.  Each pair has exactly the same L305 root, so the
  rooted quotients cancel and L305 supplies an explicit response
  column whose lift remains in the bridge ideal.
- Adding that column to L305's original rooted column gives
  `W*G_S(T+T*)W=KB_j*+B_jK*+M_S(Chat)`, with
  `Chat V* in I_j` and
  `||Chat||_F <= (5|p|+4|q|+17)||B_j||_F`.
  It is polynomial and vanishes identically when `B_j=0`.
- Complete-delay audits expose the normalization gain sharply: the
  raw L310 columns have norms approximately `1.02`, `1.21`, and
  `1.44`, while the physical ideal-valued columns vanish exactly.
  All 24 general/rank-chain/delayed records pass.  Dataset SHA-256:
  `785aaab8462afccce9c5693a31677c0215c7fd714bcc37ec28fa512bc37ec8c6`.
- This closes response selection, endpoint heredity, Smith
  divisibility, and complete-delay normalization termwise.  It does
  not prove the sign of the complete prepared coefficient.  The
  remaining elliptic gate is the all-grade finite-jet
  square-completion/margin rule for L306/L307.

## NEWEST (2026-07-26): L310 converts every root to a literal hereditary endpoint factor
- Polarized L280 gives the exact balanced identity
  `Phi(B_j*K+K*B_j)=KB_j*+B_jK*−M_S(C_j(K))`, with
  `C_j(K)=Q{S^jWK+sum_(ell<j)(S*)^(j−ell)VK*B_ell}`.
- The column is perpendicular, polynomial/analytic, and satisfies
  `||C_j(K)||<=j||K||`.  Combining it with L305 turns every rooted
  word `p(ES^jF)q` into the literal factor
  `K_pq B_j*+B_jK_pq*` plus a bounded response.
- The literal factor vanishes on `ker B_j*` and has an explicit
  bounded multiplier, so its L292 Smith valuations are automatic.
  This closes the hereditary endpoint quotient required by L309,
  term by term and without a pseudoinverse.
- Scope guard: the converting column need not lie in L305's bridge
  ideal or vanish on a complete delay.  Three exact-structure
  numerical cases have `B_1=0`, factor/response zero, but
  `||C_1(K)||=||K||=1`.  It is a nonzero endpoint-null direction and
  cannot be inserted blindly without reviving L284--L288's
  later-successor normalization debt.
- the new historical gate at that checkpoint is therefore narrower: normalize the complete sum
  of polarized columns compatibly with every delay flag, or prove
  their L301/L307 successors remain hereditary with sufficient even
  margins.  L311 subsequently closes this debt for every physical
  root by an ideal-valued relative-commutator normalization.  Dataset
  SHA-256:
  `376523b537a27d4f30f65705e8e6ab536220bf65256901294d931777a7064ae1`.

## NEWEST (2026-07-26): L309 disproves global prepared odd-response parity
- L309 repeats L308's rational scalar-copy colligation, but now
  carries out the valid triangular preparation: the half-scale
  grade-one direction, cubic cancellation, the direct grade-two
  direction, L303's quartic lower neutralizer, and then the complete
  quintic lower neutralizer.
- The first two transfers are `B_1=−4/5`, `B_2=9/25`.  The exact
  cubic upper coefficient is zero and the complete quartic pair is
  `(-81/1250,-2606/625)`, exactly matching L303 with the direct
  grade-two face included.
- Before fifth neutralization the lower quintic is `3384/3125`.
  Adding `C_5=−(1/2)V(3384/3125)` makes the lower coefficient
  exactly zero, but leaves physical upper
  `66384/15625=4.248576>0`.
- The scalar response space is zero, so this is an exact obstruction
  to the proposed statement that every prepared odd coefficient is a
  global response.  Exact rational formal Schur calculation and an
  independent floating reconstruction agree below `2.1e-14`.
- This is not a finite-flag obstruction: `B_1` is already nonzero, so
  no first surviving flag remains in the scalar example.  The
  corrected all-grade target is bounded hereditary factorization of
  every lower-neutral odd retained class through the cumulative
  active transfer row; only its surviving-flag compression must be a
  response.  L292's Smith valuations remain the rank-changing test.
  Dataset SHA-256:
  `e7aa72f5f20dd92ac9b79af64b15106bb792abf0925ef6c3ca7845619cf4761a`.

## NEWEST (2026-07-26): L308 disproves one-shot odd parity for the raw copy energy
- For the unprepared weighted L298 correction,
  `X=sum c^(2k)R_k`, `C=sum c^(2k)F_k`, L300's cubic copy trace is
  exactly zero, but the raw quintic trace is not identically zero.
- An exact scalar-copy certificate uses
  `S=[[0,3/5,−4/5],[0,4/5,3/5],[0,0,0]]`,
  `V=e_1`, `W=e_3`.  It is a stable partial isometry with
  `B_1=−4/5`, and reduced word algebra gives
  `[c^5]tr q_raw=33264/15625>0`.
- Every homogeneous endpoint response has trace zero.  In copy
  dimension one that makes the response map identically zero, so the
  quintic value cannot be hidden in L305's response part.
- The preparation is essential: cancel the cubic through L300/L307,
  carry L301's quartic successor and L303's compulsory lower
  neutralization, and only then form the next odd successor.  L309
  subsequently proves that even this prepared odd coefficient need
  not be globally response-null; its viable replacement is
  hereditary cumulative-transfer support.  L306's raw `q` cannot be
  factored directly as an even Gram series.
- The exact audit has zero cubic cyclic classes, eight nonzero
  quintic cyclic classes, and an independent floating matrix assembly
  on the rational input gives `2.128896`.  Dataset SHA-256:
  `0b755aec840af01085b7edf84078a3299174e025a2ef70d26360829ab8f603ec`.

## NEWEST (2026-07-26): L307 makes the copy reduction recursive and exposes its Gram reserve
- For an affine correction
  `X−S*XS=R+VC*+CV*`, exact expansion gives the moving residual
  `R+N(X,C)`, where `N` is the same six-term L299 successor.  If the
  preceding residual is `−R`, it cancels before the successor is
  evaluated.  Thus L306 applies at every recursive correction, not
  only to the first homogeneous L298 direction.
- The initial-copy boundary closes exactly:
  `Psi_S(VC*+CV*)=C*V+V*C`.  Hence the complete copy value is
  `Psi_S(R)+q(X,C)`, and L306's homogeneous fixed-base copy term is
  identically zero.
- At a real partial scale `0<=theta<=1`,
  `N(theta X,theta C)=theta N(X,C)+theta(1−theta)CC*` and
  `q(theta X,theta C)=theta q(X,C)+theta(1−theta)C*C`.
  Both extra terms are positive semidefinite.  The campaign's fixed
  `theta=1/2` therefore has a one-quarter correction-frame Gram term
  in both state and copy space relative to half the full-scale
  successor.
- Scope guard: if `R!=0`, the scaled affine correction cancels only
  `theta R` of the preceding `−R`, leaving `−(1−theta)R`.  L307 is
  not by itself an arbitrary-step positivity induction.  The live
  gate remains: prove
  the full-scale odd part is a bounded response and its even retained
  part is quadratically divisible by prior transfer rows, then combine
  the exact scaling term with the controlled affine forcing and earlier
  strict face margins.
- All 16 tracked affine-recursion cases pass at scales
  `0.2,0.5,0.8`, including nonvacuous grades `2,3` after a complete
  first delay and three rank-collapse scales.  Dataset SHA-256:
  `0b5f54dfb20b9fba53971ab8d03f7daad7833f4b020b7ee0c8fcf80b9dc4e385`.

## NEWEST (2026-07-26): L306 collapses the complete rooted class to one copy energy
- Put `Psi_S(Z)=V*sum_(n>=0)S^nZ(S*)^nV`.  Wold resolution and
  L305 show that the upper Stein endpoint of a bridge-ideal state
  `Z` is `Phi(Psi_S(Z))` modulo a perpendicular response whose norm
  is linear in the active bridge amplitudes and whose lift remains in
  their cumulative module.
- If `Delta=A−S`, `H=D−V`, and the prepared correction satisfies
  `X−S*XS=VC*+CV*`, the complete L299 endpoint therefore has retained
  representative
  `Phi(q)`, where
  `q=−Psi_S(Delta*XS+S*XDelta+Delta*XDelta)
     −(C*H+H*C+C*C)`.
  This is also exactly the moving full copy energy minus its fixed-base
  copy energy.
- The response-cokernel identity `Z_Y=0` makes its observability
  Gramian commute with `S,S*` and intertwine every `B_j`.  The physical
  Hermitian-gauge frame and L298/L305 prepared columns are covariant
  for those intertwinings, so all three frame products close in the
  same finite copy matrix.  Off the cokernel, each failure contains
  both an observability defect and an active bridge, giving the
  gap-free response bound.
- This is an all-series **structural reduction**, not a positivity
  theorem.  The exact remaining historical gate at that checkpoint is to prove the odd coefficients
  of `q` are response-null and its even coefficients, after compulsory
  lower neutralization, are two-ended transfer Grams dominated by
  earlier retained margins.  Do not drop the negative `C*C` cross or
  split its potentially large pieces.
- All 16 tracked mixed-grade cases pass, including nonvacuous
  grades `2,3` on the complete-first-delay stratum and three
  rank-collapse scales; 100 additional mixed-grade stress cases also
  pass.  Dataset SHA-256:
  `cc538f27b3c5395a90def31a062bad12bf2449589cecf0823b139764f1275d65`.

## NEWEST (2026-07-26): L305 makes every rooted bridge response constructive
- For every word term `T=p(ES^jF)q`, the complete upper Stein endpoint
  splits exactly as
  `Phi(B_j*W*qpV)+h.c.+M_S(C_(p,q,j))`.
- The perpendicular column is explicit, obeys
  `||C_(p,q,j)||<=(|p|+1)||B_j||`, and has
  `C_(p,q,j)V*` in the same two-sided physical bridge ideal.  This is
  rank-stable and uses no pseudoinverse or Markov-gap inverse.
- The two terms of a relative commutator
  `[ell,p(ES^jF)q]` have the same root, so the retained channels cancel
  and their difference is an explicit ideal-preserving response.
  L304's 28 terms are a finite instance.  No general assertion
  `I intersect [A,A]=[A,I]` is needed or made.
- L298's grade directions have both metric and frame lifts in the
  cumulative bridge ideal.  Consequently L299's complete six-term
  moving defect lies there.  L285 lower elimination, Stein inversion,
  and the L305 response column preserve it, closing the response/module
  part of the finite recurrence.
- The remaining algebraic gate is now only the **rooted retained
  quotient**: prove that its complete physical sum cancels at odd
  grades and is an earlier-margin-dominated two-ended transfer Gram at
  even grades.  L300 and L303--L304 establish those statements only at
  cubic and quartic order.
- All 38 unstructured, rank-changing, and completely delayed audits
  pass, including relative-root cancellation and independent
  ideal-lift reconstruction.  Dataset SHA-256:
  `b7602e5634dc372b2519f1673722ca3dce9af17b73b8c87709338dec37de920b`.

## NEWEST (2026-07-26): L304 closes rank-stable selection through the quartic
- Fix `theta=1/2` and subtract L303's retained balanced upper class
  `−(5/4)B_1B_1*−(3/4)(B_1B_1*)²` from the complete lower-neutral
  quartic coefficient.
- Put `G=ESF=VB_1*W*`.  Exact word algebra factors the remaining
  effective state trace lift into 28 relative commutators
  `alpha[ell,p(G or G*)q]`, with
  `sum|alpha|=177/16`.  Pairing against an observability Gramian
  supplies `||Z_Y||` from `[H_Y,S]` or `[H_Y,S*]`, while the bridge
  supplies `||B_1||`.  Thus the finite part is bounded by
  `(177/16)||B_1||||Z_Y||`.
- L302's one nonlocal dual-telescope substitution remains controlled
  off the cokernel:
  `|epsilon_Y|<=33 Gamma_S||N_3||||Z_Y||`, and
  `||N_3||<=(20Gamma_S+2)||B_1||`.  Hence the **complete** convenient
  quartic residual has a perpendicular response column
  `O(||B_1||)` without a Markov spectral-gap inverse.
- L300's admissible cubic column differs from the convenient
  metric-only representative by a compatible Stein direction of size
  `O(||B_1||)`.  L301 transports that difference one order forward
  with the same scaling and zero lower corner.  Therefore L303's
  lower margin and retained upper descent survive the admissible
  selection.
- The bound is uniform on a fixed repeated-block neighbourhood.
  Along every L197 analytic rank-changing arc it forces L292's Smith
  valuations, so the selected response is analytically removable.
  Cubic-through-quartic bounded selection is closed.
- This is not an all-grade induction.  the historical step at that checkpoint is now to
  identify the arbitrary-grade relative-commutator recurrence which
  keeps every response commutator inside an earlier physical bridge;
  do not compute an isolated next grade.
- The exact 28-term certificate has zero residual words.  Nineteen
  tracked cases and 79 fresh unstructured, rank-changing, delayed,
  and reducible stress cases pass.  Dataset SHA-256:
  `881bd89b39171ab2603148d54092f605683ece63d65bb61e5394e0e7347f725c`.

## NEWEST (2026-07-26): L303 closes the complete partial-retightening quartic margin
- Scale L298's grade-one retightening by a fixed
  `0<theta<1` before inserting L300's cubic response and L302's
  quartic correction.  This avoids treating exact lower-tightness as
  a compulsory intermediate invariant.
- With `A=B_1*B_1` and `L=B_1B_1*`, the complete moving Schur
  differences begin
  `[c²]Delta_-=−theta A` and `[c²]Delta_+=+4theta L`;
  the upper cubic term is an L300 Markov response and the lower cubic
  term vanishes.
- Exact word reduction including every nonlinear Schur square gives
  the quartic lower debt
  `−J_theta`, where
  `J_theta=(4theta+theta²/2)A+
  (9theta/2−3theta²/4)A²`.
  The parallel analytic correction `G_S(VJ_theta V*)` restores that
  lower coefficient and leaves the physical upper class
  `−(8theta+4theta²)L−(8theta−4theta²)L²`.
- This upper cost is supported only on the already active `B_1`
  range.  The retained quadratic gain dominates it for
  `|c|<=1/2`, while the raw lower face retains the active-range margin
  `(1−theta)c²A`.  Thus `theta=1/2` is a robust fixed choice.
  `theta=1` remains useful only as the lower-tight associated-graded
  limit; no positivity claim is made for its uncomputed higher lower
  faces.
- This closes L289's **complete pointwise Schur-graph assembly through
  quartic order**.  At the L303 checkpoint the uniformly bounded
  response through rank changes was still open; L304 subsequently
  closes it at `theta=1/2` by a relative-commutator flux estimate.
  Only the arbitrary-grade promotion remains.
- The exact degree-two-in-`theta` residuals vanish at three exact
  scales and at additional rational audit scales.  The 22 tracked
  cases and 160 fresh unstructured, rank-changing, delayed, and
  weighted stress cases pass.  Dataset SHA-256:
  `7e8cd50fbcb95d99866089c89a0d9b0b684bc32403fa153ae3c822628d8c41d0`.

## NEWEST (2026-07-26): L302 computes the quartic's exact two-ended signature
- Cancel L299's cubic with the convenient metric-only representative
  `X_3=−G_S(N_3)`.  L301 guarantees that this choice cannot change the
  quartic commutant class.
- L298's grade-one frame simplifies exactly to
  `F_1=−SWB_1/2`, and its metric Green sum has only two terms.  The
  remaining nonlocal quartic cross telescopes through the dual Stein
  solution `2(S²+(S*)²)`.
- Exact partial-isometry word reduction gives
  `tr Q_4=(15/4)||B_1||²+(1/2)tr((B_1*B_1)²)`.
  Restricting to every L206 spectral reducing block upgrades this to
  `<Y,Q_4>=(15/4)tr(YB_1B_1*)+
  (1/2)tr(Y(B_1B_1*)²)`.
- Therefore
  `Q_4−[(15/4)B_1B_1*+(1/2)(B_1B_1*)²]`
  lies in the complete Markov response range.  However, an
  adversarial two-ended check gives the exact lower corner
  `K_4=6B_1*B_1+(9/4)(B_1*B_1)²`.  Therefore the attractive
  cancellation `−G_S(N_4)` has harmful lower motion `−K_4` and must
  not be called a completed lower-tight correction.
- Adding the compulsory parallel `+G_S(VK_4V*)` restores that lower
  corner.  The resulting upper commutant class is
  `(9/4)B_1B_1*+(7/4)(B_1B_1*)²`: a positive cost, but a pure
  prior-flag factor which vanishes on `ker B_1*`.  This closes the
  pointwise two-ended homology, not the full moving metric endpoint.
- The rank-changing analytic norm of the required response column is
  not yet controlled, and the complete mixed Schur-graph endpoint has
  not been assembled.  Prove a direct L281/L292 flux/valuation bound
  for the lower-neutral response and use L289 for the graph terms; do
  not use the pointwise pseudoinverse or recompute the cubic.
- The exact certificate and 22 tracked numerical records pass, plus
  100 additional scalar/weighted stress cases.  Dataset SHA-256:
  `c8e3b4edfadcb95d8d313f20e3d1056865983df8f51ac254d1d9974d524a6816`.

## NEWEST (2026-07-26): L301 transports every first successor into the Markov range
- For any fixed-base direction
  `X−S*XS=VC*+CV*`, its first moving successor is
  `T(X,C)=−(A_1*XS+S*XA_1+D_1C*+CD_1*)`.
- Exact reduction gives L300's universal metric expression plus the
  remainder `−2{(S*)²CV*+VC*S²}`.  If
  `Z_Y=(I−VV*)H_YV`, the full dual pairing is
  `<Y,tau>=tr(XK_Y)−4 Re tr Z_Y*(S*)²C`.
- Hence
  `|<Y,tau>|<=(20||X||+4||C||)||Z_Y||`, so every successor
  annihilates the full L206 cokernel and has a perpendicular response
  preimage bounded by `10||X||+2||C||`.  No pseudoinverse or Markov
  gap occurs.
- In particular, two bounded L300 cubic preimages change the quartic
  successor only by a bounded Markov response.  The retained
  commutant/separator component of the physical quartic is therefore
  selection-independent; simultaneous cubic/quartic optimization is
  unnecessary.
- L302 subsequently computes both endpoint classes: the attractive
  upper quotient comes with a compulsory lower corner, and lower
  neutralization leaves a prior-\(B_1\) upper cost.  the historical gate at that checkpoint is
  its bounded full-graph realization, not recomputation of the
  invariant.
- Exact rational algebra plus unstructured, reducible, and repeated
  monomial audits pass.  Dataset SHA-256:
  `24ca5eade9348ee61821c24604abd8dcbb8ad1156e62cb021f89e9e193343d7c`.

## NEWEST (2026-07-26): L300 gives the moving cubic a gap-free Markov preimage
- L299's cubic frame tangents reduce exactly to
  `D_1=−2(S*)²V` and `L_1=−2S²W`.  Using L298's Stein equation
  and its additional special identity `(S*)²F_1=0` eliminates the
  frame entirely and writes the cubic forcing as one universal
  expression linear in the already bounded metric direction `R_1`.
- If `H_Y−SH_YS*=WYW*` and
  `Z_Y=(I−VV*)H_YV`, cyclic reduction gives
  `<Y,Q_3>=tr(R_1 K_Y)`, where every term of `K_Y` contains `Z_Y`
  and `||K_Y||_F<=20||Z_Y||_F`.
- L281 identifies `2||Z_Y||²` with the exact Markov Dirichlet
  energy.  Hence `Q_3` annihilates the complete L206 cokernel and has
  a perpendicular response column with
  `||C_3||<=10||R_1||<=20 Gamma_S||B_1||`.  No Markov inverse,
  rank projection, or pseudoinverse occurs; `Gamma_S` is uniformly
  bounded on each fixed repeated-block neighbourhood.
- This closes L299's full cubic range and global bounded-energy debt,
  including reducible strata.  The special identity
  `(S*)²F_1=0` is not automatic for the selected correction frame.
  Inserting `C_3` creates a **quartic successor** through its cross
  with the raw first frame; the quadratic self-cost occurs only later
  at degree six.  Compute and localize the quartic successor before
  promoting the construction to L299's all-series recursion.
- Exact rational algebra and 12 independent unstructured,
  rank-chain, and reducible audits pass.  Dataset SHA-256:
  `1c8bf694c4f0783ce9d3b38da44e046ad3a5b9a406888c1a60ebe95db8eb064d`.

## NEWEST (2026-07-26): L299 isolates and falsifies the no-extra moving correction
- If `M−A*MA=DD*` is L227's moving raw pair and `(R,F)` is L298's
  weighted fixed-base correction, put `Delta=A−S`, `H=D−V`.
  The entire corrected residual is exactly
  `N=−Delta*RS−S*RDelta−Delta*RDelta−HF*−FH*−FF*`.
- Since `Delta,H=O(c)` and `R,F=O(c²)`, `N=O(c³)`.  This explains
  why L298 closes every first completely delayed even face.
- The cubic forcing is
  `N_3=−(A_1*R_1S+S*R_1A_1+D_1F_1*+F_1D_1*)` and has zero lower
  corner.
- It is not automatically a prior-flag factor.  On an exact rational
  six-state, two-copy partial isometry with rank-one `B_1`, its upper
  Stein response is trace-zero indefinite and has surviving-flag
  compression `2985984/54865681>0` on `ker B_1*`.
- Therefore the direct all-series insertion of L298 with no further
  homology correction is false.  L298 itself remains valid; L299
  provides the compact active forcing that must now be repaired by
  L280/L296 or L289/L292.
- The rational case and nine nearby noncommuting rank-chain cases
  pass.  Dataset SHA-256:
  `caea7e9741a255ba6319e573050be51cd951889f8755bc286156885cae9622a9`.

## NEWEST (2026-07-26): L298 retightens every transfer grade explicitly
- For `A_k=B_k*B_k`, the right-orbit Gram
  `Omega(A_k)=G_S(VA_kV*)` has physical endpoint pair
  `(+A_k,+4C(A_k))`.  Its negative removes L283's raw lower face but
  initially costs `−4C(A_k)` upstairs.
- Let `D_k` be one seventh of L212's balanced column.  Its lower
  endpoint is zero and its upper endpoint is
  `+4(C(A_k)−B_kB_k*)`.
- Therefore
  `R_k=−Omega(A_k)+G_S(VD_k*+D_kV*)`, with frame column
  `F_k=−VA_k/2+D_k`, obeys the exact Stein equation
  `R_k−S*R_kS=VF_k*+F_kV*` and has physical endpoint pair
  `(-B_k*B_k,-4B_kB_k*)`.
- Adding this to L283's completely delayed raw face
  `(+B_k*B_k,-12B_kB_k*)` gives the lower-tight effective face
  `(0,-16B_kB_k*)` in **every** grade.  This extends L214--L215
  without computing another isolated coefficient.
- Arbitrary real weighted sums also close exactly for the fixed base
  colligation; weights `|c|^(2k)` give the two oriented reflected
  Grams.  This does not yet close the full moving elliptic series:
  substituting the pulled operator creates nonlinear cross terms, and
  A179/A172 remain valid warnings.
- Exact algebra proves the result.  Twenty-seven independent
  noncommuting, completely delayed, and apex audits pass through
  grade six.  Dataset SHA-256:
  `ba5a8c971e377824280ad4afdc63b029d26f108534ec506a3328db4a4c8d1c7c`.

## NEWEST (2026-07-26): L297 aligns the raw and lower-tight endpoint ledgers
- At the repeated monomial apex, the boundary metric has endpoint
  pair `(lower,upper)=(+I,−4I)`.
- L283's raw L227 canonical repair has pair `(+I,−12I)`, whereas
  L279/L225's lower-tight effective metric has pair `(0,−16I)`.
  Therefore lower retightening contributes exactly `(-I,-4I)`.
- L212's target has lower face zero and upper face `−16I` because
  `G_L=0` at the apex.  It is therefore aligned with the effective
  endpoint, not the raw upper endpoint.
- L295 remains a valid obstruction to comparing the raw upper repair
  face directly with the lower-tight L212 target, but it does not
  obstruct a correctly normalized effective comparison.
- The raw L230--L289 route may retain the favorable lower budget and
  use L296 on the surviving mixed response.  The separate
  lower-tight route is now closed on every completely delayed first
  face by L298; its live debt is the moving-series/mixed-flag
  transport.  Every future face comparison must name its endpoint
  normalization.
- Exact endpoint-ledger audits pass in multiplicities one through
  four.  Dataset SHA-256:
  `2a493ad3fec9deb6de14b0a8c65943d09eafb82ac3fea34bcf58dc78470575dc`.

## NEWEST (2026-07-26): L296 polarizes the gap-free flux estimate
- For an arbitrary copy matrix `X`, put
  `F_k(X)=sym(B_kX*)−Phi(sym(X*B_k))`; this is one eighth of L280's
  exact physical endpoint response.
- Channel adjointness gives
  `<Y,F_k(X)>=Re tr X*(YB_k−B_kPhi*(Y))`.
  L280's Dirichlet form is the sum of the squared defects on the
  right, so
  `|<Y,F_k(X)>|<=||X||sqrt(<Y,(I−Phi Phi*)Y>)`.
- Summed and flagged versions are literal.  A family `X_j` is
  controlled by `(sum_j||X_j||²)^(1/2)`, and `X=PB_k` recovers L294's
  flagged quadratic flux.
- Therefore the raw repair route does not need the mixed
  L285/L289 remainder to equal a scalar `G_k`.  Modulo favorable
  Grams and prior-flag factors whose compression is zero, it is enough
  to derive polarized response columns and prove total energy
  `<=C²||U*B_k||²`; L282 then follows with constant `C`.  L288 is in
  the zero-compression factor class and need not itself be a Markov
  response.
- Exact single, summed, and flagged audits pass on two noncommuting
  bistochastic channel families in multiplicities two and three.
  Dataset SHA-256:
  `58381726281a8102eb3a53b334b3d0ac11b7dcfaccb3eae3cc64d9750d344c48`.

## NEWEST (2026-07-26): L295 obstructs a raw-upper/lower-tight comparison
- At a repeated monomial apex, `B_L=U` is unitary and
  `Phi(K)=UKU*`, so `G_L=Phi(I)−I=0`.
- L283's raw L227/L285 canonical-repair upper endpoint is exactly
  `−12I` and retains lower endpoint `+I`.  L212's lower-tight
  canonical target is exactly `12I−28I=−16I` and has lower endpoint
  zero.  Comparing only those upper faces produces the strictly
  positive `4I` discrepancy.
- Therefore the instruction to dominate the **raw upper** face by
  L212's **lower-tight** target is false.  L297 subsequently shows
  that the missing `−4I` is exactly the upper part of the
  `(-I,-4I)` lower-retightening increment.  The correctly normalized
  completely delayed face is subsequently closed by L298; only its
  nonlinear moving-series/mixed-flag transport remains open.
- L294 itself survives unchanged and has the branch-neutral
  consequence
  `<Y,−aB_kB_k*−bG_k>_+
   <=|b|||B_k||sqrt(<Y,(I−Phi Phi*)Y>)`
  for `a>=0`.  The live repair-branch question is whether the mixed
  transport remainder, after removing L283's favorable `−12` Gram,
  has this flux form or an equally direct observability-defect
  pairing.
- Exact terminal-unitary audits pass in multiplicities one through
  four.  Dataset SHA-256:
  `d983482024c5fd0b998ea0028edafcd56d8336c01d6a4fe4aa493ca10ba0264d`.

## NEWEST (2026-07-26): L294 closes the gap-free energy estimate for the canonical flux target
- Put `A_Y=Phi*(Y)` and
  `G_k=Phi(B_k*B_k)−B_kB_k*`.  Channel adjointness gives the exact
  scalar pairing
  `<Y,G_k>=−Re tr B_k*(YB_k−B_kA_Y)`.
- L280's Dirichlet form is the sum over all grades of the squared
  defects `||YB_n−B_nA_Y||²`.  One Cauchy--Schwarz step therefore
  gives
  `|<Y,G_k>|<=||B_k|| sqrt(<Y,(I−Phi Phi*)Y>)`
  with no inverse or spectral-gap constant.
- The flagged version uses
  `G_(k,P)=P{Phi(B_k*PB_k)−B_kB_k*}P` and replaces `B_k` by
  `PB_k`; this is exactly L282's quantitative inequality.
- Consequently L212's canonical candidate
  `12B_kB_k*−28Phi(B_k*B_k)
   =−16B_kB_k*−28G_k`
  satisfies L282 with constant `28`, and L212's explicit response
  changes it exactly to the favorable Gram `−16B_kB_k*`.
- **Critical scope:** this does not prove any physical face is the
  canonical candidate.  L212 stated that identity conditionally, and
  A179 disproves raw simultaneous superposition.  L295 additionally
  disproves comparing the raw L227/L285 upper face directly with this
  lower-tight target.  L297 shows that their apex gap is precisely
  lower retightening, and L298 constructs the exact all-grade
  lower-tight complete-delay face without assuming L212's candidate
  base.  Use the general flux bound on the raw branch's mixed
  remainder, or transport L298's weighted fixed-`S` correction
  through the full moving series.
- Exact full and flagged pairings, Dirichlet identities, positive-face
  cases, and correction residuals pass in multiplicities two through
  four.  Dataset SHA-256:
  `1ea9886ad59dc1c0facfa93f43cc17b00a67db57c67a512cceca9621b5fa8110`.

## NEWEST (2026-07-26): L293 identifies the exact first Markov gap in Schur coordinates
- Right-normalize L218's transfer so its terminal unitary is `I` and
  write `Gamma_j=epsilon Delta_j+O(epsilon²)`.
- L218's reflected tangent has coefficient pairs
  `Delta_j` at degree `j` and `-Delta_j*` at degree `2L-j`.
  The two inner Parseval identities determine the Hermitian part of
  the order-two terminal coefficient; its skew part cancels from
  `I−Phi Phi*`.
- The resulting exact first face is
  `sum_j([Delta_j*,[Delta_j,H]]+
  [Delta_j,[Delta_j*,H]])`.  On Hermitian `H`, its Dirichlet form is
  `2 sum_j ||[H,Delta_j]||_F²`.
- Along any analytic arc, take the first nonzero Schur order.  The
  same formula applies to its leading jets, so a closing Markov gap
  has a fixed quadratic exponent and kernel equal to their common
  commutant.  L281's physical state-column norm has exactly twice
  this energy.
- L294 subsequently proves the exact nonlinear transfer-flux pairing
  and Cauchy--Schwarz bound.  L295 disproves comparing the raw upper
  L227/L285 face directly with L212's lower-tight target, while L297
  shows the effective apex is aligned after lower retightening.  What
  remains on the raw repair branch is to pair only the mixed transport remainder after
  separating L283's favorable Gram; L293 alone does not identify that
  remainder.
- The exact all-basis checker passes in lengths three through five,
  multiplicities two and three.  Dataset SHA-256:
  `2e5c20fabc9fed72a752a802cee12374e751327585144afb04883b93f6a17059`.

## NEWEST (2026-07-26): L292 makes bounded endpoint selection an exact valuation problem
- Along any one-variable analytic arc, write the joint earlier
  transfer row in Smith form
  `Bcal=P diag(s^nu_1,...,s^nu_rho,0) Q` and put
  `Ghat=P^-1 G P^-*` for L289's complete endpoint residual.
- A locally bounded analytic factor `G=X Bcal*+Bcal X*` exists if
  and only if: the generic-kernel block of `Ghat` is zero;
  active-active entry `(a,b)` vanishes to order at least
  `min(nu_a,nu_b)`; and active-kernel entry `(a,k)` vanishes to order
  at least `nu_a`.
- Sufficiency is constructive and triangular; every division is by
  a certified Smith power.  No pseudoinverse or rank projection is
  used.
- **New guardrail:** pointwise flag zero is not enough.  The scalar
  family `Bcal=s²,G=s` satisfies the pointwise condition, including at
  the enlarged origin flag, but forces `X=1/(2s)`.  A mixed
  `(nu_1,nu_2)=(1,3)` cross can force `s^-2`.
- By L197 curve selection, a hidden local blow-up would appear on an
  analytic arc.  A178's exact remaining selection target is therefore
  to derive these finite valuation inequalities from the physical
  arbitrary-grade endpoint recurrence, then retain the even
  finite-jet margin.  L292 does not prove those physical
  divisibilities.
- Exact positive and obstruction cases pass with zero symbolic
  residuals.  Dataset SHA-256:
  `a9c28eec243b7854bb7f06959f359a1e9bcb4a40894e6a65aeee4c5eab80136d`.

## NEWEST (2026-07-26): L291 removes nonendpoint state algebra from the flag debt
- Block elimination writes every positive-interior Schur state as
  `L_H*diag(U(H),D)L_H`.  Matching two interiors by the canonical
  positive congruence gives
  `Htilde=R*HR+diag(U(Htilde)−U(H),0)`, where
  `R=[[I,0],[K,C]]` is analytic and locally bounded.
- Thus all state motion except the literal Schur-endpoint residual is
  an endpoint-quotient-fixing graph gauge.  The gauge cannot create a
  rank-change singularity.
- For `Bcal=[B_1 ... B_r]`, a Hermitian endpoint residual `Q` has
  zero compression to the joint surviving flag exactly when
  `Q=X Bcal*+Bcal X*`.  This is an if-and-only-if statement at each
  fixed matrix, not merely L290's sufficient state-channel theorem.
- **Critical scope:** the fixed-rank converse formula uses
  `Bcal`'s Moore--Penrose inverse only to prove existence.  It may
  blow up across rank changes and must not be inserted in the metric.
  L292 now replaces it along every analytic arc by an exact
  Smith-valuation criterion; the physical valuation bounds remain
  open.
- **Strategic correction:** do not require the large transported
  metric witnesses themselves to lie in a delay word ideal.  Use
  L285/L194 to construct admissible states, then factor only L289's
  complete endpoint residual and retain the even direct-Gram/Schur
  budget.
- Ranks one through four and independently moving complex positive
  interiors pass the normal-form, flag, and converse audits.  Dataset
  SHA-256:
  `c7eea10f5d4e58b25a13fe67040d671816245dc0be8b5fafccc276a61b997bb0`.

## NEWEST (2026-07-26): L290 identifies the exact nonlinear transfer-flag invariant
- Put `D_j=ES^jF` and `B_j=W*(S*)^jV`.  Every endpoint graph from
  L289 satisfies `FJ_H=W`, hence the exact channel identity
  `D_jJ_H=VB_j*`.
- Therefore a Hermitian state motion
  `DeltaH=sum_j(A_jD_j+D_j*A_j*)` has the complete nonlinear endpoint
  motion
  `DeltaU=sum_j{(J_Htilde*A_jV)B_j*+
                 B_j(J_H*A_jV)*}`.
  This formula includes both moving graphs and every Schur cross.
- Compression to any surviving flag with `B_j*U=0` is exactly zero.
  The corresponding norm is bounded by the sum of the graph-factor
  norms times `||B_j||`, so L289's bounded analytic graphs introduce
  no additional loss.
- **Scope guard:** the relevant invariant is hereditary/one-sided:
  the channel generator must remain exposed against a graph column.
  Arbitrary words `XD_jY` need not have this property, so L290 does
  not revive the raw ideal invariance disproved by L285.
- **historical gate at that checkpoint (narrowed by L291):** factor each complete L289
  mixed-graph pairing directly in the endpoint flag ideal.  A
  hereditary state representative remains sufficient, but is not a
  required induction invariant.  L290 alone does not preserve the
  even direct Gram/Schur budgets or prove finite-jet positivity.
- Generic and rank-changing audits through defect multiplicity five
  pass.  Dataset SHA-256:
  `c0bd78fd90f4d1c7260cff29fb3a4c4d7bf8b10eb3497d6365bbc8061fef60a8`.

## NEWEST (2026-07-26): L289 makes nonlinear endpoint transport exact and the target finite-jet
- For two Hermitian upper-gap series, define the endpoint graph
  `J_H=F−P(PHP)^−1PHF`.  Exact block elimination gives
  `U(Htilde)−U(H)=J_Htilde*(Htilde−H)J_H`, with the adjoint
  orientation equal to the same Hermitian difference.
- Coefficient `n` is the triangular convolution
  `sum_(i+d+j=n) Jtilde_i* DeltaH_d J_j`.  This includes all graph
  motion and nonlinear Schur crosses without enumerating them
  separately.  L288's quartic/sextic three-term `Omega_6` is exactly
  its first nontrivial instance.
- The graph columns stay uniformly analytic and bounded near the
  repeated equality point.  Schur shorting therefore introduces no
  independent coefficient-growth obstruction.
- For fixed Crabb length `L`, choose a polynomial free row only
  through degree `2L`.  If its finite upper jet has a terminal lower
  bound `delta c^(2L)B_LB_L*`, L201 gives
  `sigma_min(B_L)>=beta>0` nearby and L194 gives a uniform
  `O(c^(2L+1))` tail.  The full endpoint is positive for
  `c<delta beta^2/M`.
- **Strategic correction:** A178 still needs an arbitrary-grade
  formula, right-ideal factorizations, Schur budgets, and uniform
  margins, but only finitely through grade `L` for each fixed
  repeated block.  Infinite summability uniformly over all grades is
  not a valid stop condition.
- Exact rational and formal-series residuals are zero.  Dataset
  SHA-256:
  `7fbf04c4b81828d756fabfe3cacd9ba05ec968287f2e583e9d35baddabc76781`.

## NEWEST (2026-07-26): L288 transports the trace-carrying sextic face through the surviving flag
- Keep L287's nonzero even trace rather than trying to erase it.  An
  explicit 141-word Hermitian witness (`l1=881/2`) and 116-word
  perpendicular correction (`l1=261`) solve the transported sextic
  Stein equation exactly.
- With L287's 40-word parallel correction, the recomputed 142-word
  sextic lift has coefficient bound `1021/2` and lies in the
  triple-delay ideal.  It vanishes when `B_1=B_2=B_3=0`.
- The raw metric corner is not enough to test L234.  After adding the
  first nonlinear quartic Schur cross, the complete sixth upper
  change is a 40-word polynomial (`l1=96`) with the exact
  factorization `A_1(ESF)+(A_1ESF)*`; `A_1` has ten words and
  `l1=24`.  Hence the motion vanishes already on `ker B_1*`, and
  L234's direct and effective sextic faces survive unchanged on
  `ker B_1* intersect ker B_2*`.
- **Scope guard:** the separately recomputed nonlinear lower Schur
  coefficient is not globally unchanged at order six.  L288 is a
  Stein/upper-flag transport certificate, not global two-endpoint
  invariance, an arbitrary-grade parity theorem, or L289's terminal
  finite-jet certificate.  Translate any final selection through
  L194's exact lower-tight chart.
- Exact residuals are zero.  The tracked audit covers unstructured,
  rank-changing multiplicities three through five, and triple-delay
  samples; an additional 35-case stress pass succeeded.  Dataset
  SHA-256:
  `aaf15d56428af913f025be1d0132d256b72a892907f6f69bdba8881781a5dfa1`.

## NEWEST (2026-07-25): L287 proves endpoint-null transport must stop at the next even face
- Recompute the sextic forcing from both L286 factor changes and both
  transported metric coefficients, subtract L284's sextic gauge, and
  perform L285's compulsory lower elimination.  The resulting
  `H_6^0` has 260 exact words (`l1=756`) and zero lower corner.
- The fixed lower parallel correction has 40 words (`l1=48`) and lies
  in the triple-delay ideal, so lower normalization is not the
  obstruction.
- The obstruction is upper trace.  On the exact rational scalar-
  defect partial isometry
  `[[0,7/9,-4/9,-4/9],[0,-4/9,1/9,-8/9],
    [0,-4/9,-8/9,1/9],[0,0,0,0]]`,
  with spectral radius `2/3` and `B_1=−4/9`, rational Stein inversion
  gives upper endpoint `151040/177147`.
- L204 proves every remaining perpendicular endpoint response has
  trace zero; at scalar defect it is identically zero.  Therefore no
  globally endpoint-null sextic completion exists.  This is an exact
  falsification, not evidence from an unsuccessful word search.
- **Strategic correction:** the recursion must cancel removable odd
  faces as in L286, but at even order retain the transported metric
  trace and combine it with L283/L234's direct transfer Gram and the
  preceding odd Schur square.  The physical even face is a budget,
  not gauge.
- Dataset SHA-256:
  `372921b75a883e3066e1d910082e5cb805d4af7b2e37d2999a9863daa3f5bb1c`.

## NEWEST (2026-07-25): L286 solves the first coupled normalized-gauge successor
- Starting from L284's normalized quartic column `C_4−N_4`, the old
  normalized quintic `C_5−N_5` must be changed because of L285's raw
  transport.  The corrected column is
  `Ctilde_5=C_5−N_5+K_5`.
- An explicit 58-word Hermitian witness `Z_5^tr` (`l1=70`) and
  48-word perpendicular correction `K_5` (`l1=127/2`) satisfy the
  complete transported Stein coboundary with zero exact residual.
  Both endpoint compressions of `Z_5^tr` vanish.
- Exact quotient reduction gives `K_5∈I_2` and
  `Ctilde_5 V*∈I_2`.  The new quintic lift has 36 words and
  coefficient bound `112`, improving L284's individually normalized
  bound `122`.  Hence it vanishes whenever `B_1=B_2=0`.
- Because the degree-four and degree-five metric changes are both
  endpoint-null, and the base has no degree-one endpoint cross, the
  full lower and upper endpoints through degree five are unchanged.
  L232's positive quartic face and L233's fifth partial-flag
  cancellation therefore survive every existing rank-changing case.
  Independent full-series matrix recomputation agrees to `4.3e−13`.
- Dataset SHA-256:
  `ba1ad08b60e250339e910d35de04d700c0b3ffd4e64ff9a7313b68dbb3c35acf`.
- **Scope guard:** this proves one genuine recursive step, not the
  arbitrary-grade recurrence.  L287 subsequently shows that the
  endpoint-null convention itself fails at the recomputed sextic
  coefficient; do not reuse L234's column unchanged or try to erase
  the physical even face.

## NEWEST (2026-07-25): L285 makes gauge transport triangular and rejects raw ideal invariance
- If `Dtilde=D+C` and `Mtilde=M+X` satisfy the same operator
  Stein-factor equation, exact coefficient subtraction gives
  `X_n−S*X_nS=R_n+C_nV*+VC_n*`, where `R_n` depends only on earlier
  gauge and metric coefficients.  The transport is triangular.
- Put `K_n=V*R_nV`.  The parallel column `−V K_n/2` removes the
  lower corner.  The remaining perpendicular column must solve
  exactly one L204 endpoint equation; L280 identifies its range with
  the Markov coboundaries.  Thus A178 and L282 are two certificate
  languages for the same coefficientwise homology, not independent
  recursions.
- The tempting stronger invariant is false.  For L284's quartic
  gauge, the raw first successor
  `R_5=−N_4D_1*−D_1N_4*−T_1*Z_4S−S*Z_4T_1` has 100 reduced words.
  Its quotient modulo the first delay ideal has 28 nonzero words
  (`l1=60`), and its `Q` corner still has 20.  Hence raw transport
  cannot simply be declared right-ideal divisible.
- Exact recurrence, Hermiticity, and lower-corner residuals are zero.
  A single-delay matrix with `B_1=0`, `B_2≠0` has nonzero raw
  successor endpoint norm `3.53e−2`, confirming that the next
  homological solve is physical rather than bookkeeping.  Dataset
  SHA-256:
  `54dc8db20c3b301f7d75978ccb126ab75741417f531ec1db4cd87eb871d7a350`.
- **historical gate at that checkpoint (updated by L291):** prove that L285's endpoint target
  has either a bounded polynomial endpoint factor with uniform
  finite-jet margins through the terminal grade or a uniformly
  bounded-energy Markov preimage.  Do not compute another isolated
  grade.  L286 subsequently supplies the first such polynomial
  preimage, through degree five only.  Raw two-sided state-ideal
  membership is not enough and, by L291, is not required.

## NEWEST (2026-07-25): L284 makes every existing preparation delay-normalized modulo an endpoint-null gauge
- If `Z=Z*`, `Y=Z−S*ZS`, and
  `EZE=FZF=EYE=QYQ=0`, then `P=QYE=N V*`,
  `Y=P+P*`, and the Stein response of `N` is exactly `Z`.  Both
  endpoint compressions vanish, so `N` is a global endpoint-null
  gauge without a response inverse.
- Explicit witnesses `Z_4,Z_5,Z_6` satisfy this certificate and
  match the delayed quotient of the stored L232--L234 columns.
  Subtracting the gauges produces endpoint-equivalent order-four,
  order-five, and order-six columns that vanish under complete delays
  of lengths one, two, and three respectively.
- Every exact Hermiticity, corner-support, forcing, perpendicularity,
  and delay-ideal residual has zero words.  Independent unstructured
  and delayed matrix audits agree.  The normalized coefficient
  l1-bounds are `25`, `122`, and `851/2`.  Dataset SHA-256:
  `9dad8b6817a2e3b56b50f7089d63252e057828d818700c4ae6d5b004ece62d17`.
- **Scope guard:** the three normalizations cannot simply be applied
  simultaneously to the stored series.  Replacing `C_4` changes the
  fifth and later metric forcing, and replacing `C_5` changes the
  sixth and later forcing.  A direct rank-chain scope guard leaves a
  fifth-flag residual `5.99e−3` after naïve simultaneous reuse.
  L285 now shows that even the raw first transport is not delay-ideal
  divisible; the correct next object is its endpoint homology target.

## NEWEST (2026-07-25): L283 closes the all-grade complete-delay matrix face
- L258 makes the edge-deleted closed-return coefficient the first
  normalized final Schur face, not merely a scalar determinant
  coefficient.
- L247's active metric reinsertion is the exact negative of L279's
  closed-return deletion response.  They cancel before taking a
  trace, leaving L278's intact face
  `(ES^kF)(S*)^k+S^kF(S*)^kE`.
- L246 conjugates that face to
  `E_1F_(k−1)+F_(k−1)E_1`.  This proves L228's former all-grade
  candidate, closing A171's covariance debt and A170's
  complete-delay face.
- Its Stein endpoint is `2B_kB_k*`.  L227's single analytic
  canonical repair therefore has exact first active faces
  `+B_k*B_k` below and `−12B_kB_k*` above on every completely delayed
  stratum.
- **Scope guard:** A172 remains load-bearing.  When earlier `B_j` are
  nonzero but singular, the unmodified canonical repair has an
  indefinite mixed cubic compression on the next partial flag.
  L197/L220 do not remove that physical mixed term by themselves.
  L230--L234 prepare the first three grades; A178's arbitrary-grade
  right-ideal recursion is still open.  L282's quantitative
  off-commutant estimate is the alternative bounded-selection route.
- Exact response/conjugacy algebra passes through grade twelve,
  numerical endpoints through grade six, and the independent full
  physical word audit through grade five.  Dataset SHA-256:
  `e45b478ae10270de59fc443ec076f8715d5dabe7654a26ff6f6df6de1842ed57`.
- **historical gate at that checkpoint:** promote the exact delayed boundary value and the
  cubic-through-sextic right-ideal certificates into an arbitrary-
  grade bounded preparation.  Do not compute another isolated grade.

## NEWEST (2026-07-25): L282 converts bounded selection to one off-commutant flux inequality
- For `K=I−Phi Phi*`, the least energy of a flagged correction is
  exactly
  `sup_(Y>=0){<Y,U*EU>−16<UYU*,K(UYU*)>}`.
  Optimizing the scale of `Y` gives the homogeneous quotient
  `(1/64)sup <Y,U*EU>_+²/<UYU*,K(UYU*)>`.
- Therefore a uniform inequality
  `<Y,U*EU>_+ <= gamma ||U*B_k|| sqrt(<UYU*,K(UYU*)>)`
  immediately yields correction energy
  `<=gamma²||U*B_k||²/64`, exactly the scale needed after L281.
- This uses Slater from L222+L279 and `K^(1/2)`, not a pseudoinverse.
  L279 proves the zero-denominator case.  L283 subsequently found the
  stronger exact complete-delay matrix face, but A172 shows that it
  does not settle mixed partial flags.  This sufficient inequality
  remains a live alternative to A178's explicit recurrence.
- The depolarizing two-copy audit compares the closed form with
  independent primal and dual SDPs.  Dataset SHA-256:
  `559e4a7ee92179c2cca7085e9942433158bbb700902988bb3a865d7451369a33`.

## NEWEST (2026-07-25): L281 square-roots the Markov response by the physical observability column
- Let `Hcal_H−S Hcal_H S*=W H W*` and
  `R_H=(I−VV*)Hcal_HV`.  Then the balanced metric column `2R_H`
  realizes the complete response exactly:
  `M_T(2P_metric^(1/2)R_H)=8(I−Phi Phi*)H`.
- Stein balance gives the exact energy identity
  `||2R_H||_F²=2<H,(I−Phi Phi*)H>`.  Its polarized form identifies
  the Markov Dirichlet pairing with the physical observability-column
  pairing.
- Therefore the closing Markov gap has already been absorbed by the
  physical state synthesis.  Do not seek a bounded right inverse or
  try to bound `H`; the only historical selection question at that checkpoint is whether the
  L279 physical endpoint admits a Poisson correction with Dirichlet
  energy bounded by its first active transfer energy along every
  L197/L220 flag.
- The shared L280/L281 checker passes the endpoint, state-energy,
  polarized-column-collapse, bistochastic, Dirichlet, and repeated-
  apex audits.  Dataset SHA-256:
  `18a8c075ceb72d50871d7ff2d1ae443e8e5732dc2ce0e3f6957dc41901dc7fda`.

## NEWEST (2026-07-25): L280 replaces rank-changing range inversion by a Markov Laplacian
- Polarizing L212 for an arbitrary copy multiplier `X` gives the exact
  response
  `8{sym(B_k X*)−Phi(sym(X*B_k))}` with an explicit state column.
- Setting `X_k=H B_k` and summing all grades yields
  `M_T(C(H))=8(I−Phi Phi*)H`.  The state column is analytic, uses no
  projection or pseudoinverse, and vanishes on the repeated monomial
  apex.
- The exact Dirichlet identity
  `<H,(I−Phi Phi*)H>=sum_k||H B_k−B_k Phi*(H)||_F²`
  identifies the fixed space with L206's self-adjoint colligation
  commutant.  Therefore
  `ran M_T=ran(I−Phi Phi*)`: this is the complete response range, not
  just L212's one coboundary direction.
- L222's semidefinite problem is now the copy-space Markov Poisson
  inequality `E+8(I−Phi Phi*)H<0`.  L279 proves its strict
  fixed-point separator condition.
- **Superseded selection formulation:** this originally left a bound
  on the long state preimage `C(H)`.  L281 replaces it by the shorter
  physical column `2R_H`, whose norm is exactly the Dirichlet energy.
  The remaining debt is bounded-energy Poisson selection.
- The exact checker passes arbitrary polarized multipliers, both
  bistochastic identities, the summed response, the Dirichlet form,
  and apex vanishing.  Dataset SHA-256:
  `18a8c075ceb72d50871d7ff2d1ae443e8e5732dc2ce0e3f6957dc41901dc7fda`.

## NEWEST (2026-07-25): L279 closes the all-grade delayed volume flux
- L277--L278 already make every coefficient below the active
  edge-deleted face zero operator-wise.
- Deleting `X_k=[c^(2k)]R` changes L258's exact closed-renewal defect
  by `P(X_k−SX_kS*)P`; the final-row return term has two zero
  background cross factors, so no loop contribution is omitted.
- L273's exact metric formula and the complete delay give
  `P(X_k−SX_kS*)P ~tr 2H_k`.  L278's intact face is independently
  `~tr 2H_k`, where
  `H_k=Q_(k+2)−(k+2)Q_1+(k+1)I`.
- Hence the full face is `~tr 4H_k`,
  `[c^(2k)]log V=4||B_k||_F²`, the full dual-Schur trace is
  `2||B_k||_F²`, and
  `tr E_(2k,eff)=−16||B_k||_F²`.
- The exact checker regenerates physical faces through grade four
  and algebraic cyclic identities through grade twelve; dataset
  SHA-256 is
  `521e3e3f1129eaa36832738fb0ab86ffd73c766437bcbf75626522c9e28276b0`.
- **What remains:** L222 pointwise feasibility is closed, but it does
  not provide a bounded analytic correction as the commutant rank
  changes.  Use L197/L220's analytic Schur flags to solve that
  selection problem before merging the repeated elliptic and normal
  faces.  Do not continue the now-closed A194 fan/volume grind.

## NEWEST (2026-07-25): L278 kills every lower even face before edge deletion
- Keep the physical metric coefficient at degree `2h` and impose only
  delays below `h`.  First-visit grouping, L244's exact coisometry,
  and L263's first-face transport identify the full frontier with the
  embedded tail frontier, shifted by two degrees.
- The universal five-word face factors exactly as
  `(ES^hF)(S*)^h+S^hF(S*)^hE`.  Adding the next complete delay kills
  the whole operator coefficient.
- Together with L277, every coefficient below the active edge-deleted
  face is now zero in arbitrary grade.  This is stronger than L276's
  lower-radial premise.
- **Sole support gate:** evaluate the active metric-deletion response
  modulo trace/cyclicity through `Q_(k+2)`.

## NEWEST (2026-07-25): L277 factors every lower odd face through the next delay
- With lower delays through `h−1`, the degree-`2h+1` unweighted
  frontier has one universal twelve-word form (`h>=2`) and an
  eight-word exceptional base at `h=1`.
- First-visit deflation embeds the preceding tail frontier and leaves
  an exhaustive eight-word new-crossing layer.  After the complete
  physical renewal is summed, its coefficients are fixed `±2`; the
  raw direct/immediate/later pieces still must not be separated.
- The complete frontier factors explicitly through
  `D_h=ES^hF` and `D_h*=F(S*)^hE`.  Adding the next complete delay
  therefore kills it identically in every grade.  Since the retained
  metric is even in `c`, every lower odd weighted closed-defect face
  vanishes as well.
- L278 subsequently proves every lower even face vanishes as well.
- **Still open:** active associated cyclic radiality through
  `Q_(k+2)`; do not compute more lower grades.

## NEWEST (2026-07-25): L276 replaces exact support by three-test cyclic forcing
- If the active full-minus-tail closed-defect face is cyclically
  radial through `Q_(k+2)`, L252 gives
  `tr Delta_k=p_k(1)n−p_k'(1)m+p_(k,k+2)||B_k||²`.
- Arbitrarily long exact axes with `B_k=0` force
  `p_k(1)=p_k'(1)=0`.  On the active monomial, L241+L247 give the
  edge-deleted value `4m` for both the full grade and its one-layer
  tail (L256 closes the grade-one endpoint), forcing
  `p_(k,k+2)=0`.
- For lower faces, no sharp index bound is required: each formal
  coefficient is a finite two-sided radial polynomial, and all of
  its radial projections are linearly independent on a sufficiently
  long shift whose exact defect is zero.  Every lower radial
  coefficient therefore vanishes operator-wise, making the active
  log determinant linear.
- Iteration ends at L256 and yields the desired edge-deleted volume
  coefficient `4||B_k||²`.
- L277 subsequently proves the lower odd cancellation.
- L278 subsequently proves the lower even face is zero, not merely
  radial.
- **Still open:** prove the active cyclic-radial envelope from
  L243/L251/L258/L272.  Do not grind the stronger literal no-extra-
  word theorem unless it falls out automatically.

## NEWEST (2026-07-25): L275 propagates coefficient one through the full return fan
- Define
  `A_(k,j)=(S*)^j S^(k+2) (S*)^(k+2−j)`,
  `B_(k,j)=S^(j+2)(S*)^(k+2)S^(k−j)`, and
  `G_k=S(S*)^(k+2)S^(k+2)S*`.  L275 proves that the active
  unweighted renewal has coefficient one on all
  `A_(k,0..k)`, `B_(k,0..k−1)`, and `G_k`.
- The load-bearing recursion is fan-restricted, not a forbidden
  whole-series equality.  Once L258 has summed the final-row loops,
  every inherited fan path is the exact tail path with its clean
  layer removed and coefficient degree lowered by two.  L274 adds
  `A_(k,k)`, `B_(k,0)`, `G_k`, and `−SG_kS*`; the tail fan supplies
  every inherited word and `+SG_kS*`, so induction closes.
- The cancellation is physical and cannot be split.  For example,
  the grade-four raw direct/immediate/later differences for
  `A_(4,0)` are `−32,+64,−32`; only their sum is zero.
- Exact quotient arithmetic verifies the fan recursion and finds no
  additional nonradial word through grade six.  Tail fan embeddings
  pass through grade twelve.
- **Still open:** L276's weaker lower radial and active cyclic-radial
  support envelopes.  Do not promote the finite exact-support audit
  or L275's coefficient theorem to a literal no-extra-word theorem.

## NEWEST (2026-07-25): L270 forces the radial moments from exact axes
- Assume L269's proposed full-minus-tail face
  `Delta_k=P_k(Q_0,...,Q_k)+R_2−R_1+G_k−SG_kS*`, with
  `G_k=S(S*)^(k+2)S^(k+2)S*`.  The deep pair has zero trace in every
  grade because its boundary flux is
  `||S^(k+2)S*V||_F²=0`; this follows from `S²S*V=0`.
- On arbitrarily long scalar monomial shifts, L241 makes the physical
  boundary metric agree with the exact zero-residual axis through the
  requested face, while L247 says active-edge deletion has zero trace
  response because `B_k=0`.  Hence
  `n p_k(1)−p_k'(1)−1=0` for every sufficiently large `n`, forcing
  `p_k(1)=0` and `p_k'(1)=−1`.
- Under a general complete matrix delay, L252 then gives
  `tr P_k=m`, the fixed `R_2−R_1` pair gives `−m`, and the deep
  divergence gives zero.  Thus the structural normal form alone
  implies exact scalar one-delay covariance; no formula for its
  irregular radial coefficients is needed.
- This is a conditional reduction, not A194: the normal form is proved
  only at grade two and audited through grade five.  L277 subsequently
  closes the lower odd vanishings; resume by proving the remaining
  support/decomposition and lower even structure directly, not by
  computing more grades.
  `proof/repeated_crabb_radial_moment_forcing.md`.

## NEWEST (2026-07-25): L271 reduces deflation to one full-face support theorem
- Write the proposed full active face as
  `D_k=U_k(Q)−R_1+R_k−4R_(k+1)+H_k`, where `H_k` consists of the
  two crossed remote words, `R_(k+2)`, and `G_k`, and the radial
  coefficients obey `u_k(1)=0`.
- For the embedded tail `T=S(I−F)=S²S*`, the complete-delay quotient
  gives exactly `Q_j(T)=Q_j+R_1−I` and `R_j(T)=R_(j+1)`.
  The zero radial symbol cancels the apparent extra `R_1`.
- Three tail remote words become the identical full remote words;
  the fourth becomes `SG_kS*`.  Subtracting the tail template from
  the full template therefore gives L269's normal form verbatim.
- The exact audit checks every generator and remote-word identity
  through grade twelve and regenerates SHA-256
  `db416897f988a131b83aa61a0a881d264053ac7ca0066c81f7e420132a5fe17c`.
  This validates the arbitrary-grade algebraic implication, not the
  open full-face premise.
- Resume by proving only the one-system support formula from
  L242--L245's first-reflection jet, L258's renewal, and L244's
  coisometric background.  Do not independently expand the tail.
  `proof/repeated_crabb_full_face_deflation.md`;
  `experiments/repeated_crabb_full_face_deflation.py`.

## NEWEST (2026-07-25): L272 closes the leakage-ideal placement gate
- With `Hhat=H L`, the shifted transfer-Hankel channel satisfies
  `L*Hhat=Hhat L` and `Hhat*Hhat=P_(K_B)=I−L_B`.  Every innermost
  paired doubled-Hardy crossing is therefore exactly a copy-free shift
  word minus one leakage insertion.
- Iterating this identity normal-orders every closed copy-dependent
  word into the two-sided ideal generated by `L_B`.  L258 Section 3,
  rather than a formal `H=0` substitution, cancels the completely
  copy-free renewal before coefficient extraction.
- For complete delay `k=r+1>=2`, one nonconstant L243 inverse-kernel
  selection first costs degree `2r+1`; two selections lie strictly
  above the active degree `2r+2`.  L245's Green columns, direct/metric
  coefficients, and L258's renewal cannot lower this valuation.
  Hence every active scalar term is
  `tr{L_B Phi(L,L*) L_B}` with copy-scalar `Phi`.  L256 supplies the
  `k=1` base.
- This closes only ideal placement/no-unpaired-transfer.  It neither
  proves lower vanishings nor bounds the unilateral boundary depth of
  `Phi`.  the unique historical theorem at that checkpoint remains L271's full-face support
  formula: prove that the copy-scalar return telescopes to the radial
  packet, fixed right-orbit endpoints, and four remote words, with
  deeper residue occurring only as L269's trace-zero divergence.
- The deterministic one-to-four-channel audit checks the shifted model
  projection, complement, intertwining, and sandwich orientations.
  Tracked SHA-256:
  `68e5a2bdc7a535f94238fa3c5f64451abea54330441b3c00e915699791342ea6`.
  `proof/repeated_crabb_hankel_leakage_ideal.md`;
  `experiments/repeated_crabb_hankel_leakage_ideal.py`.

## NEWEST (2026-07-25): L273 proves the metric fan telescope
- L219's grade-`h` metric coefficient is radial:
  `X_h=R_h−R_(h+1)+sum_(d|h)(−1)^(h/d)(Q_d−Q_(d+1))`, with top
  boundary `Q_(h+1)−R_(h+1)`.
- Below the first missing delay, formal inverse recurrence preserves
  radiality and gives top boundary
  `[Y_j]_top=−Q_(j+1)+R_(j+1)`.  The exact load-bearing reduction is
  `(Q_a−R_a)(−Q_b+R_b)=Q_min+R_min−2I` whenever
  `ES^(a+b−2)F=0` is available.
- Deleting the active coefficient gives
  `[q^k](R°)^−1=−sum_(h=1)^(k−1)X_hY_(k−h)`.  At each split, the
  first unavailable delay leaves exactly
  `Q_(h+1)R_(k−h+1)` and `R_(h+1)Q_(k−h+1)`.  Therefore the metric
  inverse has coefficient `−1` on all `2k−2` interior fan words and
  no other nonradial support.
- The observed unweighted return has a `2k+2`-word fan.  L273 proves
  that, if this unweighted support is established all-grade, the
  metric cancels every interior path and leaves exactly L271's four
  endpoints: the two crossed words, `R_(k+2)`, and `G_k`.
- This does not prove the unweighted full-fan premise or lower
  vanishings.  Those are now the sole nonradial support gate; the
  metric convolution must not be recomputed.
- The exact audit checks every inverse convolution and radial product
  through grade sixteen.  Tracked SHA-256:
  `bbf8afbe36d039f98a3c8c5df6dc257d34f4e354fe70e126265618ea89d7acff`.
  `proof/repeated_crabb_metric_fan_telescope.md`;
  `experiments/repeated_crabb_metric_fan_telescope.py`.

## NEWEST (2026-07-25): L274 fixes the four remote coefficients
- For the unweighted recursion
  `Delta_k^U=[c^(2k)]U_k−[c^(2k−2)]U_tail`, the two crossed remote
  words each split into direct, immediate-return, and later-return
  coefficients
  `1+4(−1)^k`, `−4(−1)^k`, and `0`.  Their physical sum is exactly
  one in every grade.
- The deep word `G_k=S(S*)^(k+2)S^(k+2)S*` comes only from the full
  retained direct Gram with coefficient one.  Its shifted negative
  copy comes only from subtracting the independently balanced tail
  direct Gram.  No final-row return contributes to either word.
- The proof is an all-grade run/source extraction from L125/L242's
  direct-map filtration and L273's lower inverse-metric top boundary.
  It keeps the doubled endpoint inside the complete immediate
  entry/exit product, so it does not repeat L249's false detached-port
  inference.
- L274 does not exclude additional unweighted words or prove lower
  vanishings.  Resume by classifying the remaining support; the four
  remote coefficients are banked and must not be recomputed.
- The exact direct/immediate/later-return audit passes in grades two
  through four with tracked SHA-256
  `4d360bee00b01b1f4e658c8ac4a634d78f0ea519714cd392b11dc3d22176b8d2`.
  `proof/repeated_crabb_remote_endpoint_coefficients.md`;
  `experiments/repeated_crabb_remote_endpoint_coefficients.py`.

## NEWEST (2026-07-25): L256 proves the universal relative base
- At relative grade one, deleting the active `c²` metric coefficient
  leaves `R°=I+O(c⁴)`.  Expanding L250's mass with the complete
  theta/ODE jet gives the exact cyclic pieces
  `36Q_1−12Q_3−24I` and `16Q_3−48Q_1+32I`.
- Neither piece is the answer separately.  Their sum is
  `4(Q_3−3Q_1+2I)`, whose trace is `4||B_1||_F²` by L252.
  This holds for arbitrary matrix-valued transfer data, not only
  scalar or monomial channels.
- The live theorem is now sharply isolated as **associated delay
  covariance**: after removing `k−1` clean Hardy layers, the first
  active scalar coefficient must be the universal grade-one response.
  L235 disproves whole-series tail equality, so only the associated
  coefficient may be transported.
- As part of the single-agent checkpoint, L248 and L249 were rerun
  from scratch and reproduced their tracked hashes; Git remains a
  clean linear history with no lost or conflicting artifact.
  `proof/repeated_crabb_grade_one_volume.md`.

## NEWEST (2026-07-25): L257 makes delay covariance one model-space recursion
- In the left Hardy model `K_(zB)`, all state data are now explicit:
  `F=P_0`, `E=H*P_0H`, `S*=L*`, `S=P_KL`,
  `R=H*D_RH+P_KD_LP_K`, and
  `Xi=S+c(I+P_0)L*(I+E)`.  Here
  `H_(n,j)=B_(n+j)*`.  L250's scalar mass is therefore a canonical
  functional of the matrix-inner transfer, with no state realization
  or metric root left to choose.
- Complete delay is the literal orthogonal split
  `K_(zB)=direct-sum_(j<k−1) z^j C^m direct-sum
  z^(k−1)K_(z Btilde)`.
  The all-grade gate is exactly
  `[c^(2k)]mu°_(B,k)=[c^(2k−2)]mu°_(Btilde,k−1)`.
  Once proved, iteration ends at L256 and gives `4||B_k||²`.
- The existing exact cyclic engine was refactored to evaluate an
  arbitrary partial-isometry corner.  A new audit constructs the
  original and embedded-tail masses independently and reduces their
  cyclic difference to zero in grades two through five.  The raw
  quotient differences have `36,82,159,268` words.  This is stronger
  than comparing two copies of the known target but remains finite
  evidence, not the missing proof.
- Floating continuation also shows covariance through the following
  odd coefficient and failure at the next even coefficient, matching
  L235's warning: the associated recursion is sharp and must not be
  upgraded to whole-series equality.
  `proof/repeated_crabb_left_model_volume.md`;
  `experiments/repeated_crabb_volume_delay_covariance.py`.

## NEWEST (2026-07-25): L258 resums every final-row visit into a closed return
- Write `Z=A(R°)^−1A*`, `R°=F direct-sum R_P`, and
  `Z_ret=PZP+PZF(I_F−FZF)^−1FZP`.  A direct block Schur determinant
  gives, through the active face,
  `V=det_P(I+D_ret)` with `D_ret=I_P−R_PZ_ret`.
- `D_ret` is similar to the normalized final-defect Schur residual,
  so it has exactly the same vanishing order and coefficient trace.
  Conditional on A194's still-open lower vanishings, its active trace
  is `[c^(2k)]log V`.  Expanding the sole inverse pairs every entrance
  to the final row with an exit and sums all intervening final-row loops.
- This closes the algebraic one-sided-port ambiguity: L243's later
  `B#(x)` column cannot be valued in isolation.  It does **not**
  prove the lower vanishings, the desired `+4`, remove future transfer
  coefficients from a closed path, or establish L257's one-delay
  recursion.  The remaining historical gate at that checkpoint is now the matrix-inner
  autocorrelation/Hardy-index evaluation of these closed returns.
  `proof/repeated_crabb_output_renewal_volume.md`.

## NEWEST (2026-07-25): L259 localizes every future coefficient to one leakage row
- Differentiating the model kernel gives
  `[rho]K_B(rho,x)=xI−B_1B#(x)`.  Since
  `P_(K_B)=I−T_BT_B*`, the future series `B_1B#(x)` is exactly
  the first row of the complementary causal Toeplitz leakage.
- Under complete delay `B_1=...=B_(k−1)=0`, the first nonzero
  leakage row is `(T_BT_B*)_(k,j)=B_kB_j*` for `j>=k`; it vanishes
  to the left of the diagonal and its diagonal trace is
  `||B_k||_F²`.
- This does not prove the coefficient law or its lower vanishings.
  It replaces the vague future-coefficient cancellation debt by one
  exact statement: at relative active order, L258's physical return
  must take four times the diagonal of this row and must not shift an
  off-diagonal block back into the trace.  Numerical checks already
  rule out the stronger operator covariance, so the target is
  genuinely trace-only.
  `proof/repeated_crabb_future_leakage_row.md`.

## NEWEST (2026-07-25): L260--L261 fix the allowed future-row closure
- A natural two-orientation lift of the ellipse series was tested in
  the exact grade-one cyclic algebra.  It fails as a proof device:
  both degree-one orientation sectors and all three degree-two
  sectors are individually nonzero; only their physical sums cancel.
  Theta/direct-map orientations therefore cannot be separated before
  the coisometric/Schur operation.
- The correct grouping is the whole L259 leakage row.  Since
  `L_B=T_BT_B*` is a projection, complete delay gives the exact
  Parseval closure
  `tr(P_kL_B²P_k)=tr(P_kL_BP_k)=||B_k||_F²`, or
  `sum_(j>=k)B_kB_j*B_jB_k*=B_kB_k*`.
- The live metric theorem is now sharper: show that L258's first
  physical entry/return pair supplies the ordinary Hardy adjoint row,
  with no theta/metric weight between it and L259's leakage row.
  Then L261 closes all future coefficients at once and the doubled
  remote amplitude supplies four.  Lower vanishings remain part of
  the same proof.
  `proof/repeated_crabb_orientation_split_obstruction.md`;
  `proof/repeated_crabb_leakage_row_parseval.md`.

## NEWEST (2026-07-25): L262 removes the hidden half-line metric ambiguity
- Normalize L244's half-line to the coisometry
  `U=D_R^(1/2)A_infD_R^(−1/2)` with defect column `e`.  The Wold
  synthesis `W*(z^jv)=(U*)^jev` has orthonormal columns.
- At `c=0` this synthesis is the identity Hardy basis.  Since it is
  `I+O(c)`, it is formally invertible; orthonormality therefore makes
  it formally unitary and `UW*=W*L*`.
- Consequently an actual entry/return pair through this normalized
  half-line has the ordinary unweighted Hardy norm.  There is no
  residual theta weight capable of spoiling L261's Parseval closure.
- This is not yet A194: the remaining proof must place L243's
  zero/one inverse-kernel sector and L245/L251's complete physical
  port pair between `W` and `W*` inside L258, after summing the two
  orientations required by L260.  That placement must also prove the
  lower vanishings and retain the doubled remote factor on both sides.
  `proof/repeated_crabb_half_line_wold_return.md`.

## NEWEST (2026-07-25): L263 removes the fixed-port/defect-graph ambiguity
- For a rank-`m` background Gram `dd*` and fixed Schur port `F=JJ*`,
  the first new perturbation `c^s Delta` has exact Schur face
  `T Delta T*`, where
  `T=P−Pd(J*d)^−1J*`.
- L240 has `d(0)=J_0`, so `T(0)=P`.  L262's Wold unitary is also
  `I+O(c)`.  Therefore the coefficient at its first new degree sees
  only the constant ordinary Hardy compression; motion of the theta
  defect graph, the fixed port, and the Wold coordinates is at least
  one order too late.
- This does not assert a positive Gram for the physical face (finite
  audits disprove that stronger operator claim).  It isolates the
  remaining A194 numerator: after summing both orientations, show
  that its retained **trace** is four times L261's complete
  leakage-row norm, with any other active operator terms trace-null.
  `proof/repeated_crabb_schur_graph_transport.md`.

## NEWEST (2026-07-25): L264 restores the right-half-line metric for free
- If two L250 metrics isolate the final port, agree below degree `d`,
  and differ at degree `d` by `X=PXP`, their scalar-mass response is
  `tr((F−E)X)=−tr(EX)`.  An arbitrary retained coefficient is
  therefore not free.
- The positive right wandering orbits `(S*)^jES^j`, `j>=1`, are
  orthogonal to `E`, so their retained active coefficient does have
  zero response.
- Hence replace the whole-coefficient deletion by
  `Rtri=Rcirc+q^kP([q^k]R_right)P`.  This restores L240/L244's
  complete `D_R` coefficient while leaving the active left orbit,
  final corner, and cross row deleted.
- The left orbit is not gauge-free: its `E`-trace is `||B_k||²`.
  The incorrect stronger retained-gauge claim was rejected by a
  random-matrix audit before commit.
  `proof/repeated_crabb_retained_metric_gauge.md`.

## NEWEST (2026-07-25): L265 turns the target into a boundary-flux law
- Put `E_j=(S*)^jES^j` and
  `X_k=sum_(j=0)^k(k+1−j)E_j`.  Exact telescoping gives
  `Q_(k+2)−(k+2)Q_1+(k+1)I=X_k−S*X_kS`.
- Its trace is the final-boundary flux
  `tr(FX_k)=sum_(j=1)^k(k+1−j)||B_j||²`, hence exactly
  `||B_k||²` under complete delay.
- A194 is therefore equivalent to showing that L258's complete
  active numerator is `4(X_k−S*X_kS)` modulo trace-null terms.
  This is compatible with the observed indefinite active operator;
  no positive Gram factorization is required.
- In L264's gauge the right-half-line bulk is coisometric.  The live
  calculation should use L242's universal `D_0,D_1,D_2` reflected
  jet and add L245/L251's chain flux before taking this divergence.
  `proof/repeated_crabb_triangular_stein_flux.md`.

## NEWEST (2026-07-25): L266 rigidifies Laurent and shallow shift returns
- For the leakage projection `L_B=T_BT_B*`, shift invariance of
  `BH²` and complete-delay triangular support give
  `P_k L_B Psi(L)L_B P_k=psi_0 P_kL_BP_k` for every copy-scalar
  Laurent polynomial or coefficientwise stabilized series.
- More generally, a shift word reduces to `L^a(L*)^b`.  If
  `min(a,b)<=k`, its unilateral boundary correction is supported
  strictly below the active row, so the same statement holds with the
  constant coefficient of the bilateral symbol.  Depth `k+1` can
  fail; copy-scalarity alone is not enough.
- Hence the physical return need not be the identity and need not be
  estimated term by term.  Every nonconstant forward or backward
  Hardy shift misses the first delayed diagonal exactly; the trace is
  simply `psi_0||B_k||²`.
- A nonmonomial, noncommuting Potapov-product audit passed for copy
  sizes `2,3,4`, delays `1,2,4`, and shifts through absolute degree
  four.  The exact invariant-subspace proof is independent of this
  audit.
- L269 shows how this theorem must be used physically: the fully
  assembled return is not wholly shallow.  First separate its
  depth-`k+1` trace-zero Stein divergences, then prove that the
  remainder is copy-scalar and either Laurent or of boundary depth at
  most `k`.  Compute that remainder's bilateral-symbol constant as
  four while proving all earlier coefficients vanish.
  L242's universal reflected jet and L256's grade-one result are
  checks, not substitutes for that derivation.
  `proof/repeated_crabb_toeplitz_return_rigidity.md`.

## NEWEST (2026-07-25): L267 exposes exact lossless defect transport
- For a unitary colligation
  `U=[[A,B],[C,D]]`, the Redheffer feedback
  `Phi(Z)=A+BZ(I−DZ)^−1C` satisfies the exact noncommutative identities
  `I−Phi*Phi=C*(I−Z*D*)^−1(I−Z*Z)(I−DZ)^−1C`
  and its final-defect analogue.
- Thus a genuine lossless entry/return network cannot silently insert
  an unknown positive tail metric.  It transports the tail defect
  through one explicit port column.  This is classical
  conservative-systems theory, not a novelty claim.
- A deterministic audit used sixty noncommuting feedback matrices in
  three rectangular block layouts; both identities passed with
  largest residual `1.433e−15`.  The tracked dataset hash is
  `972f9f2516914af778a637a192f3dfb0bc42bac83cda3201893be00a4e1809a3`.
- L267 does **not** identify the repeated-Crabb physical map with such
  a feedback.  L268's exact rank test subsequently rejects the direct
  single-port identification.  Analytic functional calculus does not
  generally commute with Redheffer feedback, and L251's paired channel
  alone is insufficient.  Do not promote A213's candidate
  trace-locality synthesis unless its ideal-placement and depth/row
  stop conditions are closed.
  `proof/lossless_redheffer_defect_transport.md`;
  `experiments/lossless_redheffer_defect_transport.py`.

## NEWEST (2026-07-25): L268 rejects the one-shot lossless shortcut
- On the scalar grade-two monomial, exact rational series give
  `[c^4]D_ret=diag(2,2)` on the full retained two-dimensional space.
  Every earlier coefficient vanishes.
- The independently balanced deflated tail has
  `[c^2]D_ret_tail=diag(0,4)`, with every earlier coefficient zero.
  The associated scalar traces are both four, but the face ranks are
  two and one.
- If L267 transported the tail quotient as the sole non-background
  defect, L263's graph quotient would give
  `K_full=Q(c)*K_tail Q(c)`.  Positivity of the tail's leading scalar
  face and the two-degree delay force `Q(c)=cQ_1+...`, so the full
  first face would be `4Q_1*Q_1`, of rank at most one.  This
  contradicts `2I_2`.
- Therefore the complete physical delay removal is not one
  single-port lossless feedback of the deflated quotient.  L267
  remains correct and may organize subchannels, but a larger network
  would need an extra independently evaluated defect channel.  The
  scalar theorem remains viable and is now confirmed to require the
  trace-only A213 cancellation.
- The exact checker uses no floating tolerance; tracked SHA-256:
  `fb53ef7433874e6be9887f68519bd467bf9c403e96864288b36d6aa2b4da4f65`.
  `proof/repeated_crabb_lossless_feedback_obstruction.md`;
  `experiments/repeated_crabb_lossless_feedback_obstruction.py`.

## NEWEST (2026-07-25): L269 finds the missing deep trace divergence
- Compute L258's operator-valued closed return before tracing and
  compare a full grade-two delay with its embedded grade-one tail.
  The exact difference is
  `−2I+5Q_1−R_1−3Q_2+R_2+G_2−SG_2S*`, where
  `G_2=S(S*)^4S^4S*`.
- The first five terms have zero trace by L252.  The two deep terms
  reduce to the same cyclic representative `Q_4`, so they also have
  zero trace.  This proves the grade-two scalar recursion while
  retaining the operator mismatch seen in L268.
- On the unilateral half-line, `G_2` has boundary depth three.  It is
  nonzero in the `B_1=0` operator quotient.  Therefore A213's original
  demand that the complete return have depth at most `k=2` is false.
  L266 must be applied only after deep terms are split off as
  trace-zero divergences.
- Exact grades `2..5` all collapse to
  `radial P_k+G_k−SG_kS*`, with only `7,8,9,8` words and zero cyclic
  difference.  This is much smaller than the pre-renewal expansion
  and suggests the correct all-grade recurrence, but remains finite
  evidence beyond grade two.
- The next proof must derive this split from L243/L251/L258 and prove
  the two scalar radial moment identities `P_k(1)=0`,
  `P_k'(1)=−1`.  The tracked exact dataset hash is
  `0398cee7f133b13ab4bc1f332e5d8e030f8ce86e92929e547cc9b30cce48a988`.
  `proof/repeated_crabb_closed_return_recursion.md`;
  `experiments/repeated_crabb_closed_return_recursion.py`.

## NEWEST (2026-07-25): L253 identifies the target as Toeplitz leakage
- If `T_B` is the causal Toeplitz multiplier of the matrix-inner
  transfer and `P_k` selects Hardy rows `0,...,k`, then L252's radial
  trace is exactly
  `||P_kT_BP_k||_HS²=sum_(h=1)^k(k+1−h)||B_h||²`.
  Since `T_BT_B*=I−P_(K_B)`, it is also
  `(k+1)m−tr(P_kP_(K_B)P_k)`.
- Under complete delay, the finite Toeplitz window has only one
  nonzero block, `B_k` in its bottom-left corner.  Thus the live
  target is intrinsically `4||P_kT_BP_k||_HS²=4||B_k||²`, on the
  same Hardy/model space used by L221, L236, and L244.
- What remains open is the physical identification, not the scalar
  target: move L251's closed paired channels into the Hardy frames,
  cancel the full-space coisometric background with L244, and prove
  that the residual trace is four times this finite-window leakage.
  `proof/repeated_crabb_toeplitz_window_energy.md`;
  `experiments/repeated_crabb_toeplitz_window_energy.py`.

## NEWEST (2026-07-25): L254 identifies the exact complementary channel
- Let `O_Lhat=L*O_L` be L236's left Hardy analysis map with its
  constant row removed.  The characteristic-kernel identity gives
  `O_Lhat O_Lhat*+T_BT_B*=I`.  The shift is essential: unshifted
  `O_L` realizes `K_(zB)`, whereas the Toeplitz complement uses
  `K_B`.
- Therefore L253's target is exactly the missing energy of the
  shifted left-Hardy/model channel, not merely an abstract model-space
  codimension.
- Across all grades its generating series is
  `sum_(k>=1)q^k||P_kT_BP_k||_HS²
   =(1−q)^−2 tr sum_(h>=1)q^hB_h*B_h`,
  the double Abel transform of L201's reflected Faber Gram.
- The live physical statement is now precise: show that L250's
  final-row whitening, after L244 cancels the half-line coisometric
  background and L251 keeps the doubled port paired, removes the
  `O_Lhat O_Lhat*` channel and leaves four times its `T_BT_B*`
  complement.  L254 proves the complement and generating identities,
  not that volume equality.
  `proof/repeated_crabb_hardy_faber_complement.md`;
  `experiments/repeated_crabb_toeplitz_window_energy.py`.

## NEWEST (2026-07-25): L255 puts the physical pencil in that channel
- The two missing Hardy intertwinings are
  `O_RS*=L O_R−L H P_0O_L` and
  `O_LS=L O_L−L H*P_0O_R`, with matching exact actions of `E,F`.
- Substitution gives an exact `2 x 2` representation of the full
  physical ellipse pencil.  Its diagonal background is precisely
  `Xi_inf direct-sum Xi_inf*`; its off-diagonal blocks are linear in
  the transfer Hankel cross Gram `H`; and its only diagonal
  correction is the ordered quadratic
  `c H P_0L*H*P_0`.
- Thus all colligation dependence is now in the same Hankel channel
  whose complement L254 identified, while L236's metric weights are
  fixed and diagonal in the two frames.  No free state words or
  one-sided transfer rows remain in the formulation.
- The open A194 step is the universal second channel response:
  apply the analytic functional calculus and L250 whitening to this
  block pencil, use L244 to cancel the diagonal half-line background
  and L243 to exclude higher reflection powers, then identify the
  quadratic trace as `4T_BT_B*`.
  `proof/repeated_crabb_doubled_hardy_pencil.md`;
  `experiments/repeated_crabb_hardy_two_frame.py`.

## PREVIOUS (2026-07-25): L252 reduces the target to three radial traces
- For `Q_j=(S*)^jS^j`, two elementary defect telescopes prove
  `tr Q_j=n−jm+sum_(h=1)^(j−2)(j−h−1)||B_h||²` in every grade.
  Under a complete grade-`k` delay this makes
  `4tr(Q_(k+2)−(k+2)Q_1+(k+1)I)=4||B_k||²` exactly.
- Exact cyclic word reduction of L250's edge-deleted normalized mass
  gives precisely that three-term expression through grades one to
  six.  Imposing the complete-delay ideal during multiplication keeps
  the quotient faces to `12,38,79,151,249,392` words.  Separately
  reducing the raw Stein face and final-row whitening correction shows
  that every nonradial cyclic class cancels coefficientwise between
  them; every total face has only the three radial trace classes.
  This is exact finite evidence for L244's coisometric cancellation,
  not an arbitrary-grade proof.
- The live A194 obligation is now the explicit cyclic congruence:
  reduce L251's closed paired channels, modulo the lower-delay trace
  ideal, to `4(Q_(k+2)−(k+2)Q_1+(k+1)I)`.  L244 supplies the
  word-free cancellation.  No symmetric one-sided endpoint valuation
  or positive leakage factor is being assumed.
  `proof/repeated_crabb_cyclic_radial_volume.md`;
  `experiments/repeated_crabb_cyclic_radial_volume.py`.

## PREVIOUS (2026-07-25): L251 pairs every analytic chain port
- Integrating L245's two exact Green columns before coefficient
  extraction gives a common retained operator `Q_(r,j)^f` for the two
  opposite blocks of every analytic functional calculus:
  `e_j*f(Xi)J=mu_j c^(r−j)W*Q` and `J*f(Xi)e_j=QW`.
- Every ordinary chain row has `mu_j=1`; the remote physical row alone
  has `mu_0=2`.  Thus the doubled amplitude survives the complete
  theta/Riemann map, not merely the resolvent, while remaining tied to
  the same tail channel on both orientations.
- In a scalar trace-log, block paths are closed.  The apparently
  unweighted `B#(x)` row from L243 therefore cannot be costed alone:
  it stays inside a closed product of `Q`-channels and intervening
  chain blocks.  Entrance and exit rows may differ, so this rules out
  detached one-sided bookkeeping but does not by itself close A192's
  symmetric-valuation gap.
- The paired cycle is not yet evaluated.  Insert L251 into L250, use
  L244 on the common background, and then reduce the closed tail
  cycles by matrix-inner autocorrelation.  The coefficient `+4`
  remains open.
  `proof/repeated_crabb_analytic_chain_port_balance.md`;
  `experiments/repeated_crabb_analytic_chain_port_balance.py`.

## CHECKPOINT REVALIDATION (2026-07-24): prior audit complete; frontier advanced
- The requested A164/next-actions/restart/recursion/structured-matrix/
  ledger audit was already committed as `fed9eb3` at the L234
  checkpoint and has been rechecked after L248.
- `A164` now occurs once and means L221's proved tail covariance.  The
  disproved raw L212 all-series superposition is uniquely `A179`, and
  every inbound warning points there.
- A177's sextic computation is no longer in flight: L234 proved it.
  Rolling the restart packet back to A177 would discard fourteen later
  lemmas.  The current sharply scoped gate is A194/L248's whitened
  volume coefficient `+4||B_k||_F²`; A178 remains the subsequent
  arbitrary-grade selection/convergence debt.
- The Gohberg--Semencul/displacement-rank/Schur--Levinson novelty
  classifications for L187/L193/L218 remain respectively (ii), (ii),
  and (i).  The lemma-ID gaps remain intentional and are stated at the
  top of `LEMMA_LEDGER.md`.

## NEWEST (2026-07-25): L250 removes the metric roots from A194
- Under complete delay, deleting the **whole** active coefficient
  `[c^(2k)]P_bl` makes the final defect an identity metric block with
  zero cross row through degree `2k`; the next possible contamination
  is degree `2k+2`.
- Therefore L248's normalized initial-defect mass is, through the live
  face, exactly
  `−m+tr((R−A*FA)^−1(R−A*RA))` in balanced coordinates.  No formal
  metric square root remains, while the inverse keeps the complete
  final-row/Schur whitening.
- Sylvester moves the same determinant to the output variable
  `Z=AR^−1A*`.  This is the preferred A194 interface because L244 is
  exactly the half-line identity `A_inf D_R^−1 A_inf*=D_R^−1`, and
  L243--L245 give the finite reflected state lift entering `Z`.
- The open coefficient remains `+4||B_k||²`.  Generic grades one
  through five reproduce it after all coordinate changes.  The
  whitened face itself is generally indefinite, so do **not** promote
  the doubled channel to a positive leakage Gram; use coisometric
  trace conservation and matrix-inner autocorrelation.
  `proof/repeated_crabb_edge_deleted_balanced_volume.md`;
  `experiments/repeated_crabb_edge_deleted_balanced_volume.py`.

## NEWEST (2026-07-25): L249 locks the doubled edge to the full balance
- On the exact length-three delayed monomial channel, replace only the
  physical doubled reverse edge by the multiplier `1+lambda`, while
  freezing the boundary metric.  The degree-two final-dual Schur
  residual is exactly
  `diag(1−lambda²,2lambda−2)`.
- The lower cancellation survives iff `lambda=1`, the physical
  doubled edge.  Thus the slogan “remote amplitude two, so energy
  four” is not a proof when the port is detached from the theta,
  direct-map, and metric network.
- This does not damage A194's coisometric-port route.  It prescribes
  its safe form: keep the full physical network, insert L243--L245
  into L248's initial-defect mass, apply L244's coisometry to the
  complete background, and extract the remote energy only after that
  cancellation.
  `proof/repeated_crabb_terminal_multiplier_rigidity.md`;
  `experiments/repeated_crabb_terminal_multiplier_rigidity.py`.

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
- The same volume is `det(I+D_eff)`, and
  `tr D_eff=−m+tr((I−C*FC)^−1(I−C*C))`.  This moves the open
  coefficient to the normalized mass of the **initial** defect,
  exactly the orientation in which L240 supplies its coisometric
  half-line Gram.
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
- This removes all higher delay corrections from the historical gate at that checkpoint.  It
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
  normal form, the historical gate at that checkpoint is now a single
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
  exposes a uniform right-ideal recurrence and finite-jet margins
  sufficient through the terminal transfer.  L289 proves an infinite
  summability condition is unnecessary.  L283 has since closed
  A171/L228's all-grade complete-delay anticommutator.  L284
  additionally normalizes each
  stored later column by an endpoint-null delay gauge.  The remaining
  recursion is the triangular transport of each such metric-gauge
  change into L285's endpoint homology.  Raw ideal invariance is false;
  the required right-ideal representative must be selected by solving
  that homology.  This is not another delay-face calculation.
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
  indefinite elliptic Stein slack.  The remaining historical task at that checkpoint is now
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
- The remaining historical gate at that checkpoint is no longer compatibility: it is a bounded
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
- The remaining historical issue at that checkpoint is bounded analytic selection through the rank
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
  equality manifold; its metric-flag tube is the historical frontier at that checkpoint.
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
  fallback.  The remaining historical task at that checkpoint is the nonlinear singular blow-up.
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
  At this stage the precise remaining historical gate at that checkpoint was to evaluate
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
  remaining historical gate at that checkpoint is that Weierstrass preparation and the full Stein
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
  `(I+cS^2)/(I-cS^2)`.  The remaining historical gate at that checkpoint is now solely to insert
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
  the historical calculation at that checkpoint is the elliptic deficit generated by nonzero `u`,
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
  two-copy terminal chart is now locally complete-`2`.  the next historical gate at that checkpoint is uniformizing
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
  in the flat copy core.  The remaining historical issue at that checkpoint is genuinely uniform: weighted sequences where
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
- The common flat mode creates no new irreducible equality mechanism.  The remaining historical issue at that checkpoint is
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

## Current next actions (Epoch 6, refreshed after A299)

1. **LIVE: sign the boundary-shape Schur form in the conformal
   chart.**  A292 proves that
   `2(hat(s_C)(2),...,hat(s_C)(n))` is an exact coordinate in
   `C^(n−1)` for L343's quotient.  L345 says phase covariance is
   equivalent exactly to `W^T J_phi W=0` for
   `W=J_phi^−1S_y*+iJ_phi^−1S_x*`.  L346 proves that the observed
   `tril Y_1(W)=0` is equivalent to only
   `diag Y_1(W)=0` and `Y_1(W)_(L0)=0`.  L349--L350 already close
   both endpoint diagonal entries.  Prove the remaining `n−1` Euler
   identities from L344's fixed formula.  L347 already identifies
   the corner as
   `4f'(0)hat(omega_phi s_W)(−n)/gamma`; prove this last negative
   Hardy coefficient vanishes.  L348 then turns the diagonal debt
   into root tracking.  L349--L350 close the endpoint sum and split,
   hence the repeated-root mean equation once the corner is zero.
   Derive only
   `v_j=tr(P_jG)` at each simple root from
   `Hhat^(dagger/2)(I−zeta R*)^−1 Ccal*(H_Cy−s_Cy)`,
   and use L336's exact paired endpoint polarization to derive dual
   isotropy.  Then derive a positive Gram for the negative Hermitian
   `(n−1)`-square.  Do not infer the phase identity from nearly
   paired weighted-coordinate eigenvalues or from the falsified
   diagonal-plus-rank-one inverse ansatz, do not try to cancel the
   cyclic corner with zero motion, and do not re-expand or
   re-prove the strictly negative `(n−2)^2` disk block.  L342's two
   endpoint residual squares must remain intact, and L340/L341's
   rank/nullity law remains the falsification gate for any proposed
   factor, not a second task.
   A299 additionally forbids assigning even one fixed complete copy
   of the inverse model-kernel Gram after the first pivot; keep the
   endpoint residuals and fixed port coupled.
2. **Falsify each proposed closed identity before proving it.**  The
   general jet supplies complete normal matrices and zero Schur
   blocks through dimension nine; any candidate must reproduce
   L336's endpoint sum, L340's loss Gram, L342's universal free
   Hessian and exact square gap, and the active/null singular-value
   split.  The right and left signs separately, pointwise
   `t>=delta²/2`, and the raw partial-isometry residual are already
   falsified.  Do not run another isolated finite-size SDP campaign.
3. **Promote only the proved strict quotient to a local patch.**  Use
   the fixed-model support gap, finite Blaschke compactness, and
   Gau--Wu equality uniqueness exactly as in L333.  Keep collisions
   of model zeros as separate ramified charts already exemplified by
   L192/L329.
4. **Then handle sums and global compactness.**  Extend the local
   theorem across arbitrary finite model strata and their reducing
   direct sums before asking whether every sharp sequence approaches
   that full equality set.  The CP/H-r programs remain parked, not
   concurrent frontiers.

## Superseded pre-L331 Gau--Wu actions (historical)

1. **HISTORICAL: attack the arbitrary Gau--Wu disk-model neighbourhood,
   starting with `G_a`.**  Recenter the disk-manifold and
   scalar-extremal coordinates at a fixed `0<|a|<1`.  First derive the
   exact first support/Jensen face and the transverse scalar norm
   Hessian.  Test every proposed sign against rational values of `a`
   before opening a general proof.
2. **Exploit the noncollision support gap.**  For fixed `a`,
   `spec Re(e^{-i theta}G_a)={1,-1,a cos theta}` gives a uniform
   top-support gap `1-|a|`.  Use ordinary analytic support and Riemann
   charts; do not import the ramified repeated-Crabb flag unless a
   genuine multiplicity collision occurs.
3. **Then extend across the finite Blaschke model manifold and direct
   sums.**  Keep arbitrary zeros rather than replacing the inner
   function by a monomial.  L330 proves that this enlargement is
   logically prior to global sharp-sequence classification.
4. **Only after those local patches, return to global
   classification.**  At that stage test whether every exact or
   asymptotically sharp pair approaches the full Gau--Wu disk-model
   manifold.  Keep fixed-dimension compactness and dimension-uniform
   estimates as distinct steps.

## Superseded immediate post-L329 actions (historical; corrected by L330)

1. **HISTORICAL: classify or falsify the global sharp-stratum reduction.**
   Determine whether every exact scalar ratio-two pair, or every
   asymptotically sharp sequence after affine and numerical-range
   conformal normalization, has a disk/Gau--Wu reducing model limit
   covered by L192/L329.  Test the statement first against banked
   small exact witnesses and numerical extremals; do not infer it
   from monomial equality or from `w(Phi_A(A))<=1`, which has not been
   proved in general.
2. **If the reduction survives, isolate its exact equality
   mechanism.**  Use the finite Blaschke extremal and singular-vector
   stationarity, not the stronger completely bounded similarity
   route.  For each fixed dimension, combine an equality-stratum
   classification with compactness to obtain a strict complement.
   Record separately the still-required dimension-uniform step.
3. **If a new sharp stratum appears, stop and chart it.**  Do not
   force it into the repeated-Crabb coordinates.  State the exact
   limiting operator/function/equality data and open only the minimal
   new local theorem it requires.
4. **Keep CP/H-r as ordered fallbacks.**  L319 continues to block the
   Schur-only complete-similarity shortcut.  Do not reopen a
   complete-matrix normal flag merely to strengthen L329, which is
   already closed in the scalar setting.

## Superseded pre-L329 next actions (historical; do not resume)
1. **COMPLETED BY L318: fixed-half-scale margin assembly after L317.**
   The remainder of this item is the derivation record, not a live
   instruction.
   L306 sums every rooted successor to
   `q=−Psi_S(Delta*XS+S*XDelta+Delta*XDelta)
      −(C*H+H*C+C*C)`, and L307 proves the formula recurs after a
   preceding forcing is canceled.  L308 disproves raw odd parity;
   L309 further proves that cubic/quartic preparation and complete
   quintic lower neutralization do **not** restore global odd-response
   parity.  Formulate instead, for
   `Bcal_r=[B_1 ... B_r]`, the triangular endpoint recurrence
   `U_(2r+1)=M_S(C_(2r+1))
      +X_r Bcal_r*+Bcal_r X_r*`
   with bounded analytic/Smith-divisible columns, and the even
   recurrence as the new direct Gram plus a response and a bounded
   factor through `Bcal_(r−1)`.  L310--L311 already make every
   physical root a literal hereditary factor plus an ideal-valued,
   delay-vanishing response.  Derive the remaining weighted
   square-completion from L306--L307's full matrix expression, using
   L305 only as the response/module engine and keeping the affine
   leftover in L307's scaling identity.  L312 now gives a proved
   finite-word prefix factor with a `ceil(d/2)` bound, but neither
   that result nor A259's sharper numerical conjecture may be applied
   termwise after the infinite Stein/Hardy closure: raw word length
   does not prove the required elliptic weight matching.  L313 proves
   the correct weighted factorization for every finite physical
   polynomial in `S,J`, including the complete direct-map jet.  L314
   closes the frame sector, and L315 closes the operator root
   valuation without a separate dual-Stein telescope: the bidegree
   rule `d>=2k` is invariant under all six successor terms, lower
   elimination, Stein inversion, and L311 response selection.  Do
   not redo those sectors or demand a frame/operator cancellation.
   L316 adds exact grouped channel parity, and L317 now closes the
   sharp support that parity alone could not supply.  The marked
   L313 recursion, L243/L245 port valuation, and L255/L266 closed
   leakage return prove that a pair `B_j,B_k` at elliptic degree `d`
   satisfies `j+k<=d`; equality is only L283/L298's principal direct
   Gram.  Hence the complete finite jet is already packaged as
   `Bcal(c){D_dir+cR(c)}Bcal(c)*+response`, with `R` locally bounded,
   Hermitian, polynomial, and rank-stable.

   Work directly with that grouped form.  At fixed `theta=1/2`,
   isolate each new direct Gram and complete squares against the
   `cR(c)` blocks in increasing flag order.  Charge only the prior-
   transfer Schur cost to an earlier retained margin, and use the
   extra `c` plus a sufficiently small common neighbourhood to make
   the bound strict.  Carry L307's affine leftover at every step and
   retain `−C*C` together with its positive partial-scale reserve;
   do not re-expand the square into separately estimated terms.
   Do
   not impose global odd response, confuse channel parity with
   elliptic parity, factor raw `q`, enumerate L305
   roots, launch an unrelated fifth/seventh grind, discard `−C*C`,
   split load-bearing square completions, assume A259's unproved
   endpoint-word inequality, or assume
   `I intersect [A,A]=[A,I]`.
2. **COMPLETED BY L318: fixed partial-scale finite-flag induction.**
   The remainder of this item is the derivation record, not a live
   instruction.  Combine
   the recurrence from item 1 with L283's complete-delay boundary,
   L290's hereditary endpoint factors, and L292's valuation test.
   Preserve at every active grade a fixed fraction of the lower direct
   Gram and the earlier quadratic upper gain; charge each even
   prior-bridge cost to that earlier margin.  L289 requires only a
   finite positive jet through terminal grade `L`, not infinite
   summability.  Do not demand exact lower-tightness, endpoint-null
   even gauges, raw two-sided state-ideal invariance, or use L291's
   pointwise pseudoinverse.  L305--L307 have closed response/module
   stability, affine recursion, and the partial-scale Gram identity;
   L308--L309 prove both raw and global prepared odd parity are false.
   If the hereditary recurrence is not margin-compatible, return to
   L296's polarized raw
   branch only after separating favorable Grams and bounded prior-flag
   factors; do not mix endpoint ledgers.
3. **COMPLETED BY L329: close the scalar-channel-stratum induction.**
   Work on the scalar conjecture, not the stronger complete-similarity
   normal flag.  L322 gives the sharp equality-anchor reserve, but do
   not add it to an independently chosen disk loss.  L325 replaces
   that shortcut by the actual-disk joint defect
   `Delta=min_k{<D_Xk,k>+1−sigma_X(k)}` and proves a uniform linear
   scalar gap.  L323 proves the transfer-only response is `O(delta)`.
   L324 includes every accumulated L197 quotient and proves the joint
   response is `O(delta+d)`, so its completed cost is quadratic.

   L326 closes L325's terminal zero face by Gau--Wu's classical
   reducing-model theorem, and L328 proves the summand has exactly
   one-copy dimension `L+1`.  It lies in L192's single-block tube and
   its complement has multiplicity `m−1`.  Do not rederive that model
   theorem, and do not replace its general finite Blaschke function
   by a monomial except on an explicitly nilpotent slice.

   L327 now proves the required quantitative comparison:
   `delta+d asymp Delta_k` on the same norming direction, uniformly
   along each ramified finite rank stratum.  Retain `Delta_k` until
   every direction-attached response has been paid; do not replace it
   prematurely by `min_k Delta_k`.

   Assemble the final finite curve-selected estimate.  At each L197
   layer,
   retain a fixed fraction of every accumulated square, stack all
   normal/reflected responses before completing them once, and charge
   the resulting `O((delta+d)²)` cost to L325's joint linear reserve.
   Descend L321's common eigenline when the equality-anchor leakage
   face vanishes.  Keep L318's elliptic margin separate and use
   L326/L328 to remove every exact one-copy equality summand.

   If `Lambda=0`, carry L321's common eigenline through the next
   nonzero transfer jet.  If it persists identically, L326/L328
   supply the reducing one-copy model and multiplicity descent;
   L205/L192 remain
   the explicit full-equality chart.  At a tied channel collision,
   keep L199's first
   disk/circular-normal face and L318's elliptic margin as separate
   diagonals.  Do not spend either twice.

   L319 proves that first-face positivity and exact Schur transport
   alone do not determine the stronger complete-similarity
   kernel-block response.  Do not use that shortcut.  Return to the
   endpoint-specific complete-matrix response only if the scalar
   induction leaves a term not controlled by
   L192/L197/L199/L318/L324/L325/L326/L327/L328.
4. **Then close analytic remainders and audit equality strata.**  Use
   the bounded selected columns along the finite L197/L220 flag and
   combine with the banked single-block quotient results before
   claiming a neighbourhood theorem.
5. **Retain the CP/H-r routes as fallbacks, not concurrent
   frontiers.**  Test the full block-Toeplitz CP correction if the
   local merger stalls.  Return to shifted Möbius/odd-phase H-r only
   after the repeated merger is resolved.

## Superseded pre-L279 next actions (historical; do not resume)
1. **A194/A213/A216--A225: prove the sole remaining active support envelope.**
   L267's abstract lossless identity remains valid, but L268 proves
   that a single analytic port carrying only the deflated tail defect
   cannot reproduce the physical associated face: equal traces hide
   a rank increase from one to two.  Do not search for another
   one-shot operator congruence.  A larger multiport network is useful
   only if its extra defect channel is explicitly evaluated, which is
   equivalent to the remaining trace calculation.  L272 has now
   proved that every copy-dependent active term has the single-
   sandwich form `tr{L_B Phi(L,L*) L_B}` with copy-scalar `Phi`; do
   not repeat the unpaired-`B#` or ideal-placement calculation.  It
   does not control the unilateral boundary depth of `Phi` or prove
   lower vanishings.  L273 has separately proved the complete
   metric-inverse convolution: it contributes `−1` on every one of
   the `2k−2` interior fan words and no other nonradial term.  Do not
   recompute that metric half.  L274 fixes the four remote
   coefficients, and L275 now proves coefficient one on every
   endpoint and interior word of the unweighted fan in arbitrary
   grade.  Do not recompute any fan coefficient.  L276 shows that
   literal no-extra support is unnecessary.  L277--L278 now prove
   every lower face is zero operator-wise by factoring the odd and
   intact-metric even frontiers through the next delay ideal; do not
   recompute lower grades.  Prove only that the active associated
   difference is cyclically radial through index `k+2`.
   L269 disproves the
   stronger claim that the whole
   operator return is shallow: a depth-`k+1` term survives as
   `G_k−SG_kS*`.  L271/L270 remain a valid stronger operator route,
   but L276 supersedes them as the minimal scalar stop condition:
   long axes and the exact active monomial force all three
   trace-relevant coefficients of any bounded cyclic-radial
   difference.  Do not compute individual radial coefficients or
   classify extra words beyond cyclic support.  Prove the remaining
   active L276 envelope directly from the metric-deletion response
   around L278's intact-metric face.
   L256 proves the universal relative response
   `[c²]mu°=4||B_1||_F²` for every matrix channel.  Prove that removing
   one clean Hardy layer obeys
   `[c^(2k)]mu°_(B,k)=[c^(2k−2)]mu°_(Btilde,k−1)`
   in L257's canonical split, then iterate.  L258 has already eliminated
   the first monomial/final row and resummed every visit to it as
   `Z_ret=PZP+PZF(I_F−FZF)^−1FZP`.  Evaluate the active trace of
   `D_ret=I−R_PZ_ret`; L277--L278 already close every lower
   coefficient.  Evaluate only its active face by Hardy index:
   use L251's common analytic entry/exit channel, L272's normal
   ordering, and matrix-inner autocorrelation to organize future
   coefficients.  L259
   identifies those terms exactly: after complete delay they are the
   off-diagonal blocks `B_kB_j*`, `j>k`, of the same leakage row whose
   diagonal is `B_kB_k*`.  L261 shows that if the active physical
   return supplies the ordinary Hardy adjoint row, all future blocks
   close by Parseval as
   `sum B_kB_j*B_jB_k*=B_kB_k*`.  Prove precisely that unweighted
   return-metric statement from L243--L245 and L244, with the doubled
   remote amplitude on both sides.  L262 has now proved the metric
   part: every genuine normalized half-line Wold sandwich is unitarily
   the ordinary Hardy pairing.  Do not re-prove or estimate that
   metric.  Instead show that the first nonconstant term in L258's
   closed return is exactly such a Wold sandwich after the full
   L260-mandated orientation sum, including the zero/one-kernel
   filtration and boundary/direct-map terms.  L263 now also proves
   that the fixed boundary port and analytic defect-graph/Wold motion
   cannot alter that first face: after the background graph is
   removed, compute only the constant Hardy compression of the
   physical numerator.  Identify its trace as four times the full
   leakage-row norm plus trace-null terms; do not assert a positive
   operator Gram.  This must be a trace identity—
   direct audits disprove the stronger operator covariance.  L260
   also rules out separating the two ellipse orientations.  Do not
   claim equality of the full shifted series, which L235 disproves.
   Use L255's doubled-Hardy block pencil behind that renewal.  Its
   `Xi_inf direct-sum Xi_inf*` diagonal is the L244 background.  Keep
   the terminal multiplier at its physical value `lambda=1`, pair
   opposite analytic ports as in L251, and apply L244's coisometry
   before extracting the active coefficient.  L272 has already used
   L243's reflection-power filtration to exclude multiple
   nonbackground leakage selections.
   The desired surviving quadratic trace is the shifted **left**
   complement from L254:
   `T_BT_B*=I−O_Lhat O_Lhat*`.  Under complete delay its finite Hardy
   window contains only `B_k`; equivalently L252 reduces it to
   `4tr(Q_(k+2)−(k+2)Q_1+(k+1)I)`.
   The exact grade-one calculation shows again that the raw Stein and
   whitening pieces are not meaningful separately.  Do not estimate
   the huge dual corner and Schur square, detach the doubled port
   (A195/L249), assume A192's symmetric endpoint valuation, seek a
   positive leakage-Gram factorization, or extend the grade-six audit
   instead of proving the covariance.
   Use L264's `Rtri`: the final port/cross row and active left orbit
   remain removed, while the universal right-half-line `D_R` metric
   is restored at no cost to the active scalar trace.
   In that gauge, L266 proves that any complete physical return of the
   form `P_kL_B Psi_phys(L)L_BP_k` is rigidly reduced to
   `[Psi_phys]_0||B_k||²`.  Its low-depth extension covers
   copy-scalar shift words only when every unilateral boundary
   correction lies below row `k`; depth `k+1` can fail.  Therefore
   L269 proves that deeper boundary terms do occur, paired at the
   first grade as `G_k−SG_kS*`.  Therefore prove that L243's
   zero/one-kernel sectors, L251's paired ports, and L258's renewal
   assemble in the exact L269 normal form: a radial
   `P_k(Q_0,...,Q_k)`, the fixed pair `R_2−R_1`, and
   `G_k−SG_kS*`, plus only further trace-null divergences;
   exclude copy-dependent insertions but do not falsely delete the
   deep terms.  L270 then forces the radial moments and proves that
   the full-minus-tail trace is zero, so iteration reaches L256's
   proved value four.  There is no separate symbol-constant
   calculation on this route.  Prove the lower vanishings in the same
   structural assembly.
2. **Then promote the delayed trace law, not more finite evidence.**
   Combining the `+4` volume coefficient with L247's exact metric
   contribution `−2` would give L225's total trace
   `+2||B_k||_F²`, hence the fully delayed separator trace required by
   L222.  Record separately that a bounded pathwise selection through
   rank changes would still be needed.
3. **Return to A178's all-series selection/recurrence debt.**  The
   cubic-through-sextic columns remain four separate finite
   certificates, not an induction.  Use A171/L228 and L237's exact
   arbitrary-delay reflection recurrence to seek a uniform partial
   right-ideal construction and a coefficient bound converging in
   L194's analytic chart.  Do not compute another isolated grade unless
   it exposes that recurrence.
4. **Finish the repeated normal/elliptic merger after the delayed
   trace and selection gates.** L197 closes the
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
metric; do not assume definite parity; and do not vary L245's doubled
terminal edge independently of the physical theta/metric balance
(A195).  Do not split the two formal ellipse orientations before the
theta/coisometric cancellation (A206).
Keep committing+pushing after each task (user instruction).

## Files map (handoff-ready, refreshed 2026-07-26)
Current A178 transport packet:
  proof/repeated_crabb_gauge_transport_homology.md (L285 recurrence),
  proof/repeated_crabb_first_gauge_transport.md (L286 odd successor),
  proof/repeated_crabb_sextic_transport_trace_obstruction.md
  (L287 global-null obstruction), and
  proof/repeated_crabb_second_gauge_transport.md (L288 correct
  flag-preserving even successor), and
  proof/repeated_crabb_two_graph_transport.md (L289 exact mixed-graph
  endpoint transport and finite-jet stop condition), and
  proof/repeated_crabb_graph_flag_ideal.md (L290 exact exposed-channel
  to endpoint-flag functoriality), and
  proof/repeated_crabb_endpoint_gauge_normal_form.md (L291 graph-gauge
  quotient and pointwise endpoint-factor equivalence), and
  proof/repeated_crabb_endpoint_factor_valuation.md (L292 exact
  arcwise Smith-valuation criterion and pointwise-only
  obstructions), and
  proof/repeated_crabb_schur_markov_laplacian.md (L293 exact
  Schur-jet commutator form of the alternative Markov gap), and
  proof/repeated_crabb_canonical_flux_inequality.md (L294 exact
  canonical flux/Dirichlet pairing and flagged L282 bound), and
  proof/repeated_crabb_canonical_target_scope_obstruction.md (L295
  exact obstruction to comparing the raw repair upper face directly
  with L212's lower-tight target), and
  proof/repeated_crabb_polarized_flux_inequality.md (L296 exact
  single/summed/flagged polarized response-energy bound), and
  proof/repeated_crabb_lower_retightening_normalization.md (L297 exact
  raw-to-lower-tight endpoint ledger), and
  proof/repeated_crabb_oriented_retightening_transport.md (L298 exact
  all-grade Stein-compatible oriented retightening), and
  proof/repeated_crabb_moving_retightening_defect.md (L299 exact
  six-term nonlinear defect and rational cubic flag obstruction), and
  proof/repeated_crabb_cubic_markov_flux.md (L300 exact
  frame-eliminated cubic pairing and gap-free response-column bound),
  and proof/repeated_crabb_successor_markov_transport.md (L301 exact
  universal first-successor response and selection independence),
  and proof/repeated_crabb_retightening_quartic_invariant.md (L302
  exact two-ended quartic class and prior-flag lower neutralization),
  and proof/repeated_crabb_retightening_full_quartic.md (L303 exact
  complete partial-retightening Schur-graph assembly and pointwise
  margin), and
  proof/repeated_crabb_retightening_quartic_flux.md (L304 exact
  relative-commutator flux and rank-stable analytic quartic
  selection), and
  proof/repeated_crabb_rooted_bridge_flux.md (L305 exact
  arbitrary-word rooted quotient/response split and
  ideal-preserving correction), and
  proof/repeated_crabb_copy_energy_quotient.md (L306 exact
  all-series initial-copy energy reduction and gap-free residual
  response), and
  proof/repeated_crabb_affine_copy_recurrence.md (L307 exact
  arbitrary-step affine cancellation and partial-scale Gram
  identity), and
  proof/repeated_crabb_raw_quintic_obstruction.md (L308 exact
  obstruction to one-shot raw odd parity), and
  proof/repeated_crabb_prepared_quintic_obstruction.md (L309 exact
  obstruction to global prepared odd-response parity);
  matching
  regenerators use the same
  basenames under experiments/.  Resume on the partial-retightening
  branch at fixed `theta=1/2` by constructing L306--L309's
  **hereditary prepared** affine recurrence through the finite
  transfer flag.  Response selection,
  bridge-module stability, rooted summation, affine cancellation, the
  one-quarter copy-Gram scaling term, the quartic bound, and the Smith
  gate are closed; raw/global prepared odd parity are false and
  affine positivity is not proved.  Factor every lower-neutral odd
  retained class through the cumulative active transfer row and
  isolate even direct Grams from prior-transfer costs.  Do not start
  another isolated-grade march, split the
  load-bearing terms, or assume an abstract relative-cyclic excision
  theorem.  If the prepared recurrence is not margin-compatible, return
  to the raw
  branch by
  retaining its lower budget, splitting L283's favorable `−12` Gram
  from L285/L289's mixed remainder, discarding bounded prior-flag
  factors, deriving L296 polarization columns for the surviving
  compression, and proving their total energy is
  `O(||U*B_k||²)` (or equivalently proving L292's valuations).
  L303--L304 have already assembled and bounded L289's complete
  moving-graph quartic endpoint; do not redo that bookkeeping or force
  the boundary choice `theta=1`.
  Do not
  mix those endpoint ledgers, force state witnesses into an
  ideal, use a rank-changing pseudoinverse, confuse pointwise flag
  zero with bounded divisibility, or compute grade seven.
proof/ — read in this order for the canonical gate and ordered successors:
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
