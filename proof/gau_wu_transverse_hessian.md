# The non-Crabb Gau--Wu stratum has a strict transverse Hessian

> **Campaign scope.**  This is a fixed \(3\times3\), scalar, local
> theorem.  It covers every fixed nonzero member of L330's Gau--Wu
> family and its affine-unitary orbit.  It neither proves the global
> \(3\times3\) case nor the completely bounded conjecture.

## 1. Results (L332--L333, 2026-07-26)

For \(0<q<1\), put

\[
 a={2q\over1+q^2},\qquad
 r={\sqrt2(1-q^2)\over1+q^2},\qquad
 G_a=\begin{bmatrix}0&r&-2a\\0&a&r\\0&0&0\end{bmatrix}.
 \tag{1}
\]

L330 gives \(W(G_a)=\overline{\mathbb D}\) and

\[
 f_a(G_a)=2E_{13},\qquad
 f_a(z)={z(z-a)\over1-az}.                         \tag{2}
\]

The two conclusions are:

1. **L332 (strict optimized transverse Hessian).**  After quotienting
   the exact affine-unitary Gau--Wu equality stratum and optimizing
   both moving zeros of the degree-two Blaschke product, the second
   variation of the sharp scalar norm is negative definite.  The
   equality tangent has real dimension \(13\), so the normal slice
   has dimension \(5\).
2. **L333 (non-Crabb local theorem).**  Every fixed \(G_a\) with
   \(0<|a|<1\), and every affine-unitary image of it, has a full
   operator neighbourhood on which

   \[
   \boxed{\|h(A)\|\le2\max_{z\in W(A)}|h(z)|}       \tag{3}
   \]

   for every scalar function \(h\) analytic on \(W(A)\) and
   continuous on its closure.  Equality sufficiently near \(G_a\)
   occurs only on the classical Gau--Wu disk-model stratum.

The constants are local in \(a\).  They degenerate as \(a\to0\),
where the nonnilpotent chart collides with the Crabb chart handled by
L192/L329, and as \(|a|\to1\), where the support gap closes.

## 2. The equality tangent is 13-dimensional

For complex \(b\), the Gau--Wu model is

\[
 G_b=\begin{bmatrix}
 0&\sqrt{2(1-|b|^2)}&-2\overline b\\
 0&b&\sqrt{2(1-|b|^2)}\\
 0&0&0
 \end{bmatrix}.                                   \tag{4}
\]

If \(b=\rho e^{i\psi}\) and
\(D_\psi=\operatorname{diag}(e^{i\psi},1,e^{-i\psi})\), then

\[
 G_b=e^{i\psi}D_\psi^*G_\rho D_\psi.               \tag{5}
\]

Thus the argument of \(b\) is redundant with global rotation and
unitary conjugacy.  The full nearby equality stratum has:

- eight real unitary-orbit directions, since \(G_a\) is irreducible;
- four real affine directions \(I,iI,G_a,iG_a\); and
- one real radial model direction.

This predicts dimension \(13\).  An exact determinant proves it and
simultaneously supplies a convenient complement.  Let

\[
 E(e)=\begin{bmatrix}
 0&0&0\\0&0&0\\
 e_0+ie_1&e_2+ie_3&e_4
 \end{bmatrix}.                                   \tag{6}
\]

Append these five coordinate directions to thirteen explicit
equality tangents.  In interlaced real matrix coordinates the
resulting \(18\times18\) determinant is

\[
 \boxed{
 {-2048q^2(q-1)(q+1)(q^4+1)^2\over(q^2+1)^8}\ne0.} \tag{7}
\]

Hence (6) is a genuine normal slice for every \(0<q<1\).  The earlier
numerical dimension count of \(12\) mistakenly imported the
dimension of a differently normalized disk-matrix manifold; (7)
corrects it before any proof is stacked on that count.

## 3. Smooth conformal second variation

Write

\[
 H_0(\theta)=\operatorname{Re}(e^{-i\theta}G_a),
 \qquad
 H_E(\theta)=\operatorname{Re}(e^{-i\theta}E).
\]

The eigenvalues of \(H_0(\theta)\) are

\[
 1,\quad-1,\quad m(\theta)=a\cos\theta,             \tag{8}
\]

so the top eigenvalue has the uniform gap \(1-a\).  Let \(P_+\) be
its projection and put

\[
 R_\theta={P_-(\theta)\over2}
          +{P_m(\theta)\over1-m(\theta)}.           \tag{9}
\]

For the straight physical perturbation \(G_a+\varepsilon E\), the
first two support coefficients are

\[
 s_E=\operatorname{tr}(P_+H_E),\qquad
 t_E=\operatorname{tr}(P_+H_ER_\theta H_E).        \tag{10}
\]

If \(H_s\) is the Schwarz transform of \(s_E\), put

\[
 F(w)=wH_s(w),\qquad
 \delta=\operatorname{Im}H_s(e^{i\theta})-s_E'(\theta),
\qquad
 \kappa=t_E-\frac12\delta^2,                       \tag{11}
\]

and let \(K(w)=wH_\kappa(w)\).  The inverse Riemann normalization is

\[
 T_\varepsilon
 =G_a+\varepsilon{\cal G}
       +\varepsilon^2{\cal H}+O(\varepsilon^3),
 \tag{12}
\]

\[
 {\cal G}=E-F(G_a),\qquad
 {\cal H}=-DF(G_a)[E]+(F'F)(G_a)-K(G_a).           \tag{13}
\]

These are the smooth-gap versions of the support/Riemann formulas in
L101 and the second-order calculation recorded in
`proof/general_similarity_second_order.md`.

One correction is load-bearing here: \(G_a\) is not nilpotent.
Because its minimal polynomial is \(z^2(z-a)\), evaluating an
analytic function needs its value and first derivative at zero and
its value at \(a\).  Evaluating \(DF(G_a)[E]\) needs the degree-five
Hermite interpolant, equivalently the functional calculus of

\[
 \begin{bmatrix}G_a&E\\0&G_a\end{bmatrix}.         \tag{14}
\]

Truncating the power series after degree two, as is valid at a
Crabb block, gives a false Hessian here.  The checker uses (14), not
that nilpotent shortcut.

After (1), all circle data in (10)--(13) are rational in the phase.
The Schwarz means and the values needed in (13) are the residues at
\(0,q,a\).  Thus the calculation below is exact.

## 4. Optimizing the two Blaschke zeros

Allow the two zeros to move as

\[
 \alpha_\varepsilon
 =\varepsilon(u_0+iu_1)+O(\varepsilon^2),\qquad
 \beta_\varepsilon
 =a+\varepsilon(u_2+iu_3)+O(\varepsilon^2).        \tag{15}
\]

Second-order zero accelerations do not contribute because the
zero-gradient vanishes at (2).  Expand

\[
 b_{\alpha_\varepsilon}(T_\varepsilon)
 b_{\beta_\varepsilon}(T_\varepsilon)
 =2E_{13}+\varepsilon Y_1+\varepsilon^2Y_2
  +O(\varepsilon^3).                              \tag{16}
\]

Since the top singular value in (2) is simple, its norm coefficient
is

\[
 {1\over4}\left\{
 4\operatorname{Re}(Y_2)_{13}
 +\sum_i|(Y_1)_{i3}|^2
 +\sum_{j=1}^2|(Y_1)_{1j}|^2
 \right\}.                                        \tag{17}
\]

This is a homogeneous quadratic form in \((e,u)\).  Its pure
zero-motion block is two copies of

\[
 Z_q=\begin{bmatrix}
 -{q^4+4q^2+1\over(q^2+1)^2}&-1\\
 -1&-{(q^2+1)^2\over(1-q^2)^2}
 \end{bmatrix}.                                   \tag{18}
\]

Both \(-Z_{q,11}\) and

\[
 \det Z_q={6q^2\over(1-q^2)^2}                    \tag{19}
\]

are positive.  The two zero velocities therefore have a unique
analytic optimizing branch.  Taking the Schur complement of (18)
gives

\[
 \max_u\|b_{\alpha_\varepsilon}(T_\varepsilon)
              b_{\beta_\varepsilon}(T_\varepsilon)\|
 =2+\varepsilon^2 e^TQ(q)e+O(\varepsilon^3).       \tag{20}
\]

## 5. The exact five-normal Hessian

In the order used in (6),

\[
 Q(q)=
 \begin{bmatrix}
 A&0&B&0&C\\
 0&D&0&E&0\\
 B&0&F&0&G\\
 0&E&0&H&0\\
 C&0&G&0&I
 \end{bmatrix},                                   \tag{21}
\]

where, with \(d=q^6+2q^4-4q^2+3\),

\[
\begin{aligned}
A&={q^{10}+3q^8-11q^6+17q^4-29q^2+21\over16(q^2-1)},\\
B&=-{\sqrt2 qd\over16},\qquad
C=-{q^2d\over16(q^2-1)},\\
D&={q^{10}-3q^8+5q^6-q^4-13q^2-21\over16(q^2+1)},\\
E&=-{\sqrt2 q(q^4+1)(q^4-q^2+3)\over16(q^2+1)},\\
F&={q^2(q^{12}-3q^8+3q^4-64q^2-1)
       \over8(1-q^2)^2(1+q^2)^2},\\
G&={\sqrt2 q^3(q^8+2q^6-2q^2+31)
       \over16(q^2-1)(q^2+1)^2},\\
H&={q^2(q^4+1)^2\over8(q^2-1)(q^2+1)},\\
I&={q^2(q^8+3q^6+3q^4-15q^2+16)
       \over16(q^2-1)(q^2+1)^2}.
\end{aligned}                                     \tag{22}
\]

The real block uses indices \((0,2,4)\), and the imaginary block
uses \((1,3)\).  Put \(t=q^2\) and define

\[
\begin{aligned}
p_1(t)&=t^5+3t^4-11t^3+17t^2-29t+21,\\
p_2(t)&=t^9-5t^8+12t^7+12t^6+22t^5-150t^4\\
      &\hspace{21mm}+284t^3-484t^2+337t+3,\\
p_3(t)&=t^3-3t^2+7t+3,\\
p_4(t)&=-t^5+3t^4-5t^3+t^2+13t+21.
\end{aligned}                                     \tag{23}
\]

The leading Sylvester minors of the negative blocks are

\[
\begin{aligned}
\Delta^R_1&={p_1(t)\over16(1-t)},\\
\Delta^R_2&={q^2p_2(t)\over32(1-t)^3(1+t)^2},\\
\Delta^R_3&={q^4(1+q^4)^2p_3(t)
                    \over32(1-t)^3(1+t)^2},\\
\Delta^I_1&={p_4(t)\over16(1+t)},\\
\Delta^I_2&={q^2(1+q^4)^2p_3(t)
                    \over32(1-t)(1+t)^2}.
\end{aligned}                                     \tag{24}
\]

Each polynomial in (23) is positive on \(0\le t\le1\).  A compact
exact certificate is supplied by its Bernstein coefficients:

\[
\begin{array}{c|l}
p_1&(21,76/5,111/10,38/5,21/5,2)\\
p_2&(3,364/9,580/9,1646/21,760/9,1760/21,
      1616/21,580/9,48,32)\\
p_3&(3,16/3,20/3,8)\\
p_4&(21,118/5,263/10,143/5,153/5,32).
\end{array}                                       \tag{25}
\]

Equations (24)--(25) and Sylvester's criterion prove

\[
 \boxed{Q(q)\prec0\qquad(0<q<1).}                 \tag{26}
\]

This proves L332.  It also explains the initial numerics: there are
exactly five negative modes, not five negative modes plus a soft
kernel.  Two eigenvalues become small near \(q=0\), but they are
strictly negative at every fixed nonzero model.

## 6. From the Hessian to the local theorem

Fix \(q_0\in(0,1)\).  The support gap in (8) makes the normalized
Riemann pullback \(C^2\), in fact real analytic, in a neighbourhood
of \(G_a\).  More explicitly, the support boundary
\(e^{i\theta}(h_A(\theta)+ih_A'(\theta))\) is analytic in the matrix
entries and remains strictly convex; the normalized
boundary-correspondence equation has the periodic Hilbert transform
as its linearization, so the analytic implicit-function theorem in a
fixed Hölder class gives the asserted Riemann chart.  Equations
(18)--(19) and the finite-dimensional implicit function theorem
give one analytic local maximizing branch for the two Blaschke
zeros.  L331 gives zero ambient gradient all along the equality
stratum.  Since the optimized value is identically two on that
stratum, its Hessian kills every equality tangent and every
tangent-normal cross.  Equations (7) and (26) therefore give the
Morse--Bott estimate

\[
 {\cal C}(A)\le2-c_{q_0}\,
      \operatorname{dist}(A,{\cal M}_{GW})^2       \tag{27}
\]

in a sufficiently small tube, with equality only on the local
Gau--Wu stratum \({\cal M}_{GW}\).

It remains to exclude a distant Blaschke branch.  Crouzeix's finite
extremal theorem reduces every \(3\times3\) scalar extremal to a
Blaschke product of degree at most two.  At \(G_a\), equality in such
a product invokes Gau--Wu Theorem 8.  Irreducibility makes its
reducing model the whole three-dimensional matrix.  The zeros of
the model are its eigenvalues \(\{0,0,a\}\), so the Blaschke zeros
are exactly \(\{0,a\}\), up to order and phase.  Thus the maximizing
zero pair is isolated in the compactified degree-at-most-two
parameter space.  Every other branch has a strict gap at \(G_a\),
and continuity preserves that gap nearby.

Combining this separation with (27) proves (3), hence L333.

## 7. Literature and novelty boundary

- Gau--Wu (2009) supplies the exact equality models and the reducing
  equality theorem.  Those inputs are classical.
- Crouzeix's finite Blaschke reduction is classical.
- Lewis--Overton analyze the \(12\)-dimensional manifold of
  centered \(3\times3\) disk matrices and partial smoothness of the
  numerical radius; they do not state the optimized Crouzeix
  Hessian (21) or a constant-two neighbourhood theorem.
- Greenbaum--Lewis--Overton establish first-order variational
  stationarity at classical candidate minimizers.  Li studies
  uniqueness of extremal Blaschke products, including nonuniqueness
  at some \(3\times3\) elliptic matrices.  Neither source gives
  (21)--(27).

The exact Hessian and the resulting local scalar theorem therefore
appear new, pending publication-level audit.  No novelty is claimed
for the equality classification itself.

## 8. Regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/gau_wu_transverse_hessian.py \
  --output experiments/gau_wu_transverse_hessian_s70224.jsonl
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/gau_wu_transverse_hessian_numeric.py \
  --output \
  experiments/gau_wu_transverse_hessian_numeric_s70224.jsonl
```

The checker regenerates:

1. the \(13+5\) determinant (7);
2. the support projections, reduced resolvent, Schwarz transforms,
   second normal data, and nonnilpotent Hermite calculus;
3. the joint nine-variable quadratic form in five physical
   directions and four zero velocities;
4. the zero-motion block (18), its exact Schur complement
   (21)--(22), and every Sylvester factor in (24); and
5. the positive Bernstein lists (25), followed by exact rational
   specializations at five values of \(q\).

An independent Theodorsen/Blaschke finite-difference probe at
\(a=1/2\) gave the same five normal eigenvalues to the expected
finite-difference accuracy: the maximum Hessian-entry error was
\(2.93\cdot10^{-5}\), and the maximum eigenvalue error was
\(4.01\cdot10^{-5}\).  It was used only for falsification;
the proof is the exact residue/Hermite certificate above.  The five
tracked exact records have SHA-256
`8d75b8541f19aba118bc73840b2f0c428058b693f582b7c8ab42a017b2999752`;
the independent numerical record has SHA-256
`7ea05baeef3557d64761b1e9d731f21a2dddcc670d35ee747f46c0a36c2d85b7`.
