# Candidate elliptic Newton face along the disk equality family (2026-07-23)

## 1. Status

This note records a strong numerical law and the resulting exact proof
target.  It is **not yet a lemma**.

Let `p=L+1`, let `u=(u_1,...,u_(L-1))` satisfy

\[
 u_k=\overline {u_{L-k}},                             \tag{1}
\]

and set

\[
 H(a)=\frac12I+aZ(u).
\]

Let `X(a)` be L122's exact disk matrix constructed from `H(a)`.
For real `0<c<1`,

\[
 W(X(a)+cX(a)^*)=
 \{w+c\overline w:|w|\leq1\},                        \tag{2}
\]

so its normalized Riemann pullback is

\[
 T(a,c)=\phi_c(X(a)+cX(a)^*).                        \tag{3}
\]

Write `Gamma(a,c)` for L118's locally optimized rank-one Stein
condition square minus four.  L117 and L123 give the two exact axes

\[
\begin{aligned}
\Gamma(0,c)&=-16c^{2L}+O(c^{4L}),\\
\Gamma(a,0)&=0.
\end{aligned}                                         \tag{4}
\]

The candidate joint lower face is

\[
\boxed{
\Gamma(a,c)=
-16c^{2L}
-64a^2\sum_{k=1}^{L-1}|u_k|^2c^{2k}
+\text{terms strictly above this Newton face}.}       \tag{5}
\]

The sign and the coefficient `64` are consistent in every tested
circle grade.  The live task is to prove (5), with a uniform analytic
remainder statement strong enough for the L124 merger.

## 2. Why these are the right exact curves

The matrix `X(a)` is not a linearized disk tangent.  It stays on the
exact L123 equality family for every sufficiently small `a`.  Equation
(2) is also exact because

\[
\langle(X+cX^*)x,x\rangle
=w+c\overline w,\qquad w=\langle Xx,x\rangle.
\]

Thus (3) measures only the elliptic normal to an exact equality
anchor.  There is no contamination from an approximate disk chart or
from a numerically fitted domain.

The coefficient vector in (5) has a natural Hardy-space form.  L123's
characteristic factor is

\[
 g_a(\xi)=\xi^L+2a\sum_{k=1}^{L-1}u_k\xi^k,
\]

and hence

\[
64a^2\sum_k|u_k|^2c^{2k}
=16\|g_a-\xi^L\|_{H^2(|\xi|=c)}^2.                   \tag{6}
\]

This suggests that the new face measures the failure of the
characteristic Blaschke extremal to remain monomial under elliptic
descent.

## 3. Deterministic numerical evidence

`experiments/crabb_palindromic_elliptic_face.py` uses only exact
geometric anchors:

1. construct `X(a)` from the positive Toeplitz chart;
2. form the exact ellipse (2);
3. evaluate the classical elliptic Riemann map by spectral functional
   calculus;
4. optimize L118's rank-one Stein branch; and
5. independently solve the full L21 similarity SDP.

For each size `p=3,...,7`, it tests every first offset
`1<=k<=floor(L/2)` using the phase-one vector supported at `k` and
`L-k`.  At amplitudes `.03,.06`, the change from the axis is quadratic:
the smaller/larger deficit ratio differs from `(0.03/0.06)^2=1/4`
by at most `0.0011`.  The rank-one envelope and unrestricted SDP agree
within `3.5e-9`.

At `c=.12`, the measured coefficient divided by the right side of
(5) is:

\[
\begin{array}{c|c}
\text{first offset}&\text{observed/predicted range}\\ \hline
1&0.935\text{--}0.968\\
2&0.9991\text{--}0.9995\\
3&0.999987 .
\end{array}
\]

The offset-one ratio approaches one as `c` decreases.  Its visible
finite-`c` correction is expected because (5) asserts the lowest
Newton face, not an exact formula for the whole Hessian at fixed `c`.
Offsets two and three have much later corrections.

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_palindromic_elliptic_face.py \
  --output experiments/crabb_palindromic_elliptic_face_s70223.jsonl
```

## 4. Proposed proof route

Work in L123's coefficient coordinates, avoiding matrix square roots:

\[
\begin{aligned}
K(a)&=H(a)+R^*H(a)R,\\
S(a,c)&=2K(a)^{-1}\{H(a)R+cR^*H(a)\},\\
T(a,c)&=\phi_c(S(a,c)).
\end{aligned}                                         \tag{7}
\]

Ordinary eigenvalues and Stein condition numbers in the physical
coordinates are respectively the eigenvalues and generalized
condition numbers in the `K(a)` inner product.

The required calculation should be organized as follows.

1. Expand the inverse ellipse map, then solve
   `Psi_c(T)=S(a,c)` coefficient by coefficient.  This avoids
   square roots and elliptic eigenvectors.
2. Expand the rank-one Stein recurrence together with its unique
   defect-vector minimizer from L118.
3. Use the companion form from L123.  A coefficient `u_k` cannot
   reach the endpoint defect before `k` backwards elliptic steps.
4. Show that all shorter paths cancel because `Gamma(a,0)=0`
   identically on the equality branch.
5. The two unique endpoint paths at length `k` should contribute
   `-64|u_k|^2c^(2k)`.  Distinct grades cannot mix on the first face.

This is the elliptic analogue of L120--L121's endpoint path selection,
but now at second order in the equality amplitude.

## 5. How (5) would enter the final merger

L124 writes every disk-flat coefficient as `z=u+v`, with
`||v||<=||u||`, and gives the disk-normal face

\[
-128\|u\|^2\|v\|^2.                                  \tag{8}
\]

Together, the three principal negative quantities would be

\[
c^{2L},\qquad
\sum_k|u_k|^2c^{2k},\qquad
\|u\|^2\|v\|^2.                                      \tag{9}
\]

The remaining analytic task would then be a Newton-polyhedron
absorption theorem: prove that every mixed remainder monomial lies
strictly above the lower hull generated by (9), after the already
completed strong-variable maximization.  One must not invoke that
absorption before proving the exact face and auditing sharp boundary
regimes.
