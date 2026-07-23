# The general-disk first reflection (A101, 2026-07-23)

## 1. Result and correction

The one-reflection stationarity proved by L148/L149 on the
phase-palindromic equality cone does **not** extend identically over
the whole Toeplitz disk chart.

There is an exact finite rational formula for the first reflected
derivative.  It gives a nonzero value already in size four.  However,
exact normal blow-ups indicate the structurally sufficient replacement:

\[
\boxed{
|\partial_c R(z,0)|\le C_L{\cal Q}(z),
\qquad
{\cal Q}(z)=\|z\|^4-|z^TJz|^2.
}                                                     \tag{1}
\]

Here `R` is the norm square of the corrected characteristic
Faber--Blaschke product.  Inequality (1) remains a **candidate**, not a
proved all-size statement.  If proved, every ordinary grade-one term
is absorbed by L152's disk-normal deficit after shrinking `|c|`.
The analogous marked estimate is needed for higher reflected grades.

## 2. Exact first elliptic tangent

Let `A=A(z)` and `K=K(z)` be L122's coefficient disk model and write

\[
\det(wI-A)=w g(w).
\]

The centered ellipse pullback has first tangent

\[
\boxed{
T(c)=A+cE+O(c^2),\qquad E=K^{-1}A^*K-A^3.
}                                                     \tag{2}
\]

Write

\[
g(w)=\sum_{j=0}^{L}a_jw^j,\qquad a_L=1.
\]

L148's composed natural Faber tangent is

\[
w^3g'(w)-{g'(w)-g'(0)\over w}.                      \tag{3}
\]

For the sharp raw factor, L151 adds the grade-one Hardy correction
`a_(L-1)w^(L-1)`.  Monic preparation therefore gives

\[
\boxed{
n_1=
\operatorname{rem}_g\left[
w^3g'(w)-{g'(w)-g'(0)\over w}
+a_{L-1}w^{L-1}
\right].
}                                                     \tag{4}
\]

Thus

\[
N_c=g+cn_1+O(c^2),\qquad D_c=N_c^\sharp.             \tag{5}
\]

Equations (2), (4), and rational Fréchet differentiation give

\[
\begin{aligned}
\dot B={}&
\{Dg(A)[E]+n_1(A)\}D(A)^{-1}\\
&-B(A)\{Dg^\sharp(A)[E]+n_1^\sharp(A)\}D(A)^{-1}.
                                                               \tag{6}
\end{aligned}
\]

Since `B(A)` is rank one with a simple top right line `e_L`, the
derivative is exactly

\[
\boxed{
\partial_cR(z,0)
=
{2\operatorname{Re}
  \{(B(A)e_L)^*K\dot B e_L\}
 \over e_L^*Ke_L}.
}                                                     \tag{7}
\]

This is a finite rational function of the real and imaginary parts of
the Toeplitz coefficients.  No root tracking, finite difference, or
Riemann-map quadrature enters.

## 3. Exact failure of global stationarity

For the real size-four disk point

\[
z_1=1/20,\qquad z_2=1/30,
\]

equations (2)--(7) give

\[
\boxed{
\partial_cR(z,0)
=-{45448825601074157677925589250285395456\over
5359637767128852460594708556167336924387705}
\ne0.
}                                                     \tag{8}
\]

Numerically, (8) is about `-8.48e-6`.  The small size explains why an
ordinary binary64 fit can easily mistake it for zero.

This is independent of A100's different obstruction.  A100 falsifies
equality of the primal model condition and dual norm at `c=0`; (8)
falsifies stationarity of the dual norm in the first elliptic
reflection at a general disk point.

## 4. What survives on the equality cone

If `u` is phase-palindromic, L148 proves

\[
\partial_cR(u,0)=0.                                  \tag{9}
\]

The exact formula (7) reproduces (9).  More strongly, take a fixed
phase-one equality point and an anti-palindromic normal:

\[
z=u+\epsilon v,\qquad Ju=u,\quad Jv=-v.
\]

For three rational anchors in lengths three, four, and five, exact
evaluation at
`\epsilon=1/200,1/400,1/800,1/1600` gives bounded ratios

\[
{\partial_cR(u+\epsilon v,0)\over{\cal Q}(u+\epsilon v)}
\]

converging respectively near `-4.36`, `-7.15`, and `-5.12`.
The derivative itself is not even in `epsilon`; its odd part begins at
order `epsilon^3`, while its leading normal part is quadratic.

This is the correct behavior for (1).  A merely first-order normal
zero would instead make the displayed ratio diverge like
`1/epsilon`.

## 5. Candidate proof route for the `Q` bound

The remaining all-size statement can be isolated cleanly.

1. Prove from (6), using L123's companion columns, that the first
   normal derivative of (7) vanishes at every nonzero
   phase-palindromic equality anchor.
2. Prove the apex jet
   `partial_c R(z,0)=O(||z||^4)`.  L151 already removes amplitude
   degrees below two on every raw grade; the missing assertion is the
   one-reflection cubic cancellation.
3. Apply L124's best-phase blow-up.  Exact value and first-normal
   vanishing give `O(||u||^2||v||^2)` near a fixed branch, while the
   fourth-order apex jet makes the constant uniform as `u->0`.
   Since `Q=4||u||^2||v||^2`, this proves (1).

This mirrors L152's successful uniform disk-normal argument.  The new
algebraic content is only the mixed companion recurrence in items
one and two.

For a reflected grade `k`, the desired marked extension is

\[
|\text{one-}w_k\text{ sector}|
\le C_L |c|^k {\cal Q}(z).                           \tag{10}
\]

Terms in (10) are harmless against `-a_L Q`.  The compact
two-reflection sector remains the negative L151 face.  This replaces
the false demand that the entire one-reflection sector vanish.

## 6. Regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_disk_one_reflection.py \
  --output experiments/crabb_disk_one_reflection_s70223.jsonl
```

The first line records (8) exactly.  The remaining lines use exact
rational matrix algebra before converting only the final
derivative-to-`Q` ratios to binary64.  They are an adversarial probe
of (1), not its proof.
