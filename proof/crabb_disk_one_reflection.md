# The general-disk first reflection (A101, 2026-07-23)

## 1. Result and correction

The one-reflection stationarity proved by L148/L149 on the
phase-palindromic equality cone does **not** extend identically over
the whole Toeplitz disk chart.

There is an exact finite rational formula for the first reflected
derivative.  It gives a nonzero value already in size four.  However,
an all-size endpoint recurrence proves the structurally sufficient
replacement:

\[
\boxed{
|\partial_c R(z,0)|\le C_L{\cal Q}(z),
\qquad
{\cal Q}(z)=\|z\|^4-|z^TJz|^2.
}                                                     \tag{1}
\]

Here `R` is the norm square of the corrected characteristic
Faber--Blaschke product.  Inequality (1) is proved in Section 5.  Every
ordinary grade-one term is therefore absorbed by L152's disk-normal
deficit after shrinking `|c|`.  The same recurrence gives the marked
estimate needed for higher reflected grades.

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

## 5. All-size endpoint first-jet theorem

Put `q_z=H(z)e_0` and define the endpoint-resolvent residual

\[
{\cal R}_z(\xi)
=q_z^*(\xi I-A_z)^{-1}e_L
-{g_z^\sharp(\xi)\over \xi g_z(\xi)}.                \tag{10}
\]

The missing statement is now proved:

\[
\boxed{
{\cal R}_{u+\epsilon h}(\xi)=O(\epsilon^2)
}                                                     \tag{11}
\]

at every phase-palindromic equality point `u`, for every complex disk
direction `h`.  The perturbation parameter is real, so `h` and
`conjugate(h)` vary the two Hermitian Toeplitz halves together.  These
directions span the full real disk chart.

It is enough to work in the phase-one gauge
`u_k=conjugate(u_(L-k))`; L123's diagonal gauge transports the result
back to every phase.  Write

\[
\begin{aligned}
U(\xi)&=\sum_{k=1}^{L-1}u_k\xi^k,&
g&=\xi^L+2U,&g^\sharp&=1+2U,\\
P_h(\xi)&=\sum_{k=1}^{L-1}h_k\xi^k,&
P_h^\flat(\xi)&=\sum_{k=1}^{L-1}\overline{h_{L-k}}\xi^k .
\end{aligned}                                        \tag{12}
\]

### 5.1 The tangent has only three boundary pieces

Let `T=dot H`, `dot K=T+R^*TR`, and let `C=A_u` be L123's companion.
For a row vector put

\[
\alpha(a)=\sum_{k=1}^{L-1}a_ke_{k+1}^*,
\qquad
v_h=\sum_{r=1}^{L-1}
 (h_{L-r}-\overline{h_r})e_r.
\]

Differentiating `KA=2HR` and using the companion columns gives the
inverse-free identity

\[
\boxed{
D:=K\dot A
=2TR-\dot K C
=e_0\alpha(h)+2v_h\alpha(u)
 -e_L\alpha(h^\flat).
}                                                     \tag{13}
\]

This is entrywise.  The top row is the direct tangent, the bottom row
is minus its reverse conjugate, and each interior row is the equality
row multiplied by `2(h_(L-r)-conjugate(h_r))`.  No induction or
matrix inverse remains in (13).

### 5.2 Two endpoint recurrences

Set `P=xi I-C` and

\[
\begin{aligned}
c(\xi)&=(g^\sharp+1,\xi,\xi^2,\ldots,\xi^L)^T,\\
d(\xi)&=(\xi^L,\xi^{L-1},\ldots,\xi,g^\sharp+1).
\end{aligned}
\]

Direct multiplication of the companion columns and the Toeplitz rows
gives

\[
\boxed{
Pc=\xi g e_L,\qquad dKP=\xi g e_0^*K.
}                                                     \tag{14}
\]

Consequently, as polynomial identities,

\[
\operatorname{adj}(P)e_L=c,\qquad
e_0^*K\operatorname{adj}(P)K^{-1}=d.                 \tag{15}
\]

Substitution of (13) into the two endpoint vectors telescopes to one
scalar:

\[
\boxed{
dDc=\xi\{P_hg-2g^\sharp P_h^\flat\}.
}                                                     \tag{16}
\]

Indeed `alpha(h)c=xi P_h`, `alpha(u)c=xi U`,
`alpha(h^flat)c=xi P_h^flat`, and

\[
\sum_{r=1}^{L-1}\xi^{L-r}
(h_{L-r}-\overline{h_r})=P_h-P_h^\flat.
\]

### 5.3 Characteristic and cofactor derivatives

First suppose that `g` has simple roots, which is a dense condition
on the positive equality stratum.  At a root `lambda`, the rank-one
matrix

\[
Y=\operatorname{adj}(\lambda I-C)K^{-1}
\]

has right column `YK e_L=c`, left row `e_0^*KY=d`, and pairing
`e_0^*Kc=g^\sharp`.  The last equality is also immediate from
`e_0^*K=(1/2,u_1,\ldots,u_{L-1},0)`.  L123 puts the roots of `g`
strictly inside the disk, so `g^\sharp` is nonzero there.  Hence

\[
Y={cd\over g^\sharp}.
\]

Jacobi differentiation and (16) now give

\[
\left.D_h\det(\xi I-A)\right|_{\xi=\lambda}
=-\operatorname{tr}(YD)
=2\lambda P_h^\flat(\lambda).
\]

Both sides have degree at most `L`.  Their difference vanishes on all
roots of `g`, so it is a scalar multiple of `g`.  But
`Ae_0=0` and `Ae_1=2e_0` throughout the disk chart; every
characteristic polynomial is divisible by `xi^2`.  The difference is
therefore divisible by `xi^2`, whereas a generic `g` has only one
factor `xi`.  The scalar multiple is zero.  Continuity removes the
simple-root and genericity assumptions.  Thus

\[
\boxed{
D_h\det(\xi I-A)=2\xi P_h^\flat,\qquad
D_hg=2P_h^\flat,\qquad D_hg^\sharp=2P_h.
}                                                     \tag{17}
\]

Now let `F=q^*P^{-1}e_L`.  Equations (14)--(16), together with
`dot q^*=e_0^*T`, give

\[
\begin{aligned}
\dot F
&=\dot q^*P^{-1}e_L
 +e_0^*KP^{-1}K^{-1}DP^{-1}e_L\\
&={P_h\over \xi g}
 +{P_hg-2g^\sharp P_h^\flat\over \xi g^2}\\
&={2(P_hg-g^\sharp P_h^\flat)\over \xi g^2}.          \tag{18}
\end{aligned}
\]

By (17), the last expression is exactly
`D_h(g^sharp/(xi g))`.  This proves (11).  Equivalently, the endpoint
cofactor has the especially simple derivative

\[
\boxed{
D_h\{q^*\operatorname{adj}(\xi I-A)e_L\}=2P_h(\xi).
}                                                     \tag{19}
\]

### 5.4 The `Q` bound and marked grades

L149's one-reflection Hardy calculation is obtained by applying an
analytic finite coefficient functional to (10).  Its variation causes
no extra term on the equality cone because the residual itself
vanishes there.  Equation (11) therefore says that the resulting
coefficient and its full real disk-coordinate gradient vanish on the
cone.  L155's analytic cone-division theorem gives (1), uniformly
through the singular apex; a separate cubic apex calculation is
unnecessary.

For a reflected grade `k`, the same calculation carries the explicit
mark `c^k`, so

\[
\boxed{
|\text{one-}w_k\text{ sector}|
\le C_L |c|^k {\cal Q}(z).
}                                                     \tag{20}
\]

Terms in (20) are harmless against `-a_L Q`.  The compact
two-reflection sector remains the negative L151 face.  This replaces
the false demand that the entire one-reflection sector vanish and
proves L154.

The checker also verifies a stronger ambient fact at generic complex
anchors through length ten: the gradient of the cleared determinant
form of (10) vanishes for every off-diagonal Hermitian Gram-block
variation, not just for Toeplitz variations.  That stronger ambient
identity is not needed for (11) and is not claimed all-size here.

## 6. Regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_disk_one_reflection.py \
  --output experiments/crabb_disk_one_reflection_s70223.jsonl

PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_disk_one_reflection_jets.py \
  --output experiments/crabb_disk_one_reflection_jets_s70223.jsonl
```

The first checker records (8) exactly and uses binary64 only for the
final derivative-to-`Q` diagnostic ratios.  The jet checker uses exact
Gaussian-rational algebra.  Through length ten it verifies the full
ambient cofactor gradient, the four inverse-free identities
(13)--(16), and the now-redundant subquartic apex cancellations.  These
are independent finite audits of the all-size proof above.
