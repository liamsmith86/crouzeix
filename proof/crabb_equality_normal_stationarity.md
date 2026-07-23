# First elliptic stationarity at every Crabb disk-equality anchor

## 1. Statement (L148, 2026-07-23)

Fix `L>=2` and a phase-palindromic equality coefficient vector

\[
u_j=\overline {u_{L-j}},\qquad
g(z)=z^L+2\sum_{j=1}^{L-1}u_jz^j.
\]

Let `(H,K,A,q,r)` be L123's coefficient realization:

\[
K=H+R^*HR,\quad A=2K^{-1}HR,\quad
q=He_0,\quad r=J\overline q.
\]

Thus `W(K^(1/2)AK^(-1/2))=D`, and the finite Blaschke product

\[
B(z)=\frac{g(z)}{g^\sharp(z)}
\]

satisfies

\[
B(A)=4e_0r^*,\qquad
K e_0=q,\qquad K e_L=r.                              \tag{1}
\]

Put `A_c=A+cA^\dagger_K`, apply the centered ellipse pullback, and
allow the numerator to vary through the analytic finite-Blaschke
factor selected by Weierstrass preparation.  If

\[
{\cal B}(c)=B_c(\phi_c(A_c)),
\]

then the simple top `K`-singular value is stationary:

\[
\boxed{
\frac d{dc}\|{\cal B}(c)\|_K^2\bigg|_{c=0}=0.
}                                                     \tag{2}
\]

By the touching sandwich, the same is true for L118's optimized
rank-one Stein envelope:

\[
\boxed{
\partial_c\Gamma_{\rm eq}(u,0)=0
\quad\hbox{for every nearby equality coefficient }u.
}                                                     \tag{3}
\]

This is stronger than L131/L120's
`D_u Gamma_eq(0,c)=0`: stationarity now holds all along the nonlinear
disk-equality family.  It removes the candidate L147 obstruction in
which arbitrarily many unreflected amplitudes multiply one reflected
leg.  It does not by itself prove the convergent reflected-Rees lift.

## 2. The exact first ellipse tangent

L125 gives

\[
\phi_c(z)=z-cz^3+O(c^2).
\]

The physical adjoint becomes the `K`-adjoint in coefficient
coordinates.  Therefore

\[
\boxed{
T(c):=\phi_c(A+cA^\dagger_K)
=A+cE+O(c^2),\qquad E=A^\dagger_K-A^3.
}                                                     \tag{4}
\]

For completeness, the first composed Faber polynomial can also be
written explicitly.  If `w` is the disk coordinate, then

\[
\Psi_c(w)=w+cw^3+O(c^2)
\]

and

\[
[c]\,G_c(\Psi_c(w))
=w^3g'(w)-\frac{g'(w)-g'(0)}w.                        \tag{5}
\]

Divide (5) by the monic `g`:

\[
[c]\,G_c(\Psi_c)=n_1+g\,o_1,\qquad \deg n_1<L.        \tag{6}
\]

Weierstrass preparation therefore gives

\[
N_c=g+cn_1+O(c^2),\qquad
B_c=N_c/N_c^\sharp.                                  \tag{7}
\]

The proof below does not depend on the particular remainder `n_1`;
every tangent to the local finite-Blaschke manifold has the same
first singular-value stationarity.

## 3. Inner-function variation is stationary

First hold `A` fixed and vary only `B_c`.  Since `W(A)=D`, the disk
case of Crouzeix's theorem gives

\[
\|B_c(A)\|_K\le2
\]

for every nearby finite Blaschke product.  Equality holds at `c=0`
by (1).  The top singular value in (1) is simple, so its derivative
exists and a differentiable function with a local maximum has zero
derivative.  Hence

\[
\frac d{dc}\|B_c(A)\|_K^2\bigg|_{c=0}=0.              \tag{8}
\]

This disposes of the prepared-numerator tangent without computing
`n_1`.

## 4. The companion endpoint identity

It remains to hold `B` fixed and vary the operator by `E`.  The right
top `K`-singular vector in (1) is `e_L`, because
`K^(-1)r=e_L`; the corresponding left vector is `e_0`.
Simple singular-value differentiation reduces the claim to

\[
\boxed{
\operatorname{Re}
\{q^*DB(A)[A^\dagger_K-A^3]e_L\}=0.
}                                                     \tag{9}
\]

Here is an all-size algebraic proof.  Put `d=g^\sharp` and
`y=d(A)^{-1}e_L`.  Rational differentiation and (1) reduce the scalar
in (9) to

\[
q^*Dg(A)[E]y-2r^*Dd(A)[E]y.                          \tag{10}
\]

For a polynomial `h`, set

\[
Y_0(h)=0,\qquad
Y_{m+1}(h)=A\,Y_m(h)+E A^m,
\]

so `Y_m=D(z^m)(A)[E]`.  Insert
`E=K^{-1}A^*K-A^3`, multiply the recurrence on the left by
`q^*=e_0^*K` and `r^*=e_L^*K`, and use L123's companion columns

\[
\begin{aligned}
Ae_0&=0,&Ae_1&=2e_0,\\
Ae_j&=e_{j-1}+2u_{j-1}(e_0-e_L)\quad(2\le j\le L),
\end{aligned}                                        \tag{11}
\]

together with `d(A)y=e_L`.  At every interior column, the contribution
from `K^{-1}A^*K` in `Y_m` cancels the contribution from `-A^3` in
`Y_(m+2)`.  The two uncancelled boundary columns are

\[
2u_j\,e_0^*K Y_jy
\quad\hbox{and}\quad
2\overline {u_{L-j}}\,e_L^*K Y_{L-j}y
\]

with opposite signs.  Pairing `j` with `L-j` cancels them because
`u_j=conjugate(u_(L-j))`.  The leading pair is exactly
`q^*Y_Ly-2r^*Y_0y`; its second term is zero and its first is cancelled
by the terminal term from `d(A)y=e_L`.  Thus the simultaneous
recurrence gives

\[
q^*Dg(A)[E]y=2r^*Dd(A)[E]y.                          \tag{12}
\]

Equations (10)--(12) prove (9).  The argument is sesquilinear, so it
applies to complex phase-palindromic coefficients, not only the real
checker slice.

Combining (8)--(9) proves (2).

## 5. Transfer to the optimized Stein envelope

At `c=0`, L123's rank-one Stein metric and the Blaschke lower bound
both equal four.  Near the equality point they are differentiable
and sandwich the optimal similarity square:

\[
\|B_c(T(c))\|^2
\le t_*(T(c))
\le\Gamma_{\rm eq}(u,c)+4.
\]

The difference between the two analytic outer functions is
nonnegative and zero at `c=0`, so their first derivatives agree.
Equation (2) therefore proves (3).

## 6. Consequence for L147

L147 needs the square of the reflected ideal

\[
{\mathfrak r}=(c^L,c^ku_k,c^k\overline {u_k})_k.
\]

L148 proves the geometric input that was absent from the first draft:
the coefficient of one normal/reflected leg vanishes at **every**
unreflected equality anchor.  The remaining task is algebraic rather
than variational:

1. construct the convergent Rees lift of the prepared functional
   calculus and the recentered Stein critical equations; and
2. use (3) on its exceptional divisor to conclude membership in
   `mathfrak r^2`.

No new coefficient calculation is required after that lift.

## 7. Exact regeneration

Run

```bash
.venv/bin/python -u \
  experiments/crabb_equality_normal_stationarity.py \
  --output experiments/crabb_equality_normal_stationarity_s70223.jsonl
```

The checker independently forms (4)--(7), evaluates both rational
Fréchet derivatives, and verifies:

1. `B(A)=4e_0r^*`;
2. the inner-function contribution in (8) is zero;
3. the operator contribution in (9) is zero; and
4. their total is zero.

It uses independent symbolic equality coefficients for `L=2,3,4`
and two deterministic exact rational directions for every
`2<=L<=10`, for 21 exact records.
