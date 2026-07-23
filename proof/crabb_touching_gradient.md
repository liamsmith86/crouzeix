# Upper/lower touching and the Crabb endpoint gradient (2026-07-23)

## 1. The theorem

Fix `p=L+1>=3`, `0<c<1`, and L117's elliptic-axis pullback

\[
 T_c=\phi_c(C_p+cC_p^*).
\]

Let `B_c` be the degree-`L` Chebyshev--Blaschke product used in L116,
with its unimodular phase chosen below.  Let `q_c` be the normalized
rank-one defect vector of L117's optimal metric, and put

\[
 U_c(T)=\kappa\!\left(\sum_{n\geq0}(T^*)^nq_cq_c^*T^n\right),\qquad
 R_c(T)=\|B_c(T)\|^2.                                  \tag{1}
\]

Then both functions are differentiable near `T_c`, and

\[
 \boxed{DU_c(T_c)[Z]=DR_c(T_c)[Z]\quad\hbox{for every }Z\in M_p.} \tag{2}
\]

Consequently L118's analytic envelope satisfies

\[
 \boxed{D\Gamma_p(A_c)[Y]
 =D\|B_c(T(A))\|^2_{A=A_c}[Y],\qquad A_c=C_p+cC_p^*.} \tag{3}
\]

This turns the remaining weighted-gradient problem into one scalar
endpoint functional-calculus coefficient.  With

\[
 a_c=\sqrt{k(c^{2L})},\qquad
 D_c=\operatorname{diag}(1,c^{1/2},\ldots,c^{L/2}),
\]

choose the phase of `B_c` so that

\[
 B_c(T_c)=a_cD_cJD_c^{-1}.                             \tag{4}
\]

The largest singular value in (4) is simple, with left/right singular
vectors `e_0,e_L`.  Therefore (3) is explicitly

\[
 \boxed{D\Gamma_p(A_c)[Y]
 =2a_cc^{-L/2}\operatorname{Re}
 \left\langle e_0,\,
 D_A[B_c(T(A))]_{A_c}[Y]\,e_L\right\rangle.}           \tag{5}
\]

Since `a_c c^{-L/2}->2`, L118's sufficient pure elliptic tube gate is
equivalent to proving, uniformly for unit strong directions `Y`,

\[
 \left\langle e_0,D_A[B_c(T(A))]_{A_c}[Y]e_L\right\rangle
 =o(c^L).                                               \tag{6}
\]

The exact sparse results in L118 in fact give `O(c^(L+1))` in every
tested direction.  Equation (6), not a differentiated Stein condition
number, is now the all-size target.

## 2. Proof by a touching sandwich

For every strictly stable `T` near `T_c`, L21 and von Neumann's
inequality give

\[
 R_c(T)\leq t_*(T)\leq U_c(T).                         \tag{7}
\]

Indeed, if `STS^{-1}` is a contraction, then

\[
 \|B_c(T)\|\leq\operatorname{cond}(S),
\]

while the Gramian in (1) is a feasible Stein metric.  At `T=T_c`,
L116 and L117 prove equality throughout:

\[
 R_c(T_c)=t_*(T_c)=U_c(T_c)
 ={k(c^{2L})\over c^L}.                                \tag{8}
\]

Thus the nonnegative differentiable function `U_c-R_c` has a local
minimum zero at `T_c`, which proves (2).

For completeness, differentiability has no hidden multiplicity issue.
L117's metric is diagonal with strictly increasing entries, so its
smallest and largest eigenvalues are simple.  Its Stein solution is
analytic in `T`.  Equation (4) is anti-diagonal, with singular values

\[
 a_cc^{j-L/2},\qquad 0\leq j\leq L,
\]

so its largest singular value is also simple for `0<c<1`.

L118's defect-vector stationarity and the envelope theorem identify
the derivative of `Gamma_p+4` with that of `U_c` after composing with
the analytic Riemann pullback `A\mapsto T(A)`.  This proves (3).
The standard derivative of a simple singular value squared, applied
to (4), proves (5).

## 3. Exact polynomial descent to the two-dimensional ellipse

There is a second exact form of the same scalar extremal.  Define

\[
 P_{L,c}(z)=2c^{L/2}T_L\!\left({z\over2\sqrt c}\right), \tag{9}
\]

where `T_L` is the ordinary Chebyshev polynomial.  On the elliptic
boundary `zeta+c/zeta`, one has

\[
 P_{L,c}(\zeta+c/\zeta)=\zeta^L+c^L/\zeta^L.
\]

Thus `P_{L,c}` maps the ellipse `E_c` properly onto `E_(c^L)`.
Jacobi's multiplication formula, with the centered odd phases fixed,
is equivalently

\[
 \boxed{B_c\circ\phi_c=\phi_{c^L}\circ P_{L,c}.}        \tag{10}
\]

At the matrix level, the DCT reversal identity gives

\[
 P_{L,c}(A_c)=2c^{L/2}D_cJD_c^{-1}.                   \tag{11}
\]

This matrix is a direct sum, after a permutation, of two-dimensional
blocks with off-diagonal entries `2c^j,2c^(L-j)` (and a scalar middle
entry when `L` is even).  Their numerical-range ellipses are nested
inside the endpoint block, so

\[
 W(P_{L,c}(A_c))=E_{c^L}.                              \tag{12}
\]

Hence L116's degree-`L` extremal is literally a polynomial descent of
the `p`-dimensional Crabb axis to the outer two-dimensional elliptic
block.  The remaining derivative in (6) measures the first-order
failure of this descent when the original numerical-range ellipse is
perturbed.  This gives a second interpretation of why an endpoint
path-length or Lobatto aliasing identity should close the gate.

## 4. Explicit Blaschke product and independent audit

Put `k=k(c^2)` and `K=K(k^2)`.  Up to the harmless common phase, the
real zeros of `B_c` are

\[
 \alpha_r=\sqrt{k}\,
 \operatorname{sn}\!\left(
 K-\frac{(2r+1)K}{L}\,\middle|\,k^2\right),
 \qquad 0\leq r<L.                                    \tag{13}
\]

Thus

\[
 B_c(z)=\prod_{r=0}^{L-1}{z-\alpha_r\over1-\alpha_rz}. \tag{14}
\]

The checker

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_touching_gradient.py \
  --output experiments/crabb_touching_gradient_s70223.jsonl
```

constructs (9)--(10), differentiates every rational factor directly,
and independently differentiates the fixed-defect Stein equation.
For sizes `p=3,...,10`, three values of `c`, and eight deterministic
complex matrix directions per case, the two derivatives agree to the
reported binary64 error (worst relative derivative error `9.2e-11`).
It also regenerates (11) to `1.8e-15`.  The checker uses neither an
SDP nor a conformal-map derivative.

## 5. What remains

The touching theorem is exact in every size, but it does not itself
prove (6).  It removes the metric optimization from that proof.
The remaining mechanism must be sought in:

1. the endpoint coefficient in (5);
2. the first shape derivative of the Riemann pullback;
3. Chebyshev--Lobatto/DCT cancellation for the support harmonics of a
   matrix perturbation.

This is consistent with the sparse engine's first surviving terms and
is a substantially narrower all-size identity than the original
rank-one Gramian recurrence.
