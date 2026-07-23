# Higher-order descent on the residual `p=3` Crabb quotient (2026-07-22)

## 1. Statement and scope

Let

\[
 C_3=\begin{bmatrix}0&\sqrt2&0\\0&0&\sqrt2\\0&0&0\end{bmatrix}.
\]

After quotienting the affine-unitary orbit in L66, the mode-one residual can be represented by

\[
 R_1(\zeta)=
 \begin{bmatrix}
  -\zeta/3&0&\bar\zeta\\
  0&2\zeta/3&0\\
  0&0&-\zeta/3
 \end{bmatrix}. \tag{1}
\]

The circle action

\[
 E\longmapsto e^{-i\alpha}D_\alpha^*ED_\alpha,
 \qquad D_\alpha=\operatorname{diag}(1,e^{i\alpha},e^{2i\alpha}), \tag{2}
\]

fixes `C_3` and sends `R_1(\zeta)` to `R_1(e^{-i\alpha}\zeta)`.  Thus it suffices to take
`zeta=1` and a real perturbation parameter `epsilon`.

Let `phi_epsilon` be any Riemann map from `W(C_3+epsilon R_1(1))` to the disk.  Then

\[
 \boxed{t_*(\phi_\epsilon(C_3+\epsilon R_1(1)))<4}
 \qquad(0<|\epsilon|<\epsilon_0) . \tag{3}
\]

Sections 2--3 prove this most degenerate ray at sixth order.  Section 4 treats every fixed
mixture with the residual mode-two direction at fourth order.  The combined result is raywise;
it is not yet uniform in the transition where the mode-two coordinate tends to zero with
`epsilon`.

## 2. Exact Riemann-map expansion

Normalize the inverse map `Psi_epsilon:D->W(C_3+epsilon R_1(1))` by
`Psi_epsilon(0)=0` and positive derivative.  If `h_epsilon(w)` is the simple top eigenvalue of

\[
 \operatorname{Re}(w^{-1}(C_3+\epsilon R_1(1)))), \qquad |w|=1,
\]

the normal-angle boundary is

\[
 \zeta_\epsilon(w)=w\bigl(h_\epsilon(w)-w\partial_wh_\epsilon(w)\bigr). \tag{4}
\]

Solve the characteristic equation for `h_epsilon` order by order.  A real Fourier shift
`w -> w exp(i delta_epsilon(w))` is then determined recursively by cancelling every
nonpositive Fourier mode of (4).  The remaining analytic boundary series is `Psi_epsilon`.
Its first coefficients are

\[
 \Psi_\epsilon(w)=w+\frac5{12}\epsilon w^2
 +\epsilon^2\left(\frac{73}{576}w^3-\frac{55}{576}w\right)+O(\epsilon^3). \tag{5}
\]

Writing

\[
 T_\epsilon=\phi_\epsilon(C_3+\epsilon R_1(1))
 =C_3+\sum_{j=1}^6\epsilon^jT_j+O(\epsilon^7), \tag{6}
\]

the coefficients follow recursively from

\[
 \Psi_\epsilon(T_\epsilon)=C_3+\epsilon R_1(1). \tag{7}
\]

The exact regeneration independently checks that `T_1,T_2` agree with L63's support-resolvent
derivation.  This cross-check is important: orders three through six do not rest on a numerical
Riemann-map fit.

Changing the normalization of the Riemann map only postcomposes `T_epsilon` with a disk
automorphism.  The feasible-metric definition of `t_*` is invariant under disk automorphisms
(apply the automorphism and its inverse in the same Hilbertian metric), so the chosen gauge
does not affect (3).

## 3. Rank-one Stein certificate

For a stable matrix `T` and a cyclic vector `c`, set

\[
 P(T,c)=\sum_{j\ge0}(T^*)^jcc^*T^j. \tag{8}
\]

Then

\[
 P-T^*PT=cc^*\succeq0. \tag{9}
\]

Consequently, after scaling its least eigenvalue to one, `P` is feasible in L21 and

\[
 t_*(T)\le \kappa(P):=\lambda_{\max}(P)/\lambda_{\min}(P). \tag{10}
\]

At `epsilon=0`, take `c_0=(1,0,0)^T`; then (8) is exactly
`P_0=diag(1,2,4)`.  Use the analytic defect vector

\[
 c_\epsilon=
 \begin{bmatrix}
 1\\
 -\dfrac{\sqrt2}{24}\epsilon+\dfrac{77\sqrt2}{4608}\epsilon^3\\
 -\dfrac{193}{576}\epsilon^2
 \end{bmatrix}. \tag{11}
\]

The Stein recurrence from (9) determines `P_epsilon` through order six.  Since the extreme
eigenvalues `1` and `4` of `P_0` are simple, their formal eigenvalue series are unique.  Exact
collection gives

\[
 \boxed{\kappa(P_\epsilon)
 =4-\frac{171}{4096}\epsilon^6+O(\epsilon^7).} \tag{12}
\]

For transparency, allowing the second component of `c` to start with `a_1 epsilon` gives the
second-order coefficient

\[
 \frac{288a_1^2+24\sqrt2a_1+1}{36}, \tag{13}
\]

which vanishes uniquely at `a_1=-sqrt(2)/24`.  With that choice, allowing the third component
to start with `b_2 epsilon^2` gives the fourth-order coefficient

\[
 \frac{(576b_2+193)^2}{124416}, \tag{14}
\]

which vanishes at `b_2=-193/576`.  Finally the sixth-order coefficient before choosing `a_3`
is

\[
 \frac{10616832a_3^2-354816\sqrt2a_3-49475}{1327104}; \tag{15}
\]

its minimum occurs at `a_3=77sqrt(2)/4608` and is `-171/4096`.

The exact Gramian remains positive definite for sufficiently small `epsilon` by continuity
from `P_0`.  Equations (10) and (12) therefore prove (3).  Notice that only a feasible upper
certificate is claimed; no unproved assertion that the rank-one defect is the globally optimal
metric is needed.

## 4. Quartic descent on every other residual ray

The second complex quotient coordinate has representative

\[
 R_2(\omega)=
 \begin{bmatrix}0&0&0\\ \omega&0&0\\0&\omega&0\end{bmatrix}. \tag{16}
\]

Under (2), the residual coordinates transform as

\[
 (\zeta,\omega)\longmapsto
 (e^{-i\alpha}\zeta,e^{-2i\alpha}\omega). \tag{17}
\]

Conjugating the matrix conjugates both coordinates.  The order-`j` coefficient produced by
the analytic rank-one Gramian branch is therefore a real homogeneous polynomial invariant
under (17) and conjugation.  At degree three the only possible term is a multiple of
`Re(zeta^2 conjugate(omega))`.  At degree four the only possibilities are

\[
 |\zeta|^4,\qquad |\zeta|^2|\omega|^2,\qquad |\omega|^4. \tag{18}
\]

Exact order-four calculations on the three real paths `(zeta,omega)=(1,0),(0,1),(1,1)` give

\[
 q_3=0,\qquad
 q_4(1,0)=0,\quad q_4(0,1)=-4,\quad q_4(1,1)=-\frac{63}{8}. \tag{19}
\]

The invariant classification therefore determines the full quartic:

\[
 \boxed{q_4(\zeta,\omega)
 =-4|\omega|^4-\frac{31}{8}|\zeta|^2|\omega|^2.} \tag{20}
\]

One equivariant beginning of the corresponding defect-vector branch is

\[
 c_\epsilon=
 \begin{bmatrix}
 1\\
 -\dfrac{\sqrt2}{24}\zeta\epsilon
   +\dfrac12\bar\zeta\omega\epsilon^2+O(\epsilon^3)\\
 -2\sqrt2\omega\epsilon
   -\dfrac{193}{576}\zeta^2\epsilon^2+O(\epsilon^3)
 \end{bmatrix}. \tag{21}
\]

Consequently, for each fixed pair with `omega != 0`,

\[
 t_*(\phi_\epsilon(C_3+\epsilon(R_1(\zeta)+R_2(\omega))))
 \le 4-\left(4|\omega|^4+\frac{31}{8}|\zeta|^2|\omega|^2\right)\epsilon^4
 +O(\epsilon^5)<4 \tag{22}
\]

for all sufficiently small nonzero `epsilon`.  When `omega=0`, (12) supplies the negative
sixth-order term instead.  Thus **every fixed nontrivial straight ray in the chosen canonical
representative space for the `p=3` residual quotient strictly descends**.

This raywise result is not automatically a uniform punctured-neighbourhood theorem.  In a
joint limit where `omega` tends to zero with `epsilon`, the quartic margin in (22) competes
with higher-order mixed terms and the sixth-order pure-mode margin.  A weighted blow-up of
that transition is one remaining `p=3` local debt.  L69 now constructs the exact local
affine-unitary slice; the other debt is to control couplings between its uniformly negative
second-order directions and the residual quotient.  The linear quotient count alone does not
supply this uniform estimate.

## 5. Regeneration and remaining frontier

Run

```bash
.venv/bin/python -u experiments/p3_crabb_sixth_order.py
.venv/bin/python -u experiments/p3_crabb_quartic.py
```

The scripts derive all conformal and Stein coefficients from scratch.  The sixth-order script
checks the independent L63 coefficients and verifies (13)--(15) and (12).  The quartic script
derives the three paths in (19), performs both levels of defect minimization exactly, and uses
the invariant classification (18) to verify (20).

The next `p=3` task is the weighted transition `|omega|=O(|epsilon*zeta|)` between (20) and
(12), together with the exact-slice coupling problem just noted.  Closing both
would upgrade canonical straight-ray descent to a genuine local-neighbourhood theorem.  The
higher-size residual modes and repeated-block faces remain separate debts.
