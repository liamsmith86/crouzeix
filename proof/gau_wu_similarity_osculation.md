# Exact second-order similarity/scalar osculation at Gau--Wu models

> **Status and scope.**  The finite quadratic reduction, its universal
> free curvature, and the scalar/similarity osculation identity below
> are exact and dimension-free for every fixed finite nondegenerate
> Gau--Wu disk model.  They do **not** prove the remaining sign of the
> common Hessian in arbitrary dimension.  Thus this is a unification
> and simplification of the live local gate, not a proof of the
> Crouzeix conjecture or of its completely bounded version.

## 1. The universal sharp metric

Let \(S=S_\phi\) be the \(n\)-dimensional compressed shift for
\(\phi=zf\), in the upper-triangular Takenaka basis.  Write
\[
 e_0=(1,0,\ldots,0)^T,\qquad e_L=(0,\ldots,0,1)^T,
 \quad L=n-1.
\]
Then
\[
 I-S^*S=e_0e_0^*,\qquad I-SS^*=e_Le_L^*.
\]
Put
\[
 D=\operatorname{diag}(\sqrt2,I_{n-2},1/\sqrt2),
 \qquad A=DSD^{-1}.
\]
This is the Gau--Wu disk matrix and
\(f(A)=2e_0e_L^*\).  Its sharp condition-square-four metric is
\[
 P_0=2D^{-2}=\operatorname{diag}(1,2,\ldots,2,4).
 \tag{1}
\]
Two exact rank-one identities are
\[
 P_0-A^*P_0A=e_0e_0^*                                      \tag{2}
\]
and, for
\[
 W=2D(I-e_0e_0^*)D=\operatorname{diag}(0,2,\ldots,2,1),
\]
\[
 \boxed{W-AWA^*=e_Le_L^*-4e_0e_0^*.}                       \tag{3}
\]
Thus the same sharp metric and the same strictly positive embedded
dual weight occur at every Gau--Wu model, independently of the zeros
of \(f\).

## 2. Exact reduction of the second-order metric SDP

Take the normalized inverse-Riemann jet
\[
 T_\varepsilon=A+\varepsilon G+\varepsilon^2H+o(\varepsilon^2)
 \tag{4}
\]
and expand a metric and its upper bound as
\[
 P_\varepsilon=P_0+\varepsilon Z+\varepsilon^2Y,
 \qquad t_\varepsilon=4+\varepsilon^2\eta.
 \tag{5}
\]
The first-order active equations are
\[
 Z_{00}=Z_{LL}=0,\qquad
 \Pi\Delta_1\Pi=0,                                    \tag{6}
\]
where \(\Pi=I-e_0e_0^*\) and
\[
 \Delta_1=Z-A^*ZA-Q,\qquad
 Q=G^*P_0A+A^*P_0G.                                  \tag{7}
\]
For a normalized Gau--Wu tangent, the first-order stationary trace is
\(\operatorname{tr}(WQ)=0\).  This is also the first variation of the
transferred rank-one equality in L127.

Let
\[
 x=e_0^*Z\Pi\in\mathbb C^{n-1}.
\]
Writing
\[
 A=\begin{bmatrix}0&a\\0&B\end{bmatrix},
\]
equation (6) determines the rest of \(Z\) uniquely:
\[
 Z_{\Pi\Pi}-B^*Z_{\Pi\Pi}B
 =Q_{\Pi\Pi}+a^*xB+B^*x^*a.                         \tag{8}
\]
The spectrum of \(B\) is strictly inside the disk, so (8) is an
ordinary stable Stein equation.  Taking the \(W\)-trace and using
(3) and \(\operatorname{tr}(WQ)=0\) gives \(Z_{LL}=0\)
automatically.  Hence \(x\) is the complete free first-metric
coordinate: there are exactly \(2n-2\) real variables.

Define
\[
\begin{aligned}
 R={}&H^*P_0A+A^*P_0H+G^*P_0G\\
    &\quad+G^*ZA+A^*ZG,\\
 d={}&\Pi\Delta_1e_0,\\
 \alpha(Z)={}&
 Z_{\Pi0}^*(P_0-I)_{\Pi\Pi}^{-1}Z_{\Pi0},\\
 \beta(Z)={}&
 Z_{\{0,\ldots,L-1\},L}^*
 (4I-P_0)_{\{0,\ldots,L-1\}}^{-1}
 Z_{\{0,\ldots,L-1\},L}.
\end{aligned}                                                    \tag{9}
\]
The three second-order Schur complements say
\[
\begin{aligned}
 Y_{00}&\geq\alpha(Z),\\
 \eta-Y_{LL}&\geq\beta(Z),\\
 \Pi(Y-A^*YA-R)\Pi&\succeq dd^*.
\end{aligned}                                                    \tag{10}
\]
By (3),
\[
 \operatorname{tr}W(Y-A^*YA)=Y_{LL}-4Y_{00}.          \tag{11}
\]
Consequently every feasible second-order metric obeys
\[
 \eta\geq {\cal F}_C(x)
 :=4\alpha(Z)+\beta(Z)+
 \operatorname{tr}W(R+dd^*).                         \tag{12}
\]
This bound is attained.  Set \(Y_{00}=\alpha(Z)\), choose
\(Y_{0\Pi}=0\), and solve the stable compressed Stein equation in
the last line of (10) with equality.  Equation (11) then gives its
terminal entry, and setting
\(\eta=Y_{LL}+\beta(Z)\) makes the upper Schur complement tight.
Therefore the full matrix SDP is exactly
\[
 \boxed{e_\phi(C)=\min_{x\in\mathbb C^{n-1}}{\cal F}_C(x).}
 \tag{13}
\]
No semidefinite variable remains.

## 3. Universal positive curvature in the metric variable

For \(G=H=0\), equation (8) is the homogeneous displacement equation.
After conjugating by \(D\), \(\widetilde Z=DZD\) satisfies
\[
 \Pi(\widetilde Z-S^*\widetilde ZS)\Pi=0.             \tag{14}
\]
Thus \(\widetilde Z\) is a selfadjoint truncated Toeplitz operator.
The canonical model conjugation sends \(e_0\) to \(e_L\) and fixes
selfadjoint truncated Toeplitz operators.  It therefore identifies
the norms and the reversed endpoint coordinates of
\(\widetilde Ze_0\) and \(\widetilde Ze_L\).

Write \(x=(x_{\cal M},x_L)\), where \(x_{\cal M}\in\mathbb C^{n-2}\).
The three terms in (12) are then
\[
\begin{aligned}
 4\alpha&=4\|x_{\cal M}\|^2+\frac43|x_L|^2,\\
 d^*Wd&=2\|x_{\cal M}\|^2+|x_L|^2,\\
 \beta&=2\|x_{\cal M}\|^2+\frac13|x_L|^2.
\end{aligned}                                                    \tag{15}
\]
Hence
\[
 \boxed{
 {\cal F}_0(x)=8\|x_{\cal M}\|^2+\frac83|x_L|^2.}
 \tag{16}
\]
The free curvature is strictly positive and completely independent
of the zeros or dimension.  This generalizes the \(p=3,4\)
eliminations in L62--L65 without a grade-by-grade calculation.

## 4. Exact square gap to the scalar Hessian

Let \(h\) be an arbitrary degree-preserving tangent of \(f\), and put
\[
 f_\varepsilon(T_\varepsilon)
 =2e_0e_L^*+\varepsilon F_1+\varepsilon^2F_2+\cdots.
\]
Let \({\cal J}_C(h)\) denote the coefficient of
\(\varepsilon^2\) in its ordinary operator norm.  Thus
\(4{\cal J}_C(h)\) is the coefficient in the norm square.

The first motions of the lower and upper eigenvectors of \(P_0\)
induced by \(Z\) are
\[
\begin{aligned}
 a_-&=-(P_0-I)_{\Pi\Pi}^{-1}Z_{\Pi0},\\
 a_+&=(4I-P_0)_{\{0,\ldots,L-1\}}^{-1}
       Z_{\{0,\ldots,L-1\},L}.
\end{aligned}                                                    \tag{17}
\]
Define the two matching residuals
\[
\begin{aligned}
 r_-&=(F_1e_L)_{\Pi}-2a_-,\\
 r_+&=(F_1^*e_0)_{\{0,\ldots,L-1\}}-2a_+.
\end{aligned}                                                    \tag{18}
\]
Expanding the two endpoint Schur complements and using L127 to
transfer the contraction defect gives the exact identity
\[
\boxed{\begin{aligned}
 {\cal F}_C(x)-4{\cal J}_C(h)
 ={}&r_-^*(P_0-I)_{\Pi\Pi}r_-\\
 &+r_+^*
 \left[(4I-P_0)P_0^{-1}\right]_{\{0,\ldots,L-1\}}
 r_+.
\end{aligned}}                                                   \tag{19}
\]
All second-order acceleration terms cancel in (19).  In coordinates,
the two positive weights are
\[
 \operatorname{diag}(1,\ldots,1,3),
 \qquad
 \operatorname{diag}(3,1,\ldots,1).                 \tag{20}
\]
Equation (19) is the sought canonical, data-independent factor:
the similarity/scalar gap is exactly the failure of the first
Blaschke image to carry the moving upper metric eigenvector to twice
the moving lower one, together with the adjoint failure.

For fixed \(C\), the map
\[
 (x,h)\longmapsto(r_-,r_+)
 \tag{21}
\]
is a square real-linear map of dimension \(4n-4\).  It is injective.
Indeed, at \(C=0\), a vector in its kernel makes the right side of
(19) zero.  But (16) is strictly positive for \(x\ne0\), while
L337--L338 give \({\cal J}_0(h)<0\) for every nonzero zero tangent.
Thus \(x=h=0\).  Hence (21) is bijective.

For every physical \(C\), there is therefore a unique pair
\((x_C,h_C)\) for which both residuals vanish.  Since (19) holds for
every \(x,h\), that pair simultaneously minimizes the left metric
problem and maximizes the scalar problem.  We obtain the exact
osculation theorem
\[
\boxed{
 e_\phi(C)
 =\min_x{\cal F}_C(x)
 =4\max_h{\cal J}_C(h).}                             \tag{22}
\]

This also identifies the actual second coefficient of the similarity
value, rather than only its tangent-program upper bound.  The strict
second-order lift from L62 gives
\[
 \limsup_{\varepsilon\to0}
 {t_*(T_\varepsilon)-4\over\varepsilon^2}
 \le e_\phi(C).
\]
For the matching inner tangent \(h_C\),
\[
 t_*(T_\varepsilon)
 \ge\|f_\varepsilon(T_\varepsilon)\|^2
 =4+4{\cal J}_C(h_C)\varepsilon^2+o(\varepsilon^2).
\]
Equation (22) makes the two bounds equal.  Hence the two-sided
second derivative exists along every analytic normalized curve and
is \(e_\phi(C)\).

## 5. Meaning of the theorem

Paulsen's similarity quantity is generally stronger than the scalar
Crouzeix quantity.  Equation (22) says that this extra strength has
**no quadratic cost at a finite Gau--Wu equality model**: the two
envelopes have exactly the same transverse Hessian after the natural
factor four.

This does not supply the remaining sign.  It replaces it by the
equivalent finite problem
\[
 \boxed{\min_x{\cal F}_C(x)\leq0\quad\hbox{for every physical }C.}
 \tag{23}
\]
The advantage over the previous formulation is that the variable
Hessian is the universal diagonal form (16), and the exact gap
factor (19) explains L341's two endpoint observations.  A proof of
(23) would establish both the scalar and completely bounded local
second-order inequalities at every nondegenerate Gau--Wu model.

## 6. Deterministic audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_similarity_hessian.py \
  --dimensions 3,4,5,6,7,8 --samples 2 --angle-count 512 \
  --output experiments/gau_wu_similarity_hessian_s70224.jsonl
```

The checker independently:

1. verifies (1)--(3);
2. constructs (12) from the finite Stein solve;
3. reproduces (16);
4. reconstructs both sides of the full square identity (19);
5. checks that (21) is nonsingular;
6. compares (13) with the unreduced CVXPY SDP on three directions;
7. compares (22) with the independently built scalar Hessian.

On twelve models, two in each dimension \(3,\ldots,8\), the largest
operator-norm residual in (19) is \(1.27\cdot10^{-13}\), the smallest
singular value of (21) is \(0.937\), and the universal block (16)
agrees within \(7.62\cdot10^{-15}\).  The unreduced SDP agrees within
\(3.24\cdot10^{-7}\).  The direct Schur-complement comparison in
(22) is less well conditioned near a zero collision but still agrees
within \(3.49\cdot10^{-8}\).

All twelve common Hessians are negative in the finite tests.  Eleven
are strictly separated at the conservative reporting threshold; the
soft \(n=8\) sample has largest eigenvalue
\(-6.42\cdot10^{-8}\) and is deliberately classified as numerically
unresolved.  These signs reproduce L334 and are not used as a proof.

The dataset SHA-256 is
`30f38471eea4fb3728fe3a10b154053fcc44cbbc9c13b5673f125a27fb6a8ae1`.

The finite sign tests remain evidence only; the exact statements are
the algebraic identities above.
