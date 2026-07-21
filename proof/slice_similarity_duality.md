# Similarity duality and modal reduction for the elliptic 4x4 slice

This note isolates the exact analytic content of the complete-2 route `L20`.  It proves an
SDP-dual trace criterion valid for every strictly stable matrix and then reduces the elliptic
slice to three real parameters.  The final trace inequality is **not yet proved**; the numerical
evidence recorded below is only a guide to that remaining step.

## 1. The similarity constant and its exact dual

Let (T\in\mathbb C^{n\times n}) satisfy (r(T)<1), and define

\[
 t_*(T)=\min\{t:\ \exists P=P^*,\ I\preceq P\preceq tI,\ T^*PT\preceq P\}.
\]

For a Hermitian matrix (D), write (D=D_+-D_-), where (D_+,D_-\succeq0) are its positive
and negative parts.

**Theorem (trace-ratio duality).**

\[
 \boxed{\displaystyle
 t_*(T)=\max\left\{1,\ \sup_{0\ne Z\succeq0}
 \frac{\operatorname{tr}(Z-TZT^*)_-}
      {\operatorname{tr}(Z-TZT^*)_+}\right\}.} \tag{1}
\]

In particular, (t_*(T)\le4) if and only if

\[
 \operatorname{tr}(Z-TZT^*)_-\le
 4\operatorname{tr}(Z-TZT^*)_+\qquad(Z\succeq0). \tag{2}
\]

**Proof.**  Consider the primal SDP

\[
 \min t\quad\text{subject to}\quad
 P-I\succeq0,\quad tI-P\succeq0,\quad P-T^*PT\succeq0. \tag{3}
\]

It is strictly feasible.  Indeed, strict stability makes

\[
 P_0=\sum_{j\ge0}(T^*)^jT^j
\]

convergent, with (P_0-T^*P_0T=I); a scalar enlargement of (P_0), followed by a strict
upper bound (tI), gives a Slater point.  Thus semidefinite strong duality applies.

Pair the three constraints with (X,Y,Z\succeq0), respectively.  Minimizing the Lagrangian
over (P,t) gives the dual

\[
 \max\ \operatorname{tr}X,
 \quad X-Y+Z-TZT^*=0,
 \quad \operatorname{tr}Y=1,
 \quad X,Y,Z\succeq0. \tag{4}
\]

Put (D=Z-TZT^*).  The first relation says (X=Y-D\succeq0), so (Y\succeq0), (Y\succeq D),
and

\[
 \operatorname{tr}X=1-\operatorname{tr}D. \tag{5}
\]

The least trace of a positive matrix (Y\succeq D) is (operatorname{tr}D_+): the choice
(Y=D_+) attains it, while compression to the positive spectral subspace of (D) proves the
opposite inequality.  Hence a fixed (Z) can occur in (4) precisely when
(operatorname{tr}D_+\le1).

The ray (Z\mapsto\lambda Z) scales both parts of (D).  If
(operatorname{tr}D_->\operatorname{tr}D_+), the best allowed scale has
(lambda\operatorname{tr}D_+=1), and (5) becomes
(operatorname{tr}D_-/\operatorname{tr}D_+).  Otherwise the zero ray gives the value (1).
For nonzero (Z\succeq0), the denominator is positive: if (D\preceq0), then
(Z\preceq TZT^*\preceq\cdots\preceq T^jZ(T^*)^j\to0), forcing (Z=0).  Formula (1) now
follows from strong duality.  Homogeneity also makes the supremum attainable after normalizing
(operatorname{tr}Z=1).  This proves (1) and (2).  (□)

If (3) has a solution with (t\le4), then (S=P^{1/2}) satisfies

\[
 \|STS^{-1}\|\le1,qquad \|S\|\,\|S^{-1}\|\le2. \tag{6}
\]

For (T=\phi(A)), von Neumann's inequality therefore gives the complete Crouzeix bound on the
corresponding ellipse.  Thus (2), proved uniformly for the slice below, would close that entire
matrix family without classifying its extremal Blaschke phases.

## 2. Parity reduction

Suppose a self-adjoint involution (J) anticommutes with (T): (JTJ=-T).  Averaging a feasible
metric with (JPJ) shows that the primal may be restricted to (PJ=JP).

The dual may also be restricted to (ZJ=JZ).  To see this without assuming that the nonlinear
ratio improves under averaging, average the full dual triple ((X,Y,Z)) in (4).  Positivity,
the linear constraint, traces, and the objective are all preserved.  Applying the positive-part
minimization from the proof to the averaged (Z) then yields an invariant ratio witness of the
same violating level.  Consequently (1) has the same supremum over invariant (Z).

In the parity basis write

\[
 T=\begin{bmatrix}0&B\\ C&0\end{bmatrix},\qquad
 Z=\begin{bmatrix}Z_o&0\\0&Z_e\end{bmatrix},\quad Z_o,Z_e\succeq0. \tag{7}
\]

Then the two (2\times2) blocks of (D=Z-TZT^*) are

\[
 D_o=Z_o-BZ_eB^*,\qquad D_e=Z_e-CZ_oC^*. \tag{8}
\]

Thus `L20` is exactly the assertion that the sum of the negative eigenvalues of these two small
blocks is at most four times the sum of their positive eigenvalues, for every pair
(Z_o,Z_e\succeq0).

## 3. Exact three-parameter normal form for the elliptic slice

Let

\[
 A=S_a+cS_a^T,\qquad a_1,a_2,a_3>0,\quad0<c<1.
\]

The diagonal similarity

\[
 L=\operatorname{diag}(1,\sqrt c,c,c^{3/2})
\]

turns (A) into the symmetric weighted shift (J=L^{-1}AL).  After ordering the odd indices
before the even indices,

\[
 J=\begin{bmatrix}0&R\\R^T&0\end{bmatrix},\qquad
 R=\sqrt c\begin{bmatrix}a_1&0\\a_2&a_3\end{bmatrix}. \tag{9}
\]

Scale so that the larger singular value of (R) is (e_1=1), and put (r=e_2/e_1\in(0,1)).
There are unique angles (0<u,v<\pi/2) with

\[
 R=U\begin{bmatrix}1&0\\0&r\end{bmatrix}V^T,qquad
 U=R_u,\quad V=R_v,qquad \tan v=r\tan u. \tag{10}
\]

The last relation is exactly the vanishing of the top-right entry in (9).  Conversely, every
((c,r,u)\in(0,1)^2\times(0,\pi/2)) gives positive slice weights through (9)-(10), so no slice
matrices are lost.

The ellipse has nome (q=c^2).  Let (m\in(0,1)) be determined by

\[
 e^{-\pi K(1-m)/K(m)}=c^2,
\]

put (k=\sqrt m), and define the two disk nodes

\[
 \tau_1=\sqrt k=m^{1/4},\qquad
 \tau_2=\sqrt k\,\operatorname{sn}\!\left(\frac{2K(m)}\pi\arcsin r\,\middle|\,m\right),
 \qquad \Sigma=\operatorname{diag}(\tau_1,\tau_2). \tag{11}
\]

Since (phi) is odd, the singular-vector calculus for (J), followed by the similarity (L),
gives the exact blocks of (T=\phi(A)):

\[
 \boxed{
 B=c^{-1/2}DU\Sigma V^TD^{-1},\qquad
 C=c^{1/2}DV\Sigma U^TD^{-1},\qquad D=\operatorname{diag}(1,c).} \tag{12}
\]

This is the promised three-real-parameter normal form.

## 4. Equivalent modal metric inequalities

Let (P=\operatorname{diag}(P_o,P_e)), set

\[
 S_o=DU,\quad S_e=DV,\quad
 G_o=S_o^TS_o=U^TD^2U,\quad G_e=S_e^TS_e=V^TD^2V,
\]

and use the congruence variables

\[
 H_o=S_o^TP_oS_o,\qquad H_e=S_e^TP_eS_e.
\]

Substitution of (12) shows that the (4\times4) SDP is exactly the following pair of
(2\times2) metric constraints:

\[
 \begin{aligned}
 G_o&\preceq H_o\preceq tG_o,&
 G_e&\preceq H_e\preceq tG_e,\\
 H_o&\succeq c\Sigma H_e\Sigma,&
 H_e&\succeq c^{-1}\Sigma H_o\Sigma.
 \end{aligned} \tag{13}
\]

Equations (2), (8), and (11)-(13) are three exact equivalent formulations of the remaining
slice theorem.  They expose where a proof must use the conformal coupling: arbitrary choices of
the two rotations and two diagonal entries in (12) do not satisfy the bound.

## 5. Numerical audit and remaining theorem

`experiments/slice_similarity_duality.py` verifies (12) against the independent nodal functional
calculus and compares the primal optimum with the trace ratio obtained from its SDP dual witness.
It also provides a deterministic modal grid.  On a (15\times9\times9=1215)-point grid with
(c\in[10^{-3},0.8]), every computed optimum was below (4); the largest was approximately (3.99977)
near the singular corner (c=10^{-3},r=0.97,u=0.03).  This is strong evidence but not a proof.

Several simple metric ansatzes fail before the true SDP optimum reaches four: diagonal metrics,
one-step observability Gramians, and forcing either contraction inequality in (13) to equality.
Those failures show that the remaining result is not a one-line Lyapunov estimate.

Subsequent exact progress is in `proof/slice_boundary_theorems.md` and
`proof/slice_upper_block_theorem.md`: the full bound is proved on the sharp (c\to0)
weighted-shift face, both modal blocks satisfy (\|B\|,\|C\|\le2) everywhere, and the focus node
obeys (k/c\le4/(1+c^2)^2).  `proof/slice_coupled_defects.md` then eliminates the two
contraction LMIs by a two-node Stein-kernel formula and reduces any hypothetical optimum above
four to rank-one/rank-one or rank-one/full coupled KKT faces.

The live analytic target is now the explicit cone inequality

\[
 \boxed{\operatorname{tr}D_-\le4\operatorname{tr}D_+
 \quad\text{for (8), with (B,C) constrained by (10)-(12).}} \tag{L20'}
\]

Proving `L20'` closes the elliptic 4x4 slice completely boundedly.  It remains a partial result
toward, not a resolution of, the general Crouzeix conjecture.
