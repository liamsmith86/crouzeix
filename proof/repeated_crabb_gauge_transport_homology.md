# Gauge transport is a triangular endpoint homology problem

> **Route update.**  L286 subsequently solves the first instance of
> this homology with a bounded double-delay-divisible polynomial.
> L287 then proves the next even metric coefficient cannot generally
> be endpoint-null, and L288 solves that even successor under the
> correct flag-ideal condition.  L289 packages the nonlinear endpoint
> transport in one mixed-graph identity and proves that fixed length
> requires only a finite prepared jet, not an infinitely summable
> selected-column series.  The recurrence below remains the live
> framework, but it must alternate odd homology with a retained even
> face.

## 1. Result (L285, 2026-07-25)

Let

\[
\widehat T(c)=\sum_{j\ge0}c^jT_j,\qquad
D(c)=\sum_{j\ge0}c^jD_j,\qquad
M(c)=\sum_{j\ge0}c^jM_j,
\]

where

\[
T_0=S,\qquad D_0=V,\qquad
M-\widehat T^*M\widehat T=DD^*.
\]

Change the factor and metric by

\[
\widetilde D=D+C,\qquad \widetilde M=M+X,
\qquad C_0=X_0=0.
\]

Then every coefficient \(X_n\) is governed by the exact triangular
equation

\[
\boxed{
X_n-S^*X_nS
=R_n+C_nV^*+VC_n^*,}                              \tag{1}
\]

where \(R_n\) depends only on the already chosen coefficients
\(C_1,\ldots,C_{n-1}\), \(X_1,\ldots,X_{n-1}\), and the fixed
series \(D,\widehat T\).  Explicitly,

\[
\boxed{
\begin{aligned}
R_n={}&
\sum_{\substack{i+j=n\\1\le i<n}}C_iD_j^*
+\sum_{\substack{i+j=n\\1\le j<n}}D_iC_j^*\\
&+\sum_{\substack{i+j=n\\1\le i,j<n}}C_iC_j^*
+\sum_{\substack{a+b+d=n\\d<n}}T_a^*X_dT_b .
\end{aligned}}                                    \tag{2}
\]

Thus the nonlinear propagation of an earlier gauge does not require
a new nonlinear solve at each grade.  It produces one Hermitian
forcing \(R_n\), after which the new coefficient is exactly one
instance of L204's linear endpoint homology problem.

More precisely, put

\[
K_n=V^*R_nV,\qquad
C_n^\parallel=-\frac12VK_n,\qquad
R_n^0=R_n-VK_nV^*.                                \tag{3}
\]

Then

\[
X_n^0=\mathcal G_S(R_n^0)
\]

has zero lower endpoint.  For any desired upper coefficient \(H_n\),
the remaining choice is a perpendicular column \(C_n^\perp\),
\(V^*C_n^\perp=0\), satisfying

\[
\boxed{
\mathcal M_S(C_n^\perp)
=H_n-W^*X_n^0W,}                                  \tag{4}
\]

where

\[
\mathcal M_S(C)
=W^*\mathcal G_S(VC^*+CV^*)W.
\]

Equations (1)--(4) prove that A178's explicit right-ideal recursion
and L280--L282's Markov-response route are two possible solution
methods for the **same coefficientwise equation**, not unrelated
frontiers.

There is also a decisive warning.  L284's endpoint-null quartic
gauge does **not** transport inside the first delay word ideal.
The raw order-five forcing has a nonzero exact delayed quotient.
Therefore the induction cannot merely assert that previous
endpoint-null gauges preserve ideal divisibility.  It must solve
(4), with a polynomial/right-ideal representative or with a uniform
Markov-energy estimate, at every successor.

L285 is a structural reduction and a falsification of naïve ideal
invariance.  It does not solve (4) uniformly, close A178, prove the
repeated-elliptic metric sandwich, or prove Crouzeix's conjecture.

## 2. Derivation of the triangular equation

Subtract the two Stein-factor identities:

\[
\begin{aligned}
X-\widehat T^*X\widehat T
&=(D+C)(D+C)^*-DD^*\\
&=CD^*+DC^*+CC^* .                                \tag{5}
\end{aligned}
\]

At degree \(n\), the term with \(X_n\) on the left is

\[
X_n-T_0^*X_nT_0=X_n-S^*X_nS.
\]

Every other metric term contains \(X_d\) with \(d<n\).  On the
right, the only terms containing the new column \(C_n\) are

\[
C_nD_0^*+D_0C_n^*=C_nV^*+VC_n^*.
\]

Collecting all earlier terms gives (1)--(2).

## 3. Lower elimination and endpoint homology

Since \(SV=0\), Stein inversion preserves the right-defect corner:

\[
V^*\mathcal G_S(H)V=V^*HV.                        \tag{6}
\]

The forcing added by \(C_n^\parallel\) in (3) is

\[
VC_n^{\parallel *}+C_n^\parallel V^*
=-VK_nV^*.
\]

Hence \(V^*R_n^0V=0\), and (6) gives
\(V^*X_n^0V=0\).

Adding a perpendicular column changes neither this lower endpoint
nor any earlier coefficient.  Its upper response is exactly
\(\mathcal M_S\), proving (4).  The finite Fredholm condition in
L204 characterizes solvability.  L280 identifies its range with the
Markov coboundaries

\[
\operatorname{ran}(I-\Phi\Phi^*),
\]

while L281--L282 quantify the minimum state energy.  A178 seeks a
stronger polynomial/right-ideal solution of the same equation.

## 4. The first successor already leaves the delay ideal

Apply L284's order-four normalization.  With its witness and column
\(Z_4,N_4\),

\[
C_4=-N_4,\qquad X_4=-Z_4.
\]

There are no earlier changes.  If \(D_1,T_1\) are the fixed first
factor and operator coefficients, (2) gives

\[
\boxed{
R_5=
-N_4D_1^*-D_1N_4^*
-T_1^*Z_4S-S^*Z_4T_1.}                            \tag{7}
\]

Exact partial-isometry reduction gives:

\[
\begin{array}{c|c|c}
\text{object}&\text{reduced words}&\ell^1
\\ \hline
R_5&100&104\\
R_5\pmod{\mathcal I_1}&28&60\\
QR_5Q\pmod{\mathcal I_1}&20&32.
\end{array}                                       \tag{8}
\]

Moreover,

\[
R_5=R_5^*,\qquad ER_5E=0                         \tag{9}
\]

exactly.  The nonzero quotient in (8) proves

\[
\boxed{R_5\notin\mathcal I_1.}                    \tag{10}
\]

Since even the larger two-sided delay ideal does not contain \(R_5\),
no smaller one-sided/right ideal contains it.  This is the algebraic
reason the simultaneous-reuse scope guard in L284 fails.

The matrix audit makes the distinction concrete.  On an inflated
colligation with \(B_1=0\) but \(B_2\ne0\), the raw successor metric
has upper endpoint norm

\[
3.529764545853\times10^{-2},
\]

despite its lower endpoint and upper trace vanishing to roundoff.
Thus a genuine next homological correction is required.  On the
tracked double-delay sample \(B_1=B_2=0\), that endpoint vanishes to
\(1.53\times10^{-13}\); this last observation is finite numerical
evidence only, not an all-colligation theorem.

## 5. Corrected continuation

Do not try to prove that raw gauge transport stays in a delay word
ideal; (8)--(10) disprove that statement.  The next inductive object
must instead be the target in (4):

\[
H_n-W^*\mathcal G_S(R_n^0)W.                     \tag{11}
\]

A successful continuation must prove one of:

1. (11) has a bounded polynomial preimage whose selected endpoint
   lies in the preceding transfer right ideal, with finite bounds and
   positive margins through the terminal grade; or
2. (11) has a uniformly bounded-energy Markov preimage by L282's
   off-commutant flux inequality.

This is one gate with two certificate languages.  It should not be
split into competing undocumented recursions, and computing another
isolated grade does not address it.

L286 supplies the first nontrivial positive example: after the
quartic normalization, its transported quintic correction solves
(11), preserves both endpoints, and lies in the double-delay ideal.
L287 supplies the equally important negative continuation: the next
even forcing has an exact nonzero scalar endpoint trace, so it cannot
be solved with target \(H_n=0\).  The arbitrary-grade rule must choose
the physical even target and prove its positive budget, not impose
endpoint-null transport at every order.

## 6. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_gauge_transport_homology.py \
  --output \
  experiments/repeated_crabb_gauge_transport_homology_s70225.jsonl
```

The checker independently enumerates the coefficient convolution in
(5), verifies (7), Hermiticity, and (9) with exact rational word
algebra, and records the nonzero quotients in (8).  It also evaluates
the transport on unstructured, single-delay, and double-delay
colligations.  The tracked dataset has SHA-256

```text
54dc8db20c3b301f7d75978ccb126ab75741417f531ec1db4cd87eb871d7a350
```
