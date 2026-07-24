# The repeated block-Hardy elliptic first jet

## 1. Result (L202, 2026-07-24)

Fix an L193 inverse-block-Toeplitz equality anchor and use physical
coordinates.  Let \(T\) be its disk operator and let \(P\) be its
normalized rank-\(m\) Hardy metric:

\[
P-T^*PT=VV^*,\qquad V^*V=I_m,\qquad TV=0.           \tag{1}
\]

Let \(W\) be the normalized terminal copy space.  Then

\[
\boxed{P=2I-VV^*+2WW^*.}                            \tag{2}
\]

Consequently, for

\[
S=P^{1/2}TP^{-1/2},
\]

\[
\boxed{
I-S^*S=VV^*,\qquad I-SS^*=WW^*.}                   \tag{3}
\]

Thus \(S\) is a partial isometry from \(V^\perp\) onto \(W^\perp\).
This is the concrete rank-\(m\) model operator behind L193's inner
colligation.  It also proves that the normalized metric spectrum is
exactly

\[
1\ (m\text{ times}),\qquad
2\ ((L-1)m\text{ times}),\qquad
4\ (m\text{ times})                                \tag{4}
\]

throughout the noncommuting equality manifold, not only at its apex.

For a complex ellipse direction \(\gamma\), the first normalized
ellipse-pullback tangent is

\[
E_\gamma=\gamma T^*-\overline\gamma T^3.            \tag{5}
\]

Indeed \(W(T+\varepsilon\gamma T^*)\) is the real-linear image
\(z\mapsto z+\varepsilon\gamma\overline z\) of the disk.  Applying L101
to its exact mode-two support profile gives the centered Riemann expansion

\[
\phi_{\varepsilon\gamma}(z)
=z-\varepsilon\overline\gamma z^3+O(\varepsilon^2).
\]

Evaluating this at \(T+\varepsilon\gamma T^*\) gives (5).

If

\[
\dot D=-(E_\gamma^*PT+T^*PE_\gamma)
\]

is the first variation of the Stein defect while \(P\) is held fixed,
then

\[
\boxed{(I-VV^*)\dot D(I-VV^*)=0.}                  \tag{6}
\]

Equivalently, there is a matrix column \(C\) such that

\[
\dot D=VC^*+CV^*.                                   \tag{7}
\]

Hence the rank-\(m\) Stein forcing can move to first order while the
metric itself stays \(P+O(c^2)\).  In particular the lower and upper
metric endpoints have zero ordinary first elliptic derivative at
every noncommuting block-Hardy equality anchor.

This proves physical elliptic first-jet stationarity at every equality
anchor.  It does **not** prove the all-grade prepared endpoint identity:
a path can first expose a higher reflected coefficient of L201's genuine
inner transfer, and the second-order endpoint sign is still the decisive
issue.

## 2. Exact three-eigenvalue metric

Return temporarily to coefficient coordinates.  Put

\[
D_0=E_0^*KE_0,\qquad D_L=E_L^*KE_L
\]

and define the two \(K\)-orthogonal endpoint columns

\[
{\cal V}=KE_0D_0^{-1/2},\qquad
{\cal W}=KE_LD_L^{-1/2}.                            \tag{8}
\]

At an L193 equality point, \(P_I(A-S)=0\).  Substitution of
\(A=2K^{-1}HS\), followed by block-row elimination, gives the
endpoint energy identity

\[
\boxed{
K-A^*KA
={\cal V}{\cal V}^*
-\frac12A^*{\cal V}{\cal V}^*A
-{\cal W}{\cal W}^*.}                              \tag{9}
\]

Here is an order-safe verification of the elimination.  Write
\(B=\widehat H^{-1}\).  After inserting \(KA=2HS\), multiply the
interior block rows by \(B\).  Their coefficient at \((i,j)\) is

\[
B_{i+1,j+1}-B_{ij},
\]

so every interior block vanishes exactly when \(B\) is block
Toeplitz.  The two uncancelled boundary blocks are respectively the
three terms on the right of (9).  No copy matrices are interchanged.
This is the same adjacent-principal elimination as L193 (29), now
applied to the quadratic energy rather than to \(A-S\).

Define

\[
M_0=2K-{\cal V}{\cal V}^*+2{\cal W}{\cal W}^*.     \tag{10}
\]

The terminal relation \(A^*KE_L=0\) gives
\(A^*{\cal W}=0\).  Using (9),

\[
\begin{aligned}
M_0-A^*M_0A
&=2(K-A^*KA)-{\cal V}{\cal V}^*
  +A^*{\cal V}{\cal V}^*A+2{\cal W}{\cal W}^*\\
&={\cal V}{\cal V}^*.
\end{aligned}
\]

The stable Stein equation has a unique solution, so \(M=M_0\).
Conjugating by \(K^{-1/2}\) gives (2), where

\[
V=K^{-1/2}{\cal V},\qquad W=K^{-1/2}{\cal W}.
\]

The endpoint columns are orthonormal and orthogonal because they are
distinct eigenspaces of the Hermitian matrix \(P\).

## 3. Partial-isometry colligation

Equation (1), \(PV=V\), and (2) give

\[
I-S^*S
=P^{-1/2}(P-T^*PT)P^{-1/2}
=VV^*.
\]

Thus \(S\) is a partial isometry with initial space \(V^\perp\).
Also \(T^*W=0\), directly from \(A^*KE_L=0\).  Hence
\(S^*W=0\).  The kernel of \(S^*\) has dimension \(m\), so it is
exactly \(W\mathbb C^m\), proving the second identity in (3).

## 4. Algebraic first-jet cancellation

Let

\[
Q=I-VV^*=S^*S,\qquad R=I-WW^*=SS^*.
\]

The two defect spaces are orthogonal, so their projections commute.
From (2),

\[
P=3I+Q-2R,\qquad
P^{-1}=\frac34I-\frac12Q+\frac14R.                 \tag{11}
\]

In balanced coordinates, (5) becomes

\[
\dot S
=\gamma PS^*P^{-1}-\overline\gamma S^3.            \tag{12}
\]

The first defect derivative is
\(-(\dot S^*S+S^*\dot S)\).  It is enough to show that
\(QS^*\dot SQ\) is skew-Hermitian.  The partial-isometry relations

\[
SQ=S,\quad RS=S,\quad QS^*=S^*,\quad S^*R=S^*
\]

and (11) give

\[
\begin{aligned}
QS^*\dot SQ
&=\gamma S^{*2}Q-\overline\gamma QS^2Q\\
&=\gamma(QS^2Q)^*-\overline\gamma QS^2Q.
\end{aligned}                                      \tag{13}
\]

This is skew-Hermitian, proving (6).  Since a Hermitian matrix whose
compression to \(V^\perp\) is zero is supported on the \(V\) row and
column, (7) follows.  One explicit choice is

\[
C=(I-VV^*)\dot DV+\frac12V(V^*\dot DV).             \tag{14}
\]

Replacing \(V\) by \(V+\varepsilon C+O(\varepsilon^2)\) matches the Stein
defect to first order along the real path with fixed complex phase
\(\gamma\).  L194's analytic metric chart then supplies an actual
lower/Stein-tight continuation with \(P(\varepsilon)=P+O(\varepsilon^2)\).  Its upper
endpoint therefore has no term linear in the physical elliptic
coordinate.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_elliptic_first_jet.py \
  --output \
  experiments/repeated_crabb_elliptic_first_jet_s70224.jsonl
```

The checker uses strengthened genuinely noncommuting
inverse-block-Toeplitz anchors at lengths \(2,\ldots,5\) and
multiplicities \(2,3\).  It
verifies (2)--(4), (6), and the explicit defect-row factorization
(7).  The dataset SHA-256 is
`bc4ba4a6ba3d57831c79b5169db3ed90148f2d48b11d3e8750a048bbf56aa504`.
The identities above, rather than the floating audit, prove the result.
