# One full-face support theorem would force the delayed recursion

## 1. Result (L271, 2026-07-25)

The two-system normal form left open by L269 follows algebraically
from a simpler one-system support statement.

Put

\[
Q_j=(S^*)^jS^j,\qquad R_j=S^j(S^*)^j
\]

and define the four remote words

\[
\begin{aligned}
H_k={}&S^*S^{k+2}(S^*)^{k+1}
+S^{k+1}(S^*)^{k+2}S\\
&+R_{k+2}
+G_k,\\
G_k={}&S(S^*)^{k+2}S^{k+2}S^*.
\end{aligned}                                    \tag{1}
\]

Assume that the full edge-deleted closed-return face has, in every
grade, the support form

\[
\boxed{
\mathfrak D_k
=U_k(Q_0,\ldots,Q_k)
-R_1+R_k-4R_{k+1}+H_k,}                          \tag{2}
\]

where \(U_k=\sum_{j=0}^ku_{k,j}Q_j\) and

\[
\boxed{\sum_{j=0}^ku_{k,j}=0.}                   \tag{3}
\]

Then its independently balanced deflated tail satisfies

\[
\boxed{
\begin{aligned}
\mathfrak D_{k-1}^{\rm tail}
={}&\widetilde U_{k-1}(Q_0,\ldots,Q_{k-1})\\
&-R_2+R_k-4R_{k+1}\\
&+H_k-G_k+SG_kS^* ,
\end{aligned}}                                   \tag{4}
\]

for another radial polynomial \(\widetilde U_{k-1}\).  Consequently

\[
\boxed{
\mathfrak D_k-\mathfrak D_{k-1}^{\rm tail}
=P_k(Q_0,\ldots,Q_k)
+R_2-R_1+G_k-SG_kS^*,}                           \tag{5}
\]

which is exactly L269's proposed all-grade structural normal form.
L270 then forces the two moments of \(P_k\) and proves scalar
one-delay covariance.

Thus the sole A194 gate on this route is reduced further: prove the
single full-face support theorem (2)--(3), together with its earlier
vanishings.  One need not compare two independently expanded faces.

L271 is a conditional reduction.  Direct exact word expansion proves
formula (2) in grades one and two and audits it through grade five,
but it is not yet proved in arbitrary grade.

## 2. Tail power identities

Let \(P=I-F=SS^*\) and embed the deflated tail by

\[
T=SP=S^2S^*,\qquad T^*=S(S^*)^2.                 \tag{6}
\]

Its identity is \(P=R_1\).  In the complete-delay quotient
\(ES^\ell F=0\) for \(1\leq\ell\leq k-1\), together with the
partial-isometry and orthogonal-defect relations, gives for
\(1\leq j\leq k-1\),

\[
\boxed{
\begin{aligned}
(T^*)^jT^j&=Q_j+R_1-I,\\
T^j(T^*)^j&=R_{j+1}.
\end{aligned}}                                   \tag{7}
\]

For the second identity, use
\(T^j=S^{j+1}S^*\) and
\((T^*)^j=S(S^*)^{j+1}\), then collapse the middle
\(S^*S\) by the partial-isometry relations.

For the first, the same power formulas give

\[
(T^*)^jT^j=S(S^*)^{j+1}S^{j+1}S^*.
\]

For \(j=1\), expansion of \(EF=0\) gives the first line of (7).
For \(j\ge2\), expand
\(ES^{j-1}F=0\) and its adjoint and induct on \(j\); the resulting
four-term rewrite is

\[
S^*S^{j+1}S^*
=-S^{j-1}+S^jS^*+S^*S^j,
\]

with the adjoint rewrite applied to the opposite boundary.  It
reduces the displayed word to \(Q_j+R_1-I\).  Equivalently, in this
complete-delay quotient it is the compression identity
\(P Q_jP=Q_j+P-I\).  The delay hypothesis is essential here; the
identity is not asserted for a general partial isometry.

Now apply a radial polynomial with coefficient sum zero:

\[
\begin{aligned}
\sum_{j=0}^ru_jQ_j^{\rm tail}
&=u_0R_1+\sum_{j=1}^ru_j(Q_j+R_1-I)\\
&=u_0I+\sum_{j=1}^ru_jQ_j,                       \tag{8}
\end{aligned}
\]

because \(\sum_ju_j=0\).  Deflation therefore preserves the radial
support and even the displayed coefficients; it creates no extra
\(R_1\) term.

## 3. Remote-word transport

Apply (1) at tail grade \(k-1\), replacing \(S\) by \(T\).  Direct
use of (6) and \(SS^*S=S\) gives

\[
\begin{aligned}
(T^*)T^{k+1}(T^*)^k
 &=S^*S^{k+2}(S^*)^{k+1},\\
T^k(T^*)^{k+1}T
 &=S^{k+1}(S^*)^{k+2}S,\\
T^{k+1}(T^*)^{k+1}
 &=R_{k+2},\\
T(T^*)^{k+1}T^{k+1}T^*
 &=SG_kS^*.
\end{aligned}                                    \tag{9}
\]

The right-orbit part of (2) transforms by (7):

\[
-R_1^{\rm tail}
+R_{k-1}^{\rm tail}
-4R_k^{\rm tail}
=-R_2+R_k-4R_{k+1}.                              \tag{10}
\]

Equations (8)--(10) prove (4).  Subtracting it from (2) proves
(5), with \(P_k=U_k-\widetilde U_{k-1}\).

## 4. Why the zero-sum condition is natural

On the bilateral/unitary background \(E=F=0\), every \(Q_j\), \(R_j\),
and remote word equals the identity.  The physical pencil is the
ellipse boundary parametrization and its direct map is exactly \(S\),
so the closed-return defect is zero.  The four remote words cancel
the coefficients \(-1+1-4=-4\) of the displayed right-orbit terms,
and the remaining radial symbol must obey (3).  The hard part is not
this scalar condition but proving that no other nonradial support
survives the complete L243/L251/L258 assembly.

## 5. Exact audit

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/repeated_crabb_full_face_deflation.py \
  --maximum-grade 12 \
  --output \
  experiments/repeated_crabb_full_face_deflation_s70225.jsonl
```

The checker verifies every generator identity in (7), all right-orbit
shifts, and all four remote-word identities (9) by exact word
reduction in grades two through twelve.  This is an audit of the
all-grade algebraic implication, not evidence for the still-open
premise (2).
