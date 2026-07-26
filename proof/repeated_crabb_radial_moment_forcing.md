# The radial recurrence moments are forced by the exact axes

## 1. Result (L270, 2026-07-25)

The two scalar coefficient identities anticipated after L269 do not
require another theta expansion.

Assume the all-grade one-delay operator difference has L269's
structural normal form

\[
\boxed{
\Delta_k
=P_k(Q_0,\ldots,Q_k)
+R_2-R_1
+G_k-SG_kS^*,}                                   \tag{1}
\]

where

\[
Q_j=(S^*)^jS^j,\qquad R_j=S^j(S^*)^j,
\]

\[
P_k(Q_0,\ldots,Q_k)
=\sum_{j=0}^kp_{k,j}Q_j,\qquad
p_k(x)=\sum_{j=0}^kp_{k,j}x^j,                   \tag{2}
\]

and

\[
G_k=S(S^*)^{k+2}S^{k+2}S^*.                     \tag{3}
\]

Then the exact monomial axes force

\[
\boxed{p_k(1)=0,\qquad p_k'(1)=-1.}              \tag{4}
\]

Consequently, for every complete delay
\(B_1=\cdots=B_{k-1}=0\),

\[
\boxed{\operatorname {tr}\Delta_k=0.}            \tag{5}
\]

Thus once the operator-support normal form (1) is proved from
L243/L251/L258, the one-delay scalar recursion follows automatically.
There is no remaining scalar coefficient calculation.

L270 is conditional on the structural decomposition (1), which is
proved only for \(k=2\) and exactly audited through \(k=5\) in L269.
It does not prove that decomposition in arbitrary grade.

## 2. The deep divergence has zero trace in every grade

For any finite pure partial isometry with defect frame \(V\),

\[
I-S^*S=VV^*=E,\qquad SV=0.
\]

Cyclicity gives

\[
\begin{aligned}
\operatorname {tr}(G_k-SG_kS^*)
&=\operatorname {tr}(EG_k)\\
&=\left\|S^{k+2}S^*V\right\|_F^2.                \tag{6}
\end{aligned}
\]

The orthogonal-defect hypothesis gives \(FV=0\), and hence

\[
S^2S^*V=S(I-F)V=SV=0.                            \tag{7}
\]

Therefore the norm in (6) is zero for every \(k\ge0\).  The deep
boundary term found by L269 is not merely numerically trace-null; its
vanishing flux is an elementary all-grade identity.

## 3. Force the two radial moments on longer monomial axes

Fix \(k\ge2\) and take a scalar Crabb/monomial shift of any length
\(N>k+2\).  Its transfer is \(B(z)=z^N\), so \(B_j=0\) through the
grade under consideration.  The exact all-size elliptic-axis metric
has zero Schur residual away from its rank-one defect.  L241's
all-size axis comparison says that the physical boundary metric first
departs from that exact metric only in grade \(N\).  Thus the full
physical residual vanishes through grade \(k<N\).  The independently
balanced deflated tail is the length-\(N-1\) monomial, so the same
argument makes its residual vanish through its associated grade
\(k-1<N-1\).

Finally, L247 says that deleting either active metric edge changes
the corresponding trace by a multiple of
\(\lvert B_k\rvert^2=0\) (equivalently
\(\lvert\widetilde B_{k-1}\rvert^2=0\) on the tail).  Therefore both
edge-deleted closed-return traces vanish at their associated faces,
and hence

\[
\operatorname {tr}\Delta_k=0                     \tag{8}
\]

on every such axis.  This uses only L241's exact axis comparison and
L247's proved metric split, not A194.

Let \(n=N+1\) be the state dimension.  On this shift,

\[
\operatorname {tr}Q_j
=\operatorname {tr}R_j=n-j
\qquad(0\le j\le k+2).                            \tag{9}
\]

Equations (2) and (9) give

\[
\operatorname {tr}P_k
=n\,p_k(1)-p_k'(1).                               \tag{10}
\]

The middle pair in (1) has trace

\[
\operatorname {tr}(R_2-R_1)=-1,                 \tag{11}
\]

and the deep pair has trace zero by (6)--(7).  Substituting
(10)--(11) into (8) yields

\[
n\,p_k(1)-p_k'(1)-1=0                            \tag{12}
\]

for every sufficiently large integer \(n\).  The coefficient of
\(n\) and the constant term must vanish separately, proving (4).

This argument explains why the finite radial coefficients look
arithmetically irregular while their only two relevant moments are
universal.

## 4. General delayed trace

Now return to defect multiplicity \(m\) and assume
\(B_1=\cdots=B_{k-1}=0\).  L252 gives

\[
\operatorname {tr}Q_j=n-jm
\qquad(0\le j\le k+1),                            \tag{13}
\]

and \(\operatorname {tr}R_j=\operatorname {tr}Q_j\).  By (4),

\[
\operatorname {tr}P_k
=n\,p_k(1)-m\,p_k'(1)=m.                          \tag{14}
\]

Also,

\[
\operatorname {tr}(R_2-R_1)=-m,                 \tag{15}
\]

while (6)--(7) gives zero for the last pair.  Equations
(14)--(15) prove (5).

Iterating (5), if the structural form (1) is established, transports
the full active scalar face to L256's universal grade-one value
\(4\|B_k\|_F^2\).  The only remaining A194 content on this route is
therefore the all-grade support/decomposition theorem.
