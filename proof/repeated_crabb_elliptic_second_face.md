# The oriented second elliptic face at repeated Crabb equality

## 1. Status and orientation (candidate L203, 2026-07-24)

Fix an L193 inverse-block-Toeplitz equality anchor and retain L202's
physical notation

\[
 P-T^*PT=VV^*,\qquad
 P=2I-VV^*+2WW^*.
\]

Put \(S=P^{1/2}TP^{-1/2}\).  For the real ellipse coordinate, write

\[
\begin{aligned}
 T(c)&=T+cE_1+c^2E_2+O(c^3),\\
 E_1&=T^*-T^3,\\
 E_2&=2T-(T^2T^*+TT^*T+T^*T^2)+T^5.              \tag{1}
\end{aligned}
\]

L202 identifies the only first-order partial-isometry-normal corner:

\[
 Z=W^*P^{1/2}E_1P^{-1/2}V=4B_1.                   \tag{2}
\]

Here \(Z:\mathcal V\to\mathcal W\).  Therefore the Gram on the
**upper/left** defect space is

\[
 \boxed{ZZ^*=16B_1B_1^*,}                          \tag{3}
\]

not \(Z^*Z=16B_1^*B_1\), which acts on the lower/right defect space.
The two have the same nonzero eigenvalues, so eigenvalue-only tests do
not detect this ordering error.

The candidate prepared second metric coefficient is

\[
\boxed{
 V^*XV=0,\qquad W^*XW=-ZZ^*.
}                                                   \tag{4}
\]

The accompanying linear Stein equation is solvable in every tested
noncommuting equality anchor.  This is strong evidence for (4), but a
uniform analytic solution through the rank-jumping Crabb apex has not
yet been proved.  Thus (4) is **not yet a lemma**.

## 2. Complete second Stein system

Factor the first Stein derivative in the canonical defect gauge:

\[
\begin{aligned}
 D_1&=-(E_1^*PT+T^*PE_1)=VC_1^*+C_1V^*,\\
 C_1&=(I-VV^*)D_1V+\frac12V(V^*D_1V).              \tag{5}
\end{aligned}
\]

Seek

\[
\begin{aligned}
 P(c)&=P+c^2X+O(c^3),\\
 V(c)&=V+cC_1+c^2C_2+O(c^3).
\end{aligned}
\]

The coefficient of \(c^2\) in
\(P(c)-T(c)^*P(c)T(c)=V(c)V(c)^*\) is

\[
\boxed{
X-T^*XT=F_2+VC_2^*+C_2V^*,
}                                                   \tag{6}
\]

where

\[
F_2=C_1C_1^*+E_2^*PT+T^*PE_2+E_1^*PE_1.           \tag{7}
\]

Equations (4) and (6) are a finite real-linear boundary-value
problem.  They imply that the lower endpoint stays tight to second
order and that the upper endpoint gains the negative left Gram (3).
L194 would then lift any analytic choice of its free row to an exact
lower/Stein-tight metric.

## 3. Exact scalar trace face

One part of the second face is already exact.  Define the unnormalized
L21 dual equality certificate

\[
 {\cal Z}_0=4(P^{-1}-VV^*)\succeq0.                 \tag{8}
\]

Using \(TV=0\), \(S^*W=0\), and the three-eigenvalue formula for \(P\),

\[
\boxed{
{\cal Z}_0-T{\cal Z}_0T^*=WW^*-4VV^*.
}                                                   \tag{9}
\]

Pairing (6) with \({\cal Z}_0\) kills the unknown \(C_2\) and gives

\[
\operatorname{tr}(W^*XW)-4\operatorname{tr}(V^*XV)
=\operatorname{tr}({\cal Z}_0F_2).                 \tag{10}
\]

In balanced coordinates, substitute (1), (5), and

\[
P=2I-VV^*+2WW^*,\quad
S^*S=I-VV^*,\quad SS^*=I-WW^*.
\]

A direct order-preserving reduction gives

\[
\boxed{
\operatorname{tr}({\cal Z}_0F_2)=-\operatorname{tr}(ZZ^*).
}                                                   \tag{11}
\]

Here is an explicit trace reduction.  Put
\(\dot S=P^{1/2}E_1P^{-1/2}\),
\(\ddot S=P^{1/2}E_2P^{-1/2}\), and
\[
A=(I-VV^*)[-(\dot S^*S+S^*\dot S)]V.
\]
Because \(SV=0\),

\[
A=-S^*\dot SV.
\]

The range of \(S\) is \(W^\perp\), so orthogonally decomposing
\(\dot SV\) into its \(W^\perp\) and \(W\) parts gives

\[
\|A\|_F^2+\|\dot S(I-VV^*)\|_F^2
=\|\dot S\|_F^2-\|Z\|_F^2.                         \tag{12}
\]

Consequently (11) is equivalent to

\[
\|\dot S\|_F^2+
2\operatorname{Re}\operatorname{tr}(S^*\ddot S)
=\frac34\|Z\|_F^2.                                 \tag{13}
\]

This last identity has a short partial-isometry trace proof.  Set

\[
E=VV^*,\quad F=WW^*,\quad
Q=I-E,\quad R=I-F,\quad
r=\operatorname{tr}Q,\quad
t=\|W^*S^*V\|_F^2.
\]

The relations \(S^*S=Q\), \(SS^*=R\), \(SE=0\), and \(S^*F=0\)
give

\[
J:=PS^*P^{-1}=(I+F)S^*(I+E),\qquad
\dot S=J-S^3.                                      \tag{14}
\]

Order-preserving multiplication and cyclicity of the scalar trace
then give the following five-line table:

\[
\begin{array}{c|c}
\text{term}&\text{value}\\ \hline
\|J\|_F^2&r+6m+9t\\
\operatorname{Re}\operatorname{tr}(J^*S^3)
 &\operatorname{Re}\operatorname{tr}(S^4)\\
\|S^3\|_F^2&r-2m+t\\
\operatorname{Re}\operatorname{tr}
 S^*(S^2J+SJS+JS^2)&3r+2m-t\\
\operatorname{Re}\operatorname{tr}(S^*S^5)
 &\operatorname{Re}\operatorname{tr}(S^4).
\end{array}                                        \tag{15}
\]

For example,
\(\|S^3\|_F^2=r-2m+t\) follows by applying \(S\) successively:
the first two defect losses have size \(m\), while the overlap returned
at the third step is
\(\|ESF\|_F^2=t\).  The other rows follow by expanding (14); every
factor order is retained.

Substituting (15) into

\[
\begin{aligned}
\dot S&=J-S^3,\\
\ddot S
&=2S-(S^2J+SJS+JS^2)+S^5
\end{aligned}
\]

yields

\[
\|\dot S\|_F^2+
2\operatorname{Re}\operatorname{tr}(S^*\ddot S)
=12t.
\]

Finally \(Z=4W^*S^*V\), so \(\|Z\|_F^2=16t\).  This proves
(13), hence

\[
4\operatorname{tr}\!\left[
(I-VV^*)
\left(
AA^*+\ddot S^*S+S^*\ddot S+\dot S^*\dot S
\right)
\right]
=-\|W^*\dot SV\|_F^2.                              \tag{16}
\]
This is (11).  No copy matrices are commuted.

Combining (10)--(11),

\[
\boxed{
\operatorname{tr}(W^*XW)-4\operatorname{tr}(V^*XV)
=-\|Z\|_F^2.
}                                                   \tag{17}
\]

The candidate (4) saturates this exact trace identity.  Equation (17)
explains why the coefficient is forced and why it vanishes precisely
when the grade-one transfer coefficient vanishes.

## 4. What remains

The linear map in (4),(6) changes rank at the repeated Crabb apex.  At
the apex \(B_1=0\), and its copy-valued cokernel has dimension \(m^2\);
away from a generic noncommuting anchor, only the scalar trace
compatibility (17) may remain.  Pointwise least-squares solvability
therefore does not by itself give a bounded real-analytic selection.

The next proof obligation is one of the following equivalent bridges:

1. construct \(X,C_2\) in (4),(6) by an order-safe colligation formula;
2. prove analytic divisibility of the Lyapunov--Schmidt cokernel by the
   first transfer coefficient; or
3. lift L150's orbit-complement metric to the rank-\(m\) characteristic
   colligation.

After that bridge, the same construction must be iterated on
\(\ker B_1^*\) using \(B_2,B_3,\ldots\).  L201's invertible terminal
\(B_L\) then makes the reflected flag terminate.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_elliptic_second_face.py \
  --output \
  experiments/repeated_crabb_elliptic_second_face_s70224.jsonl
```

The checker uses lengths \(2,\ldots,5\), multiplicities \(2,3\), and
genuinely noncommuting inverse-block-Toeplitz anchors.  It verifies:

1. the right transfer Gram in (2);
2. the exact trace identity (11);
3. pointwise solvability of (4),(6);
4. the zero lower endpoint; and
5. the correctly oriented upper endpoint \(-ZZ^*\).

The floating linear solve is an adversarial audit of the candidate
matrix identity, not a proof of its analytic lift.  The standard
dataset SHA-256 is
`dd6df9a351c7cd9f059cfc9050c45759446de62d62a902d2925da77ba7d62e83`.
