# An explicit analytic selection for the grade-one elliptic face

## 1. Result (L207, 2026-07-24)

L206 proves pointwise solvability of L204's oriented matrix endpoint.
The apparent rank-jump obstruction is in fact removable by one
universal polynomial formula.

Retain the balanced equality-anchor data

\[
\begin{aligned}
 S&=P^{1/2}TP^{-1/2},\\
 I-S^*S&=VV^*=:E,\qquad I-SS^*=WW^*=:F,\\
 Q&=I-E=S^*S,\qquad
 B_1=W^*S^*V.
\end{aligned}                                      \tag{1}
\]

Then the free perpendicular part of the second defect-frame motion can
be chosen as

\[
\boxed{
 \widehat C=-\frac72QSWB_1,\qquad
 C=P^{1/2}\widehat C.}                             \tag{2}
\]

It satisfies \(V^*C=0\) and, with L204's parallel column
\(C_\parallel=-VK/2\), produces a metric coefficient \(X\) obeying

\[
\boxed{
\begin{aligned}
 X-T^*XT
 &=F_2+V(C_\parallel+C)^*
       +(C_\parallel+C)V^*,\\
 V^*XV&=0,\\
 W^*XW&=-16B_1B_1^*=-ZZ^* .
\end{aligned}}                                    \tag{3}
\]

Formula (2) is polynomial in the balanced colligation data.  Those
data are real analytic in L193's inverse-block-Toeplitz coordinates,
and \(B_1=0\) at the Crabb apex.  Thus (2), and the stable Stein
solution \(X\), remain bounded and real analytic through the rank
jump.  This closes the analytic-selection debt left by L203--L206.

The identity is more general than needed: it holds for every strict
pure finite partial isometry with equal, orthogonal left and right
defects.  It does not yet treat the higher reflected grades
\(B_2,\ldots,B_L\).

## 2. Balanced second forcing

The exact balanced ellipse jets from L203 are

\[
\begin{aligned}
 J&=(I+F)S^*(I+E),\\
 \dot S&=J-S^3,\\
 \ddot S
 &=2S-(S^2J+SJS+JS^2)+S^5.                       \tag{4}
\end{aligned}
\]

The canonical first defect column and fixed second forcing are

\[
\begin{aligned}
 A&=-S^*\dot S\,V,\\
 \widehat F_2
 &=AA^*+\ddot S^*S+S^*\ddot S+\dot S^*\dot S .
\end{aligned}                                      \tag{5}
\]

Since \(P^{1/2}V=V\), L204's parallel elimination subtracts

\[
V(V^*\widehat F_2V)V^*=E\widehat F_2E.
\]

Put \(b=B_1=W^*S^*V\).  The lifted free forcing furnished by (2) is

\[
\begin{aligned}
V\widehat C^*+\widehat CV^*
 =-\frac72\{ESFS^*Q+QSFS^*E\}.                    \tag{6}
\end{aligned}
\]

Consequently the complete balanced forcing is

\[
 H=\widehat F_2-E\widehat F_2E
     +V\widehat C^*+\widehat CV^*.                 \tag{7}
\]

Use only

\[
S^*S=Q,\quad SS^*=I-F,\quad SE=0,\quad S^*F=0,
\quad EF=FE=0.                                     \tag{8}
\]

Direct order-preserving multiplication of (4)--(7) gives

\[
\begin{aligned}
H={}&-4I+4SS^*-2S^4+4(S^*)^2S^2-2(S^*)^4\\
&-\frac12S^2(S^*)^3S-S(S^*)^3S^2+2S^*S^5\\
&-\frac12S^*S^3(S^*)^2-(S^*)^2S^3S^*
  -(S^*)^3S^3+2(S^*)^5S.                          \tag{9}
\end{aligned}
\]

No copy matrices have been commuted in (9).  The coefficient \(7/2\)
in (2) is exactly what changes the raw expansion to this telescoping
form.

## 3. Endpoint functional and transfer correlations

For a state matrix \(G\), define

\[
 \Phi(G)
 =W^*{\cal G}_S(G)W
 =\sum_{j\ge0}W^*(S^*)^jGS^jW.                    \tag{10}
\]

Let

\[
B_n=W^*(S^*)^nV,\qquad B_0=W^*V=0,                \tag{11}
\]

and put

\[
\Gamma_{a,b}
=\sum_{r\ge0}(r+1)B_{a+r}B_{b+r}^*.              \tag{12}
\]

The elementary telescoping identity

\[
\begin{aligned}
B_pB_q^*
&=W^*(S^*)^p(I-S^*S)S^qW\\
&=W^*(S^*)^pS^qW
  -W^*(S^*)^{p+1}S^{q+1}W
\end{aligned}                                      \tag{13}
\]

and purity of \(S\) imply

\[
\boxed{\Phi((S^*)^aS^b)=\Gamma_{a,b}.}             \tag{14}
\]

Forward powers before an adjoint block merely reindex the orbit.  More
precisely, for \(p\le a\),

\[
\begin{aligned}
\Phi(S^p(S^*)^aS^b)&=\Gamma_{a,b+p},\\
\Phi((S^*)^bS^a(S^*)^p)&=\Gamma_{b+p,a}.           \tag{15}
\end{aligned}
\]

Equations (14)--(15) follow directly by splitting the sum in (10) at
\(j=p\), using \(W^*S=0\), and then applying
\(SS^*S=S\).  In particular, (9) gives the following completely
ordered table:

\[
\begin{array}{c|c}
\text{state word}&\Phi(\text{word})\\ \hline
I&\Gamma_{0,0}\\
SS^*&\Gamma_{1,1}\\
S^4&\Gamma_{0,4}\\
(S^*)^2S^2&\Gamma_{2,2}\\
(S^*)^4&\Gamma_{4,0}\\
S^2(S^*)^3S,\ S(S^*)^3S^2&\Gamma_{3,3}\\
S^*S^3(S^*)^2,\ (S^*)^2S^3S^*,\ (S^*)^3S^3
 &\Gamma_{3,3}\\
S^*S^5&\Gamma_{1,5}\\
(S^*)^5S&\Gamma_{5,1}.
\end{array}                                       \tag{16}
\]

Substitution into (9) yields

\[
\begin{aligned}
\Phi(H)
={}&-4\Gamma_{0,0}+4\Gamma_{1,1}
    +4\Gamma_{2,2}-4\Gamma_{3,3}\\
&-2\Gamma_{0,4}+2\Gamma_{1,5}
 -2\Gamma_{4,0}+2\Gamma_{5,1}.                    \tag{17}
\end{aligned}
\]

## 4. Exact cancellation

The diagonal part of (17) telescopes coefficientwise:

\[
-4\Gamma_{0,0}+4\Gamma_{1,1}
+4\Gamma_{2,2}-4\Gamma_{3,3}
=-4B_0B_0^*-4B_1B_1^*.                            \tag{18}
\]

The two off-diagonal differences are

\[
\begin{aligned}
-2\Gamma_{0,4}+2\Gamma_{1,5}
 &=-2\sum_{n\ge0}B_nB_{n+4}^*,\\
-2\Gamma_{4,0}+2\Gamma_{5,1}
 &=-2\sum_{n\ge0}B_{n+4}B_n^*.                    \tag{19}
\end{aligned}
\]

The genuine transfer

\[
B_H(z)=\sum_{n\ge0}B_nz^n
\]

is square matrix inner.  The fourth nonconstant Fourier coefficient of
\(B_H(e^{it})B_H(e^{it})^*=I\), and its adjoint, say exactly that both
sums in (19) vanish.  Since \(B_0=0\), equations (17)--(19) prove

\[
\boxed{\Phi(H)=-4B_1B_1^*.}                       \tag{20}
\]

This is a matrix identity; it is not obtained by scalarizing or by
commuting \(B_1\) with \(B_1^*\).

## 5. Return to physical coordinates

Let

\[
\widehat X={\cal G}_S(H),\qquad
X=P^{1/2}\widehat XP^{1/2}.                        \tag{21}
\]

Because \(SV=0\), (7) and \(V^*\widehat C=0\) give

\[
V^*\widehat XV=V^*HV=0.                           \tag{22}
\]

Also \(P^{1/2}W=2W\), while L202 gives \(Z=4B_1\).  Hence (20) gives

\[
\begin{aligned}
W^*XW
&=4W^*\widehat XW
=-16B_1B_1^*
=-ZZ^*.                                           \tag{23}
\end{aligned}
\]

Conjugating the balanced Stein equation by \(P^{1/2}\) proves the
first line of (3).  This completes the proof.

## 6. Analyticity and scope

On L193's equality chart, \(P,S,V,W\) are real analytic functions of
the positive inverse-block-Toeplitz coefficients.  Formula (2) uses
only multiplication and adjoint, so \(C\) is real analytic.  The
stable Stein inverse in (21) is real analytic as well.  No
pseudoinverse of L204's rank-changing endpoint map occurs.

At the repeated Crabb apex,

\[
B_H(z)=Uz^L,
\]

so \(B_1=0\) whenever \(L>1\), and the correction (2) vanishes
automatically.  Thus the selection is bounded through precisely the
stratum where the cokernel jumps.

L207 closes the complete **grade-one, real elliptic** selection.  The
next unresolved step is the higher reflected flag: on
\(\ker B_1^*\), derive and select the corresponding losses associated
with \(B_2,\ldots,B_L\), and then merge them with L197's later
circular-normal/Hardy residual faces.

## 7. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_elliptic_selection.py \
  --output \
  experiments/repeated_crabb_elliptic_selection_s70224.jsonl
```

The checker verifies (2)--(3) and the finite forcing expansion (9) on
all standard noncommuting equality anchors of lengths two through five
and multiplicities two and three, at two amplitudes.  It also tests 20
unstructured strict partial isometries with defect dimensions up to
three.  The latter are algebraic stress tests, not extra numerical
range assumptions.  The tracked dataset SHA-256 is
`d7b93dd5744d3cd6515994c625172503bf977bfd2c121965c00908107e4da4c5`.
