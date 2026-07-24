# The prepared sixth-order twice-delayed elliptic face

## 1. Result (L215, 2026-07-24)

Retain the partial-isometry equality data

\[
\begin{gathered}
 I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\qquad
 V^*W=0,\\
 P=2I-E+2F,\qquad T=P^{-1/2}SP^{1/2},\\
 B_j=W^*(S^*)^jV,\qquad
 {\mathfrak C}(K)=\sum_{j\ge1}B_jKB_j^*,
\qquad Q=I-E .
\end{gathered}                                      \tag{1}
\]

Assume the twice-delayed condition

\[
 B_1=B_2=0.                                        \tag{2}
\]

Continue L213--L214's exact-axis lower-tight metric gauge through
order six as specified in Section 2.  The resulting raw sixth
coefficient of the physical upper Schur complement is

\[
\boxed{
\begin{aligned}
E_{6,\rm raw}
={}&12B_3B_3^*-28{\mathfrak C}(B_3^*B_3)\\
&+2{\cal M}_T\!\left(P^{1/2}Q
\left\{
 S(S^*)^2VB_4^*B_3+
 S(S^*)^3VB_5^*B_3
\right\}\right).
\end{aligned}}                                      \tag{3}
\]

Thus the polynomial sixth-frame preparation

\[
\boxed{
\begin{aligned}
\widehat C_{6,\rm prep}
={}&-2Q\left\{
 S(S^*)^2VB_4^*B_3+
 S(S^*)^3VB_5^*B_3
\right\},\\
C_{6,\rm prep}&=P^{1/2}\widehat C_{6,\rm prep}
\end{aligned}}                                      \tag{4}
\]

gives exactly

\[
\boxed{
E_{6,\rm prep}
=12B_3B_3^*-28{\mathfrak C}(B_3^*B_3).
}                                                   \tag{5}
\]

L212's grade-three column simplifies under (2) to

\[
\widehat C_{6,\rm L212}
=-\frac72QS^3WB_3.                                  \tag{6}
\]

It changes (5) to

\[
\boxed{
E_{6,\rm final}=-16B_3B_3^*.
}                                                   \tag{7}
\]

Every column in (4),(6) is polynomial and remains analytic through
rank changes of \(B_3\).  No transfer inverse, exact kernel
projection, or pseudoinverse is used.

L215 proves the second delayed physical base and establishes the same
coercive mechanism at two consecutive grades.  It is not yet an
all-grade theorem: the higher exact-axis frame gauges must still be
organized by one generating recursion before convergence can be
claimed.

## 2. Exact-axis frame gauge through order six

The direct ellipse map through the required order is generated
exactly from L125's theta/ODE recurrence.  Write

\[
T(c)=\phi_c(T+cT^*)=\sum_{r=0}^6c^rT_r+O(c^7)
\]

and seek

\[
\begin{aligned}
P(c)&=P+\sum_{r=1}^6c^rX_r+O(c^7),\\
D(c)&=V+\sum_{r=1}^6c^rC_r+O(c^7),\\
P(c)-T(c)^*P(c)T(c)&=D(c)D(c)^*+O(c^7).             \tag{8}
\end{aligned}
\]

The zero-reflection perpendicular frame coefficients are those of

\[
\frac1{1+2c^2}
\left\{V+2\sum_{j\ge1}(-c)^j(S^*)^{2j}V\right\}.    \tag{9}
\]

The first active grade contributes the two exact-axis corrections

\[
+2c^3SWB_3-2c^5SWB_3.                              \tag{10}
\]

Consequently the complete perpendicular table is

\[
\begin{array}{c|l}
r&P^{-1/2}C_r\pmod{VM_m}\\ \hline
1&-2(S^*)^2V\\
2& 2(S^*)^4V\\
3& 4(S^*)^2V-2(S^*)^6V+2SWB_3\\
4&-4(S^*)^4V+2(S^*)^8V\\
5&-8(S^*)^2V+4(S^*)^6V-2(S^*)^{10}V-2SWB_3\\
6& 8(S^*)^4V-4(S^*)^8V+2(S^*)^{12}V .
\end{array}                                        \tag{11}
\]

As in L214, the component parallel to \(V\) is fixed recursively by
requiring the lower Schur complement of \(P(c)-I\) to vanish at every
order.  Its complementary constant block is positive definite, so
this recursion is analytic and rank-free.

The second metric coefficient is again L213's

\[
P^{-1/2}X_2P^{-1/2}=SFS^*-S^*ES.                   \tag{12}
\]

All lower Schur coefficients vanish through order six, and all upper
Schur coefficients vanish through order five.  In particular the
odd frame terms in (10) are endpoint-null at their own orders but
are essential in the nonlinear sixth coefficient.

At the repeated length-three monomial apex, \(B_3\) is unitary and
\[
SWB_3=(S^*)^2V.
\]
The third and fifth rows of (11) become respectively
\[
6(S^*)^2V,\qquad -10(S^*)^2V,
\]
which are exactly the corresponding coefficients of the all-size
elliptic-axis defect frame.  Thus (10) is not a minimum-norm
convenience; it is the terminal reflection required by the exact
axis normalization.

## 3. Ordered sixth-face reduction

At order \(r\), write the Stein equation as

\[
X_r-T^*X_rT
=G_r+VC_r^*+C_rV^*,                                \tag{13}
\]

where \(G_r\) contains all lower operator, metric, and frame
coefficients.  Substitute the exact Riemann coefficients, (11), and
the lower-parallel recursion into (13).  Before simplifying, retain
the sixth-order upper Schur square

\[
[c^6]\,
(W^*(P(c)-4I)J)
\{J(P(c)-4I)J\}^{-1}
(J(P(c)-4I)W),
\qquad J=I-WW^*.                                   \tag{14}
\]

Use only the ordered partial-isometry relations

\[
\begin{gathered}
S^*S=I-E,\quad SS^*=I-F,\quad EF=FE=0,\\
SE=0,\quad ES^*=0,\quad S^*F=0,\quad FS=0,\\
FS^*E=0,\qquad F(S^*)^2E=0,                        \tag{15}
\end{gathered}
\]

where the last line is exactly (2).  After Stein summation, group
every surviving state word by its defect endpoints.  The ordered
ledger is

\[
\begin{array}{c|c}
\text{surviving group}&\text{upper endpoint}\\ \hline
\text{oriented grade-three normal}
 &12B_3B_3^*\\
\text{right-defect transfer channel}
 &-28{\mathfrak C}(B_3^*B_3)\\
\text{future rows four and five}
 &2{\cal M}_T(P^{1/2}Q\{
 S(S^*)^2VB_4^*B_3+
 S(S^*)^3VB_5^*B_3\})\\
\text{nonzero transfer Fourier shifts}
 &0 .
\end{array}                                        \tag{16}
\]

The last row vanishes by matrix innerness:

\[
\sum_{j\ge0}B_{j+d}B_j^*=0\qquad(d\ne0).           \tag{17}
\]

The lower-parallel words cancel against the lower Schur recursion.
All \(X_2,X_3,X_4\) cross words cancel only after (14) is included.
The two uncancelled future-row words are exactly the last term of
(16).  This proves (3) without commuting copy matrices.

## 4. Preparation, L212, and the apex

Linearity of \({\cal M}_T\) makes (4) cancel the last row of (3),
proving (5).  L212 gives

\[
{\cal M}_T(C_{6,\rm L212})
=28\{{\mathfrak C}(B_3^*B_3)-B_3B_3^*\},
\]

because both lower-contamination terms contain \(B_1\) or \(B_2\).
Adding this identity to (5) proves (7).

At the repeated length-three apex,

\[
B_4=B_5=0,\qquad QS^3WB_3=0.
\]

Hence (4) and (6) both vanish, while the exact-axis raw face is
already \(-16I_m\).

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_grade_three_face.py \
  --output \
  experiments/repeated_crabb_grade_three_face_s70224.jsonl
```

The checker regenerates the direct ellipse map from the exact
theta/ODE recurrence and constructs every Stein, lower Schur, and
upper Schur coefficient through order six.  Its deterministic cases
include unstructured twice-delayed colligations, independently gauged
heterogeneous delay sums, rank-zero \(B_3\) faces, and repeated
length-three apices.  It checks (3), (5), and (7) independently.
The ordered reduction (13)--(17), not the floating audit, is the
proof.  The tracked dataset SHA-256 is
`d5502d4ffd4dc6f1f39cda9a1476dac91eaefb3b6686f72d82290215790ace13`.
