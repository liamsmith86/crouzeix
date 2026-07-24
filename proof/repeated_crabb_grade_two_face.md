# The prepared fourth-order delayed elliptic face

## 1. Result (L214, 2026-07-24)

Retain L202--L213's balanced equality data

\[
\begin{gathered}
 I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\qquad
 V^*W=0,\\
 P=2I-E+2F,\qquad T=P^{-1/2}SP^{1/2},\\
 B_j=W^*(S^*)^jV,\qquad
 {\mathfrak C}(K)=\sum_{j\ge1}B_jKB_j^* .
\end{gathered}                                      \tag{1}
\]

Assume only the delayed condition

\[
 B_1=0.                                             \tag{2}
\]

Use L213's axis-compatible second metric gauge and continue its
lower-tight Stein recursion through order four as specified below.
If \(E_{4,\rm raw}\) is the resulting fourth coefficient of the
**upper Schur complement**, then

\[
\boxed{
\begin{aligned}
E_{4,\rm raw}
={}&4{\mathfrak C}(B_2^*B_2)-20B_2B_2^*\\
&-2{\cal M}_T\!\left(
 P^{1/2}QS^*VB_3^*B_2
\right),
\qquad Q=I-E.                                      \tag{3}
\end{aligned}}
\]

The extra fourth-frame column

\[
\boxed{
\widehat C_{4,\rm prep}
=Q\{4S^2WB_2+2S^*VB_3^*B_2\},\qquad
C_{4,\rm prep}=P^{1/2}\widehat C_{4,\rm prep}
}                                                   \tag{4}
\]

changes (3) exactly to

\[
\boxed{
E_{4,\rm prep}
=12B_2B_2^*-28{\mathfrak C}(B_2^*B_2).
}                                                   \tag{5}
\]

Finally L212's grade-two column, simplified by \(B_1=0\), is

\[
\widehat C_{4,\rm L212}
=-\frac72QS^2WB_2.                                  \tag{6}
\]

It changes (5) to the coercive oriented face

\[
\boxed{
E_{4,\rm final}=-16B_2B_2^*.
}                                                   \tag{7}
\]

No inverse of \(B_2\), kernel projection, or pseudoinverse occurs.
All frame columns are polynomial in the colligation data, so the
selection is real analytic through rank changes and through the
repeated Crabb apex.

L214 closes the first genuinely delayed physical face, including the
metric Schur term that L208--L213 had left open.  It does **not** yet
prove the corresponding formula for \(B_3,B_4,\ldots\), nor the
convergent all-grade finite-\(c\) certificate.

## 2. Complete fourth-order normalization

For the real ellipse parameter, the exact direct pullback through
order four is

\[
\begin{aligned}
\phi_c(w)
={}&w-cw^3+c^2(2w+w^5)-c^3(3w^3+w^7)\\
&+c^4(w+5w^5+w^9)+O(c^5).                          \tag{8}
\end{aligned}
\]

Put

\[
T(c)=\phi_c(T+cT^*)=\sum_{r=0}^4c^rT_r+O(c^5).
\]

Seek a physical metric and rank-\(m\) Stein frame

\[
\begin{aligned}
P(c)&=P+\sum_{r=1}^4c^rX_r+O(c^5),\\
D(c)&=V+\sum_{r=1}^4c^rC_r+O(c^5),\\
P(c)-T(c)^*P(c)T(c)&=D(c)D(c)^*+O(c^5).             \tag{9}
\end{aligned}
\]

Before the fourth-order preparation (4), take the perpendicular
balanced frame coefficients

\[
\begin{array}{c|c}
r&P^{-1/2}C_r\pmod{V M_m}\\ \hline
1&-2(S^*)^2V\\
2& 2(S^*)^4V\\
3& 4(S^*)^2V-2(S^*)^6V\\
4&-4(S^*)^4V+2(S^*)^8V .
\end{array}                                        \tag{10}
\]

These are the first four perpendicular coefficients of

\[
\frac1{1+2c^2}
\left\{V+2\sum_{j\ge1}(-c)^j(S^*)^{2j}V\right\}.    \tag{11}
\]

The component parallel to \(V\) is fixed recursively by requiring
the lower Schur complement of \(P(c)-I\) to vanish coefficientwise.
This is an ordinary analytic recursion: on \(V^\perp\), the constant
lower complementary block is

\[
Q(P-I)Q,
\]

which is positive definite.  Thus only its fixed inverse is used;
no transfer rank enters.

At order two the recursion gives exactly L213:

\[
\boxed{
P^{-1/2}X_2P^{-1/2}=SFS^*-S^*ES.
}                                                   \tag{12}
\]

The lower Schur coefficients vanish through order four.  The upper
Schur coefficients vanish through order three.  Because the
complementary block of \(P-4I\) on \(W^\perp\) is invertible, its
fourth coefficient is

\[
\begin{aligned}
E_{4,\rm raw}
={}&W^*X_4W\\
&-(W^*X_2J)
\{J(P-4I)J\}^{-1}
(JX_2W),\qquad J=I-WW^*,                           \tag{13}
\end{aligned}
\]

where the inverse is restricted to \(W^\perp\).  Formula (13), not
the bare compression \(W^*X_4W\), is the physical endpoint.

## 3. Ordered reduction of the raw face

For completeness, here is the finite reduction that gives (3).
At coefficient \(r\), let \(G_r\) denote all terms in (9) not
containing \(X_r\) or the two outer copies of \(C_r\).  Then

\[
X_r-T^*X_rT
=G_r+VC_r^*+C_rV^*.                                \tag{14}
\]

Insert (8), (10), and the lower parallel components into (14).
Use (12) before expanding the \(r=3,4\) equations.  Retain every
factor order and use only

\[
\begin{gathered}
S^*S=I-E,\quad SS^*=I-F,\quad EF=FE=0,\\
SE=0,\quad ES^*=0,\quad S^*F=0,\quad FS=0,\\
FS^*E=WB_1V^*=0.                                   \tag{15}
\end{gathered}
\]

Stein-sum the fourth forcing, insert the Schur square (13), and group
the surviving words by their two defect endpoints.  The complete
ordered ledger is

\[
\begin{array}{c|c}
\text{surviving group}&\text{upper endpoint}\\ \hline
\text{diagonal right-defect orbit}
 &4{\mathfrak C}(B_2^*B_2)\\
\text{oriented grade-two normal}
 &-20B_2B_2^*\\
\text{first future-row correlation}
 &-2{\cal M}_T(P^{1/2}QS^*VB_3^*B_2)\\
\text{nonzero Fourier shifts}
 &0 .
\end{array}                                        \tag{16}
\]

The last line follows coefficientwise from the inner transfer
identity

\[
\sum_{j\ge0}B_{j+r}B_j^*=0\qquad(r\ne0).           \tag{17}
\]

The lower-parallel terms cancel against the lower Schur recursion;
the \(X_2\)-cross terms cancel against the explicit square in (13).
Nothing else remains in (16), proving (3).  This calculation is
finite at the state-word level; the only infinite notation is the
convergent Stein/transfer channel in its first row.

## 4. The two polynomial corrections

L212 with \(k=2\) and \(B_1=0\) gives

\[
-\frac72{\cal M}_T(P^{1/2}QS^2WB_2)
=28\{{\mathfrak C}(B_2^*B_2)-B_2B_2^*\}.
\]

Therefore

\[
\boxed{
{\cal M}_T(P^{1/2}QS^2WB_2)
=8\{B_2B_2^*-{\mathfrak C}(B_2^*B_2)\}.
}                                                   \tag{18}
\]

Applying (18) to the first term in (4), while the second term in
(4) cancels the last line of (3), gives

\[
\begin{aligned}
E_{4,\rm raw}+{\cal M}_T(C_{4,\rm prep})
&=4{\mathfrak C}-20B_2B_2^*
  +32(B_2B_2^*-{\mathfrak C})\\
&=12B_2B_2^*-28{\mathfrak C},
\end{aligned}
\]

which is (5).  Applying (6) adds

\[
28\{{\mathfrak C}(B_2^*B_2)-B_2B_2^*\}
\]

and proves (7).

At the repeated length-two apex, \(B_2\) is unitary, \(B_3=0\), and

\[
QS^2WB_2=0.
\]

Hence both corrections vanish, while (3) already equals
\(-16I_m\), agreeing with the exact elliptic-axis theorem.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_grade_two_face.py \
  --output \
  experiments/repeated_crabb_grade_two_face_s70224.jsonl
```

The deterministic audit uses:

1. general delayed partial isometries constrained only by \(B_1=0\);
2. independently inflated unstructured colligations;
3. independently gauged heterogeneous delay sums; and
4. repeated grade-two Crabb apices.

It constructs every operator, metric, defect-frame, and Schur
coefficient through order four.  It checks (3), (5), and (7)
separately, as well as L213's exact second metric and vanishing of
all earlier endpoint coefficients.  The ordered reduction
(14)--(18), not the floating audit, proves the result.
The tracked dataset SHA-256 is
`13728d490d218728b6c5fe367ceacf495c63575d139c9bec1f5429736fafc65d`.
