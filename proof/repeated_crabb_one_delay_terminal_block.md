# One delay is an exact rank-\(m\) terminal-crossing extension

## 1. Result (L235, 2026-07-24)

Retain the balanced repeated-Crabb colligation

\[
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\qquad EF=0,
\]

and assume the first complete delay

\[
B_1=W^*S^*V=0.                                    \tag{1}
\]

Remove the first left wandering layer \(W\mathbb C^m\).  Let
\(J:{\cal H}_1\to{\cal H}\) include its orthogonal complement and set

\[
S_1=J^*SJ,\qquad \widetilde V=J^*V,\qquad
\widetilde W=J^*SW,
\]

\[
\widetilde E=\widetilde V\widetilde V^*,\qquad
\widetilde F=\widetilde W\widetilde W^*.
\]

Then, relative to \(W\mathbb C^m\oplus J{\cal H}_1\),

\[
\boxed{
S=\begin{bmatrix}0&0\\\widetilde W&S_1\end{bmatrix},\quad
E=\begin{bmatrix}0&0\\0&\widetilde E\end{bmatrix},\quad
F=\begin{bmatrix}I&0\\0&0\end{bmatrix}.}           \tag{2}
\]

The tail \(S_1\) is again a partial isometry with orthogonal defects
\(\widetilde E,\widetilde F\), and

\[
\boxed{\widetilde B_j
=\widetilde W^*(S_1^*)^j\widetilde V=B_{j+1}.}     \tag{3}
\]

More importantly, L228's requested grouping of terminal crossings
has an exact resolvent form.  Put

\[
\Xi_S(c)=S+c(I+F)S^*(I+E)
\]

and

\[
\Xi_-(c)=S_1+cS_1^*(I+\widetilde E).
\]

Then

\[
\boxed{
\Xi_S(c)=
\begin{bmatrix}
0&2c\widetilde W^*\\
\widetilde W&\Xi_-(c)
\end{bmatrix}.}                                    \tag{4}
\]

Consequently, whenever the displayed resolvents exist,

\[
\boxed{
J^*(zI-\Xi_S(c))^{-1}J
=\left\{
zI-\Xi_-(c)-\frac{2c}{z}\widetilde F
\right\}^{-1}.}                                    \tag{5}
\]

The insertion \(2cz^{-1}\widetilde F\) is exactly one round trip through the
removed terminal layer.  Thus

\[
\begin{aligned}
J^*(zI-\Xi_S(c))^{-1}J
={}&R_-(z)
+\sum_{r\ge1}R_-(z)
\left\{\frac{2c}{z}\widetilde F R_-(z)\right\}^{r},\\
R_-(z)&=(zI-\Xi_-(c))^{-1},                         \tag{6}
\end{aligned}
\]

as a convergent Neumann series on a sufficiently large contour, or
as a formal series.  Formula (6) groups zero, one, two, and later
terminal crossings **before** expanding the scalar Riemann map.
This is the structural reduction requested in L228.

L235 does not yet prove L228's associated-graded coefficient
identity.  It replaces the uncontrolled free-word expansion by the
single rank-\(m\) insertion in (6), together with the exact
boundary-metric blocks below.

## 2. Proof of the tail colligation

The range of \(S\) is \(F^\perp\), so the top block row in (2)
vanishes.  Since \(EF=0\), \(V\) lies in the retained space.  Also

\[
(SW)^*(SW)=W^*S^*SW=W^*(I-E)W=I,
\]

so \(\widetilde W\) is an isometric frame.  Direct compression gives

\[
S_1^*S_1=I-\widetilde E,\qquad
S_1S_1^*=I-\widetilde F.
\]

Assumption (1) is precisely

\[
\widetilde W^*\widetilde V=0,
\]

and hence \(\widetilde E\widetilde F=0\).  Finally,

\[
\begin{aligned}
\widetilde W^*(S_1^*)^j\widetilde V
&=(S^jSW)^*V\\
&=W^*(S^*)^{j+1}V=B_{j+1},
\end{aligned}
\]

which proves (3).

## 3. The terminal-crossing resolvent

Taking the adjoint of (2), multiplying by the two defect weights,
and using \(\widetilde W^*\widetilde E=0\) gives

\[
(I+F)S^*(I+E)
=\begin{bmatrix}
0&2\widetilde W^*\\
0&S_1^*(I+\widetilde E)
\end{bmatrix}.
\]

This proves (4).  Notice that \(\Xi_-\) is deliberately not the
balanced tail pencil.  The latter is

\[
\Xi_{S_1}(c)
=\Xi_-(c)+c\widetilde F S_1^*(I+\widetilde E).    \tag{7}
\]

This missing left-balance term is one reason a naïve whole-series
tail identity is false.

Take the Schur complement of the upper block \(zI_m\) in
\(zI-\Xi_S(c)\).  The product of the two off-diagonal blocks is

\[
(-\widetilde W)(zI_m)^{-1}(-2c\widetilde W^*)
=\frac{2c}{z}\widetilde F.
\]

The lower-right block of the inverse is therefore (5).  Expanding
that inverse proves (6).  For every function analytic on a contour
around the pencils, holomorphic functional calculus now gives

\[
J^*f(\Xi_S(c))J
=\frac{1}{2\pi i}\int_\Gamma
\left[
f(z)R_-(z)
+\sum_{r\ge1}f(z)R_-(z)
\left\{\frac{2c}{z}\widetilde F R_-(z)\right\}^{r}
\right]\,dz.                                       \tag{8}
\]

Taking \(f=\phi_c\) makes (8) the desired zero/one/two-crossing
organization of the direct ellipse map.

## 4. Exact boundary-metric blocks

Put \(q=c^2\), and let

\[
[q^h]P_{\rm bl}^S
=S^hF(S^*)^h
+ \sum_{d\mid h}(-1)^{h/d}(S^*)^dES^d.            \tag{9}
\]

Define

\[
g_d=(S_1^*)^d\widetilde V.
\]

L221's orbit splitting at one delay is the elementary identity

\[
(S^*)^dV=WB_d+Jg_d.                                \tag{10}
\]

Moreover,

\[
S^hW=JS_1^{h-1}\widetilde W.
\]

Substitution in (9) gives the complete block coefficient

\[
\boxed{
[q^h]P_{\rm bl}^S=
\begin{bmatrix}
\displaystyle\sum_{d\mid h}\epsilon_{h,d}B_dB_d^*
&
\displaystyle\sum_{d\mid h}\epsilon_{h,d}B_dg_d^*
\\[2mm]
\displaystyle\sum_{d\mid h}\epsilon_{h,d}g_dB_d^*
&
\displaystyle
S_1^{h-1}\widetilde F(S_1^*)^{h-1}
+\sum_{d\mid h}\epsilon_{h,d}g_dg_d^*
\end{bmatrix},
}                                                  \tag{11}
\]

where \(\epsilon_{h,d}=(-1)^{h/d}\).

Thus, under the deeper delay \(B_1=\cdots=B_\ell=0\), the removed
corner and cross row of the boundary metric have no nonconstant
coefficient below \(q^{\ell+1}\).  Equations (6) and (11) isolate all
remaining terms in the one-delay recursion: terminal excursions of
the ellipse pencil and the first nonzero cross row of the metric.

## 5. A false stronger shortcut

It is tempting to replace L228's coefficient recursion by

\[
J^*{\cal K}_S(c)J
\stackrel{?}{=}c^2{\cal K}_{S_1}(c)                 \tag{12}
\]

as a whole formal series.  This is false.  On eight deterministic
complete first delays, the required leading equality

\[
J^*[c^4]{\cal K}_S J=[c^2]{\cal K}_{S_1}
\]

holds to at worst \(1.5\times10^{-13}\), while the next even
whole-series comparison has Frobenius gaps between \(0.866\) and
\(2.016\).

The failure is predicted by (7): zero-crossing propagation uses
\(\Xi_-\), not the independently balanced tail pencil, and (6)
supplies additional terminal excursions.  Therefore the legitimate
goal remains the associated-graded statement in L228.  One must show
that, after the tail delay relations are imposed, the grouped
crossing contributions in (8), together with (11), cancel through
the first active face.

## 6. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_one_delay_block.py \
  --output \
  experiments/repeated_crabb_one_delay_block_s70224.jsonl
```

The eight records cover unstructured and inflated complete first
delays, defect multiplicities two through four, and noncommuting
tails.  They audit (2)--(5), (11), the transfer shift, the leading
residual covariance, and the failure of (12).  All structural
identities are proved above; the floating audit guards their
implementation and the falsification claim.  The tracked SHA-256 is
`17d2e1c1325862ee315d043e60f7d3f235552f2b5d997e0e4f83770e46914257`.

## 7. Next gate

Insert L125's scalar expansion only after (8).  At tail delay
\(\widetilde B_1=\cdots=\widetilde B_{\ell-1}=0\),
reduce separately:

1. the zero-crossing \(\Xi_-\) contribution;
2. one insertion \(2cz^{-1}\widetilde F\);
3. two insertions, including the Schur square;
4. the first metric cross row in (11); and
5. the remaining Neumann tail, kept as a closed resolvent remainder.

The target is that their coefficient of \(c^{2\ell+2}\) is the
embedded coefficient \(c^{2\ell}\) of the tail residual.  Do not
assume that three or more round trips are degree-irrelevant: prove
that filtration if it is true, or sum their closed remainder if it
is not.  Proving these grouped statements would establish L228 for
arbitrary grade without another fixed-grade jet.
