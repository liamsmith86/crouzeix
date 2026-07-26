# A first defect-Schur trace is one logarithmic volume coefficient

> **Closure note (2026-07-25).**  L279 proves this note's formerly
> open complete-delay coefficient:
> \([c^{2k}]\log\mathcal V=4\|B_k\|_F^2\) in every grade.

## 1. Result (L248, 2026-07-24)

Let \(H(\varepsilon)=H(\varepsilon)^*\) be an analytic matrix series,
let \(F\) be a fixed projection, and put \(P=I-F\).  Assume

\[
F H(0)F=F
\]

and let

\[
K=PHP-PHF(FHF)^{-1}FHP                         \tag{1}
\]

be the Schur residual, acting on \(P\mathcal H\).  Then the following
formal identity is exact:

\[
\boxed{
\log\det(H+P)-\log\det(FHF|_{F\mathcal H})
=\operatorname {tr}\log(I_{P\mathcal H}+K).}       \tag{2}
\]

In particular, if \(K=O(\varepsilon^d)\), then

\[
\boxed{
\operatorname {tr}[\varepsilon^d]K
=[\varepsilon^d]\{
\log\det(H+P)-\log\det(FHF|_{F\mathcal H})\}.}      \tag{3}
\]

For a normalized final defect

\[
H=I-CC^*,\qquad C(0)=S,
\]

where \(S\) is a partial isometry with final defect \(F\), write
\(p=\operatorname {rank}P\) and

\[
D_F=F+\frac12P.
\]

The scalar on the right of (2) has the exact defect-volume form

\[
\boxed{
\mathcal V_C
=2^p\,
\frac{\det(I-C^*D_FC)}
     {\det(I-C^*FC)},}                             \tag{4}
\]

normalized by \(\mathcal V_S=1\).  More explicitly, if

\[
X=FC,\qquad Y=PC,\qquad M=I-X^*X,
\]

then

\[
\boxed{
\mathcal V_C
=2^p\det_{P\mathcal H}
\left(I-\frac12Y M^{-1}Y^*\right).}               \tag{5}
\]

There is a second exact form which interfaces directly with L240's
initial-defect factor.  Put

\[
G=YM^{-1/2},\qquad
\mathfrak D=I_{P\mathcal H}-GG^*,\qquad
H_R=I-C^*C.
\]

Then

\[
\boxed{
\mathcal V_C=\det(I+\mathfrak D),\qquad
\operatorname {tr}\mathfrak D
=-m+\operatorname {tr}(M^{-1}H_R),}               \tag{6}
\]

where \(m=\operatorname {rank}F\).  Thus the volume is the defect of
the main output channel after its final row has been whitened.  Its
first trace face is the excess normalized mass of the **initial**
defect, which is precisely the orientation in which L240 gives an
explicit half-line Gram.

Thus the final Schur square is not an additional object at trace
level.  It is exactly the whitening factor \(M^{-1}\) in (5).
One-sided final-row terms, including the unweighted
\(\mathcal B^\sharp(x)\) terms exposed in L243, must be assembled
inside this whitening before a coefficient is assigned an endpoint
cost.

Applied after L246's metric normalization, L247's remaining delayed
trace target becomes the scalar volume law

\[
\boxed{
[c^{2k}]\log\mathcal V_{\widetilde A^{\lnot X_k}}
=4\|B_k\|_F^2,}                                   \tag{7}
\]

together with the earlier vanishing needed to make this the first
Schur face.  Here
\(\widetilde A^{\lnot X_k}
=(R^{\lnot X_k})^{1/2}A(R^{\lnot X_k})^{-1/2}\).
Equation (7) is still open.  L248 replaces the cancellation of two
huge matrix series by one scalar, final-row-whitened determinant; it
does not assume the desired value \(4\).

## 2. Schur determinant proof

Relative to \(F\mathcal H\oplus P\mathcal H\), write

\[
H=\begin{bmatrix}A&B\\B^*&C\end{bmatrix}.
\]

Then

\[
H+P=\begin{bmatrix}A&B\\B^*&C+I\end{bmatrix}.
\]

Taking the determinant by the \(A\) block gives

\[
\det(H+P)=\det(A)\det(I+C-B^*A^{-1}B)
=\det(A)\det(I+K).
\]

All constant terms in the displayed determinants are invertible, so
formal logarithms exist.  This proves (2).  If \(K\) first appears in
degree \(d\), then

\[
\operatorname {tr}\log(I+K)
=\operatorname {tr}K+O(\varepsilon^{2d}),
\]

which proves (3).

## 3. Final-row whitening

For \(H=I-CC^*\),

\[
H+P=2I-F-CC^*.
\]

Since

\[
(2I-F)^{-1}=F+\frac12P=D_F,
\]

the matrix determinant lemma gives

\[
\det(H+P)
=2^p\det(I-C^*D_FC).                              \tag{7}
\]

Sylvester's determinant identity gives

\[
\det(FHF|_{F\mathcal H})
=\det(I-C^*FC).                                   \tag{8}
\]

Equations (7)--(8) prove (4).

Finally,

\[
I-C^*D_FC
=M-\frac12Y^*Y.
\]

Factor out \(M\), apply Sylvester once more, and obtain

\[
\frac{\det(I-C^*D_FC)}{\det M}
=\det_{P\mathcal H}
\left(I-\frac12YM^{-1}Y^*\right).
\]

This proves (5).  At \(C=S\), \(X=0\) and
\(YY^*=P\), so the last determinant is \(2^{-p}\) and
\(\mathcal V_S=1\).

Moreover,

\[
I-\frac12GG^*=\frac12(I+\mathfrak D),
\]

so (5) immediately gives the determinant identity in (6).  For the
trace, use

\[
Y^*Y=C^*C-X^*X=C^*C+M-I.
\]

Cyclicity then gives

\[
\begin{aligned}
\operatorname {tr}\mathfrak D
&=p-\operatorname {tr}(M^{-1}Y^*Y)\\
&=p-\dim\mathcal H
  +\operatorname {tr}\{M^{-1}(I-C^*C)\}\\
&=-m+\operatorname {tr}(M^{-1}H_R),
\end{aligned}
\]

which proves all of (6).

For an analytic metric \(R=I+O(\varepsilon)\), set
\(\widetilde A=R^{1/2}AR^{-1/2}\).  L246's quotient-congruence
argument says that the first fixed-frame defect face is unchanged by
this analytic coordinate change.  Therefore (3)--(6) may be applied
to \(\widetilde A\) when evaluating the first face of
\(R^{-1}-AR^{-1}A^*\).

## 4. Consequence for the live delayed flux

Direct coefficient decompositions of L247's metric-edge-deleted dual
slack are badly conditioned algebraically even when evaluated in
binary64.  In generic complete delays \(k=1,\ldots,6\), the raw corner
and Schur-square traces, divided by \(\|B_k\|_F^2\), were respectively

\[
\begin{array}{c|rrrrrr}
k&1&2&3&4&5&6\\ \hline
\text{corner}
&9.83&-13.53&116.69&-156.33&797.29&-4428.44\\
\text{square}
&5.83&-17.53&112.69&-160.33&793.29&-4432.44
\end{array}
\]

while their difference was \(4\) in every case.  These values are
diagnostic only.  Equations (2) and (5)--(6) explain exactly why trying to
bound either large term separately is the wrong proof architecture:
the invariant quantity is the whitened volume.

The next calculation should insert L243--L245's zero/one reflection
sectors directly into \(X=F\widetilde A\) and \(Y=P\widetilde A\),
not into the raw dual corner and Schur square separately.  The target
is to show that the first variation of the whitened initial-defect
mass in (6) is the squared doubled remote channel \(2B_k\); its trace
would then be (7).

## 5. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/associated_defect_volume_trace.py \
  --output \
  experiments/associated_defect_volume_trace_s70224.jsonl
```

The checker builds generic analytic unitary left/right paths through
partial isometries and inserts an unrelated first non-isometric
coefficient in degrees one through six.  It independently compares:

1. the first Schur trace with (2);
2. the defect determinant ratio (4); and
3. the final-row-whitened determinant (5); and
4. the initial-defect mass identity (6).

The algebra above proves the identities; the floating audit checks
orientations, normalizations, and formal-series implementation.
The tracked dataset SHA-256 is
`cf49efe0025b9336a84068e3c92802b5ec23a6bb40d34d4d15028771c7b5ecaa`.
