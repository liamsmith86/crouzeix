# Final-row whitening is an exact closed-return renewal

> **Closure note (2026-07-25).**  L277--L278 prove all lower
> coefficients vanish, and L279 evaluates this renewal's active
> edge-deleted trace as \(4\|B_k\|_F^2\) in every grade.

## 1. Result (L258, 2026-07-25)

Retain L250's edge-deleted balanced volume through a completely
delayed grade \(k\).  Put

\[
P=I-F,\qquad R_P=PR^\circ P,\qquad
Z=A(R^\circ)^{-1}A^* .
\]

L250 proves that, through degree \(2k\),

\[
R^\circ=F\oplus R_P
\]

and

\[
\mathcal V
=2^p\frac{\det(I-D^\circ Z)}{\det(I-FZ)},\qquad
D^\circ=F+\frac12R_P,\quad p=\operatorname{rank}P. \tag{1}
\]

Define the exact final-row renewal of the output Gram by

\[
\boxed{
\begin{aligned}
Z_{\rm ret}
={}&PZP\\
&+PZF\,(I_F-FZF)^{-1}FZP ,
\end{aligned}}                                    \tag{2}
\]

acting on \(P\mathcal H\), and put

\[
\boxed{\mathfrak D_{\rm ret}=I_P-R_PZ_{\rm ret}.} \tag{3}
\]

Then, coefficientwise through the live face,

\[
\boxed{\mathcal V=\det_{P\mathcal H}
(I_P+\mathfrak D_{\rm ret}).}                     \tag{4}
\]

Moreover, \(\mathfrak D_{\rm ret}\) is similar to the normalized
final-defect Schur residual.  Consequently, **if** the delayed
first-face vanishings required in L248/A194 hold, then

\[
\mathfrak D_{\rm ret}=O(c^{2k})
\]

and

\[
\boxed{
[c^{2k}]\log\mathcal V
=\operatorname{tr}[c^{2k}]\mathfrak D_{\rm ret}.} \tag{5}
\]

Equations (2)--(5) make the cancellation hidden in L248's
\(M^{-1}\) completely explicit: every passage into the final row is
paired with a passage back to the retained output, and every number
of intervening final-row returns is summed before a coefficient is
assigned.

This closes the algebraic “one-sided port” ambiguity.  It does not
yet prove L257's delay covariance: a closed return can still contain
future transfer coefficients, so the remaining theorem is the
Hardy-index/autocorrelation evaluation of (3), not another Schur
square.

## 2. Block determinant proof

Relative to \(F\mathcal H\oplus P\mathcal H\), write

\[
Z=
\begin{bmatrix}
Z_{FF}&Z_{FP}\\
Z_{PF}&Z_{PP}
\end{bmatrix}.
\]

Since \(FZ\) has zero lower block row,

\[
I-FZ=
\begin{bmatrix}
I_F-Z_{FF}&-Z_{FP}\\
0&I_P
\end{bmatrix},
\]

and hence

\[
\det(I-FZ)=\det(I_F-Z_{FF}).                       \tag{6}
\]

The numerator in (1) is

\[
I-D^\circ Z=
\begin{bmatrix}
I_F-Z_{FF}&-Z_{FP}\\
-\frac12R_PZ_{PF}&I_P-\frac12R_PZ_{PP}
\end{bmatrix}.                                    \tag{7}
\]

Take the Schur determinant of (7) through its upper-left block.
Its second factor is

\[
I_P-\frac12R_P
\left\{
Z_{PP}+Z_{PF}(I_F-Z_{FF})^{-1}Z_{FP}
\right\}
=I_P-\frac12R_PZ_{\rm ret}.                       \tag{8}
\]

The first factor cancels (6).  Since

\[
I_P-\frac12R_PZ_{\rm ret}
=\frac12(I_P+\mathfrak D_{\rm ret}),
\]

the factor \(2^p\) in (1) cancels \(2^{-p}\), proving (4).

For completeness, the normalized final defect, in the same split,
has Schur residual

\[
\begin{aligned}
K_F
&=I_P-R_P^{1/2}Z_{\rm ret}R_P^{1/2}\\
&=R_P^{-1/2}\mathfrak D_{\rm ret}R_P^{1/2}.       \tag{9}
\end{aligned}
\]

Thus \(\mathfrak D_{\rm ret}\) has exactly the same vanishing order
and coefficient trace as the normalized residual.  Under the open
first-face hypothesis above, equation (5) follows from

\[
\operatorname{tr}\log(I+\mathfrak D_{\rm ret})
=\operatorname{tr}\mathfrak D_{\rm ret}
+O(\mathfrak D_{\rm ret}^2).
\]

All inverses have invertible constant terms:
\(Z_{FF}(0)=0\), so (2) is also an identity of formal series.

## 3. Renewal interpretation

Expand the only inverse in (2):

\[
\boxed{
Z_{\rm ret}
=PZP+\sum_{n\geq0}
PZF(FZF)^nFZP.}                                   \tag{10}
\]

Every summand after \(PZP\) is an ordered closed path:

1. enter the final row through \(FZP\);
2. make \(n\) additional final-row returns through \(FZF\); and
3. exit through \(PZF\).

There is no commutation or scalarization in (10).  L251 now applies
literally to each entry and exit: both orientations use the same
analytic tail channel, including the physical factor two at the
remote row.  Thus an isolated \(\mathcal B^\sharp(x)\) term from
L243 is not a term of the scalar volume; it occurs only inside one
of the closed products in (10).

At the half-line background, L244 gives

\[
R_PZ_{\rm ret}=I_P,
\]

so \(\mathfrak D_{\rm ret}=0\).  Therefore the word-free part of
(3) vanishes identically before coefficient extraction.

## 4. Exact interface with the live recursion

In L257's one-delay model split

\[
\mathcal K_{zB}
=\mathbb C^m\oplus z\mathcal K_{z\widetilde B},
\]

the first summand is exactly the \(F\)-row eliminated in (2).
Hence L257's open recursion is equivalent to

\[
\boxed{
\operatorname{tr}[c^{2k}]
\mathfrak D_{\rm ret}(B,k)
=\operatorname{tr}[c^{2k-2}]
\mathfrak D_{\rm ret}(\widetilde B,k-1).}         \tag{11}
\]

The remaining proof must first prove the lower vanishings and then
evaluate the active closed returns in (10), using L243's zero/one
inverse-kernel filtration and matrix-inner autocorrelation.  L254
predicts the result:

\[
\operatorname{tr}[c^{2k}]\mathfrak D_{\rm ret}
=4\|P_kT_BP_k\|_{\rm HS}^2
=4\|B_k\|_F^2
\]

under complete delay.  L258 proves the renewal reduction, not this
last Hardy-index identity.

L264 subsequently permits \(R^\circ\) to be replaced, for this active
scalar coefficient, by

\[
R^\triangleright
=R^\circ+c^{2k}P([c^{2k}]R_{\rm right})P.
\]

The final row and cross blocks remain isolated, so every formula
above still applies through the target degree, while L240/L244's
complete right-half-line metric is restored.  The active left orbit
remains deleted because its initial-defect trace is not gauge-free.
This is the preferred setting for the remaining closed-path
evaluation.
