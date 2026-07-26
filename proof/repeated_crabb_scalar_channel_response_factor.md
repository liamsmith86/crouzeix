# Two-sided scalar-channel defect and response factorization

> **Campaign scope.**  This is a local factorization theorem on the
> scalarized zero-Jensen branch of the repeated disk chart.  It does
> not by itself complete the full curve-selection induction.

## 1. Result (L323, 2026-07-26)

Let

\[
B(z)=\sum_{n\ge1}B_nz^n
\]

be one square matrix-inner transfer from L201.  For unit defect
coordinates \(u,v\), put

\[
b(z)=u^*B(z)v,\qquad
\delta(u,v)=1-\sum_n|u^*B_nv|^2.                  \tag{1}
\]

Choose orthogonal complements and write

\[
B(z)=
\begin{bmatrix}
b(z)&r(z)\\
c(z)&D(z)
\end{bmatrix}:
\mathbb Cv\oplus v^\perp\longrightarrow
\mathbb Cu\oplus u^\perp.                          \tag{2}
\]

Then matrix-inner unitarity gives the exact two-sided identities

\[
\boxed{
\|r\|_{H^2}^2=\|c\|_{H^2}^2
=\delta(u,v),}                                     \tag{3}
\]

and, pointwise almost everywhere,

\[
\boxed{
1-|b(\zeta)|^2
=\|r(\zeta)\|^2
=\|c(\zeta)\|^2\ge0.}                              \tag{4}
\]

Consequently every nonconstant Fourier coefficient of the scalar
defect is bounded by its total mass:

\[
\boxed{
\left|\widehat{(1-|b|^2)}(k)\right|
\le\delta(u,v)\qquad(k\in\mathbb Z).}              \tag{5}
\]

Every bounded Hardy cross pairing is controlled at the same order.
For example, if \(Q(\zeta)\) is a bounded compatible matrix symbol,

\[
\boxed{
\left|
\frac1{2\pi}\int_{\mathbb T}
r(\zeta)Q(\zeta)c(\zeta)\,|d\zeta|
\right|
\le\|Q\|_\infty\delta(u,v).}                       \tag{6}
\]

At a maximizing pair for L320,

\[
\delta(u,v)=1-\sigma(B).                           \tag{7}
\]

Finally, on any fixed finite matrix-Schur chart near a split scalar
channel, every real-analytic, defect-frame-gauge-invariant scalar
response coefficient on a fixed curve-selected extremal branch which

1. vanishes on the block-diagonal inner locus, and
2. is evaluated on L199's scalar zero-Jensen winner,

obeys

\[
\boxed{
|{\cal R}(B)|\le C\{1-\sigma(B)\}.}                 \tag{8}
\]

The constant is local and may depend on the fixed Crabb length and
copy multiplicity.  The complete square against any fixed positive
normal curvature \(\gamma\) therefore costs only

\[
\frac{|{\cal R}(B)|^2}{\gamma}
=O\bigl((1-\sigma(B))^2\bigr),                     \tag{9}
\]

which is absorbed by L322's
\(\frac43(1-\sigma(B))\) reserve after one neighbourhood shrink.

Equation (8) applies directly to the pure scalar circular-normal
response at an L193 anchor.  For a joint disk/normal response, apply
it to the **difference from L192's already controlled split-channel
response**.  That difference is analytic, physical and hence
defect-frame-gauge invariant, and vanishes on every exactly split
channel.  Thus the formerly missing later transport correction need
not be computed coefficient by coefficient.  What remains is to
assemble this estimate through the finite L197/L321 channel flag and
the separately retained L318 elliptic margin.

## 2. Exact channel defect identities

On the circle, \(B(\zeta)\) is unitary.  Its first column and first row
in (2) have norm one:

\[
|b(\zeta)|^2+\|c(\zeta)\|^2=1,\qquad
|b(\zeta)|^2+\|r(\zeta)\|^2=1.                    \tag{10}
\]

This proves (4).  Integrating and applying Hardy Parseval gives

\[
\begin{aligned}
\|c\|_{H^2}^2
&=1-\sum_n|u^*B_nv|^2,\\
\|r\|_{H^2}^2
&=1-\sum_n|u^*B_nv|^2,
\end{aligned}
\]

which proves (3).

Put \(h=1-|b|^2\).  It is nonnegative and has mean \(\delta\).
Therefore

\[
|\widehat h(k)|
\le\frac1{2\pi}\int_{\mathbb T}h
=\delta,
\]

proving (5).  Equation (6) is ordinary \(L^2\)
Cauchy--Schwarz combined with (3).

No approximate-inner or spectral-factor theorem is needed: the
opposite off-channel row has exactly the same mass as the original
off-channel column.

## 3. Why a physical scalar response is quadratic in leakage

Fix a split transfer

\[
B_0=g\oplus B'
\]

and choose frames so its scalar channel is the first coordinate.  In
L218's finite matrix-Schur chart, the block-diagonal inner functions
form the subchart in which every Schur parameter preserves that
coordinate line.

The constant defect-frame changes

\[
J_{\rm out}=1\oplus(-I),\qquad
J_{\rm in}=1\oplus(-I)                             \tag{11}
\]

replace

\[
B\longmapsto J_{\rm out}BJ_{\rm in}
=\begin{bmatrix}b&-r\\-c&D\end{bmatrix}.           \tag{12}
\]

They are only changes of the left and right defect bases.  The
physical contraction \(C\), metric \(P\), operator \(T\), scalar
function norm, and every scalar prepared response are unchanged.
Thus a physical scalar response is even in the joint cross block
\((r,c)\).

The ordered Schur defect roots are analytic.  Projecting every Schur
parameter to its block diagonal gives a local analytic retraction
onto the split-inner subchart.  Its diagonal error is quadratic in
the cross parameters, because

\[
(I-\Gamma\Gamma^*)^{1/2},\qquad
(I-\Gamma^*\Gamma)^{1/2}
\]

have no diagonal term linear in an off-diagonal parameter.
Consequently a response which vanishes on the split subchart has:

1. no term involving only block-diagonal coordinates;
2. no term linear in \((r,c)\), by (12); and
3. every remaining term at least quadratic in the cross blocks.

All norms are equivalent on this fixed finite rational-inner chart.
Equation (3) therefore turns the Taylor estimate into

\[
|{\cal R}(B)|
\le C\{\|r\|_{H^2}^2+\|c\|_{H^2}^2\}
=2C\delta.                                         \tag{13}
\]

Absorb the harmless factor two into \(C\), proving (8).

This argument is local but uniform through channel-frame rotation:
compactness of the two unit spheres gives one constant on a
sufficiently small fixed repeated-block neighbourhood.  Along a
rank-changing analytic arc, choose a maximizing channel frame after
the usual finite ramification.  If the leading cross block vanishes,
L321 descends to the next common-eigenline jet rather than dividing
by a singular frame.

## 4. Normal completion

On L199's zero-Jensen winner, a scalar real/imaginary circular-normal
pair has a positive pure curvature \(\gamma>0\).  After subtracting
L192's strictly controlled split disk/normal face, the new pure
anchor response in upper-gap orientation has the form

\[
-\gamma|y|^2
+2\operatorname{Re}\{\overline y\,{\cal R}(B)\}
+\text{higher terms}.                              \tag{14}
\]

Completing the displayed square costs

\[
\frac{|{\cal R}(B)|^2}{\gamma}
\le\frac{C^2}{\gamma}(1-\sigma)^2.                 \tag{15}
\]

There are only finitely many normal modes and finitely many L197
flag layers for fixed \(L,m\).  Shrink until their summed coefficient
times \(1-\sigma\) is at most \(2/3\).  L322 retains another
\(2/3(1-\sigma)\), while the unspent half of each L199 curvature
absorbs both the \(O(1-\sigma)\) perturbation of L192's strict joint
face and the analytic terms carrying an extra neighbourhood factor.

This closes the response-size obstruction, not yet the bookkeeping
of the complete flag/elliptic assembly.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_scalar_channel_response_factor.py \
  --output \
  experiments/repeated_crabb_scalar_channel_response_factor_s70224.jsonl
```

The checker reconstructs exact \(2\times2\) matrix-inner Schur paths
and audits:

1. both identities in (3) for coordinate and dense channel frames;
2. the positive Fourier-defect bound (5);
3. shifted cross-pairing instances of (6);
4. invariance under the defect-frame sign gauge (12); and
5. exact zero leakage on diagonal split-channel paths.

The checker audits the exact Hardy identities behind (8).  The
analytic factorization itself is the matrix-Schur/gauge argument in
Section 3, not a numerical inference.  The tracked dataset regenerates
byte for byte with SHA-256

```text
beca54b375f93da923b105d4a5a1070d1c51a0da6e9f7bba8cee526f235474e2
```
