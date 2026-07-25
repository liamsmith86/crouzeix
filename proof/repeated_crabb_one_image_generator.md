# Candidate one-image generator for every delayed elliptic face

## 1. Status and candidate result (2026-07-24)

This note records a **candidate all-grade formula**, not yet a proved
lemma.  It incorporates L217's early \(B_5\) reflection and has passed
the complete physical jet through grades one to six, including
noncommuting copy matrices and rank-changing active rows.

Retain

\[
\begin{gathered}
 I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\qquad
 V^*W=0,\\
 P=2I-E+2F,\qquad T=P^{-1/2}SP^{1/2},\\
 B_n=W^*(S^*)^nV,\qquad Q=I-E .
\end{gathered}                                      \tag{1}
\]

Let

\[
 a(c)=\frac1{\vartheta_3(c^2)}
\]

and define the balanced perpendicular zero/one-image frame

\[
\boxed{
\begin{aligned}
d_\perp(c)
=2a(c)\bigg\{&
\sum_{j\ge1}\frac{(-c)^j}{1+c^{4j}}(S^*)^{2j}V\\
&+\sum_{r,j\ge1}
\frac{(-1)^{r+j}c^{2r+j}}{1+c^{4(r+j)}}
QS^rF(S^*)^{r+2j}V
\bigg\}.                                           \tag{2}
\end{aligned}}
\]

The image word has the equivalent ordered form

\[
S^rF(S^*)^{r+2j}V=S^rWB_{r+2j}.                   \tag{3}
\]

Its component parallel to \(V\) is fixed, as in L214--L215, by
requiring the lower Schur complement of the metric minus \(I\) to
vanish coefficientwise.

On the fully delayed stratum

\[
B_1=\cdots=B_{k-1}=0,                              \tag{4}
\]

the experiments support the following identities.

1. Through order \(2k-1\), the balanced metric is the explicit
   boundary-layer series

   \[
   \boxed{
   \widehat P_{\rm bl}(c)
   =I-\sum_{j\ge1}\frac{c^{2j}}{1+c^{2j}}
       (S^*)^jES^j
     +\sum_{j\ge1}c^{2j}S^jF(S^*)^j.}             \tag{5}
   \]

   Equivalently,

   \[
   [c^{2n}]\widehat P_{\rm bl}
   =S^nF(S^*)^n+
     \sum_{d\mid n}(-1)^{n/d}(S^*)^dES^d,          \tag{6}
   \]

   and every odd coefficient is zero.

2. Every lower endpoint through order \(2k\), and every upper
   endpoint below order \(2k\), vanishes.

3. Put

   \[
   L_k=B_kB_k^*,\qquad
   R_k=B_k^*B_k,\qquad
   {\mathfrak C}(K)=\sum_{n\ge1}B_nKB_n^* .
   \]

   The raw upper face selected by (2) is

   \[
   \boxed{
   E_{2k,\rm raw}
   =12L_k-28{\mathfrak C}(R_k)
    -4\,{\bf1}_{2\mid k}\,
       {\cal M}_T(P^{1/2}QS^kWB_k).}               \tag{7}
   \]

4. Hence the parity-only preparation

   \[
   \widehat C_{2k,\rm fold}
   =4\,{\bf1}_{2\mid k}\,QS^kWB_k                 \tag{8}
   \]

   gives the universal base

   \[
   E_{2k,\rm prep}
   =12L_k-28{\mathfrak C}(R_k).                    \tag{9}
   \]

   L212 then gives

   \[
   \boxed{E_{2k,\rm final}=-16B_kB_k^*.}           \tag{10}
   \]

Equations (2), (5), and (7) are the candidate content.  Equations
(8)--(10) follow from them and proved L212.

## 2. Why the generator has this form

For nome \(q=c^2\), the Jacobi Fourier expansion used by the exact
Crabb-axis theorem is

\[
\frac{k'}{\operatorname{dn}(v,k)}
=\frac{\pi}{2K}
\frac{2\pi}{K}\sum_{n\ge1}
\frac{(-1)^nq^n}{1+q^{2n}}
\cos\frac{n\pi v}{K}.                              \tag{11}
\]

After DCT sampling and physical balancing, the positive frequency
\(2j\) becomes

\[
\frac{2(-c)^j}{\vartheta_3(c^2)(1+c^{4j})}
(S^*)^{2j}V,                                      \tag{12}
\]

which is the first line of (2).

If that frequency crosses the terminal defect after \(r\) reverse
steps, its source frequency is \(2(r+j)\).  DCT reflection, physical
balancing, and the ordered terminal word give respectively

\[
(-1)^{r+j},\qquad c^{2r+j},\qquad
S^rF(S^*)^{r+2j}V.
\]

This is the second line of (2).  Its first terms are

\[
\begin{array}{c|l}
\text{frame order}&\text{terminal images}\\ \hline
3&+2SWB_3\\
4&-2SWB_5\\
5&-4SWB_3-2S^2WB_4+2SWB_7.
\end{array}                                        \tag{13}
\]

Thus L217's previously missing \(-2SWB_5\) is not exceptional: it
is the first member of a triangular image series.

At the repeated length-five apex, (2) plus the lower parallel
recursion reproduces L217's exact frame through order ten:

\[
\begin{aligned}
d_5(c)
={}&(1-2c^2+4c^4-8c^6+14c^8-30c^{10})V\\
&+(-2c+4c^3-6c^5+14c^7-26c^9)(S^*)^2V\\
&+(2c^2-6c^4+12c^6-24c^8+42c^{10})(S^*)^4V
+O(c^{11}).
\end{aligned}                                      \tag{14}
\]

This independently explains every coefficient in L217 rather than
fitting only its first obstruction.

## 3. Evidence through grade six

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_one_image_generator.py \
  --output \
  experiments/repeated_crabb_one_image_generator_s70224.jsonl
```

The 30 deterministic records cover:

1. general grade-one unstructured partial isometries;
2. fully delayed unstructured partial isometries at grades two to six;
3. independently gauged heterogeneous shifts whose \(B_k\) changes
   rank;
4. rank-zero \(B_k\) strata; and
5. repeated monomial apices.

Every record independently constructs the direct Riemann map and the
complete lower-tight Stein metric through order twelve.  It checks
(5) below the face, both endpoint hierarchies, the parity formula
(7), the preparation (9), and the final loss (10).  The largest
binary64 residual in the standard run is below \(3\times10^{-10}\).
The tracked dataset SHA-256 is
`88fd0cddb9483e553869dfcd6516bcc3727e130f6259dc24b2bef593056e2167`.

Additional exploratory holdouts included unrelated random inflated
colligations at grades four to six and gave the same identities.  No
coefficient in (2), (5), or (7) was fitted on the tracked cases.

## 4. Exact proof debt

The numerical result is unusually coherent but is not an all-grade
proof.  Two ordered statements remain.

1. **Zero/one-reflection identity.**  In the partial-isometry word
   algebra, prove that (2) and its lower-parallel completion generate
   (5) modulo terms with two terminal crossings.  Under (4), this
   proves every coefficient below \(c^{2k}\).  The scalar identity is
   (11); the work is to retain the state and copy-factor order.

2. **Two-reflection quotient.**  Show that the first surviving
   two-crossing coefficient is L207's grade-one face on the
   L209-deflated colligation.  L216 already transports its linear
   endpoint response.  The sole scalar boundary collision is the
   even-\(k\) central fold (8).  This should prove (7) without a
   fixed-grade expansion.

This route is now preferable to a grade-seven calculation.  A
failure of either ordered statement must identify an additional
two-image word; another blind high-order jet would not resolve that
structural question.
