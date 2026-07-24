# The all-size full-disk cubic response (L182, 2026-07-24)

## 1. Result

Retain L176--L181's notation.  Put \(p=L+1\), \(n=L-1\), and
write the nonconstant Toeplitz direction as

\[
 z=(z_0,\ldots,z_{n-1}),\qquad
 W=z\wedge J\overline z.
\]

For every intrinsic nonterminal anti-diagonal, define

\[
 S_t=\sum_{0\le i<t-i}(t-2i)W_{i,t-i},
 \qquad 1\le t\le n-2.                              \tag{1}
\]

Along L176's exactly recentered disk path

\[
 H(s)=I/2+sZ(z)+s^2B(W),
\]

the complete complex cubic response in circular-normal mode \(k\) is
zero unless \(3\le k\le L-3\).  In the active range it is

\[
\boxed{\begin{aligned}
G^{(3)}_{L,k}
={16(4k-1)\over L^2}\bigg\{&
\sum_{\substack{a+t=L+k-3\\k\le a\le L-2}}
 {a-t-1\over t+2}\,z_aS_t\\
&+\sum_{\substack{a+t=L-k-3\\a\ge0,\ t\ge1}}
 {a-t-1\over t+2}\,
 \overline{z_aS_t}\bigg\}.
\end{aligned}}                                      \tag{2}
\]

Thus the apparently dense cubic tensors in L180--L181 are one
triangular interval-flux transform in every size.  This proves the
response half of the arbitrary-size sixth-order Schur face.  Positivity
of the base deficit minus the completed response gain is still open.

## 2. Ordered anti-diagonal identity

The intrinsic definition (1) has the equivalent ordered form

\[
\boxed{
S_t=\sum_{\substack{i+j=t\\0\le i,j<n}}
(j-i)z_i\overline{z_{n-1-j}}.}                      \tag{3}
\]

Indeed, expand \(W_{ij}\) in (1) and interchange \(i,j\) in its
second term.  Formula (3) is useful because every summand has one
oriented endpoint path; no exterior-algebra relation is needed in the
calculation below.

## 3. Cubic path recurrence

Use the coefficient gauge and the same series as L177:

\[
\begin{aligned}
K^{-1}&=N_0+sN_1+s^2N_2+s^3N_3+O(s^4),\\
A&=A_0+sA_1+s^2A_2+s^3A_3+O(s^4).
\end{aligned}                                       \tag{4}
\]

Run the paired Horner recurrences for the characteristic polynomial,
its reversal, and their Fréchet derivatives after subtracting the
first inverse-Riemann correction.  Fix one intrinsic coordinate
\(W_{i,t-i}\), and let \(r^\pm_i(a,t)\) be its normalized coefficient
in the two possible degree-three characters

\[
z_aW_{i,t-i},\qquad
\overline{z_aW_{i,t-i}}.
\]

Circle character gives, respectively,

\[
a+t=n+k-2,\qquad a+t=n-k-2.                         \tag{5}
\]

Substitution of L176's pulse

\[
b_i^{(t)}(u)
=2{\bf1}_{\{i+1,\ldots,t-i\}}(u)
 -{2(t-2i)\over t+2}                                \tag{6}
\]

in the degree-three recurrence leaves the following scalar
first-difference system:

\[
\begin{aligned}
r^\pm_0(a,t)&={t(a-t-1)\over t+2},\\
r^\pm_{i+1}(a,t)-r^\pm_i(a,t)
&=-{2(a-t-1)\over t+2}.                             \tag{7}
\end{aligned}
\]

Here the interval endpoints in (6) give the difference equation, and
its removed mean gives the denominator \(t+2\).  The two endpoint
weights of the Crabb metric are already included in the initial value.
Every term not satisfying (5) has a different circle character and
vanishes in the \(k\)-th polarization.

Solving (7) gives

\[
\boxed{
r^\pm_i(a,t)={(a-t-1)(t-2i)\over t+2}.}             \tag{8}
\]

This is the only index calculation required.  Summing (8) against
\(W_{i,t-i}\) produces \(S_t\).

For completeness, the fixed complex Riesz representative of the
\(k\)-th support polarization has entries

\[
(N_k)_{j+k-1,j}=
\begin{cases}
1/L,&0\le j<L-k+1,\\
1/(2L),&j=L-k+1,
\end{cases}                                         \tag{9}
\]

with all other entries zero.  The endpoint/defect pairing in the
Horner recurrence contributes the remaining common factor

\[
{16(4k-1)\over L^2}.                                \tag{10}
\]

Equations (5), (8), and (10) are exactly (2).  In the plus
orientation the intrinsic range is \(a=k,\ldots,n-1\), hence
\(t=n-2,\ldots,k-1\).  In the conjugate orientation it is
\(a=0,\ldots,n-k-3\), hence \(t=n-k-2,\ldots,1\).
If \(k>L-3\), neither range contains a valid intrinsic
anti-diagonal, proving inactivity.

## 4. Checks against the earlier formulas

For the highest mode \(k=L-3\), (2) has only two terms:

\[
G^{(3)}_{L,k}
={16(4k-1)\over L^2}
\left(-{z_kS_k\over k+2}
      +{z_{k+1}S_{k-1}\over k+1}\right),            \tag{11}
\]

which is L178 equation (6).

For \(L=7,k=3\), (2) becomes

\[
{176\over147}
\left(3z_5W_{02}-4z_3W_{04}-2z_3W_{13}
      -2\overline z_0\,\overline W_{01}\right),
\]

exactly L181's independently polarized formula.

For \(L=8\), the three normalized coefficient rows are

\[
\begin{array}{c|rrrr|rr}
k&z_kS_{L-3}&z_{k+1}S_{L-4}&z_{k+2}S_{L-5}
&z_{k+3}S_{L-6}&\overline z_0\overline S_{L-k-3}
&\overline z_1\overline S_{L-k-4}\\ \hline
3&-6/7&-1/3&2/5&3/2&-3/2&-2/3\\
4&-4/7&0&4/5&0&-4/3&0\\
5&-2/7&1/3&0&0&0&0
\end{array}                                         \tag{12}
\]

after dividing by \(8(4k-1)/L^2\).  These are the unique full-rank
flux coefficients reconstructed from L181's 364-point tensor.

## 5. Exact regeneration

Run

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_cubic_response_exact.py \
  --minimum-length 6 --maximum-length 11 \
  --direction-count 2 --workers 8 \
  --output \
  experiments/crabb_full_disk_cubic_response_exact_s70224.jsonl
```

The checker does not use (7)'s recurrence implementation.  It runs the
complete characteristic/reversed-Horner/inverse-Riemann engine and
compares both real polarizations with (2), including every inactive
mode.  It uses exact SymPy arithmetic.

As a separate discovery audit, fourteen unrelated exact complex
directions in length nine give full-rank flux systems of sizes
\(11,10,9,8\) for modes \(3,4,5,6\).  Solving those systems without
supplying (2) reproduces every coefficient exactly.  This recovery is
not used in the proof.

## 6. Subsequent closure

L182 reduced the remaining arbitrary-size debt to

\[
D_{6,L}(z)
-\sum_{k=3}^{L-3}{|G^{(3)}_{L,k}(z)|^2\over4b_{L,k}}
\ge0,                                                \tag{13}
\]

with \(G^{(3)}_{L,k}\) explicitly given by (2).

L183--L184 subsequently close this inequality by a shorter route.
The disk-chart Hardy residual has a cubic skew coefficient
\({\cal R}_L\), the base deficit is \(8\|{\cal R}_L\|^2\), and the
two grade-\(\pm k\) response rows are disjoint weighted
anti-diagonal projections of \({\cal R}_L\).  Cauchy--Schwarz gives
exactly the flux-curvature gain, while L173's null lift supplies the
remaining strict margin.  See
`proof/crabb_full_disk_sixth_hardy_factor.md`.
