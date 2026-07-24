# The sixth kernel and eighth-order lift (L185--L186, 2026-07-24)

## 1. Results

Retain L183--L184's notation, put \(n=L-1\), and write the
nonconstant Toeplitz direction as

\[
z=(z_0,\ldots,z_{n-1}).
\]

Define the low Wronskian coefficients

\[
S_t=\sum_{i+j=t}(j-i)z_i\overline {z_{n-1-j}},
\qquad 1\le t\le n-2.                              \tag{1}
\]

The cubic Hardy residual has an explicit all-index formula.  In
particular,

\[
\boxed{{\cal R}_L(z)=0\quad\Longleftrightarrow\quad
S_1=\cdots=S_{n-2}=0.}                              \tag{2}
\]

The zero set in (2) has exactly two kinds of points:

1. \(z=\omega J\overline z\), \(|\omega|=1\), the known
   phase-palindromic exact equality cone;
2. \(z\) is supported on one reversal pair
   \(\{j,n-1-j\}\).

This is L185.  Thus the sixth-order kernel found in L184 has no
unclassified component.

On a strict reversal pair, put

\[
z_j=a,\qquad z_{n-1-j}=b,\qquad
d=n-1-2j>0.
\]

The first nonzero Hardy residual is fourth order and has only the two
entries

\[
\boxed{
({\cal R}^{(4)}_L)_{j,n-1-j}
={8d\over L}(|a|^2-|b|^2)^2,
\qquad
({\cal R}^{(4)}_L)_{n-1-j,j}
=-({\cal R}^{(4)}_L)_{j,n-1-j}.}                   \tag{3}
\]

Consequently the canonical eighth-order base deficit is

\[
\boxed{
D_{8,L}
={512d^2\over L^2}(|a|^2-|b|^2)^4.}                \tag{4}
\]

Every true circular-normal response vanishes through quartic order on
this pair.  Therefore (4), with no completed-square subtraction, is
the complete eighth-order Schur face.  It is positive unless
\(|a|=|b|\), precisely the intersection with the phase-palindromic
equality cone.  This is L186.

L185--L186 close the homogeneous kernel fallback.  They do not yet
prove a uniform nonlinear tube through the singular union.

## 2. Closed cubic residual

Adopt the convention \(S_t=0\) outside \(1\le t\le n-2\).  For
\(0\le r<s\le n-1\), direct evaluation of L183 equation (17) gives

\[
\boxed{\begin{aligned}
({\cal R}_L)_{rs}=4\bigg(&
 {z_{n-1-r}S_{n-2-s}\over n-s}
-{z_{n-1-s}S_{n-2-r}\over n-r}\\
&+{\overline z_s\,\overline S_{r-1}\over r+1}
-{\overline z_r\,\overline S_{s-1}\over s+1}
\bigg).
\end{aligned}}                                      \tag{5}
\]

This is the entrywise identity implicit in L184's shorter Hardy
factor argument.  It is obtained by inserting L176's interval pulse
in the order-three inverse recurrence.  On a fixed pulse
anti-diagonal, the two endpoints contribute \(4/(t+2)\); reversing
the Fourier and coefficient indices gives the last two terms.

Formula (5) proves at once that every \(S_t=0\) makes
\({\cal R}_L=0\).  The converse is less formal and is the content of
the next section.

## 3. The residual cannot hide a nonzero Wronskian

Let \(C=J\circ\) conjugation on \(\mathbb C^n\), and define

\[
u_r=z_{n-1-r},\qquad
v_r={S_{n-2-r}\over n-r}.
\]

The out-of-range convention makes \(v_{n-2}=v_{n-1}=0\).  Equation
(5) is exactly

\[
\boxed{{\cal R}_L=4\{u\wedge v-Cu\wedge Cv\}.}       \tag{6}
\]

Suppose \({\cal R}_L=0\).  If \(u\) and \(Cu\) are dependent, then
\(z\) and \(J\overline z\) are dependent, so every \(S_t=0\).

Assume they are independent.  If \(u\wedge v=0\), then
\(v\in\operatorname{span}\{u,Cu\}\).  Otherwise (6) equates two
nonzero decomposable two-forms, so their two-planes coincide and are
\(C\)-invariant.  Again

\[
v=\alpha u+\beta Cu                              \tag{7}
\]

for some scalars \(\alpha,\beta\).

Put

\[
p_i=z_i,\qquad q_i=\overline {z_{n-1-i}}.
\]

Reversing (7) gives

\[
{S_{i-1}\over i+1}=\alpha p_i+\beta q_i,
\qquad 0\le i\le n-1,                              \tag{8}
\]

where \(S_{-1}=S_0=0\).  If \(\alpha=\beta=0\), (8)
already proves the claim.  Otherwise the \(i=0,1\) equations say that
\((p_0,q_0)\) and \((p_1,q_1)\) lie on one line.  Hence

\[
S_1=p_0q_1-p_1q_0=0.
\]

Inductively, suppose \(S_1,\ldots,S_{m-1}=0\).  Equation (8) puts
each \((p_i,q_i)\), \(0\le i\le m\), on the same line.  Pairing the
ordered terms in (1) gives

\[
S_m=\sum_{\substack{i<j\\i+j=m}}
(j-i)(p_iq_j-p_jq_i)=0.
\]

This proves all the equations in (2).

## 4. Classification by a two-critical-point rational map

Introduce

\[
p(x)=\sum_{i=0}^{n-1}z_ix^i,\qquad
q(x)=\sum_{i=0}^{n-1}\overline {z_{n-1-i}}x^i.
\]

The coefficient of \(x^{t-1}\) in

\[
p q'-p'q                                               \tag{9}
\]

is \(S_t\).  The reciprocal-conjugate relation between \(p\) and
\(q\) reflects the coefficients above the middle.  Therefore (2)
is equivalent to

\[
p q'-p'q=cx^{n-2}                                   \tag{10}
\]

for a real scalar \(c\).

If \(c=0\), the two polynomials are proportional.  Applying the
reciprocal-conjugate involution twice shows that the proportionality
constant is unimodular; this is the phase-palindromic cone.

Suppose \(c\ne0\).  Remove the greatest common divisor of \(p,q\).
Equation (10) makes that divisor a monomial.  For the remaining
coprime pair \(p_0,q_0\), the rational map \(p_0/q_0\) has no critical
point except zero and infinity.  Riemann--Hurwitz forces both to be
totally ramified.  Hence

\[
{p_0(x)\over q_0(x)}
={A x^d+B\over C x^d+D}
\]

for some \(d\ge1\).  Thus \(p,q\) are supported on two common
exponents \(m,m+d\).  Their reciprocal-conjugate support has to be the
same pair, so

\[
2m+d=n-1.
\]

The support is exactly one reversal pair.  This completes L185.

## 5. The pair calculation

On one reversal pair, every low coefficient (1) vanishes.  Insert the
resulting terminal diagonal L176 correction into the order-four
version of L183 equation (17).  The residual still vanishes through
order three.  The only surviving endpoint paths run between the two
pair indices, and their signed count is \(d=n-1-2j\).  This gives
(3).

The fourth residual has Frobenius norm

\[
\|{\cal R}^{(4)}_L\|_F^2
=2\left({8d\over L}\right)^2(|a|^2-|b|^2)^4.
\]

The same endpoint-minimization argument as L183 equation (19), now
one order later, gives

\[
D_{8,L}=4\|{\cal R}^{(4)}_L\|_F^2,
\]

which is (4).

Finally, the residual in (3) lies on the central anti-diagonal:

\[
(j+1)+(n-j)=L.
\]

The order-four characteristic/Riemann coefficient comparison pairs a
true normal of mode \(k\ne0\) only with anti-diagonals \(L-k\) and
\(L+k\).  Hence every true circular-normal response is zero through
degree four.  This proves L186.

## 6. Exact regeneration

Run

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_kernel_eighth.py \
  --symbolic-maximum-length 8 \
  --exact-maximum-length 10 \
  --output \
  experiments/crabb_full_disk_kernel_eighth_s70224.jsonl
```

The unrestricted-symbol audit checks (5) and every symbolic pair
formula through length eight.  The independent endpoint and full
characteristic/Riemann audit checks every reversal pair and every true
normal through length ten.  Clean reruns must be byte-identical.

L187--L189 subsequently carry this out in
`proof/crabb_full_disk_leading_residual_tube.md`.  A transverse
weighted correction can cancel the displayed fourth residual, so a
finite eighth-order stopping argument is not uniform.  The full Hardy
residual is instead a submersion onto a curved equality manifold, and
its first eventual nonzero coefficient is absorbed in every order.
