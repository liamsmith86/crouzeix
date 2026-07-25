# Every shifted leakage return misses the first delayed diagonal

## 1. Result (L266, 2026-07-25)

Let

\[
B(z)=\sum_{j\geq1}B_jz^j
\]

be square matrix inner, let \(L\) be the forward shift on
\(H^2(\mathbb C^m)\), and put

\[
{\cal L}_B={\cal T}_B{\cal T}_B^*.
\]

Assume complete delay

\[
B_1=\cdots=B_{k-1}=0.
\]

Let \(P_k\) now denote the projection onto the **individual** Hardy
row \(z^k\mathbb C^m\), as in L259--L261.  For every scalar Laurent
polynomial

\[
\Psi(L)=\sum_{h=-a}^{b}\psi_hL^h,
\qquad
L^{-h}:=(L^*)^h\quad(h>0),
\]

one has

\[
\boxed{
P_k{\cal L}_B\Psi(L){\cal L}_BP_k
=\psi_0P_k{\cal L}_BP_k.}                         \tag{1}
\]

Taking the copy-space trace and using L259 gives

\[
\boxed{
\operatorname {tr}
\{P_k{\cal L}_B\Psi(L){\cal L}_BP_k\}
=\psi_0\|B_k\|_F^2.}                              \tag{2}
\]

The same result holds coefficientwise for a scalar formal Laurent
series whenever the requested matrix coefficient is a finite
stabilized sum.

There is a second form needed for products of half-line operators.
Let \(\Phi(L,L^*)\) be a copy-scalar noncommutative polynomial.  Reduce
each of its words using \(L^*L=I\), so that it has the unique form
\(L^a(L^*)^b\).  Assume

\[
\min(a,b)\leq k                                      \tag{3}
\]

for every reduced word.  If \(\sigma_0(\Phi)\) denotes the constant
coefficient of its **bilateral symbol** (replace \(L\) by \(z\) and
\(L^*\) by \(z^{-1}\)), then

\[
\boxed{
P_k{\cal L}_B\Phi(L,L^*){\cal L}_BP_k
=\sigma_0(\Phi)P_k{\cal L}_BP_k.}                 \tag{4}
\]

Condition (3) is a boundary-depth condition, not cosmetic.  Without
it, unilateral-shift products can contain a finite boundary
projection reaching row \(k\), and (4) can fail.

Thus A194 does not require the intervening physical return to be
literally the identity.  It is enough to prove that, after L264's
right-half-line normalization and L258's complete closure, its first
active Hardy return is copy-scalar Laurent/Toeplitz—or, more
generally, a copy-scalar shift polynomial of boundary depth at most
\(k\)—and that its bilateral-symbol constant is \(4\).  Every
nonconstant shift and every shallow unilateral boundary correction
then misses the first delayed diagonal exactly.

L266 does not prove those physical assertions.  In particular,
copy-scalarity alone does not exclude a deep unilateral boundary
term.  It removes the future-row ambiguity only after the form and
depth of the physical return have been established.

## 2. Shift invariance of the leakage range

Matrix innerness makes \({\cal T}_B\) an isometry, and

\[
\operatorname {ran}{\cal L}_B
=BH^2(\mathbb C^m)
\]

is invariant under \(L\).  Therefore

\[
(I-{\cal L}_B)L{\cal L}_B=0,
\]

or equivalently

\[
\boxed{
{\cal L}_BL^h{\cal L}_B=L^h{\cal L}_B
\qquad(h\geq0).}                                  \tag{5}
\]

Taking adjoints gives

\[
\boxed{
{\cal L}_B(L^*)^h{\cal L}_B
={\cal L}_B(L^*)^h
\qquad(h\geq0).}                                  \tag{6}
\]

No commutation of matrix coefficients is used.

## 3. Delayed triangular support

L259 gives

\[
P_i{\cal L}_BP_k=0\qquad(i<k),                    \tag{7}
\]

and, by self-adjointness,

\[
P_k{\cal L}_BP_i=0\qquad(i<k).                    \tag{8}
\]

For \(h>0\), equations (5) and (7) imply

\[
\begin{aligned}
P_k{\cal L}_BL^h{\cal L}_BP_k
&=P_kL^h{\cal L}_BP_k\\
&=0,
\end{aligned}                                    \tag{9}
\]

because \({\cal L}_BP_k\) has support only in rows at least \(k\),
and \(L^h\) moves that support strictly above row \(k\).

Likewise, equations (6) and (8) give

\[
\begin{aligned}
P_k{\cal L}_B(L^*)^h{\cal L}_BP_k
&=P_k{\cal L}_B(L^*)^hP_k\\
&=0,
\end{aligned}                                    \tag{10}
\]

because \((L^*)^hP_k\) lies in row \(k-h<k\) (or is zero).

For \(h=0\), projection idempotence gives

\[
P_k{\cal L}_B^2P_k
=P_k{\cal L}_BP_k.                                \tag{11}
\]

Linearity of (9)--(11) proves (1).  L259's diagonal identity

\[
P_k{\cal L}_BP_k=B_kB_k^*
\]

then proves (2).

## 4. Low-depth words and the unilateral boundary

Repeatedly cancel every adjacent \(L^*L\) in a word.  The result has
no \(L^*\) to the left of an \(L\), hence is \(L^a(L^*)^b\).  Put

\[
Q_r=P_0+\cdots+P_{r-1}.
\]

The unilateral-shift identities are

\[
L^a(L^*)^b=
\begin{cases}
L^{a-b}(I-Q_b),&a\geq b,\\
(I-Q_a)(L^*)^{b-a},&a<b.
\end{cases}                                      \tag{12}
\]

If \(\min(a,b)\leq k\), then the boundary term in (12) is annihilated
by (7)--(8).  The remaining Laurent shift is annihilated by
(9)--(10) unless \(a=b\), in which case (11) remains.  Applying this
word by word proves (4).

The depth threshold is sharp in this formulation.  For example,

\[
L^{k+1}(L^*)^{k+1}=I-Q_{k+1}
\]

contains \(P_k\), so its boundary correction need not vanish in the
leakage sandwich.

## 5. Interface with the physical calculation

L262 conjugates the normalized right-half-line background to the
ordinary shift with copy-scalar coefficients.  L264 restores that
complete right-half-line metric through the active face.  Therefore a
return built solely from the closed half-line background belongs to
the copy-scalar unilateral shift algebra.  It is **not** automatically
a Laurent/Toeplitz operator: mixed products of \(L\) and \(L^*\) can
carry the boundary terms in (12).

What remains to be proved is that the complete physical numerator at
the first reflected face contains no additional copy-dependent
operator between the two leakage factors and no unilateral boundary
term deep enough to reach row \(k\).  L243 limits this check to the
zero/one inverse-kernel sectors; L251 supplies common paired analytic
ports; and L258 resums every final-row loop before the return is read.
If their combination gives either

\[
P_k{\cal L}_B\Psi_{\rm phys}(L){\cal L}_BP_k,
\]

or the low-depth extension (4), then L266 reduces all future
coefficients to the single scalar
\(\sigma_0(\Psi_{\rm phys})\).  L242's delay-independent
three-coefficient jet and L256's complete grade-one calculation
suggest \(\sigma_0(\Psi_{\rm phys})=4\), but both the boundary-depth
statement and that identification remain live proof obligations.

## 6. Independent finite-inner audit

The proof was also attacked on non-monomial matrix-inner polynomials.
For independent random noncommuting projections \(P,Q\), use the
degree-two Potapov product

\[
C(z)=(I-P+zP)(I-Q+zQ),\qquad B(z)=z^kC(z).
\]

Its first delayed leakage row is generically off diagonal, so this is
not the trivial monomial case.  Direct finite Toeplitz multiplication
for copy sizes \(2,3,4\), delays \(1,2,4\), eight independent pairs
per case, and shifts \(h=\pm1,\ldots,\pm4\) gave

\[
P_k{\cal L}_BL^h{\cal L}_BP_k=0
\]

to below \(2\cdot10^{-12}\).  Every shift word through length seven
whose reduced boundary depth obeyed (3) reproduced (4) below
\(10^{-11}\).  A depth-\(k+1\) word also produced the predicted
nonzero boundary obstruction.  This audit exercises noncommuting
transfer coefficients and nonzero future leakage blocks; the exact
proof in Sections 2--4, not the numerical test, establishes the
result.
