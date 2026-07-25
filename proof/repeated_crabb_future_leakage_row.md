# The unweighted future channel is one Toeplitz-leakage row

## 1. Result (L259, 2026-07-25)

Let

\[
B(z)=\sum_{j\geq1}B_jz^j
\]

be a square matrix-inner transfer, let \({\cal T}_B\) be
multiplication by \(B\) on \(H^2(\mathbb C^m)\), and put

\[
K_B(\rho,x)
=\frac{I-B(\rho)B^\sharp(x)}{1-\rho x},
\qquad
B^\sharp(x)=\sum_{j\geq1}x^jB_j^*.
\]

Here \(P_j\) denotes the orthogonal projection onto the individual
Hardy row \(z^j\mathbb C^m\) (not the cumulative window projection
used in L253).

L254 identifies \(K_B\) as the reproducing kernel of
\({\cal K}_B\), so its coefficient matrix is
\(P_{{\cal K}_B}=I-{\cal T}_B{\cal T}_B^*\).
Differentiating at the first Hardy row gives

\[
\boxed{
[\rho]K_B(\rho,x)=xI-B_1B^\sharp(x).}             \tag{1}
\]

Equivalently, the apparently unweighted future series in L243 is
exactly one row of the complementary Toeplitz leakage:

\[
\boxed{
P_1{\cal T}_B{\cal T}_B^*P_j
=B_1B_j^*,\qquad j\geq1.}                         \tag{2}
\]

More generally, if

\[
B_1=\cdots=B_{k-1}=0,
\]

then the first nonzero leakage row is

\[
\boxed{
P_k{\cal T}_B{\cal T}_B^*P_j
=B_kB_j^*,\qquad j\geq k,}                        \tag{3}
\]

and it vanishes in columns \(j<k\).  Its diagonal trace is precisely

\[
\boxed{
\operatorname {tr}
(P_k{\cal T}_B{\cal T}_B^*P_k)
=\|B_k\|_F^2.}                                    \tag{4}
\]

Thus L243's future \(B^\sharp\)-coefficients are not independent
endpoint terms.  They are the off-diagonal blocks of the same causal
leakage row whose diagonal is the desired active energy.

L259 does **not** prove A194.  L258's closed-return renewal can still
multiply an off-diagonal block in (3) by an intervening Hardy shift
and return it to the diagonal.  The remaining theorem is now the
sharper diagonal-preservation statement

\[
\operatorname {tr}[c^{2k}]\mathfrak D_{\rm ret}
=4\operatorname {tr}
(P_k{\cal T}_B{\cal T}_B^*P_k),                  \tag{5}
\]

together with the lower vanishings.  In other words, the outstanding
autocorrelation calculation must show that L243--L245 and L258's
physical return kernel take the diagonal of (3), rather than an
uncontrolled shifted sum of its future columns.

L261 gives the more natural equivalent endpoint for that
calculation.  Since \({\cal T}_B{\cal T}_B^*\) is a projection, the
ordinary Hardy norm of the **whole** row in (3) equals its diagonal
in (4).  Thus the physical return may retain all future columns,
provided it pairs them with the unweighted adjoint row; matrix-inner
Parseval then produces (5).  L260 shows that this pairing cannot be
checked separately in the two formal ellipse orientations.

## 2. Kernel proof

Since \(B(0)=0\),

\[
B(\rho)=\rho B_1+O(\rho^2).
\]

Therefore

\[
\begin{aligned}
K_B(\rho,x)
&=\{I-\rho B_1B^\sharp(x)+O(\rho^2)\}
  \{I+\rho x+O(\rho^2)\}\\
&=I+\rho\{xI-B_1B^\sharp(x)\}+O(\rho^2),
\end{aligned}
\]

which proves (1), including the noncommutative order.

On the other hand, L254 gives

\[
P_{{\cal K}_B}=I-{\cal T}_B{\cal T}_B^*.
\]

The coefficient of \(\rho x^j\) in \(K_B\) is the
\((1,j)\) block of this projection.  Equation (1) therefore says

\[
P_1P_{{\cal K}_B}P_j
=\delta_{1j}I-B_1B_j^*.
\]

Subtracting from the corresponding block of the identity proves
(2).

## 3. Delayed row

The lower-triangular Toeplitz matrix of \({\cal T}_B\) has blocks

\[
[{\cal T}_B]_{n,\ell}
=
\begin{cases}
B_{n-\ell},&n>\ell,\\
0,&n\leq\ell.
\end{cases}
\]

Consequently

\[
[{\cal T}_B{\cal T}_B^*]_{k,j}
=\sum_{\ell=0}^{\min(k,j)-1}
B_{k-\ell}B_{j-\ell}^*.                           \tag{6}
\]

Under complete delay, every summand with \(\ell\geq1\) contains
\(B_{k-\ell}=0\).  The \(\ell=0\) term is \(B_kB_j^*\), and it is
also zero for \(j<k\).  This proves (3).  Setting \(j=k\) and taking
the copy-space trace proves (4).

Equation (6) also records why complete delay is essential.  Without
it, a leakage row is a full causal autocorrelation sum, not one
isolated transfer product.

## 4. Exact interface with L243 and L258

L243's inverse delayed model kernel contains
\(B(\rho)B^\sharp(x)\) in exactly the order used above.  At the first
tail coefficient, (1) rewrites it as

\[
B_1B^\sharp(x)=xI-[\rho]K_B(\rho,x).              \tag{7}
\]

The scalar \(xI\) belongs to L244's half-line background.  The
remaining kernel row is the model projection, and its complement is
the leakage row (2).  Hence the later \(B_j^*\) terms that looked
premature in an isolated state column are precisely off-diagonal
leakage blocks.

L258 ensures that these blocks occur only inside a closed return.
L266 now proves the exact endpoint once this return is shown to be
copy-scalar Laurent/Toeplitz:

\[
P_k{\cal L}_B\Psi(L){\cal L}_BP_k
=[\Psi]_0P_k{\cal L}_BP_k.
\]

It also covers copy-scalar shift words whose unilateral boundary depth
is at most \(k\).  Depth \(k+1\) can fail, so Wold copy-scalarity by
itself is insufficient.  Thus no separate cancellation of future
blocks is required once the return's form and depth are proved.  What
is still open is that physical form/depth assertion, the calculation
of its bilateral-symbol constant as \(4\), and all lower vanishings.
Those facts would turn (7) into (5), after which L256/L257 close the
delayed volume law.
