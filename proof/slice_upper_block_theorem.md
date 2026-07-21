# The upper modal block is also bounded by two

This note closes the conformal-displacement debt left in
`slice_boundary_theorems.md`.  It proves, with a finite outward-rounded
algebraic interval certificate for one scalar theta lemma, that the upper
off-diagonal block in the elliptic-slice normal form always has norm at most
two.  Coupled SDP witnesses remain a separate obstruction: the two block-norm
bounds alone do **not** yet prove the complete similarity theorem.

## 1. Statement and notation

Retain the parameters of the boundary note:

\[
 x=\tan u,\qquad r\le p\le1,\qquad
 h=\frac{k}{c},\qquad d_x^2=(1+x^2)(1+r^2x^2),
\]

and

\[
 B=\frac{\sqrt h}{d_x}
 \begin{bmatrix}
  1+prx^2&x(r-p)/c\\
  cx(1-pr)&p+rx^2
 \end{bmatrix}. \tag{1}
\]

It is convenient to put

\[
 \ell=\frac2{\sqrt h}=2\sqrt{\frac c k}. \tag{2}
\]

The nome is \(q=c^2\), and Jacobi's product gives the exact expression

\[
 \ell=\prod_{n\ge1}
 \left(\frac{1+q^{2n-1}}{1+q^{2n}}\right)^2. \tag{3}
\]

In particular, the nome bound from the boundary note is equivalent to
\(\ell\ge1+c^2\).

**Theorem.**  For every nondegenerate elliptic-slice configuration,

\[
 \boxed{\|B\|\le2.} \tag{4}
\]

## 2. The matrix inequality collapses to one scalar barrier

Set \(t=x^2\).  Direct expansion has the same reciprocal-quadratic structure
as for the lower block:

\[
 d_x^2\det(4I-B^TB)=E_0(1+r^2t^2)+E_1t,
 \qquad E_0=(4-h)(4-hp^2)>0. \tag{5}
\]

Since \(1+r^2t^2\ge2rt\), nonnegativity at \(t=1/r\) implies
nonnegativity for every \(t\ge0\).  At that reciprocal midpoint,

\[
 B_*=\frac{\sqrt h}{1+r}
 \begin{bmatrix}
  \sqrt r(1+p)&-(p-r)/c\\
  c(1-pr)&\sqrt r(1+p)
 \end{bmatrix}. \tag{6}
\]

The determinant of \(4I-B_*^TB_*\) factors as

\[
 \det(4I-B_*^TB_*)=
 \frac{16d_-d_+}{c^2\ell^4(1+r)^2}
 (p-p_-)(p_+-p), \tag{7}
\]

where

\[
 d_\pm=\ell(1+c^2r)\pm c(1+r) \tag{8}
\]

and

\[
 p_\pm=\ell\,
 \frac{c^2\pm c\ell(1+r)+r}
      {c^2\ell r\pm c(1+r)+\ell}. \tag{9}
\]

The factor \(d_-\) is positive.  Indeed, use \(\ell\ge1+c^2\) and check the
two endpoints in \(r\) of

\[
 (1+c^2)(1+c^2r)-c(1+r)>0.
\]

The plus root is increasing as a function of \(r\).  If \(c\ell\ge1\), then
\(p_+(r)\ge p_+(0)=c\ell\ge1\ge p\), and there is nothing more to prove.
Suppose henceforth that \(c\ell<1\).  Then

\[
 p_--r=
 -\frac{c(1+r)[\ell^2-r-c\ell(1-r)]}{d_-}\le0, \tag{10}
\]

because the bracket is
\((\ell^2-1)+(1-r)(1-c\ell)\).  Since \(p\ge r\), (7) is therefore
nonnegative exactly when \(p\le p_+(r)\); here one also uses

\[
 p_+-r=\frac{c(1+r)[\ell^2-r+c\ell(1-r)]}{d_+}>0. \tag{10a}
\]

Solving \(p=p_+(r)\) for \(r\) gives the Möbius barrier

\[
 R(p)=\frac{(c+\ell)(p-c\ell)}{(1+c\ell)(\ell-cp)}. \tag{11}
\]

Let

\[
 H(p)=\sin\!\left(sF(\arcsin p\mid m)\right),
 \qquad s=\frac{\pi}{2K(m)},\quad m=k^2. \tag{12}
\]

This is the normalized inverse ellipse map, so \(r=H(p)\).  Thus the entire
upper-block theorem has become

\[
 H(p)\ge R(p),\qquad c\ell\le p\le1. \tag{13}
\]

## 3. A cubic minorant for the inverse ellipse map

The centered inverse ellipse map has nonnegative odd Taylor coefficients
(Kanas--Sugawa).  Expanding (12) through degree three gives

\[
 H(p)=sp+a_3p^3+\sum_{j\ge2}a_{2j+1}p^{2j+1},
 \qquad
 a_3=\frac{s(1+m-s^2)}6,
 \qquad a_{2j+1}\ge0. \tag{14}
\]

Consequently \(H(p)\ge L(p):=sp+a_3p^3\) on \([0,1]\).

Define

\[
 K_0=R'(0)=\frac{(c+\ell)(1-c^2)}{\ell(1+c\ell)}. \tag{15}
\]

The following three scalar inequalities are the only theta-function input:

\[
 \begin{aligned}
  K_0&\ge s, &&\tag{I1}\\
  K_0(1+2c/\ell)&\ge s+3a_3, &&\tag{I2}\\
  s+a_3&\ge R(1). &&\tag{I3}
 \end{aligned}
\]

Assuming them momentarily, use \((1-y)^{-2}\ge1+2y\) to obtain

\[
 R'(p)=\frac{K_0}{(1-cp/\ell)^2}
 \ge K_0(1+2cp/\ell). \tag{16}
\]

The difference between the last expression and
\(L'(p)=s+3a_3p^2\) is a concave quadratic in \(p\).  Its endpoint values are
nonnegative by (I1) and (I2), so \(R'(p)\ge L'(p)\) throughout \([0,1]\).
Hence \(L-R\) is decreasing, and (I3) gives

\[
 H(p)\ge L(p)\ge R(p). \tag{17}
\]

## 4. Finite certificate for the three scalar inequalities

For completeness, this section records an independently reproducible finite
certificate rather than treating (I1)--(I3) as numerical plots.  The checker is
`experiments/slice_upper_block_certificate.py`.

For \(q=c^2\), set

\[
 \theta_N=1+2\sum_{n=1}^Nq^{n^2},\qquad
 P_N=\prod_{n=1}^N
 \left(\frac{1+q^{2n-1}}{1+q^{2n}}\right)^2. \tag{18}
\]

Two geometric-tail estimates give rigorous algebraic intervals:

\[
 \theta_N\le\theta_3(q)\le
 \theta_N+\frac{2q^{(N+1)^2}}{1-q^{2N+3}}, \tag{19}
\]

and, with \(X_N=2q^{2N+1}/(1-q^2)\),

\[
 P_N\le\ell\le P_Ne^{X_N}\le\frac{P_N}{1-X_N}. \tag{20}
\]

Here (20) follows by applying \(\log(1+x)\le x\) to the omitted product and
then \(e^x\le(1-x)^{-1}\).  Also

\[
 s=\theta_3(q)^{-2},\qquad m=\left(\frac{4c}{\ell^2}\right)^2, \tag{21}
\]

so (18)--(21) enclose every quantity in (I1)--(I3) using algebraic interval
operations only.

The certificate has three pieces.

The substitutions used below preserve the required inequality directions:
\(K_0\), \(K_0(1+2c/\ell)\), and \(R(1)\) all decrease with \(\ell\), while
\(s(3+m-s^2)/2\) and \(s(7+m-s^2)/6\) increase with both \(s\in(0,1]\)
and \(m\in(0,1)\).  Differentiating these elementary rational expressions
verifies the assertions directly.

1. On \(0<c\le1/50\), one product factor and explicit theta tails reduce
   (I1)--(I3) to rational functions.  After removing their respective zeros
   \(c^2,c,c^3\), exact coefficient domination gives normalized numerator
   lower bounds \(2.9799759964\), \(1.859816\), and \(22.0711486195\);
   the corresponding denominator magnitudes are at least \(1,1,3\).
   SymPy performs these checks over exact rationals.
2. On \([1/50,13/20]\), (19)--(21) with \(N=12\) and 50-decimal
   outward-rounded intervals certify 13,500 boxes.  The smallest returned
   lower endpoints are

   \[
   1.1557181\times10^{-3},\qquad
   3.6322322\times10^{-2},\qquad
   3.6504258\times10^{-5}. \tag{22}
   \]

3. For \(13/20\le c<1\), already \(cP_2\ge1\).  After clearing the positive
   denominator, \(cP_2-1=(1-c)Q(c)/(\text{positive})\), where \(Q\) has
   nonnegative derivative coefficients and \(Q(13/20)>0\).  Thus
   \(c\ell\ge cP_2\ge1\), which is the trivial root regime handled above.

The checker also reconstructs (7) symbolically, so the matrix factorization
and the scalar certificate are audited in one run.  Its current output is:

```text
determinant factorization: exact
interval boxes certified: 13500
  I1: lower endpoint >= 0.00115571812773
  I2: lower endpoint >= 0.036322322537
  I3: lower endpoint >= 3.65042584897e-05
```

This proves (I1)--(I3), hence (17) and the nonnegativity in (7).  Finally,
\(|\det B|=hp<4\), because the sharp nome bound gives
\(h\le4/(1+c^2)^2<4\).  Thus the two singular values of \(B\) cannot both be
larger than two.  Nonnegativity of \(\det(4I-B^TB)\) therefore puts both of
them below two and proves (4).

## 5. What remains

Together with the boundary note, both modal blocks now obey

\[
 \boxed{\|B\|\le2,\qquad\|C\|\le2.} \tag{23}
\]

This closes every one-block dual phase.  It does not control dual witnesses
with both parity blocks nonzero: the SDP optimum is numerically strictly larger
than \(\max(\|B\|^2,\|C\|^2)\) in several interior configurations.  The live
problem is therefore the genuinely coupled trace inequality in (8) of
`slice_similarity_duality.md`.
