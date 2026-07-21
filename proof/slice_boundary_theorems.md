# Boundary and one-block theorems for the elliptic-slice similarity problem

This note proves two pieces of `L20'` from `slice_similarity_duality.md`: a sharp nome bound
needed to control the modal blocks, and the complete similarity theorem on the singular
(c\to0) face.  It also proves that one of the two off-diagonal blocks of
(T=\phi(A)) always has norm at most two.  These results do not yet control coupled dual
witnesses in the interior of the three-parameter slice.

## 1. A sharp elementary bound for the focus node

Let (q=c^2\in(0,1)), let (k\in(0,1)) be the elliptic modulus with nome (q), and put

\[
 h=\frac{k}{c}.
\]

**Lemma (nome bound).**

\[
 \boxed{h\le \frac{4}{(1+c^2)^2}.} \tag{1}
\]

**Proof.**  Jacobi's product for the elliptic modulus is

\[
 k=4\sqrt q\,\prod_{n\ge1}
 \left(\frac{1+q^{2n}}{1+q^{2n-1}}\right)^4. \tag{2}
\]

The elementary inequality

\[
 (1+q^{2n})^2\le(1+q^{2n-1})(1+q^{2n+1}) \tag{3}
\]

has difference (q^{2n-1}(1-q)^2\ge0).  Multiplying (3) for (1\le n\le N) and cancelling
the repeated odd factors gives

\[
 \left[\prod_{n=1}^N\frac{1+q^{2n}}{1+q^{2n-1}}\right]^2
 \le\frac{1+q^{2N+1}}{1+q}.
\]

Letting (N\to\infty) in (2) yields (k\le4\sqrt q/(1+q)^2).  Since (sqrt q=c), this is
exactly (1).  (□)

The constant is asymptotically sharp as (c\downarrow0): (h\to4).

## 2. Explicit entries of the modal blocks

Use the exact parameters of `slice_similarity_duality.md`.  Write

\[
 x=\tan u,\qquad p=\frac{\tau_2}{\tau_1},\qquad
 d_x=\sqrt{(1+x^2)(1+r^2x^2)}.
\]

The normalized singular-value transform is

\[
 \frac1{\tau_1}U\Sigma V^T
 =\frac1{d_x}
 \begin{bmatrix}
 1+prx^2&x(r-p)\\
 x(1-pr)&p+rx^2
 \end{bmatrix}. \tag{4}
\]

The inverse ellipse map is odd with nonnegative Taylor coefficients (Kanas--Sugawa), or more
elementarily is convex on the positive radius by the Schwarz--Jack theorem.  Its normalized
secant relation therefore gives

\[
 0<r\le p\le1. \tag{5}
\]

Substitution in (12) of the modal note gives

\[
 B=\frac{\sqrt h}{d_x}
 \begin{bmatrix}
 1+prx^2&x(r-p)/c\\
 cx(1-pr)&p+rx^2
 \end{bmatrix}, \tag{6}
\]

\[
 C=\frac{c\sqrt h}{d_x}
 \begin{bmatrix}
 1+prx^2&x(1-pr)/c\\
 cx(r-p)&p+rx^2
 \end{bmatrix}. \tag{7}
\]

## 3. The lower modal block always has norm at most two

**Theorem.**  For every nondegenerate elliptic-slice configuration,

\[
 \boxed{\|C\|\le2.} \tag{8}
\]

**Proof.**  For a (2\times2) matrix (M), the condition (|M|\le2) is equivalent to
(det(4I-M^*M)\ge0), provided both singular values cannot exceed two.  Applied to (7), direct
expansion gives

\[
 d_x^2\det(4I-C^TC)=E_0(1+r^2x^4)+E_1x^2, \tag{9}
\]

where (E_0=(4-c^2h)(4-c^2hp^2)>0).  Hence, by
(1+r^2x^4\ge2rx^2), it suffices to check (9) at (x^2=1/r).  At that point

\[
 C_*=\frac{\sqrt h}{1+r}
 \begin{bmatrix}
 c\sqrt r(1+p)&1-pr\\
 -c^2(p-r)&c\sqrt r(1+p)
 \end{bmatrix}. \tag{10}
\]

Because of (5),

\[
 c^2(p-r)\le1-r\le1-pr. \tag{11}
\]

The norm of the off-diagonal part of (10) is therefore (1-pr), and the triangle inequality
reduces the proof to

\[
 \frac{\sqrt h}{1+r}\,[c\sqrt r(1+p)+1-pr]\le2. \tag{12}
\]

The bracket is affine in (p), with coefficient (sqrt r(c-\sqrt r)).

* If (r\ge c^2), its maximum on (r\le p\le1) occurs at (p=r).  After division by (1+r)
  it is (c\sqrt r+1-r\le1).  The weaker consequence (h\le4) of (1) proves (12).
* If (r\le c^2), its maximum occurs at (p=1).  With (y=\sqrt r), the resulting quotient is

  \[
  \frac{1+2cy-y^2}{1+y^2}\le\sqrt{1+c^2}, \tag{13}
  \]

  the largest eigenvalue of the real symmetric matrix
  (egin{bmatrix}1&c\\c&-1end{bmatrix}).  Equation (1) again proves (12).

Thus (9) is nonnegative for every (x\ge0).  Finally,
(|\det C|=ckp<1), so both singular values cannot simultaneously exceed two.  Therefore
(|C|\le2).  (□)

The transpose-position block (B) also satisfies (|B|\le2) on every numerical test, but the
same proof leaves one genuine conformal distortion term, ((p-r)/c).  It requires more than the
soft bounds (5), so it remains explicitly unclaimed here.

## 4. The complete similarity bound on the (c\to0) face

Normalize the bidiagonal block as in the modal note and write

\[
 R=\begin{bmatrix}a&0\\b&dend{bmatrix},\qquad \|R\|=1,qquad a,b,d\ge0.
\]

The small-nome limits (	au_1\sim2\sqrt c), (	au_2\sim2r\sqrt c) in (12) give

\[
 B_0=2\begin{bmatrix}a&0\\0&dend{bmatrix},\qquad
 C_0=2\begin{bmatrix}0&b\\0&0end{bmatrix}. \tag{14}
\]

After a permutation, (T_0=egin{bmatrix}0&B_0\\C_0&0end{bmatrix}) is the nilpotent scalar
weighted shift with successive weights (2a,2b,2d).

Since the top singular value of (R) is one,

\[
 \det(I-R^TR)=0
 \quad\Longrightarrow\quad
 b^2=(1-a^2)(1-d^2). \tag{15}
\]

Every consecutive product of the three shift weights is at most two:

\[
 2a,2b,2d\le2,\qquad
 4ab,4bd\le2,
\]

because (a\sqrt{1-a^2}\le1/2) and (d\sqrt{1-d^2}\le1/2), while

\[
 8abd=2[,2a\sqrt{1-a^2},][,2d\sqrt{1-d^2},]\le2. \tag{16}
\]

For completeness, if (W) is a nilpotent scalar weighted shift and (M) is the largest
consecutive weight product (including the empty product (1)), set

\[
 x_j=\max\{1,w_j,w_jw_{j+1},\ldots,w_j\cdots w_{n-1}\}.
\]

Then (x_j\ge w_jx_{j+1}), so the diagonal similarity (X^{-1}WX) is a contraction, while
(1\le x_j\le M).  Equations (15)-(16) give (M\le2).  Consequently

\[
 \boxed{t_*(T_0)\le4.} \tag{17}
\]

Equality includes (a=d=1/\sqrt2,b=1/2), which produces the Crabb weights
((\sqrt2,1,\sqrt2)).  Thus the constant four is already sharp on this boundary, explaining the
near-four values of the modal SDP as (c\downarrow0).

## 5. Remaining obstruction

Equations (1), (8), and (17) remove the lower-block and singular-boundary obstructions.  The
interior proof still has to control either the upper block's conformal displacement
((p-r)/c) or genuinely coupled dual witnesses with both (Z_o,Z_e\ne0).  Numerically, the
latter are exactly the phases where (t_*(T)) is larger than (max(\|B\|^2,\|C\|^2)).
