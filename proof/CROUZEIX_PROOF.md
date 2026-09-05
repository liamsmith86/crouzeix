# Crouzeix's inequality for every finite matrix and scalar polynomial

**Attribution.** This is a self-contained finite-dimensional reconstruction of
the all-powers proof of Emiel Lorist and Felix Schwenninger,
[A solution to Crouzeix's conjecture, v2, 17 August 2026](https://arxiv.org/html/2608.03841v2).
The decisive dilation lemma is theirs, not a new campaign discovery. Here the
application uses polygonal contours and a direct resolvent bound, avoiding the
boundary regularity and bounded Cauchy-transform machinery of their more
general operator theorem. Publication status and other recent results are in
`literature_refresh_20260905.md`.

## Theorem

For every positive integer d, every A in C^(d×d), and every complex polynomial p,

\[
 \|p(A)\|_2\le 2\max_{z\in W(A)}|p(z)|,
 \qquad W(A)=\{x^*Ax:x^*x=1\}.
\]

All norms below are Hilbert-space operator norms unless applied to vectors.
Inner products are linear in the first argument. No dimension, polynomial-degree,
diagonalizability, nondegeneracy, or extremality restriction is imposed.

## 1. The all-powers lemma

Let T act on a finite-dimensional complex Hilbert space H. Suppose V:H→K is an
isometry into a Hilbert space K, Q:K→K is a contraction, and, for every integer j≥1,

\[
 E_j=2V^*Q^{*j}V-T^{*j},\qquad [E_j,T]=0,
 \qquad M:=\sup_{j\ge1}\|E_j\|<\infty.                 \tag{1}
\]

Then ||T||≤2.

**Proof.** Put k=||T||. The case k≤1 needs no argument. For k>1 choose a unit
right singular vector x with T*Tx=k²x. From (1), ||T^j||≤2+M, so

\[
 m_j=\operatorname{Re}(x^*E_jT^jx)
 \quad\hbox{satisfies}\quad |m_j|\le M(2+M).          \tag{2}
\]

Set u=Q*VTx−kVx, D=||u||², w_j=T*^j x, and b_j=V*Q*^j u.
Commutativity in (1), followed by T*Tx=k²x, gives the exact identity

\[
\begin{aligned}
 km_j-m_{j+1}
 &=\operatorname{Re}\{x^*T^j(kE_j-E_{j+1}T)x\}\\
 &=k(k-1)\|w_j\|^2-2\operatorname{Re}(w_j^*b_j)\\
 &=k(k-1)\left\|w_j-\frac{b_j}{k(k-1)}\right\|^2
       -\frac{\|b_j\|^2}{k(k-1)}
 \ge -\frac{D}{k(k-1)}.                             \tag{3}
\end{aligned}
\]

For clarity, the only noncommutative cancellation in the second line is
T*^(j+1)Tx=k²T*^j x. The dilation part of (kE_j−E_(j+1)T)x is −2b_j.
The last inequality uses ||V*Q*^j||≤1.

Multiply (3) by k^(−j) and sum j=1,…,N. All intermediate m_j cancel:

\[
 m_1-k^{-N}m_{N+1}
 \ge -\frac{D}{k(k-1)^2}(1-k^{-N}).                 \tag{4}
\]

The boundedness in (2), not an assumed decay of T^j, permits N→∞. Thus
k m_1≥−D/(k−1)². Independently, expansion of u and contractivity give

\[
 D\le 2k^2-2k\operatorname{Re}(x^*V^*Q^*VTx)
   =2k^2-k^3-km_1.                                 \tag{5}
\]

Here 2V*Q*V=E_1+T* and x*T*Tx=k². Combining (4)–(5) yields

\[
 D\left(1-\frac1{(k-1)^2}\right)\le k^2(2-k).      \tag{6}
\]

If k>2 the left side is nonnegative and the right side strictly negative.
This contradiction proves the lemma. The bound M is needed only to control
the vanishing tail in (4). Although D also depends on the operators, the
contradiction works for every D≥0, so the conclusion is dimension-free. ∎

## 2. Elementary numerical-range and contour facts

W(A) is compact, by continuity on the finite-dimensional unit sphere. Every
eigenvalue belongs to W(A), by evaluating a unit eigenvector.

Here is also a finite-dimensional proof of the convexity needed below. Given
two unit vectors, compress A to their span, which has dimension at most two.
In dimension two, write a unit vector as (a,b). Its expectation is an affine
real-linear function of
(2 Re(conj(a)b), 2 Im(conj(a)b), |a|²−|b|²), which ranges over the unit sphere
in R³. A real-linear map R³→R² has a nonzero kernel, so its image of the unit
sphere equals its image of the unit ball: add a kernel vector to a minimum-norm
preimage to reach the sphere. That image is convex. The numerical range of
the compression contains the two chosen expectations and their joining segment,
and is contained in W(A). Thus W(A) is convex, including the point/segment cases.

For any nonempty compact convex K⊂C and ε>0, there is a bounded convex polygon P with
K⊂int(P) and P⊂{z:dist(z,K)<ε}. One explicit existence argument is as follows.
Choose a closed square B with K in its interior. For every z in the compact set
B\{dist(·,K)<ε}, let y be its nearest point in K and put
ν=(z−y)/|z−y|. The nearest-point inequality gives
Re(conj(ν)(w−y))≤0 for all w∈K. The closed half-plane
Re(conj(ν)(w−y))≤ε/2 contains K strictly and excludes an open neighborhood of z.
A finite subcover of those excluded neighborhoods exists. Intersect the
corresponding finitely many half-planes with B to obtain P. If the compact set
is empty, take P=B. In either case K lies in the interior and P lies in the
ε-neighborhood. This construction needs no smoothness or positive area of K.

Let Γ=∂P be counterclockwise oriented. On each open edge let ν be its outward
unit normal and ds its arc length, so dσ=iν ds. Vertices form a measure-zero set
and never enter the pointwise kernel statements. With R(σ)=(σI−A)^(−1), the
polynomial contour identity is

\[
 \frac1{2\pi i}\int_\Gamma h(\sigma)R(\sigma)\,d\sigma=h(A)
 \quad\hbox{for every polynomial }h.                \tag{7}
\]

Indeed all poles of R are eigenvalues inside Γ. Deform to a large circle
entrywise and expand R(σ)=Σ_(r≥0) A^r σ^(−r−1) there; the uniformly convergent
series and scalar Cauchy coefficient formula give (7). This proves the needed
functional calculus directly even for defective A.

## 3. A positive contour density

Fix such a polygon containing W(A) strictly and define almost everywhere on Γ

\[
 P_A(\sigma)=\frac1{2\pi}
       \{\nu R(\sigma)+\overline\nu R(\sigma)^*\}.  \tag{8}
\]

This is a bounded positive matrix density. To check the sign, put
B=σI−A. Direct multiplication gives

\[
 \operatorname{Re}(\nu B^{-1})
 =B^{-*}\operatorname{Re}(\overline\nu B)B^{-1}.    \tag{9}
\]

For every unit v, the middle quadratic form equals
Re(conj(ν)(σ−v*Av))>0 by the supporting line of that edge and W(A)⊂int(P).
This proves positivity without any scalar-to-operator norm inference.
Equation (7) with h=1 and its adjoint yield the exact mass

\[
 \int_\Gamma P_A(\sigma)\,ds=2I.                  \tag{10}
\]

## 4. Apply the lemma using direct commuting contour errors

Let f be any polynomial with max_P |f|≤1. Take
K=L²(Γ,ds;C^d), Vx(σ)=(P_A(σ)/2)^(1/2)x, and Qg(σ)=f(σ)g(σ).
The measurable positive square root exists in finite dimensions. Equation (10)
makes V an isometry; the boundary bound makes Q a contraction. Set T=f(A).
For every j≥1, multiplication and (7) show

\[
\begin{aligned}
 2V^*Q^{*j}V-T^{*j}
 &=\int_\Gamma\overline{f(\sigma)}^{\,j}P_A(\sigma)\,ds
       -(f^j)(A)^*\\
 &=\frac1{2\pi i}\int_\Gamma
       \overline{f(\sigma)}^{\,j}R(\sigma)\,d\sigma
 =:E_j.                                             \tag{11}
\end{aligned}
\]

The equality T^j=(f^j)(A) is ordinary polynomial multiplication; no degree
reduction is made. The conjugated boundary values in (11) need not be analytic.
Each R(σ) commutes with A and therefore with f(A), so E_j commutes with T.
Also, uniformly in j,

\[
 \|E_j\|\le\frac{\operatorname{length}(\Gamma)}{2\pi}
                     \max_{\sigma\in\Gamma}\|R(\sigma)\|<\infty. \tag{12}
\]

Compactness of Γ and absence of eigenvalues on it justify finiteness. We do not
use Crouzeix's conjecture, any earlier spectral-set bound, or a boundedness
theorem for a scalar Cauchy transform to obtain (12). The all-powers lemma now
proves ||f(A)||≤2. Normalizing any nonzero polynomial p by max_P |p| gives

\[
 \|p(A)\|\le2\max_{z\in P}|p(z)|.                 \tag{13}
\]

The denominator is positive because P has interior and a nonzero polynomial
cannot vanish there. The zero polynomial is immediate.

## 5. Remove the exterior polygon and check sharpness

Choose the polygons above with ε→0. Uniform continuity of p on a fixed compact
neighborhood of W(A) gives max_P |p|→max_(W(A)) |p|. Apply this to (13). The
contour bound (12) may diverge as ε→0; this does not matter: the lemma is applied
separately for each fixed polygon and its conclusion is always exactly two.
No limit is exchanged with the infinite power sum or a matrix-dependent bound.

This proves the stated theorem even when W(A) is a point or a segment, and even
when max_(W(A)) |p|=0. In the latter case (13) tends to zero without division
by that limiting supremum. Constant polynomials and dimension one are included.

For sharpness take A=[[0,2],[0,0]] and p(z)=z. Its norm is two. For a unit vector
(a,b), x*Ax=2conj(a)b, whose modulus is at most one and attains one when
|a|=|b|=1/√2. Thus max_(W(A)) |z|=1 and equality holds. ∎

## Verification boundary

The proof uses only finite-dimensional spectral theory, elementary compactness
and convex separation as proved above, scalar Cauchy's theorem, and basic L²
operator constructions. No lemma from the historical local Hessian program is
load-bearing. Tests and independent audits are recorded separately; numerical
checks are not substituted for any step here.

This is the scalar theorem required by `goal.txt`. Matrix-valued polynomials
do not automatically satisfy the commutation in (11), so no general complete
2-spectral-set conclusion is claimed.
