# The arbitrary-size elliptic Crabb axis (2026-07-22)

## 1. Status

For the `p x p` Crabb block `C_p`, put

\[
 A_c=C_p+cC_p^*,\qquad 0<c<1,\qquad L=p-1.
\]

The optimal L21 similarity square for the disk pullback is conjecturally

\[
 \boxed{t_*(\phi_c(A_c))={k(c^{2L})\over c^L}<4,}       \tag{1}
\]

where `k(q)=(theta_2(q)/theta_3(q))^2` is Jacobi's elliptic
modulus as a function of the nome.

This note proves the lower bound in (1), proves the strict scalar inequality,
and gives explicit matching upper metrics for `p=3,4`.  The
matching diagonal upper certificate is numerically confirmed for every
`p=3,...,10` on four ellipse parameters.

The all-size upper bound is **not yet proved**.  Kenan Li's 2021 thesis
proposes exactly this certificate, but page 46 says that the identities
needed for its contraction property were proved only for `p=2,...,6` and
tested numerically beyond.  No later closure was found in the targeted
literature audit.  This is therefore a sharply isolated proof gap rather
than an imported all-dimensional theorem.

## 2. Exact elliptic functional calculus

Let

\[
 D=\operatorname{diag}(1,c^{1/2},\ldots,c^{L/2}),\qquad
 H_c=\sqrt c\,(C_p+C_p^*).
\]

Then

\[
 A_c=D H_cD^{-1},\qquad
 \sigma(H_c)=\{2\sqrt c\cos(j\pi/L):0\le j\le L\}.       \tag{2}
\]

The numerical range of `A_c` is the ellipse with semiaxes `1+c` and
`1-c` and foci `+-2sqrt(c)`.  If

\[
 q=c^2,\quad k=k(q),\quad K=K(k^2),
\]

its centered Riemann map has nodal values

\[
 \tau_j=\sqrt{k}\,
   \operatorname{sn}\!\left(K(1-2j/L)\mid k^2\right).   \tag{3}
\]

Writing `U` for the orthogonal DCT-I eigenvector matrix of `H_c` and
`F=U diag(tau_j)U^T`, functional calculus gives

\[
 T_c:=\phi_c(A_c)=D F D^{-1}.                           \tag{4}
\]

This supplies a dimension-independent exact construction, with no boundary
quadrature or fitted conformal map.

## 3. Exact Chebyshev--Blaschke lower bound

The normalized degree-`L` Chebyshev--Blaschke product for the interval
`[-sqrt(k(q)),sqrt(k(q))]` maps the elliptic Lobatto points (3) to

\[
 B_L(\tau_j)=(-1)^j\sqrt{k(q^L)}                       \tag{5}
\]

up to one common sign.  This follows from the standard Jacobi-`cd`
construction and its multiplication/nesting identity.  Ng--Tsang prove its
least-deviation property; Wang's elliptic-isogeny formula gives (5)
directly.  Müller-Hermes--Szehr, Lemma 2, records the equivalent extremal
value `sqrt(k(q^L))`.

The DCT alternation identity is

\[
 U\operatorname{diag}((-1)^j)U^T=\pm J,                \tag{6}
\]

where `J` reverses coordinates.  Equations (4)--(6) therefore give

\[
 B_L(T_c)=\pm\sqrt{k(q^L)}\,D J D^{-1}.                 \tag{7}
\]

The largest singular value of `D J D^{-1}` is `c^{-L/2}`, so

\[
 \|B_L(T_c)\|^2={k(c^{2L})\over c^L}.                  \tag{8}
\]

If `S T_c S^{-1}` is a contraction, von Neumann's inequality applied to
`B_L` gives `||B_L(T_c)||<=cond(S)`.  Taking the infimum over similarities
proves

\[
 \boxed{t_*(T_c)\ge {k(c^{2L})\over c^L}.}              \tag{9}
\]

Thus any matching upper metric is automatically optimal.

## 4. The scalar bound is strictly below four

Jacobi's product is

\[
 k(Q)=4\sqrt Q\prod_{m\ge1}
  \left({1+Q^{2m}\over1+Q^{2m-1}}\right)^4.             \tag{10}
\]

With `Q=c^(2L)`, every factor is strictly less than one.  Hence

\[
 {k(c^{2L})\over c^L}<4,\qquad
 {k(c^{2L})\over c^L}=4-16c^{2L}+O(c^{4L}).             \tag{11}
\]

The axis therefore has a positive complete-similarity gap for every fixed
`c>0`, although the gap becomes extremely flat as either `c->0` or `p`
increases.

## 5. Matching upper bound in the first two live sizes

For `p=3`, odd functional calculus makes

\[
 T_c=\begin{pmatrix}0&u&0\\ l&0&u\\0&l&0\end{pmatrix},
 \quad u^2={k(c^2)\over2c},\quad l=cu.
\]

Choose `P=diag(1,r,r^2)`, where the smaller root of

\[
 l^2r^2-r+u^2=0                                       \tag{12}
\]

is used.  Then `P-T_c^*P T_c` is positive semidefinite of rank one.
The descending Landen identity gives

\[
 r^2={k(c^4)\over c^2}.
\]

Thus (9) is attained.  For `p=4`, put

\[
 \eta=\operatorname{sn}(K/3\mid k^2),\qquad
 P=\operatorname{diag}\left(
  1,{k\eta\over c},{k^2\eta^3\over c^2},
  {k^3\eta^4\over c^3}\right).                         \tag{13}
\]

The Jacobi addition formula at `K/3+2K/3=K` gives
`k(c^6)=k^3 eta^4`.  Directly splitting
`P^(1/2)T_cP^(-1/2)` into its two parity blocks gives singular values
`1,1,1,k^2 eta^2`; moreover the diagonal entries in (13) increase from one
to `k(c^6)/c^3`.  Hence (9) is attained here as well.  This is the explicit
finite-size computation in Li's thesis §4.4, reproduced in L116's
normalization.  L20 separately proves the condition-two statement for the
much larger arbitrary-weight `4 x 4` elliptic slice.

## 6. Numerical regeneration and the remaining identity

Run

```bash
.venv/bin/python -u experiments/crabb_elliptic_axis.py \
  --output experiments/crabb_elliptic_axis_s70222.jsonl
```

For `p=3,...,10` and `c in {0.05,0.15,0.4,0.6}`, the unrestricted L21 SDP,
the diagonal-metric SDP, and (1) agree to maximum absolute error
`1.66e-5` (the unrestricted error is below `5.7e-8`; the larger figure is a
nearly singular diagonal-solver case).  Identity (6) holds to `2.8e-15`.

The diagonal optimizer `P=diag(p_0,...,p_L)` consistently has

\[
 p_0=1,\quad p_L=t_*,\quad p_jp_{L-j}=t_*,             \tag{14}
\]

and `P-T_c^*P T_c` is rank one, within solver accuracy.  Set

\[
 W=DPD=\operatorname{diag}(w_j),\qquad w_j=c^jp_j.
\]

Then the desired all-size certificate reduces exactly to

\[
 W-FWF=bb^T\succeq0,\qquad
 w_jw_{L-j}=k(c^{2L}).                                 \tag{15}
\]

In the nodal basis, with `G=U^TWU` and `beta=U^Tb`,

\[
 G_{ij}={\beta_i\beta_j\over1-\tau_i\tau_j}.            \tag{16}
\]

Equation (16) is a finite Szego-kernel Gramian.  The best next attack is to
derive (15)--(16) from the Chebyshev--Blaschke model/Clark weights or an
elliptic Gauss--Lobatto quadrature formula.  This targets the conceptual
source of Li's identities instead of proving separate even/odd index
formulas by brute force.

## 7. Scope

Closing (15) for every `p` would prove the complete Crouzeix bound on the
entire pure elliptic Crabb axis and give the exact optimal constant (1).
It would not by itself prove a neighbourhood theorem near `C_p`: L115 still
requires uniform absorption of transverse coercive directions around this
axis.  Nor would it prove the general Crouzeix conjecture.
