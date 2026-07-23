# The corrected Faber--Blaschke dual square (A99, 2026-07-23)

## 1. Exact interpolating inner function

Fix a noncentral real one-grade direction `1<=k<L/2`.  Let

\[
P_0=2,\qquad P_1=z,\qquad P_m=zP_{m-1}-cP_{m-2}
\]

be the ellipse Dickson polynomials.  For a real correction parameter
`lambda`, define

\[
\boxed{
G_{a,c,\lambda}(z)
=P_L(z)+2a\{P_k(z)+P_{L-k}(z)\}
+\lambda a c^kP_{L-k}(z).
}                                                     \tag{1}
\]

For sufficiently small `a,c`, let `zeta_1,...,zeta_L` be its roots,
put `w_j=phi_c(zeta_j)`, and form

\[
B_{a,c,\lambda}(w)
=\omega\prod_{j=1}^L
\frac{w-w_j}{1-\overline{w_j}w}.                     \tag{2}
\]

The harmless unit phase `omega` is arbitrary.  Equation (2) is a
finite Blaschke product whenever the mapped roots lie in the disk,
which holds locally by continuity from the two axes.

The construction has the two exact specializations needed for a sharp
dual:

1. At `c=0`, the correction vanishes and

   \[
   G_{a,0,\lambda}
   =z^L+2a(z^k+z^{L-k}),
   \]

   the characteristic factor of L123's disk equality matrix.  Hence
   (2) is its characteristic Blaschke product.
2. At `a=0`, `G=P_L`; its roots are the elliptic Chebyshev nodes, so
   (2) is the degree-`L` Chebyshev--Blaschke product used in
   L116--L120.

Thus (2) interpolates the correct inner extremal on **both** equality
axes.  The fixed Chebyshev product fails this test: it already loses
order `a^2` on the disk equality family.

For

\[
T(a,c)=\phi_c(X(a)+cX(a)^*)
\]

on the exact phase-palindromic disk chart, inner-function monotonicity
gives

\[
t_*(T(a,c))\ge t_*(B_{a,c,\lambda}(T(a,c)))
\ge \|B_{a,c,\lambda}(T(a,c))\|^2.                  \tag{3}
\]

The open issue is to choose `lambda` sharply and pair this dual with
L138--L139's primal defect square.

## 2. The correction `lambda=2` is forced (L140)

L131 gives the exact exterior-boundary identity

\[
P_j(\zeta+c/\zeta)=\zeta^j+c^j\zeta^{-j}.
\]

For one phase-palindromic pair, the first new negative Fourier term
strictly below the central collision is therefore

\[
2a u_kc^k\zeta^{-k}.                                \tag{4}
\]

At the monomial degree-`L` inner function, write a numerator variation
as

\[
\delta N(w)=\alpha w^{L-k}.
\]

Since

\[
B=\frac{N}{N^\sharp},\qquad
\delta N^\sharp(w)=\overline\alpha w^k,
\]

the logarithmic boundary tangent has the paired frequencies

\[
\frac{\delta B}{B}
=\alpha w^{-k}-\overline\alpha w^k.                 \tag{5}
\]

Thus the unique degree-`L` inner tangent that absorbs (4) has
`alpha=2a u_kc^k` (with the conjugate placed according to the phase
gauge).  In the Faber polynomial this is exactly

\[
\boxed{
2a\overline{u_k}c^kP_{L-k}.
}                                                     \tag{6}
\]

For the real phase-one direction, (6) is `lambda=2` in (1).
Uniqueness is just Fourier orthogonality: no other numerator
coefficient produces the mode `w^(-k)`.

This proves the **location** of the scalar dual square without a root
calculation or fit.  It does not yet prove the curvature of the
operator-norm loss.

## 3. The mismatch curvature is exact (L141)

L140 identifies the unique reflected numerator coefficient.  The
quadratic cost of missing that coefficient can also be computed without
tracking roots.  In the associated Hardy model put

\[
 \epsilon=(\lambda-2)ac^k,\qquad
 N_\epsilon(w)=w^L+\epsilon w^{L-k}.
\]

Then

\[
 B_\epsilon(w)
 ={N_\epsilon(w)\over N_\epsilon^\sharp(w)}
 ={w^L+\epsilon w^{L-k}\over1+\overline\epsilon w^k}.
                                                               \tag{7}
\]

Let `C=C_(L+1)` be the canonical Crabb shift, whose first and last
weights are `sqrt(2)` and whose interior weights are one.  Since
`C^(L+1)=0`, (7) gives the exact finite identity

\[
 B_\epsilon(C)
 =(1-|\epsilon|^2)C^L+\epsilon C^{L-k}.              \tag{8}
\]

Indeed, all later terms in the inverse denominator contain
`C^(L+k)` or a higher power.  Now

\[
 C^L=2e_0e_L^*,
\]

while `C^(L-k)` has entries `sqrt(2),1,...,1,sqrt(2)` on
the `(L-k)`th superdiagonal.  Only columns `L-k` and `L` can couple
at the top singular value.  On those columns the relevant block of
`B_\epsilon(C)^*B_\epsilon(C)` is

\[
\begin{pmatrix}
2r&2\sqrt2\,\overline\epsilon(1-r)\\
2\sqrt2\,\epsilon(1-r)&4(1-r)^2+2r
\end{pmatrix},
\qquad r=|\epsilon|^2.                                \tag{9}
\]

All other eigenvalues are at most `2r`.  Expanding the larger
eigenvalue of (9) at `r=0` therefore yields

\[
 \boxed{\|B_\epsilon(C)\|^2=4-4|\epsilon|^2
        +O(|\epsilon|^4).}                            \tag{10}
\]

Consequently the part of the associated dual loss caused by a
coefficient mismatch is exactly

\[
 \boxed{4|\lambda-2|^2a^2c^{2k}.}                    \tag{11}
\]

This proves the quadratic term in `lambda` and independently confirms
that L140's correction is the unique minimizer.  It does **not** compute
the loss at the minimizer: the corrected family itself still carries
the base curvature described next.

## 4. The remaining base curvature

The deterministic probe finds the dimension- and grade-independent
limit

\[
\boxed{
\frac{
4-\|B_{a,c,\lambda}(T(a,c))\|^2
}{
16a^2c^{2k}
}
\longrightarrow
4+\frac{|\lambda-2|^2}{4}
}                                                     \tag{12}
\]

as `(a,c)->(0,0)` in the noncentral one-grade regime.

Consequences:

- the natural Faber-root choice `lambda=0` has loss
  `80a^2c^(2k)`;
- the explicit correction

  \[
  \boxed{\lambda=2}                                  \tag{13}
  \]

  improves the loss to `64a^2c^(2k)`;
- (13) is exactly a correction of the high reflected coefficient:

  \[
  G_{a,c,2}
  =P_L+2aP_k+2a(1+c^k)P_{L-k}.                      \tag{14}
  \]

The scalar square (12) is the dual counterpart of A98's primal target.
Its coefficient two matches L139's defect outer factor and its
minimum is

\[
4-64a^2c^{2k}+o(a^2c^{2k}),                         \tag{15}
\]

the same face proved centrally by L130 and found by every exact
noncentral Hessian audit.

This is a significant sharpening: the missing correction is no longer
an unspecified zero displacement or a many-parameter Blaschke
optimization.  It is the single explicit term (13)--(14).

Here L141 proves the `|\lambda-2|^2/4` part of (12).  The remaining
unproved scalar statement is only

\[
\boxed{
4-\|B_{a,c,2}(T(a,c))\|^2
=64a^2c^{2k}+o(a^2c^{2k}).
}                                                     \tag{16}
\]

## 5. Numerical adversarial checks

The persisted scan uses the coprime/noncentral cases

\[
(L,k)=(5,2),(7,3),(9,4)
\]

at `c=.18`, symmetrizes amplitudes `a=+/- .01`, and checks
`lambda=-2,0,1,2,3,4`.  After division by `16c^(2k)`, the six limiting
predictions are

\[
8,\ 5,\ 4.25,\ 4,\ 4.25,\ 5.                       \tag{17}
\]

All 18 finite-`c` values lie within `.17` of (17).  In every case
`lambda=2` improves the natural product, and all mapped zeros remain
strictly inside the disk.  Separate derivative-free optimization over
every Faber coefficient returns a dominant coefficient
`lambda_(L-k)=1.99...`; the remaining fitted coefficients affect only
higher finite-`c` terms.

The scan also falsifies two tempting shortcuts:

1. the fixed Chebyshev--Blaschke product has an unweighted
   order-`a^2` loss and cannot see A84's face;
2. the uncorrected root interpolation is valid and close, but its
   leading coefficient is 80 rather than 64.

## 6. The central collision has exact base curvature (L142)

When `L=2k`, reflection lands on the original grade.  Put `r=c^k`.
The corrected polynomial is then

\[
 G^{\rm cen}_{a,c}
 =P_{2k,c}+2a(1+c^k)P_{k,c}.                         \tag{18}
\]

The Dickson product identity

\[
 P_{2k,c}=P_{k,c}^2-2c^k
\]

factors (18) exactly as

\[
G^{\rm cen}_{a,c}
=Q_{a,r}(P_{k,c}),\qquad
Q_{a,r}(y)=y^2+2a(1+r)y-2r.                          \tag{19}
\]

Let `B_(k,c)` be the degree-`k` Chebyshev--Blaschke product and
`b_(a,r)` the degree-two root Blaschke product defined by `Q_(a,r)`.
The proper-map identity

\[
B_{k,c}\circ\phi_c=\phi_r\circ P_{k,c}
\]

and equality of the zero multisets give, up to phase,

\[
\boxed{
B^{\rm cen}_{a,c}=b_{a,r}\circ B_{k,c}.
}                                                     \tag{20}
\]

L126 reduces `B_(k,c)(T)` orthogonally to the physical size-three
central pencil.  In coefficient coordinates its matrix and coordinate
metric are

\[
\begin{aligned}
S&=
\begin{pmatrix}0&2&0\\r&0&1\\0&2r&0\end{pmatrix}
+a\begin{pmatrix}-2r&0&2\\0&0&0\\2r&0&-2\end{pmatrix},\\
K&=\begin{pmatrix}1/2&a&0\\a&1&a\\0&a&1/2\end{pmatrix}.
\end{aligned}                                        \tag{21}
\]

If `y_+`, `y_-` are the roots of `Q_(a,r)` and
`w_\pm=\phi_r(y_\pm)`, put `s=w_++w_-` and `p=w_+w_-`.
Then on this block

\[
b_{a,r}(\phi_r(S))
=
\{\phi_r(S)^2-s\phi_r(S)+pI\}
\{I-s\phi_r(S)+p\phi_r(S)^2\}^{-1}.                  \tag{22}
\]

L125's recurrence supplies the complete scalar map through total
parameter degree four:

\[
\phi_r(z)
=(1+2r^2+r^4)z-(r+3r^3)z^3
+(r^2+5r^4)z^5-r^3z^7+r^4z^9+O(r^5).               \tag{23}
\]

Newton reduction of `y_\pm` by (19), followed by simple generalized
singular-value perturbation for the pencil

\[
b_{a,r}(\phi_r(S))^*K\,b_{a,r}(\phi_r(S))-\mu K,
\]

gives zero at total degrees one, two, and three, and

\[
\boxed{
\mu_{\max}
=4-64a^2r^2-16r^4+O_{\rm total}(5).
}                                                     \tag{24}
\]

The calculation is exact over `Q[a,r]`.  Its pure-axis term independently
reproduces L117, while its mixed term is the desired curvature.  L129's
uniform inner gap makes the size-three block the active block of (20).
Consequently, for every central grade,

\[
\boxed{
[a^2c^{2k}]
\|B^{\rm cen}_{a,c}(T(a,c))\|^2=-64.
}                                                     \tag{25}
\]

Thus the corrected base curvature is now a theorem on all central
collisions.  The remaining base-curvature problem in (16) is purely
the noncentral localization from a grade `k` in length `L>2k` to this
central model.

## 7. The noncentral base curvature is also exact (L143)

Fix `1<=k<L/2` and use the corrected family (14).  Write

\[
B(a,c)=B_0(c)+aB_1(c)+a^2B_2(c)+O(a^3)
\]

for its value at the exact matrix pencil in L123 coefficient
coordinates, and write

\[
K(a)=K_0+aK_1
\]

for the coordinate Gramian.  Then

\[
\boxed{
[a^2c^j]\|B(a,c)\|_K^2=0\quad(j<2k),\qquad
[a^2c^{2k}]\|B(a,c)\|_K^2=-64.
}                                                     \tag{26}
\]

This is an all-size coefficient identity, including unequal residue
chains and the possible terminal alias.

### 7.1 Coefficient preparation

Put

\[
g_{a,c}(w)=G_{a,c,2}(\Psi_c(w)).
\]

Its local Weierstrass factorization is

\[
g_{a,c}(w)=N_{a,c}(w)O_{a,c}(w),                    \tag{27}
\]

where `N` is monic of degree `L` and `O(0)` is a unit.  The finite
Blaschke product is `N/N^sharp`.  This avoids individual roots
completely.

There is a coefficientwise proof of (27) and (26).  Let `Pi_<` and
`Pi_>=` select scalar powers below `L` and at least `L`.  If
`N_(r,s)`, `O_(r,s)`, and `g_(r,s)` denote the coefficients of
`a^r c^s`, put

\[
\begin{aligned}
R_{r,s}
&=g_{r,s}
 -\mathop{\sum_{(u,v)+(x,y)=(r,s)}}_
 {(u,v),(x,y)\ne(0,0),(r,s)}
 N_{u,v}O_{x,y},\\
N_{r,s}&=\Pi_<R_{r,s},\qquad
O_{r,s}=w^{-L}\Pi_{\ge}R_{r,s}.                     \tag{28}
\end{aligned}
\]

Starting with `N_(0,0)=w^L`, `O_(0,0)=1`, total-degree induction
makes (28) the unique preparation.  Substitute the differentiated
Dickson recurrence, L125's mutually inverse `Psi_c` and `phi_c`
recurrences, and L123's companion pencil.  Simultaneous induction on
`s` gives the following complete data that can reach the top
generalized singular pair through `c^(2k)`:

\[
\boxed{
[c^k]B_1
=4\left\{
\sum_{i=k}^{L}e_ie_{i-k}^*
+e_ke_{L-2k}^*
-e_0e_{L-k}^*
\right\}.
}                                                     \tag{29}
\]

Put `d=L-k>k`.  The only endpoint alias at or below the target is

\[
\begin{aligned}
[c^j]B_1e_L
&=-2e_d &&\text{if }j=d\le2k,\\
&=0 &&\text{otherwise for }j\le2k,\\
[c^j]e_0^*B_2e_L
&=4 &&\text{if }j=d<2k,\\
&=-16+4\,{\bf1}_{d=2k} &&\text{if }j=2k,\\
&=0 &&\text{otherwise for }j\le2k.
                                                               \tag{30}
\end{aligned}
\]

The same induction gives

\[
\boxed{
\operatorname{val}_c
\{(H_1-4K_1)e_L\}>k,
}                                                     \tag{31}
\]

where

\[
\begin{aligned}
H_1={}&B_1^*K_0B_0+B_0^*K_0B_1+B_0^*K_1B_0.
\end{aligned}
\]

For clarity about the all-size step in (28)--(31): every scalar-map
coefficient is retained symbolically.  Terms containing a nonleading
coefficient of `Psi_c` pair with the corresponding inverse coefficient
of `phi_c` and cancel by `phi_c(Psi_c(w))=w`.  The remaining edge sums
are the Catalan convolution from L125.  The Dickson recurrence then
moves an index by one, so before the terminal fold the induction is
translation invariant.  A fold can reach the top column only at
`d=L-k`; it produces exactly the two alias entries in (30).  This is
why unequal residue-chain lengths do not enter (29), and why `L=3k`
is the sole resonance on the target face.

### 7.2 Top singular calculation

Through this order the axis image has

\[
B_0e_L=2e_0,\qquad K_{0,LL}=K_{0,00}=\frac12.        \tag{32}
\]

Equation (31) says the first generalized singular-vector correction
has valuation greater than `k`.  Its Schur square therefore has
valuation greater than `2k` and cannot affect (26).

It remains to take the top Rayleigh coefficient.  In the convention
`B=B_0+aB_1+a^2B_2+...`, its quadratic Gram coefficient is

\[
\begin{aligned}
H_2={}&B_2^*K_0B_0+B_0^*K_0B_2+B_1^*K_0B_1\\
 &+B_1^*K_1B_0+B_0^*K_1B_1.
\end{aligned}
\]

From (30), the two
`B_0`--`B_2` cross terms give

\[
-32c^{2k}+8{\bf1}_{d\le2k}c^d.                      \tag{33}
\]

When the alias is present, `K_{1,0d}=K_{1,d0}=1` and
`B_1e_L=-2c^de_d`; the two `B_0`--`B_1` metric cross terms give

\[
-8{\bf1}_{d\le2k}c^d.                               \tag{34}
\]

Thus the alias cancels identically, including at the resonance
`d=2k`, and

\[
[c^j](H_2)_{LL}=0\ (j<2k),\qquad
[c^{2k}](H_2)_{LL}=-32.                              \tag{35}
\]

Division by `K_(0,LL)=1/2` proves (26).  Diagonal unitary gauging
replaces `a^2` by `|a|^2`, so the result holds for an arbitrary
one-grade phase.

Combining L142 and L143 proves the base term (16) for every central or
noncentral one-grade direction.  Together with L140--L141, the full
one-grade dual parabola (12) is now exact, not numerical.

## 8. Distinct grades polarize diagonally (L144)

Choose one representative `1<=k<=L/2` from each reflected pair and
write the corrected multi-grade factor as

\[
\boxed{
\begin{aligned}
G_{a,c,u}
=P_L
2a\sum_{k<L/2}
\{u_kP_k+\overline{u_k}(1+c^k)P_{L-k}\}\\
\quad
2a\,u_{L/2}(1+c^{L/2})P_{L/2}
\qquad(L\ {\rm even}),
\end{aligned}}                                      \tag{36}
\]

where the central coefficient is real in the fixed phase gauge.
Map the roots as before and call the resulting inner function
`B_(a,c,u)`.  Then its complete diagonal principal face is

\[
\boxed{
\|B_{a,c,u}(T(a,c,u))\|^2
=4-64a^2\sum_{k\le L/2}|u_k|^2c^{2k}
+\text{terms strictly above the diagonal face}.
}                                                     \tag{37}
\]

Equivalently, for distinct grades `k` and `ell`,

\[
\boxed{
[a^2c^{k+\ell}]
\left(
\|B_{u_k+u_\ell}(T)\|^2
-\|B_{u_k}(T)\|^2
-\|B_{u_\ell}(T)\|^2+4
\right)=0.
}                                                     \tag{38}
\]

To prove this, polarize the coefficient recurrence (28) with
independent variables `u_k`, `overline(u_k)`.  Before a terminal fold,
the first faces (29) occupy distinct shift diagonals.  Their
`K_0` products therefore have zero top Fourier grade unless the two
grades coincide.  Equation (31) separately gives first singular-vector
valuations greater than `k` and `ell`, so a mixed Schur product has
valuation strictly greater than `k+ell`.

Only a folded endpoint could defeat that orthogonality.  For `k<ell`,
it can land on the mixed face only in one of the two resonances

\[
L=2k+\ell,\qquad L=k+2\ell.                          \tag{39}
\]

In either case the polarized version of (30) contributes `+8` through
the quadratic numerator and `-8` through the corresponding
off-diagonal entry of `K_1`.  Thus it cancels exactly, just as in
(33)--(34).  A central grade is the limiting double fold and obeys
the same recurrence, consistent with L142.

Finally, a complex cross monomial carries nonzero circle grade
`k-ell`; the endpoint functional selects grade zero.  Hence both its
real and imaginary polarizations vanish, proving (38) over complex
coefficients.  The diagonal terms are L142--L143, which proves (37).

## 9. A complementary rank-one Stein defect (L145)

The corrected inner function also constructs the missing primal
direction.  Write `B=N/D`, let `x,y` be its simple top right and left
singular vectors,

\[
B(T)x=s\,y,
\]

and recall that its model space is

\[
K_B=\{p/D:\deg p<L\}.                                \tag{40}
\]

Choose a nonzero vector `q` by

\[
\boxed{
q\perp D(T)^{-1}
\operatorname{span}\{x,Tx,\ldots,T^{L-1}x\}.
}                                                     \tag{41}
\]

The orbit has rank `L` at both equality axes, so (41), with one scalar
normalization, defines an analytic one-dimensional complement locally.
Let

\[
P_T(q)=\sum_{n\ge0}(T^*)^nqq^*T^n.
\]

For an orthonormal basis `f_0,...,f_(L-1)` of `K_B`, L127's model
kernel identity is

\[
P_T(q)-B(T)^*P_T(q)B(T)
=\sum_{j=0}^{L-1}f_j(T)^*qq^*f_j(T).                 \tag{42}
\]

Equation (41) makes the right side annihilate `x`.  Hence

\[
\boxed{
x^*P_T(q)x=s^2y^*P_T(q)y.
}                                                     \tag{43}
\]

This is exact complementary slackness, not a fitted defect ansatz.

For a one-grade direction, differentiate (41) in the coefficient
gauge.  L143's right-singular correction has valuation greater than
`k`; the orbit equations are triangular at the disk endpoint.  Their
unique normalized solution is the homogeneous defect tangent.  Insert
it in the Stein recurrence and use the two simple endpoint generalized
eigenvalues.  Coefficient induction gives

\[
\boxed{
[a^2c^j]\{
\kappa(P_T(q))-\|B(T)\|^2
\}=0,\qquad 0\le j\le2k.
}                                                     \tag{44}
\]

Thus L142--L143's dual coefficient `-64` is attained by a feasible
rank-one Stein metric for every central and noncentral one-grade
direction.  Formula (44) also explains the defect corrections in A98:
they are the coefficients of the model-space complement (41), not an
ad hoc coordinate transport.

## 10. Positivity forces all mixed primal faces (L146)

Let `Q_pr(c)` be the Hermitian amplitude Hessian of the locally
optimized rank-one Stein envelope, and let `Q_du(c)` be the Hessian of
the corrected dual inner lower bound.  For every fixed small `c`,

\[
\boxed{
Q_{\rm pr}(c)-Q_{\rm du}(c)\succeq0,
}                                                     \tag{45}
\]

because the primal envelope is an upper bound for `t_*` and the
Blaschke norm is a lower bound, with equality and equal first
derivative on the axis.

L125's grade filtration, applied to the Stein and preparation
recurrences, makes the `(k,ell)` entry divisible by `c^(k+ell)`.
With

\[
D_c=\operatorname{diag}(c^k)_{k\le L/2},
\]

the associated limit

\[
E_0=\lim_{c\to0}
D_c^{-1}\{Q_{\rm pr}(c)-Q_{\rm du}(c)\}D_c^{-1}
\succeq0                                             \tag{46}
\]

therefore exists.  L145 says every diagonal entry of `E_0` is zero:
primal and dual agree in every one-grade direction.  A positive
semidefinite Hermitian matrix with zero diagonal is zero, since each
`2x2` principal minor forces its off-diagonal entry to vanish.
Consequently

\[
\boxed{
Q_{\rm pr,face}=Q_{\rm du,face}
=-64\sum_{k\le L/2}|u_k|^2c^{2k}.
}                                                     \tag{47}
\]

This proves A98's completed-square **minimum and all mixed
coefficients** without guessing a raw defect-coordinate reversal.
The explicit minimizer is obtained recursively from (41) and L139's
positive defect Hessian.  What remains for a local theorem is no
longer the Hessian sign: it is the uniform analytic remainder needed
to dominate terms above the Newton face.

## 11. Proof obligations

The complete primal/dual principal face and complementary defect are
now proved coefficientwise.  Uniform remainders remain.  The
maintainable route is:

1. Inner--outer factor

   \[
   G_{a,c,\lambda}(\Psi_c(w))
   =B_{a,c,\lambda}(w)O_{a,c,\lambda}(w)
   \]

   and expand the outer factor in the Faber/Hardy filtration.
2. Establish a uniform analytic remainder in the amplitude/ellipse
   tube, rather than only the coefficientwise formal statement.

If these steps close, (3) and (12)--(16) give the sharp lower face while the
matching L138--L139 defect vector gives the upper face.  That would
prove A84's all-grade principal coefficient rather than merely
another sufficient estimate.

## 12. Reproduction

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_faber_blaschke_dual.py \
  --output experiments/crabb_faber_blaschke_dual_s70223.jsonl
```

The script constructs the exact disk/ellipse anchors, maps the roots
with the elliptic Riemann formula, evaluates the resulting rational
inner function by matrix Horner rules, and uses only the largest
singular value.  It does not call the similarity SDP or the rank-one
envelope optimizer.

The exact central calculation is

```bash
.venv/bin/python -u \
  experiments/crabb_central_faber_blaschke.py \
  --output experiments/crabb_central_faber_blaschke_s70223.jsonl
```

It works over a symbolic polynomial ring, checks the Dickson
composition, lifts the simple generalized singular pair through total
degree four, and asserts the coefficient `-64`.

The exact noncentral preparation and generalized singular calculation
is

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_faber_blaschke_formal.py \
  --output experiments/crabb_faber_blaschke_formal_s70223.jsonl
```

It retains rational `c`-series through `2k`, prepares the numerator by
(28), checks (29) and the endpoint alias, and verifies that the Schur
term vanishes on the face while the direct Rayleigh coefficient is
`-64`.

Mixed-grade polarization, including both resonances in (39), is

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_faber_blaschke_mixed.py \
  --output experiments/crabb_faber_blaschke_mixed_s70223.jsonl
```

It subtracts the two one-grade quadratic series from the combined
series and asserts exact zero through twice the larger grade.

The model-space complementary defect and matching primal condition
Hessian are regenerated by

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_faber_blaschke_complementarity.py \
  --output \
  experiments/crabb_faber_blaschke_complementarity_s70223.jsonl
```

It solves (41) over exact rational series, verifies the differentiated
orthogonality equations, and checks (44) in central and noncentral
grades.
