# The operator-valued Hardy equality manifold at repeated Crabb blocks

## 1. Result (L193, 2026-07-24)

Fix a Crabb length \(L\ge2\) and a copy multiplicity \(m\ge1\).
Let \(S=S_L\otimes I_m\) be the block level shift on
\(\mathbb C^{L+1}\otimes\mathbb C^m\).  Let \(H\) be positive
definite on the first \(L\) levels and extend it by zero on the final
level.  Put

\[
 K=H+S^*HS,\qquad A=2K^{-1}HS.                      \tag{1}
\]

Write \(E_j:\mathbb C^m\to\mathbb C^{L+1}\otimes\mathbb C^m\)
for the \(j\)-th level embedding and set

\[
 D=E_0^*HE_0,\qquad
 Q=HE_0D^{-1/2}=KE_0D^{-1/2}.                       \tag{2}
\]

Let \(M\) be the rank-\(m\) Stein Gramian

\[
 M-A^*MA=QQ^*.                                      \tag{3}
\]

There is a finite block Hardy residual

\[
 \Psi_m(H)\in
 M_{(L-1)m}(\mathbb C)                              \tag{4}
\]

with the following properties near
\(H_0=I_L/2\otimes I_m\).

1. **Exact complete-\(2\) upper equality.**

   \[
   \boxed{K\preceq M\preceq4K}
   \tag{5}
   \]

   locally on the lower side, and

   \[
   ME_0=KE_0,\qquad
   \Psi_m(H)=0\Longleftrightarrow ME_L=4KE_L.
   \tag{6}
   \]

   Hence every residual-zero point has a feasible contraction
   similarity metric of condition square exactly four.

2. **Exact inverse-block-Toeplitz classification.**  Let
   \(\widehat H\) be the restriction of \(H\) to its first \(L\)
   levels.  Then, locally,

   \[
   \boxed{\Psi_m(H)=0
   \Longleftrightarrow \widehat H^{-1}
   \text{ is Hermitian block Toeplitz}.}             \tag{7}
   \]

   Thus the equality/certificate manifold is not merely supplied by
   an implicit-function argument: it is the inverse image of the open
   positive block-Toeplitz cone under matrix inversion.

   At the repeated Crabb point,

   \[
   \boxed{
   (D\Psi_m E)_{r,c}
   =E_{c+1,L-1-r}-E_{c,L-2-r},}
   \tag{8}
   \]

   where every displayed entry is an \(m\times m\) block.  After
   reversing the residual block rows, (7) is the Hermitian block
   diagonal-difference map.  Its real rank is

   \[
   ((L-1)m)^2,                                      \tag{9}
   \]

   and its kernel is exactly the Hermitian block-Toeplitz space,
   of real dimension

   \[
   (2L-1)m^2.                                       \tag{10}
   \]

   Thus \(\Psi_m^{-1}(0)\) is a real-analytic manifold through the
   fixed repeated Crabb block.

3. **New noncommuting equality/certificate directions.**  The block
   Toeplitz coefficients in (7) are arbitrary copy matrices, not a
   commuting normal tuple.  The analytic equality graph therefore
   contains genuinely noncommuting off-diagonal copy directions.
   These are not direct sums of scalar equality anchors in a copy
   unitary gauge.

L193 is an exact repeated-block equality **upper certificate**.  It
does not yet prove a full neighbourhood theorem: away from the
residual-zero manifold, the upper endpoint is an \(m\times m\)
matrix and its successive kernels require an operator-valued metric
flag.

## 2. The block disk factorization

The scalar L122 factorization uses no commutativity.  For
\(|w|=1\),

\[
\begin{aligned}
K-\overline wHS-wS^*H
&=(I-wS^*)H(I-\overline wS)\succeq0.                \tag{11}
\end{aligned}
\]

The right side kills

\[
 f(w)\otimes v,\qquad
 f(w)=(1,w,\ldots,w^L)^T,\quad v\in\mathbb C^m.
 \tag{12}
\]

Consequently

\[
W(K^{1/2}AK^{-1/2})=\overline{\mathbb D},            \tag{13}
\]

with an \(m\)-dimensional top support space at every boundary angle.
This is the natural block-valued disk chart through
\(I_m\otimes C_{L+1}\).

## 3. Rank-\(m\) kernel observability

Put

\[
 T=K^{1/2}AK^{-1/2},\qquad
 V=K^{1/2}E_0D^{-1/2}.                              \tag{14}
\]

Then

\[
TV=0,\qquad V^*V=I_m.                               \tag{15}
\]

Berger's theorem gives an isometry \(J\) and a unitary \(U\) such
that

\[
T^j=2J^*U^jJ,\qquad j\ge1.                          \tag{16}
\]

For nonzero integers \(r\),

\[
(JV)^*U^r(JV)
={1\over2}V^*T^rV=0,                               \tag{17}
\]

with the negative powers following by adjoint.  Thus the subspaces

\[
U^{-j}JV\mathbb C^m,\qquad j\ge0,                   \tag{18}
\]

are mutually orthogonal.  Operator-valued Bessel gives

\[
\sum_{j\ge0}(T^*)^jVV^*T^j\preceq4I.               \tag{19}
\]

Conjugating (18) by \(K^{1/2}\) gives \(M\preceq4K\),
because the Stein forcing becomes

\[
K^{1/2}VV^*K^{1/2}=QQ^*.                            \tag{20}
\]

Also \(A^jE_0=0\) for \(j\ge1\), and

\[
Q^*E_0=D^{1/2}.
\]

Therefore

\[
ME_0=QQ^*E_0=HE_0=KE_0,                            \tag{21}
\]

which is the first identity in (6).

At the repeated Crabb point, the remaining generalized eigenvalues
of \((M,K)\) equal two.  They remain uniformly separated from one
and four nearby.  Thus the exact level-zero eigenvalue in (20) is
the lower endpoint.  This proves the local lower half of (5);
alternatively it follows by continuity after restricting to the
orthogonal complement of \(E_0\mathbb C^m\).

The normalization \(D^{-1/2}\) in (2) is essential.  Using the raw
forcing \(HE_0E_0^*H\) would put the endpoint matrix \(D\), rather
than the identity, in (20), and would introduce the spurious factor
\(\kappa(D)\) into the condition number.

## 4. Operator-valued Hardy defect

The concrete dilation map is again

\[
({\cal J}_Hx)(w)
=H^{1/2}(I-\overline wS)(I-\overline wA)^{-1}x.
 \tag{22}
\]

The Fourier moment proof of L183 is block-valued without change:

\[
\int_{\mathbb T}{\cal J}_H^*{\cal J}_H\,dm=K,\qquad
\int_{\mathbb T}w^j{\cal J}_H^*{\cal J}_H\,dm
={1\over2}KA^j.                                     \tag{23}
\]

The constant orbit is \(H^{1/2}E_0\mathbb C^m\).
Its coefficient projection is

\[
Q_H=I-E_0D^{-1}E_0^*H.                              \tag{24}
\]

If \(E_0^*Kx=0\), operator-valued Bessel equality gives

\[
\boxed{
x^*(4K-M)x
=4\left\|
H^{1/2}Q_H(I-\overline wS)
(I-\overline wA)^{-1}x
\right\|_{L^2(\mathbb T)}^2.}
 \tag{25}
\]

Take \(x=E_Lv\).  Define \(\Psi_m(H)\) by the first \(L-1\)
negative Fourier block coefficients of

\[
Q_H(I-tS)(I-tA)^{-1}E_L.                            \tag{26}
\]

The inverse-Gram reflection identity in the next section shows that
these first \(L-1\) block coefficients vanish if and only if the
whole interior orbit-complement residual vanishes.  Equations
(25)--(26) therefore give

\[
\Psi_m(H)=0
\Longleftrightarrow
E_L^*(4K-M)E_L=0.                                   \tag{27}
\]

Since \(4K-M\succeq0\), zero compression implies
\((4K-M)E_L=0\).  This is the second identity in (6).

## 5. Exact inverse-Gram reflection identity

Let \(B=\widehat H^{-1}\).  Write \(R_-\) and \(R_+\) for extraction
of levels \(0,\ldots,L-2\) and \(1,\ldots,L-1\), respectively, from
the \(L\)-level coefficient space.  Define

\[
B_-=R_-BR_-^*,\qquad B_+=R_+BR_+^*,\qquad
\Delta_B=B_+-B_-,\qquad \Sigma_B=B_++B_- .          \tag{28}
\]

Let \(P_I\) extract levels \(1,\ldots,L-1\) and let \(P_2\) extract
levels \(2,\ldots,L\) from the complete \(L+1\)-level space.  Direct
block-row elimination in \(KA=2HS\) gives the exact identity

\[
\boxed{
P_I(A-S)=-\Delta_B\Sigma_B^{-1}P_2.}                \tag{29}
\]

Here is the elimination explicitly.  Fix \(y\), put \(x=Ay\), and
write

\[
u=(x_0,\ldots,x_{L-1}),\quad
v=(x_1,\ldots,x_L),\quad
z=(y_1,\ldots,y_L),\quad p=2z-u.
\]

The first, interior, and last block rows of \(Kx=2HSy\) say that
there is a vector \(c\in\mathbb C^{(L-1)m}\) with

\[
\widehat Hv=(c,0),\qquad \widehat Hp=(0,c).          \tag{30a}
\]

Let \(t=P_2y\) and \(d=P_I(A-S)y\).  The leading coordinates of
\(v\) are \(t+d\), while the trailing coordinates of \(p\) are
\(t-d\).  Applying \(B=\widehat H^{-1}\) to (30a) and putting
\(r=c/2\) gives

\[
t+d=2B_-r,\qquad t-d=2B_+r.                         \tag{30b}
\]

Their sum gives \(t=\Sigma_Br\), and their difference gives
\(d=-\Delta_Br\), proving (29).  This derivation preserves the order
of every copy-space product and uses no commutativity.

Form the terminal Krylov matrix and its tail restriction

\[
{\cal V}=[E_L,AE_L,\ldots,A^{L-2}E_L],\qquad
{\cal Q}=P_2{\cal V}.                               \tag{31}
\]

Up to the harmless block transpose used to arrange its coefficients,
the finite residual is

\[
\Psi_m(H)^T=P_I(A-S){\cal V}
=-\Delta_B\Sigma_B^{-1}{\cal Q}.                   \tag{32}
\]

At \(H_0\), \({\cal Q}\) is a weighted block reversal and is
invertible.  It remains invertible nearby, while
\(\Sigma_B\succ0\).  Therefore

\[
\Psi_m(H)=0\Longleftrightarrow\Delta_B=0.           \tag{33}
\]

The last condition is precisely

\[
B_{i+1,j+1}=B_{ij}\quad(0\le i,j\le L-2),           \tag{34}
\]

so \(B\) is Hermitian block Toeplitz.  Conversely, block Toeplitz
\(B\) makes the entire interior row \(P_I(A-S)\) vanish, not only the
finite terminal residual.  This proves (7) and supplies an explicit
analytic parameterization by the positive block-Toeplitz cone.

This exact factorization also repairs a subtle codomain issue.  The
row reversal of the raw nonlinear residual is Hermitian only to first
order, not identically.  Thus one must not apply the real
implicit-function theorem to that raw matrix as though its nonlinear
codomain were Hermitian.  Equation (32), rather than that invalid
shortcut, proves the equality-manifold statement.

## 6. Differential and dimension

Let \(H=H_0+sE\).  Since \(B_0=2I\),

\[
\dot B=-4E.                                        \tag{35}
\]

Differentiating (29), or using the two endpoint paths in L187, gives
(8) entry by entry.  After residual block-row reversal it is the
Hermitian block diagonal-difference map.  It is surjective: integrate
each block diagonal difference and choose its mean to be zero.  Its
kernel consists exactly of constant block diagonals,

\[
E_{a,b}=Z_{b-a},\qquad Z_{-j}=Z_j^*.                \tag{36}
\]

The diagonal block \(Z_0\) is Hermitian and every positive offset
block is arbitrary complex.  Hence

\[
\dim_{\mathbb R}\ker D\Psi_m
=m^2+2(L-1)m^2=(2L-1)m^2,                           \tag{37}
\]

while

\[
(Lm)^2-(2L-1)m^2=((L-1)m)^2.                       \tag{38}
\]

## 7. Why this changes the repeated-block strategy

The old repeated-\(C_3\) campaign first diagonalized normal copy
coefficients and then controlled nonnormal terminal blocks by a
copy-space metric flag.  L193 shows that the natural equality
certificate stratum for general \(p\) is larger:

\[
H_{\rm eq}(Z_0,Z_1,\ldots,Z_{L-1})
=\operatorname{Toep}(Z_0,\ldots,Z_{L-1})^{-1},
\tag{39}
\]

where the \(Z_j\) need not commute or be normal.  The nonlinear
zero-residual correction absorbs their noncommutativity inside the
block disk chart.

Thus a general repeated-block proof should not try to force all flat
copy data onto direct sums of scalar equality anchors.  The correct
plan is:

1. use the block residual as the transverse coordinate;
2. use the kernel of its operator-valued Gram square as the next
   copy-space flag;
3. promote a stable kernel directly into the smaller-multiplicity
   block equality manifold (39); and
4. add block-valued elliptic and circular-normal fibres.

This may replace the \(p=3\)-specific scalar-support terminal
classification by induction on copy multiplicity.

## 8. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/crabb_block_hardy_equality.py \
  --output \
  experiments/crabb_block_hardy_equality_s70224.jsonl
```

The checker:

1. reconstructs the full real differential on every Hermitian
   coordinate and verifies (9)--(10);
2. verifies (29) at deterministic non-Toeplitz inputs;
3. constructs \(H\) as the inverse of a positive block-Toeplitz
   matrix with noncommuting coefficient blocks; and
4. verifies residual and inverse-Toeplitz errors at roundoff,
   positive \(H\), invertible terminal Krylov tail, rank-\(m\)
   generalized endpoint eigenspaces at exactly one and four,
   \(M\preceq4K\), and (6).

The standard run covers lengths three and four at multiplicities two
and three.  Equations (11)--(38) are the all-size proof.
