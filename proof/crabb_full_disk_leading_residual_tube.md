# The full-disk equality manifold and first-residual tube (L187--L189, 2026-07-24)

## 1. Results

Put \(p=L+1\), let \(S\) be the \(p\)-square unweighted shift, and
extend a positive \(L\)-square Hermitian matrix \(H\) by a final zero
row and column.  As in L122, set

\[
 K=H+S^*HS,\qquad A=2K^{-1}HS,\qquad h=H_{00},
 \qquad q=He_0=Ke_0.                                \tag{1}
\]

Let \(M-A^*MA=qq^*\), and define

\[
 Q_H=I-e_0{e_0^*H\over h},\qquad
 {\cal Y}_H(t)=Q_H(I-tS)(I-tA)^{-1}e_L.             \tag{2}
\]

The components of \({\cal Y}_H\) in \(e_1,\ldots,e_{L-1}\)
are determined by the finite \((L-1)\)-square matrix

\[
 \Psi(H)_{r,c}
 =[t^{r+1}]\,e_{c+1}^*{\cal Y}_H(t),
 \quad0\le r,c\le L-2.                              \tag{3}
\]

The new conclusions are:

1. **L187 (full-disk equality manifold).**  Near \(H_0=I/2\),
   \(\Psi^{-1}(0)\) is a real-analytic manifold of real codimension
   \((L-1)^2\).  Its tangent space at \(H_0\) is exactly the Hermitian
   Toeplitz space.  After removing the irrelevant positive scalar
   in \(H\), its physical dimension is \(2L-2\).

   Every point of this manifold has

   \[
   \lambda_{\min}(M,K)=h,\qquad
   \lambda_{\max}(M,K)=4h,                           \tag{4}
   \]

   and its characteristic finite Blaschke product has norm exactly
   two.  Thus these are genuine scalar-Crouzeix and similarity-square
   equality points, not merely zeros of a chosen upper certificate.
   The phase-palindromic family in L123 is an explicit lower-dimensional
   section of this larger curved manifold.

2. **L188 (universal first-residual face).**  Let \(H(s)\) be any
   real-analytic full-disk path through \(H_0\), and suppose

   \[
   \Psi(H(s))=s^mF+O(s^{m+1}),\qquad F\ne0.          \tag{5}
   \]

   For

   \[
   D_{2m}=-2[s^{2m}]
   \{\lambda_{\max}(M,K)-4\lambda_{\min}(M,K)\},
   \]

   one has the order-independent factor

   \[
   \boxed{D_{2m}=4\|F\|_F^2.}                       \tag{6}
   \]

   Put \(F^-=(F-F^T)/2\).  The leading complex response against the
   real/imaginary circular normal of mode \(k\) is

   \[
   \boxed{
   G^{(m)}_{L,k}
   ={4(4k-1)\over L^2}
   \sum_{\substack{1\le a<b\le L-1\\a+b=L-k}}
   (b-a)F^-_{ab}.}                                  \tag{7}
   \]

   The corresponding sum on \(a+b=L+k\) is its conjugate.  Therefore

   \[
   \boxed{
   D_{2m}\ge
   16\sum_{k=3}^{L-3}
   {|P_{L,k}(F^-)|^2\over\binom{L-k}{3}},}           \tag{8}
   \]

   where \(P_{L,k}\) is the sum in (7).  The right side is exactly
   the completed gain using only L173's flux curvature.  L173's
   positive null lift makes the actual Schur face strict whenever
   \(F\ne0\).

3. **L189 (nonlinear full-disk/circular-normal tube).**  In L173's
   transported true-circular-normal fibres over L122's **general**
   Hermitian disk chart, the explicit rank-one Stein certificate has
   condition square at most four in a full neighbourhood of the
   Crabb point.  Equality in this tube is precisely

   \[
   H\in\Psi^{-1}(0),\qquad y=0.                     \tag{9}
   \]

   This replaces the order-by-order sixth/eighth/tenth kernel ladder
   by one analytic first-residual argument.

L189 is the requested nonlinear merger over the complete general-\(H\)
disk base and its coercive circular-normal fibres.  It does not by
itself merge the remaining elliptic/marked compact charts into the
whole local operator quotient.

## 2. The finite residual and its differential

Expanding (2) gives, in the interior coordinates,

\[
 \Psi(H)_{r,c}
 =e_{c+1}^*(A^{r+1}-SA^r)e_L.                      \tag{10}
\]

The remaining coefficients obey the characteristic recurrence, so
the first \(L-1\) rows vanish if and only if the complete
orbit-complement Hardy residual vanishes.

Let \(H=I/2+sE\), where \(E\) is an arbitrary Hermitian \(L\)-square
matrix.  Linearizing (1), or simply using (10), gives

\[
 \boxed{
 (D\Psi_{H_0}E)_{r,c}
 =E_{c+1,L-1-r}-E_{c,L-2-r}.}                      \tag{11}
\]

Let \(J\) reverse \(L-1\) coordinates.  Then

\[
 (JD\Psi_{H_0}E)_{i,c}
 =E_{c+1,i+1}-E_{c,i}.                              \tag{12}
\]

For Hermitian \(E\), the right side is Hermitian in \(i,c\).  Conversely,
given any Hermitian matrix \(G\), solve (12) diagonal by diagonal:

\[
 E_{i,i+d}=E_{0,d}+\sum_{r=0}^{i-1}G_{r+d,r},
 \qquad0\le i<L-d.                                  \tag{13}
\]

The free constants \(E_{0,d}\) are exactly the Hermitian Toeplitz
kernel.  Choosing each diagonal mean to be zero gives a unique
right inverse.  Hence

\[
 \operatorname{rank}_{\mathbb R}D\Psi=(L-1)^2,
 \qquad
 \ker D\Psi=\{\hbox{Hermitian Toeplitz matrices}\}.  \tag{14}
\]

The analytic implicit-function theorem proves the manifold statement
in L187.  It also gives a convergent equality recentering over every
small Toeplitz coordinate, not only the phase-palindromic ones.

## 3. Why residual zero is genuine equality

L183's Hardy defect identity says that, for \(x\perp_K e_0\),

\[
 x^*(4hK-M)x
 =4h\|H^{1/2}Q_H(I-tS)(I-tA)^{-1}x\|_{H^2}^2.       \tag{15}
\]

Also

\[
 Me_0=hKe_0.                                        \tag{16}
\]

Since \(e_L\perp_K e_0\), equations (3), (10), and the finite
recurrence show

\[
 \Psi(H)=0
 \quad\Longleftrightarrow\quad
 Me_L=4hKe_L.                                       \tag{17}
\]

Near \(H_0\), the two endpoint eigenvalues remain simple.  Equations
(15)--(17) prove (4).

For completeness, this upper equality has a matching scalar lower
function.  Since \(Ae_0=0\) and \(A\) has rank \(L\), write

\[
 \det(\xi I-A)=\xi g(\xi),\qquad
 B(\xi)=g(\xi)/g^\sharp(\xi).                       \tag{18}
\]

The Stein identity makes \(A\) a strict contraction in the \(M\)
metric; local controllability excludes a unitary summand, so \(B\)
is a finite Blaschke product.  The standard defect-one model identity
for the contraction \(M^{1/2}AM^{-1/2}\) says that \(B(A)\) is the
rank-one partial isometry from \(\ker(A^*M)\) to \(\ker A\).

Now \(A^*Ke_L=0\), because \(KA=2HS\) and \(He_L=0\).  On
\(\Psi=0\), (17) therefore gives

\[
 \ker(A^*M)=\mathbb Ce_L,\qquad
 \ker A=\mathbb Ce_0.                               \tag{19}
\]

The \(M\)-partial isometry preserves the \(M\)-norm from \(e_L\) to
\(e_0\).  Using (16)--(17),

\[
 {\|B(A)e_L\|_K^2\over\|e_L\|_K^2}=4.               \tag{20}
\]

The similarity upper bound is two, so (20) proves
\(\|B(A)\|_K=2\) and \(t_*(A)=4\).

## 4. Ambient stationarity on the whole manifold

L162's stationarity argument does not require Toeplitz palindromy once
it is expressed through (2).  Here is the invariant form.

Keep \(K,q\) fixed and solve

\[
 M_E-(A+E)^*M_E(A+E)=qq^*.                          \tag{21}
\]

At a point of \(\Psi^{-1}(0)\), the upper condition ratio in (21)
and the norm square of the Blaschke product (18) touch at four, with
the same simple endpoint vectors \(e_0,e_L\).  Their first
derivatives therefore agree.

The Cauchy formula for \(DB(A)[E]\) gives an endpoint boundary
integral.  The first Schwarz/Riemann correction \(h_E\) is obtained
from the same support null vector

\[
 f(w)=(1,w,\ldots,w^L)^T
\]

and the density \(f(w)^*Kf(w)\).  Subtracting the two boundary
integrands leaves precisely the orbit-complement transfer in (2).
Thus

\[
 D\kappa(A)[E-h_E(A)]
 =2\operatorname{Re}\langle{\cal Y}_H,{\cal L}_E\rangle_{H^2}
                                                               \tag{22}
\]

for an analytic transfer vector \({\cal L}_E\).  On \(\Psi=0\), the
finite recurrence makes \({\cal Y}_H=0\), and hence

\[
 \boxed{D\kappa(A)[E-h_E(A)]=0\quad\hbox{for every }E.}          \tag{23}
\]

Equation (22) is the coordinate-free version of L162 equations
(5)--(11): L162's explicit phase-palindromic kernel identity is the
special case in which the vanishing transfer can be written as a
single rational polynomial quotient.

## 5. The universal endpoint factor

Suppose (5) holds.  Analytic minimization of (15) near the simple
top vector writes \(x(s)=e_L+u(s)\).  At \(s=0\), changing \(x\) in
an \(e_L\)-transverse direction produces only the constant Hardy
coordinates, whereas \(F\) occupies the negative coordinates.
Consequently the least-squares correction satisfies \(u=O(s^m)\)
and is orthogonal to \(F\) in the leading coefficient.

At the Crabb point,

\[
 h={1\over2},\qquad e_L^*Ke_L={1\over2},\qquad
 H|_{\operatorname{span}\{e_1,\ldots,e_{L-1}\}}
 ={1\over2}I.
\]

Substitution in (15) gives

\[
 [s^{2m}]\{4\lambda_{\min}-\lambda_{\max}\}
 =2\|F\|_F^2.                                       \tag{24}
\]

The definition of \(D_{2m}\) proves (6).  Notice that this argument
does not assume \(F\) is skew.  The generic first- and second-order
residuals are not skew.

## 6. Division of the normal response by the residual

Let \({\cal G}_{L,k}(H)\) be the complex pulled ambient derivative
in circular mode \(k\).  It is analytic near \(H_0\).  Equation (23)
shows that it vanishes on the analytic submanifold \(\Psi^{-1}(0)\).
Since (14) makes \(\Psi\) a submersion, the analytic Hadamard lemma
gives

\[
 {\cal G}_{L,k}(H)={\cal C}_{L,k}(H)[\Psi(H)]        \tag{25}
\]

for an analytic family of real-linear functionals
\({\cal C}_{L,k}\).

It remains only to calculate \({\cal C}_{L,k}(H_0)\).  Insert an
arbitrary Hermitian coefficient \(E\) in the first-order
Stein-endpoint derivative (equivalently, by the touching identity,
the characteristic/reversed-Horner and Schwarz recurrences).  This
is the unrestricted first-order version of L177's upper/lower
response match.  Using (11), the two endpoint paths cancel their
transpose-symmetric part.  The surviving circle-grade row is

\[
 {\cal C}_{L,k}(H_0)[F]
 ={4(4k-1)\over L^2}
 \sum_{\substack{a<b\\a+b=L-k}}
 (b-a){F_{ab}-F_{ba}\over2}.                        \tag{26}
\]

The disk symmetry \(JF=(JF)^*\) makes the row on \(a+b=L+k\) the
complex conjugate.  Equations (5), (25), and (26) prove (7) in
**every** order \(m\).  L182's cubic identity and the weighted-pair
quartic identity are its \(m=3,4\) special cases.

## 7. Cauchy absorption in every order

The transpose-skew projection is contractive:

\[
 \|F^-\|_F\le\|F\|_F,\qquad
 \|F^-\|_F^2=2\sum_{a<b}|F^-_{ab}|^2.               \tag{27}
\]

For \(r=L-k\),

\[
 \sum_{\substack{a<b\\a+b=L-k}}(b-a)^2
 =\binom r3.                                        \tag{28}
\]

The grade \(L-k\) and \(L+k\) rows are disjoint and have equal
projection magnitude.  Cauchy--Schwarz, (6), and (27)--(28) give

\[
\begin{aligned}
 D_{2m}
 &=4\|F\|_F^2\\
 &\ge16\sum_{k=3}^{L-3}
 {|P_{L,k}(F^-)|^2\over\binom{L-k}{3}},
\end{aligned}
\]

which is (8).  From L173,

\[
 b_{L,k}^{\rm flux}
 ={(4k-1)^2\over4L^4}\binom{L-k}{3},
\qquad b_{L,k}>b_{L,k}^{\rm flux}.                  \tag{29}
\]

Equations (7)--(9) show that the right side of (8) is the complete
flux-only normal gain.  If a response is nonzero, the strict
inequality in (29) leaves margin.  If every response is zero, the
nonzero base square in (6) leaves margin.  Thus the actual
first-residual Schur face is strict for every \(F\ne0\).

## 8. The nonlinear tube

Let \(\Gamma(H,y)\) be the analytic condition-square excess of the
rank-one Stein branch in the transported true-circular-normal fibre
\(y\).  L65/L173 give a uniformly negative \(y\)-Hessian near the
origin.  Equations (4), (23) give

\[
 \Gamma(H,0)=0,\qquad D_y\Gamma(H,0)=0
 \quad(H\in\Psi^{-1}(0)).                            \tag{30}
\]

Assume, for contradiction, that positive values of \(\Gamma\)
approach the origin.  Real-analytic curve selection supplies an
analytic arc \((H(s),y(s))\) of positive values.

If \(\Psi(H(s))\) first appears in order \(m\) and \(y\) first
appears in order \(r\), then:

* \(r<m\): the negative normal quadratic is the leading term;
* \(r>m\): the negative base square (6) is the leading term;
* \(r=m\): completing the normal square gives the strict face in
  Section 7.

All three contradict positivity.  If \(\Psi(H(s))\equiv0\), then
(30) and the negative fibre Hessian give strict descent unless
\(y\equiv0\), in which case the arc is exact equality.  This also
contradicts positive values.  Hence no positive arc exists, and
curve selection proves L189.

This proof is fixed-dimension; it does not claim dimension-uniform
neighbourhood constants.

## 9. Exact and independent regeneration

Run

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_equality_manifold.py \
  --minimum-length 3 \
  --maximum-differential-length 7 \
  --formal-length 6 \
  --formal-order 5 \
  --maximum-numerical-length 6 \
  --output \
  experiments/crabb_full_disk_equality_manifold_s70224.jsonl

PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_leading_residual.py \
  --output \
  experiments/crabb_full_disk_leading_residual_s70224.jsonl

PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_pair_blowup.py \
  --minimum-length 4 --maximum-length 8 \
  --output \
  experiments/crabb_full_disk_pair_blowup_s70224.jsonl
```

The first checker:

1. verifies (11), the \(J\)-Hermitian symmetry, and the exact real
   rank through length seven;
2. constructs the zero-Toeplitz equality graph exactly through order
   five and checks all residuals, endpoint deficits through order ten,
   and every true-circular-normal response through order five; and
3. independently solves the nonlinear residual equations, obtaining
   non-Toeplitz equality points whose condition square and
   characteristic-Blaschke norm square are both four and whose full
   pulled ambient derivative is at roundoff.

The second checker tests (6)--(8) on first residual orders
\(m=1,2,3,4,5\).  The \(m=5\) record is the adversarial weighted pair
whose fourth residual was canceled exactly.  Thus the test suite
specifically guards against mistaking an eighth-order face for a
uniform theorem.  The third checker separately regenerates every
reversal pair through length eight against a dense exact transverse
direction and verifies its fourth residual, endpoint factor, and all
normal projections.
