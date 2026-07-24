# The Hardy factor and all-size sixth face (L183--L184, 2026-07-24)

## 1. Results

This note closes the arbitrary-size sixth-order inequality left by
L181--L182.

First, let \(T\) be a strict finite-dimensional numerical-radius
contraction and let \(v\in\ker T\).  Its observability Gramian

\[
 P_v=\sum_{j\ge0}(T^*)^jvv^*T^j
\]

obeys the sharp kernel bound

\[
\boxed{P_v\preceq4\|v\|^2I.}                        \tag{1}
\]

This is L183.  It is a direct consequence of Berger's unitary
\(2\)-dilation and is recorded because its equality defect, not only
the inequality, is the missing full-disk factor.

For L176's recentered disk path, let

\[
 D_{6,L}(z)=-2[s^6]\{\lambda_+(s)-4\lambda_-(s)\}.
\]

There is an explicit complex skew-symmetric
\((L-1)\times(L-1)\) matrix \({\cal R}_L(z)\), homogeneous cubic in
\(z,\overline z\), such that

\[
\boxed{D_{6,L}(z)=8\|{\cal R}_L(z)\|_{\wedge^2}^2.}  \tag{2}
\]

If \(T_{L,k}(z)\) denotes the braces in L182, so that

\[
G^{(3)}_{L,k}={16(4k-1)\over L^2}T_{L,k},
\]

then the two grade-\(\pm k\) anti-diagonals of \({\cal R}_L\) satisfy

\[
\begin{aligned}
4T_{L,k}
 &=\sum_{\substack{1\le a<b\le L-1\\a+b=L-k}}
   (b-a)({\cal R}_L)_{ab},\\
4\overline {T_{L,k}}
 &=\sum_{\substack{1\le a<b\le L-1\\a+b=L+k}}
   (b-a)({\cal R}_L)_{ab}.                           \tag{3}
\end{aligned}
\]

Consequently

\[
\boxed{
D_{6,L}\ge
256\sum_{k=3}^{L-3}
 {|T_{L,k}|^2\over\binom{L-k}{3}}.}                 \tag{4}
\]

The right side is exactly the completed response gain obtained by
using only L173's flux curvature.  L173's actual curvature has an
additional positive null-lift term.  Hence

\[
\boxed{
D_{6,L}
-\sum_{k=3}^{L-3}{|G^{(3)}_{L,k}|^2\over4b_{L,k}}
\ge0}                                                \tag{5}
\]

in every size and for every complex direction.  This is L184.

For \(L=6\), (4) is exactly the formerly conjectural sharp inequality
\(P_6\ge9|C_6|^2\).  Thus L180's weaker constant was sufficient but
not sharp; the sharp complex statement now follows without a large
SOS.

The sixth-order kernel is exactly

\[
\boxed{{\cal R}_L(z)=0.}                             \tag{6}
\]

L185--L186 subsequently classify this kernel and prove its required
eighth-order fallback; see
`proof/crabb_full_disk_kernel_eighth.md`.  Equation (5) is still a
homogeneous face, not yet the nonlinear full circular-range tube.

## 2. Kernel observability from the \(2\)-dilation

Berger's theorem gives a Hilbert space \({\cal K}\), an isometry
\(J:{\cal H}\to{\cal K}\), and a unitary \(U\) such that

\[
 T^j=2J^*U^jJ,\qquad j\ge1.                          \tag{7}
\]

Since \(Tv=0\), for every nonzero integer \(r\),

\[
\langle U^rJv,Jv\rangle
={1\over2}\langle T^rv,v\rangle=0
\]

for \(r>0\), and the negative case follows by conjugation.  Thus
\(\{U^{-j}Jv:j\ge0\}\) is an orthogonal family, each vector having
norm \(\|v\|\).

For \(j\ge1\),

\[
\langle v,T^jx\rangle
=2\langle U^{-j}Jv,Jx\rangle,
\]

whereas the \(j=0\) coefficient is
\(\langle Jv,Jx\rangle\).  Bessel's inequality gives

\[
\begin{aligned}
\langle P_vx,x\rangle
&=|\langle Jv,Jx\rangle|^2
 +4\sum_{j\ge1}|\langle U^{-j}Jv,Jx\rangle|^2\\
&\le4\|v\|^2\|x\|^2.
\end{aligned}
\]

This proves (1).  A radial limit removes strict stability.

The equality defect is also explicit.  If \(\Pi_v\) is the
orthogonal projection onto the orbit in the preceding paragraph and
\(x\perp v\), then

\[
\boxed{
4\|v\|^2\|x\|^2-\langle P_vx,x\rangle
=4\|v\|^2\|(I-\Pi_v)Jx\|^2.}                        \tag{8}
\]

## 3. The explicit disk-chart Hardy model

Use coefficient gauge.  For an arbitrary positive \(L\times L\)
matrix \(H\), extended by a final zero row and column, put

\[
K=H+S^*HS,\qquad A=2K^{-1}HS,\qquad q=He_0=Ke_0,
\]

where \(S\) is the unweighted shift.  L122 gives
\(W(K^{1/2}AK^{-1/2})=\overline{\mathbb D}\).

Let \(h=H_{00}\).  The physical kernel vector

\[
v=K^{-1/2}q=K^{1/2}e_0
\]

has norm squared \(h\) and is killed by
\(K^{1/2}AK^{-1/2}\).  Applying (1) gives the global upper sandwich

\[
\boxed{M\preceq4hK,}                                \tag{9}
\]

where \(M-A^*MA=qq^*\).  Also

\[
Me_0=hKe_0,                                         \tag{10}
\]

so \(h\) is an exact generalized eigenvalue.  Near the Crabb point it
is the simple lower endpoint.

The disk factorization supplies a concrete version of \(J\):

\[
({\cal J}_Hx)(w)
=H^{1/2}(I-\overline wS)(I-\overline wA)^{-1}x,
\qquad |w|=1.                                       \tag{11}
\]

Indeed,

\[
2K-\overline wKA-wA^*K
=2(I-wS^*)H(I-\overline wS),                        \tag{12}
\]

and Fourier integration of the corresponding Poisson kernel gives

\[
\int_{\mathbb T}\|{\cal J}_Hx\|^2\,dm=x^*Kx,
\qquad
\int_{\mathbb T}w^j{\cal J}_H^*{\cal J}_H\,dm
={1\over2}KA^j.                                     \tag{13}
\]

Moreover

\[
{\cal J}_He_0=H^{1/2}e_0
\]

is constant.  Thus the projection in (8) removes, Fourier coefficient
by Fourier coefficient, the \(H^{1/2}e_0\) component.

In coefficient coordinates this projection is

\[
Q_H=I-e_0{e_0^*H\over h}.                           \tag{14}
\]

For every \(x\) with \(e_0^*Kx=0\), equations (8) and (11)--(14) give

\[
x^*(4hK-M)x
=4h\left\|
H^{1/2}Q_H(I-\overline wS)
(I-\overline wA)^{-1}x
\right\|_{L^2(\mathbb T)}^2.                        \tag{15}
\]

## 4. The cubic skew residual

Now use

\[
H(s)=I/2+sZ(z)+s^2B(W),\qquad W=z\wedge J\overline z.
\]

Write \(A(s)=\sum_{d\ge0}s^dA_d\) and expand the expression in
(15) at \(x=e_L\).  Put

\[
{\cal Y}(s,w)
=Q_{H(s)}(I-\overline wS)
 (I-\overline wA(s))^{-1}e_L.                       \tag{16}
\]

The coefficient recurrence is finite.  With

\[
C_0(w)=\sum_{j=0}^{L}\overline w^{\,j}A_0^j,
\]

the inverse coefficients are

\[
C_d=C_0\sum_{\ell=1}^d\overline wA_\ell C_{d-\ell}.
                                                               \tag{17}
\]

Expand \(H_{00}^{-1}\) in (14), insert L176's interval pulse in
\(A_2,A_3\), and compare the negative Fourier coefficients in
coordinates \(1,\ldots,L-1\).  Direct cancellation gives

\[
\begin{aligned}
[s]{\cal Y}&=[s^2]{\cal Y}=0,\\
[s^3]{\cal Y}
&=\sum_{a,b=1}^{L-1}
 ({\cal R}_L)_{ab}\,\overline w^{\,a}e_b,\\
{\cal R}_L^T&=-{\cal R}_L.                          \tag{18}
\end{aligned}
\]

Terms in the orbit coordinate \(e_0\) have already been removed, and
the constant \(e_L\) term is killed by \(H^{1/2}\).

The cancellation in (18) is the Hardy version of L176--L177.  Its
all-index check uses only (17):

1. the two occurrences of \(Z\) on opposite endpoint paths give the
   two terms of \(W\);
2. L176's interval part cancels the order-two negative coefficient;
3. its removed mean gives the denominator \(t+2\);
4. interchanging the Fourier index and coefficient index reverses the
   remaining endpoint path, giving skew symmetry.

At \(s=0\), changing \(x=e_L+u\) changes the residual only by constant
functions in coordinates \(1,\ldots,L-1\).  These are orthogonal to
the negative modes in (18).  Hence analytic maximization of the top
generalized eigenvalue cannot reduce the order-three residual.
Since

\[
h(0)=e_L^*K(0)e_L={1\over2},
\]

equation (15) yields

\[
[s^6]\{\lambda_+-4h\}
=-4\|{\cal R}_L\|_{\wedge^2}^2.                     \tag{19}
\]

Equations (10) and (19) prove (2).

## 5. Identification of the response rows

Continue the coefficient comparison in (17), now summing one
circle-grade anti-diagonal.  If \(T_{L,k}\) is L182's flux coordinate,
the two orientations give exactly (3).  Equivalently, the degree-three
response is the pair of weighted Fourier rows of the Hardy residual.

The squared row norm on either anti-diagonal is

\[
\sum_{\substack{1\le a<b\\a+b=L-k}}(b-a)^2
=\binom{L-k}{3}.                                    \tag{20}
\]

Reflection \(a\mapsto L-a\) gives the same norm on \(a+b=L+k\).
The two anti-diagonals are disjoint, as are the pairs belonging to
different \(k\).  Applying Cauchy--Schwarz to both equations in (3)
therefore gives

\[
\|{\cal R}_L\|_{\wedge^2}^2
\ge32\sum_{k=3}^{L-3}
 {|T_{L,k}|^2\over\binom{L-k}{3}}.                  \tag{21}
\]

Multiplication by eight proves (4).

## 6. The actual Schur face

Put \(r=L-k\).  L173 splits the exact positive mode curvature as

\[
\begin{aligned}
b_{L,k}
&=b^{\rm flux}_{L,k}+b^{\rm null}_{L,k},\\
b^{\rm flux}_{L,k}
&={(4k-1)^2\over4L^4}\binom r3,\\
b^{\rm null}_{L,k}
&={2k(k-1)(k-2)\over3L^4}
  \left(r+{1\over4}\right)^2>0.                    \tag{22}
\end{aligned}
\]

Using L182,

\[
{|G^{(3)}_{L,k}|^2\over4b^{\rm flux}_{L,k}}
=256{|T_{L,k}|^2\over\binom r3}.                    \tag{23}
\]

Equations (4), (22), and (23) prove (5).

If \({\cal R}_L=0\), equations (2)--(3) make the complete face zero.
Conversely, if \({\cal R}_L\ne0\), either a response row is nonzero,
in which case the positive null lift makes (22) strict, or all response
rows vanish, in which case (2) remains strictly positive.  This proves
the kernel statement (6).

On the terminal two-coefficient edge, (21) is equality.  Equation
(22) then recovers L178's strict actual ratio and explains why its
`169` term survives.

## 7. Exact regeneration

Run

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_sixth_hardy_factor.py \
  --symbolic-maximum-length 8 \
  --exact-maximum-length 10 \
  --terminal-maximum-length 12 \
  --output \
  experiments/crabb_full_disk_sixth_hardy_factor_s70224.jsonl
```

The checker independently builds (16)--(18).  It verifies the
unrestricted-symbol transfer identities through length eight, the
endpoint identity (2) on dense exact directions through length ten,
and sharp equality in (4) on the terminal edge through length twelve.
Clean reruns must be byte-identical.

L184 closes the all-size sixth-order algebraic face.  L185--L186
subsequently solve its cubic kernel and eighth-order face.  The live
task is now the singular nonlinear blow-up.
