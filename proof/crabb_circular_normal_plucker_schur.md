# The circular-normal Plücker--Schur face (L173, 2026-07-24)

## 1. Result and corrected mechanism

Put `p=L+1` and write the nonconstant Toeplitz disk coordinate as

\[
h=(h_1,\ldots,h_{L-1}),\qquad
{\cal Q}(h)=\|h\|^4-|h^TJh|^2 .
\]

Let \(R_k^{\rm R},R_k^{\rm I}\) be the two raw support-Riesz
representatives of circular-normal Fourier mode \(k\), as in L115,
with \(3\le k\le L+1\).  The quadratic disk-to-normal gradient does
**not** vanish in every size.  Instead it factors through

\[
W=h\wedge J\overline h.                              \tag{1}
\]

Its complex response is

\[
\boxed{
G_{L,k}(h)=
{8(4k-1)\over L^2}
\sum_{\substack{1\le i<j\le L-1\\i+j=L+k}}
(j-i)\{h_i\overline{h_{L-j}}
       -h_j\overline{h_{L-i}}\}.
}                                                     \tag{2}
\]

The real and imaginary parts pair with \(R_k^{\rm R}\) and
\(R_k^{\rm I}\), respectively.  The sum is empty unless
\(3\le k\le L-3\).  Thus the cancellation seen through `L=5` was
real but low-dimensional; the first surviving term is `L=6,k=3`.

Despite this correction, the exact Schur complement is favorable:

\[
\boxed{\lambda_{\max}<512.}                          \tag{3}
\]

L152 and L157 give the same leading disk deficit for the optimized
upper envelope and the prepared characteristic-Blaschke lower
certificate:

\[
-32{\cal Q}(h)=-128\|\mathfrak p\|^2,                \tag{4}
\]

where \(\mathfrak p\) is the real Plücker vector below, (3) says that
completing every circular-normal square consumes strictly less than
the available disk deficit.  This proves the sharp homogeneous apex
face required by the full circular-range tubular lift.

It does **not** by itself prove the uniform nonlinear tube away from
the apex; that analytic patching step is recorded in Section 6.

## 2. Quadratic response calculation

At the Crabb point, in coefficient gauge,

\[
K_0=\operatorname{diag}(1/2,1,\ldots,1,1/2).
\]

The real mode-\(k\) support Riesz representative is supported on
matrix grades \(k+1\) and \(1-k\), with entries

\[
(R_k^{\rm R})_{ij}={(K_0)_{ii}\over2L};
                                                               \tag{5}
\]

the imaginary representative has the same magnitudes and opposite
signs on the two grades.

Insert \(H=I/2+xZ(h)\) into L122's exact disk model

\[
K=\widehat H+S^*\widehat H S,\qquad
A=2K^{-1}\widehat H S.
\]

Expand \(K^{-1}\) through degree two by

\[
N_0=K_0^{-1},\quad N_1=-N_0K_1N_0,\quad
N_2=N_0K_1N_0K_1N_0.                                \tag{6}
\]

Use these coefficients in the characteristic-polynomial Horner
recurrence, its reversed polynomial, and the first Schwarz/Riemann
correction.  Pair the resulting frozen Blaschke-norm derivative with
(5).  Circle-grade selection leaves only index pairs

\[
i+j=L+k.
\]

For such a pair, the two endpoint paths and the reversed-polynomial
path combine to

\[
{8(4k-1)\over L^2}(j-i)
\{h_i\overline{h_{L-j}}-h_j\overline{h_{L-i}}\}.     \tag{7}
\]

Summing (7) proves (2).  This is a finite coefficient calculation:
there is no root selection and no numerical premise.  The exact
series checker independently performs (6), the support pull, both
Horner recurrences, and the moving metric pairing before comparing
with (2).

Equivalently, putting \(r=L-k\) and collecting both ends of every
pair, the same calculation is the weighted lag-\(k\) correlation

\[
G_{L,k}(h)
={8(4k-1)\over L^2}
\sum_{a=1}^{r-1}(r-2a)h_{k+a}\overline{h_a}.         \tag{7a}
\]

This form makes the all-index path count in (7) directly checkable:
the coefficient decreases by two when both Toeplitz legs move one
step along the Crabb chain, and reversal changes its sign at the
opposite endpoint.

For its first nonzero deterministic audit, (2) gives exactly

\[
L=6,\quad k=3,\qquad G_{6,3}=-{253\over45000},
\]

the coefficient that falsified the proposed universal cancellation.

## 3. Plücker norm and mode orthogonality

Choose a unitary Takagi factor \(U\) with

\[
U^TU=J,\qquad y=Uh=a+ib.
\]

Then

\[
U(J\overline h)=\overline y,\qquad
y\wedge\overline y=-2i\mathfrak p,
\]

where

\[
\mathfrak p_{ij}=a_i b_j-a_j b_i.
\]

Consequently

\[
\boxed{{\cal Q}(h)=\|W\|^2=4\sum_{i<j}\mathfrak p_{ij}^2.}       \tag{8}
\]

Formula (2) occupies only exterior anti-diagonal \(i+j=L+k\).
Different modes are therefore orthogonal.  The real and imaginary
polarizations of one nonzero circle character are orthogonal and
have equal norm.  With \(r=L-k\),

\[
\sum_{\substack{i+j=L+k\\i<j}}(j-i)^2
=\binom r3,
\]

so either real Plücker row \(m_{L,k}\) has

\[
\boxed{
\|m_{L,k}\|^2
={128(4k-1)^2\over L^4}\binom r3 .
}                                                     \tag{9}
\]

This also proves directly that only \(r\ge3\), equivalently
\(k\le L-3\), can couple.

## 4. Exact L65 curvature

For an active mode put \(r=L-k\ge3\).  In L65's reduced coordinates,
the representative (5) becomes

\[
u={1\over L}(a,1,\ldots,1,a),\qquad
a={3+2\sqrt2\over4}.                                 \tag{10}
\]

Its positive curvature is

\[
b_{L,k}=u^TK_{r,k}u,\qquad
K_{r,k}=(B_r+c_kvv^T)^{-1},                          \tag{11}
\]

where

\[
\begin{aligned}
B_r&={I-H_r\over4},&
c_k&={3\over2k(k-1)(k-2)},\\
v&=((k+2)/\sqrt2,1,\ldots,1,(k+2)/\sqrt2)^T .
\end{aligned}
\]

The null vector of \(B_r\) is

\[
q=(1/\sqrt2,1,\ldots,1,1/\sqrt2)^T,\qquad q^Tv=L.
\]

Put

\[
\beta={q^Tu\over q^Tv}
={r-1+3\sqrt2/4\over L^2},\qquad w=u-\beta v.
\]

The rank-one lift of a one-dimensional nullspace gives

\[
b_{L,k}=w^TB_r^+w+{\beta^2\over c_k}.                \tag{12}
\]

Moreover \(B_r=\frac18D^TD\), where \(D\) is the weighted path
difference in the \(q\)-coordinates.  Its cumulative flux is

\[
\sum_{i=0}^j q_iw_i
={(2j-r+2)(4k+4-3\sqrt2)\over8L^2}.
\]

Using

\[
\sum_{j=0}^{r-2}(2j-r+2)^2={r(r-1)(r-2)\over3}
\]

in (12) yields

\[
\boxed{
b_{L,k}={1\over L^4}\left[
{(4k+4-3\sqrt2)^2r(r-1)(r-2)\over24}
+{2k(k-1)(k-2)\over3}
 \left(r-1+{3\sqrt2\over4}\right)^2
\right].
}                                                     \tag{13}
\]

## 5. Strict Schur absorption

The homogeneous joint face in one real polarization is

\[
-b_{L,k}t^2+(m_{L,k}\cdot\mathfrak p)t.
\]

Completing the square adds
\(\frac14(m\cdot\mathfrak p)^2/b\).  Because all mode/polarization
rows are orthogonal, the eigenvalues of the complete gain matrix are
exactly

\[
{\|m_{L,k}\|^2\over b_{L,k}}.
\]

Now

\[
4k+4-3\sqrt2>4k-1
\quad\Longleftrightarrow\quad 5>3\sqrt2.
\]

The first, flux term of (13) alone therefore gives

\[
b_{L,k}>
{(4k-1)^2r(r-1)(r-2)\over24L^4}
={\|m_{L,k}\|^2\over512}.                            \tag{14}
\]

This proves (3).  Inactive modes have zero coupling and strictly
positive L65 curvature.  Combining (4), (9), and (14), the completed
normal gain is strictly smaller than \(32{\cal Q}\).  Hence the
leading disk/circular-normal form is negative away from the joint
equality set.

L157 also makes the optimized-upper/prepared-dual gap \(O({\cal
Q}^2)\), and A111's gradient comparison is \(O({\cal Q})\).
With a circular normal assigned its forced apex weight two, neither
term changes this degree-four joint face.  Thus (14) applies to the
leading optimized certificate needed for the upper bound, not only
to an unrelated lower test.

## 6. Remaining analytic lift

L173 removes the sharp algebraic obstruction at the stratified Crabb
apex.  To obtain the full circular-range tube one still has to:

1. use L162's exact ambient stationarity on every positive
   phase-palindromic equality anchor;
2. transport L65's coercive normal frame analytically over the local
   Lewis--Overton disk manifold;
3. apply a compact weighted blow-up so the strict margin (14)
   persists between the apex and positive anchors; and
4. merge that estimate with L160's grade-one elliptic face and
   L163/L172's higher-grade zero rows.

No new coefficient identity is suggested by this list; it is the
uniform tubular/Morse--Bott patching step.

## 7. Exact and adversarial regeneration

Run

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_circular_normal_quadratic_exact.py \
  --minimum-length 3 --maximum-length 7 \
  --output \
  experiments/crabb_circular_normal_quadratic_exact_s70223.jsonl

PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_circular_normal_plucker_schur.py \
  --minimum-length 3 --maximum-length 12 \
  --resolution 1024 \
  --output \
  experiments/crabb_circular_normal_plucker_schur_s70223.jsonl
```

The first checker uses exact Gaussian-rational series and verifies
(2), including its first nonzero case.  The second independently
extracts the quadratic response by symmetric Richardson
extrapolation, evaluates L65's full finite program on every raw mode
representative, and compares with (9), (13), orthogonality, and the
strict threshold (14).  The earlier finite nonlinear Schur probe is
retained as a separate pre-asymptotic guard.
