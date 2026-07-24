# Jensen scalarization of repeated Crabb circular normals

## 1. Result (L199, 2026-07-24)

Fix a Crabb length \(L\) and copy multiplicity \(m\).  Combine:

- L193--L197's block-disk Hardy residual and finite copy flag;
- L61's first-support Jensen reduction at a repeated Crabb block; and
- L173/L188's scalar true-circular-normal absorption.

At the first repeated residual face, and at any later face for which one
has independently exposed a freshly recentered raw Hardy residual, one has
the following alternative.

1. The first support-compression Jensen gap is positive.  Then L61 gives
   strict descent one order before the disk residual/normal Schur face.
2. The Jensen gap is zero.  On every maximal common-winner copy space,
   every true circular-normal coefficient is scalar and has no cross block
   to the losing space.  On that winner, the scalar L188 inequality applied
   state by state gives the full copy-matrix endpoint inequality.

Consequently true circular normals are absorbed on the first raw repeated
Hardy face.  If that joint face has a kernel, every circular normal vanishes
there and the kernel reduces every leading residual block.

This is an independent support-rank audit of L198 and identifies the exact
zero-Jensen mechanism which avoids a false noncommutative tensorization.  It
does **not** yet close the complete L197 flag: after Schur orthogonalization,
a later residual need not retain the raw reflection formula below, and the
analytic normal critical graph can mix earlier active ranges into the later
quotient.  The flag lift and the elliptic soft mode both remain.

## 2. The first Jensen gate

For a repeated Crabb block, let

\[
 B_E(\theta)
 =V(\theta)^*\operatorname{Re}(e^{-i\theta}E)V(\theta)
 \in M_m^{\rm h}                                      \tag{1}
\]

be the first perturbation compressed to the repeated top-support
eigenspace.  L61 proves

\[
 J_E=2L\,\mathop{\rm mean}_\theta\lambda_{\max}B_E(\theta)
     -\lambda_{\max}M_E\ge0,                         \tag{2}
\]

and the first similarity-square derivative is \(-4J_E\).  Equality holds
exactly when \(B_E(\theta)\) has one common maximizing copy vector for
every \(\theta\).

Promote all common maximizing vectors to a maximal winner space \(W\).
Then, pointwise,

\[
 P_WB_E(\theta)P_W=\lambda(\theta)P_W,\qquad
 P_{W^\perp}B_E(\theta)P_W=0.                        \tag{3}
\]

The losing mean compression is strict; if its gap collapses, its common
kernel is promoted into \(W\), exactly as in L82/L86.

## 3. Fourier independence scalarizes the true normals

L115's true circular normals are the real and imaginary support-Riesz
representatives in distinct Fourier modes \(3,\ldots,L+1\).  In the
transported tubular coordinate, write their copy coefficients as
\(Y_{k,\mathrm R},Y_{k,\mathrm I}\).  The block-disk tangent has zero
support modes \(2,\ldots,L+1\), while affine and elliptic coordinates
occupy the lower, disjoint modes.

Take the Fourier coefficients of (3).  Independence of
\(\cos(k\theta),\sin(k\theta)\) gives, coefficient by coefficient,

\[
\begin{aligned}
 P_WY_{k,\bullet}P_W&=y_{k,\bullet}P_W,\\
 P_{W^\perp}Y_{k,\bullet}P_W&=0,
 \qquad \bullet\in\{\mathrm R,\mathrm I\}.            \tag{4}
\end{aligned}
\]

Hermiticity gives the opposite cross block.  Thus the complex normal
coordinate on \(W\) is one scalar
\(y_k=y_{k,\mathrm R}+iy_{k,\mathrm I}\), not an arbitrary copy matrix.

This observation is load-bearing.  A blind replacement of the scalar
normal amplitude by a noncommuting matrix is neither required nor
justified: such a coefficient is already outside the zero-Jensen face and
is strict at the earlier order (2).

## 4. The block residual and its scalar states

Reindex the first Hardy residual from L195 as

\[
 F=(F_{ab})_{1\le a,b\le L-1},\qquad F_{ab}\in M_{\dim W}. \tag{5}
\]

Because its first nonzero coefficient is the differential of the
Hermitian inverse-Gram residual at the Crabb point,

\[
 \boxed{F_{ab}=F_{L-b,L-a}^*.}                       \tag{6}
\]

Equivalently, \(JF\) is block Hermitian.  Define the block-valued L188
response

\[
 G_{L,k}(F)
 ={4(4k-1)\over L^2}
 \sum_{\substack{a<b\\a+b=L-k}}
 (b-a){F_{ab}-F_{ba}\over2},
 \quad 3\le k\le L-3.                               \tag{7}
\]

For a unit copy vector \(v\), put

\[
 f^v_{ab}=\langle F_{ab}v,v\rangle.                  \tag{8}
\]

Equation (6) says \(Jf^v=(Jf^v)^*\), and linearity gives

\[
 \langle G_{L,k}(F)v,v\rangle=G_{L,k}(f^v).          \tag{9}
\]

Moreover,

\[
 \sum_{a,b}|f^v_{ab}|^2
 \le \sum_{a,b}\|F_{ab}v\|^2
 =\left\langle\sum_{a,b}F_{ab}^*F_{ab}v,v\right\rangle.
                                                               \tag{10}
\]

## 5. Statewise L188 proves the matrix inequality

Let \(b_{L,k}>0\) be L173's curvature in complex mode \(k\).  On a
zero-Jensen winner, the leading disk/normal endpoint is

\[
\begin{aligned}
 {\cal Q}(F,y)={}&-4\sum_{a,b}F_{ab}^*F_{ab}\\
 &+\sum_{k=3}^{L-3}
 \left\{-b_{L,k}|y_k|^2I+
 {\,\overline y_kG_{L,k}(F)
       +y_kG_{L,k}(F)^*\over2}\right\},              \tag{11}
\end{aligned}
\]

plus strictly negative squares from response-inactive true-normal modes.
The cross row in (11) is just the scalar L188 coefficient calculation:
the copy factor is \(y_kI_W\), so no ordering ambiguity occurs.

Compress (11) to a unit vector \(v\).  By (9)--(10),

\[
\begin{aligned}
 \langle{\cal Q}(F,y)v,v\rangle
 \le{}&-4\sum_{a,b}|f^v_{ab}|^2\\
 &+\sum_k\{-b_{L,k}|y_k|^2+
          \operatorname{Re}(\overline y_kG_{L,k}(f^v))\}.
                                                               \tag{12}
\end{aligned}
\]

The right side is exactly L188's scalar joint face, including L173's
strict positive null lift.  Therefore

\[
 \boxed{{\cal Q}(F,y)\preceq0.}                      \tag{13}
\]

This proof is a vector-state amplification, not a claim that the scalar
normal Hessian is completely positive for arbitrary matrix amplitudes.

## 6. Exact kernel and compatibility with L197

L188 is strict for every nonzero scalar residual after completing the true
normal squares.  Equality in (10) also requires every \(F_{ab}v\) to be
parallel to \(v\).  Hence equality in (12) has the sharper consequences

\[
 {\cal Q}(F,y)v=0
 \Longrightarrow
 y_k=0\ (\hbox{all }k),\qquad F_{ab}v=0\ (\hbox{all }a,b). \tag{14}
\]

Using (6) on the reflected blocks also gives

\[
 F_{ab}^*v=0\qquad(\hbox{all }a,b).                  \tag{15}
\]

Thus the first joint kernel reduces every raw leading residual block and
contains no surviving true-normal coefficient.  L196 promotes precisely
this first reducing kernel.

At a later L197 member, the same statewise argument applies if the quotient
has first been identified with the copy-preserving compression of a
freshly recentered raw residual satisfying (6).  That identification is
not automatic: L197's residual is least-squares orthogonalized against
earlier active Hardy ranges, while the transported analytic normal graph
may contain response terms through those ranges.  Proving that these terms
vanish or are already paid for is the remaining circular-normal flag debt.
The separate elliptic soft mode also remains.

## 7. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_circular_jensen.py \
  --output \
  experiments/repeated_crabb_circular_jensen_s70224.jsonl
```

The checker:

1. constructs generic block residuals satisfying (6);
2. tests (12) on coordinate and dense copy states;
3. constructs the flux-Cauchy equality family whose ratio approaches the
   sharp scalar threshold as \(L\) grows;
4. verifies Fourier coefficientwise winner scalarization and zero
   winner/loser cross blocks; and
5. prescribes nontrivial residual kernels and verifies that every response
   and its adjoint annihilate them.

The standard dataset covers lengths \(6,\ldots,14\), multiplicities
\(2,3,4\), generic/sharp/kernel families, and all active response modes.
Equations (1)--(15), rather than the floating audit, prove the result.

The tracked dataset regenerates byte for byte with SHA-256

```text
9ece0fde6a29e4d5d8a58038e33e30b3c8d857713052303496a2916cd9210651
```
