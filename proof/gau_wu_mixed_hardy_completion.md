# The Gau--Wu mixed zero motion has an exact Hardy completion

> **Scope.**  This note eliminates every first-order Blaschke-zero
> velocity from L336 in model-space coordinates.  It leaves one
> purely physical quadratic inequality.  That inequality is still
> open, so this is not the arbitrary-degree Hessian theorem.

## 1. Inner tangents in one model-space coordinate

Let

\[
 \phi=zf,\qquad {\cal K}={\cal K}_\phi,\qquad
 S=S(\phi),\qquad p=f,\quad q=1,
\]

and let \(J_f\) denote the canonical conjugation on
\({\cal K}_f\).  Modulo the irrelevant constant phase tangent, every
degree-preserving finite-Blaschke tangent \(h=\dot f\) has the unique form

\[
\boxed{
 h=h_k:=\phi k-J_fk,\qquad k\in{\cal K}_f.}        \tag{1}
\]

Indeed, on the circle,

\[
 \overline f h_k=zk-\overline z\,\overline k
\]

is purely imaginary, so (1) is tangent to the inner manifold.
Conversely, if \(\overline f h\) is purely imaginary, its strictly
positive Fourier coefficients determine \(k\); the remaining
constant is precisely the phase tangent.

For a finite Blaschke product

\[
 f(z)=\prod_{j=1}^d {z-a_j\over1-\overline a_jz}
\]

and zero velocities \(\dot a_j=v_j\), differentiation gives the
concrete coordinate

\[
\boxed{
 k(z)=\sum_{j=1}^d{\overline v_j\over
                         1-\overline a_jz}.}      \tag{2}
\]

Thus, at distinct zeros, the velocity-to-\(k\) map is an isomorphism
onto \({\cal K}_f\).  Repeated zeros use the usual confluent kernel
derivatives.

## 2. Endpoint action

Since \(\phi(S)=0\), model functional calculus and the two canonical
conjugations give

\[
\boxed{
 h_k(S)q=-J_fk,\qquad h_k(S)^*p=-zk.}              \tag{3}
\]

In particular both vectors have norm \(\|k\|\).  L337 now becomes

\[
 Q(0,k):=r_2(0,k)+\ell_2(0,k)=2\|k\|^2.           \tag{4}
\]

For the kernel coordinate (2),

\[
 \|k\|^2
 =\sum_{i,j=1}^d
 {v_i\overline v_j\over1-a_i\overline a_j}.       \tag{5}
\]

This is the exact Szegő Gram that replaces the numerical zero block
in L334.

## 3. The mixed functional

Let the conformally normalized operator jet be

\[
 T_\varepsilon=S+\varepsilon C+\varepsilon^2D+o(\varepsilon^2).
\]

With the inner function fixed at \(f\), put

\[
\begin{aligned}
 U_C&=Df(S)[C]q,\\
 V_C&=Df(S)[C]^*p,\\
 Z_C&=D\phi(S)[C].
\end{aligned}                                     \tag{6}
\]

The mixed second transition is

\[
 \Lambda_C(k)
 :=\langle p,Dh_k(S)[C]q\rangle
 =\langle p,Z_Ck\rangle
  -\langle p,D(J_fk)(S)[C]q\rangle.               \tag{7}
\]

Expanding the two endpoint defects in L336 and using (3) gives the
exact joint quadratic form

\[
\boxed{
\begin{aligned}
 Q(C,k)
 ={}&Q(C,0)+2\|k\|^2+{\cal M}_C(k),\\
 {\cal M}_C(k)
 ={}&-4\operatorname{Re}\Lambda_C(k)\\
 &+2\operatorname{Re}\{
   \langle U_C,J_fk\rangle+
   \langle V_C,zk\rangle\}.
\end{aligned}}                                    \tag{8}
\]

Only first-order operator data enter \({\cal M}_C\).  A second zero
acceleration cannot enter because the zero gradient vanishes at the
extremal.

## 4. Exact elimination of every zero velocity

Equip \({\cal K}_f\) with its underlying real Hilbert inner product.
There is a unique Riesz vector \(\kappa_C\in{\cal K}_f\) satisfying

\[
\boxed{
\begin{aligned}
 \operatorname{Re}\langle\kappa_C,k\rangle
 ={}&\operatorname{Re}\Lambda_C(k)\\
 &-\frac12\operatorname{Re}\{
   \langle U_C,J_fk\rangle+
   \langle V_C,zk\rangle\}
\end{aligned}}                                    \tag{9}
\]

for every \(k\in{\cal K}_f\).  Equations (8)--(9) complete the
square:

\[
\boxed{
 Q(C,k)
 =2\|k-\kappa_C\|^2+
   \left\{Q(C,0)-2\|\kappa_C\|^2\right\}.}         \tag{10}
\]

Consequently L336's sufficient endpoint-flux gate is equivalent to
the purely physical inequality

\[
\boxed{
 Q(C,0)\geq2\|\kappa_C\|^2
 \quad\text{for every physical normal direction }C.} \tag{11}
\]

The optimizing zero tangent is exactly \(k=\kappa_C\).  At distinct
zeros, (9) is just one Szegő-Gram solve using (5); no opaque
\(2d\)-by-\(2d\) real Schur complement remains.

## 5. What remains

Equation (11) is now the sole content of the arbitrary finite
Gau--Wu Hessian gate.  The next step is to insert the L335 support
Gram into \(Q(C,0)\), express \(\kappa_C\) through the contraction and
characteristic residuals

\[
 S^*C+C^*S,\qquad D\phi(S)[C],
\]

and prove the Riesz response has norm at most the retained physical
energy.  Pure zero positivity cannot prove (11); its role was to
make the completion exact.

## 6. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_mixed_hardy_completion.py \
  --output experiments/gau_wu_mixed_hardy_completion_s70224.jsonl
```

The checker compares (2)--(5) with independent matrix and Fourier
evaluations, and compares (7)--(10) against the complete joint
Blaschke jet using block-matrix Fréchet calculus for \(Dh_k(S)[C]\).
The 12-record dataset has SHA-256
`e600a466d6ad1fd78e692f9c370a8a73b6e6fe2101e582f98e232991de83d67c`;
the largest residual is \(1.21\cdot10^{-13}\).
