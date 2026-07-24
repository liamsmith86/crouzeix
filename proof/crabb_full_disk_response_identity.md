# The full-disk ambient response identity (L177, 2026-07-24)

## 1. Result

Retain L176's notation.  Let \(g_2(z)\) be the quadratic
Toeplitz-disk coefficient of the ambient gradient of the normalized
characteristic-Blaschke square at the Crabb point.  Let

\[
 E_z=D_{I/2}X[B(z\wedge J\overline z)]
\]

be L176's explicit physical disk-tangent correction.  Then

\[
\boxed{g_2(z)=2{\cal C}_pE_z}                         \tag{1}
\]

as real covectors on the complete ambient matrix space.  Here
\({\cal C}_p=-e_p\) is L65's positive curvature.

L176 gives

\[
\langle E_z,{\cal C}_pE_z\rangle=32{\cal Q}(z).
\]

Consequently L175's complete homogeneous full-H face is

\[
\begin{aligned}
{\cal F}_4(z,Y)
&=-32{\cal Q}(z)+g_2(z)[Y]
  -\langle Y,{\cal C}_pY\rangle\\
&=-\langle Y-E_z,{\cal C}_p(Y-E_z)\rangle\le0.        \tag{2}
\end{aligned}
\]

Thus the full strong completion cancels, but never exceeds, the
Toeplitz quartic.  Equality occurs precisely at the L176 correction
modulo \(\ker{\cal C}_p\).  This proves the leading full-disk
recentring conjecture L175.

Equation (2) is a homogeneous apex theorem.  The nonlinear
full-circular-range tube and its elliptic/compact merger remain
downstream.

## 2. Intrinsic Plücker coordinates

Put \(n=L-1\) and

\[
 W=z\wedge J\overline z.
\]

Fix a paired L65 circle mode \(1\le k\le L-3\), and put

\[
 t=L-2-k,\qquad m=t+2=L-k.
\]

For \(0\le i<t-i\), let

\[
 \chi_i={\bf1}_{\{i+1,\ldots,t-i\}}\in{\mathbb R}^m.
\]

Use L65's vectors

\[
\begin{aligned}
q_0&=(1/\sqrt2,1,\ldots,1,1/\sqrt2)^T,\\
q_k&=((k+2)/\sqrt2,1,\ldots,1,(k+2)/\sqrt2)^T.
\end{aligned}
\]

The two elementary sums are

\[
q_k^Tq_0=L,\qquad q_k^T\chi_i=t-2i.                  \tag{3}
\]

Define

\[
\boxed{
x_i=16\left(\chi_i-{t-2i\over L}q_0\right).}          \tag{4}
\]

L176 equations (9)--(13) say exactly that if \(s_i\) is the reduced
L65 vector of its interval-pulse correction, then

\[
{\cal K}_{m,k}s_i=x_i.                               \tag{5}
\]

## 3. Complete paired response formula

For an arbitrary physical ambient direction \(Y\), let
\(\tau_k(Y)\in{\mathbb C}^m\) be L65 equations (12)--(14)'s reduced
vector.  The complete quadratic response in mode \(k\) is

\[
\boxed{
g_2^{(k)}(z)[Y]
=2\operatorname{Re}
\left\langle
\sum_{0\le i<t-i}W_{i,t-i}x_i,\,
\tau_k(Y)
\right\rangle.}                                      \tag{6}
\]

Every nonzero circle mode \(k>L-3\) has zero quadratic response.

Writing

\[
X=\sum_{0\le i<t-i}W_{i,t-i}x_i
\]

and using L65's raw variables \(u_j,v_j\), equation (6) is the explicit
coefficient table

\[
\begin{aligned}
\frac12g_2^{(k)}[Y]
=\operatorname{Re}\{&
\overline X_0(u_0+\eta_kv_0+\sqrt2v_1)\\
&+\sum_{a=1}^{m-2}\overline X_a(u_a+v_{a+1})\\
&+\overline X_{m-1}
 (u_{m-1}+\sqrt2v_m+\eta_kv_{m+1})\},
                                                               \tag{6a}
\end{aligned}
\]

where \(\eta_1=1/\sqrt2\) and \(\eta_k=1\) otherwise.

To prove (6), work first in L123's coefficient gauge.  Along the
Toeplitz disk ray \(H=I/2+sZ(z)\), write

\[
\begin{aligned}
K^{-1}&=N_0+sN_1+s^2N_2+O(s^3),\\
A&=A_0+sA_1+s^2A_2+O(s^3),
\end{aligned}
\]

where

\[
N_1=-N_0K_1N_0,\qquad
N_2=N_0K_1N_0K_1N_0.                                \tag{7}
\]

For a frozen raw direction \(Y\):

1. expand its support quotient through degree two;
2. subtract the first Schwarz/Riemann correction;
3. run paired Horner recurrences for the characteristic polynomial
   \(g\), its reversal \(g^\sharp\), and their Frechet derivatives;
4. pair the result with the Crabb endpoint defect and endpoint vector.

Circle character leaves only products

\[
W_{i,t-i}
=z_i\overline{z_{n-1-(t-i)}}
 -z_{t-i}\overline{z_{n-1-i}}.                       \tag{8}
\]

Retaining the raw L65 variables instead of pairing immediately with
L173's support Riesz representative gives the following coefficient
vector:

\[
16\chi_i-16{t-2i\over L}q_0.                         \tag{9}
\]

The interval term counts the two characteristic endpoint paths whose
indices lie between \(i+1\) and \(t-i\).  The second term is the zero
Fourier mode of the support quotient; its denominator is the Crabb
boundary mass \(L\), and its numerator is the interval length
\(t-2i\).  All remaining numerator, reversed-polynomial, and moving
metric paths cancel in adjacent pairs.  Thus (9) is (4), and L65's
aggregator gives the raw table (6a), proving (6).

This is the unprojected version of L173's response calculation.
Pairing (6) with L173's support representatives reproduces L173
equation (2), including the first nonzero `L=6,k=3` value.

## 4. Terminal grade-zero response

The terminal exterior anti-diagonal \(t=n-1\) produces L176's real
diagonal correction

\[
d=\operatorname{diag}B(W)\in{\mathbb R}^L.
\]

Let \(\beta=(\sqrt2,1,\ldots,1,\sqrt2)^T\) be the Crabb weights.
The derivative of L122's diagonal-H disk chart maps \(d\) to physical
weight variations by

\[
({\cal A}_Ld)_j
=\beta_j\left[
2d_j-\frac12\left(
{d_{j-1}+d_j\over\kappa_j}
+{d_j+d_{j+1}\over\kappa_{j+1}}
\right)\right],                                      \tag{10}
\]

where missing \(d_{-1},d_L\) are zero,
\(\kappa_0=\kappa_L=1/2\), and
\(\kappa_j=1\) otherwise.

Let \(\Delta\) be the first-difference matrix.  There is a unique real
covector \(q(W)\) satisfying

\[
\boxed{
{\cal A}_L^Tq=4\Delta^T\Delta d,\qquad
q^T\beta=0.}                                         \tag{11}
\]

For a physical grade-zero direction with first-superdiagonal weight
variation \(y\),

\[
\boxed{g_2^{(0)}(z)[Y]=2q(W)^T\operatorname{Re}y.}    \tag{12}
\]

The degree-two characteristic/Riemann recurrence from Section 3,
now in grade zero, gives the left side of (11) entry by entry.
The zero Fourier coefficient gives the radial constraint
\(q^T\beta=0\).  Since \({\cal A}_L\) has kernel
\({\mathbb R}{\bf1}\), these identities determine \(q\) uniquely and
prove (12).  Imaginary weight variations are infinitesimal diagonal
unitary orbits and have zero response.

## 5. Identification with L65 curvature

For paired modes, L176 gives the reduced vector

\[
s(W)=\sum_{0\le i<t-i}W_{i,t-i}s_i.
\]

Using (5), the polarization of L65 equation (18) is

\[
\begin{aligned}
2\langle E_z,{\cal C}_pY\rangle
&=2\operatorname{Re}
  \langle{\cal K}_{m,k}s(W),\tau_k(Y)\rangle\\
&=2\operatorname{Re}
  \left\langle\sum_iW_{i,t-i}x_i,\tau_k(Y)\right\rangle.
\end{aligned}
\]

This is exactly (6).

In grade zero, L176 equation (14) gives

\[
\langle E(d),{\cal C}_pE(d')\rangle
=4\langle\Delta d,\Delta d'\rangle.                  \tag{13}
\]

Every real physical weight variation is a disk-chart tangent
\({\cal A}_Ld'\) plus a radial multiple of \(\beta\).  Equations
(11)--(13) show

\[
\langle E(d),{\cal C}_pY\rangle
=q^T\operatorname{Re}y.
\]

This is (12).  Section 3's character selection gives zero response
on every inactive mode.  Hence (1) holds on every ambient direction.

Finally, A111/L157 make the optimized-upper minus prepared-dual
ambient gradient \(O({\cal Q})\).  Since \(g_2\) is quadratic in \(z\)
whereas \({\cal Q}\) is quartic, the optimized L118 upper envelope has
the same response (1).  Substitution of L176's energy proves (2).

## 6. Exact regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_full_disk_response_exact.py \
  --minimum-length 3 --maximum-length 6 \
  --output experiments/crabb_full_disk_response_exact_s70224.jsonl
```

The checker reuses L173's independent exact characteristic/Riemann
series engine.  It compares (6) and (12) with every real and imaginary
matrix unit, verifies zero on all inactive modes, and uses exact SymPy
arithmetic.  A second run is byte-identical.
