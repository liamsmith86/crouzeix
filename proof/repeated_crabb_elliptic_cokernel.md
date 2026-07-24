# The reduced cokernel equation for the repeated elliptic face

## 1. Result (L204, 2026-07-24)

L203 reduces the grade-one elliptic problem at an L193 equality
anchor to

\[
 X-T^*XT=F_2+VC_2^*+C_2V^*,                       \tag{1}
\]

with \(TV=0\), \(V^*V=I_m\), and desired endpoints

\[
 V^*XV=0,\qquad W^*XW=-ZZ^*.                       \tag{2}
\]

Here \(Z=4B_1\) is L202's defect-to-defect normal corner.  The apparent
large linear system in \(X,C_2\) has an exact smaller form.

Let

\[
 {\cal G}_T(H)=\sum_{j\ge0}(T^*)^jHT^j,             \tag{3}
\]

the stable Stein inverse, and put

\[
\begin{aligned}
 K&=V^*F_2V,\\
 C_\parallel&=-\frac12VK,\\
 F_0&=F_2-VKV^*,\\
 X_0&={\cal G}_T(F_0).                              \tag{4}
\end{aligned}
\]

Then \(V^*X_0V=0\).  Every remaining effective freedom is one column
\(C\) satisfying \(V^*C=0\), through the real-linear endpoint map

\[
 \boxed{
 {\cal M}_T(C)
 =W^*{\cal G}_T(VC^*+CV^*)W.}                      \tag{5}
\]

Consequently (1)--(2) are solvable exactly when

\[
 \boxed{
 {\cal M}_T(C)=D_T,\qquad
 D_T=-ZZ^*-W^*X_0W.}                               \tag{6}
\]

This is the finite Lyapunov--Schmidt equation hidden in L203.

Its adjoint is also explicit.  For \(Y=Y^*\in M_m\), let

\[
 {\cal Z}_Y
 =\sum_{j\ge0}T^jWYW^*(T^*)^j,
\qquad
 {\cal Z}_Y-T{\cal Z}_YT^*=WYW^*.                  \tag{7}
\]

With \(Q=I-VV^*\),

\[
 \boxed{
 \langle Y,{\cal M}_T(C)\rangle_{\rm HS}
 =2\operatorname{Re}\operatorname{tr}(C^*{\cal Z}_YV),
 \qquad
 {\cal M}_T^*(Y)=2Q{\cal Z}_YV.}                   \tag{8}
\]

Thus the sole remaining grade-one compatibility statement is

\[
 \boxed{
 \operatorname{tr}(YD_T)=0
 \quad\text{whenever}\quad
 Q{\cal Z}_YV=0.}                                  \tag{9}
\]

L206 proves (9) for every cokernel direction by identifying the
balanced dual Gramian with a self-adjoint colligation commutant and
applying L203 on its reducing spectral blocks.  The tracked
noncommutative calculations additionally show

\[
 \|D_T\|=O(\|H-H_0\|^2),\qquad
 \|C\|=O(\|H-H_0\|)                                \tag{10}
\]

along every tested line into the equality manifold.  This is exactly
the bounded divisibility required at the rank-jumping Crabb apex.
Equation (10) remains numerical evidence, not an analytic theorem.

## 2. Elimination of the lower endpoint

Decompose

\[
 C_2=VA+C,\qquad V^*C=0.                            \tag{11}
\]

Since \(TV=0\), (3) gives

\[
 V^*{\cal G}_T(H)V=V^*HV.                          \tag{12}
\]

Substitute (11) into (1) and apply \(V^*(\cdot)V\).  The lower
condition in (2) is

\[
 K+A+A^*=0.                                        \tag{13}
\]

Only the Hermitian part of \(A\) occurs in (1); its skew-Hermitian
part cancels from \(VC_2^*+C_2V^*\).  Hence the canonical choice is

\[
 A=-K/2,
\]

which gives precisely \(C_\parallel,F_0,X_0\) in (4).  Adding the
free perpendicular column \(C\) changes the metric by

\[
 X_C={\cal G}_T(VC^*+CV^*).                        \tag{14}
\]

Equation (12) makes \(V^*X_CV=0\), and its upper compression is (5).
This proves (6), including the converse.

## 3. Exact cokernel

The Stein maps

\[
 {\cal L}_T(X)=X-T^*XT,\qquad
 {\cal L}_T^\vee(Z)=Z-TZT^*
\]

are adjoints for the Hilbert--Schmidt pairing.  Therefore

\[
\begin{aligned}
\operatorname{tr}\{Y{\cal M}_T(C)\}
&=\operatorname{tr}\{
 WYW^*{\cal G}_T(VC^*+CV^*)\}\\
&=\operatorname{tr}\{
 {\cal Z}_Y(VC^*+CV^*)\}\\
&=2\operatorname{Re}\operatorname{tr}(C^*{\cal Z}_YV).
\end{aligned}                                      \tag{15}
\]

Restriction to \(V^*C=0\) inserts \(Q\) in the adjoint and proves
(8).  The finite-dimensional Fredholm alternative now gives (9).

At the Crabb apex, (5) is identically zero.  Indeed the terminal
orbit reaches \(V\) only at its last step, where \(V^*C=0\) kills the
contribution.  Thus the cokernel jumps from all of \(M_m^{\rm sa}\)
at the apex to a smaller space off it.  Applying an ordinary
constant-rank pseudoinverse across the apex would therefore be
invalid.

The identity matrix always lies in the cokernel of (5).  Equivalently,
every homogeneous endpoint correction has trace zero.  Pairing (14)
with L203's dual equality matrix proves this without an orbit
calculation.  L203's trace face is exactly

\[
 \operatorname{tr}D_T=0,                           \tag{16}
\]

the scalar instance of (9).

## 4. Numerical rank-jump audit

The deterministic equality lines use

\[
 \widehat H(s)^{-1}=2I+sN
\]

with noncommuting block-Toeplitz \(N\).  For

\[
 s\in\{0.03,0.1,0.3,1,3,8\},
\]

lengths \(2,\ldots,5\), and multiplicities \(2,3\), the checker finds:

1. (8) on complete real bases;
2. \(D_T\in\operatorname{ran}{\cal M}_T\);
3. the reconstructed full equation (1);
4. both endpoints (2);
5. \(\|D_T\|/s^2\) bounded; and
6. the minimum-frame solution \(\|C\|/s\) bounded.

For the deterministic multiplicity-two lines,
\(\operatorname{rank}{\cal M}_T=3=m^2-1\).  The multiplicity-three
lines retain one additional reducing cokernel direction and have
rank seven.  In both cases all cokernel conditions pass.  This
distinction is useful: trace compatibility alone is not the full
statement on reducible equality strata.

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_elliptic_cokernel.py \
  --output \
  experiments/repeated_crabb_elliptic_cokernel_s70224.jsonl
```

The floating range solve does not prove (9) or (10) analytically.
It identifies the exact object that the matrix-inner/model-complement
argument must control.  The tracked dataset SHA-256 is
`e722a548f36326f7e6d7277b784f269c640b8828d5971c72324a64c88b3025b2`.

## 5. Next gate

L206 closes the Fredholm compatibility in (9), and L207 subsequently
constructs the bounded real-analytic solution explicitly:

\[
 P^{-1/2}C=-\frac72(I-VV^*)SWB_1.
\]

Thus the rank jump does not require a pseudoinverse or a
Schur/Levinson recursion.

Either route must retain the left/right orientation: the target in
(6) is the left Gram \(-ZZ^*=-16B_1B_1^*\).  After grade one, the
same construction must be iterated on its kernel using
\(B_2,\ldots,B_L\); L201's invertible \(B_L\) supplies the terminal
coercivity.
