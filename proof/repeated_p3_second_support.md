# Second support reduction on the repeated-`C3` cross quotient (2026-07-22)

## 1. Effective support problem

Use the multiplicity-two model `A_0=C_3\oplus C_3` and an L74 cross
representative

\[
 X=\sum_{j=0}^2\bar\alpha_jX_j,
 \qquad Y=\sum_{j=0}^2\alpha_jY_j.                         \tag{1}
\]

The perturbation has blocks `E_21=X`, `E_12=Y`, and zero diagonal blocks.
Equation (3) of L74 says that its compression to the repeated top-support
eigenspace vanishes identically.  Thus the first support variation is zero.

For `q=e^{i theta}`, let

\[
 H(q)=\operatorname{Re}(q^{-1}C_3),\qquad
 r(q)=\frac12(q^{-1},\sqrt2,q)^T.
\]

The eigenvalues of `H(q)` are `1,0,-1`.  Its reduced resolvent below the top
eigenvalue is the exact polynomial

\[
 R(q)=(I-H(q))^\dagger
 =I-\frac14H(q)-\frac34H(q)^2.                             \tag{2}
\]

Degenerate Hermitian perturbation theory says that the second support
coefficient is the largest eigenvalue of the effective matrix

\[
 Q(q)=P_{\rm top}H_E(q)R(q)H_E(q)P_{\rm top}.             \tag{3}
\]

Because `H_E` switches the two copies, (3) is diagonal in copy space.  Call its
two entries `q_1(q),q_2(q)`.

## 2. Closed form

Exact multiplication gives

\[
 \frac{q_1+q_2}{2}=M(q),\qquad q_1-q_2=\Delta(q),          \tag{4}
\]

where

\[
\begin{aligned}
M(q)={}&\frac5{128}|\alpha_0|^2+\frac14|\alpha_1|^2
       +\frac5{72}|\alpha_2|^2\\
 &+\frac1{16}\operatorname{Re}
       (\alpha_0\bar\alpha_2q^2),                         \tag{5}\\
\Delta(q)={}&\frac{\sqrt2}{12}\operatorname{Re}(qd),\\
d={}&3\alpha_0\bar\alpha_1+4\alpha_1\bar\alpha_2.       \tag{6}
\end{aligned}
\]

Consequently the one-sided second support function is

\[
 \boxed{\kappa(q)=M(q)+\frac{\sqrt2}{24}
             |\operatorname{Re}(qd)|.}                    \tag{7}
\]

When `d` is nonzero, the absolute cosine in (7) is the complete nonsmoothness
created by the unresolved repeated top eigenvalue.  When `d=0`, the second
effective branches do not split.

## 3. Conformal collapse

Despite the absolute value, (7) is `pi`-periodic in the boundary angle.  Its
first Fourier coefficient therefore vanishes.  Its zeroth coefficient is

\[
 \widehat\kappa(0)=
 \frac5{128}|\alpha_0|^2+\frac14|\alpha_1|^2
 +\frac5{72}|\alpha_2|^2+\frac{\sqrt2}{12\pi}|d|.         \tag{8}
\]

The second Schwarz coefficient of the disk-to-domain map is

\[
 K(w)=\widehat\kappa(0)w
      +2\sum_{k\ge1}\widehat\kappa(k)w^{k+1}.             \tag{9}
\]

Since `C_3^3=0` and `kappa_hat(1)=0`, every nonconstant boundary mode disappears
under functional calculus:

\[
 \boxed{K(C_3)=\widehat\kappa(0)C_3.}                     \tag{10}
\]

The first support coefficient and hence the first Schwarz correction are zero.
Therefore the Riemann-pulled repeated operator has the finite expansion

\[
 T_\epsilon=A_0+\epsilon E
 -\epsilon^2\widehat\kappa(0)A_0+o(\epsilon^2).            \tag{11}
\]

This reduces the L74 second-order similarity question to a `6 x 6` finite Stein
problem with the single scalar (8); no boundary discretization or infinite
Fourier tail remains.

Numerical orientation only: pure generator 1 has similarity-square coefficient
approximately `-8`, while pure generators 0 and 2 are second-order flat; generic
mixtures tested so far decrease.  These observations are not used in (2)--(11)
and the universal second-order sign is not yet proved.

## 4. Regeneration

Run

```bash
.venv/bin/python -u experiments/repeated_p3_second_support.py
```

The script constructs (2)--(4) from the L74 matrices and verifies (5)--(6)
exactly.  The average of the absolute cosine and the parity conclusion leading
to (8)--(10) are elementary scalar identities.
