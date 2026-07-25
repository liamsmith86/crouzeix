# The delayed resolvent is a two-defect scattering problem

## 1. Result (L238, 2026-07-24)

Use L237 to remove a complete delay of length \(r\).  Relabel the
retained pure partial isometry and its balanced defects as

\[
T,\qquad I-T^*T=VV^*=:E,\qquad I-TT^*=WW^*=:F,
\qquad EF=0.
\]

Its transfer is

\[
\mathcal B(\rho)
=W^*(I-\rho T^*)^{-1}V
=\sum_{j\geq1}\rho^jB_j,\qquad
B_j=W^*(T^*)^jV.                                  \tag{1}
\]

L237 says that the retained resolvent denominator is

\[
D_r(z,c)
=zI-T-cT^*(I+E)-\gamma_r(z,c)F.                  \tag{2}
\]

Put

\[
z=\zeta+\rho,\qquad \rho=\frac c\zeta,\qquad
\delta_r=\gamma_r(z,c)-\rho.                      \tag{3}
\]

Then the denominator has the exact two-defect factorization

\[
\boxed{
D_r
=(\zeta I-T)(I-\rho T^*)-cT^*E-\delta_rF.}        \tag{4}
\]

Let

\[
G=(\zeta I-T)(I-\rho T^*),\qquad
R_0=G^{-1}
=(I-\rho T^*)^{-1}(\zeta I-T)^{-1},               \tag{5}
\]

\[
U=\begin{bmatrix}cT^*V&\delta_rW\end{bmatrix},
\qquad
C^*=\begin{bmatrix}V^*\\W^*\end{bmatrix}.
\]

Thus \(D_r=G-UC^*\).  Every interaction between the two endpoints is
now contained in the \(2m\)-dimensional Woodbury scattering matrix

\[
\mathcal M=C^*R_0U.                               \tag{6}
\]

Its blocks are

\[
\boxed{
\mathcal M=
\begin{bmatrix}
c\zeta^{-2}I&
\delta_r\,V^*(\zeta I-T)^{-1}W\\[2mm]
(z/\zeta)\mathcal B(\rho)&
\delta_r\,W^*R_0W
\end{bmatrix}.}                                   \tag{7}
\]

In particular, the lower-left block is exactly the reflected transfer
with no reordering of its matrix coefficients.  The upper-right block
is the oppositely oriented transfer:

\[
V^*(\zeta I-T)^{-1}W
=\sum_{j\geq1}\zeta^{-j-1}B_j^*.                  \tag{8}
\]

Whenever the inverses exist, the full retained resolvent is

\[
\boxed{
D_r^{-1}
=R_0+R_0U(I-\mathcal M)^{-1}C^*R_0.}              \tag{9}
\]

Equations (4), (7), and (9) are exact, not associated-graded
approximations.  They reduce the remaining L228 calculation to scalar
Joukowski/theta coefficients acting on a fixed \(2m\times2m\) endpoint
matrix.  The delay length enters through the single scalar
\(\delta_r\), whose first reflected coefficient is one by L237.

L238 does **not** prove L228.  It does not yet evaluate the full L125
theta/ODE contour operation, combine it with L219's boundary metric,
or take the upper right-defect Schur complement.  Those operations
must still be shown to turn the first \(\delta_r\)-dependent transfer
cell into

\[
E_1F_r+F_rE_1.
\]

## 2. Bulk factorization

Since \(TT^*=I-F\) and \(\zeta\rho=c\),

\[
\begin{aligned}
(\zeta I-T)(I-\rho T^*)
&=\zeta I-T-cT^*+\rho TT^*\\
&=(\zeta+\rho)I-T-cT^*-\rho F\\
&=zI-T-cT^*-\rho F.
\end{aligned}
\]

Subtracting \(cT^*E+(\gamma_r-\rho)F\) gives (2), proving
(4).  The factorization separates a completely factored bulk from
two rank-\(m\) endpoint corrections.

## 3. The scattering blocks

First, \(TV=0\), \(T^*W=0\), and \(V^*W=0\).  Therefore

\[
V^*(I-\rho T^*)^{-1}=V^*
\]

and

\[
(\zeta I-T)^{-1}T^*V
=\zeta^{-1}T^*V+\zeta^{-2}V.
\]

It follows that

\[
cV^*R_0T^*V=c\zeta^{-2}I,
\]

which is the upper-left block of (7).  The same first identity gives
the upper-right block directly:

\[
\delta_rV^*R_0W
=\delta_rV^*(\zeta I-T)^{-1}W.
\]

For the lower-left block, use the two-term resolvent identity above:

\[
\begin{aligned}
cW^*R_0T^*V
={}&c\zeta^{-1}
 \sum_{j\geq0}\rho^jB_{j+1}
 +c\zeta^{-2}\mathcal B(\rho)\\
={}&\left(\frac c{\zeta\rho}
          +\frac c{\zeta^2}\right)\mathcal B(\rho)\\
={}&\frac z\zeta\,\mathcal B(\rho).
\end{aligned}
\]

Here \(B_0=0\), so
\(\sum_{j\geq0}\rho^jB_{j+1}=\mathcal B(\rho)/\rho\).
The lower-right block is its defining compression
\(\delta_rW^*R_0W\).  This proves (7), including every multiplication
order.

Finally, (9) is the finite-rank Woodbury identity applied to
\(D_r=G-UC^*\).

## 4. What remains for the all-grade theorem

The exact reduction gives a bounded-size route that does not grow with
the delay:

1. express the full scalar direct map from L125 on the Joukowski
   contour \(z=\zeta+c/\zeta\);
2. substitute (9), retaining the complete scalar theta/ODE
   coefficients;
3. subtract the half-line term \(\delta_r=0\);
4. extract the first coefficient linear in
   \(\delta_r=c^r z^{-(2r-1)}+\cdots\);
5. combine it with L236's two Hardy-frame metric and take the upper
   right-defect Schur complement.

The transfer blocks in (7) show that this last coefficient can only
use the correctly oriented active Hankel cell.  What is still open is
its scalar multiplier and the cancellation with the Schur square.
No further finite-grade free-word expansion is warranted.

## 5. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_two_defect_scattering.py \
  --output \
  experiments/repeated_crabb_two_defect_scattering_s70224.jsonl
```

The six deterministic records use delays one through six, defect
multiplicities one through three, and noncommuting tails.  They audit
(4), all three simplified blocks in (7), (9), and equality with the
retained block of the original full resolvent.  The largest observed
error is below \(1.7\cdot10^{-15}\).  The tracked SHA-256 is
`c12108da1611f02e54fdce2862e7d459ee11b8010367209bca164d928ac30d30`.
