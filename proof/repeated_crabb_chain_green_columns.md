# Every exiting-chain Green column is a scalar continuant feature

## 1. Result (L245, 2026-07-24)

Retain L237's clean delay pencil

\[
L_r=N_r+c(I+F_0)N_r^*
\]

on levels \(0,\ldots,r-1\), and write

\[
H_r(z,c)=zI-L_r.
\]

Let \(e_j\) include the \(j\)-th copy-space level and put

\[
k_r=c(I+F_0)e_{r-1}.
\]

With L237's continuants

\[
\Delta_0=1,\quad \Delta_1=z,\quad
\Delta_2=z^2-2c,\quad
\Delta_j=z\Delta_{j-1}-c\Delta_{j-2},
\]

the two chain vectors needed by the **full** block inverse are

\[
\boxed{
e_{r-1}^*H_r^{-1}e_j
=\frac{\Delta_j}{\Delta_r}
\qquad(0\le j<r),}                                \tag{1}
\]

and

\[
\boxed{
e_j^*H_r^{-1}k_r
=
\begin{cases}
2c^r/\Delta_r,&j=0,\\
c^{r-j}\Delta_j/\Delta_r,&1\le j<r.
\end{cases}}                                      \tag{2}
\]

The formulas are scalar multiples of \(I_m\) at arbitrary defect
multiplicity.  Thus every amplitude entering or exiting the removed
chain is explicit; no unknown state column remains after L243.

On the Joukowski contour

\[
z=\zeta+c/\zeta,\qquad
x=\zeta^{-1},\qquad t=c/\zeta^2,
\]

L239's \(\Delta_j=\zeta^j(1+t^j)\) gives

\[
\boxed{
e_{r-1}^*H_r^{-1}e_j
=
\begin{cases}
x^r/(1+t^r),&j=0,\\
x^{r-j}(1+t^j)/(1+t^r),&1\le j<r,
\end{cases}}                                      \tag{3}
\]

and

\[
\boxed{
e_j^*H_r^{-1}k_r
=
\begin{cases}
2\rho^r/(1+t^r),&j=0,\\
\rho^{r-j}(1+t^j)/(1+t^r),&1\le j<r.
\end{cases}}                                      \tag{4}
\]

In particular their pairing is exactly the self-energy

\[
\boxed{
e_{r-1}^*H_r^{-1}k_r
=\rho+\zeta\,\frac{t^r(1-t)}{1+t^r}
=\gamma_r.}                                      \tag{5}
\]

Equations (3)--(4) have the same finite monomial-feature numerators
as L221/L243's model-kernel splitting.  All delay dependence is in
explicit scalar powers and the single denominator \(1+t^r\).

L245 does not yet perform the theta/ODE contour and Stein
calculation.  Its role is to close the “exiting-chain columns” item
left open in L242: the full finite-state resolvent is now obtained
by inserting (1)--(2) and L243's retained inverse into the standard
block formula below.

## 2. Full block inverse

With L237's tail notation, the complete resolvent denominator is

\[
\begin{bmatrix}
H_r&-k_rW^*\\
-We_{r-1}^*&D_0
\end{bmatrix},
\qquad
D_0=zI-\Xi_{r,-}.
\]

Its tail Schur complement is

\[
D_r=D_0-\gamma_rWW^*.
\]

Block inversion gives exactly

\[
\boxed{
\begin{bmatrix}
H_r^{-1}
+(H_r^{-1}k_r)W^*D_r^{-1}W
 (e_{r-1}^*H_r^{-1})
&
(H_r^{-1}k_r)W^*D_r^{-1}\\[1mm]
D_r^{-1}W(e_{r-1}^*H_r^{-1})
&
D_r^{-1}
\end{bmatrix}.}                                   \tag{6}
\]

L243 supplies \(D_r^{-1}\) with only the zero and first
nonconstant inverse-model-kernel sectors relevant through the
target order.  Equations (3)--(4) supply every remaining exterior
chain factor.  Consequently L228's unresolved computation is now a
finite Laurent-residue/Stein assembly with explicit inputs.

## 3. Continuant proof

The leading principal minors of \(H_r\) are precisely
\(\Delta_0,\ldots,\Delta_r\).  Cramer's rule along the last row gives

\[
e_{r-1}^*H_r^{-1}e_j=\Delta_j/\Delta_r,
\]

because every subdiagonal product is one.  Along the last column,
the reverse-edge product from level \(j\) to level \(r-1\) is
\(c^{r-1-j}\), with an additional factor two only when it crosses
the first edge.  Multiplication by
\(k_r=c(I+F_0)e_{r-1}\) gives (2), including the exceptional
\(r=1\) case.  Substitution of L239's continuant formula proves
(3)--(5).  Equation (6) is the ordinary Schur block inverse.

## 4. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_chain_green_columns.py \
  --output \
  experiments/repeated_crabb_chain_green_columns_s70224.jsonl
```

Delays one through eight are audited at nine complex Joukowski
points each.  The checker independently compares the matrix inverse
with (1)--(4) and verifies (5).  The tracked SHA-256 is
`17bb9477faec46b3ec81e92d32782bb377b308a9b8c96cdd6229d840021604a2`.
