# The delayed scattering Schur complement has the boundary-metric weight

## 1. Result (L239, 2026-07-24)

Retain L238's deflated pure partial isometry

\[
T,\qquad I-T^*T=VV^*=:E,\qquad I-TT^*=WW^*=:F,
\qquad EF=0,
\]

and its transfer

\[
\mathcal B(u)=W^*(I-uT^*)^{-1}V
=\sum_{j\geq1}u^jB_j.
\]

On the Joukowski contour put

\[
z=\zeta+\frac c\zeta,\qquad
x=\frac1\zeta,\qquad \rho=\frac c\zeta,\qquad
t=\rho x=\frac c{\zeta^2}.
\]

L238's terminal reflection is

\[
\delta_r(z,c)=\gamma_r(z,c)-\rho.
\]

Then its normalized endpoint weight is not merely asymptotic:

\[
\boxed{
\lambda_r(t):=\frac{\delta_r}{\zeta(1-t)}
=\frac{t^r}{1+t^r}.}                              \tag{1}
\]

Thus the complete length-\(r\) terminal chain produces exactly the
same rational weight family \(u^r/(1+u^r)\) that occurs in L219's
boundary metric.  In particular,

\[
\lambda_r(t)=t^r+O(t^{2r})
\]

with leading coefficient one for every \(r\).

There is also an exact transfer-kernel collapse of L238's remaining
endpoint block.  Define the coefficientwise-adjoint transfer

\[
\mathcal B^\sharp(x)
=V^*(I-xT)^{-1}W
=\sum_{j\geq1}x^jB_j^*.                           \tag{2}
\]

For L238's factored bulk inverse

\[
R_0=(I-\rho T^*)^{-1}(\zeta I-T)^{-1},
\]

one has

\[
\boxed{
W^*R_0W
=\frac{x}{1-t}
 \{I-t\mathcal B(\rho)\mathcal B^\sharp(x)\}.}     \tag{3}
\]

Consequently L238's \(2m\)-dimensional scattering matrix is

\[
\mathcal M=
\begin{bmatrix}
tI&
\delta_r x\mathcal B^\sharp(x)\\[1mm]
(1+t)\mathcal B(\rho)&
\displaystyle
\delta_r\frac{x}{1-t}
\{I-t\mathcal B(\rho)\mathcal B^\sharp(x)\}
\end{bmatrix}.                                    \tag{4}
\]

Schur-eliminate its upper-left block from \(I-\mathcal M\).  All
apparently separate bulk and transfer terms combine to

\[
\boxed{
\operatorname {Schur}_{11}(I-\mathcal M)
=I-\lambda_r(t)
 \{I+\mathcal B(\rho)\mathcal B^\sharp(x)\}.}      \tag{5}
\]

Equations (1), (3), and (5) are exact and retain the order
\(\mathcal B(\rho)\mathcal B^\sharp(x)\).

L239 is a stronger scalar/transfer normal form, but it is **not yet
L228**.  The holomorphic functional-calculus contour still has to be
combined with the state lift in L238's Woodbury formula, L236's two
Hardy frames, L219's metric, and the right-defect Stein Schur square.
The important gain is that the delay chain and its next corrections
have disappeared completely: every delay now enters through the one
known scalar function (1).

## 2. The terminal continuants become \(1+t^r\)

L237 defines

\[
\Delta_1=z,\qquad \Delta_2=z^2-2c,\qquad
\Delta_j=z\Delta_{j-1}-c\Delta_{j-2}.
\]

On the Joukowski contour, \(z=\zeta(1+t)\) and
\(c=t\zeta^2\).  Directly,

\[
\Delta_1=\zeta(1+t),\qquad
\Delta_2=\zeta^2(1+t^2).
\]

The recurrence then gives, by induction,

\[
\boxed{\Delta_j=\zeta^j(1+t^j)\qquad(j\geq1).}     \tag{6}
\]

For \(r\geq2\), L237's formula yields

\[
\gamma_r
=c\frac{\Delta_{r-1}}{\Delta_r}
=t\zeta\frac{1+t^{r-1}}{1+t^r}.                   \tag{7}
\]

The same formula holds at \(r=1\), because its numerator becomes
\(1+t^0=2\) and L237 has \(\gamma_1=2c/z\).  Since
\(\rho=t\zeta\),

\[
\delta_r
=\zeta\frac{t^r(1-t)}{1+t^r}.
\]

Division by \(\zeta(1-t)\) proves (1).  This exact identity also
explains the cancellation of the apparent next coefficient
\((2r-2)t^{r+1}\) in the fixed-\(z\) expansion: changing to the
Joukowski coordinate absorbs it.

## 3. The lower-right bulk block is a model kernel

Apply the standard characteristic-kernel telescoping identity to
\(T^*\), whose defect frames are \(W,V\).  Its characteristic
function in these frames is

\[
\Theta_{T^*}(u)=u\mathcal B^\sharp(u).
\]

With independent formal variables \(\rho,x\), the identity reads

\[
\boxed{
W^*(I-\rho T^*)^{-1}(I-xT)^{-1}W
=\frac{
 I-\rho x\mathcal B(\rho)\mathcal B^\sharp(x)
}{1-\rho x}.}                                     \tag{8}
\]

It can also be verified directly by expanding both resolvents and
telescoping with \(I-TT^*=WW^*\), \(T^*W=0\), and \(V^*W=0\).

Since

\[
(\zeta I-T)^{-1}=x(I-xT)^{-1},
\]

equation (8), with \(\rho x=t\), proves (3).

## 4. Endpoint Schur elimination

Substitute (2)--(3) into L238's scattering blocks.  This gives (4).
The upper-left block of \(I-\mathcal M\) is the scalar matrix

\[
(1-t)I.
\]

Its lower Schur complement is

\[
\begin{aligned}
I
&-\delta_r\frac{x}{1-t}
 \{I-t\mathcal B(\rho)\mathcal B^\sharp(x)\}\\
&-(1+t)\mathcal B(\rho)
  \frac1{1-t}
  \delta_r x\mathcal B^\sharp(x).
\end{aligned}
\]

The two ordered transfer products have total coefficient

\[
-\frac{\delta_rx}{1-t}\{-t+(1+t)\}
=-\frac{\delta_rx}{1-t}.
\]

Therefore the last display is

\[
I-\frac{\delta_rx}{1-t}
 \{I+\mathcal B(\rho)\mathcal B^\sharp(x)\},
\]

which is (5) by (1).

## 5. Why this is the correct next interface

L237 supplied only the first coefficient of the terminal reflection.
Equation (1) now sums the entire reflection on the natural contour.
L238 left one unsimplified state compression \(W^*R_0W\); equation
(3) turns it into the exact model kernel of the same transfer used by
L236.  Finally, (5) shows that right-endpoint elimination introduces
no new delay-dependent matrix:

\[
\text{delay}
\quad\longleftrightarrow\quad
\lambda_r(t)=\frac{t^r}{1+t^r}.
\]

The remaining L228 calculation should therefore be done as one
Hardy/model-kernel contour identity.  At the first active
anti-diagonal it must:

1. use the coefficient \(1\) of \(t^r\) in (1);
2. retain the state-lift off-diagonal blocks from L238, rather than
   only the retained direct block;
3. pair them with L236's active cell \((1,r)\); and
4. include the right-defect Stein Schur square.

This avoids both false shortcuts recorded in L235.

## 6. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_scattering_schur_collapse.py \
  --output \
  experiments/repeated_crabb_scattering_schur_collapse_s70224.jsonl
```

The six deterministic records use delays one through six, defect
multiplicities one through three, and noncommuting retained tails.
They audit (1), (3), the full scattering matrix (4), and its Schur
collapse (5).  The largest observed error is below
\(3.1\cdot10^{-15}\).  The tracked SHA-256 is
`0dce5fe39244a9cfc5a64859c2568442b936735731a04440d47f918061e5dd70`.
