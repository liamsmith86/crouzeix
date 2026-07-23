# Exact local slice at the `3x3` Crabb block (2026-07-22)

## 1. Affine-unitary orbit and its normal space

Let

\[
 A=C_3=\begin{bmatrix}0&\sqrt2&0\\0&0&\sqrt2\\0&0&0\end{bmatrix}.
\]

The exact value-preserving transformations are unitary conjugation, complex translation, and
nonzero complex scaling.  Rotation overlaps the unitary orbit of `A`, so the real tangent space
is

\[
 \mathcal O=\{[A,K]+\beta I+rA:K^*=-K,\ \beta\in\mathbb C,\ r\in\mathbb R\}. \tag{1}
\]

As in L66, its real dimension is `11`.  With the real Hilbert--Schmidt inner product, define

\[
 E(z,w,s,v)=
 \begin{bmatrix}
  -z/3&s&\bar z\\
  w&2z/3&-s\\
  v&w&-z/3
 \end{bmatrix},
 \qquad z,w,v\in\mathbb C,\quad s\in\mathbb R. \tag{2}
\]

Direct substitution in (1) gives

\[
 \operatorname{Re}\operatorname{tr}(E(z,w,s,v)^*O)=0
 \qquad(O\in\mathcal O). \tag{3}
\]

The space in (2) has real dimension seven, and `11+7=18`, the dimension of all complex `3x3`
matrices.  Hence

\[
 \boxed{\mathcal N:=\{E(z,w,s,v)\}=\mathcal O^\perp.} \tag{4}
\]

The four real coordinates `(z,w)` are exactly the L66 residual quotient representatives.
The remaining three coordinates `(s,v)` are transverse directions with negative L63 curvature.

## 2. Why this is an exact local reduction

Choose a linear complement to the infinitesimal stabilizer in the parameters `(K,beta,r)`.
The derivative at the identity of

\[
 (U,\beta,r;B)\longmapsto
 \Pi_{\mathcal O}\bigl(rU^*(B-\beta I)U-A\bigr) \tag{5}
\]

maps that parameter complement isomorphically onto `mathcal O` in (1).  The implicit-function theorem
therefore supplies, for every `B` sufficiently close to `A`, affine-unitary parameters for which
the transformed perturbation lies in `mathcal N`.  Since all transformations used in (5)
preserve the Riemann-pulled similarity square, the complete local problem may be studied on

\[
 A+E(z,w,s,v). \tag{6}
\]

This does not by itself give uniform estimates; it removes the orbit directions exactly and
identifies all remaining couplings which an estimate must control.

## 3. Exact second variation on the slice

L63 states, for an arbitrary perturbation,

\[
 e_3(E)=-2\bigl(\operatorname{Re}(E_{01}-E_{12})\bigr)^2
         -\frac{21}{4}|E_{20}|^2. \tag{7}
\]

For (2), `E_01-E_12=2s` and `E_20=v`, hence

\[
 \boxed{e_3(E(z,w,s,v))=-8s^2-\frac{21}{4}|v|^2.} \tag{8}
\]

Thus the Hessian is uniformly negative in exactly three slice directions and has kernel exactly
the four-real-dimensional `(z,w)` space treated raywise by L67--L68.  There are no omitted flat
directions after the exact symmetries are removed.

## 4. Remaining weighted problem

Equations (4) and (8) turn the `p=3` neighbourhood question into a finite weighted expansion.
The straight residual results are not alone uniform because the negative variables can shrink
with the flat ones and enter higher-order cross terms.  A candidate sharp chart begins with

\[
 A+\epsilon R_1(z)+\epsilon^2R_2(w)
   +\epsilon^3\bigl(sY_0+vY_3\bigr), \tag{9}
\]

provided the formal expansion confirms and eliminates every lower-weight center shift.  Here
`Y_0` is the real superdiagonal-difference direction and `Y_3` is the bottom-left direction
from (2).

L70 now derives the complete leading form and proves it nonpositive.  It has one unique weighted
zero, so the remaining target is the analytic center through that point rather than a generic
coupling estimate.  A punctured-neighbourhood theorem still requires controlling that center.

## 5. Audit

`experiments/p3_crabb_local_slice.py` independently reconstructs the orbit and the seven displayed
generators.  It checks orbit/slice/combined ranks `11/7/18`, Hilbert--Schmidt orthogonality, and
that the full L65 quadratic form restricted to the slice has exactly three negative eigenvalues.

Run

```bash
.venv/bin/python -u experiments/p3_crabb_local_slice.py
```
