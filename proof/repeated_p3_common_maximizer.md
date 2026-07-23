# Common-maximizer cross quotient for repeated `C3` blocks (2026-07-22)

## 1. Setting

Let

\[
 A_0=I_m\otimes C_3,
 \qquad
 C_3=\begin{bmatrix}0&\sqrt2&0\\0&0&\sqrt2\\0&0&0\end{bmatrix},
\]

and write a perturbation `E` in `3 x 3` copy blocks `E_ab`.  The normalized top
support vector of `Re(q^{-1}C_3)`, `|q|=1`, is

\[
 r(q)=\frac12(q^{-1},\sqrt2,q)^T.                          \tag{1}
\]

L61 says that a repeated-block direction has zero first-order Jensen gap exactly
when the compressed Hermitian perturbations have one copy vector `y` that is a
top eigenvector for every `q`.  After a constant copy-unitary change, take
`y=e_1`.

This note classifies the cross-copy linear equations in that condition modulo
the exact unitary orbit.  The order relation saying that `y` is the **top**
eigenvector, rather than merely an eigenvector, remains an inequality on the
diagonal copy compression and is not used here.

## 2. One selected/cross copy pair

Fix `xi perpendicular to y` and put

\[
 X=E_{\xi y},\qquad Y=E_{y\xi}.
\]

The off-diagonal entry of the top-support compression is, up to a positive
constant,

\[
 r(q)^*\{q^{-1}X+qY^*\}r(q).                              \tag{2}
\]

Thus `y` is a common eigenvector only if, and for the cross entries if,

\[
 r(q)^*\{q^{-1}X+qY^*\}r(q)=0\quad(|q|=1).                \tag{3}
\]

Expanding (3) as a Laurent polynomial gives a real-linear map from the 36 real
entries of `(X,Y)` of rank 14.  Hence its kernel has real dimension 22.

## 3. Remove exact unitary mixing

An infinitesimal skew-Hermitian copy mixing with cross block `K` changes the
pair by

\[
 (X,Y)=\bigl(C_3K-KC_3,
 -(C_3K^*-K^*C_3)\bigr).                                  \tag{4}
\]

Every such pair satisfies (3), as it must because unitary similarity preserves
the disk.  The map from the 18 real entries of `K` in (4) has rank 16; its
kernel consists of the two-real-dimensional scalar commutant.  Therefore the
cross quotient has real dimension

\[
 22-16=6.                                                   \tag{5}
\]

The following three complex representatives span it.  For arbitrary
`alpha_0,alpha_1,alpha_2 in C`, put

\[
 X=\sum_{j=0}^2\bar\alpha_jX_j,
 \qquad Y=\sum_{j=0}^2\alpha_jY_j,                         \tag{6}
\]

where

\[
\begin{array}{c|c|c}
j&X_j&Y_j\\ \hline
0&\operatorname{diag}(-3/4,1/4,-3/4)&E_{02}\\[2mm]
1&-E_{01}+E_{12}&-E_{01}+E_{12}\\[2mm]
2&-(4/3)E_{02}&\operatorname{diag}(1,-1/3,1).
\end{array}                                                \tag{7}
\]

Each pair (6)--(7) satisfies both (3) and real Hilbert--Schmidt orthogonality
to (4).  The six real columns have rank six.  Since the combined support and
orbit-orthogonality constraints have rank 30 in 36 variables, these columns are
the complete quotient, not just examples.

## 4. Multiplicity and consequence

For general multiplicity `m`, equation (3) decouples for each vector in
`y^perp`.  Modulo cross-copy unitary mixing, the exceptional first-order face
therefore contains

\[
 3(m-1)\quad\text{complex cross parameters}.              \tag{8}
\]

This is substantially smaller than an arbitrary pair of cross blocks, but it
is not empty.  Consequently the repeated-block problem does **not** reduce at
first order to a single diagonal `C3` perturbation.  L75--L76 determine the
second-order conformal/similarity sign for one pure selected/cross-copy pair:
it is nonpositive, with a two-complex-dimensional equality plane.  L77 proves
strict cubic descent on that pure-pair plane, and L78 proves nonpositive second
order for simultaneous star coupling across all of `y^perp`.  Higher order on
the resulting multiplicity equality set, diagonal perturbations, and the
nonsmooth case where the diagonal support compression has no strict top gap
remain.

## 5. Regeneration

Run

```bash
.venv/bin/python -u experiments/repeated_p3_common_maximizer.py
```

The script builds every Fourier constraint from (1)--(3), builds the entire
unitary-orbit map (4), and verifies the ranks `14`, `16`, and `30` and all six
canonical representatives exactly in SymPy.
