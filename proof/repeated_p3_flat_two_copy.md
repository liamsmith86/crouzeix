# Complete two-copy flat core (2026-07-22)

## 1. Second-order classification

Set the common `w` mode to zero first.  By a copy unitary, write an arbitrary
two-copy matrix in Schur form

\[
 Z=\begin{bmatrix}d_0&a\\0&d_1\end{bmatrix}.                           \tag{1}
\]

If `a=0`, then `Z` is normal and L88 gives the exact complete-2 theorem.
Assume `a!=0`.  A common eigenvector of L88's support family must reduce both
`Z^2` and `(Z*)^2`.  In two dimensions this requires `Z^2` to be normal.  But

\[
 Z^2=\begin{bmatrix}d_0^2&a(d_0+d_1)\\0&d_1^2\end{bmatrix},             \tag{2}
\]

and an upper-triangular matrix is normal only when its off-diagonal entry
vanishes.  Hence the nonnormal second-order equality condition is exactly

\[
 d_0+d_1=\operatorname{tr}Z=0.                                        \tag{3}
\]

If (3) fails, the matrix-Jensen term is strictly negative.  If (3) holds,
write

\[
 Z=\begin{bmatrix}d&a\\0&-d\end{bmatrix}.                              \tag{4}
\]

Then

\[
 Z^2=d^2I,qquad ZZ^*+Z^*Z=(2|d|^2+|a|^2)I,                            \tag{5}
\]

so the entire second effective support is scalar.

## 2. Third effective support

For the physical flat perturbation `E_Z`, exact reduced-resolvent multiplication
gives a traceless Hermitian third effective matrix.  Its upper cross entry is

\[
 \frac{a}{128q}\left[
 3d^2q^4+(2|a|^2+4|d|^2)q^2-7\bar d^2
 \right].                                                              \tag{6}
\]

If `a!=0`, (6) cannot vanish identically: its middle coefficient contains the
strictly positive scalar `2|a|^2+4|d|^2`.  Hence the top third support branch
has positive mean `m_3`.  The third matrix changes sign under `q->-q`; because
it is traceless Hermitian, its top eigenvalue is unchanged.  That eigenvalue is
therefore `pi`-periodic and its first Fourier coefficient vanishes.

The complete tight third metric, including the opposite diagonal loops in
(4), propagates exactly as in L77 and yields

\[
 \boxed{\mathcal E_3=-16m_3I_2\prec0.}                                 \tag{7}
\]

Thus every nonnormal trace-zero two-copy flat direction descends strictly at
cubic order.  L77 is the special case `d=0`.

The middle Fourier coefficient in (6) gives a quantitative version.  Since
the top eigenvalue of a traceless Hermitian `2 x 2` matrix dominates the
modulus of either off-diagonal entry, and the mean modulus of a Fourier
function dominates every Fourier coefficient,

\[
 m_3\ge
 \frac{|a|}{128}\left(2|a|^2+4|d|^2\right).             \tag{8}
\]

Consequently

\[
 \boxed{
 \lambda_{\max}(\mathcal E_3)
 \le-\frac{|a|}{4}\left(|a|^2+2|d|^2\right).}           \tag{9}
\]

This coercivity is uniform on normalized trace-zero blocks away from the
normal face `a=0`.  It is the terminal estimate needed by L93's future
weighted compactness argument.

## 3. Restore the common mode

If `w!=0`, L91 says a common top vector must reduce `Z`.  A nonnormal `2 x 2`
matrix is unitarily irreducible, so its matrix-Jensen gap is strict.  If `Z`
is reducible, it is normal and L88 applies.

Combining all cases closes the complete two-copy L87 flat core:

- normal `Z`: exact complete-2 bound from L88;
- nonnormal `Z,w!=0`: strict second-order Jensen gap from L91;
- nonnormal `Z,w=0,tr Z!=0`: strict second-order Jensen gap from (2);
- nonnormal `Z,w=0,tr Z=0`: strict cubic descent from (7).

The remaining weighted issue is uniformity as `w`, `tr Z`, and the cubic
nonnormal margin tend to zero together; there is no unclassified fixed
two-copy direction.  L93 subsequently proves that these one-/two-dimensional
blocks are the only irreducible terminal equality blocks in any copy
multiplicity.

## 4. Regeneration

Run

```bash
.venv/bin/python -u experiments/repeated_p3_flat_two_copy_third.py
.venv/bin/python -u experiments/repeated_p3_third_sign.py
```

The first checker derives (5)--(7) with arbitrary complex `d,a` and builds the
full third metric through the shared exact propagation module.  The second
independently regenerates the L77 `d=0` boundary case after the same refactor.
