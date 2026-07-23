# Three-copy nonnormal flat stratum (2026-07-22)

## 1. Schur normal form

Take the L88 copy matrix with zero diagonal and three copies.  After phase
choices its general strictly upper form is

\[
 N=\begin{bmatrix}0&a&b\\0&0&c\\0&0&0\end{bmatrix}.                   \tag{1}
\]

L88 gives

\[
 128Q_N(q)=5(NN^*+N^*N)-3\{q^2N^2+q^{-2}(N^*)^2\}.                    \tag{2}
\]

Since

\[
 N^2=acE_{02},                                                         \tag{3}
\]

there are only two equality mechanisms.

## 2. Exact common-top classification

Suppose first that `ac!=0` and `x` is a common eigenvector of every `Q_N(q)`.
The `q^2` and `q^-2` coefficients in (2) force

\[
 E_{02}x=E_{20}x=0,
\]

so `x` is the middle coordinate `e_1`.  The constant coefficient acts by

\[
 (NN^*+N^*N)e_1
 =\bar c b\,e_0+(|a|^2+|c|^2)e_1+\bar b a\,e_2.                       \tag{4}
\]

Thus `e_1` is a common eigenvector exactly when `b=0`.

On that face, its eigenvalue in the scale of (2) is
`5(|a|^2+|c|^2)`.  The defect on `span{e_0,e_2}` has diagonal

\[
 5|c|^2,\qquad5|a|^2
\]

and determinant

\[
 16|ac|^2>0.                                                          \tag{5}
\]

Hence `e_1` is the unique common top branch for every `q`.

If `ac=0`, then `N^2=0` and (2) is independent of `q`, so it automatically
has a fixed top eigenspace.  We have proved the exact alternative

\[
 \boxed{\text{the L88 matrix-Jensen endpoint is flat only if }abc=0.} \tag{6}
\]

When `abc!=0`, there is no common eigenvector and the matrix-Jensen term is
strictly negative.

## 3. Existing theorems close both equality faces

### Square-zero face

If `ac=0`, a nonzero `3 x 3` square-zero matrix has rank at most one.  Its
range lies in its kernel, and unitary singular-value reduction gives one
`2 x 2` nilpotent copy block plus a zero copy.  In the L88 tensor form this is
exactly L77's pure generator-zero pair, together with an uncoupled `C_3`.
L77 supplies strict cubic descent on the pair.  The extra `C_3` block is
already contained in the pair's enlarged numerical-range domain and does not
weaken the complete-2 estimate.

### Length-two path face

If `b=0` and `ac!=0`, the selected middle copy is coupled independently to
the two outer copies.  In its selected-copy basis the two flat coefficient
vectors are independent, so this is precisely L80's rank-two star.  The
strict complement margin (5) and L80's improved metric make the full second
endpoint negative definite.

Therefore every nonzero zero-diagonal three-copy Schur direction is closed:

- `abc!=0`: strict at canonical second order;
- `b=0, ac!=0`: strict at optimized second order;
- `ac=0`: strict at cubic order.

## 4. Scope and later resolution

This theorem closes the smallest genuinely nonnormal copy-matrix stratum.
L90--L92 subsequently insert the Schur diagonal/common `w` mode and close the
two-copy terminal block.  L93 then replaces the proposed larger-matrix
path/star induction by a dimension-free metric-kernel flag and a
two-dimensional Clifford terminal classification.  Thus no larger fixed
Schur stratum remains; only the uniform gap-collapse problem does.

## 5. Regeneration

Run

```bash
.venv/bin/python -u experiments/repeated_p3_flat_three_copy.py
```

The checker derives (3)--(5) with independent formal adjoints and verifies
that both square-zero faces lose every angular Fourier mode.
