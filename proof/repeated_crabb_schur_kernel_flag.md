# The ordered Schur kernel is an exact finite model-space flag

## 1. Result (L220, 2026-07-24)

Use L218's recursion

\[
F_j(z)=\Gamma_j+
D_{\Gamma_j^*}\,zF_{j+1}(z)
\bigl(I+\Gamma_j^*zF_{j+1}(z)\bigr)^{-1}D_{\Gamma_j},
                                                               \tag{1}
\]

with \(F_L=U\) unitary.  Put

\[
\begin{aligned}
K_j(z,w)
&=\frac{I-F_j(z)^*F_j(w)}{1-\overline zw},\\
A_j(w)
&=\bigl(I+\Gamma_j^*wF_{j+1}(w)\bigr)^{-1}
  D_{\Gamma_j}.
\end{aligned}                                                \tag{2}
\]

Then every Schur step has the exact ordered kernel splitting

\[
\boxed{
\begin{aligned}
K_j(z,w)
={}&A_j(z)^*A_j(w)\\
&+\overline zw\,
  A_j(z)^*K_{j+1}(z,w)A_j(w).
\end{aligned}}                                                \tag{3}
\]

Define

\[
\Phi_j(w)=w^jA_j(w)A_{j-1}(w)\cdots A_0(w),
\qquad 0\le j<L.                                      \tag{4}
\]

Iteration of (3) gives

\[
\boxed{
K_0(z,w)=\sum_{j=0}^{L-1}\Phi_j(z)^*\Phi_j(w).
}                                                               \tag{5}
\]

Near the repeated Crabb transfer, the feature columns in (5) span
\(\mathbb C^{Lm}\).  Consequently (5) is an exact isometric
identification of the \(Lm\)-dimensional model space of \(F_0\) with
\(L\) orthogonal copy-coordinate layers.

At the monomial apex all \(\Gamma_j=0\), so

\[
A_j=I,\qquad \Phi_j(w)=w^jI_m.                       \tag{6}
\]

Thus the ordinary Hardy coefficient grades are the associated
graded of an exact noncommutative Schur-orthogonal flag.  Away from
the apex, the ordered factors in (4) give the unique triangular
whitening of those grades.

L220 is a transfer/model-space theorem.  It does not yet say that
the grouped physical Riemann/metric endpoint is the negative norm of
these feature coordinates.  It reduces that remaining assertion to
an identification problem rather than another grade-by-grade
orthogonalization.

## 2. One-step derivation

L218's one-step identity, with

\[
Z(w)=wF_{j+1}(w),
\]

is

\[
\begin{aligned}
I-F_j(z)^*F_j(w)
={}&D_{\Gamma_j}
 (I+Z(z)^*\Gamma_j)^{-1}\\
&\quad\cdot[I-Z(z)^*Z(w)]\\
&\quad\cdot(I+\Gamma_j^*Z(w))^{-1}
 D_{\Gamma_j}.                                      \tag{7}
\end{aligned}
\]

The middle factor splits without any commutation:

\[
\begin{aligned}
I-Z(z)^*Z(w)
={}&I-\overline zwF_{j+1}(z)^*F_{j+1}(w)\\
={}&(1-\overline zw)I\\
&+\overline zw
  [I-F_{j+1}(z)^*F_{j+1}(w)].
\end{aligned}                                                \tag{8}
\]

Divide (7) by \(1-\overline zw\), substitute (2), and retain the
factor order.  This gives (3).

## 3. Finite iteration and isometry

Because \(F_L=U\) is constant unitary,

\[
K_L(z,w)=0.                                           \tag{9}
\]

Substitute (3) successively.  At the second layer the outer factors
combine as

\[
A_0(z)^*A_1(z)^*
=\bigl(A_1(z)A_0(z)\bigr)^*,
\]

and the same ordered multiplication persists at every layer.
This proves (4)--(5).

Let

\[
\Phi(w)=
\begin{bmatrix}
\Phi_0(w)\\ \vdots\\ \Phi_{L-1}(w)
\end{bmatrix}.
\]

Equation (5) reads \(K_0(z,w)=\Phi(z)^*\Phi(w)\).
At the apex, evaluating at \(L\) distinct points gives a block
Vandermonde matrix tensored with \(I_m\), hence full rank \(Lm\).
The same finite evaluation remains full rank in a neighbourhood.

L218 gives

\[
\operatorname{wind}\det F_0=Lm,
\]

so the finite model space with kernel \(K_0\) also has dimension
\(Lm\).  The spanning feature realization is therefore minimal:
every model-space vector has unique feature coordinates, and its
model norm is their Euclidean copy-space norm.  This is the precise
sense in which the \(L\) layers are orthogonal.

## 4. Consequence for the elliptic campaign

The candidate one-image frame is written in raw state-orbit
coordinates.  L217 shows that those coordinates mix future transfer
coefficients before their nominal grades.  Formula (4) supplies the
canonical replacement:

1. raw Hardy grades are first passed through the ordered factors
   \(A_j\);
2. the resulting \(\Phi_j\)-coordinates are exactly orthogonal; and
3. at the apex they reduce to the familiar monomials used in the
   scalar L190 compact face.

The next load-bearing identity should therefore be stated in these
coordinates:

> the first two-reflection physical upper endpoint is \(-16\) times
> the Gram of the corresponding \(\Phi\)-feature coefficients.

If that identity is established, negativity and all mixed-grade
orthogonality follow from (5), while L216 transports a delayed first
active layer to the grade-one physical calculation.

One tempting stronger statement is false: simply insert every
gradewise L212 face column into one raw frame series.  Two generic
noncommuting probes have the correct grade-one endpoint at order two
to \(4.5\times10^{-14}\), but acquire order-three residuals of norms
\(7.37\) and \(1.64\).  Their lower endpoints remain zero through
order six.  Thus this is a failure of the proposed upper identity,
not of lower tightening.  Later layers must be formed by the
Schur-orthogonal quotient (3), not by superposing raw face
representatives.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_schur_kernel_flag.py \
  --output \
  experiments/repeated_crabb_schur_kernel_flag_s70224.jsonl
```

The 16 records cover lengths two through five, multiplicities two
and three, repeated apices, and noncommuting
inverse-block-Toeplitz equality anchors.  They verify:

1. the direct kernel equals the ordered feature sum (5);
2. finite feature evaluations span all \(Lm\) coordinates; and
3. the apex features are exactly the monomials (6).

Equations (7)--(9), not the floating audit, prove L220.  The tracked
dataset SHA-256 is
`40b9611043413ea447fc729cc97e84242a658879bd346563307d0e5922d7dfb4`.

Regenerate the separate raw-superposition falsification with

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_raw_face_superposition.py \
  --output \
  experiments/repeated_crabb_raw_face_superposition_s70224.jsonl
```

Its tracked SHA-256 is
`8f3644128a5803fdeb9f23a82dc6f619eb92984aea8c550dc6abc5afab9d4874`.
