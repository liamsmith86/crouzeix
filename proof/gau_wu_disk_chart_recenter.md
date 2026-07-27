# Recentring the disk chart at every Gau--Wu equality model

> **Status and scope.**  The polynomial-frame recentering, the
> inverse-Toeplitz identification, the strict disk-tangent block, and
> the boundary-shape quotient below are exact for every finite
> nondegenerate Gau--Wu model.  The sign of the final
> \(2n-2\)-dimensional quotient form remains open.  Thus this is a
> reduction of L342's live Hessian gate, not a local or global proof
> of the Crouzeix conjecture.

## 1. A polynomial support frame at an arbitrary model

Write
\[
 \phi=zf,\qquad S=S_\phi,\qquad
 A=XSX^{-1},
 \quad X=\sqrt2\oplus I\oplus {1\over\sqrt2},
 \tag{1}
\]
and let \(p\in\ker S\), \(q\in\ker S^*\).  As in L339, put
\[
 L=1\oplus\sqrt2 I\oplus1,\qquad
 y_\zeta=L(I-\overline\zeta S)^{-1}q.              \tag{2}
\]
Then \(y_\zeta\) spans the kernel of
\(I-\operatorname {Re}(\overline\zeta A)\).

Let
\[
 \delta(t)=\det(I-tS)=\sum_j\delta_jt^j
\]
and define
\[
 r_j=\sum_{k=0}^j\delta_kS^{j-k}q,\qquad0\le j<n.
 \tag{3}
\]
The adjugate identity gives
\[
 \operatorname {adj}(I-tS)q=\sum_{j=0}^{n-1}t^jr_j.
 \tag{4}
\]
Consequently the matrix
\[
 \boxed{C=[Lr_{n-1},Lr_{n-2},\ldots,Lr_0]}         \tag{5}
\]
satisfies, for
\({\bf f}_\zeta=(1,\zeta,\ldots,\zeta^{n-1})^T\),
\[
 C{\bf f}_\zeta
 =\zeta^{n-1}\delta(\overline\zeta)y_\zeta.        \tag{6}
\]
The scalar factor never vanishes on the circle.  The defect vector
\(q\) is cyclic for the finite completely nonunitary model \(S\), so
\(C\) is invertible.

This constructs the Lewis--Overton polynomial support frame directly
from the model resolvent; no boundary eigenvector phase selection is
needed.

## 2. Exact inverse-Toeplitz chart coordinates

Put
\[
 K=C^*C,\qquad A_c=C^{-1}AC,
\]
and let \(R\) be the unweighted superdiagonal shift.  Congruencing
the support slack by \(C\) gives the positive Laurent pencil
\[
 Q_\zeta
 =K-\overline\zeta\,{C^*AC\over2}
       -\zeta\,{C^*A^*C\over2},\qquad
 Q_\zeta{\bf f}_\zeta=0.                           \tag{7}
\]

The coefficient equations in (7) give a unique Hermitian matrix
\(\widehat H\), with final row and column zero, such that
\[
\boxed{
 K=\widehat H+R^*\widehat HR,\qquad
 KA_c=2\widehat HR.}                               \tag{8}
\]
Indeed, if \(G=C^*AC/2\), then
\(\widehat H=GR^*\), and the coefficient recurrence in (7) gives
both identities.  Moreover
\[
 Q_\zeta
 =(I-\zeta R^*)\widehat H(I-\overline\zeta R).
 \tag{9}
\]
Since \(I-\overline\zeta R\) is invertible and \(Q_\zeta\succeq0\),
\(\widehat H\succeq0\).  Its only kernel is the final coordinate.
Writing \(H\) for its positive leading \((n-1)\)-square block,
(8) is precisely L122's general disk chart.

This point is on its inverse-Toeplitz equality locus.  To see this
without a local-near-Crabb assumption, use L187's exact factorization
\[
 \Psi(H)=-(B_+-B_-)(B_++B_-)^{-1}{\cal Q},
 \qquad B=H^{-1}.                                  \tag{10}
\]
Here \({\cal Q}\) is invertible because the terminal vector is cyclic
for the companion matrix \(A_c\).  More explicitly, L183's canonical
Stein metric \(M\) obeys
\[
 hK\preceq M\preceq4hK.
\]
The characteristic Gau--Wu Blaschke product sends the terminal
coefficient vector to a vector of \(K\)-norm exactly twice as large.
Its \(M\)-contractivity sandwiches the squared norms between
\(4h\|e_L\|_K^2\) and the same number.  Every inequality is therefore
an equality, in particular
\((4hK-M)e_L=0\).  L183's Hardy defect then gives
\(\Psi(H)=0\).  Equation (10) therefore yields
\[
 \boxed{H^{-1}\ \hbox{is Hermitian Toeplitz}.}     \tag{11}
\]
Thus every finite nondegenerate Gau--Wu equality model, not only a
model near the Crabb collision, is an exact point of the same chart
used in L187--L189.

## 3. The strictly negative disk-tangent block

Vary \(H\) through all positive Hermitian matrices in the chart,
fixing its irrelevant positive scale.  These are genuine disk
matrices, so Okubo--Ando gives \(t_*\le4\) along every such curve.
At the inverse-Toeplitz anchor, the equality submanifold has tangent
\[
 \dot B_+=\dot B_-,
 \tag{12}
\]
of physical dimension \(2n-4\).

Differentiating (10) at \(B_+=B_-\) shows that the kernel of
\(D\Psi\) is exactly (12); cyclicity keeps \({\cal Q}\) invertible.
L183's exact Hardy defect is the squared norm of this residual.
Hence, after quotienting the equality tangent and scale, every
nonzero general-\(H\) disk direction has a strictly negative
second-order similarity coefficient.  By L342 it has the same
strict scalar coefficient up to the factor four.

The resulting subspace \({\cal D}_\phi\) of L342's physical normal
space has dimension
\[
\boxed{
 \dim_{\mathbb R}{\cal D}_\phi
 =\{(n-1)^2-1\}-(2n-4)=(n-2)^2.}                  \tag{13}
\]
The common Hessian is negative definite on
\({\cal D}_\phi\).

## 4. The canonical boundary-shape quotient

For an ambient first direction \(E\), let \(s_E(\zeta)\) be its first
support variation.  Put
\[
 \omega_\phi(\zeta)
 =|\delta(\overline\zeta)|^2D_f(\zeta)
 ={1\over2}\|C{\bf f}_\zeta\|^2.                  \tag{14}
\]
From (6),
\[
 \omega_\phi(\zeta)s_E(\zeta)
 ={1\over4}
 \left\langle C{\bf f}_\zeta,
  \{\overline\zeta E+\zeta E^*\}C{\bf f}_\zeta
 \right\rangle                                    \tag{15}
\]
is a real trigonometric polynomial of degree at most \(n\).
The complete polynomial support frame makes the ambient map onto the
\(2n+1\) real coefficients of that polynomial.

Translations and positive scaling occupy exactly the three support
functions
\[
 \omega_\phi(\zeta),\qquad
 \omega_\phi(\zeta)\cos\theta,\qquad
 \omega_\phi(\zeta)\sin\theta.                     \tag{16}
\]
Unitary, rotational, and model-zero equality tangents have zero
support variation.  Therefore the physical, affine-quotiented shape
map
\[
 \boxed{
 \Sigma_\phi(E)
 =[\omega_\phi s_E]\ \bmod\
 \omega_\phi\operatorname {span}_{\mathbb R}
 \{1,\cos\theta,\sin\theta\}}
 \tag{17}
\]
has rank \(2n-2\).

Every disk-chart tangent has constant disk support before taking the
physical slice, so \({\cal D}_\phi\subseteq\ker\Sigma_\phi\).
The physical normal dimension is
\[
 (n-1)^2+1=(n-2)^2+(2n-2).
 \tag{18}
\]
Equations (13), (17), and (18) prove
\[
\boxed{\ker\Sigma_\phi={\cal D}_\phi.}             \tag{19}
\]
Thus no quadratic-size ambiguity remains: the directions not already
controlled by the global disk theorem are exactly \(2n-2\) boundary
shape coordinates.

Let \({\cal H}_\phi\) be L342's common optimized Hessian.  Since its
restriction to \({\cal D}_\phi\) is negative definite, define the
finite Schur form
\[
 \boxed{
 \widehat{\cal H}_\phi(\sigma)
 =\max_{\Sigma_\phi(C)=\sigma}{\cal H}_\phi(C).}   \tag{20}
\]
The maximum is unique along each affine fibre, and
\[
 {\cal H}_\phi\preceq0
 \quad\Longleftrightarrow\quad
 \widehat{\cal H}_\phi\preceq0.                   \tag{21}
\]
This is the revised sole local gate.  It has real dimension
\(2n-2\), exactly L342's boundary-row dimension.  L339's support
port is therefore not merely a low-rank numerical hint: it is the
canonical quotient left after the known disk theorem has removed
the entire \((n-2)^2\)-dimensional bulk.

## 5. Deterministic audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_disk_chart_recenter.py
```

The checker independently constructs (5) from the adjugate
coefficients, verifies (6), (8), and (11), differentiates the full
Hermitian disk chart, and compares its projected tangent with the
kernel of (16).  It also records the two sign checks without using
them as proofs: the disk restriction is uniformly strict in the
sample, and the open shape Schur form remains negative.

The tracked audit covers one generic model in each dimension
\(4,\ldots,8\).  The largest polynomial-null, chart, and
inverse-Toeplitz errors are respectively
\(2.12\cdot10^{-15}\), \(1.10\cdot10^{-14}\), and
\(2.23\cdot10^{-14}\); the weighted-support Fourier tail beyond
degree \(n\) is at most \(1.48\cdot10^{-15}\).
The disk/shape kernel projectors agree within
\(2.27\cdot10^{-8}\).  The disk block's largest eigenvalue is at most
\(-4.07\), while the open quotient's softest tested eigenvalue is
\(-2.16\cdot10^{-9}\).  That final numerical sign is evidence only.

The dataset SHA-256 is
`01dabed0b2aad2b584210928396051bd74a585cb5832fec49b3a83036ee10103`.
