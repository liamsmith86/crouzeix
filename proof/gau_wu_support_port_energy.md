# The Gau--Wu support reserve is an unweighted port energy

> **Scope.**  This note rewrites L335's positive weighted support
> Gram as the ordinary \(L^2\)-energy of the canonical Ando support
> factor.  It supplies the correct port for L338's remaining
> contraction estimate, but does not prove that estimate.

## 1. The support factor and its null state

Let

\[
 \phi=zf,\qquad S=S(\phi),\qquad
 p=f\in\ker S,\quad q=1\in\ker S^*,
\]

and decompose

\[
 {\cal K}_\phi=\mathbb Cp\oplus{\cal M}\oplus\mathbb Cq.
\]

The model shift is a defect-one partial isometry:

\[
 S^*S=I-pp^*,\qquad SS^*=I-qq^*.
\]

Hence

\[
 U=S+qp^*
\]

is unitary.  Put

\[
\begin{aligned}
 D&=0\oplus {1\over\sqrt2}I_{\cal M}\oplus1,\\
 E&=1\oplus {1\over\sqrt2}I_{\cal M}\oplus0,\\
 L&=1\oplus\sqrt2 I_{\cal M}\oplus1.
\end{aligned}                                     \tag{1}
\]

For the Gau--Wu matrix

\[
 A=XSX^{-1},\qquad
 X=\sqrt2\oplus I_{\cal M}\oplus{1\over\sqrt2},
\]

one has

\[
 A=2EUD.                                           \tag{2}
\]

For \(\zeta\in\mathbb T\), define the support pencil

\[
 F_\zeta=UD-\zeta E.                               \tag{3}
\]

Direct multiplication gives the exact Ando factorization

\[
\boxed{
 F_\zeta^*F_\zeta
 =I-\operatorname{Re}(\overline\zeta A).}         \tag{4}
\]

Let

\[
 k_\zeta^\phi=(I-\overline\zeta S)^{-1}q,\qquad
 y_\zeta=Lk_\zeta^\phi.                            \tag{5}
\]

Then

\[
\boxed{
 F_\zeta y_\zeta=0,\qquad
 \|y_\zeta\|^2=2D_f(\zeta),}                       \tag{6}
\]

where

\[
 D_f(\zeta)={\zeta f'(\zeta)\over f(\zeta)}.
\]

Thus the normalized top-support vector of
\(\operatorname{Re}(\overline\zeta A)\) is

\[
 v_\zeta={y_\zeta\over\sqrt{2D_f(\zeta)}}.         \tag{7}
\]

## 2. Proof of the null-state norm

The standard boundary model-kernel identity gives

\[
 \|k_\zeta^\phi\|^2
 =D_\phi(\zeta)=1+D_f(\zeta).                     \tag{8}
\]

Reproduction at \(0\) and at \(p=f\) gives

\[
 |\langle q,k_\zeta^\phi\rangle|=1,\qquad
 |\langle p,k_\zeta^\phi\rangle|=|f(\zeta)|=1.
\]

The middle component therefore has squared norm \(D_f(\zeta)-1\).
Applying \(L\) proves

\[
 \|Lk_\zeta^\phi\|^2
 =1+2(D_f(\zeta)-1)+1=2D_f(\zeta).
\]

The equation \(F_\zeta y_\zeta=0\) follows either directly from
\((I-\overline\zeta S)k_\zeta^\phi=q\) in the three summands, or from
the defect-one unitary completion \(U=S+qp^*\).  Equations (4) and
(6) then prove (7).

## 3. L335 becomes an ordinary port norm

For a physical direction \(C\), write

\[
 H_C(\zeta)
 ={\,\overline\zeta C+\zeta C^*\over2}.
\]

Let \((F_\zeta^*)^\dagger\) be the Moore--Penrose inverse on the
orthogonal complement of \(y_\zeta\), and define the minimal support
response

\[
\boxed{
 \gamma_C(\zeta)
 =(F_\zeta^*)^\dagger H_C(\zeta)y_\zeta.}         \tag{9}
\]

Since \(F_\zeta^*F_\zeta\) is the support slack in (4),

\[
\begin{aligned}
 \|\gamma_C(\zeta)\|^2
 &=\langle H_Cy_\zeta,
   (F_\zeta^*F_\zeta)^\dagger H_Cy_\zeta\rangle\\
 &=2D_f(\zeta)\,
   \langle H_Cv_\zeta,
   (I-\operatorname{Re}(\overline\zeta A))^\dagger
   H_Cv_\zeta\rangle.                             \tag{10}
\end{aligned}
\]

The last scalar is L335's second-support coefficient \(t_C(\zeta)\);
the pseudoinverse automatically removes the top-vector component.
Consequently L335's reserve is exactly

\[
\boxed{
 {\cal P}_\phi(C)
 =2\int_{\mathbb T}D_f(\zeta)t_C(\zeta)\,dm(\zeta)
 =\int_{\mathbb T}\|\gamma_C(\zeta)\|^2\,dm(\zeta).} \tag{11}
\]

There is no longer a separate Poisson weight: it is the norm of the
unnormalized support null state.

## 4. Revised form of the live estimate

L338 reduced the full endpoint-flux sign to

\[
 Q(C,0)\geq2\|\kappa_C\|^2.
\]

The second-support part of \(Q(C,0)\) is
\(2{\cal P}_\phi(C)\).  Equation (11) makes that contribution

\[
 2\|\gamma_C\|_{L^2(\mathbb T)}^2.
\]

The remaining task is therefore a lossless-port estimate: express
L338's Riesz response \(\kappa_C\), together with the other
first-order endpoint loss, as the output of the unitary colligation
\(U=S+qp^*\), and prove that its energy is bounded by the physical
input energy plus \(2\|\gamma_C\|_{L^2}^2\).  This note identifies
the input port exactly; it does not assume the contraction.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_support_port_energy.py \
  --output experiments/gau_wu_support_port_energy_s70224.jsonl
```

The checker independently verifies the Ando factor, the boundary
null state and its norm, and the fully polarized equality (11).
The 12-record dataset has SHA-256
`cae577fc8c856279b734f6e6a7481d0eefc42f531af5ccea80e4863a2371a2d3`;
the largest port-Gram residual is \(1.60\cdot10^{-14}\).
