# Independent audit of the August 2026 scalar proof

**Audit date:** 2026-09-05.  **Auditor:** independently assigned physical-factor
agent, redirected before continuing the older local-Hessian program.

**Primary source:** E. Lorist and F. L. Schwenninger, *A solution to
Crouzeix's conjecture*, arXiv:2608.03841v1, 2026-08-04,
<https://arxiv.org/html/2608.03841v1>.
The revised v2, dated 2026-08-17, was also inspected after the independent
reconstruction: <https://arxiv.org/html/2608.03841v2>.  Its reorganized
power-lemma argument agrees with the finite-dimensional derivation below.

**Conclusion of this mathematical audit:** the finite-dimensional scalar
argument checks.  Below is a reconstruction with the limit justification
made explicit and with the function-space Cauchy-transform assertion
removed entirely: a direct resolvent integral suffices.  This is an audit
and reconstruction of the cited authors' proof, not a claim of independent
discovery.  It does not establish the completely bounded conjecture or
resolve the old local-factor proof debts by their proposed methods.

## 1. The power lemma, checked without omitted limit steps

Let `H` be finite dimensional, `V:H→K` an isometry, and `Q` a contraction.
Suppose

\[
 E_j=2V^*Q^{*j}V-T^{*j}\quad(j\geq1),\qquad
 [E_j,T]=0,\qquad M:=\sup_j\|E_j\|<\infty.
\]

Write `κ=||T||`.  Only `κ>1` needs consideration.  Choose a unit right
singular vector `x`, so `T*Tx=κ²x`.  Use the inner product linear in its
first argument.  Set

\[
 S_j^*=2V^*Q^{*j}V,\quad
 m_j=\operatorname{Re}\langle E_jT^jx,x\rangle,\quad
 y_j=(S_{j+1}^*T-\kappa S_j^*)x.
\]

Commutation gives exactly

\[
\begin{aligned}
\kappa m_j-m_{j+1}
 &=\operatorname{Re}\langle
 T^j(\kappa E_j-E_{j+1}T)x,x\rangle\\
 &=(\kappa^2-\kappa)\|T^{*j}x\|^2
   -\operatorname{Re}\langle y_j,T^{*j}x\rangle.
\end{aligned}
\tag{A1}
\]

The potentially delicate noncommutative substitution is legitimate:
`T^{*(j+1)}Tx=T^{*j}(T*Tx)=κ²T^{*j}x`.
There is no step moving `T` through `T*`.

Put `b=κ(κ−1)>0` and

\[
 D=\|Q^*VTx-\kappa Vx\|^2\geq0.
\]

Completing the square in (A1), then using
`y_j=2V*Q^{*j}(Q*VTx−κVx)`, yields

\[
 \kappa m_j-m_{j+1}
 \geq-\frac{\|y_j\|^2}{4b}
 \geq-\frac{D}{b}.
\tag{A2}
\]

Iterate the first-order scalar recurrence, retaining its terminal term:

\[
 m_1\geq \kappa^{-N}m_{N+1}
       -\frac{D}{\kappa(\kappa-1)}
          \sum_{j=1}^N\kappa^{-j}.
\tag{A3}
\]

The terminal term does vanish.  The hypotheses give
`||T^j||=||T^{*j}||≤2+M`, and therefore
`|m_j|≤M(2+M)` for every `j`.  Thus (A3) gives

\[
 m_1\geq-\frac{D}{\kappa(\kappa-1)^2}.
\tag{A4}
\]

Independently, expand `D`, use contractivity and the isometry, and insert
`2V*Q*V=E_1+T*`:

\[
\begin{aligned}
 D&\leq2\kappa^2
     -2\kappa\operatorname{Re}\langle V^*Q^*VTx,x\rangle\\
  &=2\kappa^2-\kappa m_1-\kappa^3.
\end{aligned}
\tag{A5}
\]

Combining (A4) and (A5),

\[
 D\left(1-\frac1{(\kappa-1)^2}\right)
 \leq\kappa^2(2-\kappa).
\tag{A6}
\]

For `κ>2` the left side is nonnegative and the right side strictly
negative.  Hence `||T||≤2`.

### Hypothesis stress tests

* **Uniformity is essential.** With `Q=0`, `V=I`, and `T=3I`, every
  `E_j=−3^jI` commutes with `T`, but uniformity fails and the conclusion
  is false.
* **Commutation is essential.** With `Q=0`, `V=I`, and
  `T=[[0,3],[0,0]]`, the sequence is bounded and zero for `j≥2`, but
  `E_1=−T*` does not commute with `T`; again the conclusion is false.
* **One common contraction matters.** The proof uses the factorization
  of `y_j` through powers of the *same* `Q`.  Independent bounds on
  each `S_j` do not justify that factorization.
* **Equality is not a contradiction.** At `κ=2`, (A6) has both its
  scalar factor and its right side zero.
* **The lemma is scalar in the later application.** Matrix-amplified
  analytic functions need not commute with the correction integrals.

## 2. Resolvent construction: no abstract Cauchy-transform bound needed

Let `A` be any finite complex matrix.  Choose a bounded convex domain
`Ω` with smooth boundary and `W(A)⊂Ω`.  Write `σ` for its positively
oriented boundary point, `ν(σ)` for its outward unit complex normal, and
`ds` for arclength.  Then `dσ=iν ds`.  Define

\[
 B_\sigma=\sigma I-A,\qquad
 P(\sigma)=\frac1{2\pi}
       \{\nu B_\sigma^{-1}+\overline\nu B_\sigma^{-*}\}.
\tag{A7}
\]

The spectrum lies in `W(A)` (use an eigenvector), so all the inverses
exist on the boundary.

### Positivity, with the adjoints in the correct order

Convex support gives
`Re(conj(ν)(σ−z))>0` for every `z∈W(A)`.  Consequently

\[
 B_\sigma^*\operatorname{Re}(\nu B_\sigma^{-1})B_\sigma
 =\operatorname{Re}(\overline\nu B_\sigma)\succ0.
\tag{A8}
\]

Thus `P(σ)≻0`.  Cauchy's integral formula gives

\[
 \int_{\partial\Omega}P(\sigma)\,ds=2I.
\tag{A9}
\]

Let `f` be a polynomial normalized by `max_{closure Ω}|f|=1`, put
`T=f(A)`, and take

\[
 K=L^2(\partial\Omega,ds;H),\qquad
 (Vx)(\sigma)=2^{-1/2}P(\sigma)^{1/2}x,\qquad
 (Qg)(\sigma)=f(\sigma)g(\sigma).
\tag{A10}
\]

Equation (A9) proves `V*V=I`; normalization proves `||Q||≤1`.

Direct contour calculation, using only the ordinary polynomial
resolvent formula, gives

\[
\begin{aligned}
 2V^*Q^{*j}V
 &=\int\overline{f(\sigma)}^{\,j}P(\sigma)\,ds\\
 &=T^{*j}+\frac1{2\pi i}
       \int_{\partial\Omega}
       \overline{f(\sigma)}^{\,j}(\sigma I-A)^{-1}\,d\sigma.
\end{aligned}
\tag{A11}
\]

Define `E_j` to be the last integral.  Its integrand need not be
holomorphic as a function of the boundary variable.  This causes no
problem: it is an ordinary continuous contour integral, not an
application of Cauchy's theorem to `conj(f)^j`.

Each resolvent commutes with `f(A)`, hence so does its integral `E_j`.
Moreover,

\[
 \|E_j\|\leq
 \frac{\operatorname{length}(\partial\Omega)}{2\pi}
 \max_{\sigma\in\partial\Omega}\|(\sigma I-A)^{-1}\|=:M_{A,\Omega}.
\tag{A12}
\]

This number is finite, may depend on `A`, dimension and `Ω`, and is
independent of `j`.  No Crouzeix estimate is used to obtain it.  The
power lemma therefore gives

\[
 \|p(A)\|\leq2\max_{\overline\Omega}|p|
\tag{A13}
\]

for every polynomial `p`; the zero polynomial is handled separately.
Dependence of the intermediate bound (A12) on dimension does not
survive into (A13), and no uniform bound as the contour shrinks is
required: the lemma is applied separately on each fixed contour.

## 3. Smooth outer approximation, including degenerate numerical ranges

Here is an explicit smoothing argument, so no smoothness of `W(A)`
is assumed.  Let `C=W(A)`, a compact convex subset of the plane.  For
`δ>0`, take a nonnegative smooth mollifier `ρ_δ` of integral one
supported in the radius-`δ` disk, and set

\[
 F_\delta(z)=\int d(z-y,C)^2\rho_\delta(y)\,dy,
 \qquad\Omega_\delta=\{z:F_\delta(z)<2\delta^2\}.
\tag{A14}
\]

Squared distance to a nonempty closed convex set is convex.
Convolution preserves convexity and makes `F_δ` smooth.  It is
coercive since `C` is bounded.  For `z∈C`, `F_δ(z)≤δ²`, while
`d(z,C)≥3δ` implies `F_δ(z)≥4δ²`.  Consequently

\[
 C\subset\Omega_\delta\subset\{z:d(z,C)<3\delta\}.
\tag{A15}
\]

The level `F_δ=2δ²` is regular: a vanishing gradient of a differentiable
convex function is a global minimizer, whereas here the minimum is at
most `δ²`.  Thus `Ω_δ` is bounded, convex and smoothly bounded.  This
construction applies unchanged when `C` is a line segment or singleton.

For every fixed polynomial `p`, continuity on a fixed compact
neighborhood of `C` and (A15) give

\[
 \lim_{\delta\downarrow0}\max_{\overline{\Omega_\delta}}|p|
 =\max_C|p|.
\tag{A16}
\]

Apply (A13) first, then take this scalar limit.  The result is the
requested inequality for every finite dimension and every polynomial.
If the limiting maximum is zero, the same limit gives `p(A)=0`;
no division by that limiting maximum occurs.

### Convexity need not be imported as a black box

Compactness of `W(A)` follows from compactness of the unit sphere.
For convexity, any two numerical-range values come from the compression
of `A` to the span of their two unit vectors.  For a two-dimensional
compression, write the matrix in the Pauli basis.  Its numerical range
is the image of the real unit sphere `S²⊂R³` under an affine real-linear
map into `R²` (the pure-state Bloch representation).  A linear map
`R³→R²` has nonzero kernel; its image of the unit sphere equals its
image of the closed unit ball, by adding a suitable kernel vector to
the minimum-norm preimage.  That image is convex.  The one-dimensional
compression is a singleton.  Therefore the segment between the
original two numerical-range values belongs to `W(A)`.

## 4. Adversarial coverage table

| Potential gap | Audit result |
| --- | --- |
| Noncommuting rearrangement in the recurrence | Only `[E_j,T]=0` is used; the singular-vector substitution is explicitly ordered in (A1). |
| Tail term could grow too quickly | `||T^j||≤2+M` follows immediately from the defining correction identity; (A3)'s tail vanishes. |
| Circular boundedness of `g(A)` | Direct contour estimate (A12), without any spectral-set inequality. |
| Cauchy transform of continuous boundary data might not have the asserted regularity | Removed: only the finite operator integral in (A11) is needed. |
| Hidden positivity/complete positivity theorem | The concrete positive kernel (A8) and `L²` construction suffice. |
| Dimension dependence | The auxiliary bound `M` controls the terminal term on each fixed contour. The nonnegative scalar `D` survives in (A6), whose sign contradiction is valid for every `D≥0`; the final constant is exactly two in all dimensions. |
| Nonnormal, defective, repeated eigenvalues | Resolvents and singular vectors remain available; no diagonalization or simple-spectrum assumption. |
| Flat portions, corners, empty interior of `W(A)` | Explicit smooth outer approximation (A14)–(A16). |
| Norm normalization fails when `max_{W(A)}|p|=0` | Normalize on outer domains and pass to the limit instead. |
| Infinite-dimensional reduction in the source | Not needed for `goal.txt`, whose matrices are finite dimensional. No audit claim about that extra scope here. |
| Complete spectral-set conclusion | Not claimed; scalar commutation is the exact limitation. |
| Old local-Hessian factor debts | Superseded for the scalar goal, not separately proved by this argument. |

## 5. Reproducible algebra diagnostic

`experiments/global_power_lemma_exact_audit_20260905.py` checks the
recurrence, square-completion identity, contraction estimates and finite
iteration with exact rational arithmetic on a nonnormal example with
nonzero commuting corrections.  It also verifies counterexamples when
either boundedness or commutation is removed.  These checks are
diagnostics; the general argument is Sections 1–3, not extrapolation
from finitely many powers.

The diagnostic's exact nonnormal example is

\[
 Q=\begin{pmatrix}1/3&4/5\\0&1/3\end{pmatrix},\quad
 T=\begin{pmatrix}1/3&8/5\\0&1/3\end{pmatrix},\quad
 V=I,\quad x=\frac{(1,5)^T}{\sqrt{26}},\quad\kappa=5/3.
\]

Here `I−Q*Q` is positive definite, and the binomial formula for a
square-zero off-diagonal matrix proves `E_j=3^{-j}I` for **all** `j`.
Thus the test is not merely approximately admissible.

## 6. Comparison with the campaign's final polynomial/polygon proof

After finishing Sections 1–5 independently, I read
`proof/CROUZEIX_PROOF.md` in full.  Its Sections 1–5 give a complete proof
of the actual `goal.txt` scalar matrix inequality; no local Hessian result
or unproved claim is used.  I checked the following independently:

1. Its equations (3)–(6) have the same signs, geometric-series coefficient,
   and terminal bound as (A1)–(A6) above.  In particular multiplying the
   recurrence by `k^{-j}` produces exactly the displayed finite telescoping
   formula, not a missing factor of `k`.
2. Its finite supporting-half-plane construction gives a full-dimensional
   polygon containing the compact convex set strictly, including the
   singleton and segment cases.  The finite exclusion cover proves
   containment in the prescribed open neighborhood.  This is a valid
   alternative to my independent smooth-distance construction.
3. Positivity holds on every open edge.  Ignoring vertices in the arclength
   integral is legitimate; the density stays bounded since the whole
   contour is separated from the finite spectrum.
4. The normal identity `dσ=iν ds` and the adjoint half of the resolvent
   formula produce precisely its equation (11), with no conjugation or
   factor-of-two error.
5. Its all-power boundedness estimate and commutation are direct and
   noncircular.  Contour shrinking is performed only after proving the
   uniform final constant on each fixed polygon.
6. Its final sharpness computation has supremum one and matrix norm two,
   as required; equality is not mislabeled as a counterexample.

**Explicit audit judgment:** the inspected final proof establishes the
finite-dimensional scalar theorem in the stated full generality.  I found
no missing hypothesis, invalid noncommutative step, unproved all-dimension
reduction, or uncontrolled boundary/limit operation.  The stronger
matrix-valued polynomial statement remains outside its claim.

Final audited `proof/CROUZEIX_PROOF.md` SHA-256:
`2cdacb8a05bd4a57c542ee4f8761b78f95f09dee7285a99a76f35fe5441fd500`.
This includes the final precision edits requiring the compact convex set
to be nonempty and distinguishing the vanishing-tail bound `M` from the
surviving nonnegative error square `D`.

Executed diagnostic command:

```bash
.venv/bin/python experiments/global_power_lemma_exact_audit_20260905.py \
  --output experiments/global_power_lemma_exact_audit_20260905.json
.venv/bin/python -m py_compile \
  experiments/global_power_lemma_exact_audit_20260905.py
```

Both succeeded.  The JSON reports all checks passed with exact rational
arithmetic, including eight nontrivial recurrence steps.
