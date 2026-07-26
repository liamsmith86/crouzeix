# A nearby scalar equality model removes exactly one Crabb copy

> **Campaign scope.**  L326 proves that a Gau--Wu equality summand
> near \(C_p\otimes I_m\) has dimension \(pr\).  This note proves the
> sharper local fact \(r=1\).  The summand may still be a general
> finite-Blaschke disk model; it is not asserted to be nilpotent or
> monomial.

## 1. Result (L328, 2026-07-26)

Fix \(p=L+1\ge3\) and \(m\ge1\).  After shrinking the neighbourhood
in L326, every scalar equality

\[
 \|f(T)x\|=2,\qquad
 w(T)\le1,\qquad
 T\ \hbox{near }C_p\otimes I_m                    \tag{1}
\]

has a Gau--Wu reducing model summand of dimension exactly \(p\):

\[
 T\cong A_f\oplus T',
\qquad
 A_f=X_\phi S(\phi)X_\phi^{-1},
\qquad
 \dim H(\phi)=p.                                   \tag{2}
\]

Consequently \(A_f\) is a single scalar disk-equality block near
\(C_p\), while

\[
 T'\ \hbox{is near }C_p\otimes I_{m-1}.            \tag{3}
\]

Thus every exact terminal equality lowers copy multiplicity by
**one**, not merely by an unspecified positive number of copies.
L192's single-block tube applies to the first summand, and induction
applies to the complement.

## 2. Exclude a larger irreducible model by compactness

Suppose the conclusion were false.  L326 then gives a sequence

\[
 T_n\longrightarrow C_p\otimes I_m
\]

whose Gau--Wu summands \(A_n\) have dimension

\[
 d=pr,\qquad r\ge2,                                \tag{4}
\]

after taking a subsequence.  Let \(Q_n\) be their reducing
projections.  Passing to another subsequence and using L326's local
commutant argument gives

\[
 Q_n\longrightarrow Q,\qquad
 (C_p\otimes I_m)|_{\operatorname{ran}Q}
 \cong C_p\otimes I_r.                             \tag{5}
\]

After near-identity unitary alignment of the subspaces,

\[
 A_n\longrightarrow C_p\otimes I_r                \tag{6}
\]

up to unitary equivalence.

## 3. The Gau--Wu matrix limit is one long Crabb chain

Write \(\phi_n=zf_n\).  It is a finite Blaschke product of degree
\(d\), and the eigenvalues of \(S(\phi_n)\), hence of \(A_n\), are
the zeros of \(\phi_n\), counted with multiplicity.  Since every
eigenvalue of \(T_n\) tends to zero, all \(d\) zeros of \(\phi_n\)
tend to zero.

Gau--Wu Corollary 3 gives an explicit upper-triangular matrix for
\(A_n\).  If its ordered zeros are \(a_1,\ldots,a_d\), its basic
entries are

\[
 b_{ij}=(-1)^{j-i-1}
 \overline a_{i+1}\cdots\overline a_{j-1}
 \sqrt{(1-|a_i|^2)(1-|a_j|^2)}.                   \tag{7}
\]

The endpoint-adjacent entries are multiplied by \(\sqrt2\), the
interior adjacent entries are not, and the corner has factor \(2\).
As every \(a_i\to0\),

\[
\begin{aligned}
 b_{i,i+1}&\longrightarrow1,\\
 b_{ij}&\longrightarrow0\qquad(j\ge i+2).
\end{aligned}
\]

Therefore the explicit Gau--Wu matrices converge to the single
length-\(d\) Crabb shift:

\[
 A_n\longrightarrow C_d=C_{pr}                   \tag{8}
\]

up to unitary equivalence.  Equations (6) and (8) would force

\[
 C_{pr}\cong C_p\otimes I_r.                       \tag{9}
\]

## 4. An exact power invariant gives the contradiction

Take the power \(p-1\).  The repeated short chains traverse both
Crabb endpoint weights, whereas a segment of length \(p-1<pr-1\)
inside the long chain can traverse at most one.  Hence

\[
\begin{aligned}
\|(C_p\otimes I_r)^{p-1}\|^2&=4,\\
\|C_{pr}^{p-1}\|^2&=2
\qquad(r\ge2).                                    \tag{10}
\end{aligned}
\]

Unitary equivalence preserves power norms, so (10) contradicts (9).
This proves \(r=1\).

The argument does not use nilpotency of the nearby \(T_n\).  It uses
nilpotency only for the **limit matrices** obtained after all
Blaschke zeros tend to zero.  Thus it is compatible with L326's
scope correction: the nearby \(f_n\) remain general finite Blaschke
products of degree \(p-1\).

## 5. Consequence for the local induction

At an exact zero face, L326 supplies the reducing model.  L328 now
places that model inside L192's fixed single-\(C_p\) tube and leaves
exactly the multiplicity-\((m-1)\) problem on the complement.
Repeating terminates after at most \(m\) steps.

There is no need for a new local theorem around a hypothetical
irreducible \(pr\)-dimensional disk model approaching \(r\) repeated
short chains: (10) proves that such a model cannot occur.

## 6. Audit

Run

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/repeated_crabb_model_size_rigidity.py \
  --output \
  experiments/repeated_crabb_model_size_rigidity_s70224.jsonl
```

The exact checker verifies (10) for \(3\le p\le7\) and
\(r\in\{2,3\}\).  The all-\(p,r\) proof is the consecutive-weight
argument in Section 4; the finite computation is an adversarial
audit of its orientation and exponent.
