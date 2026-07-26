# Lossless feedback transports both defect operators exactly

## 1. Result (L267, 2026-07-25)

Let

\[
{\cal U}=
\begin{bmatrix}A&B\\ C&D\end{bmatrix}:
{\cal E}\oplus{\cal X}\longrightarrow
{\cal E}\oplus{\cal X}
\]

be unitary.  Let \(Z\) be any operator on \({\cal X}\) for which the
formal inverses below exist, and form the Redheffer feedback

\[
\Phi_{\cal U}(Z)
=A+BZ(I-DZ)^{-1}C.                               \tag{1}
\]

Then its initial and final defects are exact congruences of the two
defects of \(Z\):

\[
\boxed{
\begin{aligned}
I-\Phi_{\cal U}(Z)^*\Phi_{\cal U}(Z)
={}&C^*(I-Z^*D^*)^{-1}\\
&\quad (I-Z^*Z)(I-DZ)^{-1}C,
\end{aligned}}                                   \tag{2}
\]

\[
\boxed{
\begin{aligned}
I-\Phi_{\cal U}(Z)\Phi_{\cal U}(Z)^*
={}&B(I-ZD)^{-1}\\
&\quad (I-ZZ^*)(I-D^*Z^*)^{-1}B^*.
\end{aligned}}                                   \tag{3}
\]

No commutativity is used.  The identities hold for bounded
contractions and coefficientwise for formal operator series with
invertible constant denominators.

Thus a genuinely lossless scalar entry/return network cannot insert
an unknown positive weight into a tail defect.  It only applies the
explicit entrance column and its adjoint.  This is the abstract
network identity sought in A194.

L267 does **not** close A194.  L268 subsequently proves that the full
one-delay closed-return quotient cannot be obtained by transporting
the deflated quotient as the sole non-background defect through one
analytic lossless port: their exact first faces have ranks two and
one, despite equal traces.  Functional calculus does not commute with
Redheffer feedback in general, and L251's common channel is
insufficient.  L267 may organize individual subchannels, but the
complete physical trace still requires an additional cancellation.

## 2. Initial-defect proof

Put

\[
R=(I-DZ)^{-1}.
\]

The resolvent identity gives

\[
R=I+DZR,\qquad ZR=(I-ZD)^{-1}Z.                  \tag{4}
\]

Unitarity of \({\cal U}\) gives

\[
\begin{aligned}
I-A^*A&=C^*C,& A^*B&=-C^*D,\\
I-B^*B&=D^*D,& B^*A&=-D^*C.
\end{aligned}}                                   \tag{5}
\]

Expand the left side of (2) using (1):

\[
\begin{aligned}
I-\Phi^*\Phi
={}&I-A^*A
-A^*BZR-C^*R^*Z^*B^*A\\
&-C^*R^*Z^*B^*BZR C.
\end{aligned}}                                   \tag{6}
\]

Substitute (5) into (6), factor \(C^*\) and \(C\), and use (4):

\[
\begin{aligned}
I-\Phi^*\Phi
={}&C^*\{I+DZR+R^*Z^*D^*\\
&\qquad-R^*Z^*(I-D^*D)ZR\}C\\
={}&C^*R^*(I-Z^*Z)RC.
\end{aligned}
\]

Since \(R^*=(I-Z^*D^*)^{-1}\), this is (2).

## 3. Final-defect proof

Apply (2) to the adjoint unitary

\[
{\cal U}^*=
\begin{bmatrix}A^*&C^*\\ B^*&D^*\end{bmatrix}
\]

and the feedback variable \(Z^*\).  The resulting feedback is
\(\Phi_{\cal U}(Z)^*\).  Taking adjoints and using

\[
Z(I-DZ)^{-1}=(I-ZD)^{-1}Z
\]

gives exactly (3).

## 4. Interface with the delayed campaign

The desired use is an associated-order factorization

\[
\widetilde A_{\rm full}
=\Phi_{{\cal U}_{r,c}}(\widetilde A_{\rm tail})
+O(c^{2k+1}),                                    \tag{7}
\]

where \({\cal U}_{r,c}\) is a universal copy-scalar lossless port and
the error is too late for the active trace.  Equations (2)--(3) would
then transport L244's half-line cancellation and make L261's Hardy
return automatic.  L245 fixes the remote amplitude and L241 fixes
the resulting scalar normalization.

Equation (7), with the deflated quotient as the sole feedback defect,
is false by L268's exact rank obstruction.  A larger lossless network
could evade that obstruction only by adding another independent
defect channel.  Such a channel would itself have to be evaluated and
cancelled, so it would not reduce A194 to L267.  The live route is
therefore A213's trace-ideal calculation, not a search for another
one-shot port.

The identity itself is classical in conservative systems/Redheffer
theory.  Related shorted-defect formulas appear in Y. Arlinskii,
*Schur parameters, Toeplitz matrices, and Krein shorted operators*,
arXiv:1109.4020.

## 5. Independent numerical audit

`experiments/lossless_redheffer_defect_transport.py` generates
deterministic noncommuting feedback matrices for external/state
dimensions `(1,2)`, `(2,3)`, and `(3,2)`.  For each size it samples
twenty complex feedbacks and evaluates both sides of (2) and (3)
independently.  The largest Frobenius residual was
`1.433e-15`.

The tracked output is
`experiments/lossless_redheffer_defect_transport_s70225.jsonl`, with
SHA-256
`972f9f2516914af778a637a192f3dfb0bc42bac83cda3201893be00a4e1809a3`.
The checker passes Ruff and `py_compile`.  This audit tests the
noncommutative algebra of L267; it does not assume the physical
factorization (7), whose direct single-port form is disproved by
L268.
