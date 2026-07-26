# Joint scalar-channel and disk-residual response factorization

> **Campaign scope.**  This is a curvewise, fixed-dimension statement
> on the zero-Jensen branch near one repeated Crabb block.  It repairs
> an overbroad downstream use of L323: transfer leakage alone does not
> control a response which also contains off-channel disk residuals.
> The result does not yet assert the full repeated-block neighbourhood
> theorem.

## 1. Result (L324, 2026-07-26)

Fix a Crabb length \(L\), copy multiplicity \(m\), and a
real-analytic arc in L194's repeated disk/normal chart.  After finite
ramification, choose analytic input/output scalar-channel frames
\(v(s),u(s)\) for L201's equality-anchor transfer.  Put

\[
 \delta(s)
 =1-\sum_n|u(s)^*B_n(s)v(s)|^2.                  \tag{1}
\]

Let \(d(s)\ge0\) be the **accumulated scalar disk loss** on the same
channel after L197's exact Schur congruences.  Equivalently, if
\(\eta_0,\ldots,\eta_{r-1}\) are the active analytic Hardy quotient
columns seen by the transported channel vector, then, up to uniformly
positive analytic weights,

\[
 \boxed{d(s)\asymp\sum_{j=0}^{r-1}\|\eta_j(s)\|^2.} \tag{2}
\]

The terminal zero block is omitted from (2).  It contains no
unpaid Hardy quotient column.

Let \({\cal R}\) be any scalar prepared disk/normal or
disk/reflected response coefficient on L199's zero-Jensen winner.
Subtract the response of the channel-preserving scalar model
controlled by L192.  Then

\[
 \boxed{
 |{\cal R}_{\rm full}(s)-{\cal R}_{\rm split}(s)|
 \le C\{\delta(s)+d(s)\}.}                         \tag{3}
\]

Here \(C\) is uniform on one sufficiently small fixed chart.  In
particular, completing finitely many such responses against L192's
retained scalar normal/reflected curvatures costs

\[
 \boxed{O\bigl((\delta+d)^2\bigr).}                \tag{4}
\]

After shrinking,

\[
 O\bigl((\delta+d)^2\bigr)
 \le \varepsilon_0(\delta+d).                    \tag{5}
\]

L325 subsequently supplies one joint linear reserve for these two
losses on the same actual-disk norming direction.  Choose
\(\varepsilon_0\) below its local constant; that reserve then
absorbs every later **scalar** transport correction after shrinking.
The abstract negative Schur quotient in L319 cannot occur as an
unpaid scalar response: its actual channel-breaking part must carry
either transfer leakage or one of the accumulated disk square
factors.

L324 is deliberately an accumulated-loss theorem.  Replacing \(d\)
by only the newest L197 quotient is invalid, because a response can
pass through an earlier active Hardy range.  It also does not claim
that the complete matrix endpoint is positive; L319 still blocks
that stronger inference.

## 2. Analytic square coordinates for the disk loss

Along the chosen arc, L197 gives an analytic congruence

\[
 S_{\rm disk}(s)\sim
 \operatorname{diag}\bigl(
 s^{2\nu_0}A_0(s),\ldots,
 s^{2\nu_{r-1}}A_{r-1}(s),0\bigr),
 \qquad A_j(0)\succ0.                              \tag{6}
\]

Take the positive analytic square root of each \(A_j\) and define

\[
 \eta_j(s)=s^{\nu_j}A_j(s)^{1/2}w_j(s),            \tag{7}
\]

where \(w_j\) is the channel vector transported through the preceding
triangular congruences.  The scalar compression of (6) is exactly
\(\sum_j\|\eta_j\|^2\).  Changing among the physical Hardy quotient
columns, L197's triangular columns, and (7) uses analytic invertible
maps.  Their norms are uniformly equivalent on a compact
neighbourhood, proving (2).

The first square is L195's

\[
 4\sum_{a,b}\|F_{ab}v\|^2.                         \tag{8}
\]

L196's inverse diagonal-difference formula proves that its kernel
removes both the leading channel column and row of the transverse
inverse-Gram coordinate.  Every later L197 step is the exact
least-squares quotient after the earlier Hardy ranges have been
removed.  Consequently the **collection** \((\eta_j)_j\), not merely
its last member, contains every off-channel disk variable which can
enter a later response on the selected scalar state.

This is the endpoint-specific information absent from L319's two
arbitrary Hermitian germs.

## 3. The joint defect-frame gauge

Use the channel frames to decompose the equality transfer as

\[
 B=\begin{bmatrix}b&r\\c&D\end{bmatrix}.           \tag{9}
\]

L323 proves

\[
 \|r\|_{H^2}^2=\|c\|_{H^2}^2=\delta.              \tag{10}
\]

Decompose every accumulated Hardy quotient column into its scalar,
complement, and cross-channel blocks.  The channel-preserving scalar
blocks are precisely the data retained in
\({\cal R}_{\rm split}\).  The remaining blocks have norm bounded by
a fixed multiple of

\[
 \left(\sum_j\|\eta_j\|^2\right)^{1/2}.            \tag{11}
\]

Now reverse both complementary defect frames.  This is only the
coordinate change

\[
 J_{\rm out}=1\oplus(-I),\qquad
 J_{\rm in}=1\oplus(-I).                           \tag{12}
\]

It changes the sign of:

1. both transfer cross blocks \((r,c)\); and
2. every off-channel accumulated disk quotient block.

It leaves the physical operator, scalar function, norming
functional, metric, and response unchanged.  Therefore the Taylor
series of
\({\cal R}_{\rm full}-{\cal R}_{\rm split}\) contains:

1. no constant term in these joint cross variables; and
2. no term of total cross degree one.

Every surviving monomial has total cross degree at least two.  On a
fixed finite chart, Taylor's theorem and Cauchy--Schwarz give

\[
\begin{aligned}
 |{\cal R}_{\rm full}-{\cal R}_{\rm split}|
 &\le C_0\left(
   \|r\|_{H^2}+\|c\|_{H^2}
   +\Bigl(\sum_j\|\eta_j\|^2\Bigr)^{1/2}
   \right)^2\\
 &\le C_1\left(
   \delta+\sum_j\|\eta_j\|^2\right).
                                                               \tag{13}
\end{aligned}
\]

Equations (2), (10), and (13) prove (3).

The mixed monomial
\(\sqrt\delta\,\sqrt d\) is allowed and is important.  It is charged
by

\[
 2\sqrt{\delta d}\le\delta+d,                      \tag{14}
\]

not falsely assigned to transfer leakage alone.  This is the precise
correction to the overly broad application sentence after L323.

## 4. Finite square completion

Let \(y_1,\ldots,y_N\) be the finitely many scalar true-normal and
marked reflected coordinates remaining on the zero-Jensen winner.
L192 supplies a positive definite retained curvature matrix
\(\Gamma\succeq\gamma I\) for their split-channel joint face.  After
retaining half of it, the new response has the form

\[
 -\frac12 y^*\Gamma y
+2\operatorname{Re}\langle y,R\rangle,
\qquad
\|R\|\le C(\delta+d).                              \tag{15}
\]

Completing one vector square gives

\[
 -\frac12 y^*\Gamma y
+2\operatorname{Re}\langle y,R\rangle
\le \frac2\gamma\|R\|^2
\le C_2(\delta+d)^2.                              \tag{16}
\]

Because \(0\le\delta\le1\) and \(d\to0\), shrink until

\[
 C_2(\delta+d)\le\varepsilon_0.                   \tag{17}
\]

Then (5) follows.
On the equality manifold, L322 retains at least
\((4/3)\delta\).  Off it, use L325's joint actual-disk reserve rather
than adding that number to a separately selected disk square.  L197
still identifies the accumulated disk part of L325's joint defect.

No normal curvature, channel reserve, or earlier disk square is
spent twice: all response coefficients are stacked before the one
completion (16).

## 5. What remains

L324 removes the scalar version of L319's later-response obstruction.
The remaining assembly is finite but must still be written
explicitly:

1. use L61/L199 to enter a maximal zero-Jensen winner;
2. apply (3)--(17) simultaneously with L325's joint actual-disk gap;
3. retain L318's elliptic margin in the marked reflected coordinates;
4. if \(\delta\) vanishes to the current order, follow L321's common
   eigenline; and
5. if it vanishes identically, use L205 to split a scalar channel and
   descend in copy multiplicity, with L192 as the base.

The only acceptable proof of the next theorem is this finite
curve-selected induction.  Neither abstract ideal membership nor
L197 positivity alone supplies it.

## 6. Audit

Run

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/repeated_crabb_joint_channel_disk_factor.py \
  --output \
  experiments/repeated_crabb_joint_channel_disk_factor_s70224.jsonl
```

The exact checker enumerates joint-sign-invariant response monomials,
verifies that subtraction of the split locus removes degrees zero and
one, and checks the rational Young/square-completion bounds with
multiple finite flags.  It audits the algebra in Sections 3--4; the
endpoint identification in Section 2 is L195--L197's exact Hardy
quotient construction.

The standard dataset has 15 records and SHA-256

```text
c06279887606fcff4954bdc059b8b7fd8a118fbdff4ea6b4e25c7d795eb35472
```
