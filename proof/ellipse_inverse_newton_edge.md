# The inverse ellipse map has a Catalan Newton edge (2026-07-23)

## 1. Statement

For real `c` near zero, let

\[
 E_c=\{w+c\overline w:|w|<1\}
\]

and let `Psi_c` be the inverse of the centered normalized Riemann map
from `E_c` to the disk.  Write

\[
 \Psi_c(z)=\sum_{n\geq0}b_n(c)z^{2n+1}.
\]

Then, for every fixed `n`,

\[
 \boxed{b_n(c)=C_n c^n+O(c^{n+2}),}                  \tag{1}
\]

where

\[
 C_n={1\over n+1}{2n\choose n}
\]

is the `n`th Catalan number.  More precisely,

\[
 b_n(c)=c^n h_n(c^2)
\]

for a real-analytic function `h_n` with `h_n(0)=C_n`.  Thus the
lowest bivariate edge is

\[
 \boxed{
 \Psi_c(z)=
 z+cz^3+2c^2z^5+5c^3z^7+14c^4z^9+\cdots
 }
                                                               \tag{2}
\]

and every correction to the displayed coefficient at circle grade
`n` contains at least two additional powers of `c`.

This is a scalar conformal-map lemma.  It supplies the graded
propagation rule anticipated in A84, but does **not** by itself prove
the optimized Stein deficit asserted there.

## 2. Exact differential equation

Put `q=c^2`, let `k=k(q)` be Jacobi's elliptic modulus, and write

\[
 \alpha={\pi\over2K(k^2)}.
\]

The classical inverse ellipse map is

\[
 \Psi_c(z)=2\sqrt c\,
 \sin\!\left[
   \alpha F\!\left(\arcsin {z\over\sqrt k}\,\middle|\,k^2\right)
 \right].                                                   \tag{3}
\]

Differentiating (3), then eliminating the sine and cosine, gives

\[
 (\Psi_c')^2(k-z^2)(1-kz^2)
 =\alpha^2(4c-\Psi_c^2).                                  \tag{4}
\]

There is no choice of elliptic-function convention left in (4);
the normalization is fixed by

\[
 \Psi_c'(0)=2\alpha\sqrt{c/k}>0.                           \tag{5}
\]

Jacobi's theta identities make all relevant parameter dependence
regular:

\[
\begin{aligned}
a(c)&=\sum_{j\geq0}c^{2j(j+1)},\\
t(c)&=1+2\sum_{j\geq1}c^{2j^2},\\
{k\over c}&={4a(c)^2\over t(c)^2},\qquad
\alpha={1\over t(c)^2},\\
b_0(c)&={1\over a(c)t(c)}.
\end{aligned}                                             \tag{6}
\]

In particular, `k/c`, `alpha`, and `b_0` are analytic functions of
`c^2`, with values `4`, `1`, and `1` at zero.

## 3. Newton blow-up

Set

\[
 \Psi_c(z)=z\,h(c,z^2),\qquad
 h(c,x)=H(c,cx).
\]

If `y=cx` and `D=H+2yH_y`, multiplying (4) by `c` gives

\[
 D^2\left[
 ck-(1+k^2)y+{k\over c}y^2
 \right]
 =\alpha^2(4c^2-yH^2).                                  \tag{7}
\]

Every coefficient in (7) is analytic in `(c^2,y)`.  Starting from
`H(c,0)=b_0(c)`, comparison of successive powers of `y` uniquely
determines a formal series

\[
 H(c,y)=\sum_{n\geq0}h_n(c^2)y^n.                       \tag{8}
\]

The same coefficient recursion follows from the analytic solution
(3), so the formal series converges in a neighbourhood of the origin.
Substitution into `h(c,x)=H(c,cx)` proves the divisibility and parity
claim in (1).

At `c=0`, equation (7), after cancelling its common factor `-y`,
becomes

\[
 (H_0+2yH_0')^2(1-4y)=H_0^2,\qquad H_0(0)=1.           \tag{9}
\]

The normalized branch satisfies

\[
 H_0+2yH_0'={H_0\over\sqrt{1-4y}}.
\]

Its unique solution is

\[
 H_0(y)={2\over1+\sqrt{1-4y}}
        ={1-\sqrt{1-4y}\over2y}
        =\sum_{n\geq0}C_ny^n.                          \tag{10}
\]

Equations (8)--(10) prove (1)--(2).

## 4. Relevance to the Crabb elliptic face

In L123's companion coordinates, a coefficient at offset `k` must
propagate through `k` odd functional-calculus steps before it can
reach the endpoint defect.  Equation (1) makes that filtration exact:
the first available scalar coefficient is `C_k c^k`, and no term of
lower elliptic weight exists.

This proves the order-selection part of A84's proposed path argument.
What remains is substantial:

1. linearize `T(a,c)=phi_c(S(a,c))` in the equality amplitude;
2. solve the optimized rank-one Stein recurrence to the same grade;
3. prove that its two endpoint paths combine to
   `-64|u_k|^2c^(2k)`; and
4. control interactions and analytic remainders uniformly.

The Catalan coefficient need not survive separately in the final
answer: metric optimization and the direct `cA^\dagger` term can
cancel intermediate Catalan paths.

## 5. Exact regeneration

Run

```bash
.venv/bin/python -u experiments/ellipse_inverse_newton_edge.py \
  --output experiments/ellipse_inverse_newton_edge_s70223.jsonl
```

The script uses rational theta series and solves (4) coefficient by
coefficient.  It checks the differential-equation residual, the
`c^n` divisibility, the parity gap, and the Catalan leading
coefficients through the requested degree.
