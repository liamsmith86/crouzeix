# sym4 program (Epoch 5) — first beyond-literature target

Family: A 4×4 zero-diag tridiagonal real, weights (a₁,b₁,a₂,b₂,a₃,b₃); σ(A) = {±e₁, ±e₂};
A⁴ = pA² − qI; symmetries A ~ −A (D = diag(1,−1,1,−1)) and A real.

## Proved en route (2026-07-20)
**Sym3 ellipse fact — analytic proof:** for the 3×3 checkerboard pattern, det(xH+yS) ≡ 0
(rows 1,3 of xH+yS have only the middle entry ⟹ proportional-degenerate), so the Kippenhahn
determinant factors t·(t² + c₂(x,y)): dual = line (⟷ point 0 = center eigenvalue) ∪ conic
(⟷ ellipse, foci ±e). Hence W(A) = that ellipse exactly. [Elementary + standard Kippenhahn.]
Sym4 analog: det(xH+yS+tI) is EVEN in t: t⁴ + c₂(x,y)t² + c₄(x,y) — a conic in t²; the
boundary-generating curve is a "bi-conic"; expect the ζ = z²-image structure to be the hull of a
tractable curve (work out — the ζ-domain is NOT an ellipse, per probe residuals 1e-4–1e-3).

## Graded collapse (derived; verify via sweep)
Sublattice split odd {1,3} / even {2,4}: A = [[0, A₊],[A₋, 0]], A₊ = [[a₁,0],[b₂,a₃]],
A₋ = [[b₁,a₂],[0,b₃]]; N_odd = A₊A₋, N_even = A₋A₊ (2×2, eigenvalues {e₁², e₂²} each).
Extremal f₀ odd (zeros {0, ±α} — probe-confirmed in the real-pair sector): f₀ = z·F(z²)-type,
f₀(A) = [[0, A₊F(N_even)],[A₋F(N_odd), 0]] ⟹ **K = max(‖A₊F(N_even)‖, ‖A₋F(N_odd)‖)** — two
coupled GKL-type 2×2 problems; F(e_j²) = B(τ_j)/e_j, τ_j = φ(e_j).
ρ: f₀g₀ even ⟹ (f₀g₀)(A) block-diagonal; x₀ in one sublattice:
**ρ = Re[Σ_j B(τ_j)·g₀(e_j)·q_j]**, q_j = ⟨Q_j x₀,x₀⟩ (block spectral projections), Σq_j = 1.
Sym3 analogs to hunt in sweep data: stationarity law for α (midpoint-of-{τ₁²,τ₂²}-analog?),
q_j-law (π₀ = 1/2-analog), sign law for B(τ_j)g₀(e_j).

## Sector caveats
- Imaginary-eigenvalue-pair sector (probe case C) behaves differently (boundary-zero artifacts);
  treat separately.
- Degenerate corners: e₂ → 0 (collapses toward sym3-like), e₁ → e₂ (double pair).

## Data: experiments/sym4_sweep.py → sym4_sweep_s61.jsonl (running).
