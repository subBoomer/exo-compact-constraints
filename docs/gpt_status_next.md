# Project Status Brief — Exotic Compact Object Constraints

**Repository**: 30 source modules, 7 simulations, 10 figures
**Stack**: numpy, scipy, matplotlib, astropy, sympy

---

## Current Results

### Comparison Table (all objects at M=1 geometric units)

| Object | d_sh/M | χ² | Exclusion |
|---|---|---|---|
| Schwarzschild | 10.39 | 0.00 | — |
| Kerr (a=0.7) | 10.05 | 0.14 | — |
| MT Wormhole (b0=2M) | **4.00** | **75.80** | **Excluded** |
| RotatingWH (Teo, b0=2M) | **5.44** | **22.74** | **Excluded** |
| Boson Star (C=0.12) | 10.39* | 0.00* | *nonphysical shadow |
| Gravastar (r1=3M) | 10.39 | 0.00 | — |
| ECO (R=0.5) | 10.05 | 0.30 | — |

### Throat-radius exclusion map (from Phase 4)
Wormholes with b0/M < 4.5 are excluded by EHT at >95% CL.
EHT-consistent range: b0/M ∈ [4.5, 5.8].

### LIGO ringdown constraints (from Phase A)
- δf/f < 0.05, δτ/τ < 0.2 (GW150914, 90% CI)
- Only wormholes show large deviations (MT WH: δω/ω ≈ 0.56)

---

## Physics Bugs Fixed (with your help)

| # | Bug | Fix |
|---|---|---|
| 1 | Phi-omega rescaling sign | `omega_phys = omega_g * exp(-Phi_inf)` [was +] |
| 2 | Energy density 1/2 factor | Complex U(1) scalar convention |
| 3 | Gauge omega > m_s | Required for sigma_2 < 0 at center |
| 4 | Collapse stopping condition | Compactified coordinate x = r/(1+r) |
| 5 | Christoffel RK4 divergence | Hamiltonian formulation for geodesics |
| 6 | Teo shadow formula | `b_c = e * b0` (confirmed by ray tracing) |

---

## Remaining Issues

### 1. Boson star EKG solver (experimental)
- Compactified-coordinate solver works for 10/11 sigma_c values
- alpha(0) < 1 (physical redshift) for all solutions ✓
- M_max ≈ 1.1 vs expected Kaup 0.633 — factor ~2 from fixed-step RK4
- Remaining issue: secondary (excited) branch, not ground state
- Needs adaptive mesh or better integrator

### 2. Teo wormhole ray tracing (Hamiltonian)
- Works for Schwarzschild (verified: boundary at 10.2M vs 10.39M, error < 2%)
- Works for Teo (verified: b_c = e*b0, agrees with analytic prediction)
- BUT: 2D image-plane scan is slow (~2s per ray, need ~50×50 = 2500 rays)
- Equatorial scan is fast (20 rays, ~5s) but gives only approximate boundary

### 3. Boson star production metric (tanh approximation)
- Calibrated tanh profile uses p = 2+5C
- Shadows are nonphysical — flagged in docstring
- Need full EKG solution for accurate shadows

---

## Next Steps — Priority Order

### Step 1: Accelerate Hamiltonian ray tracer (highest impact)
**Goal**: Make 2D image-plane scans feasible (target: 50×50 pixels in < 1 minute)

**Options**:
- (a) Cache metric values — pre-compute g_μν on a grid and interpolate instead of evaluating at each step
- (b) Use lower precision for non-critical rays — RK4 with h=1.0 instead of 0.5
- (c) Parallelize over alpha-beta grid
- (d) Implement logarithmic radial coordinate (fewer steps in weak field)

**Question**: Which approach is best for this use case?

### Step 2: 2D Teo wormhole shadow boundary
**Goal**: Full (α, β) shadow image for the Teo wormhole with rotation.

**Plan**:
- 50×50 scan of the image plane (α ∈ [-10, 10], β ∈ [-10, 10])
- Binary classification map (escape vs throat-cross)
- Edge extraction for the shadow boundary
- Compare with Kerr shadows at same spin

### Step 3: Wormhole parameter scan with correct Teo shadow
**Goal**: Update the exclusion map with the correct `b_c = e * b0` formula.

**Plan**:
- Re-run throat-radius scan with `b_c = e * b0` for Teo wormhole
- Compute χ² against EHT for the corrected shadow
- Find the b0 range that's consistent with EHT (currently b0/M ∈ [4.5, 5.8] under old approx)

### Step 4: Publication pipeline
**Goal**: Paper-ready figures and tables.

**Plan**:
- Consistent color scheme across all figures
- Legend with parameter values (spin, throat radius, compactness)
- Multi-panel summary figure (shadow overlay + ringdown + exclusion)
- LaTeX table auto-generation (already partially done)

### Step 5: Boson star EKG solver refinement (optional)
**Goal**: Converge to Kaup sequence with M_max ≈ 0.633.

**Plan**:
- Replace fixed-step RK4 with scipy's solve_ivp (BDF/Radau for stiff ODEs)
- Or implement a proper shooting-on-mass approach
- This is a standalone research project — lower priority

---

## Questions for GPT

1. **Ray tracing acceleration**: We need to do 50×50 = 2500 pixel scans. Currently ~2s per pixel. Best strategy?
   - Cache metric on a grid?
   - Use logarithmic r coordinate?
   - Or is there a smarter boundary-finding algorithm that doesn't need full scans?

2. **Teo shadow with spin**: The Teo wormhole with frame-dragging (a > 0) produces an asymmetric shadow. Is the Hamiltonian formulation sufficient, or does the g_tp coupling cause issues near the throat?

3. **Publication figure**: If we had to produce ONE impactful figure from this framework, what would you recommend? The wormhole exclusion map with the corrected Teo shadow?

4. **EKG solver convergence**: Is the Kaup limit factor-2 error worth fixing, or is the tanh approximation good enough for a first paper?
