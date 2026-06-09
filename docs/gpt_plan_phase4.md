# Phase 4 Plan: Parameter Space Exploration + Exclusion Contours

## Current state (after Phase 3 GPT review)

All 7 metrics implemented and verified:
- Schwarzschild, Kerr — GR baselines
- MT Wormhole, Teo Wormhole — throat-radius shadows
- BosonStar — tanh approximation (shadows nonphysical, flagged)
- Gravastar — de Sitter interior + tanh shell
- ECO — Kerr + quantum corrections

## Phase 4 objective

Systematically scan parameter space to answer:

1. **What wormhole throat radius is consistent with EHT?**
2. **How does shadow asymmetry evolve with spin?**
3. **What (mass, spin) combinations are excluded by current constraints?**
4. **Are there degenerate regions where different objects give the same observables?**

---

## Calculations to perform

### Scan 1: Wormhole throat radius

Parameter: b0 ∈ [1.5M, 10M], step 0.25M
Objects: MT Wormhole, RotatingWormhole (Teo)
Observables: d_sh, chi^2 against EHT (42 ± 3 μas for M87*)

Output: b0 vs chi^2 plot, threshold b0 where chi^2 < 2

### Scan 2: Spin

Parameter: a ∈ [0, 0.99], step 0.05
Objects: Kerr, ECO (R=0.5), RotatingWormhole (Teo, b0=5M)
Observables: d_sh/M, asymmetry A, ISCO/M

Output: Spin vs d_sh + A overlay plot

### Scan 3: Combined exclusion

Parameters: M ∈ [1, 100] Msun, a ∈ [0, 0.99]
Objects: Kerr, ECO (R=0.3, 0.7)
Observables: chi^2 composite (shadow + ringdown + tidal)

Output: 2D exclusion contour plot

### Scan 4: Degeneracy map

For each pair of objects, find parameter combinations where
|observable_A - observable_B| < measurement uncertainty

---

## Questions for GPT

1. **Wormhole throat scan**: Our current MT wormhole uses b(r) = b0²/r and Phi=0. For a wormhole with ADM mass M, what's the relationship between b0 and M? Can b0 be arbitrarily larger than M, or is there a physical upper bound?

2. **Exclusion contours**: We plan to use chi^2 = Σ[(O_model - O_obs)/sigma_obs]² against EHT (shadow diameter + asymmetry) and LIGO (ringdown frequency + tidal deformability). Are these the right observables for joint constraints? Should we add any others (e.g., photon ring autocorrelation)?

3. **Degeneracy detection**: Are there known parameter combinations where different object classes produce identical shadow + ringdown signatures? For example, a Kerr BH with spin a=0.9 vs an ECO with R=0.3 at the same mass?

4. **Mass scaling**: All geometric observables (d_sh/M, ISCO/M, precession) are dimensionless in M. The only M-dependent observable is f_RD ∝ 1/M. For the exclusion contours, should we:
   a) Fix M=10 Msun and vary a (current approach)
   b) Vary both M and a independently
   c) Use M87* mass (6.5e9 Msun) for shadow constraints and 10 Msun for ringdown

5. **EHT measurement uncertainties**: We currently use d_sh = 42 ± 3 μas for M87*. But the actual EHT constraint is asymmetric and depends on spin prior. Should we use:
   a) Symmetric 42 ± 3 μas (current)
   b) Asymmetric 39-45 μas (EHT 68% CI)
   c) The full likelihood approximation from EHT Paper VI

6. **Most actionable question**: Based on our framework so far, what single calculation would be most valuable for a first paper figure?

---

## Implementation plan

```
simulations/
  scan_throat_radius.py       # Scan 1
  scan_spin.py                 # Scan 2  
  scan_exclusion.py            # Scan 3
  scan_degeneracy.py           # Scan 4

visualizations/
  exclusion_throat_radius.png
  exclusion_spin.png
  exclusion_mass_spin.png
  degeneracy_map.png

docs/
  exclusion_summary.md
```
