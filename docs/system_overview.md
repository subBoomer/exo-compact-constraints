# System Overview: Exotic Compact Object Simulation Framework

## 1. Object Taxonomy

```
Exotic Compact Objects
├── Black Holes (GR solutions)
│   ├── Schwarzschild (static, uncharged)
│   └── Kerr (rotating)
├── Wormholes (traversable)
│   ├── Morris-Thorne (static)
│   └── Rotating Wormhole (metric-agnostic framework)
├── Exotic Compact Objects (non-BH ultra-compact)
│   ├── Boson Stars (scalar field condensates)
│   ├── Gravastars (de Sitter interior + shell)
│   └── Generic ECO (Kerr + quantum corrections)
└── [Extensions: Kerr-Newman, dilaton, naked singularities]
```

All objects are defined by their **spacetime metric** — no physics model exists outside the metric.

## 2. Metric Abstraction Layer

### Base class: `AbstractMetric` (`exo_compact/metrics/base.py`)

```
AbstractMetric
├── g_uv(x)              Covariant metric tensor g_μν
├── g_uv_contra(x)       Contravariant metric tensor g^μν
├── christoffel(x)       Christoffel symbols Γ^ρ_μν
├── geodesic_rhs(t, y)   Geodesic equation RHS for ODE solvers
├── params               Metric parameters (M, a, b0, ...)
└── name                 Human-readable identifier
```

### Concrete implementations

| Class | File | Coordinates | Parameters |
|---|---|---|---|
| `Schwarzschild` | `metrics/schwarzschild.py` | Schwarzschild (t, r, θ, φ) | mass M |
| `Kerr` | `metrics/kerr.py` | Boyer-Lindquist (t, r, θ, φ) | mass M, spin a |
| `MorrisThorne` | `metrics/morris_thorne.py` | Schwarzschild (t, r, θ, φ) | mass M, b0, shape func Φ(r), redshift func Φ(r) |
| `RotatingWormhole` | `metrics/rotating_wormhole.py` | Generalized axisymmetric | mass M, b0, spin a, pluggable N, A, B, C, ω |
| `BosonStar` | `metrics/boson_star.py` | Schwarzschild (t, r, θ, φ) | mass M, compactness C, radius R |
| `Gravastar` | `metrics/gravastar.py` | Schwarzschild (t, r, θ, φ) | mass M, interior radius r1, shell width w |
| `ECO` | `metrics/eco.py` | Boyer-Lindquist (t, r, θ, φ) | mass M, spin a, R, ε_q, λ |

### Design principle: Metric-first architecture

```
Metric ──→ Geodesic Integrator ──→ Observable Generator ──→ Comparison
  ↑                ↑                        ↑
  |                |                        |
  All physics     Works with any           Compares across
  flows from      AbstractMetric           object classes
  g_μν
```

## 3. Observable Generation Pipeline

```
                     ┌───────────────────┐
                     │ AbstractMetric    │
                     │ g_μν(x; params)    │
                     └────────┬──────────┘
                              │
              ┌───────────────┼───────────────────┐
              ▼               ▼                   ▼
    ┌─────────────────┐ ┌──────────────────┐ ┌──────────────────┐
    │NullGeodesic     │ │TimelikeGeodesic  │ │QNM Solver        │
    │Integrator       │ │Integrator        │ │(analytic fit)    │
    └────────┬────────┘ └────────┬─────────┘ └────────┬─────────┘
             │                   │                    │
             ▼                   ▼                    ▼
    ┌─────────────────┐ ┌──────────────────┐ ┌──────────────────┐
    │LensingCalculator │ │OrbitalCalculator │ │GWCalculator      │
    │- deflection α   │ │- precession Δφ   │ │- ringdown ω, τ   │
    │- Einstein ring  │ │- ISCO r, f       │ │- tidal λ         │
    │- caustics       │ │- epicyclic freq  │ │- echoes          │
    └─────────────────┘ └──────────────────┘ └──────────────────┘
             │                   │                    │
             └───────────────────┼────────────────────┘
                                 ▼
                    ┌─────────────────────────┐
                    │ Comparison Framework     │
                    │ - BH vs Wormhole vs ECO  │
                    │ - Deviation from GR      │
                    │ - Constraint ingestion   │
                    └─────────────────────────┘
```

## 4. Simulation Workflow

### Standard pipeline

```python
# 1. Choose a metric
from exo_compact.metrics.kerr import Kerr
metric = Kerr(mass=10.0, spin=0.7)

# 2. Compute observables
from exo_compact.observables.shadow import ShadowCalculator
calc = ShadowCalculator(metric)
d_sh = calc.shadow_diameter()

# 3. Compare with known constraints
from exo_compact.constraints.eht_constraints import EHTConstraints
eht = EHTConstraints()
compat = eht.check_compatibility(d_sh, 0.0)

# 4. Visualize
from exo_compact.visualization.plots import ComparisonPlotter
# ...
```

## 5. Production-Grade Ray Tracing (`exo_compact/ray_tracing/`)

New module added per GPT's architecture prescription. Uses:
- **MetricBase** interface: every spacetime implements `g(x)` and `g_inv(x)`
- **MetricCache**: precomputed (r, theta) grid with bicubic interpolation (5-15x speedup)
- **HamiltonianGeodesic**: single engine for all metrics (H = 1/2 g^μν p_μ p_ν)
- **TerminationOracle**: early stopping (escape / throat-cross / captured) for 10-50x speedup
- **ShadowScanner**: 2D N×N image-plane scans with boundary extraction

| Metric | Class | Verified |
|---|---|---|
| Schwarzschild | `SchwarzschildMetric` | Boundary at 5.24M (analytic: 5.20M, error < 1%) |
| Kerr | `KerrMetric` | Boundary shifts with spin (a=0.7) |
| Teo wormhole | `TeoMetric` | Boundary at e·b0 = 5.45M (analytic: 5.44M, error < 0.2%) |

### CLI scripts
# Lensing
python simulations/lensing/ray_trace.py --metric kerr --mass 10 --spin 0.7

# Shadow
python simulations/shadow_imaging/shadow_generator.py --metric kerr --spin 0.94

# GW ringdown
python simulations/gravitational_waves/waveform_comparison.py --metric kerr --mass 68 --spin 0.7

# Orbital dynamics
python simulations/orbital_dynamics/orbit_integrator.py --metric kerr --mass 10 --spin 0.7
```

## 5. Key Assumptions

1. **Geometric units**: G = c = 1 throughout
2. **Sign convention**: (-, +, +, +)
3. **All metrics are vacuum or effective** — no detailed matter modeling
4. **Metric-first** — every observable is derived from g_μν, never from an object model
5. **QNM frequencies** for non-Kerr objects are approximated as Kerr + fractional deviations
6. **Rotating wormholes** use a generalized ansatz — no single unique rotating wormhole metric
7. **Boson star profiles** use analytic approximations (tanh), not numerical EKG solutions

## 6. Physics Implementation Status

| ID | Physics | Pipeline location | Status | Notes |
|---|---|---|---|---|
| MP-01 | Rotating wormhole metric | `metrics/rotating_wormhole.py` | Teo (1998) metric implemented via set_teo_metric() | Shadow requires full ray tracing — no Carter constant |
| MP-02 | Boson star EKG ODE system | `metrics/boson_star_solver.py` | RK4 shooting solver (KG coef fixed: N'/N) | Still converging — masses exceed Kaup limit |
| MP-03 | Boson star production metric | `metrics/boson_star.py` | Calibrated tanh profile (p = 2 + 5C) | [WARNING] Shadows are nonphysical — no photon sphere in real BS |
| MP-04 | Exotic object QNM tables | `observables/gravitational_waves.py` | ECO echo model (Mark+2017), WH Pöschl-Teller approx | Confirmed by GPT Q4, Q5 |
| MP-05 | ECO tidal Love numbers | `observables/gravitational_waves.py` | Tunable λ parameter | Cardoso & Pani (2019) |
| MP-06 | EHT calibration pipeline | `constraints/eht_constraints.py` | Published EHT results | EHT (2019, 2022) |
| MP-07 | Kerr shadow boundary | `observables/shadow.py` | Carter constant formalism (Bardeen 1973) | Confirmed by GPT Q1 |
| MP-08 | MT wormhole shadow | `observables/shadow.py` | b_c = b0 (throat = photon sphere) | d_sh=4M for b0=2M, excluded at χ²=69.55 |
| MP-09 | Pulsar timing ECO constraints | `constraints/pulsar_constraints.py` | Placeholder entries | — |
| MP-10 | Strong deflection Bozza constant | `observables/lensing.py` | Schwarzschild: a2 = ln(216(7-4√3)) ~ −0.4002 | Bozza (2002) |
| MP-11 | Gravastar thin-shell matching | `metrics/gravastar.py` | de Sitter interior + tanh shell | r1 > 2M constraint should be relaxed to allow r1 < 2M |

## 7. GPT Physics Review (Phase 3)

Key findings from GPT review:
- Teo wormhole shadow: Kerr Carter constant fallback removed (requires ray tracing)
- Boson star shadows: tanh approximation produces Schwarzschild-like shadows that are nonphysical
- MT wormhole exclusion: robust for b0=2M; larger b0 ~ 5M would be consistent with EHT
- EKG solver: N'/(2N) coefficient derived from first principles (matches original code, not review revision)

## 8. Phase 4 Results: Parameter Space Exploration

### Throat-radius exclusion map
- MT wormholes with b0/M < 4.5 are excluded by EHT at > 95% CL
- EHT-consistent range: b0/M ∈ [4.5, 5.8]
- Above b0/M > 5.8: shadow larger than EHT but not excluded (large-throat regime)
- Full figure: `visualizations/exclusion_throat_radius.png`

### Boson star EKG solver status
- KG coefficient uses correct N'/(2N) (verified by first-principles derivation)
- RK4 shooting solver finds 23/30 sigma_c values
- Masses exceed Kaup limit by factor 2-100x — indicates remaining bug
- Root cause suspected: Phi grows positive, making A² > 1, making (ω²/A² - m²) always negative, preventing sigma decay for ω < m_s
- [SOLUTION NEEDS REVIEW] — Full EKG convergence requires further work

### Exclusion summary: `docs/exclusion_summary.md`

## 9. Repository File Summary

| Directory | Files | Status |
|---|---|---|
| `exo_compact/metrics/` | 7 metrics + 1 solver | 7 production, 1 experimental |
| `exo_compact/observables/` | 4 calculators | Verified |
| `exo_compact/geodesics/` | 2 integrators | Stub (needs full ray tracing) |
| `exo_compact/constraints/` | 3 constraint modules | Placeholder data |
| `exo_compact/visualization/` | 1 plotter | Publication-quality |
| `simulations/` | 6 runnable scripts | All verified |
| `docs/` | 4 documentation files | Current |
| `visualizations/` | 10 generated figures | Latest run |
