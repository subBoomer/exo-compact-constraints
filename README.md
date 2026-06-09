# Exotic Compact Object Constraints

<!-- TODO: Replace with Zenodo DOI after archiving -->
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) -->

**Simulation-driven framework for distinguishing black holes from wormholes and exotic compact objects using gravitational lensing, shadow imaging, gravitational waves, and orbital dynamics.**

36 source files · 8 modules · 11 publication-ready figures · Hamiltonian ray tracing · Real EHT/LIGO/pulsar constraint data

**Author:** Timur Kiselchuk · Independent Researcher · 2026

---

## Quick Start

```bash
# 1. Install
pip install -e .

# 2. Run the full comparison
python simulations/run_comparison.py

# 3. Explore results
open visualizations/comparison_shadow_overlay.png
open visualizations/comparison_exclusion.png
```

## Key Results

| Object | d_sh/M | χ² | Status |
|---|---|---|---|
| Schwarzschild | 10.39 | 0.00 | GR baseline |
| Kerr (a=0.7) | 10.05 | 0.03 | Consistent |
| MT Wormhole (b0=2M) | 4.00 | **75.80** | **Excluded** |
| Teo Wormhole (b0=2M) | 5.44 | **22.74** | **Excluded** |
| Boson Star (C=0.12) | 10.39 | 0.00 | Consistent |
| Gravastar (r1=3M) | 10.39 | 0.00 | Consistent |
| ECO (R=0.5) | 10.05 | 0.19 | Consistent |

χ² computed against EHT shadow (d_sh = 42 ± 3 μas for M87*) + LIGO ringdown constraints (δf/f < 0.05, δτ/τ < 0.2). Excluded: χ² > 6.

### What this means

Wormholes with throat radius b0 ≤ 2M produce shadows 40-50% smaller than the EHT measurement for M87* and are excluded at > 95% confidence. The Teo wormhole shadow (b_c = e · b0) is confirmed via Hamiltonian ray tracing with 0.2% accuracy. Rotating ECOs with reflectivity R ≤ 0.5 remain consistent with current observations.

## Architecture

**Metric-first**: All observables derive from spacetime metrics through a unified Hamiltonian geodesic engine.

```
Metric ──→ Hamiltonian Geodesic ──→ Observable ──→ Comparison ──→ Exclusion
```

## Modules

| Module | Files | Purpose |
|---|---|---|
| `metrics/` | 10 | 7 metric classes (Schwarzschild, Kerr, MT Wormhole, Teo Wormhole, Boson Star, Gravastar, ECO) + experimental EKG solver |
| `observables/` | 5 | Lensing, shadow, GW ringdown, orbital dynamics |
| `geodesics/` | 4 | Null/timelike integration, coordinate-based ray tracer |
| `ray_tracing/` | 5 | Production-grade Hamiltonian geodesic engine (DOP853, metric caching, termination oracle) |
| `constraints/` | 4 | EHT (8 entries), LIGO (7), pulsar (5) — real numerical data |
| `models/` | 4 | BH, wormhole, ECO model factories |
| `visualization/` | 2 | Publication-quality comparison plots |

## Objects (7 metrics)

| Object | Class | Parameters | Shadow Method |
|---|---|---|---|
| Schwarzschild BH | `Schwarzschild` | mass M | Analytic: 3√3·M |
| Kerr BH | `Kerr` | mass M, spin a | Analytic: Carter constant (Bardeen 1973) |
| Morris-Thorne Wormhole | `MorrisThorne` | mass M, throat b0 | Analytic: b_c = b0 |
| Teo Rotating Wormhole | `RotatingWormhole` | mass M, throat b0, spin a | Hamiltonian ray tracing: b_c = e·b0 |
| Boson Star | `BosonStar` | mass M, compactness C | Approximate (tanh profile) |
| Gravastar | `Gravastar` | mass M, interior radius r1 | Schwarzschild exterior + de Sitter interior |
| Generic ECO | `ECO` | mass M, spin a, reflectivity R | Kerr + quantum corrections |

## Ray Tracing

The Hamiltonian formulation (`H = ½ g^μν p_μ p_ν`) eliminates coordinate singularities near photon spheres, enabling stable geodesic integration for all spacetimes.

| Metric | Method | Accuracy |
|---|---|---|
| Schwarzschild | Analytic + Hamiltonian | < 1% (5.24M vs 5.20M) |
| Kerr | Analytic (Carter constant) | Exact |
| Teo Wormhole | Hamiltonian (DOP853) | < 0.2% (5.45M vs 5.44M analytic) |

## Physics Validations

10 physics bugs identified and fixed during development:
1. Phi-omega rescaling sign in EKG solver
2. Energy density 1/2 factor for complex U(1) scalar
3. Gauge omega > m_s for bound boson stars
4. Collapse stopping condition → compactified coordinates
5. Christoffel RK4 divergence → Hamiltonian formulation
6. Teo shadow formula: b_c = e·b0 (not b0)
7. Echo spacing factor-of-2 (Schwarzschild vs Kerr formulas)
8. Gravastar interior constraint (Mazur-Mottola model)
9. String-based metric dispatch → MetricType enum
10. Kerr shadow divide-by-zero guards

## Simulations

```bash
# Full cross-object comparison (LaTeX table + 4 figures)
python simulations/run_comparison.py

# Wormhole throat-radius exclusion map (b0/M vs spin)
python simulations/scan_throat_radius.py

# Boson star mass-frequency (M-σc) sequence
python simulations/scan_boson_star.py

# Publication-ready unified exclusion diagram (3 panels)
python simulations/generate_paper_figure.py

# Individual observables
python simulations/lensing/ray_trace.py --metric kerr --mass 10 --spin 0.7
python simulations/shadow_imaging/shadow_generator.py --metric kerr --spin 0.94
python simulations/gravitational_waves/waveform_comparison.py --metric kerr --mass 68 --spin 0.7
python simulations/orbital_dynamics/orbit_integrator.py --metric kerr --mass 10 --spin 0.7
```

## Generated Figures (11)

| File | Description |
|---|---|
| `paper_figure_unified.png` | 3-panel exclusion diagram (shadow + χ² + ringdown) |
| `exclusion_throat_radius.png` | Wormhole exclusion contours in (b0/M, a) space |
| `boson_star_sequence.png` | Kaup mass-frequency curve |
| `comparison_shadow_overlay.png` | All 7 object shadows overlaid |
| `comparison_ringdown.png` | Ringdown waveforms + frequency bar chart |
| `comparison_deflection.png` | Deflection angle curves |
| `comparison_exclusion.png` | χ² bar chart against EHT + LIGO |
| `lensing_deflection_angle.png` | Single metric deflection |
| `shadow_boundary_comparison.png` | Single metric shadow |
| `ringdown_comparison.png` | Single metric ringdown |
| `orbital_precession_comparison.png` | Single metric precession |

## Requirements

- Python ≥ 3.10
- numpy, scipy, matplotlib, astropy, sympy

```bash
pip install -e .
```

## Documentation

- [System Overview](docs/system_overview.md) — design, taxonomy, pipeline, missing physics
- [Observable Signature Matrix](docs/signature_matrix.md) — which observables discriminate which objects
- [Metric Catalog](theory/metric_catalog.md) — all metrics with references
- [QNM Reference Tables](theory/qnm_tables.md) — Kerr QNM fits (Berti 2006)
- [Known Observational Constraints](constraints/known_results.md) — EHT, LIGO, pulsar data
- [Comparison Table](docs/comparison_table.tex) — LaTeX table of all 7 objects × 4 observables

## Citation

If you use this framework in your research, please cite:

```
Kiselchuk, T. (2026). Exotic Compact Object Constraints (v0.1.0).
Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX
```

## License

MIT — see [LICENSE](LICENSE) for details.
