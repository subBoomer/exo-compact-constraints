# Final Project Audit — Briefing for GPT

## Project Status

**Repository**: 36 source files, 8 modules, 9 simulations, 11 figures
**Framework**: 7 metric classes, Hamiltonian ray tracing, constraint ingestion, comparison pipeline

### Key Results

| Object | d_sh/M | χ² | Status |
|---|---|---|---|
| Schwarzschild | 10.39 | 0.00 | GR baseline |
| Kerr (a=0.7) | 10.05 | 0.14 | Consistent |
| MT Wormhole (b0=2M) | 4.00 | 75.80 | Excluded |
| Teo Wormhole (b0=2M) | 5.44 | 22.74 | Excluded |
| Boson Star (C=0.12) | 10.39* | 0.00* | *nonphysical shadow |
| Gravastar (r1=3M) | 10.39 | 0.00 | Consistent |
| ECO (R=0.5) | 10.05 | 0.30 | Consistent |

### Physics Bugs Fixed (with your help)
1. Phi-omega rescaling sign
2. Energy density 1/2 factor
3. Gauge omega > m_s
4. Collapse stopping condition (compactified coordinates)
5. Christoffel RK4 divergence (Hamiltonian formulation)
6. Teo shadow formula (b_c = e × b0, verified by ray tracing)

---

## Audit Findings — Issues Requiring Your Input

A comprehensive codebase audit found 35 issues (2 critical, 7 high, 12 medium, 14 low). The most important ones for your review:

### Critical: Kerr shadow divide-by-zero (shadow.py:107-113)

`_kerr_xi_eta` divides by `a_star` and `(1.0 - r_norm)` with no guards:
```python
xi = (r_norm²(r_norm-3) + a_star²(r_norm+1)) / (a_star * (1-r_norm))
eta = r_norm³(4a_star² - r_norm(r_norm-3)²) / (a_star² * (1-r_norm)²)
```
When `a=0` (Schwarzschild limit) or `r_norm=1` (extremal prograde orbit), these produce NaN/Inf. The callers guard against `a=0` by routing to Schwarzschild, but the function itself is unsafe for standalone use.

**Question 1**: Should we add guards (`if abs(a_star) < 1e-15: raise ValueError`) or compute the analytic limit (a→0, r_norm→1) directly?

### High: Echo spacing factor-of-2 discontinuity (gravitational_waves.py:187-190)

The echo spacing formula has two branches:
- Schwarzschild (a=0): `Δt = 4M × ln(1/ε)`
- Kerr (a small): `Δt = 4M × (1 + 1/√(1-a²)) × ln(1/ε) → 8M × ln(1/ε)` as a→0

The Kerr formula gives **twice** the Schwarzschild value in the a→0 limit. This is a factor-of-2 inconsistency — either the Schwarzschild formula should be `8M`, or the Kerr prefactor should be `2` instead of `4`.

**Question 2**: Which is physically correct? The Mark et al. (2017) cavity model — does the echo spacing have a factor-of-2 difference between Schwarzschild and the a→0 Kerr limit?

### High: Gravastar interior constraint is backwards (gravastar.py:48-49)

Current: `if r1 <= 2*M: raise MetricError("Interior radius must be > 2M")`

Per Mazur & Mottola (2001), the de Sitter interior should be INSIDE 2M (the shell straddles the would-be horizon at r=2M). The constraint should be `r1 < 2M`. But changing this would break existing code that creates Gravastar instances with r1 > 2M.

**Question 3**: Should we relax to `r1 < 2M` (Mazur-Mottola model) or keep r1 > 2M (safe approximation)? Or add a flag?

### High: String-based metric dispatch is fragile (cross-cutting)

Throughout the observables layer, metric types are identified by `self.metric.name` string comparisons:
```python
if name == "MorrisThorne": ...
if name == "Kerr" and a == 0.0: ...
```
Any rename silently breaks dispatching. There are 20+ such comparisons across 5 files.

**Question 4**: Best practice for replacing string-based dispatch? Options:
- (a) Add an `isinstance`-like method to AbstractMetric (e.g., `metric.is_type("Kerr")`)
- (b) Use an enum with a `metric_type` property
- (c) Keep strings but centralize in a registry module

### Medium: Exclusion chi² uses GR shadow for all metrics (comparison.py:183)

The chi² against EHT uses `d_sh_GR = 10.3923` (Schwarzschild value) for ALL objects. For Kerr with a=0.7, the expected shadow is ~9.8M at 90° inclination, so the comparison is biased against spinning GR objects.

**Question 5**: Should the chi² reference be spin-dependent for Kerr (using the analytic Bardeen formula for the expected shadow at each spin)?

### Medium: Periastron precession uses hardcoded 6M ISCO (orbital.py:37)

The check `if p <= 6.0 * M` uses the Schwarzschild ISCO for all metrics, incorrectly rejecting valid Kerr orbits between `M` (prograde extremal) and `6M`. Should use `self.isco_radius()`.

**Question 6**: Simple fix — just replace `6.0` with `self.isco_radius()/M`. Any side effects?

### Medium: Null initial condition fallback is silent (ray_tracing/engine.py:94-97)

When `p_r` can't be solved from the null condition, the code silently falls back to `p_r = -1.0`, creating a non-null initial momentum. This produces incorrect geodesics without warning.

**Question 7**: Should this raise a `ValueError` instead of silently continuing?

---

## Minor Issues (for reference)

| Issue | File | Description |
|---|---|---|
| ECO shadow falls through to Schwarzschild default | shadow.py:52-56 | Missing ECO branch in critical_impact_parameter |
| Gravastar _Lambda returns length not 1/length² | gravastar.py:66-71 | Dimensional error in de Sitter matching |
| Echo_spacing no a_star clamp near extremal | gravitational_waves.py:190 | Diverges at a=1 |
| Boundary extraction misses diagonal transitions | ray_tracing/shadow.py:67-76 | Only checks horizontal/vertical |
| Two separate metric hierarchies | metrics/base.py vs ray_tracing/base.py | AbstractMetric and MetricBase are parallel but not unified |
| Duplicate physics constants (G, c, Msun) | shadow.py, gravitational_waves.py, orbital.py | Should be in a shared constants module |

---

## Priority for Your Response

**HIGH**: Questions 1-4 (critical/high bugs needing physics input)
**MEDIUM**: Questions 5-7 (design decisions)
**LOW**: Minor issues (noted for reference)
