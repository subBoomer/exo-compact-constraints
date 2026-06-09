# Phase 3 Implementation Results — Briefing for Physics Review

We implemented your previous answers (Q1–Q8) in the exo-compact-constraints framework.
Below are the specific implementations and the numerical issues encountered. Please verify correctness.

---

## 1. Boson Star EKG Solver — Needs Review

**File:** `exo_compact/metrics/boson_star_solver.py`

We implemented the ODE system from your Q3 answer:

```python
# ODE system (state: m, Phi, sigma, u)
def ekg_rhs_fixed(r, y, omega, m_s):
    m, Phi, sigma, u = y
    N = 1.0 - 2.0 * m / r
    inv_A2 = exp(-2 * Phi)

    rho = (omega**2 * inv_A2 + m_s**2) * sigma**2 + u**2 * N
    pr  = (omega**2 * inv_A2 - m_s**2) * sigma**2 - u**2 * N

    dm   = 4 * pi * r**2 * rho
    dPhi = (m + 4 * pi * r**3 * pr) / (r * (r - 2*m))
    Np   = -2 * dm / r + 2 * m / r**2
    du   = -(2/r + dPhi + 0.5 * Np/N) * u - (inv_A2 * omega**2 - m_s**2) * sigma / N
```

**Shooting bisection:** binary search on omega in (0.001, 0.999*m_s).

**Problem:** For ALL omega values tested, the integration collapses when N = 1 - 2m/r → 0 at small radii (r ~ 0.5–9). This happens for every sigma_c in [0.01, 0.3] and every omega. The mass at collapse is typically 0.1–4 (code units, m_s=1), which far exceeds the Kaup limit M_max ≈ 0.633.

We checked and corrected:
1. σ'²/B² term: was u²/N, corrected to u² * N (since B² = 1/N)
2. B'/B = -N'/(2N) coefficient: correctly implemented
3. Bisection direction: corrected (collapse → omega too low → raise lower bound)

**Question:** Is there still a sign error in rho or p_r? We used:
- rho = (ω²/A² + m_s²)σ² + N(σ')²
- p_r = (ω²/A² - m_s²)σ² - N(σ')²

Should the σ'² term in rho be σ'²/B² (which equals Nσ'²) or B²σ'² (which equals σ'²/N)?

---

## 2. Boson Star Production Metric — Calibrated Tanh

**File:** `exo_compact/metrics/boson_star.py`

Since the numerical solver isn't converging, we use a calibrated tanh mass profile:

```python
m(r) = M * tanh(r/R)^p    where p = 2 + 5*C
R = M / C (radius from compactness)
```

The redshift function uses:
- Interior: Phi(r) ≈ ln(1 - 0.5*m(r)/r) (Newtonian)
- Exterior: Phi(r) ≈ 0.5*ln(1 - 2M/r) (Schwarzschild)
- Smooth interpolation via tanh²

**Questions:**
1. Is p = 2 + 5*C a reasonable heuristic for mini-boson-star profiles?
2. The redshift interpolation — should it use the TOV equation dPhi/dr = (m + 4πr³p)/(r(r-2m)) instead?
3. For compactness C = 0.12, R = M/C ≈ 8.3M. External observers see Schwarzschild exterior — is the shadow truly indistinguishable from Schwarzschild?

---

## 3. Teo Rotating Wormhole — Metric Implementation

**File:** `exo_compact/metrics/rotating_wormhole.py`, method `set_teo_metric()`

```python
def set_teo_metric(self):
    r0 = self.throat_radius
    J = self.spin * r0          # a = J/M^2, M = r0
    N(r) = exp(-r0/r)           # lapse
    A(r) = 1/sqrt(1 - b0/r)     # radial metric (b0 = r0)
    B(r,theta) = r               # angular
    C(r,theta) = r * sin(theta)  # angular
    omega(r) = 2*J/r^3           # frame dragging
```

Metric: ds² = -N² dt² + A² dr² + B² dθ² + C² (dφ - ω dt)²

**Questions:**
1. Is this the correct Teo metric form (PRD 58, 024014)?
2. We assume M = r0 (the throat radius). Is this correct for the ADM mass?
3. Can |a| = J/M² truly exceed 1, as some references claim?
4. Currently, the shadow uses the Kerr boundary formula (which is wrong — Teo lacks Carter separability). Should we remove this fallback and require full ray tracing?

---

## 4. Gravastar — Shell Implementation

**File:** `exo_compact/metrics/gravastar.py`

Interior: de Sitter g_tt = -(1 - r²/Λ²) with Λ = r₁ (interior radius).
Exterior: Schwarzschild.
Shell: tanh interpolation between de Sitter and Schwarzschild over width w.

Current constraint: r₁ > 2M (must be outside Schwarzschild radius).

**Question:** From your Q7 answer, the Mazur-Mottola construction places the shell at r₁ ≈ 2M(1-ε), r₂ ≈ 2M(1+ε) with ε ~ 10⁻¹⁰-10⁻⁵. This means r₁ < 2M. Should we relax the r₁ > 2M check to allow sub-2M interior radii?

---

## 5. Comparison Results — Check Physical Reasonableness

All 7 metrics compared at M=1 geometric units unless specified.

| Object | d_sh/M | Asymmetry | f_RD [Hz] (10 Msun) | ISCO/M | χ² | 
|---|---|---|---|---|---|
| Schwarzschild | 10.39 | 0 | 2478 | 6.00 | 0.00 |
| Kerr (a=0.7) | 10.05 | 0.017 | 2152 | 3.39 | 0.14 |
| MT Wormhole (b0=2M) | 4.00 | 0 | 2478 | 6.00 | 69.55 |
| RotatingWH (Teo, a=0.7) | 10.05 | 0.035 | 2152 | 3.39 | 0.23 |
| Boson Star (C=0.12) | 10.39 | 0 | 2478 | 6.00 | 0.00 |
| Gravastar (r1=3M) | 10.39 | 0 | 2478 | 6.00 | 0.00 |
| ECO (R=0.5, ε=0.01) | 10.05 | 0.035 | 2152 | 3.39 | 0.23 |

**Questions:**
1. The MT wormhole (b0=2M) is excluded at χ²=69.55. Is this conclusion robust? A wormhole with b0 closer to 3M would have a larger shadow. What throat radius would be consistent with EHT?
2. The RotatingWH shows d_sh=10.05 (same as Kerr) — but this uses the Kerr boundary formula, which is wrong for Teo. How much would the Teo shadow differ from Kerr's?
3. BosonStar and Gravastar show χ²=0 — is this because their shadows really are indistinguishable from Schwarzschild for these parameters?
4. Kerr a=0.7: d_sh=10.05M instead of 10.39M. Is this ~3% reduction the correct magnitude for a=0.7?

---

## Priority for your review

**HIGH: Question 1 (EKG solver) + Question 5.1 (wormhole exclusion robustness)**
**MEDIUM: Questions 2, 3, 5.3, 5.4**
**LOW: Question 4 (gravastar shell)**
