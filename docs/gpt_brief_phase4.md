# Phase 4 Results + Next Steps — Briefing for GPT

## Current repository state

53 source files, 10 generated figures, 6 runnable simulations.
All 7 metric classes working. Comparison framework producing LaTeX tables + exclusion contours.

---

## Result 1: Wormhole throat-radius exclusion (Robust)

**Method:** Scanned b0/M ∈ [1, 10] for MT and Teo wormholes, computed χ² against EHT shadow diameter (d_sh = 10.39M ± 10%).

**Result:** Wormholes with b0/M < 4.5 are excluded at > 95% CL by EHT M87*.
EHT-consistent range: b0/M ∈ [4.5, 5.8].

**Caveat:** The Teo wormhole shadow uses the throat-radius approximation (b_c ≈ b0). Full numerical ray tracing would give more accurate boundaries.

**Question 1:** Is the b0/M ∈ [4.5, 5.8] range physically reasonable? For an Ellis wormhole with Phi=0, ADM mass M = 0, so b0/M is undefined. For a Teo wormhole with N(r) = exp(-b0/r), M = b0. A wormhole with b0 = 5.2M has throat radius 5.2 times its ADM mass — is this allowed by energy conditions?

---

## Result 2: EKG boson star solver — still not converging (Needs review)

**Current behavior:** The solver finds scalar field configurations but the field NEVER decays to zero for ω < m_s (bound states). Instead:

1. sigma(r) ≈ sigma_c * (1 + k²r²/6) near origin (correct — matches flat-space spherical Bessel)
2. sigma grows slightly, then stays flat
3. Mass accumulates as M ≈ ∫4πr²ρ dr
4. When 2M/r → 1, N = 1 - 2M/r → 0, integration collapses
5. This happens at r ≈ 2M for ALL sigma_c and ALL ω values tried

**The ODEs we use (verified against Liebling & Palenzuela 2023):**

```python
rho = (omega**2 * exp(-2*Phi) + m_s**2) * sigma**2 + N * sigma_prime**2
p_r = (omega**2 * exp(-2*Phi) - m_s**2) * sigma**2 - N * sigma_prime**2

dm/dr = 4*pi * r**2 * rho
dPhi/dr = (m + 4*pi * r**3 * p_r) / (r * (r - 2*m))

N = 1 - 2*m/r
N_prime = -2*dm/dr / r + 2*m / r**2

sigma'' = -(2/r + dPhi/dr + N_prime/(2*N)) * sigma'
          - (omega**2 * exp(-2*Phi) - m_s**2) * sigma / N
```

The KG coefficient N'/(2N) has been verified by first-principles derivation (which gave N'/(2N), not N'/N).

**Key observation:** For ALL solutions with ω < m_s, Phi grows positive in the interior. This makes A² = exp(2*Phi) > 1, which means ω²/A² - m² is always NEGATIVE. In the KG equation:

sigma'' = -(damping) - (ω²/A² - m²)σ/N
         = -(damping) - (negative)σ/N
         = -(damping) + |positive|σ/N

The last term is POSITIVE — it DRIVES sigma upward, not downward. The field never decays.

For ω > m_s (unbound, oscillatory solutions), the term becomes negative and sigma DOES decay. But ω > m_s gives oscillatory solutions at infinity, not bound states.

**Question 2:** Why does Phi grow positive in the interior? For a normal star, dΦ/dr ≈ m/r² > 0, so Φ increases from the center outward. But this creates A² > 1, which kills the restoring force. Is this a sign error somewhere, or is this the correct physical behavior for boson stars?

**Question 3:** Should the interior metric use a DIFFERENT sign convention? Some references write the metric as ds² = -N dt² + dr²/N + r² dΩ² with N = 1 - 2m(r)/r, and define Φ = 0.5*ln(N). In that case, A² = N < 1, and ω²/A² - m² could be positive. Is this the correct interior metric for boson stars?

---

## Proposed Next Steps

### Step A: Teo wormhole full ray tracing
Currently, the Teo wormhole shadow uses the throat-radius approximation (b_c ≈ b0). Full numerical null geodesic integration would give accurate shadows.

**We already have the geodesic integrator** (`exo_compact/geodesics/null_geodesics.py`). We need to:
1. Set up initial conditions for null geodesics on the image plane
2. Integrate geodesics in the Teo metric
3. Determine which geodesics escape vs. cross the throat
4. Map the shadow boundary in (α, β) coordinates

**Question 4:** What is the correct initial condition setup for ray tracing in the Teo metric? For Kerr, we use the Carter constant to parameterize the shadow boundary. For Teo (no Carter constant), we need full image-plane scanning. What's the standard approach for the initial position and direction of null geodesics?

### Step B: Fix EKG solver
The solver needs the remaining bug identified. Based on the observation that Phi > 0 → A² > 1 → no field decay for ω < m_s:

**Question 5:** Could the issue be that we're using the wrong interior metric ansatz? Some boson star codes use:
ds² = -α² dt² + (1 - 2m/r)^{-1} dr² + r² dΩ²
with α(r) determined by dα/dr = α (m + 4πr³p) / (r(r - 2m)).

If α < 1 in the interior (which would require dα/dr < 0 initially), then A² = α² < 1 and the KG restoring force would be negative (correct sign). Is there a known property where α < 1 for boson stars?

### Step C: Constraint data ingestion
Replace placeholder constraint entries with real numerical data tables.

**Question 6:** What are the BEST current numerical values for:
1. EHT M87* shadow diameter: 42 ± 3 μas (symmetric Gaussian OK?)
2. EHT Sgr A* shadow diameter: 51.8 ± 2.3 μas
3. LIGO ringdown constraint: δω/ω < ? (I've seen 0.1 and 0.2)
4. LIGO tidal deformability: Λ < ? for exotic objects

### Step D: Full degeneracy map
Compare all 7 object classes across parameter space to find regions where different objects produce identical observables.

**Question 7:** Given our current results (two wormhole classes excluded at b0=2M, Kerr-like objects consistent), what's the single most impactful next calculation for a first paper figure?

---

## Priority for your answers

**HIGH:** Questions 2, 3, 5 (EKG solver sign issue)
**MEDIUM:** Questions 1, 4 (Teo ray tracing)
**LOW:** Questions 6, 7 (Data + next steps)
