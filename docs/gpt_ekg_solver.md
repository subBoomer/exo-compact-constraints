# EKG Solver: Remaining Convergence Issue

## Current state

Two bugs fixed in previous iterations:
1. **Phi-omega rescaling sign**: Changed `omega * exp(+Phi_inf)` → `omega * exp(-Phi_inf)`. Now `omega_phys < m_s` and `alpha(0) = exp(Phi(0)) < 1` (physical redshift at center).
2. **Energy density 1/2 factor**: Added `0.5*` to rho and p_r for complex U(1) scalar convention.

The solver finds bound-state solutions with correct asymptotic normalization but **ALL solutions are on the unstable (upper) branch** — mass decreases with sigma_c instead of increasing.

---

## The ODE system (current, verified correct)

```python
def ekg_rhs_fixed(r, y, omega, m_s):
    m, Phi, sigma, u = y      # u = dsigma/dr
    N = 1.0 - 2.0*m / r

    inv_A2 = exp(-2*Phi)

    rho = 0.5 * ((omega**2 * inv_A2 + m_s**2) * sigma**2 + u**2 * N)
    p_r = 0.5 * ((omega**2 * inv_A2 - m_s**2) * sigma**2 - u**2 * N)

    dm   = 4*pi * r**2 * rho
    dPhi = (m + 4*pi * r**3 * p_r) / (r * (r - 2*m))
    Np   = -2*dm/r + 2*m/r**2

    du   = -(2/r + dPhi + 0.5*Np/N) * u
    du  -= (omega**2 * inv_A2 - m_s**2) * sigma / N

    return [dm, dPhi, u, du]
```

---

## Numerical trace for sigma_c = 0.14, omega = 0.85

```
 step     r       m       Phi    sigma        u        N     2m/r
    0  0.001  0.0000   0.0000  0.1400  1.30e-05  1.0000  0.0000
   40  0.401  0.0049   0.0042  0.1406  3.73e-03  0.9755  0.0245
   80  0.801  0.0395   0.0165  0.1425  8.52e-03  0.9016  0.0984
  120  1.201  0.1342   0.0367  0.1466  1.57e-02  0.7765  0.2235
  160  1.601  0.3255   0.0633  0.1545  2.87e-02  0.5933  0.4067
  200  2.001  0.6631   0.0820  0.1724  8.76e-02  0.3372  0.6628
  220  2.201  1.0797  -1.2566  0.2213  1.35e+01 -1.10e-02  1.0110
                                    ^^^^^^  ^^^^^^^^^
                                    sigma growing   u blows up
```

Key observations:
- sigma **never decreases** — it grows from 0.14 to 0.22
- The integration always collapses at 2M/r ≈ 1 (when N → 0)
- u (sigma') blows up at the collapse because N'/N → -∞ creates anti-damping
- The field never enters the Yukawa decay regime

---

## Mathematical root cause

The KG equation is:
```
sigma'' = -(2/r + Phi' + N'/(2N)) * sigma' - (omega^2/A^2 - m_s^2) * sigma / N
```

For bound states (omega < m_s): `omega^2/A^2 - m_s^2 < 0`, so the last term is POSITIVE:

```
sigma'' = -(damping) + |omega^2/A^2 - m_s^2| * sigma / N
```

For sigma to turn over (sigma'' < 0 at sigma' = 0), we need:
```
omega^2/A^2 - m_s^2 > 0   →   omega > A * m_s
```

At the turning point: sigma' = 0, so sigma'' = -(omega^2/A^2 - m_s^2) * sigma / N

For sigma'' < 0 (the turnover): we need `omega > A * m_s`.

Our solutions have `omega/m_s ≈ 0.84` and `alpha(0) = A(0) ≈ 0.84`. Since A increases outward (Phi grows positive from center), A > alpha(0) everywhere except at r=0. So `omega < A * m_s` holds everywhere, and sigma can NEVER turn over.

**In other words: for all our solutions, the lapse A(r) is always larger than omega/m_s, making the KG source term always positive and sigma always increase. The field never decays because the condition for turning is never met.**

---

## The puzzle

Published boson star solutions (Liebling & Palenzuela 2023, Kaup 1968) show M(sigma_c) curves with:
- Stable branch: M increases with sigma_c from 0 to M_max ≈ 0.633
- At sigma_c = 0.14: M ≈ 0.55, omega/m_s ≈ 0.89, alpha(0) ≈ 0.78

For our solver at sigma_c = 0.14:
- M ≈ 0.66, omega/m_s ≈ 0.84, alpha(0) ≈ 0.84

The difference: **our alpha(0) is too large** (0.84 vs 0.78), which means the redshift at center is too weak. With alpha(0) = 0.78, `omega > A(0) * m_s` would be `0.84 > 0.78`, easily satisfied.

## Our alpha(0) = exp(-Phi_inf)

Phi_inf is taken from the last saved integration point (before the collapse). If the collapse happens at too small r, Phi_inf is too small, alpha(0) is too large.

**Question: Could this be because our r_max is too small, so Phi hasn't converged to its asymptotic value?** Our integration typically collapses at r ≈ 2-3, where Phi is still evolving rapidly. If we could integrate further (past the mode where 2M/r approaches 1), would Phi converge to a larger asymptotic value, making alpha(0) smaller?

Or is there a more fundamental issue with how we handle the boundary conditions?

---

## Specific questions

1. Our integration uses `H = 0.01` with RK4. Could the step size be causing numerical dissipation that prevents sigma from reaching its turnover point?

2. The shooting bisection converges on the unstable branch. Should we instead shoot on a DIFFERENT observable — like the asymptotic decay rate `d/dr(ln sigma) → -k` — to find the stable branch?

3. Could the near-origin expansion for Phi_2 be wrong when rho_c has the 1/2 factor? Currently:
   ```python
   Phi_2 = 2*pi * ((omega**2 - m_s**2) * sigma_c**2 + rho_c / 3.0)
   ```
   where `rho_c = 0.5 * (omega**2 + m_s**2) * sigma_c**2`

4. Should we enforce `alpha(0) < omega/m_s` as an ADDITIONAL constraint during shooting, to ensure the field can turn over?

5. The solve_boson_star function currently does simple bisection on omega. For the stable branch, should we use a different algorithm (e.g., Newton-Raphson on `sigma(r_max)` or on the asymptotic decay rate)?

---

## What we need

A specific, implementable suggestion for:
1. What shooting criterion to use (instead of simply checking if sigma < 1e-10)
2. How to ensure the solver converges to the stable branch
3. A diagnostic check to verify we're on the correct branch

Target: reproduce the M(sigma_c) sequence with M_max ≈ 0.633 at the correct sigma_c.

---

## Files

The full solver is at: `exo_compact/metrics/boson_star_solver.py` (222 lines)
