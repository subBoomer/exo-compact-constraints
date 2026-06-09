# Metric Catalog

## 1. Schwarzschild (1916)

**Coordinates**: Schwarzschild (t, r, θ, φ)

**Metric**:
```
ds² = -(1 - 2M/r) dt² + (1 - 2M/r)⁻¹ dr² + r² dθ² + r² sin²θ dφ²
```

**Parameters**: M (mass)

**Horizon**: r_h = 2M

**Photon sphere**: r_ph = 3M

**ISCO**: r_ISCO = 6M

**Source**: K. Schwarzschild, Sitzungsber. Preuss. Akad. Wiss. (1916) 189.

---

## 2. Kerr (1963)

**Coordinates**: Boyer-Lindquist (t, r, θ, φ)

**Metric**:
```
ds² = -(1 - 2Mr/Σ) dt² - (4Mar sin²θ/Σ) dt dφ + (Σ/Δ) dr²
      + Σ dθ² + (r² + a² + 2Ma²r sin²θ/Σ) sin²θ dφ²
```
where Σ = r² + a²cos²θ, Δ = r² - 2Mr + a²

**Parameters**: M (mass), a (spin)

**Horizons**: r± = M ± √(M² - a²)

**Photon sphere**: depends on a and orbit direction

**ISCO**: depends on a (Bardeen, Press & Teukolsky 1972)

**Source**: R.P. Kerr, PRL 11 (1963) 237.

---

## 3. Morris-Thorne Wormhole (1988)

**Coordinates**: Schwarzschild-like (t, r, θ, φ)

**Metric**:
```
ds² = -e^{2Φ(r)} dt² + dr²/(1 - b(r)/r) + r² dθ² + r² sin²θ dφ²
```

**Parameters**: b0 (throat radius), M (mass), Φ(r) (redshift function), b(r) (shape function)

**Throat**: r = b0, where b(b0) = b0

**Flaring-out condition**: b'(r) < b(r)/r

**Source**: M.S. Morris & K.S. Thorne, Am. J. Phys. 56 (1988) 395.

---

## 4. Rotating Wormhole (Generalized Ansatz)

**Coordinates**: Generalized stationary axisymmetric (t, r, θ, φ)

**Metric**:
```
ds² = -N² dt² + A² dr² + B² dθ² + C² (dφ - ω dt)²
```

**Parameters**: M, b0, a, pluggable N(r,θ), A(r,θ), B(r,θ), C(r,θ), ω(r,θ)

**Known special cases**:
- Teo wormhole (PRD 58, 1998, 024014)
- Kashargin & Sushkov (Grav. Cosmol. 14, 2008, 80)

**Source**: This framework; Teo (1998) for a specific rotating wormhole.

---

## 5. Boson Star (Effective Metric)

**Coordinates**: Schwarzschild-like (t, r, θ, φ)

**Metric**:
```
ds² = -e^{2Φ(r)} dt² + (1 - 2m(r)/r)⁻¹ dr² + r² dθ² + r² sin²θ dφ²
```

**Parameters**: M, compactness C = M/R

**m(r) profile**: [ASSUMPTION] M tanh³(r/R)

**Source**: Kaup (1968), Ruffini & Bonazzola (1969), Liebling & Palenzuela (2023).

---

## 6. Gravastar (Mazur & Mottola 2001)

**Three regions**:
- Interior (r < r1): de Sitter, p = -ρ, ds² ≈ -(1 - r²/Λ²) dt² + ...
- Shell (r1 < r < r2): thin shell of stiff matter
- Exterior (r > r2): Schwarzschild

**Parameters**: M, r1 (interior radius), w (shell width)

**Source**: P.O. Mazur & E. Mottola, PRD 64 (2001) 104022.

---

## 7. Generic ECO (Effective Kerr + Corrections)

**Metric**: g_uv(ECO) = g_uv(Kerr) + δg_uv

**Corrections**:
- Reflection coefficient R ∈ [0, 1] (for echoes)
- Quantum correction ε_q (horizon-scale modifications)
- Tidal deformability λ (Love number, λ_BH = 0)

**Parameters**: M, a, R, ε_q, λ

**Source**: Cardoso & Pani (2019), Living Rev. Rel. 22, 4.
