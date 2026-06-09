# Observable Signature Matrix

Mapping: Spacetime model → predicted observable

## Legend

- ✓ = Well-known (analytic formula exists)
- ∼ = Approximate (analytic fit or perturbative)
- ? = Unknown / [REQUIRES INPUT]
- — = Not applicable

## Signature Matrix

| Observable | Schwarzschild | Kerr | Morris-Thorne WH | Rotating WH | Boson Star | Gravastar | ECO (generic) |
|---|---|---|---|---|---|---|---|
| **Lensing** | | | | | | | |
| Weak deflection α(b) | ✓ α=4M/b | ∼ α=4M/b ± 4Ma/b² | ✓ α≈4M/b + b₀²/2b² | ? | ∼ ≈Schwarzschild | ∼ ≈Schwarzschild | ∼ ≈Kerr+δε |
| Strong deflection | ✓ α=−ln(b/bc−1) | ∼ numerical | ? | ? | ∼ shifted | ∼ shifted | ∼ Kerr+δ |
| Einstein ring θ_E | ✓ θ_E²=4M D_ls/D_l D_s | ∼ small correction | ∼ ≈Schwarzschild | ? | ∼ ≈Schwarzschild | ∼ ≈Schwarzschild | ∼ ≈Kerr |
| **Shadow** | | | | | | | |
| Critical impact param b_c | ✓ 3√3 M ≈ 5.196M | ✓ depends on a, i | ? | ? | ∼ ≳ 3√3 M | ∼ > 2M | ∼ ≈Kerr + ε_q |
| Shadow diameter d_sh | ✓ 2⋅3√3 M ≈ 10.39M | ✓ depends on a, i | ? | ? | ∼ ≳ 10.39M | ∼ ≳ 10.39M | ∼ ≈Kerr |
| Shadow asymmetry A | ✓ 0 | ✓ up to ∼0.1 | ✓ 0 (spherical) | ? | ✓ 0 (spherical) | ✓ 0 (spherical) | ∼ small |
| **GW Ringdown** | | | | | | | |
| Dominant QNM ω (l=2,m=2) | ✓ ω=0.374/M | ✓ ω=0.374+0.324(1−a)^0.31 | ? | ? | ? | ? | ∼ Kerr+δω(R) |
| Damping time τ | ✓ τ=10.85M | ✓ τ≈10.85exp[0.76(1−a)^0.18] | ? | ? | ? | ? | ∼ Kerr+δτ(R) |
| Ringdown frequency f_RD | ✓ | ✓ | ? | ? | ∼ shifted | ∼ shifted | ∼ Kerr+δ |
| Tidal deformability λ | ✓ 0 | ✓ 0 | — | — | ∼ > 0 | ∼ large | ∼ > 0 |
| **Orbital** | | | | | | | |
| Periastron precession Δφ | ✓ 6πM/p | ∼ 6πM/p + LT | ? | ? | ∼ modified | ∼ ≈Schwarzschild | ∼ ≈Kerr |
| ISCO radius r_ISCO | ✓ 6M | ✓ M(3+Z₂∓...) | ? | ? | ∼ > 6M | ∼ > 6M | ∼ ≈Kerr+ε |
| Lense-Thirring Ω_LT | ✓ 0 | ✓ 2J/r³ | ✓ 0 | ? | ✓ 0 | ✓ 0 | ✓ 2J/r³+δ |
| Epicyclic frequencies | ✓ | ✓ | ? | ? | ∼ modified | ∼ ≈Schwarzschild | ∼ ≈Kerr |

## Notes

1. All Schwarzschild entries are exact within GR.
2. Kerr entries use the analytic fits of Bardeen (1973) and Berti et al. (2006).
3. Morris-Thorne entries marked "∼" assume weak-field observations far from throat.
4. Boson star and gravastar entries marked "∼" use analytic approximations to the metric.
5. ECO entries depend on the specific quantum correction model.
6. Rotating wormhole entries marked "?" require a specific metric choice.

## Key Discriminants

The most promising observables for distinguishing object classes:

| Observable | BH vs WH | BH vs BS | BH vs Gravastar | BH vs ECO |
|---|---|---|---|---|
| Shadow diameter d_sh | Strong if WH throat > 3M | Moderate | Weak (similar) | Weak |
| Shadow asymmetry A | Strong (WH spin?) | Weak (spherical) | Weak (spherical) | Weak (≈Kerr) |
| Ringdown frequency | Strong (modified cavity) | Moderate (shifted) | Strong (echoes) | Strong (echoes) |
| Ringdown damping τ | Strong | Moderate | Strong | Strong |
| Tidal deformability | — | Strong (> 0) | Strong (Love #) | Moderate |
| Orbital precession | Moderate | Weak | Weak | Weak |
| ISCO radius | Strong (throat cutoff) | Moderate (>6M) | Weak (≈6M) | Weak |

## References

- Bardeen (1973), "Kerr metric black holes", in *Black Holes* (DeWitt & DeWitt)
- Berti, Cardoso & Will (2006), PRD 73, 064030
- Bozza (2002), PRD 66, 103001
- Cardoso & Pani (2019), Living Rev. Rel. 22, 4
- Darwin (1959), MNRAS 120, 240
- Mazur & Mottola (2001), PRD 64, 104022
- Morris & Thorne (1988), Am. J. Phys. 56, 395
- Teo (1998), PRD 58, 024014
- EHT Collaboration (2019, 2022), ApJL
