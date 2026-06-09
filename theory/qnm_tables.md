# Quasinormal Mode Reference Tables

## 1. Kerr QNM Frequencies

Fitting formula for (l=2, m=2, n=0) mode:
```
ω = [0.4434 + 0.3236(1 - a)^0.3105] / M
τ = 10.853 * exp[0.7626(1 - a)^0.1775] * M
```

Source: Berti, Cardoso & Will (2006), PRD 73, 064030, Tables VIII-X.

### Selected (l=2, m=2, n=0) values

| a | Mω | M/τ | f (10 M_sun) [Hz] | τ [ms] |
|---|---|---|---|---|
| 0.00 | 0.3737 | 0.0890 | 241 | 0.56 |
| 0.10 | 0.3775 | 0.0900 | 243 | 0.55 |
| 0.20 | 0.3826 | 0.0914 | 247 | 0.54 |
| 0.30 | 0.3895 | 0.0932 | 251 | 0.53 |
| 0.40 | 0.3987 | 0.0957 | 257 | 0.51 |
| 0.50 | 0.4110 | 0.0990 | 265 | 0.50 |
| 0.60 | 0.4276 | 0.1037 | 276 | 0.48 |
| 0.70 | 0.4508 | 0.1107 | 291 | 0.45 |
| 0.80 | 0.4848 | 0.1218 | 313 | 0.41 |
| 0.90 | 0.5389 | 0.1405 | 348 | 0.35 |
| 0.98 | 0.6326 | 0.1733 | 408 | 0.28 |

## 2. Exotic Object QNM Shifts

[REQUIRES INPUT] QNM frequencies for non-Kerr objects.

Expected deviations (from literature survey):
- **ECOs with surface**: ω_ECO ≈ ω_Kerr + δω(R, ε_q), where δω depends on reflectivity
- **Wormholes**: QNM spectrum shifted due to modified effective potential and boundary conditions
- **Boson stars**: QNM frequencies depend on compactness; typically lower than BH equivalent

## 3. References

- Berti, E., Cardoso, V., Will, C.M. "Gravitational-wave spectroscopy of massive black holes with the space interferometer LISA." PRD 73, 064030 (2006).
- Kokkotas, K.D., Schmidt, B.G. "Quasi-normal modes of stars and black holes." Living Rev. Rel. 2, 2 (1999).
- Cardoso, V., Pani, P. "Testing the nature of dark compact objects: a status report." Living Rev. Rel. 22, 4 (2019).
