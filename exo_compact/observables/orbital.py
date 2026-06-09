import numpy as np
from typing import Dict, Any
from ..metrics.base import AbstractMetric


class OrbitalCalculator:
    """
    Orbital dynamics observables from spacetime metrics.

    Computes:
    - Periastron precession (Δφ per orbit)
    - Lense-Thirring precession (frame dragging)
    - ISCO radius and frequency
    - Orbital frequencies (Ω_r, Ω_θ, Ω_φ)
    - Epicyclic frequencies (for QPO models)
    """

    def __init__(self, metric: AbstractMetric):
        self.metric = metric

    def periastron_precession(self, semi_major: float, eccentricity: float = 0.0) -> float:
        """
        Periastron precession per orbit.

        For Schwarzschild (1st post-Newtonian):
        Δφ = 6πM / (a(1 - e²))

        For Kerr: additional Lense-Thirring term.

        Source: Einstein (1915), PR 47, 898 (Mercury perihelion);
        Weinberg (1972), Gravitation and Cosmology.
        """
        M = self.metric.params.get("mass", 1.0)
        spin = self.metric.params.get("spin", 0.0)

        p = semi_major * (1.0 - eccentricity ** 2)
        r_isco = self.isco_radius()
        if p <= r_isco:
            raise ValueError(f"Orbit with p={p:.2f} is below ISCO ({r_isco:.2f})")

        delta_phi_schwarz = 6.0 * np.pi * M / p

        if spin != 0.0:
            delta_phi_lt = 4.0 * np.pi * spin / (p ** 1.5)
            return delta_phi_schwarz + delta_phi_lt

        return delta_phi_schwarz

    def lense_thirring_precession(self, radius: float, inclination: float = np.pi / 2) -> float:
        """
        Lense-Thirring precession (frame dragging) rate.

        Ω_LT = 2J / r³ (per orbit)

        Source: Lense & Thirring (1918), Phys. Z. 19, 156.
        """
        M = self.metric.params.get("mass", 1.0)
        spin = self.metric.params.get("spin", 0.0)

        J = spin * M
        return 2.0 * J / radius ** 3

    def isco_radius(self) -> float:
        """
        Innermost Stable Circular Orbit radius.

        For Schwarzschild: r_ISCO = 6M
        For Kerr (prograde): r_ISCO = M * (3 + Z₂ ∓ sqrt((3 - Z₁)(3 + Z₁ + 2Z₂)))
        where Z₁ = 1 + (1 - a²/M²)^{1/3} * ((1 + a/M)^{1/3} + (1 - a/M)^{1/3})
        Source: Bardeen, Press & Teukolsky (1972), ApJ 178, 347.
        """
        M = self.metric.params.get("mass", 1.0)
        spin = self.metric.params.get("spin", 0.0)

        if spin == 0.0:
            return 6.0 * M

        a_star = spin / M
        z1 = 1.0 + (1.0 - a_star ** 2) ** (1.0 / 3.0) * (
            (1.0 + a_star) ** (1.0 / 3.0) + (1.0 - a_star) ** (1.0 / 3.0)
        )
        z2 = np.sqrt(3.0 * a_star ** 2 + z1 ** 2)

        r_isco_pro = M * (3.0 + z2 - np.sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2)))
        return float(r_isco_pro)

    def isco_frequency(self) -> float:
        """
        Orbital frequency at ISCO in Hz.

        Ω_ISCO = M^{1/2} / (r_ISCO^{3/2} + a M^{1/2})
        """
        M = self.metric.params.get("mass", 1.0)
        spin = self.metric.params.get("spin", 0.0)

        M_sun_kg = 1.98847e30
        G = 6.67430e-11
        c = 299792458.0
        mass_kg = M * M_sun_kg
        M_in_meters = G * mass_kg / c ** 2

        r_isco = self.isco_radius()

        omega = np.sqrt(M) / (r_isco ** 1.5 + spin * np.sqrt(M))
        freq_hz = omega / (M_in_meters / c) / (2.0 * np.pi)

        return float(freq_hz)

    def epicyclic_frequencies(self, radius: float) -> Dict[str, float]:
        """
        Radial and vertical epicyclic frequencies for nearly circular orbits.

        For Schwarzschild:
        Ω_r = Ω_φ * sqrt(1 - 6M/r)
        Ω_θ = Ω_φ

        Source: Kato, Fukue & Mineshige (2008), Black-Hole Accretion Disks.
        """
        M = self.metric.params.get("mass", 1.0)

        omega_phi = np.sqrt(M / radius ** 3)
        omega_r = omega_phi * np.sqrt(1.0 - 6.0 * M / radius)
        omega_theta = omega_phi

        return {
            "omega_phi": float(omega_phi),
            "omega_r": float(omega_r),
            "omega_theta": float(omega_theta),
        }

    def compare_precession(
        self, radius_range: np.ndarray
    ) -> Dict[str, np.ndarray]:
        """
        Compare orbital precession for current metric vs Schwarzschild.
        Filters out radii below ISCO.
        """
        from ..metrics.schwarzschild import Schwarzschild
        from .orbital import OrbitalCalculator as BaseCalc

        M = self.metric.params.get("mass", 1.0)
        r_min = 6.0 * M * 1.05
        valid = radius_range > r_min
        if not np.all(valid):
            import warnings
            warnings.warn(f"Filtering {np.sum(~valid)} radii below ~{r_min:.1f}M")
        r_valid = radius_range[valid]

        schwarz_metric = Schwarzschild(mass=M)
        schwarz_calc = BaseCalc(schwarz_metric)

        delta_phi_metric = np.array([
            self.periastron_precession(r) for r in r_valid
        ])
        delta_phi_schwarz = np.array([
            schwarz_calc.periastron_precession(r) for r in r_valid
        ])

        return {
            "radius": r_valid,
            "delta_phi_metric": delta_phi_metric,
            "delta_phi_schwarzschild": delta_phi_schwarz,
        }
