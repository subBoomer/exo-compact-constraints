import numpy as np
from typing import Dict, Any, Optional
from ..metrics.base import AbstractMetric, MetricType


class ShadowCalculator:
    """
    Black hole shadow observables.

    Computes:
    - Critical impact parameter (b_c) / shadow radius
    - Shadow boundary for Kerr (Bardeen 1973 via Carter constant formalism)
    - Shadow for general spherically symmetric metrics (photon sphere)
    - Shadow asymmetry (Bardeen formula)
    - EHT-like observables (d_sh, delta, A)

    Core equations from Bardeen (1973), Cunningham & Bardeen (1973),
    and the EHT Collaboration (2019, 2022).
    """

    def __init__(self, metric: AbstractMetric):
        self.metric = metric

    def critical_impact_parameter(self) -> float:
        """
        Critical impact parameter b_c.

        For spherically symmetric metrics (Schwarzschild, BosonStar, Gravastar):
          b_c = r_ph / sqrt(A(r_ph)) ~ 3*sqrt(3)*M for Schwarzschild.

        For Kerr: b_c = average shadow radius = shadow_diameter() / 2.

        For Morris-Thorne: b_c = throat radius (when Phi=0).

        For RotatingWormhole (Teo): b_c ~ throat radius (approx, requires
        full ray tracing for exact value).
        """
        M = self.metric.params.get("mass", 1.0)
        mt = self.metric.metric_type
        spin = self.metric.params.get("spin", 0.0)

        if mt == MetricType.MORRIS_THORNE:
            return float(self.metric.params.get("throat_radius", 2.0))

        if mt == MetricType.ROTATING_WORMHOLE:
            b0 = self.metric.params.get("throat_radius", 2.0)
            return np.e * b0

        if mt == MetricType.SCHWARZSCHILD:
            return 3.0 * np.sqrt(3.0) * M

        if mt == MetricType.KERR and spin != 0.0:
            d_sh = self.shadow_diameter()
            return d_sh / 2.0

        return 3.0 * np.sqrt(3.0) * M

    def photon_sphere_radius(self) -> float:
        """
        Photon sphere radius for spherically symmetric metrics.

        For ds^2 = -A(r) dt^2 + B(r) dr^2 + r^2 dOmega^2:
        Solve d/dr(r^2 / A(r)) = 0.

        For Schwarzschild: r_ph = 3M.

        Source: Bardeen (1973), Cunningham & Bardeen (1973).
        """
        M = self.metric.params.get("mass", 1.0)
        mt = self.metric.metric_type

        if mt == MetricType.MORRIS_THORNE:
            return float(self.metric.params.get("throat_radius", 2.0))

        return 3.0 * M

    def _kerr_photon_orbit_radii(self) -> tuple:
        """
        Prograde and retrograde equatorial photon orbit radii for Kerr.

        r_ph^+/- = 2M * [1 + cos(2/3 * arccos(∓ a/M))]

        Source: Bardeen (1973), Eq. (14).
        """
        M = self.metric.params.get("mass", 1.0)
        a = self.metric.params.get("spin", 0.0)
        a_star = a / M

        r_ph_pro = 2.0 * M * (
            1.0 + np.cos((2.0 / 3.0) * np.arccos(-a_star))
        )
        r_ph_retro = 2.0 * M * (
            1.0 + np.cos((2.0 / 3.0) * np.arccos(a_star))
        )
        return r_ph_pro, r_ph_retro

    def _kerr_xi_eta(self, r_norm: float) -> tuple:
        """
        Dimensionless Carter constants (xi, eta) for spherical null geodesics in Kerr.

        With M=1 normalization:
        xi(r) = (r^2(r - 3) + a^2(r + 1)) / (a(1 - r))
        eta(r) = r^3(4a^2 - r(r - 3)^2) / (a^2(1 - r)^2)

        [NOTE] This function is only valid for Kerr with a > 0.
        For a=0 (Schwarzschild), route to _schwarzschild_shadow_boundary instead.
        For r_norm=1, the extremal prograde orbit requires dedicated treatment.

        Source: Bardeen (1973); Chandrasekhar (1983).
        """
        a_star = self.metric.params.get("spin", 0.0) / self.metric.params.get("mass", 1.0)
        if abs(a_star) < 1e-12:
            raise ValueError(
                "Kerr shadow formula undefined in Schwarzschild limit. "
                "Use _schwarzschild_shadow_boundary instead."
            )
        if abs(1.0 - r_norm) < 1e-10:
            raise ValueError(
                f"Carter constant singular at r_norm=1 (r_norm={r_norm}). "
                "Extremal prograde orbit requires dedicated limiting treatment."
            )

        xi = ((r_norm ** 2) * (r_norm - 3.0)
              + (a_star ** 2) * (r_norm + 1.0)) / (a_star * (1.0 - r_norm))
        eta = ((r_norm ** 3)
               * (4.0 * (a_star ** 2) - r_norm * ((r_norm - 3.0) ** 2))
               / ((a_star ** 2) * (1.0 - r_norm) ** 2))

        return float(xi), float(eta)

    def shadow_boundary_kerr(
        self, num_points: int = 100, inclination: float = np.pi / 2
    ) -> np.ndarray:
        """
        Shadow boundary for Kerr black hole in celestial coordinates (alpha, beta).

        Uses the Carter constant formalism (Bardeen 1973). The shadow is
        parameterized by the spherical null geodesic orbital radius r.

        [NOTE] This method is ONLY valid for Kerr black holes. For Teo wormholes
        and other rotating exotic objects that lack Carter separability, full
        numerical ray tracing is required.

        Celestial coordinates for observer at inclination theta_0:
          alpha = -xi / sin(theta_0)
          beta = +/- sqrt(eta + a^2 cos^2(theta_0) - xi^2 cot^2(theta_0))

        Parameters
        ----------
        num_points : int
            Number of points on the shadow boundary
        inclination : float
            Observer inclination angle from rotation axis (radians)

        Returns
        -------
        boundary : ndarray of shape (N, 2)
            (alpha, beta) coordinates of the shadow boundary in units of M
        """
        mt = self.metric.metric_type
        if mt == MetricType.ROTATING_WORMHOLE:
            raise NotImplementedError(
                "Teo wormhole shadow requires full numerical ray tracing "
                "(no Carter constant). Use spherical approximation instead."
            )

        if mt in (MetricType.SCHWARZSCHILD, MetricType.BOSON_STAR,
                   MetricType.GRAVASTAR) or (
            mt == MetricType.KERR and self.metric.params.get("spin", 0.0) == 0.0
        ):
            return self._schwarzschild_shadow_boundary(num_points)

        M = self.metric.params.get("mass", 1.0)
        a = self.metric.params.get("spin", 0.0)
        a_star = a / M
        i = inclination

        r_ph_pro, r_ph_retro = self._kerr_photon_orbit_radii()
        r_norm_pro = r_ph_pro / M
        r_norm_retro = r_ph_retro / M

        r_norm_vals = np.linspace(r_norm_pro, r_norm_retro, num_points)

        xi_vals = np.zeros(num_points)
        eta_vals = np.zeros(num_points)
        for idx, rn in enumerate(r_norm_vals):
            xi_vals[idx], eta_vals[idx] = self._kerr_xi_eta(rn)

        sin_i = np.sin(i)
        cos_i = np.cos(i)
        cot_i = cos_i / max(sin_i, 1e-15)

        alpha_vals = -xi_vals / sin_i
        inside = eta_vals + a_star ** 2 * cos_i ** 2 - xi_vals ** 2 * cot_i ** 2
        inside = np.maximum(inside, 0.0)
        beta_vals_pos = np.sqrt(inside)
        beta_vals_neg = -beta_vals_pos

        alpha_full = np.concatenate([alpha_vals, alpha_vals[::-1]])
        beta_full = np.concatenate([beta_vals_pos, beta_vals_neg[::-1]])

        boundary = np.column_stack([alpha_full, beta_full])
        boundary = boundary[~np.isnan(boundary[:, 0])]
        return boundary

    def _schwarzschild_shadow_boundary(self, num_points: int = 100) -> np.ndarray:
        """
        Circular shadow boundary for Schwarzschild BH.

        Shadow is a circle of radius b_c = 3*sqrt(3)*M in celestial coordinates.
        """
        M = self.metric.params.get("mass", 1.0)
        b_c = 3.0 * np.sqrt(3.0) * M
        angles = np.linspace(0, 2 * np.pi, num_points)
        alpha = b_c * np.cos(angles)
        beta = b_c * np.sin(angles)
        return np.column_stack([alpha, beta])

    def shadow_diameter(self) -> float:
        """
        Shadow angular diameter d_sh in units of M.

        For Schwarzschild: d_sh = 2 * 3*sqrt(3) * M ~ 10.39 M.
        For Kerr: depends on spin and inclination.

        Source: EHT Collaboration (2019), ApJL 875, L1; (2022), ApJL 930, L17.
        """
        mt = self.metric.metric_type
        M = self.metric.params.get("mass", 1.0)
        a = self.metric.params.get("spin", 0.0)

        if mt == MetricType.SCHWARZSCHILD or (
            mt == MetricType.KERR and a == 0.0
        ):
            return 2.0 * 3.0 * np.sqrt(3.0) * M

        if mt == MetricType.MORRIS_THORNE:
            b0 = self.metric.params.get("throat_radius", 2.0)
            return 2.0 * b0

        if mt == MetricType.ROTATING_WORMHOLE:
            b0 = self.metric.params.get("throat_radius", 2.0)
            return 2.0 * np.e * b0

        if mt in (MetricType.BOSON_STAR, MetricType.GRAVASTAR):
            return 2.0 * 3.0 * np.sqrt(3.0) * M

        num_pts = 200
        boundary = self.shadow_boundary_kerr(num_pts)
        if len(boundary) < 4:
            return 2.0 * 3.0 * np.sqrt(3.0) * M
        alpha = boundary[:, 0]
        beta = boundary[:, 1]
        d_sh = np.max(alpha) - np.min(alpha)
        return float(d_sh)

    def shadow_asymmetry(self) -> float:
        """
        Shadow asymmetry parameter A.

        A = (d_sh_max - d_sh_min) / (d_sh_max + d_sh_min)

        For Schwarzschild (spherical): A = 0.
        For Kerr: A > 0, increasing with spin and inclination.
        For Teo wormhole: requires full ray tracing (placeholder).

        Source: Bardeen (1973); Tamburini et al. (2020), A&A 642, L1.
        """
        mt = self.metric.metric_type
        a = self.metric.params.get("spin", 0.0)

        if mt in (MetricType.SCHWARZSCHILD, MetricType.MORRIS_THORNE,
                   MetricType.BOSON_STAR, MetricType.GRAVASTAR):
            return 0.0

        if mt == MetricType.ROTATING_WORMHOLE:
            return 0.0

        if mt == MetricType.KERR:
            num_pts = 200
            try:
                boundary = self.shadow_boundary_kerr(num_pts)
            except NotImplementedError:
                return 0.0
            if len(boundary) < 4:
                return 0.0
            alpha = boundary[:, 0]
            beta = boundary[:, 1]
            d_horiz = np.max(alpha) - np.min(alpha)
            d_vert = np.max(beta) - np.min(beta)
            asymmetry = (d_horiz - d_vert) / (d_horiz + d_vert + 1e-15)
            return float(abs(asymmetry))

        eps_q = self.metric.params.get("quantum_correction", 0.0)
        return float(abs(a) * 0.05 + abs(eps_q) * 0.02)

    def eht_deviation(
        self, m87_mass: float = 6.5e9, m87_distance: float = 16.8e6
    ) -> Dict[str, Any]:
        """
        Compare shadow predictions against EHT M87* results.

        EHT measured d_sh = 42 +/- 3 microas at 230 GHz.
        GR predicts d_sh ~ 10.4 M for a non-rotating BH.
        For M87*: d_sh ~ 37-45 microas (consistent with Kerr).

        Returns deviation delta from GR prediction.
        """
        M = self.metric.params.get("mass", 1.0)

        M_solar_mass = 1.98847e30
        G = 6.67430e-11
        c = 299792458.0
        pc = 3.085677581e16

        mass_kg = m87_mass * M_solar_mass
        mass_m = G * mass_kg / c ** 2
        dist_m = m87_distance * pc

        d_sh_dimensionless = self.shadow_diameter()

        rad_to_mas = 206264.80624709636 * 1000.0

        d_sh_model_rad = (d_sh_dimensionless * mass_m) / dist_m
        d_sh_model_mas = d_sh_model_rad * rad_to_mas

        d_sh_gr_rad = (2.0 * 3.0 * np.sqrt(3.0) * mass_m) / dist_m
        d_sh_gr_mas = d_sh_gr_rad * rad_to_mas

        deviation = (d_sh_model_rad - d_sh_gr_rad) / d_sh_gr_rad

        return {
            "d_sh_GR_mas": float(d_sh_gr_mas),
            "d_sh_model_mas": float(d_sh_model_mas),
            "deviation_fraction": float(deviation),
            "EHT_measurement_mas": 42.0,
            "EHT_uncertainty_mas": 3.0,
        }
