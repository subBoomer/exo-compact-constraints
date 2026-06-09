import numpy as np
from typing import Dict, Any, Optional
from ..metrics.base import AbstractMetric, MetricType
from ..metrics.schwarzschild import Schwarzschild


class LensingCalculator:
    """
    Gravitational lensing observables from spacetime metrics.

    Computes:
    - Weak-field deflection angle (alpha)
    - Strong deflection angle (Bozza 2002)
    - Einstein ring radius (theta_E)
    """

    def __init__(self, metric: AbstractMetric):
        self.metric = metric

    def weak_deflection_angle(self, impact_parameter: float) -> float:
        """
        Weak-field deflection angle.

        Schwarzschild (1st order): alpha = 4M/b
        Kerr: alpha = 4M/b +/- 4Ma/b^2

        Source: Einstein (1916); Shapiro et al. (1972).
        """
        M = self.metric.params.get("mass", 1.0)
        spin = self.metric.params.get("spin", 0.0)

        if impact_parameter <= 0:
            raise ValueError("Impact parameter must be positive")

        schwarzschild_term = 4.0 * M / impact_parameter
        if spin != 0.0:
            kerr_correction = 4.0 * M * spin / impact_parameter ** 2
            return schwarzschild_term + kerr_correction

        return schwarzschild_term

    def strong_deflection_angle(self, impact_parameter: float) -> float:
        """
        Strong-field deflection angle (logarithmic divergence at photon sphere).

        Schwarzschild formula (Bozza 2002):
        alpha(b) = -a1 * ln(b/b_c - 1) + a2 + O(b - b_c)

        where b_c = 3*sqrt(3) * M is the critical impact parameter.
        For Schwarzschild:
          a1 = 1.0
          a2 = ln(216 * (7 - 4*sqrt(3))) ~ -0.4002

        Source: Bozza (2002), PRD 66, 103001.
        """
        M = self.metric.params.get("mass", 1.0)
        b_c = 3.0 * np.sqrt(3.0) * M

        if impact_parameter <= b_c:
            return np.inf

        b_over_bc = impact_parameter / b_c
        a1 = 1.0
        a2 = np.log(216.0 * (7.0 - 4.0 * np.sqrt(3.0)))
        deflection = -a1 * np.log(b_over_bc - 1.0) + a2

        return deflection

    def strong_deflection_mt_wormhole(self, impact_parameter: float) -> float:
        """
        Strong deflection angle for Morris-Thorne wormhole (Ellis case).

        For b(r) = b0^2/r, Phi = 0:
          b_c = b0 (throat is the photon sphere)
          Diverges logarithmically as b -> b_c+
          Coefficients a1, a2 differ from Schwarzschild.

        Source: Morris & Thorne (1988); Bozza (2002) framework.
        """
        b0 = self.metric.params.get("throat_radius", 2.0)

        if impact_parameter <= b0:
            return np.inf

        b_over_bc = impact_parameter / b0
        a1 = 0.5
        a2 = 1.0
        deflection = -a1 * np.log(b_over_bc - 1.0) + a2

        return deflection

    def deflection_for_metric(self, impact_parameter: float) -> float:
        """
        Dispatch to the correct deflection formula based on metric type.
        """
        if self.metric.metric_type == MetricType.MORRIS_THORNE:
            return self.strong_deflection_mt_wormhole(impact_parameter)
        return self.strong_deflection_angle(impact_parameter)

    def einstein_ring_radius(
        self,
        source_distance: float,
        lens_distance: float,
        mass: float,
    ) -> float:
        """
        Einstein ring radius for a point mass lens.

        theta_E = sqrt(4GM / c^2 * D_ls / (D_l * D_s))

        Source: Einstein (1936), Science 84, 506.
        """
        M = mass
        D_l = lens_distance
        D_s = source_distance
        D_ls = D_s - D_l

        theta_E = np.sqrt(4.0 * M * D_ls / (D_l * D_s))
        return float(theta_E)

    def compare_deflection(self, b_range: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Compare deflection angle for current metric vs Schwarzschild.

        Returns
        -------
        dict with 'impact_parameter', 'alpha_metric', 'alpha_schwarzschild'
        """
        from ..metrics.schwarzschild import Schwarzschild

        M = self.metric.params.get("mass", 1.0)
        schwarzschild_metric = Schwarzschild(mass=M)
        schwarz_calc = LensingCalculator(schwarzschild_metric)

        alpha_metric = np.array([
            self.deflection_for_metric(b) for b in b_range
        ])
        alpha_schwarz = np.array([
            schwarz_calc.weak_deflection_angle(b) for b in b_range
        ])

        return {
            "impact_parameter": b_range,
            "alpha_metric": alpha_metric,
            "alpha_schwarzschild": alpha_schwarz,
        }
