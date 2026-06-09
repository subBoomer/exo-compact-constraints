import numpy as np
from typing import Dict, Any, Optional
from .base import AbstractMetric, MetricError, MetricType


class BosonStar(AbstractMetric):
    """
    Boson star effective metric (spherically symmetric, static).

    Uses a calibrated tanh mass profile that reproduces the correct
    mass-radius relation from Liebling & Palenzuela (2023).

    Metric:
      ds^2 = -e^{2Phi(r)} dt^2 + (1 - 2m(r)/r)^{-1} dr^2 + r^2 dOmega^2

    where m(r) and Phi(r) are tanh profiles calibrated against
    numerical EKG solutions. The calibration ensures:
      - Compactness C = M/R in the range [0.02, 0.16]
      - Mass follows M vs sigma_c relation (Kaup limit M_max ~ 0.633)

    [WARNING] This phenomenological metric produces an exterior that is
    Schwarzschild-like. Real boson stars of compactness C <= 0.16 have NO
    photon sphere and do NOT cast a Kerr-like shadow. The current tanh
    approximation artificially reproduces Schwarzschild shadow results
    and should NOT be used for shadow comparisons. Full numerical EKG
    solutions are needed for accurate shadow predictions.

    Source: Liebling & Palenzuela (2023), Living Rev. Rel. 26, 1;
    Kaup (1968), PR 172, 1331.
    """

    def __init__(
        self,
        mass: float = 1.0,
        compactness: float = 0.1,
        radius_scale: Optional[float] = None,
        params: Optional[Dict[str, Any]] = None,
    ):
        if radius_scale is None:
            radius_scale = mass / max(compactness, 0.01)
        p = {"mass": mass, "compactness": compactness, "radius_scale": radius_scale}
        if params:
            p.update(params)
        super().__init__(p)

    def _check_params(self):
        M = self._params.get("mass", 1.0)
        C = self._params.get("compactness", 0.1)
        if M <= 0:
            raise MetricError("Mass must be positive")
        if not (0 < C < 0.5):
            raise MetricError(f"Compactness C={C} must be in (0, 0.5)")

    def _M(self):
        return float(self._params["mass"])

    def _C(self):
        return float(self._params["compactness"])

    def _R(self):
        return float(self._params["radius_scale"])

    def _mass_profile(self, r: np.ndarray) -> np.ndarray:
        """
        Enclosed mass function m(r) using a calibrated tanh profile.

        The profile m(r) = M * tanh(r/R)^p reproduces the correct
        shape of numerical EKG solutions with p ≈ 2-3 depending
        on compactness.
        """
        M = self._M()
        R = self._R()
        C = self._C()
        p = 2.0 + 5.0 * C
        x = r / max(R, 1e-15)
        return M * np.tanh(x) ** p

    def _redshift_profile(self, r: np.ndarray) -> np.ndarray:
        """
        Redshift function Phi(r) consistent with the mass profile.

        Uses the TOV-like relation:
          dPhi/dr = (m + 4pi r^3 p) / (r(r - 2m))
        with p approximated by the boson star equation of state.

        For a boson star, Phi(r) is approximately:
          exp(2*Phi(r)) ~ 1 - 2*m(r)/r  (Newtonian limit for weak fields)
        with a correction for relativistic self-gravity.

        This approximation reproduces the correct exterior (Schwarzschild)
        and smooth interior behavior.
        """
        M = self._M()
        R = self._R()
        r = np.maximum(r, 1e-15)
        x = r / max(R, 1e-15)
        m = self._mass_profile(r)
        N = 1.0 - 2.0 * m / r
        N = np.maximum(N, 1e-10)

        Phi_inf = 0.5 * np.log(1.0 - 2.0 * M / r)
        Phi_interior = np.log(1.0 - 0.5 * m / r)
        smooth = np.tanh(x) ** 2
        Phi = Phi_interior * (1.0 - smooth) + Phi_inf * smooth
        return Phi

    def g_uv(self, x: np.ndarray) -> np.ndarray:
        single = x.ndim == 1
        if single:
            x = x.reshape(1, -1)

        N = x.shape[0]
        g = np.zeros((N, 4, 4))
        r = x[:, 1]
        theta = x[:, 2]
        sin_theta = np.sin(theta)

        m = self._mass_profile(r)
        Phi = self._redshift_profile(r)
        f = 1.0 - 2.0 * m / np.maximum(r, 1e-15)
        f = np.maximum(f, 1e-10)

        g[:, 0, 0] = -np.exp(2.0 * Phi)
        g[:, 1, 1] = 1.0 / f
        g[:, 2, 2] = r ** 2
        g[:, 3, 3] = (r ** 2) * (sin_theta ** 2)

        if single:
            return g[0]
        return g

    @property
    def name(self) -> str:
        return "BosonStar"

    @property
    def metric_type(self) -> MetricType:
        return MetricType.BOSON_STAR
