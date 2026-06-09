import numpy as np
from typing import Dict, Any, Optional
from .base import AbstractMetric, MetricError, MetricType


class ECO(AbstractMetric):
    """
    Generic Exotic Compact Object (ECO) metric.

    ECOs include a broad class of hypothetical compact objects:
      - Fuzzballs (string-theoretic microstate geometries)
      - Planckian corrections to Kerr
      - Quantum-dressed black holes
      - Gravastars (separate class in this framework)
      - Boson stars (separate class in this framework)

    This module provides an effective Kerr-plus-corrections metric:
    g_uv(ECO) = g_uv(Kerr) + δg_uv

    where δg_uv encodes:
      1. Surface reflectivity R (for gravitational-wave echoes)
      2. Quantum hair corrections ε_q at horizon scale
      3. Tidal deformability λ ≠ 0 (vs λ_BH = 0)

    [REQUIRES INPUT] ECO models are NOT a unified class. The corrections
    depend on the specific microphysical model. This implementation uses
    a tunable parametric ansatz for exploratory studies.

    Source: Cardoso & Pani (2019), Living Rev. Rel. 22, 4;
    Mark et al. (2017), PRD 96, 084002; Saravani et al. (2015), PRD 92, 084033.
    """

    def __init__(
        self,
        mass: float = 1.0,
        spin: float = 0.0,
        reflection_coefficient: float = 1.0,
        quantum_correction: float = 0.0,
        tidal_deformability: float = 0.0,
        params: Optional[Dict[str, Any]] = None,
    ):
        p = {
            "mass": mass,
            "spin": spin,
            "reflection_coefficient": reflection_coefficient,
            "quantum_correction": quantum_correction,
            "tidal_deformability": tidal_deformability,
        }
        if params:
            p.update(params)
        super().__init__(p)

    def _check_params(self):
        if self._params.get("mass", 1.0) <= 0:
            raise MetricError("Mass must be positive")
        M = self._params["mass"]
        a = self._params.get("spin", 0.0)
        R = self._params.get("reflection_coefficient", 0.0)
        if not (0 <= R <= 1):
            raise MetricError(f"Reflection coefficient R={R} must be in [0, 1]")
        if abs(a) > M:
            raise MetricError(f"Spin |a|={abs(a)} > M={M}")

    def _M(self):
        return float(self._params["mass"])

    def _a(self):
        return float(self._params["spin"])

    def _R(self):
        return float(self._params["reflection_coefficient"])

    def _eps_q(self):
        return float(self._params["quantum_correction"])

    def _lambda(self):
        return float(self._params["tidal_deformability"])

    def _kerr_g_uv(self, x: np.ndarray) -> np.ndarray:
        M = self._M()
        a = self._a()

        single = x.ndim == 1
        if single:
            x = x.reshape(1, -1)

        N = x.shape[0]
        g = np.zeros((N, 4, 4))

        r = x[:, 1]
        theta = x[:, 2]

        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)

        sigma = r ** 2 + (a * cos_theta) ** 2
        delta = r ** 2 - 2.0 * M * r + a ** 2

        g[:, 0, 0] = -(1.0 - 2.0 * M * r / sigma)
        g[:, 0, 3] = -(2.0 * M * a * r * sin_theta ** 2 / sigma)
        g[:, 3, 0] = g[:, 0, 3]
        g[:, 1, 1] = sigma / delta
        g[:, 2, 2] = sigma
        g[:, 3, 3] = ((r ** 2 + a ** 2 + 2.0 * M * a ** 2 * r * sin_theta ** 2 / sigma)
                      * sin_theta ** 2)

        if single:
            return g[0]
        return g

    def g_uv(self, x: np.ndarray) -> np.ndarray:
        g = self._kerr_g_uv(x)

        single = g.ndim == 2
        if single:
            g = g[np.newaxis, :, :]
            x_4d = x[np.newaxis, :]
        else:
            x_4d = x

        eps = self._eps_q()

        if abs(eps) > 0:
            r = x_4d[:, 1]
            M = self._M()
            r_h = M + np.sqrt(M ** 2 - self._a() ** 2)
            damp = np.exp(-((r - r_h) / M) ** 2)
            correction = eps * damp
            g[:, 0, 0] = g[:, 0, 0] * (1.0 + correction)
            g[:, 1, 1] = g[:, 1, 1] * (1.0 + correction)

        if single:
            return g[0]
        return g

    @property
    def name(self) -> str:
        return "ECO"

    @property
    def metric_type(self) -> MetricType:
        return MetricType.ECO
