import numpy as np
from typing import Dict, Any, Optional
from .base import AbstractMetric, MetricError, MetricType


class Kerr(AbstractMetric):
    """
    Kerr black hole metric in Boyer-Lindquist coordinates.

    ds² = -(1 - 2Mr/Σ) dt² - (4Mar sin²θ/Σ) dt dφ + (Σ/Δ) dr²
          + Σ dθ² + (r² + a² + 2Ma²r sin²θ/Σ) sin²θ dφ²

    where Σ = r² + a² cos²θ, Δ = r² - 2Mr + a²

    Source: Kerr (1963), Boyer & Lindquist (1967).
    """

    def __init__(self, mass: float = 1.0, spin: float = 0.0,
                 params: Optional[Dict[str, Any]] = None):
        p = {"mass": mass, "spin": spin}
        if params:
            p.update(params)
        super().__init__(p)

    def _check_params(self):
        if self._params.get("mass", 1.0) <= 0:
            raise MetricError("Mass must be positive")
        M = self._params["mass"]
        a = self._params.get("spin", 0.0)
        if abs(a) > M:
            raise MetricError(f"Spin |a| = {abs(a)} > M = {M} violates cosmic censorship")

    def _M(self):
        return float(self._params["mass"])

    def _a(self):
        return float(self._params["spin"])

    def g_uv(self, x: np.ndarray) -> np.ndarray:
        M = self._M()
        a = self._a()

        single = x.ndim == 1
        if single:
            x = x.reshape(1, -1)

        N = x.shape[0]
        g = np.zeros((N, 4, 4))

        t = x[:, 0]
        r = x[:, 1]
        theta = x[:, 2]
        phi = x[:, 3]

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

    @property
    def name(self) -> str:
        return "Kerr"

    @property
    def metric_type(self) -> MetricType:
        return MetricType.KERR
