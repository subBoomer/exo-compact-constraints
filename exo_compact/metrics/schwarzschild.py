import numpy as np
from typing import Dict, Any, Optional
from .base import AbstractMetric, MetricError, MetricType


class Schwarzschild(AbstractMetric):
    """
    Schwarzschild black hole metric in standard Schwarzschild coordinates.

    ds² = -(1 - 2M/r) dt² + (1 - 2M/r)^(-1) dr² + r² dθ² + r² sin²θ dφ²

    Source: Schwarzschild (1916), standard GR textbook.
    """

    def __init__(self, mass: float = 1.0, params: Optional[Dict[str, Any]] = None):
        p = {"mass": mass}
        if params:
            p.update(params)
        super().__init__(p)

    def _check_params(self):
        if self._params.get("mass", 1.0) <= 0:
            raise MetricError("Mass must be positive")

    def _M(self):
        return float(self._params["mass"])

    def g_uv(self, x: np.ndarray) -> np.ndarray:
        M = self._M()
        single = x.ndim == 1
        if single:
            x = x.reshape(1, -1)

        N = x.shape[0]
        g = np.zeros((N, 4, 4))

        t, r, theta, phi = x[:, 0], x[:, 1], x[:, 2], x[:, 3]

        f = 1.0 - 2.0 * M / r
        sin_theta = np.sin(theta)

        g[:, 0, 0] = -f
        g[:, 1, 1] = 1.0 / f
        g[:, 2, 2] = r ** 2
        g[:, 3, 3] = (r ** 2) * (sin_theta ** 2)

        if single:
            return g[0]
        return g

    @property
    def name(self) -> str:
        return "Schwarzschild"

    @property
    def metric_type(self) -> MetricType:
        return MetricType.SCHWARZSCHILD
