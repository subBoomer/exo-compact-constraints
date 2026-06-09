import numpy as np
from typing import Dict, Any, Optional, Callable
from .base import AbstractMetric, MetricError, MetricType


class MorrisThorne(AbstractMetric):
    """
    Morris-Thorne traversable wormhole metric in Schwarzschild-like coordinates.

    ds² = -e^{2Φ(r)} dt² + dr²/(1 - b(r)/r) + r² (dθ² + sin²θ dφ²)

    Default parameters: Schwarzschild-like exterior with flaring-out throat.

    Source: Morris & Thorne (1988), Am. J. Phys. 56, 395.
    """

    def __init__(
        self,
        throat_radius: float = 2.0,
        mass: float = 1.0,
        shape_function: Optional[Callable] = None,
        redshift_function: Optional[Callable] = None,
        params: Optional[Dict[str, Any]] = None
    ):
        p = {
            "throat_radius": throat_radius,
            "mass": mass,
        }
        if params:
            p.update(params)

        self._shape_func = shape_function or self._default_shape
        self._redshift_func = redshift_function or self._default_redshift

        super().__init__(p)

    def _check_params(self):
        b0 = self._params.get("throat_radius", 2.0)
        M = self._params.get("mass", 1.0)
        if b0 <= 0:
            raise MetricError("Throat radius must be positive")
        if M <= 0:
            raise MetricError("Mass must be positive")

    def _default_shape(self, r: np.ndarray) -> np.ndarray:
        """
        Default shape function: b(r) = b0²/r (flaring-out, asymptotically flat).
        Satisfies b(r)/r <= 1 at throat r=b0, flaring-out condition b'(r) < b(r)/r.
        """
        b0 = float(self._params["throat_radius"])
        return np.where(r >= b0, b0 ** 2 / r, 0.0)

    def _default_redshift(self, r: np.ndarray) -> np.ndarray:
        """
        Default redshift function: Φ(r) = 0 (zero tidal force).
        Redshift function governs gravitational time dilation.
        """
        return np.zeros_like(r)

    def _throat(self):
        return float(self._params["throat_radius"])

    def g_uv(self, x: np.ndarray) -> np.ndarray:
        b0 = self._throat()
        single = x.ndim == 1
        if single:
            x = x.reshape(1, -1)

        N = x.shape[0]
        g = np.zeros((N, 4, 4))

        r = x[:, 1]
        theta = x[:, 2]

        Phi = self._redshift_func(r)
        b_val = self._shape_func(r)

        f = 1.0 - b_val / np.maximum(r, 1e-15)
        f = np.maximum(f, 1e-15)
        sin_theta = np.sin(theta)

        g[:, 0, 0] = -np.exp(2.0 * Phi)
        g[:, 1, 1] = 1.0 / f
        g[:, 2, 2] = r ** 2
        g[:, 3, 3] = (r ** 2) * (sin_theta ** 2)

        if single:
            return g[0]
        return g

    @property
    def name(self) -> str:
        return "MorrisThorne"

    @property
    def metric_type(self) -> MetricType:
        return MetricType.MORRIS_THORNE
