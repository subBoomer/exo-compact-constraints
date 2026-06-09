import numpy as np
from typing import Dict, Any, Optional, Callable
from .base import AbstractMetric, MetricError, MetricType


class RotatingWormhole(AbstractMetric):
    """
    Parametric, metric-agnostic rotating wormhole framework.

    Uses the general stationary axisymmetric metric ansatz:
    ds² = -N² dt² + A² dr² + B² dθ² + C² (dφ - ω dt)²

    where N, A, B, C, ω are functions of (r, θ) determined by
    pluggable component functions.

    Known special cases:
      - Teo wormhole (Teo 1998, PRD 58, 024014)
      - Kashargin & Sushkov (2008, grav. & cosm.)
      - Generalised rotating thin-shell wormholes

    The default implementation uses a Teo-like ansatz:
      N(r) = exp(-M/r)
      A(r) = 1 / sqrt(1 - b(r)/r)
      B(r,θ) = r * f(r,θ)
      C(r,θ) = r * f(r,θ) * sin(θ)
      ω(r,θ) = 2J / r³ (Lense-Thirring-like frame dragging, a = J/M)

    [REQUIRES INPUT] The rotating wormhole metric does NOT follow
    from GR + exotic matter uniquely. Different shape/redshift/rotation
    combinations produce different observables. This module provides the
    parametric framework; specific models must be chosen by the user or
    systematically enumerated.

    Source: Teo (1998), PRD 58, 024014; Morris & Thorne (1988).
    """

    def __init__(
        self,
        mass: float = 1.0,
        throat_radius: float = 2.0,
        spin: float = 0.0,
        shape_function: Optional[Callable] = None,
        rotation_profile: Optional[Callable] = None,
        component_N: Optional[Callable] = None,
        component_A: Optional[Callable] = None,
        component_B: Optional[Callable] = None,
        component_C: Optional[Callable] = None,
        component_omega: Optional[Callable] = None,
        params: Optional[Dict[str, Any]] = None,
    ):
        p = {
            "mass": mass,
            "throat_radius": throat_radius,
            "spin": spin,
        }
        if params:
            p.update(params)

        self._component_N = component_N or self._default_N
        self._component_A = component_A or self._default_A
        self._component_B = component_B or self._default_B
        self._component_C = component_C or self._default_C
        self._component_omega = component_omega or self._default_omega

        super().__init__(p)

    def _check_params(self):
        M = self._params.get("mass", 1.0)
        b0 = self._params.get("throat_radius", 2.0)
        if M <= 0:
            raise MetricError("Mass must be positive")
        if b0 <= 0:
            raise MetricError("Throat radius must be positive")

    def _M(self):
        return float(self._params["mass"])

    def _b0(self):
        return float(self._params["throat_radius"])

    def _a(self):
        return float(self._params["spin"])

    def _default_shape(self, r: np.ndarray) -> np.ndarray:
        b0 = self._b0()
        return np.where(r >= b0, b0 ** 2 / r, 0.0)

    def _default_N(self, r: np.ndarray, theta: np.ndarray) -> np.ndarray:
        """
        Default lapse function: N(r) = exp(-M/r) (asymptotically Schwarzschild).
        [ASSUMPTION] Arbitrary choice; does NOT derive from any specific
        energy-momentum tensor.
        """
        M = self._M()
        return np.exp(-M / np.maximum(r, 1e-15))

    def _default_A(self, r: np.ndarray, theta: np.ndarray) -> np.ndarray:
        """
        Default radial metric function: A² = 1/(1 - b(r)/r).
        Uses Morris-Thorne shape function.
        Note: b(r) is the shape function, not throat radius.
        """
        b_val = self._default_shape(r)
        f = 1.0 - b_val / np.maximum(r, 1e-15)
        f = np.maximum(f, 1e-15)
        return 1.0 / np.sqrt(f)

    def _default_B(self, r: np.ndarray, theta: np.ndarray) -> np.ndarray:
        return r.copy()

    def _default_C(self, r: np.ndarray, theta: np.ndarray) -> np.ndarray:
        sin_theta = np.sin(theta)
        return max(r) * np.ones_like(r) * sin_theta

    def _default_omega(self, r: np.ndarray, theta: np.ndarray) -> np.ndarray:
        """
        Default frame-dragging: ω(r) ≈ 2J/r³ (Kerr-like asymptotic falloff).
        [ASSUMPTION] Uses Lense-Thirring falloff; not a solution of
        the Einstein equations for wormhole matter.
        """
        J = self._a() * self._M()
        return 2.0 * J / np.maximum(r, 1e-15) ** 3

    def set_teo_metric(self):
        """
        Configure the metric as a Teo (1998) rotating wormhole.

        Teo metric components:
          N(r) = exp(-r0/r), r0 = throat radius
          b(r) = b0 (constant throat radius)
          K(r,theta) = 1
          omega(r) = 2J/r^3 (frame dragging)

        ADM mass: M = r0
        Spin parameter: a = J/M^2
        |a| can exceed 1 (unlike Kerr).

        Source: Teo (1998), PRD 58, 024014.
        """
        r0 = self._b0()
        J = self._a() * r0

        self._component_N = lambda r, theta: np.exp(-r0 / np.maximum(r, 1e-15))
        self._component_A = lambda r, theta: 1.0 / np.sqrt(
            1.0 - self._b0() / np.maximum(r, 1e-15))
        self._component_B = lambda r, theta: r.copy()
        self._component_C = lambda r, theta: r * np.sin(theta)
        self._component_omega = lambda r, theta: 2.0 * J / np.maximum(r, 1e-15) ** 3

    def g_uv(self, x: np.ndarray) -> np.ndarray:
        single = x.ndim == 1
        if single:
            x = x.reshape(1, -1)

        N_pts = x.shape[0]
        g = np.zeros((N_pts, 4, 4))

        r = x[:, 1]
        theta = x[:, 2]

        N_val = self._component_N(r, theta)
        A_val = self._component_A(r, theta)
        B_val = self._component_B(r, theta)
        C_val = self._component_C(r, theta)
        omega_val = self._component_omega(r, theta)

        g_tt = -N_val ** 2 + C_val ** 2 * omega_val ** 2
        g_tp = -C_val ** 2 * omega_val
        g_pp = C_val ** 2

        g[:, 0, 0] = g_tt
        g[:, 0, 3] = g_tp
        g[:, 3, 0] = g_tp
        g[:, 1, 1] = A_val ** 2
        g[:, 2, 2] = B_val ** 2
        g[:, 3, 3] = g_pp

        if single:
            return g[0]
        return g

    @property
    def name(self) -> str:
        return "RotatingWormhole"

    @property
    def metric_type(self) -> MetricType:
        return MetricType.ROTATING_WORMHOLE
