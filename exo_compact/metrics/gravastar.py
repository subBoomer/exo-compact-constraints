import numpy as np
from typing import Dict, Any, Optional
from .base import AbstractMetric, MetricError, MetricType


class Gravastar(AbstractMetric):
    """
    Gravitational Vacuum Star (Gravastar) metric.

    Gravastars are ultra-compact objects with a de Sitter interior,
    a thin shell of stiff matter, and a Schwarzschild exterior.

    Three-region structure:
      Region I (interior):  r < r1,  de Sitter: p = -ρ, ds² = -(1 - r²/Λ²) dt² + ...
      Region II (shell):    r1 < r < r2, thin shell with surface density
      Region III (exterior): r > r2, Schwarzschild

    [ASSUMPTION] This implementation uses a simplified two-parameter
    model with a smoothed transition between the interior and exterior.
    The thin-shell discontinuity is approximated by a hyperbolic tangent
    interpolation.

    Source: Mazur & Mottola (2001), PRD 64, 104022;
    Visser & Wiltshire (2004), CQG 21, 1135.
    """

    def __init__(
        self,
        mass: float = 1.0,
        interior_radius: float = 2.1,
        shell_width: float = 0.1,
        mode: str = "safe",
        params: Optional[Dict[str, Any]] = None,
    ):
        p = {
            "mass": mass,
            "interior_radius": interior_radius,
            "shell_width": shell_width,
            "mode": mode,
        }
        if params:
            p.update(params)
        super().__init__(p)

    def _check_params(self):
        if self._params.get("mass", 1.0) <= 0:
            raise MetricError("Mass must be positive")
        r1 = self._params.get("interior_radius", 2.0)
        M = self._params["mass"]
        mode = self._params.get("mode", "safe")
        if mode == "physical":
            if not (0 < r1 < 2.0 * M):
                raise MetricError(
                    f"Mazur-Mottola physical mode requires 0 < r1 < 2M. "
                    f"Got r1={r1}, 2M={2*M}"
                )
        else:
            if r1 <= 2.0 * M:
                raise MetricError(
                    f"Safe mode requires r1 > 2M. Got r1={r1}, 2M={2*M}. "
                    f"Use mode='physical' for Mazur-Mottola gravastar."
                )

    def _M(self):
        return float(self._params["mass"])

    def _r1(self):
        return float(self._params["interior_radius"])

    def _w(self):
        return float(self._params["shell_width"])

    def _r2(self):
        return self._r1() + self._w()

    def _Lambda(self):
        """
        Effective cosmological constant for de Sitter interior.

        For matching of de Sitter (interior) to Schwarzschild (exterior):
        Lambda = 3 / r_shell^2 where r_shell is the shell radius.
        This ensures g_tt continuity across the shell.

        Source: Mazur & Mottola (2001), PRD 64, 104022, Eq. (15).
        """
        r_shell = self._r1()
        return 3.0 / max(r_shell ** 2, 1e-15)

    def g_uv(self, x: np.ndarray) -> np.ndarray:
        M = self._M()
        r1 = self._r1()
        r2 = self._r2()
        L = self._Lambda()

        single = x.ndim == 1
        if single:
            x = x.reshape(1, -1)

        N = x.shape[0]
        g = np.zeros((N, 4, 4))

        r = x[:, 1]
        theta = x[:, 2]

        sin_theta = np.sin(theta)

        def interpolate_f(ri):
            if ri < r1:
                return 1.0 - ri ** 2 / L ** 2
            elif ri < r2:
                x = (ri - r1) / max(self._w(), 1e-15)
                w_de = 1.0 - ri ** 2 / L ** 2
                w_sch = 1.0 - 2.0 * M / ri
                smooth = 0.5 * (1.0 + np.tanh(5.0 * (x - 0.5)))
                return w_de * (1.0 - smooth) + w_sch * smooth
            else:
                return 1.0 - 2.0 * M / ri

        f = np.array([interpolate_f(ri) for ri in r])
        f = np.maximum(f, 1e-15)

        g[:, 0, 0] = -f
        g[:, 1, 1] = 1.0 / f
        g[:, 2, 2] = r ** 2
        g[:, 3, 3] = (r ** 2) * (sin_theta ** 2)

        if single:
            return g[0]
        return g

    @property
    def name(self) -> str:
        return "Gravastar"

    @property
    def metric_type(self) -> MetricType:
        return MetricType.GRAVASTAR
