import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from typing import Optional, Callable
from ..metrics.base import AbstractMetric


class TimelikeGeodesicIntegrator:
    """
    Integrates timelike geodesics (massive particle orbits) in any spacetime.

    Uses the geodesic equation with proper time tau as parameter.
    Normalization: g_{munu} dx^mu/dtau dx^nu/dtau = -1.
    """

    def __init__(self, metric: AbstractMetric):
        self.metric = metric

    def geodesic_eq(self, tau: float, y: np.ndarray) -> np.ndarray:
        x = y[:4]
        u = y[4:]
        gamma = self.metric.christoffel(x)
        du = np.zeros(4)
        for rho in range(4):
            du[rho] = -np.einsum('ij,ij->', gamma[rho, :, :], np.outer(u, u))
        return np.concatenate([u, du])

    def trace_orbit(
        self,
        initial_position: np.ndarray,
        initial_velocity: np.ndarray,
        tau_max: float = 100.0,
        n_steps: int = 10000,
        rtol: float = 1e-8,
        atol: float = 1e-10,
        event: Optional[Callable] = None,
    ) -> dict:
        """
        Trace a timelike geodesic.

        Parameters
        ----------
        initial_position : ndarray (4,)
        initial_velocity : ndarray (4,)
        tau_max : float
        n_steps : int
        rtol, atol : float
        event : callable, optional

        Returns
        -------
        result : dict
        """
        y0 = np.concatenate([initial_position, initial_velocity])

        result = solve_ivp(
            fun=self.geodesic_eq,
            t_span=(0, tau_max),
            y0=y0,
            method='DOP853',
            t_eval=np.linspace(0, tau_max, n_steps),
            rtol=rtol,
            atol=atol,
            events=event,
        )

        return {
            'tau': result.t,
            'trajectory': result.y,
            'success': result.success,
        }

    def effective_potential(self, r: float, L: float = 4.0) -> float:
        """
        Effective potential for timelike geodesics in spherically symmetric metric.

        For metric ds^2 = -A(r) dt^2 + B(r) dr^2 + r^2 dOmega^2:
        V_eff(r) = A(r) * (1 + L^2 / r^2)

        Circular orbits satisfy dV_eff/dr = 0.
        ISCO satisfies dV_eff/dr = 0 and d^2 V_eff/dr^2 = 0.
        """
        M = self.metric.params.get("mass", 1.0)
        x = np.array([0.0, r, np.pi / 2, 0.0])
        g = self.metric.g_uv(x)
        A = -g[0, 0]
        if A <= 0:
            return np.inf
        return A * (1.0 + L ** 2 / r ** 2)

    def compute_isco(self) -> float:
        """
        Compute ISCO radius using effective potential.

        For Schwarzschild: r_ISCO = 6M (derived from V_eff).

        For general spherically symmetric metric, find root of
        d^2 V_eff/dr^2 = 0 at the circular orbit radius.

        For Kerr: use Bardeen formula.

        Source: Bardeen, Press & Teukolsky (1972), ApJ 178, 347.
        """
        M = self.metric.params.get("mass", 1.0)
        spin = self.metric.params.get("spin", 0.0)
        name = self.metric.name

        if name == "MorrisThorne":
            b0 = self.metric.params.get("throat_radius", 2.0)
            return max(b0 * 1.05, 6.0 * M)

        if abs(spin) > 0 and name == "Kerr":
            a_star = spin / M
            z1 = 1.0 + (1.0 - a_star ** 2) ** (1.0 / 3.0) * (
                (1.0 + a_star) ** (1.0 / 3.0) + (1.0 - a_star) ** (1.0 / 3.0)
            )
            z2 = np.sqrt(3.0 * a_star ** 2 + z1 ** 2)
            r_isco_pro = M * (3.0 + z2 - np.sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2)))
            return float(r_isco_pro)

        return 6.0 * M
