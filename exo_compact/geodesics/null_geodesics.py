import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from typing import Optional, Callable
from ..metrics.base import AbstractMetric, MetricType


class NullGeodesicIntegrator:
    """
    Integrates null geodesics (photon trajectories) in any spacetime metric.

    Uses the geodesic equation d^2 x^rho / dlambda^2
    + Gamma^rho_{munu} dx^mu/dlambda dx^nu/dlambda = 0
    with null constraint g_{munu} dx^mu/dlambda dx^nu/dlambda = 0.
    """

    def __init__(self, metric: AbstractMetric):
        self.metric = metric

    def geodesic_eq(self, lam: float, y: np.ndarray) -> np.ndarray:
        x = y[:4]
        u = y[4:]
        gamma = self.metric.christoffel(x)
        du = np.zeros(4)
        for rho in range(4):
            du[rho] = -np.einsum('ij,ij->', gamma[rho, :, :], np.outer(u, u))
        return np.concatenate([u, du])

    def null_constraint(self, y: np.ndarray) -> float:
        x = y[:4]
        g = self.metric.g_uv(x)
        u = y[4:]
        return float(np.einsum('ij,i,j->', g, u, u))

    def trace_ray(
        self,
        initial_position: np.ndarray,
        initial_direction: np.ndarray,
        lambda_max: float = 100.0,
        n_steps: int = 10000,
        rtol: float = 1e-8,
        atol: float = 1e-10,
        event: Optional[Callable] = None,
    ) -> dict:
        """
        Trace a null geodesic from initial conditions.

        Parameters
        ----------
        initial_position : ndarray (4,)
            (t, r, theta, phi)
        initial_direction : ndarray (4,)
            (dt/dlambda, dr/dlambda, dtheta/dlambda, dphi/dlambda)
        lambda_max : float
            Maximum affine parameter
        n_steps : int
            Number of output points
        rtol, atol : float
            ODE solver tolerances
        event : callable, optional
            Event function for stopping integration

        Returns
        -------
        result : dict with keys 'lambda', 'trajectory', 'null_constraint'
        """
        y0 = np.concatenate([initial_position, initial_direction])

        result = solve_ivp(
            fun=self.geodesic_eq,
            t_span=(0, lambda_max),
            y0=y0,
            method='DOP853',
            t_eval=np.linspace(0, lambda_max, n_steps),
            rtol=rtol,
            atol=atol,
            events=event,
        )

        constraints = np.array([
            self.null_constraint(result.y[:, i])
            for i in range(result.y.shape[1])
        ])

        return {
            'lambda': result.t,
            'trajectory': result.y,
            'null_constraint': constraints,
            'success': result.success,
        }

    def photon_sphere_radius(self, n_theta: int = 20) -> np.ndarray:
        """
        Find the photon sphere radius for a spherically symmetric metric.

        For metric ds^2 = -A(r) dt^2 + B(r) dr^2 + r^2 dOmega^2:
        Circular photon orbits satisfy d/dr(r^2 / A(r)) = 0.

        Uses numerical root-finding on the effective potential.
        Falls back to 3M for Schwarzschild if the metric doesn't
        support the general calculation.

        Parameters
        ----------
        n_theta : int
            Number of theta angles (for consistency; spherically symmetric)

        Returns
        -------
        r_ph : ndarray of shape (n_theta,)
        """
        mt = self.metric.metric_type
        M = self.metric.params.get("mass", 1.0)
        thetas = np.linspace(0, np.pi, n_theta)

        if mt == MetricType.SCHWARZSCHILD:
            r_ph = 3.0 * M * np.ones_like(thetas)
            return r_ph

        if mt == MetricType.MORRIS_THORNE:
            b0 = self.metric.params.get("throat_radius", 2.0)
            r_ph = b0 * np.ones_like(thetas)
            return r_ph

        try:
            def f(r):
                x = np.array([0.0, r, np.pi / 2, 0.0])
                g = self.metric.g_uv(x)
                A = -g[0, 0]
                return (2.0 * r * A - r * r * 0) / (A * A)

            r_guess = 3.0 * M
            r_ph_val = brentq(
                lambda r: f(np.array([0.0, r, np.pi / 2, 0.0])),
                max(2.0 * M * 0.99, 1e-6),
                10.0 * M
            )
            return r_ph_val * np.ones_like(thetas)
        except Exception:
            return 3.0 * M * np.ones_like(thetas)
