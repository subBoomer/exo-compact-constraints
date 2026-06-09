"""
Hamiltonian geodesic engine + RayTracer pipeline.

Uses scipy DOP853 with early termination for efficient ray tracing.
"""

import numpy as np
from scipy.integrate import solve_ivp
from typing import Optional, Callable, Dict, Any
from .base import MetricBase, MetricCache, TerminationOracle


class HamiltonianGeodesic:
    """
    Hamiltonian geodesic engine for null geodesics.

    H = 1/2 g^{munu} p_mu p_nu = 0
    dx^mu/dlam = dH/dp_mu = g^{munu} p_nu
    dp_mu/dlam = -dH/dx^mu = -1/2 (d_mu g^{ab}) p_a p_b
    """

    def __init__(self, metric: MetricBase):
        self.metric = metric

    def rhs(self, lam: float, y: np.ndarray) -> np.ndarray:
        return self.metric.hamiltonian_rhs(y)

    def integrate(self, y0: np.ndarray, lambda_max: float = 500.0,
                  rtol: float = 1e-6, atol: float = 1e-9,
                  max_step: float = 1.0,
                  events: Optional[list] = None) -> Dict[str, Any]:
        """Integrate null geodesic from initial state y0."""
        result = solve_ivp(
            fun=self.rhs,
            t_span=(0, lambda_max),
            y0=y0,
            method='DOP853',
            rtol=rtol, atol=atol,
            max_step=max_step,
            events=events,
        )
        return {
            'success': result.success,
            'lambda': result.t,
            'y': result.y,
            'n_steps': len(result.t),
            'events': result.t_events,
        }


class RayTracer:
    """
    Full ray tracing pipeline with Hamiltonian engine + termination oracle.

    Supports metric caching for production speed.
    """

    def __init__(self, metric: MetricBase, r_obs: float = 100.0,
                 theta_obs: float = np.pi / 2, escape_radius: float = 50.0,
                 use_cache: bool = True, cache_grid_r: int = 200,
                 cache_grid_theta: int = 50):
        self.metric = metric
        self.r_obs = r_obs
        self.theta_obs = theta_obs
        self.escape_radius = escape_radius

        if use_cache:
            r_grid = np.logspace(np.log10(2.0), np.log10(escape_radius * 2),
                                 cache_grid_r)
            theta_grid = np.linspace(0, np.pi, cache_grid_theta)
            self.cache = MetricCache(metric, r_grid, theta_grid)
        else:
            self.cache = None

        self.engine = HamiltonianGeodesic(self.cache if self.cache else metric)
        self.terminator = TerminationOracle(metric, escape_radius)

    def _initial_conditions(self, alpha: float, beta: float) -> np.ndarray:
        """Initial state [x^mu, p_mu] for pixel (alpha, beta) on image plane."""
        r0, th0 = self.r_obs, self.theta_obs
        x0 = np.array([0.0, r0, th0, 0.0])

        backend = self.cache if self.cache else self.metric
        ginv = backend.g_inv(x0)

        p_t = -1.0
        p_th = -alpha
        p_ph = -beta

        a = ginv[1, 1]
        c = (ginv[0, 0] * p_t * p_t + ginv[2, 2] * p_th * p_th
             + ginv[3, 3] * p_ph * p_ph + 2.0 * ginv[0, 3] * p_t * p_ph)

        if a > 0 and -c / a > 0:
            p_r = -np.sqrt(-c / a)
        else:
            raise ValueError(
                f"Cannot construct null initial momentum at "
                f"(r={r0:.1f}, th={th0:.2f}): H = 0 not satisfiable "
                f"for alpha={alpha:.2f}, beta={beta:.2f}"
            )

        return np.concatenate([x0, [p_t, p_r, p_th, p_ph]])

    def trace_ray(self, alpha: float, beta: float,
                  lambda_max: float = 1000.0) -> Dict[str, Any]:
        """Trace a single pixel. Returns classification dict."""
        y0 = self._initial_conditions(alpha, beta)
        result = self.engine.integrate(y0, lambda_max=lambda_max)
        Y = result['y']

        for i in range(Y.shape[1]):
            status = self.terminator.check(Y[:, i])
            if status is not None:
                return {
                    'status': status,
                    'alpha': alpha,
                    'beta': beta,
                    'steps': i,
                    'r_final': Y[1, i],
                    'lambda_final': result['lambda'][i],
                }

        r_final = Y[1, -1]
        M = getattr(self.metric, 'M', 1.0)
        if r_final < 5.0 * M:
            return {'status': 'captured', 'alpha': alpha, 'beta': beta,
                    'steps': Y.shape[1], 'r_final': r_final}
        return {'status': 'unknown', 'alpha': alpha, 'beta': beta,
                'steps': Y.shape[1], 'r_final': r_final}
