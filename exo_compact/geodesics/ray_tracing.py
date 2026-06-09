"""
Hamiltonian geodesic ray tracer for stationary metrics.

Uses the Hamiltonian formulation (H = 1/2 g^mu^nu p_mu p_nu = 0) instead
of the Christoffel-based geodesic equation. This avoids coordinate singularity
artifacts near photon spheres and horizons where g_rr diverges.

State: (x^mu, p_mu) where p_mu are COVARIANT momentum components.
Equations: dx^mu/dlam = g^mu^nu p_nu,  dp_mu/dlam = -1/2 (d_mu g^alpha^beta) p_alpha p_beta

Reference: Carter (1968), PR 174, 1559; Gralla, Lupsasca, Marrone (2019), PRD 100.
"""

import numpy as np
from typing import Dict, Any, Tuple


class RayTracer:
    """
    Hamiltonian geodesic ray tracer for null geodesics.

    Solves dx^mu/dlam = g^mu^nu p_nu, dp_mu/dlam = -1/2 (d_mu g^alpha^beta) p_alpha p_beta
    using fixed-step RK4. The Hamiltonian H = 1/2 g^mu^nu p_mu p_nu is conserved.
    """

    def __init__(self, metric, r_obs=200.0, theta_obs=np.pi / 2, dx=1e-6):
        self.metric = metric
        self.r_obs = r_obs
        self.theta_obs = theta_obs
        self.b0 = metric.params.get("throat_radius", 2.0)
        self.dx = dx

    def _inverse_metric(self, x):
        """Inverse metric g^mu^nu at position x via matrix inversion."""
        g = self.metric.g_uv(x)
        return np.linalg.inv(g)

    def _metric_derivatives(self, x):
        """Partial derivatives d_mu g_ab at x via central finite differences."""
        dg = np.zeros((4, 4, 4))
        for mu in range(4):
            xp = x.copy()
            xm = x.copy()
            xp[mu] += self.dx
            xm[mu] -= self.dx
            gp = self.metric.g_uv(xp)
            gm = self.metric.g_uv(xm)
            dg[mu] = (gp - gm) / (2.0 * self.dx)
        return dg

    def _inverse_metric_derivatives(self, x):
        """Partial derivatives d_mu g^ab via dg^ab = -g^ac g^bd dg_cd."""
        g_inv = self._inverse_metric(x)
        dg = self._metric_derivatives(x)
        dg_inv = np.zeros((4, 4, 4))
        for mu in range(4):
            dg_inv[mu] = -np.einsum('ac,bd,cd->ab', g_inv, g_inv, dg[mu])
        return dg_inv

    def hamiltonian_rhs(self, y):
        """RHS of Hamiltonian equations. y = [x^mu, p_mu]."""
        x = y[:4]
        p = y[4:]

        g_inv = self._inverse_metric(x)
        dx = np.einsum('ab,b->a', g_inv, p)

        dg_inv = self._inverse_metric_derivatives(x)
        dp = np.zeros(4)
        for mu in range(4):
            dp[mu] = -0.5 * np.einsum('ab,a,b->', dg_inv[mu], p, p)

        return np.concatenate([dx, dp])

    def rk4_step(self, y, h):
        k1 = self.hamiltonian_rhs(y)
        k2 = self.hamiltonian_rhs(y + 0.5 * h * k1)
        k3 = self.hamiltonian_rhs(y + 0.5 * h * k2)
        k4 = self.hamiltonian_rhs(y + h * k3)
        return y + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    def _initial_conditions(self, alpha, beta):
        """
        Initial state [x^mu, p_mu] for pixel (alpha, beta).

        Sets p_t = -1 (unit energy), p_theta = -alpha, p_phi = -beta
        (impact parameter mapping), then solves for p_r from H = 0.
        The affine parameter is rescaled so that dr/dlam ~ O(1).
        """
        x0 = np.array([0.0, self.r_obs, self.theta_obs, 0.0])
        ginv = self._inverse_metric(x0)

        p_t = -1.0
        p_theta = -alpha
        p_phi = -beta

        a = ginv[1, 1]
        c = (ginv[0, 0] * p_t * p_t + ginv[2, 2] * p_theta * p_theta
             + ginv[3, 3] * p_phi * p_phi + 2.0 * ginv[0, 3] * p_t * p_phi)

        if a > 0 and -c / a > 0:
            p_r = -np.sqrt(-c / a)
        elif a != 0:
            p_r = -np.sqrt(max(-c / a, 0.0))
        else:
            p_r = -1.0

        p = np.array([p_t, p_r, p_theta, p_phi])
        H = 0.5 * np.einsum('ab,a,b->', ginv, p, p)

        if abs(H) > 1e-12:
            p /= np.sqrt(2.0 * abs(H))

        return np.concatenate([x0, p])

    def _hamiltonian(self, y):
        """H = 1/2 g^mu^nu p_mu p_nu. Should be ~0 for null geodesics."""
        x, p = y[:4], y[4:]
        g_inv = self._inverse_metric(x)
        return 0.5 * np.einsum('ab,a,b->', g_inv, p, p)

    def trace_pixel(self, alpha, beta, lambda_max=500.0, h=0.5):
        """Trace null geodesic for pixel (alpha, beta). Returns classification dict."""
        is_bh = self.metric.name in ("Schwarzschild", "Kerr")
        if is_bh:
            r_horizon = 2.5 * self.metric.params.get("mass", 1.0)
        else:
            r_horizon = self.b0 * 1.05

        y = self._initial_conditions(alpha, beta)
        r_start = self.r_obs
        r_current = r_start
        went_in = False
        n_steps = int(lambda_max / h)

        for step in range(n_steps):
            y_new = self.rk4_step(y, h)
            if not np.all(np.isfinite(y_new)):
                break

            y = y_new
            r_new = y[1]

            if r_new < r_horizon:
                return {"status": "throat_cross", "steps": step, "r_min": r_new}
            if r_new < r_current:
                went_in = True
            if went_in and r_new > r_current and r_new > 0.9 * r_start:
                return {"status": "escape", "steps": step, "r_min": r_current}
            if r_new > 1.5 * r_start:
                return {"status": "escape", "steps": step, "r_min": r_current}

            r_current = r_new

        if r_current < r_horizon * 2:
            return {"status": "throat_cross", "steps": n_steps, "r_min": r_current}
        return {"status": "trapped", "steps": n_steps, "r_min": r_current}

    def scan_equatorial(self, alpha_range, n_alpha=40, beta=0.0,
                        lambda_max=500.0, h=0.5, verbose=True):
        """Scan equatorial plane to find shadow boundary."""
        alphas = np.linspace(alpha_range[0], alpha_range[1], n_alpha)
        classifications = []
        r_mins = []

        for a in alphas:
            result = self.trace_pixel(a, beta, lambda_max=lambda_max, h=h)
            classifications.append(result["status"])
            r_mins.append(result["r_min"])

        if verbose:
            n_escape = sum(1 for c in classifications if c == "escape")
            n_cross = sum(1 for c in classifications if c == "throat_cross")
            n_trapped = sum(1 for c in classifications if c == "trapped")
            print(f"    {n_escape} escape, {n_cross} throat-cross, {n_trapped} trapped")

        boundary_alphas = []
        for i in range(len(alphas) - 1):
            if classifications[i] != classifications[i + 1]:
                boundary_alphas.append((alphas[i] + alphas[i + 1]) / 2)

        return {
            "alpha": alphas,
            "classifications": classifications,
            "r_min": np.array(r_mins),
            "boundary_alpha": boundary_alphas,
        }

    def teo_shadow_radius(self, n_alpha=40, lambda_max=500.0, h=0.5, verbose=True):
        """Compute Teo equatorial shadow radius via ray tracing."""
        scan = self.scan_equatorial(
            (1.0, 8.0), n_alpha=n_alpha, beta=0.0,
            lambda_max=lambda_max, h=h, verbose=verbose
        )
        bounds = scan["boundary_alpha"]
        if len(bounds) >= 1:
            return float(np.mean(bounds))
        return 2.0 * self.b0
