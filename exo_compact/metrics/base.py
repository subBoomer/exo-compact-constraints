from abc import ABC, abstractmethod
from enum import Enum, auto
import numpy as np
from typing import Dict, Any, Optional


class MetricType(Enum):
    SCHWARZSCHILD = auto()
    KERR = auto()
    MORRIS_THORNE = auto()
    ROTATING_WORMHOLE = auto()
    BOSON_STAR = auto()
    GRAVASTAR = auto()
    ECO = auto()


class MetricError(Exception):
    pass


class AbstractMetric(ABC):
    """
    Metric-first architecture: all observables flow from spacetime geometry.

    Every compact object is defined by its metric tensor g_μν(x; params).
    Concrete subclasses implement g_uv() and g_uv_contra() in their
    chosen coordinate system. Higher-order geometric quantities
    (Christoffel symbols, Ricci, geodesic equations) are derived
    automatically via finite-differencing or sympy, with caching.

    Coordinate convention: x = (t, r, theta, phi)
    Sign convention: (-, +, +, +)
    Units: G = c = 1 (geometric units)
    """

    def __init__(self, params: Optional[Dict[str, Any]] = None):
        self._params = params or {}
        self._check_params()

    @abstractmethod
    def _check_params(self):
        pass

    @abstractmethod
    def g_uv(self, x: np.ndarray) -> np.ndarray:
        """
        Covariant metric tensor g_μν at event x.

        Parameters
        ----------
        x : ndarray of shape (4,) or (N, 4)
            Coordinates (t, r, theta, phi)

        Returns
        -------
        g : ndarray of shape (4, 4) or (N, 4, 4)
        """
        pass

    def g_uv_contra(self, x: np.ndarray) -> np.ndarray:
        """
        Contravariant metric tensor g^μν at event x.
        Default: matrix inverse of g_uv.
        """
        g = self.g_uv(x)
        single_point = g.ndim == 2
        if single_point:
            g = g[np.newaxis, :, :]
        result = np.array([np.linalg.inv(g[i]) for i in range(g.shape[0])])
        return result[0] if single_point else result

    def christoffel(self, x: np.ndarray, dx: float = 1e-6) -> np.ndarray:
        """
        Christoffel symbols Γ^ρ_μν via finite-difference derivatives of g_uv.

        Parameters
        ----------
        x : ndarray of shape (4,)
        dx : float, step size for numerical differentiation

        Returns
        -------
        Gamma : ndarray of shape (4, 4, 4)
            Gamma[rho, mu, nu] = Γ^ρ_μν
        """
        g = self.g_uv(x)
        g_inv = np.linalg.inv(g)

        Gamma = np.zeros((4, 4, 4))
        for mu in range(4):
            x_plus = x.copy()
            x_minus = x.copy()
            x_plus[mu] += dx
            x_minus[mu] -= dx
            g_plus = self.g_uv(x_plus)
            g_minus = self.g_uv(x_minus)

            for nu in range(4):
                for sigma in range(4):
                    dg_nu_sigma = (g_plus[nu, sigma] - g_minus[nu, sigma]) / (2 * dx)
                    for rho in range(4):
                        g_inv_nu = g_inv[rho, :]
                        Gamma[rho, mu, nu] += 0.5 * g_inv_nu[sigma] * dg_nu_sigma

        for mu in range(4):
            for nu in range(4):
                for sigma in range(4):
                    x_plus_nu = x.copy()
                    x_minus_nu = x.copy()
                    x_plus_nu[nu] += dx
                    x_minus_nu[nu] -= dx
                    g_plus_nu = self.g_uv(x_plus_nu)
                    g_minus_nu = self.g_uv(x_minus_nu)

                    x_plus_sigma = x.copy()
                    x_minus_sigma = x.copy()
                    x_plus_sigma[sigma] += dx
                    x_minus_sigma[sigma] -= dx
                    g_plus_sigma = self.g_uv(x_plus_sigma)
                    g_minus_sigma = self.g_uv(x_minus_sigma)

                    dg_mu_nu_sigma = (g_plus_nu[mu, sigma] - g_minus_nu[mu, sigma]) / (2 * dx)
                    dg_mu_sigma_nu = (g_plus_sigma[mu, nu] - g_minus_sigma[mu, nu]) / (2 * dx)

                    for rho in range(4):
                        g_inv_rho = g_inv[rho, :]
                        Gamma[rho, mu, nu] += 0.5 * g_inv_rho[mu] * (dg_mu_nu_sigma + dg_mu_sigma_nu)

        return Gamma

    def geodesic_rhs(self, tau: float, y: np.ndarray) -> np.ndarray:
        """
        Right-hand side of the geodesic equation d²x^ρ/dτ² + Γ^ρ_μν dx^μ/dτ dx^ν/dτ = 0.

        Parameters
        ----------
        tau : float, affine parameter (not used, included for ODE solver interface)
        y : ndarray of shape (8,)
            First 4: x^μ (t, r, theta, phi)
            Last 4: dx^μ/dτ

        Returns
        -------
        dydt : ndarray of shape (8,)
        """
        x = y[:4]
        u = y[4:]
        gamma = self.christoffel(x)
        du = np.zeros(4)
        for rho in range(4):
            du[rho] = -np.einsum('mn,mn->', gamma[rho, :, :], np.outer(u, u))
        return np.concatenate([u, du])

    @property
    def params(self) -> Dict[str, Any]:
        return dict(self._params)

    @property
    @abstractmethod
    def name(self) -> str:
        pass

    @property
    @abstractmethod
    def metric_type(self) -> MetricType:
        pass

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self._params})"
