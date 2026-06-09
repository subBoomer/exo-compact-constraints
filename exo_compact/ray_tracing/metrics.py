"""
Concrete metric implementations for the ray tracing framework.

Each metric implements the MetricBase interface used by the
Hamiltonian geodesic engine and caching layer.
"""

import numpy as np
from .base import MetricBase


class SchwarzschildMetric(MetricBase):
    """Schwarzschild BH: g_{muu} in standard coordinates. M = mass."""

    def __init__(self, mass: float = 1.0):
        self.M = mass
        self.a = 0.0

    def g(self, x: np.ndarray) -> np.ndarray:
        M = self.M
        r = max(x[1], 1e-10)
        th = x[2]
        f = 1.0 - 2.0 * M / r
        g = np.diag([-f, 1.0 / f, r * r, (r * np.sin(th)) ** 2])
        return g


class KerrMetric(MetricBase):
    """Kerr BH: g_{muu} in Boyer-Lindquist coordinates. M = mass, a = spin."""

    def __init__(self, mass: float = 1.0, spin: float = 0.0):
        self.M = mass
        self.a = min(spin, mass * 0.999)

    def g(self, x: np.ndarray) -> np.ndarray:
        M, a = self.M, self.a
        r = max(x[1], 1e-10)
        th = x[2]
        s, c = np.sin(th), np.cos(th)
        sigma = r * r + (a * c) ** 2
        delta = r * r - 2.0 * M * r + a * a
        g = np.zeros((4, 4))
        g[0, 0] = -(1.0 - 2.0 * M * r / sigma)
        g[0, 3] = -(2.0 * M * a * r * s * s / sigma)
        g[3, 0] = g[0, 3]
        g[1, 1] = sigma / delta
        g[2, 2] = sigma
        g[3, 3] = ((r * r + a * a + 2.0 * M * a * a * r * s * s / sigma) * s * s)
        return g


class TeoMetric(MetricBase):
    """
    Teo rotating wormhole (Teo 1998, PRD 58, 024014).

    Metric: ds^2 = -N^2 dt^2 + (1-b0/r)^{-1} dr^2
            + r^2 (dth^2 + sin^2(th) (dph - omega dt)^2)
    N(r) = exp(-b0/r), omega(r) = 2J/r^3
    ADM mass: M = b0, spin: a = J/M^2
    """

    def __init__(self, throat_radius: float = 2.0, spin: float = 0.0):
        self.throat_radius = throat_radius
        self.M = throat_radius
        self.a = spin
        self.J = spin * throat_radius

    def g(self, x: np.ndarray) -> np.ndarray:
        b0 = self.throat_radius
        J = self.J
        r = max(x[1], 1.01 * b0)
        th = x[2]
        s = np.sin(th)

        N_val = np.exp(-b0 / r)
        omega_val = 2.0 * J / r ** 3
        A_val = 1.0 / np.sqrt(max(1.0 - b0 / r, 1e-10))
        C_val = r * s

        g = np.zeros((4, 4))
        g[0, 0] = -N_val ** 2 + C_val ** 2 * omega_val ** 2
        g[0, 3] = -C_val ** 2 * omega_val
        g[3, 0] = g[0, 3]
        g[1, 1] = A_val ** 2
        g[2, 2] = r ** 2
        g[3, 3] = C_val ** 2
        return g
