"""
ShadowScanner: 2D image-plane scans for shadow boundary extraction.

Supports NxN pixel scans with progress reporting.
Returns binary shadow maps and boundary contours.
"""

import numpy as np
from typing import Tuple, Optional
from .engine import RayTracer


class ShadowScanner:
    """
    Image-plane scanner for shadow boundary extraction.

    Scans (alpha, beta) grid and classifies each pixel.
    Outputs binary maps and boundary contours.
    """

    def __init__(self, tracer: RayTracer):
        self.tracer = tracer

    def scan(self, alpha_range: Tuple[float, float],
             beta_range: Tuple[float, float],
             n_alpha: int = 50, n_beta: int = 50,
             lambda_max: float = 1000.0,
             verbose: bool = True) -> dict:
        """
        Full 2D image-plane scan.

        Parameters
        ----------
        alpha_range : (min, max) — horizontal impact parameter
        beta_range : (min, max) — vertical impact parameter
        n_alpha, n_beta : grid resolution
        lambda_max : max affine parameter
        verbose : print progress

        Returns
        -------
        dict with keys:
            alpha_grid, beta_grid — coordinate arrays
            classified — (n_alpha, n_beta) int array (0=escape, 1=captured, -1=unknown)
            boundary — (N, 2) array of boundary points
        """
        alphas = np.linspace(alpha_range[0], alpha_range[1], n_alpha)
        betas = np.linspace(beta_range[0], beta_range[1], n_beta)
        classified = np.zeros((n_alpha, n_beta), dtype=int)

        for i, a in enumerate(alphas):
            for j, b in enumerate(betas):
                result = self.tracer.trace_ray(a, b, lambda_max=lambda_max)
                s = result['status']
                if s in ('escape',):
                    classified[i, j] = 0
                elif s in ('captured', 'throat_cross'):
                    classified[i, j] = 1
                else:
                    classified[i, j] = -1

            if verbose and i % 10 == 0:
                n_done = (i + 1) * n_beta
                n_total = n_alpha * n_beta
                print(f'    row {i + 1}/{n_alpha} ({n_done}/{n_total} pixels)')

        boundary = []
        for i in range(n_alpha - 1):
            for j in range(n_beta - 1):
                if classified[i, j] != classified[i + 1, j]:
                    a = (alphas[i] + alphas[i + 1]) / 2
                    boundary.append([a, betas[j]])
                if classified[i, j] != classified[i, j + 1]:
                    b = (betas[j] + betas[j + 1]) / 2
                    boundary.append([alphas[i], b])

        return {
            'alpha_grid': alphas,
            'beta_grid': betas,
            'classified': classified,
            'boundary': np.array(boundary) if boundary else np.zeros((0, 2)),
        }

    def equatorial_scan(self, alpha_range: Tuple[float, float],
                        n_alpha: int = 30, lambda_max: float = 1000.0,
                        verbose: bool = True) -> dict:
        """Equatorial scan (beta=0) for fast boundary estimation."""
        return self.scan(alpha_range, (0.0, 0.1), n_alpha, 2,
                         lambda_max, verbose)

    def shadow_radius(self, alpha_range: Tuple[float, float] = (1.0, 10.0),
                      n_alpha: int = 40, lambda_max: float = 1000.0,
                      verbose: bool = True) -> float:
        """Equatorial shadow radius estimate."""
        scan = self.equatorial_scan(alpha_range, n_alpha, lambda_max, verbose)
        boundary = scan['boundary']
        if len(boundary) > 0:
            return float(np.mean(boundary[:, 0]))
        M = getattr(self.tracer.metric, 'M', 1.0)
        return 3.0 * np.sqrt(3.0) * M
