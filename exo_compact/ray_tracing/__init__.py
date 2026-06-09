"""
Production-grade GR ray tracing architecture.

Design (per GPT prescription):
- MetricBase interface: all metrics implement g(x) and g_inv(x)
- MetricCache: precomputes metric on (r, theta) grid for 5-15x speedup
- HamiltonianGeodesic: single engine for all spacetimes
- TerminationOracle: early stopping for 10-50x speedup
- RayTracer: full pipeline with adaptive DOP853
- ShadowScanner: 2D image-plane scans
"""

from .base import MetricBase, MetricCache, TerminationOracle
from .engine import HamiltonianGeodesic, RayTracer
from .shadow import ShadowScanner
from .metrics import KerrMetric, TeoMetric, SchwarzschildMetric

__all__ = [
    "MetricBase", "MetricCache", "TerminationOracle",
    "HamiltonianGeodesic", "RayTracer",
    "ShadowScanner",
    "KerrMetric", "TeoMetric", "SchwarzschildMetric",
]
