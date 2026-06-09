from typing import Dict, Any, Optional, Callable
from ..metrics.morris_thorne import MorrisThorne
from ..metrics.rotating_wormhole import RotatingWormhole
from ..metrics.base import AbstractMetric


class WormholeModel:
    """
    Wormhole model factory.

    Supports Morris-Thorne (static) and RotatingWormhole (metric-agnostic).
    """

    def __init__(self, wh_type: str = "morris_thorne", **params):
        self.wh_type = wh_type.lower()
        self.params = params
        self.metric = self._build_metric()

    def _build_metric(self) -> AbstractMetric:
        if self.wh_type == "morris_thorne":
            return MorrisThorne(**self.params)
        elif self.wh_type == "rotating":
            return RotatingWormhole(**self.params)
        else:
            raise ValueError(f"Unknown wormhole type: {self.wh_type}")

    def __repr__(self) -> str:
        return f"WormholeModel({self.wh_type}, throat={self.params.get('throat_radius', 2.0)})"
