from typing import Dict, Any, Optional
from ..metrics.schwarzschild import Schwarzschild
from ..metrics.kerr import Kerr
from ..metrics.base import AbstractMetric


class BlackHoleModel:
    """
    Black hole model factory. Instantiates the correct metric
    based on the requested BH type.
    """

    def __init__(self, bh_type: str = "schwarzschild", **params):
        self.bh_type = bh_type.lower()
        self.params = params
        self.metric = self._build_metric()

    def _build_metric(self) -> AbstractMetric:
        if self.bh_type == "schwarzschild":
            return Schwarzschild(**self.params)
        elif self.bh_type == "kerr":
            return Kerr(**self.params)
        else:
            raise ValueError(f"Unknown black hole type: {self.bh_type}")

    def __repr__(self) -> str:
        return f"BlackHoleModel({self.bh_type}, M={self.params.get('mass', 1.0)})"
