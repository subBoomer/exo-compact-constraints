from typing import Dict, Any, Optional
from ..metrics.boson_star import BosonStar
from ..metrics.gravastar import Gravastar
from ..metrics.eco import ECO
from ..metrics.base import AbstractMetric


class ECOModel:
    """
    Exotic Compact Object model factory.

    Supports: boson_star, gravastar, eco (generic ECO with Kerr + corrections).
    """

    ECO_TYPES = {
        "boson_star": BosonStar,
        "gravastar": Gravastar,
        "eco": ECO,
    }

    def __init__(self, eco_type: str = "eco", **params):
        if eco_type.lower() not in self.ECO_TYPES:
            raise ValueError(
                f"Unknown ECO type: {eco_type}. "
                f"Available: {list(self.ECO_TYPES.keys())}"
            )
        self.eco_type = eco_type.lower()
        self.params = params
        self.metric = self._build_metric()

    def _build_metric(self) -> AbstractMetric:
        cls = self.ECO_TYPES[self.eco_type]
        return cls(**self.params)

    def __repr__(self) -> str:
        return f"ECOModel({self.eco_type})"
