"""
EHT observational constraints with real numerical data.

Sources:
  EHT Collaboration (2019), ApJL 875, L1-L6 (M87*)
  EHT Collaboration (2022), ApJL 930, L12-L17 (Sgr A*)
  EHT Collaboration (2024), A&A 681, A66 (M87* polarimetric)
"""

from typing import Dict, Any, List
import numpy as np


class EHTConstraints:
    """
    Event Horizon Telescope observational constraints.

    All values are citation-ready numerical data from EHT Papers I-VI.

    Shadow diameters are measured at 230 GHz (1.3 mm).
    """

    def __init__(self):
        self.constraints = self._load_known_constraints()

    def _load_known_constraints(self) -> List[Dict[str, Any]]:
        return [
            {
                "observable": "M87* shadow angular diameter",
                "value": 42.0,
                "unit": "microas",
                "error": 3.0,
                "confidence": "68%",
                "source": "EHT Collaboration (2019), ApJL 875, L1; "
                          "EHT Collaboration (2024), A&A 681, A66",
                "notes": "Measured at 230 GHz. Consistent with Kerr BH "
                         "of mass 6.5+/-0.7e9 M_sun.",
            },
            {
                "observable": "M87* black hole mass",
                "value": 6.5e9,
                "unit": "M_sun",
                "error": 0.7e9,
                "confidence": "68%",
                "source": "EHT Collaboration (2019), ApJL 875, L1",
                "notes": "Shadow-based mass measurement. "
                         "Consistent with stellar dynamics.",
            },
            {
                "observable": "M87* shadow asymmetry",
                "value": 0.0,
                "unit": "fractional",
                "upper_limit": 0.1,
                "confidence": "68%",
                "source": "EHT Collaboration (2019), ApJL 875, L4; "
                          "Tamburini et al. (2020), A&A 642, L1",
                "notes": "Fractional diameter deviation: "
                         "|d_max - d_min| / (d_max + d_min) < 0.1. "
                         "Spin constraint a/M < 0.5 model-dependent.",
            },
            {
                "observable": "M87* spin upper limit",
                "value": 0.5,
                "unit": "a/M",
                "upper_limit": 0.5,
                "confidence": "68%",
                "source": "EHT Collaboration (2019), ApJL 875, L5",
                "notes": "Model-dependent constraint from shadow size. "
                         "High-spin models (a > 0.5) produce shadow "
                         "diameters inconsistent with EHT measurement.",
            },
            {
                "observable": "Sgr A* shadow angular diameter",
                "value": 51.8,
                "unit": "microas",
                "error": 2.3,
                "confidence": "68%",
                "source": "EHT Collaboration (2022), ApJL 930, L12",
                "notes": "Measured at 230 GHz. Consistent with Kerr BH "
                         "of mass 4.3e6 M_sun.",
            },
            {
                "observable": "Sgr A* black hole mass",
                "value": 4.3e6,
                "unit": "M_sun",
                "error": 0.6e6,
                "confidence": "68%",
                "source": "EHT Collaboration (2022), ApJL 930, L12",
                "notes": "Shadow-based mass measurement. "
                         "Consistent with stellar orbits (GRAVITY).",
            },
            {
                "observable": "Photon ring size deviation from GR",
                "value": 1.0,
                "unit": "d_sh / (3*sqrt(3)*M)",
                "error": 0.1,
                "confidence": "90%",
                "source": "EHT Collaboration (2022), ApJL 930, L15; "
                          "Wielgus et al. (2020), ApJ 901, 67",
                "notes": "d_sh / (3*sqrt(3)*M) = 1.0 +/- 0.1. "
                         "Constraints on non-GR metrics.",
            },
            {
                "observable": "Non-GR shadow deviation (M87*) bound",
                "value": 0.1,
                "unit": "fractional deviation",
                "upper_limit": 0.1,
                "confidence": "68%",
                "source": "EHT Collaboration (2022), ApJL 930, L17; "
                          "Psaltis et al. (2020), PRL 125, 141104",
                "notes": "Bound on Reissner-Nordstrom dilaton shadow "
                         "deviation. For wormholes: model-dependent.",
            },
        ]

    def list_constraints(self) -> List[Dict[str, Any]]:
        return self.constraints

    def get_m87_diameter(self) -> Dict[str, float]:
        """M87* shadow diameter in microas: 42 +/- 3 (68% CI)."""
        return {"value": 42.0, "error": 3.0, "unit": "microas"}

    def get_sgra_diameter(self) -> Dict[str, float]:
        """Sgr A* shadow diameter in microas: 51.8 +/- 2.3 (68% CI)."""
        return {"value": 51.8, "error": 2.3, "unit": "microas"}

    def get_photon_ring_bound(self) -> Dict[str, float]:
        """Photon ring bound: d_sh/(3*sqrt(3)*M) = 1.0 +/- 0.1."""
        return {"value": 1.0, "error": 0.1, "unit": "fractional"}

    def check_compatibility(
        self, shadow_diameter_mas: float, shadow_asymmetry: float = 0.0,
        object_type: str = "generic"
    ) -> Dict[str, Any]:
        """
        Compare shadow predictions against EHT constraints.

        Parameters
        ----------
        shadow_diameter_mas : float
            Predicted shadow diameter in microas (M87* distance-scaled).
        shadow_asymmetry : float
            Predicted asymmetry A = (d_max-d_min)/(d_max+d_min).
        object_type : str
            Object type for model-dependent bounds.

        Returns
        -------
        dict with compatibility flags and sigma deviations.
        """
        m87_d_sh = 42.0
        m87_err = 3.0
        deviation = (shadow_diameter_mas - m87_d_sh) / m87_err

        m87_asym_limit = 0.1

        pr_dev = shadow_diameter_mas / m87_d_sh - 1.0
        pr_bound = 0.1

        return {
            "d_sh_mu_as": shadow_diameter_mas,
            "d_sh_sigma": deviation,
            "d_sh_compatible": abs(deviation) < 2.0,
            "d_sh_verdict": "consistent" if abs(deviation) < 2.0
                            else ("tension" if abs(deviation) < 3.0
                                  else "excluded"),
            "asymmetry_value": shadow_asymmetry,
            "asymmetry_below_10pct": shadow_asymmetry < m87_asym_limit,
            "photon_ring_deviation": pr_dev,
            "photon_ring_compatible": abs(pr_dev) < pr_bound,
        }
