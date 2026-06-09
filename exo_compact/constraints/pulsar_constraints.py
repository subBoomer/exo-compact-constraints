"""
Pulsar timing constraints on compact objects with real numerical data.

Sources:
  PSR B1913+16: Weisberg & Huang (2016), ApJ 829, 55
  Double pulsar: Kramer et al. (2021), PRD 103, 102002
  NANOGrav 15yr: Agazie et al. (2023), ApJL 951, L1
  EPTA DR2: Antoniadis et al. (2023), A&A 678, A48
"""

from typing import Dict, Any, List


class PulsarConstraints:
    """
    Pulsar timing constraints on compact objects and gravity theories.

    Tracks:
    - Binary pulsar orbital decay (strong-field GR tests)
    - Double pulsar precession (frame-dragging tests)
    - NS maximum mass (compactness bounds)
    - PTA gravitational wave background
    """

    def __init__(self):
        self.constraints = self._load_known_constraints()

    def _load_known_constraints(self) -> List[Dict[str, Any]]:
        return [
            {
                "observable": "Orbital decay PSR B1913+16",
                "value": -2.46e-12,
                "unit": "s^-1 (period derivative / period)",
                "error": 0.06e-12,
                "gr_agreement": 0.001,
                "confidence": "99.99%",
                "source": "Weisberg & Huang (2016), ApJ 829, 55",
                "notes": "First indirect GW detection. Orbital period decay "
                         "agrees with GR quadrupole formula to 0.1%. "
                         "Constrains dipole radiation from alternative gravity.",
            },
            {
                "observable": "Orbital precession PSR J0737-3039A/B",
                "value": 16.90,
                "unit": "deg/yr",
                "error": 0.01,
                "confidence": "99.99%",
                "source": "Kramer et al. (2021), PRD 103, 102002",
                "notes": "Double pulsar system. Precession rate consistent "
                         "with GR. Constrains strong-field frame-dragging.",
            },
            {
                "observable": "Neutron star maximum mass",
                "value": 2.35,
                "unit": "M_sun",
                "error": 0.17,
                "confidence": "68%",
                "source": "Romani et al. (2022), ApJL 935, L17 (PSR J0952-0607); "
                          "Fonseca et al. (2021), ApJL 915, L12 (PSR J0740+6620)",
                "notes": "Highest precision NS mass measurement. "
                         "M = 2.35 +/- 0.17 Msun for the black widow pulsar. "
                         "PSR J0740+6620: M = 2.08 +/- 0.07 Msun.",
            },
            {
                "observable": "Gravitational wave background amplitude",
                "value": 2.4e-15,
                "unit": "characteristic strain at f=1/yr",
                "range": "2.0-4.0e-15",
                "confidence": "~50%",
                "source": "NANOGrav (2023), ApJL 951, L1; "
                          "EPTA (2023), A&A 678, A48",
                "notes": "Candidate GWB detection from PTA data. "
                         "Preferred amplitude h_c = 2.4e-15 at f=1/yr. "
                         "Constrains supermassive BH binary merger rate.",
            },
            {
                "observable": "Strong equivalence principle bound",
                "value": 1e-3,
                "unit": "Nordtvedt parameter eta",
                "upper_limit": 1e-3,
                "confidence": "95%",
                "source": "Lunar Laser Ranging; "
                          "Hofmann & Muller (2018), CQG 35, 035015",
                "notes": "eta < 10^-3 from LLR. "
                         "Pulsar timing gives comparable bound.",
            },
        ]

    def list_constraints(self) -> List[Dict[str, Any]]:
        return self.constraints

    def get_max_ns_mass(self) -> Dict[str, float]:
        """Maximum neutron star mass from pulsar timing."""
        return {"value": 2.35, "error": 0.17, "unit": "M_sun"}

    def check_ns_mass_compatibility(self, mass_model: float) -> Dict[str, Any]:
        """Check if object mass exceeds maximum NS mass."""
        max_mass = 2.35
        mass_err = 0.17
        deviation = (mass_model - max_mass) / mass_err
        return {
            "model_mass_Msun": mass_model,
            "max_ns_mass_Msun": max_mass,
            "exceeds_max_NS": mass_model > (max_mass + 2 * mass_err),
            "sigma_above_NS_max": deviation,
        }
