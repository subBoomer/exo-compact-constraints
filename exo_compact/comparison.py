"""
Cross-object comparison framework for all metric classes.
"""

import numpy as np
from typing import Dict, Any, List, Optional, Tuple
from .metrics.schwarzschild import Schwarzschild
from .metrics.kerr import Kerr
from .metrics.morris_thorne import MorrisThorne
from .metrics.rotating_wormhole import RotatingWormhole
from .metrics.boson_star import BosonStar
from .metrics.gravastar import Gravastar
from .metrics.eco import ECO
from .metrics.base import AbstractMetric, MetricType
from .observables.shadow import ShadowCalculator
from .observables.lensing import LensingCalculator
from .observables.gravitational_waves import GravitationalWaveCalculator
from .observables.orbital import OrbitalCalculator


class ComparisonSuite:
    """
    Runs all metrics through all observables and produces structured comparison.
    """

    DEFAULT_CONFIGS = {
        "Schwarzschild": {"mass": 1.0},
        "Kerr": {"mass": 1.0, "spin": 0.7},
        "MorrisThorne": {"mass": 1.0, "throat_radius": 2.0},
        "RotatingWormhole": {"mass": 2.0, "throat_radius": 2.0, "spin": 0.0},
        "BosonStar": {"mass": 1.0, "compactness": 0.12, "radius_scale": 8.0},
        "Gravastar": {"mass": 1.0, "interior_radius": 3.0, "shell_width": 0.3},
        "ECO": {"mass": 1.0, "spin": 0.7, "reflection_coefficient": 0.5,
                "quantum_correction": 0.01},
    }

    def __init__(self, configs: Optional[Dict[str, Dict]] = None):
        self.configs = configs or self.DEFAULT_CONFIGS
        self.metrics: Dict[str, AbstractMetric] = {}
        self.results: Dict[str, Dict[str, Any]] = {}
        self._build_metrics()

    def _build_metrics(self):
        factories = {
            "Schwarzschild": Schwarzschild,
            "Kerr": Kerr,
            "MorrisThorne": MorrisThorne,
            "RotatingWormhole": RotatingWormhole,
            "BosonStar": BosonStar,
            "Gravastar": Gravastar,
            "ECO": ECO,
        }
        for name, params in self.configs.items():
            if name in factories:
                try:
                    metric = factories[name](**params)
                    if name == "RotatingWormhole" and params.get("spin", 0) != 0:
                        metric.set_teo_metric()
                    self.metrics[name] = metric
                except Exception as e:
                    print(f"  Warning: {name} failed: {e}")

    def compute_all(self):
        for name, metric in self.metrics.items():
            self.results[name] = self._compute_single(metric, name)
        return self.results

    def _compute_single(self, metric, name):
        r = {"name": name, "params": dict(metric.params)}
        sc = ShadowCalculator(metric)
        lc = LensingCalculator(metric)
        gc = GravitationalWaveCalculator(metric)
        oc = OrbitalCalculator(metric)
        M = metric.params.get("mass", 1.0)

        try:
            b_c = sc.critical_impact_parameter()
            d_sh = sc.shadow_diameter()
            asym = sc.shadow_asymmetry()
            r["shadow"] = {"b_c/M": b_c / M, "d_sh/M": d_sh / M, "asymmetry": asym}
        except Exception as e:
            r["shadow"] = {"error": str(e)}

        try:
            alpha_w = lc.weak_deflection_angle(max(10.0 * M, M * 5))
            alpha_s = lc.deflection_for_metric(max(6.0 * M, M * 5 + 0.1))
            r["lensing"] = {"alpha_weak": alpha_w, "alpha_strong": alpha_s}
        except Exception as e:
            r["lensing"] = {"error": str(e)}

        try:
            f_rd = gc.ringdown_frequency(2, 2, 0)
            tau_rd = gc.ringdown_damping_time(2, 2, 0)
            qnm_dev = gc.qnm_deviation_from_kerr()
            tidal = gc.tidal_deformability()
            f_rd_10 = f_rd * M / 10.0
            tau_rd_10 = tau_rd * 10.0 / M
            r["gw"] = {
                "f_RD/M": f_rd, "tau_RD/M": tau_rd,
                "f_RD_10Msun": f_rd_10, "tau_RD_10Msun": tau_rd_10,
                "delta_omega": qnm_dev.get("delta_omega", 0),
                "delta_tau": qnm_dev.get("delta_tau", 0),
                "tidal_lambda": tidal,
            }
        except Exception as e:
            r["gw"] = {"error": str(e)}

        try:
            r_isco = oc.isco_radius()
            dphi = oc.periastron_precession(max(10.0 * M, r_isco * 1.5), 0.0)
            epi = oc.epicyclic_frequencies(max(10.0 * M, r_isco * 1.5))
            r["orbital"] = {
                "r_ISCO/M": r_isco / M, "dphi": dphi,
                "omega_phi": epi["omega_phi"],
            }
        except Exception as e:
            r["orbital"] = {"error": str(e)}

        return r

    def _val(self, d, key, default="---"):
        if "error" in d:
            return default
        v = d.get(key, default)
        if isinstance(v, float):
            if abs(v) < 1e-3 or abs(v) >= 1e6:
                return f"{v:.2e}"
            return f"{v:.4f}"
        return str(v)

    def to_latex_table(self) -> str:
        rows_data = [
            ("$d_{\\rm sh}/M$", "d_sh/M"),
            ("Asymmetry $A$", "asymmetry"),
            ("$f_{\\rm RD}$ [10M$_\\odot$] [Hz]", "f_RD_10Msun"),
            ("$\\tau_{\\rm RD}$ [10M$_\\odot$] [s]", "tau_RD_10Msun"),
            ("$\\delta\\omega/\\omega$", "delta_omega"),
            ("$\\lambda$ (tidal)", "tidal_lambda"),
            ("$r_{\\rm ISCO}/M$", "r_ISCO/M"),
            ("$\\Delta\\phi$ [rad]", "dphi"),
            ("$\\alpha(10M)$ [rad]", "alpha_weak"),
        ]
        obs_order = ["shadow", "shadow", "gw", "gw", "gw", "gw", "orbital", "orbital", "lensing"]

        lines = []
        lines.append(r"\begin{table}[t]")
        lines.append(r"\centering")
        lines.append(r"\caption{Cross-object observable comparison. "
                      r"Geometric units $G=c=1$.}")
        lines.append(r"\label{tab:comparison}")
        lines.append(r"\begin{tabular}{l|cccccccc}")
        lines.append(r"\toprule")
        lines.append(r"Observable & Schw & Kerr & MT WH & Rot WH & "
                      r"Boson & Grava & ECO \\")
        lines.append(r"\midrule")
        for (label, key), obs in zip(rows_data, obs_order):
            vals = [self._val(r.get(obs, {}), key) for r in self.results.values()]
            lines.append(f"  {label} & " + " & ".join(vals) + r" \\")
        lines.append(r"\bottomrule")
        lines.append(r"\end{tabular}")
        lines.append(r"\end{table}")
        return "\n".join(lines)

    def exclusion_chi2(self) -> Dict[str, float]:
        """
        Compute chi^2 against EHT + LIGO constraints.

        Shadow reference: uses the GR-predicted diameter for each metric
        (spin-dependent for Kerr, 10.39M for Schwarzschild). The EHT
        constraint for M87* is d_sh ≃ 10.39M with 10% tolerance.
        For Kerr at spin a, the GR reference is the analytic Kerr shadow
        diameter at that spin (equatorial observer).

        Ringdown reference: delta_omega and delta_tau from Kerr.
        """
        chi2_scores = {}
        for name, r in self.results.items():
            score = 0.0
            metric = self.metrics.get(name)

            sh = r.get("shadow", {})
            if "d_sh/M" in sh and metric is not None:
                d_sh_norm = sh["d_sh/M"]
                mt = metric.metric_type
                M = metric.params.get("mass", 1.0)
                spin = metric.params.get("spin", 0.0)

                if mt == MetricType.SCHWARZSCHILD:
                    d_sh_ref = 10.3923
                elif mt == MetricType.KERR:
                    from .observables.shadow import ShadowCalculator
                    sc = ShadowCalculator(metric)
                    d_sh_ref = sc.shadow_diameter() / M
                else:
                    from .metrics.kerr import Kerr
                    gr = Kerr(mass=M, spin=min(abs(spin), M * 0.99))
                    from .observables.shadow import ShadowCalculator
                    sc = ShadowCalculator(gr)
                    d_sh_ref = sc.shadow_diameter() / M

                score += ((d_sh_norm - d_sh_ref) / (0.1 * d_sh_ref)) ** 2

            if "asymmetry" in sh:
                score += (sh["asymmetry"] / 0.1) ** 2

            gw_d = r.get("gw", {})
            if "delta_omega" in gw_d:
                score += (gw_d["delta_omega"] / 0.1) ** 2
            if "delta_tau" in gw_d:
                score += (gw_d["delta_tau"] / 0.2) ** 2
            if "tidal_lambda" in gw_d:
                score += (gw_d["tidal_lambda"] / 10.0) ** 2

            chi2_scores[name] = score
        return chi2_scores
