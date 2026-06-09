#!/usr/bin/env python3
"""
Generate the unified exclusion diagram (GPT-recommended 3-panel figure).

Panel A: Shadow morphology — all 7 objects with ray-traced boundaries
Panel B: Exclusion bar chart — chi^2 against EHT + LIGO constraints
Panel C: Ringdown consistency — delta f/f vs f_RD

Usage:
    python simulations/generate_paper_figure.py

Output:
    visualizations/paper_figure_unified.png (300 DPI, 18x6 inches)
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
from exo_compact.ray_tracing import TeoMetric, RayTracer, ShadowScanner
from exo_compact.metrics import (
    MorrisThorne, RotatingWormhole, BosonStar, Gravastar, ECO
)
from exo_compact.observables.shadow import ShadowCalculator
from exo_compact.observables.gravitational_waves import GravitationalWaveCalculator
from exo_compact.comparison import ComparisonSuite
from exo_compact.visualization.plots import ComparisonPlotter


def get_shadow_boundary(metric, name):
    """Get shadow boundary. Uses ray tracing only for Teo; analytic for others."""
    if name == "Teo wormhole":
        m = TeoMetric(throat_radius=2.0, spin=0.0)
        t = RayTracer(m, use_cache=False)
        s = ShadowScanner(t)
        scan = s.equatorial_scan((4.0, 7.0), n_alpha=15)
        bounds = scan["boundary"]
        if len(bounds) > 0:
            bc = np.mean(bounds[:, 0])
            ang = np.linspace(0, 2*np.pi, 200)
            return np.column_stack([bc*np.cos(ang), bc*np.sin(ang)])

    sc = ShadowCalculator(metric)
    b_c = sc.critical_impact_parameter()
    M = metric.params.get("mass", 1.0)
    ang = np.linspace(0, 2*np.pi, 200)
    return np.column_stack([b_c/M*np.cos(ang), b_c/M*np.sin(ang)])


def main():
    print("=" * 60)
    print("UNIFIED EXCLUSION DIAGRAM")
    print("=" * 60)

    suite = ComparisonSuite()
    results = suite.compute_all()
    chi2 = suite.exclusion_chi2()

    print("\nComputing shadow boundaries...")
    shadow_bounds = {}
    shadow_order = ["Schwarzschild", "Kerr", "MorrisThorne",
                    "Teo wormhole", "BosonStar", "Gravastar", "ECO"]
    for name in shadow_order:
        if name == "Teo wormhole":
            m = TeoMetric(throat_radius=2.0, spin=0.0)
        elif name in suite.metrics:
            m = suite.metrics[name]
        else:
            continue
        print(f"  {name}...", end=" ")
        b = get_shadow_boundary(m, name)
        if len(b) > 2:
            shadow_bounds[name] = b
            print(f"done ({len(b)} pts)")
        else:
            print("skipped")

    print("\nComputing ringdown data...")
    ringdown_data = {}
    for name, metric in suite.metrics.items():
        gc = GravitationalWaveCalculator(metric)
        f_rd = gc.ringdown_frequency(2, 2, 0) * metric.params.get("mass", 1.0) / 10.0
        qnm = gc.qnm_deviation_from_kerr()
        kerr_metric = suite.metrics.get("Kerr")
        if kerr_metric:
            gc_kerr = GravitationalWaveCalculator(kerr_metric)
            f_kerr = gc_kerr.ringdown_frequency(2, 2, 0) * kerr_metric.params.get("mass", 1.0) / 10.0
            delta_omega = (f_rd - f_kerr) / f_kerr
        else:
            delta_omega = qnm.get("delta_omega", 0)
        ringdown_data[name] = {"f_RD": f_rd, "delta_omega": delta_omega}

    print("\nGenerating figure...")
    plotter = ComparisonPlotter()
    plotter.plot_unified_exclusion(
        shadow_bounds, chi2, ringdown_data,
        save_path="visualizations/paper_figure_unified.png"
    )
    print("\nSaved: visualizations/paper_figure_unified.png")
    print("=" * 60)


if __name__ == "__main__":
    main()
