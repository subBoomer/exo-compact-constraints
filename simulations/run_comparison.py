#!/usr/bin/env python3
"""
Master comparison runner.

Orchestrates all metrics through all observables and generates:
1. All-object shadow overlay figure
2. All-object ringdown comparison figure
3. All-object deflection overlay figure
4. Exclusion contour (chi^2 against EHT + LIGO)
5. LaTeX comparison table

Usage:
    python simulations/run_comparison.py

Outputs:
    visualizations/comparison_shadow_overlay.png
    visualizations/comparison_ringdown.png
    visualizations/comparison_deflection.png
    visualizations/comparison_exclusion.png
    docs/comparison_table.tex
"""

import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
from exo_compact.comparison import ComparisonSuite
from exo_compact.visualization.plots import ComparisonPlotter
from exo_compact.observables.shadow import ShadowCalculator
from exo_compact.observables.lensing import LensingCalculator
from exo_compact.observables.gravitational_waves import GravitationalWaveCalculator
from exo_compact.observables.orbital import OrbitalCalculator


def main():
    print("=" * 60)
    print("EXO-COMPACT CONSTRAINTS: Cross-Object Comparison")
    print("=" * 60)

    suite = ComparisonSuite()
    print(f"\nMetrics loaded: {len(suite.metrics)}")
    for name in suite.metrics:
        print(f"  - {name}")

    print("\nComputing all observables...")
    results = suite.compute_all()
    print(f"Results computed for {len(results)} objects")

    print("\n" + "-" * 60)
    print("OBSERVABLE SUMMARY")
    print("-" * 60)
    for name, r in results.items():
        print(f"\n  {name}:")
        for obs_type in ("shadow", "gw", "orbital", "lensing"):
            if obs_type in r:
                vals = r[obs_type]
                if "error" in vals:
                    print(f"    {obs_type}: ERROR {vals['error']}")
                else:
                    items = " | ".join(f"{k}={v:.4f}" if isinstance(v, float)
                                       else f"{k}={v}" for k, v in vals.items()
                                       if isinstance(v, (int, float)))
                    print(f"    {obs_type}: {items}")

    print("\n" + "-" * 60)
    print("LATEX TABLE")
    print("-" * 60)
    latex = suite.to_latex_table()
    print(latex)

    table_path = os.path.join("docs", "comparison_table.tex")
    os.makedirs("docs", exist_ok=True)
    with open(table_path, "w") as f:
        f.write(latex)
    print(f"\nSaved: {table_path}")

    print("\n" + "-" * 60)
    print("CHI-SQUARED EXCLUSION")
    print("-" * 60)
    chi2 = suite.exclusion_chi2()
    for name, score in sorted(chi2.items(), key=lambda x: x[1]):
        label = "GOOD" if score < 2 else "TENSION" if score < 6 else "EXCLUDED"
        print(f"  {name:20s}: chi2 = {score:.2f}  [{label}]")

    print("\n" + "-" * 60)
    print("GENERATING FIGURES")
    print("-" * 60)

    plotter = ComparisonPlotter()

    # --- Shadow overlay ---
    print("  Shadow overlay...")
    shadow_boundaries = {}
    for name, metric in suite.metrics.items():
        try:
            sc = ShadowCalculator(metric)
            spin = metric.params.get("spin", 0.0)
            b0 = metric.params.get("throat_radius", None)
            if name in ("Schwarzschild",) or (name == "Kerr" and spin == 0.0):
                b = sc._schwarzschild_shadow_boundary(200)
            elif name in ("MorrisThorne", "RotatingWormhole"):
                r0 = b0 if b0 else metric.params.get("mass", 1.0) * 2
                ang = np.linspace(0, 2 * np.pi, 200)
                b = np.column_stack([r0 * np.cos(ang), r0 * np.sin(ang)])
                if name == "RotatingWormhole":
                    name += " (throat approx)"
            else:
                b = sc.shadow_boundary_kerr(200, np.pi / 2)
            if len(b) > 2:
                shadow_boundaries[name] = b
        except Exception as e:
            print(f"    {name}: shadow error {e}")

    plotter.plot_shadow_overlay(shadow_boundaries,
                                "visualizations/comparison_shadow_overlay.png")

    # --- Ringdown overlay ---
    print("  Ringdown overlay...")
    ringdown_data = {}
    for name, metric in suite.metrics.items():
        try:
            gc = GravitationalWaveCalculator(metric)
            tau = gc.ringdown_damping_time(2, 2, 0)
            t = np.linspace(0, max(tau * 3, 0.01), 1000)
            h = gc.ringdown_waveform(t, 2, 2, 0)

            name_label = name
            if name == "Kerr":
                name_label = f"Kerr (a={metric.params.get('spin', 0)})"
            elif name == "ECO":
                name_label = f"ECO (R={metric.params.get('reflection_coefficient', 0)})"
            elif name == "MorrisThorne":
                name_label = f"MT WH (b0={metric.params.get('throat_radius', 0)})"

            f_rd = gc.ringdown_frequency(2, 2, 0)
            ringdown_data[name_label] = {"t": t, "h": h, "f_RD": f_rd}
        except Exception as e:
            print(f"    {name}: ringdown error {e}")

    plotter.plot_ringdown_overlay(ringdown_data,
                                  "visualizations/comparison_ringdown.png")

    # --- Deflection overlay ---
    print("  Deflection overlay...")
    deflection_data = {}
    for name, metric in suite.metrics.items():
        try:
            lc = LensingCalculator(metric)
            M = metric.params.get("mass", 1.0)
            b_range = np.linspace(max(M * 3, 5), 50, 100)
            alpha = np.array([lc.deflection_for_metric(b) for b in b_range])
            finite = np.isfinite(alpha)
            deflection_data[name] = {
                "b": b_range[finite], "alpha": alpha[finite]
            }
        except Exception as e:
            print(f"    {name}: deflection error {e}")

    plotter.plot_deflection_overlay(deflection_data,
                                    "visualizations/comparison_deflection.png")

    # --- Exclusion contours ---
    print("  Exclusion contours...")
    plotter.plot_exclusion_contours(chi2,
                                    "visualizations/comparison_exclusion.png")

    print("\n" + "=" * 60)
    print("COMPARISON COMPLETE")
    print("=" * 60)
    print("\nGenerated files:")
    for f in ["comparison_shadow_overlay.png", "comparison_ringdown.png",
              "comparison_deflection.png", "comparison_exclusion.png"]:
        path = f"visualizations/{f}"
        if os.path.exists(path):
            print(f"  {path}")
    print(f"  docs/comparison_table.tex")


if __name__ == "__main__":
    main()
