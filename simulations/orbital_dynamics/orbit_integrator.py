#!/usr/bin/env python3
"""
Orbital dynamics integrator.

Computes orbital precession and ISCO for a chosen metric.
Generates precession comparison plots.

Usage:
    python orbit_integrator.py --metric schwarzschild --mass 10 --radius 20

Output:
    visualizations/orbital_precession_comparison.png
"""

import argparse
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

from exo_compact.metrics.schwarzschild import Schwarzschild
from exo_compact.metrics.kerr import Kerr
from exo_compact.observables.orbital import OrbitalCalculator
from exo_compact.visualization.plots import ComparisonPlotter


def parse_args():
    parser = argparse.ArgumentParser(description="Orbital dynamics integrator")
    parser.add_argument("--metric", choices=["schwarzschild", "kerr"],
                        default="schwarzschild")
    parser.add_argument("--mass", type=float, default=10.0)
    parser.add_argument("--spin", type=float, default=0.0)
    parser.add_argument("--r_min", type=float, default=7.0)
    parser.add_argument("--r_max", type=float, default=100.0)
    parser.add_argument("--n_points", type=int, default=50)
    parser.add_argument("--output", type=str,
                        default="visualizations/orbital_precession_comparison.png")
    return parser.parse_args()


def main():
    args = parse_args()

    if args.metric == "schwarzschild":
        metric = Schwarzschild(mass=args.mass)
    else:
        metric = Kerr(mass=args.mass, spin=args.spin)

    calc = OrbitalCalculator(metric)
    r_isco = calc.isco_radius()

    print(f"Metric: {metric.name}")
    print(f"  M = {args.mass} M_sun, a = {args.spin}")
    print(f"  ISCO radius r_ISCO = {r_isco:.4f} M")

    r_range = np.linspace(max(args.r_min, r_isco * 1.1), args.r_max, args.n_points)
    comparison = calc.compare_precession(r_range)

    plotter = ComparisonPlotter()
    plotter.metric_name = f"{metric.name} (a={args.spin})"
    plotter.plot_orbital_precession(comparison, save_path=args.output)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
