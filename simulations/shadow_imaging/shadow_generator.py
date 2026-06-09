#!/usr/bin/env python3
"""
Black hole shadow generator.

Computes shadow boundary and diameter for a chosen metric.
Generates shadow comparison plots.

Usage:
    python shadow_generator.py --metric schwarzschild --mass 4e6

Output:
    visualizations/shadow_boundary_comparison.png
"""

import argparse
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

from exo_compact.metrics.schwarzschild import Schwarzschild
from exo_compact.metrics.kerr import Kerr
from exo_compact.observables.shadow import ShadowCalculator
from exo_compact.visualization.plots import ComparisonPlotter


def parse_args():
    parser = argparse.ArgumentParser(description="Shadow generator")
    parser.add_argument("--metric", choices=["schwarzschild", "kerr"],
                        default="schwarzschild")
    parser.add_argument("--mass", type=float, default=4e6)
    parser.add_argument("--spin", type=float, default=0.0)
    parser.add_argument("--inclination", type=float, default=90.0)
    parser.add_argument("--output", type=str,
                        default="visualizations/shadow_boundary_comparison.png")
    return parser.parse_args()


def main():
    args = parse_args()

    if args.metric == "schwarzschild":
        metric = Schwarzschild(mass=args.mass)
    else:
        metric = Kerr(mass=args.mass, spin=args.spin)

    calc = ShadowCalculator(metric)
    b_c = calc.critical_impact_parameter()
    d_sh = calc.shadow_diameter()

    print(f"Metric: {metric.name}")
    print(f"  M = {args.mass}")
    print(f"  Critical impact parameter b_c = {b_c:.4f} [{args.mass:.2e} M]")
    print(f"  Shadow diameter d_sh = {d_sh:.4f} [{args.mass:.2e} M]")

    boundary = calc.shadow_boundary_kerr(
        inclination=np.radians(args.inclination)
    )

    bc_ref = 3.0 * np.sqrt(3.0)
    comparison = {
        metric.name: boundary,
        "Schwarzschild": np.column_stack([
            np.linspace(-bc_ref, bc_ref, 100),
            np.sqrt(bc_ref ** 2 - np.linspace(-bc_ref, bc_ref, 100) ** 2)
        ]),
    }

    plotter = ComparisonPlotter()
    plotter.plot_shadow_boundary(comparison, save_path=args.output)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
