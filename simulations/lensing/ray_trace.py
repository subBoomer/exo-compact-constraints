#!/usr/bin/env python3
"""
Lensing ray tracer simulation.

Computes deflection angles and Einstein rings for a chosen metric.
Generates comparison plots.

Usage:
    python ray_trace.py --metric schwarzschild --mass 4e6 --num_rays 100

Output:
    visualizations/lensing_deflection_angle.png
"""

import argparse
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

from exo_compact.metrics.schwarzschild import Schwarzschild
from exo_compact.metrics.kerr import Kerr
from exo_compact.metrics.morris_thorne import MorrisThorne
from exo_compact.observables.lensing import LensingCalculator
from exo_compact.visualization.plots import ComparisonPlotter


def parse_args():
    parser = argparse.ArgumentParser(description="Lensing ray tracer")
    parser.add_argument("--metric", choices=["schwarzschild", "kerr", "morris_thorne"],
                        default="schwarzschild")
    parser.add_argument("--mass", type=float, default=4e6)
    parser.add_argument("--spin", type=float, default=0.0)
    parser.add_argument("--throat", type=float, default=2.0)
    parser.add_argument("--num_rays", type=int, default=100)
    parser.add_argument("--b_min", type=float, default=1.0)
    parser.add_argument("--b_max", type=float, default=100.0)
    parser.add_argument("--output", type=str,
                        default="visualizations/lensing_deflection_angle.png")
    return parser.parse_args()


def main():
    args = parse_args()

    if args.metric == "schwarzschild":
        metric = Schwarzschild(mass=args.mass)
    elif args.metric == "kerr":
        metric = Kerr(mass=args.mass, spin=args.spin)
    elif args.metric == "morris_thorne":
        metric = MorrisThorne(mass=args.mass, throat_radius=args.throat)

    calc = LensingCalculator(metric)
    b_range = np.linspace(args.b_min, args.b_max, args.num_rays)
    comparison = calc.compare_deflection(b_range)

    plotter = ComparisonPlotter()
    plotter.plot_deflection_angle(
        comparison, metric_name=metric.name, save_path=args.output
    )
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
