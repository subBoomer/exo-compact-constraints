#!/usr/bin/env python3
"""
Gravitational wave waveform comparison.

Computes ringdown frequencies and damping times for a chosen metric.
Generates ringdown waveform comparison plots.

Usage:
    python waveform_comparison.py --metric kerr --mass 10 --spin 0.7

Output:
    visualizations/ringdown_comparison.png
"""

import argparse
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

from exo_compact.metrics.schwarzschild import Schwarzschild
from exo_compact.metrics.kerr import Kerr
from exo_compact.observables.gravitational_waves import GravitationalWaveCalculator
from exo_compact.visualization.plots import ComparisonPlotter


def parse_args():
    parser = argparse.ArgumentParser(description="GW waveform comparison")
    parser.add_argument("--metric", choices=["schwarzschild", "kerr"],
                        default="kerr")
    parser.add_argument("--mass", type=float, default=10.0)
    parser.add_argument("--spin", type=float, default=0.0)
    parser.add_argument("--l", type=int, default=2)
    parser.add_argument("--m", type=int, default=2)
    parser.add_argument("--output", type=str,
                        default="visualizations/ringdown_comparison.png")
    return parser.parse_args()


def main():
    args = parse_args()

    if args.metric == "schwarzschild":
        metric = Schwarzschild(mass=args.mass)
    else:
        metric = Kerr(mass=args.mass, spin=args.spin)

    calc = GravitationalWaveCalculator(metric)
    freq = calc.ringdown_frequency(l=args.l, m=args.m)
    tau = calc.ringdown_damping_time(l=args.l, m=args.m)

    print(f"Metric: {metric.name}")
    print(f"  M = {args.mass} M_sun, a = {args.spin}")
    print(f"  Ringdown frequency f_RD = {freq:.2f} Hz (l={args.l}, m={args.m})")
    print(f"  Damping time (tau_RD) = {tau:.4f} s")

    t = np.linspace(0, tau * 3, 1000)
    h = calc.ringdown_waveform(t, l=args.l, m=args.m)

    name = f"{metric.name} (a={args.spin})"
    schwarz_metric = Schwarzschild(mass=args.mass)
    schwarz_calc = GravitationalWaveCalculator(schwarz_metric)
    h_schwarz = schwarz_calc.ringdown_waveform(t, l=args.l, m=args.m)

    waveforms = {
        name: {"t": t, "h": h},
        "Schwarzschild": {"t": t, "h": h_schwarz},
    }

    plotter = ComparisonPlotter()
    plotter.plot_ringdown_waveform(waveforms, save_path=args.output)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
