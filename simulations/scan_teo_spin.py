#!/usr/bin/env python3
"""
Teo wormhole shadow scans at multiple spin values.

Uses the production-grade Hamiltonian ray tracer to find
the shadow boundary for Teo wormholes at a=0, 0.5, 0.9.

Usage:
    python simulations/scan_teo_spin.py

Output:
    visualizations/teo_shadow_spin_comparison.png
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import matplotlib.pyplot as plt
from exo_compact.ray_tracing import TeoMetric, RayTracer, ShadowScanner


def find_shadow_radius(metric, n_alpha=20):
    """Find equatorial shadow radius via ray tracing."""
    tracer = RayTracer(metric, r_obs=100.0, use_cache=False)
    scanner = ShadowScanner(tracer)
    scan = scanner.equatorial_scan((3.0, 7.0), n_alpha=n_alpha, lambda_max=500.0)
    bounds = scan["boundary"]
    if len(bounds) > 0:
        return float(np.mean(bounds[:, 0]))
    return None


def main():
    print("=" * 60)
    print("TEO WORMHOLE SPIN SCAN")
    print("=" * 60)

    spins = [0.0, 0.5, 0.9]
    results = {}

    for a in spins:
        print(f"\nSpin a={a}...", flush=True)
        metric = TeoMetric(throat_radius=2.0, spin=a)
        bc = find_shadow_radius(metric, n_alpha=12)
        results[a] = bc
        if bc:
            print(f"  Shadow radius: {bc:.4f}M")
        else:
            print(f"  Shadow radius: not found")

    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)
    for a, bc in results.items():
        print(f"  Teo a={a:.1f}: b_c = {bc:.4f}M")
    print()

    print("Generating comparison figure...")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    ax = axes[0]
    for i, (a, bc) in enumerate(results.items()):
        ang = np.linspace(0, 2 * np.pi, 200)
        ax.plot(bc * np.cos(ang), bc * np.sin(ang), "-",
                label=f"Teo a={a:.1f}", linewidth=2, color=colors[i])
    bc_schw = 3.0 * np.sqrt(3.0)
    ax.plot(bc_schw * np.cos(ang), bc_schw * np.sin(ang), "--",
            label="Schwarzschild", linewidth=2, color="gray", alpha=0.7)
    ax.set_xlabel(r"$\alpha$ [$M$]", fontsize=13)
    ax.set_ylabel(r"$\beta$ [$M$]", fontsize=13)
    ax.set_title("Teo Wormhole Shadows (Ray-Traced)", fontsize=14)
    ax.set_aspect("equal")
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    spin_vals = list(results.keys())
    radius_vals = [results[s] for s in spin_vals]
    ax.plot(spin_vals, radius_vals, "o-", color="blue", markersize=8, linewidth=2)
    ax.axhline(y=np.e * 2, color="red", linestyle="--", alpha=0.6,
               label=f"Analytic (a=0): e*b0 = {np.e*2:.2f}M")
    ax.axhline(y=3*np.sqrt(3), color="gray", linestyle=":", alpha=0.6,
               label=f"Schwarzschild: {3*np.sqrt(3):.2f}M")
    ax.set_xlabel("Spin parameter a", fontsize=13)
    ax.set_ylabel(r"Shadow radius $b_c$ [$M$]", fontsize=13)
    ax.set_title("Shadow Radius vs Spin", fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    save_path = "visualizations/teo_shadow_spin_comparison.png"
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {save_path}")


if __name__ == "__main__":
    main()
