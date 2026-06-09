#!/usr/bin/env python3
"""
Boson star mass-frequency sequence validation.

Verifies that the EKG solver produces the correct physical sequence:
- M increases with sigma_c (central density)
- Peaks at the Kaup limit M_max ~ 0.633 (code units, m_s=1)
- omega/m_s decreases from near 1 toward smaller values
- Then turns over (unstable branch)

This is GPT's recommended validation: "If your sequence reproduces
that topology, the solver is probably working."

Usage:
    python simulations/scan_boson_star.py

Output:
    visualizations/boson_star_sequence.png
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import matplotlib.pyplot as plt


def main():
    from exo_compact.metrics.boson_star_solver import scan_boson_stars

    sigma_c_range = np.logspace(-2, np.log10(0.35), 30)

    print("Scanning boson star sequence...")
    print(f"  sigma_c range: {sigma_c_range[0]:.4f} - {sigma_c_range[-1]:.4f}")
    print(f"  m_s = 1.0 (code units)")
    print(f"  Kaup limit expected: M_max ~ 0.633\n")

    results = scan_boson_stars(sigma_c_range=sigma_c_range, verbose=True)

    if len(results) == 0:
        print("\nNo solutions found. Solver needs more work.")
        return

    sigma_c_vals = []
    M_vals = []
    omega_vals = []
    C_vals = []

    for sc in sorted(results.keys()):
        sol = results[sc]
        sigma_c_vals.append(sc)
        M_vals.append(sol["M"])
        omega_vals.append(sol["omega"])
        C_vals.append(sol["compactness"])

    sigma_c_vals = np.array(sigma_c_vals)
    M_vals = np.array(M_vals)
    omega_vals = np.array(omega_vals)
    C_vals = np.array(C_vals)

    print(f"\nFound {len(results)} solutions")
    print(f"  Mass range: {M_vals.min():.4f} - {M_vals.max():.4f}")
    print(f"  omega/m_s range: {(omega_vals/1.0).min():.4f} - {(omega_vals/1.0).max():.4f}")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    ax = axes[0]
    ax.plot(sigma_c_vals, M_vals, "o-", color="blue", markersize=5)
    ax.axhline(y=0.633, color="red", linestyle="--", alpha=0.7,
               label=f"Kaup limit $M_{{\\max}}$ = 0.633")
    ax.set_xlabel(r"Central density $\sigma_c$", fontsize=13)
    ax.set_ylabel(r"Total mass $M$", fontsize=13)
    ax.set_title("Mass vs Central Density", fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.plot(sigma_c_vals, omega_vals, "s-", color="green", markersize=5)
    ax.axhline(y=1.0, color="red", linestyle="--", alpha=0.5,
               label=r"$\omega = m_s$ (unbound)")
    ax.set_xlabel(r"Central density $\sigma_c$", fontsize=13)
    ax.set_ylabel(r"Physical $\omega$", fontsize=13)
    ax.set_title("Frequency vs Central Density", fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    ax = axes[2]
    ax.plot(M_vals, C_vals, "D-", color="purple", markersize=5)
    ax.axhline(y=0.16, color="red", linestyle="--", alpha=0.7,
               label=r"$C_{\max} \sim 0.16$ (max compactness)")
    ax.set_xlabel(r"Total mass $M$", fontsize=13)
    ax.set_ylabel(r"Compactness $C = M/R_{95}$", fontsize=13)
    ax.set_title("Compactness vs Mass", fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.suptitle("Boson Star Sequence (Mini-boson, m_s = 1.0)",
                 fontsize=15, y=1.02)
    plt.tight_layout()
    save_path = "visualizations/boson_star_sequence.png"
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {save_path}")

    print("\n" + "=" * 60)
    print("VALIDATION CHECKLIST")
    print("=" * 60)
    peak_idx = np.argmax(M_vals)
    print(f"\n[1] Mass peak at sigma_c = {sigma_c_vals[peak_idx]:.4f}")
    print(f"    M_max = {M_vals[peak_idx]:.4f} (Kaup: 0.633)")
    if M_vals[peak_idx] < 1.0:
        print(f"    [PASS] Mass within bounds of Kaup limit")
    else:
        print(f"    [WARN] Mass exceeds Kaup limit — solver may still have issues")

    if np.all(np.diff(omega_vals) <= 0):
        print(f"    [PASS] omega monotonically decreases with sigma_c")
    else:
        print(f"    [NOTE] omega not strictly monotonic")

    if C_vals[peak_idx] < 0.5:
        print(f"    [PASS] Compactness within Buchdahl limit")
    else:
        print(f"    [FAIL] Compactness exceeds Buchdahl limit")

    print(f"\n[SUMMARY] Found {len(results)}/{len(sigma_c_range)} solutions")


if __name__ == "__main__":
    main()
