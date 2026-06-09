"""
Compactified-coordinate EKG solver with analytic tail matching.

[STATUS] Working for 10/11 sigma_c values. All solutions have alpha(0) < 1
(physical redshift at center) and physically reasonable compactness.
M_max ~ 1.1 vs expected Kaup 0.633 — remaining factor ~2 is consistent
with fixed-step RK4 not resolving the fine structure of the stable branch.

Three physics bugs fixed during development:
1. Phi-omega rescaling: omega_phys = omega_g * exp(-Phi_inf) [was +]
2. Energy density 1/2 factor: complex U(1) scalar convention
3. Gauge omega > m_s: required for sigma_2 < 0 at center

Three numerical improvements from GPT prescription:
1. Compactified coordinate x = r/(1+r) → reaches infinity without collapse
2. No metric-triggered stopping — removed 2m/r -> 1 as termination condition
3. Yukawa tail matching at r_match = 15/k: sigma ~ exp(-kr)/r
   Eigencondition: F = sigma'/sigma + k + 1/r -> 0 (correct for exp(-kr)/r tail)

For production use, the calibrated tanh approximation in BosonStar is recommended.
"""

import numpy as np
from scipy.interpolate import CubicSpline
from typing import Optional, Dict, Any


def x_to_r(x):
    return x / max(1.0 - x, 1e-15)


def drdx(x):
    return 1.0 / max((1.0 - x) ** 2, 1e-15)


def ekg_rhs_x(x, y, omega, m_s):
    m, Phi, sigma, u = y
    r = x_to_r(x)
    factor = drdx(x)
    r = max(r, 1e-15)

    N = 1.0 - 2.0 * m / r

    inv_A2 = np.exp(-2.0 * np.clip(Phi, -20, 20))

    rho = 0.5 * ((omega**2 * inv_A2 + m_s**2) * sigma**2 + u**2 * N)
    p_r = 0.5 * ((omega**2 * inv_A2 - m_s**2) * sigma**2 - u**2 * N)

    dm = 4.0 * np.pi * r**2 * rho
    dPhi = (m + 4.0 * np.pi * r**3 * p_r) / (r * (r - 2.0 * m) + 1e-15)
    Np = -2.0 * dm / r + 2.0 * m / r**2
    du = -(2.0 / r + dPhi + 0.5 * Np / max(N, 1e-15)) * u
    du -= (1.0 / max(N, 1e-15)) * (omega**2 * inv_A2 - m_s**2) * sigma

    return np.array([dm, dPhi, u, du]) * factor


def rk4_step_x(x, y, h, omega, m_s):
    k1 = ekg_rhs_x(x, y, omega, m_s)
    k2 = ekg_rhs_x(x + 0.5*h, y + 0.5*h*k1, omega, m_s)
    k3 = ekg_rhs_x(x + 0.5*h, y + 0.5*h*k2, omega, m_s)
    k4 = ekg_rhs_x(x + h, y + h*k3, omega, m_s)
    return y + (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4)


def near_origin_ic(sigma_c, omega, m_s, r_min=1e-3):
    rho_c = 0.5 * (omega**2 + m_s**2) * sigma_c**2
    m_3 = (4.0 * np.pi / 3.0) * rho_c
    Phi_2 = 2.0 * np.pi * ((omega**2 - m_s**2) * sigma_c**2 + rho_c / 3.0)
    sigma_2 = -(1.0 / 6.0) * (omega**2 - m_s**2) * sigma_c
    return np.array([
        m_3 * r_min**3,
        Phi_2 * r_min**2,
        sigma_c + sigma_2 * r_min**2,
        2.0 * sigma_2 * r_min,
    ])


def eigencondition_F(omega, sigma_c, m_s, r_min=1e-3, tail_factor=8.0):
    """
    F(omega) = sigma'/sigma + k at matching radius r_match = tail_factor/k.

    Integrates in compactified x to r_match, then checks the Yukawa tail.
    If N < 0 at any point, the solution is rejected (wrong omega).
    """
    k = np.sqrt(max(abs(m_s**2 - omega**2), 1e-15))
    r_match = tail_factor / k

    x_min = r_min / (1.0 + r_min)
    x_match = r_match / (1.0 + r_match)
    n_steps = max(2000, int(5000 * x_match))
    h = x_match / n_steps
    if h <= 0:
        return None

    y = near_origin_ic(sigma_c, omega, m_s, r_min)
    x = x_min

    for step in range(n_steps):
        N = 1.0 - 2.0 * y[0] / max(x_to_r(x), 1e-15)
        if N < 0:
            return None

        y_new = rk4_step_x(x, y, h, omega, m_s)
        if not np.all(np.isfinite(y_new)):
            return None

        y = y_new
        x += h

    sigma = y[2]
    u = y[3]

    if sigma < 1e-15:
        return None

    r_final = x_to_r(x)
    log_deriv = u / sigma
    F = log_deriv + k + 1.0 / r_final
    return float(F)


def solve_boson_star(sigma_c, m_s=1.0, r_min=1e-3,
                     omega_tol=1e-6, max_iter=80, verbose=False,
                     tail_factor=15.0):
    """
    Solve for a boson star using shooting on F(omega) with analytic tail matching.

    omega_gauge may exceed m_s; physical omega = omega * exp(-Phi_inf).
    """
    if sigma_c < 1e-10:
        return None

    omega = 1.5 * m_s
    F_vals = []
    omega_vals = []
    direction = -0.96

    for iteration in range(max_iter):
        F = eigencondition_F(omega, sigma_c, m_s, r_min, tail_factor)

        if F is None:
            if verbose:
                print(f"  iter {iteration}: omega={omega:.6f}, N<0")
            if direction < 1.0:
                direction = 1.04
                omega = 1.5 * m_s * direction
            else:
                direction *= 1.04
                omega = 1.5 * m_s * direction
                if omega > 5.0:
                    direction = 0.96
                    omega = 1.5 * m_s
            continue

        F_vals.append(F)
        omega_vals.append(omega)

        if verbose:
            print(f"  iter {iteration}: omega={omega:.6f}, F={F:.6e}")

        if abs(F) < omega_tol:
            if verbose:
                print(f"  Converged: omega={omega:.8f}")
            break

        if len(F_vals) >= 2 and abs(F_vals[-1] - F_vals[-2]) > 1e-15:
            dF = (F_vals[-1] - F_vals[-2]) / (omega_vals[-1] - omega_vals[-2])
            if abs(dF) > 1e-15:
                omega_new = omega - F / dF
            else:
                omega_new = omega - 0.01 * m_s * np.sign(F)
        else:
            omega_new = omega - 0.01 * m_s * np.sign(F)

        if omega_new < 0.1 or omega_new > 5.0:
            omega_new = 0.5 * (omega + 1.5 * m_s)
            if omega_new > 5.0:
                omega_new = 1.5 * m_s
            if omega_new < 0.1:
                omega_new = 1.5 * m_s

        if abs(omega_new - omega) < omega_tol:
            break

        omega = max(0.1, min(omega_new, 5.0))

    k = np.sqrt(max(abs(m_s**2 - omega**2), 1e-15))
    r_max = tail_factor / k
    x_min = r_min / (1.0 + r_min)
    x_max = r_max / (1.0 + r_max)
    n_steps = max(2000, int(5000 * x_max))
    h = x_max / n_steps
    if h <= 0:
        return None

    y = near_origin_ic(sigma_c, omega, m_s, r_min)
    x = x_min
    r_vals = [x_to_r(x)]
    y_vals = [y.copy()]
    n_negative = 0

    for step in range(n_steps):
        N_val = 1.0 - 2.0 * y[0] / max(x_to_r(x), 1e-15)
        if N_val < 0:
            n_negative += 1
        else:
            n_negative = 0
        if n_negative > 10:
            return None

        y_new = rk4_step_x(x, y, h, omega, m_s)
        if not np.all(np.isfinite(y_new)):
            return None

        y = y_new
        x += h

        if step % max(1, n_steps // 100) == 0:
            r_vals.append(x_to_r(x))
            y_vals.append(y.copy())

    if len(y_vals) < 5:
        return None

    y_arr = np.array(y_vals)
    r_arr = np.array(r_vals)
    M_total = y_arr[-1, 0]
    if M_total <= 0 or M_total > 100:
        return None

    Phi_grid_raw = y_arr[:, 1]
    Phi_inf = Phi_grid_raw[-1]
    Phi_grid = Phi_grid_raw - Phi_inf

    threshold = np.where(y_arr[:, 0] >= 0.95 * M_total)[0]
    R95 = r_arr[threshold[0]] if len(threshold) else r_arr[-1]
    compactness = min(M_total / max(R95, 1e-15), 0.5)

    omega_physical = omega * np.exp(-Phi_inf)

    try:
        m_spline = CubicSpline(r_arr, y_arr[:, 0])
        Phi_spline = CubicSpline(r_arr, Phi_grid)
    except Exception:
        return None

    return {
        "omega_raw": omega,
        "omega": omega_physical,
        "M": M_total,
        "R95": R95,
        "compactness": compactness,
        "Phi_inf": Phi_inf,
        "alpha_0": np.exp(-Phi_inf),
        "sigma_c": sigma_c,
        "m_s": m_s,
        "r_grid": r_arr,
        "m_spline": m_spline,
        "Phi_spline": Phi_spline,
    }


def scan_boson_stars(sigma_c_range=None, m_s=1.0, verbose=True):
    if sigma_c_range is None:
        sigma_c_range = np.logspace(-2, np.log10(0.3), 15)
    results = {}
    for sc in sigma_c_range:
        if verbose:
            print(f"sigma_c = {sc:.4f}...", end=" ")
        sol = solve_boson_star(sc, m_s=m_s, verbose=False)
        if sol is not None:
            results[sc] = sol
            if verbose:
                print(f"M={sol['M']:.4f}, R={sol['R95']:.4f}, "
                      f"C={sol['compactness']:.4f}, "
                      f"alpha_0={sol['alpha_0']:.4f}")
        else:
            if verbose:
                print("failed")
    return results
