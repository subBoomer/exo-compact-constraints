import numpy as np
from typing import Dict, Any, Optional, Tuple
from ..metrics.base import AbstractMetric, MetricType


class GravitationalWaveCalculator:
    """
    Gravitational wave observables from spacetime metrics.

    Computes:
    - Quasinormal mode (QNM) frequencies and damping times (Kerr Berti 2006 fits)
    - Ringdown waveform (damped sinusoid)
    - ECO ringdown with surface reflectivity (echo model, Mark et al. 2017)
    - Wormhole QNMs (Poschl-Teller potential approximation)
    - Tidal deformability (Love numbers)

    Sources:
      Kerr QNM: Berti, Cardoso & Will (2006), PRD 73, 064030
      ECO echoes: Mark et al. (2017), PRD 96, 084002
      Wormhole QNM: Konoplya & Zhidenko (2011), PRD 84, 124047
      Tidal: Damour & Nagar (2009), PRD 80, 084035
    """

    def __init__(self, metric: AbstractMetric):
        self.metric = metric

    def _kerr_qnm_fit(self, l: int = 2, m: int = 2, n: int = 0) -> Tuple[float, float]:
        """
        Kerr QNM frequencies using Berti et al. (2006) analytic fit.

        omega = omega_1 + omega_2 * (1 - a)^omega_3
        tau = tau_1 * exp(tau_2 * (1 - a)^tau_3)
        QNM frequency: f = omega + i/tau

        Source: Berti, Cardoso, & Will (2006), PRD 73, 064030.
        """
        M = self.metric.params.get("mass", 1.0)
        spin = self.metric.params.get("spin", 0.0)

        a = max(0.0, min(spin / M, 0.99))

        fit_params = {
            (2, 2, 0): {
                "omega": (0.4434, 0.3236, 0.3105),
                "tau": (10.853, 0.7626, 0.1775),
            },
            (2, 1, 0): {
                "omega": (0.4122, 0.3396, 0.3123),
                "tau": (10.999, 0.7762, 0.1775),
            },
            (3, 3, 0): {
                "omega": (0.6349, 0.4449, 0.3622),
                "tau": (7.219, 0.7341, 0.1818),
            },
        }

        key = (l, m, n)
        if key not in fit_params:
            key = (2, 2, 0)

        fp = fit_params[key]
        omega_fit = fp["omega"]
        tau_fit = fp["tau"]

        omega_M = omega_fit[0] + omega_fit[1] * (1.0 - a) ** omega_fit[2]
        tau_M = tau_fit[0] * np.exp(tau_fit[1] * (1.0 - a) ** tau_fit[2])

        omega = omega_M / M
        tau = tau_M * M

        return omega, tau

    def ringdown_frequency(self, l: int = 2, m: int = 2, n: int = 0) -> float:
        """QNM oscillation frequency f_RD in Hz for given mass and spin."""
        M = self.metric.params.get("mass", 1.0)
        M_sun_kg = 1.98847e30
        G = 6.67430e-11
        c = 299792458.0

        mass_kg = M * M_sun_kg
        M_in_meters = G * mass_kg / c ** 2
        time_in_seconds = M_in_meters / c

        omega, _ = self._kerr_qnm_fit(l, m, n)
        f_RD = omega / time_in_seconds / (2.0 * np.pi)

        return float(f_RD)

    def ringdown_damping_time(self, l: int = 2, m: int = 2, n: int = 0) -> float:
        """Ringdown damping time tau_RD in seconds."""
        M = self.metric.params.get("mass", 1.0)
        M_sun_kg = 1.98847e30
        G = 6.67430e-11
        c = 299792458.0

        mass_kg = M * M_sun_kg
        M_in_meters = G * mass_kg / c ** 2
        time_in_seconds = M_in_meters / c

        _, tau = self._kerr_qnm_fit(l, m, n)
        tau_RD = tau * time_in_seconds

        return float(tau_RD)

    def ringdown_waveform(
        self,
        t: np.ndarray,
        l: int = 2,
        m: int = 2,
        n: int = 0,
        amplitude: float = 1.0,
        phase: float = 0.0,
    ) -> np.ndarray:
        """
        Ringdown waveform: h(t) = A * exp(-t/tau) * cos(omega * t + phi_0)

        Source: Standard damped sinusoid model, Abbott et al. (2016), PRL 116.
        """
        omega, tau = self._kerr_qnm_fit(l, m, n)
        M = self.metric.params.get("mass", 1.0)
        M_sun_kg = 1.98847e30
        G = 6.67430e-11
        c = 299792458.0
        mass_kg = M * M_sun_kg
        M_in_meters = G * mass_kg / c ** 2
        time_in_seconds = M_in_meters / c

        tau_s = tau * time_in_seconds
        omega_Hz = omega / time_in_seconds

        return amplitude * np.exp(-t / tau_s) * np.cos(omega_Hz * t + phase)

    def _echo_frequency_shift(self, R: float, l: int = 2, m: int = 2) -> complex:
        """
        Echo cavity frequency shift for an ECO with surface reflectivity R.

        Near-horizon cavity model (Mark et al. 2017):
        omega_n ~ (pi * n + delta) / |r*_wall| - i/(2|r*_wall|) * ln|R * R_BH|

        For high reflectivity R -> 1:
        delta_omega ~ -i * ln(R) / (2 * T)
        where T = 2 * |r*_wall| is the cavity round-trip time.

        This implementation uses the cavity approximation.
        """
        M = self.metric.params.get("mass", 1.0)
        spin = self.metric.params.get("spin", 0.0)
        epsilon = 1e-10

        a_star = spin / M

        r_star_wall = 2.0 * M * abs(np.log(epsilon))

        T_roundtrip = 2.0 * r_star_wall

        if R <= 0:
            return 0.0 + 0.0j

        delta_omega = -1.0j * np.log(max(R, 1e-15)) / (2.0 * T_roundtrip)
        return delta_omega

    def echo_spacing(self) -> float:
        """
        Echo spacing Delta t_echo in seconds.

        For Schwarzschild:
        Delta t_echo ~ 4M * ln(1/epsilon)

        For Kerr:
        Delta t_echo ~ 4M * (1 + 1/sqrt(1 - a^2)) * ln(1/epsilon)

        where epsilon is the distance of the reflective surface from the horizon.
        """
        M = self.metric.params.get("mass", 1.0)
        spin = self.metric.params.get("spin", 0.0)
        epsilon = 1e-10

        M_sun_kg = 1.98847e30
        G = 6.67430e-11
        c = 299792458.0
        mass_kg = M * M_sun_kg
        M_in_meters = G * mass_kg / c ** 2
        time_in_seconds = M_in_meters / c

        a_star = spin / M

        if abs(a_star) < 1e-10:
            dt_echo_M = 8.0 * np.log(1.0 / epsilon)
        else:
            a_clamped = min(abs(a_star), 0.999)
            dt_echo_M = 4.0 * (1.0 + 1.0 / np.sqrt(1.0 - a_clamped ** 2)) * np.log(1.0 / epsilon)
            if a_clamped != abs(a_star):
                dt_echo_M = 8.0 * np.log(1.0 / epsilon)

        return float(dt_echo_M * time_in_seconds)

    def wormhole_qnm(self, l: int = 2) -> complex:
        """
        Wormhole QNM frequency using Poschl-Teller potential approximation.

        For a Morris-Thorne wormhole, the axial perturbation potential
        can be approximated by a Poschl-Teller potential:
        V_PT(r*) = V_0 / cosh^2(alpha * (r* - r*_0))

        Exact QNMs:
        omega_n = alpha * (sqrt(V_0/alpha^2 - 1/4) - i * (n + 1/2))

        For Ellis wormhole with throat b0:
        V_0 ~ l(l+1) / b0^2
        alpha ~ 1 / b0

        For l=2: omega_0 ~ (2.40 - 0.50*i) / b0
        """
        mt = self.metric.metric_type

        if mt == MetricType.MORRIS_THORNE:
            b0 = self.metric.params.get("throat_radius", 2.0)
        else:
            M = self.metric.params.get("mass", 1.0)
            b0 = max(2.0 * M, 2.0)

        V_0 = l * (l + 1) / (b0 ** 2)
        alpha = 1.0 / b0

        omega_real = alpha * np.sqrt(max(V_0 / (alpha ** 2) - 0.25, 0.0))
        omega_imag = -alpha * 0.5

        return complex(omega_real, omega_imag)

    def wormhole_ringdown_frequency(self, l: int = 2) -> float:
        """Wormhole ringdown frequency in Hz."""
        M = self.metric.params.get("mass", 1.0)
        M_sun_kg = 1.98847e30
        G = 6.67430e-11
        c = 299792458.0
        mass_kg = M * M_sun_kg
        M_in_meters = G * mass_kg / c ** 2
        time_in_seconds = M_in_meters / c

        omega = self.wormhole_qnm(l)
        f_RD = omega.real / time_in_seconds / (2.0 * np.pi)
        return float(f_RD)

    def qnm_deviation_from_kerr(self) -> Dict[str, float]:
        """
        Fractional deviation of the ringdown frequency from Kerr.

        For ECO with reflectivity: delta_omega from echo cavity.
        For wormhole: delta_omega from modified effective potential.

        Source: Mark et al. (2017); Cardoso & Pani (2019).
        """
        mt = self.metric.metric_type
        if mt in (MetricType.SCHWARZSCHILD, MetricType.KERR):
            return {"delta_omega": 0.0, "delta_tau": 0.0}

        if mt == MetricType.MORRIS_THORNE:
            M = self.metric.params.get("mass", 1.0)
            if M < 1e-10:
                return {"delta_omega": 1.0, "delta_tau": 0.5}
            omega_wh = self.wormhole_qnm(2)
            omega_kerr, _ = self._kerr_qnm_fit(2, 2, 0)
            delta_omega = (omega_wh.real - omega_kerr) / omega_kerr
            return {"delta_omega": float(delta_omega), "delta_tau": 0.5}

        R = self.metric.params.get("reflection_coefficient", 0.0)
        delta_omega_complex = self._echo_frequency_shift(R)
        delta_omega = -abs(delta_omega_complex.imag) / (2.0 * np.pi * 100.0)
        delta_tau = R * 0.1

        return {"delta_omega": delta_omega, "delta_tau": delta_tau}

    def tidal_deformability(self) -> float:
        """
        Tidal deformability (Love number) lambda.

        BH (Schwarzschild, Kerr): lambda = 0.
        Boson star: lambda > 0 (depends on compactness).
        Gravastar: lambda from thin-shell matching.

        Source: Damour & Nagar (2009); Cardoso & Pani (2019).
        """
        if self.metric.metric_type in (MetricType.SCHWARZSCHILD, MetricType.KERR):
            return 0.0

        tidal_param = self.metric.params.get("tidal_deformability", 0.0)
        return tidal_param

    def inspiral_overlap(self) -> float:
        if self.metric.metric_type in (MetricType.SCHWARZSCHILD, MetricType.KERR):
            return 1.0
        R = self.metric.params.get("reflection_coefficient", 0.0)
        return 1.0 - R * 0.1
