import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, Any, Optional, List
from ..metrics.base import AbstractMetric


class ComparisonPlotter:
    """
    Publication-quality figure generation for metric comparisons.
    Supports multi-object overlay plots for research papers.
    """

    def __init__(self, style: str = "default"):
        self.style = style
        plt.style.use("seaborn-v0_8-colorblind")
        self.colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    def _style_ax(self, ax, xlabel: str, ylabel: str, title: str):
        ax.set_xlabel(xlabel, fontsize=14)
        ax.set_ylabel(ylabel, fontsize=14)
        ax.set_title(title, fontsize=16)
        ax.legend(fontsize=11, loc="best")
        ax.grid(True, alpha=0.3)

    def plot_shadow_boundary(
        self,
        boundaries: Dict[str, np.ndarray],
        save_path: Optional[str] = None,
    ):
        fig, ax = plt.subplots(figsize=(8, 8))
        for i, (label, boundary) in enumerate(boundaries.items()):
            c = self.colors[i % len(self.colors)]
            ax.plot(boundary[:, 0], boundary[:, 1], "-",
                    label=label, linewidth=2, color=c)
        self._style_ax(ax, r"$\alpha$ [$M$]", r"$\beta$ [$M$]",
                       "Shadow Boundary Comparison")
        ax.set_aspect("equal")
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()

    def plot_shadow_overlay(self, boundaries: Dict[str, np.ndarray],
                            save_path: Optional[str] = None):
        """Shadow overlay with all object classes.
        First boundary in dict is used for inset zoom."""
        fig, (ax_main, ax_zoom) = plt.subplots(1, 2, figsize=(14, 7),
                                                gridspec_kw={"width_ratios": [1, 0.7]})

        for i, (label, boundary) in enumerate(boundaries.items()):
            c = self.colors[i % len(self.colors)]
            ax_main.plot(boundary[:, 0], boundary[:, 1], "-",
                         label=label, linewidth=2, color=c)
            ax_zoom.plot(boundary[:, 0], boundary[:, 1], "-",
                         linewidth=2, color=c)

        self._style_ax(ax_main, r"$\alpha$ [$M$]", r"$\beta$ [$M$]",
                       "All-Object Shadow Comparison")
        ax_main.set_aspect("equal")

        ax_zoom.set_xlim(-5.5, 5.5)
        ax_zoom.set_ylim(-5.5, 5.5)
        self._style_ax(ax_zoom, r"$\alpha$ [$M$]", r"$\beta$ [$M$]",
                       "Zoom: Photon Ring Region")
        ax_zoom.set_aspect("equal")

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()

    def plot_deflection_angle(
        self, comparison_data: Dict[str, np.ndarray],
        metric_name: str = "Model",
        save_path: Optional[str] = None,
    ):
        fig, ax = plt.subplots(figsize=(8, 6))
        b = comparison_data["impact_parameter"]
        ax.plot(b, comparison_data["alpha_metric"], "-",
                label=metric_name, linewidth=2)
        if "alpha_schwarzschild" in comparison_data:
            ax.plot(b, comparison_data["alpha_schwarzschild"], "--",
                    label="Schwarzschild", linewidth=2, alpha=0.7)
        self._style_ax(ax, "Impact parameter b [M]",
                       "Deflection angle [rad]",
                       "Gravitational Lensing Comparison")
        ax.set_yscale("log")
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()

    def plot_deflection_overlay(self, deflection_data: Dict[str, Dict],
                                save_path: Optional[str] = None):
        fig, ax = plt.subplots(figsize=(8, 6))
        for i, (label, dd) in enumerate(deflection_data.items()):
            c = self.colors[i % len(self.colors)]
            ax.plot(dd["b"], dd["alpha"], "-",
                    label=label, linewidth=2, color=c)
        self._style_ax(ax, "Impact parameter b [M]",
                       "Deflection angle [rad]",
                       "Deflection Angle: All Objects")
        ax.set_yscale("log")
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()

    def plot_ringdown_waveform(
        self, waveforms: Dict[str, Dict[str, np.ndarray]],
        save_path: Optional[str] = None,
    ):
        fig, ax = plt.subplots(figsize=(10, 6))
        for label, wf in waveforms.items():
            ax.plot(wf["t"], wf["h"], "-", label=label, linewidth=1.5)
        self._style_ax(ax, "Time t [s]", "Strain h(t)",
                       "Ringdown Waveform Comparison")
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()

    def plot_ringdown_overlay(self, ringdown_data: Dict[str, Dict],
                              save_path: Optional[str] = None):
        fig, (ax_wf, ax_qnm) = plt.subplots(1, 2, figsize=(14, 6))

        for i, (label, rd) in enumerate(ringdown_data.items()):
            c = self.colors[i % len(self.colors)]
            if "t" in rd and "h" in rd:
                ax_wf.plot(rd["t"], rd["h"], "-", label=label,
                          linewidth=1.5, color=c)

        self._style_ax(ax_wf, "Time t [s]", "Strain h(t)",
                       "Ringdown Waveforms")

        labels = []
        freqs = []
        for i, (label, rd) in enumerate(ringdown_data.items()):
            if "f_RD" in rd:
                labels.append(label)
                freqs.append(rd["f_RD"])
                c = self.colors[i % len(self.colors)]
                ax_qnm.bar(i, rd["f_RD"], color=c, alpha=0.7)

        ax_qnm.set_xticks(range(len(labels)))
        ax_qnm.set_xticklabels(labels, rotation=45, ha="right", fontsize=10)
        self._style_ax(ax_qnm, "", "f_RD [Hz]",
                       "Ringdown Frequency Comparison")

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()

    def plot_orbital_precession(
        self, precession_data: Dict[str, np.ndarray],
        save_path: Optional[str] = None,
    ):
        fig, ax = plt.subplots(figsize=(8, 6))
        r = precession_data["radius"]
        ax.plot(r, precession_data["delta_phi_metric"], "-",
                label=getattr(self, "metric_name", "Model"), linewidth=2)
        if "delta_phi_schwarzschild" in precession_data:
            ax.plot(r, precession_data["delta_phi_schwarzschild"], "--",
                    label="Schwarzschild", linewidth=2, alpha=0.7)
        self._style_ax(ax, "Orbital radius r [M]",
                       "Precession per orbit [rad]",
                       "Periastron Precession Comparison")
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()

    def plot_exclusion_contours(self, chi2_data: Dict[str, float],
                                save_path: Optional[str] = None):
        fig, ax = plt.subplots(figsize=(10, 6))
        labels = list(chi2_data.keys())
        scores = list(chi2_data.values())
        colors = [self.colors[i % len(self.colors)] for i in range(len(labels))]

        bars = ax.bar(range(len(labels)), scores, color=colors, alpha=0.7)
        ax.axhline(y=2.0, color="red", linestyle="--", alpha=0.5,
                   label=r"$\chi^2=2$ (good fit)")
        ax.axhline(y=6.0, color="orange", linestyle="--", alpha=0.5,
                   label=r"$\chi^2=6$ (tension)")

        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=11)
        for bar, score in zip(bars, scores):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                    f"{score:.1f}", ha="center", va="bottom", fontsize=10)

        self._style_ax(ax, "", r"$\chi^2$",
                       "Exclusion: Agreement with Observations")
        ax.set_yscale("log")
        ax.legend(fontsize=11)
        ax.set_ylim(bottom=0.1)

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()

    def plot_unified_exclusion(self, shadow_boundaries, chi2_data,
                                ringdown_data, save_path=None):
        """
        GPT-recommended unified exclusion diagram (3 panels).

        Panel A: Shadow morphology (all objects)
        Panel B: Exclusion contours (chi2)
        Panel C: Ringdown consistency (f_RD deviation)

        Single most impactful figure for an ECO paper.
        """
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        ax = axes[0]
        for i, (label, boundary) in enumerate(shadow_boundaries.items()):
            c = colors[i % len(colors)]
            ax.plot(boundary[:, 0], boundary[:, 1], "-",
                    label=label, linewidth=2, color=c)
        ax.set_xlabel(r"$\alpha$ [$M$]", fontsize=14)
        ax.set_ylabel(r"$\beta$ [$M$]", fontsize=14)
        ax.set_title("A: Shadow Morphology", fontsize=15, fontweight="bold")
        ax.set_aspect("equal")
        ax.legend(fontsize=9, loc="upper right", framealpha=0.8)
        ax.grid(True, alpha=0.3)

        ax = axes[1]
        labels = list(chi2_data.keys())
        scores = list(chi2_data.values())
        bar_colors = [colors[i % len(colors)] for i in range(len(labels))]
        bars = ax.bar(range(len(labels)), scores, color=bar_colors, alpha=0.8,
                       edgecolor="black", linewidth=0.5)
        ax.axhline(y=2.0, color="green", linestyle="--", alpha=0.6,
                   linewidth=1.5, label=r"$\chi^2=2$ (consistent)")
        ax.axhline(y=6.0, color="red", linestyle="--", alpha=0.6,
                   linewidth=1.5, label=r"$\chi^2=6$ (excluded)")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
        ax.set_ylabel(r"$\chi^2$", fontsize=13)
        ax.set_title("B: Observational Exclusion", fontsize=15, fontweight="bold")
        ax.legend(fontsize=10)
        ax.set_yscale("log")
        ax.set_ylim(bottom=0.3, top=max(scores) * 3)
        for bar, score in zip(bars, scores):
            lbl = "EXCL" if score > 6 else ("OK" if score < 2 else "TEN")
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.2,
                    lbl, ha="center", va="bottom", fontsize=7, fontweight="bold")

        ax = axes[2]
        markers = ["o", "s", "D", "^", "v", "<", ">"]
        for i, (label, rd) in enumerate(ringdown_data.items()):
            c = colors[i % len(colors)]
            m = markers[i % len(markers)]
            if "f_RD" in rd and "delta_omega" in rd:
                ax.scatter(rd["f_RD"], rd["delta_omega"], marker=m, s=100,
                           color=c, label=label, edgecolor="black",
                           linewidth=0.5, zorder=5)
        ax.axhline(y=0.0, color="black", linestyle="-", alpha=0.3)
        ax.axhline(y=0.05, color="red", linestyle="--", alpha=0.5,
                   label=r"LIGO bound ($\delta f/f < 0.05$)")
        ax.axhline(y=-0.05, color="red", linestyle="--", alpha=0.5)
        ax.set_xlabel(r"$f_{\rm RD}$ [Hz] (10 $M_\odot$)", fontsize=13)
        ax.set_ylabel(r"$\delta f / f_{\rm Kerr}$", fontsize=13)
        ax.set_title("C: Ringdown Consistency", fontsize=15, fontweight="bold")
        ax.legend(fontsize=9, loc="upper right")
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.subplots_adjust(wspace=0.3)
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()
