"""Figure 2 right: actual events from the C++ Xie2024 source sampler.

No separate Python spiral model or old four-arm/bar overlays are used.
The five Xie components, their unequal weights, broadening, and central
component are all contained in the generated event positions.
"""

import matplotlib.pyplot as plt
import numpy as np

import gryphon_plots as gp
from figure2_common import HALOS, RUN, galactocentric_xy, horizon_kpc, load_reference


def plot_spirals():
    params = load_reference()
    events = gp.load(f"{RUN}/events_000000.txt")
    x, y = galactocentric_xy(events, params)
    r_sun = params["sunkpc"]
    radius = params["rgkpc"]

    fig, ax = gp.new_figure(figsize=(9.5, 9.0))
    gp.set_axes(ax, xlabel=r"$x$ [kpc]", ylabel=r"$y$ [kpc]",
                xlim=(-radius, radius), ylim=(-radius, radius))
    ax.set_aspect("equal")
    ax.grid(False, which="both")
    ax.set_xticks(np.arange(-20, 21, 10))
    ax.set_yticks(np.arange(-20, 21, 10))

    for halo, color in reversed(HALOS):
        horizon = horizon_kpc(halo)
        ax.add_patch(plt.Circle((r_sun, 0), horizon, facecolor=color,
                                edgecolor="none", alpha=0.08, zorder=1))
        ax.add_patch(plt.Circle((r_sun, 0), horizon, fill=False,
                                edgecolor=color, linewidth=2, zorder=3))
        ax.text(r_sun + halo - 0.35, -halo, rf"$H={halo:g}$",
                color=color, fontsize=19, ha="right", va="center", zorder=5,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))

    ax.scatter(x, y, s=3.0, c="0.25", alpha=0.50, edgecolors="none",
               rasterized=True, zorder=2)
    ax.scatter(0, 0, marker="+", s=130, linewidths=2.5, c="black",
               label="Galactic centre", zorder=6)
    ax.scatter(r_sun, 0, marker="*", s=260, c="gold", edgecolors="black",
               linewidths=1, label="Sun", zorder=7)
    ax.text(0.04, 0.96, "Xie et al. (2024)", fontsize=22,
            transform=ax.transAxes, va="top")
    ax.text(0.04, 0.90, rf"$T={params['maxtimemyr']:g}$ Myr", fontsize=22,
            transform=ax.transAxes, va="top")
    ax.legend(loc="lower left", fontsize=19, handletextpad=0.3)
    print(f"Xie2024 snapshot: {len(events)} events; Sun at ({r_sun:g}, 0) kpc")
    return gp.savefig(fig, "gryphon_spirals")


if __name__ == "__main__":
    gp.style.use()
    plot_spirals()
