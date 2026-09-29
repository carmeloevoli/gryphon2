"""Pulsar positions of one realization on the Galactic plane -> gryphon_spirals.pdf

Source: `inspectEvents` on the Steiman-Cameron et al. (2010) spiral pattern,
integrated over a short enough time for one realization to stay plottable.
Event positions are Sun-centred, hence the -R_sun shift applied to the arm
model; R_sun and the integration time come from the run itself.
"""

import math

import matplotlib.pyplot as plt
import numpy as np

import gryphon_plots as gp

RUN = "events_1Myr"
EVENTS_FILE = f"{RUN}/events_000000.txt"

# Steiman-Cameron et al. (2010): pitch alpha, scale a [kpc], theta_0 [deg]
ARMS = [
    (0.242, 0.246, 13.6, "tab:green", "Sagittarius-Carina"),
    (0.279, 0.608, 15.6, "tab:red", "Scutum-Crux"),
    (0.249, 0.449, 13.5, "tab:blue", "Norma-Cygnus"),
    (0.240, 0.378, 13.5, "tab:orange", "Perseus"),
]

# Illustrative distances travelled by diffusion, drawn around the Sun
HORIZONS_KPC = [(2, "tab:red", 3), (4, "tab:orange", 2), (8, "tab:olive", 1)]


def plot_arm(ax, r_sun, alpha, a, theta_0_deg, color, label, rmin=3.0, rmax=20.0, npts=2000):
    r = np.linspace(rmin, rmax, npts)
    phi = np.log(r / a) / alpha + np.deg2rad(theta_0_deg)
    ax.plot(
        r * np.cos(-math.pi / 2.0 + phi) - r_sun,
        r * np.sin(-math.pi / 2.0 + phi),
        color=color,
        label=label,
        lw=2.05,
        zorder=6,
    )


def plot_bar(ax, r_sun, length_kpc=3.0, angle_deg=25.0):
    theta = np.deg2rad(angle_deg)
    x = np.array([length_kpc, -length_kpc]) * np.cos(theta) - r_sun
    y = np.array([-length_kpc, length_kpc]) * np.sin(theta)
    ax.plot(x, y, lw=6, color="k", label="MW bar", zorder=6)


def plot_sun(ax):
    ax.scatter(0.0, 0.0, s=59, label="Sun", color="darkgray", zorder=6)
    ax.hlines(0, -30, 30, color="darkgray", linestyle="--", lw=1.5, zorder=4)
    ax.vlines(0, -30, 30, color="darkgray", linestyle="--", lw=1.5, zorder=4)


def plot_spirals(figname="gryphon_spirals"):
    events = gp.load(EVENTS_FILE)
    params = gp.load_params(RUN)
    r_sun = params["sunkpc"]

    fig, ax = gp.new_figure(figsize=(10.0, 9.0))
    gp.set_axes(ax, xlabel="x [kpc]", ylabel="y [kpc]", xlim=[-20, 20], ylim=[-20, 20])

    for alpha, a, theta_0_deg, color, label in ARMS:
        plot_arm(ax, r_sun, alpha, a, theta_0_deg, color, label)

    plot_bar(ax, r_sun)
    plot_sun(ax)

    ax.scatter(events["x"], events["y"], c="tab:gray", s=4.5, alpha=0.75,
               edgecolors="none", zorder=5)
    ax.text(8, 15, rf"T = {params['maxtimemyr']:.0f} Myr", fontsize=28)

    for radius, color, zorder in HORIZONS_KPC:
        ax.add_patch(plt.Circle((0, 0), radius, color=color, alpha=0.2, zorder=zorder))

    ax.legend(loc="lower right", fontsize=14)
    gp.savefig(fig, figname)


if __name__ == "__main__":
    gp.style.use()
    plot_spirals()
