"""Burst injection spectrum of a single young pulsar -> gryphon_spectrum.pdf

Source: `inspectRandomPulsars` on the fiducial model, with the random birth
properties disabled, i.e. the reference pulsar sitting at the centre of the
birth distributions. Its P0 and B_* are read back from the run.
"""

import numpy as np

import gryphon_plots as gp

RUN = "youngpulsars_H5"
SPECTRUM_FILE = f"{RUN}/injection_000000.txt"


def plot_injection_spectrum(figname="gryphon_spectrum"):
    table = gp.load(SPECTRUM_FILE)
    E, Q = table.columns("E", "Q(E)")

    params = gp.load_params(RUN)
    p0_msec = params["youngpulsarsp0ms"]
    b_star_gauss = params["youngpulsarsb0gauss"]
    E2Q = np.power(E, 2.0) * Q

    fig, ax = gp.new_figure(figsize=(10.5, 8.0))
    gp.set_axes(
        ax,
        xlabel="E [GeV]",
        ylabel=r"E$^2$ Q(E) [GeV]",
        xscale="log",
        xlim=[E.min(), E.max()],
        yscale="log",
        ylim=[1e-3 * E2Q.max(), 3.0 * E2Q.max()],
    )

    ax.plot(E, E2Q, color="tab:blue", lw=4, zorder=3)

    ax.text(
        0.05,
        0.12,
        rf"$P_0 = {p0_msec:.0f}$ ms",
        transform=ax.transAxes,
        fontsize=26,
    )
    ax.text(
        0.05,
        0.03,
        rf"$B_* = {b_star_gauss / 1e12:.1f} \times 10^{{12}}$ G",
        transform=ax.transAxes,
        fontsize=26,
    )

    gp.savefig(fig, figname)


if __name__ == "__main__":
    gp.style.use()
    plot_injection_spectrum()
