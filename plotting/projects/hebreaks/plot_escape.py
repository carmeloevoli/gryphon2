"""Figure 2 left: escape times from the paper's pure-diffusion run."""

import numpy as np

import gryphon_plots as gp
from figure2_common import HALOS, RUN, escape_time_myr, load_reference


def plot_escape():
    params = load_reference()
    table = gp.load(f"{RUN}/propagation_000000.txt")
    energy, reference_time = table.columns("E", "t_diff")
    # The C++ diagnostic is the input, not a legacy curve rescaled by eye.
    np.testing.assert_allclose(
        reference_time, escape_time_myr(energy, params["hkpc"], params),
        rtol=2e-6, atol=0,
    )

    fig, ax = gp.new_figure(figsize=(9.5, 9.0))
    gp.set_axes(ax, xlabel=r"$E$ [GeV]", ylabel=r"$\tau_{\rm esc}$ [Myr]",
                xscale="log", yscale="log", xlim=(1e3, 1e6), ylim=(0.1, 20))
    ax.set_box_aspect(1)
    ax.grid(False, which="both")
    for halo, color in HALOS:
        ax.plot(energy, reference_time * halo / params["hkpc"],
                color=color, lw=3, label=rf"$H={halo:g}$ kpc")
    ax.legend(loc="upper right", fontsize=23)
    ax.text(0.06, 0.10, r"$\tau_{\rm esc}=H^2/(2D)$", transform=ax.transAxes,
            fontsize=25)
    ax.text(0.06, 0.04, rf"$D_0/H={params['d0h']:g}$ kpc Myr$^{{-1}}$",
            transform=ax.transAxes, fontsize=21)
    return gp.savefig(fig, "gryphon_escape")


if __name__ == "__main__":
    gp.style.use()
    plot_escape()
