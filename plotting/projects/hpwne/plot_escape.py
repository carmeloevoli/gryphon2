"""Diffusive escape time versus proton energy -> gryphon_escape.pdf

Source: `inspectPropagation` run on the fiducial model, so the curve is the one
of the population that produces the flux. Since D scales linearly with H in
gryphon (D0/H is the input parameter), t_diff = H^2 / 2D is itself linear in H,
and the other halo sizes are exact rescalings of the same curve.
"""

import gryphon_plots as gp

RUN = "youngpulsars_H5"
PROPAGATION = f"{RUN}/propagation_000000.txt"

HALO_SIZES = [
    (2.0, "tab:orange", 1.45),
    (4.0, "tab:red", 2.9),
    (8.0, "tab:purple", 5.8),
]


def plot_escape(figname="gryphon_escape"):
    table = gp.load(PROPAGATION)
    E, t_diff = table.columns("E", "t_diff")
    reference_h_kpc = gp.load_params(RUN)["hkpc"]

    fig, ax = gp.new_figure(figsize=(10.5, 8.0))
    gp.set_axes(
        ax,
        xlabel="E [GeV]",
        ylabel="escape timescale [Myr]",
        xscale="log",
        xlim=[1e3, 1e6],
        yscale="log",
        ylim=[0.1, 15],
    )

    for halo_kpc, color, label_y in HALO_SIZES:
        scaling = halo_kpc / reference_h_kpc
        ax.plot(E, scaling * t_diff, "-", color=color, lw=4)
        ax.text(2e3, label_y, rf"$H = {halo_kpc:.0f}$ kpc", rotation=-22, fontsize=20, color=color)

    ax.text(2e3, 0.2, r"$\tau \sim \frac{H^2}{2D}$", fontsize=40)

    gp.savefig(fig, figname)


if __name__ == "__main__":
    gp.style.use()
    plot_escape()
