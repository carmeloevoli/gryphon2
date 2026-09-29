"""Birth-property distributions -> gryphon_Erot_distribution.pdf, gryphon_Emax_distribution.pdf

Source: `inspectRandomPulsars` on the fiducial model, which samples its pulsar
birth distributions and dumps one row per realization:

    P0 [ms]  B0 [G]  rotational energy [erg]  Emax [GeV]  tau0 [kyr]

Both quantities span decades, so they are histogrammed per dex: the y axis is
dN/dlog10(x), which integrates to one over the plotted range.
"""

import numpy as np

import gryphon_plots as gp

PULSARS_FILE = "youngpulsars_H5/pulsars_000000.txt"
BINS_PER_DEX = 20


def plot_distribution(values, xlabel, unit, color, figname):
    """One normalized log-binned distribution, with its median marked."""
    lo = np.floor(np.log10(values.min()))
    hi = np.ceil(np.log10(values.max()))
    nbins = int((hi - lo) * BINS_PER_DEX)
    edges = np.logspace(lo, hi, nbins + 1)

    counts, _ = np.histogram(values, bins=edges)
    density = counts / (len(values) / BINS_PER_DEX)  # dN / dlog10
    median = np.median(values)

    fig, ax = gp.new_figure(figsize=(10.5, 8.0))
    gp.set_axes(
        ax,
        xlabel=xlabel,
        ylabel=r"dN/dlog$_{10}$",
        xscale="log",
        xlim=[10**lo, 10**hi],
        ylim=[0.0, 1.35 * density.max()],
    )

    ax.stairs(density, edges, fill=True, color=color, alpha=0.5, zorder=2)
    ax.stairs(density, edges, color=color, lw=2.5, zorder=3)
    # stop the median line short of the caption line
    ax.vlines(median, 0.0, 1.1 * density.max(), color="k", linestyle="--", lw=3, zorder=4)

    exponent = int(np.floor(np.log10(median)))
    ax.text(
        0.04,
        0.92,
        rf"median $= {median / 10 ** exponent:.1f} \times 10^{{{exponent}}}$ {unit}",
        transform=ax.transAxes,
        fontsize=22,
    )

    gp.savefig(fig, figname)
    return median


def plot_erot(pulsars):
    return plot_distribution(
        pulsars["rotational energy"],
        xlabel=r"E$_{\rm rot}$ [erg]",
        unit="erg",
        color="tab:blue",
        figname="gryphon_Erot_distribution",
    )


def plot_emax(pulsars):
    return plot_distribution(
        pulsars["Emax"],
        xlabel=r"E$_{\rm max}$ [GeV]",
        unit="GeV",
        color="tab:red",
        figname="gryphon_Emax_distribution",
    )


if __name__ == "__main__":
    gp.style.use()
    catalog = gp.load(PULSARS_FILE)
    print(f"{len(catalog)} sampled pulsars")
    print(f"median E_rot : {plot_erot(catalog):.3g} erg")
    print(f"median E_max : {plot_emax(catalog):.3g} GeV")
