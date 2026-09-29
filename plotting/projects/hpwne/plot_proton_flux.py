"""Proton flux from a Galactic population of young pulsars -> gryphon_proton_youngpulsars.pdf

Source: `runYoungPulsars`, one file per seed. The simulation is run with
efficiency = 1, so the CR efficiency xi is applied here as a plain rescaling of
the flux -- keep it explicit rather than folded into an anonymous factor.
"""

import numpy as np

import gryphon_plots as gp

SLOPE = 2.7  # flux is plotted as E^SLOPE * I
BAND = (5.0, 95.0)  # percentiles of the realization ensemble, as defined in the paper
NORMALIZATION_ENERGY = 3e6  # GeV; where the delta = 0 model is matched to the data
REFERENCE_DATA = "LHAASO_QGSJET-II-04_H_totalEnergy.txt"

MODELS = [
    # stem, xi, label, style
    (
        "youngpulsars_H5",
        0.1,
        r"$P_0 = 60$ ms, $\xi = 0.1$",
        dict(color="tab:blue", ls="-", lw=4, zorder=4),
    ),
    (
        "youngpulsars_H5_100msec",
        1.0,
        r"$P_0 = 100$ ms, $\xi = 1$",
        dict(color="tab:gray", ls="--", lw=4, zorder=3),
    ),
    (
        # delta = 0 with the slower birth spin: xi is set by matching the
        # measured proton flux at 3 PeV (see NORMALIZATION_ENERGY below).
        "youngpulsars_H5_130msec_flatD",
        0.15,
        r"$P_0 = 130$ ms, $\xi = 0.15$, $\delta = 0$",
        dict(color="tab:red", ls="-", lw=4, zorder=5),
    ),
]

DATASETS = [
    ("DAMPE_H_kineticEnergy.txt", "s", "tab:red", "DAMPE"),
    ("CALET_H_kineticEnergy.txt", "D", "tab:orange", "CALET"),
    ("ISS-CREAM_H_kineticEnergy.txt", "o", "tab:purple", "ISS-CREAM"),
    ("GRAPES-3_H_totalEnergy.txt", "v", "tab:brown", "GRAPES-3"),
    ("LHAASO_QGSJET-II-04_H_totalEnergy.txt", "^", "tab:green", "LHAASO"),
]

EMIN, EMAX = 1e4, 1e8


def plot_proton_flux(figname="gryphon_proton_youngpulsars"):
    fig, ax = gp.new_figure(figsize=(12.0, 8.5))
    gp.set_axes(
        ax,
        xlabel="E [GeV]",
        ylabel=r"E$^{2.7}$ I [GeV$^{1.7}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]",
        xscale="log",
        xlim=[EMIN, EMAX],
        yscale="log",
        ylim=[1e2, 3e4],
    )

    E_data, I_data, _, _ = gp.get_data(REFERENCE_DATA)
    I_reference = np.interp(NORMALIZATION_ENERGY, E_data, I_data)

    for stem, xi, label, style in MODELS:
        runs = gp.load_ensemble(stem)
        E = runs.x
        scaling = xi * E**SLOPE

        median = scaling * runs.median("I")
        lo, hi = runs.band("I", *BAND)

        ax.plot(E, median, label=label, **style)
        ax.fill_between(
            E, scaling * lo, scaling * hi, color=style["color"], alpha=0.25, lw=0, zorder=2
        )

        # The efficiency that would put this model on the data at the
        # normalization energy, so the value used above stays checkable.
        required = I_reference / np.interp(NORMALIZATION_ENERGY, E, runs.median("I"))
        print(
            f"{stem}: {len(runs)} realizations, xi = {xi} "
            f"(xi to match {NORMALIZATION_ENERGY:.0e} GeV: {required:.3f})"
        )

    for filename, fmt, color, label in DATASETS:
        gp.plot_data(
            ax, filename, slope=SLOPE, fmt=fmt, color=color, label=label,
            emin=EMIN, emax=EMAX, zorder=5,
        )

    ax.legend(fontsize=18, loc="lower left")
    gp.savefig(fig, figname)


if __name__ == "__main__":
    gp.style.use()
    plot_proton_flux()
