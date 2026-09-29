"""Figure 5: fixed, unselected current-code realizations, not tail statistics."""

from __future__ import annotations

import hashlib
import json

import matplotlib.pyplot as plt
import numpy as np

import gryphon_plots as gp

SEEDS = tuple(range(10))
SLOPE = 2.7
TARGET = 1e4  # mean of E^2.7 I_median over the common logarithmic energy grid
DATASETS = (
    ("CALET_H_kineticEnergy.txt", "^", "CALET", 12),
    ("DAMPE_H_kineticEnergy.txt", "s", "DAMPE", 10),
    ("CREAM_H_kineticEnergy.txt", "p", "CREAM", 11),
)
EXPECTED = dict(spiralmodel="Xie2024", transportmodel="PureDiffusion",
                injectionmodel="GalacticRandom", d0h=0.42, e0=1000.,
                delta=0.36, ddelta=-1., injslope=2.34, injemax=0.,
                efficiency=0.1, snrateyr=0.02, maxtimemyr=100.,
                emin=100., emax=1e6, esize=64., sunkpc=8.5,
                discsizepc=50., rgkpc=20., varyenergy="false",
                varyslope="false", syntheticassociations="false", pid="H")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_params(params, halo):
    for key, expected in dict(EXPECTED, hkpc=float(halo)).items():
        value = params.get(key)
        equal = (isinstance(value, (int, float)) and
                 np.isclose(value, expected, rtol=1e-10, atol=0)) if isinstance(
                     expected, float) else value == expected
        if not equal:
            raise ValueError(f"H={halo}: expected {key}={expected}, got {value}")


def scaled_spectra(energy, flux):
    """One scalar per halo, never per realization or energy-dependent detrending."""
    energy, flux = np.asarray(energy), np.asarray(flux)
    if (energy.ndim != 1 or flux.ndim != 2 or flux.shape[1] != energy.size
            or not np.all(np.isfinite(energy)) or not np.all(energy > 0)
            or not np.all(np.diff(energy) > 0)
            or not np.all(np.isfinite(flux)) or not np.all(flux > 0)):
        raise ValueError("expected finite positive spectra on an increasing energy grid")
    weighted = flux * energy**SLOPE
    median = np.median(weighted, axis=0)
    factor = TARGET / np.mean(median)
    return factor * weighted, factor * median, float(factor)


def load_model(halo, proj):
    directory = proj.run(f"figure5_H{halo}")
    manifest = json.loads((directory / "figure5_manifest.json").read_text())
    if manifest["seeds"] != list(SEEDS):
        raise ValueError("Figure 5 requires fixed seeds 0--9")
    params = gp.load_params(directory)
    validate_params(params, halo)
    energy, rows, records = None, [], []
    for seed in SEEDS:
        path = directory / f"flux_{seed:06d}.txt"
        saved = json.loads((directory / f"complete_{seed:06d}.json").read_text())
        checksum = digest(path)
        if saved != {"flux_sha256": checksum}:
            raise ValueError(f"{path}: checksum mismatch")
        table = gp.load(path)
        if int(table.meta["seed"]) != seed:
            raise ValueError(f"{path}: seed mismatch")
        if energy is None:
            energy = table["E"]
        elif not np.array_equal(energy, table["E"]):
            raise ValueError("inconsistent energy grids")
        rows.append(table["I"])
        records.append(dict(seed=seed, flux_sha256=checksum, path=str(path)))
    if not np.allclose(energy, np.geomspace(100., 1e6, 64), rtol=1e-8, atol=0):
        raise ValueError("unexpected Figure 5 energy grid")
    curves, median, factor = scaled_spectra(energy, rows)
    summary = dict(halo_kpc=halo, manifest=manifest, params=params, files=records,
                   display_flux_factor=factor, effective_efficiency=0.1 * factor,
                   normalization="one factor per halo: mean(E^2.7 I_median) = 10000",
                   plotted_range=[float(curves.min()), float(curves.max())])
    return energy, curves, median, summary


def main():
    gp.style.use()
    proj = gp.project()
    models = {halo: load_model(halo, proj) for halo in (2, 4)}
    # Shared axes; extend them if needed so no realization is hidden.
    low = min(6000., *(0.95 * model[1].min() for model in models.values()))
    high = max(16000., *(1.05 * model[1].max() for model in models.values()))
    for halo, (energy, curves, median, summary) in models.items():
        fig, ax = gp.new_figure(figsize=(12.5, 8.5))
        gp.set_axes(ax, xlabel=r"$E$ [GeV]",
                    ylabel=r"$E^{2.7} I$ [GeV$^{1.7}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]",
                    xscale="log", xlim=(100., 1e6), ylim=(low, high))
        ax.grid(False, which="both")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(3, 3), useMathText=True)
        # Experimental points are contextual; no fits, energy shifts or data rescaling.
        for filename, marker, label, zorder in DATASETS:
            gp.plot_data(ax, filename, slope=SLOPE, fmt=marker, color="silver",
                         label=label, zorder=zorder, emin=100., emax=1e6)
        for curve, color in zip(curves, plt.get_cmap("tab10").colors):
            ax.plot(energy, curve, color=color, lw=2.6, zorder=20)
        ax.plot(energy, median, color="black", ls="--", lw=3.5, zorder=21,
                label="Median of 10 realizations")
        ax.text(0.97, 0.96, rf"$H={halo}$ kpc", ha="right", va="top",
                transform=ax.transAxes, fontsize=25)
        ax.legend(loc="upper left", fontsize=18)
        gp.savefig(fig, f"gryphon_realizations_H{halo}")
        print(f"H={halo}: display factor {summary['display_flux_factor']:.6g}; "
              f"effective efficiency {summary['effective_efficiency']:.6g}")
    summary = dict(seeds=list(SEEDS), number_per_halo=len(SEEDS),
                   median="pointwise median of the ten displayed realizations only",
                   selection="no selection by flux or spectral features",
                   data_errors="statistical and systematic errors in quadrature",
                   data={name: digest(proj.kiss / name) for name, *_ in DATASETS},
                   models={str(h): m[3] for h, m in models.items()})
    (proj.figures / "figure5_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
