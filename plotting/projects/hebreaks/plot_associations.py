"""Matched synthetic-association comparison with finite-ensemble uncertainty."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

import gryphon_plots as gp
from feature_statistics import binomial_interval, excursions, summarize

MODELS = ("independent", "mixed", "clustered")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_model(directory: Path):
    completed = sorted(directory.glob("complete_*.json"))
    if not completed:
        raise ValueError(f"{directory}: no completed benchmark realizations")
    spectra, populations, seeds = [], [], []
    energy = None
    for marker in completed:
        seed = int(marker.stem.split("_")[-1])
        flux_file = directory / f"flux_{seed:06d}.txt"
        pop_file = directory / f"population_{seed:06d}.txt"
        if json.loads(marker.read_text()) != {"flux_sha256": digest(flux_file),
                                              "population_sha256": digest(pop_file)}:
            raise ValueError(f"{directory}: seed {seed} failed checksum verification")
        table, pop = gp.load(flux_file), gp.load(pop_file)
        if energy is None:
            energy = table["E"]
        elif not np.array_equal(table["E"], energy):
            raise ValueError(f"{directory}: inconsistent energy grids")
        spectra.append(table["I"])
        populations.append([pop[name][0] for name in
                            ("parents", "field_SN", "clustered_SN", "retained_SN", "expected_SN")])
        seeds.append(seed)
    return energy, np.stack(spectra), np.array(populations), seeds


def main() -> None:
    gp.style.use()
    import matplotlib.pyplot as plt

    proj = gp.project()
    fig, axes = plt.subplots(1, 2, figsize=(16, 7.7))
    thresholds = ((0.16, 0.19), (0.20, 0.35))
    summary = {"status": "pilot; rare tails require larger ensembles",
               "fit_range_gev": [1000., 1e6], "search_range_gev": [1e4, 1e5],
               "running_slope_window_bins": 5,
               "intervals": "pointwise exact 95% Clopper-Pearson; not simultaneous bands",
               "models": {}}
    reference_params = reference_seeds = reference_energy = None
    rows = []
    for model, color in zip(MODELS, ("tab:blue", "tab:orange", "tab:red")):
        directory = proj.run(f"associations_{model}")
        energy, flux, population, seeds = load_model(directory)
        params = gp.load_params(directory)
        fraction = params["associationfraction"]
        label = rf"{model.capitalize()}: $f_{{\rm a}}={fraction:g}$"
        matched = {key: value for key, value in params.items()
                   if key not in ("simname", "seed", "associationfraction")}
        if reference_params is None:
            reference_params, reference_seeds, reference_energy = matched, seeds, energy
        elif matched != reference_params or seeds != reference_seeds or not np.array_equal(energy, reference_energy):
            raise ValueError("comparison requires matched parameters, seed ranges, and energy grids")
        residual, slope = excursions(energy, flux)
        all_sn = population[:, 1] + population[:, 2]
        result = {"seeds": seeds, "parameters": params,
                  "provenance": json.loads((directory / "benchmark_manifest.json").read_text()),
                  "population": {"mean_explosions": float(all_sn.mean()),
                                 "expected_explosions": float(population[:, 4].mean()),
                                 "mean_clustered_fraction": float(population[:, 2].sum() / all_sn.sum()),
                                 "mean_retained_explosions": float(population[:, 3].mean()),
                                 "standard_error_mean_explosions": float(all_sn.std(ddof=1) / np.sqrt(len(seeds)))},
                  "residual": summarize(residual, thresholds[0]),
                  "slope_excursion": summarize(slope, thresholds[1]),
                  "mean_flux_1tev": float(flux[:, 0].mean()),
                  "standard_error_mean_flux_1tev": float(flux[:, 0].std(ddof=1) / np.sqrt(len(seeds)))}
        summary["models"][model] = result
        for seed, a, slope_range in zip(seeds, residual, slope):
            rows.append((model, seed, a, slope_range))
        for ax, values, band in zip(axes, (residual, slope), thresholds):
            x = np.geomspace(1e-3, 0.8, 280)
            counts = np.count_nonzero(values[:, None] >= x, axis=0)
            lo, hi = binomial_interval(counts, len(values))
            ax.fill_between(x, np.maximum(lo, 1e-8), hi, color=color, alpha=0.10, linewidth=0)
            nonzero = counts > 0
            ax.plot(x, np.where(nonzero, counts / len(values), np.nan), color=color, label=label)
    n = len(reference_seeds)
    for ax, band, xlabel in zip(axes, thresholds,
                               (r"Maximum residual $A$", r"Running-slope excursion $\Delta\gamma_{\rm run}$")):
        ax.axvspan(*band, color="0.65", alpha=0.30, label="Data-motivated range")
        ax.axhline(1. / n, color="0.45", linestyle=":", linewidth=1.2)
        ax.set(xscale="log", yscale="log", xlim=(1e-3, 0.8), ylim=(max(1e-6, 0.3 / n), 1.05),
               xlabel=xlabel, ylabel=r"Survival probability $P(\geq x)$")
        ax.grid(False, which="both")
        ax.tick_params(labelsize=23)
    axes[0].legend(fontsize=18, loc="upper right")
    members = reference_params["associationmembers"]
    halo = reference_params["hkpc"]
    axes[1].text(0.97, 0.95, rf"Pilot: {n} realizations/model" + "\n" +
                rf"$H={halo:g}$ kpc, $N_{{\rm a}}={members:g}$" +
                "\n" + r"$10<E<100$ TeV", transform=axes[1].transAxes,
                va="top", ha="right", fontsize=18)
    proj.figures.mkdir(parents=True, exist_ok=True)
    fig.savefig(proj.figures / "gryphon_association_statistics.pdf", bbox_inches="tight")
    fig.savefig(proj.figures / "gryphon_association_statistics.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    (proj.figures / "association_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (proj.figures / "association_statistics.txt").open("w") as out:
        out.write("# model seed residual slope_excursion\n")
        for model, seed, residual, slope in rows:
            out.write(f"{model} {seed} {residual:.10e} {slope:.10e}\n")
    print(json.dumps({model: {"A95": data["residual"]["percentiles"]["95"],
                              "slope95": data["slope_excursion"]["percentiles"]["95"],
                              "events": data["population"]}
                      for model, data in summary["models"].items()}, indent=2))


if __name__ == "__main__":
    main()
