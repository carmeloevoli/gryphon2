"""Verified 100 GeV--1 PeV ensembles with primary and control search windows.

Writes a separate analysis directory. Never updates the manuscript in place.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
import gryphon_plots as gp
from feature_statistics import excursions, running_indices, summarize
from plot_halos import THRESHOLDS, digest, load_model, survival
from plot_variations import continuous_index_spectrum

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from run_results import NAMES, expected_parameters, validate_parameters, atomic_json

GRID = (100., 1e6, 65)
WINDOWS = {
    "primary_1_100_TeV": ((100., 1e6), (1e3, 1e5)),
    "secondary_10_100_TeV": ((100., 1e6), (1e4, 1e5)),
    # Separates the changed reference fit from the changed search window.
    # Grid spacing and history still differ from the archived calculation.
    "legacy_reference_10_100_TeV": ((1e3, 1e6), (1e4, 1e5)),
}
GROUPS = {
    "halos": (("halos_H2", "$H=2$ kpc", "tab:orange"),
              ("halos_H4", "$H=4$ kpc", "tab:red"),
              ("halos_H8", "$H=8$ kpc", "tab:purple")),
    "variations": (("halos_H4", "Identical", "tab:red"),
                   ("variations_energy", "Energy scatter", "tab:blue"),
                   ("variations_index", "Index scatter", "tab:green")),
    "low_rate": (("halos_H4", "Reference rate", "tab:orange"),
                 ("low_rate", "One tenth the rate", "tab:blue")),
    "associations": (("associations_independent", "$f_a=0$", "tab:blue"),
                     ("associations_mixed", "$f_a=0.5$", "tab:orange"),
                     ("associations_clustered", "$f_a=1$", "tab:red")),
}


def manifest_name(name):
    if name.startswith("halos_"):
        return "halos_manifest.json"
    if name.startswith("variations_"):
        return "variations_manifest.json"
    return "low_rate_manifest.json" if name == "low_rate" else "benchmark_manifest.json"


def statistic_summary(a, d):
    return {"residual": summarize(a, THRESHOLDS[0]),
            "slope_excursion": summarize(d, THRESHOLDS[1]),
            "joint_lower_thresholds": summarize(((a >= .16) & (d >= .20)).astype(float), (.5,))}


def analyse(runs, number):
    if number < 2:
        raise ValueError("at least two realizations are required")
    summary = {
        "status": "complete requested prefix; amplitude diagnostics, not data-fit probabilities",
        "number_of_realizations_per_model": number,
        "production_target_per_model": 10000,
        "full_production_sample": number >= 10000,
        "grid_gev": list(GRID), "history_myr": 200.,
        "running_slope_window_bins": 5, "running_slope_window_dex": .25,
        "fit_weighting": "equal weights in log flux on a logarithmic energy grid",
        "selection": "consecutive seeds from zero, no feature selection",
        "benchmark_caveat": "historical observed amplitudes; unmatched estimators, not observed values of these statistics",
        "intervals": "pointwise exact two-sided 95% binomial; zero counts also have one-sided 95% limits",
        "windows": {key: {"fit_range_gev": list(fit), "search_range_gev": list(search)}
                    for key, (fit, search) in WINDOWS.items()},
        "source_sha256": {p.name: digest(p) for p in
                          (Path(__file__), Path(__file__).with_name("feature_statistics.py"),
                           Path(__file__).with_name("plot_halos.py"),
                           Path(__file__).with_name("plot_variations.py"), ROOT / "scripts/run_results.py")},
        "models": {},
    }
    values, spectra = {}, {}
    libraries = executable = reference_params = None
    for name in NAMES:
        directory = runs / name
        associated = name.startswith("associations_")
        energy, flux, provenance, hashes = load_model(
            directory, number, manifest_name(name), ("population",) if associated else (), grid=GRID)
        validate_parameters(provenance["config_text"], expected_parameters(name, "wide"))
        validate_parameters((directory / "params.ini").read_text(), expected_parameters(name, "wide"))
        params = gp.load_params(directory)
        # All resolved settings except the deliberate variations must agree.
        ignored = {"simname", "seed", "hkpc", "snrateyr", "efficiency", "varyenergy", "varyslope",
                   "injslopesigma", "syntheticassociations", "associationfraction"}
        shared = {k: v for k, v in params.items() if k not in ignored}
        if reference_params is None:
            reference_params = shared
        elif shared != reference_params:
            raise ValueError(f"{name}: unplanned difference in resolved parameters")
        if libraries is None:
            libraries = provenance["library_sha256"]
            executable = provenance["executable_sha256"]
        if libraries != provenance["library_sha256"]:
            raise ValueError("all nine models must use the same physics library")
        if not associated and executable != provenance["executable_sha256"]:
            raise ValueError("independent models must use the same executable")
        data = {"parameters": params, "provenance": provenance,
                "flux_sha256_by_seed": hashes, "windows": {}}
        if associated:
            counts, population_hashes = [], []
            for seed in range(number):
                path = directory / f"population_{seed:06d}.txt"
                table = gp.load(path)
                if table.meta["simname"] != name or int(table.meta["seed"]) != seed or len(table) != 1:
                    raise ValueError("population identity/row count mismatch")
                row = [float(table[k][0]) for k in
                       ("parents", "field_SN", "clustered_SN", "retained_SN", "expected_SN")]
                if (not np.all(np.isfinite(row)) or min(row) < 0 or row[-1] != 4e6 or
                        row[3] > row[1] + row[2]):
                    raise ValueError("invalid 200 Myr association population counts")
                counts.append(row)
                population_hashes.append(digest(path))
            counts = np.asarray(counts)
            total = counts[:, 1] + counts[:, 2]
            data["population"] = {"expected_explosions": 4e6,
                                  "mean_explosions": float(total.mean()),
                                  "standard_error": float(total.std(ddof=1) / np.sqrt(number)),
                                  "sha256_by_seed": population_hashes}
        values[name] = {}
        for key, (fit, search) in WINDOWS.items():
            a, d = excursions(energy, flux, fit_range=fit, search_range=search)
            values[name][key] = (a, d)
            result = statistic_summary(a, d)
            result["sample_size_checks"] = {
                str(n): statistic_summary(a[:n], d[:n])
                for n in sorted({max(1, number // 4), max(1, number // 2), number})}
            ma, md = excursions(energy, np.median(flux, axis=0), fit_range=fit, search_range=search)
            result["median_spectrum_statistics"] = {"residual": float(ma[0]), "slope_excursion": float(md[0])}
            data["windows"][key] = result
        # A nested search must never reduce either statistic when the fit is fixed.
        for broad, narrow in zip(values[name]["primary_1_100_TeV"], values[name]["secondary_10_100_TeV"]):
            if np.any(broad + 1e-12 < narrow):
                raise ValueError("nested-window invariant failed")
        centres, indices = running_indices(energy, flux)
        data["energy_gev"] = energy.tolist()
        data["flux_percentiles"] = {str(q): np.percentile(flux, q, axis=0).tolist() for q in (5, 50, 95)}
        data["index_energy_gev"] = centres.tolist()
        data["index_percentiles"] = {str(q): np.percentile(indices, q, axis=0).tolist() for q in (5, 50, 95)}
        summary["models"][name] = data
        spectra[name] = (energy, flux, centres, indices)
    summary["continuous_index_limit"] = {}
    for key, (fit, search) in WINDOWS.items():
        a, d = excursions(energy, continuous_index_spectrum(energy), fit_range=fit, search_range=search)
        summary["continuous_index_limit"][key] = {"residual": float(a[0]), "slope_excursion": float(d[0]),
                                                  "subtracted_from_realizations": False}
    return summary, values, spectra


def plot_survival(values, group, window, number, output):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(16, 7.7))
    xmax = max(.8, 1.1 * max(float(v.max()) for name, _, _ in GROUPS[group]
                           for v in values[name][window]))
    for name, label, color in GROUPS[group]:
        for ax, v in zip(axes, values[name][window]):
            x = np.unique(np.r_[1e-3, v[v >= 1e-3], xmax, *THRESHOLDS[0], *THRESHOLDS[1]])
            empirical, lo, hi = survival(v, x)
            ax.fill_between(x, np.maximum(lo, 1e-9), hi, step="pre", color=color, alpha=.12, linewidth=0)
            ax.step(x, empirical, where="pre", color=color, label=label)
    for ax, band, label in zip(axes, THRESHOLDS, (r"Maximum residual $A$", r"Running-slope excursion $\Delta\gamma_{\rm run}$")):
        ax.axvspan(*band, color=".7", alpha=.3, label="Amplitude benchmark")
        ax.axhline(1 / number, linestyle=":", color=".4", linewidth=1.2)
        ax.set(xscale="log", yscale="log", xlim=(1e-3, xmax), ylim=(.25 / number, 1.05),
               xlabel=label, ylabel=r"Survival probability $P(\geq x)$")
        ax.grid(False, which="both")
    axes[0].legend(loc="lower left", fontsize=17, frameon=True, facecolor="white", framealpha=1)
    fit, search = WINDOWS[window]
    detail = f"{number:,} realizations/model\nSearch: {search[0]/1e3:g}--{search[1]/1e3:g} TeV\n"
    detail += f"Reference fit: {fit[0]/1e3:g}--{fit[1]/1e3:g} TeV"
    if number < 10000:
        detail += "\nPilot sample"
    axes[1].text(.97, .96, detail, transform=axes[1].transAxes, ha="right", va="top", fontsize=16)
    for ext in ("pdf", "png"):
        fig.savefig(output / f"gryphon_{group}_statistics.{ext}", dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_spectra(spectra, number, output):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(16, 7.7))
    for name, label, color in GROUPS["variations"]:
        energy, flux, centres, indices = spectra[name]
        pivot = int(np.argmin(abs(energy - 1e3)))
        if not np.isclose(energy[pivot], 1e3):
            raise ValueError("wide grid must include the 1 TeV normalization pivot")
        factor = 1 / (np.median(flux[:, pivot]) * energy[pivot]**2.7)
        for ax, x, y in ((axes[0], energy, flux * energy**2.7 * factor), (axes[1], centres, indices)):
            lo, median, hi = np.percentile(y, [5, 50, 95], axis=0)
            ax.fill_between(x, lo, hi, color=color, alpha=.15, linewidth=0)
            ax.plot(x, median, color=color, label=label)
    for ax in axes:
        ax.set(xscale="log", xlim=(100, 1e6), xlabel="$E$ [GeV]")
        ax.grid(False, which="both")
        ax.axvspan(1e3, 1e5, color=".5", alpha=.06)
    axes[0].set_ylabel(r"Normalized $E^{2.7}J(E)$")
    axes[1].set_ylabel(r"Running index $\gamma_{\rm run}$")
    axes[0].legend(fontsize=18)
    axes[1].text(.03, .03, f"{number:,} realizations/model\nMedian and central 90\\%" +
                 ("\nPilot sample" if number < 10000 else ""), transform=axes[1].transAxes, fontsize=17)
    for ext in ("pdf", "png"):
        fig.savefig(output / f"gryphon_variations_spectra.{ext}", dpi=140, bbox_inches="tight")
    plt.close(fig)


def write_results(summary, values, spectra, output):
    output.mkdir(parents=True, exist_ok=True)
    for window in WINDOWS:
        folder = output / window
        folder.mkdir(exist_ok=True)
        rows = []
        with (folder / "statistics.csv").open("w") as stream:
            writer = csv.writer(stream)
            writer.writerow(("model", "seed", "residual", "slope_excursion"))
            for name in NAMES:
                writer.writerows((name, seed, a, d) for seed, (a, d) in enumerate(zip(*values[name][window])))
        for name, data in summary["models"].items():
            stats = data["windows"][window]
            a, d = stats["residual"], stats["slope_excursion"]
            rows.append((name, a["percentiles"]["95"], d["percentiles"]["95"],
                         *[a["thresholds"][str(t)]["exceedances"] for t in THRESHOLDS[0]],
                         *[d["thresholds"][str(t)]["exceedances"] for t in THRESHOLDS[1]]))
        with (folder / "summary.csv").open("w") as stream:
            writer = csv.writer(stream)
            writer.writerow(("model", "A95", "slope95", "n_A16", "n_A19", "n_slope20", "n_slope35"))
            writer.writerows(rows)
        tex = [r"\begin{tabular}{lrrrrrr}", r"\hline", r"Model & $A_{95}$ & $\Delta\gamma_{95}$ & $n_{A\geq.16}$ & $n_{A\geq.19}$ & $n_{\Delta\gamma\geq.20}$ & $n_{\Delta\gamma\geq.35}$ \\", r"\hline"]
        for name, a, d, *counts in rows:
            tex.append(name.replace("_", r"\_") + f" & {100*a:.2f}\\% & {d:.4f} & " +
                       " & ".join(map(str, counts)) + r" \\")
        tex += [r"\hline", r"\end{tabular}"]
        (folder / "summary_table.tex").write_text("\n".join(tex) + "\n")
        for group in GROUPS:
            plot_survival(values, group, window, summary["number_of_realizations_per_model"], folder)
    plot_spectra(spectra, summary["number_of_realizations_per_model"], output)
    atomic_json(output / "summary.json", summary)
    print(f"Verified and analysed {summary['number_of_realizations_per_model']} seeds in all nine models: {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=Path, default=ROOT / "runs/hebreaks_wide_10k")
    parser.add_argument("--seeds", type=int, default=10000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summary, values, spectra = analyse(args.runs.resolve(), args.seeds)
    gp.style.use()
    write_results(summary, values, spectra, args.output.resolve())


if __name__ == "__main__":
    main()
