"""Current-model low-rate and association ensembles, with verified provenance."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import gryphon_plots as gp
from feature_statistics import excursions, summarize
from plot_halos import THRESHOLDS, digest, load_model, survival, validate_parameters as validate_halo

GROUPS = {
    "low_rate": (("reference", "halos_H4", "halos_manifest.json", "tab:orange",
                  r"$\mathcal{R}=2\ \mathrm{century}^{-1}$"),
                 ("rare", "low_rate", "low_rate_manifest.json", "tab:blue",
                  r"$\mathcal{R}_{\rm HE}=0.2\ \mathrm{century}^{-1}$")),
    "association": tuple((name, f"associations_{name}", "benchmark_manifest.json", color,
                           rf"{name.capitalize()}: $f_{{\rm a}}={fraction:g}$")
                         for name, color, fraction in (("independent", "tab:blue", 0),
                                                      ("mixed", "tab:orange", .5),
                                                      ("clustered", "tab:red", 1))),
}


def validate_parameters(params, group, name):
    reference = dict(params)
    if group == "low_rate" and name == "rare":
        if params.get("snrateyr") != .002 or params.get("efficiency") != 1.:
            raise ValueError("low-rate model requires 0.002/year and efficiency=1")
        reference.update(snrateyr=.02, efficiency=.1)
    if group == "association":
        expected = {"syntheticassociations": "true", "associationmembers": 100.,
                    "associationradiuspc": 30., "associationvelocitykms": 3.,
                    "associationfraction": {"independent": 0., "mixed": .5, "clustered": 1.}[name]}
        for key, value in expected.items():
            if params.get(key) != value:
                raise ValueError(f"{name}: incorrect {key}")
        reference["syntheticassociations"] = "false"
    validate_halo(reference, 4)


def shared_parameters(params, group):
    ignored = {"simname", "seed"}
    ignored.update(("snrateyr", "efficiency") if group == "low_rate" else ("associationfraction",))
    return {k: v for k, v in params.items() if k not in ignored}


def population_data(directory, number):
    rows, hashes = [], []
    for seed in range(number):
        path = directory / f"population_{seed:06d}.txt"
        table = gp.load(path)
        if int(table.meta["seed"]) != seed or table.meta["simname"] != directory.name:
            raise ValueError("population seed/model mismatch")
        if len(table) != 1:
            raise ValueError("population table requires exactly one row")
        rows.append([table[key][0] for key in
                     ("parents", "field_SN", "clustered_SN", "retained_SN", "expected_SN")])
        hashes.append(digest(path))
    population = np.asarray(rows)
    if (not np.all(np.isfinite(population)) or np.any(population < 0) or
            not np.all(population[:, 4] == 2e6) or
            np.any(population[:, 3] > population[:, 1] + population[:, 2])):
        raise ValueError("invalid population normalization/counts")
    count = population[:, 1] + population[:, 2]
    return {"mean_explosions": float(count.mean()), "expected_explosions": 2e6,
            "standard_error_mean_explosions": float(count.std(ddof=1)/np.sqrt(number)),
            "mean_clustered_fraction": float(population[:, 2].sum()/count.sum()),
            "mean_retained_explosions": float(population[:, 3].mean())}, hashes


def analyse(proj, group, number):
    summary = {"status": "current cutoff-free model; not a data-fit p-value",
               "number_of_realizations_per_model": number,
               "fit_range_gev": [1e3, 1e6], "search_range_gev": [1e4, 1e5],
               "running_slope_window_bins": 5,
               "intervals": "pointwise exact two-sided 95% Clopper-Pearson; not simultaneous",
               "zero_count_limits": "one-sided 95%: 1 - 0.05**(1/N)",
               "selection": "consecutive seeds 0,...,N-1, not selected for spectral features",
               "analysis_source_sha256": {name: digest(Path(__file__).with_name(name)) for name in
                                           ("plot_population_results.py", "plot_halos.py", "feature_statistics.py")},
               "models": {}}
    reference = reference_grid = reference_build = None
    values, rows = {}, []
    for name, run, manifest, color, label in GROUPS[group]:
        directory = proj.run(run)
        extra = ("population",) if group == "association" else ()
        energy, flux, provenance, hashes = load_model(directory, number, manifest, extra)
        params = gp.load_params(directory)
        validate_parameters(params, group, name)
        matched = shared_parameters(params, group)
        build = {k: provenance[k] for k in ("executable_sha256", "library_sha256")}
        if reference is None:
            reference, reference_grid, reference_build = matched, energy, build
        elif matched != reference or not np.array_equal(energy, reference_grid) or build != reference_build:
            raise ValueError("population comparison requires matched physics, grids, and builds")
        residual, slope = excursions(energy, flux)
        data = {"seeds": list(range(number)), "parameters": params, "provenance": provenance,
                "flux_sha256_by_seed": hashes,
                "residual": summarize(residual, THRESHOLDS[0]),
                "slope_excursion": summarize(slope, THRESHOLDS[1]),
                "sample_size_checks": {}}
        for n in sorted({max(1, number//4), max(1, number//2), number}):
            data["sample_size_checks"][str(n)] = {
                "residual": summarize(residual[:n], THRESHOLDS[0]),
                "slope_excursion": summarize(slope[:n], THRESHOLDS[1])}
        if group == "association":
            data["population"], data["population_sha256_by_seed"] = population_data(directory, number)
        summary["models"][name] = data
        values[name] = (residual, slope)
        rows.extend((name, seed, a, d) for seed, (a, d) in enumerate(zip(residual, slope)))
    # Association executable differs, but propagation/injection library must
    # still be the same build as the common H=4 reference used in A/B/D.
    halo_manifest = json.loads((proj.run("halos_H4") / "halos_manifest.json").read_text())
    if reference_build["library_sha256"] != halo_manifest["library_sha256"]:
        raise ValueError("population and halo ensembles use different physics libraries")
    return values, summary, rows


def plot(values, group, number, output):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(16, 7.7))
    xmax = max(.8, 1.1*max(v.max() for pair in values.values() for v in pair))
    for name, _, _, color, label in GROUPS[group]:
        for ax, data in zip(axes, values[name]):
            x = np.unique(np.r_[1e-3, data[data >= 1e-3], xmax, *THRESHOLDS[0], *THRESHOLDS[1]])
            empirical, lower, upper = survival(data, x)
            ax.fill_between(x, np.maximum(lower, 1e-9), upper, step="pre",
                            color=color, alpha=.10, linewidth=0, zorder=1)
            ax.step(x, empirical, where="pre", color=color, label=label, zorder=3)
    for ax, band, xlabel in zip(axes, THRESHOLDS,
                               (r"Maximum residual $A$", r"Running-slope excursion $\Delta\gamma_{\rm run}$")):
        ax.axvspan(*band, color=".65", alpha=.30, zorder=0, label="Data range")
        ax.axhline(1/number, color=".45", linestyle=":", linewidth=1.2, zorder=2)
        ax.set(xscale="log", yscale="log", xlim=(1e-3, xmax), ylim=(.25/number, 1.05),
               xlabel=xlabel, ylabel=r"Survival probability $P(\geq x)$")
        ax.grid(False, which="both")
        ax.tick_params(labelsize=23)
    axes[0].legend(fontsize=18, loc="lower left", frameon=True,
                   facecolor="white", edgecolor="none", framealpha=1)
    detail = "$H=4$ kpc, Xie" + (", $N_{\\rm a}=100$" if group == "association" else "")
    axes[1].text(.97, .95, f"{number:,} realizations/model\n" + detail + "\n" +
                 r"$10<E<100$ TeV", transform=axes[1].transAxes, ha="right", va="top", fontsize=16)
    for extension in ("pdf", "png"):
        fig.savefig(output / f"gryphon_{group}_statistics.{extension}", bbox_inches="tight", dpi=140)
    plt.close(fig)


def write_table(path, group, models):
    rows = ["% Generated by plot_population_results.py.", r"\begin{tabular}{lcccccc}",
            r"\hline\hline", r"Population & $A_{95}$ & $(\Delta\gamma_{\rm run})_{95}$ & "
            r"$n_{A\geq0.16}$ & $n_{A\geq0.19}$ & $n_{\Delta\gamma\geq0.20}$ & $n_{\Delta\gamma\geq0.35}$ \\",
            r"\hline"]
    labels = {"reference": "Reference", "rare": "Low rate", "independent": "$f_{\\rm a}=0$",
              "mixed": "$f_{\\rm a}=0.5$", "clustered": "$f_{\\rm a}=1$"}
    for name in models:
        a, d = (models[name][key] for key in ("residual", "slope_excursion"))
        counts = [a["thresholds"][str(t)]["exceedances"] for t in THRESHOLDS[0]]
        counts += [d["thresholds"][str(t)]["exceedances"] for t in THRESHOLDS[1]]
        rows.append(f"{labels[name]} & {100*a['percentiles']['95']:.2f}\\% & {d['percentiles']['95']:.4f} & " +
                    " & ".join(map(str, counts)) + r" \\")
    rows.extend((r"\hline\hline", r"\end{tabular}"))
    path.write_text("\n".join(rows) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, default=10000)
    args = parser.parse_args()
    gp.style.use()
    proj = gp.project()
    proj.figures.mkdir(parents=True, exist_ok=True)
    for group in GROUPS:
        values, summary, rows = analyse(proj, group, args.seeds)
        plot(values, group, args.seeds, proj.figures)
        (proj.figures / f"{group}_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        with (proj.figures / f"{group}_statistics.txt").open("w") as out:
            out.write("# model seed residual slope_excursion\n")
            for name, seed, a, d in rows:
                out.write(f"{name} {seed} {a:.10e} {d:.10e}\n")
        write_table(proj.figures / f"{group}_table.tex", group, summary["models"])
        print(json.dumps({name: {k: v for k, v in m.items() if k in
                                 ("residual", "slope_excursion", "population")}
                          for name, m in summary["models"].items()}, indent=2), flush=True)


if __name__ == "__main__":
    main()
