"""A single, unselected matched-catalogue illustration of Eq. (13)."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import gryphon_plots as gp

from plot_halos import digest
from plot_variations import MODELS, LINESTYLES, validate_parameters, shared_parameters

ROOT = Path(__file__).resolve().parents[3]


def validate_table(table, original):
    if int(table.meta["seed"]) != 0:
        raise ValueError("illustration requires a priori seed 0")
    energy, flux, neff, largest = table.columns("E", "I", "N_eff", "q_max")
    count = int(table.meta["sourcecount"])
    if len(energy) != 48 or not np.allclose(energy, np.geomspace(1e3, 1e6, 48), rtol=1e-12):
        raise ValueError("unexpected diagnostic energy grid")
    if (not np.all(np.isfinite(table.values)) or np.any(flux <= 0) or
            np.any(neff < 1.-1e-12) or np.any(neff > count*(1+1e-12)) or
            np.any(largest <= 0) or np.any(largest > 1+1e-12)):
        raise ValueError("invalid participation weights/counts")
    # sum q_i^2 lies between q_max^2 and q_max for positive weights.
    if np.any(1/neff < largest**2*(1-1e-10)) or np.any(1/neff > largest*(1+1e-10)):
        raise ValueError("participation number inconsistent with largest source weight")
    if (not np.allclose(energy, original["E"], rtol=1e-8) or
            not np.allclose(flux, original["I"], rtol=1e-8, atol=0)):
        raise ValueError("diagnostic does not reproduce the original seed-0 spectrum")
    return float(np.max(np.abs(flux / original["I"] - 1)))


def analyse(root):
    curves, models = {}, {}
    catalogues, counts, builds, params_shared = [], [], [], []
    for name, (run, _, _) in MODELS.items():
        directory = root / run
        manifest = json.loads((directory / "participation_manifest.json").read_text())
        if manifest["seed"] != 0 or digest(directory / "participation_000000.txt") != json.loads(
                (directory / "complete_000000.json").read_text())["participation_sha256"]:
            raise ValueError("changed seed/diagnostic output")
        if hashlib.sha256(manifest["config_text"].encode()).hexdigest() != manifest["config_sha256"]:
            raise ValueError("changed diagnostic config manifest")
        reference = Path(manifest["reference_directory"])
        original = reference / "flux_000000.txt"
        if digest(original) != manifest["reference_flux_sha256"]:
            raise ValueError("original reference spectrum changed")
        table = gp.load(directory / "participation_000000.txt")
        if table.meta["simname"] != run:
            raise ValueError("diagnostic model header mismatch")
        error = validate_table(table, gp.load(original))
        params = gp.load_params(directory)
        validate_parameters(params, name)
        if shared_parameters(params) != shared_parameters(gp.load_params(reference)):
            raise ValueError("diagnostic and original physics differ")
        params_shared.append(shared_parameters(params))
        builds.append({k: manifest[k] for k in ("executable_sha256", "library_sha256")})
        catalogues.append(table.meta["cataloguefnv1a64"])
        counts.append(int(table.meta["sourcecount"]))
        curves[name] = table["N_eff"]
        models[name] = {"provenance": manifest, "parameters": params,
                        "participation_sha256": digest(table.path),
                        "maximum_relative_flux_difference_from_overnight": error,
                        "N_eff": table["N_eff"].tolist(), "q_max": table["q_max"].tolist()}
    for values in (catalogues, counts, builds, params_shared):
        if any(value != values[0] for value in values):
            raise ValueError("comparison requires identical catalogues and matched physics/builds")
    energy = table["E"]
    for name, neff in curves.items():
        models[name]["ratio_to_identical"] = (neff/curves["identical"]).tolist()
        models[name]["N_eff_at_1_10_100_1000_TeV_log_interpolated"] = np.exp(
            np.interp(np.log([1e3, 1e4, 1e5, 1e6]), np.log(energy), np.log(neff))).tolist()
    summary = {"seed": 0, "selection": "single fixed catalogue, not selected by spectral features",
               "source_count": counts[0], "catalogue_fnv1a64": catalogues[0],
               "energy_gev": energy.tolist(), "models": models,
               "interpretation": "one-realization illustration, not ensemble quantiles or feature probabilities",
               "analysis_source_sha256": digest(Path(__file__))}
    return energy, curves, summary


def plot(energy, curves, output):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FuncFormatter
    fig, axes = plt.subplots(1, 2, figsize=(16, 7.7))
    for name, (_, label, color) in MODELS.items():
        axes[0].plot(energy, curves[name], color=color, ls=LINESTYLES[name], label=label)
        axes[1].plot(energy, curves[name]/curves["identical"], color=color,
                     ls=LINESTYLES[name], label=label)
    for ax in axes:
        ax.axvspan(1e4, 1e5, color=".65", alpha=.18, zorder=0)
        ax.set(xscale="log", yscale="log", xlim=(1e3, 1e6), xlabel=r"$E$ [GeV]")
        ax.grid(False, which="both")
        ax.tick_params(labelsize=23)
    axes[0].set_ylabel(r"Effective source number $N_{\rm eff}$")
    axes[1].set_ylabel(r"$N_{\rm eff}/N_{\rm eff,identical}$")
    axes[0].legend(loc="best", fontsize=19, frameon=True, facecolor="white", edgecolor="none", framealpha=1)
    axes[1].yaxis.set_major_locator(FixedLocator([.1, .2, .5, 1, 2, 5, 10]))
    axes[1].yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    axes[1].text(.96, .96, "$H=4$ kpc, Xie\nSame catalogue, seed 0\nOne realization",
                 ha="right", va="top", transform=axes[1].transAxes, fontsize=18,
                 bbox=dict(facecolor="white", edgecolor="none", pad=3))
    for extension in ("pdf", "png"):
        fig.savefig(output / f"gryphon_effective_sources.{extension}", bbox_inches="tight", dpi=140)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=Path, default=ROOT / "runs/hebreaks_effective_sources")
    args = parser.parse_args()
    gp.style.use()
    output = gp.project().figures
    output.mkdir(parents=True, exist_ok=True)
    energy, curves, summary = analyse(args.runs)
    plot(energy, curves, output)
    (output / "effective_sources_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    for name, model in summary["models"].items():
        print(name, model["N_eff_at_1_10_100_1000_TeV_log_interpolated"],
              "flux check", model["maximum_relative_flux_difference_from_overnight"])


if __name__ == "__main__":
    main()
