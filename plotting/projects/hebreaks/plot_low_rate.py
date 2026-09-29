"""Restyle the historical low-rate comparison without rerunning simulations."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import re
import tarfile

import numpy as np

import gryphon_plots as gp
from feature_statistics import binomial_interval, excursions, summarize

THRESHOLDS = ((0.16, 0.19), (0.20, 0.35))
MODELS = (("baseline", "tab:orange", r"$\mathcal{R}=2\ \mathrm{century}^{-1}$"),
          ("rare", "tab:blue", r"$\mathcal{R}_{\rm HE}=0.2\ \mathrm{century}^{-1}$"))


def load_legacy(summary: dict):
    """Read the original ensembles, never substitute the current H=4 sample."""
    spectra = {}
    energy = None

    def add(rows, seed, raw, name):
        nonlocal energy
        if seed in rows:
            raise ValueError(f"duplicate realization: {name}")
        table = np.loadtxt(io.BytesIO(raw))
        if table.shape != (48, 2):
            raise ValueError(f"{name}: expected 48 energies and two columns")
        if energy is None:
            energy = table[:, 0]
        elif not np.array_equal(table[:, 0], energy):
            raise ValueError(f"{name}: inconsistent energy grid")
        rows[seed] = table[:, 1]

    rows = {}
    archive = Path(summary["baseline_archive"]).expanduser()
    pattern = re.compile(r"done/test_fixed4_(\d+)\.txt")
    # Stream the archive without extracting or modifying any source files.
    with tarfile.open(archive, "r|gz") as stream:
        for member in stream:
            match = pattern.fullmatch(member.name)
            if member.isfile() and match:
                add(rows, int(match[1]), stream.extractfile(member).read(), member.name)
    size = summary["residual"]["baseline"]["number_of_realizations"]
    if set(rows) != set(range(size)):
        raise ValueError("historical baseline must contain exactly seeds 0,...,N-1")
    spectra["baseline"] = np.stack([rows[i] for i in range(size)])

    rows = {}
    directory = Path(summary["rare_directory"]).expanduser()
    for path in directory.glob("test_raresn_*.txt"):
        seed = int(path.stem.rsplit("_", 1)[1])
        add(rows, seed, path.read_bytes(), str(path))
    size = summary["residual"]["rare"]["number_of_realizations"]
    if set(rows) != set(range(1, size + 1)):
        raise ValueError("historical low-rate sample must contain exactly seeds 1,...,N")
    spectra["rare"] = np.stack([rows[i] for i in range(1, size + 1)])
    return energy, spectra


def check_legacy(values: np.ndarray, expected: dict) -> None:
    """A styling change must reproduce the archived numerical results."""
    if len(values) != expected["number_of_realizations"]:
        raise ValueError("historical sample size changed")
    for q, value in expected["percentiles"].items():
        if not np.isclose(np.percentile(values, float(q)), value, rtol=1e-8, atol=1e-12):
            raise ValueError(f"historical percentile {q} changed")
    if not np.isclose(values.max(), expected["maximum"], rtol=1e-8, atol=1e-12):
        raise ValueError("historical maximum changed")
    for threshold, record in expected["thresholds"].items():
        if np.count_nonzero(values >= float(threshold)) != record["exceedances"]:
            raise ValueError(f"historical exceedance count at {threshold} changed")


def survival(values: np.ndarray, thresholds: np.ndarray):
    """P(X >= x), with memory use linear rather than quadratic in sample size."""
    ordered = np.sort(values)
    counts = len(ordered) - np.searchsorted(ordered, thresholds, side="left")
    lower, upper = binomial_interval(counts, len(ordered))
    empirical = np.where(counts > 0, counts / len(ordered), np.nan)
    return empirical, lower, upper


def main() -> None:
    proj = gp.project()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path,
                        default=proj.paper.parent / "modelstats/low_rate_summary.json")
    args = parser.parse_args()
    historical = json.loads(args.summary.read_text())
    print("Reading the historical baseline and low-rate spectra...", flush=True)
    energy, spectra = load_legacy(historical)
    values, results = {}, {}
    for model, _, _ in MODELS:
        values[model] = excursions(energy, spectra[model])
        results[model] = {}
        for key, data, band in zip(("residual", "slope_excursion"), values[model], THRESHOLDS):
            check_legacy(data, historical[key][model])
            results[model][key] = summarize(data, band)
    print("All archived percentiles, maxima and exceedance counts reproduced.", flush=True)

    gp.style.use()
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(16, 7.7))
    xmax = max(0.8, 1.1 * max(v.max() for pair in values.values() for v in pair))
    nmax = max(len(pair[0]) for pair in values.values())
    for model, color, label in MODELS:
        for ax, data, band in zip(axes, values[model], THRESHOLDS):
            # Resolve every tail step; a dense log grid suffices for the body.
            tail = np.sort(data)[-200:]
            x = np.unique(np.r_[np.geomspace(1e-3, xmax, 600), tail, band])
            empirical, lower, upper = survival(data, x)
            ax.fill_between(x, np.maximum(lower, 1e-9), upper, step="pre",
                            color=color, alpha=0.10, linewidth=0, zorder=1)
            ax.step(x, empirical, where="pre", color=color, label=label, zorder=3)
            ax.axhline(1 / len(data), color=color, linestyle=":", linewidth=1.2,
                       alpha=0.8, zorder=2)
    for ax, band, xlabel in zip(axes, THRESHOLDS,
                               (r"Maximum residual $A$", r"Running-slope excursion $\Delta\gamma_{\rm run}$")):
        ax.axvspan(*band, color="0.65", alpha=0.30, zorder=0, label="Data-motivated range")
        ax.set(xscale="log", yscale="log", xlim=(1e-3, xmax), ylim=(0.3 / nmax, 1.05),
               xlabel=xlabel, ylabel=r"Survival probability $P(\geq x)$")
        ax.grid(False, which="both")
        ax.tick_params(labelsize=23)
    axes[0].legend(fontsize=18, loc="center left", bbox_to_anchor=(0.01, 0.49),
                   labelspacing=0.35)
    nref, nrare = (len(values[name][0]) for name in ("baseline", "rare"))
    axes[1].text(0.97, 0.95, "Legacy ensembles\n" +
                 f"Reference: {nref:,}\nLow-rate: {nrare:,}\n" +
                 "$H=4$ kpc\n" + r"$10<E<100$ TeV",
                 transform=axes[1].transAxes, va="top", ha="right", fontsize=16)
    proj.figures.mkdir(parents=True, exist_ok=True)
    fig.savefig(proj.figures / "gryphon_low_rate_statistics.pdf", bbox_inches="tight")
    fig.savefig(proj.figures / "gryphon_low_rate_statistics.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    audit = {
        "status": "restyled historical ensembles; no new simulations",
        "historical_summary": str(args.summary.resolve()),
        "historical_summary_sha256": hashlib.sha256(args.summary.read_bytes()).hexdigest(),
        "baseline_archive": historical["baseline_archive"],
        "rare_directory": historical["rare_directory"],
        "verification": "all archived percentiles, maxima, sizes and exceedance counts reproduced",
        "intervals": "pointwise exact two-sided 95% Clopper-Pearson; not simultaneous",
        "models": results,
    }
    (proj.figures / "low_rate_replot_summary.json").write_text(json.dumps(audit, indent=2) + "\n")


if __name__ == "__main__":
    main()
