#!/usr/bin/env python3
"""Plot saved runMSP tables. No source sampling, injection or propagation in Python."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]


def read_pairs(path):
    result = {}
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            key, value = line.split("=", 1)
            result[key.strip()] = value.strip()
    return result


def read_ams(path):
    data = np.loadtxt(path)
    r, flux = data[:, 0], data[:, 1]
    mass = 0.00051099895  # GeV
    energy = np.sqrt(r*r + mass*mass) - mass
    jacobian = (energy + mass)/r
    errors = np.array([np.hypot(data[:, 2], data[:, 4]), np.hypot(data[:, 3], data[:, 5])])
    return energy, flux*jacobian, errors*jacobian


def plot(energy, spectra, population, horizon, params, data, dest):
    plt.rcParams.update({"font.family": "serif", "font.size": 14,
                         "axes.labelsize": 16, "axes.grid": False,
                         "axes.linewidth": 1.2, "lines.linewidth": 2,
                         "xtick.direction": "in", "ytick.direction": "in",
                         "xtick.top": True, "ytick.right": True,
                         "mathtext.fontset": "dejavuserif", "pdf.fonttype": 42})
    blue, orange = "#00246a", "#d97900"
    low, median, high = np.quantile(spectra, [.16, .5, .84], axis=0)
    eta = float(params["efficiency"])
    fig, ax = plt.subplots(figsize=(8.7, 5.2), layout="constrained")
    ax.fill_between(energy, energy**3*low, energy**3*high, color=blue, alpha=.2,
                    label="16–84% synthetic-population range")
    ax.loglog(energy, energy**3*median, color=blue, label=rf"Jelly MSP median, $\eta_{{\rm pair}}={eta:g}$")
    ax.loglog(energy, energy**3*spectra[0], color=orange, lw=1.4, label="First realization (fixed seed)")
    ax.loglog(energy, energy**3*median*.1, "--", color=blue,
              label=rf"Same median $\times\,0.1$ ($\eta_{{\rm pair}}={.1*eta:g}$)")
    e, f, err = data
    keep = (e >= energy.min()) & (e <= energy.max())
    ax.errorbar(e[keep], (e**3*f)[keep], yerr=err[:, keep]*e[keep]**3,
                fmt="o", ms=3.4, color="black", elinewidth=1, capsize=2,
                label="AMS-02: total observed positrons", zorder=5)
    ax.set(xlim=(energy.min(), energy.max()), ylim=(.3, max(150, 1.6*np.max(energy**3*high))),
           xlabel=r"Positron energy $E$ [GeV]",
           ylabel=r"$E^3 J_{e^+}$ [GeV$^2$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]")
    ax.legend(fontsize=10, loc="lower left", frameon=False)
    fig.savefig(dest/"msp_synthetic_spectrum.pdf")
    plt.close(fig)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4.5), layout="constrained")
    rsun, radius = float(params["sunkpc"]), float(params["rgkpc"])
    axs[0].scatter(population[:, 0]+rsun, population[:, 1], s=.5, alpha=.2, color=blue, rasterized=True)
    axs[0].plot(rsun, 0, "*", color=orange, ms=14, label="Sun")
    axs[0].set(xlim=(-radius, radius), ylim=(-radius, radius), aspect="equal", xlabel="$x$ [kpc]", ylabel="$y$ [kpc]")
    axs[0].legend(frameon=False)
    bins = np.linspace(26, 36, 60)
    axs[1].hist(np.log10(population[:, 6]), bins=bins,
                weights=np.full(len(population), 1/(len(population)*np.diff(bins)[0])),
                color=blue, histtype="step", lw=2)
    axs[1].set(xlabel=r"$\log_{10}(L_{\rm sd}/{\rm erg\,s^{-1}})$", ylabel="Probability density [dex$^{-1}$]")
    fig.savefig(dest/"msp_synthetic_population.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.7, 4.8), layout="constrained")
    for e in np.unique(horizon[:, 0]):
        rows = horizon[np.isclose(horizon[:, 0], e)]
        ax.plot(rows[:, 1], rows[:, 2], label=f"{e:g} GeV")
    ax.set(xlim=(0, 12), ylim=(0, 1.02), xlabel="Distance from the Sun [kpc]",
           ylabel="Cumulative fraction of the MSP flux")
    ax.legend(frameon=False, title="First realization")
    fig.savefig(dest/"msp_synthetic_horizon.pdf")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=Path, default=ROOT/"runs/msp_ensemble")
    parser.add_argument("--output", type=Path, default=Path(__file__).parent/"output")
    parser.add_argument("--data", type=Path, default=ROOT.parent/"KISS-CosmicRayDataBase/kiss_tables/AMS-02_e+_rigidity.txt")
    parser.add_argument("--publish", type=Path, help="msp.slides root; copy plots and generated numbers there")
    args = parser.parse_args()
    manifest = json.loads((args.runs/"manifest.json").read_text())
    if manifest["engine"] != "gryphon2/runMSP" or manifest["status"] != "complete":
        parser.error("Need a completed native runMSP ensemble")
    if args.publish and not (args.publish/"synthetic_msp.tex").is_file():
        parser.error("--publish must point to the existing msp.slides deck")
    fluxes, totals, local, outside, checks, table_hashes = [], [], [], [], [], []
    common, energy = None, None
    for seed in manifest["seeds"]:
        directory = args.runs/f"seed_{seed:06d}"/manifest["model"]
        suffix = f"_{seed:06d}.txt"
        params = read_pairs(directory/f"params{suffix}")
        comparison = {k: v for k, v in params.items() if k != "seed"}
        if common is not None and common != comparison:
            parser.error(f"Inconsistent native parameters at seed {seed}")
        common = comparison
        flux_path = directory/f"flux{suffix}"
        values = np.loadtxt(flux_path, ndmin=2)
        if not np.all(np.isfinite(values)) or np.any(values[:, 1] <= 0):
            parser.error(f"Invalid spectrum at seed {seed}")
        if energy is not None and not np.array_equal(energy, values[:, 0]):
            parser.error("Native energy grids differ")
        energy = values[:, 0]
        fluxes.append(values[:, 1])
        summary = read_pairs(directory/f"summary{suffix}")
        totals.append(float(summary["total_spin_down_erg_s"]))
        local.append(int(summary["sources_within_1kpc"]))
        outside.append(int(summary["outside_halo"]))
        checks.append(np.max(np.loadtxt(directory/f"checks{suffix}")[:, 1:], axis=0))
        table_hashes.append(hashlib.sha256(flux_path.read_bytes()).hexdigest())
    if args.publish:
        expected = {"mspsources": 49000, "mspperiodminms": 1.5, "msplog10b": 8,
                    "mspsigmalog10b": .2, "mspinertiagcm2": 1e45, "mspemingev": 10,
                    "injslope": 1, "injemax": 400, "efficiency": 1, "a": 1.9, "b": 5,
                    "r1kpc": 0, "rgkpc": 20, "sunkpc": 8.5, "discsizepc": 1000,
                    "hkpc": 5, "e0": 10, "d0h": .04, "delta": .56, "ddelta": -1,
                    "bmug": 1, "urad": .475, "emin": 10, "emax": 2000, "esize": 81,
                    "mspquadrature": 256, "mspdistances": 1601}
        if params["spiralmodel"] != "Jelly" or any(
                not np.isclose(float(params[k]), v, rtol=1e-10, atol=0) for k, v in expected.items()):
            parser.error("Slide text describes jelly.ini; changed parameters require a text revision before publishing")
    spectra = np.array(fluxes)
    check_max = np.max(checks, axis=0)
    if not np.all(np.isfinite(check_max)) or np.max(check_max) > .005:
        parser.error("Native numerical checks exceed tolerance")
    first_seed = manifest["seeds"][0]
    first_dir = args.runs/f"seed_{first_seed:06d}"/manifest["model"]
    population = np.loadtxt(first_dir/f"population_{first_seed:06d}.txt")
    horizon = np.loadtxt(first_dir/f"horizon_{first_seed:06d}.txt")
    contributions = np.loadtxt(first_dir/f"contributions_{first_seed:06d}.txt")
    data = read_ams(args.data)
    q16, median, q84 = np.quantile(spectra, [.16, .5, .84], axis=0)
    result = dict(engine=manifest["engine"], manifest=manifest, parameters=common,
                  flux_sha256=table_hashes, data_path=str(args.data.resolve()),
                  data_sha256=hashlib.sha256(args.data.read_bytes()).hexdigest(),
                  checks=dict(table_vs_direct=float(check_max[0]), refinement=float(check_max[1])),
                  median_total_spin_down_erg_s=float(np.median(totals)),
                  median_sources_within_1kpc=float(np.median(local)),
                  first_seed=first_seed, first_outside_halo=outside[0], diagnostics={})
    macros = dict(MspRealizations=str(len(spectra)), MspSourceCount=f"{int(params['mspsources']):,}",
                  MspFirstSeed=str(first_seed), MspLastSeed=str(manifest["seeds"][-1]),
                  MspTotalPower=f"{np.median(totals)/1e38:.2f}", MspLocalCount=f"{np.median(local):.0f}",
                  MspGridError=f"{100*check_max[1]:.3f}", MspDirectError=f"{100*check_max[0]:.3f}")
    for e, suffix in {10: "Ten", 100: "Hundred", 300: "ThreeHundred", 1000: "Thousand"}.items():
        idx = np.flatnonzero(np.isclose(energy, e))
        if len(idx) != 1:
            parser.error(f"Need exactly one {e} GeV diagnostic row")
        idx = idx[0]
        row = contributions[np.isclose(contributions[:, 0], e)][0]
        values = dict(q16=float(q16[idx]), median=float(median[idx]), q84=float(q84[idx]),
                      neff=float(row[1]), brightest_fraction=float(row[2]), d50_kpc=float(row[3]), d90_kpc=float(row[4]))
        macros[f"MspFlux{suffix}"] = f"{e**3*median[idx]:.2f}"
        macros[f"MspRadius{suffix}"] = f"{row[4]:.2f}"
        macros[f"MspNeff{suffix}"] = f"{row[1]:.0f}"
        if data[0].min() <= e <= data[0].max():
            ratio = median[idx]/np.exp(np.interp(np.log(e), np.log(data[0]), np.log(data[1])))
            values["median_to_ams_ratio"] = float(ratio)
            macros[f"MspRatio{suffix}"] = f"{ratio:.1f}"
            macros[f"MspEfficiency{suffix}"] = f"{100/ratio:.1f}"
        result["diagnostics"][str(e)] = values
    args.output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output/"ensemble.npz", energy_GeV=energy, flux=spectra)
    np.savetxt(args.output/"spectrum.csv", np.column_stack([energy, q16, median, q84, spectra[0]]),
               delimiter=",", header="energy_GeV,q16_flux,median_flux,q84_flux,first_flux; flux: GeV^-1 m^-2 s^-1 sr^-1")
    (args.output/"summary.json").write_text(json.dumps(result, indent=2)+"\n")
    (args.output/"numbers.tex").write_text("% Generated from gryphon2/runMSP outputs; do not edit.\n"+
        "\n".join("\\newcommand{\\"+k+"}{"+v+"}" for k, v in macros.items())+"\n")
    plot(energy, spectra, population, horizon, params, data, args.output)
    if args.publish:
        for name in ("spectrum", "population", "horizon"):
            file = f"msp_synthetic_{name}.pdf"
            shutil.copy2(args.output/file, args.publish/"figures"/file)
        native = args.publish/"calculation/native_output"
        native.mkdir(parents=True, exist_ok=True)
        for file in ("numbers.tex", "summary.json", "spectrum.csv", "ensemble.npz"):
            shutil.copy2(args.output/file, native/file)
    print(json.dumps(dict(realizations=len(spectra), checks=result["checks"], diagnostics=result["diagnostics"]), indent=2))


if __name__ == "__main__":
    main()
