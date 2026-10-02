"""Record a nested 100/200/300 Myr convergence check at H=8 kpc."""
import argparse
import json
from pathlib import Path
import numpy as np
import gryphon_plots as gp
from feature_statistics import excursions
from plot_halos import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--build-dir", type=Path, required=True)
    args = parser.parse_args()
    args.directory = args.directory.resolve()
    args.build_dir = args.build_dir.resolve()
    params = gp.load_params(args.directory)
    expected = {"hkpc": 8., "emin": 100., "emax": 1e6, "esize": 65.,
                "maxtimemyr": 300., "injslope": 2.34, "varyslope": "false",
                "varyenergy": "false", "injemax": 0.}
    if any(params.get(k) != v for k, v in expected.items()):
        raise ValueError("expected the dedicated 300 Myr H=8 diagnostic")
    files = sorted(args.directory.glob("history_*.txt"))
    if not files:
        raise ValueError("no history diagnostic tables")
    result = {"description": "nested age cuts of the same catalogue and injections; H=8, identical sources",
              "parameters": params, "max_relative_flux_tolerance": 1e-5,
              "max_statistic_change_tolerance": 1e-5,
              "build_sha256": {p.name: digest(p) for p in
                               [args.build_dir / "inspectHistoryConvergence", *args.build_dir.glob("libgryphon.*")]},
              "seeds": {}}
    for path in files:
        x = np.loadtxt(path)
        if (x.shape != (65, 4) or not np.all(np.isfinite(x)) or np.any(x <= 0) or
                not np.allclose(x[:, 0], np.geomspace(100, 1e6, 65), rtol=1e-12)):
            raise ValueError("invalid history table")
        if np.any(np.diff(x[:, 1:], axis=1) < -1e-12 * x[:, -1, None]):
            raise ValueError("nested positive contributions must be monotone")
        stats = [excursions(x[:, 0], x[:, i], fit_range=(100, 1e6), search_range=(1e3, 1e5))
                 for i in (1, 2, 3)]
        flux_change = float(np.max(abs(x[:, 3] - x[:, 2]) / x[:, 3]))
        stat_change = max(float(abs(a[0] - b[0])) for a, b in zip(stats[1], stats[2]))
        result["seeds"][path.stem] = {
            "sha256": digest(path), "100_to_300_max_relative_flux_change": float(np.max(abs(x[:, 3]-x[:, 1])/x[:, 3])),
            "200_to_300_max_relative_flux_change": flux_change,
            "200_to_300_max_absolute_statistic_change": stat_change,
            "statistics_by_history_myr": {str(t): {"residual": float(a[0]), "slope_excursion": float(d[0])}
                                          for t, (a, d) in zip((100, 200, 300), stats)}}
        if flux_change > 1e-5 or stat_change > 1e-5:
            raise ValueError("200 Myr did not pass the history convergence tolerance")
    result["passed"] = True
    (args.directory / "convergence.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
