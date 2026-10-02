# Extended proton-spectrum ensemble

This production is separate from the archived 1 TeV--1 PeV ensembles. All
nine scenarios retain their previous physical parameters and consecutive,
unselected seeds, with three changes:

- 65 logarithmic energy points from **100 GeV to 1 PeV** (16 intervals per decade).
- **200 Myr** of explosion history for every model, preserving matched histories
  across halos and injection variants. Association births extend to 248 Myr;
  the expected association-run explosion count is 4,000,000.
- Primary feature search **1--100 TeV**, fitting a reference power law over
  **100 GeV--1 PeV**. Five centered points span 0.25 decades, close to the
  previous 0.2553 decades.

The index-scatter ensemble uses parent width **0.30**, conditioned on gamma>2.
This agrees with the previous run metadata and Results/Appendix; the current
manuscript's model-section value 0.15 is inconsistent and must be reconciled
when the new results are incorporated.

## Convergence and validation

Build separately, without replacing any previous production binaries:

```sh
cmake -S . -B build-hebreaks-wide -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-14 \
  -DBUILD_EXPERIMENTS=OFF -DBUILD_DIAGNOSTICS=ON -DENABLE_TESTING=ON
cmake --build build-hebreaks-wide --target runAnisotropy runAssociations \
  inspectHistoryConvergence gryphon_tests -j 8
build-hebreaks-wide/inspectHistoryConvergence configs/hebreaks/wide/history_H8.ini \
  --outdir runs/hebreaks_history_check
.venv/bin/python plotting/projects/hebreaks/check_wide_history.py \
  runs/hebreaks_history_check/history_H8 --build-dir build-hebreaks-wide
```

Adjust the compiler path for the machine. Install the analysis dependencies
with `python3.13 -m venv .venv` and `.venv/bin/python -m pip install -e ./plotting`.

The diagnostic sums the identical 300 Myr catalogue with nested 100, 200, and
300 Myr age cuts, without drawing different sources or injection properties.
For seed 0 and H=8 kpc, the largest relative flux changes compared with 300 Myr
are 8.51e-4 (100 Myr) and 1.73e-6 (200 Myr). The checker also verifies that
the absolute changes in the two primary statistics are below 1e-5. This is
a history-truncation diagnostic, not a convergence test of rare-tail counts.
Its JSON records parameters, table checksum, and build fingerprints.

## Production and automatic analysis

From the repository root:

```sh
.venv/bin/python scripts/run_results.py --profile wide \
  --build-dir build-hebreaks-wide --dry-run
mkdir -p runs
nohup .venv/bin/python -u scripts/run_wide_workflow.py \
  >> runs/hebreaks_wide_workflow.log 2>&1 < /dev/null &
```

The workflow first completes and analyses an unselected prefix of 100 seeds
per model, then extends the same ensembles to **10,000 per model**. Final
analysis runs only after every requested spectrum has passed verification.
Six simultaneous models use two flux threads each; seeds within a model
remain serial. Pilot timings on this machine imply a multi-day calculation.
The production child prevents idle sleep on macOS. Keep the machine powered
and available; interruption is safely resumable with the same command.

```sh
cat runs/hebreaks_wide_10k/workflow_status.json
.venv/bin/python scripts/run_results.py --profile wide --status
tail -f runs/hebreaks_wide_workflow.log
```

`workflow_status.json` distinguishes simulation, analysis, stopped/failed, and
complete states. `status.json` contains per-model counts and measured times
for the current production stage. Send SIGTERM to the workflow PID for a
graceful stop. Re-running the same command verifies and skips completed seeds.
Do not rebuild the dedicated binaries or edit simulation/analysis sources
during this run. Frozen config/build fingerprints and a workflow source
manifest prevent mixing changed calculations. There is no import from legacy
ensembles, and no existing manuscript content is overwritten.

## Analysis products

`runs/hebreaks_wide_10k/analysis_preview_100/` is explicitly a pilot analysis;
`runs/hebreaks_wide_10k/analysis/` contains the completed requested ensemble.
Each has a provenance-rich `summary.json`, the spectrum/index-quantile figure,
and these three subdirectories:

| Analysis | Reference fit | Search |
| --- | --- | --- |
| `primary_1_100_TeV` | 0.1--1000 TeV | 1--100 TeV |
| `secondary_10_100_TeV` | 0.1--1000 TeV | 10--100 TeV |
| `legacy_reference_10_100_TeV` | 1--1000 TeV | 10--100 TeV |

The first two isolate the search-window change with an identical fitted
reference. The third separates the effect of changing that reference fit;
it still uses the new grid and 200 Myr histories, not the original raw runs.
All three contain four survival-function figures (halo, injection, low rate,
associations), per-seed CSV statistics, a compact CSV, and a LaTeX table.
Summaries include exact binomial intervals, zero-count limits, joint lower
threshold counts, fixed-prefix convergence checks, median-spectrum diagnostics,
and the analytic continuous-source index-mixture contribution. No smooth
curvature is subtracted from the spectra.

The gray amplitude bands remain illustrative benchmarks from the old
observational analysis. They are **not** fitted observational values of the
new estimators, nor data-model p-values. This energy-window extension does not
replace a matched observation/simulation comparison.

To regenerate analysis without rerunning simulations:

```sh
.venv/bin/python plotting/projects/hebreaks/analyse_wide.py \
  --runs runs/hebreaks_wide_10k --seeds 10000 \
  --output runs/hebreaks_wide_10k/analysis
```

The existing paper figures and prose retain their previous statistics until
the completed wider results are reviewed and incorporated together.
