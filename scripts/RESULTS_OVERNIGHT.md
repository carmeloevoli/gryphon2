# Overnight proton-paper results

For the new **100 GeV--1 PeV** production and **1--100 TeV** analysis, use
[`run_wide_workflow.py`](run_wide_workflow.py) and the
[wide-run guide](../configs/hebreaks/wide/README.md). The commands below retain
the original energy range and remain available to reproduce that analysis.

From the `gryphon2` repository root, run:

```sh
python3 scripts/run_results.py --dry-run
nohup python3 -u scripts/run_results.py >> runs/hebreaks_results_10k.log 2>&1 < /dev/null &
```

The second command launches production in the background and survives closing
the terminal. The runner starts `caffeinate` automatically on macOS to prevent
idle sleep while it runs. Keep the Mac plugged in and the laptop lid open;
`caffeinate` does not guarantee operation with the lid closed or through logout.
Alternatively run `python3 -u scripts/run_results.py` in a terminal kept open.
Python 3.10+ and the existing C++ executables suffice; no Python packages need
installing and **the script does not rebuild the code**.

## What runs

The default is **10,000 total realizations per scenario**, seeds 0--9999,
not 10,000 additional realizations and not a tenfold increase in the number
of supernovae per realization. There are nine scenarios:

| Results subsection | Output directories | Cases |
| --- | --- | --- |
| A: identical sources | `halos_H2`, `halos_H4`, `halos_H8` | H = 2, 4, 8 kpc |
| B: source variations | `variations_energy`, `variations_index` | 0.54 dex energy scatter; parent index width 0.30 conditioned on gamma > 2 |
| D: rare accelerators | `low_rate` | 0.2 per century, energy per accelerator multiplied by ten |
| E: associations | `associations_independent`, `associations_mixed`, `associations_clustered` | f_a = 0, 0.5, 1, with 100 progenitors per association |

Subsection C fits the experimental data and needs no additional independent
Monte Carlo ensemble. Figures in the model section are not rerun.

All spectra use the current Xie distribution, D(1 TeV)/H = 0.42 kpc/Myr,
delta = 0.36, no injection cutoff, normalization above 10 GeV, 100 Myr of
explosion history, and 48 logarithmic energies from 1 TeV to 1 PeV. The
association runs retain the 3--48 Myr delay distribution and stellar spreading,
including in their independent control. Births extend to 148 Myr lookback;
the expected number of explosions before causal selection is 2,000,000.
The ordinary H=4 sample is not substituted for the association control.

The new low-rate and association configurations are in `configs/hebreaks/results/`.
The earlier 30 Myr association configs and all legacy simulations remain untouched.
The default output root is `runs/hebreaks_results_10k/`, with one directory
per model. Old and new physics are never mixed. Verified current-model
halo/variation outputs are **copied**, not hard-linked, into the new root;
at present this reuses 5,000 spectra and leaves 85,000 new realizations to run.

## Progress and completion

```sh
tail -f runs/hebreaks_results_10k.log
python3 scripts/run_results.py --status
```

The status report contains each model's completed/remaining count, elapsed time,
and measured seconds per new realization. Production is complete only when
`state` is `complete` and all nine counts are 10,000. A final `COMPLETE` line
also appears in the log. The status timestamp matters: after a crash or power
loss, the saved report may lag behind the individual completion markers.

Default parallelism is six models with two C++ threads each on this 16-core
machine. Seeds within a model stay serial to protect its shared `params.ini`.
Short batches rotate among models so every scenario advances, rather than
finishing the halos before starting the associations. On slower machines the
default model concurrency is reduced automatically; use `--jobs` and `--threads`
to tune it. A full 90,000-spectrum target is substantial: **do not assume it
will finish within one night**. Partial ensembles remain usable and resumable.
There is no default time limit. To stop after roughly eight production hours:

```sh
python3 -u scripts/run_results.py --max-hours 8
```

The limit stops scheduling work and lets current seeds finish. It is not a
wall-clock guarantee including input verification. A time-limited/stopped run
returns a nonzero status (130), not a false success.

## Restart safely

Run the **same launch command** again. Verified complete seeds are skipped;
interrupted seeds without completion markers are recomputed. A lock prevents
two instances from writing into this output root simultaneously. To stop a
background run cleanly, send `kill -TERM PID`, using the `pid` from `--status`,
then wait for `STOPPED` before restarting. In foreground, use Ctrl-C.

Each model retains its frozen `input.ini`, resolved `params.ini`, build/config
manifest, flux files, logs, and per-seed SHA-256 completion markers. Association
markers also cover the population-count table. The runner checks the energy
grid, positive finite fluxes, model/seed headers, resolved parameters, and
population normalization before marking a seed complete. It checks executable
and `libgryphon` hashes before and after each new seed.

**Do not rebuild or replace the C++ executables/libraries while running or
between restarts of the same ensemble.** A different build or configuration
is rejected, not silently appended. For an intentionally changed build use
`--no-reuse --outdir runs/ANOTHER_NEW_NAME`. Corrupted completed outputs also
cause an error rather than silently changing the ensemble.

`--seeds N` changes the total target and may extend an existing run.
`--models NAME ...` is available for diagnostic/subset runs; it is not needed
for the complete overnight production. `--dry-run` performs read-only
build/config checks and reports eligible imports; full completion checksums
are verified at startup, not during that quick preview.

## Tomorrow's analysis

The runner only generates and validates simulation outputs. It does not modify
the manuscript, publish figures, or automatically replace 1,000-realization
statistics with incomplete 10,000-realization results. `results_plan.json`
records the run definitions and `status.json` records progress.

The plotting project now defaults to `runs/hebreaks_results_10k` and 10,000
realizations per scenario. From the repository root:

```sh
make -C plotting/projects/hebreaks test publish-results PYTHON="$PWD/.venv/bin/python"
make -C /Users/carmeloevoli/Work/overleaf/hebreaks.paper
```

This verifies every input checksum and the matched model parameters, then
updates all five Monte Carlo result PDFs and four numerical tables/summaries.
Use `results` instead of `publish-results` for a non-publishing preview.
The current low-rate and association path is `plot_population_results.py`;
**do not use the historical `plot_low_rate.py` or `plot_associations.py`**.
Text and captions must still be reviewed when changing the ensemble size.

Tests (do not launch production):

```sh
python3 -m unittest discover -s scripts -p 'test_run_*.py'
```
