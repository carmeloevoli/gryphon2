# Synthetic stellar associations

For the complete **10,000-realization per scenario** results production, use
[`scripts/run_results.py`](../../scripts/run_results.py), documented in
[the overnight guide](../../scripts/RESULTS_OVERNIGHT.md). Its new low-rate
and association configs in `results/` all use 100 Myr of explosion history.
The older 30 Myr benchmark configurations described below remain unchanged.

The `figure5_H2.ini` and `figure5_H4.ini` configs regenerate the identical-source
proton spectra in Figure 5, independently of the association benchmark below.
Use `scripts/run_figure5.py` with `runAnisotropy`; the reproduction commands,
fixed seeds and display normalization are documented in
[the plotting README](../../plotting/projects/hebreaks/README.md#figure-5-individual-proton-realizations).

The separate `figure2_xie2024.ini` config is a 1 Myr, independent-source
snapshot for Figure 2, with the same Xie2024 geometry and D0/H=0.42 transport
normalization as the reference model. It is used by `inspectEvents` and
`inspectPropagation`, not by the association-ensemble runner. Reproduction
commands are in [the plotting documentation](../../plotting/projects/hebreaks/README.md).

The separate `snapshot_independent.ini` and `snapshot_clustered.ini` configs
provide the matched 10 Myr association maps in the model section. They are
also consumed by `inspectEvents`, not the ensemble runner; reproduction
commands and parent-ID conventions are in the same plotting documentation.
These event/escape diagnostics do not evaluate the injection spectrum; their
historical cutoff entries have no effect on Figures 2 or 3.

This is a controlled clustering experiment, **not an OB-association catalogue
reconstruction**. `runAssociations` uses the usual transport/injection machinery.
Existing generators retain their behaviour unless `syntheticassociations=true`.

## Population

The three supplied configs differ only in `associationfraction` (0, 0.5, 1) and
their names. They all have a Galactic SN rate of 2/century, the same energy per
SN, Xie2024 centres, H=4 kpc, and 30 Myr of explosion history. Association births
are Poisson at f R / N, with **exactly N=100 core-collapse progenitors per parent**.
This is not a total stellar membership or a measured richness distribution.

Rather than inventing a mass-lifetime law, sample the normalized **single-star**
delay-time fit in Zapartas et al. (2017), A&A 601 A29, Appendix A, Eq. (10):
https://doi.org/10.1051/0004-6361/201629685
The IMF is already integrated into this fit. Delays are 3--48 Myr. The sampler
integrates the piecewise density analytically and interpolates a table of 8193
inverse-CDF quantiles (maximum CDF discretization <=1/8192).

Parents extend to T+48 Myr lookback; an event has age birth_age - delay and is
retained if 0 < age < T. This is essential for stationarity. The expectation
before the common light-cone cut is R*T=600,000 SNe for every supplied scenario.
Realized counts fluctuate; do not renormalize each realization to that count.

Offsets are isotropic Gaussians with one-dimensional standard deviation
sqrt(r_a^2 + (sigma_v * delay)^2), r_a=30 pc and sigma_v=3 km/s. These are
illustrative parameters, not catalogue-derived measurements. Positions are
conditioned on the existing Galactic-volume and 1 pc observer exclusions.
The independent component draws the **same** delays and offsets but a new centre
per event, matching the one-point source density of the clustered case.
An independent control without this spreading would confound clustering with
changes in the mean source distribution. Temporal stationarity is also matched.

The parent centres are stationary; no Galactic orbital motion, binary-delayed
SNe, runaways, mass-dependent CR efficiency, Type Ia component, or collective
wind/superbubble acceleration is included. Each SN remains a distinct burst.

## Run and extend

From the repository root:

```sh
cmake --build build --target runAssociations gryphon_tests -j 6
python3 scripts/run_associations.py --seeds 256 --outdir runs/hebreaks_associations_uncut
```

The runner uses up to three concurrent models, each with three energy-bin
threads; `--jobs` and `--threads` control these. Seeds are serial within each
model to avoid concurrent writes to params.ini. Output goes to
`runs/hebreaks_associations_uncut/associations_{independent,mixed,clustered}/`
with the command above (the runner's default root remains
`runs/hebreaks_associations/`).
The flux headers identify each seed; params.ini records the last executed seed.
Per-seed completion records checksum the spectrum and population-count table.
The manifest records config text, executable hash, and shared-library hash. Re-running the command
skips verified complete seeds; interrupted seeds are rerun. A changed config
or executable requires a new `--outdir` to avoid combining different ensembles.

Increase `--seeds` to extend the ensemble with the same build and settings.
Change N, radius, or velocity in **all three** copied configs for a matched
sensitivity test, and supply `--config-dir` and a new `--outdir`. Association
richness, a richness distribution, delay prescriptions, and Galactic motion
must be tested before interpreting this benchmark as a Galactic population.

## Figures and statistics

Python requires numpy, matplotlib, scipy. Install the plotting package in a
virtual environment, or use the project Makefile with a suitable interpreter:

```sh
make -C plotting/projects/hebreaks test PYTHON="$PWD/.venv/bin/python"
make -C plotting/projects/hebreaks PYTHON="$PWD/.venv/bin/python"
make -C plotting/projects/hebreaks publish PYTHON="$PWD/.venv/bin/python"
```

Set `GRYPHON_RUNS` to the absolute path of `runs/hebreaks_associations_uncut`
to plot the new ensemble; `GRYPHON_PAPER` changes the publish
destination. Publishing copies only the association figure and statistics.
The figure follows gryphon.mplstyle, uses large panels, and has no grid.
The residual is max|J/J_PL-1| in 10--100 TeV, fitting unweighted log J over
1 TeV--1 PeV. Running slopes use five centred bins and exclude two edge bins.
Exact pointwise 95% binomial intervals accompany the survival curves; these are
not simultaneous confidence bands. Zero counts have a separate one-sided 95%
upper limit in association_summary.json. No significance for a data-model fit
is claimed: the experimental fit uses different sampling/uncertainties.

The configs match the paper and **current numerical runModels driver**:
D0/H=0.42 kpc/Myr, Emax -> infinity (`injemax=0`), and exact CR-energy
normalization above 10 GeV. All association cases have fixed gamma=2.34>2,
so the uncut energy integral converges. When enabled,
source-energy scatter is 0.54 dex (standard deviation of log10 energy);
the association benchmark keeps it disabled in all three cases.
The existing pilot used Emax=10^12 GeV and normalization above 1 GeV.
Any later 10 PeV runs likewise differ from the new no-cutoff prescription.
Its statistics and figures must be regenerated in a **new output directory**,
not relabelled or mixed with new realizations. No simulation output is changed
by updating these configs. This comparison also remains separate from the
legacy Steiman2010 ensemble.

The initial 256-realization/model run is a **pilot**, not a measurement of
10^-4--10^-5 tail probabilities. With no exceedances, its one-sided 95% limit
is about 0.0116, not zero.
