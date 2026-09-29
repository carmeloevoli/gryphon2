# Synthetic stellar associations

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
