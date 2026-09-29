# hpwne — hadronic emission from young pulsars

Figures for `~/Work/overleaf/hpwne.paper`.

| Figure (in the paper)              | Script                       | Run                            | Table            |
| ---------------------------------- | ---------------------------- | ------------------------------ | ---------------- |
| `gryphon_spectrum.pdf`             | `plot_injection_spectrum.py` | `youngpulsars_H5`              | `injection_*`    |
| `gryphon_spirals.pdf`              | `plot_spirals.py`            | `events_1Myr`                  | `events_*`       |
| `gryphon_Erot_distribution.pdf`    | `plot_distributions.py`      | `youngpulsars_H5`              | `pulsars_*`      |
| `gryphon_Emax_distribution.pdf`    | `plot_distributions.py`      | `youngpulsars_H5`              | `pulsars_*`      |
| `gryphon_proton_youngpulsars.pdf`  | `plot_proton_flux.py`        | three `youngpulsars_H5*` runs  | `flux_*`         |
| `gryphon_crab_sed.pdf`             | `plot_crab.py`               | none (analytic)                | `data/lhaaso_*`  |

`plot_escape.py` (`gryphon_escape.pdf`, from `propagation_*`) is still built but
is no longer used in the paper — the escape-time panel was dropped from Fig. 2.

`plot_crab.py` is the odd one out: it reads no run at all. The Crab section of
the paper is a single-source calculation — the injection model of Sec. II A
evaluated for the Crab pulsar's own P and Pdot, then folded with the Kelner et
al. (2006) pion-decay yield — so the whole thing lives in the script, and the
numbers quoted in the text (P0, E_M, the residence times, and the limits on xi)
are printed on every run.

Its one input is `data/lhaaso_crab_sed.txt`, the measured WCDA and KM2A points
of the Crab copied verbatim from `LHAASOSED.txt` in the LHAASO public data
release (`publishScience.zip`, linked from
<https://english.ihep.cas.cn/lhaaso/pdl/202110/t20211026_286779.html>), the
release accompanying Cao et al., Science 373, 425 (2021), which holds the points
of Figs. 3 and 4 of that paper. These are the measurement itself, not the
published log-parabola fit. The distinction matters: above 300 TeV the points
sit above the fit with large asymmetric errors, so using the fit as the
comparison would tighten the limits on xi by about an order of magnitude and
overstate the case.

```sh
make            # all figures into ./figures
make publish    # copy them into the paper
```

Every table lives in `runs/hpwne/<run>/` next to the `params.ini` that produced
it, and carries a provenance header. Three of the five paper figures share the
single `youngpulsars_H5` run, so the injection spectrum and the birth-property
distributions are guaranteed to describe the same population as the predicted
flux.

## Regenerating the runs

The models are defined in [`configs/hpwne/`](../../../configs/hpwne). From the
build directory:

```sh
# flux ensembles (100 realizations each)
for c in youngpulsars_H5 youngpulsars_H5_100msec youngpulsars_H5_flatD \
         youngpulsars_H5_{100,130}msec_flatD; do
    ../scripts/run_ensemble.sh ./runYoungPulsars ../configs/hpwne/$c.ini 100 ../runs/hpwne
done

# diagnostics of the fiducial model
./inspectRandomPulsars ../configs/hpwne/youngpulsars_H5.ini --outdir ../runs/hpwne
./inspectPropagation   ../configs/hpwne/youngpulsars_H5.ini --outdir ../runs/hpwne
./inspectEvents        ../configs/hpwne/events_1Myr.ini     --outdir ../runs/hpwne
```

## Model parameters

All runs use H = 5 kpc, and `runYoungPulsars` is run with `efficiency = 1`, so
the CR efficiency ξ is applied in `plot_proton_flux.py` as an explicit rescaling
of the flux:

- `youngpulsars_H5` — fiducial, P₀ = 60 ms, D ∝ E^0.36, plotted with ξ = 0.1
- `youngpulsars_H5_100msec` — P₀ = 100 ms, plotted with ξ = 1
- `youngpulsars_H5_130msec_flatD` — P₀ = 130 ms with δ = 0, plotted with ξ = 0.15,
  fixed by matching the LHAASO proton flux at 3 PeV (the script recomputes and
  prints this value on every run — currently 0.147 — so the rounded ξ in the
  label stays checkable)

`youngpulsars_H5_flatD` (60 ms) and `_100msec_flatD` are also kept: they are the
δ = 0 cases at other birth periods, quoted in the text but not plotted.

The band is the 5th–95th percentile envelope over the 100 realizations, matching
the definition in the paper.

Geometry constants are read from `params.ini` rather than restated in the
scripts: `plot_spirals.py` takes R_sun and the integration time from the run,
`plot_escape.py` takes the reference halo size, and `plot_injection_spectrum.py`
takes P₀ and B_*.

## Known gaps

- The original scripts for these figures were not all recoverable, so axis
  ranges and annotations are reconstructions: they reproduce the physics but
  may differ cosmetically from the PDFs currently in the paper. Check before
  `make publish`.
- The spiral snapshot and the birth-property histograms are new random draws
  (the diagnostics are now seeded from their config rather than from whatever
  state the old tools happened to be in). The distributions are unchanged —
  median E_rot is 6.1 × 10⁴⁸ erg, as quoted in the paper — but the individual
  realization plotted in `gryphon_spirals.pdf` is not the same one.
- The flat-D model (`youngpulsars_H5_flatD`) reads the user's "D = 0.42
  kpc^2/Myr constant in energy" as delta = 0 at the same normalization, i.e.
  D/H = 0.42 kpc/Myr at all energies. That is the reading that reproduces the
  factor ~10 drop in D/H at 1 PeV the paper's text anticipates.
