# Synthetic MSPs in gryphon2

The C++ application `apps/runMSP.cpp` samples a present-day MSP population using
the native `galaxy::makeGalaxy` position sampler, draws spin parameters, and sums
the stationary positron flux with the native diffusion/loss kernel. Python only
launches repeated C++ runs and plots their saved outputs.

## Build and run

From the gryphon2 root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target runMSP gryphon_tests -j 4
./build/gryphon_tests --gtest_filter='MSP.*'
./build/runMSP configs/msp/jelly.ini --seed 1 --outdir runs/msp
```

One run contains 49,000 sources, takes a few seconds on the development host,
and includes convergence checks. Change the seed for another population.
An existing flux table is never overwritten: use a new seed or output directory.
Outputs are under `<outdir>/<configuration-stem>/`:

- `population_<seed>.txt`: **heliocentric** x, y, z in kpc, period in seconds,
  field in gauss, Pdot, and spin-down power in erg/s.
- `flux_<seed>.txt`: energy in GeV and **positron-only** intensity in
  GeV^-1 m^-2 s^-1 sr^-1. The equal electron component has the same intensity.
- `contributions_<seed>.txt`: effective source count `(sum J)^2/sum J^2`,
  brightest-source fraction, and distances enclosing 50% and 90% of the flux.
- `horizon_<seed>.txt`: distance-ordered cumulative flux, sampled every 50 sources.
- `checks_<seed>.txt`: interpolation/direct-kernel and resolution comparisons.
- `summary_<seed>.txt`: totals and completion marker, written after all checks.
- `params_<seed>.txt`: resolved, re-runnable parameters for this seed.
  `params.ini` also records the latest run in the same directory.

A failed run can leave partial tables but no completed summary. Do not use its
flux as a result. The batch runner resumes only complete, matching runs.

To run and plot an ensemble (100 realizations by default):

```sh
python3 scripts/run_msp.py --dry-run
python3 scripts/run_msp.py --realizations 100 --jobs 4
.venv/bin/python plotting/projects/msp/plot_spectrum.py
```

The batch runner uses seeds 20260924..20261023 and separate per-seed directories
inside `runs/msp_ensemble`. It resumes matching completed runs and refuses to
mix changed config, executable, or shared-library hashes. For changed code or
parameters, choose a new `--output`. Python 3 launches the executable; only the
plotter needs NumPy and Matplotlib. Its `--runs` selects the ensemble directory.
Add `--publish /path/to/msp.slides` to copy the three plots and computed numbers
into the slides. Publication checks the parameters against the documented slide
benchmark. The pointwise 16--84% band describes population scatter at fixed model
parameters, not a numerical error or parameter-uncertainty band.

## Spatial model and population

`spiralmodel=Jelly` selects **GalaxyJelly**, an axisymmetric spatial distribution.
Its native surface density is

    Sigma(R) proportional to [(R+R1)/(Rsun+R1)]^a exp[-b (R-Rsun)/(Rsun+R1)]
    p(R) proportional to R Sigma(R), 0 <= R <= Rg
    z ~ Normal(0, h^2), azimuth uniform in [0, 2 pi)

The benchmark uses a=1.9, b=5, R1=0, Rg=20 kpc, Rsun=8.5 kpc and
**h=1 kpc (Gaussian standard deviation)**. This vertical width is an explicit
MSP benchmark, not the thin SN disk default. Jelly excludes positions within
1 pc of the observer. There is no upper heliocentric distance selection.

The fixed count N=49,000 and spin distributions follow the MSP2_base benchmark
in [Faucher-Giguere & Loeb (2010), section 2.7 and Table 2](https://arxiv.org/abs/0904.3102v5):

    p(P) = Pmin/P^2, P >= 1.5 ms (no upper period cutoff)
    log10(B/G) ~ Normal(8, 0.2^2)
    Pdot = (B / 3.2e19 G)^2 / P, with P in seconds
    Lsd = 4 pi^2 I Pdot / P^3, I = 1e45 g cm^2

The field convention follows the [ATNF documentation](https://heasarc.gsfc.nasa.gov/W3Browse/all/atnfpulsar.html).
Period, field and position are independent. The count is their 2 million
simulated MSPs times their radio normalization 0.0245 (radio beaming 0.5 and
luminosity floor 0.1 mJy kpc^2). **We replace their spatial distribution with
Jelly without re-fitting that normalization**. This hybrid is an illustrative
benchmark, not a reproduction or new fit of their entire radio population model.
There is no further beaming correction to the escaping pairs.

No supernova times, birth rate, ages, recycling evolution, or orbits enter this
extant-population calculation. `snrateyr`, `timestep`, `maxtimemyr`, the legacy
nuclear `pid`, and the other injection models' parameters do not set this flux.

## Injection and transport

Each MSP emits steadily at its present spin-down power:

    Q+(E) = A (E/E0)^(-alpha) exp(-E/Ec), E >= Emin
    integral_Emin^infinity E Q+(E) dE = eta_pair Lsd / 2

The defaults are alpha=1, Ec=400 GeV, E0=Emin=10 GeV, eta_pair=1. The exact
exponential-cutoff energy normalization includes the lower limit. Efficiency
refers to **both charges combined**, and is deliberately optimistic. This is an
injection rate, not the burst-integrated particle spectrum used by the SN/PWN
factory; those factories reject `injectionmodel=MSP`.

`kernel::ContinuousLosses` integrates the stationary solution

    J+(E) = c/[4 pi b(E)] sum_i integral_E^infinity Q+,i(Es) G_i(lambda^2) dEs
    lambda^2 = 4 integral_E^Es D(E')/b(E') dE'

It reuses `DiffusionLossesKernel::lambda2`, `::b`, and the existing slab boundary
prescription, rather than a separate Python propagator. All internal quantities
use gryphon2 CGS units. The shared diffusion-length expression now uses `expm1`
to avoid subtractive cancellation when Es is very close to E.

The supplied config gives D(E)=0.2 (E/10 GeV)^0.56 kpc^2/Myr, H=5 kpc,
B=1 microgauss and U_rad=0.475 eV/cm^3. The native Thomson loss formula yields
b(10 GeV)=0.1605408 GeV/Myr. Free escape applies at z=+-H, with no radial escape
boundary. Sources outside the slab have zero flux without renormalizing N.
Positions and luminosities are constant over the cooling history. There is no
Klein-Nishina correction, source confinement, solar modulation, or positron
background. The calculated spectrum is the MSP contribution alone.

## Numerics

The injection-energy integral uses 256 Gauss-Legendre nodes in
v=ln[(Es-E)/E], from -32 to ln(40 Ec/E). A free-space response per unit spin-down
power is tabulated at 1,601 log-distance nodes and evaluated at every real and
image-source distance. The lower distance is 0.1 pc; the upper limit covers all
images for R<=Rg. The signed sum uses compensated summation. The image order is
chosen from the maximum diffusion length with a 1e-14 Gaussian-tail tolerance
(orders -8..8 for the benchmark). No source is clipped to the table edges.

Every full-population spectrum is checked with **both** twice the quadrature
and twice the distance resolution. The nearest 16 sources are additionally
checked against direct energy integration using `utils::halo_function`, whose
image/eigenfunction treatment is independent of the cached image interpolation.
Either error exceeding 0.5% aborts the run. The default run's maximum spectrum
change on refinement is about 0.004%.

The eight `MSP.*` unit tests cover energy normalization and charge split,
sampling moments and reproducibility, the native spatial distribution,
cached/direct responses, agreement with the time-integrated native burst
kernel, efficiency scaling, config round-tripping, and invalid inputs.
Same-seed reproducibility refers to the same C++ standard-library implementation.
