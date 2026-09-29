# LHAASO proton anisotropy benchmarks

`source_knee.ini` and `transport_knee.ini` share the same source population,
low-energy injection slope, and LHAASO knee parameters. The former places the
smooth break in the source spectrum; the latter places it in the diffusion
coefficient. Both grids extend to 100 PeV so the high-energy separation and the
limit of the diffusion approximation are visible.

From `build/`, run an ensemble with:

```sh
for seed in $(seq 0 999); do
    GRYPHON_CR_THREADS=4 ./runAnisotropy \
        ../configs/lhaasoanisotropy/source_knee.ini \
        --seed "$seed" --outdir ../runs/lhaasoanisotropy
    GRYPHON_CR_THREADS=4 ./runAnisotropy \
        ../configs/lhaasoanisotropy/transport_knee.ini \
        --seed "$seed" --outdir ../runs/lhaasoanisotropy
done
```

Each `flux_<seed>.txt` contains the proton intensity and all three components
of the dimensionless dipole vector. The components use the simulation frame;
they are not yet the first harmonic in right ascension measured by LHAASO.

Build and publish the current comparison figure with:

```sh
cd plotting/projects/lhaasoanisotropy
make
make publish
```
