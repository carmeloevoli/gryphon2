# Simulation configurations

One directory per project, one file per model. The file name is the simulation
name, which is also the name of the output directory:

```sh
cd build
./runYoungPulsars ../configs/hpwne/youngpulsars_H4.ini --seed 0 --outdir ../runs/hpwne
```

produces

```
runs/hpwne/youngpulsars_H4/
    params.ini          every parameter as actually used, re-runnable
    flux_000000.txt     one table per seed, with a provenance header
```

Run an ensemble by looping over seeds:

```sh
for seed in $(seq 0 99); do
    ./runYoungPulsars ../configs/hpwne/youngpulsars_H4.ini --seed $seed --outdir ../runs/hpwne
done
```

`params.ini` is a valid configuration file, so a run can always be reproduced
from its own output directory. Any key may be written as `key = value` or
`key value`; see [../data/example_input.txt](../data/example_input.txt) for the
full list and the expected units.

## Diagnostics

The `inspect*` tools take the same command line and write into the same run
directories, so pointing one at a model config puts its diagnostics next to
that model's results:

```sh
./inspectPropagation ../configs/hpwne/youngpulsars_H4.ini --outdir ../runs/hpwne
# -> runs/hpwne/youngpulsars_H4/propagation_000000.txt
```

What a diagnostic *sweeps* stays in its source (the three slopes of
`inspectInjectionSpectrum`, the ten ages of `inspectSingleSource`, the sample
counts); what describes the *model* comes from the configuration.
`configs/diagnostics/` holds standalone configurations for the tools that are
not tied to a project.
