# 🦅 GRYPHON

**Green-function cosmic-ray propagation with stochastic sources**

### What is GRYPHON?

A fast, semi-analytical framework for Galactic cosmic-ray transport.
Built on Green functions. Designed for discrete, stochastic sources.

### Features

- Green-function solver
- Flexible injection models
- tochastic source populations
- Nuclei + electrons

### Quick start

```sh
git clone https://github.com/carmeloevoli/gryphon2.git
cd gryphon2
mkdir build
cd build
cmake ..
make -j
make tests
```

### Plots

Figures live with the code, one directory per project:

```sh
cd plotting/projects/hpwne
make            # build the figures
make publish    # copy them into the paper
```

See [plotting/README.md](plotting/README.md).

For the proton paper's nine current-model ensembles (10,000 realizations each),
see the [overnight results runner](scripts/RESULTS_OVERNIGHT.md).

### Dependencies and requirements

Required:
- CMake (≥ 3.1)
- GCC (`g++`) or Clang (`clang++`) with C++17 support
- GNU Scientific Library (GSL)

### Status

🚧 Work in progress

### Citation

Coming soon

### Contact

Carmelo Evoli
