# gryphon plotting

Python figures for gryphon projects. One shared library, one directory per
project, so that switching project means changing directory and nothing else.

```
plotting/
  src/gryphon_plots/     the shared library (style, IO, project paths, CR data)
  projects/<name>/       one directory per paper/project
      project.toml       where runs, figures, paper and data live
      plot_*.py          one script per figure
      Makefile           `make`, `make publish`, `make clean`
      figures/           generated PDFs (not tracked)
```

## Running

From a project directory:

```sh
make                      # build every figure into ./figures
make plot_escape.py       # build a single figure
make publish              # copy ./figures/*.pdf into the paper directory
```

The `Makefile` puts `plotting/src` on `PYTHONPATH`, so nothing has to be
installed. If you prefer a real install (to run scripts directly, or from an
editor), use `pip install -e plotting/` in a virtual environment.

## Writing a plot script

```python
import gryphon_plots as gp

gp.style.use()
table = gp.load("inspect_pure_diffusion.txt")     # resolved in the runs dir
fig, ax = gp.new_figure()
gp.set_axes(ax, xlabel="E [GeV]", ylabel="t [Myr]", xscale="log", yscale="log")
ax.plot(table["E"], table["t_diff"])              # columns addressed by name
gp.savefig(fig, "gryphon_escape")                 # written to the figures dir
```

Constants that belong to the simulation are read back from it rather than
restated:

```python
params = gp.load_params("youngpulsars_H4")   # the run's own params.ini
r_sun = params["sunkpc"]
```

Rules that keep this maintainable:

- **No absolute paths in scripts.** Every location comes from `project.toml`,
  which `gp.project()` finds by walking up from the script. Any entry can be
  overridden from the environment (`GRYPHON_RUNS`, `GRYPHON_FIGURES`,
  `GRYPHON_PAPER`, `GRYPHON_KISS`), so a one-off run against a different set of
  outputs needs no edit.
- **Columns by name, not by index.** `gp.load` parses the header of the output
  file, so a change of column order in the C++ code raises a `KeyError`
  instead of silently producing a wrong figure.
- **One script, one figure, named as in the paper.** `plot_escape.py` produces
  `gryphon_escape.pdf`, which is what `\includegraphics` refers to.
- **Publishing is explicit.** Scripts write to the project's local `figures/`;
  `make publish` is the only thing that touches the paper directory.

## Adding a project

Copy `projects/hpwne/project.toml` and `Makefile` into a new directory, set
`name`, `runs` and `paper`, and add scripts. Nothing else is shared state.

## Output format

Config-driven runs are written as one directory per model, and every table
opens with a provenance header:

```
runs/hpwne/youngpulsars_H4/
    params.ini          every parameter as used, re-runnable
    flux_000000.txt     one table per seed
```

```
# gryphon = 2.0
# git = 81e7805f18... (dirty)
# date = 2026-08-14T17:30:55Z
# simname = youngpulsars_H4
# seed = 0
# config = ../configs/hpwne/youngpulsars_H4.ini
# params = params.ini
# columns: E [GeV] | I [GeV^-1 m^-2 s^-1 sr^-1]
```

Everything above `columns:` lands in `table.meta`, so a figure can state the
code version and the parameters it came from. `gp.load` also still reads the
older `# E [GeV] - I [...]` header, and `gp.load_ensemble` accepts both the
directory layout and the flat `<stem>_<seed>.txt` files that preceded it.
