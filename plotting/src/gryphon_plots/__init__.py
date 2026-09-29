"""Shared plotting library for gryphon projects.

Typical use from a project script:

    import gryphon_plots as gp

    gp.style.use()
    table = gp.load("inspect_pure_diffusion.txt")
    fig, ax = gp.new_figure()
    gp.set_axes(ax, xlabel="E [GeV]", ylabel="t [Myr]", xscale="log", yscale="log")
    ax.plot(table["E"], table["t_diff"])
    gp.savefig(fig, "gryphon_escape")

`gp.project()` returns the enclosing project (see project.toml); every path used
above is resolved through it, so scripts contain no absolute paths.
"""

from . import data, figure, io, style
from .data import get_data, plot_data
from .figure import new_figure, publish, savefig, set_axes
from .io import Ensemble, Table, load, load_ensemble, load_params
from .project import Project, find_project_dir, project

__all__ = [
    "Ensemble",
    "Project",
    "Table",
    "data",
    "figure",
    "find_project_dir",
    "get_data",
    "io",
    "load",
    "load_ensemble",
    "load_params",
    "new_figure",
    "plot_data",
    "project",
    "publish",
    "savefig",
    "set_axes",
    "style",
]
