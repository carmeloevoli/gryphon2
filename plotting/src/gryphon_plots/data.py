"""Experimental cosmic-ray data from the KISS database.

The tables are expected as

    E  y  err_sta_lo  err_sta_up  err_sys_lo  err_sys_up

with the directory declared as `kiss` in project.toml.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from .project import Project, project


def kiss_path(filename: str, proj: Project | None = None):
    proj = proj or project()
    if proj.kiss is None:
        raise RuntimeError(f"{proj} declares no 'kiss' directory in project.toml")
    path = proj.kiss / filename
    if not path.is_file():
        raise FileNotFoundError(f"missing KISS table: {path}")
    return path


def get_data(
    filename: str,
    slope: float = 0.0,
    norm: float = 1.0,
    emin: float = 0.0,
    emax: float = np.inf,
    add_systematics: bool = True,
    proj: Project | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (E, E^slope * y, err_lo, err_up), scaled and restricted in energy."""
    E, y, sta_lo, sta_up, sys_lo, sys_up = np.loadtxt(
        kiss_path(filename, proj), usecols=range(6), unpack=True
    )

    if add_systematics:
        err_lo = np.sqrt(sta_lo**2 + sys_lo**2)
        err_up = np.sqrt(sta_up**2 + sys_up**2)
    else:
        err_lo, err_up = sta_lo, sta_up

    scaling = norm * np.power(E / norm, slope)
    y, err_lo, err_up = scaling * y, scaling * err_lo, scaling * err_up

    keep = (E > emin) & (E < emax)
    return E[keep], y[keep], err_lo[keep], err_up[keep]


def plot_data(
    ax: plt.Axes,
    filename: str,
    slope: float = 0.0,
    norm: float = 1.0,
    fmt: str = "o",
    color: str = "tab:gray",
    label: str | None = None,
    zorder: int = 1,
    emin: float = 0.0,
    emax: float = np.inf,
    add_systematics: bool = True,
    proj: Project | None = None,
) -> None:
    """Plot a KISS table as error bars, scaled by E^slope."""
    E, y, err_lo, err_up = get_data(
        filename, slope, norm, emin, emax, add_systematics, proj
    )
    ax.errorbar(
        E,
        y,
        yerr=[err_lo, err_up],
        fmt=fmt,
        color=color,
        markeredgecolor=color,
        label=label,
        capsize=4.2,
        markersize=8.5,
        elinewidth=2.0,
        capthick=2.0,
        zorder=zorder,
    )
