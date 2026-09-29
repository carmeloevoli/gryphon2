"""Figure helpers: consistent axes, and one place that decides where PDFs go."""

from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib.pyplot as plt

from .project import Project, project


def new_figure(figsize: tuple[float, float] = (10.5, 8.0)):
    """A figure/axes pair with the gryphon default size."""
    return plt.subplots(figsize=figsize)


def set_axes(
    ax: plt.Axes,
    xlabel: str,
    ylabel: str,
    xscale: str = "linear",
    yscale: str = "linear",
    xlim: tuple | None = None,
    ylim: tuple | None = None,
) -> None:
    """Label and scale an axes object."""
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)


def savefig(
    fig: plt.Figure,
    name: str,
    proj: Project | None = None,
    dpi: int = 300,
    close: bool = True,
) -> Path:
    """Save into the project's figures directory. `name` needs no extension."""
    proj = proj or project()
    filename = name if Path(name).suffix else f"{name}.pdf"
    path = proj.figures / filename
    path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0.1)
    print(f"wrote {path}")
    if close:
        plt.close(fig)
    return path


def publish(proj: Project | None = None, pattern: str = "*.pdf") -> list[Path]:
    """Copy the project's figures into the paper directory declared in project.toml."""
    proj = proj or project()
    if proj.paper is None:
        raise RuntimeError(f"{proj} declares no 'paper' directory in project.toml")
    if not proj.paper.is_dir():
        raise FileNotFoundError(f"paper directory does not exist: {proj.paper}")

    copied = []
    for source in sorted(proj.figures.glob(pattern)):
        target = proj.paper / source.name
        action = "overwrite" if target.exists() else "new"
        shutil.copy2(source, target)
        print(f"{action:>9}  {target}")
        copied.append(target)

    if not copied:
        print(f"nothing to publish: no {pattern} in {proj.figures}")
    return copied
