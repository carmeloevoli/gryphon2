"""Project discovery and path resolution.

A project is a directory containing a ``project.toml``:

    name    = "hpwne"
    runs    = "../../../build/output"   # where the C++ code writes its tables
    figures = "figures"                 # where plots are written
    paper   = "~/Work/overleaf/hpwne.paper/figures"
    kiss    = "~/Work/codes/KISS-CosmicRayDataBase/kiss_tables"

Relative paths are resolved against the project directory, ``~`` is expanded.
Every entry can be overridden from the environment with GRYPHON_RUNS,
GRYPHON_FIGURES, GRYPHON_PAPER, GRYPHON_KISS, and the project itself can be
selected with GRYPHON_PROJECT.
"""

from __future__ import annotations

import os
import sys
import tomllib
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

PROJECT_FILE = "project.toml"


@dataclass(frozen=True)
class Project:
    """Resolved paths of a single project."""

    name: str
    root: Path
    runs: tuple[Path, ...]
    figures: Path
    paper: Path | None = None
    kiss: Path | None = None

    def run(self, filename: str) -> Path:
        """Absolute path of a simulation output, searched in each runs directory."""
        for directory in self.runs:
            candidate = directory / filename
            if candidate.exists():
                return candidate

        searched = ", ".join(str(d) for d in self.runs)
        raise FileNotFoundError(f"no '{filename}' in {searched}")

    def __str__(self) -> str:
        return f"project '{self.name}' at {self.root}"


def _resolve(value: str | os.PathLike, root: Path) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (root / path).resolve()


def _starting_points(start: str | os.PathLike | None) -> list[Path]:
    """Where to begin the search: the explicit argument, then $GRYPHON_PROJECT,
    then the current directory, then the directory of the running script."""
    if start is not None:
        return [Path(start)]

    env = os.environ.get("GRYPHON_PROJECT")
    if env:
        return [Path(env)]

    starts = [Path.cwd()]
    main_file = getattr(sys.modules.get("__main__"), "__file__", None)
    if main_file is not None:
        starts.append(Path(main_file))
    return starts


def find_project_dir(start: str | os.PathLike | None = None) -> Path:
    """Walk up looking for project.toml, so scripts work from any directory."""
    tried = []
    for entry in _starting_points(start):
        here = Path(entry).expanduser().resolve()
        if here.is_file():
            here = here.parent
        tried.append(here)

        for candidate in (here, *here.parents):
            if (candidate / PROJECT_FILE).is_file():
                return candidate

    searched = ", ".join(str(p) for p in tried)
    raise FileNotFoundError(
        f"no {PROJECT_FILE} above {searched}; "
        "run from inside a project directory or set GRYPHON_PROJECT"
    )


@lru_cache(maxsize=None)
def _load(root: Path) -> Project:
    with open(root / PROJECT_FILE, "rb") as f:
        cfg = tomllib.load(f)

    def raw(key: str, default=None):
        override = os.environ.get(f"GRYPHON_{key.upper()}")
        return override.split(os.pathsep) if override else cfg.get(key, default)

    def entry(key: str, default: str | None = None) -> Path | None:
        value = raw(key, default)
        if value is None:
            return None
        if isinstance(value, list):
            value = value[0]
        return _resolve(value, root)

    def entries(key: str, default) -> tuple[Path, ...]:
        value = raw(key, default)
        if isinstance(value, str):
            value = [value]
        return tuple(_resolve(v, root) for v in value)

    return Project(
        name=cfg.get("name", root.name),
        root=root,
        # `runs` may list several directories: the first one holding a file wins,
        # which keeps figures working while outputs move between layouts.
        runs=entries("runs", ["output"]),
        figures=entry("figures", "figures"),
        paper=entry("paper"),
        kiss=entry("kiss"),
    )


def project(start: str | os.PathLike | None = None) -> Project:
    """Return the project containing `start` (default: the current directory)."""
    return _load(find_project_dir(start))
