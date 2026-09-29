"""Reading gryphon output tables.

Tables are plain text with a commented header. Two header dialects are
understood, so plots keep working across the output-format migration:

  legacy (current C++ output)
      # E [GeV] - I [GeV^-1 m^-2 s^-1 sr^-1]

  contract (planned: provenance + explicit column list)
      # gryphon 2.0  git=81e7805  date=2026-08-14T18:59
      # project = hpwne
      # seed = 3
      # columns: E[GeV] I[GeV-1 m-2 s-1 sr-1]

Columns are addressed by name so that a change of column order in the C++ code
cannot silently corrupt a figure.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .project import Project, project

_COLUMN_TOKEN = re.compile(r"^(?P<name>.*?)\s*(?:\[(?P<unit>[^\]]*)\])?$")
_META_LINE = re.compile(r"^(?P<key>[A-Za-z_][\w .-]*?)\s*=\s*(?P<value>.*)$")


def _canonical(name: str) -> str:
    """Loose column-name key: 'rotational energy' == 'rotational_energy' == 'RotationalEnergy'."""
    return re.sub(r"[^a-z0-9]", "", name.lower())


def _parse_column_spec(spec: str) -> tuple[str, str]:
    match = _COLUMN_TOKEN.match(spec.strip())
    if match is None:  # pragma: no cover - the regex matches anything
        return spec.strip(), ""
    return match.group("name").strip(), (match.group("unit") or "").strip()


@dataclass
class Table:
    """A gryphon output table: named columns plus whatever the header recorded."""

    path: Path
    values: np.ndarray
    names: list[str] = field(default_factory=list)
    units: dict[str, str] = field(default_factory=dict)
    meta: dict[str, str] = field(default_factory=dict)

    def __len__(self) -> int:
        return self.values.shape[0]

    @property
    def ncols(self) -> int:
        return self.values.shape[1]

    def index(self, key: str | int) -> int:
        if isinstance(key, int):
            return key
        wanted = _canonical(key)
        for i, name in enumerate(self.names):
            if _canonical(name) == wanted:
                return i
        raise KeyError(f"{self.path.name}: no column '{key}' in {self.names}")

    def __getitem__(self, key: str | int) -> np.ndarray:
        return self.values[:, self.index(key)]

    def unit(self, key: str | int) -> str:
        return self.units.get(self.names[self.index(key)], "")

    def columns(self, *keys: str | int) -> tuple[np.ndarray, ...]:
        """Unpack several columns at once: `E, I = table.columns('E', 'I')`."""
        return tuple(self[key] for key in keys)


def _parse_header(lines: list[str]) -> tuple[list[str], dict[str, str], dict[str, str]]:
    names: list[str] = []
    units: dict[str, str] = {}
    meta: dict[str, str] = {}

    for raw in lines:
        line = raw.lstrip("#").strip()
        if not line:
            continue

        if line.lower().startswith("columns:"):
            spec_list = line.split(":", 1)[1]
            # '|' separates columns whose units contain spaces
            specs = spec_list.split("|") if "|" in spec_list else spec_list.split()
        elif " - " in line and "=" not in line:
            specs = line.split(" - ")
        else:
            match = _META_LINE.match(line)
            if match:
                meta[_canonical(match.group("key"))] = match.group("value").strip()
            continue

        for spec in specs:
            name, unit = _parse_column_spec(spec)
            if not name:
                continue
            names.append(name)
            units[name] = unit

    return names, units, meta


def load(path: str | os.PathLike, proj: Project | None = None) -> Table:
    """Load a table. Relative paths are resolved against the project's runs directory."""
    path = Path(path).expanduser()
    if not path.is_absolute():
        path = (proj or project()).run(str(path))
    if not path.is_file():
        raise FileNotFoundError(f"missing gryphon output: {path}")

    header: list[str] = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.startswith("#"):
                break
            header.append(line)

    names, units, meta = _parse_header(header)
    values = np.atleast_2d(np.loadtxt(path))

    # An unnamed or partially named header must not shift column addressing.
    if len(names) != values.shape[1]:
        names = names[: values.shape[1]]
        names += [f"c{i}" for i in range(len(names), values.shape[1])]

    return Table(path=path, values=values, names=names, units=units, meta=meta)


def load_params(path: str | os.PathLike, proj: Project | None = None) -> dict[str, float | str]:
    """Read the `params.ini` written next to a run.

    `path` may be the file itself or the run directory holding it. Values that
    look numeric are converted, so a plot can take a geometry constant from the
    run that produced its data instead of restating it.
    """
    path = Path(path).expanduser()
    if not path.is_absolute():
        path = (proj or project()).run(str(path))
    if path.is_dir():
        path = path / "params.ini"
    if not path.is_file():
        raise FileNotFoundError(f"missing parameter file: {path}")

    params: dict[str, float | str] = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line or "=" not in line:
                continue
            key, value = (part.strip() for part in line.split("=", 1))
            try:
                params[key] = float(value)
            except ValueError:
                params[key] = value

    return params


@dataclass
class Ensemble:
    """A set of runs of the same model differing only by seed."""

    stem: str
    seeds: list[int]
    tables: list[Table]

    def __len__(self) -> int:
        return len(self.tables)

    @property
    def x(self) -> np.ndarray:
        return self.tables[0][0]

    def stack(self, column: str | int = 1) -> np.ndarray:
        """Realizations as a (nseeds, npoints) array."""
        rows = [t[column] for t in self.tables]
        widths = {len(r) for r in rows}
        if len(widths) != 1:
            raise ValueError(f"{self.stem}: runs have inconsistent lengths {sorted(widths)}")
        return np.vstack(rows)

    def median(self, column: str | int = 1) -> np.ndarray:
        return np.median(self.stack(column), axis=0)

    def band(self, column: str | int = 1, lo: float = 5.0, hi: float = 95.0) -> tuple:
        """Percentile envelope over realizations (default: 5th-95th, as in the paper)."""
        data = self.stack(column)
        return np.percentile(data, lo, axis=0), np.percentile(data, hi, axis=0)


def _seeded_files(directory: Path, prefix: str) -> list[tuple[int, Path]]:
    """Files named `<prefix>_<seed>.txt` in `directory`, as (seed, path).

    The seed pattern is anchored, so `youngpulsars_H4` does not swallow the
    runs of `youngpulsars_H4_100msec`.
    """
    if not directory.is_dir():
        return []

    pattern = re.compile(rf"^{re.escape(prefix)}_(\d+)\.txt$")
    found = [
        (int(match.group(1)), path)
        for path in directory.iterdir()
        if (match := pattern.match(path.name))
    ]
    found.sort()
    return found


def load_ensemble(stem: str, table: str = "flux", proj: Project | None = None) -> Ensemble:
    """Load every seed of one model.

    Two layouts are accepted, in this order:

        <runs>/<stem>/<table>_<seed>.txt   one directory per model (current)
        <runs>/<stem>_<seed>.txt           flat files (pre-`runs/` output)
    """
    proj = proj or project()

    for directory in proj.runs:
        found = _seeded_files(directory / stem, table) or _seeded_files(directory, stem)
        if found:
            return Ensemble(
                stem=stem,
                seeds=[seed for seed, _ in found],
                tables=[load(path, proj) for _, path in found],
            )

    searched = ", ".join(str(d) for d in proj.runs)
    raise FileNotFoundError(f"no runs of '{stem}' in {searched}")
