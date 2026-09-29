"""The one and only gryphon matplotlib style."""

from __future__ import annotations

import os
from pathlib import Path

STYLE_FILE = Path(__file__).with_name("gryphon.mplstyle")


def use(backend: str | None = None) -> None:
    """Activate the gryphon style.

    Scripts save figures rather than show them, so a non-interactive backend is
    selected by default; override with `backend` or $GRYPHON_MPL_BACKEND
    (e.g. "MacOSX") to work interactively.
    """
    import matplotlib

    matplotlib.use(backend or os.environ.get("GRYPHON_MPL_BACKEND", "Agg"))

    import matplotlib.pyplot as plt

    plt.style.use(STYLE_FILE)
