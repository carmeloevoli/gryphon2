"""Matched ten-Myr source maps, using C++ event/parent diagnostics only."""

import json

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

import gryphon_plots as gp


RUNS = ("snapshot_independent_10myr", "snapshot_clustered_10myr")
WINDOW_MYR = 10.0
HALF_WIDTH = 2.0  # kpc, chosen before sampling; same square in both panels
N_HIGHLIGHT = 6
COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")


def window_mask(table, half_width=HALF_WIDTH):
    """Face-on, Sun-centred square; include all heights and all selected ages."""
    return (np.abs(table["x"]) <= half_width) & (np.abs(table["y"]) <= half_width)


def highlighted_parents(table, r_sun, half_width=HALF_WIDTH, number=N_HIGHLIGHT):
    """Nearest parent centres with >=2 visible events, not richest groups.

    Deterministic distance ordering (ID breaks ties); no seed search or
    ranking by spectral effect. Centres must themselves be inside the view.
    """
    visible = window_mask(table, half_width)
    ids, counts = np.unique(table["parent_id"][visible].astype(int), return_counts=True)
    candidates = []
    for parent, count in zip(ids, counts):
        if parent <= 0 or count < 2:
            continue
        index = np.flatnonzero(table["parent_id"] == parent)[0]
        cx = table["centre_x"][index] - r_sun
        cy, cz = table["centre_y"][index], table["centre_z"][index]
        if abs(cx) <= half_width and abs(cy) <= half_width:
            candidates.append((cx**2 + cy**2 + cz**2, int(parent)))
    return [parent for _, parent in sorted(candidates)[:number]]


def validate_pair(params):
    """Do not silently compare another window, rate, geometry or population."""
    ignored = {"simname", "associationfraction"}
    left, right = params
    if {k: v for k, v in left.items() if k not in ignored} != {
        k: v for k, v in right.items() if k not in ignored
    }:
        raise ValueError("snapshot parameters differ beyond association fraction/name")
    for index, p in enumerate(params):
        for key, expected in (("syntheticassociations", "true"), ("spiralmodel", "Xie2024")):
            if str(p[key]).lower() != expected.lower():
                raise ValueError(f"requires {key}={expected}")
        for key, expected in (("associationfraction", index), ("seed", 0),
                              ("maxtimemyr", WINDOW_MYR), ("snrateyr", 0.02),
                              ("associationmembers", 100), ("associationradiuspc", 30),
                              ("associationvelocitykms", 3), ("sunkpc", 8.5)):
            if not np.isclose(p[key], expected, rtol=1e-12, atol=0):
                raise ValueError(f"requires {key}={expected}")


def validate_events(table, params):
    if not np.all(np.isfinite(table.values)):
        raise ValueError("nonfinite event metadata")
    age, delay, birth = table["age"], table["delay"], table["birth_age"]
    if not np.all((age > 0) & (age < params["maxtimemyr"])):
        raise ValueError("events outside the requested explosion-age window")
    if not np.all((delay >= 3) & (delay <= 48)):
        raise ValueError("delays outside the implemented single-star distribution")
    if not np.allclose(birth - delay, age, rtol=0, atol=1e-9):
        raise ValueError("inconsistent parent and explosion ages")
    parents = table["parent_id"]
    if not np.all(parents == parents.astype(int)):
        raise ValueError("noninteger parent IDs")
    if params["associationfraction"] == 0 and np.any(parents != 0):
        raise ValueError("independent control contains grouped events")
    if params["associationfraction"] == 1 and np.any(parents <= 0):
        raise ValueError("clustered snapshot contains field events")
    for parent in np.unique(parents[parents > 0]):
        members = parents == parent
        for key in ("birth_age", "centre_x", "centre_y", "centre_z"):
            if np.ptp(table[key][members]) > 1e-9:
                raise ValueError(f"inconsistent {key} within parent {parent}")


def plot_snapshot():
    params = [gp.load_params(run) for run in RUNS]
    validate_pair(params)
    tables = [gp.load(f"{run}/event_origins_000000.txt") for run in RUNS]
    for table, p in zip(tables, params):
        validate_events(table, p)
    selected = highlighted_parents(tables[1], params[1]["sunkpc"])
    fig, axes = plt.subplots(1, 2, figsize=(15.6, 8.0), sharex=True, sharey=True)
    summary = {"window_myr": params[0]["maxtimemyr"], "half_width_kpc": HALF_WIDTH, "seed": 0,
               "highlighted_parent_ids": selected, "cases": {}}
    for i, (ax, table, p, run) in enumerate(zip(axes, tables, params, RUNS)):
        mask = window_mask(table)
        ax.scatter(table["x"][mask], table["y"][mask], s=16, c="0.38",
                   alpha=0.75, linewidths=0, zorder=2)
        ax.scatter(0, 0, marker="*", s=260, c="gold", edgecolors="black",
                   linewidths=1.0, zorder=8)
        ax.set(xlim=(-HALF_WIDTH, HALF_WIDTH), ylim=(-HALF_WIDTH, HALF_WIDTH),
               xlabel=r"$x-x_\odot$ [kpc]")
        ax.set_aspect("equal")
        ax.set_xticks(np.arange(-2, 3))
        ax.set_yticks(np.arange(-2, 3))
        ax.tick_params(labelsize=23)
        ax.grid(False, which="both")
        ax.set_title((r"Independent: $f_{\rm a}=0$", r"Clustered: $f_{\rm a}=1$")[i],
                     fontsize=27, pad=18)
        ax.text(0.04, 0.96, rf"$0<a<{p['maxtimemyr']:g}$ Myr", transform=ax.transAxes, va="top",
                fontsize=22, bbox=dict(facecolor="white", edgecolor="none", alpha=.85, pad=2))
        ax.text(0.96, 0.04, rf"${mask.sum()}$ explosions", transform=ax.transAxes,
                ha="right", fontsize=19,
                bbox=dict(facecolor="white", edgecolor="none", alpha=.85, pad=2))
        population = gp.load(f"{run}/population_000000.txt")
        raw_count = int(population["field_explosions"][0] + population["clustered_explosions"][0])
        if int(population["retained_explosions"][0]) != len(table):
            raise ValueError("population count disagrees with metadata")
        summary["cases"][run] = {
            "f_a": p["associationfraction"], "all_events_before_causal_cut": raw_count,
            "all_retained_events": len(table), "events_in_view": int(mask.sum()),
            "parents_in_view": int(len(np.unique(table["parent_id"][mask]))) if i else None,
        }
    axes[0].set_ylabel(r"$y$ [kpc]")
    table = tables[1]
    mask = window_mask(table)
    for parent, color in zip(selected, COLORS):
        members = mask & (table["parent_id"] == parent)
        index = np.flatnonzero(members)[0]
        axes[1].scatter(table["x"][members], table["y"][members], s=16,
                        color=color, linewidths=0, zorder=4)
        axes[1].scatter(table["centre_x"][index] - params[1]["sunkpc"],
                        table["centre_y"][index], marker="+", s=85,
                        linewidths=1.6, color=color, zorder=5)
    handles = [Line2D([], [], marker="*", color="none", markerfacecolor="gold",
                      markeredgecolor="black", markersize=14, label="Sun")]
    axes[0].legend(handles=handles, loc="lower left", fontsize=19,
                   handlelength=1, handletextpad=.3, frameon=True,
                   facecolor="white", edgecolor="none", framealpha=.9)
    fig.tight_layout(w_pad=1.4)
    output = gp.savefig(fig, "gryphon_association_snapshot")
    # Machine-readable provenance/counts; not ensemble or rare-tail statistics.
    gp.project().figures.joinpath("association_snapshot.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return output


if __name__ == "__main__":
    gp.style.use()
    plot_snapshot()
