"""Run-backed geometry and transport for the two panels of Figure 2."""

import numpy as np

import gryphon_plots as gp

RUN = "figure2_xie2024"
HALOS = ((2.0, "tab:orange"), (4.0, "tab:red"), (8.0, "tab:purple"))


def load_reference():
    """Reject stale diagnostics from another source/transport prescription."""
    params = gp.load_params(RUN)
    for key, expected in (("spiralmodel", "Xie2024"),
                          ("transportmodel", "PureDiffusion"),
                          ("syntheticassociations", "false")):
        if str(params[key]).lower() != expected.lower():
            raise ValueError(f"Figure 2 requires {key}={expected}")
    for key, expected in (("d0h", 0.42), ("e0", 1000.0), ("delta", 0.36),
                          ("ddelta", -1.0), ("snrateyr", 0.02),
                          ("maxtimemyr", 1.0)):
        if not np.isclose(params[key], expected, rtol=1e-10, atol=0):
            raise ValueError(f"Figure 2 requires {key}={expected}")
    return params


def escape_time_myr(energy_gev, halo_kpc, params):
    """H^2/(2D), with D/H in kpc/Myr and energies in GeV."""
    diffusion = params["d0h"] * halo_kpc * (
        np.asarray(energy_gev) / params["e0"]
    ) ** params["delta"]
    return halo_kpc**2 / (2.0 * diffusion)


def horizon_kpc(halo_kpc):
    return np.sqrt(2.0) * np.asarray(halo_kpc)


def galactocentric_xy(events, params):
    """inspectEvents stores Sun-centred positions; Gryphon's Sun is at (+R,0)."""
    return events["x"] + params["sunkpc"], events["y"]
