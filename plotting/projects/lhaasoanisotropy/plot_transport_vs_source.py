"""Spectrum-degenerate proton-knee models and their dipole predictions."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import gryphon_plots as gp


# The IAU equatorial-to-Galactic rotation matrix. For row vectors the inverse
# transformation is multiplication by this matrix. gryphon's +x axis points
# away from the Galactic centre, while the conventional Galactic x axis points
# towards it.
EQUATORIAL_TO_GALACTIC = np.array(
    [
        [-0.0548755604, -0.8734370902, -0.4838350155],
        [0.4941094279, -0.4448296300, 0.7469822445],
        [-0.8676661490, -0.1980763734, 0.4559837762],
    ]
)


def lhaaso_exposure_factor(latitude_deg=29.3577, zenith_max_deg=40.0):
    """Coefficient mapping an equatorial dipole to the RA first harmonic."""
    latitude = np.deg2rad(latitude_deg)
    zenith_max = np.deg2rad(zenith_max_deg)
    declination = np.linspace(-np.pi / 2, np.pi / 2, 200001)
    xi = (
        np.cos(zenith_max) - np.sin(latitude) * np.sin(declination)
    ) / (np.cos(latitude) * np.cos(declination))
    alpha_max = np.where(
        xi <= -1,
        np.pi,
        np.where(xi >= 1, 0, np.arccos(np.clip(xi, -1, 1))),
    )
    exposure = (
        np.cos(latitude) * np.cos(declination) * np.sin(alpha_max)
        + alpha_max * np.sin(latitude) * np.sin(declination)
    )
    normalization = np.trapezoid(exposure * np.cos(declination), declination)
    harmonic = np.trapezoid(exposure * np.cos(declination) ** 2, declination)
    return harmonic / normalization


def lhaaso_first_harmonic(ensemble):
    """Small-anisotropy RA amplitude and phase for every realization."""
    simulation = np.stack(
        [
            ensemble.stack("dipole_x"),
            ensemble.stack("dipole_y"),
            ensemble.stack("dipole_z"),
        ],
        axis=-1,
    )
    galactic = simulation.copy()
    galactic[..., 0] *= -1
    equatorial = galactic @ EQUATORIAL_TO_GALACTIC
    amplitude = lhaaso_exposure_factor() * np.hypot(
        equatorial[..., 0], equatorial[..., 1]
    )
    phase = np.rad2deg(np.arctan2(equatorial[..., 1], equatorial[..., 0]))
    return amplitude, phase


def lhaaso_data():
    path = Path(__file__).resolve().parents[3] / "data" / "lhaaso_proton_2025.txt"
    return np.loadtxt(path).T


def shape_normalization(energy_pev, flux_pev, data_energy, data_flux):
    mask = (data_energy >= 0.3) & (data_energy <= 3.0)
    model = np.interp(np.log(data_energy[mask]), np.log(energy_pev), np.log(flux_pev))
    return np.exp(np.mean(np.log(data_flux[mask]) - model))


gp.style.use()
plt.rcParams.update(
    {
        "font.size": 8,
        "font.weight": "normal",
        "axes.labelsize": 10,
        "axes.labelpad": 4,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.major.pad": 3,
        "ytick.major.pad": 3,
        "lines.linewidth": 1.6,
    }
)
source = gp.load_ensemble("source_knee")
transport = gp.load_ensemble("transport_knee")
data_e, data_f, data_stat, data_syst = lhaaso_data()

fig, (ax_flux, ax_dipole) = plt.subplots(
    1, 2, figsize=(7.0, 2.6), sharex=True, constrained_layout=True
)

for ensemble, label, color in (
    (source, "source knee", "tab:blue"),
    (transport, "transport knee", "tab:red"),
):
    energy = ensemble.x / 1e6
    fluxes = ensemble.stack("I") * 1e6
    median_flux = np.median(fluxes, axis=0)
    scale = shape_normalization(energy, median_flux, data_e, data_f)
    lo_flux, hi_flux = np.percentile(fluxes * scale, [16, 84], axis=0)
    ax_flux.plot(energy, energy**2.75 * median_flux * scale, color=color, label=label)
    ax_flux.fill_between(
        energy, energy**2.75 * lo_flux, energy**2.75 * hi_flux, color=color, alpha=0.2
    )

    harmonic, _ = lhaaso_first_harmonic(ensemble)
    harmonic /= harmonic[:, :1]
    lo_d, med_d, hi_d = np.percentile(harmonic, [16, 50, 84], axis=0)
    ax_dipole.plot(energy, med_d, color=color, label=label)
    ax_dipole.fill_between(energy, lo_d, hi_d, color=color, alpha=0.2)

errors = np.sqrt(data_stat**2 + data_syst**2)
ax_flux.errorbar(
    data_e, data_e**2.75 * data_f, yerr=data_e**2.75 * errors,
    fmt="o", color="black", ms=3, lw=0.8, label="LHAASO protons",
)
for axis in (ax_flux, ax_dipole):
    axis.axvspan(22.3, 100.0, color="0.92", zorder=-10)
ax_flux.set_xscale("log")
ax_flux.set_yscale("log")
ax_flux.set_xlim(0.1, 100.0)
ax_flux.set_xlabel(r"$E$ [PeV]")
ax_flux.set_ylabel(r"$E^{2.75} I_p$ [arb. units]")
ax_flux.legend(frameon=False, loc="lower left")

ax_dipole.axvline(3.3, color="0.5", ls="--", lw=1)
ax_dipole.set_xscale("log")
ax_dipole.set_yscale("log")
ax_dipole.set_xlim(0.1, 100.0)
ax_dipole.set_xlabel(r"$E$ [PeV]")
ax_dipole.set_ylabel(r"$r_1(E)/r_1(0.1\,\mathrm{PeV})$")

gp.savefig(fig, "transport_vs_source")
