"""Hadronic gamma-ray emission of the Crab Nebula -> gryphon_crab_sed.pdf

A single-source consistency test of the spin-down injection model, in the spirit
of the lepto-hadronic Crab modelling of Peng et al. (2022): the same injection
prescription used for the Galactic population is applied to the Crab pulsar, the
pion-decay emission of the resulting protons is computed, and the result is
compared with the LHAASO measurement.

This figure needs no simulation run. Everything is analytic:

  * the Crab birth properties follow from the measured P, Pdot and the historical
    age under the same magnetic-dipole assumption (n = 3) used in the paper;
  * the proton content follows from the E^-1 spin-down fluence, truncated below
    the potential drop available today -- the model injects nothing below it yet;
  * the pion-decay yield uses the Kelner et al. (2006) parametrization, which is
    accurate above 0.1 TeV and for E_gamma/E_p > 1e-3. The protons sit at tens of
    PeV, so the first condition is always met; the second fails only below about
    30 TeV, far from the peak and from the energy at which the limits below are
    set (0.79 PeV, where E_gamma/E_p ~ 0.03). There the curves rest on the
    Feynman-scaling asymptote F_gamma -> B_gamma/x and are indicative only.

Two limiting assumptions bracket how much of the injected proton energy is still
inside the nebula and available as a target for p-p collisions:

  calorimetric   every proton ever injected is still there. Not physical at these
                 energies -- it is the absolute upper bound on the signal.
  escape limited steady state, W_p = xi * L_sd * t_res, with t_res between the
                 light-crossing time R/c (free streaming, the shortest residence
                 time possible) and the Bohm value at the energy injected today.

The printed output gives, for each case, the efficiency xi that would saturate
the measured flux -- the number quoted in the paper.
"""

import pathlib

import numpy as np

import gryphon_plots as gp

# --- constants, cgs ------------------------------------------------------
C_LIGHT = 2.99792458e10
E_CHARGE = 4.80320425e-10  # esu
M_SUN = 1.98892e33
GEV = 1.602176634e-3  # erg
PC = 3.0856775814913673e18
KPC = 1e3 * PC
YEAR = 3.15576e7
MBARN = 1e-27
M_PROTON = 1.67262192e-24

# --- neutron star, as in Sec. II A of the paper --------------------------
R_STAR = 1.0e6
M_STAR = 1.4 * M_SUN
I_STAR = 0.4 * M_STAR * R_STAR**2

# --- Crab observables ----------------------------------------------------
PERIOD = 33.7e-3  # s
PERIOD_DOT = 4.21e-13  # s/s
AGE = 972.0 * YEAR  # SN 1054
DISTANCE = 2.0 * KPC
R_NEBULA = 2.0 * PC
M_EJECTA = 4.5 * M_SUN
B_NEBULA = 110e-6  # G, the one-zone value of Cao et al. (2021)

# Target density: the ejecta mass spread over the nebular volume. Peng et al.
# treat it as a free parameter and need 0.2-7 cm^-3 to fit the >200 TeV points.
N_TARGET = M_EJECTA / (M_PROTON * (4.0 / 3.0) * np.pi * R_NEBULA**3)

# --- efficiency carried over from the population model -------------------
XI = 0.15  # the most favourable knee model of Sec. III B

EMIN, EMAX = 1e4, 2e7  # GeV, gamma-ray energy range of the figure


class Crab:
    """Birth properties of the Crab pulsar under the paper's dipole model."""

    def __init__(self):
        self.omega = 2.0 * np.pi / PERIOD
        omega_dot = 2.0 * np.pi * PERIOD_DOT / PERIOD**2
        self.spin_down_luminosity = I_STAR * self.omega * omega_dot
        self.characteristic_age = PERIOD / (2.0 * PERIOD_DOT)

        # B_* from L_sd = 2 B^2 R^6 Omega^4 / 3c^3, Eq. (3) of the paper
        self.b_star = np.sqrt(
            3.0 * C_LIGHT**3 * self.spin_down_luminosity
            / (2.0 * R_STAR**6 * self.omega**4)
        )

        # (P0/P)^2 = 1 - T_age / tau_c for n = 3
        self.spin_ratio_squared = 1.0 - AGE / self.characteristic_age
        self.initial_period = PERIOD * np.sqrt(self.spin_ratio_squared)
        self.omega_0 = 2.0 * np.pi / self.initial_period
        self.rotational_energy = 0.5 * I_STAR * self.omega_0**2

        # Eqs. (5) and (9): the wind potential drop, then and now
        self.energy_max = self._potential_energy(self.omega_0)
        self.energy_today = self._potential_energy(self.omega)

        # rotational energy already radiated away
        self.energy_released = (1.0 - self.spin_ratio_squared) * self.rotational_energy

    def _potential_energy(self, omega):
        """Z e Delta Phi, in GeV, for Z = 1."""
        return E_CHARGE * self.b_star * R_STAR**3 * omega**2 / (2.0 * C_LIGHT**2) / GEV

    def summary(self):
        return "\n".join(
            [
                f"  L_sd        = {self.spin_down_luminosity:.2e} erg/s",
                f"  tau_c       = {self.characteristic_age / YEAR:.0f} yr",
                f"  B_*         = {self.b_star / 1e12:.2f} x 10^12 G",
                f"  P_0         = {self.initial_period * 1e3:.1f} ms",
                f"  E_rot       = {self.rotational_energy:.2e} erg"
                f"  ({self.energy_released:.2e} erg already released)",
                f"  E_M         = {self.energy_max / 1e6:.0f} PeV",
                f"  E_w(today)  = {self.energy_today / 1e6:.0f} PeV",
                f"  n_H         = {N_TARGET:.1f} cm^-3 ({M_EJECTA / M_SUN:.1f} Msun"
                f" in {R_NEBULA / PC:.1f} pc)",
            ]
        )


def inelastic_cross_section(energy):
    """p-p inelastic cross section, Kelner et al. (2006) Eq. (79). E in GeV."""
    threshold = 1.22
    logarithm = np.log(energy / 1.0e3)
    sigma = (34.3 + 1.88 * logarithm + 0.25 * logarithm**2) * MBARN
    return np.where(
        energy > threshold, sigma * (1.0 - (threshold / energy) ** 4) ** 2, 0.0
    )


def photon_yield(x, energy_proton):
    """F_gamma(x, E_p) of Kelner et al. (2006) Eq. (58), x = E_gamma / E_p."""
    logarithm = np.log(energy_proton / 1.0e3)
    b = 1.30 + 0.14 * logarithm + 0.011 * logarithm**2
    beta = 1.0 / (1.79 + 0.11 * logarithm + 0.008 * logarithm**2)
    k = 1.0 / (0.801 + 0.049 * logarithm + 0.014 * logarithm**2)

    x_beta = x**beta
    numerator = 1.0 - x_beta
    denominator = 1.0 + k * x_beta * numerator
    correction = (
        1.0 / np.log(x)
        - 4.0 * beta * x_beta / numerator
        - 4.0 * k * beta * x_beta * (1.0 - 2.0 * x_beta) / denominator
    )
    return b * np.log(x) / x * (numerator / denominator) ** 4 * correction


def pion_decay_sed(energy_gamma, energy_proton, spectrum, n_target):
    """E^2 dN/dE [erg cm^-2 s^-1] from protons `spectrum` [GeV^-1] at `energy_proton`."""
    log_proton = np.log(energy_proton)
    sigma = inelastic_cross_section(energy_proton)

    flux = np.empty_like(energy_gamma)
    for i, e_gamma in enumerate(energy_gamma):
        x = np.clip(e_gamma / energy_proton, 1e-12, 1.0 - 1e-9)
        integrand = np.where(
            e_gamma < energy_proton,
            sigma * spectrum * photon_yield(x, energy_proton),
            0.0,
        )
        flux[i] = C_LIGHT * n_target * np.trapezoid(integrand, log_proton)

    return flux * energy_gamma**2 * GEV / (4.0 * np.pi * DISTANCE**2)


def calorimetric_protons(crab, xi, n_points=400):
    """Every proton injected so far: dN/dE = xi E_rot / (E_M E), truncated below E_w(today).

    The truncation is the exact statement of Eq. (7): at time t the wind injects
    at E_w(t) only, so energies below the potential drop available today have not
    been populated yet. Integrating E dN/dE over the band returns the rotational
    energy already released, as it must.
    """
    energy = np.logspace(
        np.log10(crab.energy_today), np.log10(crab.energy_max), n_points
    )
    spectrum = xi * (crab.rotational_energy / GEV) / (crab.energy_max * energy)
    return energy, spectrum


def escape_limited_protons(crab, xi, residence_time, n_points=400):
    """Steady state: W_p = xi L_sd t_res, all of it at the energy injected today.

    Over a residence time this short the potential drop moves by less than a per
    cent, so the retained protons are monoenergetic to any relevant accuracy;
    they are spread over a narrow band here only to reuse the same integrator.
    """
    energy_budget = xi * crab.spin_down_luminosity * residence_time / GEV  # GeV
    lo, hi = 0.99 * crab.energy_today, 1.01 * crab.energy_today
    energy = np.logspace(np.log10(lo), np.log10(hi), n_points)
    spectrum = np.full_like(energy, 1.0) / energy  # dN/dE ∝ 1/E over the band
    spectrum *= energy_budget / np.trapezoid(spectrum * energy, energy)
    return energy, spectrum


def bohm_residence_time(energy, radius=R_NEBULA, field=B_NEBULA):
    """R^2 / 6D with D = r_L c / 3. E in GeV."""
    larmor = energy * GEV / (E_CHARGE * field)
    return radius**2 / (6.0 * larmor * C_LIGHT / 3.0)


def lhaaso_crab():
    """The measured LHAASO SED of the Crab, as released by the collaboration.

    Returns (detector, E [GeV], E^2 dN/dE, err_lo, err_up, kind) with the fluxes
    converted from TeV cm^-2 s^-1 to erg cm^-2 s^-1. `kind` is "point" except for
    the 1.6-2.5 PeV bin, which held no event and is published as a 90% CL upper
    limit; for that row the limit sits in `err_up` and the flux is the excess.
    """
    path = pathlib.Path(__file__).parent / "data" / "lhaaso_crab_sed.txt"
    table = np.genfromtxt(path, dtype=None, encoding="utf-8")
    detector = np.array([row[0] for row in table])
    kind = np.array([row[5] for row in table])
    energy = np.array([row[1] for row in table]) * 1e3  # TeV -> GeV
    flux, err_lo, err_up = (
        np.array([row[i] for row in table]) * 1e3 * GEV for i in (2, 3, 4)
    )
    return detector, energy, flux, err_lo, err_up, kind


def saturating_efficiency(energy_gamma, model_sed, xi_model, sigma=1.0):
    """The xi at which the model would first exceed the measurement.

    The hadronic component has to fit underneath a spectrum that is already
    accounted for leptonically, so the criterion is that it stay below every
    measured flux -- and below the published 90% CL limit in the bin where no
    photon was detected. `sigma` sets how much headroom the measurement is given:
    sigma = 1 allows the component to reach the upper end of each error bar and
    is the conservative choice; sigma = 0 holds it under the central values.
    Returns the efficiency and the energy of the point that sets it.
    """
    _, energy, flux, _, err_up, kind = lhaaso_crab()
    ceiling = np.where(kind == "limit", err_up, flux + sigma * err_up)

    model = np.interp(energy, energy_gamma, model_sed)
    inside = (energy >= energy_gamma[0]) & (energy <= energy_gamma[-1])
    ratio = np.where(inside, ceiling / model, np.inf)
    return xi_model * ratio.min(), energy[np.argmin(ratio)]


def plot_crab(figname="gryphon_crab_sed"):
    crab = Crab()
    print("Crab pulsar, magnetic-dipole reconstruction:")
    print(crab.summary())

    energy_gamma = np.logspace(np.log10(EMIN), np.log10(EMAX), 140)

    cases = {}
    energy, spectrum = calorimetric_protons(crab, XI)
    cases["calorimetric"] = pion_decay_sed(energy_gamma, energy, spectrum, N_TARGET)

    t_ballistic = R_NEBULA / C_LIGHT
    t_bohm = bohm_residence_time(crab.energy_today)
    for name, t_res in [("ballistic", t_ballistic), ("bohm", t_bohm)]:
        energy, spectrum = escape_limited_protons(crab, XI, t_res)
        cases[name] = pion_decay_sed(energy_gamma, energy, spectrum, N_TARGET)

    print(f"\nresidence times: R/c = {t_ballistic / YEAR:.1f} yr, "
          f"Bohm at {crab.energy_today / 1e6:.0f} PeV = {t_bohm / YEAR:.0f} yr")
    print(f"retained fraction of the injected energy: "
          f"{crab.spin_down_luminosity * t_ballistic / crab.energy_released:.1e} "
          f"to {crab.spin_down_luminosity * t_bohm / crab.energy_released:.1e}")

    print("\nefficiency saturating the LHAASO flux, as xi (n_H / 5 cm^-3):")
    print(f"  {'':13s} {'+1 sigma':>22s}   {'central value':>22s}")
    for name, sed in cases.items():
        row = [saturating_efficiency(energy_gamma, sed, XI, s) for s in (1.0, 0.0)]
        cells = "   ".join(
            f"{xi * N_TARGET / 5.0:8.2e} at {e / 1e6:5.2f} PeV" for xi, e in row
        )
        print(f"  {name:13s} {cells}")

    fraction = {
        name: np.interp(7.94e5, energy_gamma, sed) for name, sed in cases.items()
    }
    _, energy_d, flux_d, _, _, _ = lhaaso_crab()
    measured = np.interp(7.94e5, energy_d, flux_d)
    print("\nhadronic share of the measured flux at 0.79 PeV, at xi = 0.15:")
    for name, value in fraction.items():
        print(f"  {name:13s} {value / measured:6.2f}")

    # --- figure ----------------------------------------------------------
    fig, ax = gp.new_figure(figsize=(12.0, 8.5))
    gp.set_axes(
        ax,
        xlabel=r"E$_\gamma$ [GeV]",
        ylabel=r"E$_\gamma^2$ dN/dE$_\gamma$ [erg cm$^{-2}$ s$^{-1}$]",
        xscale="log",
        xlim=[EMIN, EMAX],
        yscale="log",
        ylim=[2e-16, 3e-11],
    )

    detector, energy, flux, err_lo, err_up, kind = lhaaso_crab()
    for name, marker in [("WCDA", "s"), ("KM2A", "o")]:
        shown = (detector == name) & (kind == "point") & (energy >= EMIN)
        if not shown.any():
            continue
        ax.errorbar(
            energy[shown], flux[shown],
            yerr=[err_lo[shown], err_up[shown]],
            fmt=marker, color="k", markeredgecolor="k", label=f"LHAASO-{name}",
            capsize=4.2, markersize=8.5, elinewidth=2.0, capthick=2.0, zorder=6,
        )

    for e, limit in zip(energy[kind == "limit"], err_up[kind == "limit"]):
        ax.errorbar(
            e, limit, yerr=0.35 * limit, uplims=True,
            fmt="none", color="k", elinewidth=2.0, capthick=2.0, zorder=6,
        )
        ax.hlines(limit, e / 1.26, e * 1.26, color="k", lw=2.0, zorder=6)

    ax.plot(
        energy_gamma, cases["calorimetric"], color="tab:blue", lw=4, zorder=5,
        label=r"$\xi = 0.15$, full retention",
    )
    ax.fill_between(
        energy_gamma, cases["ballistic"], cases["bohm"],
        color="tab:red", alpha=0.3, lw=0, zorder=3,
    )
    ax.plot(energy_gamma, cases["bohm"], color="tab:red", lw=3, zorder=4)
    ax.plot(
        energy_gamma, cases["ballistic"], color="tab:red", ls="--", lw=3, zorder=4,
        label=r"$\xi = 0.15$, escape limited",
    )

    ax.set_title(
        rf"$P_0 = {crab.initial_period * 1e3:.0f}$ ms, "
        rf"$E_M = {crab.energy_max / 1e6:.0f}$ PeV, "
        rf"$n_{{\rm H}} = {N_TARGET:.1f}$ cm$^{{-3}}$",
        fontsize=24,
    )
    ax.legend(fontsize=19, loc="lower center")
    gp.savefig(fig, figname)


if __name__ == "__main__":
    gp.style.use()
    plot_crab()
