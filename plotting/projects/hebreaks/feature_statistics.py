"""The manuscript's ideal-spectrum excursion statistics, not data-fit p-values."""

from __future__ import annotations

import numpy as np
from scipy.stats import beta


def _validate_spectra(energy: np.ndarray, flux: np.ndarray):
    energy = np.asarray(energy, dtype=float)
    flux = np.atleast_2d(np.asarray(flux, dtype=float))
    if energy.ndim != 1 or flux.shape[1] != energy.size or energy.size < 5:
        raise ValueError("need a one-dimensional energy grid and matching spectra with >=5 bins")
    if not np.all(np.isfinite(energy)) or not np.all(np.isfinite(flux)):
        raise ValueError("nonfinite spectrum")
    if np.any(energy <= 0) or np.any(np.diff(energy) <= 0) or np.any(flux <= 0):
        raise ValueError("energies must increase and energies/fluxes must be positive")
    return energy, flux


def running_indices(energy: np.ndarray, flux: np.ndarray):
    """Five-bin centered indices for every spectrum; omit both two-bin edges."""
    energy, flux = _validate_spectra(energy, flux)
    x, y = np.log(energy / 1e3), np.log(flux)
    slopes = []
    for i in range(2, energy.size - 2):
        dx = x[i - 2:i + 3] - np.mean(x[i - 2:i + 3])
        window = y[:, i - 2:i + 3]
        dy = window - window.mean(axis=1, keepdims=True)
        slopes.append(-(dy @ dx) / (dx @ dx))
    return energy[2:-2], np.stack(slopes, axis=1)


def excursions(energy: np.ndarray, flux: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    energy, flux = _validate_spectra(energy, flux)
    fit = (energy >= 1e3) & (energy <= 1e6)
    search = (energy >= 1e4) & (energy <= 1e5)
    if fit.sum() < 5 or search.sum() < 2:
        raise ValueError("insufficient fit/search coverage")
    x, y = np.log(energy / 1e3), np.log(flux)
    design = np.column_stack((np.ones(fit.sum()), x[fit]))
    coefficients = np.linalg.lstsq(design, y[:, fit].T, rcond=None)[0]
    fitted = coefficients[0, :, None] + coefficients[1, :, None] * x
    residual = np.max(np.abs(np.expm1(y[:, search] - fitted[:, search])), axis=1)
    centres, slopes = running_indices(energy, flux)
    mask = (centres >= 1e4) & (centres <= 1e5)
    if mask.sum() < 2:
        raise ValueError("insufficient interior bins for slope statistic")
    slope_excursion = np.ptp(slopes[:, mask], axis=1)
    return residual, slope_excursion


def binomial_interval(count, size: int, confidence: float = 0.95):
    """Exact two-sided Clopper-Pearson interval, including zero/all successes."""
    k = np.asarray(count)
    if size < 1 or np.any(k < 0) or np.any(k > size):
        raise ValueError("invalid binomial counts")
    tail = (1. - confidence) / 2.
    lower = np.where(k == 0, 0., beta.ppf(tail, k, size - k + 1))
    upper = np.where(k == size, 1., beta.ppf(1. - tail, k + 1, size - k))
    return lower, upper


def summarize(values: np.ndarray, thresholds: tuple[float, ...]) -> dict:
    values = np.asarray(values)
    n = len(values)
    if n < 1 or not np.all(np.isfinite(values)):
        raise ValueError("empty/nonfinite statistic")
    result = {"number_of_realizations": n,
              "percentiles": {str(q): float(np.percentile(values, q)) for q in (50, 90, 95, 99)},
              "maximum": float(values.max()), "thresholds": {}}
    for threshold in thresholds:
        count = int(np.count_nonzero(values >= threshold))
        lo, hi = binomial_interval(count, n)
        entry = {"exceedances": count, "probability": count / n,
                 "interval95_two_sided": [float(lo), float(hi)]}
        if count == 0:
            entry["upper95_one_sided"] = float(-np.expm1(np.log(0.05) / n))
        result["thresholds"][str(threshold)] = entry
    return result
