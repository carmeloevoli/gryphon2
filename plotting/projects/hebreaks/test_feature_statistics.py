import unittest

import numpy as np

from feature_statistics import binomial_interval, excursions, summarize


class FeatureStatisticsTests(unittest.TestCase):
    def setUp(self):
        self.energy = np.geomspace(1e3, 1e6, 48)

    def test_power_law_and_normalization_invariance(self):
        flux = self.energy ** -2.7
        residual, slope = excursions(self.energy, np.stack((flux, 1e8 * flux)))
        np.testing.assert_allclose(residual, 0, atol=2e-13)
        np.testing.assert_allclose(slope, 0, atol=2e-12)

    def test_log_parabola_running_slope(self):
        x = np.log(self.energy / 1e3)
        kappa = 0.08
        flux = np.exp(-2.7 * x + 0.5 * kappa * x ** 2)
        _, slope = excursions(self.energy, flux)
        selected = x[(self.energy >= 1e4) & (self.energy <= 1e5)]
        self.assertAlmostEqual(slope[0], kappa * np.ptp(selected), places=11)
        residual2, slope2 = excursions(self.energy, 1000 * flux)
        residual, _ = excursions(self.energy, flux)
        np.testing.assert_allclose(residual2, residual, atol=1e-12)
        np.testing.assert_allclose(slope2, slope, atol=1e-12)

    def test_zero_exceedances_are_an_upper_limit(self):
        result = summarize(np.zeros(256), (0.16,))["thresholds"]["0.16"]
        self.assertEqual(result["exceedances"], 0)
        self.assertAlmostEqual(result["upper95_one_sided"], 1 - 0.05 ** (1 / 256))
        lo, hi = binomial_interval(0, 256)
        self.assertEqual(lo, 0)
        self.assertGreater(hi, result["upper95_one_sided"])

    def test_edge_flux_changes_do_not_enter_running_index_search(self):
        flux = self.energy ** -2.7
        flux[:2] *= 20
        flux[-2:] *= 0.1
        residual, slope = excursions(self.energy, flux)
        self.assertGreater(residual[0], 0.01)  # The global reference fit does change.
        self.assertLess(slope[0], 2e-12)  # No boundary derivative leaks into the search.

    def test_zero_count_limit_for_the_thousand_realization_first_pass(self):
        result = summarize(np.zeros(1000), (0.20,))["thresholds"]["0.2"]
        self.assertAlmostEqual(result["upper95_one_sided"], 0.0029912495450953314)

    def test_invalid_spectra_fail(self):
        with self.assertRaises(ValueError):
            excursions(self.energy, np.zeros(48))
        with self.assertRaises(ValueError):
            excursions(self.energy[::-1], np.ones(48))

    def test_wide_window_matches_analytic_curvature_and_nested_search(self):
        energy = np.geomspace(100., 1e6, 65)
        x = np.log(energy / 1e3)
        kappa = .08
        flux = np.exp(-2.7*x + .5*kappa*x*x)
        a, d = excursions(energy, flux, fit_range=(100., 1e6), search_range=(1e3, 1e5))
        a2, d2 = excursions(energy, flux, fit_range=(100., 1e6), search_range=(1e4, 1e5))
        self.assertAlmostEqual(d[0], kappa * np.log(100.), places=11)
        self.assertAlmostEqual(d2[0], kappa * np.log(10.), places=11)
        self.assertGreaterEqual(a[0], a2[0])
        self.assertGreaterEqual(d[0], d2[0])
        aa, dd = excursions(energy, flux * 1e35, fit_range=(100., 1e6), search_range=(1e3, 1e5))
        np.testing.assert_allclose(aa, a, atol=1e-12)
        np.testing.assert_allclose(dd, d, atol=1e-12)

    def test_wide_window_rejects_missing_coverage_or_slope_padding(self):
        energy = np.geomspace(100., 1e6, 65)
        flux = energy**-2.7
        with self.assertRaisesRegex(ValueError, "cover"):
            excursions(energy[1:], flux[1:], fit_range=(100., 1e6), search_range=(1e3, 1e5))
        with self.assertRaisesRegex(ValueError, "padding"):
            excursions(energy, flux, fit_range=(100., 1e6), search_range=(100., 1e5))
        with self.assertRaisesRegex(ValueError, "contained"):
            excursions(energy, flux, fit_range=(100., 1e6), search_range=(1e3, 2e6))


if __name__ == "__main__":
    unittest.main()
