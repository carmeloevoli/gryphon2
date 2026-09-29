"""Checks of the Figure 5 parameters and display-only normalization."""

from pathlib import Path
import unittest

import numpy as np
import gryphon_plots as gp

from plot_realizations import EXPECTED, TARGET, scaled_spectra, validate_params


class RealizationTests(unittest.TestCase):
    def test_common_factor_preserves_relative_amplitudes_and_shapes(self):
        energy = np.geomspace(100., 1e6, 64)
        flux = np.array([energy**-2.7, 2 * energy**-2.6, 0.5 * energy**-2.8])
        curves, median, factor = scaled_spectra(energy, flux)
        np.testing.assert_allclose(curves, factor * flux * energy**2.7)
        np.testing.assert_allclose(curves[1] / curves[0], flux[1] / flux[0])
        self.assertAlmostEqual(np.mean(median), TARGET)

    def test_invalid_spectra_rejected(self):
        for flux in ([[1., -1.]], [[1., np.nan]], [[1., np.inf]]):
            with self.assertRaises(ValueError):
                scaled_spectra([100., 1000.], flux)

    def test_wrong_cutoff_or_scatter_rejected(self):
        params = dict(EXPECTED, hkpc=2.)
        validate_params(params, 2)
        for key, value in (("injemax", 1e6), ("varyslope", "true"), ("d0h", 0.2)):
            with self.assertRaises(ValueError):
                validate_params(dict(params, **{key: value}), 2)

    def test_configs_only_differ_in_halo_and_name(self):
        root = Path(__file__).resolve().parents[3]
        params = [gp.load_params(root / f"configs/hebreaks/figure5_H{h}.ini")
                  for h in (2, 4)]
        for h, p in zip((2, 4), params):
            validate_params(p, h)
            p.pop("hkpc")
            p.pop("simname")
        self.assertEqual(*params)


if __name__ == "__main__":
    unittest.main()
