from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.integrate import quad
from scipy.stats import norm

import gryphon_plots as gp
from feature_statistics import excursions, running_indices, summarize
from plot_variations import (MODELS, continuous_index_spectrum, scaled_spectra,
                             shared_parameters, validate_parameters, write_table)


class VariationTests(unittest.TestCase):
    def test_configs_change_only_the_requested_injection_properties(self):
        root = Path(__file__).resolve().parents[3]
        reference = None
        for model, (run, _, _) in MODELS.items():
            params = gp.load_params(root / f"configs/hebreaks/{run}.ini")
            params["timestep"] = 1.
            validate_parameters(params, model)
            if reference is None:
                reference = shared_parameters(params)
            self.assertEqual(shared_parameters(params), reference)
            if model == "index":
                params["injslopesigma"] = 0.15
                with self.assertRaisesRegex(ValueError, "injslopesigma"):
                    validate_parameters(params, model)
            params["varyenergy"] = "true" if params["varyenergy"] == "false" else "false"
            with self.assertRaisesRegex(ValueError, "varyenergy"):
                validate_parameters(params, model)

    def test_display_normalization_is_one_scalar_and_statistics_invariant(self):
        energy = np.geomspace(1e3, 1e6, 48)
        flux = np.array([energy**-2.7, 3*energy**-2.6, 0.5*energy**-2.8])
        curves, factor = scaled_spectra(energy, flux)
        self.assertAlmostEqual(np.median(curves[:, 0]), 1.)
        np.testing.assert_allclose(curves / (flux * energy**2.7), factor)
        np.testing.assert_allclose(curves[1] / curves[0], flux[1] / flux[0])
        np.testing.assert_allclose(excursions(energy, flux*factor), excursions(energy, flux), atol=2e-13)

    def test_shared_running_index_matches_slope_statistic(self):
        energy = np.geomspace(1e3, 1e6, 48)
        flux = energy**-2.7 * np.exp(0.02*np.log(energy/1e3)**2)
        centres, indices = running_indices(energy, flux)
        np.testing.assert_array_equal(centres, energy[2:-2])
        np.testing.assert_allclose(indices[0], 2.7-0.04*np.log(centres/1e3), atol=1e-12)
        mask = (centres >= 1e4) & (centres <= 1e5)
        np.testing.assert_allclose(np.ptp(indices[:, mask], axis=1), excursions(energy, flux)[1])

    def test_continuous_index_limit_is_energy_normalized_and_truncated(self):
        energy = np.geomspace(1e3, 1e6, 9)
        mu, sigma = 2.34, 0.30
        numerical = []
        for e in energy:
            y = np.log(e/10.)
            # Remove exp(-2*y) during quadrature for good relative precision.
            integral = quad(lambda u: u*norm.pdf((u+2-mu)/sigma)/sigma*np.exp(-u*y),
                            0., np.inf, epsabs=1e-12, epsrel=1e-11)[0]
            numerical.append(integral / norm.cdf((mu-2)/sigma) * np.exp(-2*y) * (e/1e3)**-0.36)
        np.testing.assert_allclose(continuous_index_spectrum(energy), numerical, rtol=1e-10, atol=0.)

    def test_generated_table_has_complete_alignment(self):
        models = {name: {"residual": summarize(np.array([0.01, 0.2]), (0.16, 0.19)),
                         "slope_excursion": summarize(np.array([0.02, 0.3]), (0.20, 0.35))}
                  for name in MODELS}
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "table.tex"
            write_table(path, models)
            text = path.read_text()
        self.assertIn(r"\begin{tabular}{lcccccc}", text)
        self.assertTrue(text.endswith("\\end{tabular}\n"))
        for _, label, _ in MODELS.values():
            self.assertIn(f"{label} & 19.05\\% & 0.2860 & 1 & 1 & 1 & 0", text)


if __name__ == "__main__":
    unittest.main()
