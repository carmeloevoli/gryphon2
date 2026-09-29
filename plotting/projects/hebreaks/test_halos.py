import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

import gryphon_plots as gp
from feature_statistics import summarize
from plot_halos import HALOS, load_model, survival, validate_parameters, write_table


class HaloStatisticsTests(unittest.TestCase):
    def test_dedicated_configs_match_the_model(self):
        root = Path(__file__).resolve().parents[3]
        for halo in HALOS:
            params = gp.load_params(root / f"configs/hebreaks/halos_H{halo}.ini")
            # Explicit default used by the production code; full params.ini is
            # checked by the plotter after each simulation.
            params["timestep"] = 1.
            validate_parameters(params, halo)
            params["injemax"] = 1e7
            with self.assertRaisesRegex(ValueError, "injemax"):
                validate_parameters(params, halo)

    def test_survival_includes_ties_and_retains_zero_count_uncertainty(self):
        values = np.array([0.01, 0.02, 0.02, 0.1])
        empirical, lower, upper = survival(values, np.array([0.01, 0.02, 0.1, 0.2]))
        np.testing.assert_array_equal(empirical[:3], [1., 0.75, 0.25])
        self.assertTrue(np.isnan(empirical[-1]))
        self.assertEqual(lower[-1], 0.)
        self.assertGreater(upper[-1], 0.)
        self.assertTrue(np.all(np.diff(upper) <= 0))

    def test_loader_requires_completed_matching_checksum_and_seed(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "halos_H2"
            directory.mkdir()
            config = "simname = halos_H2\n"
            (directory / "halos_manifest.json").write_text(json.dumps({
                "config_text": config,
                "config_sha256": hashlib.sha256(config.encode()).hexdigest(),
            }))
            with self.assertRaisesRegex(ValueError, "not complete"):
                load_model(directory, 1)
            flux = directory / "flux_000000.txt"
            energy = np.geomspace(1e3, 1e6, 48)
            np.savetxt(flux, np.column_stack((energy, energy ** -2.7)),
                       header="seed = 0\nsimname = halos_H2\ncolumns: E [GeV] | I")
            marker = directory / "complete_000000.json"
            marker.write_text(json.dumps({"flux_sha256": hashlib.sha256(flux.read_bytes()).hexdigest()}))
            grid, spectra, _, _ = load_model(directory, 1)
            np.testing.assert_allclose(grid, energy)
            np.testing.assert_allclose(spectra[0], energy ** -2.7)
            population = directory / "population_000000.txt"
            with self.assertRaisesRegex(ValueError, "missing population"):
                load_model(directory, 1, extra_stems=("population",))
            population.write_text("population metadata\n")
            marker.write_text(json.dumps({
                "flux_sha256": hashlib.sha256(flux.read_bytes()).hexdigest(),
                "population_sha256": hashlib.sha256(population.read_bytes()).hexdigest(),
            }))
            load_model(directory, 1, extra_stems=("population",))
            population.write_text("changed population metadata\n")
            with self.assertRaisesRegex(ValueError, "checksum"):
                load_model(directory, 1, extra_stems=("population",))
            marker.write_text(json.dumps({"flux_sha256": hashlib.sha256(flux.read_bytes()).hexdigest()}))
            flux.write_text(flux.read_text().replace("seed = 0", "seed = 1"))
            with self.assertRaisesRegex(ValueError, "checksum"):
                load_model(directory, 1)
            marker.write_text(json.dumps({"flux_sha256": hashlib.sha256(flux.read_bytes()).hexdigest()}))
            with self.assertRaisesRegex(ValueError, "header mismatch"):
                load_model(directory, 1)

    def test_generated_table_contains_the_complete_revtex_alignment(self):
        models = {str(h): {"residual": summarize(np.array([0.01, 0.2]), (0.16, 0.19)),
                           "slope_excursion": summarize(np.array([0.02, 0.3]), (0.2, 0.35))}
                  for h in HALOS}
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "table.tex"
            write_table(path, models)
            text = path.read_text()
        self.assertIn(r"\begin{tabular}{ccccccc}", text)
        self.assertTrue(text.endswith("\\end{tabular}\n"))
        for halo in HALOS:
            self.assertIn(f"{halo} & 19.05\\% & 0.2860 & 1 & 1 & 1 & 0", text)


if __name__ == "__main__":
    unittest.main()
