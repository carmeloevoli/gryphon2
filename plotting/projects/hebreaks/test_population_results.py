from pathlib import Path
import tempfile
import unittest

import numpy as np
import gryphon_plots as gp

from feature_statistics import summarize
from plot_population_results import validate_parameters, shared_parameters, write_table


class PopulationResultsTests(unittest.TestCase):
    def params(self, name):
        root = Path(__file__).resolve().parents[3]
        params = gp.load_params(root / "configs/hebreaks/halos_H4.ini")
        params.update(gp.load_params(root / f"configs/hebreaks/{name}.ini"))
        # The overnight runner explicitly extends the original association
        # configs to 100 Myr; plotting validates the resolved production file.
        params.update(timestep=1., maxtimemyr=100.)
        return params

    def test_association_controls_and_current_physics(self):
        matched = []
        for name in ("independent", "mixed", "clustered"):
            params = self.params(f"associations_{name}")
            validate_parameters(params, "association", name)
            matched.append(shared_parameters(params, "association"))
            for key, value in (("maxtimemyr", 30.), ("injemax", 1e12),
                               ("associationmembers", 10.), ("efficiency", 1.)):
                with self.subTest(model=name, key=key):
                    with self.assertRaises(ValueError):
                        validate_parameters(dict(params, **{key: value}), "association", name)
        self.assertEqual(matched[0], matched[1])
        self.assertEqual(matched[0], matched[2])

    def test_rate_and_energy_are_both_required(self):
        params = self.params("halos_H4")
        validate_parameters(params, "low_rate", "reference")
        rare = dict(params, snrateyr=.002, efficiency=1.)
        validate_parameters(rare, "low_rate", "rare")
        self.assertEqual(shared_parameters(params, "low_rate"),
                         shared_parameters(rare, "low_rate"))
        for incorrect in (params, dict(rare, efficiency=.1), dict(rare, snrateyr=.02)):
            with self.assertRaises(ValueError):
                validate_parameters(incorrect, "low_rate", "rare")

    def test_table_contains_counts_not_probabilities(self):
        models = {"reference": {"residual": summarize(np.array([.01, .2]), (.16, .19)),
                                 "slope_excursion": summarize(np.array([.02, .3]), (.2, .35))}}
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "table.tex"
            write_table(path, "low_rate", models)
            result = path.read_text()
        self.assertIn("Reference & 19.05\\% & 0.2860 & 1 & 1 & 1 & 0", result)
        self.assertTrue(result.endswith("\\end{tabular}\n"))


if __name__ == "__main__":
    unittest.main()
