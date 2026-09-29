from pathlib import Path
import unittest

import numpy as np
from gryphon_plots.io import Table
from plot_effective_sources import validate_table


class EffectiveSourcesTests(unittest.TestCase):
    def fixture(self):
        energy = np.geomspace(1e3, 1e6, 48)
        flux = energy**-2.7
        table = Table(Path("diagnostic.txt"), np.column_stack((energy, flux, np.full(48, 2.), np.full(48, .5))),
                      names=["E", "I", "N_eff", "q_max"], meta={"seed": "0", "sourcecount": "3"})
        original = Table(Path("original.txt"), np.column_stack((energy, flux)), names=["E", "I"])
        return table, original

    def test_matched_total_and_bounded_participation(self):
        table, original = self.fixture()
        self.assertEqual(validate_table(table, original), 0.)

    def test_reject_different_seed_or_total(self):
        table, original = self.fixture()
        table.meta["seed"] = "1"
        with self.assertRaisesRegex(ValueError, "seed 0"):
            validate_table(table, original)
        table.meta["seed"] = "0"
        table.values[:, 1] *= 1.01
        with self.assertRaisesRegex(ValueError, "original seed-0"):
            validate_table(table, original)

    def test_reject_impossible_effective_count(self):
        for value in (.5, 4., np.nan):
            table, original = self.fixture()
            table.values[:, 2] = value
            with self.assertRaisesRegex(ValueError, "invalid participation"):
                validate_table(table, original)
        table, original = self.fixture()
        table.values[:, 3] = .9
        with self.assertRaisesRegex(ValueError, "largest source"):
            validate_table(table, original)


if __name__ == "__main__":
    unittest.main()
