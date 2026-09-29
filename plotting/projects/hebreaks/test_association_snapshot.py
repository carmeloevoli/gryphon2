"""Selection and consistency tests without requiring simulation outputs."""

from copy import deepcopy
from pathlib import Path
import unittest

import numpy as np

import gryphon_plots as gp
from plot_association_snapshot import (
    highlighted_parents, validate_events, validate_pair, window_mask,
)


def table(rows):
    names = ["age", "x", "y", "z", "parent_id", "birth_age", "delay",
             "centre_x", "centre_y", "centre_z"]
    return gp.Table(Path("test.txt"), np.asarray(rows, dtype=float), names)


def event(parent, x, y, centre_x, age=.5, delay=10.):
    return [age, x, y, 0., parent, age + delay, delay, centre_x, 0., 0.]


class SnapshotTests(unittest.TestCase):
    def test_same_projected_window_includes_edges(self):
        data = {"x": np.array([-2., 0., 2., 2.01, 0]),
                "y": np.array([0., 2., -2., 0., -2.01])}
        np.testing.assert_array_equal(window_mask(data), [True, True, True, False, False])

    def test_parent_selection_by_distance_not_richness(self):
        data = table([event(7, .1, 0, 8.6), event(7, .2, 0, 8.6),
                      event(3, 1, 0, 9.5), event(3, 1.1, 0, 9.5), event(3, 1.2, 0, 9.5),
                      event(2, 0, 0, 8.5),  # only one visible event
                      event(2, 3, 0, 8.5),
                      event(0, .1, 0, 8.5), event(0, .1, 0, 8.5),  # field
                      event(4, 1.9, 0, 11), event(4, 1.9, 0, 11)])  # centre outside
        self.assertEqual(highlighted_parents(data, 8.5), [7, 3])
        self.assertEqual(highlighted_parents(data, 8.5, number=1), [7])
        data.values = data.values[::-1]
        self.assertEqual(highlighted_parents(data, 8.5), [7, 3])

    def test_parameters_must_match(self):
        base = dict(simname="a", associationfraction=0, syntheticassociations="true",
                    spiralmodel="Xie2024", seed=0, maxtimemyr=10, snrateyr=.02,
                    associationmembers=100, associationradiuspc=30,
                    associationvelocitykms=3, sunkpc=8.5)
        pair = [base, dict(base, simname="b", associationfraction=1)]
        validate_pair(pair)
        bad = deepcopy(pair)
        bad[1]["associationradiuspc"] = 0
        with self.assertRaises(ValueError):
            validate_pair(bad)
        for wrong_window in (1, 30):
            bad = deepcopy(pair)
            for p in bad:
                p["maxtimemyr"] = wrong_window
            with self.subTest(window=wrong_window), self.assertRaises(ValueError):
                validate_pair(bad)

    def test_physical_parent_metadata(self):
        p = dict(maxtimemyr=10, associationfraction=1)
        valid = table([event(1, .1, 0, 8.5), event(1, .2, 0, 8.5)])
        validate_events(valid, p)
        # A retained explosion may have a parent much older than the window.
        old = table([event(1, .1, 0, 8.5, age=9.5, delay=48)])
        validate_events(old, p)
        for column, value in (("age", 10.1), ("birth_age", 1), ("delay", 50),
                              ("parent_id", .5), ("parent_id", 0),
                              ("centre_x", 9), ("x", np.nan)):
            bad = deepcopy(valid)
            bad.values[0, bad.index(column)] = value
            with self.subTest(column=column, value=value), self.assertRaises(ValueError):
                validate_events(bad, p)


if __name__ == "__main__":
    unittest.main()
