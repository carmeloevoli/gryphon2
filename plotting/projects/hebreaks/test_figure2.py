"""Analytic and coordinate checks for Figure 2 (no simulation files needed)."""

import unittest

import numpy as np

from figure2_common import escape_time_myr, galactocentric_xy, horizon_kpc


class Figure2Tests(unittest.TestCase):
    def setUp(self):
        self.params = dict(d0h=0.42, e0=1000.0, delta=0.36, sunkpc=8.5)

    def test_escape_time_normalization(self):
        np.testing.assert_allclose(
            [escape_time_myr(1000, h, self.params) for h in (2, 4, 8)],
            [2.380952380952381, 4.761904761904762, 9.523809523809524],
        )

    def test_escape_energy_and_halo_scaling(self):
        energy = np.geomspace(1e3, 1e6, 100)
        np.testing.assert_allclose(escape_time_myr(energy, 8, self.params),
                                   2 * escape_time_myr(energy, 4, self.params))
        self.assertAlmostEqual(
            escape_time_myr(1e6, 4, self.params) / escape_time_myr(1e3, 4, self.params),
            1000**-0.36,
        )

    def test_horizon_consistent_with_escape_time(self):
        for h in (2, 4, 8):
            diffusion = self.params["d0h"] * h
            self.assertAlmostEqual(horizon_kpc(h),
                                   np.sqrt(4 * diffusion * escape_time_myr(1000, h, self.params)))

    def test_sun_centred_to_galactocentric_coordinates(self):
        events = {"x": np.array([-8.5, 0, 2]), "y": np.array([0, 0, 3])}
        x, y = galactocentric_xy(events, self.params)
        np.testing.assert_allclose(x, [0, 8.5, 10.5])
        np.testing.assert_allclose(y, events["y"])


if __name__ == "__main__":
    unittest.main()
