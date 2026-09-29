"""Regression checks for the historical Figure 10 restyling."""

import unittest

import numpy as np

from feature_statistics import summarize
from plot_low_rate import check_legacy, survival


class LowRateTests(unittest.TestCase):
    def test_survival_includes_ties_and_preserves_zero_count_limits(self):
        data = np.array([0.01, 0.02, 0.02, 0.10])
        x = np.array([0.001, 0.02, 0.10, 0.11])
        empirical, lower, upper = survival(data, x)
        np.testing.assert_allclose(empirical[:3], [1, 0.75, 0.25])
        self.assertTrue(np.isnan(empirical[-1]))
        self.assertEqual(lower[-1], 0)
        self.assertGreater(upper[-1], 0)
        self.assertEqual(upper[0], 1)

    def test_replot_rejects_changed_statistics(self):
        data = np.array([0.01, 0.02, 0.08, 0.19])
        expected = summarize(data, (0.16, 0.19))
        check_legacy(data, expected)
        with self.assertRaisesRegex(ValueError, "sample size"):
            check_legacy(data[:-1], expected)
        with self.assertRaisesRegex(ValueError, "percentile"):
            check_legacy(data * 1.01, expected)
        expected["thresholds"]["0.16"]["exceedances"] = 0
        with self.assertRaisesRegex(ValueError, "exceedance count"):
            check_legacy(data, expected)


if __name__ == "__main__":
    unittest.main()
