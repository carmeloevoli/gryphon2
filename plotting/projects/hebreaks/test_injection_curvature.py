"""Verify the paper's exact uncut, energy-normalized truncated-Gaussian average."""

import unittest

import numpy as np
from scipy.integrate import quad
from scipy.stats import norm


class InjectionCurvatureTests(unittest.TestCase):
    def test_truncated_average_and_curvature_against_quadrature(self):
        mu, sigma = 2.34, 0.30
        a = (mu - 2) / sigma
        for y in (0., np.log(10), np.log(1e3), np.log(1e5)):
            with self.subTest(y=y):
                z = a - sigma * y
                b = norm.pdf(z) + z * norm.cdf(z)
                # Divide out exp(-2y) to keep quadrature well scaled. Integrate
                # u=gamma-2>0; the irrelevant E_CR/E_min^2 factor is also omitted.
                def weight(u):
                    return u * norm.pdf((u + 2 - mu) / sigma) / sigma * np.exp(-u * y)
                moments = [quad(lambda u: u**n * weight(u), 0., np.inf,
                                epsabs=1e-12, epsrel=1e-11)[0] for n in range(3)]
                numerical = moments[0] / norm.cdf(a)
                analytic = sigma * b / norm.cdf(a) * np.exp(
                    -(mu - 2) * y + 0.5 * sigma**2 * y**2)
                self.assertAlmostEqual(analytic / numerical, 1., places=10)
                variance = moments[2] / moments[0] - (moments[1] / moments[0])**2
                curvature = sigma**2 * (1 + norm.pdf(z) / b - (norm.cdf(z) / b)**2)
                self.assertGreater(curvature, 0.)
                self.assertAlmostEqual(curvature, variance, places=11)
                # Independent numerical second derivative of the log average.
                def log_average(offset):
                    zz = a - sigma * offset
                    bb = norm.pdf(zz) + zz * norm.cdf(zz)
                    return np.log(bb) - mu * offset + 0.5 * sigma**2 * offset**2
                step = 1e-3
                difference = (log_average(y + step) - 2 * log_average(y)
                              + log_average(y - step)) / step**2
                self.assertAlmostEqual(curvature, difference, places=7)


if __name__ == "__main__":
    unittest.main()
