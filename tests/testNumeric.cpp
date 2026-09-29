#include <algorithm>
#include <cmath>

#include "gryphon.h"
#include "gtest/gtest.h"

namespace gryphon {

TEST(HaloFunction, ReturnsZeroForInvalidParameters) {
  EXPECT_DOUBLE_EQ(utils::halo_function(0.0, 1.0, 0.0, 0.0), 0.0);
  EXPECT_DOUBLE_EQ(utils::halo_function(1.0, 0.0, 0.0, 0.0), 0.0);
}

TEST(HaloFunction, IsSymmetricInZAndZs) {
  const double H = 4.0 * cgs::kpc;
  const double l2 = std::pow(0.7 * H, 2);
  const double lhs = utils::halo_function(l2, H, 0.4 * H, -0.15 * H);
  const double rhs = utils::halo_function(l2, H, -0.15 * H, 0.4 * H);
  const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});

  EXPECT_NEAR(lhs, rhs, 1e-12 * scale);
}

TEST(HaloFunction, StaysFiniteAndNonNegativeForLargeDiffusionLength) {
  const double H = 4.0 * cgs::kpc;
  const double value = utils::halo_function(std::pow(20.0 * H, 2), H, 0.3 * H, -0.1 * H);

  EXPECT_TRUE(std::isfinite(value));
  EXPECT_GE(value, 0.0);
}

TEST(HaloFunction, DecreasesAwayFromSource) {
  const double H = 4.0 * cgs::kpc;
  const double l2 = std::pow(0.8 * H, 2);
  const double onSource = utils::halo_function(l2, H, 0.0, 0.0);
  const double offSource = utils::halo_function(l2, H, 0.7 * H, 0.0);

  EXPECT_TRUE(std::isfinite(onSource));
  EXPECT_TRUE(std::isfinite(offSource));
  EXPECT_GE(onSource, offSource);
}

namespace {

// Brute-force image sum, the definition of the halo function. It is accurate in
// double precision as long as lambda/H stays below ~5, which is enough to pin
// down the eigenmode branch just above the switch at lambda/H = 3.
double haloFunctionByImages(double l2, double H, double z, double zs, size_t images = 4000) {
  auto gaussian = [l2](double x) { return std::exp(-x * x / l2); };

  double value = gaussian(z - zs);
  for (size_t n = 1; n <= images; ++n) {
    const double sign = (n % 2 == 0) ? 1. : -1.;
    const double nn = static_cast<double>(n);
    value += sign * (gaussian(z - (sign * zs + 2. * nn * H)) +
                     gaussian(z - (sign * zs - 2. * nn * H)));
  }
  return value;
}

}  // namespace

TEST(HaloFunction, EigenmodeBranchMatchesImageSum) {
  const double H = 4.0 * cgs::kpc;
  const double zs = 0.05 * cgs::kpc;

  // Above lambda/H = 3 the implementation switches to the eigenmode expansion;
  // it must reproduce the image sum it replaces.
  for (const double ratio : {3.2, 3.5, 4.0, 4.5}) {
    const double l2 = std::pow(ratio * H, 2);
    const double expected = haloFunctionByImages(l2, H, 0., zs);
    const double actual = utils::halo_function(l2, H, 0., zs);

    EXPECT_NEAR(actual, expected, 1e-6 * expected) << "at lambda/H = " << ratio;
  }
}

TEST(HaloFunction, IsContinuousAcrossRepresentationSwitch) {
  const double H = 4.0 * cgs::kpc;
  const double zs = 0.05 * cgs::kpc;

  const double below = utils::halo_function(std::pow(2.9999 * H, 2), H, 0., zs);
  const double above = utils::halo_function(std::pow(3.0001 * H, 2), H, 0., zs);

  // The two representations are the same function, so the switch must not be
  // visible in the result. The tolerance leaves room for the genuine variation
  // of the halo function between the two sampling points, which falls by about
  // 10% per per cent in lambda in this regime, and is still far tighter than
  // any constant mismatch between the two branches.
  EXPECT_NEAR(above, below, 5e-3 * below);
}

TEST(HaloFunction, IsSmallAtHaloBoundary) {
  const double H = 4.0 * cgs::kpc;
  const double value = utils::halo_function(std::pow(0.8 * H, 2), H, H, 0.0);

  EXPECT_GE(value, 0.0);
  EXPECT_LT(value, 1e-8);
}

TEST(HaloFunction, AnalyticDerivativeMatchesFiniteDifference) {
  const double H = 4.0 * cgs::kpc;
  const double z = 0.2 * H;
  const double zs = -0.15 * H;
  const double step = 1e-5 * H;

  for (const double ratio : {0.8, 3.5}) {
    const double l2 = std::pow(ratio * H, 2);
    const double finiteDifference =
        (utils::halo_function(l2, H, z + step, zs) -
         utils::halo_function(l2, H, z - step, zs)) /
        (2. * step);
    const double analytic = utils::halo_function_dz(l2, H, z, zs);
    const double scale = std::max(1. / H, std::abs(finiteDifference));

    EXPECT_NEAR(analytic, finiteDifference, 2e-6 * scale) << "at lambda/H = " << ratio;
  }
}

}  // namespace gryphon
