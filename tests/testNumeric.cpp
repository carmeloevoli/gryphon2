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

TEST(HaloFunction, IsSmallAtHaloBoundary) {
  const double H = 4.0 * cgs::kpc;
  const double value = utils::halo_function(std::pow(0.8 * H, 2), H, H, 0.0);

  EXPECT_GE(value, 0.0);
  EXPECT_LT(value, 1e-8);
}

}  // namespace gryphon
