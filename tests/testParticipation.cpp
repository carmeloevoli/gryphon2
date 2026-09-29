#include <limits>

#include "gryphon/core/participation.h"
#include "gtest/gtest.h"

using gryphon::core::Participation;

TEST(Participation, EmptyAndEqualSources) {
  Participation weights;
  EXPECT_DOUBLE_EQ(weights.effectiveCount(), 0.);
  weights.add(0.);
  EXPECT_DOUBLE_EQ(weights.effectiveCount(), 0.);
  for (int n = 1; n <= 100; ++n) {
    weights.add(1.);
    EXPECT_DOUBLE_EQ(weights.effectiveCount(), n);
    EXPECT_DOUBLE_EQ(weights.maxFraction(), 1. / n);
  }
}

TEST(Participation, UnequalSourcesAndExtremeRescaling) {
  for (double scale : {1e-200, 1., 1e200}) {
    Participation weights;
    for (double value : {1., 2., 3.}) weights.add(scale * value);
    EXPECT_NEAR(weights.effectiveCount(), 36. / 14., 1e-14);
    EXPECT_NEAR(weights.maxFraction(), .5, 1e-15);
  }
  Participation dominated;
  dominated.add(1.);
  dominated.add(1e-20);
  EXPECT_DOUBLE_EQ(dominated.effectiveCount(), 1.);
}

TEST(Participation, RejectsInvalidWeights) {
  Participation weights;
  EXPECT_THROW(weights.add(-1.), std::invalid_argument);
  EXPECT_THROW(weights.add(std::numeric_limits<double>::infinity()), std::invalid_argument);
  EXPECT_THROW(weights.add(std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
}
