#include <memory>

#include "gryphon.h"
#include "gtest/gtest.h"

namespace gryphon {

TEST(RNG, rangeUnity) {
  RandomNumberGenerator rng(1234);
  EXPECT_DOUBLE_EQ(rng.min(), 0.);
  EXPECT_DOUBLE_EQ(rng.max(), 1.);
}

TEST(RNG, reproducibleForSameSeed) {
  RandomNumberGenerator rng1(4242);
  RandomNumberGenerator rng2(4242);

  for (size_t i = 0; i < 32; ++i) {
    EXPECT_DOUBLE_EQ(rng1(), rng2());
  }
}

TEST(RNG, withinRange) {
  RandomNumberGenerator rng(5678);
  for (size_t i = 0; i < 100000; ++i) {
    {
      auto r = rng();
      EXPECT_GE(r, 0.0);
      EXPECT_LE(r, 1.0);
    }
    {
      auto r = rng.uniform(-0.1, 0.1);
      EXPECT_GE(r, -0.1);
      EXPECT_LE(r, 0.1);
    }
  }
}

TEST(RNG, uniformMomentsAreReasonable) {
  RandomNumberGenerator rng(12);
  const size_t N = 100000;
  double sum = 0.0;
  double sumSquared = 0.0;

  for (size_t i = 0; i < N; ++i) {
    const double value = rng();
    sum += value;
    sumSquared += std::pow(value - 0.5, 2.0);
  }

  EXPECT_NEAR(sum / static_cast<double>(N), 0.5, 0.003);
  EXPECT_NEAR(sumSquared / static_cast<double>(N - 1), 1. / 12., 0.003);
}

}  // namespace gryphon