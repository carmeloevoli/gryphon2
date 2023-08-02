#include <memory>

#include "gryphon.h"
#include "gtest/gtest.h"

namespace gryphon {

TEST(RNG, rangeUnity) {
  RandomNumberGenerator rng = utils::RNG<double>(1234);
  EXPECT_DOUBLE_EQ(rng.min(), 0.);
  EXPECT_DOUBLE_EQ(rng.max(), 1.);
}

TEST(RNG, withinRange) {
  RandomNumberGenerator rng = utils::RNG<double>(5678);
  for (size_t i = 0; i < 1000000; ++i) {
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

TEST(RNG, uniform) {
  RandomNumberGenerator rng = utils::RNG<double>(12);
  size_t N = 1000000;
  double sum = 0;
  double sumSquared = 0;

  for (size_t i = 0; i < N; ++i) {
    sum += rng();
    sumSquared += std::pow(rng() - 0.5, 2.0);
  }

  EXPECT_NEAR(sum / (double)N, 0.5, 0.001);
  EXPECT_NEAR(sumSquared / (double)(N - 1), 1. / 12., 0.001);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

}  // namespace gryphon