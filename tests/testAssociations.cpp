#include <cmath>
#include <limits>
#include <map>
#include <vector>

#include "gryphon.h"
#include "gryphon/galaxy/associations.h"
#include "gtest/gtest.h"

namespace gryphon {
namespace {
core::Input associationInput(double fraction = 1.) {
  core::Input in;
  in.set_syntheticAssociations(true);
  in.set_associationFraction(fraction);
  in.set_associationRadius(0.);
  in.set_associationVelocity(0.);
  in.set_sunRadius(0.);
  in.set_maxtime(1. * cgs::Myr);
  in.set_rate(0.001 / cgs::year);
  return in;
}
const auto fixedCentre = [](RandomNumberGenerator&) {
  return utils::Vector3d(2. * cgs::pc, 0., 0.);
};
}  // namespace

TEST(AssociationDelays, NormalizedMonotoneAndInvertible) {
  EXPECT_DOUBLE_EQ(galaxy::associationDelayCdf(0.), 0.);
  EXPECT_DOUBLE_EQ(galaxy::associationDelayCdf(100.), 1.);
  double previous = 0.;
  for (int i = 0; i <= 100; ++i) {
    const double p = i / 100.;
    const double t = galaxy::associationDelayQuantile(p);
    EXPECT_GE(t, previous);
    EXPECT_NEAR(galaxy::associationDelayCdf(t), p, 1e-12);
    previous = t;
  }
  EXPECT_GT(galaxy::associationDelayQuantile(0.5), 18.);
  EXPECT_LT(galaxy::associationDelayQuantile(0.5), 22.);
  EXPECT_THROW(galaxy::associationDelayQuantile(-0.1), std::invalid_argument);
}

TEST(AssociationInput, RejectsInvalidParameters) {
  auto in = associationInput();
  in.set_associationFraction(1.1);
  EXPECT_THROW(in.validate(), std::invalid_argument);
  in.set_associationFraction(0.5);
  in.set_associationMembers(0);
  EXPECT_THROW(in.validate(), std::invalid_argument);
  in.set_associationMembers(100);
  in.set_associationRadius(-1.);
  EXPECT_THROW(in.validate(), std::invalid_argument);
  in.set_associationRadius(0.);
  in.set_associationVelocity(std::numeric_limits<double>::infinity());
  EXPECT_THROW(in.validate(), std::invalid_argument);
}

TEST(AssociationGeneration, DeterministicAndFreshWithPhysicalAges) {
  auto in = associationInput();
  galaxy::AssociationGenerationStats firstStats, secondStats;
  RandomNumberGenerator firstRng(170), secondRng(170);
  const auto first = galaxy::generateAssociations(in, firstRng, fixedCentre, firstStats);
  const auto second = galaxy::generateAssociations(in, secondRng, fixedCentre, secondStats);
  ASSERT_GT(first.size(), 0u);
  ASSERT_EQ(first.size(), second.size());
  EXPECT_EQ(firstStats.fieldExplosions, 0u);
  EXPECT_GT(firstStats.associations, 0u);
  for (size_t i = 0; i < first.size(); ++i) {
    EXPECT_DOUBLE_EQ(first[i]->age, second[i]->age);
    EXPECT_GT(first[i]->age, 0.);
    EXPECT_LT(first[i]->age, in.max_time());
    EXPECT_DOUBLE_EQ(first[i]->pos.x, 2. * cgs::pc);
  }
}

TEST(AssociationGeneration, OldParentsPreserveRateAndUniformExplosionAges) {
  // T=1 Myr is shorter than *every* stellar delay. Omitting parent prehistory
  // would produce zero SNe. Ensemble count variance is bounded by N*R*T.
  for (const double fraction : {0., 0.5, 1.}) {
    auto in = associationInput(fraction);
    constexpr int samples = 100;
    double total = 0., meanAge = 0.;
    size_t halfCounts[2] = {0, 0};
    for (int seed = 0; seed < samples; ++seed) {
      RandomNumberGenerator rng(1000 + seed);
      galaxy::AssociationGenerationStats stats;
      const auto events = galaxy::generateAssociations(in, rng, fixedCentre, stats);
      total += stats.fieldExplosions + stats.clusteredExplosions;
      if (fraction == 0.) { EXPECT_EQ(stats.associations, 0u); }
      if (fraction == 1.) { EXPECT_EQ(stats.fieldExplosions, 0u); }
      for (const auto& event : events) {
        meanAge += event->age / in.max_time();
        ++halfCounts[event->age > 0.5 * in.max_time() ? 1 : 0];
      }
    }
    const double expected = in.sn_rate() * in.max_time() * samples;
    EXPECT_NEAR(total, expected, 0.05 * expected);
    EXPECT_NEAR(meanAge / (halfCounts[0] + halfCounts[1]), 0.5, 0.01);
    EXPECT_NEAR(double(halfCounts[0]) / (halfCounts[0] + halfCounts[1]), 0.5, 0.02);
  }
}

TEST(AssociationGeneration, ReusesCentresButDoesNotSynchronizeExplosions) {
  auto in = associationInput();
  in.set_maxtime(60. * cgs::Myr);
  size_t drawn = 0;
  auto sampler = [&drawn](RandomNumberGenerator&) {
    return utils::Vector3d((2. + ++drawn) * cgs::pc, 0., 0.);
  };
  RandomNumberGenerator rng(299);
  galaxy::AssociationGenerationStats stats;
  const auto events = galaxy::generateAssociations(in, rng, sampler, stats);
  ASSERT_GT(events.size(), 1u);
  size_t pairs = 0;
  for (size_t i = 1; i < events.size(); ++i) {
    if (events[i]->pos.x == events[i - 1]->pos.x) {
      ++pairs;
      EXPECT_NE(events[i]->age, events[i - 1]->age);
    }
  }
  EXPECT_GT(pairs, events.size() / 2);
}

TEST(AssociationGeneration, FactoryUsesOptInAndFiniteSpreading) {
  auto in = associationInput(0.5);
  in.set_spiralModel(SpiralModel::Xie2024);
  in.set_associationRadius(30. * cgs::pc);
  in.set_associationVelocity(3. * cgs::km / cgs::second);
  auto model = galaxy::makeGalaxy(in);
  RandomNumberGenerator rng(77);
  model->generate(rng, false);
  EXPECT_GT(model->association_stats().associations, 0u);
  EXPECT_EQ(model->size(), model->association_stats().retainedExplosions);
  for (const auto& event : model->get_events()) {
    EXPECT_LE(event->pos.getModule(), cgs::c_light * event->age);
    EXPECT_LE(event->pos.getModule(), in.R_g());
  }
}

TEST(AssociationDiagnostics, MetadataDoesNotChangeEventsOrRandomStream) {
  for (double fraction : {0., 0.5, 1.}) {
    auto in = associationInput(fraction);
    in.set_associationRadius(30. * cgs::pc);
    in.set_associationVelocity(3. * cgs::km / cgs::second);
    RandomNumberGenerator first(77), second(77);
    galaxy::AssociationGenerationStats firstStats, secondStats;
    std::vector<galaxy::AssociationEventOrigin> origins;
    const auto plain = galaxy::generateAssociations(in, first, fixedCentre, firstStats);
    const auto diagnostic = galaxy::generateAssociations(
        in, second, fixedCentre, secondStats, &origins);
    ASSERT_EQ(plain.size(), diagnostic.size());
    ASSERT_EQ(origins.size(), diagnostic.size());
    ASSERT_GT(origins.size(), 0u);
    EXPECT_EQ(firstStats.clusteredExplosions, secondStats.clusteredExplosions);
    EXPECT_EQ(firstStats.fieldExplosions, secondStats.fieldExplosions);
    EXPECT_DOUBLE_EQ(first(), second());
    std::map<size_t, double> parentBirths;
    for (size_t i = 0; i < plain.size(); ++i) {
      EXPECT_DOUBLE_EQ(plain[i]->age, diagnostic[i]->age);
      EXPECT_DOUBLE_EQ(plain[i]->pos.x, diagnostic[i]->pos.x);
      EXPECT_DOUBLE_EQ(plain[i]->pos.y, diagnostic[i]->pos.y);
      EXPECT_DOUBLE_EQ(plain[i]->pos.z, diagnostic[i]->pos.z);
      const auto& origin = origins[i];
      EXPECT_NEAR((origin.birthAge - origin.delay - plain[i]->age) / cgs::Myr, 0., 1e-12);
      EXPECT_GE(origin.delay / cgs::Myr, 3.);
      EXPECT_LE(origin.delay / cgs::Myr, 48.);
      EXPECT_LE(origin.birthAge / cgs::Myr, 49.);
      EXPECT_DOUBLE_EQ(origin.centre.x, 2. * cgs::pc);
      if (fraction == 0.) { EXPECT_EQ(origin.parentId, 0u); }
      if (fraction == 1.) { EXPECT_GT(origin.parentId, 0u); }
      if (origin.parentId > 0) {
        auto [parent, inserted] = parentBirths.emplace(origin.parentId, origin.birthAge);
        EXPECT_DOUBLE_EQ(parent->second, origin.birthAge);
      }
    }
    // Reusing a diagnostic buffer must not append stale origins.
    galaxy::generateAssociations(in, second, fixedCentre, secondStats, &origins);
    EXPECT_EQ(origins.size(), secondStats.retainedExplosions);
  }
}

TEST(AssociationDiagnostics, FactoryForwardsAndClearsMetadata) {
  auto in = associationInput();
  RandomNumberGenerator rng(4);
  std::vector<galaxy::AssociationEventOrigin> origins;
  auto model = galaxy::makeGalaxy(in);
  model->generate(rng, false, &origins);
  ASSERT_GT(origins.size(), 0u);
  EXPECT_EQ(origins.size(), model->size());
  in.set_syntheticAssociations(false);
  auto independent = galaxy::makeGalaxy(in);
  independent->generate(rng, false, &origins);
  EXPECT_TRUE(origins.empty());
}
}  // namespace gryphon
