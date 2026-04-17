#include <algorithm>
#include <cmath>
#include <vector>

#include "gryphon.h"
#include "gtest/gtest.h"

namespace gryphon {
namespace {

core::Input makeGalaxyTestInput(SpiralModel model) {
  core::Input in;
  in.set_spiralModel(model);
  in.set_maxtime(0.05 * cgs::Myr);
  in.set_rate(1.0 / 50.0 / cgs::year);
  return in;
}

void expectPhysicalEvents(const core::Events& events, double maxTime) {
  ASSERT_FALSE(events.empty());
  for (const auto& event : events) {
    ASSERT_NE(event, nullptr);
    EXPECT_GT(event->age, 0.0);
    EXPECT_LE(event->age, maxTime);
    EXPECT_TRUE(std::isfinite(event->pos.getX()));
    EXPECT_TRUE(std::isfinite(event->pos.getY()));
    EXPECT_TRUE(std::isfinite(event->pos.getZ()));
    EXPECT_LE(event->pos.getModule(), cgs::c_light * event->age * (1.0 + 1e-12));
  }
}

}  // namespace

TEST(GalaxyFactory, CreatesSupportedModels) {
  for (const auto model :
       {SpiralModel::Uniform, SpiralModel::Jelly, SpiralModel::Steiman2010,
        SpiralModel::Xie2024}) {
    core::Input in = makeGalaxyTestInput(model);
    auto galaxyModel = galaxy::makeGalaxy(in);
    EXPECT_NE(galaxyModel, nullptr);
  }
}

TEST(GalaxyFactory, RejectsUnsupportedModels) {
  for (const auto model : {SpiralModel::Faucher2006, SpiralModel::Vallee2008}) {
    core::Input in = makeGalaxyTestInput(model);
    EXPECT_THROW(galaxy::makeGalaxy(in), std::invalid_argument);
  }
}

TEST(GalaxyGeneration, SupportedModelsProducePhysicalEvents) {
  for (const auto model :
       {SpiralModel::Uniform, SpiralModel::Jelly, SpiralModel::Steiman2010,
        SpiralModel::Xie2024}) {
    core::Input in = makeGalaxyTestInput(model);
    RandomNumberGenerator rng(2026);
    auto galaxyModel = galaxy::makeGalaxy(in);
    galaxyModel->generate(rng, false);

    const auto& events = galaxyModel->get_events();
    expectPhysicalEvents(events, in.max_time());

    if (model == SpiralModel::Uniform) {
      for (const auto& event : events) {
        EXPECT_DOUBLE_EQ(event->pos.getZ(), 0.0);
      }
    }
  }
}

TEST(GalaxyGeneration, RegenerateWithSameSeedProducesSameSequence) {
  core::Input in = makeGalaxyTestInput(SpiralModel::Jelly);
  auto galaxyModel = galaxy::makeGalaxy(in);

  RandomNumberGenerator rng1(77);
  galaxyModel->generate(rng1, false);
  const auto firstRunSize = galaxyModel->size();
  ASSERT_GT(firstRunSize, 0u);

  std::vector<double> ages;
  std::vector<utils::Vector3d> positions;
  const auto& firstEvents = galaxyModel->get_events();
  const size_t sampleCount = std::min<size_t>(10, firstEvents.size());
  ages.reserve(sampleCount);
  positions.reserve(sampleCount);
  for (size_t i = 0; i < sampleCount; ++i) {
    ages.push_back(firstEvents[i]->age);
    positions.push_back(firstEvents[i]->pos);
  }

  RandomNumberGenerator rng2(77);
  galaxyModel->generate(rng2, false);
  const auto& secondEvents = galaxyModel->get_events();

  ASSERT_EQ(secondEvents.size(), firstRunSize);
  for (size_t i = 0; i < sampleCount; ++i) {
    EXPECT_DOUBLE_EQ(secondEvents[i]->age, ages[i]);
    EXPECT_DOUBLE_EQ(secondEvents[i]->pos.getX(), positions[i].getX());
    EXPECT_DOUBLE_EQ(secondEvents[i]->pos.getY(), positions[i].getY());
    EXPECT_DOUBLE_EQ(secondEvents[i]->pos.getZ(), positions[i].getZ());
  }
}

}  // namespace gryphon
