#include <cmath>
#include <memory>

#include "gryphon.h"
#include "gtest/gtest.h"

namespace gryphon {
namespace {

core::Input singleSourceInput() {
  core::Input input;
  input.set_simEmin(10. * cgs::TeV);
  input.set_simEmax(100. * cgs::TeV);
  input.set_simEsize(3);
  input.set_injEmax(100. * cgs::PeV);
  input.set_transportModel(TransportModel::PureDiffusion);
  input.set_injectionModel(InjectionModel::SinglePowerLaw);
  return input;
}

injection::InjectionSpectra spectraFor(const core::Input& input, size_t count) {
  const auto spectrum = std::make_shared<injection::SinglePowerLawSpectrum>(input);
  return injection::InjectionSpectra(count, spectrum);
}

}  // namespace

TEST(Anisotropy, SingleSourceMatchesBurstSolution) {
  const auto input = singleSourceInput();
  const double age = 0.1 * cgs::Myr;
  const double distance = 0.3 * cgs::kpc;
  const core::Events events{
      std::make_shared<core::Event>(age, utils::Vector3d(distance, 0., 0.))};
  auto kernel = std::make_shared<kernel::PureDiffusionKernel>(input);

  core::CosmicRays cosmicRays(input, kernel, spectraFor(input, events.size()), events);
  cosmicRays.run();

  const double expected = 3. * distance / (2. * cgs::c_light * age);
  for (const auto& dipole : cosmicRays.get_dipole()) {
    EXPECT_NEAR(dipole.x, expected, 1e-12 * expected);
    EXPECT_NEAR(dipole.y, 0., 1e-15);
    EXPECT_NEAR(dipole.z, 0., 1e-15);
  }
}

TEST(Anisotropy, OppositeEqualSourcesCancel) {
  const auto input = singleSourceInput();
  const double age = 0.1 * cgs::Myr;
  const double distance = 0.3 * cgs::kpc;
  const core::Events events{
      std::make_shared<core::Event>(age, utils::Vector3d(distance, 0., 0.)),
      std::make_shared<core::Event>(age, utils::Vector3d(-distance, 0., 0.))};
  auto kernel = std::make_shared<kernel::PureDiffusionKernel>(input);

  core::CosmicRays cosmicRays(input, kernel, spectraFor(input, events.size()), events);
  cosmicRays.run();

  for (const auto& dipole : cosmicRays.get_dipole()) {
    EXPECT_NEAR(dipole.getModule(), 0., 1e-15);
  }
}

TEST(PureDiffusionKernel, SmoothBreakChangesHighEnergySlope) {
  core::Input input;
  input.set_refEnergy(1. * cgs::TeV);
  input.set_delta(0.36);
  input.set_ddelta(1.0);
  input.set_diffusionSmoothness(0.1);
  input.set_diffusionBreakEnergy(3.3 * cgs::PeV);
  const kernel::PureDiffusionKernel kernel(input);

  auto logarithmicSlope = [&kernel](double low, double high) {
    return std::log(kernel.D(high) / kernel.D(low)) / std::log(high / low);
  };

  EXPECT_NEAR(logarithmicSlope(1. * cgs::TeV, 10. * cgs::TeV), 0.36, 1e-8);
  EXPECT_NEAR(logarithmicSlope(33. * cgs::PeV, 330. * cgs::PeV), 1.36, 1e-8);
}

TEST(SmoothBrokenPowerLaw, HasConfiguredAsymptoticSlopes) {
  core::Input input;
  input.set_injectionModel(InjectionModel::SmoothBrokenPowerLaw);
  input.set_injSlope(2.15);
  input.set_injDeltaSlope(0.99);
  input.set_injSmoothness(0.2673);
  input.set_injBreakEnergy(3.3 * cgs::PeV);
  input.set_injEmax(1e3 * cgs::PeV);
  const injection::SmoothBrokenPowerLawSpectrum spectrum(input);

  auto indexBetween = [&spectrum](double low, double high) {
    return -std::log(spectrum.get(high) / spectrum.get(low)) / std::log(high / low);
  };

  EXPECT_NEAR(indexBetween(1. * cgs::TeV, 10. * cgs::TeV), 2.15, 1e-8);
  EXPECT_NEAR(indexBetween(33. * cgs::PeV, 330. * cgs::PeV), 3.14, 3e-5);
}

}  // namespace gryphon
