#include <cmath>

#include "gryphon.h"
#include "gtest/gtest.h"

namespace gryphon {

TEST(InjectionSpectrum, SinglePowerLawIsPositiveAndDecreasesWithEnergy) {
  core::Input in;
  in.set_injSlope(2.3);
  in.set_injEmax(1e3 * cgs::TeV);
  in.set_efficiency(0.2);

  const injection::SinglePowerLawSpectrum spectrum(in);
  const double lowEnergy = spectrum.get(10.0 * cgs::GeV);
  const double highEnergy = spectrum.get(100.0 * cgs::GeV);

  EXPECT_TRUE(std::isfinite(lowEnergy));
  EXPECT_TRUE(std::isfinite(highEnergy));
  EXPECT_GT(lowEnergy, 0.0);
  EXPECT_GT(highEnergy, 0.0);
  EXPECT_GT(lowEnergy, highEnergy);
}

TEST(InjectionSpectrum, SinglePowerLawRejectsInvalidCutoff) {
  core::Input in;
  in.set_refEnergy(1.0 * cgs::TeV);
  in.set_injEmax(0.5 * cgs::TeV);

  EXPECT_THROW({ const injection::SinglePowerLawSpectrum spectrum(in); }, std::invalid_argument);
}

TEST(InjectionSpectrum, GalacticRandomUsesConfiguredValuesWhenVariationsDisabled) {
  core::Input in;
  in.set_injSlope(2.25);
  in.set_efficiency(0.15);
  in.disable_varyenergy();
  in.disable_varyslope();

  RandomNumberGenerator rng1(1234);
  RandomNumberGenerator rng2(5678);
  const injection::GalacticRandomSpectrum spectrum1(in, rng1);
  const injection::GalacticRandomSpectrum spectrum2(in, rng2);

  EXPECT_DOUBLE_EQ(spectrum1.alpha, in.injSlope());
  EXPECT_DOUBLE_EQ(spectrum1.crEnergy, in.injEfficiency() * cgs::E_SN);
  EXPECT_DOUBLE_EQ(spectrum1.get(10.0 * cgs::GeV), spectrum2.get(10.0 * cgs::GeV));
}

TEST(InjectionSpectrum, GalacticRandomIsReproducibleForSameSeed) {
  core::Input in;
  in.set_efficiency(0.2);
  in.set_injEmax(1e3 * cgs::TeV);
  in.enable_varyenergy();
  in.enable_varyslope();

  RandomNumberGenerator rng1(2024);
  RandomNumberGenerator rng2(2024);
  const injection::GalacticRandomSpectrum spectrum1(in, rng1);
  const injection::GalacticRandomSpectrum spectrum2(in, rng2);

  EXPECT_DOUBLE_EQ(spectrum1.alpha, spectrum2.alpha);
  EXPECT_DOUBLE_EQ(spectrum1.crEnergy, spectrum2.crEnergy);
  EXPECT_DOUBLE_EQ(spectrum1.Q0, spectrum2.Q0);
  EXPECT_DOUBLE_EQ(spectrum1.get(5.0 * cgs::TeV), spectrum2.get(5.0 * cgs::TeV));
}

}  // namespace gryphon
