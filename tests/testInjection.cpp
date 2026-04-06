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
  in.set_injEmax(0.5 * cgs::GeV);

  EXPECT_THROW({ const injection::SinglePowerLawSpectrum spectrum(in); }, std::invalid_argument);
}

TEST(InjectionSpectrum, SinglePowerLawNormalizesFromFixedGeVThreshold) {
  core::Input in;
  in.set_refEnergy(1.0 * cgs::TeV);
  in.set_injSlope(2.3);
  in.set_injEmax(10.0 * cgs::GeV);

  const injection::SinglePowerLawSpectrum spectrum(in);
  const double value = spectrum.get(5.0 * cgs::GeV);

  EXPECT_TRUE(std::isfinite(value));
  EXPECT_GT(value, 0.0);
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

TEST(InjectionFactory, SinglePowerLawReusesSpectrumAcrossEvents) {
  core::Input in;
  in.set_injectionModel(InjectionModel::SinglePowerLaw);
  in.set_injEmax(1e3 * cgs::TeV);

  core::Events events;
  events.emplace_back(std::make_shared<core::Event>(1.0 * cgs::Myr,
                                                    utils::Vector3d(1.0 * cgs::kpc, 0.0, 0.0)));
  events.emplace_back(std::make_shared<core::Event>(2.0 * cgs::Myr,
                                                    utils::Vector3d(2.0 * cgs::kpc, 0.0, 0.0)));

  RandomNumberGenerator rng(11);
  const auto spectra = injection::makeInjectionSpectra(in, events, rng);

  ASSERT_EQ(spectra.size(), events.size());
  ASSERT_NE(spectra[0], nullptr);
  ASSERT_NE(spectra[1], nullptr);
  EXPECT_EQ(spectra[0], spectra[1]);
}

TEST(InjectionFactory, GalacticRandomBuildsIndependentSpectraPerEvent) {
  core::Input in;
  in.set_injectionModel(InjectionModel::GalacticRandom);
  in.set_injEmax(1e3 * cgs::TeV);
  in.enable_varyenergy();
  in.enable_varyslope();

  core::Events events;
  events.emplace_back(std::make_shared<core::Event>(1.0 * cgs::Myr,
                                                    utils::Vector3d(1.0 * cgs::kpc, 0.0, 0.0)));
  events.emplace_back(std::make_shared<core::Event>(2.0 * cgs::Myr,
                                                    utils::Vector3d(2.0 * cgs::kpc, 0.0, 0.0)));

  RandomNumberGenerator rng(2024);
  const auto spectra = injection::makeInjectionSpectra(in, events, rng);

  ASSERT_EQ(spectra.size(), events.size());
  ASSERT_NE(spectra[0], nullptr);
  ASSERT_NE(spectra[1], nullptr);
  EXPECT_NE(spectra[0], spectra[1]);

  const auto* first = dynamic_cast<const injection::GalacticRandomSpectrum*>(spectra[0].get());
  const auto* second = dynamic_cast<const injection::GalacticRandomSpectrum*>(spectra[1].get());
  ASSERT_NE(first, nullptr);
  ASSERT_NE(second, nullptr);
  EXPECT_TRUE(first->alpha != second->alpha || first->crEnergy != second->crEnergy);
}

}  // namespace gryphon
