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

TEST(InjectionSpectrum, RandomEmaxUsesConfiguredValueWhenScatterDisabled) {
  core::Input in;
  in.set_injectionModel(InjectionModel::RandomEmax);
  in.set_injSlope(2.2);
  in.set_injEmax(300.0 * cgs::TeV);
  in.set_injEmaxSigmaDex(0.0);
  in.set_injEmaxMin(10.0 * cgs::TeV);
  in.set_injEmaxMax(2.0 * cgs::PeV);

  RandomNumberGenerator rng1(123);
  RandomNumberGenerator rng2(456);
  const injection::RandomEmaxSpectrum spectrum1(in, rng1);
  const injection::RandomEmaxSpectrum spectrum2(in, rng2);

  EXPECT_DOUBLE_EQ(spectrum1.alpha, in.injSlope());
  EXPECT_DOUBLE_EQ(spectrum1.crEnergy, in.injEfficiency() * cgs::E_SN);
  EXPECT_DOUBLE_EQ(spectrum1.Emax, in.injEmax());
  EXPECT_DOUBLE_EQ(spectrum1.get(10.0 * cgs::GeV), spectrum2.get(10.0 * cgs::GeV));
}

TEST(InjectionSpectrum, RandomEmaxIsReproducibleForSameSeed) {
  core::Input in;
  in.set_injectionModel(InjectionModel::RandomEmax);
  in.set_injSlope(2.1);
  in.set_injEmax(200.0 * cgs::TeV);
  in.set_injEmaxSigmaDex(0.4);
  in.set_injEmaxMin(10.0 * cgs::TeV);
  in.set_injEmaxMax(2.0 * cgs::PeV);

  RandomNumberGenerator rng1(2026);
  RandomNumberGenerator rng2(2026);
  const injection::RandomEmaxSpectrum spectrum1(in, rng1);
  const injection::RandomEmaxSpectrum spectrum2(in, rng2);

  EXPECT_DOUBLE_EQ(spectrum1.Emax, spectrum2.Emax);
  EXPECT_DOUBLE_EQ(spectrum1.Q0, spectrum2.Q0);
  EXPECT_DOUBLE_EQ(spectrum1.get(5.0 * cgs::TeV), spectrum2.get(5.0 * cgs::TeV));
}

TEST(InjectionSpectrum, RandomEmaxDrawStaysWithinConfiguredBounds) {
  core::Input in;
  in.set_injectionModel(InjectionModel::RandomEmax);
  in.set_injEmax(150.0 * cgs::TeV);
  in.set_injEmaxSigmaDex(0.6);
  in.set_injEmaxMin(20.0 * cgs::TeV);
  in.set_injEmaxMax(800.0 * cgs::TeV);

  for (unsigned long seed = 0; seed < 128; ++seed) {
    RandomNumberGenerator rng(seed);
    const injection::RandomEmaxSpectrum spectrum(in, rng);
    EXPECT_GE(spectrum.Emax, in.injEmaxMin());
    EXPECT_LE(spectrum.Emax, in.injEmaxMax());
    EXPECT_GT(spectrum.get(10.0 * cgs::GeV), 0.0);
  }
}

TEST(InjectionSpectrum, PWNUsesConfiguredValuesWhenVariationsDisabled) {
  core::Input in;
  in.set_injectionModel(InjectionModel::PWN);
  in.set_efficiency(0.12);
  in.set_pwnP0(85.0 * cgs::msec);
  in.set_pwnSigmaP0(0.0);
  in.set_pwnAlpha1(1.4);
  in.set_pwnAlpha2(2.7);
  in.set_pwnEmin(1.0 * cgs::GeV);
  in.set_pwnEbreak(100.0 * cgs::GeV);
  in.disable_varyenergy();
  in.disable_varyslope();

  RandomNumberGenerator rng1(111);
  RandomNumberGenerator rng2(222);
  const injection::PWNSpectrum spectrum1(in, rng1);
  const injection::PWNSpectrum spectrum2(in, rng2);

  EXPECT_DOUBLE_EQ(spectrum1.alphaBelowBreak, in.pwnAlpha1());
  EXPECT_DOUBLE_EQ(spectrum1.alphaAboveBreak, in.pwnAlpha2());
  EXPECT_DOUBLE_EQ(spectrum1.conversionEfficiency, in.injEfficiency());
  EXPECT_DOUBLE_EQ(spectrum1.initialPeriod, in.pwnP0());
  EXPECT_DOUBLE_EQ(spectrum1.Ebreak, in.pwnEbreak());
  EXPECT_DOUBLE_EQ(spectrum1.Emin, in.pwnEmin());
  const double B_s = std::pow(10., 12.65) * cgs::gauss;
  const double expected_emax = cgs::elementary_charge * 2. * M_PI * M_PI * B_s *
                               pow3(cgs::pulsar_radius) / cgs::c_2 / pow2(in.pwnP0());
  EXPECT_DOUBLE_EQ(spectrum1.Ecut, expected_emax);
  EXPECT_DOUBLE_EQ(spectrum1.Emax, expected_emax);
  EXPECT_DOUBLE_EQ(spectrum1.Ecut, spectrum1.Emax);
  EXPECT_GT(spectrum1.Emax, spectrum1.Ebreak);
  EXPECT_TRUE(std::isfinite(spectrum1.spinDownAge));
  EXPECT_GT(spectrum1.spinDownAge, 0.0);
  EXPECT_DOUBLE_EQ(spectrum1.rotEnergy, spectrum2.rotEnergy);
  EXPECT_DOUBLE_EQ(spectrum1.crEnergy, spectrum2.crEnergy);
  EXPECT_DOUBLE_EQ(spectrum1.get(10.0 * cgs::GeV), spectrum2.get(10.0 * cgs::GeV));
  EXPECT_DOUBLE_EQ(spectrum1.get(0.5 * cgs::GeV), 0.0);
  EXPECT_DOUBLE_EQ(spectrum1.get(2.0 * spectrum1.Emax), 0.0);
}

TEST(InjectionSpectrum, PWNIsReproducibleForSameSeed) {
  core::Input in;
  in.set_injectionModel(InjectionModel::PWN);
  in.set_efficiency(0.2);
  in.set_pwnP0(100.0 * cgs::msec);
  in.set_pwnSigmaP0(30.0 * cgs::msec);
  in.set_pwnAlpha1(1.5);
  in.set_pwnAlpha2(2.8);
  in.set_pwnEmin(1.0 * cgs::GeV);
  in.set_pwnEbreak(250.0 * cgs::GeV);
  in.enable_varyenergy();
  in.enable_varyslope();

  RandomNumberGenerator rng1(2025);
  RandomNumberGenerator rng2(2025);
  const injection::PWNSpectrum spectrum1(in, rng1);
  const injection::PWNSpectrum spectrum2(in, rng2);

  EXPECT_DOUBLE_EQ(spectrum1.alphaBelowBreak, in.pwnAlpha1());
  EXPECT_DOUBLE_EQ(spectrum1.alphaAboveBreak, in.pwnAlpha2());
  EXPECT_DOUBLE_EQ(spectrum1.alphaBelowBreak, spectrum2.alphaBelowBreak);
  EXPECT_DOUBLE_EQ(spectrum1.alphaAboveBreak, spectrum2.alphaAboveBreak);
  EXPECT_DOUBLE_EQ(spectrum1.initialPeriod, spectrum2.initialPeriod);
  EXPECT_DOUBLE_EQ(spectrum1.Ecut, spectrum2.Ecut);
  EXPECT_DOUBLE_EQ(spectrum1.Emax, spectrum2.Emax);
  EXPECT_DOUBLE_EQ(spectrum1.Ecut, spectrum1.Emax);
  EXPECT_DOUBLE_EQ(spectrum1.conversionEfficiency, spectrum2.conversionEfficiency);
  EXPECT_DOUBLE_EQ(spectrum1.rotEnergy, spectrum2.rotEnergy);
  EXPECT_DOUBLE_EQ(spectrum1.spinDownAge, spectrum2.spinDownAge);
  EXPECT_DOUBLE_EQ(spectrum1.crEnergy, spectrum2.crEnergy);
  EXPECT_DOUBLE_EQ(spectrum1.kickVelocity.getX(), spectrum2.kickVelocity.getX());
  EXPECT_DOUBLE_EQ(spectrum1.kickVelocity.getY(), spectrum2.kickVelocity.getY());
  EXPECT_DOUBLE_EQ(spectrum1.kickVelocity.getZ(), spectrum2.kickVelocity.getZ());
  EXPECT_DOUBLE_EQ(spectrum1.kickSpeed, spectrum2.kickSpeed);
  EXPECT_DOUBLE_EQ(spectrum1.Q0, spectrum2.Q0);
  EXPECT_DOUBLE_EQ(spectrum1.get(5.0 * cgs::TeV), spectrum2.get(5.0 * cgs::TeV));
}

TEST(InjectionSpectrum, PWNKeepsConversionEfficiencyFixedAcrossSeeds) {
  core::Input in;
  in.set_injectionModel(InjectionModel::PWN);
  in.set_efficiency(0.17);
  in.set_pwnP0(90.0 * cgs::msec);
  in.set_pwnSigmaP0(0.0);
  in.set_pwnAlpha1(1.5);
  in.set_pwnAlpha2(2.5);
  in.set_pwnEmin(1.0 * cgs::GeV);
  in.set_pwnEbreak(100.0 * cgs::GeV);
  in.enable_varyenergy();

  RandomNumberGenerator rng1(1);
  RandomNumberGenerator rng2(2);
  const injection::PWNSpectrum spectrum1(in, rng1);
  const injection::PWNSpectrum spectrum2(in, rng2);

  EXPECT_DOUBLE_EQ(spectrum1.conversionEfficiency, in.injEfficiency());
  EXPECT_DOUBLE_EQ(spectrum2.conversionEfficiency, in.injEfficiency());
  EXPECT_DOUBLE_EQ(spectrum1.conversionEfficiency, spectrum2.conversionEfficiency);
}

TEST(InjectionSpectrum, PWNComputesSpinDownAgeFromInitialPeriod) {
  core::Input in;
  in.set_injectionModel(InjectionModel::PWN);
  in.set_efficiency(0.1);
  in.set_pwnP0(120.0 * cgs::msec);
  in.set_pwnSigmaP0(0.0);
  in.set_pwnAlpha1(1.5);
  in.set_pwnAlpha2(2.5);
  in.set_pwnEmin(1.0 * cgs::GeV);
  in.set_pwnEbreak(100.0 * cgs::GeV);

  RandomNumberGenerator rng(7);
  const injection::PWNSpectrum spectrum(in, rng);

  const double B_s = std::pow(10., 12.65) * cgs::gauss;
  const double omega_0 = 2. * M_PI / in.pwnP0();
  const double moment_of_inertia = 0.4 * cgs::pulsar_mass * pow2(cgs::pulsar_radius);
  const double expected =
      3. * cgs::c_3 * moment_of_inertia / (pow2(B_s) * pow6(cgs::pulsar_radius) * pow2(omega_0));

  EXPECT_DOUBLE_EQ(spectrum.spinDownAge, expected);
}

TEST(InjectionSpectrum, PWNNormalizesSourceTermToEnergy) {
  core::Input in;
  in.set_injectionModel(InjectionModel::PWN);
  in.set_efficiency(0.1);
  in.set_pwnP0(120.0 * cgs::msec);
  in.set_pwnSigmaP0(0.0);
  in.set_pwnAlpha1(1.5);
  in.set_pwnAlpha2(2.5);
  in.set_pwnEmin(1.0 * cgs::GeV);
  in.set_pwnEbreak(100.0 * cgs::GeV);

  RandomNumberGenerator rng(19);
  const injection::PWNSpectrum spectrum(in, rng);

  auto integrate_energy = [&spectrum](double start, double stop) {
    if (!(start > 0.) || !(stop > start)) return 0.;

    auto integrand = [&spectrum](double logE) {
      const double E = std::exp(logE);
      return E * E * spectrum.get(E);
    };

    return utils::simpsonIntegration<double>(integrand, std::log(start), std::log(stop), 4096);
  };

  double total_energy = 0.;
  if (spectrum.Emin < spectrum.Ebreak) {
    total_energy += integrate_energy(spectrum.Emin, spectrum.Ebreak);
  }
  if (spectrum.Ebreak < spectrum.Ecut) {
    total_energy += integrate_energy(spectrum.Ebreak, spectrum.Ecut);
  }

  EXPECT_NEAR(total_energy, spectrum.crEnergy, spectrum.crEnergy * 1e-6);
}

TEST(InjectionSpectrum, PWNUsesBrokenPowerLawShape) {
  core::Input in;
  in.set_injectionModel(InjectionModel::PWN);
  in.set_efficiency(0.1);
  in.set_pwnP0(100.0 * cgs::msec);
  in.set_pwnSigmaP0(0.0);
  in.set_pwnAlpha1(1.0);
  in.set_pwnAlpha2(3.0);
  in.set_pwnEmin(1.0 * cgs::GeV);
  in.set_pwnEbreak(100.0 * cgs::GeV);
  in.disable_varyenergy();

  RandomNumberGenerator rng(17);
  const injection::PWNSpectrum spectrum(in, rng);

  const double below_ratio = spectrum.get(10.0 * cgs::GeV) / spectrum.get(100.0 * cgs::GeV);
  const double above_ratio = spectrum.get(100.0 * cgs::GeV) / spectrum.get(1000.0 * cgs::GeV);
  const double expected_below_ratio =
      10.0 * std::exp((100.0 * cgs::GeV - 10.0 * cgs::GeV) / spectrum.Ecut);
  const double expected_above_ratio =
      1000.0 * std::exp((1000.0 * cgs::GeV - 100.0 * cgs::GeV) / spectrum.Ecut);

  EXPECT_NEAR(below_ratio, expected_below_ratio, 1e-3);
  EXPECT_NEAR(above_ratio, expected_above_ratio, 1e-3);
}

TEST(InjectionSpectrum, PWNUsesFixedInitialPeriodWhenRandomDrawDisabled) {
  core::Input in;
  in.set_injectionModel(InjectionModel::PWN);
  in.set_pwnAlpha1(1.5);
  in.set_pwnAlpha2(2.5);
  in.set_pwnEmin(1.0 * cgs::GeV);
  in.set_pwnEbreak(100.0 * cgs::GeV);
  in.set_pwnP0(20.0 * cgs::msec);
  in.set_pwnSigmaP0(40.0 * cgs::msec);
  in.set_pwnRandomInitialPeriod(false);

  RandomNumberGenerator rng1(13);
  RandomNumberGenerator rng2(29);
  const injection::PWNSpectrum spectrum1(in, rng1);
  const injection::PWNSpectrum spectrum2(in, rng2);

  EXPECT_DOUBLE_EQ(spectrum1.initialPeriod, in.pwnP0());
  EXPECT_DOUBLE_EQ(spectrum2.initialPeriod, in.pwnP0());
}

TEST(InjectionSpectrum, PWNDrawsStrictlyPositiveInitialPeriod) {
  core::Input in;
  in.set_injectionModel(InjectionModel::PWN);
  in.set_pwnAlpha1(1.5);
  in.set_pwnAlpha2(2.5);
  in.set_pwnEmin(1.0 * cgs::GeV);
  in.set_pwnEbreak(100.0 * cgs::GeV);
  in.set_pwnP0(20.0 * cgs::msec);
  in.set_pwnSigmaP0(40.0 * cgs::msec);

  for (unsigned long seed = 0; seed < 128; ++seed) {
    RandomNumberGenerator rng(seed);
    const injection::PWNSpectrum spectrum(in, rng);
    EXPECT_GT(spectrum.initialPeriod, 0.0);
  }
}

TEST(InjectionSpectrum, PWNDrawsFiniteKickVelocityWithinExpectedRange) {
  core::Input in;
  in.set_injectionModel(InjectionModel::PWN);
  in.set_pwnAlpha1(1.5);
  in.set_pwnAlpha2(2.5);
  in.set_pwnEmin(1.0 * cgs::GeV);
  in.set_pwnEbreak(100.0 * cgs::GeV);

  RandomNumberGenerator rng(91);
  const injection::PWNSpectrum spectrum(in, rng);

  EXPECT_TRUE(std::isfinite(spectrum.kickVelocity.getX()));
  EXPECT_TRUE(std::isfinite(spectrum.kickVelocity.getY()));
  EXPECT_TRUE(std::isfinite(spectrum.kickVelocity.getZ()));
  EXPECT_TRUE(std::isfinite(spectrum.kickSpeed));
  EXPECT_GE(spectrum.kickVelocity.getX(), 0.0);
  EXPECT_GE(spectrum.kickVelocity.getY(), 0.0);
  EXPECT_GE(spectrum.kickVelocity.getZ(), 0.0);
  EXPECT_LE(spectrum.kickVelocity.getX(), 3000. * cgs::km / cgs::sec);
  EXPECT_LE(spectrum.kickVelocity.getY(), 3000. * cgs::km / cgs::sec);
  EXPECT_LE(spectrum.kickVelocity.getZ(), 3000. * cgs::km / cgs::sec);
}

TEST(InjectionSpectrum, YoungPulsarsUsesConfiguredValuesWhenRandomDrawsDisabled) {
  core::Input in;
  in.set_injectionModel(InjectionModel::YoungPulsars);
  in.set_pid(core::H);
  in.set_efficiency(0.13);
  in.set_youngPulsarsP0(60.0 * cgs::msec);
  in.set_youngPulsarsSigmaP0(10.0 * cgs::msec);
  in.set_youngPulsarsRandomInitialPeriod(false);
  in.set_youngPulsarsB0(2.5e12 * cgs::gauss);
  in.set_youngPulsarsSigmaLog10B(0.5);
  in.set_youngPulsarsRandomMagneticField(false);

  RandomNumberGenerator rng1(111);
  RandomNumberGenerator rng2(222);
  const injection::YoungPulsarsSpectrum spectrum1(in, rng1);
  const injection::YoungPulsarsSpectrum spectrum2(in, rng2);

  const double moment_of_inertia = 0.4 * cgs::pulsar_mass * pow2(cgs::pulsar_radius);
  const double omega_0 = 2. * M_PI / in.youngPulsarsP0();
  const double expected_rot_energy = 0.5 * moment_of_inertia * pow2(omega_0);
  const double expected_emax =
      cgs::elementary_charge * 2. * M_PI * M_PI * in.youngPulsarsB0() *
      pow3(cgs::pulsar_radius) / cgs::c_2 / pow2(in.youngPulsarsP0());
  const double expected_tau0 =
      3. * cgs::c_3 * moment_of_inertia /
      (pow2(in.youngPulsarsB0()) * pow6(cgs::pulsar_radius) * pow2(omega_0));

  EXPECT_DOUBLE_EQ(spectrum1.initialPeriod, in.youngPulsarsP0());
  EXPECT_DOUBLE_EQ(spectrum1.surfaceMagneticField, in.youngPulsarsB0());
  EXPECT_DOUBLE_EQ(spectrum1.chargeNumber, 1.0);
  EXPECT_DOUBLE_EQ(spectrum1.conversionEfficiency, in.injEfficiency());
  EXPECT_DOUBLE_EQ(spectrum1.rotEnergy, expected_rot_energy);
  EXPECT_DOUBLE_EQ(spectrum1.crEnergy, in.injEfficiency() * expected_rot_energy);
  EXPECT_DOUBLE_EQ(spectrum1.Emax, expected_emax);
  EXPECT_DOUBLE_EQ(spectrum1.Ecut, expected_emax);
  EXPECT_DOUBLE_EQ(spectrum1.tau0, expected_tau0);
  EXPECT_DOUBLE_EQ(spectrum1.Q0, spectrum1.crEnergy / spectrum1.Emax);
  EXPECT_DOUBLE_EQ(spectrum1.get(10.0 * cgs::GeV), spectrum2.get(10.0 * cgs::GeV));
  EXPECT_DOUBLE_EQ(spectrum1.get(0.0), 0.0);
}

TEST(InjectionSpectrum, YoungPulsarsIsReproducibleForSameSeed) {
  core::Input in;
  in.set_injectionModel(InjectionModel::YoungPulsars);

  RandomNumberGenerator rng1(2026);
  RandomNumberGenerator rng2(2026);
  const injection::YoungPulsarsSpectrum spectrum1(in, rng1);
  const injection::YoungPulsarsSpectrum spectrum2(in, rng2);

  EXPECT_DOUBLE_EQ(spectrum1.initialPeriod, spectrum2.initialPeriod);
  EXPECT_DOUBLE_EQ(spectrum1.surfaceMagneticField, spectrum2.surfaceMagneticField);
  EXPECT_DOUBLE_EQ(spectrum1.Emax, spectrum2.Emax);
  EXPECT_DOUBLE_EQ(spectrum1.rotEnergy, spectrum2.rotEnergy);
  EXPECT_DOUBLE_EQ(spectrum1.tau0, spectrum2.tau0);
  EXPECT_DOUBLE_EQ(spectrum1.crEnergy, spectrum2.crEnergy);
  EXPECT_DOUBLE_EQ(spectrum1.Q0, spectrum2.Q0);
  EXPECT_DOUBLE_EQ(spectrum1.get(5.0 * cgs::TeV), spectrum2.get(5.0 * cgs::TeV));
}

TEST(InjectionSpectrum, YoungPulsarsDrawsPhysicalPulsarParameters) {
  core::Input in;
  in.set_injectionModel(InjectionModel::YoungPulsars);
  in.set_youngPulsarsP0(60.0 * cgs::msec);
  in.set_youngPulsarsSigmaP0(10.0 * cgs::msec);
  in.set_youngPulsarsB0(2.5e12 * cgs::gauss);
  in.set_youngPulsarsSigmaLog10B(0.5);

  for (unsigned long seed = 0; seed < 128; ++seed) {
    RandomNumberGenerator rng(seed);
    const injection::YoungPulsarsSpectrum spectrum(in, rng);
    EXPECT_GT(spectrum.initialPeriod, 0.0);
    EXPECT_GT(spectrum.surfaceMagneticField, 0.0);
    EXPECT_GT(spectrum.Emax, 0.0);
    EXPECT_GT(spectrum.tau0, 0.0);
    EXPECT_GT(spectrum.rotEnergy, 0.0);
    EXPECT_GT(spectrum.get(1.0 * cgs::TeV), 0.0);
  }
}

TEST(InjectionSpectrum, YoungPulsarsNormalizesBurstSpectrumToRotationalEnergyFraction) {
  core::Input in;
  in.set_injectionModel(InjectionModel::YoungPulsars);
  in.set_efficiency(0.1);
  in.set_youngPulsarsP0(60.0 * cgs::msec);
  in.set_youngPulsarsRandomInitialPeriod(false);
  in.set_youngPulsarsB0(2.5e12 * cgs::gauss);
  in.set_youngPulsarsRandomMagneticField(false);

  RandomNumberGenerator rng(19);
  const injection::YoungPulsarsSpectrum spectrum(in, rng);

  const double start = 1e-6 * spectrum.Emax;
  const double stop = 100.0 * spectrum.Emax;
  auto integrand = [&spectrum](double logE) {
    const double E = std::exp(logE);
    return E * E * spectrum.get(E);
  };

  const double total_energy =
      utils::simpsonIntegration<double>(integrand, std::log(start), std::log(stop), 4096);
  const double expected_energy =
      spectrum.crEnergy * (std::exp(-start / spectrum.Emax) - std::exp(-stop / spectrum.Emax));

  EXPECT_NEAR(total_energy, expected_energy, expected_energy * 1e-6);
}

TEST(InjectionSpectrum, YoungPulsarsUsesInverseEnergyBurstShape) {
  core::Input in;
  in.set_injectionModel(InjectionModel::YoungPulsars);
  in.set_youngPulsarsP0(60.0 * cgs::msec);
  in.set_youngPulsarsRandomInitialPeriod(false);
  in.set_youngPulsarsB0(2.5e12 * cgs::gauss);
  in.set_youngPulsarsRandomMagneticField(false);

  RandomNumberGenerator rng(17);
  const injection::YoungPulsarsSpectrum spectrum(in, rng);

  const double low_energy = 0.1 * spectrum.Emax;
  const double high_energy = 0.5 * spectrum.Emax;
  const double ratio = spectrum.get(low_energy) / spectrum.get(high_energy);
  const double expected_ratio =
      (high_energy / low_energy) * std::exp((high_energy - low_energy) / spectrum.Emax);

  EXPECT_NEAR(ratio, expected_ratio, expected_ratio * 1e-12);
}

TEST(InjectionFactory, SinglePowerLawReusesSpectrumAcrossEvents) {
  core::Input in;
  in.set_injectionModel(InjectionModel::SinglePowerLaw);
  in.set_injEmax(1e3 * cgs::TeV);

  core::Events events;
  events.emplace_back(
      std::make_shared<core::Event>(1.0 * cgs::Myr, utils::Vector3d(1.0 * cgs::kpc, 0.0, 0.0)));
  events.emplace_back(
      std::make_shared<core::Event>(2.0 * cgs::Myr, utils::Vector3d(2.0 * cgs::kpc, 0.0, 0.0)));

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
  events.emplace_back(
      std::make_shared<core::Event>(1.0 * cgs::Myr, utils::Vector3d(1.0 * cgs::kpc, 0.0, 0.0)));
  events.emplace_back(
      std::make_shared<core::Event>(2.0 * cgs::Myr, utils::Vector3d(2.0 * cgs::kpc, 0.0, 0.0)));

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

TEST(InjectionFactory, RandomEmaxBuildsIndependentSpectraPerEventWhenScatterEnabled) {
  core::Input in;
  in.set_injectionModel(InjectionModel::RandomEmax);
  in.set_injEmax(150.0 * cgs::TeV);
  in.set_injEmaxSigmaDex(0.5);
  in.set_injEmaxMin(10.0 * cgs::TeV);
  in.set_injEmaxMax(2.0 * cgs::PeV);

  core::Events events;
  events.emplace_back(
      std::make_shared<core::Event>(1.0 * cgs::Myr, utils::Vector3d(1.0 * cgs::kpc, 0.0, 0.0)));
  events.emplace_back(
      std::make_shared<core::Event>(2.0 * cgs::Myr, utils::Vector3d(2.0 * cgs::kpc, 0.0, 0.0)));

  RandomNumberGenerator rng(6060);
  const auto spectra = injection::makeInjectionSpectra(in, events, rng);

  ASSERT_EQ(spectra.size(), events.size());
  ASSERT_NE(spectra[0], nullptr);
  ASSERT_NE(spectra[1], nullptr);
  EXPECT_NE(spectra[0], spectra[1]);

  const auto* first = dynamic_cast<const injection::RandomEmaxSpectrum*>(spectra[0].get());
  const auto* second = dynamic_cast<const injection::RandomEmaxSpectrum*>(spectra[1].get());
  ASSERT_NE(first, nullptr);
  ASSERT_NE(second, nullptr);
  EXPECT_NE(first->Emax, second->Emax);
}

TEST(InjectionFactory, RandomEmaxReusesSpectrumAcrossEventsWhenScatterDisabled) {
  core::Input in;
  in.set_injectionModel(InjectionModel::RandomEmax);
  in.set_injEmax(200.0 * cgs::TeV);
  in.set_injEmaxSigmaDex(0.0);
  in.set_injEmaxMin(10.0 * cgs::TeV);
  in.set_injEmaxMax(2.0 * cgs::PeV);

  core::Events events;
  events.emplace_back(
      std::make_shared<core::Event>(1.0 * cgs::Myr, utils::Vector3d(1.0 * cgs::kpc, 0.0, 0.0)));
  events.emplace_back(
      std::make_shared<core::Event>(2.0 * cgs::Myr, utils::Vector3d(2.0 * cgs::kpc, 0.0, 0.0)));

  RandomNumberGenerator rng(7070);
  const auto spectra = injection::makeInjectionSpectra(in, events, rng);

  ASSERT_EQ(spectra.size(), events.size());
  ASSERT_NE(spectra[0], nullptr);
  ASSERT_NE(spectra[1], nullptr);
  EXPECT_EQ(spectra[0], spectra[1]);
}

TEST(InjectionFactory, PWNBuildsIndependentSpectraPerEvent) {
  core::Input in;
  in.set_injectionModel(InjectionModel::PWN);
  in.set_pwnAlpha1(1.6);
  in.set_pwnAlpha2(2.7);
  in.set_pwnEmin(1.0 * cgs::GeV);
  in.set_pwnEbreak(150.0 * cgs::GeV);
  in.set_pwnP0(100.0 * cgs::msec);
  in.set_pwnSigmaP0(25.0 * cgs::msec);
  in.disable_varyenergy();
  in.disable_varyslope();

  core::Events events;
  events.emplace_back(
      std::make_shared<core::Event>(1.0 * cgs::Myr, utils::Vector3d(1.0 * cgs::kpc, 0.0, 0.0)));
  events.emplace_back(
      std::make_shared<core::Event>(2.0 * cgs::Myr, utils::Vector3d(2.0 * cgs::kpc, 0.0, 0.0)));

  RandomNumberGenerator rng(3030);
  const auto spectra = injection::makeInjectionSpectra(in, events, rng);

  ASSERT_EQ(spectra.size(), events.size());
  ASSERT_NE(spectra[0], nullptr);
  ASSERT_NE(spectra[1], nullptr);
  EXPECT_NE(spectra[0], spectra[1]);

  const auto* first = dynamic_cast<const injection::PWNSpectrum*>(spectra[0].get());
  const auto* second = dynamic_cast<const injection::PWNSpectrum*>(spectra[1].get());
  ASSERT_NE(first, nullptr);
  ASSERT_NE(second, nullptr);
  EXPECT_TRUE(first->crEnergy != second->crEnergy || first->initialPeriod != second->initialPeriod);
}

TEST(InjectionFactory, PWNBuildsIndependentSpectraEvenWhenSpectrumShapeIsFixed) {
  core::Input in;
  in.set_injectionModel(InjectionModel::PWN);
  in.set_pwnAlpha1(1.8);
  in.set_pwnAlpha2(2.9);
  in.set_pwnEmin(1.0 * cgs::GeV);
  in.set_pwnEbreak(200.0 * cgs::GeV);
  in.set_pwnP0(100.0 * cgs::msec);
  in.set_pwnSigmaP0(0.0);
  in.disable_varyenergy();
  in.enable_varyslope();

  core::Events events;
  events.emplace_back(
      std::make_shared<core::Event>(1.0 * cgs::Myr, utils::Vector3d(1.0 * cgs::kpc, 0.0, 0.0)));
  events.emplace_back(
      std::make_shared<core::Event>(2.0 * cgs::Myr, utils::Vector3d(2.0 * cgs::kpc, 0.0, 0.0)));

  RandomNumberGenerator rng(4040);
  const auto spectra = injection::makeInjectionSpectra(in, events, rng);

  ASSERT_EQ(spectra.size(), events.size());
  ASSERT_NE(spectra[0], nullptr);
  ASSERT_NE(spectra[1], nullptr);
  EXPECT_NE(spectra[0], spectra[1]);

  const auto* first = dynamic_cast<const injection::PWNSpectrum*>(spectra[0].get());
  const auto* second = dynamic_cast<const injection::PWNSpectrum*>(spectra[1].get());
  ASSERT_NE(first, nullptr);
  ASSERT_NE(second, nullptr);
  EXPECT_DOUBLE_EQ(first->alphaBelowBreak, in.pwnAlpha1());
  EXPECT_DOUBLE_EQ(first->alphaAboveBreak, in.pwnAlpha2());
  EXPECT_DOUBLE_EQ(first->conversionEfficiency, in.injEfficiency());
  EXPECT_DOUBLE_EQ(second->conversionEfficiency, in.injEfficiency());
  EXPECT_TRUE(first->kickVelocity.getX() != second->kickVelocity.getX() ||
              first->kickVelocity.getY() != second->kickVelocity.getY() ||
              first->kickVelocity.getZ() != second->kickVelocity.getZ());
}

TEST(InjectionFactory, YoungPulsarsBuildsIndependentSpectraPerEvent) {
  core::Input in;
  in.set_injectionModel(InjectionModel::YoungPulsars);

  core::Events events;
  events.emplace_back(
      std::make_shared<core::Event>(1.0 * cgs::Myr, utils::Vector3d(1.0 * cgs::kpc, 0.0, 0.0)));
  events.emplace_back(
      std::make_shared<core::Event>(2.0 * cgs::Myr, utils::Vector3d(2.0 * cgs::kpc, 0.0, 0.0)));

  RandomNumberGenerator rng(9090);
  const auto spectra = injection::makeInjectionSpectra(in, events, rng);

  ASSERT_EQ(spectra.size(), events.size());
  ASSERT_NE(spectra[0], nullptr);
  ASSERT_NE(spectra[1], nullptr);
  EXPECT_NE(spectra[0], spectra[1]);

  const auto* first = dynamic_cast<const injection::YoungPulsarsSpectrum*>(spectra[0].get());
  const auto* second = dynamic_cast<const injection::YoungPulsarsSpectrum*>(spectra[1].get());
  ASSERT_NE(first, nullptr);
  ASSERT_NE(second, nullptr);
  EXPECT_TRUE(first->initialPeriod != second->initialPeriod ||
              first->surfaceMagneticField != second->surfaceMagneticField);
}

}  // namespace gryphon
