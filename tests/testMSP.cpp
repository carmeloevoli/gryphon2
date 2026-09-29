#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numeric>
#include <unistd.h>

#include "gryphon.h"
#include "gryphon/injection/MSP.h"
#include "gryphon/kernel/ContinuousLosses.h"
#include "gtest/gtest.h"

namespace gryphon {
namespace {
core::Input benchmark() {
  core::Input in;
  in.set_injectionModel(InjectionModel::MSP);
  in.set_transportModel(TransportModel::DiffusionLosses);
  in.set_spiralModel(SpiralModel::Jelly);
  in.set_discsize(cgs::kpc);
  in.set_halosize(5.*cgs::kpc);
  in.set_simEmin(10.*cgs::GeV);
  in.set_simEmax(1000.*cgs::GeV);
  in.set_simEsize(3);
  in.set_refEnergy(10.*cgs::GeV);
  in.set_D0_over_H(.04*cgs::kpc/cgs::Myr);
  in.set_delta(.56);
  in.set_Urad(.475*cgs::eV/cgs::cm3);
  in.set_injSlope(1.);
  in.set_injEmax(400.*cgs::GeV);
  in.set_efficiency(1.);
  return in;
}
}  // namespace

TEST(MSP, EnergyNormalizationAndChargeSplit) {
  auto in = benchmark();
  for (double alpha : {1., 2., 2.7}) {
    in.set_injSlope(alpha);
    injection::MSPInjection source(in);
    const double luminosity = 1e34*cgs::erg/cgs::second;
    const double power = utils::QAGIntegration<double>([&](double loge) {
      const double e = std::exp(loge)*cgs::GeV;
      return e*e*source.rate(e, luminosity);
    }, std::log(10.), std::log(400.*80.), 1000, 1e-9);
    EXPECT_NEAR(power/luminosity, .5, 1e-9);
    EXPECT_EQ(source.rate(9.*cgs::GeV, luminosity), 0.);
  }
}

TEST(MSP, SpinSamplingAndReproducibility) {
  auto in = benchmark();
  injection::MSPInjection spectrum(in);
  RandomNumberGenerator rng(432), same(432);
  std::vector<double> periods;
  double logB = 0., logB2 = 0., power = 0.;
  constexpr size_t count = 100000;
  for (size_t i = 0; i < count; ++i) {
    const auto s = spectrum.sample(rng);
    const auto t = spectrum.sample(same);
    EXPECT_DOUBLE_EQ(s.period, t.period);
    EXPECT_DOUBLE_EQ(s.luminosity, t.luminosity);
    EXPECT_GE(s.period, in.mspPeriodMin());
    EXPECT_NEAR(3.2e19*std::sqrt(s.period*s.periodDot)/s.field, 1., 1e-12);
    periods.push_back(s.period);
    logB += std::log10(s.field);
    logB2 += pow2(std::log10(s.field));
    power += s.luminosity;
  }
  std::sort(periods.begin(), periods.end());
  EXPECT_NEAR(periods[count/2]/cgs::msec, 3., .05);
  EXPECT_NEAR(logB/count, 8., .005);
  EXPECT_NEAR(logB2/count-pow2(logB/count), .04, .002);
  const double b2 = std::pow(10., 2.*in.mspLog10B()) * std::exp(2.*pow2(std::log(10.)*in.mspSigmaLog10B()));
  const double expected = 4.*M_PI*M_PI*in.mspInertia()/pow2(3.2e19)*b2/(5.*pow4(in.mspPeriodMin()));
  EXPECT_NEAR(power/count/expected, 1., .04);
}

TEST(MSP, JellySamplingUsesNativeProfileAndObserverFrame) {
  const auto in = benchmark();
  auto galaxy = galaxy::makeGalaxy(in);
  RandomNumberGenerator rng(9876);
  double xmean = 0., z2 = 0., rmean = 0.;
  constexpr size_t count = 20000;
  for (size_t i = 0; i < count; ++i) {
    const auto p = galaxy->sample_position(rng);
    const double r = std::hypot(p.x+in.R_sun(), p.y);
    EXPECT_LE(r, in.R_g()*(1.+1e-12));
    EXPECT_GE(p.getModule(), cgs::pc);
    xmean += p.x;
    z2 += pow2(p.z);
    rmean += r;
  }
  EXPECT_NEAR(xmean/count/cgs::kpc, -in.R_sun()/cgs::kpc, .12);
  EXPECT_NEAR(z2/count/pow2(in.h()), 1., .04);
  core::SourceProfile profile(in);
  const auto integral = [&](int power) {
    return utils::QAGIntegration<double>([&](double r) {
      return std::pow(r, power)*profile.get(r*cgs::kpc);
    }, 0., in.R_g()/cgs::kpc, 1000, 1e-8);
  };
  EXPECT_NEAR(rmean/count/cgs::kpc, integral(2)/integral(1), .08);
}

TEST(MSP, CachedResponseMatchesExistingHaloFunction) {
  const auto in = benchmark();
  kernel::ContinuousLosses solver(in);
  for (size_t e = 0; e < solver.energies().size(); ++e) {
    for (auto p : {utils::Vector3d(.001*cgs::kpc, 0., 0.),
                   utils::Vector3d(.3*cgs::kpc, .4*cgs::kpc, .1*cgs::kpc),
                   utils::Vector3d(4.*cgs::kpc, 0., cgs::kpc),
                   utils::Vector3d(.2*cgs::kpc, 0., 4.9*cgs::kpc)}) {
      const double direct = solver.directResponse(solver.energies()[e], p);
      EXPECT_NEAR(solver.response(e, p)/direct, 1., .003);
    }
    EXPECT_EQ(solver.response(e, utils::Vector3d(cgs::kpc, 0., 5.*cgs::kpc)), 0.);
    EXPECT_EQ(solver.response(e, utils::Vector3d(cgs::kpc, 0., 6.*cgs::kpc)), 0.);
  }
}

TEST(MSP, StationaryIntegralMatchesTimeIntegralOfBurstKernel) {
  const auto in = benchmark();
  injection::MSPInjection rate(in);
  kernel::ContinuousLosses solver(in);
  kernel::DiffusionLossesKernel burst(in);
  const utils::Vector3d pos(.3*cgs::kpc, 0., .1*cgs::kpc);
  const double e = 100.*cgs::GeV, loss = burst.energyLossTimescale(e);
  const double integrated = utils::QAGIntegration<double>([&](double logt) {
    const double t = std::exp(logt)*cgs::Myr;
    return t*burst.flux(e, t, pos, [&](double es) { return rate.rate(es, 1.); });
  }, std::log(1e-8), std::log(loss/cgs::Myr), 1000, 1e-7);
  EXPECT_NEAR(integrated/solver.directResponse(e, pos), 1., 2e-6);
}

TEST(MSP, EfficiencyScalingAndZero) {
  auto in = benchmark();
  kernel::ContinuousLosses full(in);
  in.set_efficiency(.1);
  kernel::ContinuousLosses tenth(in);
  in.set_efficiency(0.);
  kernel::ContinuousLosses zero(in);
  const utils::Vector3d pos(.3*cgs::kpc, 0., .1*cgs::kpc);
  for (size_t e = 0; e < full.energies().size(); ++e) {
    EXPECT_NEAR(tenth.response(e, pos)/full.response(e, pos), .1, 1e-12);
    EXPECT_EQ(zero.response(e, pos), 0.);
  }
}

TEST(MSP, ConfigurationRoundTripAndBurstSafety) {
  auto in = benchmark();
  in.set_mspSources(12345);
  in.set_mspPeriodMin(2.*cgs::msec);
  in.set_mspLog10B(8.2);
  in.set_mspSigmaLog10B(.3);
  in.set_mspInertia(1.2e45);
  in.set_mspEmin(5.*cgs::GeV);
  in.set_mspQuadrature(512);
  in.set_mspDistances(3201);
  char path[] = "/tmp/gryphon_msp_XXXXXX";
  const int fd = mkstemp(path);
  ASSERT_GE(fd, 0);
  close(fd);
  in.write_params_file(path);
  const core::Input recovered(path);
  EXPECT_EQ(recovered.mspSources(), 12345ul);
  EXPECT_DOUBLE_EQ(recovered.mspPeriodMin(), in.mspPeriodMin());
  EXPECT_DOUBLE_EQ(recovered.mspLog10B(), in.mspLog10B());
  EXPECT_DOUBLE_EQ(recovered.mspSigmaLog10B(), in.mspSigmaLog10B());
  EXPECT_NEAR(recovered.mspInertia()/in.mspInertia(), 1., 1e-12);
  EXPECT_DOUBLE_EQ(recovered.mspEmin(), in.mspEmin());
  EXPECT_EQ(recovered.mspQuadrature(), in.mspQuadrature());
  EXPECT_EQ(recovered.mspDistances(), in.mspDistances());
  EXPECT_EQ(std::remove(path), 0);
  RandomNumberGenerator rng(1);
  EXPECT_THROW(injection::makeInjectionSpectrum(in, rng), std::invalid_argument);
  EXPECT_THROW(injection::makeInjectionSpectra(in, {}, rng), std::invalid_argument);
}

TEST(MSP, RejectsInvalidOrUnsupportedConfigurations) {
  auto good = benchmark();
  EXPECT_NO_THROW(good.validate());
  for (double value : {-1., std::numeric_limits<double>::infinity(), std::nan("")}) {
    auto in = good;
    in.set_mspInertia(value);
    EXPECT_THROW(in.validate(), std::invalid_argument);
    in = good;
    in.set_efficiency(value);
    EXPECT_THROW(in.validate(), std::invalid_argument);
  }
  auto in = good;
  in.set_mspSources(0);
  EXPECT_THROW(in.validate(), std::invalid_argument);
  in = good;
  in.set_transportModel(TransportModel::PureDiffusion);
  EXPECT_THROW(in.validate(), std::invalid_argument);
  in = good;
  in.set_ddelta(.2);
  EXPECT_THROW(in.validate(), std::invalid_argument);
}
}  // namespace gryphon
