#include "gryphon/injection/YoungPulsars.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace gryphon {

namespace {

double drawPositiveGaussian(double mean, double sigma, const char* quantity,
                            RandomNumberGenerator& rng) {
  if (sigma < 0.) {
    throw std::invalid_argument(
        std::string("YoungPulsarsSpectrum requires non-negative sigma for ") + quantity);
  }
  if (mean <= 0.) {
    throw std::invalid_argument(std::string("YoungPulsarsSpectrum requires positive mean for ") +
                                quantity);
  }
  if (sigma == 0.) return mean;

  for (size_t draw = 0; draw < 1024; ++draw) {
    const double candidate = rng.normal(mean, sigma);
    if (candidate > 0.) return candidate;
  }

  throw std::runtime_error(std::string("Failed to draw positive ") + quantity +
                           " for YoungPulsarsSpectrum");
}

double drawLogNormalMagneticField(double mean_field, double sigma_log10,
                                  RandomNumberGenerator& rng) {
  if (mean_field <= 0.) {
    throw std::invalid_argument("YoungPulsarsSpectrum requires mean B_* > 0");
  }
  if (sigma_log10 < 0.) {
    throw std::invalid_argument("YoungPulsarsSpectrum requires sigmaLog10B >= 0");
  }
  if (sigma_log10 == 0.) return mean_field;

  const double log10_mean = std::log10(mean_field / cgs::gauss);
  const double sampled_log10 = rng.normal(log10_mean, sigma_log10);
  return std::pow(10., sampled_log10) * cgs::gauss;
}

double pulsarMomentOfInertia() { return 0.4 * cgs::pulsar_mass * pow2(cgs::pulsar_radius); }

double rotationalEnergy(double period) {
  const double omega_0 = 2. * M_PI / period;
  return 0.5 * pulsarMomentOfInertia() * pow2(omega_0);
}

double spinDownAge(double period, double surface_magnetic_field) {
  const double omega_0 = 2. * M_PI / period;
  return 3. * cgs::c_3 * pulsarMomentOfInertia() /
         (pow2(surface_magnetic_field) * pow6(cgs::pulsar_radius) * pow2(omega_0));
}

}  // namespace

namespace injection {

YoungPulsarsSpectrum::YoungPulsarsSpectrum(const core::Input& in, RandomNumberGenerator& rng)
    : InjectionSpectrum(in) {
  m_initialPeriod =
      in.youngPulsarsRandomInitialPeriod()
          ? drawPositiveGaussian(in.youngPulsarsP0(), in.youngPulsarsSigmaP0(), "P0", rng)
          : in.youngPulsarsP0();
  m_surfaceMagneticField =
      in.youngPulsarsRandomMagneticField()
          ? drawLogNormalMagneticField(in.youngPulsarsB0(), in.youngPulsarsSigmaLog10B(), rng)
          : in.youngPulsarsB0();
  m_chargeNumber = std::max(1, in.pid().get_Z());
  m_conversionEfficiency = in.injEfficiency();
  m_rotEnergy = rotationalEnergy(m_initialPeriod);
  m_crenergy = m_conversionEfficiency * m_rotEnergy;
  m_Emax = maximumPotentialDropEnergy();
  m_tau0 = spinDownAge(m_initialPeriod, m_surfaceMagneticField);
  m_Q0 = m_crenergy / m_Emax;
}

double YoungPulsarsSpectrum::maximumPotentialDropEnergy() const {
  return m_chargeNumber * cgs::elementary_charge * 2. * M_PI * M_PI * m_surfaceMagneticField *
         pow3(cgs::pulsar_radius) / cgs::c_2 / pow2(m_initialPeriod);
}

double YoungPulsarsSpectrum::get(double E) const {
  if (!(E > 0.)) return 0.;
  return m_Q0 / E * std::exp(-E / m_Emax);
}

}  // namespace injection
}  // namespace gryphon
