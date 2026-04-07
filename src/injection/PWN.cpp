#include "gryphon/injection/PWN.h"

#include <gsl/gsl_sf_erf.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "gryphon/utils/numeric.h"

namespace gryphon {

namespace {

constexpr double kKickWeight1 = 0.90;
constexpr double kKickSigma1 = 160. * cgs::km / cgs::sec;
constexpr double kKickSigma2 = 780. * cgs::km / cgs::sec;
constexpr double kKickVelocityMax = 3000. * cgs::km / cgs::sec;
constexpr double kSurfaceMagneticField = std::pow(10., 12.65) * cgs::gauss;

double pickInitialPeriod(double mean_period, double sigma_period, RandomNumberGenerator& rng) {
  if (sigma_period < 0.) {
    throw std::invalid_argument("PWNSpectrum requires sigmaP0 >= 0");
  }
  if (mean_period <= 0.) {
    throw std::invalid_argument("PWNSpectrum requires mean P0 > 0");
  }
  if (sigma_period == 0.) return mean_period;

  while (true) {
    const double sampled_period = rng.normal(mean_period, sigma_period);
    if (sampled_period > 0.) return sampled_period;
  }
}

double pulsarMomentOfInertia() { return 0.4 * cgs::pulsar_mass * pow2(cgs::pulsar_radius); }

double rotationalEnergy(double period) {
  const double omega_0 = 2. * M_PI / period;
  return 0.5 * pulsarMomentOfInertia() * pow2(omega_0);
}

double spin_down_age(double period) {
  const double omega_0 = 2. * M_PI / period;
  return 3. * cgs::c_3 * pulsarMomentOfInertia() /
         (pow2(kSurfaceMagneticField) * pow6(cgs::pulsar_radius) * pow2(omega_0));
}

double maximum_potential_drop_energy(double period) {
  return cgs::elementary_charge * 2. * M_PI * M_PI * kSurfaceMagneticField *
         pow3(cgs::pulsar_radius) / cgs::c_2 / pow2(period);
}

double pickVelocityComponent(double w_1, double sigma_1, double sigma_2, double v_max, double rnd) {
  const double mixture_ratio = (1. - w_1) / w_1;
  const double normalization =
      gsl_sf_erf(v_max / M_SQRT2 / sigma_1) + mixture_ratio * gsl_sf_erf(v_max / M_SQRT2 / sigma_2);
  const double target = rnd * normalization;

  auto cdf_minus_target = [=](double v) {
    return gsl_sf_erf(v / M_SQRT2 / sigma_1) + mixture_ratio * gsl_sf_erf(v / M_SQRT2 / sigma_2) -
           target;
  };

  return utils::rootFinder<double>(cdf_minus_target, 0., v_max, 100, 1e-7);
}

utils::Vector3d drawPulsarKickVelocity(RandomNumberGenerator& rng) {
  const double v_x =
      pickVelocityComponent(kKickWeight1, kKickSigma1, kKickSigma2, kKickVelocityMax, rng());
  const double v_y =
      pickVelocityComponent(kKickWeight1, kKickSigma1, kKickSigma2, kKickVelocityMax, rng());
  const double v_z =
      pickVelocityComponent(kKickWeight1, kKickSigma1, kKickSigma2, kKickVelocityMax, rng());
  return utils::Vector3d(v_x, v_y, v_z);
}

}  // namespace

namespace injection {

PWNSpectrum::PWNSpectrum(const core::Input& in, RandomNumberGenerator& rng)
    : InjectionSpectrum(in) {
  m_alphaBelowBreak = in.pwnAlpha1();
  m_alphaAboveBreak = in.pwnAlpha2();
  m_Emin = in.pwnEmin();
  m_Ebreak = in.pwnEbreak();
  m_initialPeriod = pickInitialPeriod(in.pwnP0(), in.pwnSigmaP0(), rng);
  m_Ecut = maximum_potential_drop_energy(m_initialPeriod);
  m_conversionEfficiency = in.injEfficiency();
  m_rotEnergy = rotationalEnergy(m_initialPeriod);
  m_spinDownAge = spin_down_age(m_initialPeriod);
  m_crenergy = m_conversionEfficiency * m_rotEnergy;
  m_kickVelocity = drawPulsarKickVelocity(rng);
  m_kickSpeed = m_kickVelocity.getModule();
  m_Q0 = source_normalization();
}

double PWNSpectrum::source_normalization() const {
  if (!(m_Emin > 0.)) {
    throw std::invalid_argument("PWNSpectrum requires Emin > 0");
  }
  if (!(m_Ecut > 0.)) {
    throw std::invalid_argument("PWNSpectrum requires Ecut > 0");
  }
  if (!(m_Emin < m_Ebreak && m_Ebreak < m_Ecut)) {
    throw std::invalid_argument("PWNSpectrum requires Emin < Ebreak < Ecut");
  }

  auto integrate_energy_weighted_shape = [this](double start, double stop) {
    if (!(start > 0.) || !(stop > start)) return 0.;

    auto integrand = [this](double logE) {
      const double E = std::exp(logE);
      return E * E * spectralShape(E);
    };

    return utils::simpsonIntegration<double>(integrand, std::log(start), std::log(stop), 4096);
  };

  double integral = 0.;
  if (m_Emin < m_Ebreak) {
    integral += integrate_energy_weighted_shape(m_Emin, m_Ebreak);
  }
  if (m_Ebreak < m_Ecut) {
    integral += integrate_energy_weighted_shape(m_Ebreak, m_Ecut);
  }

  if (!(integral > 0.) || !std::isfinite(integral)) {
    throw std::invalid_argument("PWNSpectrum normalization integral must be positive and finite");
  }

  return m_crenergy / integral;
}

double PWNSpectrum::spectralShape(double E) const {
  if (E < m_Emin || E > m_Ecut) return 0.;

  const double alpha = (E < m_Ebreak) ? m_alphaBelowBreak : m_alphaAboveBreak;
  return std::pow(E / m_Ebreak, -alpha) * std::exp(-E / m_Ecut);
}

double PWNSpectrum::get(double E) const { return m_Q0 * spectralShape(E); }

}  // namespace injection
}  // namespace gryphon
