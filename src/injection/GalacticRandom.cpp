#include "gryphon/injection/GalacticRandom.h"

#include <stdexcept>

#include "gryphon/utils/numeric.h"

namespace gryphon {

namespace {

double pickSnEnergy(RandomNumberGenerator& rng) {
  auto logEnergy = rng.normal(std::log10(3e50), 0.54);  // adimensional
  return std::max(std::pow(10., logEnergy) * cgs::erg, 0.);
}

double pickSlope(double mean, double sdev, RandomNumberGenerator& rng) {
  auto slope = rng.normal(mean, sdev);
  return std::max(0.01, slope);
}

}  // namespace

namespace injection {

GalacticRandomSpectrum::GalacticRandomSpectrum(const core::Input& in, RandomNumberGenerator& rng)
    : InjectionSpectrum(in) {
  m_alpha = (in.doVarySlope()) ? pickSlope(in.injSlope(), in.injSlopeSigma(), rng) : in.injSlope();
  m_crenergy = in.injEfficiency() * ((in.doVaryEnergy()) ? pickSnEnergy(rng) : cgs::E_SN);
  m_Emax = in.injEmax();
  m_Q0 = source_normalization();
}

double GalacticRandomSpectrum::source_normalization() const {
  const double prefactor = m_crenergy / pow2(m_E0);
  const double dalpha = m_alpha - 2.;
  const double tol = 1e-12;

  if (m_Emax <= 0.) {
    throw std::invalid_argument("GalacticRandomSpectrum requires injEmax > E0 for normalization");
  }

  const double ratio = m_Emax / m_E0;
  if (ratio <= 1.) {
    throw std::invalid_argument("GalacticRandomSpectrum requires injEmax > E0 for normalization");
  }

  if (std::abs(dalpha) <= tol) {
    return prefactor / std::log(ratio);
  }

  if (m_alpha < 2.) {
    const double exponent = 2. - m_alpha;
    return prefactor * exponent / (std::pow(ratio, exponent) - 1.);
  }

  return prefactor * dalpha / (1. - std::pow(ratio, 2. - m_alpha));
}

double GalacticRandomSpectrum::get(double E) const {
  auto value = m_Q0 * std::pow(E / m_E0, -m_alpha);
  if (m_Emax > 0.) value *= std::exp(-(E / m_Emax));
  return value;
}

}  // namespace injection
}  // namespace gryphon