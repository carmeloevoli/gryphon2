#include "gryphon/injection/SinglePowerLaw.h"

#include <stdexcept>

#include "gryphon/utils/numeric.h"

namespace gryphon {

namespace injection {

SinglePowerLawSpectrum::SinglePowerLawSpectrum(const core::Input& in)
    : InjectionSpectrum(in), m_E0(in.E_0()) {
  m_alpha = in.injSlope();
  m_crenergy = in.injEfficiency() * cgs::E_SN;
  m_Emax = in.injEmax();
  m_Q0 = source_normalization();
}

double SinglePowerLawSpectrum::source_normalization() const {
  const double prefactor = m_crenergy / pow2(m_E0);
  const double dalpha = m_alpha - 2.;
  const double tol = 1e-12;

  if (m_Emax <= 0.) {
    throw std::invalid_argument("SinglePowerLawSpectrum requires injEmax > E0 for normalization");
  }

  const double ratio = m_Emax / m_E0;
  if (ratio <= 1.) {
    throw std::invalid_argument("SinglePowerLawSpectrum requires injEmax > E0 for normalization");
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

double SinglePowerLawSpectrum::get(double E) const {
  auto value = m_Q0 * std::pow(E / m_E0, -m_alpha);
  if (m_Emax > 0.) value *= std::exp(-(E / m_Emax));
  return value;
}

}  // namespace injection
}  // namespace gryphon