#include "gryphon/injection/SinglePowerLaw.h"

#include <stdexcept>

#include "gryphon/utils/numeric.h"

namespace gryphon {

namespace injection {

SinglePowerLawSpectrum::SinglePowerLawSpectrum(const core::Input& in)
    : InjectionSpectrum(in), m_E0(in.E_0()), m_Emin(1. * cgs::GeV) {
  m_alpha = in.injSlope();
  m_crenergy = in.injEfficiency() * cgs::E_SN;
  m_Emax = in.injEmax();
  m_Q0 = source_normalization();
}

double SinglePowerLawSpectrum::source_normalization() const {
  const double prefactor = m_crenergy / pow2(m_E0);
  const double dalpha = m_alpha - 2.;
  const double tol = 1e-12;

  if (m_Emax <= m_Emin) {
    throw std::invalid_argument(
        "SinglePowerLawSpectrum requires injEmax > 1 GeV for normalization");
  }

  const double min_ratio = m_Emin / m_E0;
  const double max_ratio = m_Emax / m_E0;

  if (std::abs(dalpha) <= tol) {
    return prefactor / std::log(max_ratio / min_ratio);
  }

  const double exponent = 2. - m_alpha;
  return prefactor * exponent / (std::pow(max_ratio, exponent) - std::pow(min_ratio, exponent));
}

double SinglePowerLawSpectrum::get(double E) const {
  auto value = m_Q0 * std::pow(E / m_E0, -m_alpha);
  if (m_Emax > 0.) value *= std::exp(-(E / m_Emax));
  return value;
}

}  // namespace injection
}  // namespace gryphon