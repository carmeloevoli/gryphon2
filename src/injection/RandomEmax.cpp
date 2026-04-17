#include "gryphon/injection/RandomEmax.h"

#include <cmath>
#include <stdexcept>

#include "gryphon/utils/numeric.h"

namespace gryphon {

namespace {

double drawEmax(const core::Input& in, RandomNumberGenerator& rng) {
  // Interpret injEmax as the median cutoff of the underlying log-normal draw.
  const double characteristic_emax = in.injEmax();
  const double sigma_dex = in.injEmaxSigmaDex();
  const double emax_min = in.injEmaxMin();
  const double emax_max = in.injEmaxMax();

  if (sigma_dex <= 0.) return characteristic_emax;

  const double sigma_ln = sigma_dex * std::log(10.);
  const double mean_ln = std::log(characteristic_emax);

  for (size_t draw = 0; draw < 1024; ++draw) {
    const double candidate = std::exp(rng.normal(mean_ln, sigma_ln));
    if (candidate >= emax_min && candidate <= emax_max) return candidate;
  }

  throw std::runtime_error("Failed to draw RandomEmax cutoff within configured bounds");
}

}  // namespace

namespace injection {

RandomEmaxSpectrum::RandomEmaxSpectrum(const core::Input& in, RandomNumberGenerator& rng)
    : InjectionSpectrum(in), m_E0(in.E_0()), m_Emin(1. * cgs::GeV) {
  m_alpha = in.injSlope();
  m_crenergy = in.injEfficiency() * cgs::E_SN;
  m_Emax = drawEmax(in, rng);
  m_Q0 = source_normalization();
}

double RandomEmaxSpectrum::source_normalization() const {
  const double prefactor = m_crenergy / pow2(m_E0);
  const double dalpha = m_alpha - 2.;
  const double tol = 1e-12;

  if (m_Emax <= m_Emin) {
    throw std::invalid_argument("RandomEmaxSpectrum requires Emax > 1 GeV for normalization");
  }

  const double min_ratio = m_Emin / m_E0;
  const double max_ratio = m_Emax / m_E0;

  if (std::abs(dalpha) <= tol) {
    return prefactor / std::log(max_ratio / min_ratio);
  }

  const double exponent = 2. - m_alpha;
  return prefactor * exponent / (std::pow(max_ratio, exponent) - std::pow(min_ratio, exponent));
}

double RandomEmaxSpectrum::get(double E) const {
  if (E < m_Emin || E >= m_Emax) return 0.;

  const double value = m_Q0 * std::pow(E / m_E0, -m_alpha);
  return value;
}

}  // namespace injection
}  // namespace gryphon
