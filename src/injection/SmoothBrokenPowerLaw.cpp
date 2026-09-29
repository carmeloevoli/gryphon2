#include "gryphon/injection/SmoothBrokenPowerLaw.h"

#include <cmath>
#include <stdexcept>

#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace injection {

SmoothBrokenPowerLawSpectrum::SmoothBrokenPowerLawSpectrum(const core::Input& input)
    : InjectionSpectrum(input),
      m_E0(input.E_0()),
      m_Emin(1. * cgs::GeV),
      m_Emax(input.injEmax()),
      m_Ebreak(input.injBreakEnergy()),
      m_alpha(input.injSlope()),
      m_deltaAlpha(input.injDeltaSlope()),
      m_smoothness(input.injSmoothness()),
      m_crEnergy(input.injEfficiency() * cgs::E_SN),
      m_Q0(sourceNormalization()) {}

double SmoothBrokenPowerLawSpectrum::shape(double E) const {
  if (E < m_Emin || E >= m_Emax) return 0.;
  const double transition =
      std::pow(E / m_Ebreak, m_deltaAlpha / m_smoothness);
  return std::pow(E / m_E0, -m_alpha) * std::pow(1. + transition, -m_smoothness);
}

double SmoothBrokenPowerLawSpectrum::sourceNormalization() const {
  if (!(m_Emax > m_Emin)) {
    throw std::invalid_argument(
        "SmoothBrokenPowerLawSpectrum requires injEmax > 1 GeV for normalization");
  }

  auto energyPerLogBin = [this](double logEnergy) {
    const double energy = std::exp(logEnergy);
    return energy * energy * shape(energy);
  };
  const double integral = utils::simpsonIntegration<double>(
      energyPerLogBin, std::log(m_Emin), std::log(m_Emax), 4096);
  if (!(integral > 0.) || !std::isfinite(integral)) {
    throw std::runtime_error("Could not normalize SmoothBrokenPowerLawSpectrum");
  }
  return m_crEnergy / integral;
}

double SmoothBrokenPowerLawSpectrum::get(double E) const { return m_Q0 * shape(E); }

}  // namespace injection
}  // namespace gryphon
