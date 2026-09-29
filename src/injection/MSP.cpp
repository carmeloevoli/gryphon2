#include "gryphon/injection/MSP.h"

#include <cmath>
#include <stdexcept>

#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace injection {

MSPInjection::MSPInjection(const core::Input& input) : m_input(input) {
  input.validate();
  if (input.injectionModel() != InjectionModel::MSP)
    throw std::invalid_argument("MSPInjection requires injectionmodel=MSP");
  const double lo = input.mspEmin()/input.injEmax();
  // Exact exponential normalization to a negligible exp(-80) tail. Logarithmic
  // integration remains well behaved for steep indices and small Emin/Ecut.
  const double integral = utils::QAGIntegration<double>([&](double logx) {
    const double x = std::exp(logx);
    return std::exp((2.-input.injSlope())*logx-x);
  }, std::log(lo), std::log(lo+80.), 1000, 1e-10);
  m_norm = std::pow(input.E_0(), input.injSlope()) *
           std::pow(input.injEmax(), 2.-input.injSlope()) * integral;
  if (!(m_norm > 0.) || !std::isfinite(m_norm))
    throw std::runtime_error("Nonfinite MSP injection normalization");
}

MSPSpin MSPInjection::sample(RandomNumberGenerator& rng) const {
  const double period = m_input.mspPeriodMin()/(1.-rng());
  const double logb = m_input.mspSigmaLog10B() == 0. ? m_input.mspLog10B() :
      rng.normal(m_input.mspLog10B(), m_input.mspSigmaLog10B());
  const double field = std::pow(10., logb)*cgs::gauss;
  const double pdot = pow2(field/(3.2e19*cgs::gauss))*cgs::second/period;
  const double luminosity = 4.*M_PI*M_PI*m_input.mspInertia()*pdot/pow3(period);
  if (!std::isfinite(luminosity) || !(luminosity > 0.))
    throw std::runtime_error("Nonfinite sampled MSP power; check field dispersion");
  return {period, field, pdot, luminosity};
}

double MSPInjection::rate(double energy, double spinDownPower) const {
  if (energy < m_input.mspEmin()) return 0.;
  return .5*m_input.injEfficiency()*spinDownPower/m_norm *
         std::pow(energy/m_input.E_0(), -m_input.injSlope()) *
         std::exp(-energy/m_input.injEmax());
}

}  // namespace injection
}  // namespace gryphon
