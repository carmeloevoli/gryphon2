#include "gryphon/kernel/DiffusionLossesKernel.h"

#include <cassert>
#include <cmath>
#include <limits>

#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace kernel {

double DiffusionLossesKernel::lambda2(double E, double Es) const {
  assert(Es >= E);
  const auto x = E / m_E0;
  return 4. * m_D0 * m_E0 / m_b0 / (1. - m_delta) *
         std::pow(x, m_delta - 1.) *
         (-std::expm1((m_delta - 1.) * std::log(Es / E)));
}

double DiffusionLossesKernel::tau(double E, double Es) const {
  assert(Es >= E);
  return m_E0 / m_b0 * (m_E0 / E - m_E0 / Es);
}

double DiffusionLossesKernel::Estar(double E, double dt) const {
  assert(dt >= 0.);
  if (dt == 0.) return E;
  auto value = m_E0 / (m_E0 / E - dt * m_b0 / m_E0);
  return (value < 0.) ? 1e10 * E : value;
}

FluxContribution DiffusionLossesKernel::contribution(
    double E, double dt, const utils::Vector3d& pos,
    const InjectionSpectrum& injection) const {
  if (dt <= 0.) return {};

  const auto tLoss = tau(E, std::numeric_limits<double>::infinity());
  if (dt >= tLoss) return {};

  const auto Es = Estar(E, dt);
  if (!std::isfinite(Es) || Es <= E) return {};

  const auto lambda2Value = lambda2(E, Es);
  if (lambda2Value <= 0.) return {};

  const auto d2 = pow2(pos.x) + pow2(pos.y);
  auto prefactor = injection(Es) / std::pow(M_PI * lambda2Value, 1.5);
  prefactor *= b(Es) / b(E);
  prefactor *= std::exp(-d2 / lambda2Value);
  prefactor *= cgs::c_light / 4. / M_PI;

  const double halo = utils::halo_function(lambda2Value, m_H, 0., pos.z);
  const double halo_gradient = utils::halo_function_dz(lambda2Value, m_H, 0., pos.z);

  FluxContribution result;
  result.flux = prefactor * halo;
  result.gradient.x = result.flux * 2. * pos.x / lambda2Value;
  result.gradient.y = result.flux * 2. * pos.y / lambda2Value;
  result.gradient.z = prefactor * halo_gradient;
  return result;
}

}  // namespace kernel
}  // namespace gryphon
