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
  const auto xs = Es / m_E0;
  return 4. * m_D0 * m_E0 / m_b0 / (1. - m_delta) *
         (std::pow(x, m_delta - 1.) - std::pow(xs, m_delta - 1.));
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

double DiffusionLossesKernel::flux(double E, double dt, const utils::Vector3d& pos,
                                   const InjectionSpectrum& injection) const {
  if (dt <= 0.) return 0.;

  const auto tLoss = tau(E, std::numeric_limits<double>::infinity());
  if (dt >= tLoss) return 0.;

  const auto Es = Estar(E, dt);
  if (!std::isfinite(Es) || Es <= E) return 0.;

  const auto lambda2Value = lambda2(E, Es);
  if (lambda2Value <= 0.) return 0.;

  const auto d2 = pos.getModuleSquared();
  auto value = injection(Es) / std::pow(M_PI * lambda2Value, 1.5);
  value *= b(Es) / b(E);
  value *= std::exp(-d2 / lambda2Value);
  value *= utils::halo_function(lambda2Value, m_H, 0., pos.z);

  return cgs::c_light / 4. / M_PI * value;
}

}  // namespace kernel
}  // namespace gryphon
