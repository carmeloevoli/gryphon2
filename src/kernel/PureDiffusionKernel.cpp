#include "gryphon/kernel/PureDiffusionKernel.h"

#include <cmath>

#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace kernel {

namespace {

double softplus(double value) {
  return std::max(0., value) + std::log1p(std::exp(-std::abs(value)));
}

}  // namespace

double PureDiffusionKernel::D(double E) const {
  const double log_energy_ratio = std::log(E / m_E0);
  double log_diffusion = std::log(m_D0) + m_delta * log_energy_ratio;

  // ddelta = -1 is retained as the legacy "no break" sentinel.  Otherwise
  // delta changes smoothly by ddelta around Eb, while D(E0) remains D0.
  if (m_ddelta != -1.) {
    const double transition = m_ddelta / m_s;
    log_diffusion +=
        m_s * (softplus(transition * std::log(E / m_Eb)) -
               softplus(transition * std::log(m_E0 / m_Eb)));
  }
  return std::exp(log_diffusion);
}

FluxContribution PureDiffusionKernel::contribution(
    double E, double dt, const utils::Vector3d& pos,
    const InjectionSpectrum& injection) const {
  if (dt < cgs::t_ST) return {};

  const auto lambda2 = 4. * D(E) * dt;
  if (!(lambda2 > 0.)) return {};

  auto prefactor = injection(E);
  prefactor /= std::pow(M_PI * lambda2, 1.5);

  const auto d2 = pow2(pos.x) + pow2(pos.y);
  prefactor *= std::exp(-(d2 / lambda2));
  prefactor *= cgs::c_light / 4. / M_PI;

  const double halo = utils::halo_function(lambda2, m_H, 0., pos.z);
  const double halo_gradient = utils::halo_function_dz(lambda2, m_H, 0., pos.z);

  FluxContribution result;
  result.flux = prefactor * halo;
  result.gradient.x = result.flux * 2. * pos.x / lambda2;
  result.gradient.y = result.flux * 2. * pos.y / lambda2;
  result.gradient.z = prefactor * halo_gradient;

  return result;
}

}  // namespace kernel
}  // namespace gryphon
