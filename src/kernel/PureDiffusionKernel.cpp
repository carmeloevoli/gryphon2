#include "gryphon/kernel/PureDiffusionKernel.h"

#include <cmath>

#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace kernel {

double PureDiffusionKernel::flux(double E, double dt, const utils::Vector3d& pos,
                                 const InjectionSpectrum& injection) const {
  if (dt < cgs::t_ST) return 0.;

  const auto lambda2 = 4. * D(E) * dt;
  auto value = injection(E);
  value /= std::pow(M_PI * lambda2, 1.5);

  const auto d2 = pow2(pos.x) + pow2(pos.y);
  value *= std::exp(-(d2 / lambda2));
  value *= utils::halo_function(lambda2, m_H, 0., pos.z);

  return cgs::c_light / 4. / M_PI * value;
}

}  // namespace kernel
}  // namespace gryphon
