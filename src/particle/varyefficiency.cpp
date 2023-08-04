#include "gryphon/particle/varyefficiency.h"

#include "gryphon/utils/logging.h"
#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace particle {

VaryEfficiencyParticle::VaryEfficiencyParticle(const core::Input& in, RandomNumberGenerator& rng)
    : Particle(in), m_rng(rng), m_varEfficienvy(in.injEfficiencyVariance) {
  m_Q0 = source_normalization();
}

double VaryEfficiencyParticle::source_normalization() const {
  double value = m_epsilon * cgs::E_SN / pow2(m_E0);
  value *= (m_alpha - 2.);
  if (m_Emax > 0.) value /= 1. - std::pow(m_Emax / m_E0, 2. - m_alpha);
  return value;
}

double VaryEfficiencyParticle::Q(double E) const {
  auto value = m_Q0 * std::pow(E / m_E0, -m_alpha);
  if (m_Emax > 0.) value *= std::exp(-(E / m_Emax));
  return m_rng.normal(1., m_varEfficienvy) * value;
}

double VaryEfficiencyParticle::get(double E, double dt, utils::Vector3d pos) const {
  const auto lambda_2 = 4. * m_D.get(E) * dt;
  auto value = Q(E);
  value /= std::pow(M_PI * lambda_2, 1.5);
  const auto d2 = pow2(pos.x) + pow2(pos.y);
  value *= std::exp(-(d2 / lambda_2));
  value *= utils::halo_function(lambda_2, m_H, 0., pos.z);
  return cgs::c_light / 4. / M_PI * value;
}

}  // namespace particle
}  // namespace gryphon