#include "gryphon/particle/hpwne.h"

// #include "gryphon/core/pid.h"
// #include "gryphon/utils/logging.h"
#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace particle {

HPWNE::HPWNE(const core::Input& in)
    : Particle(in.pid, std::make_shared<core::Event>()), m_prop(core::Propagation(in)) {
  //   m_alpha = in.injSlope;
  //   m_crenergy = in.injEfficiency * cgs::E_SN;
  //   m_Emax = in.injEmax;
  m_H = in.H;
  //   m_Q0 = source_normalization();
}

// double pickSnEnergy(RandomNumberGenerator& rng) {
//   auto logEnergy = rng.normal(std::log10(3e50), 0.54);  // adimensional
//   return std::max(std::pow(10., logEnergy) * cgs::erg, 0.);
// }

// double pickSlope(double mean, double sdev, RandomNumberGenerator& rng) {
//   auto slope = rng.normal(mean, sdev);
//   return std::max(2.01, slope);
// }

HPWNE::HPWNE(const core::Input& in, const std::shared_ptr<core::Event>& event,
             RandomNumberGenerator& rng)
    : Particle(in.pid, event), m_prop(core::Propagation(in)) {
  //   m_alpha = (in.doVarySlope) ? pickSlope(in.injSlope, in.injSlopeSigma, rng) : in.injSlope;
  //   m_crenergy = in.injEfficiency * ((in.doVaryEnergy) ? pickSnEnergy(rng) : cgs::E_SN);
  //   m_Emax = in.injEmax;
  //   m_H = in.H;
  //   m_Q0 = source_normalization();
}

// double FixedSpectrumParticle::source_normalization() const {
//   double value = m_crenergy / pow2(m_E0);
//   value *= (m_alpha - 2.);
//   if (m_Emax > 0.) value /= 1. - std::pow(m_Emax / m_E0, 2. - m_alpha);
//   return value;
// }

double HPWNE::Q(double E) const {
  auto value = 0;
  return value;
}

double HPWNE::get(double E, double dt, utils::Vector3d pos) const {
  const auto lambda_2 = 4. * m_prop.D(E) * dt;
  auto value = 1;
  value /= std::pow(M_PI * lambda_2, 1.5);
  const auto d2 = pow2(pos.x) + pow2(pos.y);
  value *= std::exp(-(d2 / lambda_2));
  value *= utils::halo_function(lambda_2, m_H, 0., pos.z);
  return cgs::c_light / 4. / M_PI * value;
}

double HPWNE::get(double E) const { return get(E, m_age, m_pos); }

}  // namespace particle
}  // namespace gryphon