#include "gryphon/particle/msp.h"

#include "gryphon/utils/logging.h"
#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace particle {

MSP::MSP(const core::Input& in) : Particle(in.pid, std::make_shared<core::Event>()) {}

MSP::MSP(const core::Input& in, const std::shared_ptr<core::Event>& event,
         RandomNumberGenerator& rng)
    : Particle(in.pid, event) {}

double MSP::Q(double E) const {
  if (E < cgs::MeV) return 0.;
  const auto x = m_Ec / m_E0;
  const auto I = (std::pow(x, 2. - m_alpha) - 1.) / (2. - m_alpha);
  const auto Q0 = (x >= 1) ? m_Luminosity / pow2(m_E0) / I : 0;
  return Q0 * std::pow(E / m_E0, -m_alpha) * std::exp(-E / m_Ec);
}

double MSP::b(double E) const { return m_b0 * std::pow(E / m_E0, 2.); }

double MSP::lambda2(double E, double E_s) const {
  assert(E_s >= E);
  const auto x = E / m_E0;
  const auto x_s = E_s / m_E0;
  return 4. * m_D0 * m_E0 / m_b0 / (1. - m_delta) *
         (std::pow(x, m_delta - 1.) - std::pow(x_s, m_delta - 1.));
}

double MSP::tau(double E, double E_s) const {
  assert(E_s >= E);
  return m_E0 / m_b0 * (m_E0 / E - m_E0 / E_s);
}

double MSP::Estar(double E, double dt) const {
  assert(dt >= 0);
  if (dt == 0.) return E;
  auto value = m_E0 / (m_E0 / E - dt * m_b0 / m_E0);
  return (value < 0.) ? 1e4 * E : value;
}

// double MSP::get(double E, double dt, utils::Vector3d pos) const {
//   if (dt <= 0.) return 0.;
//   if (dt >= pow2(m_E0) / m_b0 / E) return 0.;
//   const auto E_s = Estar(E, dt);
//   const auto lambda_2 = lambda2(E, E_s);
//   const auto d2 = pos.getModuleSquared();
//   auto value = Q(E_s) / std::pow(M_PI * lambda_2, 1.5);
//   value *= b(E_s) / b(E);
//   value *= std::exp(-d2 / lambda_2);
//   return cgs::c_light / 4. / M_PI * value;
// }

double MSP::get(double E, double dt, utils::Vector3d pos) const {
  if (dt <= 0.) return 0.;
  if (dt >= pow2(m_E0) / m_b0 / E) return 0.;
  auto I = utils::QAGIntegration<double>(
      [this, E, pos](double t) {
        const auto E_s = Estar(E, t);
        const auto lambda_2 = lambda2(E, E_s);
        auto value = Q(E_s) / std::pow(M_PI * lambda_2, 1.5);
        value *= b(E_s) / b(E);
        const auto d2 = pos.getModuleSquared();
        value *= std::exp(-(d2 / lambda_2));
        return value;
      },
      0., dt, 1000, 1e-2);
  return cgs::c_light / 4. / M_PI * I;
}

double MSP::get(double E) const { return get(E, m_age, m_pos); }

}  // namespace particle
}  // namespace gryphon