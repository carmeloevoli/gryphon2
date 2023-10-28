#include "gryphon/particle/msp.h"

#include "gryphon/utils/logging.h"
#include "gryphon/utils/numeric.h"

// #define BURSTED

namespace gryphon {
namespace particle {

MSP::MSP(const core::Input& in) : Particle(in.pid, std::make_shared<core::Event>()) {}

MSP::MSP(const core::Input& in, const std::shared_ptr<core::Event>& event,
         RandomNumberGenerator& rng)
    : Particle(in.pid, event) {}

double MSP::Q(double E) const {
  if (E < cgs::MeV || E > 5. * m_Ec) return 0.;
  const auto x = m_Ec / m_E0;
  assert(x > 1);
  const auto I = (std::pow(x, 2. - m_alpha) - 1.) / (2. - m_alpha);
  const auto Q0 = m_Luminosity / pow2(m_E0) / I;
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
  return (value < 0.) ? 1e10 * E : value;
}

#ifdef BURSTED
double MSP::get(double E, double dt, utils::Vector3d pos) const {
  if (dt <= 0.) return 0.;
  if (dt >= pow2(m_E0) / m_b0 / E) return 0.;
  double value = 0;
  const auto E_s = Estar(E, dt);
  if (E_s < 100. * cgs::PeV) {
    const auto lambda_2 = lambda2(E, E_s);
    const auto d2 = pos.getModuleSquared();
    value = Q(E_s) / std::pow(M_PI * lambda_2, 1.5);
    value *= b(E_s) / b(E);
    value *= std::exp(-d2 / lambda_2);
  }
  return cgs::c_light / 4. / M_PI * value;
}
#else
double MSP::get(double E, double dt, utils::Vector3d pos) const {
  if (dt <= 0.) return 0.;
  const auto dt_max = pow2(m_E0) / m_b0 / E;
  const auto d2 = pos.getModuleSquared();
  auto I = utils::QAGIntegration<double>(
      [this, E, d2, dt](double t) {
        const auto E_s = Estar(E, dt - t);
        double value = 0.;
        if (E_s < 100. * cgs::PeV) {
          const auto lambda_2 = lambda2(E, E_s);
          value = Q(E_s) / std::pow(M_PI * lambda_2, 1.5);
          value *= b(E_s) / b(E);
          if (d2 > 100. * lambda_2) return 0.;
          value *= std::exp(-(d2 / lambda_2));
        }
        return value;
      },
      std::max(0., dt - dt_max), dt, 1000, 1e-3);
  return cgs::c_light / 4. / M_PI * I;
}
#endif

double MSP::get(double E) const { return get(E, m_age, m_pos); }

}  // namespace particle
}  // namespace gryphon