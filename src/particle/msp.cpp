#include "gryphon/particle/msp.h"

#include "gryphon/utils/logging.h"
#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace particle {

MSP::MSP(const core::Input& in) : Particle(in.pid, std::make_shared<core::Event>()) {
  LOGD << escapeTimescale() / cgs::Myr;
  LOGD << lossesTimescale() / cgs::Myr;
}

MSP::MSP(const core::Input& in, const std::shared_ptr<core::Event>& event)
    : Particle(in.pid, event) {}

double MSP::Q(double E) const {
  const auto x = m_Ec / m_E0;
  const auto I = (std::pow(x, 2. - m_alpha) - 1.) / (2. - m_alpha);
  const auto Q0 = (x >= 1) ? m_Luminosity / pow2(m_E0) / I : 0;
  return (E >= m_E0) ? Q0 * std::pow(E / m_E0, -m_alpha) * std::exp(-E / m_Ec) : 0;
}

double MSP::b(double E) const { return m_b0 * std::pow(E / m_E0, 2.); }

double MSP::lambda2(double E, double E_s) const {
  assert(E_s >= E);
  const auto x = E / m_E0;
  const auto x_s = E_s / m_E0;
  return 4. * m_D0 * m_E0 / m_b0 / (1. - m_delta) *
         (std::pow(x, m_delta - 1.) - std::pow(x_s, m_delta - 1.));
}

double MSP::get(double E, double dt, utils::Vector3d pos) const {
  if (dt <= 0.) return 0.;
  auto I = utils::QAGIntegration<double>(
      [&](double t) {
        const auto E_s = 1.;  // TBD
        return b(E_s) / std::pow(M_PI * lambda2(E, E_s), 1.5) * Q(E_s);
      },
      0., dt, 1000, 1e-6);
  double value = I / b(E);
  return value;
}

double MSP::get(double E) const { return get(E, m_age, m_pos); }

}  // namespace particle
}  // namespace gryphon