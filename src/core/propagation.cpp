#include "gryphon/core/propagation.h"

#include <cassert>

#include "gryphon/utils/logging.h"

namespace gryphon {
namespace core {

Propagation::Propagation(const Input& in) {
  m_E0 = in.E_0;
  m_D0 = in.D0_over_H * in.H;
  m_H = in.H;
  m_delta = in.delta;
  auto U_CMB = 0.25 * cgs::eV / cgs::cm3;
  auto U_B = pow2(in.B_field) / 8. / M_PI;
  auto gamma2 = pow2(m_E0 / cgs::electron_mass_c2);
  m_b0 = 4. / 3. * cgs::c_light * cgs::sigma_th * (U_CMB + U_B) * gamma2;
}

double Propagation::lambda2(double E, double E_s) const {
  assert(E_s >= E);
  const auto x = E / m_E0;
  const auto x_s = E_s / m_E0;
  return 4. * m_D0 * m_E0 / m_b0 / (1. - m_delta) *
         (std::pow(x, m_delta - 1.) - std::pow(x_s, m_delta - 1.));
}

double Propagation::tau(double E, double E_s) const {
  assert(E_s >= E);
  return m_E0 / m_b0 * (m_E0 / E - m_E0 / E_s);
}

double Propagation::Estar(double E, double dt) const {
  assert(dt >= 0);
  if (dt == 0.) return E;
  auto value = m_E0 / (m_E0 / E - dt * m_b0 / m_E0);
  return (value < 0.) ? 1e10 * E : value;
}

}  // namespace core
}  // namespace gryphon