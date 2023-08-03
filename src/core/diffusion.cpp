#include "gryphon/core/diffusion.h"

#include <cmath>

namespace gryphon {
namespace core {

DiffusionCoefficient::DiffusionCoefficient(const Input& input)
    : m_H(input.H),
      m_D0(input.D0_over_H * input.H),
      m_E0(input.E_0),
      m_delta(input.delta),
      m_ddelta(input.ddelta),
      m_s(input.s),
      m_Eb(input.E_b) {}

double DiffusionCoefficient::get(const double& E) const {
  auto value = m_D0 * std::pow(E / m_E0, m_delta);
  if (m_ddelta > 0.) value /= std::pow(1. + std::pow(E / m_Eb, m_ddelta / m_s), m_s);
  return value;
}

double DiffusionCoefficient::escapeTimescale(const double& E) const {
  return pow2(m_H) / 2. / get(E);
}

}  // namespace core
}  // namespace gryphon