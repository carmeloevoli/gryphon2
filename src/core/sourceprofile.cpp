#include "gryphon/core/sourceprofile.h"

#include <cmath>

namespace gryphon {
namespace core {

SourceProfile::SourceProfile(const Input& input)
    : m_a(input.a), m_b(input.b), m_R1(input.R_1), m_Rsun(input.R_sun) {}

double SourceProfile::get(const double& r) const {
  double value = std::pow((r + m_R1) / (m_Rsun + m_R1), m_a);
  value *= std::exp(-m_b * (r - m_Rsun) / (m_Rsun + m_R1));
  return value;
}

}  // namespace core
}  // namespace gryphon
