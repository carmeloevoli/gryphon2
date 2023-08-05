#include "gryphon/core/sourceprofile.h"

#include <cmath>

#include "gryphon/utils/logging.h"
#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace core {

SourceProfile::SourceProfile(const Input& input)
    : m_a(input.a), m_b(input.b), m_R1(input.R_1), m_Rsun(input.R_sun), m_Rg(input.R_g) {
  auto Imax = integrateProfile(m_Rg);
  LOGD << "normalization factor : " << Imax / pow2(cgs::kpc);
  m_profileIntegral.cacheTable([this, Imax](double r) { return integrateProfile(r) / Imax; },
                               {0, m_Rg});
}

double SourceProfile::get(const double& r) const {
  double value = std::pow((r + m_R1) / (m_Rsun + m_R1), m_a);
  value *= std::exp(-m_b * (r - m_Rsun) / (m_Rsun + m_R1));
  return value;
}

double SourceProfile::integrateProfile(const double& r) const {
  return utils::QAGIntegration<double>([&](double x) { return x * get(x); }, 0., r, 1000, 1e-6);
}

double SourceProfile::pick_radius(const double& rnd) const {
  return m_profileIntegral.getInverse(rnd);
}

}  // namespace core
}  // namespace gryphon
