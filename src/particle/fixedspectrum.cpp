#include "gryphon/particle/fixedspectrum.h"

#include "gryphon/utils/logging.h"

namespace gryphon {
namespace particle {

double FixedSpectrumParticle::compute_normalization() const {
  double value = m_epsilon * cgs::E_SN / pow2(m_E0);
  value *= (m_alpha - 2.);
  if (m_Emax > 0.) value /= 1. - std::pow(m_Emax / m_E0, 2. - m_alpha);
  return value;
}

double FixedSpectrumParticle::Q(double E) const {
  auto value = m_Q0 * std::pow(E / m_E0, -m_alpha);
  if (m_Emax > 0.) value *= std::exp(-(E / m_Emax));
  return value;
}

double G1D(double x, double sigma_squared) { return std::exp(-(x * x) / sigma_squared); }

double FixedSpectrumParticle::summation_on_z(double lambda_2, double L, double z_obs, double z_s,
                                             int max_n) const {
  double value = 0;
  for (int n = -max_n; n <= max_n; n += 1) {
    int sign = (n % 2 == 0) ? 1 : -1;
    double z_n = 2.0 * static_cast<double>(n) * L + static_cast<double>(sign) * z_s;
    value += static_cast<double>(sign) * G1D(z_obs - z_n, lambda_2);
  }
  return value;
}

double FixedSpectrumParticle::get(double E, double dt, utils::Vector3d pos) const {
  const auto lambda_2 = 4. * m_D.get(E) * dt;
  auto value = Q(E);
  value /= std::pow(M_PI * lambda_2, 1.5);
  const auto d2 = pow2(pos.x) + pow2(pos.y);
  value *= std::exp(-(d2 / lambda_2));
  value *= summation_on_z(lambda_2, m_H, 0., pos.z, 10);
  // const auto n_bar = gas_density * 2. * h_disk / m_H;
  // const auto spallation_rate = n_bar * cgs::c_light * sigma_in(E, m_pid);
  // value *= std::exp(-spallation_rate * dt);
  return cgs::c_light / 4. / M_PI * value;
}

}  // namespace particle
}  // namespace gryphon