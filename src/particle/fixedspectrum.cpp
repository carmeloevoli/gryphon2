#include "gryphon/particle/fixedspectrum.h"

#include <gsl/gsl_sf_gamma.h>

namespace gryphon {
namespace particle {

double FixedSpectrumParticle::compute_normalization() const {
  const double x_c = m_E0 / m_Emax;
  const double total_energy = m_epsilon * cgs::E_SN;
  auto value = total_energy * std::pow(x_c, 2. - m_alpha);
  value /= pow2(m_E0) * gsl_sf_gamma_inc(2. - m_alpha, x_c);
  return value;
}

double FixedSpectrumParticle::get(double E, double dt, utils::Vector3d pos) const {
  const auto lambda_2 = 4. * m_D.get(E) * dt;
  auto value = m_Q0 * std::pow(E / m_E0, -m_alpha);
  if (m_Emax > 0.) value *= std::exp(-(E / m_Emax));
  value /= std::pow(M_PI * lambda_2, 1.5);
  const auto d2 = pow2(pos.x) + pow2(pos.y);
  value *= std::exp(-(d2 / lambda_2));
  // value *= Utilities::get_summation_on_z(lambda_2, m_H, r_obs.z, r_source.z, 10);
  // const auto n_bar = gas_density * 2. * h_disk / m_H;
  // const auto spallation_rate = n_bar * cgs::c_light * sigma_in(E, m_pid);
  // value *= std::exp(-spallation_rate * dt);
  return value;
}

}  // namespace particle
}  // namespace gryphon