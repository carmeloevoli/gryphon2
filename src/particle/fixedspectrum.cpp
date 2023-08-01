#include "gryphon/particle/fixedspectrum.h"

// class Particle {
//  public:
//   Particle(double slope, double Emax, const DiffusionCoefficient& D) {}
//   virtual ~Particle() = default;

//   double get(double E, double dt, Vector3d pos) {
//     const auto lambda_2 = 4. * m_D(E) * dt;
//     auto value = m_Q.get(E);
//     value /= std::pow(M_PI * lambda_2, 1.5);
//     const auto d2 = pow2(pos.x) + pow2(pos.y);
//     value *= std::exp(-(d2 / lambda_2));
//     // value *= Utilities::get_summation_on_z(lambda_2, m_H, r_obs.z, r_source.z, 10);
//     // const auto n_bar = gas_density * 2. * h_disk / m_H;
//     // const auto spallation_rate = n_bar * cgs::c_light * sigma_in(E, m_pid);
//     // value *= std::exp(-spallation_rate * dt);
//     return value;
//   }
// }

// #endif