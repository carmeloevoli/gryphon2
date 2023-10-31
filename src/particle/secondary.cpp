#include "gryphon/particle/secondary.h"

#include <cmath>

#include "gryphon/core/xsecs.h"
#include "gryphon/utils/numeric.h"

double T2beta(const double& T) {
  double beta = std::sqrt(T * (T + 2. * cgs::proton_mass_c2));
  beta /= T + cgs::proton_mass_c2;
  return beta;
}

namespace gryphon {
namespace particle {

SecondaryPositrons::SecondaryPositrons(const core::Input& in)
    : Particle(in.pid, std::make_shared<core::Event>()) {
  m_E0 = in.E_0;
  m_H = in.H;
  m_D0 = in.D0_over_H * in.H;
  m_delta = in.delta;
  m_Rg = in.R_g;
}

double SecondaryPositrons::Q(double E) const { return Q_pp(E) + Q_Hep(E); }

double SecondaryPositrons::Q_pp(double E) const {
  auto I = utils::QAGIntegration<double>(
      [this, E](double lnEprime) {
        const auto Eprime = std::exp(lnEprime);
        auto beta = T2beta(Eprime);
        auto value = beta * H_lis.get(Eprime) * posXsecs.get(Orusa2022::Channel::pp, Eprime, E);
        return Eprime * value;
      },
      std::log(E), std::log(cgs::PeV), 1000, 1e-3);
  return 4. * M_PI * I;
}

double SecondaryPositrons::Q_Hep(double E) const {
  auto I = utils::QAGIntegration<double>(
      [this, E](double lnEprime) {
        const auto Eprime = std::exp(lnEprime);
        auto beta = T2beta(Eprime);
        auto value = beta * He_lis.get(Eprime) * posXsecs.get(Orusa2022::Channel::Hep, Eprime, E);
        return Eprime * value;
      },
      std::log(E), std::log(cgs::PeV), 1000, 1e-3);
  return 4. * M_PI * I;
}

double SecondaryPositrons::b(double E) const { return m_b0 * std::pow(E / m_E0, 2.); }

double SecondaryPositrons::lambda2(double E, double E_s) const {
  assert(E_s >= E);
  const auto x = E / m_E0;
  const auto x_s = E_s / m_E0;
  return 4. * m_D0 * m_E0 / m_b0 / (1. - m_delta) *
         (std::pow(x, m_delta - 1.) - std::pow(x_s, m_delta - 1.));
}

double SecondaryPositrons::tau(double E, double E_s) const {
  assert(E_s >= E);
  return m_E0 / m_b0 * (m_E0 / E - m_E0 / E_s);
}

double SecondaryPositrons::Estar(double E, double dt) const {
  assert(dt >= 0);
  if (dt == 0.) return E;
  auto value = m_E0 / (m_E0 / E - dt * m_b0 / m_E0);
  return (value < 0.) ? 1e10 * E : value;
}

double SecondaryPositrons::get(double E) const {
  const auto dt_max = pow2(m_E0) / m_b0 / E;
  const auto factor = cgs::c_light / 4. / std::pow(M_PI, 1.5) * m_mu / m_ismMass;
  auto I = utils::QAGIntegration<double>(
      [this, E, dt_max](double t) {
        const auto E_s = Estar(E, dt_max - t);
        double value = 0.;
        if (E_s < 100. * cgs::PeV) {
          const auto lambda_2 = lambda2(E, E_s);
          value = Q(E_s) / std::sqrt(lambda_2);
          value *= b(E_s) / b(E);
          value *= 1. - std::exp(-pow2(m_Rg) / lambda_2);
          value *= utils::halo_function(lambda_2, m_H, 0., 0.);
        }
        return value;
      },
      0., dt_max, 1000, 1e-3);
  return factor * I;
}

}  // namespace particle
}  // namespace gryphon