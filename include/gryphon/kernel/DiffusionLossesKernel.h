#ifndef GRYPHON_KERNEL_DIFFUSIONLOSSESKERNEL_H
#define GRYPHON_KERNEL_DIFFUSIONLOSSESKERNEL_H

#include <cmath>

#include "gryphon/core/input.h"
#include "gryphon/kernel/greenkernel.h"

namespace gryphon {
namespace kernel {

class DiffusionLossesKernel final : public GreenKernel {
 public:
  explicit DiffusionLossesKernel(const core::Input& in)
      : m_D0(in.D0_over_H() * in.H()),
        m_E0(in.E_0()),
        m_delta(in.delta()),
        m_H(in.H()),
        m_b0(4. / 3. * cgs::c_light * cgs::sigma_th *
             (in.U_rad() + pow2(in.B_field()) / 8. / M_PI) *
             pow2(in.E_0() / cgs::electron_mass_c2)) {}

  double flux(double E, double dt, const utils::Vector3d& pos,
              const InjectionSpectrum& injection) const override;

  double diffusionTimescale(double E) const override { return pow2(m_H) / (2. * D(E)); }

  double energyLossTimescale(double E) const override { return E / b(E); }

 public:
  inline double D(double E) const { return m_D0 * std::pow(E / m_E0, m_delta); }
  inline double b(double E) const { return m_b0 * pow2(E / m_E0); }
  double tau(double E, double Es) const;
  double Estar(double E, double dt) const;
  double lambda2(double E, double Es) const;

 private:
  double m_D0;
  double m_E0;
  double m_delta;
  double m_H;
  double m_b0;
};

}  // namespace kernel
}  // namespace gryphon

#endif  // GRYPHON_KERNEL_DIFFUSIONLOSSESKERNEL_H
