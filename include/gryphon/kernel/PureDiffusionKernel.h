#ifndef GRYPHON_KERNEL_PUREDIFFUSIONKERNEL_H
#define GRYPHON_KERNEL_PUREDIFFUSIONKERNEL_H

#include <cmath>

#include "gryphon/core/input.h"
#include "gryphon/kernel/greenkernel.h"

namespace gryphon {
namespace kernel {

class PureDiffusionKernel final : public GreenKernel {
 public:
  explicit PureDiffusionKernel(const core::Input& in)
      : m_D0(in.D0_over_H() * in.H()),
        m_E0(in.E_0()),
        m_delta(in.delta()),
        m_ddelta(in.ddelta()),
        m_s(in.s()),
        m_Eb(in.E_b()),
        m_H(in.H()) {}

  FluxContribution contribution(double E, double dt, const utils::Vector3d& pos,
                                const InjectionSpectrum& injection) const override;

  double D(double E) const;
  double diffusionCoefficient(double E) const override { return D(E); }

  double diffusionTimescale(double E) const override { return pow2(m_H) / (2. * D(E)); }

  double energyLossTimescale(double E) const override {
    return std::numeric_limits<double>::infinity();
  }

 private:
  double m_D0;
  double m_E0;
  double m_delta;
  double m_ddelta;
  double m_s;
  double m_Eb;
  double m_H;
};

}  // namespace kernel
}  // namespace gryphon

#endif  // GRYPHON_KERNEL_PUREDIFFUSIONKERNEL_H
