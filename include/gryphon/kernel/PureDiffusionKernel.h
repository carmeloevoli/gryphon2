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
      : m_D0(in.D0_over_H * in.H), m_E0(in.E_0), m_delta(in.delta), m_H(in.H) {}

  double flux(double E, double dt, const utils::Vector3d& pos,
              const InjectionSpectrum& injection) const override;

  inline double D(double E) const { return m_D0 * std::pow(E / m_E0, m_delta); }

 private:
  double m_D0;
  double m_E0;
  double m_delta;
  double m_H;
};

}  // namespace kernel
}  // namespace gryphon

#endif  // GRYPHON_KERNEL_PUREDIFFUSIONKERNEL_H
