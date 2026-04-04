#ifndef GRYPHON_KERNEL_GREENKERNEL_H
#define GRYPHON_KERNEL_GREENKERNEL_H

#include <functional>

#include "gryphon/core/cgs.h"
#include "gryphon/utils/vector3.h"

namespace gryphon {
namespace kernel {

using InjectionSpectrum = std::function<double(double)>;

class GreenKernel {
 public:
  virtual ~GreenKernel() = default;

  virtual double flux(double E, double dt, const utils::Vector3d& pos,
                      const InjectionSpectrum& injection) const = 0;
};

}  // namespace kernel
}  // namespace gryphon

#endif  // GRYPHON_KERNEL_GREENKERNEL_H
