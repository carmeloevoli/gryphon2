#ifndef GRYPHON_KERNEL_GREENKERNEL_H
#define GRYPHON_KERNEL_GREENKERNEL_H

#include <functional>

#include "gryphon/core/cgs.h"
#include "gryphon/utils/vector3.h"

namespace gryphon {
namespace kernel {

using InjectionSpectrum = std::function<double(double)>;

// Contribution of one source at the observer.  The gradient is taken with
// respect to the observer position in the same Cartesian frame as Event::pos.
// Since flux is proportional to density, grad(flux)/flux = grad(n)/n.
struct FluxContribution {
  double flux = 0.;
  utils::Vector3d gradient;
};

class GreenKernel {
 public:
  virtual ~GreenKernel() = default;

  virtual FluxContribution contribution(double E, double dt, const utils::Vector3d& pos,
                                        const InjectionSpectrum& injection) const = 0;

  double flux(double E, double dt, const utils::Vector3d& pos,
              const InjectionSpectrum& injection) const {
    return contribution(E, dt, pos, injection).flux;
  }

  virtual double diffusionCoefficient(double E) const = 0;

  virtual double diffusionTimescale(double E) const = 0;

  virtual double energyLossTimescale(double E) const = 0;
};

}  // namespace kernel
}  // namespace gryphon

#endif  // GRYPHON_KERNEL_GREENKERNEL_H
