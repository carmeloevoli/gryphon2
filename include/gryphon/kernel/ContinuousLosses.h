#ifndef GRYPHON_KERNEL_CONTINUOUSLOSSES_H
#define GRYPHON_KERNEL_CONTINUOUSLOSSES_H

#include <vector>

#include "gryphon/injection/MSP.h"
#include "gryphon/kernel/DiffusionLossesKernel.h"

namespace gryphon {
namespace kernel {

// Steady point sources. The free-space energy integrals are cached, then the
// slab images are summed at each actual source position. All units are CGS.
class ContinuousLosses {
 public:
  explicit ContinuousLosses(const core::Input& input);
  // Result per unit spin-down power (1 erg/s); multiply by source luminosity.
  double response(size_t energyIndex, const utils::Vector3d& position) const;
  // Independent direct energy integral using gryphon's halo_function.
  double directResponse(double energy, const utils::Vector3d& position,
                        size_t quadrature = 512) const;
  const std::vector<double>& energies() const { return m_energies; }
  int imageOrder() const { return m_images; }

 private:
  core::Input m_input;
  DiffusionLossesKernel m_kernel;
  injection::MSPInjection m_injection;
  std::vector<double> m_energies;
  std::vector<double> m_logDistances;
  std::vector<std::vector<double>> m_response;
  int m_images;
  double interpolate(size_t energyIndex, double distance) const;
};

}  // namespace kernel
}  // namespace gryphon
#endif
