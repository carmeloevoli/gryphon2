#ifndef GRYPHON_CORE_COSMICRAYS_H
#define GRYPHON_CORE_COSMICRAYS_H

#include <memory>
#include <vector>

#include "gryphon/core/event.h"
#include "gryphon/core/input.h"
#include "gryphon/injection/InjectionSpectrum.h"
#include "gryphon/kernel/greenkernel.h"

namespace gryphon {
namespace core {

// using EventInjectionSpectrum = std::function<double(double, const Event&)>;

class CosmicRays {
 public:
  CosmicRays(const core::Input& input, const std::shared_ptr<const kernel::GreenKernel> kernel,
             injection::InjectionSpectra injections,
             const Events& events);

  void run();

  const std::vector<double>& get_energyAxis() const { return m_E; }
  const std::vector<double>& get_flux() const { return m_I; }

 private:
  std::vector<double> m_E;
  std::vector<double> m_I;
  std::shared_ptr<const kernel::GreenKernel> m_kernel;
  injection::InjectionSpectra m_injections;
  Events m_events;
};

}  // namespace core
}  // namespace gryphon

#endif  // GRYPHON_CORE_COSMICRAYS_H
