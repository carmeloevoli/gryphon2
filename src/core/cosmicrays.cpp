#include "gryphon/core/cosmicrays.h"

#include <stdexcept>

#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace core {

CosmicRays::CosmicRays(const core::Input& input,
                       const std::shared_ptr<const kernel::GreenKernel> kernel,
                       const std::shared_ptr<const injection::InjectionSpectrum> injection,
                       const Events& events)
    : m_kernel(kernel), m_injection(injection), m_events(events) {
  if (!m_kernel) throw std::invalid_argument("GreenKernel cannot be null");
  if (!m_injection) throw std::invalid_argument("InjectionSpectrum cannot be null");
  m_E = utils::LogAxis<double>(input.E_min(), input.E_max(), input.E_size());
  m_I.assign(m_E.size(), 0.);
}

void CosmicRays::run() {
  std::fill(m_I.begin(), m_I.end(), 0.);

  for (size_t i = 0; i < m_E.size(); ++i) {
    const auto E = m_E[i];
    double flux = 0.;

    for (const auto& event : m_events) {
      if (!event) continue;
      flux += m_kernel->flux(E, event->age, event->pos,
                             [this](double Eprime) { return m_injection->get(Eprime); });
    }

    m_I[i] = flux;
  }
}

}  // namespace core
}  // namespace gryphon
