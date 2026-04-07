#include "gryphon/core/cosmicrays.h"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <thread>

#include "gryphon/utils/logging.h"
#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace core {

namespace {

size_t parse_positive_env_var(const char* name, size_t fallback) {
  const char* value = std::getenv(name);
  if (value == nullptr || *value == '\0') return fallback;

  errno = 0;
  char* end = nullptr;
  const unsigned long long parsed = std::strtoull(value, &end, 10);
  if (errno != 0 || end == value || (end != nullptr && *end != '\0') || parsed == 0 ||
      parsed > std::numeric_limits<size_t>::max()) {
    return fallback;
  }

  return static_cast<size_t>(parsed);
}

}  // namespace

CosmicRays::CosmicRays(const core::Input& input,
                       const std::shared_ptr<const kernel::GreenKernel> kernel,
                       injection::InjectionSpectra injections, const Events& events)
    : m_kernel(kernel), m_injections(std::move(injections)), m_events(events) {
  if (!m_kernel) throw std::invalid_argument("GreenKernel cannot be null");
  if (m_injections.size() != m_events.size()) {
    throw std::invalid_argument("Number of event spectra must match number of events");
  }
  for (size_t i = 0; i < m_events.size(); ++i) {
    if (m_events[i] && !m_injections[i]) {
      throw std::invalid_argument("InjectionSpectrum cannot be null for a non-null event");
    }
  }
  m_E = utils::LogAxis<double>(input.E_min(), input.E_max(), input.E_size());
  m_I.assign(m_E.size(), 0.);
}

void CosmicRays::run() {
  std::fill(m_I.begin(), m_I.end(), 0.);

  const size_t n_bins = m_E.size();
  if (n_bins == 0) return;

  auto compute_bin = [this](size_t i) {
    const auto E = m_E[i];
    double flux = 0.;

    for (size_t event_index = 0; event_index < m_events.size(); ++event_index) {
      const auto& event = m_events[event_index];
      if (!event) continue;
      const auto* injection = m_injections[event_index].get();
      if (injection == nullptr) continue;
      flux += m_kernel->flux(E, event->age, event->pos,
                             [injection](double Eprime) { return injection->get(Eprime); });
    }

    m_I[i] = flux;
  };

  const auto detected_threads = std::thread::hardware_concurrency();
  const size_t hardware_threads = detected_threads == 0 ? 1 : detected_threads;
  const size_t requested_threads = parse_positive_env_var("GRYPHON_CR_THREADS", hardware_threads);
  const size_t worker_count = std::min(n_bins, std::max<size_t>(1, requested_threads));
  const bool run_in_parallel = worker_count > 1;
  const size_t actual_worker_count = run_in_parallel ? worker_count : 1;

  LOGI << "CosmicRays::run using " << actual_worker_count << " worker thread"
       << (actual_worker_count == 1 ? "" : "s")
       << (run_in_parallel ? " in parallel" : " (serial execution)")
       << " [requested=" << requested_threads << ", hardware=" << hardware_threads
       << ", bins=" << n_bins << "]";

  if (!run_in_parallel) {
    for (size_t i = 0; i < n_bins; ++i) {
      compute_bin(i);
    }
    return;
  }

  std::vector<std::thread> workers;
  workers.reserve(worker_count);

  const size_t base_chunk = n_bins / worker_count;
  const size_t extra = n_bins % worker_count;
  size_t begin = 0;

  for (size_t w = 0; w < worker_count; ++w) {
    const size_t chunk_size = base_chunk + (w < extra ? 1 : 0);
    const size_t end = begin + chunk_size;
    workers.emplace_back([begin, end, &compute_bin]() {
      for (size_t i = begin; i < end; ++i) {
        compute_bin(i);
      }
    });
    begin = end;
  }

  for (auto& worker : workers) {
    worker.join();
  }
}

}  // namespace core
}  // namespace gryphon
