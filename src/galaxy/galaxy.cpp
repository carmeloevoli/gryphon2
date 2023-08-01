#include "gryphon/galaxy/galaxy.h"

#include "gryphon/utils/logging.h"
#include "gryphon/utils/progressbar.h"
#include "gryphon/utils/timer.h"

namespace gryphon {
namespace galaxy {

Galaxy::Galaxy(const core::Input& input)
    : m_rate(input.sn_rate), m_dt(input.time_step), m_tObs(input.max_time), m_radius(input.R_g) {
  m_sun.setX(input.R_sun);
}

void Galaxy::generate(RandomNumberGenerator& rng, bool show_bar) {
  if (m_dt * m_rate > 1) throw std::runtime_error("time step > 1 / rate");
  // utils::Timer timer("generate Galaxy");
  auto bar = std::make_shared<utils::ProgressBar>(m_tObs / m_dt);
  bar->start("generate Galaxy");
  double t = 0;
  while (t < m_tObs) {
    const auto r = rng();
    if (r < m_dt * m_rate) {
      const auto pos = get_position(rng);
      if (pos.getDistanceTo(m_sun) < cgs::c_light * (m_tObs - t)) {
        auto sn = std::make_shared<core::Event>(m_tObs - t, pos - m_sun);  // m_rng, m_input);
        m_events.emplace_back(sn);
      }
    }
    t += m_dt;
    if (show_bar) bar->update();
  }
  LOGD << "expected events : " << m_tObs * m_rate;
  LOGD << "fraction of simulated events : " << m_events.size() / (m_tObs * m_rate);
  LOGD << "memory occupied : " << m_events.capacity() * sizeof core::Event() / 1e9 << " Gb";
}

}  // namespace galaxy
}  // namespace gryphon