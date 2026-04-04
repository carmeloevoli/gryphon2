#include "gryphon/galaxy/jelly.h"

#include <stdexcept>

#include "gryphon/utils/logging.h"

namespace gryphon {
namespace galaxy {

namespace {
constexpr size_t kMaxSamplingTries = 1000000;
}

utils::Vector3d GalaxyJelly::get_position(RandomNumberGenerator& rng) const {
  utils::Vector3d pos;
  for (size_t tries = 0; tries < kMaxSamplingTries; ++tries) {
    const auto z = rng.normal(0., m_h);
    const auto r = m_profile.pick_radius(rng());
    const auto angle = rng.uniform(0., 2. * M_PI);
    const auto x = r * std::cos(angle);
    const auto y = r * std::sin(angle);
    pos.setXYZ(x, y, z);
    const auto distance = pos.getDistanceTo(m_sun);
    assert((pos - m_GC).getR() <= m_radius);
    if (distance >= cgs::pc) return pos;
  }

  throw std::runtime_error("GalaxyJelly::get_position failed to sample a valid event position");
}

}  // namespace galaxy
}  // namespace gryphon