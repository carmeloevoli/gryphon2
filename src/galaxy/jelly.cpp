#include "gryphon/galaxy/jelly.h"

#include "gryphon/utils/logging.h"

namespace gryphon {
namespace galaxy {

utils::Vector3d GalaxyJelly::get_position(RandomNumberGenerator& rng) const {
  double distance = 0;
  utils::Vector3d pos;
  while (distance < cgs::pc) {
    const auto z = rng.normal(0., m_h);
    const auto r = m_profile.pick_radius(rng());
    const auto angle = rng.uniform(0., 2. * M_PI);
    const auto x = r * std::cos(angle);
    const auto y = r * std::sin(angle);
    pos.setXYZ(x, y, z);
    distance = pos.getDistanceTo(m_sun);
    assert((pos - m_GC).getR() <= m_radius);
  }
  return pos;
}

}  // namespace galaxy
}  // namespace gryphon