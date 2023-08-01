#include "gryphon/galaxy/uniform.h"

#include "gryphon/utils/logging.h"

namespace gryphon {
namespace galaxy {

utils::Vector3d GalaxyUniform::get_position(RandomNumberGenerator& rng) const {
  double distance = 0;
  const double z = 0;
  utils::Vector3d pos;
  while (distance < cgs::pc) {
    const auto r = std::sqrt(rng()) * m_radius;
    const auto angle = rng.uniform(0., 2. * M_PI);
    const auto x = r * std::cos(angle);
    const auto y = r * std::sin(angle);
    pos.setXYZ(x, y, z);
    distance = pos.getDistanceTo(m_sun);
    assert(pos.getDistanceTo(m_GC) <= m_radius);
  }
  return pos;
}

}  // namespace galaxy
}  // namespace gryphon