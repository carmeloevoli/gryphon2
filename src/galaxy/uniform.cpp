#include "gryphon/galaxy/uniform.h"

#include <stdexcept>

#include "gryphon/utils/logging.h"

namespace gryphon {
namespace galaxy {

namespace {
constexpr size_t kMaxSamplingTries = 1000000;
}

utils::Vector3d GalaxyUniform::get_position(RandomNumberGenerator& rng) const {
  const double z = 0;
  utils::Vector3d pos;
  for (size_t tries = 0; tries < kMaxSamplingTries; ++tries) {
    const auto r = std::sqrt(rng()) * m_radius;
    const auto angle = rng.uniform(0., 2. * M_PI);
    const auto x = r * std::cos(angle);
    const auto y = r * std::sin(angle);
    pos.setXYZ(x, y, z);
    const auto distance = pos.getDistanceTo(m_sun);
    assert(pos.getDistanceTo(m_GC) <= m_radius);
    if (distance >= cgs::pc) return pos;
  }

  throw std::runtime_error("GalaxyUniform::get_position failed to sample a valid event position");
}

}  // namespace galaxy
}  // namespace gryphon