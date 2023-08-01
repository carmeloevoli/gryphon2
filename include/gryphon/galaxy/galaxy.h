#ifndef GRYPHON_GALAXY_H
#define GRYPHON_GALAXY_H

#include "gryphon/core/event.h"
#include "gryphon/core/input.h"
#include "gryphon/utils/random.h"
#include "gryphon/utils/vector3.h"

namespace gryphon {
namespace galaxy {

class Galaxy {
 public:
  Galaxy(const core::Input& input);
  virtual ~Galaxy() = default;

  void generate(RandomNumberGenerator& rng, bool show_bar = true);

  const core::Events& get_events() const { return m_events; }

 protected:
  virtual utils::Vector3d get_position(RandomNumberGenerator& rng) const = 0;

 protected:
  double m_rate;
  double m_dt;
  double m_tObs;
  double m_radius;
  utils::Vector3d m_sun;
  core::Events m_events;
  const utils::Vector3d m_GC{0., 0., 0.};
};

}  // namespace galaxy
}  // namespace gryphon

#endif  //  GRYPHON_GALAXY_H