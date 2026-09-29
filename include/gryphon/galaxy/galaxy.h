#ifndef GRYPHON_GALAXY_H
#define GRYPHON_GALAXY_H

#include "gryphon/core/event.h"
#include "gryphon/core/input.h"
#include "gryphon/galaxy/associations.h"
#include "gryphon/utils/random.h"
#include "gryphon/utils/vector3.h"

namespace gryphon {
namespace galaxy {

class Galaxy {
 public:
  Galaxy(const core::Input& input);
  virtual ~Galaxy() = default;

  inline size_t size() const { return m_events.size(); }

  void generate(RandomNumberGenerator& rng, bool show_bar = true);
  void generate(RandomNumberGenerator& rng, bool show_bar,
                std::vector<AssociationEventOrigin>* origins);

  const core::Events& get_events() const { return m_events; }
  // Sample the configured spatial law without generating SN births or assigning
  // an age. Returned coordinates are heliocentric, like Event::pos.
  utils::Vector3d sample_position(RandomNumberGenerator& rng) const {
    return get_position(rng) - m_sun;
  }
  const AssociationGenerationStats& association_stats() const { return m_associationStats; }

 protected:
  virtual utils::Vector3d get_position(RandomNumberGenerator& rng) const = 0;

 protected:
  double m_rate;
  double m_dt;
  double m_tObs;
  double m_radius;
  double m_h;
  utils::Vector3d m_sun;
  core::Events m_events;
  core::Input m_input;
  AssociationGenerationStats m_associationStats;
  const utils::Vector3d m_GC{0., 0., 0.};
};

}  // namespace galaxy
}  // namespace gryphon

#endif  //  GRYPHON_GALAXY_H
