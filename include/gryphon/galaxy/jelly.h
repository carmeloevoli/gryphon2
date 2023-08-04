#ifndef GALAXY_JELLY_H
#define GALAXY_JELLY_H

#include "gryphon/core/sourceprofile.h"
#include "gryphon/galaxy/galaxy.h"

namespace gryphon {
namespace galaxy {

class GalaxyJelly final : public Galaxy {
 public:
  GalaxyJelly(const core::Input& input) : Galaxy(input), m_profile(input) {}

 protected:
  core::SourceProfile m_profile;

 protected:
  utils::Vector3d get_position(RandomNumberGenerator& rng) const override;
};

}  // namespace galaxy
}  // namespace gryphon

#endif