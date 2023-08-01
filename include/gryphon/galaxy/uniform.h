#ifndef GALAXY_UNIFORM_H
#define GALAXY_UNIFORM_H

#include "gryphon/galaxy/galaxy.h"

namespace gryphon {
namespace galaxy {

class GalaxyUniform : public Galaxy {
 public:
  GalaxyUniform(const core::Input& input) : Galaxy(input) {}

 protected:
  utils::Vector3d get_position(RandomNumberGenerator& rng) const override;
};

}  // namespace galaxy
}  // namespace gryphon

#endif