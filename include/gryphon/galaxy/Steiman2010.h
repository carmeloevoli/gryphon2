#ifndef GRYPHON_GALAXY_STEIMAN2010_H
#define GRYPHON_GALAXY_STEIMAN2010_H

#include <array>

#include "gryphon/core/sourceprofile.h"
#include "gryphon/galaxy/galaxy.h"

namespace gryphon {
namespace galaxy {

struct ArmParams {
  double theta_0;
  double alpha;
  double a;
};

class GalaxySteiman2010 final : public Galaxy {
 public:
  GalaxySteiman2010(const core::Input& input) : Galaxy(input), m_profile(input) { init_spirals(); }

 protected:
  core::SourceProfile m_profile;
  std::array<ArmParams, 4> m_arms;

 protected:
  utils::Vector3d get_position(RandomNumberGenerator& rng) const override;
  void init_spirals();
};

}  // namespace galaxy
}  // namespace gryphon

#endif