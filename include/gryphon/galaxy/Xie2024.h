#ifndef GRYPHON_GALAXY_XIE2024_H
#define GRYPHON_GALAXY_XIE2024_H

#include <array>

#include "gryphon/galaxy/galaxy.h"

namespace gryphon {
namespace galaxy {

class GalaxyXie2024 final : public Galaxy {
 public:
  GalaxyXie2024(const core::Input& input);

 protected:
  utils::Vector3d get_position(RandomNumberGenerator& rng) const override;

 private:
  static constexpr size_t kRadiusGridSize = 4096;

  std::array<double, kRadiusGridSize> m_radiusGrid{};
  std::array<double, kRadiusGridSize> m_radiusCdf{};
  double m_gcCdf = 0.;

  void init_radial_distribution();
  double radius_scale() const;
  double radial_density(double radius) const;
  double cdf_at_radius(double radius) const;
  double pick_radius(RandomNumberGenerator& rng, double cdf_min, double cdf_max) const;
};

}  // namespace galaxy
}  // namespace gryphon

#endif  // GRYPHON_GALAXY_XIE2024_H
