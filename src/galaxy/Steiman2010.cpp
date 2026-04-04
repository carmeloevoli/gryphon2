#include "gryphon/galaxy/Steiman2010.h"

#include <array>
#include <stdexcept>

#include "gryphon/utils/logging.h"
#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace galaxy {

namespace {

constexpr size_t kMaxSamplingTries = 1000000;
const std::array<double, 4> kArmCumulativeWeights = {
    169. / (169. + 266. + 339. + 176.),
    (169. + 266.) / (169. + 266. + 339. + 176.),
    (169. + 266. + 339.) / (169. + 266. + 339. + 176.),
    1.0,
};
const double kSpiralRadiusThreshold = 3. * cgs::kpc;
const double kBarLength = 3. * cgs::kpc;
const double kBarAngle = utils::deg2rad(-25.);
const double kBarBlurSigma = 0.5 * cgs::kpc;
const double kSpiralPhaseSigmaDeg = 15.;

size_t pickSpiralArmIndex(double rnd) {
  for (size_t i = 0; i < kArmCumulativeWeights.size(); ++i) {
    if (rnd < kArmCumulativeWeights[i]) return i;
  }
  return kArmCumulativeWeights.size() - 1;
}

double get_angle(const double& r, const ArmParams& arm) {
  return std::log(r / arm.a) / arm.alpha + arm.theta_0;
}

}  // namespace

void GalaxySteiman2010::init_spirals() {
  m_arms = {ArmParams{utils::deg2rad(13.6), 0.242, 0.246 * cgs::kpc},
            ArmParams{utils::deg2rad(15.6), 0.279, 0.608 * cgs::kpc},
            ArmParams{utils::deg2rad(13.5), 0.249, 0.449 * cgs::kpc},
            ArmParams{utils::deg2rad(13.5), 0.240, 0.378 * cgs::kpc}};
}

utils::Vector3d GalaxySteiman2010::get_position(RandomNumberGenerator& rng) const {
  utils::Vector3d pos;
  for (size_t tries = 0; tries < kMaxSamplingTries; ++tries) {
    const auto z = rng.normal(0., m_h);
    const auto r_0 = m_profile.pick_radius(rng());
    if (r_0 > kSpiralRadiusThreshold) {
      const auto& arm = m_arms[pickSpiralArmIndex(rng())];
      auto phi = get_angle(r_0, arm);
      phi += utils::deg2rad(rng.normal(0., kSpiralPhaseSigmaDeg));
      pos.setXYZ(r_0 * std::cos(-M_PI / 2. + phi), r_0 * std::sin(-M_PI / 2. + phi), z);
    } else {
      const auto d = rng.uniform(-1., 1.) * kBarLength;
      const auto blur_x = rng.normal(d * std::cos(kBarAngle), kBarBlurSigma);
      const auto blur_y = rng.normal(d * std::sin(kBarAngle), kBarBlurSigma);
      pos.setXYZ(blur_x, blur_y, z);
    }

    const auto distance = pos.getDistanceTo(m_sun);
    assert(pos.getDistanceTo(m_GC) <= 1.01 * m_radius);
    if (distance >= cgs::pc) return pos;
  }

  throw std::runtime_error(
      "GalaxySteiman2010::get_position failed to sample a valid event position");
}

}  // namespace galaxy
}  // namespace gryphon