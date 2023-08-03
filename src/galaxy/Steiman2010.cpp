#include "gryphon/galaxy/Steiman2010.h"

#include "gryphon/utils/logging.h"
#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace galaxy {

void GalaxySteiman2010::init_spirals() {
  ArmParams Norma = {utils::deg2rad(13.6), 0.242, 0.246 * cgs::kpc};
  m_arms.insert(std::pair<size_t, ArmParams>(1, Norma));

  ArmParams Carina = {utils::deg2rad(15.6), 0.279, 0.608 * cgs::kpc};
  m_arms.insert(std::pair<size_t, ArmParams>(2, Carina));

  ArmParams Perseus = {utils::deg2rad(13.5), 0.249, 0.449 * cgs::kpc};
  m_arms.insert(std::pair<size_t, ArmParams>(3, Perseus));

  ArmParams Crux = {utils::deg2rad(13.5), 0.240, 0.378 * cgs::kpc};
  m_arms.insert(std::pair<size_t, ArmParams>(4, Crux));
}

int pick_a_spiral_arm(double rnd) {
  double B_total = 169 + 266 + 339 + 176;
  if (rnd < 169 / B_total)
    return 1;
  else if (rnd < (169 + 266) / B_total)
    return 2;
  else if (rnd < (169 + 266 + 339) / B_total)
    return 3;
  else
    return 4;
}

double get_angle(const double& r, const ArmParams& arm) {
  return std::log(r / arm.a) / arm.alpha + arm.theta_0;
}

utils::Vector3d GalaxySteiman2010::get_position(RandomNumberGenerator& rng) const {
  double distance = 0;
  utils::Vector3d pos;
  while (distance < cgs::pc) {
    const auto z = rng.normal(0., m_h);
    const auto r_0 = m_profile.pick_radius(rng());
    if (r_0 > 3. * cgs::kpc) {
      const auto i = pick_a_spiral_arm(rng());
      const auto arm = m_arms.at(i);
      auto phi = get_angle(r_0, arm);
      phi += utils::deg2rad(rng.normal(0., 15.));
      pos.setXYZ(r_0 * std::cos(-M_PI / 2. + phi), r_0 * std::sin(-M_PI / 2. + phi), z);
    } else {
      const auto bar_lenght = 3. * cgs::kpc;
      const auto bar_angle = utils::deg2rad(-25);
      const auto d = rng.uniform(-1, 1) * bar_lenght;
      const auto blur_x = rng.normal(d * std::cos(bar_angle), 0.5 * cgs::kpc);
      const auto blur_y = rng.normal(d * std::sin(bar_angle), 0.5 * cgs::kpc);
      pos.setXYZ(blur_x, blur_y, z);
    }
    distance = pos.getDistanceTo(m_sun);
    assert(pos.getDistanceTo(m_GC) <= m_radius);
  }
  return pos;
}

}  // namespace galaxy
}  // namespace gryphon