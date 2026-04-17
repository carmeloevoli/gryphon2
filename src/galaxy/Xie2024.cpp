#include "gryphon/galaxy/Xie2024.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace galaxy {

namespace {

struct Xie2024ArmParams {
  double cumulative_weight;
  double theta_0;
  double tan_pitch;
  double r_0;
};

constexpr size_t kMaxSamplingTries = 1000000;
constexpr double kXieRsun = 8.3 * cgs::kpc;
constexpr double kRadialA = 9.03;
constexpr double kRadialB = 13.99;
constexpr double kRadialRpdf = 3.76 * cgs::kpc;
constexpr double kGcRadius = 3.57 * cgs::kpc;
constexpr double kArmRadialScatter = 0.07;
constexpr double kThetaBroadeningScale = 0.35;
constexpr double kArmCumulativeWeight = 84.92 / 100.;

const std::array<Xie2024ArmParams, 5> kArms = {
    // Xie et al. 2024 Table 4; coordinates follow Hou & Han 2014 Eq. 18.
    Xie2024ArmParams{19.77 / 100., utils::deg2rad(44.4),
                     std::tan(utils::deg2rad(11.43)), 3.35 * cgs::kpc},
    Xie2024ArmParams{(19.77 + 20.91) / 100., utils::deg2rad(120.0),
                     std::tan(utils::deg2rad(9.84)), 3.71 * cgs::kpc},
    Xie2024ArmParams{(19.77 + 20.91 + 21.20) / 100., utils::deg2rad(218.6),
                     std::tan(utils::deg2rad(10.38)), 3.56 * cgs::kpc},
    Xie2024ArmParams{(19.77 + 20.91 + 21.20 + 20.17) / 100., utils::deg2rad(330.3),
                     std::tan(utils::deg2rad(10.54)), 3.67 * cgs::kpc},
    Xie2024ArmParams{kArmCumulativeWeight, utils::deg2rad(55.1),
                     std::tan(utils::deg2rad(2.77)), 8.21 * cgs::kpc},
};

size_t pickArmIndex(double rnd) {
  for (size_t i = 0; i < kArms.size(); ++i) {
    if (rnd < kArms[i].cumulative_weight) return i;
  }
  return kArms.size() - 1;
}

double getArmTheta(double radius, const Xie2024ArmParams& arm, double radius_scale) {
  return std::log(radius / (arm.r_0 * radius_scale)) / arm.tan_pitch + arm.theta_0;
}

}  // namespace

GalaxyXie2024::GalaxyXie2024(const core::Input& input) : Galaxy(input) {
  init_radial_distribution();
}

double GalaxyXie2024::radius_scale() const {
  return m_sun.getX() > 0. ? m_sun.getX() / kXieRsun : 1.;
}

double GalaxyXie2024::radial_density(double radius) const {
  const auto model_radius = radius / radius_scale();
  const auto normalized_radius = (model_radius + kRadialRpdf) / (kXieRsun + kRadialRpdf);
  return std::pow(normalized_radius, kRadialA) *
         std::exp(-kRadialB * (model_radius - kXieRsun) / (kXieRsun + kRadialRpdf));
}

void GalaxyXie2024::init_radial_distribution() {
  const auto dr = m_radius / static_cast<double>(kRadiusGridSize - 1);
  double integral = 0.;
  double previous_value = 0.;

  m_radiusGrid[0] = 0.;
  m_radiusCdf[0] = 0.;
  for (size_t i = 1; i < kRadiusGridSize; ++i) {
    const auto radius = static_cast<double>(i) * dr;
    const auto value = radius * radial_density(radius);
    integral += 0.5 * (previous_value + value) * dr;
    m_radiusGrid[i] = radius;
    m_radiusCdf[i] = integral;
    previous_value = value;
  }

  if (!(integral > 0.)) {
    throw std::runtime_error("GalaxyXie2024 radial distribution has zero normalization");
  }
  for (auto& cdf : m_radiusCdf) {
    cdf /= integral;
  }

  m_gcCdf = cdf_at_radius(std::min(kGcRadius * radius_scale(), m_radius));
}

double GalaxyXie2024::cdf_at_radius(double radius) const {
  if (radius <= m_radiusGrid.front()) return m_radiusCdf.front();
  if (radius >= m_radiusGrid.back()) return m_radiusCdf.back();

  const auto it = std::lower_bound(m_radiusGrid.begin(), m_radiusGrid.end(), radius);
  const auto idx = static_cast<size_t>(std::distance(m_radiusGrid.begin(), it));
  const auto r0 = m_radiusGrid[idx - 1];
  const auto r1 = m_radiusGrid[idx];
  const auto c0 = m_radiusCdf[idx - 1];
  const auto c1 = m_radiusCdf[idx];
  return c0 + (c1 - c0) * (radius - r0) / (r1 - r0);
}

double GalaxyXie2024::pick_radius(RandomNumberGenerator& rng, double cdf_min,
                                  double cdf_max) const {
  if (!(cdf_max > cdf_min)) return m_radiusGrid.front();

  const auto target = cdf_min + rng() * (cdf_max - cdf_min);
  const auto it = std::lower_bound(m_radiusCdf.begin(), m_radiusCdf.end(), target);
  if (it == m_radiusCdf.begin()) return m_radiusGrid.front();
  if (it == m_radiusCdf.end()) return m_radiusGrid.back();

  const auto idx = static_cast<size_t>(std::distance(m_radiusCdf.begin(), it));
  const auto c0 = m_radiusCdf[idx - 1];
  const auto c1 = m_radiusCdf[idx];
  const auto r0 = m_radiusGrid[idx - 1];
  const auto r1 = m_radiusGrid[idx];
  return r0 + (r1 - r0) * (target - c0) / (c1 - c0);
}

utils::Vector3d GalaxyXie2024::get_position(RandomNumberGenerator& rng) const {
  utils::Vector3d pos;
  const auto scale = radius_scale();
  for (size_t tries = 0; tries < kMaxSamplingTries; ++tries) {
    const auto z = rng.normal(0., m_h);
    double radius = 0.;
    double theta = 0.;
    const auto component = rng();

    if (component >= kArmCumulativeWeight || m_gcCdf >= 1.) {
      radius = pick_radius(rng, 0., m_gcCdf);
      theta = rng.uniform(0., 2. * M_PI);
    } else {
      const auto raw_radius = pick_radius(rng, m_gcCdf, 1.);
      if (!(raw_radius > 0.)) continue;

      const auto& arm = kArms[pickArmIndex(component)];
      const auto theta_raw = getArmTheta(raw_radius, arm, scale);
      theta = theta_raw +
              rng.uniform(0., 2. * M_PI) *
                  std::exp(-kThetaBroadeningScale * raw_radius / cgs::kpc);
      radius = rng.normal(raw_radius, kArmRadialScatter * raw_radius);
      if (!(radius > 0.) || radius > m_radius) continue;
    }

    pos.setXYZ(radius * std::cos(theta), radius * std::sin(theta), z);

    const auto distance = pos.getDistanceTo(m_sun);
    if (pos.getDistanceTo(m_GC) <= m_radius && distance >= cgs::pc) return pos;
  }

  throw std::runtime_error("GalaxyXie2024::get_position failed to sample a valid event position");
}

}  // namespace galaxy
}  // namespace gryphon
