#include "gryphon/galaxy/associations.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace gryphon::galaxy {
namespace {

constexpr double kMinDelay = 3.;
constexpr double kBreakDelay = 25.;
constexpr double kMaxDelay = 48.;

// Integral of [a + b log10(t) + c log10(t)^2] dt/t. Common ln(10)
// and dimensional factors cancel in the normalized CDF. The second segment
// has ten times the prefactor of the first in the published fit.
double primitive(double t, bool late) {
  const double x = std::log10(t);
  return late ? 10. * (-4.85 * x + 6.55 * x * x / 2. - 1.92 * x * x * x / 3.)
              : -2.83 * x + 8.70 * x * x / 2. - 2.07 * x * x * x / 3.;
}

double cumulative(double t) {
  const double early = primitive(std::min(t, kBreakDelay), false) -
                       primitive(kMinDelay, false);
  return early + (t > kBreakDelay ? primitive(t, true) - primitive(kBreakDelay, true) : 0.);
}

double drawDelay(RandomNumberGenerator& rng) {
  // Uniform-quantile interpolation: <= 1/8192 absolute CDF discretization.
  static const auto quantiles = [] {
    std::array<double, 8193> values{};
    for (size_t i = 0; i < values.size(); ++i)
      values[i] = associationDelayQuantile(double(i) / (values.size() - 1));
    return values;
  }();
  const double index = rng() * (quantiles.size() - 1);
  const auto i = std::min(static_cast<size_t>(index), quantiles.size() - 2);
  return (quantiles[i] + (index - i) * (quantiles[i + 1] - quantiles[i])) * cgs::Myr;
}

double waitTime(RandomNumberGenerator& rng, double rate) {
  // 1-u is in (0,1]; use a positive open variate to avoid a zero increment.
  double u;
  do { u = rng(); } while (u == 0.);
  return -std::log(u) / rate;
}

}  // namespace

double associationDelayCdf(double delayMyr) {
  if (!std::isfinite(delayMyr)) throw std::invalid_argument("delay must be finite");
  if (delayMyr <= kMinDelay) return 0.;
  if (delayMyr >= kMaxDelay) return 1.;
  return cumulative(delayMyr) / cumulative(kMaxDelay);
}

double associationDelayQuantile(double probability) {
  if (!(probability >= 0. && probability <= 1.))
    throw std::invalid_argument("delay probability must be in [0,1]");
  if (probability == 0.) return kMinDelay;
  if (probability == 1.) return kMaxDelay;
  double lo = kMinDelay, hi = kMaxDelay;
  for (int i = 0; i < 48; ++i) {
    const double mid = 0.5 * (lo + hi);
    if (associationDelayCdf(mid) < probability) lo = mid;
    else hi = mid;
  }
  return 0.5 * (lo + hi);
}

core::Events generateAssociations(
    const core::Input& input, RandomNumberGenerator& rng,
    const std::function<utils::Vector3d(RandomNumberGenerator&)>& positionSampler,
    AssociationGenerationStats& stats) {
  return generateAssociations(input, rng, positionSampler, stats, nullptr);
}

core::Events generateAssociations(
    const core::Input& input, RandomNumberGenerator& rng,
    const std::function<utils::Vector3d(RandomNumberGenerator&)>& positionSampler,
    AssociationGenerationStats& stats, std::vector<AssociationEventOrigin>* origins) {
  input.validate();
  stats = {};
  core::Events events;
  events.reserve(static_cast<size_t>(1.1 * input.sn_rate() * input.max_time() + 1.));
  if (origins) {
    origins->clear();
    origins->reserve(events.capacity());
  }
  const utils::Vector3d sun(input.R_sun(), 0., 0.);

  const auto addEvent = [&](double age, double delay, const utils::Vector3d& centre,
                            size_t parentId, double birthAge) {
    // Ballistic relative spreading, not orbital evolution. Both components
    // use the same offsets, so f=0 and f=1 have the same one-point density.
    const double sigma = std::hypot(input.associationRadius(),
                                    input.associationVelocity() * delay);
    utils::Vector3d position = centre;
    if (sigma > 0.) {
      bool accepted = false;
      for (size_t trial = 0; trial < 100000; ++trial) {
        position = centre + utils::Vector3d(rng.normal(0., sigma), rng.normal(0., sigma),
                                           rng.normal(0., sigma));
        if (position.getModule() <= input.R_g() && position.getDistanceTo(sun) >= cgs::pc) {
          accepted = true;
          break;
        }
      }
      if (!accepted) throw std::runtime_error("could not sample association displacement");
    }
    const auto relative = position - sun;
    if (relative.getModule() < cgs::c_light * age) {
      events.emplace_back(std::make_shared<core::Event>(age, relative));
      if (origins) origins->push_back({parentId, birthAge, delay, centre});
    }
  };

  const double fieldRate = (1. - input.associationFraction()) * input.sn_rate();
  if (fieldRate > 0.) {
    for (double age = waitTime(rng, fieldRate); age < input.max_time();
         age += waitTime(rng, fieldRate)) {
      ++stats.fieldExplosions;
      const double delay = drawDelay(rng);
      addEvent(age, delay, positionSampler(rng), 0, age + delay);
    }
  }

  const double birthRate = input.associationFraction() * input.sn_rate() /
                           input.associationMembers();
  if (birthRate > 0.) {
    // Include old parents! For each delay, birth ages [delay,T+delay] have
    // length T, giving E[N_SN] = R*T without an initial-transient deficit.
    const double oldestBirth = input.max_time() + kMaxDelay * cgs::Myr;
    for (double birthAge = kMinDelay * cgs::Myr + waitTime(rng, birthRate);
         birthAge < oldestBirth; birthAge += waitTime(rng, birthRate)) {
      ++stats.associations;
      const auto centre = positionSampler(rng);
      for (ulong member = 0; member < input.associationMembers(); ++member) {
        const double delay = drawDelay(rng);
        const double age = birthAge - delay;
        if (!(age > 0. && age < input.max_time())) continue;
        ++stats.clusteredExplosions;
        addEvent(age, delay, centre, stats.associations, birthAge);
      }
    }
  }
  stats.retainedExplosions = events.size();
  return events;
}

}  // namespace gryphon::galaxy
