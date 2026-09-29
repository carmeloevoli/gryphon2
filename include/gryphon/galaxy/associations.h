#ifndef GRYPHON_GALAXY_ASSOCIATIONS_H
#define GRYPHON_GALAXY_ASSOCIATIONS_H

#include <cstddef>
#include <functional>
#include <vector>

#include "gryphon/core/event.h"
#include "gryphon/core/input.h"
#include "gryphon/utils/random.h"

namespace gryphon::galaxy {

// Normalized single-star DTD, Zapartas et al. (2017), Appendix A, Eq. (10).
// Arguments/results are in Myr, not CGS. The IMF is already integrated out.
double associationDelayCdf(double delayMyr);
double associationDelayQuantile(double probability);

struct AssociationGenerationStats {
  size_t associations = 0;        // Parents drawn, including those with no retained SN.
  size_t fieldExplosions = 0;     // Before the observer's light-cone selection.
  size_t clusteredExplosions = 0;
  size_t retainedExplosions = 0;
};

// Optional diagnostic information, in the same order as retained events.
// All times/positions are CGS; centres are Galactocentric. Parent 0 denotes
// independent field events, positive IDs identify shared association births.
struct AssociationEventOrigin {
  size_t parentId;
  double birthAge;
  double delay;
  utils::Vector3d centre;
};

// positionSampler returns Galactocentric centres from the existing spatial law.
core::Events generateAssociations(
    const core::Input& input, RandomNumberGenerator& rng,
    const std::function<utils::Vector3d(RandomNumberGenerator&)>& positionSampler,
    AssociationGenerationStats& stats);
core::Events generateAssociations(
    const core::Input& input, RandomNumberGenerator& rng,
    const std::function<utils::Vector3d(RandomNumberGenerator&)>& positionSampler,
    AssociationGenerationStats& stats, std::vector<AssociationEventOrigin>* origins);

}  // namespace gryphon::galaxy

#endif
