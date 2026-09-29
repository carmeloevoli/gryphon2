#ifndef GRYPHON_CORE_OBSERVABLES_H
#define GRYPHON_CORE_OBSERVABLES_H

#include "gryphon/core/cosmicrays.h"
#include "gryphon/core/run.h"

namespace gryphon {
namespace core {

// Write the spectrum and the complete realization-level dipole vector.  The
// vector components use the simulation frame: +x points from the Galactic
// centre through the Sun, +y lies in the Galactic plane, and +z points north.
void dumpFluxAndDipole(const RunOutput& run, const CosmicRays& cosmicRays);

}  // namespace core
}  // namespace gryphon

#endif  // GRYPHON_CORE_OBSERVABLES_H
