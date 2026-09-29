#ifndef GRYPHON_INJECTION_MSP_H
#define GRYPHON_INJECTION_MSP_H

#include "gryphon/core/input.h"
#include "gryphon/utils/random.h"

namespace gryphon {
namespace injection {

struct MSPSpin {
  double period;       // s (present-day, not birth period)
  double field;        // gauss
  double periodDot;    // s/s
  double luminosity;   // erg/s
};

// A RATE spectrum. Deliberately not a burst InjectionSpectrum: passing a
// particles/(erg s) rate to CosmicRays::run would silently lose a time integral.
class MSPInjection {
 public:
  explicit MSPInjection(const core::Input& input);
  MSPSpin sample(RandomNumberGenerator& rng) const;
  double rate(double energy, double spinDownPower) const;
  double normalization() const { return m_norm; }

 private:
  core::Input m_input;
  double m_norm;
};

}  // namespace injection
}  // namespace gryphon
#endif
