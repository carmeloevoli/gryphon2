#ifndef GRYPHON_INJECTION_SMOOTHBROKENPOWERLAW_H
#define GRYPHON_INJECTION_SMOOTHBROKENPOWERLAW_H

#include "gryphon/core/input.h"
#include "gryphon/injection/InjectionSpectrum.h"

namespace gryphon {
namespace injection {

// Q(E) proportional to E^-alpha below the break and
// E^-(alpha+deltaAlpha) above it, with a smooth transition.
class SmoothBrokenPowerLawSpectrum final : public InjectionSpectrum {
 public:
  explicit SmoothBrokenPowerLawSpectrum(const core::Input& input);

  double get(double E) const override;

 private:
  double shape(double E) const;
  double sourceNormalization() const;

  double m_E0;
  double m_Emin;
  double m_Emax;
  double m_Ebreak;
  double m_alpha;
  double m_deltaAlpha;
  double m_smoothness;
  double m_crEnergy;
  double m_Q0;
};

}  // namespace injection
}  // namespace gryphon

#endif  // GRYPHON_INJECTION_SMOOTHBROKENPOWERLAW_H
