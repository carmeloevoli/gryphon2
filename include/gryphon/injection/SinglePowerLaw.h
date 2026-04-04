#ifndef GRYPHON_PARTICLE_SINGLEPOWERLAWSPECTRUM_H
#define GRYPHON_PARTICLE_SINGLEPOWERLAWSPECTRUM_H

#include "gryphon/core/input.h"
#include "gryphon/injection/InjectionSpectrum.h"
#include "gryphon/utils/random.h"

namespace gryphon {
namespace injection {

class SinglePowerLawSpectrum final : public InjectionSpectrum {
 public:
  explicit SinglePowerLawSpectrum(const core::Input& in);
  virtual ~SinglePowerLawSpectrum() = default;

  double get(double E) const override;

  const double& alpha = m_alpha;
  const double& crEnergy = m_crenergy;

 private:
  double source_normalization() const;

 private:
  double m_E0;
  double m_alpha;
  double m_crenergy;
  double m_Emax;
  double m_Q0;
};

}  // namespace injection
}  // namespace gryphon

#endif  // GRYPHON_PARTICLE_SINGLEPOWERLAWSPECTRUM_H