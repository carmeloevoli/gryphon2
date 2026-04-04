#ifndef GRYPHON_PARTICLE_GALACTICRANDOM_H
#define GRYPHON_PARTICLE_GALACTICRANDOM_H

#include <memory>

#include "gryphon/core/input.h"
#include "gryphon/injection/InjectionSpectrum.h"
#include "gryphon/utils/random.h"

namespace gryphon {
namespace injection {

class GalacticRandomSpectrum final : public InjectionSpectrum {
 public:
  GalacticRandomSpectrum(const core::Input& in, RandomNumberGenerator& rng);
  virtual ~GalacticRandomSpectrum() = default;

  double get(double E) const override;

  const double& alpha = m_alpha;
  const double& crEnergy = m_crenergy;
  const double& Q0 = m_Q0;

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

#endif  // GRYPHON_PARTICLE_GALACTICRANDOM_H