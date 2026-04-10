#ifndef GRYPHON_PARTICLE_RANDOMEMAX_H
#define GRYPHON_PARTICLE_RANDOMEMAX_H

#include "gryphon/core/input.h"
#include "gryphon/injection/InjectionSpectrum.h"
#include "gryphon/utils/random.h"

namespace gryphon {
namespace injection {

class RandomEmaxSpectrum final : public InjectionSpectrum {
 public:
  RandomEmaxSpectrum(const core::Input& in, RandomNumberGenerator& rng);
  virtual ~RandomEmaxSpectrum() = default;

  double get(double E) const override;

  const double& alpha = m_alpha;
  const double& crEnergy = m_crenergy;
  const double& Emax = m_Emax;
  const double& Q0 = m_Q0;

 private:
  double source_normalization() const;

 private:
  double m_E0;
  double m_Emin;
  double m_alpha;
  double m_crenergy;
  double m_Emax;
  double m_Q0;
};

}  // namespace injection
}  // namespace gryphon

#endif  // GRYPHON_PARTICLE_RANDOMEMAX_H
