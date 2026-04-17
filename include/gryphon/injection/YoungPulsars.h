#ifndef GRYPHON_PARTICLE_YOUNGPULSARS_H
#define GRYPHON_PARTICLE_YOUNGPULSARS_H

#include "gryphon/core/input.h"
#include "gryphon/injection/InjectionSpectrum.h"
#include "gryphon/utils/random.h"

namespace gryphon {
namespace injection {

class YoungPulsarsSpectrum final : public InjectionSpectrum {
 public:
  YoungPulsarsSpectrum(const core::Input& in, RandomNumberGenerator& rng);
  virtual ~YoungPulsarsSpectrum() = default;

  double get(double E) const override;

  const double& initialPeriod = m_initialPeriod;
  const double& surfaceMagneticField = m_surfaceMagneticField;
  const double& chargeNumber = m_chargeNumber;
  const double& conversionEfficiency = m_conversionEfficiency;
  const double& rotEnergy = m_rotEnergy;
  const double& crEnergy = m_crenergy;
  const double& Emax = m_Emax;
  const double& Ecut = m_Emax;
  const double& tau0 = m_tau0;
  const double& Q0 = m_Q0;

 private:
  double maximumPotentialDropEnergy() const;

 private:
  double m_initialPeriod;
  double m_surfaceMagneticField;
  double m_chargeNumber;
  double m_conversionEfficiency;
  double m_rotEnergy;
  double m_crenergy;
  double m_Emax;
  double m_tau0;
  double m_Q0;
};

}  // namespace injection
}  // namespace gryphon

#endif  // GRYPHON_PARTICLE_YOUNGPULSARS_H
