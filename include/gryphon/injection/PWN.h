#ifndef GRYPHON_PARTICLE_PWNSPECTRUM_H
#define GRYPHON_PARTICLE_PWNSPECTRUM_H

#include "gryphon/core/input.h"
#include "gryphon/injection/InjectionSpectrum.h"
#include "gryphon/utils/random.h"
#include "gryphon/utils/vector3.h"

namespace gryphon {
namespace injection {

class PWNSpectrum final : public InjectionSpectrum {
 public:
  PWNSpectrum(const core::Input& in, RandomNumberGenerator& rng);
  virtual ~PWNSpectrum() = default;

  double get(double E) const override;

  const double& alphaBelowBreak = m_alphaBelowBreak;
  const double& alphaAboveBreak = m_alphaAboveBreak;
  const double& initialPeriod = m_initialPeriod;
  const double& conversionEfficiency = m_conversionEfficiency;
  const double& rotEnergy = m_rotEnergy;
  const double& spinDownAge = m_spinDownAge;
  const double& crEnergy = m_crenergy;
  const double& Ebreak = m_Ebreak;
  const double& Ecut = m_Ecut;
  const double& Emin = m_Emin;
  const double& Emax = m_Ecut;
  const utils::Vector3d& kickVelocity = m_kickVelocity;
  const double& kickSpeed = m_kickSpeed;
  const double& Q0 = m_Q0;

 private:
  double source_normalization() const;
  double spectralShape(double E) const;

 private:
  double m_Emin;
  double m_alphaBelowBreak;
  double m_alphaAboveBreak;
  double m_initialPeriod;
  double m_conversionEfficiency;
  double m_rotEnergy;
  double m_spinDownAge;
  double m_crenergy;
  double m_Ebreak;
  double m_Ecut;
  utils::Vector3d m_kickVelocity;
  double m_kickSpeed;
  double m_Q0;
};

}  // namespace injection
}  // namespace gryphon

#endif  // GRYPHON_PARTICLE_PWNSPECTRUM_H
