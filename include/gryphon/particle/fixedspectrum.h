#ifndef GRYPHON_PARTICLE_FIXEDSPECTRUM_H
#define GRYPHON_PARTICLE_FIXEDSPECTRUM_H

#include "gryphon/core/diffusion.h"
#include "gryphon/core/input.h"
#include "gryphon/particle/particle.h"
#include "gryphon/utils/random.h"

namespace gryphon {
namespace particle {

class FixedSpectrumParticle final : public Particle {
 public:
  FixedSpectrumParticle(const core::Input& in, const std::shared_ptr<core::Event>& event,
                        RandomNumberGenerator& rng);
  FixedSpectrumParticle(const core::Input& in);
  virtual ~FixedSpectrumParticle() = default;

  double Q(double E) const override;
  double get(double E, double dt, utils::Vector3d pos) const;
  double get(double E) const override;

  const double& alpha = m_alpha;
  const double& crEnergy = m_crenergy;

 protected:
  void init(const core::Input& in);
  double source_normalization() const;

 protected:
  const double m_E0 = 1. * cgs::GeV;
  double m_alpha;
  double m_crenergy;
  double m_Emax;
  double m_H;
  double m_Q0;
  core::DiffusionCoefficient m_D;
};

}  // namespace particle
}  // namespace gryphon

#endif  // GRYPHON_PARTICLE_FIXEDSPECTRUM_H