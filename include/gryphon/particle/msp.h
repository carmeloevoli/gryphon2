#ifndef GRYPHON_PARTICLE_MSP_H
#define GRYPHON_PARTICLE_MSP_H

#include "gryphon/core/event.h"
#include "gryphon/core/input.h"
#include "gryphon/particle/particle.h"
#include "gryphon/utils/random.h"

namespace gryphon {
namespace particle {

class MSP final : public Particle {
 public:
  MSP(const core::Input& in, const std::shared_ptr<core::Event>& event, RandomNumberGenerator& rng);
  MSP(const core::Input& in);
  virtual ~MSP() = default;

  double Q(double E) const override;
  double b(double E) const;
  double tau(double E, double E_s) const;
  double Estar(double E, double dt) const;
  double lambda2(double E, double E_s) const;
  double get(double E, double dt, utils::Vector3d pos) const;
  double get(double E) const override;

  inline double escapeTimescale(double H) const { return pow2(H) / 2. / m_D0; }
  inline double lossesTimescale() const { return m_E0 / m_b0; }

 protected:
  const double m_E0 = 10. * cgs::GeV;
  const double m_D0 = 0.2 * cgs::kpc2 / cgs::Myr;
  const double m_delta = 0.56;
  const double m_b0 = 0.3 * cgs::GeV / cgs::Myr;
  const double m_alpha = 1.2;
  const double m_Ec = cgs::TeV;
  const double m_Luminosity = 1.5e32 * cgs::erg / cgs::sec;
};

}  // namespace particle
}  // namespace gryphon

#endif