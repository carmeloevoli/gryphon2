#ifndef GRYPHON_CORE_INPUT_H
#define GRYPHON_CORE_INPUT_H

#include <string>

#include "gryphon/core/cgs.h"
#include "gryphon/core/pid.h"

using ulong = unsigned long int;

enum class SpiralModel {
  Uniform,      // Uniform in (x,y) and z = 0
  Jelly,        // Only profile in r and z is retained
  Steiman2010,  // Steiman-Cameron et al., ApJ, 722, 1460–1473, 2010
  Faucher2006,  // Faucher-Giguere and Kaspi, ApJ, 643, 332–355, 2006
  Vallee2008    // Vallee, AJ, 135, 1301-1310, 2008
};

enum class ParticleModel { SingleSpectrum, VaryingSlopeSpectrum, VaryingEmaxSpectrum };

namespace gryphon {
namespace core {

class Input {
 private:
  ulong _seed = 0;
  std::string _simname = "test";
  // output energy vector
  double _E_min = 10. * cgs::GeV;
  double _E_max = 10. * cgs::TeV;
  ulong _E_size = 3 * 16;
  // galaxy size
  double _H = 5. * cgs::kpc;
  double _Rg = 20. * cgs::kpc;
  double _Rsun = 8.5 * cgs::kpc;
  // diffusion coefficient
  double _D0_over_H = 0.44e28 * cgs::cm2 / cgs::sec / cgs::kpc;
  double _E_0 = cgs::GeV;
  double _delta = 0.5;
  double _ddelta = 0.2;
  double _s = 0.1;
  double _E_b = 312. * cgs::GeV;
  // source profile
  double _a = 0;
  double _b = 0;
  double _R1 = 0;
  //   // SNR spectrum
  //   double _SNR_alpha = 2.25;
  //   double _SNR_E_c = gamma_max * cgs::electron_mass_c2;  // constants.h;
  //   double _SNR_efficiency = 1e-3;
  //   // Pulsar parameters
  //   double _pulsar_fraction = 0.8;
  //   double _pulsar_mean_period = 300. * cgs::msec;
  //   double _pulsar_sdev_period = 150. * cgs::msec;
  //   double _PWN_alpha_low = 1.4;
  //   double _PWN_alpha_high = 2.3;
  //   double _PWN_E_break = 500. * cgs::GeV;
  //   double _PWN_efficiency = 0.5;
  //   // energy losses
  //   double _b_0 = 4e-7 * cgs::eV / cgs::sec;  // only for tests
  //   double _f_UV = 1;
  //   double _B_halo = cgs::microgauss;
  // simulation parameters
  double _sn_rate = 1. / 30. / cgs::year;
  double _time_step = 1. * cgs::year;
  double _max_time = 1. * cgs::Myr;
  //   double _max_radius = 10. * cgs::kpc;
  // models
  //   RunMode _runMode;
  PID::PID _pid = PID::H;
  ParticleModel _particleModel = ParticleModel::SingleSpectrum;
  SpiralModel _spiralModel = SpiralModel::Uniform;

 protected:
  //   void read_params_file(const std::string& filename);
  //   void set_simname(const std::string& inifilename);
  //   void set_outputdir(const std::string& dir);
  //   void validator();

 public:
  Input();
  Input(const std::string& filename);
  virtual ~Input() = default;
  void print();
  //   void save_on_file(const std::shared_ptr<RNG> rng);
  //   void set_seed(const size_t& seed);
  //   void set_params(const std::string& key, const double& value);

  const std::string& simname = _simname;
  const unsigned long int& seed = _seed;
  const double& E_min = _E_min;
  const double& E_max = _E_max;
  const size_t& E_size = _E_size;
  const double& H = _H;
  const double& R_g = _Rg;
  const double& R_sun = _Rsun;
  const double& D0_over_H = _D0_over_H;
  const double& E_0 = _E_0;
  const double& delta = _delta;
  const double& ddelta = _ddelta;
  const double& s = _s;
  const double& E_b = _E_b;
  const double& a = _a;
  const double& b = _b;
  const double& R_1 = _R1;

  //    const double& D_delta = _D_delta;
  //    const double& D_E_break = _D_E_break;
  //    const double& D_ddelta = _D_ddelta;
  //    const double& D_smoothness = _D_smoothness;
  //    const double& D0_over_H = _D0_over_H;
  //    const double& E_0 = _E_0;
  //    const double& SNR_alpha = _SNR_alpha;
  //    const double& SNR_efficiency = _SNR_efficiency;
  //    const double& SNR_E_c = _SNR_E_c;
  //    const double& pulsar_fraction = _pulsar_fraction;
  //    const double& pulsar_mean_period = _pulsar_mean_period;
  //    const double& pulsar_sdev_period = _pulsar_sdev_period;
  //    const double& PWN_alpha_low = _PWN_alpha_low;
  //    const double& PWN_alpha_high = _PWN_alpha_high;
  //    const double& PWN_E_break = _PWN_E_break;
  //    const double& PWN_efficiency = _PWN_efficiency;
  //    const double& b_0 = _b_0;
  //    const double& f_UV = _f_UV;
  //    const double& B_halo = _B_halo;
  const double& max_time = _max_time;
  const double& sn_rate = _sn_rate;
  const double& time_step = _time_step;
  //    const double& max_radius = _max_radius;
  //    const RunMode& runMode = _runMode;
  const ParticleModel& particleModel = _particleModel;
  const SpiralModel& spiralModel = _spiralModel;
  const PID::PID& pid = _pid;
};

}  // namespace core
}  // namespace gryphon

#endif  // GRYPHON_CORE_INPUT_H