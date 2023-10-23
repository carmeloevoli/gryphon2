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

namespace gryphon {
namespace core {

class Input {
 private:
  ulong _seed = 69;
  std::string _simname = "test";
  // output energy vector
  double _E_min = cgs::TeV;
  double _E_max = cgs::PeV;
  ulong _E_size = 3 * 16;
  // galaxy size
  double _H = 4. * cgs::kpc;
  double _h = 50. * cgs::pc;
  double _Rg = 20. * cgs::kpc;
  double _Rsun = 8.5 * cgs::kpc;
  // diffusion coefficient parameters from Schroer+, PRD 103, 2001
  double _D0_over_H = 0.42 * cgs::kpc / cgs::Myr;
  double _E_0 = cgs::TeV;
  double _delta = 0.36;
  double _ddelta = -1.;
  double _s = 0.1;
  double _E_b = 312. * cgs::GeV;
  // source profile parameters from Lorimer2006
  double _a = 1.9;
  double _b = 5.0;
  double _R1 = 0.;
  // SNR spectrum
  double _injSlope = 2.34;
  double _injSlopeSigma = 0.15;
  double _injEmax = -1;
  double _injEfficiency = 0.1;
  // simulation parameters
  double _sn_rate = 1. / 50. / cgs::year;
  double _time_step = 1. * cgs::year;
  double _max_time = 100. * cgs::Myr;
  // models
  core::PID _pid = core::H;
  bool _doVarySlope = false;
  bool _doVaryEnergy = false;
  SpiralModel _spiralModel = SpiralModel::Uniform;

 protected:
  //   void read_params_file(const std::string& filename);
  //   void set_outputdir(const std::string& dir);
  //   void validator();

 public:
  Input();
  Input(const std::string& filename);
  virtual ~Input() = default;
  void print();

  inline void set_simname(const std::string& name) { _simname = name; }
  inline void set_simEmin(const double& E_min) { _E_min = E_min; }
  inline void set_simEmax(const double& E_max) { _E_max = E_max; }
  inline void set_simEsize(const ulong& E_size) { _E_size = E_size; }
  inline void set_seed(const unsigned long int& seed) { _seed = seed; }
  inline void set_maxtime(const double& time) { _max_time = time; }
  inline void set_halosize(const double& H) { _H = H; }
  inline void set_injEmax(const double& Emax) { _injEmax = Emax; }
  inline void set_efficiency(const double& epsilon) { _injEfficiency = epsilon; }
  inline void set_rate(const double& rate) { _sn_rate = rate; }
  inline void enable_varyenergy() { _doVaryEnergy = true; }
  inline void enable_varyslope() { _doVarySlope = true; }

  const std::string& simname = _simname;
  const unsigned long int& seed = _seed;
  const double& E_min = _E_min;
  const double& E_max = _E_max;
  const size_t& E_size = _E_size;
  const double& H = _H;
  const double& h = _h;
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
  const double& injSlope = _injSlope;
  const double& injSlopeSigma = _injSlopeSigma;
  const double& injEmax = _injEmax;
  const double& injEfficiency = _injEfficiency;
  const double& max_time = _max_time;
  const double& sn_rate = _sn_rate;
  const double& time_step = _time_step;
  const core::PID& pid = _pid;
  const bool& doVarySlope = _doVarySlope;
  const bool& doVaryEnergy = _doVaryEnergy;
  const SpiralModel& spiralModel = _spiralModel;
};

}  // namespace core
}  // namespace gryphon

#endif  // GRYPHON_CORE_INPUT_H