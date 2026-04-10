#ifndef GRYPHON_CORE_INPUT_H
#define GRYPHON_CORE_INPUT_H

#include <string>
#include <utility>

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

enum class TransportModel { PureDiffusion, DiffusionLosses };

enum class InjectionModel {
  SinglePowerLaw,
  GalacticRandom,
  RandomEmax,
  PWN,
  MSP,
  SecondaryPositrons
};

namespace gryphon {
namespace core {

class Input {
 private:
  std::string _simname = "test";
  ulong _seed = 69;
  // output energy vector
  double _E_min = cgs::TeV;
  double _E_max = cgs::PeV;
  ulong _E_size = 3 * 16;
  // galaxy size
  double _H = 4. * cgs::kpc;
  double _h = 50. * cgs::pc;
  double _R_g = 20. * cgs::kpc;
  double _R_sun = 8.5 * cgs::kpc;
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
  double _injEmax = cgs::PeV;
  double _injEmaxSigmaDex = 0.;
  double _injEmaxMin = 10. * cgs::GeV;
  double _injEmaxMax = 10. * cgs::PeV;
  double _injEfficiency = 0.1;
  // PWN spectrum
  double _pwnP0 = 0.1 * cgs::second;
  double _pwnSigmaP0 = 0.;
  bool _pwnRandomInitialPeriod = true;
  double _pwnAlpha1 = 1.5;
  double _pwnAlpha2 = 2.6;
  double _pwnEbreak = 100. * cgs::GeV;
  double _pwnEmin = 1. * cgs::GeV;
  // Energy losses
  double _B_field = cgs::microgauss;
  double _U_rad = 0.25 * cgs::eV / cgs::cm3;
  // simulation parameters
  double _sn_rate = 1. / 50. / cgs::year;
  double _time_step = 1. * cgs::year;
  double _max_time = 100. * cgs::Myr;
  // models
  core::PID _pid = core::H;
  bool _doVarySlope = false;
  bool _doVaryEnergy = false;
  SpiralModel _spiralModel = SpiralModel::Uniform;
  TransportModel _transportModel = TransportModel::PureDiffusion;
  InjectionModel _injectionModel = InjectionModel::SinglePowerLaw;

 public:
  Input() = default;
  explicit Input(const std::string& filename);
  Input(const Input& other) = default;
  Input(Input&& other) noexcept = default;
  Input& operator=(const Input& other) = default;
  Input& operator=(Input&& other) noexcept = default;
  virtual ~Input() = default;
  void read_params_file(const std::string& filename);
  void validate() const;
  void print() const;

  inline void set_simname(std::string name) { _simname = std::move(name); }
  inline void set_simEmin(double E_min) noexcept { _E_min = E_min; }
  inline void set_simEmax(double E_max) noexcept { _E_max = E_max; }
  inline void set_simEsize(ulong E_size) noexcept { _E_size = E_size; }
  inline void set_seed(unsigned long int seed) noexcept { _seed = seed; }
  inline void set_refEnergy(double E_0) noexcept { _E_0 = E_0; }
  inline void set_maxtime(double time) noexcept { _max_time = time; }
  inline void set_halosize(double H) noexcept { _H = H; }
  inline void set_discsize(double h) noexcept { _h = h; }
  inline void set_galaxyRadius(double R_g) noexcept { _R_g = R_g; }
  inline void set_sunRadius(double R_sun) noexcept { _R_sun = R_sun; }
  inline void set_injSlope(double slope) noexcept { _injSlope = slope; }
  inline void set_injSlopeSigma(double sigma) noexcept { _injSlopeSigma = sigma; }
  inline void set_injEmax(double Emax) noexcept { _injEmax = Emax; }
  inline void set_injEmaxSigmaDex(double sigma_dex) noexcept { _injEmaxSigmaDex = sigma_dex; }
  inline void set_injEmaxMin(double Emin) noexcept { _injEmaxMin = Emin; }
  inline void set_injEmaxMax(double Emax) noexcept { _injEmaxMax = Emax; }
  inline void set_efficiency(double epsilon) noexcept { _injEfficiency = epsilon; }
  inline void set_pwnP0(double P0) noexcept { _pwnP0 = P0; }
  inline void set_pwnSigmaP0(double sigmaP0) noexcept { _pwnSigmaP0 = sigmaP0; }
  inline void set_pwnRandomInitialPeriod(bool doRandom) noexcept {
    _pwnRandomInitialPeriod = doRandom;
  }
  inline void set_pwnAlpha1(double alpha) noexcept { _pwnAlpha1 = alpha; }
  inline void set_pwnAlpha2(double alpha) noexcept { _pwnAlpha2 = alpha; }
  inline void set_pwnEbreak(double Ebreak) noexcept { _pwnEbreak = Ebreak; }
  inline void set_pwnEmin(double Emin) noexcept { _pwnEmin = Emin; }
  inline void set_rate(double rate) noexcept { _sn_rate = rate; }
  inline void set_D0_over_H(double D0_over_H) noexcept { _D0_over_H = D0_over_H; }
  inline void set_delta(double delta) noexcept { _delta = delta; }
  inline void set_Bfield(double B) noexcept { _B_field = B; }
  inline void set_Urad(double U) noexcept { _U_rad = U; }
  inline void set_pid(core::PID pid) noexcept { _pid = pid; }
  inline void enable_varyenergy() noexcept { _doVaryEnergy = true; }
  inline void disable_varyenergy() noexcept { _doVaryEnergy = false; }
  inline void enable_varyslope() noexcept { _doVarySlope = true; }
  inline void disable_varyslope() noexcept { _doVarySlope = false; }
  inline void set_spiralModel(SpiralModel model) noexcept { _spiralModel = model; }
  inline void set_transportModel(TransportModel model) noexcept { _transportModel = model; }
  inline void set_injectionModel(InjectionModel model) noexcept { _injectionModel = model; }

  const std::string& simname() const noexcept { return _simname; }
  unsigned long int seed() const noexcept { return _seed; }
  double simEmin() const noexcept { return _E_min; }
  double simEmax() const noexcept { return _E_max; }
  ulong simEsize() const noexcept { return _E_size; }
  double E_min() const noexcept { return _E_min; }
  double E_max() const noexcept { return _E_max; }
  ulong E_size() const noexcept { return _E_size; }
  double H() const noexcept { return _H; }
  double h() const noexcept { return _h; }
  double R_g() const noexcept { return _R_g; }
  double R_sun() const noexcept { return _R_sun; }
  double D0_over_H() const noexcept { return _D0_over_H; }
  double E_0() const noexcept { return _E_0; }
  double delta() const noexcept { return _delta; }
  double ddelta() const noexcept { return _ddelta; }
  double s() const noexcept { return _s; }
  double E_b() const noexcept { return _E_b; }
  double a() const noexcept { return _a; }
  double b() const noexcept { return _b; }
  double R_1() const noexcept { return _R1; }
  double injSlope() const noexcept { return _injSlope; }
  double injSlopeSigma() const noexcept { return _injSlopeSigma; }
  double injEmax() const noexcept { return _injEmax; }
  double injEmaxSigmaDex() const noexcept { return _injEmaxSigmaDex; }
  double injEmaxMin() const noexcept { return _injEmaxMin; }
  double injEmaxMax() const noexcept { return _injEmaxMax; }
  double injEfficiency() const noexcept { return _injEfficiency; }
  double pwnP0() const noexcept { return _pwnP0; }
  double pwnSigmaP0() const noexcept { return _pwnSigmaP0; }
  bool pwnRandomInitialPeriod() const noexcept { return _pwnRandomInitialPeriod; }
  double pwnAlpha1() const noexcept { return _pwnAlpha1; }
  double pwnAlpha2() const noexcept { return _pwnAlpha2; }
  double pwnEbreak() const noexcept { return _pwnEbreak; }
  double pwnEmin() const noexcept { return _pwnEmin; }
  double max_time() const noexcept { return _max_time; }
  double sn_rate() const noexcept { return _sn_rate; }
  double time_step() const noexcept { return _time_step; }
  core::PID pid() const noexcept { return _pid; }
  bool doVarySlope() const noexcept { return _doVarySlope; }
  bool doVaryEnergy() const noexcept { return _doVaryEnergy; }
  SpiralModel spiralModel() const noexcept { return _spiralModel; }
  TransportModel transportModel() const noexcept { return _transportModel; }
  InjectionModel injectionModel() const noexcept { return _injectionModel; }
  double B_field() const noexcept { return _B_field; }
  double U_rad() const noexcept { return _U_rad; }
};

}  // namespace core
}  // namespace gryphon

#endif  // GRYPHON_CORE_INPUT_H
