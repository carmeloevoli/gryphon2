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
  Xie2024,      // Xie et al., arXiv:2402.11428, 2024
  Faucher2006,  // Faucher-Giguere and Kaspi, ApJ, 643, 332–355, 2006
  Vallee2008    // Vallee, AJ, 135, 1301-1310, 2008
};

enum class TransportModel { PureDiffusion, DiffusionLosses };

enum class InjectionModel {
  SinglePowerLaw,
  SmoothBrokenPowerLaw,
  GalacticRandom,
  RandomEmax,
  PWN,
  YoungPulsars,
  MSP,
  SecondaryPositrons
};

namespace gryphon {
namespace core {

class Input {
 private:
  std::string _simname = "test";
  std::string _configFile;
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
  // Change of diffusion slope above E_b.  -1 retains the legacy unbroken
  // power law; any other value enables the smooth break.
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
  double _injDeltaSlope = 1.;
  double _injBreakEnergy = 3.3 * cgs::PeV;
  double _injSmoothness = 0.27;
  // GalacticRandom additionally accepts 0 for an exactly uncut power law;
  // then the parent index is >2 and varying indices are conditioned on >2.
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
  // Young pulsar burst-injection spectrum
  double _youngPulsarsP0 = 60. * cgs::msec;
  double _youngPulsarsSigmaP0 = 10. * cgs::msec;
  bool _youngPulsarsRandomInitialPeriod = true;
  double _youngPulsarsB0 = 2.5e12 * cgs::gauss;
  double _youngPulsarsSigmaLog10B = 0.5;
  bool _youngPulsarsRandomMagneticField = true;
  // Present-day MSP snapshot; stationary pair injection (not SN bursts).
  ulong _mspSources = 49000;
  double _mspPeriodMin = 1.5 * cgs::msec;
  double _mspLog10B = 8.;
  double _mspSigmaLog10B = 0.2;
  double _mspInertia = 1e45 * cgs::gram * cgs::cm2;
  double _mspEmin = 10. * cgs::GeV;
  ulong _mspQuadrature = 256;
  ulong _mspDistances = 1601;
  // Energy losses
  double _B_field = cgs::microgauss;
  double _U_rad = 0.25 * cgs::eV / cgs::cm3;
  // simulation parameters
  double _sn_rate = 1. / 50. / cgs::year;
  double _time_step = 1. * cgs::year;
  double _max_time = 100. * cgs::Myr;
  // Opt-in stationary association population; the legacy generator is unchanged.
  bool _syntheticAssociations = false;
  double _associationFraction = 1.;
  ulong _associationMembers = 100;
  double _associationRadius = 30. * cgs::pc;
  double _associationVelocity = 3. * cgs::km / cgs::second;
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
  // Dump every parameter in a form read_params_file() accepts, so that a run
  // can be reproduced from the file stored next to its output.
  void write_params_file(const std::string& filename) const;
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
  inline void set_injDeltaSlope(double delta) noexcept { _injDeltaSlope = delta; }
  inline void set_injBreakEnergy(double energy) noexcept { _injBreakEnergy = energy; }
  inline void set_injSmoothness(double smoothness) noexcept { _injSmoothness = smoothness; }
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
  inline void set_youngPulsarsP0(double P0) noexcept { _youngPulsarsP0 = P0; }
  inline void set_youngPulsarsSigmaP0(double sigmaP0) noexcept {
    _youngPulsarsSigmaP0 = sigmaP0;
  }
  inline void set_youngPulsarsRandomInitialPeriod(bool doRandom) noexcept {
    _youngPulsarsRandomInitialPeriod = doRandom;
  }
  inline void set_youngPulsarsB0(double B0) noexcept { _youngPulsarsB0 = B0; }
  inline void set_youngPulsarsSigmaLog10B(double sigmaLog10B) noexcept {
    _youngPulsarsSigmaLog10B = sigmaLog10B;
  }
  inline void set_youngPulsarsRandomMagneticField(bool doRandom) noexcept {
    _youngPulsarsRandomMagneticField = doRandom;
  }
  inline void set_rate(double rate) noexcept { _sn_rate = rate; }
  void set_mspSources(ulong value) noexcept { _mspSources = value; }
  void set_mspPeriodMin(double value) noexcept { _mspPeriodMin = value; }
  void set_mspLog10B(double value) noexcept { _mspLog10B = value; }
  void set_mspSigmaLog10B(double value) noexcept { _mspSigmaLog10B = value; }
  void set_mspInertia(double value) noexcept { _mspInertia = value; }
  void set_mspEmin(double value) noexcept { _mspEmin = value; }
  void set_mspQuadrature(ulong value) noexcept { _mspQuadrature = value; }
  void set_mspDistances(ulong value) noexcept { _mspDistances = value; }
  void set_syntheticAssociations(bool enabled) noexcept { _syntheticAssociations = enabled; }
  void set_associationFraction(double fraction) noexcept { _associationFraction = fraction; }
  void set_associationMembers(ulong members) noexcept { _associationMembers = members; }
  void set_associationRadius(double radius) noexcept { _associationRadius = radius; }
  void set_associationVelocity(double velocity) noexcept { _associationVelocity = velocity; }
  inline void set_D0_over_H(double D0_over_H) noexcept { _D0_over_H = D0_over_H; }
  inline void set_delta(double delta) noexcept { _delta = delta; }
  inline void set_ddelta(double ddelta) noexcept { _ddelta = ddelta; }
  inline void set_diffusionSmoothness(double smoothness) noexcept { _s = smoothness; }
  inline void set_diffusionBreakEnergy(double energy) noexcept { _E_b = energy; }
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
  // Path of the parameter file this input was read from, empty if built in code.
  const std::string& configFile() const noexcept { return _configFile; }
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
  double injDeltaSlope() const noexcept { return _injDeltaSlope; }
  double injBreakEnergy() const noexcept { return _injBreakEnergy; }
  double injSmoothness() const noexcept { return _injSmoothness; }
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
  double youngPulsarsP0() const noexcept { return _youngPulsarsP0; }
  double youngPulsarsSigmaP0() const noexcept { return _youngPulsarsSigmaP0; }
  bool youngPulsarsRandomInitialPeriod() const noexcept {
    return _youngPulsarsRandomInitialPeriod;
  }
  double youngPulsarsB0() const noexcept { return _youngPulsarsB0; }
  double youngPulsarsSigmaLog10B() const noexcept { return _youngPulsarsSigmaLog10B; }
  bool youngPulsarsRandomMagneticField() const noexcept {
    return _youngPulsarsRandomMagneticField;
  }
  double max_time() const noexcept { return _max_time; }
  ulong mspSources() const noexcept { return _mspSources; }
  double mspPeriodMin() const noexcept { return _mspPeriodMin; }
  double mspLog10B() const noexcept { return _mspLog10B; }
  double mspSigmaLog10B() const noexcept { return _mspSigmaLog10B; }
  double mspInertia() const noexcept { return _mspInertia; }
  double mspEmin() const noexcept { return _mspEmin; }
  ulong mspQuadrature() const noexcept { return _mspQuadrature; }
  ulong mspDistances() const noexcept { return _mspDistances; }
  double sn_rate() const noexcept { return _sn_rate; }
  bool syntheticAssociations() const noexcept { return _syntheticAssociations; }
  double associationFraction() const noexcept { return _associationFraction; }
  ulong associationMembers() const noexcept { return _associationMembers; }
  double associationRadius() const noexcept { return _associationRadius; }
  double associationVelocity() const noexcept { return _associationVelocity; }
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
