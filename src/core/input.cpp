#include "gryphon/core/input.h"

#include "gryphon/utils/logging.h"

namespace gryphon {
namespace core {

Input::Input() {  // validator();
}

Input::Input(const std::string& filename) {
  // read_params_file(filename);
  // validator();
}

void Input::print() {
  LOGD << "seed : " << _seed;
  LOGD << "E_min : " << _E_min / cgs::GeV << " GeV";
  LOGD << "E_max : " << _E_max / cgs::GeV << " GeV";
  LOGD << "E_size : " << _E_size;
  LOGD << "halo size : " << _H / cgs::kpc << " kpc";
  LOGD << "disc size : " << _h / cgs::kpc << " kpc";
  LOGD << "galaxy radius : " << _Rg / cgs::kpc << " kpc";
  LOGD << "sun position : " << _Rsun / cgs::kpc << " kpc";
  LOGD << "D_0 / H : " << _D0_over_H / (cgs::kpc / cgs::Myr) << " kpc/Myr";
  LOGD << "E_0 : " << _E_0 / cgs::GeV << " GeV";
  LOGD << "delta : " << _delta;
  LOGD << "ddelta : " << _ddelta;
  LOGD << "s : " << _s;
  LOGD << "E_b : " << _E_b / cgs::GeV << " GeV";
  LOGD << "a : " << _a;
  LOGD << "b : " << _b;
  LOGD << "R_1 : " << _R1 / cgs::kpc << " kpc";
  LOGD << "inj slope : " << _injSlope;
  LOGD << "inj slope sigma : " << _injSlopeSigma;
  if (_injEmax > 0) {
    LOGD << "inj Emax : " << _injEmax / cgs::GeV << "GeV";
  }
  LOGD << "inj efficiency : " << _injEfficiency;
  LOGD << "max time : " << _max_time / cgs::Myr << " Myr";
  LOGD << "SN rate : " << _sn_rate / (1. / cgs::year) << " yr-1";
  LOGD << "time step : " << _time_step / cgs::year << " yr";
  LOGD << "PID : " << _pid;
  LOGD << "DoVarySlope : " << std::boolalpha << _doVarySlope;
  LOGD << "DoVaryEnergy : " << std::boolalpha << _doVaryEnergy;
}

}  // namespace core
}  // namespace gryphon