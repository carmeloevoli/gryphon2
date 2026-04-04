#include "gryphon/core/input.h"

#include "gryphon/utils/logging.h"

namespace gryphon {
namespace core {

std::string spiralModelToString(SpiralModel model) {
  switch (model) {
    case SpiralModel::Uniform:
      return "Uniform";
    case SpiralModel::Jelly:
      return "Jelly";
    case SpiralModel::Steiman2010:
      return "Steiman2010";
    case SpiralModel::Faucher2006:
      return "Faucher2006";
    case SpiralModel::Vallee2008:
      return "Vallee2008";
    default:
      return "Unknown";
  }
}

std::string transportModelToString(TransportModel model) {
  switch (model) {
    case TransportModel::PureDiffusion:
      return "PureDiffusion";
    case TransportModel::DiffusionLosses:
      return "DiffusionLosses";
    default:
      return "Unknown";
  }
}

std::string injectionModelToString(InjectionModel model) {
  switch (model) {
    case InjectionModel::SinglePowerLaw:
      return "SinglePowerLaw";
    case InjectionModel::GalacticRandom:
      return "GalacticRandom ";
    case InjectionModel::MSP:
      return "MSP";
    case InjectionModel::SecondaryPositrons:
      return "SecondaryPositrons";
    default:
      return "Unknown";
  }
}

Input::Input() {  // validator();
}

Input::Input(const std::string& filename) {
  // read_params_file(filename);
  // validator();
}

void Input::print() {
  LOGD << "Simulation name : " << _simname;
  LOGD << "seed : " << _seed;
  LOGD << "E_min : " << _E_min / cgs::GeV << " GeV";
  LOGD << "E_max : " << _E_max / cgs::GeV << " GeV";
  LOGD << "E_size : " << _E_size;
  LOGD << "halo size : " << _H / cgs::kpc << " kpc";
  LOGD << "disc size : " << _h / cgs::kpc << " kpc";
  LOGD << "galaxy radius : " << _R_g / cgs::kpc << " kpc";
  LOGD << "sun position : " << _R_sun / cgs::kpc << " kpc";
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
    LOGD << "inj Emax : " << _injEmax / cgs::GeV << " GeV";
  }
  LOGD << "inj efficiency : " << _injEfficiency;
  LOGD << "B field : " << _B_field / cgs::microgauss << " muG";
  LOGD << "SN rate : " << _sn_rate / (1. / cgs::year) << " yr-1";
  LOGD << "time step : " << _time_step / cgs::year << " yr";
  LOGD << "max time : " << _max_time / cgs::Myr << " Myr";
  LOGD << "PID : " << _pid;
  LOGD << "DoVarySlope : " << std::boolalpha << _doVarySlope;
  LOGD << "DoVaryEnergy : " << std::boolalpha << _doVaryEnergy;
  LOGD << "Source profile model : " << spiralModelToString(_spiralModel);
  LOGD << "Transport model : " << transportModelToString(_transportModel);
  LOGD << "Injection model : " << injectionModelToString(_injectionModel);
}

}  // namespace core
}  // namespace gryphon