#include "gryphon/core/input.h"

#include <sstream>
#include <stdexcept>
#include <vector>

#include "gryphon/utils/logging.h"

namespace gryphon {
namespace core {

namespace {

const char* spiralModelToString(SpiralModel model) noexcept {
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

const char* transportModelToString(TransportModel model) noexcept {
  switch (model) {
    case TransportModel::PureDiffusion:
      return "PureDiffusion";
    case TransportModel::DiffusionLosses:
      return "DiffusionLosses";
    default:
      return "Unknown";
  }
}

const char* injectionModelToString(InjectionModel model) noexcept {
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

}  // namespace

Input::Input(const std::string& filename) {
  // read_params_file(filename);
  validate();
}

void Input::validate() const {
  std::vector<std::string> errors;
  auto addError = [&errors](const std::string& message) { errors.push_back(message); };

  if (_simname.empty()) addError("simname cannot be empty");

  if (!(_E_min > 0.)) addError("E_min must be > 0");
  if (!(_E_max > _E_min)) addError("E_max must be greater than E_min");
  if (_E_size < 2) addError("E_size must be >= 2");

  if (!(_H > 0.)) addError("H must be > 0");
  if (!(_h >= 0.)) addError("h must be >= 0");
  if (!(_R_g > 0.)) addError("R_g must be > 0");
  if (!(_R_sun >= 0. && _R_sun <= _R_g)) addError("R_sun must satisfy 0 <= R_sun <= R_g");

  if (!(_D0_over_H > 0.)) addError("D0_over_H must be > 0");
  if (!(_E_0 > 0.)) addError("E_0 must be > 0");
  if (!(_E_b > 0.)) addError("E_b must be > 0");

  if (!(_injSlope > 0.)) addError("injSlope must be > 0");
  if (!(_injSlopeSigma >= 0.)) addError("injSlopeSigma must be >= 0");
  if (!(_injEfficiency >= 0.)) addError("injEfficiency must be >= 0");
  if ((_injectionModel == InjectionModel::SinglePowerLaw ||
       _injectionModel == InjectionModel::GalacticRandom) &&
      !(_injEmax > _E_0)) {
    addError("injEmax must be greater than E_0 for the selected injection model");
  }

  if (!(_B_field >= 0.)) addError("B_field must be >= 0");
  if (!(_U_rad >= 0.)) addError("U_rad must be >= 0");

  if (!(_sn_rate > 0.)) addError("sn_rate must be > 0");
  if (!(_time_step > 0.)) addError("time_step must be > 0");
  if (!(_max_time > 0.)) addError("max_time must be > 0");
  if (_time_step * _sn_rate > 1.) {
    addError("time_step * sn_rate must be <= 1 for galaxy generation");
  }

  if (!errors.empty()) {
    std::ostringstream message;
    message << "Input validation failed:";
    for (const auto& error : errors) {
      message << "\n - " << error;
    }
    throw std::invalid_argument(message.str());
  }
}

void Input::print() const {
  validate();
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