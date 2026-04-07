#include "gryphon/core/input.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "gryphon/utils/io.h"
#include "gryphon/utils/logging.h"

namespace gryphon {
namespace core {

namespace {

namespace fs = std::filesystem;

[[noreturn]] void throwParseError(const std::string& filename, size_t line_number,
                                  const std::string& message) {
  throw std::invalid_argument(filename + ":" + std::to_string(line_number) + ": " + message);
}

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
    case InjectionModel::PWN:
      return "PWN";
    case InjectionModel::MSP:
      return "MSP";
    case InjectionModel::SecondaryPositrons:
      return "SecondaryPositrons";
    default:
      return "Unknown";
  }
}

}  // namespace

Input::Input(const std::string& filename) { read_params_file(filename); }

void Input::read_params_file(const std::string& filename) {
  const fs::path input_path(filename);
  std::ifstream file(input_path);
  if (!file.is_open()) {
    throw std::runtime_error("could not open input parameter file '" + filename + "'");
  }

  const auto derived_simname = input_path.stem().string();
  if (!derived_simname.empty()) _simname = derived_simname;

  std::string line;
  size_t line_number = 0;
  while (std::getline(file, line)) {
    ++line_number;
    const std::string content = utils::stripComment(line);
    if (content.empty()) continue;

    const auto key_value = utils::splitKeyValue(content);
    const std::string key = utils::trim(key_value.first);
    const std::string value = utils::trim(key_value.second);
    if (key.empty()) continue;
    if (value.empty()) {
      throwParseError(filename, line_number, "missing value for '" + key + "'");
    }

    const auto normalized_key = utils::normalizeToken(key);

    if (normalized_key == "simname" || normalized_key == "name") {
      _simname = value;
    } else if (normalized_key == "seed") {
      _seed = utils::parseUnsignedLongValue(filename, line_number, key, value);
    } else if (normalized_key == "emin" || normalized_key == "simemin" ||
               normalized_key == "simenergymin" || normalized_key == "simemingev") {
      _E_min = utils::parseDoubleValue(filename, line_number, key, value) * cgs::GeV;
    } else if (normalized_key == "emax" || normalized_key == "simemax" ||
               normalized_key == "simenergymax" || normalized_key == "simemaxgev") {
      _E_max = utils::parseDoubleValue(filename, line_number, key, value) * cgs::GeV;
    } else if (normalized_key == "esize" || normalized_key == "simesize" ||
               normalized_key == "energysize") {
      _E_size = utils::parseUnsignedLongValue(filename, line_number, key, value);
    } else if (normalized_key == "hkpc" || normalized_key == "halosize" ||
               normalized_key == "halosizekpc" || normalized_key == "halosizeinkpc") {
      _H = utils::parseDoubleValue(filename, line_number, key, value) * cgs::kpc;
    } else if (normalized_key == "discsizepc" || normalized_key == "discheightpc") {
      _h = utils::parseDoubleValue(filename, line_number, key, value) * cgs::pc;
    } else if (normalized_key == "discsizekpc" || normalized_key == "discheightkpc") {
      _h = utils::parseDoubleValue(filename, line_number, key, value) * cgs::kpc;
    } else if (normalized_key == "rgkpc" || normalized_key == "galaxyradius" ||
               normalized_key == "galaxyradiuskpc") {
      _R_g = utils::parseDoubleValue(filename, line_number, key, value) * cgs::kpc;
    } else if (normalized_key == "sunkpc" || normalized_key == "sunradius" ||
               normalized_key == "sunradiuskpc") {
      _R_sun = utils::parseDoubleValue(filename, line_number, key, value) * cgs::kpc;
    } else if (normalized_key == "d0h" || normalized_key == "d0overh") {
      _D0_over_H =
          utils::parseDoubleValue(filename, line_number, key, value) * cgs::kpc / cgs::Myr;
    } else if (normalized_key == "e0" || normalized_key == "refenergy" ||
               normalized_key == "refenergygev") {
      _E_0 = utils::parseDoubleValue(filename, line_number, key, value) * cgs::GeV;
    } else if (normalized_key == "delta") {
      _delta = utils::parseDoubleValue(filename, line_number, key, value);
    } else if (normalized_key == "ddelta") {
      _ddelta = utils::parseDoubleValue(filename, line_number, key, value);
    } else if (normalized_key == "s") {
      _s = utils::parseDoubleValue(filename, line_number, key, value);
    } else if (normalized_key == "eb" || normalized_key == "ebgev") {
      _E_b = utils::parseDoubleValue(filename, line_number, key, value) * cgs::GeV;
    } else if (normalized_key == "a") {
      _a = utils::parseDoubleValue(filename, line_number, key, value);
    } else if (normalized_key == "b") {
      _b = utils::parseDoubleValue(filename, line_number, key, value);
    } else if (normalized_key == "r1" || normalized_key == "r1kpc") {
      _R1 = utils::parseDoubleValue(filename, line_number, key, value) * cgs::kpc;
    } else if (normalized_key == "injslope" || normalized_key == "snralpha") {
      _injSlope = utils::parseDoubleValue(filename, line_number, key, value);
    } else if (normalized_key == "injslopesigma") {
      _injSlopeSigma = utils::parseDoubleValue(filename, line_number, key, value);
    } else if (normalized_key == "injemax" || normalized_key == "snrcutoff") {
      _injEmax = utils::parseDoubleValue(filename, line_number, key, value) * cgs::GeV;
    } else if (normalized_key == "efficiency" || normalized_key == "injefficiency" ||
               normalized_key == "snrefficiency") {
      _injEfficiency = utils::parseDoubleValue(filename, line_number, key, value);
    } else if (normalized_key == "p0" || normalized_key == "p0sec" ||
               normalized_key == "pwnp0" || normalized_key == "pwnp0sec") {
      _pwnP0 = utils::parseDoubleValue(filename, line_number, key, value) * cgs::second;
    } else if (normalized_key == "p0ms" || normalized_key == "pwnp0ms") {
      _pwnP0 = utils::parseDoubleValue(filename, line_number, key, value) * cgs::msec;
    } else if (normalized_key == "sigmap0" || normalized_key == "sigmap0sec" ||
               normalized_key == "pwnsigmap0" || normalized_key == "pwnsigmap0sec") {
      _pwnSigmaP0 =
          utils::parseDoubleValue(filename, line_number, key, value) * cgs::second;
    } else if (normalized_key == "sigmap0ms" || normalized_key == "pwnsigmap0ms") {
      _pwnSigmaP0 = utils::parseDoubleValue(filename, line_number, key, value) * cgs::msec;
    } else if (normalized_key == "pwnalpha1" || normalized_key == "pwnalphabelowbreak" ||
               normalized_key == "pwnalphalow") {
      _pwnAlpha1 = utils::parseDoubleValue(filename, line_number, key, value);
    } else if (normalized_key == "pwnalpha2" || normalized_key == "pwnalphaabovebreak" ||
               normalized_key == "pwnalphahigh") {
      _pwnAlpha2 = utils::parseDoubleValue(filename, line_number, key, value);
    } else if (normalized_key == "pwnebreak" || normalized_key == "pwnebreakgev" ||
               normalized_key == "pwnbreak" || normalized_key == "pwnbreakgev") {
      _pwnEbreak = utils::parseDoubleValue(filename, line_number, key, value) * cgs::GeV;
    } else if (normalized_key == "pwnemin" || normalized_key == "pwnemingev") {
      _pwnEmin = utils::parseDoubleValue(filename, line_number, key, value) * cgs::GeV;
    } else if (normalized_key == "pwnecut" || normalized_key == "pwnecutgev" ||
               normalized_key == "pwncutoff" || normalized_key == "pwncutoffgev" ||
               normalized_key == "pwnemax" || normalized_key == "pwnemaxgev") {
      throwParseError(
          filename, line_number,
          "PWN cutoff energy is derived from the potential drop and must not be set explicitly");
    } else if (normalized_key == "bmug" || normalized_key == "bfield" ||
               normalized_key == "bfieldmug") {
      _B_field = utils::parseDoubleValue(filename, line_number, key, value) * cgs::microgauss;
    } else if (normalized_key == "urad" || normalized_key == "uradevcm3") {
      _U_rad = utils::parseDoubleValue(filename, line_number, key, value) * cgs::eV / cgs::cm3;
    } else if (normalized_key == "snrate" || normalized_key == "snrateyr" ||
               normalized_key == "snrateperyear") {
      _sn_rate = utils::parseDoubleValue(filename, line_number, key, value) / cgs::year;
    } else if (normalized_key == "rate" || normalized_key == "snratepercentury") {
      _sn_rate = utils::parseDoubleValue(filename, line_number, key, value) / 100. / cgs::year;
    } else if (normalized_key == "timestep" || normalized_key == "timestepsyr" ||
               normalized_key == "timestepperyear") {
      _time_step = utils::parseDoubleValue(filename, line_number, key, value) * cgs::year;
    } else if (normalized_key == "maxtimemyr" || normalized_key == "maxtime" ||
               normalized_key == "maxtimeinmyr") {
      _max_time = utils::parseDoubleValue(filename, line_number, key, value) * cgs::Myr;
    } else if (normalized_key == "pid" || normalized_key == "particletype") {
      _pid = utils::parsePidValue(filename, line_number, value);
    } else if (normalized_key == "varyenergy" || normalized_key == "dovaryenergy") {
      _doVaryEnergy = utils::parseBoolValue(filename, line_number, key, value);
    } else if (normalized_key == "varyslope" || normalized_key == "dovaryslope") {
      _doVarySlope = utils::parseBoolValue(filename, line_number, key, value);
    } else if (normalized_key == "spiralmodel") {
      _spiralModel = utils::parseSpiralModelValue(filename, line_number, value);
    } else if (normalized_key == "transportmodel" || normalized_key == "lossesmodel") {
      _transportModel = utils::parseTransportModelValue(filename, line_number, value);
    } else if (normalized_key == "injectionmodel" || normalized_key == "snrmodel") {
      _injectionModel = utils::parseInjectionModelValue(filename, line_number, value);
    } else {
      throwParseError(filename, line_number, "unknown parameter '" + key + "'");
    }
  }

  validate();
}

void Input::validate() const {
  std::vector<std::string> errors;
  auto addError = [&errors](const std::string& message) { errors.push_back(message); };

  if (_simname.empty()) addError("simname cannot be empty");

  if (!(_E_min > 1. * cgs::GeV)) addError("E_min must be greater than 1 GeV");
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
  if (_injectionModel == InjectionModel::PWN && !(_pwnP0 > 0.)) {
    addError("PWN mean P0 must be > 0");
  }
  if (_injectionModel == InjectionModel::PWN && !(_pwnSigmaP0 >= 0.)) {
    addError("PWN sigmaP0 must be >= 0");
  }
  if (_injectionModel == InjectionModel::PWN && !(_pwnAlpha1 > 0.)) {
    addError("PWN alpha1 must be > 0");
  }
  if (_injectionModel == InjectionModel::PWN && !(_pwnAlpha2 > 0.)) {
    addError("PWN alpha2 must be > 0");
  }
  if (_injectionModel == InjectionModel::PWN && !(_pwnEmin > 0.)) {
    addError("PWN Emin must be > 0");
  }
  if (_injectionModel == InjectionModel::PWN && !(_pwnEmin < _pwnEbreak)) {
    addError("PWN Emin must be smaller than Ebreak");
  }
  if ((_injectionModel == InjectionModel::SinglePowerLaw ||
       _injectionModel == InjectionModel::GalacticRandom) &&
      !(_injEmax > 1. * cgs::GeV)) {
    addError("injEmax must be greater than 1 GeV for the selected injection model");
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
  LOGD << "inj efficiency : " << _injEfficiency;
  if (_injectionModel == InjectionModel::PWN) {
    LOGD << "PWN alpha below break : " << _pwnAlpha1;
    LOGD << "PWN alpha above break : " << _pwnAlpha2;
    LOGD << "PWN Emin : " << _pwnEmin / cgs::GeV << " GeV";
    LOGD << "PWN Ebreak : " << _pwnEbreak / cgs::GeV << " GeV";
    LOGD << "PWN Ecut (= Emax) : derived from potential drop";
    LOGD << "PWN mean P0 : " << _pwnP0 / cgs::msec << " ms";
    LOGD << "PWN sigma P0 : " << _pwnSigmaP0 / cgs::msec << " ms";
  } else {
    LOGD << "inj slope : " << _injSlope;
    LOGD << "inj slope sigma : " << _injSlopeSigma;
    if (_injEmax > 0) {
      LOGD << "inj Emax : " << _injEmax / cgs::GeV << " GeV";
    } else {
      LOGD << "inj Emax : no cutoff";
    }
  }
  LOGD << "B field : " << _B_field / cgs::microgauss << " muG";
  LOGD << "U_B : " << (_B_field * _B_field / (8. * M_PI)) / (cgs::eV / cgs::cm3) << " eV/cm3";
  LOGD << "U_rad : " << _U_rad / (cgs::eV / cgs::cm3) << " eV/cm3";
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
