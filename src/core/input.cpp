#include "gryphon/core/input.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "gryphon/utils/logging.h"

namespace gryphon {
namespace core {

namespace {

std::string trim(const std::string& value) {
  const auto begin = value.find_first_not_of(" \t\r\n");
  if (begin == std::string::npos) return "";
  const auto end = value.find_last_not_of(" \t\r\n");
  return value.substr(begin, end - begin + 1);
}

std::string stripComment(const std::string& line) {
  return trim(line.substr(0, line.find('#')));
}

std::string normalizeToken(const std::string& value) {
  std::string normalized;
  normalized.reserve(value.size());
  for (const unsigned char ch : value) {
    if (std::isalnum(ch)) normalized.push_back(static_cast<char>(std::tolower(ch)));
  }
  return normalized;
}

std::string basenameWithoutExtension(const std::string& filename) {
  const auto slash = filename.find_last_of("/\\");
  const std::string basename =
      slash == std::string::npos ? filename : filename.substr(slash + 1);
  const auto dot = basename.find_last_of('.');
  return dot == std::string::npos ? basename : basename.substr(0, dot);
}

std::pair<std::string, std::string> splitKeyValue(const std::string& line) {
  const auto eq = line.find('=');
  if (eq != std::string::npos) {
    return {trim(line.substr(0, eq)), trim(line.substr(eq + 1))};
  }

  std::istringstream input(line);
  std::string key;
  if (!(input >> key)) return {"", ""};

  std::string value;
  std::getline(input, value);
  return {trim(key), trim(value)};
}

[[noreturn]] void throwParseError(const std::string& filename, size_t line_number,
                                  const std::string& message) {
  throw std::invalid_argument(filename + ":" + std::to_string(line_number) + ": " + message);
}

double parseDoubleValue(const std::string& filename, size_t line_number, const std::string& key,
                        const std::string& value) {
  try {
    size_t parsed = 0;
    const double result = std::stod(value, &parsed);
    if (parsed != value.size()) {
      throwParseError(filename, line_number, "invalid numeric value for '" + key + "'");
    }
    return result;
  } catch (const std::exception&) {
    throwParseError(filename, line_number, "invalid numeric value for '" + key + "'");
  }
}

ulong parseUnsignedLongValue(const std::string& filename, size_t line_number, const std::string& key,
                             const std::string& value) {
  if (!value.empty() && value[0] == '-') {
    throwParseError(filename, line_number, "invalid integer value for '" + key + "'");
  }
  try {
    size_t parsed = 0;
    const unsigned long result = std::stoul(value, &parsed);
    if (parsed != value.size()) {
      throwParseError(filename, line_number, "invalid integer value for '" + key + "'");
    }
    return result;
  } catch (const std::exception&) {
    throwParseError(filename, line_number, "invalid integer value for '" + key + "'");
  }
}

bool parseBoolValue(const std::string& filename, size_t line_number, const std::string& key,
                    const std::string& value) {
  const auto normalized = normalizeToken(value);
  if (normalized == "1" || normalized == "true" || normalized == "yes" ||
      normalized == "on" || normalized == "enable" || normalized == "enabled") {
    return true;
  }
  if (normalized == "0" || normalized == "false" || normalized == "no" ||
      normalized == "off" || normalized == "disable" || normalized == "disabled") {
    return false;
  }
  throwParseError(filename, line_number, "invalid boolean value for '" + key + "'");
}

SpiralModel parseSpiralModelValue(const std::string& filename, size_t line_number,
                                  const std::string& value) {
  const auto normalized = normalizeToken(value);
  if (normalized == "uniform") return SpiralModel::Uniform;
  if (normalized == "jelly") return SpiralModel::Jelly;
  if (normalized == "steiman2010" || normalized == "steiman") return SpiralModel::Steiman2010;
  if (normalized == "faucher2006" || normalized == "faucher") return SpiralModel::Faucher2006;
  if (normalized == "vallee2008" || normalized == "vallee") return SpiralModel::Vallee2008;
  throwParseError(filename, line_number, "unknown spiral model '" + value + "'");
}

TransportModel parseTransportModelValue(const std::string& filename, size_t line_number,
                                        const std::string& value) {
  const auto normalized = normalizeToken(value);
  if (normalized == "purediffusion") return TransportModel::PureDiffusion;
  if (normalized == "diffusionlosses" || normalized == "losses") {
    return TransportModel::DiffusionLosses;
  }
  throwParseError(filename, line_number, "unknown transport model '" + value + "'");
}

InjectionModel parseInjectionModelValue(const std::string& filename, size_t line_number,
                                        const std::string& value) {
  const auto normalized = normalizeToken(value);
  if (normalized == "singlepowerlaw") return InjectionModel::SinglePowerLaw;
  if (normalized == "galacticrandom") return InjectionModel::GalacticRandom;
  if (normalized == "msp") return InjectionModel::MSP;
  if (normalized == "secondarypositrons") return InjectionModel::SecondaryPositrons;
  throwParseError(filename, line_number, "unknown injection model '" + value + "'");
}

core::PID parsePidValue(const std::string& filename, size_t line_number, const std::string& value) {
  const auto normalized = normalizeToken(value);
  if (normalized == "h" || normalized == "p" || normalized == "proton" || normalized == "11") {
    return core::H;
  }
  if (normalized == "he" || normalized == "helium" || normalized == "24" ||
      normalized == "42") {
    return core::He;
  }
  throwParseError(filename, line_number, "unknown PID '" + value + "'");
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
  read_params_file(filename);
}

void Input::read_params_file(const std::string& filename) {
  std::ifstream file(filename.c_str());
  if (!file.is_open()) {
    throw std::runtime_error("could not open input parameter file '" + filename + "'");
  }

  const auto derived_simname = basenameWithoutExtension(filename);
  if (!derived_simname.empty()) _simname = derived_simname;

  std::string line;
  size_t line_number = 0;
  while (std::getline(file, line)) {
    ++line_number;
    const std::string content = stripComment(line);
    if (content.empty()) continue;

    const auto key_value = splitKeyValue(content);
    const std::string key = trim(key_value.first);
    const std::string value = trim(key_value.second);
    if (key.empty()) continue;
    if (value.empty()) {
      throwParseError(filename, line_number, "missing value for '" + key + "'");
    }

    const auto normalized_key = normalizeToken(key);

    if (normalized_key == "simname" || normalized_key == "name") {
      _simname = value;
    } else if (normalized_key == "seed") {
      _seed = parseUnsignedLongValue(filename, line_number, key, value);
    } else if (normalized_key == "emin" || normalized_key == "simemin" ||
               normalized_key == "simenergymin" || normalized_key == "simemingev") {
      _E_min = parseDoubleValue(filename, line_number, key, value) * cgs::GeV;
    } else if (normalized_key == "emax" || normalized_key == "simemax" ||
               normalized_key == "simenergymax" || normalized_key == "simemaxgev") {
      _E_max = parseDoubleValue(filename, line_number, key, value) * cgs::GeV;
    } else if (normalized_key == "esize" || normalized_key == "simesize" ||
               normalized_key == "energysize") {
      _E_size = parseUnsignedLongValue(filename, line_number, key, value);
    } else if (normalized_key == "hkpc" || normalized_key == "halosize" ||
               normalized_key == "halosizekpc" || normalized_key == "halosizeinkpc") {
      _H = parseDoubleValue(filename, line_number, key, value) * cgs::kpc;
    } else if (normalized_key == "discsizepc" || normalized_key == "discheightpc") {
      _h = parseDoubleValue(filename, line_number, key, value) * cgs::pc;
    } else if (normalized_key == "discsizekpc" || normalized_key == "discheightkpc") {
      _h = parseDoubleValue(filename, line_number, key, value) * cgs::kpc;
    } else if (normalized_key == "rgkpc" || normalized_key == "galaxyradius" ||
               normalized_key == "galaxyradiuskpc") {
      _R_g = parseDoubleValue(filename, line_number, key, value) * cgs::kpc;
    } else if (normalized_key == "sunkpc" || normalized_key == "sunradius" ||
               normalized_key == "sunradiuskpc") {
      _R_sun = parseDoubleValue(filename, line_number, key, value) * cgs::kpc;
    } else if (normalized_key == "d0h" || normalized_key == "d0overh") {
      _D0_over_H = parseDoubleValue(filename, line_number, key, value) * cgs::kpc / cgs::Myr;
    } else if (normalized_key == "e0" || normalized_key == "refenergy" ||
               normalized_key == "refenergygev") {
      _E_0 = parseDoubleValue(filename, line_number, key, value) * cgs::GeV;
    } else if (normalized_key == "delta") {
      _delta = parseDoubleValue(filename, line_number, key, value);
    } else if (normalized_key == "ddelta") {
      _ddelta = parseDoubleValue(filename, line_number, key, value);
    } else if (normalized_key == "s") {
      _s = parseDoubleValue(filename, line_number, key, value);
    } else if (normalized_key == "eb" || normalized_key == "ebgev") {
      _E_b = parseDoubleValue(filename, line_number, key, value) * cgs::GeV;
    } else if (normalized_key == "a") {
      _a = parseDoubleValue(filename, line_number, key, value);
    } else if (normalized_key == "b") {
      _b = parseDoubleValue(filename, line_number, key, value);
    } else if (normalized_key == "r1" || normalized_key == "r1kpc") {
      _R1 = parseDoubleValue(filename, line_number, key, value) * cgs::kpc;
    } else if (normalized_key == "injslope" || normalized_key == "snralpha") {
      _injSlope = parseDoubleValue(filename, line_number, key, value);
    } else if (normalized_key == "injslopesigma") {
      _injSlopeSigma = parseDoubleValue(filename, line_number, key, value);
    } else if (normalized_key == "injemax" || normalized_key == "snrcutoff") {
      _injEmax = parseDoubleValue(filename, line_number, key, value) * cgs::GeV;
    } else if (normalized_key == "efficiency" || normalized_key == "injefficiency" ||
               normalized_key == "snrefficiency") {
      _injEfficiency = parseDoubleValue(filename, line_number, key, value);
    } else if (normalized_key == "bmug" || normalized_key == "bfield" ||
               normalized_key == "bfieldmug") {
      _B_field = parseDoubleValue(filename, line_number, key, value) * cgs::microgauss;
    } else if (normalized_key == "urad" || normalized_key == "uradevcm3") {
      _U_rad = parseDoubleValue(filename, line_number, key, value) * cgs::eV / cgs::cm3;
    } else if (normalized_key == "snrate" || normalized_key == "snrateyr" ||
               normalized_key == "snrateperyear") {
      _sn_rate = parseDoubleValue(filename, line_number, key, value) / cgs::year;
    } else if (normalized_key == "rate" || normalized_key == "snratepercentury") {
      _sn_rate = parseDoubleValue(filename, line_number, key, value) / 100. / cgs::year;
    } else if (normalized_key == "timestep" || normalized_key == "timestepsyr" ||
               normalized_key == "timestepperyear") {
      _time_step = parseDoubleValue(filename, line_number, key, value) * cgs::year;
    } else if (normalized_key == "maxtimemyr" || normalized_key == "maxtime" ||
               normalized_key == "maxtimeinmyr") {
      _max_time = parseDoubleValue(filename, line_number, key, value) * cgs::Myr;
    } else if (normalized_key == "pid" || normalized_key == "particletype") {
      _pid = parsePidValue(filename, line_number, value);
    } else if (normalized_key == "varyenergy" || normalized_key == "dovaryenergy") {
      _doVaryEnergy = parseBoolValue(filename, line_number, key, value);
    } else if (normalized_key == "varyslope" || normalized_key == "dovaryslope") {
      _doVarySlope = parseBoolValue(filename, line_number, key, value);
    } else if (normalized_key == "spiralmodel") {
      _spiralModel = parseSpiralModelValue(filename, line_number, value);
    } else if (normalized_key == "transportmodel" || normalized_key == "lossesmodel") {
      _transportModel = parseTransportModelValue(filename, line_number, value);
    } else if (normalized_key == "injectionmodel" || normalized_key == "snrmodel") {
      _injectionModel = parseInjectionModelValue(filename, line_number, value);
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
  LOGD << "inj slope : " << _injSlope;
  LOGD << "inj slope sigma : " << _injSlopeSigma;
  if (_injEmax > 0) {
    LOGD << "inj Emax : " << _injEmax / cgs::GeV << " GeV";
  } else {
    LOGD << "inj Emax : no cutoff";
  }
  LOGD << "inj efficiency : " << _injEfficiency;
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
