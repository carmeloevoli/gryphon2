#include "gryphon/utils/io.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "gryphon/utils/logging.h"

namespace gryphon {
namespace utils {

namespace {

namespace fs = std::filesystem;

[[noreturn]] void throwParseError(const std::string& filename, size_t line_number,
                                  const std::string& message) {
  throw std::invalid_argument(filename + ":" + std::to_string(line_number) + ": " + message);
}

bool isSkippableLine(const std::string& line) {
  const auto first = line.find_first_not_of(" \t\r\n");
  return first == std::string::npos || line[first] == '#';
}

std::ifstream openInputFile(const fs::path& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("could not open file '" + filename.string() + "'");
  }
  return file;
}

void ensureDirectoryExists(const fs::path& directory) {
  std::error_code error;
  const bool exists = fs::exists(directory, error);
  if (error) {
    throw std::runtime_error("could not inspect directory '" + directory.string() + "': " +
                             error.message());
  }

  if (exists) {
    if (!fs::is_directory(directory, error)) {
      if (error) {
        throw std::runtime_error("could not inspect directory '" + directory.string() + "': " +
                                 error.message());
      }
      throw std::runtime_error("output path '" + directory.string() +
                               "' exists and is not a directory");
    }
    return;
  }

  if (!fs::create_directories(directory, error) && error) {
    throw std::runtime_error("could not create directory '" + directory.string() + "': " +
                             error.message());
  }
}

}  // namespace

std::string trim(const std::string& value) {
  const auto begin = value.find_first_not_of(" \t\r\n");
  if (begin == std::string::npos) return "";
  const auto end = value.find_last_not_of(" \t\r\n");
  return value.substr(begin, end - begin + 1);
}

std::string stripComment(const std::string& line) { return trim(line.substr(0, line.find('#'))); }

std::string normalizeToken(const std::string& value) {
  std::string normalized;
  normalized.reserve(value.size());
  for (const unsigned char ch : value) {
    if (std::isalnum(ch)) normalized.push_back(static_cast<char>(std::tolower(ch)));
  }
  return normalized;
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

std::string removeExtensionIniFilename(std::string inputFilename) {
  fs::path input_path(std::move(inputFilename));
  if (input_path.extension() != ".ini") {
    throw std::invalid_argument("Wrong input filename! It must be : <filename>.ini");
  }

  input_path.replace_extension();
  return input_path.string();
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

unsigned long parseUnsignedLongValue(const std::string& filename, size_t line_number,
                                     const std::string& key, const std::string& value) {
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
  if (normalized == "1" || normalized == "true" || normalized == "yes" || normalized == "on" ||
      normalized == "enable" || normalized == "enabled") {
    return true;
  }
  if (normalized == "0" || normalized == "false" || normalized == "no" || normalized == "off" ||
      normalized == "disable" || normalized == "disabled") {
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
  if (normalized == "randomemax" || normalized == "randomcutoff") {
    return InjectionModel::RandomEmax;
  }
  if (normalized == "pwn" || normalized == "pulsarwindnebula") return InjectionModel::PWN;
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

unsigned long parseSeed(const char* arg) {
  if (arg == nullptr || *arg == '\0' || arg[0] == '-') {
    throw std::invalid_argument("Seed must be a non-negative integer");
  }

  errno = 0;
  char* end = nullptr;
  const unsigned long value = std::strtoul(arg, &end, 10);
  if (errno == ERANGE || end == arg || (end != nullptr && *end != '\0')) {
    throw std::invalid_argument("Seed must be a non-negative integer");
  }

  return value;
}

size_t countFileLines(const std::string& filename) {
  size_t counter = 0;
  std::string line;
  auto file = openInputFile(fs::path(filename));
  while (getline(file, line)) {
    if (!isSkippableLine(line)) counter++;
  }
  return counter;
}

bool fileExists(const std::string& filename) {
  std::error_code error;
  return fs::exists(fs::path(filename), error) && !error;
}

std::vector<std::string> split(std::string s, std::string delimiter) {
  std::vector<std::string> v;
  size_t pos = 0;
  std::string token;
  while ((pos = s.find(delimiter)) != std::string::npos) {
    token = s.substr(0, pos);
    v.push_back(token);
    s.erase(0, pos + delimiter.length());
  }
  v.push_back(s);
  return v;
}

std::vector<double> loadRow(std::string filePath, size_t iRow, std::string delimiter) {
  if (iRow >= countFileLines(filePath)) throw std::runtime_error("row index outside file size");
  std::vector<double> v;
  size_t count = 0;
  std::string line;
  auto file = openInputFile(fs::path(filePath));
  while (getline(file, line)) {
    if (!isSkippableLine(line)) {
      if (iRow == count) {
        auto s = split(line, delimiter);
        v.resize(s.size());
        std::transform(s.begin(), s.end(), v.begin(),
                       [](const std::string& value) { return std::stod(value); });
        break;
      }
      count++;
    }
  }
  return v;
}

std::vector<std::vector<double> > loadFileByRow(std::string filePath, std::string delimiter) {
  std::vector<std::vector<double> > rows;
  std::string line;
  auto file = openInputFile(fs::path(filePath));
  while (getline(file, line)) {
    if (!isSkippableLine(line)) {
      std::vector<double> v;
      auto s = split(line, delimiter);
      v.resize(s.size());
      std::transform(s.begin(), s.end(), v.begin(),
                     [](const std::string& value) { return std::stod(value); });
      rows.push_back(v);
    }
  }
  return rows;
}

OutputFile::OutputFile(const std::string& name) : filename(name) {
  const fs::path output_directory("output");
  ensureDirectoryExists(output_directory);

  const fs::path output_path = output_directory / name;
  filename = output_path.string();
  out.open(output_path);
  if (!out.is_open()) {
    throw std::runtime_error("could not open output file '" + filename + "'");
  }
}

OutputFile::~OutputFile() {
  if (!out.is_open()) return;

  out.flush();
  if (!out) {
    LOGE << "failed writing output file " << filename;
    return;
  }

  out.close();
  if (!out) {
    LOGE << "failed closing output file " << filename;
    return;
  }

  LOGD << "created output file " << filename;
}

}  // namespace utils
}  // namespace gryphon
