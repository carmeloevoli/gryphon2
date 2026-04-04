#include "gryphon/utils/io.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <sys/types.h>

#include "gryphon/utils/logging.h"

namespace gryphon {
namespace utils {

namespace {

bool isSkippableLine(const std::string& line) {
  const auto first = line.find_first_not_of(" \t\r\n");
  return first == std::string::npos || line[first] == '#';
}

std::ifstream openInputFile(const std::string& filename) {
  std::ifstream file(filename.c_str());
  if (!file.is_open()) throw std::runtime_error("could not open file '" + filename + "'");
  return file;
}

void ensureDirectoryExists(const std::string& directory) {
  struct stat info;
  if (stat(directory.c_str(), &info) == 0) {
    if (!S_ISDIR(info.st_mode)) {
      throw std::runtime_error("output path '" + directory + "' exists and is not a directory");
    }
    return;
  }

  if (errno != ENOENT) {
    throw std::runtime_error("could not inspect directory '" + directory + "': " +
                             std::string(std::strerror(errno)));
  }

  if (mkdir(directory.c_str(), 0755) != 0 && errno != EEXIST) {
    throw std::runtime_error("could not create directory '" + directory + "': " +
                             std::string(std::strerror(errno)));
  }
}

}  // namespace

std::string removeExtensionIniFilename(std::string inputFilename) {
  auto ext = inputFilename.substr(inputFilename.rfind('.') + 1);
  if (ext != "ini")
    throw std::invalid_argument("Wrong input filename! It must be : <filename>.ini");
  size_t lastindex = inputFilename.find_last_of(".");
  return inputFilename.substr(0, lastindex);
}

size_t countFileLines(const std::string& filename) {
  size_t counter = 0;
  std::string line;
  auto file = openInputFile(filename);
  while (getline(file, line)) {
    if (!isSkippableLine(line)) counter++;
  }
  return counter;
}

bool fileExists(const std::string& filename) {
  std::ifstream f(filename.c_str());
  return f.good();
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
  auto file = openInputFile(filePath);
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
  auto file = openInputFile(filePath);
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
  ensureDirectoryExists("output");
  out.open("output/" + name);
  if (!out.is_open()) throw std::runtime_error("could not open output file 'output/" + name + "'");
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
