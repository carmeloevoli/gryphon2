#ifndef GRYPHON_UTILS_IO_H
#define GRYPHON_UTILS_IO_H

#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "gryphon/core/input.h"

namespace gryphon {
namespace utils {

// String Utilities
std::string trim(const std::string& value);
std::string stripComment(const std::string& line);
std::string normalizeToken(const std::string& value);
std::pair<std::string, std::string> splitKeyValue(const std::string& line);
std::string removeExtensionIniFilename(std::string inputFilename);
double parseDoubleValue(const std::string& filename, size_t line_number, const std::string& key,
                        const std::string& value);
unsigned long parseUnsignedLongValue(const std::string& filename, size_t line_number,
                                     const std::string& key, const std::string& value);
bool parseBoolValue(const std::string& filename, size_t line_number, const std::string& key,
                    const std::string& value);
SpiralModel parseSpiralModelValue(const std::string& filename, size_t line_number,
                                  const std::string& value);
TransportModel parseTransportModelValue(const std::string& filename, size_t line_number,
                                        const std::string& value);
InjectionModel parseInjectionModelValue(const std::string& filename, size_t line_number,
                                        const std::string& value);
core::PID parsePidValue(const std::string& filename, size_t line_number, const std::string& value);
unsigned long parseSeed(const char* arg);

// Files
size_t countFileLines(const std::string& filename);
bool fileExists(const std::string& filename);
std::vector<std::string> split(std::string s, std::string delimiter = " ");
std::vector<double> loadRow(std::string filePath, size_t iRow, std::string delimiter = " ");
std::vector<std::vector<double> > loadFileByRow(std::string filePath, std::string delimiter = " ");

// Output file
class OutputFile {
  std::string filename;
  std::ofstream out;

 public:
  OutputFile(const std::string& name);
  ~OutputFile();

  template <typename T>
  OutputFile& operator<<(const T& value) {
    out << value;
    if (!out) throw std::runtime_error("failed writing output file " + filename);
    return *this;
  }
};

}  // namespace utils
}  // namespace gryphon

#endif  // GRYPHON_UTILS_IO_H
