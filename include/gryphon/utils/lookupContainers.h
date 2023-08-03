#ifndef GRYPHON_UTILS_LOOKUPTABLE_H
#define GRYPHON_UTILS_LOOKUPTABLE_H

#include <vector>

#include "gryphon/utils/io.h"
#include "gryphon/utils/numeric.h"
#include "gryphon/utils/progressbar.h"
#include "gryphon/utils/timer.h"

namespace gryphon {
namespace utils {

// typedef std::vector<double>::const_iterator const_iterator;

template <size_t xSize>
class LookupArray {
 public:
  LookupArray() {
    if (xSize < 2) throw std::runtime_error("x-axis size must be > 1");
    m_xAxis.reserve(xSize);
    m_array.reserve(xSize);
  }

  inline double get(double x) const { return utils::interpolate(x, m_xAxis, m_array); }
  inline double getInverse(double x) const { return utils::interpolate(x, m_array, m_xAxis); }
  inline double spline(double x) const { return utils::cspline(x, m_xAxis, m_array); }
  inline bool xIsInside(double x) const { return x >= m_xAxis.front() && x <= m_xAxis.back(); }

 public:
  void loadTable(const std::string& filePath, size_t iCol = 1) {
    if (!utils::fileExists(filePath))
      throw std::runtime_error("file data for lookup array does not exist");
    auto v = utils::loadFileByRow(filePath, ",");
    for (size_t i = 0; i < xSize; ++i) {
      auto line = v.at(i);
      m_xAxis.emplace_back(line[0]);
      m_array.emplace_back(line[iCol]);
    }
    assert(m_xAxis.size() == xSize && m_array.size() == xSize);
  }

  void cacheTable(const std::function<double(double)>& func,
                  const std::pair<double, double>& range) {
    const double dx = (range.second - range.first) / (double)(xSize - 1);
    // Progressbar init
    auto progressbar = std::make_shared<ProgressBar>(xSize);
    auto progressbar_mutex = std::make_shared<std::mutex>();
    progressbar->setMutex(progressbar_mutex);
    progressbar->start("Start caching vector");
    for (size_t i = 0; i < xSize; ++i) {
      progressbar->update();
      auto x = (double)i * dx + range.first;
      auto f_x = func(x);
      m_xAxis.emplace_back(x);
      m_array.emplace_back(f_x);
    }
    assert(m_xAxis.size() == xSize && m_array.size() == xSize);
  }

 protected:
  std::vector<double> m_xAxis;
  std::vector<double> m_array;
};

template <size_t xSize, size_t ySize>
class LookupTable {
 public:
  LookupTable() {
    if (xSize < 2) throw std::runtime_error("x-axis size must be > 1");
    if (ySize < 2) throw std::runtime_error("y-axis size must be > 1");
    m_xAxis.reserve(xSize);
    m_yAxis.reserve(ySize);
    m_table.reserve(xSize * ySize);
  }

  inline double get(double x, double y) const {
    return utils::interpolate2d(x, y, m_xAxis, m_yAxis, m_table);
  }

  bool xIsInside(double x) const { return x >= m_xAxis.front() && x <= m_xAxis.back(); }
  bool yIsInside(double y) const { return y >= m_yAxis.front() && y <= m_yAxis.back(); }

 public:
  //    void loadTable() {
  //     auto v = utils::loadFileByRow(m_filePath, ",");
  //     size_t counter = 0;
  //     for (size_t i = 0; i < xSize; ++i) {
  //       for (size_t j = 0; j < ySize; ++j) {
  //         auto line = v.at(counter);
  //         if (line.size() != 3) throw std::runtime_error("error in reading table values");
  //         if (j == 0) m_xAxis.emplace_back(line[0]);
  //         if (i == 0) m_yAxis.emplace_back(line[1]);
  //         m_table.emplace_back(line[2]);
  //         counter++;
  //       }
  //     }
  //     assert(m_xAxis.size() == xSize && m_yAxis.size() == ySize);
  //   }

  void cacheTable(const std::function<double(double, double)>& func,
                  const std::pair<double, double>& xRange,
                  const std::pair<double, double>& yRange) {
    const double dx = (xRange.second - xRange.first) / (double)(xSize - 1);
    const double dy = (yRange.second - yRange.first) / (double)(ySize - 1);
    // Progressbar init
    auto progressbar = std::make_shared<ProgressBar>(xSize * ySize);
    auto progressbar_mutex = std::make_shared<std::mutex>();
    progressbar->setMutex(progressbar_mutex);
    progressbar->start("Start caching table");

    for (size_t i = 0; i < xSize; ++i) {
      auto x = (double)i * dx + xRange.first;
      m_xAxis.emplace_back(x);
      for (size_t j = 0; j < ySize; ++j) {
        progressbar->update();
        auto y = (double)j * dy + yRange.first;
        if (i == 0) m_yAxis.emplace_back(y);
        auto f_xy = func(x, y);
        m_table.emplace_back(f_xy);
      }
    }
    assert(m_xAxis.size() == xSize && m_yAxis.size() == ySize);
    assert(m_table.size() == xSize * ySize);
  }

 protected:
  std::vector<double> m_xAxis;
  std::vector<double> m_yAxis;
  std::vector<double> m_table;
};

}  // namespace utils
}  // namespace gryphon

#endif