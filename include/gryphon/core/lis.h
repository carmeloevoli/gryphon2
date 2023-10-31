#ifndef GRYPHON_CORE_LIS_H_
#define GRYPHON_CORE_LIS_H_

#include <string>
#include <vector>

#include "gryphon/core/pid.h"

namespace gryphon {
namespace Boschini2017 {

class LIS {
 public:
  LIS(core::PID pid, const std::string& filename);
  virtual ~LIS() = default;
  double get(const double& T_n) const;

 protected:
  void read_filename(const std::string& filename);

 protected:
  int Z;
  int A;
  std::vector<double> logR_table;
  std::vector<double> logdJdR_table;
};

class H_LIS final : public LIS {
 public:
  H_LIS() : LIS(core::H, "data/proton_lis_Boschini2017.txt") {}
};

class He_LIS final : public LIS {
 public:
  He_LIS() : LIS(core::He, "data/helium_lis_Boschini2017.txt") {}
};

}  // namespace Boschini2017
}  // namespace gryphon

#endif