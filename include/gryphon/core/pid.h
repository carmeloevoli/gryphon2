#ifndef GRYPHON_CORE_PID_H
#define GRYPHON_CORE_PID_H

#include <cassert>
#include <cmath>
#include <ostream>
#include <string>

namespace gryphon {
namespace core {

class PID {
 public:
  PID() { set(0, 0); }

  PID(const int& Z, const int& A) {
    assert(A > 0);
    assert(Z <= A);
    set(Z, A);
  }

  virtual ~PID() {}

  void set(const int& Z, const int& A) {
    m_Z = Z;
    m_A = A;
    m_id = A * 1000 + Z;
  }

  int get_Z() const { return m_Z; }
  int get_A() const { return m_A; }

  double get_Z_over_A() const {
    return (m_A > 0) ? fabs(static_cast<double>(m_Z) / static_cast<double>(m_A)) : 0;
  }

  int get_id() const { return m_id; }

  bool operator==(const PID& other) const { return m_id == other.m_id; }
  bool operator!=(const PID& other) const { return m_id != other.m_id; }
  bool is_H() const { return (m_Z == 1); }
  bool is_He() const { return (m_Z == 2); }

  friend std::ostream& operator<<(std::ostream& stream, const PID& pid) {
    stream << "(" << pid.get_A() << "," << pid.get_Z() << ")";
    return stream;
  }

  std::string to_string() const {
    std::string ss;
    ss = "(" + std::to_string(m_Z) + "," + std::to_string(m_A) + ")";
    return ss;
  }

 protected:
  int m_Z;
  int m_A;
  int m_id;
};

static const PID H = PID(1, 1);
static const PID He = PID(2, 4);

}  // namespace core
}  // namespace gryphon

#endif  // SIMPROP_PID_H