/* Vector3.h
 * Adapted from CRPropa3
 * https://crpropa.desy.de/
 */

#ifndef GRYPHON_UTILS_VECTOR3_H
#define GRYPHON_UTILS_VECTOR3_H

#include <cmath>
#include <limits>
#include <ostream>

namespace gryphon {
namespace utils {

template <typename T>
class Vector3 {
 public:
  T x, y, z;

  Vector3() : x(0), y(0), z(0) {}

  // Provides implicit conversion
  template <typename U>
  Vector3(const Vector3<U> &v) : x(v.x), y(v.y), z(v.z) {}

  explicit Vector3(const double *v) : x(v[0]), y(v[1]), z(v[2]) {}

  explicit Vector3(const float *v) : x(v[0]), y(v[1]), z(v[2]) {}

  explicit Vector3(const T &X, const T &Y, const T &Z) : x(X), y(Y), z(Z) {}

  explicit Vector3(T t) : x(t), y(t), z(t) {}

  void setX(const T X) { x = X; }

  void setY(const T Y) { y = Y; }

  void setZ(const T Z) { z = Z; }

  void setXYZ(const T X, const T Y, const T Z) {
    x = X;
    y = Y;
    z = Z;
  }

  T getX() const { return (x); }

  T getY() const { return (y); }

  T getZ() const { return (z); }

  void setRThetaZ(const T r, const T theta, const T z) {
    this->x = r * sin(theta);
    this->y = r * cos(theta);
    this->z = z;
  }

  // magnitude (2-norm) of the vector
  T getModule() const { return (std::sqrt(x * x + y * y + z * z)); }

  T getModuleSquared() const { return (x * x + y * y + z * z); }

  // square of magnitude of the vector in the plane
  T getR() const {
    // LOGD << "check this";
    return (std::sqrt(x * x + y * y));
  }

  T getR2() const {
    // LOGD << "check this";
    return (x * x + y * y);
  }

  // return the azimuth angle
  T getPhi() const {
    T eps = std::numeric_limits<T>::min();
    if ((fabs(x) < eps) && (fabs(y) < eps))
      return (0.0);
    else
      return (std::atan2(y, x));
  }

  // return the unit-vector e_r
  Vector3<T> getUnitVector() const { return (*this / getModule()); }

  // return the unit-vector e_phi
  Vector3<T> getUnitVectorPhi() const { return (Vector3<T>(-sin(getPhi()), cos(getPhi()), 0)); }

  // return the angle [0, pi] between the vectors
  T getAngleTo(const Vector3<T> &v) const {
    T cosdistance = dot(v) / v.getR() / getR();
    // In some directions cosdistance is > 1 on some compilers
    // This ensures that the correct result is returned
    if (cosdistance >= 1.)
      return (0);
    else if (cosdistance <= -1.)
      return (M_PI);
    else
      return (acos(cosdistance));
  }

  // return true if the angle between the vectors is smaller than a threshold
  bool isParallelTo(const Vector3<T> &v, T maxAngle) const { return (getAngleTo(v) < maxAngle); }

  // linear distance to a given vector
  T getDistanceTo(const Vector3<T> &point) const {
    Vector3<T> d = *this - point;
    return (d.getModule());
  }

  // return the component parallel to a second vector
  // 0 if the second vector has 0 magnitude
  Vector3<T> getParallelTo(const Vector3<T> &v) const {
    T vmag = v.getModule();
    if (vmag == std::numeric_limits<T>::min()) return (Vector3<T>(0.));
    return (v * dot(v) / vmag);
  }

  // return the component perpendicular to a second vector
  // 0 if the second vector has 0 magnitude
  Vector3<T> getPerpendicularTo(const Vector3<T> &v) const {
    if (v.getR() == std::numeric_limits<T>::min()) return (Vector3<T>(0.));
    return ((*this) - getParallelTo(v));
  }

  Vector3<T> getRotated2D(T angle) const {
    T r = getR();
    T X = r * std::sin(angle);
    T Y = r * std::cos(angle);
    return (Vector3<T>(X, Y, z));
  }

  // return vector with values limited to the range [lower, upper]
  Vector3<T> clip(T lower, T upper) const {
    Vector3<T> out;
    out.x = std::max(lower, std::min(x, upper));
    out.y = std::max(lower, std::min(y, upper));
    out.z = std::max(lower, std::min(z, upper));
    return (out);
  }

  // return vector with absolute values
  Vector3<T> abs() const { return (Vector3<T>(std::abs(x), std::abs(y), std::abs(z))); }

  // return vector with floored values
  Vector3<T> floor() const { return (Vector3<T>(std::floor(x), std::floor(y), std::floor(z))); }

  // return vector with ceiled values
  Vector3<T> ceil() const { return (Vector3<T>(std::ceil(x), std::ceil(y), std::ceil(z))); }

  // minimum element
  T min() const { return (std::min(x, std::min(y, z))); }

  // maximum element
  T max() const { return (std::max(x, std::max(y, z))); }

  // dot product
  T dot(const Vector3<T> &v) const { return (x * v.x + y * v.y + z * v.z); }

  // cross product
  Vector3<T> cross(const Vector3<T> &v) const {
    return (Vector3<T>(y * v.z - v.y * z, z * v.x - v.z * x, x * v.y - v.x * y));
  }

  // returns true if all elements of the first vector are smaller than those in
  // the second vector
  bool operator<(const Vector3<T> &v) const {
    if (x > v.x)
      return (false);
    else if (x < v.x)
      return (true);
    if (y > v.y)
      return (false);
    else if (y < v.y)
      return (true);
    if (z >= v.z)
      return (false);
    else
      return (true);
  }

  // returns true if all elements of the two vectors are equal
  bool operator==(const Vector3<T> &v) const {
    if (x != v.x) return (false);
    if (y != v.y) return (false);
    if (z != v.z) return (false);
    return (true);
  }

  Vector3<T> operator+(const Vector3<T> &v) const { return (Vector3(x + v.x, y + v.y, z + v.z)); }

  Vector3<T> operator+(const T &f) const { return (Vector3(x + f, y + f, z + f)); }

  Vector3<T> operator-(const Vector3<T> &v) const { return (Vector3(x - v.x, y - v.y, z - v.z)); }

  Vector3<T> operator-(const T &f) const { return (Vector3(x - f, y - f, z - f)); }

  // element-wise multiplication
  Vector3<T> operator*(const Vector3<T> &v) const { return (Vector3(x * v.x, y * v.y, z * v.z)); }

  Vector3<T> operator*(const T &v) const { return (Vector3(x * v, y * v, z * v)); }

  // element-wise division
  Vector3<T> operator/(const Vector3<T> &v) const { return (Vector3(x / v.x, y / v.y, z / v.z)); }

  Vector3<T> operator/(const T &f) const { return (Vector3(x / f, y / f, z / f)); }

  // element-wise modulo operation
  Vector3<T> operator%(const Vector3<T> &v) const {
    return (Vector3(fmod(x, v.x), fmod(y, v.y), fmod(z, v.z)));
  }

  Vector3<T> operator%(const T &f) const { return (Vector3(fmod(x, f), fmod(y, f), fmod(z, f))); }

  Vector3<T> &operator-=(const Vector3<T> &v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return (*this);
  }

  Vector3<T> &operator-=(const T &f) {
    x -= f;
    y -= f;
    z -= f;
    return (*this);
  }

  Vector3<T> &operator+=(const Vector3<T> &v) {
    x += v.x;
    y += v.y;
    z += v.z;
    return (*this);
  }

  Vector3<T> &operator+=(const T &f) {
    x += f;
    y += f;
    z += f;
    return (*this);
  }

  // element-wise multiplication
  Vector3<T> &operator*=(const Vector3<T> &v) {
    x *= v.x;
    y *= v.y;
    z *= v.z;
    return (*this);
  }

  Vector3<T> &operator*=(const T &f) {
    x *= f;
    y *= f;
    z *= f;
    return (*this);
  }

  // element-wise division
  Vector3<T> &operator/=(const Vector3<T> &v) {
    x /= v.x;
    y /= v.y;
    z /= v.z;
    return (*this);
  }

  Vector3<T> &operator/=(const T &f) {
    x /= f;
    y /= f;
    z /= f;
    return (*this);
  }

  // element-wise modulo operation
  Vector3<T> &operator%=(const Vector3<T> &v) {
    x = fmod(x, v.x);
    y = fmod(y, v.y);
    z = fmod(z, v.z);
    return (*this);
  }

  Vector3<T> &operator%=(const T &f) {
    x = fmod(x, f);
    y = fmod(y, f);
    z = fmod(z, f);
    return (*this);
  }

  Vector3<T> &operator=(const Vector3<T> &v) {
    x = v.x;
    y = v.y;
    z = v.z;
    return (*this);
  }

  Vector3<T> &operator=(const T &f) {
    x = f;
    y = f;
    z = f;
    return (*this);
  }
};

template <typename T>
inline std::ostream &operator<<(std::ostream &out, const Vector3<T> &v) {
  out << v.x << " " << v.y << " " << v.z;
  return (out);
}

template <typename T>
inline std::istream &operator>>(std::istream &in, Vector3<T> &v) {
  in >> v.x >> v.y >> v.z;
  return (in);
}

template <typename T>
inline Vector3<T> operator*(T f, const Vector3<T> &v) {
  return (Vector3<T>(v.x * f, v.y * f, v.z * f));
}

using Vector3d = Vector3<double>;
using Vector3f = Vector3<float>;

}  // namespace utils
}  // namespace gryphon

#endif  // GRYPHON_UTILS_VECTOR3_H