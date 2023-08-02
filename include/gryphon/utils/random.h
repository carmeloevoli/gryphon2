#ifndef GRYPHON_UTILS_RANDOM_H
#define GRYPHON_UTILS_RANDOM_H

#include <random>

namespace gryphon {
namespace utils {

template <class FloatType = double,
          class = std::enable_if_t<std::is_floating_point<FloatType>::value> >
class RNG {
 public:
  typedef FloatType result_type;
  typedef std::mt19937_64 generator_type;
  typedef std::uniform_real_distribution<FloatType> uniform_distribution;

  explicit RNG(const int64_t seed) { eng = generator_type(seed); }

  // generate next random value in distribution
  result_type operator()() { return dist(eng); }
  // will always yield 0.0 for this class type
  constexpr result_type min() const { return dist.min(); }
  // will always yield 1.0 for this class type
  constexpr result_type max() const { return dist.max(); }
  // resets internal state such that next call to operator()
  // does not rely on previous call
  void reset_distribution_state() { dist.reset(); }
  // uniform distribution
  result_type uniform(double vMin, double vMax) { return dist(eng) * (vMax - vMin) + vMin; }
  // normal distibution
  result_type normal(double mean, double stdev) {
    std::normal_distribution<double> norm(mean, stdev);
    return norm(eng);
  }

 private:
  generator_type eng;
  uniform_distribution dist;
};

}  // namespace utils

using RandomNumberGenerator = gryphon::utils::RNG<double>;

}  // namespace gryphon

#endif  // GRYPHON_UTILS_RANDOM_H