#ifndef GRYPHON_UTILS_RANDOM_H
#define GRYPHON_UTILS_RANDOM_H

#include <cstdint>
#include <random>
#include <type_traits>

namespace gryphon {
namespace utils {

template <class FloatType = double,
          class = std::enable_if_t<std::is_floating_point<FloatType>::value> >
class RNG {
 public:
  using result_type = FloatType;
  using generator_type = std::mt19937_64;
  using uniform_distribution = std::uniform_real_distribution<result_type>;
  using seed_type = std::uint64_t;

  explicit RNG(seed_type seed) : eng(seed), dist(0.0, 1.0) {}

  // RNG state should be owned explicitly by the caller and passed by reference.
  // Deleting copy prevents accidental duplication of the random stream.
  RNG(const RNG&) = delete;
  RNG& operator=(const RNG&) = delete;
  RNG(RNG&&) noexcept = default;
  RNG& operator=(RNG&&) noexcept = default;

  // generate next random value in the unit interval
  result_type operator()() { return dist(eng); }
  // will always yield 0.0 for this class type
  constexpr result_type min() const { return dist.min(); }
  // will always yield 1.0 for this class type
  constexpr result_type max() const { return dist.max(); }
  // resets internal distribution state without changing the engine seed
  void reset_distribution_state() { dist.reset(); }
  void reseed(seed_type seed) {
    eng.seed(seed);
    dist.reset();
  }
  // uniform distribution
  result_type uniform(result_type vMin, result_type vMax) {
    std::uniform_real_distribution<result_type> uniformDist(vMin, vMax);
    return uniformDist(eng);
  }
  // normal distribution
  result_type normal(result_type mean, result_type stdev) {
    std::normal_distribution<result_type> norm(mean, stdev);
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