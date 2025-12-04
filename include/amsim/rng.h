#ifndef AMSIMCPP_RNG_H
#define AMSIMCPP_RNG_H

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <type_traits>
#if __cpp_lib_bitops
#include <bit>
#endif

/// Random number generation utilities
namespace amsim::rng {

/// @brief Return seed or generate one from system clock
/// @param seed Input seed (0 to auto-generate)
/// @return Seed value
inline uint64_t auto_seed(uint64_t seed) {
  if (seed != 0) return seed;
  return static_cast<uint64_t>(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());
}

/// @brief Xoshiro256** pseudorandom number generator
///
/// Fast, high-quality PRNG with 256-bit state and excellent statistical
/// properties.
class Xoshiro256ss {
 public:
  Xoshiro256ss(const Xoshiro256ss&) = delete;
  Xoshiro256ss& operator=(const Xoshiro256ss&) = delete;

  // @brief Get Xoshiro256** instance with default settings
  // @return Reference to singleton Xoshiro256** instance
  static Xoshiro256ss& get_instance() {
    thread_local static Xoshiro256ss rng;
    return rng;
  }

  /// @brief Seed the Xoshiro256** generator
  /// @param seed Seed value
  void set_seed(std::uint64_t seed) noexcept {
    std::uint64_t x = seed ? seed : 0x9e3779b97f4a7c15ULL;
    s_[0] = splitMix64Step(x);
    s_[1] = splitMix64Step(x);
    s_[2] = splitMix64Step(x);
    s_[3] = splitMix64Step(x);
    if ((s_[0] | s_[1] | s_[2] | s_[3]) == 0) s_[0] = 1;
  }

  /// @brief Generate next random value
  /// @return 64-bit random value
  std::uint64_t next() noexcept {
    const std::uint64_t result = rotl(s_[1] * 5, 7) * 9;
    const std::uint64_t t = s_[1] << 17;
    s_[2] ^= s_[0];
    s_[3] ^= s_[1];
    s_[1] ^= s_[2];
    s_[0] ^= s_[3];
    s_[2] ^= t;
    s_[3] = rotl(s_[3], 45);
    return result;
  }

 private:
  explicit Xoshiro256ss() = default;

  std::uint64_t s_[4]{};  ///< Generator state

  /// @brief Rotate left operation
  /// @param x Value to rotate
  /// @param k Number of bits to rotate
  /// @return Rotated value
  static uint64_t rotl(uint64_t x, int k) noexcept {
#if __cpp_lib_bitops
    return std::rotl(x, k);
#else
    return (x << k) | (x >> (64 - k));
#endif
  }

  /// @brief SplitMix64 generator step
  /// @param x State variable (modified in-place)
  /// @return Random value
  static std::uint64_t splitMix64Step(std::uint64_t& x) noexcept {
    std::uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }
};

/// @brief Threshold type for given bit precision
template <int BITS>
using ThrT = std::conditional_t<
    BITS == 8,
    std::uint8_t,
    std::conditional_t<BITS == 16, std::uint16_t, std::uint64_t>>;

/// @brief Generate mask for low k bits
/// @param k Number of bits
/// @return Bit mask
inline std::uint64_t lowbits_mask(unsigned k) noexcept {
  if (k == 0) return 0ULL;
  if (k >= 64) return ~0ULL;
  return (1ULL << k) - 1ULL;
}

/// @brief Convert probability to threshold value
/// @tparam BITS Precision (8, 16, or 64)
/// @param p Probability
/// @return Threshold value
template <int BITS>
inline ThrT<BITS> prob_to_thr(double p) noexcept {
  static_assert(
      BITS == 8 || BITS == 16 || BITS == 64, "BITS must be 8, 16, 64");
  if (p <= 0.0) return ThrT<BITS>(0);
  if (p >= 1.0) return ThrT<BITS>(~ThrT<BITS>(0));
  long double v = std::ldexp(static_cast<long double>(p), BITS);
  auto t = static_cast<std::uint64_t>(v);
  if constexpr (BITS < 64) {
    const std::uint64_t cap = (1ULL << BITS) - 1ULL;
    t = std::min(t, cap);
  }
  return ThrT<BITS>(t);
}

/// @brief Bernoulli 64-bit word generator
///
/// Generates 64-bit words where each bit is an independent Bernoulli trial.
///
/// @tparam BITS Precision for threshold comparisons (8, 16, or 64)
template <int BITS = 16>
struct BernoulliWord {
  static_assert(BITS == 8 || BITS == 16 || BITS == 64);
  using T = ThrT<BITS>;  ///< Threshold type

  /// @brief Construct generator
  BernoulliWord() = default;

  /// @brief Set success probability
  /// @param p Probability
  void set_prob(double p) noexcept {
    p = std::max(0.0, p);
    p = std::min(p, 1.0);
    T thresh = prob_to_thr<BITS>(p);
    for (std::size_t j = 0; j < 64; ++j) tj_[j] = thresh;
  }

  /// @brief Set per-bit probabilities (uses first value)
  /// @param f Array of probabilities
  void set_probs(const double* ptr) noexcept {
    double p;
    for (std::size_t j = 0; j < 64; ++j) {
      p = *(ptr + j);
      p = std::max(0.0, p);
      p = std::min(p, 1.0);
      T thresh = prob_to_thr<BITS>(p);
      tj_[j] = thresh;
    }
  }

  /// @brief Generate single random bit
  /// @return Random boolean
  bool coinflip() noexcept {
    return (Xoshiro256ss::get_instance().next() & 1ULL) != 0U;
  }

  /// @brief Generate 64-bit Bernoulli word
  /// @param valid_bits Number of valid bits to generate
  /// @return Word with random bits
  std::uint64_t sample(unsigned valid_bits = 64) noexcept {
    std::uint64_t w = 0;
    if constexpr (BITS == 8) {
      for (int j = 0; j < 64;) {
        std::uint64_t r = Xoshiro256ss::get_instance().next();
        for (int k = 0; k < 8 && j < 64; ++k, ++j) {
          auto rv = static_cast<std::uint8_t>(r >> 56);
          r <<= 8;
          w |= static_cast<std::uint64_t>(-(rv < tj_[j])) & (1ULL << j);
        }
      }
    } else if constexpr (BITS == 16) {
      for (int j = 0; j < 64;) {
        std::uint64_t r = Xoshiro256ss::get_instance().next();
        for (int k = 0; k < 4 && j < 64; ++k, ++j) {
          auto rv = static_cast<std::uint16_t>(r >> 48);
          r <<= 16;
          w |= static_cast<std::uint64_t>(-(rv < tj_[j])) & (1ULL << j);
        }
      }
    } else {
      for (int j = 0; j < 64; ++j) {
        std::uint64_t rv = Xoshiro256ss::get_instance().next();
        w |= static_cast<std::uint64_t>(-(rv < tj_[j])) & (1ULL << j);
      }
    }
    if (valid_bits < 64) w &= lowbits_mask(valid_bits);
    return w;
  }

 private:
  std::array<T, 64> tj_;  ///< Probability threshold
};

/// @brief Type alias for 16-bit precision Bernoulli word generator
using BW16 = BernoulliWord<16>;

/// @brief Convert 64-bit integer to uniform [0,1) value
/// @param x Random integer
/// @return Uniform value in [0,1) with 53-bit precision
inline double u01_53(const uint64_t x) noexcept {
  return ((x >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

/// @brief Standard normal generator using polar method
///
/// Generates pairs of independent standard normal variates using the
/// Box-Muller polar method.
struct NormalPolar {
  /// @brief Construct generator
  NormalPolar() = default;

  /// @brief Generate pair of standard normal variates
  /// @return Array of two independent N(0,1) values
  static std::array<double, 2> two() noexcept {
    double u;
    double v;
    double s;
    do {
      u = 2.0 * u01_53(Xoshiro256ss::get_instance().next()) - 1.0;
      v = 2.0 * u01_53(Xoshiro256ss::get_instance().next()) - 1.0;
      s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);
    const double m = std::sqrt(-2.0 * std::log(s) / s);
    return {u * m, v * m};
  }

  /// @brief Generate a single standard normal variate
  /// @return A single N(0, 1) variate
  static double single() noexcept { return two()[0]; }

  /// @brief Fill array with standard normal variates
  /// @param out Output array
  /// @param n Number of values to generate
  static void fill(double* out, std::size_t n) noexcept {
    std::size_t i = 0;
    for (; i + 1 < n; i += 2) {
      auto z = two();
      out[i] = z[0];
      out[i + 1] = z[1];
    }
    if (i < n) out[i] = two()[0];
  }

  // @brief Generate a vector of standard normal variates
  // @param n_elem Number of elements to generate
  // @return A double vector comprised of n_elem N(0, 1) variates
  static std::vector<double> sample(std::size_t n_elem) noexcept {
    std::vector<double> out(n_elem);
    fill(out.data(), n_elem);
    return out;
  }
};

/// @brief Uniform [0,a) generator
struct UniformRange {
  /// @brief Construct generator
  UniformRange() = default;

  /// @brief Sample uniform value in [0,a)
  /// @param a Upper bound
  /// @return Uniform value
  static double sample(double a) noexcept {
    return a * u01_53(Xoshiro256ss::get_instance().next());
  }

  /// @brief Fill array with uniform [0,a) values
  /// @param out Output array
  /// @param n Number of values
  /// @param a Upper bound
  static void fill(double* out, std::size_t n, double a) noexcept {
    for (std::size_t i = 0; i < n; ++i)
      out[i] = a * u01_53(Xoshiro256ss::get_instance().next());
  }
};

/// @brief Uniform integer [0,n) generator
struct UniformIntRange {
  /// @brief Construct generator
  UniformIntRange() = default;

  /// @brief Sample uniform integer in [0,n)
  /// @param n Upper bound
  /// @return Uniform integer
  static std::size_t sample(const std::size_t n) {
    if (n == 0) return 0;

    const std::size_t thresh = UINT64_MAX - (UINT64_MAX % n);
    std::size_t x;

    do {
      x = Xoshiro256ss::get_instance().next();
    } while (x >= thresh);

    return x % n;
  }
};

inline void set_seed(std::uint64_t seed) {
  Xoshiro256ss::get_instance().set_seed(seed);
}

}  // namespace amsim::rng

#endif  // AMSIMCPP_RNG_H
